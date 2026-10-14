//! Run normalisation of label-free MS1 quantities across runs.
//!
//! Each run gets a log2 offset `o_r(t)` at apex RT `t`, fitted on reference
//! cells (one log2 quantity `l` per precursor and run, typically the
//! identified features of high quality). A precursor enters the fit when it
//! has a reference cell in at least [`NormParams::min_runs`] runs (capped at
//! the number of runs, and at least 2), counting only runs with at least
//! [`NormParams::min_cells`] such cells.
//!
//! The offsets start at 0 and are refined for [`NormParams::rounds`] rounds:
//!
//! * the precursor's reference level is `ref_p = mean over its runs of
//!   (l - o)`, and each of its cells gets `z = l - ref_p`;
//! * a reference cell's offset becomes the median of the `z` of the
//!   [`NormParams::half_window`] reference cells of its run just before it
//!   and as many just after it in apex RT order, its own `z` left out.
//!
//! Any (precursor, run, RT) then gets the median of the last round's `z` of
//! the `half_window` reference cells of the run before and after that RT
//! (the reference cell itself left out when the precursor is one).
//!
//! The offset carries the run's level and its changes of ion signal along the
//! gradient, which a single factor per run cannot follow (electrospray drifts
//! of seconds to minutes). Measuring each run against the precursor's mean
//! over all its runs, own run included, keeps a drift of one run out of the
//! other runs' offsets. Like every normalisation, it assumes most precursors
//! do not change between runs: the median sets the level.
//!
//! Runs below `min_cells`, and every run of a set with fewer than two such
//! runs, get no offset.

use std::collections::HashMap;

/// One reference cell: precursor `precursor` measured in run `run`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NormCell {
    pub precursor: u32,
    pub run: u16,
    /// log2 of the quantity; must be finite.
    pub log2_quantity: f64,
    /// Apex RT in the run's own time, seconds; must be finite.
    pub rt: f64,
}

/// Parameters of [`fit`].
#[derive(Debug, Clone, PartialEq)]
pub struct NormParams {
    /// Fewest runs with a reference cell for a precursor to enter the fit
    /// (capped at the number of runs, and at least 2).
    pub min_runs: usize,
    /// Reference cells on each side of `t` in the RT-local median.
    pub half_window: usize,
    /// Fewest reference cells for a run to be normalised.
    pub min_cells: usize,
    /// Refinement rounds.
    pub rounds: usize,
}

impl Default for NormParams {
    fn default() -> Self {
        Self {
            min_runs: 3,
            half_window: 50,
            min_cells: 20,
            rounds: 5,
        }
    }
}

/// The reference cells of one run, in apex RT order, with their last `z`.
#[derive(Debug, Clone, PartialEq)]
struct RunCurve {
    rt: Vec<f64>,
    z: Vec<f64>,
    /// Position in `rt` of each reference precursor of the run.
    position: HashMap<u32, usize>,
}

impl RunCurve {
    /// Median of the `z` of the `h` cells before and after `rt`, without the
    /// reference cell of `precursor` when it is one.
    fn at(&self, precursor: u32, rt: f64, h: usize) -> Option<f64> {
        let n = self.z.len();
        let mut w: Vec<f64> = match self.position.get(&precursor) {
            Some(&k) => self.z[k.saturating_sub(h)..k]
                .iter()
                .chain(&self.z[k + 1..(k + 1 + h).min(n)])
                .copied()
                .collect(),
            None => {
                let p = self.rt.partition_point(|&x| x < rt);
                self.z[p.saturating_sub(h)..(p + h).min(n)].to_vec()
            }
        };
        median(&mut w)
    }
}

/// Fitted offsets of every run (see the module docs).
#[derive(Debug, Clone, PartialEq)]
pub struct RunNormalization {
    curves: Vec<Option<RunCurve>>,
    half_window: usize,
    /// Precursors that entered the fit.
    pub precursors: usize,
}

/// Median of `v` (mean of the two middle values for an even length); `None`
/// when empty. Reorders `v`.
fn median(v: &mut [f64]) -> Option<f64> {
    let n = v.len();
    if n == 0 {
        return None;
    }
    let cmp = |a: &f64, b: &f64| a.total_cmp(b);
    let (_, &mut hi, _) = v.select_nth_unstable_by(n / 2, cmp);
    if n % 2 == 1 {
        return Some(hi);
    }
    let lo = v[..n / 2].iter().copied().max_by(cmp).unwrap_or(hi);
    Some(0.5 * (lo + hi))
}

/// Fits the offsets of runs `0..n_runs` on `cells` (see the module docs).
/// Cells with a non-finite value, a run out of range or a repeated
/// (precursor, run) after the first are ignored.
pub fn fit(cells: &[NormCell], n_runs: usize, params: &NormParams) -> RunNormalization {
    let none = |precursors| RunNormalization {
        curves: vec![None; n_runs],
        half_window: params.half_window,
        precursors,
    };
    let min_runs = params.min_runs.min(n_runs).max(2);
    let mut by_precursor: HashMap<u32, Vec<NormCell>> = HashMap::new();
    for c in cells {
        if (c.run as usize) < n_runs && c.log2_quantity.is_finite() && c.rt.is_finite() {
            let v = by_precursor.entry(c.precursor).or_default();
            if !v.iter().any(|o| o.run == c.run) {
                v.push(*c);
            }
        }
    }
    let eligible = |v: &Vec<NormCell>| v.len() >= min_runs;
    let mut per_run = vec![0usize; n_runs];
    for c in by_precursor.values().filter(|v| eligible(v)).flatten() {
        per_run[c.run as usize] += 1;
    }
    let active: Vec<bool> = per_run.iter().map(|&n| n >= params.min_cells).collect();
    if n_runs < 2 || active.iter().filter(|&&a| a).count() < 2 {
        return none(0);
    }
    // Precursor groups over the active runs, in precursor order.
    let mut groups: Vec<Vec<NormCell>> = by_precursor
        .into_values()
        .map(|v| v.into_iter().filter(|c| active[c.run as usize]).collect())
        .filter(eligible)
        .collect();
    groups.sort_by_key(|v| v[0].precursor);

    // Per run, the (group, cell) indices of its reference cells in RT order.
    let mut order: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n_runs];
    for (gi, v) in groups.iter().enumerate() {
        for (ci, c) in v.iter().enumerate() {
            order[c.run as usize].push((gi, ci));
        }
    }
    for o in &mut order {
        o.sort_by(|&(ga, ca), &(gb, cb)| {
            let (a, b) = (&groups[ga][ca], &groups[gb][cb]);
            a.rt.total_cmp(&b.rt).then(a.precursor.cmp(&b.precursor))
        });
    }
    let mut offset: Vec<Vec<f64>> = groups.iter().map(|v| vec![0.0; v.len()]).collect();
    let mut z: Vec<Vec<f64>> = offset.clone();
    let h = params.half_window;
    for _ in 0..params.rounds.max(1) {
        for (gi, v) in groups.iter().enumerate() {
            let level = v
                .iter()
                .zip(&offset[gi])
                .map(|(c, o)| c.log2_quantity - o)
                .sum::<f64>()
                / v.len() as f64;
            for (ci, c) in v.iter().enumerate() {
                z[gi][ci] = c.log2_quantity - level;
            }
        }
        for o in &order {
            let zs: Vec<f64> = o.iter().map(|&(g, c)| z[g][c]).collect();
            for (k, &(g, c)) in o.iter().enumerate() {
                let mut w: Vec<f64> = zs[k.saturating_sub(h)..k]
                    .iter()
                    .chain(&zs[k + 1..(k + 1 + h).min(zs.len())])
                    .copied()
                    .collect();
                offset[g][c] = median(&mut w).unwrap_or(0.0);
            }
        }
    }
    let curves = order
        .iter()
        .enumerate()
        .map(|(run, o)| {
            active[run].then(|| RunCurve {
                rt: o.iter().map(|&(g, c)| groups[g][c].rt).collect(),
                z: o.iter().map(|&(g, c)| z[g][c]).collect(),
                position: o
                    .iter()
                    .enumerate()
                    .map(|(k, &(g, c))| (groups[g][c].precursor, k))
                    .collect(),
            })
        })
        .collect();
    RunNormalization {
        curves,
        half_window: h,
        precursors: groups.len(),
    }
}

impl RunNormalization {
    /// The log2 offset of run `run` for precursor `precursor` at apex RT `rt`;
    /// `None` for a run that is not normalised or a non-finite `rt`.
    pub fn offset(&self, run: u16, precursor: u32, rt: f64) -> Option<f64> {
        let c = self.curves.get(run as usize)?.as_ref()?;
        if !rt.is_finite() {
            return None;
        }
        c.at(precursor, rt, self.half_window)
    }

    /// `quantity * 2^-offset`; `None` when the run is not normalised.
    pub fn normalize(&self, run: u16, precursor: u32, rt: f64, quantity: f64) -> Option<f64> {
        self.offset(run, precursor, rt)
            .map(|o| quantity * (-o).exp2())
    }

    /// Whether run `run` is normalised.
    pub fn is_normalized(&self, run: u16) -> bool {
        matches!(self.curves.get(run as usize), Some(Some(_)))
    }

    /// Median and standard deviation of run `run`'s offsets over its
    /// reference cells, and their number; `None` for a run that is not
    /// normalised.
    pub fn summary(&self, run: u16) -> Option<(f64, f64, usize)> {
        let c = self.curves.get(run as usize)?.as_ref()?;
        let mut v: Vec<f64> = c
            .position
            .iter()
            .filter_map(|(&p, &k)| c.at(p, c.rt[k], self.half_window))
            .collect();
        let n = v.len();
        let mean = v.iter().sum::<f64>() / n as f64;
        let sd = (v.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n as f64).sqrt();
        Some((median(&mut v)?, sd, n))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Deterministic uniform numbers in [0, 1).
    struct Lcg(u64);
    impl Lcg {
        fn unit(&mut self) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (self.0 >> 11) as f64 / (1u64 << 53) as f64
        }
        fn normal(&mut self) -> f64 {
            let (u, v) = (self.unit().max(1e-12), self.unit());
            (-2.0 * u.ln()).sqrt() * (2.0 * std::f64::consts::PI * v).cos()
        }
    }

    /// `n` precursors eluting over 0..3000 s in `runs` runs, abundance
    /// 20 + 5 * N(0,1) log2, plus `gain(run, rt)` and `noise` * N(0,1).
    fn cells(
        n: u32,
        runs: u16,
        noise: f64,
        gain: impl Fn(u16, f64) -> f64,
    ) -> (Vec<NormCell>, Vec<f64>) {
        let mut g = Lcg(42);
        let mut out = Vec::new();
        let mut truth = Vec::new();
        for p in 0..n {
            let a = 20.0 + 5.0 * g.normal();
            let rt = 3000.0 * g.unit();
            truth.push(a);
            for run in 0..runs {
                out.push(NormCell {
                    precursor: p,
                    run,
                    log2_quantity: a + gain(run, rt) + noise * g.normal(),
                    rt,
                });
            }
        }
        (out, truth)
    }

    fn corrected(norm: &RunNormalization, c: &NormCell) -> f64 {
        c.log2_quantity - norm.offset(c.run, c.precursor, c.rt).unwrap()
    }

    #[test]
    fn constant_run_factors_are_removed() {
        let shift = [0.0, 0.6, -0.4];
        let (cells, _) = cells(2000, 3, 0.05, |r, _| shift[r as usize]);
        let norm = fit(&cells, 3, &NormParams::default());
        let mean = shift.iter().sum::<f64>() / 3.0;
        for r in 0..3u16 {
            let o = norm.offset(r, u32::MAX, 1500.0).unwrap();
            assert!(
                (o - (shift[r as usize] - mean)).abs() < 0.02,
                "run {r}: {o}"
            );
        }
        // Each precursor reads the same in every run after normalisation.
        for c in cells.iter().filter(|c| c.precursor < 50) {
            let first = corrected(&norm, &cells[3 * c.precursor as usize]);
            assert!((corrected(&norm, c) - first).abs() < 0.3);
        }
    }

    #[test]
    fn rt_dependent_gain_is_followed() {
        // Run 1 loses up to 0.8 log2 in the middle of the gradient.
        let gain = |r: u16, t: f64| {
            if r == 1 {
                -0.8 * (std::f64::consts::PI * t / 3000.0).sin()
            } else {
                0.0
            }
        };
        let (cells, _) = cells(6000, 3, 0.1, gain);
        let norm = fit(&cells, 3, &NormParams::default());
        // Per precursor, the run-1 minus run-0 offset is the gain at its RT.
        let mut worst: f64 = 0.0;
        let mut total = 0.0;
        for c in cells.iter().filter(|c| c.run == 1) {
            let o1 = norm.offset(1, c.precursor, c.rt).unwrap();
            let o0 = norm.offset(0, c.precursor, c.rt).unwrap();
            let e = (o1 - o0 - gain(1, c.rt)).abs();
            worst = worst.max(e);
            total += e;
        }
        assert!(total / 6000.0 < 0.03, "mean error {}", total / 6000.0);
        assert!(worst < 0.12, "max error {worst}");
        let mid =
            norm.offset(1, u32::MAX, 1500.0).unwrap() - norm.offset(0, u32::MAX, 1500.0).unwrap();
        assert!((mid + 0.8).abs() < 0.05, "{mid}");
    }

    #[test]
    fn a_minority_of_changed_precursors_does_not_move_the_level() {
        // A fifth of the precursors are 4x up in run 2: the median ignores them.
        let (mut cells, _) = cells(3000, 3, 0.05, |_, _| 0.0);
        for c in cells
            .iter_mut()
            .filter(|c| c.run == 2 && c.precursor % 5 == 0)
        {
            c.log2_quantity += 2.0;
        }
        let norm = fit(&cells, 3, &NormParams::default());
        // A mean would put run 2 0.4 above runs 0 and 1 (0.2 * 2); the
        // median keeps the runs within a small fraction of that.
        let o: Vec<f64> = (0..3u16)
            .map(|r| norm.offset(r, u32::MAX, 1500.0).unwrap())
            .collect();
        for r in 1..3 {
            assert!((o[r] - o[0]).abs() < 0.08, "{o:?}");
        }
    }

    #[test]
    fn a_reference_cell_is_left_out_of_its_own_offset() {
        // Five precursors in two runs, run 1 reads +1 except precursor 2 (+5).
        let mut cells = Vec::new();
        for p in 0..5u32 {
            for run in 0..2u16 {
                let up = if run == 1 {
                    if p == 2 {
                        5.0
                    } else {
                        1.0
                    }
                } else {
                    0.0
                };
                cells.push(NormCell {
                    precursor: p,
                    run,
                    log2_quantity: 10.0 + up,
                    rt: 100.0 * p as f64,
                });
            }
        }
        let params = NormParams {
            min_cells: 3,
            ..NormParams::default()
        };
        let norm = fit(&cells, 2, &params);
        // The window covers the whole run: precursor 2's own +5 never enters
        // its own offset, so it keeps its excess of 4 over the others.
        let d = |p: u32| {
            let c1 = cells[2 * p as usize + 1];
            let c0 = cells[2 * p as usize];
            corrected(&norm, &c1) - corrected(&norm, &c0)
        };
        assert!((d(2) - 4.0).abs() < 1e-9, "{}", d(2));
        assert!(d(0).abs() < 1e-9, "{}", d(0));
    }

    #[test]
    fn runs_without_enough_cells_are_not_normalised() {
        let (mut cells, _) = cells(100, 3, 0.05, |_, _| 0.0);
        // Run 2 keeps 10 cells: fewer than min_cells.
        cells.retain(|c| c.run != 2 || c.precursor < 10);
        let params = NormParams {
            min_runs: 2,
            ..NormParams::default()
        };
        let norm = fit(&cells, 3, &params);
        assert!(norm.is_normalized(0) && norm.is_normalized(1));
        assert!(!norm.is_normalized(2));
        assert_eq!(norm.normalize(2, 0, 10.0, 1e6), None);
        assert_eq!(norm.precursors, 100);
        // With min_runs 3 only 10 precursors qualify: no run has 20 cells.
        let norm = fit(&cells, 3, &NormParams::default());
        assert!((0..3).all(|r| !norm.is_normalized(r)));
        // A single run is never normalised.
        let one: Vec<NormCell> = cells.iter().filter(|c| c.run == 0).copied().collect();
        let norm = fit(&one, 1, &NormParams::default());
        assert!(!norm.is_normalized(0));
        assert_eq!(norm.normalize(0, 0, 10.0, 1e6), None);
    }

    #[test]
    fn two_runs_use_precursors_seen_in_both() {
        let (cells, _) = cells(500, 2, 0.05, |r, _| if r == 1 { 1.0 } else { 0.0 });
        let norm = fit(&cells, 2, &NormParams::default());
        assert_eq!(norm.precursors, 500);
        let d =
            norm.offset(1, u32::MAX, 1500.0).unwrap() - norm.offset(0, u32::MAX, 1500.0).unwrap();
        assert!((d - 1.0).abs() < 0.02, "{d}");
    }

    #[test]
    fn invalid_and_repeated_cells_are_ignored() {
        let (mut cells, _) = cells(200, 3, 0.05, |_, _| 0.0);
        let clean = fit(&cells, 3, &NormParams::default());
        cells.push(NormCell {
            precursor: 0,
            run: 1,
            log2_quantity: 99.0,
            rt: 5.0,
        });
        cells.push(NormCell {
            precursor: 1,
            run: 7,
            log2_quantity: 1.0,
            rt: 5.0,
        });
        cells.push(NormCell {
            precursor: 1,
            run: 0,
            log2_quantity: f64::NAN,
            rt: 5.0,
        });
        assert_eq!(fit(&cells, 3, &NormParams::default()), clean);
    }

    #[test]
    fn median_of_odd_and_even_lengths() {
        assert_eq!(median(&mut []), None);
        assert_eq!(median(&mut [3.0, 1.0, 2.0]), Some(2.0));
        assert_eq!(median(&mut [4.0, 1.0, 3.0, 2.0]), Some(2.5));
    }
}
