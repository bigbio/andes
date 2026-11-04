//! Same-run claims on an MS1 feature: which other probed target of the run
//! sits on the same feature.
//!
//! A target row (identified or transferred, fitted) is *claimed* by another
//! target row of the same run when the two have the same charge, neutral
//! masses within the run's learned mass window of each other and apexes
//! within its learned apex window ([`Windows`]). The claim is *co-apex* when the apexes are
//! closer than the sum of the two fitted widths: the two probes then measure
//! one MS1 feature, whose quantity is inherently shared between the two
//! compounds (isobaric residue combinations, permutations) or belongs to the
//! claimant alone (a transfer into a run that does not contain the analyte).
//!
//! The claim is a fact about the run's probes, not a verdict: nothing here
//! removes or re-weights a measurement. The best claimant is the identified
//! one when there is any, then the one with the nearest apex.
//!
//! The mass and apex windows are learned per run ([`learn_windows`]) from
//! the agreement of the run's own identified peptides seen at two charge
//! states (one compound measured twice: the difference of their measured
//! mass errors and of their fitted apexes) against the density of chance
//! partners: the window is the contour where the same-compound density falls
//! to the chance density. A run with too few such pairs has no claims.

/// Extent of the pair search for the chance density: relative neutral mass,
/// parts per million, and apex difference, seconds. This is the extent of the
/// search, not a claim window.
pub const SEARCH_PPM: f64 = 20.0;
pub const SEARCH_S: f64 = 60.0;
/// Fewest same-compound pairs a run needs for its windows.
pub const MIN_PAIRS: usize = 30;
/// Quantile of the deviation from the median that equals one standard
/// deviation of a normal law.
const ONE_SD_QUANTILE: f64 = 0.6827;

/// One fitted target row of a run, as the claim search sees it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Probe {
    /// Precursor id (distinguishes the row from its own claimants).
    pub precursor: u32,
    pub charge: u8,
    /// Neutral mass, Da.
    pub mass: f64,
    /// Fitted apex and width, seconds.
    pub t0: f64,
    pub w: f64,
    /// The precursor is identified in this run.
    pub identified: bool,
    /// Id of the peptidoform (the charge states of one peptidoform share it).
    pub peptide: u32,
    /// Measured mass error of the fit, ppm.
    pub mass_err_ppm: f64,
}

/// The claim windows of a run and what they were learned from.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Windows {
    /// Half-width of the neutral-mass window, ppm.
    pub ppm: f64,
    /// Half-width of the apex window, seconds.
    pub seconds: f64,
    /// Robust standard deviation of the mass-error and apex differences of
    /// one compound seen at two charge states.
    pub sd_ppm: f64,
    pub sd_s: f64,
    /// Same-compound pairs behind the standard deviations.
    pub pairs: usize,
}

/// The claims on one row.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Claim {
    /// Best claimant: identified first, then nearest apex; `None` without a
    /// claimant.
    pub claimant: Option<u32>,
    /// The next two claimants in the same order.
    pub others: [Option<u32>; 2],
    /// Apex of the best claimant minus the row's apex, seconds.
    pub dt_s: Option<f32>,
    /// Smallest absolute apex difference over all claimants, seconds.
    pub min_abs_dt_s: Option<f32>,
    /// The best claimant is identified in the run.
    pub identified: bool,
    /// The best claimant is co-apex.
    pub co_apex: bool,
    /// Claimants within the window.
    pub n: u16,
    /// Co-apex claimants.
    pub n_co_apex: u16,
    /// Co-apex claimants identified in the run.
    pub n_identified_co_apex: u16,
}

impl Claim {
    /// The claimants, best first (at most three).
    pub fn claimants(&self) -> impl Iterator<Item = u32> + '_ {
        self.claimant
            .into_iter()
            .chain(self.others.iter().flatten().copied())
    }
}

/// The value below which `q` of the sorted `v` lie, interpolated linearly.
fn quantile(v: &[f64], q: f64) -> f64 {
    let pos = q * (v.len() - 1) as f64;
    let (lo, hi) = (pos.floor() as usize, pos.ceil() as usize);
    v[lo] + (v[hi] - v[lo]) * (pos - lo as f64)
}

/// Robust standard deviation of `x`: the one-sigma quantile of the deviation
/// from the median; when over a third of the values sit exactly on the
/// median, the root mean square of the deviations clipped at `clip`.
fn robust_sd(x: &mut [f64], clip: f64) -> f64 {
    x.sort_by(f64::total_cmp);
    let median = quantile(x, 0.5);
    let mut dev: Vec<f64> = x.iter().map(|v| (v - median).abs()).collect();
    dev.sort_by(f64::total_cmp);
    let v = quantile(&dev, ONE_SD_QUANTILE);
    if v > 0.0 {
        return v;
    }
    let ms = dev.iter().map(|d| d.min(clip).powi(2)).sum::<f64>() / dev.len() as f64;
    ms.sqrt()
}

/// The claim windows of a run from its fitted `probes`; `None` when the run
/// has fewer than [`MIN_PAIRS`] same-compound pairs or a degenerate spread.
pub fn learn_windows(probes: &[Probe]) -> Option<Windows> {
    let fitted = |p: &Probe| p.mass.is_finite() && p.t0.is_finite() && p.w.is_finite();
    // Same compound, two charge states: lower charge minus higher charge.
    let mut by_peptide: std::collections::HashMap<u32, Vec<&Probe>> = Default::default();
    for p in probes
        .iter()
        .filter(|p| p.identified && fitted(p) && p.mass_err_ppm.is_finite())
    {
        by_peptide.entry(p.peptide).or_default().push(p);
    }
    let (mut de, mut dt) = (Vec::new(), Vec::new());
    for group in by_peptide.values() {
        for (i, a) in group.iter().enumerate() {
            for b in &group[i + 1..] {
                let (lo, hi) = if a.charge <= b.charge { (a, b) } else { (b, a) };
                if lo.charge == hi.charge {
                    continue;
                }
                de.push(lo.mass_err_ppm - hi.mass_err_ppm);
                dt.push(lo.t0 - hi.t0);
            }
        }
    }
    let pairs = de.len();
    if pairs < MIN_PAIRS {
        return None;
    }
    let sd_ppm = robust_sd(&mut de, SEARCH_PPM);
    let sd_s = robust_sd(&mut dt, SEARCH_S);
    if !(sd_ppm > 0.0 && sd_s > 0.0) {
        return None;
    }
    // Chance partners of the identified rows: the outer part of the search
    // rectangle, where no same-compound pair reaches.
    let order = sorted_by_mass(probes);
    let (mut anchors, mut outer) = (0usize, 0usize);
    for me in probes.iter().filter(|p| p.identified && fitted(p)) {
        anchors += 1;
        let tol = me.mass * SEARCH_PPM * 1e-6;
        let lo = order
            .partition_point(|&j| (probes[j].charge, probes[j].mass) < (me.charge, me.mass - tol));
        for &j in &order[lo..] {
            let p = &probes[j];
            if p.charge != me.charge || p.mass > me.mass + tol {
                break;
            }
            let dt = (p.t0 - me.t0).abs();
            if p.precursor == me.precursor || dt > SEARCH_S {
                continue;
            }
            let dppm = ((p.mass - me.mass) / me.mass * 1e6).abs();
            if dppm > SEARCH_PPM {
                continue;
            }
            if dppm > SEARCH_PPM / 2.0 || dt > SEARCH_S / 2.0 {
                outer += 1;
            }
        }
    }
    if anchors == 0 {
        return None;
    }
    // Chance partners per anchor per ppm per second over the outer
    // three quarters of the rectangle (+-20 ppm x +-60 s).
    let u = (outer as f64 + 1.0) / (anchors as f64 * 3.0 * SEARCH_PPM * SEARCH_S);
    let f0 = 1.0 / (2.0 * std::f64::consts::PI * sd_ppm * sd_s);
    let l = (2.0 * (f0 / u).ln()).max(0.0).sqrt();
    Some(Windows {
        ppm: (l * sd_ppm).min(SEARCH_PPM),
        seconds: (l * sd_s).min(SEARCH_S),
        sd_ppm,
        sd_s,
        pairs,
    })
}

/// Indices of the fitted probes ordered by (charge, mass).
fn sorted_by_mass(probes: &[Probe]) -> Vec<usize> {
    let mut order: Vec<usize> = (0..probes.len())
        .filter(|&i| {
            let p = &probes[i];
            p.mass.is_finite() && p.t0.is_finite() && p.w.is_finite()
        })
        .collect();
    order.sort_unstable_by(|&a, &b| {
        (probes[a].charge, probes[a].mass)
            .partial_cmp(&(probes[b].charge, probes[b].mass))
            .expect("finite masses")
    });
    order
}

/// The claims on the rows `rows` (indices into `probes`) from all `probes`
/// within `windows`; one [`Claim`] per index of `rows`, in order. Rows or
/// probes without a finite apex, width or mass take no part. Without
/// windows (a run with too few same-compound pairs) no row has a claim.
pub fn assign(probes: &[Probe], rows: &[usize], windows: Option<Windows>) -> Vec<Claim> {
    let Some(win) = windows else {
        return vec![Claim::default(); rows.len()];
    };
    let order = sorted_by_mass(probes);
    rows.iter()
        .map(|&i| {
            let me = &probes[i];
            if !(me.mass.is_finite() && me.t0.is_finite() && me.w.is_finite()) {
                return Claim::default();
            }
            let tol = me.mass * win.ppm * 1e-6;
            let lo = order.partition_point(|&j| {
                (probes[j].charge, probes[j].mass) < (me.charge, me.mass - tol)
            });
            let mut claim = Claim::default();
            // The three best claimants: (identified, apex distance, precursor).
            let mut best: Vec<(bool, f64, u32, usize)> = Vec::with_capacity(4);
            let mut min_dt = f64::INFINITY;
            for &j in &order[lo..] {
                let p = &probes[j];
                if p.charge != me.charge || p.mass > me.mass + tol {
                    break;
                }
                let dt = p.t0 - me.t0;
                if p.precursor == me.precursor || dt.abs() > win.seconds {
                    continue;
                }
                claim.n = claim.n.saturating_add(1);
                min_dt = min_dt.min(dt.abs());
                if dt.abs() <= me.w + p.w {
                    claim.n_co_apex = claim.n_co_apex.saturating_add(1);
                    claim.n_identified_co_apex += u16::from(p.identified);
                }
                best.push((p.identified, dt.abs(), p.precursor, j));
                best.sort_by(|a, b| b.0.cmp(&a.0).then(a.1.total_cmp(&b.1)).then(a.2.cmp(&b.2)));
                best.truncate(3);
            }
            if let Some(&(identified, _, id, at)) = best.first() {
                let p = &probes[at];
                let dt = p.t0 - me.t0;
                claim.claimant = Some(id);
                claim.others = [best.get(1).map(|b| b.2), best.get(2).map(|b| b.2)];
                claim.dt_s = Some(dt as f32);
                claim.min_abs_dt_s = Some(min_dt as f32);
                claim.identified = identified;
                claim.co_apex = dt.abs() <= me.w + p.w;
            }
            claim
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn probe(precursor: u32, mass: f64, t0: f64, identified: bool) -> Probe {
        Probe {
            precursor,
            charge: 2,
            mass,
            t0,
            w: 4.0,
            identified,
            peptide: precursor,
            mass_err_ppm: 0.0,
        }
    }

    fn windows() -> Option<Windows> {
        Some(Windows {
            ppm: 2.0,
            seconds: 60.0,
            sd_ppm: 1.0,
            sd_s: 1.0,
            pairs: 100,
        })
    }

    #[test]
    fn co_apex_identified_claimant_is_best() {
        let probes = [
            probe(0, 1000.0, 100.0, false),
            probe(1, 1000.0 + 1.0e-3, 130.0, true),
            probe(2, 1000.0 + 1.0e-3, 102.0, false),
            probe(3, 1000.0 + 0.1, 100.0, true),
        ];
        let c = assign(&probes, &[0], windows());
        assert_eq!(c[0].claimant, Some(1));
        assert_eq!(c[0].others, [Some(2), None]);
        assert_eq!(c[0].n, 2);
        assert_eq!(c[0].n_co_apex, 1);
        assert_eq!(c[0].n_identified_co_apex, 0);
        assert!(c[0].identified && !c[0].co_apex);
        assert_eq!(c[0].dt_s, Some(30.0));
        assert_eq!(c[0].min_abs_dt_s, Some(2.0));
        assert_eq!(c[0].claimants().collect::<Vec<_>>(), [1, 2]);
    }

    #[test]
    fn nearest_apex_among_equals_and_window() {
        let probes = [
            probe(0, 1000.0, 100.0, false),
            probe(1, 1000.0, 105.0, false),
            probe(2, 1000.0, 99.0, false),
            probe(3, 1000.0, 400.0, true),
        ];
        let c = assign(&probes, &[0], windows());
        assert_eq!(c[0].claimant, Some(2));
        assert!(c[0].co_apex && !c[0].identified);
        assert_eq!(c[0].n, 2);
        assert_eq!(c[0].n_co_apex, 2);
    }

    #[test]
    fn the_windows_set_the_claimants() {
        let probes = [
            probe(0, 1000.0, 100.0, false),
            probe(1, 1000.0 * (1.0 + 3.0e-6), 100.0, true),
            probe(2, 1000.0, 108.0, true),
        ];
        let narrow = Windows {
            ppm: 1.0,
            seconds: 5.0,
            ..windows().unwrap()
        };
        assert_eq!(assign(&probes, &[0], Some(narrow))[0], Claim::default());
        let wide = Windows {
            ppm: 4.0,
            seconds: 10.0,
            ..narrow
        };
        assert_eq!(assign(&probes, &[0], Some(wide))[0].n, 2);
        assert_eq!(assign(&probes, &[0], None)[0], Claim::default());
    }

    #[test]
    fn charge_mass_and_self_are_excluded() {
        let mut other_charge = probe(1, 1000.0, 100.0, true);
        other_charge.charge = 3;
        let probes = [
            probe(0, 1000.0, 100.0, false),
            other_charge,
            probe(2, 1000.0 * (1.0 + 3.0e-6), 100.0, true),
        ];
        let c = assign(&probes, &[0], windows());
        assert_eq!(c[0], Claim::default());
    }

    #[test]
    fn unfitted_probes_take_no_part() {
        let mut nan = probe(1, 1000.0, f64::NAN, true);
        nan.w = f64::NAN;
        let probes = [probe(0, 1000.0, 100.0, false), nan];
        assert_eq!(
            assign(&probes, &[0, 1], windows()),
            vec![Claim::default(); 2]
        );
    }

    /// Peptides seen at two charge states with a small spread of the mass
    /// error and apex, and unrelated partners spread over the run.
    fn run_probes(n: usize) -> Vec<Probe> {
        let mut out = Vec::new();
        let mut state = 12345u64;
        let mut next = move || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        for k in 0..n {
            let t = 100.0 + 5000.0 * next();
            let m = 800.0 + 1500.0 * next();
            let (e1, e2) = (next() - 0.5, next() - 0.5);
            let (d1, d2) = (next() - 0.5, next() - 0.5);
            for (c, e, d, z) in [(0u32, e1, d1, 2u8), (1, e2, d2, 3)] {
                out.push(Probe {
                    precursor: (2 * k) as u32 + c,
                    charge: z,
                    mass: m * f64::from(z) / 2.0 + m * 1e-6 * e,
                    t0: t + d,
                    w: 4.0,
                    identified: true,
                    peptide: k as u32,
                    mass_err_ppm: e,
                });
            }
        }
        out
    }

    #[test]
    fn windows_follow_the_same_compound_spread() {
        let w = learn_windows(&run_probes(400)).expect("enough pairs");
        // Uniform differences of width 1 have a standard deviation near 0.41
        // (triangular law, quantile 0.6827 of |d|).
        assert!(w.sd_ppm > 0.2 && w.sd_ppm < 0.6, "{w:?}");
        assert!(w.sd_s > 0.2 && w.sd_s < 0.6, "{w:?}");
        assert!(w.ppm > 3.0 * w.sd_ppm && w.ppm <= SEARCH_PPM, "{w:?}");
        assert!(w.seconds > 3.0 * w.sd_s && w.seconds <= SEARCH_S, "{w:?}");
        assert_eq!(w.pairs, 400);
        // Twice the mass-error spread widens the mass window.
        let mut loose = run_probes(400);
        for p in &mut loose {
            p.mass_err_ppm *= 2.0;
        }
        let w2 = learn_windows(&loose).expect("enough pairs");
        assert!(w2.sd_ppm > 1.5 * w.sd_ppm && w2.ppm > w.ppm, "{w2:?} {w:?}");
    }

    #[test]
    fn a_run_without_pairs_has_no_windows() {
        let probes: Vec<Probe> = (0..40)
            .map(|i| probe(i, 1000.0 + f64::from(i), 100.0, true))
            .collect();
        assert_eq!(learn_windows(&probes), None);
    }

    #[test]
    fn exact_agreement_gives_a_root_mean_square_spread() {
        let mut x = vec![0.0; 80];
        x.extend([3.0, -3.0, 6.0]);
        // 68 % quantile of the deviation is 0: the clipped root mean square.
        let sd = robust_sd(&mut x, 60.0);
        assert!((sd - ((9.0 + 9.0 + 36.0) / 83.0f64).sqrt()).abs() < 1e-12);
    }
}
