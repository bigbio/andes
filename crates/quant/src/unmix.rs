use model::mass::ISOTOPE;
use model::stats::{median_mad, Even, MAD_TO_SD};

use crate::nnls::nnls;

pub const MAX_INTERFERENCE: f64 = 0.8;
pub const MIN_FIT_SCANS: usize = 100;
const ENVELOPE_PPM: f64 = 10.0;
const ENVELOPE_ISOTOPES: std::ops::RangeInclusive<i32> = -1..=8;
const DENSITY_HALF_WIDTH: f64 = 5.0;
const BISQUARE_K: f64 = 4.685;
const ROBUST_ROUNDS: usize = 8;
const N_COEF: usize = 5;
const MIN_MAD: f64 = 1e-6;
const SLOPE_ROUNDS: usize = 8;

#[derive(Debug, Clone, PartialEq)]
pub struct ScanInput {
    pub charge: u8,
    pub reporters: Vec<f32>,
    pub pic: f64,
    pub window: Option<WindowRecord>,
    pub protein: Option<u32>,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WindowRecord {
    pub total: f64,
    pub precursor: f64,
    pub hidden: f64,
}

impl WindowRecord {
    pub fn precursor_fraction(&self) -> f64 {
        if self.total > 0.0 {
            (self.precursor / self.total).clamp(0.0, 1.0)
        } else {
            1.0
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Survey<'a> {
    pub rt: f64,
    pub peaks: &'a [(f64, f32)],
    pub detection_limit: f32,
}

fn window_sums(peaks: &[(f64, f32)], lo: f64, hi: f64, mono_mz: f64, charge: u8) -> (f64, f64) {
    let start = peaks.partition_point(|p| p.0 < lo);
    let end = start + peaks[start..].partition_point(|p| p.0 <= hi);
    let window = &peaks[start..end];
    let total: f64 = window.iter().map(|p| p.1 as f64).sum();
    let z = charge.max(1) as f64;
    let mut envelope = 0.0;
    for i in ENVELOPE_ISOTOPES {
        let expected = mono_mz + i as f64 * ISOTOPE / z;
        if expected < lo || expected > hi {
            continue;
        }
        let tol = expected * ENVELOPE_PPM * 1e-6;
        let a = window.partition_point(|p| p.0 < expected - tol);
        let best = window[a..]
            .iter()
            .take_while(|p| p.0 <= expected + tol)
            .map(|p| p.1 as f64)
            .fold(0.0, f64::max);
        envelope += best;
    }
    (total, envelope.min(total))
}

fn centroid_density(peaks: &[(f64, f32)], centre: f64) -> f64 {
    let lo = peaks.partition_point(|p| p.0 < centre - DENSITY_HALF_WIDTH);
    let hi = peaks.partition_point(|p| p.0 <= centre + DENSITY_HALF_WIDTH);
    (hi - lo) as f64 / (2.0 * DENSITY_HALF_WIDTH)
}

pub fn window_record(
    prev: Option<Survey<'_>>,
    next: Option<Survey<'_>>,
    rt: f64,
    lo: f64,
    hi: f64,
    mono_mz: f64,
    charge: u8,
) -> Option<WindowRecord> {
    if hi.partial_cmp(&lo) != Some(std::cmp::Ordering::Greater) {
        return None;
    }
    let weight_prev = match (prev, next) {
        (Some(p), Some(n)) if n.rt > p.rt => ((n.rt - rt) / (n.rt - p.rt)).clamp(0.0, 1.0),
        (Some(_), Some(_)) => 0.5,
        (Some(_), None) => 1.0,
        (None, Some(_)) => 0.0,
        (None, None) => return None,
    };
    let centre = 0.5 * (lo + hi);
    let mut rec = WindowRecord {
        total: 0.0,
        precursor: 0.0,
        hidden: 0.0,
    };
    for (scan, w) in [(prev, weight_prev), (next, 1.0 - weight_prev)] {
        let Some(s) = scan else { continue };
        if w <= 0.0 {
            continue;
        }
        let (total, precursor) = window_sums(s.peaks, lo, hi, mono_mz, charge);
        rec.total += w * total;
        rec.precursor += w * precursor;
        rec.hidden += w * s.detection_limit as f64 * centroid_density(s.peaks, centre) * (hi - lo);
    }
    Some(rec)
}

#[derive(Debug, Clone, PartialEq)]
pub struct InterferenceModel {
    pub precursor: [f64; 3],
    pub co_isolated: f64,
    pub hidden: f64,
    pub n: usize,
    pub r2_log: f64,
}

fn charge_class(z: u8) -> usize {
    match z {
        0..=2 => 0,
        3 => 1,
        _ => 2,
    }
}

fn design_row(s: &ScanInput) -> Option<[f64; N_COEF]> {
    let w = s.window?;
    if !(s.pic.is_finite() && s.pic > 0.0) {
        return None;
    }
    let f = w.precursor_fraction();
    let mut row = [0.0; N_COEF];
    row[charge_class(s.charge)] = s.pic * f;
    row[3] = s.pic * (1.0 - f);
    row[4] = if w.hidden.is_finite() { w.hidden } else { 0.0 };
    Some(row)
}

fn summed(r: &[f32]) -> f64 {
    r.iter().map(|&v| (v as f64).max(0.0)).sum()
}

impl InterferenceModel {
    pub fn fit(scans: &[ScanInput]) -> Option<Self> {
        let (rows, y): (Vec<[f64; N_COEF]>, Vec<f64>) = scans
            .iter()
            .filter_map(|s| {
                let total = summed(&s.reporters);
                (total > 0.0).then_some(())?;
                Some((design_row(s)?, total))
            })
            .unzip();
        let m = rows.len();
        if m < MIN_FIT_SCANS {
            return None;
        }
        let mut weight = vec![1.0f64; m];
        let mut beta = vec![0.0f64; N_COEF];
        let mut a = vec![0.0f64; m * N_COEF];
        let mut b = vec![0.0f64; m];
        for _ in 0..ROBUST_ROUNDS {
            for i in 0..m {
                let s = weight[i] / y[i];
                for j in 0..N_COEF {
                    a[i * N_COEF + j] = rows[i][j] * s;
                }
                b[i] = weight[i];
            }
            beta = nnls(&a, m, N_COEF, &b, 50 * N_COEF);
            let resid: Vec<f64> = (0..m)
                .map(|i| (y[i] - dot(&rows[i], &beta)) / y[i])
                .collect();
            let (_, raw_mad) = median_mad(&resid, Even::Upper).unwrap_or((0.0, 0.0));
            let mad = (MAD_TO_SD * raw_mad).max(MIN_MAD);
            for (w, r) in weight.iter_mut().zip(&resid) {
                let u = r / (BISQUARE_K * mad);
                *w = if u.abs() < 1.0 { 1.0 - u * u } else { 0.0 };
            }
        }
        if beta.iter().all(|&v| v <= 0.0) {
            return None;
        }
        let logs: Vec<(f64, f64)> = rows
            .iter()
            .zip(&y)
            .map(|(r, &t)| (t.ln(), dot(r, &beta).max(f64::MIN_POSITIVE).ln()))
            .collect();
        let mean = logs.iter().map(|l| l.0).sum::<f64>() / m as f64;
        let ss_tot: f64 = logs.iter().map(|l| (l.0 - mean).powi(2)).sum();
        let ss_res: f64 = logs.iter().map(|l| (l.0 - l.1).powi(2)).sum();
        Some(Self {
            precursor: [beta[0], beta[1], beta[2]],
            co_isolated: beta[3],
            hidden: beta[4],
            n: m,
            r2_log: if ss_tot > 0.0 {
                1.0 - ss_res / ss_tot
            } else {
                0.0
            },
        })
    }

    pub fn interference(&self, scan: &ScanInput) -> Option<f64> {
        let row = design_row(scan)?;
        let beta = self.coefficients();
        let predicted = dot(&row, &beta);
        if predicted <= 0.0 {
            return None;
        }
        let own = row[charge_class(scan.charge)] * beta[charge_class(scan.charge)];
        Some((1.0 - own / predicted).clamp(0.0, MAX_INTERFERENCE))
    }

    fn coefficients(&self) -> [f64; N_COEF] {
        [
            self.precursor[0],
            self.precursor[1],
            self.precursor[2],
            self.co_isolated,
            self.hidden,
        ]
    }
}

fn dot(row: &[f64; N_COEF], beta: &[f64]) -> f64 {
    row.iter().zip(beta).map(|(a, b)| a * b).sum()
}

pub fn background_profile(scans: &[ScanInput]) -> Option<Vec<f64>> {
    let n = scans.first()?.reporters.len();
    let mut sum = vec![0.0f64; n];
    for s in scans.iter().filter(|s| s.reporters.len() == n) {
        for (acc, &v) in sum.iter_mut().zip(&s.reporters) {
            *acc += (v as f64).max(0.0);
        }
    }
    let total: f64 = sum.iter().sum();
    (total > 0.0).then(|| sum.iter().map(|v| v / total).collect())
}

#[derive(Debug, Clone, PartialEq)]
pub struct Unmixed {
    pub interference: f32,
    pub reporters: Vec<f32>,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Strength {
    Fixed(f64),
    Learned,
}

fn raw_shares(reporters: &[f32]) -> Option<Vec<f64>> {
    let total = summed(reporters);
    (total > 0.0).then(|| {
        reporters
            .iter()
            .map(|&v| (v as f64).max(0.0) / total)
            .collect()
    })
}

pub fn learn_strength(scans: &[ScanInput], eil: &[Option<f64>], profile: &[f64]) -> Option<f64> {
    let n = profile.len();
    let mut groups: std::collections::BTreeMap<u32, Vec<usize>> = std::collections::BTreeMap::new();
    for (i, s) in scans.iter().enumerate() {
        if let (Some(g), Some(_)) = (s.protein, eil[i]) {
            if s.reporters.len() == n && summed(&s.reporters) > 0.0 {
                groups.entry(g).or_default().push(i);
            }
        }
    }
    struct Term {
        q: Vec<f64>,
        r: Vec<f64>,
        e_ref: f64,
        e: f64,
        w: f64,
    }
    let mut terms: Vec<Term> = Vec::new();
    for members in groups.values().filter(|m| m.len() >= 2) {
        let q: Vec<Vec<f64>> = members
            .iter()
            .map(|&i| raw_shares(&scans[i].reporters).unwrap_or_else(|| vec![0.0; n]))
            .collect();
        let e: Vec<f64> = members.iter().map(|&i| eil[i].unwrap_or(0.0)).collect();
        let c: Vec<f64> = e.iter().map(|e| (1.0 - e).powi(2)).collect();
        let mass: f64 = c.iter().sum();
        let mut sum_q = vec![0.0; n];
        for (qk, ck) in q.iter().zip(&c) {
            for (a, &v) in sum_q.iter_mut().zip(qk) {
                *a += ck * v;
            }
        }
        let sum_e: f64 = c.iter().zip(&e).map(|(c, e)| c * e).sum();
        for (k, qk) in q.into_iter().enumerate() {
            let c_ref = mass - c[k];
            if c_ref <= 0.0 {
                continue;
            }
            terms.push(Term {
                r: sum_q
                    .iter()
                    .zip(&qk)
                    .map(|(&a, &v)| (a - c[k] * v) / c_ref)
                    .collect(),
                q: qk,
                e_ref: (sum_e - c[k] * e[k]) / c_ref,
                e: e[k],
                w: c_ref / (c_ref + 1.0),
            });
        }
    }
    let mut lambda = 1.0f64;
    let mut info = (0.0f64, 0.0f64);
    for _ in 0..SLOPE_ROUNDS {
        let (mut num, mut den) = (0.0f64, 0.0f64);
        info = (0.0, 0.0);
        for t in &terms {
            let keep = 1.0 - lambda * t.e_ref;
            let (mut xy, mut xx) = (0.0, 0.0);
            for ((&b, &r), &q) in profile.iter().zip(&t.r).zip(&t.q) {
                let clean = (r - lambda * t.e_ref * b) / keep;
                let x = t.e * (b - clean);
                xy += x * (q - clean);
                xx += x * x;
            }
            num += t.w * xy;
            den += t.w * xx;
            info.0 += t.w * xx;
            info.1 += (t.w * xx) * (t.w * xx);
        }
        if den <= 0.0 {
            return None;
        }
        lambda = (num / den).clamp(0.0, 1.0);
    }
    (info.0 * info.0 / info.1 >= MIN_FIT_SCANS as f64).then_some(lambda)
}

#[derive(Debug, Clone, PartialEq)]
pub struct UnmixRun {
    pub model: Option<InterferenceModel>,
    pub strength: f64,
    pub scans: Vec<Option<Unmixed>>,
}

pub fn unmix_run(scans: &[ScanInput], strength: Strength) -> UnmixRun {
    let none = |model| UnmixRun {
        model,
        strength: 0.0,
        scans: vec![None; scans.len()],
    };
    let Some(model) = InterferenceModel::fit(scans) else {
        return none(None);
    };
    let Some(profile) = background_profile(scans) else {
        return none(None);
    };
    let eil: Vec<Option<f64>> = scans.iter().map(|s| model.interference(s)).collect();
    let lambda = match strength {
        Strength::Fixed(l) => l.clamp(0.0, 1.0),
        Strength::Learned => learn_strength(scans, &eil, &profile).unwrap_or(1.0),
    };
    let out = scans
        .iter()
        .zip(&eil)
        .map(|(s, eil)| {
            let eil = (*eil)?;
            let total = summed(&s.reporters);
            let reporters = s
                .reporters
                .iter()
                .zip(&profile)
                .map(|(&v, &p)| ((v as f64) - lambda * eil * total * p).max(0.0) as f32)
                .collect();
            Some(Unmixed {
                interference: eil as f32,
                reporters,
            })
        })
        .collect();
    UnmixRun {
        model: Some(model),
        strength: lambda,
        scans: out,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn survey(rt: f64, peaks: &[(f64, f32)]) -> Survey<'_> {
        Survey {
            rt,
            peaks,
            detection_limit: 10.0,
        }
    }

    #[test]
    fn window_record_interpolates_between_the_bordering_scans() {
        let before = [
            (499.0, 5.0f32),
            (500.0, 100.0),
            (500.2, 50.0),
            (500.5017, 50.0),
        ];
        let after = [(500.0, 300.0f32), (500.2, 50.0), (500.5017, 150.0)];
        let rec = window_record(
            Some(survey(10.0, &before)),
            Some(survey(20.0, &after)),
            15.0,
            499.8,
            500.6,
            500.0,
            2,
        )
        .unwrap();
        assert!((rec.precursor - 300.0).abs() < 1e-9, "{rec:?}");
        assert!((rec.total - 350.0).abs() < 1e-9, "{rec:?}");
        assert!((rec.precursor_fraction() - 300.0 / 350.0).abs() < 1e-12);
        let dens_before = 4.0 / 10.0;
        let dens_after = 3.0 / 10.0;
        let expect = 0.5 * 10.0 * dens_before * 0.8 + 0.5 * 10.0 * dens_after * 0.8;
        assert!((rec.hidden - expect).abs() < 1e-9, "{rec:?}");
    }

    #[test]
    fn window_record_needs_a_bordering_scan_and_a_window() {
        assert_eq!(window_record(None, None, 1.0, 1.0, 2.0, 1.5, 2), None);
        let p = [(500.0, 1.0f32)];
        assert_eq!(
            window_record(Some(survey(0.0, &p)), None, 1.0, 2.0, 2.0, 1.5, 2),
            None
        );
        let one = window_record(Some(survey(0.0, &p)), None, 1.0, 499.0, 501.0, 500.0, 1).unwrap();
        assert_eq!(one.precursor_fraction(), 1.0);
    }

    fn synthetic(n: usize) -> Vec<ScanInput> {
        let (a, b, c) = (1.0, 0.5, 2.0);
        (0..n)
            .map(|i| {
                let pic = 1e5 * (1.0 + (i % 7) as f64);
                let f = 0.3 + 0.7 * ((i * 37) % 100) as f64 / 100.0;
                let hidden = 1e4 * ((i * 13) % 5) as f64;
                let total = a * pic * f + b * pic * (1.0 - f) + c * hidden;
                let reporters = vec![(total / 4.0) as f32; 4];
                ScanInput {
                    charge: 2,
                    reporters,
                    pic,
                    window: Some(WindowRecord {
                        total: 1000.0,
                        precursor: 1000.0 * f,
                        hidden,
                    }),
                    protein: None,
                }
            })
            .collect()
    }

    #[test]
    fn fit_recovers_the_generating_model() {
        let scans = synthetic(400);
        let m = InterferenceModel::fit(&scans).unwrap();
        assert!((m.precursor[0] - 1.0).abs() < 1e-6, "{m:?}");
        assert!((m.co_isolated - 0.5).abs() < 1e-6, "{m:?}");
        assert!((m.hidden - 2.0).abs() < 1e-6, "{m:?}");
        assert!(m.r2_log > 0.999, "{m:?}");
        let s = &scans[5];
        let w = s.window.unwrap();
        let f = w.precursor_fraction();
        let own = s.pic * f;
        let total = own + 0.5 * s.pic * (1.0 - f) + 2.0 * w.hidden;
        let eil = m.interference(s).unwrap();
        assert!((eil - (1.0 - own / total).min(MAX_INTERFERENCE)).abs() < 1e-6);
    }

    #[test]
    fn too_few_scans_are_not_fitted() {
        let run = unmix_run(&synthetic(MIN_FIT_SCANS - 1), Strength::Fixed(1.0));
        assert!(run.model.is_none());
        assert!(run.scans.iter().all(Option::is_none));
    }

    #[test]
    fn interference_is_removed_with_the_run_profile() {
        let mut scans = synthetic(400);
        scans[0].reporters = vec![400.0, 100.0, 100.0, 100.0];
        let run = unmix_run(&scans, Strength::Fixed(1.0));
        let (model, out) = (run.model.unwrap(), run.scans);
        let eil = model.interference(&scans[0]).unwrap();
        let profile = background_profile(&scans).unwrap();
        let u = out[0].as_ref().unwrap();
        assert!((u.interference as f64 - eil).abs() < 1e-6);
        for (k, (&v, &p)) in scans[0].reporters.iter().zip(&profile).enumerate() {
            let expect = ((v as f64) - eil * 700.0 * p).max(0.0);
            assert!((u.reporters[k] as f64 - expect).abs() < 1e-2, "{k}");
        }
        assert!(u.reporters[1] / u.reporters[0] < 0.25);
    }

    #[test]
    fn scans_without_a_window_record_are_left_out() {
        let mut scans = synthetic(400);
        scans[3].window = None;
        let out = unmix_run(&scans, Strength::Fixed(1.0)).scans;
        assert!(out[3].is_none());
        assert!(out[4].is_some());
    }

    fn proteins(share: f64) -> (Vec<ScanInput>, Vec<Option<f64>>, Vec<f64>) {
        let bg = vec![0.25; 4];
        let mut scans = Vec::new();
        let mut eil = Vec::new();
        for i in 0..300usize {
            let g = i % 30;
            let mut truth = [0.1; 4];
            truth[g % 4] = 0.7;
            let e = 0.05 + 0.7 * ((i * 53) % 100) as f64 / 100.0;
            let mix = share * e;
            let reporters = (0..4)
                .map(|k| (1000.0 * ((1.0 - mix) * truth[k] + mix * bg[k])) as f32)
                .collect();
            scans.push(ScanInput {
                charge: 2,
                reporters,
                pic: 1.0,
                window: None,
                protein: Some(g as u32),
            });
            eil.push(Some(e));
        }
        (scans, eil, bg)
    }

    #[test]
    fn learned_strength_is_full_when_the_interference_is_real() {
        let (scans, eil, bg) = proteins(1.0);
        let l = learn_strength(&scans, &eil, &bg).unwrap();
        assert!((l - 1.0).abs() < 0.02, "{l}");
    }

    #[test]
    fn learned_strength_follows_the_share_that_is_real() {
        let (scans, eil, bg) = proteins(0.5);
        let l = learn_strength(&scans, &eil, &bg).unwrap();
        assert!((l - 0.5).abs() < 0.05, "{l}");
    }

    #[test]
    fn learned_strength_is_zero_when_the_estimates_are_noise() {
        let (scans, eil, bg) = proteins(0.0);
        let l = learn_strength(&scans, &eil, &bg).unwrap();
        assert!(l < 0.05, "{l}");
    }

    #[test]
    fn learned_strength_needs_scans_that_share_proteins() {
        let (mut scans, eil, bg) = proteins(1.0);
        for (i, s) in scans.iter_mut().enumerate() {
            s.protein = Some(i as u32);
        }
        assert_eq!(learn_strength(&scans, &eil, &bg), None);
    }
}
