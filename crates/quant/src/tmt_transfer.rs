#![allow(clippy::needless_range_loop)]

use crate::ms2::Fragment;
use model::Tolerance;

pub const DEFAULT_WRONG_LOG2: f64 = 0.5;
pub const N_SCAN: usize = 9;
pub const N_XFER: usize = 8;
pub const MIN_TRAIN_ROWS: usize = 200;
pub const MIN_CALIBRATION: usize = 200;
pub const FALLBACK_PPM: f64 = 10.0;
pub const FALLBACK_RT_K: f64 = 3.0;
pub const PPM_BOUNDS: (f64, f64) = (1.0, 25.0);
pub const EMPTY_GAP_LOG10: f64 = 0.7;
pub const MIN_VALUE_ROWS: usize = 100;
const VALUE_L2: f64 = 3.0;
pub const LEAK_BINS: usize = 5;

const ROBUST_SIGMA: f64 = 1.4826;

pub fn quantile(sorted: &[f64], q: f64) -> f64 {
    if sorted.is_empty() {
        return f64::NAN;
    }
    let pos = q.clamp(0.0, 1.0) * (sorted.len() - 1) as f64;
    let lo = pos.floor() as usize;
    let hi = pos.ceil() as usize;
    sorted[lo] + (sorted[hi] - sorted[lo]) * (pos - lo as f64)
}

pub fn median_finite(v: &[f64]) -> f64 {
    let mut s: Vec<f64> = v.iter().copied().filter(|x| x.is_finite()).collect();
    s.sort_by(f64::total_cmp);
    quantile(&s, 0.5)
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PpmWindow {
    pub centre: f64,
    pub half_width: f64,
}

impl PpmWindow {
    pub fn learn(errors: &[f64]) -> Option<Self> {
        let mut e: Vec<f64> = errors
            .iter()
            .copied()
            .filter(|x| x.is_finite() && x.abs() <= 50.0)
            .collect();
        if e.len() < MIN_CALIBRATION {
            return None;
        }
        e.sort_by(f64::total_cmp);
        let centre = quantile(&e, 0.5);
        let mut dev: Vec<f64> = e.iter().map(|x| (x - centre).abs()).collect();
        dev.sort_by(f64::total_cmp);
        let hw = (3.0 * ROBUST_SIGMA * quantile(&dev, 0.5)).clamp(PPM_BOUNDS.0, PPM_BOUNDS.1);
        Some(Self {
            centre,
            half_width: hw,
        })
    }

    pub fn fallback() -> Self {
        Self {
            centre: 0.0,
            half_width: FALLBACK_PPM,
        }
    }

    pub fn contains(&self, ppm: f64) -> bool {
        (ppm - self.centre).abs() <= self.half_width
    }

    pub fn z(&self, ppm: f64) -> f64 {
        (ppm - self.centre).abs() / self.half_width
    }
}

pub fn learn_rt_multiplier(z: &[f64]) -> (f64, bool) {
    let mut a: Vec<f64> = z
        .iter()
        .filter(|x| x.is_finite())
        .map(|x| x.abs())
        .collect();
    if a.len() < MIN_CALIBRATION {
        return (FALLBACK_RT_K, false);
    }
    a.sort_by(f64::total_cmp);
    ((3.0 * ROBUST_SIGMA * quantile(&a, 0.5)).max(1.0), true)
}

pub fn log2_deviation(a: &[f32], b: &[f32]) -> Option<f64> {
    let shared: Vec<usize> = (0..a.len().min(b.len()))
        .filter(|&i| a[i] > 0.0 && b[i] > 0.0)
        .collect();
    if shared.len() < 3 {
        return None;
    }
    let sa: f64 = shared.iter().map(|&i| a[i] as f64).sum();
    let sb: f64 = shared.iter().map(|&i| b[i] as f64).sum();
    let ss: f64 = shared
        .iter()
        .map(|&i| {
            let r = (a[i] as f64 / sa) / (b[i] as f64 / sb);
            let l = r.log2();
            l * l
        })
        .sum();
    Some((ss / shared.len() as f64).sqrt())
}

#[derive(Debug, Clone, PartialEq)]
pub struct RunProfile {
    pub loading: Vec<f64>,
    pub empty: Vec<usize>,
    pub bin_edges: Vec<f64>,
    pub p99: Vec<f64>,
}

fn proportions(v: &[f32]) -> Option<Vec<f64>> {
    let s: f64 = v.iter().map(|&x| x.max(0.0) as f64).sum();
    (s > 0.0).then(|| v.iter().map(|&x| x.max(0.0) as f64 / s).collect())
}

pub fn log_signal(corrected: &[f32]) -> f64 {
    (1.0 + corrected.iter().map(|&x| x.max(0.0) as f64).sum::<f64>()).log10()
}

impl RunProfile {
    pub fn from_identified(corrected: &[Vec<f32>]) -> Self {
        let n = corrected.first().map_or(0, Vec::len);
        let props: Vec<Vec<f64>> = corrected.iter().filter_map(|c| proportions(c)).collect();
        let mut loading = vec![0.0; n];
        for (c, l) in loading.iter_mut().enumerate() {
            let col: Vec<f64> = props.iter().map(|p| p[c]).collect();
            *l = median_finite(&col);
        }
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&a, &b| loading[a].total_cmp(&loading[b]));
        let logs: Vec<f64> = order
            .iter()
            .map(|&i| (loading[i].max(0.0) + 1e-6).log10())
            .collect();
        let mut best = (0usize, 0.0f64);
        for k in 1..n {
            let gap = logs[k] - logs[k - 1];
            if gap > best.1 {
                best = (k, gap);
            }
        }
        let empty: Vec<usize> = if n >= 2 && best.1 >= EMPTY_GAP_LOG10 {
            let mut e: Vec<usize> = order[..best.0].to_vec();
            e.sort_unstable();
            e
        } else {
            Vec::new()
        };
        let mut profile = Self {
            loading,
            empty,
            bin_edges: Vec::new(),
            p99: Vec::new(),
        };
        if profile.empty.is_empty() {
            return profile;
        }
        let mut rows: Vec<(f64, f64)> = corrected
            .iter()
            .filter_map(|c| profile.empty_prop(c).map(|e| (log_signal(c), e)))
            .collect();
        if rows.len() < MIN_CALIBRATION {
            profile.empty.clear();
            return profile;
        }
        rows.sort_by(|a, b| a.0.total_cmp(&b.0));
        let per = rows.len() / LEAK_BINS;
        let mut edges = Vec::new();
        let mut p99 = Vec::new();
        for b in 0..LEAK_BINS {
            let lo = b * per;
            let hi = if b + 1 == LEAK_BINS {
                rows.len()
            } else {
                (b + 1) * per
            };
            let mut e: Vec<f64> = rows[lo..hi].iter().map(|r| r.1).collect();
            e.sort_by(f64::total_cmp);
            p99.push(quantile(&e, 0.99));
            if b + 1 < LEAK_BINS {
                edges.push(rows[hi].0);
            }
        }
        profile.bin_edges = edges;
        profile.p99 = p99;
        profile
    }

    pub fn empty_prop(&self, corrected: &[f32]) -> Option<f64> {
        if self.empty.is_empty() {
            return None;
        }
        let p = proportions(corrected)?;
        Some(self.empty.iter().filter_map(|&c| p.get(c)).sum())
    }

    pub fn leak(&self, corrected: &[f32]) -> bool {
        let Some(e) = self.empty_prop(corrected) else {
            return false;
        };
        if self.p99.is_empty() {
            return false;
        }
        let s = log_signal(corrected);
        let bin = self.bin_edges.iter().take_while(|&&x| s >= x).count();
        e > self.p99[bin.min(self.p99.len() - 1)]
    }

    pub fn loading_cosine(&self, other: &Self) -> f64 {
        let dot: f64 = self
            .loading
            .iter()
            .zip(&other.loading)
            .map(|(a, b)| a * b)
            .sum();
        let na = self.loading.iter().map(|a| a * a).sum::<f64>().sqrt();
        let nb = other.loading.iter().map(|a| a * a).sum::<f64>().sqrt();
        if na > 0.0 && nb > 0.0 {
            dot / (na * nb)
        } else {
            f64::NAN
        }
    }
}

pub fn fragment_evidence(
    ladder: &[Fragment],
    peaks: &[(f32, f32)],
    tol: Tolerance,
    min_mz: f64,
) -> (f64, f64) {
    let mut used: Vec<(f32, bool)> = Vec::new();
    for &(mz, inten) in peaks {
        if (mz as f64) < min_mz || inten.is_nan() || inten <= 0.0 {
            continue;
        }
        let d = tol.as_da(mz as f64);
        let at = ladder.partition_point(|f| f.mz < mz as f64 - d);
        let hit = ladder[at..]
            .iter()
            .take_while(|f| f.mz <= mz as f64 + d)
            .next()
            .is_some();
        used.push((inten, hit));
    }
    let total: f64 = used.iter().map(|u| u.0 as f64).sum();
    if total <= 0.0 {
        return (0.0, 0.0);
    }
    let explained: f64 = used.iter().filter(|u| u.1).map(|u| u.0 as f64).sum();
    used.sort_by(|a, b| b.0.total_cmp(&a.0));
    let top = used.len().min(10);
    let top_hit = used[..top].iter().filter(|u| u.1).count();
    (explained / total, top_hit as f64 / top as f64)
}

#[allow(clippy::too_many_arguments)]
pub fn scan_features(
    corrected: &[f32],
    purity: Option<f32>,
    sn_sum: Option<f32>,
    profile: &RunProfile,
    frag_current: f64,
    frag_top10: f64,
    charge: u8,
    ms2_current: f64,
) -> [f64; N_SCAN] {
    let n = corrected.len().max(1);
    let detected = corrected.iter().filter(|&&x| x > 0.0).count() as f64 / n as f64;
    [
        log_signal(corrected),
        detected,
        purity.map_or(f64::NAN, f64::from),
        sn_sum.map_or(f64::NAN, |s| (1.0 + s.max(0.0) as f64).log10()),
        profile.empty_prop(corrected).unwrap_or(0.0),
        frag_current,
        frag_top10,
        f64::from(charge),
        (1.0 + ms2_current.max(0.0)).log10(),
    ]
}

#[derive(Debug, Clone, PartialEq)]
pub struct Logit {
    fill: Vec<f64>,
    mean: Vec<f64>,
    sd: Vec<f64>,
    w: Vec<f64>,
    b: f64,
}

fn sigmoid(z: f64) -> f64 {
    1.0 / (1.0 + (-z.clamp(-30.0, 30.0)).exp())
}

fn solve(mut a: Vec<f64>, mut b: Vec<f64>, n: usize) -> Option<Vec<f64>> {
    for col in 0..n {
        let piv =
            (col..n).max_by(|&i, &j| a[i * n + col].abs().total_cmp(&a[j * n + col].abs()))?;
        if a[piv * n + col].abs() < 1e-12 {
            return None;
        }
        if piv != col {
            for k in 0..n {
                a.swap(piv * n + k, col * n + k);
            }
            b.swap(piv, col);
        }
        for r in col + 1..n {
            let f = a[r * n + col] / a[col * n + col];
            if f != 0.0 {
                for k in col..n {
                    a[r * n + k] -= f * a[col * n + k];
                }
                b[r] -= f * b[col];
            }
        }
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let s: f64 = (i + 1..n).map(|k| a[i * n + k] * x[k]).sum();
        x[i] = (b[i] - s) / a[i * n + i];
    }
    Some(x)
}

impl Logit {
    pub fn fit(x: &[Vec<f64>], y: &[f64], weights: &[f64], l2: f64) -> Option<Self> {
        Self::fit_min(x, y, weights, l2, MIN_TRAIN_ROWS)
    }

    pub fn fit_min(
        x: &[Vec<f64>],
        y: &[f64],
        weights: &[f64],
        l2: f64,
        min_rows: usize,
    ) -> Option<Self> {
        let n = x.len();
        let pos = y.iter().filter(|&&v| v > 0.5).count();
        let min_class = min_rows / 10;
        if n < min_rows || pos < min_class || n - pos < min_class || y.len() != n {
            return None;
        }
        let d = x[0].len();
        let mut fill = vec![0.0; d];
        let mut mean = vec![0.0; d];
        let mut sd = vec![1.0; d];
        for j in 0..d {
            let col: Vec<f64> = x.iter().map(|r| r[j]).collect();
            fill[j] = median_finite(&col);
            if !fill[j].is_finite() {
                fill[j] = 0.0;
            }
            let vals: Vec<f64> = col
                .iter()
                .map(|&v| if v.is_finite() { v } else { fill[j] })
                .collect();
            mean[j] = vals.iter().sum::<f64>() / n as f64;
            let var = vals.iter().map(|v| (v - mean[j]).powi(2)).sum::<f64>() / n as f64;
            sd[j] = if var > 1e-12 { var.sqrt() } else { 1.0 };
        }
        let mut model = Self {
            fill,
            mean,
            sd,
            w: vec![0.0; d],
            b: 0.0,
        };
        let z: Vec<Vec<f64>> = x.iter().map(|r| model.standardise(r)).collect();
        let p = d + 1;
        let mut theta = vec![0.0; p];
        let wsum: f64 = weights.iter().sum();
        let base =
            (y.iter().zip(weights).map(|(a, w)| a * w).sum::<f64>() / wsum).clamp(1e-4, 1.0 - 1e-4);
        theta[d] = (base / (1.0 - base)).ln();
        for _ in 0..60 {
            let mut grad = vec![0.0; p];
            let mut hess = vec![0.0; p * p];
            for i in 0..n {
                let mut eta = theta[d];
                for j in 0..d {
                    eta += theta[j] * z[i][j];
                }
                let pr = sigmoid(eta);
                let wt = weights[i];
                let r = wt * (pr - y[i]);
                let s = wt * pr * (1.0 - pr);
                for j in 0..d {
                    grad[j] += r * z[i][j];
                    for k in j..d {
                        hess[j * p + k] += s * z[i][j] * z[i][k];
                    }
                    hess[j * p + d] += s * z[i][j];
                }
                grad[d] += r;
                hess[d * p + d] += s;
            }
            for j in 0..p {
                for k in 0..j {
                    hess[j * p + k] = hess[k * p + j];
                }
            }
            for j in 0..d {
                grad[j] += l2 * theta[j];
                hess[j * p + j] += l2;
            }
            hess[d * p + d] += 1e-9;
            let step = solve(hess, grad, p)?;
            let mut max = 0.0f64;
            for j in 0..p {
                theta[j] -= step[j];
                max = max.max(step[j].abs());
            }
            if max < 1e-6 {
                break;
            }
        }
        model.w = theta[..d].to_vec();
        model.b = theta[d];
        Some(model)
    }

    fn standardise(&self, x: &[f64]) -> Vec<f64> {
        x.iter()
            .enumerate()
            .map(|(j, &v)| {
                let v = if v.is_finite() { v } else { self.fill[j] };
                (v - self.mean[j]) / self.sd[j]
            })
            .collect()
    }

    pub fn predict(&self, x: &[f64]) -> f64 {
        let z = self.standardise(x);
        sigmoid(self.b + self.w.iter().zip(&z).map(|(w, v)| w * v).sum::<f64>())
    }
}

#[derive(Debug, Default)]
pub struct Training {
    pub floor_x: Vec<Vec<f64>>,
    pub floor_y: Vec<f64>,
    pub value_x: Vec<Vec<f64>>,
    pub value_y: Vec<f64>,
    pub id_x: Vec<Vec<f64>>,
    pub id_y: Vec<f64>,
}

#[derive(Debug, Clone)]
pub struct TransferRisk {
    floor: Logit,
    value: Option<Logit>,
    identity: Logit,
}

#[derive(Debug, Clone, PartialEq)]
pub enum RiskError {
    TooFewRows { model: &'static str, rows: usize },
}

impl std::fmt::Display for RiskError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RiskError::TooFewRows { model, rows } => write!(
                f,
                "the {model} model has too few training rows ({rows}, at least {MIN_TRAIN_ROWS} with both classes needed)"
            ),
        }
    }
}

impl TransferRisk {
    pub fn fit(t: &Training) -> Result<Self, RiskError> {
        let ones = |n: usize| vec![1.0; n];
        let floor = Logit::fit(&t.floor_x, &t.floor_y, &ones(t.floor_y.len()), 1.0).ok_or(
            RiskError::TooFewRows {
                model: "noise-floor",
                rows: t.floor_y.len(),
            },
        )?;
        let identity = Logit::fit(&t.id_x, &t.id_y, &ones(t.id_y.len()), 1.0).ok_or(
            RiskError::TooFewRows {
                model: "identity-risk",
                rows: t.id_y.len(),
            },
        )?;
        let (mut vx, mut vy) = (Vec::new(), Vec::new());
        for (x, y) in t.value_x.iter().zip(&t.value_y) {
            if identity.predict(x) < 0.5 {
                vx.push(x.clone());
                vy.push(*y);
            }
        }
        let value = Logit::fit_min(&vx, &vy, &ones(vy.len()), VALUE_L2, MIN_VALUE_ROWS);
        Ok(Self {
            floor,
            value,
            identity,
        })
    }

    pub fn value_risk(&self, scan: &[f64; N_SCAN], xfer: &[f64; N_XFER]) -> f64 {
        let mut full = Vec::with_capacity(N_SCAN + N_XFER);
        full.extend_from_slice(scan);
        full.extend_from_slice(xfer);
        match &self.value {
            Some(v) => v.predict(&full),
            None => self.floor.predict(scan),
        }
    }

    pub fn value_fitted(&self) -> bool {
        self.value.is_some()
    }

    pub fn predict(&self, scan: &[f64; N_SCAN], xfer: &[f64; N_XFER]) -> (f64, f64) {
        let (p_id, g, f) = self.parts(scan, xfer);
        (p_id + (1.0 - p_id) * g, f)
    }

    pub fn parts(&self, scan: &[f64; N_SCAN], xfer: &[f64; N_XFER]) -> (f64, f64, f64) {
        let mut full = Vec::with_capacity(N_SCAN + N_XFER);
        full.extend_from_slice(scan);
        full.extend_from_slice(xfer);
        let p_h = self.identity.predict(&full).clamp(1e-9, 1.0 - 1e-9);
        let p_id = (p_h / (1.0 - p_h)).min(1.0);
        let f = self.floor.predict(scan);
        let g = self.value.as_ref().map_or(f, |v| v.predict(&full));
        (p_id, g, f)
    }
}

pub fn calibration(pred: &[f64], observed: &[f64], bins: usize) -> Vec<(f64, f64, usize)> {
    let n = pred.len().min(observed.len());
    if n == 0 || bins == 0 {
        return Vec::new();
    }
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| pred[a].total_cmp(&pred[b]).then(a.cmp(&b)));
    (0..bins.min(n))
        .map(|b| {
            let lo = b * n / bins.min(n);
            let hi = (b + 1) * n / bins.min(n);
            let rows = &order[lo..hi];
            let m = rows.len() as f64;
            (
                rows.iter().map(|&i| pred[i]).sum::<f64>() / m,
                rows.iter().map(|&i| observed[i]).sum::<f64>() / m,
                rows.len(),
            )
        })
        .collect()
}

pub fn quant_q(p_wrong: &[f64], floor: &[f64]) -> Vec<f64> {
    let n = p_wrong.len();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| p_wrong[a].total_cmp(&p_wrong[b]).then(a.cmp(&b)));
    let mut run = vec![0.0; n];
    let mut sum = 0.0;
    for (k, &i) in order.iter().enumerate() {
        sum += p_wrong[i] - floor[i];
        run[k] = sum / (k + 1) as f64;
    }
    for k in (0..n.saturating_sub(1)).rev() {
        run[k] = run[k].min(run[k + 1]);
    }
    let mut out = vec![0.0; n];
    for (k, &i) in order.iter().enumerate() {
        out[i] = run[k].max(0.0);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ppm_window_is_learned_from_the_errors_and_bounded() {
        let errs: Vec<f64> = (0..1000)
            .map(|i| 1.0 + (i % 21) as f64 * 0.1 - 1.0)
            .collect();
        let w = PpmWindow::learn(&errs).unwrap();
        assert!((w.centre - 1.0).abs() < 0.2, "{w:?}");
        assert!((w.half_width - 2.2239).abs() < 0.05, "{w:?}");
        assert!(w.contains(1.5) && !w.contains(30.0));
        assert!(PpmWindow::learn(&errs[..50]).is_none());
        let mut with_outliers = errs.clone();
        with_outliers.extend(std::iter::repeat_n(500.0, 300));
        let w2 = PpmWindow::learn(&with_outliers).unwrap();
        assert!(w2.half_width < 3.0, "{w2:?}");
    }

    #[test]
    fn rt_multiplier_is_three_robust_sigmas_and_falls_back_when_few_pairs() {
        let z: Vec<f64> = (0..1000).map(|i| (i as f64 / 1000.0) * 4.0 - 2.0).collect();
        let (k, learned) = learn_rt_multiplier(&z);
        assert!(learned && (k - 4.4478).abs() < 0.05, "{k}");
        assert_eq!(learn_rt_multiplier(&z[..10]), (FALLBACK_RT_K, false));
        let tight: Vec<f64> = vec![0.1; 500];
        assert_eq!(learn_rt_multiplier(&tight).0, 1.0);
    }

    #[test]
    fn log2_deviation_ignores_scale_and_needs_three_shared_channels() {
        let a = [100.0f32, 200.0, 400.0, 0.0, 800.0];
        let b = [10.0f32, 20.0, 40.0, 5.0, 80.0];
        assert!(log2_deviation(&a, &b).unwrap() < 1e-9);
        let c = [400.0f32, 200.0, 100.0, 0.0, 800.0];
        assert!(log2_deviation(&a, &c).unwrap() > 0.5);
        assert!(log2_deviation(&[1.0, 0.0, 2.0, 0.0], &[1.0, 1.0, 2.0, 0.0]).is_none());
    }

    fn loaded_scan(carrier: f32, empty: f32) -> Vec<f32> {
        vec![carrier, 40.0, 50.0, 60.0, empty, empty]
    }

    #[test]
    fn empty_channels_and_the_identified_p99_price_the_spill() {
        let ident: Vec<Vec<f32>> = (0..1000)
            .map(|i| loaded_scan(100.0 + (i % 50) as f32 * 10.0, 0.2 + (i % 10) as f32 * 0.01))
            .collect();
        let p = RunProfile::from_identified(&ident);
        assert_eq!(p.empty, vec![4, 5]);
        assert_eq!(p.p99.len(), LEAK_BINS);
        assert!(!p.leak(&loaded_scan(300.0, 0.2)));
        assert!(p.leak(&loaded_scan(300.0, 30.0)));
        let full: Vec<Vec<f32>> = (0..300)
            .map(|i| vec![100.0 + i as f32, 90.0, 110.0, 95.0])
            .collect();
        let q = RunProfile::from_identified(&full);
        assert!(q.empty.is_empty() && !q.leak(&[100.0, 1.0, 1.0, 1.0]));
        assert!(p.loading_cosine(&p) > 0.999);
    }

    #[test]
    fn fragment_evidence_counts_ladder_hits_outside_the_reporter_region() {
        let ladder =
            crate::ms2::fragment_ladder(&[71.037, 114.043, 128.095, 113.084, 97.053, 156.101], 2);
        let hit = ladder[2].mz as f32;
        let peaks = vec![(126.13f32, 9999.0f32), (hit, 100.0), (700.5, 100.0)];
        let (cur, top) = fragment_evidence(&ladder, &peaks, Tolerance::Ppm(20.0), 160.0);
        assert!(
            (cur - 0.5).abs() < 1e-6 && (top - 0.5).abs() < 1e-6,
            "{cur} {top}"
        );
        assert_eq!(
            fragment_evidence(&ladder, &[], Tolerance::Ppm(20.0), 160.0),
            (0.0, 0.0)
        );
    }

    struct Lcg(u64);
    impl Lcg {
        fn next(&mut self) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((self.0 >> 33) as f64) / ((1u64 << 31) as f64)
        }
    }

    #[test]
    fn logit_recovers_a_separable_signal_and_fills_missing_covariates() {
        let mut g = Lcg(7);
        let (mut x, mut y) = (Vec::new(), Vec::new());
        for _ in 0..2000 {
            let a = g.next();
            let b = g.next();
            let label = if a + 0.2 * (b - 0.5) > 0.6 { 1.0 } else { 0.0 };
            x.push(vec![a, b, if g.next() < 0.1 { f64::NAN } else { b }]);
            y.push(label);
        }
        let w = vec![1.0; x.len()];
        let m = Logit::fit(&x, &y, &w, 1.0).unwrap();
        assert!(m.predict(&[0.9, 0.5, 0.5]) > 0.9);
        assert!(m.predict(&[0.1, 0.5, f64::NAN]) < 0.1);
        assert!(Logit::fit(&x[..50], &y[..50], &w[..50], 1.0).is_none());
        let one_class = vec![0.0; x.len()];
        assert!(Logit::fit(&x, &one_class, &w, 1.0).is_none());
    }

    fn training(rng: &mut Lcg) -> Training {
        let mut t = Training::default();
        for _ in 0..1500 {
            let q = rng.next();
            let mut scan = [0.0; N_SCAN];
            scan[0] = 3.0 + 2.0 * q;
            scan[1] = 0.5 + 0.5 * q;
            scan[5] = q;
            let wrong = if rng.next() < 0.02 + 0.3 * (1.0 - q) {
                1.0
            } else {
                0.0
            };
            let floor_wrong = if rng.next() < 0.01 + 0.1 * (1.0 - q) {
                1.0
            } else {
                0.0
            };
            t.floor_x.push(scan.to_vec());
            t.floor_y.push(floor_wrong);
            let mut full = scan.to_vec();
            full.extend([0.8, 0.5, 6.0, 0.5, 3.0, 2.0, 1.0, 2.0]);
            t.value_x.push(full.clone());
            t.value_y.push(wrong);
            t.id_x.push(full);
            t.id_y.push(0.0);
            let mut dec = scan.to_vec();
            dec.extend([0.2 * rng.next(), 0.2, 2.0, 0.5, 3.0, 2.0, 1.0, 2.0]);
            t.id_x.push(dec);
            t.id_y.push(1.0);
        }
        t
    }

    #[test]
    fn risk_orders_weak_scans_above_strong_ones_and_flags_decoy_like_candidates() {
        let mut rng = Lcg(3);
        let t = training(&mut rng);
        let risk = TransferRisk::fit(&t).unwrap();
        let strong = {
            let mut s = [0.0; N_SCAN];
            s[0] = 5.0;
            s[1] = 1.0;
            s[5] = 1.0;
            s
        };
        let weak = {
            let mut s = strong;
            s[0] = 3.0;
            s[1] = 0.5;
            s[5] = 0.0;
            s
        };
        let x = [0.8, 0.5, 6.0, 0.5, 3.0, 2.0, 1.0, 2.0];
        let (ps, fs) = risk.predict(&strong, &x);
        let (pw, fw) = risk.predict(&weak, &x);
        assert!(pw > ps && fw > fs, "{pw} {ps} {fw} {fs}");
        let decoy_like = [0.05, 0.2, 2.0, 0.5, 3.0, 2.0, 1.0, 2.0];
        let (pd, _) = risk.predict(&strong, &decoy_like);
        assert!(pd > ps);
        let few = Training {
            value_x: t.value_x[..20].to_vec(),
            value_y: t.value_y[..20].to_vec(),
            ..Training {
                floor_x: t.floor_x.clone(),
                floor_y: t.floor_y.clone(),
                id_x: t.id_x.clone(),
                id_y: t.id_y.clone(),
                ..Training::default()
            }
        };
        let fallback = TransferRisk::fit(&few).unwrap();
        assert!(!fallback.value_fitted());
        let (_, g, f) = fallback.parts(&strong, &x);
        assert_eq!(g, f);
        let small = Training {
            floor_x: t.floor_x[..10].to_vec(),
            floor_y: t.floor_y[..10].to_vec(),
            ..Training::default()
        };
        assert!(matches!(
            TransferRisk::fit(&small),
            Err(RiskError::TooFewRows {
                model: "noise-floor",
                ..
            })
        ));
    }

    #[test]
    fn calibration_groups_by_prediction() {
        let pred = [0.1, 0.9, 0.2, 0.8];
        let obs = [0.0, 1.0, 0.0, 1.0];
        let t = calibration(&pred, &obs, 2);
        assert_eq!(t.len(), 2);
        assert!((t[0].0 - 0.15).abs() < 1e-12 && t[0].1 == 0.0 && t[0].2 == 2);
        assert!((t[1].0 - 0.85).abs() < 1e-12 && t[1].1 == 1.0);
        assert!(calibration(&[], &[], 5).is_empty());
    }

    #[test]
    fn quant_q_is_a_monotone_running_mean_of_the_excess() {
        let p = [0.01, 0.02, 0.5, 0.03, 0.9];
        let f = [0.01, 0.01, 0.01, 0.01, 0.01];
        let q = quant_q(&p, &f);
        assert!(q[0].abs() < 1e-12);
        assert!((q[1] - 0.005).abs() < 1e-12);
        assert!((q[3] - 0.01).abs() < 1e-12);
        assert!(q[2] > q[3] && q[4] > q[2]);
        let q2 = quant_q(&[0.001, 0.2], &[0.05, 0.05]);
        assert_eq!(q2[0], 0.0);
        assert!(q2[1] >= q2[0]);
    }
}
