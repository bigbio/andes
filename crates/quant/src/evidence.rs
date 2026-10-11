//! Template-fit evidence for one (precursor, run) hypothesis.
//!
//! The hypothesis "the precursor elutes near `t_hat` with a Gaussian profile
//! and its theoretical isotope ratios" is fitted directly to the MS1 data, and
//! one log-likelihood ratio summarises the evidence. Identified pairs,
//! transfers and decoys all go through [`fit`], so their statistics are the
//! same kind of number.
//!
//! Data. With `W = 4 s + 2 w_prior`, the MS1 scans with RT in
//! `[t_hat - W, t_hat + W]` give, for isotope `i`, the trace `y_i[j]`: the most
//! intense centroid within `tol_ppm` of `mu + i * 1.003355 / z` (0 when none),
//! and one more trace `y_m1` at `mu - 1.003355 / z`.
//!
//! Model. `y_i[j] = b_i + A * pi_i * g(t_j; t0, w) + noise`, with
//! `g = exp(-(t - t0)^2 / (2 w^2))`, where:
//!
//! * `b_i` is the median of `y_i` over the window;
//! * `sigma_i = 1.4826 * median|y_i[j+1] - y_i[j]| / sqrt(2)`, floored at
//!   `0.5 * L`, where the detection limit `L` is the 5th percentile of all
//!   centroid intensities, at any m/z, in the window's scans (1.0 when there
//!   are none). A missing centroid means "below detection", so a trace with
//!   no centroids is uncertain to about `L`, not to 0;
//! * `pi` is the envelope scaled to sum to 1 over its 4 entries.
//!
//! Fit. Over `t0` in the scan RTs inside `[t_hat - 3 s, t_hat + 3 s]` (the
//! scan nearest `t_hat` when none) and `w` in `{0.7, 1.0, 1.4} * w_prior`, the
//! weighted least-squares amplitude over the used isotopes is
//! `A = max(0, num / den)`, with
//! `num = sum_i pi_i sum_j (y_i[j] - b_i) g_j / sigma_i^2` and
//! `den = sum_i pi_i^2 sum_j g_j^2 / sigma_i^2`, and the log-likelihood ratio
//! against "no peak" is `LLR = A num - A^2 den / 2`. The fit maximises
//! `LLR - (t0 - t_hat)^2 / (2 s^2)`; the first maximum in (t0, w) order wins.
//!
//! Interference. At the optimum, each isotope `i` in `{1, 2, 3}` gets its own
//! amplitude `A_i = sum_j (y_i[j] - b_i) g_j / (pi_i sum_j g_j^2)` with standard
//! error `se_i = sigma_i / (pi_i sqrt(sum_j g_j^2))`. If the largest
//! `(A_i - A) / sqrt(se_i^2 + A_se^2)` exceeds 3, that one isotope is dropped
//! and `A` is refitted at the same `(t0, w)`. The test is one-sided because
//! interference only adds signal. Isotope 0 is never dropped.
//!
//! The per-probe arithmetic, [`fit_traces`], takes the extracted traces as
//! slices, so a batch extraction can reuse it and reproduce [`fit`] exactly.

use model::mass::ISOTOPE;

use crate::decoys::Probe;
use crate::ms1_index::Ms1RunIndex;

/// Isotopes the envelope carries (M+0 to M+3).
pub const MAX_ISOTOPES: usize = 4;
/// A window with fewer MS1 scans than this is not fitted.
pub const MIN_SCANS: usize = 5;
/// Width grid, as multiples of `w_prior`.
const WIDTH_GRID: [f64; 3] = [0.7, 1.0, 1.4];
/// Consistency factor of the median absolute deviation for a Gaussian.
const MAD_SCALE: f64 = 1.4826;
/// Quantile of the window's centroid intensities taken as the detection limit.
const DETECTION_QUANTILE: f64 = 0.05;
/// Detection limit when the window has no centroids, intensity units.
const DEFAULT_DETECTION_LIMIT: f64 = 1.0;
/// Noise floor as a fraction of the detection limit.
const NOISE_FLOOR_OF_LIMIT: f64 = 0.5;
/// One-sided z threshold of the interference test.
const EXCLUSION_Z: f64 = 3.0;
/// Floor, seconds, on `s` and `w_prior` (a zero or non-finite value would
/// divide by zero in the RT prior or the elution profile).
const MIN_SCALE_S: f64 = 1e-3;

/// Extraction settings.
#[derive(Debug, Clone, PartialEq)]
pub struct EvidenceParams {
    /// Centroid tolerance: `|m/z - target| <= target * tol_ppm * 1e-6`.
    pub tol_ppm: f64,
    /// Isotopes extracted and fitted, 1 to [`MAX_ISOTOPES`] (clamped).
    pub n_isotopes: usize,
}

impl Default for EvidenceParams {
    fn default() -> Self {
        Self {
            tol_ppm: 10.0,
            n_isotopes: MAX_ISOTOPES,
        }
    }
}

/// The fitted hypothesis and its summary statistics. Every field is finite
/// for finite input.
#[derive(Debug, Clone, PartialEq)]
pub struct Evidence {
    /// `A`, the total amplitude of the isotope pattern at the apex.
    pub amplitude: f64,
    /// `1 / sqrt(den)` over the used isotopes.
    pub amplitude_se: f64,
    /// `A * w * sqrt(2 pi)`: the area of the elution profile of the whole
    /// pattern (`pi` sums to 1).
    pub quantity: f64,
    /// `amplitude_se * w * sqrt(2 pi)`.
    pub quantity_se: f64,
    /// Log-likelihood ratio against "no peak" at the optimum (0 when `A = 0`).
    pub llr: f64,
    /// `(t0 - t_hat) / s`.
    pub rt_z: f64,
    /// `w / w_prior`.
    pub width_ratio: f64,
    /// Cosine between the per-isotope heights `A_i * pi_i` (each isotope's
    /// own projection on `g`) and `pi` over the used isotopes; 1 when the
    /// observed ratios match the theoretical ones, 0 when either vector is
    /// zero.
    pub envelope_cos: f64,
    /// Minimum over isotopes 0 and 1 of the Pearson correlation between the
    /// trace and `g` over `[t0 - 2w, t0 + 2w]`; 0 with fewer than 3 scans
    /// there, and 0 for a constant trace.
    pub coherence: f64,
    /// The M-1 trace's amplitude projected on `g` (clamped at 0), divided by
    /// `pi_0 * max(A, amplitude_se)`; near 0 when `mu` is monoisotopic. The
    /// floor at one standard error keeps it finite when `A = 0`.
    pub m1_ratio: f64,
    /// Intensity-weighted ppm error of the isotope-0 centroids over
    /// `[t0 - w, t0 + w]` (0 when there are none).
    pub mass_err_ppm: f64,
    /// Scans inside `[t0 - 2w, t0 + 2w]`.
    pub n_points: u32,
    /// Bit `i` set when isotope `i` is in the final fit.
    pub isotopes_used: u8,
    /// Fitted apex, seconds.
    pub t0: f64,
    /// Fitted Gaussian width (standard deviation), seconds.
    pub w: f64,
}

/// The RT prior of a probe, with the floors of [`fit`] applied.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Prior {
    pub t_hat: f64,
    pub s: f64,
    pub w_prior: f64,
}

impl Prior {
    /// `s` and `w_prior` floored at 1 ms (a NaN becomes the floor).
    pub fn of(probe: &Probe) -> Self {
        Self {
            t_hat: probe.t_hat,
            s: probe.s.max(MIN_SCALE_S),
            w_prior: probe.w_prior.max(MIN_SCALE_S),
        }
    }

    /// Half-width `W = 4 s + 2 w_prior` of the extraction window.
    pub fn half_window(&self) -> f64 {
        4.0 * self.s + 2.0 * self.w_prior
    }
}

/// Reusable buffers for [`fit_with`]: the extracted traces and the work space
/// of [`fit_traces`].
#[derive(Debug, Default)]
pub struct FitScratch {
    rt: Vec<f64>,
    y: [Vec<f32>; MAX_ISOTOPES],
    m1: Vec<f32>,
    mz0: Vec<f64>,
    intensities: Vec<f32>,
    work: Workspace,
}

/// Work space of [`fit_traces`]; contents are overwritten on every call.
#[derive(Debug, Default)]
pub struct Workspace {
    sort: Vec<f32>,
    z: Vec<f64>,
    g: Vec<f64>,
}

/// Fit `probe` (monoisotopic m/z, charge `charge`, isotope fractions
/// `envelope`) to `ms1`, which must hold one FAIMS compensation-voltage group.
/// `None` only when the extraction window holds fewer than [`MIN_SCANS`] scans.
pub fn fit(
    ms1: &Ms1RunIndex,
    probe: &Probe,
    charge: u8,
    envelope: &[f32; 4],
    p: &EvidenceParams,
) -> Option<Evidence> {
    fit_with(ms1, probe, charge, envelope, p, &mut FitScratch::default())
}

/// [`fit`] with caller-owned buffers, for fitting many probes without
/// reallocating.
pub fn fit_with(
    ms1: &Ms1RunIndex,
    probe: &Probe,
    charge: u8,
    envelope: &[f32; 4],
    p: &EvidenceParams,
    scratch: &mut FitScratch,
) -> Option<Evidence> {
    debug_assert!(charge > 0, "charge must be positive");
    let prior = Prior::of(probe);
    let half = prior.half_window();
    let scans = ms1.scan_window(prior.t_hat - half, prior.t_hat + half);
    if scans.len() < MIN_SCANS {
        return None;
    }
    let n_iso = p.n_isotopes.clamp(1, MAX_ISOTOPES);
    let detection_limit = detection_limit(ms1, scans.clone(), &mut scratch.intensities);
    traces(
        ms1,
        scans,
        probe.mz,
        charge.max(1),
        n_iso,
        p.tol_ppm,
        scratch,
    );
    let FitScratch {
        rt,
        y,
        m1,
        mz0,
        work,
        ..
    } = scratch;
    let mut ys: [&[f32]; MAX_ISOTOPES] = [&[]; MAX_ISOTOPES];
    for (slot, trace) in ys.iter_mut().zip(y.iter()).take(n_iso) {
        *slot = trace;
    }
    Some(fit_traces(
        rt,
        &ys[..n_iso],
        m1,
        mz0,
        probe.mz,
        envelope,
        &prior,
        detection_limit,
        work,
    ))
}

/// The 5th percentile (nearest rank) of every centroid intensity in `scans`,
/// or 1.0 when they hold no centroid.
fn detection_limit(ms1: &Ms1RunIndex, scans: std::ops::Range<usize>, buf: &mut Vec<f32>) -> f64 {
    buf.clear();
    for j in scans {
        buf.extend(ms1.peaks(j).iter().map(|p| p.1));
    }
    if buf.is_empty() {
        return DEFAULT_DETECTION_LIMIT;
    }
    let rank = (DETECTION_QUANTILE * (buf.len() - 1) as f64).round() as usize;
    let (_, v, _) = buf.select_nth_unstable_by(rank, f32::total_cmp);
    *v as f64
}

/// Fill the scratch traces over `scans`: RT, the isotope traces 0..`n_iso`,
/// the M-1 trace and the isotope-0 centroid m/z (0.0 when none).
fn traces(
    ms1: &Ms1RunIndex,
    scans: std::ops::Range<usize>,
    mu: f64,
    charge: u8,
    n_iso: usize,
    tol_ppm: f64,
    out: &mut FitScratch,
) {
    let step = ISOTOPE / charge as f64;
    let slice = |target: f64| {
        let tol = target * tol_ppm * 1e-6;
        (target - tol, target + tol)
    };
    out.rt.clear();
    out.rt.extend(scans.clone().map(|j| ms1.rt(j)));
    let (lo, hi) = slice(mu);
    out.mz0.clear();
    out.y[0].clear();
    for j in scans.clone() {
        let (mz, inten) = ms1.apex_centroid(j, lo, hi).unwrap_or((0.0, 0.0));
        out.mz0.push(mz);
        out.y[0].push(inten);
    }
    for i in 1..n_iso {
        let (lo, hi) = slice(mu + i as f64 * step);
        out.y[i].clear();
        out.y[i].extend(scans.clone().map(|j| ms1.max_in(j, lo, hi)));
    }
    let (lo, hi) = slice(mu - step);
    out.m1.clear();
    out.m1.extend(scans.map(|j| ms1.max_in(j, lo, hi)));
}

/// The template fit on extracted traces.
///
/// `rt` is ascending and non-empty; `traces[i]` is isotope `i` (1 to
/// [`MAX_ISOTOPES`] traces), `m1` the M-1 trace and `mz0` the isotope-0
/// centroid m/z (0.0 when none), all of `rt`'s length. `mu` is the probe's
/// monoisotopic m/z; `detection_limit` is `L` (see the module docs; a
/// non-finite or non-positive value is replaced by 1.0). A degenerate envelope (non-finite entries, a
/// non-positive sum or `pi_0 = 0`) is treated as `[1, 0, 0, 0]`; isotopes with
/// `pi_i = 0` are not used.
#[allow(clippy::too_many_arguments)]
pub fn fit_traces(
    rt: &[f64],
    traces: &[&[f32]],
    m1: &[f32],
    mz0: &[f64],
    mu: f64,
    envelope: &[f32; 4],
    prior: &Prior,
    detection_limit: f64,
    work: &mut Workspace,
) -> Evidence {
    debug_assert!(!rt.is_empty() && !traces.is_empty() && traces.len() <= MAX_ISOTOPES);
    let n_iso = traces.len().min(MAX_ISOTOPES);
    let traces = &traces[..n_iso];
    let pi = normalised(envelope);
    let limit = if detection_limit.is_finite() && detection_limit > 0.0 {
        detection_limit
    } else {
        DEFAULT_DETECTION_LIMIT
    };
    let floor = NOISE_FLOOR_OF_LIMIT * limit;
    let mut b = [0.0f64; MAX_ISOTOPES];
    let mut sigma = [1.0f64; MAX_ISOTOPES];
    let mut mask = 0u8;
    for i in 0..n_iso {
        (b[i], sigma[i]) = background_and_noise(traces[i], floor, &mut work.sort);
        if pi[i] > 0.0 {
            mask |= 1 << i;
        }
    }
    let (t0, w) = best_alignment(rt, traces, &pi, &b, &sigma, mask, prior, work);
    let g = &mut work.g;
    g.clear();
    g.extend(rt.iter().map(|&t| gaussian(t, t0, w)));
    let sgg: f64 = g.iter().map(|v| v * v).sum();
    let mut proj = [0.0f64; MAX_ISOTOPES];
    for i in 0..n_iso {
        proj[i] = project(traces[i], b[i], g);
    }
    let mask = exclude_interfered(&proj, sgg, &pi, &sigma, mask);
    let (m1_b, _) = background_and_noise(m1, floor, &mut work.sort);
    let m1_proj = project(m1, m1_b, g);
    features(
        rt, traces, mz0, mu, prior, &pi, &sigma, &proj, sgg, m1_proj, mask, t0, w, g,
    )
}

/// `envelope / sum(envelope)`, or `[1, 0, 0, 0]` when that is not a valid
/// distribution with `pi_0 > 0`.
fn normalised(envelope: &[f32; 4]) -> [f64; MAX_ISOTOPES] {
    let e = envelope.map(|v| v as f64);
    let sum: f64 = e.iter().sum();
    if !sum.is_finite() || sum <= 0.0 || e.iter().any(|&v| v < 0.0) || e[0] <= 0.0 {
        return [1.0, 0.0, 0.0, 0.0];
    }
    e.map(|v| v / sum)
}

fn gaussian(t: f64, t0: f64, w: f64) -> f64 {
    let d = (t - t0) / w;
    (-0.5 * d * d).exp()
}

/// `sum_j (y_j - b) g_j`.
fn project(y: &[f32], b: f64, g: &[f64]) -> f64 {
    y.iter().zip(g).map(|(&v, &gj)| (v as f64 - b) * gj).sum()
}

/// The trace's background (median over the window) and noise
/// (`1.4826 * median|successive difference| / sqrt(2)`, floored at `floor`).
fn background_and_noise(y: &[f32], floor: f64, sort: &mut Vec<f32>) -> (f64, f64) {
    if y.is_empty() {
        return (0.0, floor);
    }
    sort.clear();
    sort.extend_from_slice(y);
    sort.sort_unstable_by(f32::total_cmp);
    let b = median_sorted(sort);
    sort.clear();
    sort.extend(y.windows(2).map(|w| (w[1] - w[0]).abs()));
    sort.sort_unstable_by(f32::total_cmp);
    let mad = median_sorted(sort);
    let sigma = (MAD_SCALE * mad / std::f64::consts::SQRT_2).max(floor);
    (b, sigma)
}

/// Median of an ascending slice (mean of the middle two for an even
/// length; 0 when empty).
fn median_sorted(v: &[f32]) -> f64 {
    match v.len() {
        0 => 0.0,
        n if n % 2 == 1 => v[n / 2] as f64,
        n => 0.5 * (v[n / 2 - 1] as f64 + v[n / 2] as f64),
    }
}

/// The (t0, w) maximising `LLR - (t0 - t_hat)^2 / (2 s^2)` over the grid, with
/// every isotope in `mask`.
#[allow(clippy::too_many_arguments)]
fn best_alignment(
    rt: &[f64],
    traces: &[&[f32]],
    pi: &[f64; MAX_ISOTOPES],
    b: &[f64; MAX_ISOTOPES],
    sigma: &[f64; MAX_ISOTOPES],
    mask: u8,
    prior: &Prior,
    work: &mut Workspace,
) -> (f64, f64) {
    // z_j = sum_i pi_i (y_i[j] - b_i) / sigma_i^2, so num = sum_j z_j g_j and
    // den = d * sum_j g_j^2.
    let z = &mut work.z;
    z.clear();
    z.resize(rt.len(), 0.0);
    let mut d = 0.0;
    for (i, trace) in traces.iter().enumerate() {
        if mask & (1 << i) == 0 {
            continue;
        }
        let c = pi[i] / (sigma[i] * sigma[i]);
        d += pi[i] * c;
        for (zj, &v) in z.iter_mut().zip(trace.iter()) {
            *zj += c * (v as f64 - b[i]);
        }
    }
    let (t_hat, s) = (prior.t_hat, prior.s);
    let lo = rt.partition_point(|&t| t < t_hat - 3.0 * s);
    let hi = rt.partition_point(|&t| t <= t_hat + 3.0 * s);
    let candidates = if lo < hi {
        lo..hi
    } else {
        let k = nearest(rt, t_hat);
        k..k + 1
    };
    let mut best = (f64::NEG_INFINITY, rt[candidates.start], prior.w_prior);
    for k in candidates {
        let t0 = rt[k];
        let penalty = (t0 - t_hat) * (t0 - t_hat) / (2.0 * s * s);
        for f in WIDTH_GRID {
            let w = f * prior.w_prior;
            let (mut num, mut sgg) = (0.0, 0.0);
            for (&t, &zj) in rt.iter().zip(z.iter()) {
                let g = gaussian(t, t0, w);
                num += zj * g;
                sgg += g * g;
            }
            let den = d * sgg;
            let a = (num / den).max(0.0);
            let objective = a * num - 0.5 * a * a * den - penalty;
            if objective > best.0 {
                best = (objective, t0, w);
            }
        }
    }
    (best.1, best.2)
}

/// Index of the RT nearest `t` (the earlier one on a tie); `rt` non-empty.
fn nearest(rt: &[f64], t: f64) -> usize {
    let i = rt.partition_point(|&x| x < t);
    if i == 0 {
        0
    } else if i == rt.len() || t - rt[i - 1] <= rt[i] - t {
        i - 1
    } else {
        i
    }
}

/// `num` and `den` over the isotopes in `mask`, from the projections
/// `proj[i] = sum_j (y_i[j] - b_i) g_j` and `sgg = sum_j g_j^2`.
fn num_den(
    proj: &[f64; MAX_ISOTOPES],
    sgg: f64,
    pi: &[f64; MAX_ISOTOPES],
    sigma: &[f64; MAX_ISOTOPES],
    mask: u8,
) -> (f64, f64) {
    let (mut num, mut den) = (0.0, 0.0);
    for i in 0..MAX_ISOTOPES {
        if mask & (1 << i) != 0 {
            let inv = 1.0 / (sigma[i] * sigma[i]);
            num += pi[i] * proj[i] * inv;
            den += pi[i] * pi[i] * sgg * inv;
        }
    }
    (num, den)
}

/// Drop the one isotope of {1, 2, 3} whose own amplitude exceeds the joint
/// amplitude by more than [`EXCLUSION_Z`] combined standard errors, if any.
fn exclude_interfered(
    proj: &[f64; MAX_ISOTOPES],
    sgg: f64,
    pi: &[f64; MAX_ISOTOPES],
    sigma: &[f64; MAX_ISOTOPES],
    mask: u8,
) -> u8 {
    let (num, den) = num_den(proj, sgg, pi, sigma, mask);
    let a = (num / den).max(0.0);
    let a_se2 = 1.0 / den;
    let mut worst: Option<(usize, f64)> = None;
    for i in 1..MAX_ISOTOPES {
        if mask & (1 << i) == 0 {
            continue;
        }
        let a_i = proj[i] / (pi[i] * sgg);
        let se_i = sigma[i] / (pi[i] * sgg.sqrt());
        let score = (a_i - a) / (se_i * se_i + a_se2).sqrt();
        if worst.is_none_or(|(_, s)| score > s) {
            worst = Some((i, score));
        }
    }
    match worst {
        Some((i, score)) if score > EXCLUSION_Z => mask & !(1 << i),
        _ => mask,
    }
}

/// The [`Evidence`] at the final (t0, w) and isotope set.
#[allow(clippy::too_many_arguments)]
fn features(
    rt: &[f64],
    traces: &[&[f32]],
    mz0: &[f64],
    mu: f64,
    prior: &Prior,
    pi: &[f64; MAX_ISOTOPES],
    sigma: &[f64; MAX_ISOTOPES],
    proj: &[f64; MAX_ISOTOPES],
    sgg: f64,
    m1_proj: f64,
    mask: u8,
    t0: f64,
    w: f64,
    g: &[f64],
) -> Evidence {
    let (num, den) = num_den(proj, sgg, pi, sigma, mask);
    let amplitude = (num / den).max(0.0);
    let amplitude_se = 1.0 / den.sqrt();
    let llr = if amplitude > 0.0 {
        amplitude * num - 0.5 * amplitude * amplitude * den
    } else {
        0.0
    };
    let area = w * (2.0 * std::f64::consts::PI).sqrt();

    let (mut dot, mut hh, mut pp) = (0.0, 0.0, 0.0);
    for i in 0..MAX_ISOTOPES {
        if mask & (1 << i) != 0 {
            let h = proj[i] / sgg;
            dot += h * pi[i];
            hh += h * h;
            pp += pi[i] * pi[i];
        }
    }
    let envelope_cos = if hh > 0.0 && pp > 0.0 {
        dot / (hh * pp).sqrt()
    } else {
        0.0
    };

    let lo = rt.partition_point(|&t| t < t0 - 2.0 * w);
    let hi = rt.partition_point(|&t| t <= t0 + 2.0 * w);
    let n_points = hi - lo;
    let coherence = if n_points < 3 {
        0.0
    } else {
        traces
            .iter()
            .take(2)
            .map(|y| pearson(&y[lo..hi], &g[lo..hi]))
            .fold(f64::INFINITY, f64::min)
    };

    let m1_ratio = m1_proj.max(0.0) / sgg / (pi[0] * amplitude.max(amplitude_se));

    let lo1 = rt.partition_point(|&t| t < t0 - w);
    let hi1 = rt.partition_point(|&t| t <= t0 + w);
    let (mut sw, mut swe) = (0.0, 0.0);
    for (&mz, &v) in mz0[lo1..hi1].iter().zip(&traces[0][lo1..hi1]) {
        if mz > 0.0 && v > 0.0 && mu > 0.0 {
            sw += v as f64;
            swe += v as f64 * (mz - mu) / mu * 1e6;
        }
    }
    let mass_err_ppm = if sw > 0.0 { swe / sw } else { 0.0 };

    Evidence {
        amplitude,
        amplitude_se,
        quantity: amplitude * area,
        quantity_se: amplitude_se * area,
        llr,
        rt_z: (t0 - prior.t_hat) / prior.s,
        width_ratio: w / prior.w_prior,
        envelope_cos,
        coherence,
        m1_ratio,
        mass_err_ppm,
        n_points: n_points as u32,
        isotopes_used: mask,
        t0,
        w,
    }
}

/// Pearson correlation of `y` and `g`; 0 when either is constant.
fn pearson(y: &[f32], g: &[f64]) -> f64 {
    let n = y.len() as f64;
    let my = y.iter().map(|&v| v as f64).sum::<f64>() / n;
    let mg = g.iter().sum::<f64>() / n;
    let (mut sxy, mut sxx, mut syy) = (0.0, 0.0, 0.0);
    for (&v, &gj) in y.iter().zip(g) {
        let (dy, dg) = (v as f64 - my, gj - mg);
        sxy += dy * dg;
        sxx += dy * dy;
        syy += dg * dg;
    }
    if sxx > 0.0 && syy > 0.0 {
        (sxy / (sxx * syy).sqrt()).clamp(-1.0, 1.0)
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::decoys::{Hypothesis, PEPTIDE_MASS_SLOPE};
    use crate::ms1_index::Ms1Scan;
    use model::mass::ISOTOPE;

    /// (m/z, charge, isotope fractions, apex RT, Gaussian width, amplitude).
    /// The amplitude is `A` of the model: isotope `i` peaks at `A * pi_i`.
    type Feature = (f64, u8, [f32; 4], f64, f64, f64);

    const ENV: [f32; 4] = [0.5, 0.3, 0.15, 0.05];
    const MZ: f64 = 650.3271;
    const AMP: f64 = 1.0e6;
    /// 1% of the most intense isotope's apex (`AMP * ENV[0]`).
    const NOISE: f32 = 5.0e3;

    struct Rng(u64);

    impl Rng {
        fn next_u64(&mut self) -> u64 {
            self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut x = self.0;
            x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            x ^ (x >> 31)
        }

        fn uniform(&mut self) -> f64 {
            ((self.next_u64() >> 11) as f64 + 0.5) / (1u64 << 53) as f64
        }

        fn normal(&mut self) -> f64 {
            let (u, v) = (self.uniform(), self.uniform());
            (-2.0 * u.ln()).sqrt() * (2.0 * std::f64::consts::PI * v).cos()
        }
    }

    /// MS1 scans every `scan_dt` seconds from 0 to past the last feature.
    /// Every feature puts a centroid at its M-1 position and at each of its 4
    /// isotopes in every scan: a baseline of `5 * noise`, plus the Gaussian
    /// elution profile, plus Gaussian noise of standard deviation `noise`.
    /// Non-positive intensities are dropped, as centroiding would.
    fn synth_scans(features: &[Feature], noise: f32, scan_dt: f64, seed: u64) -> Vec<Ms1Scan> {
        let noise = noise as f64;
        let end = features
            .iter()
            .map(|f| f.3 + 10.0 * f.4)
            .fold(60.0, f64::max)
            + 60.0;
        let mut rng = Rng(seed);
        let n = (end / scan_dt) as usize + 1;
        (0..n)
            .map(|k| {
                let t = k as f64 * scan_dt;
                let mut peaks = Vec::with_capacity(features.len() * 5);
                for &(mz, z, env, apex, width, amp) in features {
                    let g = (-(t - apex) * (t - apex) / (2.0 * width * width)).exp();
                    for iso in -1i32..4 {
                        let signal = if iso >= 0 {
                            amp * env[iso as usize] as f64 * g
                        } else {
                            0.0
                        };
                        let inten = 5.0 * noise + signal + noise * rng.normal();
                        if inten > 0.0 {
                            peaks.push((mz + iso as f64 * ISOTOPE / z as f64, inten as f32));
                        }
                    }
                }
                Ms1Scan {
                    rt: t,
                    peaks,
                    ..Ms1Scan::default()
                }
            })
            .collect()
    }

    fn synth_run(features: &[Feature], noise: f32, scan_dt: f64) -> Ms1RunIndex {
        Ms1RunIndex::new(synth_scans(features, noise, scan_dt, 1))
    }

    fn probe(mz: f64, t_hat: f64, s: f64, w_prior: f64) -> Probe {
        Probe {
            kind: Hypothesis::Target,
            mz,
            t_hat,
            s,
            w_prior,
        }
    }

    fn params() -> EvidenceParams {
        EvidenceParams {
            tol_ppm: 10.0,
            n_isotopes: 4,
        }
    }

    fn assert_finite(e: &Evidence) {
        for (name, v) in [
            ("amplitude", e.amplitude),
            ("amplitude_se", e.amplitude_se),
            ("quantity", e.quantity),
            ("quantity_se", e.quantity_se),
            ("llr", e.llr),
            ("rt_z", e.rt_z),
            ("width_ratio", e.width_ratio),
            ("envelope_cos", e.envelope_cos),
            ("coherence", e.coherence),
            ("m1_ratio", e.m1_ratio),
            ("mass_err_ppm", e.mass_err_ppm),
            ("t0", e.t0),
            ("w", e.w),
        ] {
            assert!(v.is_finite(), "{name} = {v} is not finite");
        }
    }

    #[test]
    fn recovers_amplitude_apex_and_width() {
        let dt = 1.0;
        let apex = 300.4;
        let run = synth_run(&[(MZ, 2, ENV, apex, 5.0, AMP)], NOISE, dt);
        let e = fit(&run, &probe(MZ, 303.0, 6.0, 5.0), 2, &ENV, &params()).unwrap();
        eprintln!(
            "recover: A/A_true-1 = {:.4}, t0-apex = {:.3} s, width_ratio = {}, llr = {:.0}, \
             envelope_cos = {:.4}, coherence = {:.4}, m1_ratio = {:.4}, used = {:#06b}",
            e.amplitude / AMP - 1.0,
            e.t0 - apex,
            e.width_ratio,
            e.llr,
            e.envelope_cos,
            e.coherence,
            e.m1_ratio,
            e.isotopes_used
        );
        assert_finite(&e);
        assert!((e.amplitude / AMP - 1.0).abs() < 0.05);
        assert!((e.t0 - apex).abs() < dt);
        assert!((0.7..=1.4).contains(&e.width_ratio));
        assert_eq!(e.isotopes_used, 0b1111);
        assert!(e.envelope_cos > 0.99 && e.coherence > 0.9);
    }

    /// Upper bound on the LLR of a noise-only window: 1.5 times the largest
    /// null LLR measured over 100 noise seeds of the absent-precursor run at
    /// 1 s and 2 s scan spacing (2.75).
    const NULL_LLR_BOUND: f64 = 4.1;

    #[test]
    fn absent_precursor_has_small_llr() {
        let run = synth_run(&[(MZ, 2, ENV, 300.0, 5.0, 0.0)], NOISE, 1.0);
        let e = fit(&run, &probe(MZ, 300.0, 6.0, 5.0), 2, &ENV, &params()).unwrap();
        eprintln!(
            "absent: llr = {:.2}, amplitude = {:.0} ({:.4} of AMP, {:.2} SE)",
            e.llr,
            e.amplitude,
            e.amplitude / AMP,
            e.amplitude / e.amplitude_se
        );
        assert_finite(&e);
        assert!(e.llr < NULL_LLR_BOUND);
        assert!(e.amplitude < 0.05 * AMP);
    }

    #[test]
    fn interfered_isotope_is_excluded_and_quantity_unbiased() {
        let w = 5.0;
        let iso2 = MZ + 2.0 * ISOTOPE / 2.0;
        let interferer = (
            iso2 * (1.0 + 3e-6),
            2,
            [1.0, 0.0, 0.0, 0.0],
            300.0,
            w,
            5.0 * AMP * ENV[2] as f64,
        );
        let run = synth_run(&[(MZ, 2, ENV, 300.0, w, AMP), interferer], NOISE, 1.0);
        let e = fit(&run, &probe(MZ, 300.0, 6.0, w), 2, &ENV, &params()).unwrap();
        let truth = AMP * w * (2.0 * std::f64::consts::PI).sqrt();
        eprintln!(
            "interfered: used = {:#06b}, quantity/true-1 = {:.4}, w = {}",
            e.isotopes_used,
            e.quantity / truth - 1.0,
            e.w
        );
        assert_finite(&e);
        assert_eq!(e.isotopes_used & 0b0100, 0);
        assert_eq!(e.isotopes_used & 0b0001, 1);
        assert!((e.quantity / truth - 1.0).abs() < 0.1);
    }

    #[test]
    fn rt_prior_prefers_peak_near_prediction() {
        let dt = 1.0;
        let run = synth_run(
            &[(MZ, 2, ENV, 300.0, 1.0, AMP), (MZ, 2, ENV, 305.0, 1.0, AMP)],
            0.0,
            dt,
        );
        let near = fit(&run, &probe(MZ, 300.0, 3.0, 1.0), 2, &ENV, &params()).unwrap();
        let far = fit(&run, &probe(MZ, 305.0, 3.0, 1.0), 2, &ENV, &params()).unwrap();
        eprintln!(
            "rt prior: t0 = {} (t_hat 300), {} (t_hat 305)",
            near.t0, far.t0
        );
        assert!((near.t0 - 300.0).abs() < dt);
        assert!((far.t0 - 305.0).abs() < dt);
    }

    #[test]
    fn m1_ratio_flags_non_monoisotopic_hypothesis() {
        let run = synth_run(&[(MZ, 2, ENV, 300.0, 5.0, AMP)], NOISE, 1.0);
        let mono = fit(&run, &probe(MZ, 300.0, 6.0, 5.0), 2, &ENV, &params()).unwrap();
        let plus1 = fit(
            &run,
            &probe(MZ + ISOTOPE / 2.0, 300.0, 6.0, 5.0),
            2,
            &ENV,
            &params(),
        )
        .unwrap();
        eprintln!(
            "m1_ratio: mono = {:.4}, +1 = {:.4}",
            mono.m1_ratio, plus1.m1_ratio
        );
        assert_finite(&mono);
        assert_finite(&plus1);
        assert!(mono.m1_ratio < 0.1);
        assert!(plus1.m1_ratio > 0.5);
    }

    #[test]
    fn decoy_probe_in_empty_region_scores_below_target() {
        let decoy_mz = MZ + 7.0 * PEPTIDE_MASS_SLOPE / 2.0;
        let run = synth_run(
            &[
                (MZ, 2, ENV, 300.0, 5.0, AMP),
                (decoy_mz, 2, ENV, 300.0, 5.0, 0.0),
            ],
            NOISE,
            1.0,
        );
        let t = fit(&run, &probe(MZ, 300.0, 6.0, 5.0), 2, &ENV, &params()).unwrap();
        let mut dp = probe(decoy_mz, 300.0, 6.0, 5.0);
        dp.kind = Hypothesis::MassShiftDecoy(7);
        let d = fit(&run, &dp, 2, &ENV, &params()).unwrap();
        eprintln!("decoy: llr target = {:.0}, llr decoy = {:.2}", t.llr, d.llr);
        assert_finite(&d);
        assert!(t.llr > 20.0 * d.llr.max(1.0));
    }

    #[test]
    fn extract_uses_donor_faims_group() {
        let at = |cv: f32, offset: f64, amp: f64, seed: u64| {
            synth_scans(&[(MZ, 2, ENV, 300.0, 5.0, amp)], NOISE, 2.0, seed)
                .into_iter()
                .map(move |mut s| {
                    s.rt += offset;
                    s.faims_cv = Some(cv);
                    s
                })
        };
        let scans: Vec<Ms1Scan> = at(-45.0, 0.0, AMP, 1)
            .chain(at(-60.0, 1.0, 0.0, 2))
            .collect();
        let groups = Ms1RunIndex::new(scans).split_by_faims();
        let group = |cv: f32| &groups.iter().find(|g| g.0 == cv).unwrap().1;
        let pr = probe(MZ, 300.0, 6.0, 5.0);
        let on = fit(group(-45.0), &pr, 2, &ENV, &params()).unwrap();
        let off = fit(group(-60.0), &pr, 2, &ENV, &params()).unwrap();
        eprintln!(
            "faims: -45 llr = {:.0}, A/A_true = {:.4}; -60 llr = {:.2}",
            on.llr,
            on.amplitude / AMP,
            off.llr
        );
        assert!((on.amplitude / AMP - 1.0).abs() < 0.1);
        assert!(off.llr < NULL_LLR_BOUND);
    }

    #[test]
    fn short_window_is_none_and_empty_signal_is_zero() {
        let scans = |n: usize, peaks: Vec<(f64, f32)>| {
            Ms1RunIndex::new(
                (0..n)
                    .map(|k| Ms1Scan {
                        rt: k as f64,
                        peaks: peaks.clone(),
                        ..Ms1Scan::default()
                    })
                    .collect(),
            )
        };
        // W = 4 * 0.25 + 2 * 0.25 = 1.5 s holds 4 scans at 1 s spacing.
        let pr = probe(MZ, 10.5, 0.25, 0.25);
        assert_eq!(fit(&scans(30, vec![]), &pr, 2, &ENV, &params()), None);
        let wide = probe(MZ, 10.0, 3.0, 2.0);
        let e = fit(&scans(30, vec![(400.0, 1e5)]), &wide, 2, &ENV, &params()).unwrap();
        assert_finite(&e);
        assert_eq!(e.amplitude, 0.0);
        assert_eq!(e.llr, 0.0);
        assert_eq!(e.coherence, 0.0);
        assert!(e.n_points > 0 && e.isotopes_used & 1 == 1);
    }

    #[test]
    fn degenerate_probe_and_envelope_stay_finite() {
        let run = synth_run(&[(MZ, 2, ENV, 300.0, 5.0, AMP)], NOISE, 1.0);
        let wide = probe(MZ, 300.0, 6.0, 5.0);
        assert_eq!(
            fit(&run, &probe(MZ, 300.0, 0.0, 0.0), 2, &ENV, &params()),
            None
        );
        for (pr, env) in [
            (probe(MZ, 300.0, 0.0, 5.0), ENV),
            (probe(MZ, 300.0, f64::NAN, 5.0), ENV),
            (wide, [0.0; 4]),
            (wide, [1.0, 0.0, 0.0, 0.0]),
            (wide, [f32::NAN, 1.0, 0.0, 0.0]),
        ] {
            let e = fit(&run, &pr, 2, &env, &params()).unwrap();
            assert_finite(&e);
        }
        let one = EvidenceParams {
            tol_ppm: 10.0,
            n_isotopes: 1,
        };
        let e = fit(&run, &wide, 2, &ENV, &one).unwrap();
        assert_finite(&e);
        assert_eq!(e.isotopes_used, 1);
    }

    #[test]
    fn fit_is_deterministic() {
        let run = synth_run(&[(MZ, 2, ENV, 300.0, 5.0, AMP)], NOISE, 1.0);
        let pr = probe(MZ, 301.0, 6.0, 5.0);
        let a = fit(&run, &pr, 2, &ENV, &params());
        assert!(a.is_some());
        assert_eq!(a, fit(&run, &pr, 2, &ENV, &params()));
    }

    #[test]
    fn missing_isotope_trace_does_not_collapse_amplitude() {
        let scans = synth_scans(&[(MZ, 2, ENV, 300.0, 5.0, AMP)], NOISE, 1.0, 5);
        let m3 = MZ + 3.0 * ISOTOPE / 2.0;
        let without_m3: Vec<Ms1Scan> = scans
            .iter()
            .cloned()
            .map(|mut s| {
                s.peaks.retain(|p| (p.0 - m3).abs() > 0.01);
                s
            })
            .collect();
        let pr = probe(MZ, 300.0, 6.0, 5.0);
        let intact = fit(&Ms1RunIndex::new(scans), &pr, 2, &ENV, &params()).unwrap();
        let missing = fit(&Ms1RunIndex::new(without_m3), &pr, 2, &ENV, &params()).unwrap();
        eprintln!(
            "missing M+3: A/A_true = {:.4} (intact {:.4}), llr = {:.0} (intact {:.0}), used = {:#06b}",
            missing.amplitude / AMP,
            intact.amplitude / AMP,
            missing.llr,
            intact.llr,
            missing.isotopes_used
        );
        assert_finite(&missing);
        assert!((missing.amplitude / AMP - 1.0).abs() < 0.15);
        assert!(missing.llr > 0.5 * intact.llr && missing.llr < 2.0 * intact.llr);
    }

    #[test]
    fn timing_10k_probes_on_3000_scan_run() {
        let mut rng = Rng(7);
        let features: Vec<Feature> = (0..300)
            .map(|_| {
                (
                    400.0 + 800.0 * rng.uniform(),
                    2,
                    ENV,
                    60.0 + 2800.0 * rng.uniform(),
                    5.0,
                    AMP * (0.01 + rng.uniform()),
                )
            })
            .chain([(900.0, 2, ENV, 2900.0, 5.0, 0.0)])
            .collect();
        let mut scans = synth_scans(&features, NOISE, 1.0, 3);
        scans.truncate(3000);
        let run = Ms1RunIndex::new(scans);
        let probes: Vec<Probe> = (0..10_000)
            .map(|k| {
                let f = &features[k % features.len()];
                let mz = if k % 2 == 0 {
                    f.0
                } else {
                    400.0 + 800.0 * rng.uniform()
                };
                probe(mz, f.3 + 4.0 * (rng.uniform() - 0.5), 6.0, 5.0)
            })
            .collect();
        let start = std::time::Instant::now();
        let mut fitted = 0usize;
        for pr in &probes {
            if fit(&run, pr, 2, &ENV, &params()).is_some() {
                fitted += 1;
            }
        }
        eprintln!(
            "timing: {} probes ({} fitted) on {} scans ({} centroids) in {:.3} s",
            probes.len(),
            fitted,
            run.len(),
            run.total_peaks(),
            start.elapsed().as_secs_f64()
        );
        assert_eq!(fitted, probes.len());
    }
}
