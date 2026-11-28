//! MS1 feature of a precursor with a wide isotope envelope, such as a cross-linked peptide
//! pair: the elution peak is found on the isotope traces taken together, weighted by the
//! theoretical envelope, and the quantity is a fit of the isotope areas to that envelope that
//! leaves out isotopes another ion adds signal to. The width of a typical elution peak is
//! estimated from the autocorrelation of those traces.

use model::mass::ISOTOPE;

use crate::lfq::{FeatureQuant, LfqParams, LfqTarget};
use crate::ms1::Ms1RunIndex;
use crate::xic::{cosine, find_peak, gaussian_smooth, integrate, savitzky_golay5};

/// An isotope whose area is above this multiple of the area the envelope expects is left out
/// of the fit: interference only adds signal.
const INTERFERENCE_EXCESS: f64 = 1.5;

/// The fewest isotopes a fit keeps, whatever their excess.
const MIN_FIT_ISOTOPES: usize = 2;

/// The elution width is read where the autocorrelation of a pattern trace falls to this
/// share of its value at a lag of one grid step.
const ACF_HALF: f64 = 0.5;

/// The fewest grid steps of a pattern trace that an elution width is estimated from.
const MIN_ACF_POINTS: usize = 8;

/// The pattern signal of one scan: the least-squares amplitude of the envelope in the
/// isotope intensities, scaled by the squared cosine between them, so that a scan where one
/// isotope trace alone carries signal counts for little.
fn pattern_signal(traces: &[Vec<f32>], envelope: &[f64], scan: usize) -> f32 {
    let (mut dot, mut obs2, mut env2) = (0.0f64, 0.0f64, 0.0f64);
    for (trace, &e) in traces.iter().zip(envelope) {
        let v = f64::from(trace[scan]);
        dot += v * e;
        obs2 += v * v;
        env2 += e * e;
    }
    if obs2 <= 0.0 || env2 <= 0.0 {
        return 0.0;
    }
    let amplitude = dot / env2;
    let cos2 = dot * dot / (obs2 * env2);
    (amplitude * cos2) as f32
}

fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}

/// The weighted median of `values`; `None` when the weights sum to nothing.
fn weighted_median(values: &mut [(f64, f64)]) -> Option<f64> {
    let total: f64 = values.iter().map(|v| v.1).sum();
    if !positive(total) {
        return None;
    }
    values.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut cum = 0.0;
    for &(v, w) in values.iter() {
        cum += w;
        if cum >= 0.5 * total {
            return Some(v);
        }
    }
    values.last().map(|v| v.0)
}

/// The amplitude of `envelope` in the isotope `areas`: a least-squares fit over the isotopes
/// whose area is not in excess of what a robust first estimate (the weighted median of the
/// area-to-envelope ratios) expects. With an envelope that sums to one, the amplitude is the
/// summed area of a clean feature.
pub fn envelope_amplitude(areas: &[f64], envelope: &[f64]) -> f64 {
    let mut ratios: Vec<(f64, f64)> = areas
        .iter()
        .zip(envelope)
        .filter(|(_, &e)| e > 0.0)
        .map(|(&a, &e)| (a / e, e))
        .collect();
    let Some(first) = weighted_median(&mut ratios) else {
        return 0.0;
    };
    let mut kept: Vec<(f64, f64)> = areas
        .iter()
        .zip(envelope)
        .filter(|(&a, &e)| e > 0.0 && a <= INTERFERENCE_EXCESS * first * e)
        .map(|(&a, &e)| (a, e))
        .collect();
    if kept.len() < MIN_FIT_ISOTOPES {
        kept = areas
            .iter()
            .zip(envelope)
            .filter(|(_, &e)| e > 0.0)
            .map(|(&a, &e)| (a, e))
            .collect();
        kept.sort_by(|x, y| (x.0 / x.1).total_cmp(&(y.0 / y.1)));
        kept.truncate(MIN_FIT_ISOTOPES);
    }
    let env2: f64 = kept.iter().map(|(_, e)| e * e).sum();
    if env2 <= 0.0 {
        return 0.0;
    }
    kept.iter().map(|(a, e)| a * e).sum::<f64>() / env2
}

struct Extraction {
    rt: Vec<f64>,
    traces: Vec<Vec<f32>>,
    signal: Vec<f32>,
}

/// The isotope traces of `target` at `mono_mz` over its extraction window and their pattern
/// signal; `None` when the window holds no scan.
fn extract(
    index: &Ms1RunIndex,
    target: &LfqTarget,
    mono_mz: f64,
    params: &LfqParams,
) -> Option<Extraction> {
    let envelope = &target.envelope;
    if envelope.is_empty() {
        return None;
    }
    let z = f64::from(target.charge.max(1));
    let window = index.scan_window(
        target.anchor_rt - params.rt_window_s,
        target.anchor_rt + params.rt_window_s,
    );
    if window.is_empty() {
        return None;
    }
    let rt: Vec<f64> = window.clone().map(|i| index.rt(i)).collect();
    let traces: Vec<Vec<f32>> = (0..envelope.len())
        .map(|k| {
            let mz = mono_mz + k as f64 * ISOTOPE / z;
            index.xic(mz, params.tol.as_da(mz), window.clone())
        })
        .collect();
    let signal = (0..rt.len())
        .map(|i| pattern_signal(&traces, envelope, i))
        .collect();
    Some(Extraction { rt, traces, signal })
}

/// The retention times, in seconds, and the pattern signal of `target` over its extraction
/// window.
pub fn pattern_trace(
    index: &Ms1RunIndex,
    target: &LfqTarget,
    params: &LfqParams,
) -> Option<(Vec<f64>, Vec<f32>)> {
    extract(index, target, target.mono_mz, params).map(|e| (e.rt, e.signal))
}

/// The standard deviation, in seconds, of the elution peak in a pattern trace, from its
/// autocorrelation on a grid of `step_s` seconds.
///
/// The autocorrelation of a Gaussian peak of standard deviation `s` at lag `t` is
/// `exp(-t^2 / (4 s^2))`. Noise that is independent between scans adds to lag zero only, so
/// the lags are taken relative to lag one: with `k` the lag, in grid steps, at which the
/// autocorrelation is half of that at one step, `s = step * sqrt((k^2 - 1) / (4 ln 2))`.
/// `None` when the trace is shorter than `MIN_ACF_POINTS` grid steps, holds no signal, or its
/// autocorrelation does not fall to half within half of the trace.
pub fn elution_sd_s(rt: &[f64], signal: &[f32], step_s: f64) -> Option<f64> {
    let (&first, &last) = (rt.first()?, rt.last()?);
    if rt.len() != signal.len() || !positive(step_s) || !positive(last - first) {
        return None;
    }
    let n = ((last - first) / step_s).floor() as usize + 1;
    if n < MIN_ACF_POINTS {
        return None;
    }
    let mut y = Vec::with_capacity(n);
    let mut j = 0;
    for i in 0..n {
        let t = first + i as f64 * step_s;
        while j + 2 < rt.len() && rt[j + 1] < t {
            j += 1;
        }
        let (t0, t1) = (rt[j], rt[j + 1]);
        let w = if t1 > t0 {
            ((t - t0) / (t1 - t0)).clamp(0.0, 1.0)
        } else {
            0.0
        };
        y.push(f64::from(signal[j]) * (1.0 - w) + f64::from(signal[j + 1]) * w);
    }
    let at = |k: usize| -> f64 { y[..n - k].iter().zip(&y[k..]).map(|(a, b)| a * b).sum() };
    let one = at(1);
    if !positive(one) {
        return None;
    }
    let mut previous = 1.0;
    for k in 2..n / 2 {
        let r = at(k) / one;
        if r <= ACF_HALF {
            let crossing = (k - 1) as f64 + (previous - ACF_HALF) / (previous - r);
            let sd = step_s * ((crossing * crossing - 1.0) / (4.0 * std::f64::consts::LN_2)).sqrt();
            return positive(sd).then_some(sd);
        }
        previous = r;
    }
    None
}

/// The feature of `target` at `mono_mz` (the target's own m/z, or a shifted one for a decoy).
///
/// The apex is the local maximum of the smoothed pattern signal next to the scan where that
/// signal times a Gaussian prior around the anchor, of standard deviation `prior_sd_s`
/// seconds, is largest. Every isotope of `target.envelope` is traced; `params.n_isotopes` and
/// `params.min_mono_share` are not read. `area` is the envelope amplitude over
/// `target.isotope_share`, `apex_intensity` the smoothed pattern signal at the apex, and
/// `cosine` compares all traced isotope areas with the envelope.
pub fn quantify_envelope(
    index: &Ms1RunIndex,
    target: &LfqTarget,
    mono_mz: f64,
    params: &LfqParams,
    prior_sd_s: f64,
) -> Option<FeatureQuant> {
    let envelope = &target.envelope;
    let Extraction { rt, traces, signal } = extract(index, target, mono_mz, params)?;
    let smooth = match params.gaussian_width_s {
        Some(width) => gaussian_smooth(&rt, &signal, width),
        None => savitzky_golay5(&signal),
    };
    let sd = prior_sd_s.max(f64::MIN_POSITIVE);
    let weighted = |i: usize| {
        let dz = (rt[i] - target.anchor_rt) / sd;
        f64::from(smooth[i].max(0.0)) * (-0.5 * dz * dz).exp()
    };
    let start = (0..rt.len())
        .max_by(|&a, &b| weighted(a).total_cmp(&weighted(b)).then(b.cmp(&a)))
        .filter(|&i| weighted(i) > 0.0)?;
    let peak = find_peak(&smooth, start, params.max_apex_climb)?;
    let apex_intensity = smooth[peak.apex];
    if apex_intensity <= 0.0 {
        return None;
    }
    let isotope_areas: Vec<f64> = traces
        .iter()
        .map(|t| integrate(&rt, t, peak.left, peak.right))
        .collect();
    let amplitude = envelope_amplitude(&isotope_areas, envelope);
    if !positive(amplitude) || !positive(target.isotope_share) {
        return None;
    }
    let n_found = isotope_areas.iter().filter(|a| **a > 0.0).count();
    Some(FeatureQuant {
        apex_rt: rt[peak.apex],
        rt_start: rt[peak.left],
        rt_stop: rt[peak.right],
        area: amplitude / target.isotope_share,
        apex_intensity,
        cosine: cosine(&isotope_areas, envelope) as f32,
        isotope_areas,
        n_isotopes_found: n_found as u32,
        rt_delta_s: rt[peak.apex] - target.anchor_rt,
        n_scans: (peak.right - peak.left + 1) as u32,
    })
}

#[cfg(test)]
#[path = "envelope/tests/envelope.rs"]
mod tests;
