//! Label-free MS1 quantification of one target in one run.
//!
//! A target is a (peptidoform, charge) with an anchor retention time (its
//! best PSM) and a theoretical isotope envelope. Quantification extracts the
//! first `n_isotopes` isotope chromatograms in an RT window around the anchor,
//! finds the chromatographic peak that holds the anchor on the smoothed
//! monoisotopic trace, integrates every isotope over the same boundaries and
//! scores the envelope against theory. The same routine, run at a shifted m/z,
//! produces the decoy twin for the feature q-value.

use model::mass::ISOTOPE;
use model::Tolerance;

use crate::ms1_index::Ms1RunIndex;
use crate::xic::{cosine, find_peak, integrate, savitzky_golay5};

/// Extraction settings.
#[derive(Debug, Clone)]
pub struct LfqParams {
    /// MS1 m/z tolerance for every isotope trace.
    pub tol: Tolerance,
    /// Half-width of the RT window around the anchor, seconds.
    pub rt_window_s: f64,
    /// Isotopes extracted (monoisotopic + n−1).
    pub n_isotopes: usize,
    /// Envelope cosine below which the driver leaves a feature out of `lfq.tsv`
    /// (the feature itself is still returned and reported in `lfq_features.tsv`).
    pub min_cosine: f32,
    /// m/z shift of the decoy twin (Th); IonQuant uses +11.0, Sage +11.06.
    pub decoy_mz_shift: f64,
    /// Scans the apex may be from the anchor scan.
    pub max_apex_climb: usize,
    /// A feature whose monoisotopic area holds less than this fraction of its
    /// theoretical share of the envelope is not the target ion: the M+1..
    /// traces line up with a different, isotope-shifted species (on PXD001819
    /// a deamidated N-G peptide integrated at ~100x the true value). Median
    /// observed/theoretical share is 1.07 there; 0.9 % of features fall below 0.3.
    pub min_mono_share: f64,
}

impl Default for LfqParams {
    fn default() -> Self {
        Self {
            tol: Tolerance::Ppm(10.0),
            rt_window_s: 60.0,
            n_isotopes: 4,
            min_cosine: 0.7,
            decoy_mz_shift: 11.0,
            max_apex_climb: 12,
            min_mono_share: 0.3,
        }
    }
}

/// One thing to quantify.
#[derive(Debug, Clone)]
pub struct LfqTarget {
    /// Theoretical monoisotopic m/z of the precursor.
    pub mono_mz: f64,
    pub charge: u8,
    /// Retention time of the identifying PSM, seconds.
    pub anchor_rt: f64,
    /// Theoretical isotope envelope (first `n_isotopes`, sums to 1).
    pub envelope: Vec<f64>,
}

/// A quantified feature.
#[derive(Debug, Clone, PartialEq)]
pub struct FeatureQuant {
    pub apex_rt: f64,
    pub rt_start: f64,
    pub rt_stop: f64,
    /// Summed isotope areas (intensity · seconds).
    pub area: f64,
    /// Monoisotopic trace height at the apex scan.
    pub apex_intensity: f32,
    /// Area per isotope, monoisotopic first.
    pub isotope_areas: Vec<f64>,
    /// Cosine similarity of the isotope areas to the theoretical envelope.
    pub cosine: f32,
    /// Isotopes (from the monoisotope up, contiguous) with signal.
    pub n_isotopes_found: u32,
    /// Apex RT minus anchor RT, seconds.
    pub rt_delta_s: f64,
    /// Scans inside the boundaries.
    pub n_scans: u32,
}

/// Quantify `target` in `index`. `None` when no MS1 scan lies in the window or
/// no monoisotopic signal surrounds the anchor.
pub fn quantify(
    index: &Ms1RunIndex,
    target: &LfqTarget,
    params: &LfqParams,
) -> Option<FeatureQuant> {
    quantify_at(index, target, target.mono_mz, params)
}

/// The decoy twin: the same extraction at `mono_mz + decoy_mz_shift`.
pub fn quantify_decoy(
    index: &Ms1RunIndex,
    target: &LfqTarget,
    params: &LfqParams,
) -> Option<FeatureQuant> {
    quantify_at(
        index,
        target,
        target.mono_mz + params.decoy_mz_shift,
        params,
    )
}

fn quantify_at(
    index: &Ms1RunIndex,
    target: &LfqTarget,
    mono_mz: f64,
    params: &LfqParams,
) -> Option<FeatureQuant> {
    let n_iso = params.n_isotopes.max(1);
    let z = target.charge.max(1) as f64;
    let window = index.scan_window(
        target.anchor_rt - params.rt_window_s,
        target.anchor_rt + params.rt_window_s,
    );
    if window.is_empty() {
        return None;
    }
    let anchor_scan = index.nearest_scan(target.anchor_rt)?;
    let anchor = anchor_scan.clamp(window.start, window.end - 1) - window.start;
    let rt: Vec<f64> = window.clone().map(|i| index.rt(i)).collect();

    let traces: Vec<Vec<f32>> = (0..n_iso)
        .map(|k| {
            let mz = mono_mz + k as f64 * ISOTOPE / z;
            index.xic(mz, params.tol.as_da(mz), window.clone())
        })
        .collect();

    let smooth = savitzky_golay5(&traces[0]);
    let peak = find_peak(&smooth, anchor, params.max_apex_climb)?;
    // The apex on the smoothed trace may sit on a zero of the raw trace (a
    // missing scan); fall back to the highest raw value inside the bounds.
    let mut apex_intensity = traces[0][peak.apex];
    if apex_intensity <= 0.0 {
        apex_intensity = traces[0][peak.left..=peak.right]
            .iter()
            .copied()
            .fold(0.0f32, f32::max);
    }
    if apex_intensity <= 0.0 {
        return None;
    }

    let isotope_areas: Vec<f64> = traces
        .iter()
        .map(|t| integrate(&rt, t, peak.left, peak.right))
        .collect();
    let mut n_found = 0u32;
    for a in &isotope_areas {
        if *a > 0.0 {
            n_found += 1;
        } else {
            break;
        }
    }
    let theo: Vec<f64> = target.envelope.iter().copied().take(n_iso).collect();
    let cos = cosine(&isotope_areas, &theo) as f32;
    let area: f64 = isotope_areas.iter().sum();
    if let Some(&theo_mono) = theo.first() {
        let theo_total: f64 = theo.iter().sum();
        if theo_mono > 0.0 && area > 0.0 && theo_total > 0.0 {
            let share = (isotope_areas[0] / area) / (theo_mono / theo_total);
            if share < params.min_mono_share {
                return None;
            }
        }
    }

    Some(FeatureQuant {
        apex_rt: rt[peak.apex],
        rt_start: rt[peak.left],
        rt_stop: rt[peak.right],
        area,
        apex_intensity,
        isotope_areas,
        cosine: cos,
        n_isotopes_found: n_found,
        rt_delta_s: rt[peak.apex] - target.anchor_rt,
        n_scans: (peak.right - peak.left + 1) as u32,
    })
}

/// Feature score for the target/decoy competition: envelope similarity
/// cubed, proximity to the anchor, and relative height (Sage's hybrid rule).
/// `max_apex` is the largest apex among the run's features (normalises the
/// height term to 0..1).
pub fn feature_score(fq: &FeatureQuant, rt_window_s: f64, max_apex: f32) -> f32 {
    let cos = fq.cosine.clamp(0.0, 1.0);
    let prox = (1.0 - (fq.rt_delta_s.abs() / rt_window_s.max(1e-9)).min(1.0)) as f32;
    let height = if max_apex > 0.0 {
        (fq.apex_intensity / max_apex).clamp(0.0, 1.0).sqrt()
    } else {
        0.0
    };
    cos * cos * cos * prox.powf(1.0 / 3.0) * height
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ms1_index::Ms1Scan;

    /// A run with one peptide eluting as a Gaussian at 100 s (sigma 4 s),
    /// charge 2, mono m/z 600.3, envelope [0.5, 0.3, 0.15, 0.05], plus noise
    /// peaks far away.
    fn synthetic_run() -> (Ms1RunIndex, LfqTarget) {
        let env = [0.5, 0.3, 0.15, 0.05];
        let mono = 600.3;
        let scans = (0..200)
            .map(|i| {
                let rt = i as f64; // 1 Hz
                let g = (-0.5 * ((rt - 100.0) / 4.0).powi(2)).exp();
                let mut peaks: Vec<(f64, f32)> = vec![(400.0, 50.0), (900.0, 50.0)];
                for (k, e) in env.iter().enumerate() {
                    let inten = (1e6 * e * g) as f32;
                    if inten > 1.0 {
                        peaks.push((mono + k as f64 * ISOTOPE / 2.0, inten));
                    }
                }
                peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
                Ms1Scan { rt, peaks }
            })
            .collect();
        let target = LfqTarget {
            mono_mz: mono,
            charge: 2,
            anchor_rt: 97.0,
            envelope: env.to_vec(),
        };
        (Ms1RunIndex::new(scans), target)
    }

    #[test]
    fn quantifies_the_synthetic_peak() {
        let (index, target) = synthetic_run();
        let params = LfqParams::default();
        let fq = quantify(&index, &target, &params).expect("feature");
        assert!((fq.apex_rt - 100.0).abs() <= 1.0, "{fq:?}");
        assert!(fq.rt_start < 95.0 && fq.rt_stop > 105.0, "{fq:?}");
        assert!(fq.cosine > 0.99, "{fq:?}");
        assert_eq!(fq.n_isotopes_found, 4);
        // Gaussian area = height · sigma · sqrt(2π) ≈ 0.5e6 · 4 · 2.5066 ≈ 5.01e6 for the mono.
        assert!(
            (fq.isotope_areas[0] - 5.01e6).abs() / 5.01e6 < 0.03,
            "{fq:?}"
        );
        assert!((fq.area - 1.0e7).abs() / 1.0e7 < 0.03, "{fq:?}");
        assert!((fq.rt_delta_s - 3.0).abs() <= 1.0);
    }

    #[test]
    fn decoy_twin_finds_nothing_and_scores_zero() {
        let (index, target) = synthetic_run();
        let params = LfqParams::default();
        assert!(quantify_decoy(&index, &target, &params).is_none());
        let fq = quantify(&index, &target, &params).unwrap();
        let s = feature_score(&fq, params.rt_window_s, fq.apex_intensity);
        assert!(s > 0.9, "{s}");
    }

    #[test]
    fn a_missing_monoisotope_is_not_the_target() {
        // The species that elutes sits one isotope above the target (a
        // deamidated form, say): the target's M+1.. traces hold its envelope
        // and the target's monoisotopic trace only a trace of noise (1 %).
        let (real, _) = synthetic_run();
        let target = LfqTarget {
            mono_mz: 600.3 - ISOTOPE / 2.0,
            charge: 2,
            anchor_rt: 97.0,
            envelope: vec![0.5, 0.3, 0.15, 0.05],
        };
        let scans: Vec<Ms1Scan> = (0..real.len())
            .map(|i| {
                let mut peaks = real.peaks(i).to_vec();
                let mono = peaks
                    .iter()
                    .find(|p| (p.0 - 600.3).abs() < 1e-6)
                    .map_or(0.0, |p| p.1);
                if mono > 0.0 {
                    peaks.push((target.mono_mz, mono * 0.01));
                }
                Ms1Scan {
                    rt: real.rt(i),
                    peaks,
                }
            })
            .collect();
        let index = Ms1RunIndex::new(scans);
        let lenient = LfqParams {
            min_mono_share: 0.0,
            ..LfqParams::default()
        };
        assert!(quantify(&index, &target, &lenient).is_some());
        assert!(quantify(&index, &target, &LfqParams::default()).is_none());
    }

    #[test]
    fn outside_the_run_is_none() {
        let (index, mut target) = synthetic_run();
        target.anchor_rt = 5000.0;
        assert!(quantify(&index, &target, &LfqParams::default()).is_none());
    }
}
