//! Per-run fragment mass error, measured on confident PSMs of a sample of the
//! run's spectra.
//!
//! A high-resolution analyzer's fragment error differs by instrument and by run:
//! an Astral run measures about -2.5 ppm with a spread of under 2 ppm, an
//! Orbitrap run can sit +4 ppm off, and TOF analyzers spread wider. The search
//! scores a sample once, keeps the target PSMs at 1% target-decoy q-value,
//! re-matches their most intense b/y ions in a wide window and takes a robust
//! centre and spread of the errors. [`scoring_crate::scoring::FragmentErrorModel`]
//! then sets the PIN feature tolerance and the ppm-Gaussian term from them.

use model::peptide::Peptide;
use model::spectrum::Spectrum;
use scoring_crate::scoring::fragment_ions::predict_by_ions;
use scoring_crate::scoring::FragmentErrorModel;

use crate::match_engine::PreparedSearch;
use crate::tdc::{confident_target_indices, ScoredLabel};

/// Half-width of the window the sampled PSMs' ions are re-matched in (ppm).
/// Wider than any supported analyzer's fragment error, so the measurement is
/// not cut off by the default feature window.
pub const MEASURE_WINDOW_PPM: f64 = 50.0;
/// Most intense matched ions taken per PSM.
pub const TOP_IONS_PER_PSM: usize = 7;
/// Confident PSMs needed for a measurement.
pub const MIN_CONFIDENT_PSMS: usize = 100;
/// Target-decoy q-value of the PSMs measured.
pub const MEASURE_Q: f64 = 0.01;
/// Spectra scored for the measurement.
pub const SAMPLE_SPECTRA: usize = 3000;

/// Signed ppm errors of the `k` most intense b/y ions (charge 1) of `peptide`
/// found within [`MEASURE_WINDOW_PPM`] in `spec` (peaks sorted by m/z).
pub fn top_ion_errors_ppm(spec: &Spectrum, peptide: &Peptide, k: usize) -> Vec<f64> {
    let peaks = &spec.peaks;
    let mut matched: Vec<(f32, f64)> = Vec::new();
    for ion in predict_by_ions(peptide, 1..=1) {
        let tol = ion.mz * MEASURE_WINDOW_PPM * 1e-6;
        let lo = peaks.partition_point(|p| p.0 < ion.mz - tol);
        let best = peaks[lo..]
            .iter()
            .take_while(|p| p.0 <= ion.mz + tol)
            .min_by(|a, b| (a.0 - ion.mz).abs().total_cmp(&(b.0 - ion.mz).abs()));
        if let Some(&(mz, intensity)) = best {
            matched.push((intensity, (mz - ion.mz) / ion.mz * 1e6));
        }
    }
    matched.sort_by(|a, b| b.0.total_cmp(&a.0));
    matched.truncate(k);
    matched.into_iter().map(|(_, e)| e).collect()
}

fn median(v: &mut [f64]) -> f64 {
    v.sort_by(f64::total_cmp);
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2]
    } else {
        0.5 * (v[n / 2 - 1] + v[n / 2])
    }
}

/// Median and normal-consistent MAD of `errors`.
fn centre_and_spread(errors: &[f64]) -> (f64, f64) {
    let mut v = errors.to_vec();
    let c = median(&mut v);
    let mut dev: Vec<f64> = errors.iter().map(|e| (e - c).abs()).collect();
    (c, 1.4826 * median(&mut dev))
}

/// Robust centre and spread of pooled ion errors. One trimming pass at four
/// spreads drops chance matches from the wide window before the final estimate.
pub fn estimate(errors: &[f64], psms: usize) -> Option<FragmentErrorModel> {
    if errors.len() < MIN_CONFIDENT_PSMS {
        return None;
    }
    let (c0, s0) = centre_and_spread(errors);
    let half = (4.0 * s0).max(1.0);
    let kept: Vec<f64> = errors
        .iter()
        .copied()
        .filter(|e| (e - c0).abs() <= half)
        .collect();
    if kept.len() < MIN_CONFIDENT_PSMS {
        return None;
    }
    let (offset_ppm, sigma_ppm) = centre_and_spread(&kept);
    Some(FragmentErrorModel::new(
        offset_ppm,
        sigma_ppm,
        kept.len(),
        psms,
    ))
}

/// Why a measurement was not made.
#[derive(Debug, Clone, PartialEq)]
pub enum NotMeasured {
    /// Fewer than [`MIN_CONFIDENT_PSMS`] target PSMs at [`MEASURE_Q`].
    TooFewPsms { confident: usize, scored: usize },
    /// Too few matched ions survived the robust trim.
    TooFewIons { confident: usize },
}

impl std::fmt::Display for NotMeasured {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            NotMeasured::TooFewPsms { confident, scored } => write!(
                f,
                "{confident} of {scored} sampled spectra have a target PSM at {MEASURE_Q} q-value \
                 (need {MIN_CONFIDENT_PSMS})"
            ),
            NotMeasured::TooFewIons { confident } => {
                write!(f, "too few matched ions in {confident} confident PSMs")
            }
        }
    }
}

/// Score `sample`, keep the rank-1 target PSMs at [`MEASURE_Q`] and measure
/// their fragment error. The sample must suit `prepared`'s retrieval (one
/// fragment-index slice when the index is on).
pub fn measure_on_queues(
    prepared: &PreparedSearch<'_>,
    sample: &[Spectrum],
    queues: &[crate::TopNQueue],
) -> Result<FragmentErrorModel, NotMeasured> {
    let mut tops: Vec<(usize, u32)> = Vec::new();
    let mut labels: Vec<ScoredLabel> = Vec::new();
    for (i, q) in queues.iter().enumerate() {
        if let Some(p) = q.peek_top() {
            let idx = p.primary_candidate_idx();
            tops.push((i, idx));
            labels.push(ScoredLabel {
                score: p.rank_score,
                is_decoy: prepared.candidate_is_decoy(idx),
            });
        }
    }
    let confident = confident_target_indices(&labels, MEASURE_Q);
    if confident.len() < MIN_CONFIDENT_PSMS {
        return Err(NotMeasured::TooFewPsms {
            confident: confident.len(),
            scored: labels.len(),
        });
    }
    let mut errors: Vec<f64> = Vec::new();
    for &k in &confident {
        let (si, idx) = tops[k];
        let Some(spec) = sample.get(si) else {
            continue;
        };
        let peptide = prepared.candidate_peptide(idx);
        errors.extend(top_ion_errors_ppm(spec, &peptide, TOP_IONS_PER_PSM));
    }
    estimate(&errors, confident.len()).ok_or(NotMeasured::TooFewIons {
        confident: confident.len(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::AminoAcid;

    fn peptide(residues: &[u8]) -> Peptide {
        let aas: Vec<AminoAcid> = residues
            .iter()
            .map(|&r| AminoAcid::standard(r).unwrap())
            .collect();
        Peptide::new(aas, b'K', b'K')
    }

    /// The most intense matched ions are taken, each at its signed ppm error;
    /// ions outside the 50 ppm window are not matched.
    #[test]
    fn top_ion_errors_take_the_most_intense_matches() {
        let p = peptide(b"PEPTIDEK");
        let ions = predict_by_ions(&p, 1..=1);
        let mut peaks: Vec<(f64, f32)> = ions
            .iter()
            .enumerate()
            .map(|(i, ion)| (ion.mz * (1.0 + 3.0e-6), 100.0 + i as f32))
            .collect();
        // One ion 80 ppm off: outside the window.
        peaks[0].0 = ions[0].mz * (1.0 + 80.0e-6);
        peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
        let spec = Spectrum {
            peaks,
            ..Default::default()
        };
        let errs = top_ion_errors_ppm(&spec, &p, 7);
        assert_eq!(errs.len(), 7);
        assert!(errs.iter().all(|e| (e - 3.0).abs() < 1e-6), "{errs:?}");
        let all = top_ion_errors_ppm(&spec, &p, 100);
        assert_eq!(all.len(), ions.len() - 1);
    }

    /// Errors around -2.5 ppm (spread 1.6) with 10% chance matches spread over
    /// the whole window: the trimmed estimate recovers the centre and spread.
    #[test]
    fn estimate_is_robust_to_chance_matches() {
        let mut state = 12345u64;
        let mut uniform = || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        let mut errs = Vec::new();
        for i in 0..5000 {
            if i % 10 == 0 {
                errs.push(-50.0 + 100.0 * uniform());
            } else {
                let (u1, u2) = (uniform().max(1e-12), uniform());
                let n = (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();
                errs.push(-2.5 + 1.6 * n);
            }
        }
        let m = estimate(&errs, 700).expect("estimate");
        assert!((m.offset_ppm + 2.5).abs() < 0.15, "{m:?}");
        assert!((m.sigma_ppm - 1.6).abs() < 0.25, "{m:?}");
        assert_eq!(m.psms, 700);
        assert!(m.feature_tolerance_ppm() > 7.0 && m.feature_tolerance_ppm() < 11.0);
    }

    #[test]
    fn too_few_errors_give_no_estimate() {
        assert_eq!(estimate(&[1.0; 10], 10), None);
    }
}
