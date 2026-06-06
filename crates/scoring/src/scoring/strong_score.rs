//! Strong-score signal (S1) and null/competition denominator terms (S2).
//!
//! S1: context-intensity spectral similarity (`intensity_signal`).
//! S2: mass-competition evidence, listwise score gap, and candidate rank entropy.
//! Per-peak chance-match surprise lives in `match_engine::compute_psm_features`
//! as `ChanceMatchSurprise` (reused as the first null term).

use model::peptide::Peptide;

/// Half-width (Da) for `local_peak_density` in chance-match and competition terms.
pub const DENSITY_HW: f64 = 50.0;

use crate::intensity_model::{IntensityIonType, IntensityModel};
use crate::scoring::fragment_ions::{predict_by_ions, IonKind};
use crate::scoring::ScoredSpectrum;

/// Cosine similarity between two non-negative vectors (spectral-angle form).
/// Returns 0..=1; higher = better agreement. Empty or zero-norm → 0.
pub fn spectral_cosine_similarity(pred: &[f64], obs: &[f64]) -> f64 {
    if pred.is_empty() || pred.len() != obs.len() {
        return 0.0;
    }
    let dot: f64 = pred.iter().zip(obs).map(|(p, o)| p * o).sum::<f64>();
    let norm_p: f64 = pred.iter().map(|x| x * x).sum::<f64>().sqrt();
    let norm_o: f64 = obs.iter().map(|x| x * x).sum::<f64>().sqrt();
    if norm_p <= 0.0 || norm_o <= 0.0 {
        return 0.0;
    }
    (dot / (norm_p * norm_o)).clamp(0.0, 1.0)
}

/// Residues N- and C-side of the cleavage for b_idx or y_idx (1-based position).
pub fn flank_residues(seq: &[u8], kind: IonKind, idx: u32) -> Option<(u8, u8)> {
    let n = seq.len();
    let i = idx as usize;
    if i < 1 || i >= n {
        return None;
    }
    match kind {
        IonKind::B => Some((seq[i - 1], seq[i])),
        IonKind::Y => {
            let left = n - i;
            Some((seq[left - 1], seq[left]))
        }
    }
}

fn position_bin(idx: u32, pep_len: usize) -> i32 {
    if pep_len == 0 {
        return 0;
    }
    ((10.0 * idx as f64 / pep_len as f64).round() as i32).clamp(0, 10)
}

/// Signal numerator: spectral cosine between predicted and observed relative
/// intensities over charge-1 b/y ions. Returns 0.0 when `model` is `None`.
pub fn intensity_signal(
    model: Option<&IntensityModel>,
    scored_spec: &ScoredSpectrum<'_>,
    peptide: &Peptide,
    precursor_charge: u8,
    nce_bin: &str,
    feature_tol: f64,
    feature_tol_is_ppm: bool,
) -> f32 {
    let model = match model {
        Some(m) => m,
        None => return 0.0,
    };
    let n = peptide.length();
    if n < 2 {
        return 0.0;
    }

    let base_peak = scored_spec
        .dump_active_peaks()
        .iter()
        .map(|(_, _, intensity)| *intensity as f64)
        .fold(0.0_f64, f64::max);
    if base_peak <= 0.0 {
        return 0.0;
    }

    let seq: Vec<u8> = peptide.residues.iter().map(|aa| aa.residue).collect();
    let predicted = predict_by_ions(peptide, 1..=1);
    let mut pred_vec = Vec::with_capacity(predicted.len());
    let mut obs_vec = Vec::with_capacity(predicted.len());

    for ion in &predicted {
        let (flank_n, flank_c) = match flank_residues(&seq, ion.kind, ion.position) {
            Some(f) => f,
            None => continue,
        };
        let ion_type = match ion.kind {
            IonKind::B => IntensityIonType::B,
            IonKind::Y => IntensityIonType::Y,
        };
        let (mean_log, _) = model.predict_log_rel(
            ion_type,
            flank_n,
            flank_c,
            position_bin(ion.position, n),
            i32::from(precursor_charge),
            nce_bin,
        );
        pred_vec.push(mean_log.exp());

        let tol_da = if feature_tol_is_ppm {
            ion.mz * feature_tol / 1e6
        } else {
            feature_tol
        };
        let obs_rel = scored_spec
            .nearest_peak_full(ion.mz, tol_da)
            .map(|(_, intensity, _)| (f64::from(intensity) / base_peak).max(0.0))
            .unwrap_or(0.0);
        obs_vec.push(obs_rel);
    }

    spectral_cosine_similarity(&pred_vec, &obs_vec) as f32
}

/// Matched-ion tuple: (intensity, observed_mz, predicted_mz, is_b_ion).
pub type MatchedIon = (f32, f64, f64, bool);

/// S2 null term 2: `Σ 1/(1+competition)` over matched charge-1 ions.
/// `competition` = within-peptide alternative-mass ambiguity + local peak
/// density (peaks/Da) as a cheap global mass-crowding proxy.
pub fn mass_competition_evidence(
    scored_spec: &ScoredSpectrum<'_>,
    matched_ions: &[MatchedIon],
    theo_mz_list: &[f64],
    feature_tol: f64,
    feature_tol_is_ppm: bool,
) -> f32 {
    matched_ions
        .iter()
        .map(|&(_, obs, pred, _)| {
            let tol_da = if feature_tol_is_ppm {
                obs * feature_tol / 1e6
            } else {
                feature_tol
            };
            let ambiguity = theo_mz_list
                .iter()
                .filter(|&&theo| {
                    (theo - obs).abs() <= tol_da && (theo - pred).abs() > 1e-9
                })
                .count();
            let rho = scored_spec.local_peak_density(obs, DENSITY_HW);
            let competition = ambiguity as f64 + rho;
            1.0 / (1.0 + competition)
        })
        .sum::<f64>() as f32
}

/// S2 listwise term: RawScore gap between the top two retained candidates.
/// `scores_best_first` must be sorted descending.
pub fn listwise_score_gap(scores_best_first: &[f32]) -> f32 {
    if scores_best_first.len() < 2 {
        return 0.0;
    }
    (scores_best_first[0] - scores_best_first[1]).max(0.0)
}

/// S2 listwise term: Shannon entropy of a softmax over retained candidate RawScores.
/// Higher = more ambiguous top-K field. Returns 0 when fewer than two scores.
pub fn candidate_rank_entropy(scores: &[f32]) -> f32 {
    if scores.len() < 2 {
        return 0.0;
    }
    let max_s = scores.iter().copied().fold(f32::NEG_INFINITY, f32::max);
    let exp_scores: Vec<f64> = scores
        .iter()
        .map(|&s| f64::from(s - max_s).exp())
        .collect();
    let sum: f64 = exp_scores.iter().sum();
    if sum <= 0.0 {
        return 0.0;
    }
    exp_scores
        .iter()
        .map(|&e| {
            let p = e / sum;
            if p > 0.0 { -p * p.ln() } else { 0.0 }
        })
        .sum::<f64>() as f32
}

/// Inputs for [`fuse_strong_score`] (mirrors the S1/S2 PIN feature columns).
#[derive(Debug, Clone, Copy, Default)]
pub struct StrongScoreInputs {
    pub intensity_signal: f32,
    pub chance_match_surprise: f32,
    pub mass_competition_evidence: f32,
    pub candidate_rank_entropy: f32,
    pub listwise_score_gap: f32,
}

/// Fuse S1 signal with S2 null terms: `strong = signal − null`.
///
/// Null is the log-domain coincidental-match cost (higher ⇒ more coincidental):
/// - chance: `−chance_match_surprise` (high surprise ⇒ low coincidence cost)
/// - competition: `−mass_competition_evidence`
/// - listwise: `candidate_rank_entropy − listwise_score_gap`
pub fn fuse_strong_score(f: &StrongScoreInputs) -> f32 {
    let null = -f.chance_match_surprise
        - f.mass_competition_evidence
        + f.candidate_rank_entropy
        - f.listwise_score_gap;
    f.intensity_signal - null
}

/// Minimum scored candidates before per-spectrum z-score calibration (mirrors Tailor).
pub const STRONG_CAL_MIN_CANDIDATES: u32 = 100;

/// Welford online mean/variance for per-spectrum null statistics.
#[derive(Debug, Clone, Copy, Default)]
pub struct OnlineStats {
    n: u64,
    mean: f64,
    m2: f64,
}

impl OnlineStats {
    pub fn push(&mut self, x: f32) {
        let x = f64::from(x);
        self.n += 1;
        let delta = x - self.mean;
        self.mean += delta / self.n as f64;
        let delta2 = x - self.mean;
        self.m2 += delta * delta2;
    }

    pub fn count(&self) -> u64 {
        self.n
    }

    pub fn mean(&self) -> f64 {
        self.mean
    }

    pub fn population_stdev(&self) -> f64 {
        if self.n < 2 {
            return 1.0;
        }
        (self.m2 / self.n as f64).sqrt().max(1e-6)
    }
}

/// Z-score `score` against a per-spectrum null pool; returns `score` uncalibrated
/// when `n < STRONG_CAL_MIN_CANDIDATES`.
pub fn strong_score_zscore(score: f32, null: &OnlineStats) -> f32 {
    if null.n < STRONG_CAL_MIN_CANDIDATES as u64 {
        return score;
    }
    ((f64::from(score) - null.mean()) / null.population_stdev()) as f32
}

/// Leave-one-out z-score of `this` against other retained `strong_score` values.
/// Used when top-N retains multiple candidates; falls back to `this` when alone.
pub fn strong_score_calibrated_loo(retained_strong: &[f32], this: f32) -> f32 {
    if retained_strong.len() < 2 {
        return this;
    }
    let n = retained_strong.len();
    let this_d = f64::from(this);
    let sum: f64 = retained_strong.iter().map(|&s| f64::from(s)).sum();
    let sum_sq: f64 = retained_strong.iter().map(|&s| f64::from(s) * f64::from(s)).sum();
    // Leave exactly ONE copy of `this` out, so mean and variance are over the
    // same subset. The prior `filter(s != this)` dropped EVERY tied value,
    // making the variance disagree with the mean when ties were present.
    let n_others = (n - 1) as f64;
    let mean_others = (sum - this_d) / n_others;
    let var_others = ((sum_sq - this_d * this_d) / n_others - mean_others * mean_others).max(0.0);
    let sigma = var_others.sqrt().max(1e-6);
    ((this_d - mean_others) / sigma) as f32
}

/// S4 calibration: prefer LOO among retained strong scores when top-N ≥ 2; otherwise
/// z-score against the per-spectrum scored-candidate null (`pin_score` pool).
pub fn strong_score_calibrated(
    strong: f32,
    retained_strong: &[f32],
    pin_null: &OnlineStats,
) -> f32 {
    if retained_strong.len() >= 2 {
        strong_score_calibrated_loo(retained_strong, strong)
    } else {
        strong_score_zscore(strong, pin_null)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::intensity_model::IntensityModel;
    use arrow::array::{Float64Array, Int32Array, Int64Array, StringArray};
    use arrow::datatypes::{DataType, Field, Schema};
    use arrow::record_batch::RecordBatch;
    use model::amino_acid::AminoAcid;
    use model::peptide::Peptide;
    use model::spectrum::Spectrum;
    use parquet::arrow::ArrowWriter;
    use std::fs::File;
    use std::path::Path;
    use tempfile::NamedTempFile;

    fn write_model(path: &Path, rows: &[(&str, &str, &str, i32, i32, &str, i64, f64, f64)]) {
        let schema = Schema::new(vec![
            Field::new("ion_type", DataType::Utf8, false),
            Field::new("flank_n", DataType::Utf8, false),
            Field::new("flank_c", DataType::Utf8, false),
            Field::new("pos_bin", DataType::Int32, false),
            Field::new("charge", DataType::Int32, false),
            Field::new("nce_bin", DataType::Utf8, false),
            Field::new("count", DataType::Int64, false),
            Field::new("mean_log_rel", DataType::Float64, false),
            Field::new("var_log_rel", DataType::Float64, false),
        ]);
        let ion: Vec<_> = rows.iter().map(|r| r.0).collect();
        let fn_: Vec<_> = rows.iter().map(|r| r.1).collect();
        let fc: Vec<_> = rows.iter().map(|r| r.2).collect();
        let pb: Vec<_> = rows.iter().map(|r| r.3).collect();
        let ch: Vec<_> = rows.iter().map(|r| r.4).collect();
        let nce: Vec<_> = rows.iter().map(|r| r.5).collect();
        let cnt: Vec<_> = rows.iter().map(|r| r.6).collect();
        let mean: Vec<_> = rows.iter().map(|r| r.7).collect();
        let var: Vec<_> = rows.iter().map(|r| r.8).collect();
        let batch = RecordBatch::try_new(
            std::sync::Arc::new(schema),
            vec![
                std::sync::Arc::new(StringArray::from(ion)),
                std::sync::Arc::new(StringArray::from(fn_)),
                std::sync::Arc::new(StringArray::from(fc)),
                std::sync::Arc::new(Int32Array::from(pb)),
                std::sync::Arc::new(Int32Array::from(ch)),
                std::sync::Arc::new(StringArray::from(nce)),
                std::sync::Arc::new(Int64Array::from(cnt)),
                std::sync::Arc::new(Float64Array::from(mean)),
                std::sync::Arc::new(Float64Array::from(var)),
            ],
        )
        .unwrap();
        let file = File::create(path).unwrap();
        let mut writer = ArrowWriter::try_new(file, batch.schema(), None).unwrap();
        writer.write(&batch).unwrap();
        writer.close().unwrap();
    }

    fn pep(seq: &[u8]) -> Peptide {
        let residues: Vec<AminoAcid> = seq
            .iter()
            .map(|&r| AminoAcid::standard(r).unwrap())
            .collect();
        Peptide::new(residues, b'K', b'R')
    }

    #[test]
    fn spectral_cosine_perfect_and_orthogonal() {
        assert!((spectral_cosine_similarity(&[1.0, 2.0], &[1.0, 2.0]) - 1.0).abs() < 1e-9);
        assert!(spectral_cosine_similarity(&[1.0, 0.0], &[0.0, 1.0]) < 1e-9);
        assert_eq!(spectral_cosine_similarity(&[], &[]), 0.0);
    }

    #[test]
    fn flank_residues_b_and_y() {
        let seq = b"ARCDE";
        assert_eq!(flank_residues(seq, IonKind::B, 2), Some((b'R', b'C')));
        assert_eq!(flank_residues(seq, IonKind::Y, 2), Some((b'C', b'D')));
    }

    #[test]
    fn intensity_signal_zero_without_model() {
        let spec = Spectrum {
            peaks: vec![(500.0, 1000.0)],
            ..Default::default()
        };
        let ss = ScoredSpectrum::new_without_filtering(&spec);
        let signal = intensity_signal(None, &ss, &pep(b"AR"), 2, "unknown", 20.0, true);
        assert_eq!(signal, 0.0);
    }

    #[test]
    fn intensity_signal_higher_when_observed_matches_bright_prediction() {
        let tmp = NamedTempFile::new().unwrap();
        write_model(
            tmp.path(),
            &[
                ("y", "R", "C", 5, 2, "unknown", 100, -0.2, 0.1),
                ("b", "A", "R", 1, 2, "unknown", 100, -2.5, 0.1),
            ],
        );
        let model = IntensityModel::load(tmp.path()).unwrap();
        let peptide = pep(b"ARCDE");
        let predicted = predict_by_ions(&peptide, 1..=1);
        let y3 = predicted
            .iter()
            .find(|p| p.kind == IonKind::Y && p.position == 3)
            .expect("y3");
        let mz = y3.mz;
        let spec = Spectrum {
            peaks: vec![(200.0, 10.0), (mz, 1000.0), (mz + 0.001, 50.0)],
            ..Default::default()
        };
        let ss = ScoredSpectrum::new_without_filtering(&spec);
        let good = intensity_signal(Some(&model), &ss, &peptide, 2, "unknown", 20.0, true);

        let wrong_pep = pep(b"FGHIK");
        let bad = intensity_signal(Some(&model), &ss, &wrong_pep, 2, "unknown", 20.0, true);
        assert!(good > bad, "good={good} bad={bad}");
        assert!(good > 0.1);
    }

    #[test]
    fn mass_competition_lower_in_crowded_region() {
        let peptide = pep(b"ARCDE");
        let predicted = predict_by_ions(&peptide, 1..=1);
        let y3 = predicted
            .iter()
            .find(|p| p.kind == IonKind::Y && p.position == 3)
            .unwrap();
        let sparse = Spectrum {
            peaks: vec![(y3.mz, 1000.0)],
            ..Default::default()
        };
        let crowded = Spectrum {
            peaks: (0..20)
                .map(|i| (y3.mz - 0.5 + i as f64 * 0.05, 50.0 + i as f32))
                .collect(),
            ..Default::default()
        };
        let ss_sparse = ScoredSpectrum::new_without_filtering(&sparse);
        let ss_crowded = ScoredSpectrum::new_without_filtering(&crowded);
        let mut theo: Vec<f64> = Vec::new();
        for p in &predicted {
            theo.push(p.mz);
        }
        let matched = vec![(1000.0_f32, y3.mz, y3.mz, false)];
        let sparse_ev = mass_competition_evidence(&ss_sparse, &matched, &theo, 20.0, true);
        let crowded_ev = mass_competition_evidence(&ss_crowded, &matched, &theo, 20.0, true);
        assert!(sparse_ev > crowded_ev);
    }

    #[test]
    fn candidate_rank_entropy_uniform_high_dominant_low() {
        let uniform = candidate_rank_entropy(&[10.0, 10.0, 10.0]);
        let dominant = candidate_rank_entropy(&[100.0, 1.0, 1.0]);
        assert!(uniform > dominant);
        assert_eq!(candidate_rank_entropy(&[5.0]), 0.0);
    }

    #[test]
    fn strong_score_zscore_and_loo() {
        let mut null = OnlineStats::default();
        for v in [10.0_f32, 12.0, 11.0, 13.0] {
            null.push(v);
        }
        // Below minimum candidate count → passthrough.
        assert_eq!(strong_score_zscore(5.0, &null), 5.0);
        for _ in 0..100 {
            null.push(10.0);
        }
        let z = strong_score_zscore(15.0, &null);
        assert!(z > 0.0);
        let loo = strong_score_calibrated_loo(&[1.0, 3.0, 5.0], 5.0);
        assert!(loo > 0.0);
        assert_eq!(strong_score_calibrated_loo(&[2.0], 2.0), 2.0);
    }

    #[test]
    fn fuse_strong_score_increases_with_surprise_and_evidence() {
        let base = StrongScoreInputs {
            intensity_signal: 0.5,
            chance_match_surprise: 0.0,
            mass_competition_evidence: 0.0,
            candidate_rank_entropy: 0.0,
            listwise_score_gap: 0.0,
        };
        let low = fuse_strong_score(&base);
        let high = fuse_strong_score(&StrongScoreInputs {
            chance_match_surprise: 3.0,
            mass_competition_evidence: 2.0,
            listwise_score_gap: 1.0,
            ..base
        });
        assert!(high > low);
        let ambiguous = fuse_strong_score(&StrongScoreInputs {
            candidate_rank_entropy: 2.0,
            ..base
        });
        assert!(ambiguous < low);
    }

    #[test]
    fn listwise_score_gap_basic() {
        assert_eq!(listwise_score_gap(&[10.0, 7.0, 3.0]), 3.0);
        assert_eq!(listwise_score_gap(&[5.0]), 0.0);
    }

    #[test]
    fn missing_observed_ions_reduce_signal() {
        let tmp = NamedTempFile::new().unwrap();
        write_model(
            tmp.path(),
            &[("y", "R", "C", 5, 2, "unknown", 100, 0.0, 0.01)],
        );
        let model = IntensityModel::load(tmp.path()).unwrap();
        let peptide = pep(b"ARCDE");
        let predicted = predict_by_ions(&peptide, 1..=1);
        let y3 = predicted
            .iter()
            .find(|p| p.kind == IonKind::Y && p.position == 3)
            .unwrap();
        let spec = Spectrum {
            peaks: vec![(y3.mz, 1000.0)],
            ..Default::default()
        };
        let ss = ScoredSpectrum::new_without_filtering(&spec);
        let partial = intensity_signal(Some(&model), &ss, &peptide, 2, "unknown", 20.0, true);
        let empty_spec = Spectrum {
            peaks: vec![(400.0, 1000.0)],
            ..Default::default()
        };
        let ss_empty = ScoredSpectrum::new_without_filtering(&empty_spec);
        let none = intensity_signal(Some(&model), &ss_empty, &peptide, 2, "unknown", 20.0, true);
        assert!(partial > none);
    }
}
