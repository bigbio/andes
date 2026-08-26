//! Strong-score signal (S1) and null/competition denominator terms (S2).
//!
//! S1: context-intensity spectral similarity (`intensity_signal`).
//! S2: mass-competition evidence, listwise score gap, and candidate rank entropy.
//! Per-peak chance-match surprise lives in `match_engine::compute_psm_features`
//! as `ChanceMatchSurprise` (reused as the first null term).

use model::peptide::Peptide;

/// Half-width (Da) for `local_peak_density` in chance-match and competition terms.
pub const DENSITY_HW: f64 = 50.0;

/// Finite clamp (finding 3.4) for a GBDT regressor's predicted log-relative
/// fragment intensity before `exp`. The frag-intensity model emits a *log*
/// relative intensity; a corrupt/degenerate tree (extreme leaf values) can
/// return a magnitude large enough that `exp` overflows to `+inf`, which then
/// poisons every downstream sum (`sum_p`, `explained`, the chance LLR). Real
/// calibrated outputs sit within a few units of zero, so this bound is never
/// engaged on the happy path — `exp(±50)` (≈5e21 / ≈2e-22) is already far past
/// any physical relative intensity. Clamping keeps the derived features finite
/// instead of `inf`/`NaN`.
const LOG_INTENSITY_CLAMP: f64 = 50.0;

/// Convert a GBDT regressor's predicted log-relative intensity to a finite,
/// non-negative linear intensity. Clamps the log value to
/// `[-LOG_INTENSITY_CLAMP, LOG_INTENSITY_CLAMP]` (a no-op for real outputs) and
/// substitutes `0.0` for a non-finite prediction so callers never propagate
/// `inf`/`NaN`.
/// Predict the linear relative intensity of every `predict_by_ions(peptide, 1..=2)`
/// fragment ONCE, in that iteration order.
///
/// `intensity_signal` and `frag_llr_battery` are called back-to-back on the same
/// `(peptide, precursor_charge, frag_model)` and each used to walk the identical ion
/// list, extract byte-identical features, and run the identical 300-tree GBDT — so the
/// single most expensive thing in the search (`Tree::eval`, ~50% of wall) was being done
/// twice per candidate for the same numbers. Callers compute this once and hand the slice
/// to both. Values are positionally aligned with `predict_by_ions(peptide, 1..=2)`; a
/// caller passing a slice built any other way will silently mis-attribute intensities.
pub fn predict_frag_intensities(
    frag_model: &GbdtPeakModel,
    peptide: &Peptide,
    precursor_charge: u8,
) -> Vec<f64> {
    let ions = predict_by_ions(peptide, 1..=2);
    // Build every feature row first, then walk the ensemble ONCE with trees on the
    // outside (see `GbdtPeakModel::predict_value_batch`). Bit-identical to the
    // per-ion call it replaces.
    let feats: Vec<[f32; crate::frag_features::N_FRAG_FEATURES]> = ions
        .iter()
        .map(|ion| {
            extract_frag_features(
                peptide,
                ion.kind,
                ion.position,
                precursor_charge,
                ion.charge,
                0.0,
            )
        })
        .collect();
    let rows: Vec<&[f32]> = feats.iter().map(|f| f.as_slice()).collect();
    let mut raw = vec![0.0f32; rows.len()];
    frag_model.predict_value_batch(&rows, &mut raw);
    raw.iter()
        .map(|&log_rel| {
            let v = f64::from(log_rel);
            if !v.is_finite() {
                return 0.0;
            }
            v.clamp(-LOG_INTENSITY_CLAMP, LOG_INTENSITY_CLAMP).exp()
        })
        .collect()
}

#[inline]
fn predicted_linear_intensity(g: &GbdtPeakModel, feats: &[f32]) -> f64 {
    let log_rel = f64::from(g.predict_value(feats));
    if !log_rel.is_finite() {
        return 0.0;
    }
    log_rel.clamp(-LOG_INTENSITY_CLAMP, LOG_INTENSITY_CLAMP).exp()
}

use crate::frag_features::extract_frag_features;
use crate::ion_features::extract_ion_features;
use crate::gbdt_eval::GbdtPeakModel;
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
        IonKind::B | IonKind::C => Some((seq[i - 1], seq[i])),
        IonKind::Y | IonKind::Z => {
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
/// intensities over charge-1 b/y ions.
///
/// - If `frag_model` is `Some`, it drives the predicted intensities (charge=1,
///   nce=0.0 to match training in `frag_dataset.rs`).
/// - Otherwise falls back to the coarse `IntensityModel` lookup table
///   (existing behaviour, byte-identical when `frag_model` is `None`).
/// - Returns 0.0 when both `model` and `frag_model` are `None` (no predictor).
#[allow(clippy::too_many_arguments)]
pub fn intensity_signal(
    model: Option<&IntensityModel>,
    frag_model: Option<&GbdtPeakModel>,
    scored_spec: &ScoredSpectrum<'_>,
    peptide: &Peptide,
    precursor_charge: u8,
    nce_bin: &str,
    feature_tol: f64,
    feature_tol_is_ppm: bool,
    // Predictions from `predict_frag_intensities`, positionally aligned with
    // `predict_by_ions(peptide, 1..=2)`. `None` recomputes them inline.
    precomputed: Option<&[f64]>,
) -> f32 {
    // Require at least one predictor.
    if model.is_none() && frag_model.is_none() {
        return 0.0;
    }
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
    // Predict 1+ AND 2+ fragments: high-charge precursors put real signal in 2+
    // ions, which a 1+-only prediction is blind to. Training (frag_dataset) uses
    // the SAME 1..=2 range and the same frag-charge feature.
    let predicted = predict_by_ions(peptide, 1..=2);
    let mut pred_vec = Vec::with_capacity(predicted.len());
    let mut obs_vec = Vec::with_capacity(predicted.len());

    for (ion_idx, ion) in predicted.iter().enumerate() {
        // Predicted LINEAR relative intensity, kept finite (finding 3.4): the
        // GBDT/intensity-model paths emit a log-relative value; an extreme
        // prediction would otherwise `exp` to `+inf` and corrupt the cosine.
        let pred_intensity = if let Some(g) = frag_model {
            // Frag-intensity regressor path (precursor charge, nce=0.0 matches training).
            // Reuse the caller's predictions when supplied — same ion order, same
            // features, same model, so this is the identical value without a second
            // 300-tree walk.
            match precomputed.and_then(|p| p.get(ion_idx).copied()) {
                Some(v) => v,
                None => {
                    let feats = extract_frag_features(peptide, ion.kind, ion.position, precursor_charge, ion.charge, 0.0);
                    predicted_linear_intensity(g, &feats)
                }
            }
        } else {
            // Existing coarse table path (fallback; model is Some here).
            let (flank_n, flank_c) = match flank_residues(&seq, ion.kind, ion.position) {
                Some(f) => f,
                None => continue,
            };
            let ion_type = match ion.kind {
                IonKind::B | IonKind::C => IntensityIonType::B,
                IonKind::Y | IonKind::Z => IntensityIonType::Y,
            };
            // model is guaranteed Some when frag_model is None (checked above).
            #[allow(clippy::unwrap_used)]
            let (mean_log, _) = model.unwrap().predict_log_rel(
                ion_type,
                flank_n,
                flank_c,
                position_bin(ion.position, n),
                i32::from(precursor_charge),
                nce_bin,
            );
            // Preserve the historical f64→f32→f64 round-trip exactly so the
            // happy-path cosine is byte-identical; only add finiteness/clamp.
            let log_rel = f64::from(mean_log as f32);
            if log_rel.is_finite() {
                log_rel.clamp(-LOG_INTENSITY_CLAMP, LOG_INTENSITY_CLAMP).exp()
            } else {
                0.0
            }
        };
        pred_vec.push(pred_intensity);

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

/// Frag-intensity LLR battery: three additive discriminative PIN features
/// derived from the frag-intensity GBDT's per-fragment predicted intensity,
/// deployed as likelihood-ratio signals instead of a single cosine (the cosine
/// normalizes both vectors and smears target/decoy separation, so a more
/// accurate intensity model gives no PSM lift — confirmed 3× on Astral).
///
/// Returns `(explained, chance_llr, topk_observed)`:
/// - `explained` = `Σ(matched·pred) / Σpred` — of the intensity the model
///   predicts, the fraction actually observed. Asymmetric (NOT normalized by the
///   observed vector like cosine) → true peptides explain their predicted-bright
///   ions; decoys don't.
/// - `chance_llr` = `Σ matched·pred·max(0, −ln p_chance)` — predicted-intensity-
///   weighted chance-match surprise. Fuses the prediction with the local-noise
///   denominator (rewards a predicted-bright match in a sparse region); the
///   cosine ignores chance entirely.
/// - `topk_observed` = fraction of the top-`FRAG_TOPK` predicted-most-intense
///   ions that are observed — rank-based, intensity-scale-free, robust.
///
/// All `0.0` when `frag_model` is `None` (neutral additive features). Enumerates
/// `1..=2` to match the serve cosine + the trainer.
pub fn frag_llr_battery(
    frag_model: Option<&GbdtPeakModel>,
    scored_spec: &ScoredSpectrum<'_>,
    peptide: &Peptide,
    precursor_charge: u8,
    feature_tol: f64,
    feature_tol_is_ppm: bool,
    // Predictions from `predict_frag_intensities`, positionally aligned with
    // `predict_by_ions(peptide, 1..=2)`. `None` recomputes them inline.
    precomputed: Option<&[f64]>,
) -> (f32, f32, f32) {
    const FRAG_TOPK: usize = 6;
    let g = match frag_model {
        Some(g) => g,
        None => return (0.0, 0.0, 0.0),
    };
    if peptide.length() < 2 {
        return (0.0, 0.0, 0.0);
    }
    let predicted = predict_by_ions(peptide, 1..=2);
    let mut sum_p = 0.0f64;
    let mut sum_matched_p = 0.0f64;
    let mut sum_chance_llr = 0.0f64;
    // (predicted intensity, matched?) for the top-K observed-fraction feature.
    let mut pred_matched: Vec<(f64, bool)> = Vec::with_capacity(predicted.len());

    for (ion_idx, ion) in predicted.iter().enumerate() {
        // Reuse the caller's predictions when supplied (see `predict_frag_intensities`).
        let p = match precomputed.and_then(|pp| pp.get(ion_idx).copied()) {
            Some(v) => v,
            None => {
                let feats = extract_frag_features(
                    peptide,
                    ion.kind,
                    ion.position,
                    precursor_charge,
                    ion.charge,
                    0.0,
                );
                predicted_linear_intensity(g, &feats)
            }
        };
        let tol_da = if feature_tol_is_ppm {
            ion.mz * feature_tol / 1e6
        } else {
            feature_tol
        };
        let matched = scored_spec.nearest_peak_full(ion.mz, tol_da).is_some();
        sum_p += p;
        if matched {
            sum_matched_p += p;
            let rho = scored_spec.local_peak_density(ion.mz, DENSITY_HW);
            // p_chance = probability of a random peak in the match window.
            let p_chance = (rho * 2.0 * tol_da).clamp(1e-12, 1.0);
            let surprise = (-p_chance.ln()).max(0.0);
            sum_chance_llr += p * surprise;
        }
        pred_matched.push((p, matched));
    }

    let explained = if sum_p > 0.0 {
        (sum_matched_p / sum_p) as f32
    } else {
        0.0
    };

    // Fraction of the top-K predicted-most-intense ions that were observed.
    pred_matched.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
    let k = pred_matched.len().min(FRAG_TOPK);
    let topk_observed = if k > 0 {
        pred_matched[..k].iter().filter(|(_, m)| *m).count() as f32 / k as f32
    } else {
        0.0
    };

    (explained, sum_chance_llr as f32, topk_observed)
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
    // Leave exactly ONE copy of `this` out, so mean and variance are over the
    // same subset.
    let n_others = (n - 1) as f64;
    let mean_others = (sum - this_d) / n_others;
    // Two-pass variance about the mean. The previous `E[x^2] - E[x]^2` form
    // cancels catastrophically when the retained scores are large and tightly
    // clustered: both terms are big and nearly equal, so the difference keeps
    // only the low-order bits and can land at or below zero. Measured on a real
    // 55,171-row PIN: 53 rows blew past |1e6| and inflated the column's stdev
    // 72x over the inliers, which after Percolator's standardisation leaves the
    // feature carrying almost nothing but those spikes.
    let mut ss = 0.0f64;
    let mut skipped_self = false;
    for &v in retained_strong {
        let d = f64::from(v);
        if !skipped_self && d == this_d {
            skipped_self = true;
            continue;
        }
        let dev = d - mean_others;
        ss += dev * dev;
    }
    let var_others = (ss / n_others).max(0.0);
    // A degenerate spectrum - every retained score effectively identical - carries
    // NO discriminative information about this candidate, and the honest encoding
    // of that is 0.0, not a huge z-score. The previous code floored sigma at an
    // absolute 1e-6 and divided a non-zero numerator by it, manufacturing values
    // in the millions out of nothing. Judge degeneracy RELATIVE to the score
    // scale so the test is not itself scale-blind.
    let scale = mean_others.abs().max(1.0);
    let sigma = var_others.sqrt();
    if sigma <= 1e-9 * scale {
        return 0.0;
    }
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

/// Decoy-aware rich-ion LLR: `Σ` over matched b/y ions of the per-ion
/// `predict_logit(P(signal | features))` from the rich-ion GBDT. The model is
/// trained on target-matched (positive) vs decoy-matched (negative) ions, so a
/// true peptide's matched ions score high and a decoy's chance-matched ions
/// score low — the SUM separates target from decoy (assignment-aware, NOT the
/// isotropic per-peak quality that sank the intensity cosine). Returns 0.0 when
/// no rich-ion model is loaded (neutral additive feature; RankScore untouched).
/// Enumerates `1..=2` to match the trainer (`ion_dataset`).
pub fn rich_ion_llr(
    model: Option<&GbdtPeakModel>,
    scored_spec: &ScoredSpectrum<'_>,
    peptide: &Peptide,
    precursor_charge: u8,
    feature_tol: f64,
    feature_tol_is_ppm: bool,
) -> f32 {
    let g = match model {
        Some(g) => g,
        None => return 0.0,
    };
    if peptide.length() < 2 {
        return 0.0;
    }
    // Collect the feature rows for MATCHED ions first, then walk the rich-ion
    // ensemble once with trees on the outside. Bit-identical: the rows are gathered
    // in the same ion order the per-ion loop used, and the f64 sum below adds the
    // per-row logits in that same order.
    let mut rows_buf: Vec<Vec<f32>> = Vec::new();
    for ion in predict_by_ions(peptide, 1..=2) {
        let tol_da = if feature_tol_is_ppm {
            ion.mz * feature_tol / 1e6
        } else {
            feature_tol
        };
        if scored_spec.nearest_peak_full(ion.mz, tol_da).is_some() {
            rows_buf.push(extract_ion_features(
                peptide,
                scored_spec,
                ion.kind,
                ion.position,
                precursor_charge,
                ion.charge,
                feature_tol,
                feature_tol_is_ppm,
            ).to_vec());
        }
    }
    if rows_buf.is_empty() {
        return 0.0;
    }
    let rows: Vec<&[f32]> = rows_buf.iter().map(|r| r.as_slice()).collect();
    let mut logits = vec![0.0f32; rows.len()];
    g.predict_logit_batch(&rows, &mut logits);
    let mut sum = 0.0f64;
    for l in &logits {
        sum += f64::from(*l);
    }
    sum as f32
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

    /// (ion_type, flank_n, flank_c, pos_bin, charge, nce_bin, count, mean, var)
    type ModelRow<'a> = (&'a str, &'a str, &'a str, i32, i32, &'a str, i64, f64, f64);

    fn write_model(path: &Path, rows: &[ModelRow]) {
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
        let signal = intensity_signal(None, None, &ss, &pep(b"AR"), 2, "unknown", 20.0, true, None);
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
        let good = intensity_signal(Some(&model), None, &ss, &peptide, 2, "unknown", 20.0, true, None);

        let wrong_pep = pep(b"FGHIK");
        let bad = intensity_signal(Some(&model), None, &ss, &wrong_pep, 2, "unknown", 20.0, true, None);
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
        let partial = intensity_signal(Some(&model), None, &ss, &peptide, 2, "unknown", 20.0, true, None);
        let empty_spec = Spectrum {
            peaks: vec![(400.0, 1000.0)],
            ..Default::default()
        };
        let ss_empty = ScoredSpectrum::new_without_filtering(&empty_spec);
        let none = intensity_signal(Some(&model), None, &ss_empty, &peptide, 2, "unknown", 20.0, true, None);
        assert!(partial > none);
    }

    /// Build a constant-leaf GbdtPeakModel (apply_sigmoid=false) that returns
    /// `leaf_value` for every input. Uses the public struct fields + `to_bytes`
    /// / `from_bytes` round-trip to honour the validation path.
    fn const_leaf_gbdt(leaf_value: f32) -> crate::gbdt_eval::GbdtPeakModel {
        use crate::gbdt_eval::{GbdtPeakModel, Tree};
        let m = GbdtPeakModel {
            n_features: crate::frag_features::N_FRAG_FEATURES as u32,
            apply_sigmoid: false,
            trees: vec![Tree {
                feature: vec![-1],            // single leaf node
                threshold: vec![0.0],
                left: vec![-1],
                right: vec![-1],
                value: vec![leaf_value],
                default_left: vec![1] }],
            iso_x: vec![],
            iso_y: vec![],
        };
        // Round-trip through bytes to exercise the validator.
        GbdtPeakModel::from_bytes(&m.to_bytes()).expect("const-leaf gbdt round-trip")
    }

    #[test]
    fn predicted_linear_intensity_is_finite_for_extreme_leaf() {
        // An extreme positive leaf would `exp` to +inf without the clamp.
        let g_big = const_leaf_gbdt(1.0e6);
        let v = super::predicted_linear_intensity(&g_big, &[0.0f32; crate::frag_features::N_FRAG_FEATURES]);
        assert!(v.is_finite(), "clamped exp must be finite, got {v}");
        assert!(v > 0.0);
        // An extreme negative leaf clamps toward 0, still finite.
        let g_small = const_leaf_gbdt(-1.0e6);
        let v2 = super::predicted_linear_intensity(&g_small, &[0.0f32; crate::frag_features::N_FRAG_FEATURES]);
        assert!(v2.is_finite() && v2 >= 0.0, "got {v2}");
        // A normal leaf (0.0) is unchanged: exp(0) == 1.
        let g0 = const_leaf_gbdt(0.0);
        let v3 = super::predicted_linear_intensity(&g0, &[0.0f32; crate::frag_features::N_FRAG_FEATURES]);
        assert!((v3 - 1.0).abs() < 1e-9, "exp(0) must be 1.0, got {v3}");
    }

    #[test]
    fn frag_llr_battery_finite_with_extreme_model() {
        // A frag model with an extreme leaf must not produce NaN/inf features.
        let g = const_leaf_gbdt(1.0e6);
        let peptide = pep(b"ARCDE");
        let spec = Spectrum { peaks: vec![(300.0, 1000.0)], ..Default::default() };
        let ss = ScoredSpectrum::new_without_filtering(&spec);
        let (explained, chance, topk) = frag_llr_battery(Some(&g), &ss, &peptide, 2, 20.0, true, None);
        for (name, v) in [("explained", explained), ("chance", chance), ("topk", topk)] {
            assert!(v.is_finite(), "{name} must be finite, got {v}");
        }
    }

    #[test]
    fn intensity_signal_uses_frag_model_when_present() {
        // A constant-leaf frag model returning 0.0 (=> exp 1.0 for every ion).
        // Build ARCDE with a peak at y3 so the observed vector is non-zero.
        let g = const_leaf_gbdt(0.0); // predict_value returns 0.0 => exp(0) = 1.0

        let peptide = pep(b"ARCDE");
        let predicted = predict_by_ions(&peptide, 1..=1);
        let y3 = predicted
            .iter()
            .find(|p| p.kind == IonKind::Y && p.position == 3)
            .expect("y3");
        let spec = Spectrum {
            peaks: vec![(y3.mz, 1000.0)],
            ..Default::default()
        };
        let ss = ScoredSpectrum::new_without_filtering(&spec);

        // frag_model Some, table model None => should still compute cosine > 0.
        let sig = intensity_signal(None, Some(&g), &ss, &peptide, 2, "unknown", 20.0, true, None);
        assert!(sig > 0.0, "expected signal > 0 with frag model, got {sig}");
    }

    #[test]
    fn intensity_signal_falls_back_to_table_when_no_frag_model() {
        // With frag_model = None the result must equal the table-only path.
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

        let with_table = intensity_signal(Some(&model), None, &ss, &peptide, 2, "unknown", 20.0, true, None);
        // Same call was already exercised in `missing_observed_ions_reduce_signal`;
        // just verify it's > 0 (table path active, not early-return).
        assert!(with_table > 0.0, "table fallback expected > 0, got {with_table}");
    }

    #[test]
    fn rich_ion_llr_zero_without_model() {
        // No rich-ion model loaded ⇒ neutral additive feature (0.0), RankScore untouched.
        let peptide = pep(b"ARCDE");
        let predicted = predict_by_ions(&peptide, 1..=1);
        let y3 = predicted.iter().find(|p| p.kind == IonKind::Y && p.position == 3).unwrap();
        let spec = Spectrum { peaks: vec![(y3.mz, 1000.0)], ..Default::default() };
        let ss = ScoredSpectrum::new_without_filtering(&spec);
        assert_eq!(rich_ion_llr(None, &ss, &peptide, 2, 20.0, true), 0.0);
    }

    #[test]
    fn intensity_signal_zero_when_no_predictor_at_all() {
        let spec = Spectrum {
            peaks: vec![(500.0, 1000.0)],
            ..Default::default()
        };
        let ss = ScoredSpectrum::new_without_filtering(&spec);
        // Both None => must return 0.0 exactly.
        let sig = intensity_signal(None, None, &ss, &pep(b"ARCDE"), 2, "unknown", 20.0, true, None);
        assert_eq!(sig, 0.0, "no predictor must yield 0.0");
    }
}
