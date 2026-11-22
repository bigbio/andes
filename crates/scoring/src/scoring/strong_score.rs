use model::peptide::Peptide;

pub const DENSITY_HW: f64 = 50.0;

const LOG_INTENSITY_CLAMP: f64 = 50.0;

pub fn predict_frag_intensities(
    frag_model: &GbdtPeakModel,
    peptide: &Peptide,
    precursor_charge: u8,
) -> Vec<f64> {
    frag_intensities_from_raw(&predict_frag_raw(frag_model, peptide, precursor_charge))
}

/// The model's raw (log relative) value per b/y ion of charge 1 and 2 of `peptide`, in the
/// order of [`predict_by_ions`]; [`frag_intensities_from_raw`] turns them into the
/// intensities of [`predict_frag_intensities`].
pub fn predict_frag_raw(
    frag_model: &GbdtPeakModel,
    peptide: &Peptide,
    precursor_charge: u8,
) -> Vec<f32> {
    let ions = predict_by_ions(peptide, 1..=2);
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
    predict_frag_values_memo(frag_model, &feats)
}

/// Linear intensities from the raw values of [`predict_frag_raw`].
pub fn frag_intensities_from_raw(raw: &[f32]) -> Vec<f64> {
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

type FragRowKey = [u32; crate::frag_features::N_FRAG_FEATURES];

const FRAG_MEMO_CAP: usize = 1 << 17;

thread_local! {
    static FRAG_MEMO: std::cell::RefCell<(usize, rustc_hash::FxHashMap<FragRowKey, f32>)> =
        std::cell::RefCell::new((0, rustc_hash::FxHashMap::default()));
}

fn predict_frag_values_memo(
    frag_model: &GbdtPeakModel,
    feats: &[[f32; crate::frag_features::N_FRAG_FEATURES]],
) -> Vec<f32> {
    let model_id = frag_model.trees.as_ptr() as usize ^ frag_model.trees.len().rotate_left(48);
    let key =
        |f: &[f32; crate::frag_features::N_FRAG_FEATURES]| -> FragRowKey { f.map(f32::to_bits) };
    FRAG_MEMO.with(|cell| {
        let mut memo = cell.borrow_mut();
        if memo.0 != model_id {
            memo.0 = model_id;
            memo.1.clear();
        }
        let mut out = vec![0.0f32; feats.len()];
        let mut miss_idx: Vec<usize> = Vec::new();
        for (i, f) in feats.iter().enumerate() {
            match memo.1.get(&key(f)) {
                Some(&v) => out[i] = v,
                None => miss_idx.push(i),
            }
        }
        if !miss_idx.is_empty() {
            let rows: Vec<&[f32]> = miss_idx.iter().map(|&i| feats[i].as_slice()).collect();
            let mut raw = vec![0.0f32; rows.len()];
            frag_model.predict_value_batch(&rows, &mut raw);
            if memo.1.len() + miss_idx.len() > FRAG_MEMO_CAP {
                memo.1.clear();
            }
            for (&i, &v) in miss_idx.iter().zip(&raw) {
                out[i] = v;
                memo.1.insert(key(&feats[i]), v);
            }
        }
        out
    })
}

#[inline]
fn predicted_linear_intensity(g: &GbdtPeakModel, feats: &[f32]) -> f64 {
    let log_rel = f64::from(g.predict_value(feats));
    if !log_rel.is_finite() {
        return 0.0;
    }
    log_rel
        .clamp(-LOG_INTENSITY_CLAMP, LOG_INTENSITY_CLAMP)
        .exp()
}

use crate::frag_features::extract_frag_features;
use crate::gbdt_eval::GbdtPeakModel;
use crate::ion_features::extract_ion_features;
use crate::scoring::fragment_ions::{predict_by_ions, IonKind, PredictedIon};
use crate::scoring::ScoredSpectrum;

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

pub fn intensity_signal(
    frag_model: Option<&GbdtPeakModel>,
    scored_spec: &ScoredSpectrum<'_>,
    peptide: &Peptide,
    precursor_charge: u8,
    feature_tol: f64,
    feature_tol_is_ppm: bool,
    precomputed: Option<&[f64]>,
) -> f32 {
    let Some(g) = frag_model else {
        return 0.0;
    };
    let n = peptide.length();
    if n < 2 {
        return 0.0;
    }

    let base_peak = scored_spec.max_active_intensity();
    if base_peak <= 0.0 {
        return 0.0;
    }

    let predicted = predict_by_ions(peptide, 1..=2);
    let mut pred_vec = Vec::with_capacity(predicted.len());
    let mut obs_vec = Vec::with_capacity(predicted.len());

    for (ion_idx, ion) in predicted.iter().enumerate() {
        let pred_intensity = match precomputed.and_then(|p| p.get(ion_idx).copied()) {
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

pub fn frag_llr_battery(
    frag_model: Option<&GbdtPeakModel>,
    scored_spec: &ScoredSpectrum<'_>,
    peptide: &Peptide,
    precursor_charge: u8,
    feature_tol: f64,
    feature_tol_is_ppm: bool,
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
    let mut pred_matched: Vec<(f64, bool)> = Vec::with_capacity(predicted.len());

    for (ion_idx, ion) in predicted.iter().enumerate() {
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

    pred_matched.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
    let k = pred_matched.len().min(FRAG_TOPK);
    let topk_observed = if k > 0 {
        pred_matched[..k].iter().filter(|(_, m)| *m).count() as f32 / k as f32
    } else {
        0.0
    };

    (explained, sum_chance_llr as f32, topk_observed)
}

/// Neutral-loss class of a theoretical fragment ion: none, water, ammonia, or a
/// modification-specific loss identified by its `PredictedIon::loss_class`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FragmentLoss {
    None,
    Water,
    Ammonia,
    Modification(u8),
}

/// Identity of one theoretical fragment ion of a peptide: series, position and loss class.
/// Two different ions can share an m/z bit for bit; identity, not m/z, tells them apart.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FragmentIon {
    pub kind: IonKind,
    pub position: u32,
    pub loss: FragmentLoss,
}

impl FragmentIon {
    pub fn of_predicted(p: &PredictedIon) -> Self {
        Self {
            kind: p.kind,
            position: p.position,
            loss: if p.loss_class == 0 {
                FragmentLoss::None
            } else {
                FragmentLoss::Modification(p.loss_class)
            },
        }
    }
}

pub type MatchedIon = (f32, f64, f64, bool, FragmentIon);

pub type TheoreticalIons = Vec<(f64, FragmentIon)>;

/// Charge-1 theoretical ions of a peptide with the water and ammonia loss of each,
/// every entry keeping the identity of the ion it comes from.
pub fn theoretical_fragment_ions(charge1_predicted: &[PredictedIon]) -> TheoreticalIons {
    let mut out = TheoreticalIons::with_capacity(3 * charge1_predicted.len());
    for p in charge1_predicted {
        let ion = FragmentIon::of_predicted(p);
        out.push((p.mz, ion));
        out.push((
            p.mz - model::mass::H2O,
            FragmentIon {
                loss: FragmentLoss::Water,
                ..ion
            },
        ));
        out.push((
            p.mz - model::mass::NH3,
            FragmentIon {
                loss: FragmentLoss::Ammonia,
                ..ion
            },
        ));
    }
    out
}

/// Number of theoretical ions other than `matched` that lie within the feature tolerance
/// of the observed peak at `obs`.
pub fn fragment_ambiguity(
    theoretical: &[(f64, FragmentIon)],
    obs: f64,
    feature_tol: f64,
    feature_tol_is_ppm: bool,
    matched: FragmentIon,
) -> usize {
    let tol_da = if feature_tol_is_ppm {
        obs * feature_tol / 1e6
    } else {
        feature_tol
    };
    theoretical
        .iter()
        .filter(|&&(theo, ion)| (theo - obs).abs() <= tol_da && ion != matched)
        .count()
}

pub fn mass_competition_evidence(
    scored_spec: &ScoredSpectrum<'_>,
    matched_ions: &[MatchedIon],
    theoretical: &[(f64, FragmentIon)],
    feature_tol: f64,
    feature_tol_is_ppm: bool,
) -> f32 {
    matched_ions
        .iter()
        .map(|&(_, obs, _, _, ion)| {
            let ambiguity =
                fragment_ambiguity(theoretical, obs, feature_tol, feature_tol_is_ppm, ion);
            let rho = scored_spec.local_peak_density(obs, DENSITY_HW);
            let competition = ambiguity as f64 + rho;
            1.0 / (1.0 + competition)
        })
        .sum::<f64>() as f32
}

pub fn listwise_score_gap(scores_best_first: &[f32]) -> f32 {
    if scores_best_first.len() < 2 {
        return 0.0;
    }
    (scores_best_first[0] - scores_best_first[1]).max(0.0)
}

pub fn candidate_rank_entropy(scores: &[f32]) -> f32 {
    if scores.len() < 2 {
        return 0.0;
    }
    let max_s = scores.iter().copied().fold(f32::NEG_INFINITY, f32::max);
    let exp_scores: Vec<f64> = scores.iter().map(|&s| f64::from(s - max_s).exp()).collect();
    let sum: f64 = exp_scores.iter().sum();
    if sum <= 0.0 {
        return 0.0;
    }
    exp_scores
        .iter()
        .map(|&e| {
            let p = e / sum;
            if p > 0.0 {
                -p * p.ln()
            } else {
                0.0
            }
        })
        .sum::<f64>() as f32
}

#[derive(Debug, Clone, Copy, Default)]
pub struct StrongScoreInputs {
    pub intensity_signal: f32,
    pub chance_match_surprise: f32,
    pub mass_competition_evidence: f32,
    pub candidate_rank_entropy: f32,
    pub listwise_score_gap: f32,
}

pub fn fuse_strong_score(f: &StrongScoreInputs) -> f32 {
    let null = -f.chance_match_surprise - f.mass_competition_evidence + f.candidate_rank_entropy
        - f.listwise_score_gap;
    f.intensity_signal - null
}

pub const STRONG_CAL_MIN_CANDIDATES: u32 = 100;

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

pub fn strong_score_zscore(score: f32, null: &OnlineStats) -> f32 {
    if null.n < STRONG_CAL_MIN_CANDIDATES as u64 {
        return score;
    }
    ((f64::from(score) - null.mean()) / null.population_stdev()) as f32
}

/// Fewest other retained scores the leave-one-out z-score of `strong_score_calibrated_loo`
/// estimates its spread from; with fewer, the PSM gets 0 (no calibrated evidence).
pub const STRONG_CAL_MIN_LOO_OTHERS: usize = 4;

pub fn strong_score_calibrated_loo(retained_strong: &[f32], this: f32) -> f32 {
    if retained_strong.len() < 2 {
        return this;
    }
    let n = retained_strong.len();
    if n - 1 < STRONG_CAL_MIN_LOO_OTHERS {
        return 0.0;
    }
    let this_d = f64::from(this);
    let sum: f64 = retained_strong.iter().map(|&s| f64::from(s)).sum();
    let n_others = (n - 1) as f64;
    let mean_others = (sum - this_d) / n_others;
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
    let scale = mean_others.abs().max(1.0);
    let sigma = var_others.sqrt();
    if sigma <= 1e-9 * scale {
        return 0.0;
    }
    ((this_d - mean_others) / sigma) as f32
}

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
    let mut rows_buf: Vec<[f32; crate::ion_features::N_ION_FEATURES]> = Vec::new();
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
            ));
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
    use model::amino_acid::AminoAcid;
    use model::peptide::Peptide;
    use model::spectrum::Spectrum;

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
    fn intensity_signal_zero_without_model() {
        let spec = Spectrum {
            peaks: vec![(500.0, 1000.0)],
            ..Default::default()
        };
        let ss = ScoredSpectrum::new_without_filtering(&spec);
        let signal = intensity_signal(None, &ss, &pep(b"AR"), 2, 20.0, true, None);
        assert_eq!(signal, 0.0);
    }

    #[test]
    fn a_coincident_ion_of_another_series_counts_as_ambiguity() {
        let peptide = pep(b"GAKAG");
        let predicted = predict_by_ions(&peptide, 1..=1);
        let theoretical = theoretical_fragment_ions(&predicted);
        let b2 = predicted
            .iter()
            .find(|p| p.kind == IonKind::B && p.position == 2)
            .unwrap();
        let y2 = predicted
            .iter()
            .find(|p| p.kind == IonKind::Y && p.position == 2)
            .unwrap();
        let y2_water_loss = y2.mz - model::mass::H2O;
        assert!((y2_water_loss - b2.mz).abs() < 1e-9);
        let n = fragment_ambiguity(
            &theoretical,
            b2.mz,
            20.0,
            true,
            FragmentIon::of_predicted(b2),
        );
        assert_eq!(n, 1);
    }

    #[test]
    fn an_unambiguous_peak_has_no_competitor_and_full_evidence() {
        let peptide = pep(b"ARCDE");
        let predicted = predict_by_ions(&peptide, 1..=1);
        let theoretical = theoretical_fragment_ions(&predicted);
        let y3 = predicted
            .iter()
            .find(|p| p.kind == IonKind::Y && p.position == 3)
            .unwrap();
        let n = fragment_ambiguity(
            &theoretical,
            y3.mz,
            20.0,
            true,
            FragmentIon::of_predicted(y3),
        );
        assert_eq!(n, 0);
        let sparse = Spectrum {
            peaks: vec![(y3.mz, 1000.0)],
            ..Default::default()
        };
        let ss = ScoredSpectrum::new_without_filtering(&sparse);
        let matched = vec![(
            1000.0_f32,
            y3.mz,
            y3.mz,
            false,
            FragmentIon::of_predicted(y3),
        )];
        let ev = mass_competition_evidence(&ss, &matched, &theoretical, 20.0, true);
        let expected = 1.0 / (1.0 + 1.0 / (2.0 * DENSITY_HW));
        assert!((f64::from(ev) - expected).abs() < 1e-6);
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
        let theo = theoretical_fragment_ions(&predicted);
        let matched = vec![(
            1000.0_f32,
            y3.mz,
            y3.mz,
            false,
            FragmentIon::of_predicted(y3),
        )];
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
        assert_eq!(strong_score_zscore(5.0, &null), 5.0);
        for _ in 0..100 {
            null.push(10.0);
        }
        let z = strong_score_zscore(15.0, &null);
        assert!(z > 0.0);
        let loo = strong_score_calibrated_loo(&[1.0, 3.0, 5.0, 2.0, 4.0], 5.0);
        assert!(loo > 0.0);
        assert_eq!(strong_score_calibrated_loo(&[2.0], 2.0), 2.0);
    }

    #[test]
    fn a_leave_one_out_z_score_needs_four_other_scores() {
        let pool = [22.766172_f32, 32.432117, 22.761606];
        assert_eq!(strong_score_calibrated_loo(&pool, 32.432117), 0.0);
        let pool = [22.766172_f32, 32.432117, 22.761606, 22.7601, 22.7655];
        let z = strong_score_calibrated_loo(&pool, 32.432117);
        assert!(z > 0.0 && z.is_finite());
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

    fn const_leaf_gbdt(leaf_value: f32) -> crate::gbdt_eval::GbdtPeakModel {
        use crate::gbdt_eval::{GbdtPeakModel, Tree};
        let m = GbdtPeakModel {
            packed: Default::default(),
            n_features: crate::frag_features::N_FRAG_FEATURES as u32,
            feature_names: crate::gbdt_eval::indexed_feature_names(
                crate::frag_features::N_FRAG_FEATURES as u32,
            ),
            apply_sigmoid: false,
            trees: vec![Tree {
                feature: vec![-1],
                threshold: vec![0.0],
                left: vec![-1],
                right: vec![-1],
                value: vec![leaf_value],
                default_left: vec![1],
            }],
            iso_x: vec![],
            iso_y: vec![],
        };
        GbdtPeakModel::from_bytes(&m.to_bytes()).expect("const-leaf gbdt round-trip")
    }

    #[test]
    fn predicted_linear_intensity_is_finite_for_extreme_leaf() {
        let g_big = const_leaf_gbdt(1.0e6);
        let v = super::predicted_linear_intensity(
            &g_big,
            &[0.0f32; crate::frag_features::N_FRAG_FEATURES],
        );
        assert!(v.is_finite(), "clamped exp must be finite, got {v}");
        assert!(v > 0.0);
        let g_small = const_leaf_gbdt(-1.0e6);
        let v2 = super::predicted_linear_intensity(
            &g_small,
            &[0.0f32; crate::frag_features::N_FRAG_FEATURES],
        );
        assert!(v2.is_finite() && v2 >= 0.0, "got {v2}");
        let g0 = const_leaf_gbdt(0.0);
        let v3 = super::predicted_linear_intensity(
            &g0,
            &[0.0f32; crate::frag_features::N_FRAG_FEATURES],
        );
        assert!((v3 - 1.0).abs() < 1e-9, "exp(0) must be 1.0, got {v3}");
    }

    #[test]
    fn frag_llr_battery_finite_with_extreme_model() {
        let g = const_leaf_gbdt(1.0e6);
        let peptide = pep(b"ARCDE");
        let spec = Spectrum {
            peaks: vec![(300.0, 1000.0)],
            ..Default::default()
        };
        let ss = ScoredSpectrum::new_without_filtering(&spec);
        let (explained, chance, topk) =
            frag_llr_battery(Some(&g), &ss, &peptide, 2, 20.0, true, None);
        for (name, v) in [("explained", explained), ("chance", chance), ("topk", topk)] {
            assert!(v.is_finite(), "{name} must be finite, got {v}");
        }
    }

    #[test]
    fn intensity_signal_uses_frag_model_when_present() {
        let g = const_leaf_gbdt(0.0);

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

        let sig = intensity_signal(Some(&g), &ss, &peptide, 2, 20.0, true, None);
        assert!(sig > 0.0, "expected signal > 0 with frag model, got {sig}");
    }

    #[test]
    fn rich_ion_llr_zero_without_model() {
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
        assert_eq!(rich_ion_llr(None, &ss, &peptide, 2, 20.0, true), 0.0);
    }

    #[test]
    fn intensity_signal_zero_when_no_predictor_at_all() {
        let spec = Spectrum {
            peaks: vec![(500.0, 1000.0)],
            ..Default::default()
        };
        let ss = ScoredSpectrum::new_without_filtering(&spec);
        let sig = intensity_signal(None, &ss, &pep(b"ARCDE"), 2, 20.0, true, None);
        assert_eq!(sig, 0.0, "no predictor must yield 0.0");
    }
}
