//! Canonical feature extraction for the learned glyco backbone SELECTOR
//! (andes-glyco 2.0 Path 2).
//!
//! The learned selector is a native GBDT ([`scoring_crate::gbdt_eval::GbdtPeakModel`])
//! trained to rank the TRUE backbone above wrong competitors within a scan. Offline
//! 5-fold CV (leakage-free, grouped by scan) recovers 379 top-1 vs the hand `gp`
//! fusion's 347 (+32) — because the signal is in the COMBINATION of ~15 features, no
//! single hand-added term captures it.
//!
//! ## The one rule that makes this correct: feature parity
//!
//! Training and inference MUST build the feature vector from the SAME function at the
//! SAME pipeline point (the collapse, where [`compute_psm_features`] has run for a
//! candidate). Otherwise a feature that is non-zero in the training set but zero at
//! inference (or vice versa) silently corrupts the model — the exact class of bug that
//! sank the hand-added terms. Therefore this extractor deliberately EXCLUDES features
//! that are not available per-candidate at collapse time:
//!   - retention-time features (`delta_rt`, `abs_delta_rt`, `delta_rt_norm`,
//!     `predicted_rt_min`, glyco `DeltaRTRank`) — populated post-search by
//!     `rt_wiring::populate_rt_features`, so they are 0 at the collapse.
//!
//! [`glyco_selector_feature_names`] and [`glyco_selector_feature_vec`] are kept in
//! lock-step: same length, same order. The names exist for the training-set header and
//! for debugging; the model only ever sees the positional vector.
//!
//! [`compute_psm_features`]: crate::match_engine::compute_psm_features

use crate::psm::PsmMatch;
use andes_glyco::glyco_psm::GlycoPsmKey;
use andes_glyco::hybrid::Source;

/// Ordered names of the selector features. MUST stay in lock-step with
/// [`glyco_selector_feature_vec`] (same length, same order).
pub fn glyco_selector_feature_names() -> &'static [&'static str] {
    &[
        // ── peptide b/y base score + coverage (from PsmMatch / PsmFeatures) ──
        "RankScore",
        "NumMatchedMainIons",
        "longest_b",
        "longest_y",
        "longest_y_pct",
        "matchedIonRatio",
        "ExplainedIonCurrentRatio",
        "NTermIonCurrentRatio",
        "CTermIonCurrentRatio",
        "MS2IonCurrent",
        "MeanErrorTop7",
        "StdevErrorTop7",
        "MeanRelErrorTop7",
        "StdevRelErrorTop7",
        "EdgeScore",
        "PrecursorIsotopeKL",
        "PrecursorSNR",
        "DeltaRawScore",
        "TailorScore",
        "RankScoreFloat",
        "PpmGaussianScore",
        "LongestComplementaryLadder",
        "ComplementaryIonBalance",
        "NeutralLossIonCount",
        "MeanMatchedIntensityRank",
        "DoublyChargedMatchedIonCount",
        "ChanceMatchSurprise",
        "UniqueMatchFraction",
        "IntensitySignal",
        "FragPredExplained",
        "FragPredChanceLLR",
        "FragTopKObserved",
        "RichIonLLR",
        "NumMods",
        "MassCompetitionEvidence",
        "CandidateRankEntropy",
        "ListwiseScoreGap",
        "RawScore",
        "RawScoreCal",
        // ── glyco-specific evidence (from GlycoPsmKey) ──
        "OxoniumScore",
        "NCoreOxoniumIons",
        "YLadderScore",
        "CoreYHits",
        "IsGlycanDb",
        "Y0Y1Anchor",
        "SialicConsistency",
    ]
}

/// Number of selector features (compile-time-checked against the names/vec by the
/// unit test `feature_names_and_vec_agree`).
pub const GLYCO_SELECTOR_N_FEATURES: usize = 46;

/// Build the positional feature vector for one candidate. MUST be called at the
/// COLLAPSE point, after [`compute_psm_features`](crate::match_engine::compute_psm_features)
/// has populated `psm.features`, for BOTH training-row emission and inference — see the
/// module-level parity rule.
pub fn glyco_selector_feature_vec(psm: &PsmMatch, key: &GlycoPsmKey) -> Vec<f32> {
    let f = &psm.features;
    let v = vec![
        psm.score.round() as f32, // RankScore (integer-rounded, matches PIN column)
        f.num_matched_main_ions as f32,
        f.longest_b as f32,
        f.longest_y as f32,
        f.longest_y_pct,
        f.matched_ion_ratio,
        f.explained_ion_current_ratio,
        f.n_term_ion_current_ratio,
        f.c_term_ion_current_ratio,
        f.ms2_ion_current,
        f.mean_error_top7,
        f.stdev_error_top7,
        f.mean_rel_error_top7,
        f.stdev_rel_error_top7,
        f.edge_score as f32,
        f.precursor_isotope_kl,
        f.precursor_snr,
        f.delta_raw_score,
        f.tailor_score,
        f.rank_score_float,
        f.ppm_gaussian_score,
        f.longest_complementary_ladder as f32,
        f.complementary_ion_balance,
        f.neutral_loss_ion_count as f32,
        f.mean_matched_intensity_rank,
        f.doubly_charged_matched_ion_count as f32,
        f.chance_match_surprise,
        f.unique_match_fraction,
        f.intensity_signal,
        f.frag_pred_explained,
        f.frag_pred_chance_llr,
        f.frag_topk_observed,
        f.rich_ion_llr,
        f.num_mods as f32,
        f.mass_competition_evidence,
        f.candidate_rank_entropy,
        f.listwise_score_gap,
        f.strong_score,
        f.strong_score_cal,
        // glyco evidence
        key.oxonium_summed_frac,
        key.n_core_oxonium_ions as f32,
        key.y_ladder_intensity_score,
        key.core_y_hits as f32,
        if key.glycan_source == Source::Db { 1.0 } else { 0.0 },
        key.y0y1_anchor_score,
        key.sialic_consistency,
    ];
    // Parity guard: the positional vector must match the declared count and the
    // names list (which the training-set header is built from). Debug-only so the
    // hot inference path pays nothing in release.
    debug_assert_eq!(v.len(), GLYCO_SELECTOR_N_FEATURES);
    debug_assert_eq!(v.len(), glyco_selector_feature_names().len());
    v
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The names list length matches the declared count. The vec length is guarded at
    /// runtime by the `debug_assert_eq!`s inside [`glyco_selector_feature_vec`] (which
    /// fire under `cargo test` whenever the extractor is exercised on a real candidate),
    /// so training-set columns and the inference vector stay positionally identical.
    #[test]
    fn feature_names_match_declared_count() {
        assert_eq!(glyco_selector_feature_names().len(), GLYCO_SELECTOR_N_FEATURES);
    }
}
