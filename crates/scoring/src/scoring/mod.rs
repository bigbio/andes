//! Rank-based PSM scoring using the loaded Param model.

pub mod fragment_ions;
pub mod psm_score;
pub mod rank_scorer;
pub mod scored_spectrum;
pub mod strong_score;

pub use fragment_ions::{predict_by_ions, predict_cz_ions, PredictedIon};
pub use psm_score::{
    cz_hyperscore_psm, cz_matched_intensity_frac, cz_structure_features, hyperscore_psm,
    hyperscore_psm_with_matches,
    psm_edge_existence_facts,
    psm_edge_score, score_psm,
    score_psm_float,
    init_cz_settings, CzSettings,
};
pub use rank_scorer::RankScorer;
pub use scored_spectrum::{init_scoring_settings, IonMatchFact, ScoredSpectrum, ScoringSettings};
pub use strong_score::{
    candidate_rank_entropy, frag_llr_battery, fuse_strong_score, intensity_signal, listwise_score_gap,
    predict_frag_intensities,
    rich_ion_llr,
    mass_competition_evidence, strong_score_calibrated, strong_score_calibrated_loo,
    strong_score_zscore, OnlineStats, StrongScoreInputs, DENSITY_HW, STRONG_CAL_MIN_CANDIDATES,
};
