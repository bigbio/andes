//! Phase 5 scoring: replaces Phase 4e's mass-error placeholder with
//! rank-based PSM scoring using the Param model loaded in Phase 2.

pub mod fragment_ions;
pub mod rank_scorer;
pub mod scored_spectrum;

pub use fragment_ions::{predict_by_ions, PredictedIon};
pub use rank_scorer::RankScorer;
pub use scored_spectrum::ScoredSpectrum;
