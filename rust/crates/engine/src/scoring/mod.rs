//! Phase 5 scoring: replaces Phase 4e's mass-error placeholder with
//! rank-based PSM scoring using the Param model loaded in Phase 2.
//!
//! Submodules:
//! - `scored_spectrum`: per-spectrum precomputed state (peak ranks).
//! - More to come in Tasks 2-7 (fragment_ions, rank_scorer, psm_score).

pub mod scored_spectrum;

pub use scored_spectrum::ScoredSpectrum;
