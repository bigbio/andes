//! The search engine: candidate peptide generation, the fragment index, precursor matching,
//! open and chimeric search, decoy handling and PSM scoring and refinement. It turns spectra from
//! `input` and a protein database into scored PSMs that `rescore` then calibrates.

pub mod candidates;
mod chimeric;
pub mod crosslinks;
pub mod decoy;
pub mod dia;
pub mod engine;
pub mod fragments;
pub mod glycopeptides;
pub mod index;
pub mod ipf;
pub mod isomer;
pub mod mimics;
pub mod open;
pub mod params;
pub mod precursor;
pub mod psm;
pub mod refinement;
pub mod sites;
pub mod tdc;

pub use candidates::enumerate_candidates;
pub use decoy::{decoy_db, normalize_decoy_prefix, DEFAULT_DECOY_PREFIX};
pub use precursor::prepass::{
    apply_shift_for_mode, apply_tightened_precursor_tolerance, build_spec_keys,
    learn_calibration_stats, prepass_search_params, run_by_file, CalibrationStats,
    FileCalibration, SpecKey,
};
pub use engine::{
    coisolation_pass, match_spectra, PreparedSearch, ScorerLoader,
};
pub use precursor::calibration::{
    adjusted_observed_neutral_mass, robust_sigma_ppm, tightened_tolerance_ppm, PrecursorCalMode,
};
pub use precursor::matching::{matches_precursor, MassError};
pub use psm::{PsmFeatures, PsmMatch, TopNQueue};
pub use refinement::config::{RefineConfig, RefineMod};
pub use index::SearchIndex;
pub use params::{ScoreMode, SearchParams};
pub use tdc::{confident_target_indices, ScoredLabel};
