//! Peptide database search engine: candidate enumeration, precursor matching,
//! scoring, and PSM aggregation.
//!
//! Contains candidate generation, suffix array, search index, precursor
//! matching, PSM structures, and the match engine.
//! Depends on `model` and `scoring` crates.

pub mod candidate_gen;
mod chimeric_features;
pub(crate) mod coisolation;
pub mod decoy;
pub mod distinct_peptide;
pub mod match_engine;
pub mod mass_calibrator;
pub mod precursor_cal;
pub mod precursor_matching;
pub mod psm;
pub mod sa_walk;
pub mod search_index;
pub mod search_params;
pub(crate) mod fragment_index;
pub(crate) mod sage_index;
pub(crate) mod shared_fragment;
pub mod suffix_array;

// Convenience re-exports.
pub use candidate_gen::enumerate_candidates;
pub use decoy::{reverse_db, target_plus_decoy, DEFAULT_DECOY_PREFIX};
pub use match_engine::{match_spectra, run_pass2_coisolation, PreparedSearch};
pub use mass_calibrator::{
    apply_shift_for_mode, apply_tightened_precursor_tolerance, build_spec_keys,
    learn_calibration_stats, prepass_search_params, CalibrationStats, SpecKey,
};
pub use precursor_cal::{
    PrecursorCalMode, adjusted_observed_neutral_mass, robust_sigma_ppm, tightened_tolerance_ppm,
};
pub use precursor_matching::{matches_precursor, MassError};
pub use psm::{PsmFeatures, PsmMatch, TopNQueue};
pub use search_index::SearchIndex;
pub use search_params::SearchParams;
pub use suffix_array::SuffixArray;
