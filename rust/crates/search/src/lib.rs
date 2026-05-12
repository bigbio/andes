//! Search sub-system for MS-GF+ Rust port.
//!
//! Contains candidate generation, suffix array, search index, precursor
//! matching, PSM structures, and the match engine.
//! Depends on `model` and `scoring` crates.

pub mod candidate_gen;
pub mod decoy;
pub mod distinct_peptide;
pub mod match_engine;
pub mod precursor_matching;
pub mod psm;
pub mod sa_walk;
pub mod search_index;
pub mod search_params;
pub mod suffix_array;

// Convenience re-exports.
pub use candidate_gen::enumerate_candidates;
pub use decoy::{reverse_db, target_plus_decoy, DEFAULT_DECOY_PREFIX};
pub use match_engine::{match_spectra, PreparedSearch};
pub use precursor_matching::{matches_precursor, MassError};
pub use psm::{PsmFeatures, PsmMatch, TopNQueue};
pub use search_index::SearchIndex;
pub use search_params::SearchParams;
pub use suffix_array::SuffixArray;
