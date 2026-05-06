//! Domain model for MS-GF+ Rust port.
//!
//! Phase 1 milestone: amino acids, modifications, peptides, enzymes,
//! tolerances. Pure CPU + types. No I/O. See
//! `docs/superpowers/2026-05-03-phase1-engine-domain-model-design.md`.

pub mod aa_set;
pub mod activation;
pub mod amino_acid;
pub mod candidate_gen;
pub mod compact_fasta;
pub mod decoy;
pub mod enzyme;
pub mod gf;
pub mod instrument;
pub mod mass;
pub mod match_engine;
pub mod modification;
pub mod output;
pub mod param_model;
pub mod peptide;
pub mod precursor_matching;
pub mod protein;
pub mod protocol;
pub mod psm;
pub mod scoring;
pub mod search_index;
pub mod search_params;
pub mod spectrum;
pub mod suffix_array;
pub mod tolerance;

#[cfg(test)]
pub(crate) mod testutil;

// Convenience re-exports for the most-used types. Downstream crates
// (input, cli, integration tests) prefer `use engine::Peptide` over the
// qualified path. Internal plumbing (GF, Param sub-types, scoring
// internals) is intentionally NOT re-exported here — use the qualified
// `engine::module::Type` path for those.
pub use aa_set::{AaSetError, AminoAcidSet, AminoAcidSetBuilder};
pub use activation::ActivationMethod;
pub use amino_acid::AminoAcid;
pub use candidate_gen::enumerate_candidates;
pub use compact_fasta::{CompactFastaError, CompactFastaSequence, ProteinAnnotation};
pub use decoy::{reverse_db, target_plus_decoy, DEFAULT_DECOY_PREFIX};
pub use enzyme::Enzyme;
pub use instrument::InstrumentType;
pub use mass::{nominal_from, H2O, PROTON};
pub use match_engine::match_spectra;
pub use modification::{ModLocation, ModParseError, Modification, ResidueSpec};
pub use param_model::{Param, ParamParseError};
pub use peptide::Peptide;
pub use precursor_matching::{matches_precursor, MassError};
pub use protein::{Protein, ProteinDb};
pub use protocol::Protocol;
pub use psm::{PsmFeatures, PsmMatch, TopNQueue};
pub use scoring::{RankScorer, ScoredSpectrum};
pub use search_index::SearchIndex;
pub use search_params::SearchParams;
pub use spectrum::Spectrum;
pub use suffix_array::SuffixArray;
pub use tolerance::{PrecursorTolerance, Tolerance};
