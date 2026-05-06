//! Domain model for MS-GF+ Rust port.
//!
//! Phase 1 milestone: amino acids, modifications, peptides, enzymes,
//! tolerances. Pure CPU + types. No I/O. See
//! `docs/superpowers/2026-05-03-phase1-engine-domain-model-design.md`.

// Model modules now live in the `model` crate — re-export for compatibility.
pub use model::aa_set;
pub use model::activation;
pub use model::amino_acid;
pub use model::compact_fasta;
pub use model::enzyme;
pub use model::instrument;
pub use model::mass;
pub use model::modification;
pub use model::peptide;
pub use model::protein;
pub use model::protocol;
pub use model::spectrum;
pub use model::tolerance;

pub mod candidate_gen;
pub mod decoy;
pub mod gf;
pub mod match_engine;
pub mod output;
pub mod param_model;
pub mod precursor_matching;
pub mod psm;
pub mod scoring;
pub mod search_index;
pub mod search_params;
pub mod suffix_array;

#[cfg(test)]
pub(crate) mod testutil;

// Convenience re-exports for the most-used types. Downstream crates
// (input, cli, integration tests) prefer `use engine::Peptide` over the
// qualified path. Internal plumbing (GF, Param sub-types, scoring
// internals) is intentionally NOT re-exported here — use the qualified
// `engine::module::Type` path for those.
pub use model::{
    AaSetError, AminoAcidSet, AminoAcidSetBuilder,
    ActivationMethod,
    AminoAcid,
    CompactFastaError, CompactFastaSequence, ProteinAnnotation,
    Enzyme,
    InstrumentType,
    nominal_from, H2O, PROTON,
    ModLocation, ModParseError, Modification, ResidueSpec,
    Peptide,
    Protein, ProteinDb,
    Protocol,
    Spectrum,
    PrecursorTolerance, Tolerance,
};
pub use candidate_gen::enumerate_candidates;
pub use decoy::{reverse_db, target_plus_decoy, DEFAULT_DECOY_PREFIX};
pub use match_engine::match_spectra;
pub use param_model::{Param, ParamParseError};
pub use precursor_matching::{matches_precursor, MassError};
pub use psm::{PsmFeatures, PsmMatch, TopNQueue};
pub use scoring::{RankScorer, ScoredSpectrum};
pub use search_index::SearchIndex;
pub use search_params::SearchParams;
pub use suffix_array::SuffixArray;
