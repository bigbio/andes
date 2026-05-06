//! MS-GF+ search engine — now a facade that re-exports from focused crates.
//!
//! Model types → `model` crate
//! Scoring / GF / Param types → `scoring` crate
//! Remaining: candidate_gen, search_index, match_engine, output, etc.

// Model modules re-exported from the `model` crate.
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

// Scoring modules re-exported from the `scoring` crate.
// Note: the `scoring` sub-module is re-exported as `scoring` here;
// the `gf` and `param_model` sub-modules come from the scoring crate root.
pub use scoring_crate::gf;
pub use scoring_crate::param_model;
pub use scoring_crate::scoring;

pub mod candidate_gen;
pub mod decoy;
pub mod match_engine;
pub mod output;
pub mod precursor_matching;
pub mod psm;
pub mod search_index;
pub mod search_params;
pub mod suffix_array;

#[cfg(test)]
pub(crate) mod testutil;

// Convenience re-exports for the most-used types.
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
pub use scoring_crate::{Param, ParamParseError, RankScorer, ScoredSpectrum};
pub use candidate_gen::enumerate_candidates;
pub use decoy::{reverse_db, target_plus_decoy, DEFAULT_DECOY_PREFIX};
pub use match_engine::match_spectra;
pub use precursor_matching::{matches_precursor, MassError};
pub use psm::{PsmFeatures, PsmMatch, TopNQueue};
pub use search_index::SearchIndex;
pub use search_params::SearchParams;
pub use suffix_array::SuffixArray;
