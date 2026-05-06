//! MS-GF+ search engine — now a facade that re-exports from focused crates.
//!
//! Model types → `model` crate
//! Scoring / GF / Param types → `scoring` crate
//! Search types → `search` crate
//! Output → still here (moves in Phase 4)

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
pub use scoring_crate::gf;
pub use scoring_crate::param_model;
pub use scoring_crate::scoring;

// Search modules re-exported from the `search` crate.
pub use search::candidate_gen;
pub use search::decoy;
pub use search::match_engine;
pub use search::precursor_matching;
pub use search::psm;
pub use search::search_index;
pub use search::search_params;
pub use search::suffix_array;

pub mod output;

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
pub use search::{
    enumerate_candidates, reverse_db, target_plus_decoy, DEFAULT_DECOY_PREFIX,
    match_spectra, matches_precursor, MassError,
    PsmFeatures, PsmMatch, TopNQueue,
    SearchIndex, SearchParams, SuffixArray,
};
