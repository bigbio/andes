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

// Convenience re-exports for the most-used types. Downstream crates
// (input/output/cli) prefer `use engine::Peptide` over the qualified path.
pub use aa_set::{AaSetError, AminoAcidSet, AminoAcidSetBuilder};
pub use activation::ActivationMethod;
pub use instrument::InstrumentType;
pub use amino_acid::AminoAcid;
pub use compact_fasta::{CompactFastaError, CompactFastaSequence, ProteinAnnotation};
pub use enzyme::Enzyme;
pub use mass::{nominal_from, C, H, H2O, N, O, PROTON, S, INTEGER_MASS_SCALER};
pub use modification::{ModLocation, ModParseError, Modification, ResidueSpec};
pub use param_model::{
    FragmentOffsetFrequency, IonType, Param, ParamParseError, Partition,
    PrecursorOffsetFrequency, SpecDataType,
};
pub use candidate_gen::{enumerate_candidates, Candidate};
pub use peptide::{Peptide, PeptideParseError};
pub use precursor_matching::{matches_precursor, MassError};
pub use decoy::{reverse_db, target_plus_decoy, DEFAULT_DECOY_PREFIX};
pub use protein::{Protein, ProteinDb};
pub use protocol::Protocol;
pub use match_engine::match_spectra;
pub use psm::{PsmMatch, TopNQueue};
pub use search_index::{SearchIndex, SearchIndexError};
pub use gf::{ScoreBound, ScoreDist, GeneratingFunction, GfError, PrimitiveAaGraph, GeneratingFunctionGroup};
pub use scoring::{score_psm, RankScorer, ScoredSpectrum};
pub use search_params::SearchParams;
pub use spectrum::Spectrum;
pub use suffix_array::{SuffixArray, SuffixArrayError};
pub use tolerance::{PrecursorTolerance, Tolerance};
