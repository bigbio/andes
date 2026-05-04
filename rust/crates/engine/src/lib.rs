//! Domain model for MS-GF+ Rust port.
//!
//! Phase 1 milestone: amino acids, modifications, peptides, enzymes,
//! tolerances. Pure CPU + types. No I/O. See
//! `docs/superpowers/2026-05-03-phase1-engine-domain-model-design.md`.

pub mod aa_set;
pub mod activation;
pub mod amino_acid;
pub mod enzyme;
pub mod instrument;
pub mod mass;
pub mod modification;
pub mod param_model;
pub mod peptide;
pub mod protocol;
pub mod spectrum;
pub mod tolerance;

// Convenience re-exports for the most-used types. Downstream crates
// (input/output/cli) prefer `use engine::Peptide` over the qualified path.
pub use aa_set::{AaSetError, AminoAcidSet, AminoAcidSetBuilder};
pub use activation::ActivationMethod;
pub use instrument::InstrumentType;
pub use amino_acid::AminoAcid;
pub use enzyme::Enzyme;
pub use mass::{nominal_from, C, H, H2O, N, O, PROTON, S, INTEGER_MASS_SCALER};
pub use modification::{ModLocation, ModParseError, Modification, ResidueSpec};
pub use param_model::{
    FragmentOffsetFrequency, IonType, Param, ParamParseError, Partition,
    PrecursorOffsetFrequency, SpecDataType,
};
pub use peptide::{Peptide, PeptideParseError};
pub use protocol::Protocol;
pub use spectrum::Spectrum;
pub use tolerance::{PrecursorTolerance, Tolerance};
