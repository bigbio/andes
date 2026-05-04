//! Domain model for MS-GF+ Rust port.
//!
//! Phase 1 milestone: amino acids, modifications, peptides, enzymes,
//! tolerances. Pure CPU + types. No I/O. See
//! `docs/superpowers/2026-05-03-phase1-engine-domain-model-design.md`.

pub mod aa_set;
pub mod amino_acid;
pub mod enzyme;
pub mod mass;
pub mod modification;
pub mod peptide;
pub mod tolerance;

// Convenience re-exports for the most-used types. Downstream crates
// (input/output/cli) prefer `use engine::Peptide` over the qualified path.
pub use aa_set::{AaSetError, AminoAcidSet, AminoAcidSetBuilder};
pub use amino_acid::AminoAcid;
pub use enzyme::Enzyme;
pub use mass::{nominal_from, C, H, H2O, N, O, PROTON, S, INTEGER_MASS_SCALER};
pub use modification::{ModLocation, ModParseError, Modification, ResidueSpec};
pub use peptide::{Peptide, PeptideParseError};
pub use tolerance::{PrecursorTolerance, Tolerance};
