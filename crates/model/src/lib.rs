//! Domain model for MS-GF+ Rust port.
//!
//! Pure types: amino acids, modifications, peptides, enzymes,
//! tolerances, spectra, proteins, masses, activation, instrument,
//! protocol, compact FASTA. No I/O, no scoring.

pub mod aa_set;
pub mod activation;
pub mod amino_acid;
pub mod compact_fasta;
pub mod enzyme;
pub mod instrument;
pub mod mass;
pub mod modification;
pub mod peptide;
pub mod protein;
pub mod protocol;
pub mod spectrum;
pub mod tolerance;

// Convenience re-exports for the most-used types.
pub use aa_set::{AaSetError, AminoAcidSet, AminoAcidSetBuilder};
pub use activation::ActivationMethod;
pub use amino_acid::AminoAcid;
pub use compact_fasta::{CompactFastaError, CompactFastaSequence, ProteinAnnotation};
pub use enzyme::Enzyme;
pub use instrument::InstrumentType;
pub use mass::{nominal_from, H2O, PROTON};
pub use modification::{ModLocation, ModParseError, Modification, ResidueSpec};
pub use peptide::Peptide;
pub use protein::{Protein, ProteinDb};
pub use protocol::Protocol;
pub use spectrum::Spectrum;
pub use tolerance::{PrecursorTolerance, Tolerance};
