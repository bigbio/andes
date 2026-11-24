//! Shared domain types used by every other crate: amino acids, modifications, peptides,
//! proteins, spectra, tolerances, instruments, enzymes and FDR helpers. It sits at the bottom of
//! the dependency graph and holds no pipeline stage of its own.
pub mod activation;
pub mod amino_acid;
pub mod annotation;
pub mod crosslinkers;
pub mod enzyme;
pub mod fdr;
pub mod glycans;
pub mod hash;
pub mod instrument;
pub mod isomers;
pub mod isotope;
pub mod mass;
pub mod modification;
pub mod peptide;
pub mod protein;
pub mod protocol;
pub mod reporter;
pub mod scan;
pub mod spectrum;
pub mod stats;
pub mod tolerance;
pub mod unimod;

pub use amino_acid::set::{AaSetError, AminoAcidSet, AminoAcidSetBuilder, Terminus};
pub use activation::ActivationMethod;
pub use amino_acid::AminoAcid;
pub use annotation::{Annotation, Edit, Processing};
pub use enzyme::Enzyme;
pub use instrument::InstrumentType;
pub use isomers::IsomerClass;
pub use mass::{nominal_from, H2O, PROTON};
pub use modification::spec::{ModDecl, ModSpecError};
pub use modification::{ModLocation, Modification, ResidueSpec};
pub use peptide::Peptide;
pub use protein::{
    ChainEnds, EntryKind, EntryOrigin, ModSite, PriorTier, Protein, ProteinAnnotation, ProteinDb,
    SiteAnchor, SiteClaim, Substitution, VariantRegion,
};
pub use protocol::Protocol;
pub use scan::{Ms1Scan, NoiseProfile, ProductScan};
pub use spectrum::Spectrum;
pub use tolerance::{PrecursorTolerance, Tolerance};
pub use unimod::{UnimodDb, UnimodMod, UnimodResolveError};
