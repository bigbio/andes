//! Cross-linking search support: cross-linker definitions, peptide-pair candidates, decoys,
//! doublet handling, error control and roll-up of cross-linked PSMs to sites. It extends the
//! search stage for cross-linked samples and feeds its own FDR step.

pub mod crosslinker;
pub mod decoy;
pub mod doublet;
pub mod fdr;
pub mod pair;
pub use model::crosslinkers as reagents;
pub mod rollup;
pub mod sites;
pub mod spec;
