//! Input-side readers for MS-GF+ Rust port: spectrum file formats
//! (MGF in Phase 3a; mzML in Phase 3b) and `.fasta` (Phase 4a).

pub mod fasta;
pub mod mgf;
pub mod mzml;

pub use model::{Protein, ProteinDb, Spectrum};
pub use fasta::{FastaParseError, FastaReader};
pub use mgf::{MgfParseError, MgfReader};
pub use mzml::{MzMLParseError, MzMLReader};
