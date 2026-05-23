//! Input-side readers for MS-GF+ Rust port: MGF and mzML spectrum files
//! and `.fasta` protein databases.

pub mod fasta;
pub mod mgf;
pub mod mzml;

pub use model::{InstrumentType, Protein, ProteinDb, Spectrum};
pub use fasta::{FastaParseError, FastaReader};
pub use mgf::{MgfParseError, MgfReader};
pub use mzml::{detect_instrument_type, MzMLParseError, MzMLReader};
