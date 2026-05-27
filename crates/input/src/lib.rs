//! Input readers: MGF, mzML, FASTA.

pub mod fasta;
pub mod mgf;
pub mod mzml;

pub use model::{InstrumentType, Protein, ProteinDb, Spectrum};
pub use fasta::{FastaParseError, FastaReader};
pub use mgf::{MgfParseError, MgfReader};
pub use mzml::{detect_instrument_type, MzMLParseError, MzMLReader};
