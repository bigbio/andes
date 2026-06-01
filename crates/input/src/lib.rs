//! Input readers: MGF, mzML, FASTA, and (optional) Bruker timsTOF `.d`.

pub mod fasta;
pub mod mgf;
pub mod mzml;
#[cfg(feature = "timstof")]
pub mod timstof;

pub use model::{InstrumentType, Protein, ProteinDb, Spectrum};
pub use fasta::{FastaParseError, FastaReader};
pub use mgf::{MgfParseError, MgfReader};
pub use mzml::{detect_instrument_type, Ms1Link, MzMLParseError, MzMLReader};
#[cfg(feature = "timstof")]
pub use timstof::{TimsTofParseError, TimsTofReader};
