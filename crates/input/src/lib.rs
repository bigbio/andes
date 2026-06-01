//! Input readers: MGF, mzML, FASTA, and (optional) Thermo `.raw`.

pub mod fasta;
pub mod mgf;
pub mod mzml;
#[cfg(feature = "thermo")]
pub mod thermo;

pub use model::{InstrumentType, Protein, ProteinDb, Spectrum};
pub use fasta::{FastaParseError, FastaReader};
pub use mgf::{MgfParseError, MgfReader};
pub use mzml::{detect_instrument_type, Ms1Link, MzMLParseError, MzMLReader};
#[cfg(feature = "thermo")]
pub use thermo::{ThermoParseError, ThermoRawReader};
