//! Input readers: MGF, mzML, FASTA, and (optional) Thermo `.raw` / Bruker timsTOF `.d`.

pub mod fasta;
pub mod isobaric;
pub mod mgf;
pub mod mzml;
#[cfg(feature = "thermo")]
pub mod thermo;
#[cfg(feature = "timstof")]
pub mod timstof;

pub use model::{InstrumentType, Protein, ProteinDb, Spectrum};
pub use fasta::{FastaParseError, FastaReader};
pub use isobaric::{detect_isobaric, IsobaricLabel};
pub use mgf::{MgfParseError, MgfReader};
pub use mzml::{detect_instrument_type, Ms1Link, MzMLParseError, MzMLReader};
#[cfg(feature = "thermo")]
pub use thermo::{ThermoParseError, ThermoRawReader};
#[cfg(feature = "timstof")]
pub use timstof::{TimsTofParseError, TimsTofReader};
