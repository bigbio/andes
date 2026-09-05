//! Input readers: MGF, mzML, FASTA, and (optional) Thermo `.raw` / Bruker timsTOF `.d`.

pub mod fasta;
pub mod isobaric;
pub mod mgf;
pub mod mzml;
#[cfg(feature = "thermo")]
pub mod thermo;
#[cfg(feature = "timstof")]
pub mod timstof;

pub use fasta::{FastaParseError, FastaReader};
pub use isobaric::{detect_isobaric, IsobaricLabel};
pub use mgf::{MgfParseError, MgfReader};
pub use model::{InstrumentType, Protein, ProteinDb, Spectrum};
pub use mzml::{detect_instrument_type, Ms1Link, MzMLParseError, MzMLReader};
#[cfg(feature = "thermo")]
pub use thermo::{ThermoParseError, ThermoRawReader};
#[cfg(feature = "timstof")]
pub use timstof::{TimsTofParseError, TimsTofReader};

/// Open a spectrum/FASTA file for buffered reading, **transparently
/// gzip-decompressing** when the path ends in `.gz` (case-insensitive).
///
/// Returns a boxed `BufRead + Send` so the generic MGF / mzML / FASTA readers —
/// and the streaming parser thread — consume gzipped and plain inputs uniformly.
/// Callers that dispatch on file extension should strip a trailing `.gz` first to
/// detect the underlying format (e.g. `spectra.mgf.gz` → MGF).
pub fn open_buf_maybe_gz(
    path: &std::path::Path,
) -> std::io::Result<Box<dyn std::io::BufRead + Send>> {
    let file = std::fs::File::open(path)?;
    let is_gz = path
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.eq_ignore_ascii_case("gz"))
        .unwrap_or(false);
    if is_gz {
        Ok(Box::new(std::io::BufReader::new(
            flate2::read::GzDecoder::new(std::io::BufReader::new(file)),
        )))
    } else {
        Ok(Box::new(std::io::BufReader::new(file)))
    }
}
