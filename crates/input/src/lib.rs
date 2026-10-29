//! Input readers: MGF, mzML, FASTA, UniProt XML, and (optional) Thermo `.raw` /
//! Bruker timsTOF `.d`.

pub mod acq_meta;
pub mod annotation;
pub mod annotations_file;
pub mod database;
pub mod fasta;
pub mod isobaric;
pub mod mgf;
pub mod mzml;
pub mod peff;
#[cfg(feature = "thermo")]
pub mod thermo;
#[cfg(feature = "timstof")]
pub mod timstof;
#[cfg(feature = "timstof")]
mod timstof_ms1;
pub mod uniprot;
mod uniprot_ptm;

pub use acq_meta::{
    infer_fragment_resolution, sample_mzml_ms2_meta, Ms2MetaSample, ResolutionEvidence,
};
pub use annotation::{
    add_site_knockoffs, build_entry, knockoff_annotations, KnockoffOptions, SiteKnockoffStats,
};
pub use annotations_file::{load_annotations, sequence_md5, AnnotationLoadStats, AnnotationsError};
pub use database::{load_protein_database, DatabaseError, DatabaseFormat, LoadedDatabase};
pub use fasta::{FastaParseError, FastaReader};
pub use isobaric::{detect_isobaric, IsobaricLabel};
pub use mgf::{MgfParseError, MgfReader};
pub use model::{InstrumentType, Protein, ProteinDb, Spectrum};
pub use mzml::{
    detect_instrument_type, read_spectra_parallel, Ms1Link, MzMLParseError, MzMLReader,
};
pub use peff::{load_peff, PeffError, PeffLoadStats};
#[cfg(feature = "thermo")]
pub use thermo::{ThermoParseError, ThermoRawReader};
#[cfg(feature = "timstof")]
pub use timstof::{TimsTofParseError, TimsTofReader};
pub use uniprot::{
    load_uniprot_xml, EvidenceLevel, UniProtLoadStats, UniProtOptions, UniProtParseError,
    UniProtXmlReader,
};

/// Open a spectrum/FASTA file for buffered reading, **transparently
/// gzip-decompressing** when the path ends in `.gz` (case-insensitive).
///
/// Returns a boxed `BufRead + Send` so the generic MGF / mzML / FASTA readers —
/// and the streaming parser thread — consume gzipped and plain inputs uniformly.
/// Callers that dispatch on file extension should strip a trailing `.gz` first to
/// detect the underlying format (e.g. `spectra.mgf.gz` → MGF). Every gzip member
/// is read, so concatenated (`cat a.gz b.gz`) and bgzip files decode in full.
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
            flate2::read::MultiGzDecoder::new(std::io::BufReader::new(file)),
        )))
    } else {
        Ok(Box::new(std::io::BufReader::new(file)))
    }
}

#[cfg(test)]
mod tests {
    use std::io::{Read, Write};

    #[test]
    fn multi_member_gzip_is_read_in_full() {
        let mut bytes = Vec::new();
        for part in ["BEGIN IONS\n", "END IONS\n"] {
            let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
            enc.write_all(part.as_bytes()).unwrap();
            bytes.extend(enc.finish().unwrap());
        }
        let dir = std::env::temp_dir().join(format!("andes-gz-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("two_members.mgf.gz");
        std::fs::write(&path, &bytes).unwrap();
        let mut text = String::new();
        super::open_buf_maybe_gz(&path)
            .unwrap()
            .read_to_string(&mut text)
            .unwrap();
        std::fs::remove_dir_all(&dir).ok();
        assert_eq!(text, "BEGIN IONS\nEND IONS\n");
    }
}
