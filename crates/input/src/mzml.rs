//! mzML input: a streaming reader over the PSI mzML schema, a parallel reader that cuts one
//! file into byte pieces, and instrument detection from the run header.

mod activation;
mod binary;
mod centroid;
mod cv;
mod instrument;
mod parallel;
mod reader;

pub(crate) use activation::log_ethcd_once;
pub use activation::{ethcd_as_etd_policy, filter_has_supplemental_collision, init_ethcd_as_etd};
pub use instrument::detect_instrument_type;
pub use parallel::{read_ms1_linked_parallel, read_spectra_parallel, Ms1LinkedRead};
pub use reader::{Ms1Link, MzMLReader, MAX_SKIPPED_ERRORS};

#[derive(Debug, thiserror::Error)]
pub enum MzMLParseError {
    #[error("XML parse error: {0}")]
    Xml(#[from] quick_xml::Error),

    #[error("base64 decode error: {0}")]
    Base64(#[from] base64::DecodeError),

    #[error("zlib decode error: {0}")]
    Zlib(std::io::Error),

    #[error("mzML structure: {0}")]
    Structure(String),

    #[error("mismatched binary array lengths: m/z {mz_len} vs intensity {int_len}")]
    LengthMismatch { mz_len: usize, int_len: usize },

    #[error(
        "unsupported binary encoding {accession} ({name}); re-convert without it \
         (e.g. msconvert without --numpress*/--zstd)"
    )]
    UnsupportedCompression {
        accession: &'static str,
        name: &'static str,
    },

    #[error("truncated binary array: decoded {decoded_len} bytes is not a multiple of the {elem_bytes}-byte element width")]
    TruncatedBinary {
        decoded_len: usize,
        elem_bytes: usize,
    },
}

impl From<std::io::Error> for MzMLParseError {
    fn from(e: std::io::Error) -> Self {
        MzMLParseError::Zlib(e)
    }
}

fn attr_str(e: &quick_xml::events::BytesStart<'_>, name: &[u8]) -> Option<String> {
    e.attributes()
        .filter_map(|a| a.ok())
        .find(|a| a.key.local_name().as_ref() == name)
        .and_then(|a| {
            std::str::from_utf8(a.value.as_ref())
                .ok()
                .map(str::to_owned)
        })
}

#[cfg(test)]
#[path = "mzml/tests/mzml.rs"]
mod tests;
