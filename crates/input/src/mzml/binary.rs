//! Decoding of one mzML `<binaryDataArray>`: base64, optional zlib, 32- or 64-bit
//! little-endian floats.

use base64::{engine::general_purpose::STANDARD, Engine as _};
use byteorder::{LittleEndian, ReadBytesExt};
use flate2::read::ZlibDecoder;

use super::MzMLParseError;

#[derive(Debug)]
pub(super) struct BinaryArrayCtx {
    pub(super) is_mz: bool,
    pub(super) is_intensity: bool,
    pub(super) is_noise_mz: bool,
    pub(super) is_noise: bool,
    pub(super) precision_bits: u8,
    pub(super) zlib: bool,
    pub(super) unsupported: Option<(&'static str, &'static str)>,
    pub(super) b64_text: String,
}

impl BinaryArrayCtx {
    fn new() -> Self {
        BinaryArrayCtx {
            is_mz: false,
            is_intensity: false,
            is_noise_mz: false,
            is_noise: false,
            precision_bits: 64,
            zlib: false,
            unsupported: None,
            b64_text: String::new(),
        }
    }
}

pub(super) fn decode_binary_array(ctx: &BinaryArrayCtx) -> Result<Vec<f64>, MzMLParseError> {
    if let Some((accession, name)) = ctx.unsupported {
        return Err(MzMLParseError::UnsupportedCompression { accession, name });
    }
    let trimmed = ctx.b64_text.trim();
    if trimmed.is_empty() {
        return Ok(Vec::new());
    }

    let raw = STANDARD.decode(trimmed)?;

    let bytes: Vec<u8> = if ctx.zlib {
        let mut decoder = ZlibDecoder::new(&raw[..]);
        let mut out = Vec::with_capacity(raw.len() * 2);
        std::io::Read::read_to_end(&mut decoder, &mut out).map_err(MzMLParseError::Zlib)?;
        out
    } else {
        raw
    };

    let elem_bytes = if ctx.precision_bits == 64 {
        8usize
    } else {
        4usize
    };
    if bytes.len() % elem_bytes != 0 {
        return Err(MzMLParseError::TruncatedBinary {
            decoded_len: bytes.len(),
            elem_bytes,
        });
    }

    let mut cur = std::io::Cursor::new(&bytes);
    let mut out: Vec<f64> = Vec::with_capacity(bytes.len() / elem_bytes);

    if ctx.precision_bits == 64 {
        while let Ok(v) = cur.read_f64::<LittleEndian>() {
            out.push(v);
        }
    } else {
        while let Ok(v) = cur.read_f32::<LittleEndian>() {
            out.push(v as f64);
        }
    }

    Ok(out)
}
