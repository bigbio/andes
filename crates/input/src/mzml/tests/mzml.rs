//! Tests of the mzML reader, one file per module of `mzml`, and the XML builders they share.

use std::io::Cursor;

use base64::{engine::general_purpose::STANDARD, Engine as _};
use byteorder::LittleEndian;

use model::Spectrum;

use crate::mzml::MzMLReader;

#[path = "activation.rs"]
mod activation;
#[path = "binary.rs"]
mod binary;
#[path = "centroid.rs"]
mod centroid;
#[path = "cv.rs"]
mod cv;
#[path = "instrument.rs"]
mod instrument;
#[path = "parallel.rs"]
mod parallel;
#[path = "reader.rs"]
mod reader;

fn collect_ok(xml: &str) -> Vec<Spectrum> {
    MzMLReader::new(Cursor::new(xml))
        .map(|r| r.expect("parse error"))
        .collect()
}

fn wrap_spectra(spectra: &str) -> String {
    format!(
        r#"<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml">
  <run>
    <spectrumList count="1" defaultDataProcessingRef="dp">
      {spectra}
    </spectrumList>
  </run>
</mzML>"#
    )
}

fn encode_f64_b64(vals: &[f64]) -> String {
    use byteorder::WriteBytesExt;
    let mut buf: Vec<u8> = Vec::with_capacity(vals.len() * 8);
    for &v in vals {
        buf.write_f64::<LittleEndian>(v).unwrap();
    }
    STANDARD.encode(&buf)
}

fn encode_f64_zlib_b64(vals: &[f64]) -> String {
    use byteorder::WriteBytesExt;
    use flate2::{write::ZlibEncoder, Compression};
    use std::io::Write;

    let mut raw: Vec<u8> = Vec::with_capacity(vals.len() * 8);
    for &v in vals {
        raw.write_f64::<LittleEndian>(v).unwrap();
    }
    let mut enc = ZlibEncoder::new(Vec::new(), Compression::default());
    enc.write_all(&raw).unwrap();
    STANDARD.encode(enc.finish().unwrap())
}

fn bda_block(cv_array: &str, compression_cv: &str, b64: &str) -> String {
    format!(
        r#"<binaryDataArray>
              <cvParam accession="MS:1000523" name="64-bit float" value=""/>
              <cvParam accession="{compression_cv}" name="compression" value=""/>
              <cvParam accession="{cv_array}" name="" value=""/>
              <binary>{b64}</binary>
            </binaryDataArray>"#
    )
}

fn bda_plain(cv_array: &str, b64: &str) -> String {
    bda_block(cv_array, "MS:1000576", b64)
}

fn bda_zlib(cv_array: &str, b64: &str) -> String {
    bda_block(cv_array, "MS:1000574", b64)
}

fn ms2_spectrum_xml(
    id: &str,
    mz_bda: &str,
    int_bda: &str,
    precursor_mz: f64,
    charge: Option<i32>,
) -> String {
    let charge_param = match charge {
        Some(z) => {
            format!(r#"<cvParam accession="MS:1000041" name="charge state" value="{z}"/>"#)
        }
        None => String::new(),
    };
    format!(
        r#"<spectrum index="0" id="{id}" defaultArrayLength="2">
              <cvParam accession="MS:1000511" name="ms level" value="2"/>
              <scanList count="1">
                <scan>
                  <cvParam accession="MS:1000016" name="scan start time" value="1.5"
                           unitAccession="UO:0000031" unitName="minute"/>
                </scan>
              </scanList>
              <precursorList count="1">
                <precursor>
                  <selectedIonList count="1">
                    <selectedIon>
                      <cvParam accession="MS:1000744" name="selected ion m/z"
                               value="{precursor_mz}"/>
                      {charge_param}
                    </selectedIon>
                  </selectedIonList>
                </precursor>
              </precursorList>
              <binaryDataArrayList count="2">
                {mz_bda}
                {int_bda}
              </binaryDataArrayList>
            </spectrum>"#
    )
}

fn malformed_spectrum() -> String {
    ms2_spectrum_xml(
        "scan=bad",
        &bda_plain("MS:1000514", "!!!not-base64!!!"),
        &bda_plain("MS:1000515", "!!!not-base64!!!"),
        500.5,
        Some(2),
    )
}
