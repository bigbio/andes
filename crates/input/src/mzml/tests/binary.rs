//! Tests of binary array decoding: zlib, 64-bit, refused numpress, truncated arrays.

use std::io::Cursor;

use base64::{engine::general_purpose::STANDARD, Engine as _};

use super::*;
use crate::mzml::{MzMLParseError, MzMLReader};

#[test]
fn decodes_zlib_compressed_peaks() {
    let mz_vals = [150.0_f64, 300.0, 450.0];
    let int_vals = [2000.0_f64, 1000.0, 500.0];

    let spec = format!(
        r#"<spectrum index="0" id="scan=7" defaultArrayLength="3">
              <cvParam accession="MS:1000511" name="ms level" value="2"/>
              <scanList count="1"><scan/></scanList>
              <precursorList count="1">
                <precursor>
                  <selectedIonList count="1">
                    <selectedIon>
                      <cvParam accession="MS:1000744" name="selected ion m/z" value="500.0"/>
                    </selectedIon>
                  </selectedIonList>
                </precursor>
              </precursorList>
              <binaryDataArrayList count="2">
                {mz}
                {int}
              </binaryDataArrayList>
            </spectrum>"#,
        mz = bda_zlib("MS:1000514", &encode_f64_zlib_b64(&mz_vals)),
        int = bda_zlib("MS:1000515", &encode_f64_zlib_b64(&int_vals)),
    );
    let spectra = collect_ok(&wrap_spectra(&spec));

    assert_eq!(spectra.len(), 1);
    let peaks = &spectra[0].peaks;
    assert_eq!(peaks.len(), 3);
    assert!((peaks[0].0 - 150.0).abs() < 1e-6, "first m/z");
    assert!((peaks[1].0 - 300.0).abs() < 1e-6, "second m/z");
    assert!((peaks[2].0 - 450.0).abs() < 1e-6, "third m/z");
}

#[test]
fn decodes_uncompressed_64bit_peaks() {
    let mz_b64 = encode_f64_b64(&[200.0, 400.0]);
    let int_b64 = encode_f64_b64(&[5000.0, 2500.0]);

    let spec = ms2_spectrum_xml(
        "scan=3",
        &bda_plain("MS:1000514", &mz_b64),
        &bda_plain("MS:1000515", &int_b64),
        600.0,
        None,
    );
    let spectra = collect_ok(&wrap_spectra(&spec));

    assert_eq!(spectra.len(), 1);
    let peaks = &spectra[0].peaks;
    assert_eq!(peaks.len(), 2);
    assert!((peaks[0].0 - 200.0).abs() < 1e-6);
    assert!((peaks[1].0 - 400.0).abs() < 1e-6);
    assert!((peaks[0].1 - 5000.0_f32).abs() < 1.0);
    assert!((peaks[1].1 - 2500.0_f32).abs() < 1.0);
}

#[test]
fn numpress_array_is_an_error_not_garbage_peaks() {
    let b64 = encode_f64_b64(&[100.0]);
    let xml = ms2_spectrum_xml(
        "scan=1",
        &bda_block("MS:1000514", "MS:1002312", &b64),
        &bda_plain("MS:1000515", &encode_f64_b64(&[1000.0])),
        500.5,
        Some(2),
    );
    let mut reader = MzMLReader::new(Cursor::new(wrap_spectra(&xml)));
    assert_eq!(reader.by_ref().filter_map(Result::ok).count(), 0);
    assert_eq!(reader.skipped_count(), 1);
    let err = &reader.skipped_errors()[0];
    assert!(
        err.contains("MS:1002312") && err.contains("Numpress"),
        "{err}"
    );
}

#[test]
fn truncation_zlib_decodes_like_zlib() {
    let mz: Vec<u8> = [100.0f64, 200.0]
        .iter()
        .flat_map(|v| v.to_le_bytes())
        .collect();
    let mut enc = flate2::write::ZlibEncoder::new(Vec::new(), flate2::Compression::default());
    std::io::Write::write_all(&mut enc, &mz).unwrap();
    let b64 = STANDARD.encode(enc.finish().unwrap());
    let xml = ms2_spectrum_xml(
        "scan=1",
        &bda_block("MS:1000514", "MS:1003088", &b64),
        &bda_plain("MS:1000515", &encode_f64_b64(&[1.0, 2.0])),
        500.5,
        Some(2),
    );
    let spectra = collect_ok(&wrap_spectra(&xml));
    assert_eq!(spectra[0].peaks, vec![(100.0, 1.0), (200.0, 2.0)]);
}

#[test]
fn truncated_binary_array_is_rejected() {
    let truncated_b64 = STANDARD.encode([0u8; 12]);
    let spec = ms2_spectrum_xml(
        "scan=1",
        &bda_plain("MS:1000514", &truncated_b64),
        &bda_plain("MS:1000515", &encode_f64_b64(&[1000.0, 2000.0])),
        500.5,
        Some(2),
    );
    let mut reader = MzMLReader::new(Cursor::new(wrap_spectra(&spec))).strict();
    let first = reader.next().expect("should yield a result");
    assert!(
        matches!(
            first,
            Err(MzMLParseError::TruncatedBinary { elem_bytes: 8, .. })
        ),
        "truncated f64 array must be rejected, got {first:?}"
    );
}
