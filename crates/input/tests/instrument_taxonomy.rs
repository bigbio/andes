//! Taxonomy tests for the two new explicit instrument classes:
//! `InstrumentType::OrbitrapAstral` and `InstrumentType::TimsTOF`.
//!
//! These tests do NOT require real fixture files — they exercise the enum
//! variants and detection helpers directly.

use std::io::Cursor;

use input::mzml::detect_instrument_type;
use model::InstrumentType;

// ── Orbitrap Astral — mzML instrument-model cvParam ─────────────────────────

/// PSI-MS accession `MS:1003378` (Orbitrap Astral instrument model) must map
/// to `InstrumentType::OrbitrapAstral`, NOT the generic `QExactive`.
#[test]
fn instrument_taxonomy_astral_accession_ms1003378() {
    // Build a minimal mzML with the Orbitrap Astral instrument-model cvParam
    // at the IC level (not as an analyzer cvParam).
    let xml = wrap_with_model_cv_ic("IC1", "MS:1003378", "Orbitrap Astral");
    let result = detect_instrument_type(Cursor::new(xml));
    assert_eq!(
        result,
        Some(InstrumentType::OrbitrapAstral),
        "MS:1003378 must map to OrbitrapAstral"
    );
}

/// Name-substring fallback: an unknown accession whose `name` attribute
/// contains "astral" (case-insensitive) must also map to `OrbitrapAstral`.
#[test]
fn instrument_taxonomy_astral_name_fallback_case_insensitive() {
    // Use a placeholder accession that isn't MS:1003378 to exercise the
    // name-string fallback path.
    let xml = wrap_with_model_cv_ic("IC1", "MS:9999999", "Orbitrap ASTRAL HF");
    let result = detect_instrument_type(Cursor::new(xml));
    assert_eq!(
        result,
        Some(InstrumentType::OrbitrapAstral),
        "instrument name containing 'astral' must map to OrbitrapAstral"
    );
}

/// `OrbitrapAstral` is high-resolution.
#[test]
fn instrument_taxonomy_astral_is_high_resolution() {
    assert!(
        InstrumentType::OrbitrapAstral.is_high_resolution(),
        "OrbitrapAstral must be high-resolution"
    );
}

/// Round-trip through `name()` / `from_name()`.
#[test]
fn instrument_taxonomy_astral_name_round_trip() {
    let v = InstrumentType::OrbitrapAstral;
    assert_eq!(InstrumentType::from_name(v.name()), Some(v));
}

// ── TimsTOF ──────────────────────────────────────────────────────────────────

/// `TimsTOF` variant exists and is high-resolution.
#[test]
fn instrument_taxonomy_timstof_is_high_resolution() {
    assert!(
        InstrumentType::TimsTOF.is_high_resolution(),
        "TimsTOF must be high-resolution"
    );
}

/// Round-trip through `name()` / `from_name()`.
#[test]
fn instrument_taxonomy_timstof_name_round_trip() {
    let v = InstrumentType::TimsTOF;
    assert_eq!(InstrumentType::from_name(v.name()), Some(v));
}

// ── Regression: existing variants unchanged ───────────────────────────────────

/// Existing `QExactive` still maps from a plain Orbitrap analyzer cvParam
/// (`MS:1000484`) — `OrbitrapAstral` must not shadow the generic Orbitrap path.
#[test]
fn instrument_taxonomy_generic_orbitrap_still_qexactive() {
    let xml = wrap_with_analyzer_cv_ic("IC1", "MS:1000484");
    let result = detect_instrument_type(Cursor::new(xml));
    assert_eq!(result, Some(InstrumentType::QExactive));
}

// ── Helpers ───────────────────────────────────────────────────────────────────

/// Build a minimal mzML with one IC that carries an *instrument-model* cvParam
/// (at IC level, not inside `<analyzer>`). Mirrors the `detect_instrument_qexactive_model_cv_param`
/// test in mzml.rs.
fn wrap_with_model_cv_ic(ic_id: &str, accession: &str, name: &str) -> String {
    let mz_b64 = encode_f64_b64(&[100.0]);
    let int_b64 = encode_f64_b64(&[1000.0]);
    format!(
        r#"<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml">
  <instrumentConfigurationList count="1">
    <instrumentConfiguration id="{ic_id}">
      <cvParam accession="{accession}" name="{name}" value=""/>
      <componentList count="2">
        <source order="1"><cvParam accession="MS:1000398" name="nanoelectrospray" value=""/></source>
        <analyzer order="2"/>
      </componentList>
    </instrumentConfiguration>
  </instrumentConfigurationList>
  <run id="r" defaultInstrumentConfigurationRef="{ic_id}">
    <spectrumList count="1" defaultDataProcessingRef="dp">
      <spectrum index="0" id="scan=1" defaultArrayLength="1">
        <cvParam accession="MS:1000511" name="ms level" value="2"/>
        <scanList count="1"><scan instrumentConfigurationRef="{ic_id}"/></scanList>
        <precursorList count="1">
          <precursor>
            <selectedIonList count="1">
              <selectedIon>
                <cvParam accession="MS:1000744" name="selected ion m/z" value="500.5"/>
              </selectedIon>
            </selectedIonList>
          </precursor>
        </precursorList>
        <binaryDataArrayList count="2">
          <binaryDataArray encodedLength="16">
            <cvParam accession="MS:1000514" name="m/z array" value=""/>
            <cvParam accession="MS:1000521" name="32-bit float" value=""/>
            <binary>{mz_b64}</binary>
          </binaryDataArray>
          <binaryDataArray encodedLength="16">
            <cvParam accession="MS:1000515" name="intensity array" value=""/>
            <cvParam accession="MS:1000521" name="32-bit float" value=""/>
            <binary>{int_b64}</binary>
          </binaryDataArray>
        </binaryDataArrayList>
      </spectrum>
    </spectrumList>
  </run>
</mzML>"#
    )
}

/// Build a minimal mzML with one IC that carries an *analyzer* cvParam (inside
/// `<analyzer>`). Used for the regression test.
fn wrap_with_analyzer_cv_ic(ic_id: &str, analyzer_cv: &str) -> String {
    let mz_b64 = encode_f64_b64(&[100.0]);
    let int_b64 = encode_f64_b64(&[1000.0]);
    format!(
        r#"<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml">
  <instrumentConfigurationList count="1">
    <instrumentConfiguration id="{ic_id}">
      <componentList count="2">
        <source order="1"><cvParam accession="MS:1000398" name="nanoelectrospray" value=""/></source>
        <analyzer order="2">
          <cvParam accession="{analyzer_cv}" name="" value=""/>
        </analyzer>
      </componentList>
    </instrumentConfiguration>
  </instrumentConfigurationList>
  <run id="r" defaultInstrumentConfigurationRef="{ic_id}">
    <spectrumList count="1" defaultDataProcessingRef="dp">
      <spectrum index="0" id="scan=1" defaultArrayLength="1">
        <cvParam accession="MS:1000511" name="ms level" value="2"/>
        <scanList count="1"><scan instrumentConfigurationRef="{ic_id}"/></scanList>
        <precursorList count="1">
          <precursor>
            <selectedIonList count="1">
              <selectedIon>
                <cvParam accession="MS:1000744" name="selected ion m/z" value="500.5"/>
              </selectedIon>
            </selectedIonList>
          </precursor>
        </precursorList>
        <binaryDataArrayList count="2">
          <binaryDataArray encodedLength="16">
            <cvParam accession="MS:1000514" name="m/z array" value=""/>
            <cvParam accession="MS:1000521" name="32-bit float" value=""/>
            <binary>{mz_b64}</binary>
          </binaryDataArray>
          <binaryDataArray encodedLength="16">
            <cvParam accession="MS:1000515" name="intensity array" value=""/>
            <cvParam accession="MS:1000521" name="32-bit float" value=""/>
            <binary>{int_b64}</binary>
          </binaryDataArray>
        </binaryDataArrayList>
      </spectrum>
    </spectrumList>
  </run>
</mzML>"#
    )
}

fn encode_f64_b64(vals: &[f64]) -> String {
    use base64::Engine;
    let mut bytes = Vec::with_capacity(vals.len() * 8);
    for &v in vals {
        bytes.extend_from_slice(&v.to_le_bytes());
    }
    base64::engine::general_purpose::STANDARD.encode(&bytes)
}
