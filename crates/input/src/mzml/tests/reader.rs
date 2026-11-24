//! Tests of the streaming reader: spectrum fields, precursor and isolation window, peak
//! filtering, MS1/MS3/noise capture, chunked reading and recovery from malformed spectra.

use std::io::Cursor;

use model::scan::ScanContext;
use model::Spectrum;

use super::*;
use crate::mzml::reader::{extract_scan_from_id, scan_number};
use crate::mzml::MzMLReader;

#[test]
fn parses_minimal_mzml_with_one_ms2_spectrum() {
    let mz_b64 = encode_f64_b64(&[100.0, 200.0]);
    let int_b64 = encode_f64_b64(&[1000.0, 500.0]);

    let spec = ms2_spectrum_xml(
        "scan=1",
        &bda_plain("MS:1000514", &mz_b64),
        &bda_plain("MS:1000515", &int_b64),
        500.5,
        Some(2),
    );
    let spectra = collect_ok(&wrap_spectra(&spec));

    assert_eq!(spectra.len(), 1, "expected exactly one MS2 spectrum");
    assert_eq!(spectra[0].peaks.len(), 2, "expected two peaks");
}

#[test]
fn filters_out_ms1_spectra() {
    let mz_b64 = encode_f64_b64(&[100.0]);
    let int_b64 = encode_f64_b64(&[100.0]);

    let ms1 = format!(
        r#"<spectrum index="0" id="scan=1" defaultArrayLength="1">
              <cvParam accession="MS:1000511" name="ms level" value="1"/>
              <scanList count="1"><scan/></scanList>
              <binaryDataArrayList count="2">
                {mz}
                {int}
              </binaryDataArrayList>
            </spectrum>"#,
        mz = bda_plain("MS:1000514", &mz_b64),
        int = bda_plain("MS:1000515", &int_b64),
    );

    let ms2_mz_b64 = encode_f64_b64(&[200.0, 300.0]);
    let ms2_int_b64 = encode_f64_b64(&[800.0, 400.0]);
    let ms2 = ms2_spectrum_xml(
        "scan=2",
        &bda_plain("MS:1000514", &ms2_mz_b64),
        &bda_plain("MS:1000515", &ms2_int_b64),
        500.0,
        Some(2),
    );

    let xml = format!(
        r#"<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml">
  <run>
    <spectrumList count="2" defaultDataProcessingRef="dp">
      {ms1}
      {ms2}
    </spectrumList>
  </run>
</mzML>"#
    );

    let spectra = collect_ok(&xml);
    assert_eq!(spectra.len(), 1, "only the MS2 should be emitted");
    assert_eq!(spectra[0].scan, Some(2));
}

#[test]
fn extracts_scan_number_from_id_attr() {
    let mz_b64 = encode_f64_b64(&[100.0]);
    let int_b64 = encode_f64_b64(&[1000.0]);

    let spec = ms2_spectrum_xml(
        "controllerType=0 controllerNumber=1 scan=1234",
        &bda_plain("MS:1000514", &mz_b64),
        &bda_plain("MS:1000515", &int_b64),
        500.0,
        None,
    );
    let spectra = collect_ok(&wrap_spectra(&spec));

    assert_eq!(spectra.len(), 1);
    assert_eq!(spectra[0].scan, Some(1234));
}

#[test]
fn extracts_precursor_mz_and_charge() {
    let mz_b64 = encode_f64_b64(&[100.0]);
    let int_b64 = encode_f64_b64(&[1000.0]);

    let spec = ms2_spectrum_xml(
        "scan=10",
        &bda_plain("MS:1000514", &mz_b64),
        &bda_plain("MS:1000515", &int_b64),
        500.5,
        Some(2),
    );
    let spectra = collect_ok(&wrap_spectra(&spec));

    assert_eq!(spectra.len(), 1);
    assert!((spectra[0].precursor_mz - 500.5).abs() < 1e-6);
    assert_eq!(spectra[0].precursor_charge, Some(2));
}

#[test]
fn first_selected_ion_gives_mz_and_charge_together() {
    let mz_b64 = encode_f64_b64(&[100.0]);
    let int_b64 = encode_f64_b64(&[1000.0]);
    let ions = |second: &str| {
        format!(
            r#"<spectrum index="0" id="scan=7" defaultArrayLength="1">
              <cvParam accession="MS:1000511" name="ms level" value="2"/>
              <scanList count="1">
                <scan>
                  <cvParam accession="MS:1000016" name="scan start time" value="1.5"
                           unitAccession="UO:0000031" unitName="minute"/>
                </scan>
              </scanList>
              <precursorList count="1">
                <precursor>
                  <selectedIonList count="2">
                    <selectedIon>
                      <cvParam accession="MS:1000744" name="selected ion m/z" value="500.25"/>
                      <cvParam accession="MS:1000041" name="charge state" value="2"/>
                      <cvParam accession="MS:1000042" name="peak intensity" value="1000"/>
                    </selectedIon>
                    <selectedIon>{second}</selectedIon>
                  </selectedIonList>
                </precursor>
              </precursorList>
              <binaryDataArrayList count="2">
                {}
                {}
              </binaryDataArrayList>
            </spectrum>"#,
            bda_plain("MS:1000514", &mz_b64),
            bda_plain("MS:1000515", &int_b64),
        )
    };
    for second in [
        r#"<cvParam accession="MS:1000744" name="selected ion m/z" value="733.70"/>"#,
        r#"<cvParam accession="MS:1000744" name="selected ion m/z" value="733.70"/>
               <cvParam accession="MS:1000041" name="charge state" value="3"/>
               <cvParam accession="MS:1000042" name="peak intensity" value="50"/>"#,
        r#"<cvParam accession="MS:1000041" name="charge state" value="4"/>"#,
    ] {
        let spectra = collect_ok(&wrap_spectra(&ions(second)));
        assert_eq!(spectra.len(), 1);
        assert!((spectra[0].precursor_mz - 500.25).abs() < 1e-9, "{second}");
        assert_eq!(spectra[0].precursor_charge, Some(2), "{second}");
        assert_eq!(spectra[0].precursor_intensity, Some(1000.0), "{second}");
    }
}

#[test]
fn out_of_range_charge_state_is_treated_as_absent() {
    let mz_b64 = encode_f64_b64(&[100.0]);
    let int_b64 = encode_f64_b64(&[1000.0]);
    for z in [300, -2] {
        let spec = ms2_spectrum_xml(
            "scan=10",
            &bda_plain("MS:1000514", &mz_b64),
            &bda_plain("MS:1000515", &int_b64),
            500.5,
            Some(z),
        );
        let spectra = collect_ok(&wrap_spectra(&spec));
        assert_eq!(spectra[0].precursor_charge, None, "charge {z}");
    }
}

#[test]
fn charge_state_zero_is_treated_as_absent() {
    let mz_b64 = encode_f64_b64(&[100.0]);
    let int_b64 = encode_f64_b64(&[1000.0]);
    let spec = ms2_spectrum_xml(
        "scan=10",
        &bda_plain("MS:1000514", &mz_b64),
        &bda_plain("MS:1000515", &int_b64),
        500.5,
        Some(0),
    );
    let spectra = collect_ok(&wrap_spectra(&spec));
    assert_eq!(spectra.len(), 1);
    assert_eq!(spectra[0].precursor_charge, None);
}

#[test]
fn non_finite_and_nonpositive_peaks_are_filtered() {
    let mz_b64 = encode_f64_b64(&[f64::NAN, 100.0, -50.0, f64::INFINITY, 200.0]);
    let int_b64 = encode_f64_b64(&[10.0, 1000.0, 10.0, 10.0, 500.0]);
    let spec = ms2_spectrum_xml(
        "scan=11",
        &bda_plain("MS:1000514", &mz_b64),
        &bda_plain("MS:1000515", &int_b64),
        500.5,
        Some(2),
    );
    let spectra = collect_ok(&wrap_spectra(&spec));
    assert_eq!(spectra.len(), 1);
    let mzs: Vec<f64> = spectra[0].peaks.iter().map(|&(mz, _)| mz).collect();
    assert_eq!(
        mzs,
        vec![100.0, 200.0],
        "only finite positive-m/z peaks survive"
    );
}

#[test]
fn nonpositive_precursor_mz_is_skipped() {
    let mz_b64 = encode_f64_b64(&[100.0]);
    let int_b64 = encode_f64_b64(&[1000.0]);
    let spec = ms2_spectrum_xml(
        "scan=12",
        &bda_plain("MS:1000514", &mz_b64),
        &bda_plain("MS:1000515", &int_b64),
        0.0,
        Some(2),
    );
    let spectra = collect_ok(&wrap_spectra(&spec));
    assert!(
        spectra.is_empty(),
        "spectrum with precursor m/z 0 must be skipped"
    );
}

#[test]
fn peaks_sorted_ascending_by_mz() {
    let mz_b64 = encode_f64_b64(&[300.0, 100.0, 200.0]);
    let int_b64 = encode_f64_b64(&[3.0, 1.0, 2.0]);

    let spec = format!(
        r#"<spectrum index="0" id="scan=5" defaultArrayLength="3">
              <cvParam accession="MS:1000511" name="ms level" value="2"/>
              <scanList count="1"><scan/></scanList>
              <precursorList count="1">
                <precursor>
                  <selectedIonList count="1">
                    <selectedIon>
                      <cvParam accession="MS:1000744" name="selected ion m/z" value="600.0"/>
                    </selectedIon>
                  </selectedIonList>
                </precursor>
              </precursorList>
              <binaryDataArrayList count="2">
                {mz}
                {int}
              </binaryDataArrayList>
            </spectrum>"#,
        mz = bda_plain("MS:1000514", &mz_b64),
        int = bda_plain("MS:1000515", &int_b64),
    );
    let spectra = collect_ok(&wrap_spectra(&spec));

    assert_eq!(spectra.len(), 1);
    let mzs: Vec<f64> = spectra[0].peaks.iter().map(|p| p.0).collect();
    assert_eq!(mzs, vec![100.0, 200.0, 300.0]);
}

#[test]
fn parses_real_test_fixture() {
    let fixture =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../test-fixtures/tiny.pwiz.mzML");

    let file = std::fs::File::open(&fixture).expect("failed to open tiny.pwiz.mzML");
    let spectra: Vec<Spectrum> = MzMLReader::new(std::io::BufReader::new(file))
        .map(|r| r.expect("parse error"))
        .collect();

    assert!(!spectra.is_empty(), "expected at least one MS2 spectrum");
    let s = &spectra[0];
    assert!(!s.peaks.is_empty(), "MS2 spectrum should have peaks");
    assert!(s.precursor_mz > 0.0, "precursor m/z should be positive");
}

#[test]
fn extract_scan_plain() {
    assert_eq!(extract_scan_from_id("scan=1234"), Some(1234));
}

#[test]
fn extract_scan_pwiz_format() {
    assert_eq!(
        extract_scan_from_id("controllerType=0 controllerNumber=1 scan=42"),
        Some(42)
    );
}

#[test]
fn thermo_trailer_monoisotopic_overrides_selected_ion_mz() {
    let mz_b64 = encode_f64_b64(&[100.0]);
    let int_b64 = encode_f64_b64(&[1000.0]);
    let xml = format!(
        r#"<spectrum index="0" id="scan=1" defaultArrayLength="1">
              <cvParam accession="MS:1000511" name="ms level" value="2"/>
              <scanList count="1">
                <scan>
                  <cvParam accession="MS:1000016" name="scan start time"
                           value="1.5" unitAccession="UO:0000031"
                           unitName="minute"/>
                  <userParam name="[Thermo Trailer Extra]Monoisotopic M/Z:"
                             type="xsd:float" value="625.0037"/>
                </scan>
              </scanList>
              <precursorList count="1">
                <precursor>
                  <selectedIonList count="1">
                    <selectedIon>
                      <cvParam accession="MS:1000744" name="selected ion m/z"
                               value="625.338134765625"/>
                      <cvParam accession="MS:1000041" name="charge state" value="3"/>
                    </selectedIon>
                  </selectedIonList>
                </precursor>
              </precursorList>
              <binaryDataArrayList count="2">
                {}
                {}
              </binaryDataArrayList>
            </spectrum>"#,
        bda_plain("MS:1000514", &mz_b64),
        bda_plain("MS:1000515", &int_b64),
    );
    let spectra = collect_ok(&wrap_spectra(&xml));
    assert_eq!(spectra.len(), 1);
    assert!(
        (spectra[0].precursor_mz - 625.0037).abs() < 1e-6,
        "expected Thermo trailer monoisotopic m/z (625.0037), got {}",
        spectra[0].precursor_mz
    );
}

#[test]
fn precursor_mz_falls_back_to_selected_ion_without_trailer() {
    let mz_b64 = encode_f64_b64(&[100.0]);
    let int_b64 = encode_f64_b64(&[1000.0]);
    let spec = ms2_spectrum_xml(
        "scan=42",
        &bda_plain("MS:1000514", &mz_b64),
        &bda_plain("MS:1000515", &int_b64),
        500.5,
        Some(2),
    );
    let spectra = collect_ok(&wrap_spectra(&spec));
    assert_eq!(spectra.len(), 1);
    assert!((spectra[0].precursor_mz - 500.5).abs() < 1e-6);
}

#[test]
fn zero_thermo_trailer_does_not_override() {
    let mz_b64 = encode_f64_b64(&[100.0]);
    let int_b64 = encode_f64_b64(&[1000.0]);
    let xml = format!(
        r#"<spectrum index="0" id="scan=1" defaultArrayLength="1">
              <cvParam accession="MS:1000511" name="ms level" value="2"/>
              <scanList count="1">
                <scan>
                  <userParam name="[Thermo Trailer Extra]Monoisotopic M/Z:"
                             type="xsd:float" value="0"/>
                </scan>
              </scanList>
              <precursorList count="1">
                <precursor>
                  <selectedIonList count="1">
                    <selectedIon>
                      <cvParam accession="MS:1000744" name="selected ion m/z"
                               value="700.25"/>
                    </selectedIon>
                  </selectedIonList>
                </precursor>
              </precursorList>
              <binaryDataArrayList count="2">
                {}
                {}
              </binaryDataArrayList>
            </spectrum>"#,
        bda_plain("MS:1000514", &mz_b64),
        bda_plain("MS:1000515", &int_b64),
    );
    let spectra = collect_ok(&wrap_spectra(&xml));
    assert_eq!(spectra.len(), 1);
    assert!(
        (spectra[0].precursor_mz - 700.25).abs() < 1e-6,
        "zero-trailer must fall back to selectedIon m/z; got {}",
        spectra[0].precursor_mz
    );
}

#[test]
fn extract_scan_missing() {
    assert_eq!(extract_scan_from_id("spectrum=1"), None);
    assert_eq!(extract_scan_from_id("index=101"), Some(101));
    assert_eq!(extract_scan_from_id("scanId=7"), Some(7));
    assert_eq!(
        extract_scan_from_id("controllerType=0 controllerNumber=1 scan=42 index=9"),
        Some(42)
    );
    assert_eq!(
        scan_number("sample=1 period=1 cycle=4 experiment=2", Some(16)),
        Some(17)
    );
    assert_eq!(scan_number("index=101", Some(0)), Some(101));
}

#[test]
fn parses_isolation_window_offsets() {
    let mz_b64 = encode_f64_b64(&[100.0]);
    let int_b64 = encode_f64_b64(&[1000.0]);
    let xml = format!(
        r#"<spectrum index="0" id="scan=1" defaultArrayLength="1">
              <cvParam accession="MS:1000511" name="ms level" value="2"/>
              <scanList count="1"><scan/></scanList>
              <precursorList count="1">
                <precursor>
                  <isolationWindow>
                    <cvParam accession="MS:1000827" name="isolation window target m/z"
                             value="500.5"/>
                    <cvParam accession="MS:1000828" name="isolation window lower offset"
                             value="1.5"/>
                    <cvParam accession="MS:1000829" name="isolation window upper offset"
                             value="1.5"/>
                  </isolationWindow>
                  <selectedIonList count="1">
                    <selectedIon>
                      <cvParam accession="MS:1000744" name="selected ion m/z"
                               value="500.5"/>
                    </selectedIon>
                  </selectedIonList>
                </precursor>
              </precursorList>
              <binaryDataArrayList count="2">
                {mz}
                {int}
              </binaryDataArrayList>
            </spectrum>"#,
        mz = bda_plain("MS:1000514", &mz_b64),
        int = bda_plain("MS:1000515", &int_b64),
    );
    let spectra = collect_ok(&wrap_spectra(&xml));
    assert_eq!(spectra.len(), 1);
    assert_eq!(spectra[0].isolation_lower_offset, Some(1.5));
    assert_eq!(spectra[0].isolation_upper_offset, Some(1.5));
    assert_eq!(spectra[0].isolation_target_mz, Some(500.5));
    assert!((spectra[0].precursor_mz - 500.5).abs() < 1e-6);
}

#[test]
fn product_isolation_window_does_not_overwrite_precursor_window() {
    let mz_b64 = encode_f64_b64(&[100.0]);
    let int_b64 = encode_f64_b64(&[1000.0]);
    let xml = format!(
        r#"<spectrum index="0" id="scan=1" defaultArrayLength="1">
              <cvParam accession="MS:1000511" name="ms level" value="2"/>
              <scanList count="1"><scan/></scanList>
              <precursorList count="1">
                <precursor>
                  <isolationWindow>
                    <cvParam accession="MS:1000827" name="isolation window target m/z" value="500.5"/>
                    <cvParam accession="MS:1000828" name="isolation window lower offset" value="0.7"/>
                    <cvParam accession="MS:1000829" name="isolation window upper offset" value="0.7"/>
                  </isolationWindow>
                  <selectedIonList count="1">
                    <selectedIon>
                      <cvParam accession="MS:1000744" name="selected ion m/z" value="500.5"/>
                    </selectedIon>
                  </selectedIonList>
                </precursor>
              </precursorList>
              <productList count="1">
                <product>
                  <isolationWindow>
                    <cvParam accession="MS:1000827" name="isolation window target m/z" value="800.0"/>
                    <cvParam accession="MS:1000828" name="isolation window lower offset" value="400.0"/>
                    <cvParam accession="MS:1000829" name="isolation window upper offset" value="400.0"/>
                  </isolationWindow>
                </product>
              </productList>
              <binaryDataArrayList count="2">
                {mz}
                {int}
              </binaryDataArrayList>
            </spectrum>"#,
        mz = bda_plain("MS:1000514", &mz_b64),
        int = bda_plain("MS:1000515", &int_b64),
    );
    let spectra = collect_ok(&wrap_spectra(&xml));
    assert_eq!(spectra.len(), 1);
    assert_eq!(spectra[0].isolation_lower_offset, Some(0.7));
    assert_eq!(spectra[0].isolation_upper_offset, Some(0.7));
    assert_eq!(spectra[0].isolation_target_mz, Some(500.5));
    assert_eq!(spectra[0].peaks.len(), 1);
}

#[test]
fn missing_isolation_window_yields_none() {
    let mz_b64 = encode_f64_b64(&[100.0]);
    let int_b64 = encode_f64_b64(&[1000.0]);
    let spec = ms2_spectrum_xml(
        "scan=1",
        &bda_plain("MS:1000514", &mz_b64),
        &bda_plain("MS:1000515", &int_b64),
        500.5,
        Some(2),
    );
    let spectra = collect_ok(&wrap_spectra(&spec));
    assert_eq!(spectra.len(), 1);
    assert_eq!(spectra[0].isolation_lower_offset, None);
    assert_eq!(spectra[0].isolation_upper_offset, None);
    assert_eq!(spectra[0].isolation_target_mz, None);
}

#[test]
fn ms1_capture_links_two_ms2_to_preceding_ms1() {
    let ms1_mz_b64 = encode_f64_b64(&[111.0, 222.0]);
    let ms1_int_b64 = encode_f64_b64(&[10.0, 20.0]);
    let ms1 = format!(
        r#"<spectrum index="0" id="scan=1" defaultArrayLength="2">
              <cvParam accession="MS:1000511" name="ms level" value="1"/>
              <scanList count="1"><scan/></scanList>
              <binaryDataArrayList count="2">
                {mz}
                {int}
              </binaryDataArrayList>
            </spectrum>"#,
        mz = bda_plain("MS:1000514", &ms1_mz_b64),
        int = bda_plain("MS:1000515", &ms1_int_b64),
    );

    let ms2a = ms2_spectrum_xml(
        "scan=2",
        &bda_plain("MS:1000514", &encode_f64_b64(&[300.0, 400.0])),
        &bda_plain("MS:1000515", &encode_f64_b64(&[1.0, 2.0])),
        500.0,
        Some(2),
    );
    let ms2b = ms2_spectrum_xml(
        "scan=3",
        &bda_plain("MS:1000514", &encode_f64_b64(&[500.0])),
        &bda_plain("MS:1000515", &encode_f64_b64(&[3.0])),
        600.0,
        Some(2),
    );

    let xml = format!(
        r#"<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml">
  <run>
    <spectrumList count="3" defaultDataProcessingRef="dp">
      {ms1}
      {ms2a}
      {ms2b}
    </spectrumList>
  </run>
</mzML>"#
    );

    let (spectra, link) = MzMLReader::new(Cursor::new(xml))
        .with_ms1_capture(true)
        .read_with_ms1()
        .expect("read_with_ms1 failed");

    assert_eq!(spectra.len(), 2, "expected exactly two MS2 spectra");
    assert_eq!(spectra[0].scan, Some(2));
    assert_eq!(spectra[1].scan, Some(3));

    assert_eq!(link.ms1_peaks.len(), 1, "expected one captured MS1");
    let ms1_peaks = &link.ms1_peaks[0];
    assert_eq!(ms1_peaks.len(), 2);
    assert!((ms1_peaks[0].0 - 111.0).abs() < 1e-6);
    assert!((ms1_peaks[1].0 - 222.0).abs() < 1e-6);
    assert!((ms1_peaks[0].1 - 10.0_f32).abs() < 1e-3);
    assert!((ms1_peaks[1].1 - 20.0_f32).abs() < 1e-3);

    assert_eq!(link.ms2_to_ms1, vec![Some(0), Some(0)]);
}

#[test]
fn run_ms1_capture_keeps_rt_and_peaks_without_emitting_ms1() {
    let ms1 = |id: &str, rt_min: f64, mz: f64| {
        format!(
            r#"<spectrum index="0" id="{id}" defaultArrayLength="1">
              <cvParam accession="MS:1000511" name="ms level" value="1"/>
              <scanList count="1"><scan>
                <cvParam accession="MS:1000016" name="scan start time" value="{rt_min}"
                         unitAccession="UO:0000031" unitName="minute"/>
              </scan></scanList>
              <binaryDataArrayList count="2">{mzb}{intb}</binaryDataArrayList>
            </spectrum>"#,
            mzb = bda_plain("MS:1000514", &encode_f64_b64(&[mz])),
            intb = bda_plain("MS:1000515", &encode_f64_b64(&[7.0])),
        )
    };
    let ms2 = ms2_spectrum_xml(
        "scan=2",
        &bda_plain("MS:1000514", &encode_f64_b64(&[300.0])),
        &bda_plain("MS:1000515", &encode_f64_b64(&[1.0])),
        350.0,
        Some(2),
    );
    let xml = wrap_spectra(&format!(
        "{}{}{}",
        ms1("scan=1", 1.0, 400.0),
        ms2,
        ms1("scan=3", 1.5, 401.0)
    ));
    let mut reader = MzMLReader::new(Cursor::new(xml)).with_run_ms1_capture(true);
    let spectra: Vec<Spectrum> = reader.by_ref().map(|r| r.expect("parse")).collect();
    assert_eq!(spectra.len(), 1, "only the MS2 is emitted");
    let ms1s = reader.take_run_ms1();
    assert_eq!(ms1s.len(), 2);
    assert_eq!(ms1s[0].rt, 60.0);
    assert_eq!(ms1s[1].rt, 90.0);
    assert_eq!(ms1s[1].peaks, vec![(401.0, 7.0)]);
    assert!(reader.take_run_ms1().is_empty(), "take drains");
}

#[test]
fn scan_context_records_survey_reference_scan_window_and_faims_voltage() {
    let spectrum = |id: &str, level: u32, extra: &str| {
        format!(
            r#"<spectrum index="0" id="{id}" defaultArrayLength="1">
              <cvParam accession="MS:1000511" name="ms level" value="{level}"/>
              <scanList count="1"><scan>
                <cvParam accession="MS:1000016" name="scan start time" value="1.0"
                         unitAccession="UO:0000031" unitName="minute"/>
                <cvParam accession="MS:1001581" name="FAIMS compensation voltage" value="-45"/>
                <scanWindowList count="1"><scanWindow>
                  <cvParam accession="MS:1000501" name="scan window lower limit" value="148.0"/>
                  <cvParam accession="MS:1000500" name="scan window upper limit" value="1866.0"/>
                </scanWindow></scanWindowList>
              </scan></scanList>
              {extra}
              <binaryDataArrayList count="2">{mzb}{intb}</binaryDataArrayList>
            </spectrum>"#,
            mzb = bda_plain("MS:1000514", &encode_f64_b64(&[400.0])),
            intb = bda_plain("MS:1000515", &encode_f64_b64(&[7.0])),
        )
    };
    let precursor = r#"<precursorList count="1"><precursor spectrumRef="scan=1">
                  <selectedIonList count="1"><selectedIon>
                    <cvParam accession="MS:1000744" name="selected ion m/z" value="500.5"/>
                    <cvParam accession="MS:1000041" name="charge state" value="2"/>
                  </selectedIon></selectedIonList>
                </precursor></precursorList>"#;
    let xml = wrap_spectra(&format!(
        "{}{}",
        spectrum("scan=1", 1, ""),
        spectrum("scan=2", 2, precursor)
    ));
    let mut reader = MzMLReader::new(Cursor::new(xml.clone()))
        .with_run_ms1_capture(true)
        .with_scan_context_capture(true);
    let spectra: Vec<Spectrum> = reader.by_ref().map(|r| r.expect("parse")).collect();
    assert_eq!(spectra.len(), 1);
    let ms1 = reader.take_run_ms1();
    assert_eq!(ms1[0].id, "scan=1");
    assert_eq!(ms1[0].faims_cv, Some(-45.0));
    let contexts = reader.take_scan_contexts();
    assert_eq!(
        contexts,
        vec![ScanContext {
            id: "scan=2".into(),
            parent_id: Some("scan=1".into()),
            scan_window: Some((148.0, 1866.0)),
            faims_cv: Some(-45.0),
            noise: None,
            ion_mobility: None,
        }]
    );
    let mut plain = MzMLReader::new(Cursor::new(xml));
    assert_eq!(plain.by_ref().count(), 1);
    assert!(plain.take_scan_contexts().is_empty());
}

#[test]
fn product_scan_capture_links_ms3_to_its_ms2_and_hides_it() {
    let ms2 = ms2_spectrum_xml(
        "controllerType=0 controllerNumber=1 scan=10",
        &bda_plain("MS:1000514", &encode_f64_b64(&[300.0])),
        &bda_plain("MS:1000515", &encode_f64_b64(&[1.0])),
        350.0,
        Some(2),
    );
    let ms3 = format!(
        r#"<spectrum index="1" id="controllerType=0 controllerNumber=1 scan=11" defaultArrayLength="2">
              <cvParam accession="MS:1000511" name="ms level" value="3"/>
              <scanList count="1"><scan>
                <cvParam accession="MS:1000016" name="scan start time" value="2.0"
                         unitAccession="UO:0000031" unitName="minute"/>
              </scan></scanList>
              <precursorList count="2">
                <precursor spectrumRef="controllerType=0 controllerNumber=1 scan=10">
                  <selectedIonList count="1"><selectedIon>
                    <cvParam accession="MS:1000744" name="selected ion m/z" value="500.5"/>
                  </selectedIon></selectedIonList>
                  <activation><cvParam accession="MS:1000422" name="HCD" value=""/></activation>
                </precursor>
                <precursor spectrumRef="controllerType=0 controllerNumber=1 scan=10">
                  <selectedIonList count="1"><selectedIon>
                    <cvParam accession="MS:1000744" name="selected ion m/z" value="612.2"/>
                  </selectedIon></selectedIonList>
                </precursor>
              </precursorList>
              <binaryDataArrayList count="2">{mzb}{intb}</binaryDataArrayList>
            </spectrum>"#,
        mzb = bda_plain("MS:1000514", &encode_f64_b64(&[126.1277, 127.1248])),
        intb = bda_plain("MS:1000515", &encode_f64_b64(&[1000.0, 2000.0])),
    );
    let xml = wrap_spectra(&format!("{ms2}{ms3}"));
    let mut plain = MzMLReader::new(Cursor::new(xml.clone()));
    let n = plain
        .by_ref()
        .collect::<Result<Vec<_>, _>>()
        .expect("parse")
        .len();
    assert_eq!(n, 1);
    assert!(plain.take_product_scans().is_empty());
    let mut reader = MzMLReader::new(Cursor::new(xml)).with_product_scan_capture(true);
    let spectra: Vec<Spectrum> = reader.by_ref().map(|r| r.expect("parse")).collect();
    assert_eq!(spectra.len(), 1);
    assert_eq!(spectra[0].scan, Some(10));
    let ms3s = reader.take_product_scans();
    assert_eq!(ms3s.len(), 1);
    let p = &ms3s[0];
    assert_eq!(p.scan, Some(11));
    assert_eq!(p.ms_level, 3);
    assert_eq!(
        p.parent_id.as_deref(),
        Some("controllerType=0 controllerNumber=1 scan=10")
    );
    assert_eq!(p.rt_seconds, Some(120.0));
    assert_eq!(
        p.precursor_mz, 500.5,
        "the first selected ion gives the precursor m/z, as the first \
             precursor gives the parent reference"
    );
    assert_eq!(p.peaks.len(), 2);
    assert_eq!(p.peaks[1], (127.1248, 2000.0));
}

#[test]
fn noise_capture_reads_the_sampled_noise_arrays() {
    let noise_arrays = format!(
        "{}{}",
        bda_plain("MS:1002743", &encode_f64_b64(&[120.0, 130.0, 200.0, 900.0])),
        bda_plain("MS:1002744", &encode_f64_b64(&[100.0, 200.0, 300.0, 50.0])),
    );
    let ms2 = format!(
        r#"<spectrum index="0" id="scan=10" defaultArrayLength="1">
              <cvParam accession="MS:1000511" name="ms level" value="2"/>
              <scanList count="1"><scan/></scanList>
              <precursorList count="1"><precursor>
                <selectedIonList count="1"><selectedIon>
                  <cvParam accession="MS:1000744" name="selected ion m/z" value="500.5"/>
                  <cvParam accession="MS:1000041" name="charge state" value="2"/>
                </selectedIon></selectedIonList>
              </precursor></precursorList>
              <binaryDataArrayList count="4">{mz}{int}{noise_arrays}</binaryDataArrayList>
            </spectrum>"#,
        mz = bda_plain("MS:1000514", &encode_f64_b64(&[126.1277])),
        int = bda_plain("MS:1000515", &encode_f64_b64(&[1000.0])),
    );
    let ms3 = format!(
        r#"<spectrum index="1" id="scan=11" defaultArrayLength="1">
              <cvParam accession="MS:1000511" name="ms level" value="3"/>
              <scanList count="1"><scan/></scanList>
              <precursorList count="1"><precursor spectrumRef="scan=10">
                <selectedIonList count="1"><selectedIon>
                  <cvParam accession="MS:1000744" name="selected ion m/z" value="500.5"/>
                </selectedIon></selectedIonList>
              </precursor></precursorList>
              <binaryDataArrayList count="4">{mz}{int}{noise_arrays}</binaryDataArrayList>
            </spectrum>"#,
        mz = bda_plain("MS:1000514", &encode_f64_b64(&[126.1277])),
        int = bda_plain("MS:1000515", &encode_f64_b64(&[1000.0])),
    );
    let xml = wrap_spectra(&format!("{ms2}{ms3}"));
    let mut reader = MzMLReader::new(Cursor::new(xml.clone()))
        .with_product_scan_capture(true)
        .with_scan_context_capture(true)
        .with_noise_capture(true);
    let spectra: Vec<Spectrum> = reader.by_ref().map(|r| r.expect("parse")).collect();
    assert_eq!(spectra.len(), 1);
    assert_eq!(spectra[0].peaks, vec![(126.1277, 1000.0)]);
    let contexts = reader.take_scan_contexts();
    let ms3s = reader.take_product_scans();
    for noise in [&contexts[0].noise, &ms3s[0].noise] {
        let noise = noise.as_ref().expect("noise captured");
        assert_eq!(noise.len(), 3);
        assert!((noise.at(125.0) - 150.0).abs() < 1e-4);
    }
    let mut plain = MzMLReader::new(Cursor::new(xml))
        .with_product_scan_capture(true)
        .with_scan_context_capture(true);
    let _: Vec<Spectrum> = plain.by_ref().map(|r| r.expect("parse")).collect();
    assert!(plain.take_scan_contexts()[0].noise.is_none());
    assert!(plain.take_product_scans()[0].noise.is_none());
}

#[test]
fn chunked_ms1_linkage_matches_batch_across_boundaries() {
    let ms1a = format!(
        r#"<spectrum index="0" id="scan=1" defaultArrayLength="2">
              <cvParam accession="MS:1000511" name="ms level" value="1"/>
              <scanList count="1"><scan/></scanList>
              <binaryDataArrayList count="2">{mz}{int}</binaryDataArrayList>
            </spectrum>"#,
        mz = bda_plain("MS:1000514", &encode_f64_b64(&[111.0, 222.0])),
        int = bda_plain("MS:1000515", &encode_f64_b64(&[10.0, 20.0])),
    );
    let ms1b = format!(
        r#"<spectrum index="4" id="scan=5" defaultArrayLength="1">
              <cvParam accession="MS:1000511" name="ms level" value="1"/>
              <scanList count="1"><scan/></scanList>
              <binaryDataArrayList count="2">{mz}{int}</binaryDataArrayList>
            </spectrum>"#,
        mz = bda_plain("MS:1000514", &encode_f64_b64(&[333.0])),
        int = bda_plain("MS:1000515", &encode_f64_b64(&[30.0])),
    );
    let mk_ms2 = |id: &str, mz: f64| {
        ms2_spectrum_xml(
            id,
            &bda_plain("MS:1000514", &encode_f64_b64(&[mz])),
            &bda_plain("MS:1000515", &encode_f64_b64(&[1.0])),
            mz + 50.0,
            Some(2),
        )
    };
    let xml = format!(
        r#"<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml"><run><spectrumList count="6" defaultDataProcessingRef="dp">
{a}{m1}{m2}{m3}{b}{m4}
</spectrumList></run></mzML>"#,
        a = ms1a,
        m1 = mk_ms2("scan=2", 300.0),
        m2 = mk_ms2("scan=3", 400.0),
        m3 = mk_ms2("scan=4", 410.0),
        b = ms1b,
        m4 = mk_ms2("scan=6", 600.0),
    );

    let (b_spectra, b_link) = MzMLReader::new(Cursor::new(xml.clone()))
        .with_ms1_capture(true)
        .read_with_ms1()
        .expect("batch read");
    let batch_resolved: Vec<Option<&Vec<(f64, f32)>>> = b_link
        .ms2_to_ms1
        .iter()
        .map(|o| o.map(|i| &b_link.ms1_peaks[i]))
        .collect();

    let mut chunked_spectra: Vec<Spectrum> = Vec::new();
    let mut chunked_resolved: Vec<Option<Vec<(f64, f32)>>> = Vec::new();
    let (errc, _errs) = MzMLReader::new(Cursor::new(xml))
        .with_ms1_capture(true)
        .read_with_ms1_chunked(2, usize::MAX, |chunk, link| {
            for (i, s) in chunk.into_iter().enumerate() {
                chunked_resolved.push(link.ms2_to_ms1[i].map(|j| link.ms1_peaks[j].clone()));
                chunked_spectra.push(s);
            }
        });

    assert_eq!(errc, 0, "clean input → no parse errors");
    assert_eq!(
        chunked_spectra.len(),
        b_spectra.len(),
        "same MS2 count as batch"
    );
    assert_eq!(chunked_resolved.len(), batch_resolved.len());
    for (i, (cr, br)) in chunked_resolved
        .iter()
        .zip(batch_resolved.iter())
        .enumerate()
    {
        match (cr, br) {
            (Some(c), Some(b)) => {
                assert_eq!(c, *b, "MS2 #{i} linked to a different MS1 than batch")
            }
            (None, None) => {}
            _ => panic!("MS2 #{i} linkage presence differs from batch: {cr:?} vs {br:?}"),
        }
    }
    assert_eq!(chunked_resolved[2].as_ref().map(|v| v.len()), Some(2));
    assert_eq!(chunked_resolved[3].as_ref().map(|v| v.len()), Some(1));
}

#[test]
fn chunked_read_keeps_run_ms1_and_scan_contexts_for_quantification() {
    let ms1 = |id: &str, rt_min: f64, mz: f64| {
        format!(
            r#"<spectrum index="0" id="{id}" defaultArrayLength="1">
              <cvParam accession="MS:1000511" name="ms level" value="1"/>
              <scanList count="1"><scan>
                <cvParam accession="MS:1000016" name="scan start time" value="{rt_min}"
                         unitAccession="UO:0000031" unitName="minute"/>
              </scan></scanList>
              <binaryDataArrayList count="2">{mzb}{intb}</binaryDataArrayList>
            </spectrum>"#,
            mzb = bda_plain("MS:1000514", &encode_f64_b64(&[mz])),
            intb = bda_plain("MS:1000515", &encode_f64_b64(&[7.0])),
        )
    };
    let ms2 = |id: &str, mz: f64| {
        ms2_spectrum_xml(
            id,
            &bda_plain("MS:1000514", &encode_f64_b64(&[mz])),
            &bda_plain("MS:1000515", &encode_f64_b64(&[1.0])),
            mz + 50.0,
            Some(2),
        )
    };
    let xml = wrap_spectra(&format!(
        "{}{}{}{}{}{}",
        ms1("scan=1", 1.0, 400.0),
        ms2("scan=2", 300.0),
        ms2("scan=3", 310.0),
        ms1("scan=4", 1.5, 401.0),
        ms2("scan=5", 320.0),
        ms1("scan=6", 2.0, 402.0),
    ));
    let mut reader = MzMLReader::new(Cursor::new(xml.clone()))
        .with_run_ms1_capture(true)
        .with_scan_context_capture(true);
    let mut linked: Vec<Option<Vec<(f64, f32)>>> = Vec::new();
    let mut n_ms2 = 0usize;
    let (errc, _) = reader.read_with_ms1_chunked(2, usize::MAX, |chunk, link| {
        n_ms2 += chunk.len();
        for i in 0..chunk.len() {
            linked.push(link.ms2_to_ms1[i].map(|j| link.ms1_peaks[j].clone()));
        }
    });
    assert_eq!(errc, 0);
    assert_eq!(n_ms2, 3, "only MS2 is emitted");
    assert_eq!(
        linked,
        vec![
            Some(vec![(400.0, 7.0)]),
            Some(vec![(400.0, 7.0)]),
            Some(vec![(401.0, 7.0)])
        ]
    );
    let run = reader.take_run_ms1();
    assert_eq!(
        run.iter()
            .map(|s| (s.id.as_str(), s.rt))
            .collect::<Vec<_>>(),
        vec![("scan=1", 60.0), ("scan=4", 90.0), ("scan=6", 120.0)],
        "every MS1 of the run, the trailing one included"
    );
    assert_eq!(run[2].peaks, vec![(402.0, 7.0)]);
    let contexts = reader.take_scan_contexts();
    assert_eq!(
        contexts.iter().map(|c| c.id.as_str()).collect::<Vec<_>>(),
        vec!["scan=2", "scan=3", "scan=5"]
    );
    let mut plain = MzMLReader::new(Cursor::new(xml))
        .with_run_ms1_capture(true)
        .with_scan_context_capture(true);
    assert_eq!(plain.by_ref().count(), 3);
    assert_eq!(plain.take_run_ms1(), run);
    assert_eq!(plain.take_scan_contexts(), contexts);
}

#[test]
fn iterator_resyncs_past_a_malformed_spectrum() {
    let good1 = ms2_spectrum_xml(
        "scan=1",
        &bda_plain("MS:1000514", &encode_f64_b64(&[300.0])),
        &bda_plain("MS:1000515", &encode_f64_b64(&[1.0])),
        350.0,
        Some(2),
    );
    let bad = ms2_spectrum_xml(
        "scan=2",
        &bda_plain("MS:1000514", "@@@not-valid-base64@@@"),
        &bda_plain("MS:1000515", &encode_f64_b64(&[1.0])),
        450.0,
        Some(2),
    );
    let good2 = ms2_spectrum_xml(
        "scan=3",
        &bda_plain("MS:1000514", &encode_f64_b64(&[500.0])),
        &bda_plain("MS:1000515", &encode_f64_b64(&[1.0])),
        550.0,
        Some(2),
    );
    let xml = format!(
        r#"<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml"><run><spectrumList count="3" defaultDataProcessingRef="dp">
{good1}{bad}{good2}
</spectrumList></run></mzML>"#
    );

    let got = collect_ok(&xml);
    assert_eq!(
        got.len(),
        2,
        "both good spectra must survive the resync (no truncation)"
    );
    assert_eq!(got[0].scan, Some(1));
    assert_eq!(
        got[1].scan,
        Some(3),
        "the post-error spectrum must still be parsed"
    );
}

#[test]
fn chunked_resyncs_past_a_malformed_spectrum() {
    let good1 = ms2_spectrum_xml(
        "scan=1",
        &bda_plain("MS:1000514", &encode_f64_b64(&[300.0])),
        &bda_plain("MS:1000515", &encode_f64_b64(&[1.0])),
        350.0,
        Some(2),
    );
    let bad = ms2_spectrum_xml(
        "scan=2",
        &bda_plain("MS:1000514", "@@@not-valid-base64@@@"),
        &bda_plain("MS:1000515", &encode_f64_b64(&[1.0])),
        450.0,
        Some(2),
    );
    let good2 = ms2_spectrum_xml(
        "scan=3",
        &bda_plain("MS:1000514", &encode_f64_b64(&[500.0])),
        &bda_plain("MS:1000515", &encode_f64_b64(&[1.0])),
        550.0,
        Some(2),
    );
    let xml = format!(
        r#"<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml"><run><spectrumList count="3" defaultDataProcessingRef="dp">
{good1}{bad}{good2}
</spectrumList></run></mzML>"#
    );

    let mut got: Vec<Spectrum> = Vec::new();
    let (errc, errs) = MzMLReader::new(Cursor::new(xml))
        .with_ms1_capture(true)
        .read_with_ms1_chunked(5000, usize::MAX, |chunk, _link| got.extend(chunk));

    assert_eq!(errc, 1, "exactly one malformed spectrum should be counted");
    assert_eq!(errs.len(), 1, "the error message should be captured");
    assert_eq!(
        got.len(),
        2,
        "both good spectra must survive the resync (no truncation)"
    );
    assert_eq!(got[0].scan, Some(1));
    assert_eq!(
        got[1].scan,
        Some(3),
        "the post-error spectrum must still be parsed"
    );
}

#[test]
fn ms1_capture_ms2_before_any_ms1_links_to_none() {
    let ms2a = ms2_spectrum_xml(
        "scan=1",
        &bda_plain("MS:1000514", &encode_f64_b64(&[300.0])),
        &bda_plain("MS:1000515", &encode_f64_b64(&[1.0])),
        500.0,
        Some(2),
    );
    let ms1 = format!(
        r#"<spectrum index="1" id="scan=2" defaultArrayLength="1">
              <cvParam accession="MS:1000511" name="ms level" value="1"/>
              <scanList count="1"><scan/></scanList>
              <binaryDataArrayList count="2">
                {mz}
                {int}
              </binaryDataArrayList>
            </spectrum>"#,
        mz = bda_plain("MS:1000514", &encode_f64_b64(&[999.0])),
        int = bda_plain("MS:1000515", &encode_f64_b64(&[5.0])),
    );
    let ms2b = ms2_spectrum_xml(
        "scan=3",
        &bda_plain("MS:1000514", &encode_f64_b64(&[400.0])),
        &bda_plain("MS:1000515", &encode_f64_b64(&[2.0])),
        600.0,
        Some(2),
    );

    let xml = format!(
        r#"<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml">
  <run>
    <spectrumList count="3" defaultDataProcessingRef="dp">
      {ms2a}
      {ms1}
      {ms2b}
    </spectrumList>
  </run>
</mzML>"#
    );

    let (spectra, link) = MzMLReader::new(Cursor::new(xml))
        .with_ms1_capture(true)
        .read_with_ms1()
        .expect("read_with_ms1 failed");

    assert_eq!(spectra.len(), 2, "expected two MS2 spectra");
    assert_eq!(link.ms1_peaks.len(), 1, "one captured MS1");
    assert_eq!(link.ms2_to_ms1, vec![None, Some(0)]);
}

#[test]
fn ms1_capture_off_yields_empty_link() {
    let mz_b64 = encode_f64_b64(&[100.0, 200.0]);
    let int_b64 = encode_f64_b64(&[1000.0, 500.0]);
    let spec = ms2_spectrum_xml(
        "scan=1",
        &bda_plain("MS:1000514", &mz_b64),
        &bda_plain("MS:1000515", &int_b64),
        500.5,
        Some(2),
    );
    let (spectra, link) = MzMLReader::new(Cursor::new(wrap_spectra(&spec)))
        .read_with_ms1()
        .expect("read_with_ms1 failed");
    assert_eq!(spectra.len(), 1);
    assert!(
        link.ms1_peaks.is_empty(),
        "no MS1 captured when capture off"
    );
    assert_eq!(link.ms2_to_ms1, vec![None]);
}

#[test]
fn strict_iterator_surfaces_malformed_spectrum_as_err() {
    let xml = wrap_spectra(&malformed_spectrum());
    let mut reader = MzMLReader::new(Cursor::new(xml)).strict();
    let first = reader.next().expect("should yield a result");
    assert!(first.is_err(), "strict mode must surface the parse error");
}

#[test]
fn tolerant_iterator_skips_and_counts_malformed_spectrum() {
    let good = ms2_spectrum_xml(
        "scan=2",
        &bda_plain("MS:1000514", &encode_f64_b64(&[150.0, 300.0])),
        &bda_plain("MS:1000515", &encode_f64_b64(&[2000.0, 1000.0])),
        500.5,
        Some(2),
    );
    let xml = wrap_spectra(&format!("{}{}", malformed_spectrum(), good));
    let mut reader = MzMLReader::new(Cursor::new(xml));
    let collected: Vec<_> = reader.by_ref().filter_map(|r| r.ok()).collect();
    assert_eq!(collected.len(), 1, "tolerant mode keeps the good spectrum");
    assert!(
        reader.skipped_count() >= 1,
        "tolerant mode records the dropped scan"
    );
}
