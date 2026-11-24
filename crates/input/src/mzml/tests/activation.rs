//! Tests of activation parsing: single and combined terms, user params, EThcD/ETciD
//! detection from terms and from the filter string.

use std::io::Cursor;

use model::{ActivationMethod, Spectrum};

use super::*;
use crate::mzml::{filter_has_supplemental_collision, MzMLReader};

fn spectrum_xml_with_activation(activation_cv: Option<&str>) -> String {
    let mz_b64 = encode_f64_b64(&[100.0]);
    let int_b64 = encode_f64_b64(&[1000.0]);
    let act_block = match activation_cv {
        Some(cv) => format!(
            r#"<activation>
                     <cvParam accession="{cv}" name="" value=""/>
                   </activation>"#
        ),
        None => String::new(),
    };
    format!(
        r#"<spectrum index="0" id="scan=1" defaultArrayLength="1">
              <cvParam accession="MS:1000511" name="ms level" value="2"/>
              <scanList count="1"><scan/></scanList>
              <precursorList count="1">
                <precursor>
                  <selectedIonList count="1">
                    <selectedIon>
                      <cvParam accession="MS:1000744" name="selected ion m/z"
                               value="500.5"/>
                      <cvParam accession="MS:1000041" name="charge state" value="2"/>
                    </selectedIon>
                  </selectedIonList>
                  {act_block}
                </precursor>
              </precursorList>
              <binaryDataArrayList count="2">
                {mz}
                {int}
              </binaryDataArrayList>
            </spectrum>"#,
        mz = bda_plain("MS:1000514", &mz_b64),
        int = bda_plain("MS:1000515", &int_b64),
    )
}

#[test]
fn parses_cid_activation() {
    let spectra = collect_ok(&wrap_spectra(&spectrum_xml_with_activation(Some(
        "MS:1000133",
    ))));
    assert_eq!(spectra.len(), 1);
    assert_eq!(spectra[0].activation_method, Some(ActivationMethod::CID));
}

#[test]
fn parses_hcd_activation() {
    let spectra = collect_ok(&wrap_spectra(&spectrum_xml_with_activation(Some(
        "MS:1000422",
    ))));
    assert_eq!(spectra.len(), 1);
    assert_eq!(spectra[0].activation_method, Some(ActivationMethod::HCD));
}

#[test]
fn parses_etd_activation() {
    let spectra = collect_ok(&wrap_spectra(&spectrum_xml_with_activation(Some(
        "MS:1000598",
    ))));
    assert_eq!(spectra.len(), 1);
    assert_eq!(spectra[0].activation_method, Some(ActivationMethod::ETD));
}

#[test]
fn parses_ecd_as_etd() {
    let spectra = collect_ok(&wrap_spectra(&spectrum_xml_with_activation(Some(
        "MS:1000250",
    ))));
    assert_eq!(spectra.len(), 1);
    assert_eq!(spectra[0].activation_method, Some(ActivationMethod::ETD));
}

#[test]
fn missing_activation_block_yields_none() {
    let spectra = collect_ok(&wrap_spectra(&spectrum_xml_with_activation(None)));
    assert_eq!(spectra.len(), 1);
    assert_eq!(spectra[0].activation_method, None);
}

fn spectrum_xml_with_two_activations(cv1: &str, cv2: &str) -> String {
    let mz_b64 = encode_f64_b64(&[100.0]);
    let int_b64 = encode_f64_b64(&[1000.0]);
    format!(
        r#"<spectrum index="0" id="scan=1" defaultArrayLength="1">
              <cvParam accession="MS:1000511" name="ms level" value="2"/>
              <scanList count="1"><scan/></scanList>
              <precursorList count="1">
                <precursor>
                  <selectedIonList count="1">
                    <selectedIon>
                      <cvParam accession="MS:1000744" name="selected ion m/z" value="500.5"/>
                      <cvParam accession="MS:1000041" name="charge state" value="2"/>
                    </selectedIon>
                  </selectedIonList>
                  <activation>
                    <cvParam accession="{cv1}" name="" value=""/>
                    <cvParam accession="{cv2}" name="" value=""/>
                  </activation>
                </precursor>
              </precursorList>
              <binaryDataArrayList count="2">
                {mz}
                {int}
              </binaryDataArrayList>
            </spectrum>"#,
        mz = bda_plain("MS:1000514", &mz_b64),
        int = bda_plain("MS:1000515", &int_b64),
    )
}

#[test]
fn ethcd_etd_plus_hcd_is_labelled_ethcd() {
    let spectra = collect_ok(&wrap_spectra(&spectrum_xml_with_two_activations(
        "MS:1000598",
        "MS:1000422",
    )));
    assert_eq!(spectra.len(), 1);
    assert_eq!(spectra[0].activation_method, Some(ActivationMethod::EThcD));
}

#[test]
fn ethcd_detection_is_order_independent() {
    let spectra = collect_ok(&wrap_spectra(&spectrum_xml_with_two_activations(
        "MS:1000422",
        "MS:1000598",
    )));
    assert_eq!(spectra.len(), 1);
    assert_eq!(spectra[0].activation_method, Some(ActivationMethod::EThcD));
}

#[test]
fn supplemental_activation_terms_mark_ethcd_and_etcid() {
    for supp in ["MS:1002678", "MS:1002679"] {
        for (a, b) in [("MS:1000598", supp), (supp, "MS:1000598")] {
            let spectra = collect_ok(&wrap_spectra(&spectrum_xml_with_two_activations(a, b)));
            assert_eq!(
                spectra[0].activation_method,
                Some(ActivationMethod::EThcD),
                "{a} + {b}"
            );
        }
    }
}

#[test]
fn supplemental_collision_read_from_the_filter_chain() {
    for f in [
        "FTMS + p NSI d Full ms2 784.37@etd50.00@hcd25.00 [120.00-2000.00]",
        "ITMS + c NSI d Full ms2 652.10@etd100.00@cid35.00 [50.00-1315.00]",
        "FTMS + c NSI Full ms3 900.00@cid35.00 450.20@etd50.00@hcd20.00 [100.00-1000.00]",
        "FTMS + p NSI d Full ms2 391.1936@etd55.49 391.1936@hcd20.00 [100.0000-1184.0000]",
        "FTMS + c NSI d sa Full ms2 361.9050@etd143.26@hcd28.00 [100.0000-734.0000]",
    ] {
        assert!(filter_has_supplemental_collision(f), "{f}");
    }
    for f in [
        "FTMS + p NSI d Full ms2 784.37@etd50.00 [120.00-2000.00]",
        "FTMS + p NSI d Full ms2 784.37@hcd28.00 [120.00-2000.00]",
        "FTMS + c NSI Full ms3 900.00@cid35.00 450.20@etd50.00 [100.00-1000.00]",
        "FTMS + c NSI Full ms3 900.00@hcd35.00 450.20@etd50.00 [100.00-1000.00]",
        "FTMS + p NSI Full ms [350.00-1800.00]",
        "",
    ] {
        assert!(!filter_has_supplemental_collision(f), "{f}");
    }
}

#[test]
fn etd_term_with_supplemental_filter_string_is_labelled_ethcd() {
    let xml = spectrum_xml_with_activation(Some("MS:1000598")).replace(
        "<scanList count=\"1\"><scan/></scanList>",
        "<scanList count=\"1\"><scan><cvParam accession=\"MS:1000512\" name=\"filter string\" \
             value=\"FTMS + p NSI d Full ms2 500.50@etd55.49 500.50@hcd20.00 [100.00-1000.00]\"/>\
             </scan></scanList>",
    );
    assert!(
        xml.contains("MS:1000512"),
        "fixture must carry the filter string"
    );
    let spectra = collect_ok(&wrap_spectra(&xml));
    assert_eq!(spectra[0].activation_method, Some(ActivationMethod::EThcD));
    let xml = spectrum_xml_with_activation(Some("MS:1000598")).replace(
        "<scanList count=\"1\"><scan/></scanList>",
        "<scanList count=\"1\"><scan><cvParam accession=\"MS:1000512\" name=\"filter string\" \
             value=\"FTMS + p NSI d Full ms2 500.50@etd55.49 [100.00-1000.00]\"/></scan></scanList>",
    );
    let spectra = collect_ok(&wrap_spectra(&xml));
    assert_eq!(spectra[0].activation_method, Some(ActivationMethod::ETD));
}

#[test]
fn combined_and_child_activation_terms() {
    for (cv, want) in [
        ("MS:1002631", ActivationMethod::EThcD),
        ("MS:1003182", ActivationMethod::EThcD),
        ("MS:1002472", ActivationMethod::CID),
        ("MS:1002481", ActivationMethod::HCD),
        ("MS:1003246", ActivationMethod::UVPD),
    ] {
        let spectra = collect_ok(&wrap_spectra(&spectrum_xml_with_activation(Some(cv))));
        assert_eq!(spectra[0].activation_method, Some(want), "{cv}");
    }
    let spectra = collect_ok(&wrap_spectra(&spectrum_xml_with_activation(Some(
        "MS:1002678",
    ))));
    assert_eq!(spectra[0].activation_method, None);
}

#[test]
fn etcid_etd_plus_cid_is_labelled_ethcd() {
    let spectra = collect_ok(&wrap_spectra(&spectrum_xml_with_two_activations(
        "MS:1000598",
        "MS:1000133",
    )));
    assert_eq!(spectra.len(), 1);
    assert_eq!(spectra[0].activation_method, Some(ActivationMethod::EThcD));
}

#[test]
fn multiple_activations_first_wins() {
    let mz_b64 = encode_f64_b64(&[100.0]);
    let int_b64 = encode_f64_b64(&[1000.0]);
    let xml = format!(
        r#"<spectrum index="0" id="scan=1" defaultArrayLength="1">
              <cvParam accession="MS:1000511" name="ms level" value="3"/>
              <scanList count="1"><scan/></scanList>
              <precursorList count="2">
                <precursor>
                  <selectedIonList count="1">
                    <selectedIon>
                      <cvParam accession="MS:1000744" name="selected ion m/z" value="500.5"/>
                    </selectedIon>
                  </selectedIonList>
                  <activation>
                    <cvParam accession="MS:1000133" name="CID" value=""/>
                  </activation>
                </precursor>
                <precursor>
                  <selectedIonList count="1">
                    <selectedIon>
                      <cvParam accession="MS:1000744" name="selected ion m/z" value="350.0"/>
                    </selectedIon>
                  </selectedIonList>
                  <activation>
                    <cvParam accession="MS:1000422" name="HCD" value=""/>
                  </activation>
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

    let wrapped = format!(
        r#"<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml">
  <run>
    <spectrumList count="1" defaultDataProcessingRef="dp">
      {xml}
    </spectrumList>
  </run>
</mzML>"#
    );
    let spectra: Vec<Spectrum> = MzMLReader::new(Cursor::new(wrapped))
        .with_ms_level_range(2, 3)
        .map(|r| r.expect("parse error"))
        .collect();
    assert_eq!(spectra.len(), 1);
    assert_eq!(spectra[0].activation_method, Some(ActivationMethod::CID));
}

fn activation_of(activation_xml: &str) -> Option<ActivationMethod> {
    let mz_b64 = encode_f64_b64(&[100.0]);
    let int_b64 = encode_f64_b64(&[1000.0]);
    let xml = format!(
        r#"<spectrum index="0" id="scan=1" defaultArrayLength="1">
              <cvParam accession="MS:1000511" name="ms level" value="2"/>
              <scanList count="1"><scan/></scanList>
              <precursorList count="1">
                <precursor>
                  <selectedIonList count="1"><selectedIon>
                    <cvParam accession="MS:1000744" name="selected ion m/z" value="500.5"/>
                  </selectedIon></selectedIonList>
                  <activation>{activation_xml}</activation>
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
    let wrapped = format!(
        r#"<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml">
  <run><spectrumList count="1" defaultDataProcessingRef="dp">
      {xml}
  </spectrumList></run>
</mzML>"#
    );
    let spectra: Vec<Spectrum> = MzMLReader::new(Cursor::new(wrapped))
        .map(|r| r.expect("parse error"))
        .collect();
    assert_eq!(spectra.len(), 1);
    spectra[0].activation_method
}

#[test]
fn activation_user_param_sets_activation() {
    assert_eq!(
        activation_of(r#"<userParam name="CID" value=""/>"#),
        Some(ActivationMethod::CID)
    );
    assert_eq!(
        activation_of(r#"<userParam name="beam-type collision-induced dissociation"/>"#),
        Some(ActivationMethod::HCD)
    );
    assert_eq!(
        activation_of(r#"<userParam name="ETD"/>"#),
        Some(ActivationMethod::ETD)
    );
    assert_eq!(
        activation_of(r#"<userParam name="collision energy" value="35"/>"#),
        None
    );
    assert_eq!(
        activation_of(
            r#"<cvParam accession="MS:1000422" name="HCD" value=""/><userParam name="CID"/>"#
        ),
        Some(ActivationMethod::HCD)
    );
}

#[test]
fn separate_activation_blocks_do_not_falsely_trigger_ethcd() {
    let mz_b64 = encode_f64_b64(&[100.0]);
    let int_b64 = encode_f64_b64(&[1000.0]);
    let xml = format!(
        r#"<spectrum index="0" id="scan=1" defaultArrayLength="1">
              <cvParam accession="MS:1000511" name="ms level" value="3"/>
              <scanList count="1"><scan/></scanList>
              <precursorList count="2">
                <precursor>
                  <selectedIonList count="1"><selectedIon>
                    <cvParam accession="MS:1000744" name="selected ion m/z" value="500.5"/>
                  </selectedIon></selectedIonList>
                  <activation><cvParam accession="MS:1000133" name="CID" value=""/></activation>
                </precursor>
                <precursor>
                  <selectedIonList count="1"><selectedIon>
                    <cvParam accession="MS:1000744" name="selected ion m/z" value="350.0"/>
                  </selectedIon></selectedIonList>
                  <activation><cvParam accession="MS:1000598" name="ETD" value=""/></activation>
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
    let wrapped = format!(
        r#"<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml">
  <run><spectrumList count="1" defaultDataProcessingRef="dp">
      {xml}
  </spectrumList></run>
</mzML>"#
    );
    let spectra: Vec<Spectrum> = MzMLReader::new(Cursor::new(wrapped))
        .with_ms_level_range(2, 3)
        .map(|r| r.expect("parse error"))
        .collect();
    assert_eq!(spectra.len(), 1);
    assert_ne!(
        spectra[0].activation_method,
        Some(ActivationMethod::EThcD),
        "separate CID + ETD blocks must not be mis-detected as EThcD"
    );
    assert_eq!(spectra[0].activation_method, Some(ActivationMethod::ETD));
}

#[test]
fn etd_plus_supplemental_collisional_is_labelled_ethcd() {
    let mz_b64 = encode_f64_b64(&[100.0]);
    let int_b64 = encode_f64_b64(&[1000.0]);
    let xml = format!(
        r#"<spectrum index="0" id="scan=1" defaultArrayLength="1">
              <cvParam accession="MS:1000511" name="ms level" value="2"/>
              <scanList count="1"><scan/></scanList>
              <precursorList count="1">
                <precursor>
                  <selectedIonList count="1">
                    <selectedIon>
                      <cvParam accession="MS:1000744" name="selected ion m/z" value="500.5"/>
                    </selectedIon>
                  </selectedIonList>
                  <activation>
                    <cvParam accession="MS:1000133" name="CID" value=""/>
                    <cvParam accession="MS:1000598" name="ETD" value=""/>
                  </activation>
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
    assert_eq!(spectra[0].activation_method, Some(ActivationMethod::EThcD));
}
