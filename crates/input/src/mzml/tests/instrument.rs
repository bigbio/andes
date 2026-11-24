//! Tests of instrument detection from the mzML header.

use std::io::Cursor;

use model::InstrumentType;

use super::*;
use crate::mzml::detect_instrument_type;

fn wrap_with_instrument_configs(
    instrument_configs: &str,
    default_ic_ref: &str,
    spectra_xml: &str,
) -> String {
    format!(
        r#"<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml">
  <instrumentConfigurationList count="1">
    {instrument_configs}
  </instrumentConfigurationList>
  <run id="r" defaultInstrumentConfigurationRef="{default_ic_ref}">
    <spectrumList count="1" defaultDataProcessingRef="dp">
      {spectra_xml}
    </spectrumList>
  </run>
</mzML>"#
    )
}

fn ic_block(id: &str, analyzer_cv: &str) -> String {
    format!(
        r#"<instrumentConfiguration id="{id}">
              <componentList count="3">
                <source order="1">
                  <cvParam accession="MS:1000398" name="nanoelectrospray" value=""/>
                </source>
                <analyzer order="2">
                  <cvParam accession="{analyzer_cv}" name="" value=""/>
                </analyzer>
                <detector order="3">
                  <cvParam accession="MS:1000624" name="inductive detector" value=""/>
                </detector>
              </componentList>
            </instrumentConfiguration>"#
    )
}

fn ms2_spectrum_with_ic_ref(ic_ref: &str) -> String {
    let mz_b64 = encode_f64_b64(&[100.0]);
    let int_b64 = encode_f64_b64(&[1000.0]);
    format!(
        r#"<spectrum index="0" id="scan=1" defaultArrayLength="1">
              <cvParam accession="MS:1000511" name="ms level" value="2"/>
              <scanList count="1">
                <scan instrumentConfigurationRef="{ic_ref}"/>
              </scanList>
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
                {mz}
                {int}
              </binaryDataArrayList>
            </spectrum>"#,
        mz = bda_plain("MS:1000514", &mz_b64),
        int = bda_plain("MS:1000515", &int_b64),
    )
}

#[test]
fn detect_instrument_orbitrap_analyzer_to_qexactive() {
    let xml = wrap_with_instrument_configs(
        &ic_block("IC1", "MS:1000484"),
        "IC1",
        &ms2_spectrum_with_ic_ref("IC1"),
    );
    let result = detect_instrument_type(Cursor::new(xml));
    assert_eq!(result, Some(InstrumentType::QExactive));
}

#[test]
fn detect_instrument_ion_trap_analyzer_to_lowres() {
    let xml = wrap_with_instrument_configs(
        &ic_block("IC1", "MS:1000291"),
        "IC1",
        &ms2_spectrum_with_ic_ref("IC1"),
    );
    let result = detect_instrument_type(Cursor::new(xml));
    assert_eq!(result, Some(InstrumentType::LowRes));
}

#[test]
fn detect_instrument_quad_ion_trap_to_lowres() {
    let xml = wrap_with_instrument_configs(
        &ic_block("IC1", "MS:1000082"),
        "IC1",
        &ms2_spectrum_with_ic_ref("IC1"),
    );
    let result = detect_instrument_type(Cursor::new(xml));
    assert_eq!(result, Some(InstrumentType::LowRes));
}

#[test]
fn detect_instrument_fticr_to_highres() {
    let xml = wrap_with_instrument_configs(
        &ic_block("IC1", "MS:1000079"),
        "IC1",
        &ms2_spectrum_with_ic_ref("IC1"),
    );
    let result = detect_instrument_type(Cursor::new(xml));
    assert_eq!(result, Some(InstrumentType::HighRes));
}

#[test]
fn detect_instrument_tof_analyzer() {
    let xml = wrap_with_instrument_configs(
        &ic_block("IC1", "MS:1000084"),
        "IC1",
        &ms2_spectrum_with_ic_ref("IC1"),
    );
    let result = detect_instrument_type(Cursor::new(xml));
    assert_eq!(result, Some(InstrumentType::TOF));
}

#[test]
fn detect_instrument_ms2_referenced_ic_wins_pxd001819_pattern() {
    let ics = format!(
        "{}\n{}",
        ic_block("IC1", "MS:1000484"),
        ic_block("IC2", "MS:1000264"),
    );
    let xml = wrap_with_instrument_configs(&ics, "IC1", &ms2_spectrum_with_ic_ref("IC2"));
    let result = detect_instrument_type(Cursor::new(xml));
    assert_eq!(result, Some(InstrumentType::LowRes));
}

#[test]
fn detect_instrument_falls_back_to_default_ic_when_scan_lacks_ref() {
    let mz_b64 = encode_f64_b64(&[100.0]);
    let int_b64 = encode_f64_b64(&[1000.0]);
    let spec = format!(
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
    let xml = wrap_with_instrument_configs(&ic_block("IC1", "MS:1000484"), "IC1", &spec);
    let result = detect_instrument_type(Cursor::new(xml));
    assert_eq!(result, Some(InstrumentType::QExactive));
}

#[test]
fn detect_instrument_returns_none_when_no_ic_info() {
    let mz_b64 = encode_f64_b64(&[100.0]);
    let int_b64 = encode_f64_b64(&[1000.0]);
    let spec = ms2_spectrum_xml(
        "scan=1",
        &bda_plain("MS:1000514", &mz_b64),
        &bda_plain("MS:1000515", &int_b64),
        500.5,
        None,
    );
    let xml = wrap_spectra(&spec);
    let result = detect_instrument_type(Cursor::new(xml));
    assert_eq!(result, None);
}

#[test]
fn detect_instrument_qexactive_model_cv_param() {
    let ic = r#"<instrumentConfiguration id="IC1">
            <cvParam accession="MS:1001911" name="Q Exactive" value=""/>
            <componentList count="3">
              <source order="1">
                <cvParam accession="MS:1000398" name="nanoelectrospray" value=""/>
              </source>
              <analyzer order="2"/>
              <detector order="3"/>
            </componentList>
          </instrumentConfiguration>"#;
    let xml = wrap_with_instrument_configs(ic, "IC1", &ms2_spectrum_with_ic_ref("IC1"));
    let result = detect_instrument_type(Cursor::new(xml));
    assert_eq!(result, Some(InstrumentType::QExactive));
}

#[test]
fn detect_instrument_astral_model_in_referenceable_param_group() {
    let xml = r#"<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml">
  <referenceableParamGroupList count="1">
    <referenceableParamGroup id="commonInstrumentParams">
      <cvParam accession="MS:1003378" name="Orbitrap Astral" value=""/>
      <cvParam accession="MS:1000529" name="instrument serial number" value="FSN123"/>
    </referenceableParamGroup>
  </referenceableParamGroupList>
  <instrumentConfigurationList count="1">
    <instrumentConfiguration id="IC1">
      <referenceableParamGroupRef ref="commonInstrumentParams"/>
      <componentList count="3">
        <source order="1">
          <cvParam accession="MS:1000398" name="nanoelectrospray" value=""/>
        </source>
        <analyzer order="2">
          <cvParam accession="MS:1000484" name="orbitrap" value=""/>
        </analyzer>
        <detector order="3"/>
      </componentList>
    </instrumentConfiguration>
  </instrumentConfigurationList>
  <run id="r" defaultInstrumentConfigurationRef="IC1">
    <spectrumList count="1" defaultDataProcessingRef="dp">
      {spectrum}
    </spectrumList>
  </run>
</mzML>"#
        .replace("{spectrum}", &ms2_spectrum_with_ic_ref("IC1"));
    let result = detect_instrument_type(Cursor::new(xml));
    assert_eq!(result, Some(InstrumentType::OrbitrapAstral));
}

#[test]
fn detect_instrument_fusion_ion_trap_ms2_is_low_res() {
    let doc = |ms2_ic: &str| {
        r#"<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml">
  <referenceableParamGroupList count="1">
    <referenceableParamGroup id="common">
      <cvParam accession="MS:1002416" name="Orbitrap Fusion" value=""/>
    </referenceableParamGroup>
  </referenceableParamGroupList>
  <instrumentConfigurationList count="2">
    <instrumentConfiguration id="IC1">
      <referenceableParamGroupRef ref="common"/>
      <componentList count="1">
        <analyzer order="2"><cvParam accession="MS:1000484" name="orbitrap" value=""/></analyzer>
      </componentList>
    </instrumentConfiguration>
    <instrumentConfiguration id="IC2">
      <referenceableParamGroupRef ref="common"/>
      <componentList count="1">
        <analyzer order="2"><cvParam accession="MS:1000083" name="radial ejection linear ion trap" value=""/></analyzer>
      </componentList>
    </instrumentConfiguration>
  </instrumentConfigurationList>
  <run id="r" defaultInstrumentConfigurationRef="IC1">
    <spectrumList count="1" defaultDataProcessingRef="dp">
      {spectrum}
    </spectrumList>
  </run>
</mzML>"#
            .replace("{spectrum}", &ms2_spectrum_with_ic_ref(ms2_ic))
    };
    assert_eq!(
        detect_instrument_type(Cursor::new(doc("IC2"))),
        Some(InstrumentType::LowRes)
    );
    assert_eq!(
        detect_instrument_type(Cursor::new(doc("IC1"))),
        Some(InstrumentType::QExactive)
    );
}
