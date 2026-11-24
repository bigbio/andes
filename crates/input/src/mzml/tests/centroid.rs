//! Tests of on-the-fly centroiding of profile spectra.

use super::*;
use crate::mzml::centroid::centroid_profile;

#[test]
fn centroid_profile_collapses_one_peak_to_intensity_weighted_centroid() {
    let prof = [
        (100.0, 1.0f32),
        (100.1, 4.0),
        (100.2, 9.0),
        (100.3, 4.0),
        (100.4, 1.0),
    ];
    let c = centroid_profile(&prof);
    assert_eq!(c.len(), 1, "one peak ⇒ one centroid");
    assert!(
        (c[0].0 - 100.2).abs() < 1e-9,
        "centroid m/z {} != 100.2",
        c[0].0
    );
    assert_eq!(c[0].1, 9.0, "centroid intensity is the apex");
}

#[test]
fn centroid_profile_separates_two_peaks() {
    let prof = [
        (200.0, 1.0f32),
        (200.1, 5.0),
        (200.2, 1.0),
        (200.3, 0.0),
        (200.4, 2.0),
        (200.5, 8.0),
        (200.6, 2.0),
    ];
    let c = centroid_profile(&prof);
    assert_eq!(c.len(), 2, "two peaks ⇒ two centroids");
    assert!((c[0].0 - 200.1).abs() < 1e-9);
    assert!((c[1].0 - 200.5).abs() < 1e-9);
}

fn profile_spectrum_xml(profile: bool) -> String {
    let mz_b64 = encode_f64_b64(&[150.0, 150.1, 150.2]);
    let int_b64 = encode_f64_b64(&[2.0, 9.0, 2.0]);
    let rep = if profile {
        r#"<cvParam accession="MS:1000128" name="profile spectrum" value=""/>"#
    } else {
        r#"<cvParam accession="MS:1000127" name="centroid spectrum" value=""/>"#
    };
    format!(
        r#"<spectrum index="0" id="scan=1" defaultArrayLength="3">
              <cvParam accession="MS:1000511" name="ms level" value="2"/>
              {rep}
              <scanList count="1"><scan/></scanList>
              <precursorList count="1"><precursor><selectedIonList count="1"><selectedIon>
                <cvParam accession="MS:1000744" name="selected ion m/z" value="500.5"/>
                <cvParam accession="MS:1000041" name="charge state" value="2"/>
              </selectedIon></selectedIonList></precursor></precursorList>
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
fn profile_mode_spectrum_is_centroided_on_parse() {
    let prof = collect_ok(&wrap_spectra(&profile_spectrum_xml(true)));
    assert_eq!(prof.len(), 1);
    assert_eq!(prof[0].peaks.len(), 1, "profile 3-point peak ⇒ 1 centroid");
    assert!((prof[0].peaks[0].0 - 150.1).abs() < 1e-9);

    let cent = collect_ok(&wrap_spectra(&profile_spectrum_xml(false)));
    assert_eq!(
        cent[0].peaks.len(),
        3,
        "centroided spectrum is not re-centroided"
    );
}
