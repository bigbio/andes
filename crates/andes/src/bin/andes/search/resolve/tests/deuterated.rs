//! Tests of the deuterated-tag warning.

use super::checks::deuterated_tag_warning;
use model::Tolerance;
use quant::Plex;

#[test]
fn warns_only_for_deuterated_kits_under_a_narrow_window() {
    let w = deuterated_tag_warning(Some(Plex::Tmt35), Tolerance::Ppm(5.0)).unwrap();
    assert!(w.contains("8.8 ppm") && w.contains("5ppm"), "{w}");
    assert!(deuterated_tag_warning(Some(Plex::Tmt32), Tolerance::Da(0.005)).is_some());
    assert!(deuterated_tag_warning(Some(Plex::Tmt35), Tolerance::Ppm(10.0)).is_none());
    assert!(deuterated_tag_warning(Some(Plex::Tmt35), Tolerance::Da(0.02)).is_none());
    assert!(deuterated_tag_warning(Some(Plex::Tmt18), Tolerance::Ppm(5.0)).is_none());
    assert!(deuterated_tag_warning(None, Tolerance::Ppm(5.0)).is_none());
}
