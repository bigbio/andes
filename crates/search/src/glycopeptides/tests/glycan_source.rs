use super::{glycan_source, GlycanSource};

#[test]
fn hcd_with_oxonium_uses_the_ion_index() {
    assert_eq!(glycan_source(false, true, false), GlycanSource::IonIndex);
}

#[test]
fn hcd_without_oxonium_gets_nothing() {
    assert_eq!(glycan_source(false, false, false), GlycanSource::None);
}

#[test]
fn etd_without_oxonium_falls_back_to_the_full_list() {
    assert_eq!(glycan_source(false, false, true), GlycanSource::FullList);
}

#[test]
fn etd_with_oxonium_behaves_like_any_other_scan() {
    assert_eq!(glycan_source(false, true, true), GlycanSource::IonIndex);
    assert_eq!(glycan_source(true, true, true), GlycanSource::FullList);
}

#[test]
fn full_glycan_db_uses_the_whole_list_when_oxonium_fired() {
    assert_eq!(glycan_source(true, true, false), GlycanSource::FullList);
}

#[test]
fn full_glycan_db_without_oxonium_on_hcd_still_gets_nothing() {
    assert_eq!(glycan_source(true, false, false), GlycanSource::None);
}
