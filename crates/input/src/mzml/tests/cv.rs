//! Tests of the controlled-vocabulary helpers.

use crate::mzml::cv::is_minute_unit;

#[test]
fn minute_unit_by_accession_legacy_term_or_name() {
    assert!(is_minute_unit("UO:0000031", "minute"));
    assert!(is_minute_unit("MS:1000038", ""));
    assert!(is_minute_unit("", "minute"));
    assert!(is_minute_unit("", "Min"));
    assert!(!is_minute_unit("UO:0000010", "second"));
    assert!(!is_minute_unit("", "second"));
    assert!(!is_minute_unit("", ""));
}
