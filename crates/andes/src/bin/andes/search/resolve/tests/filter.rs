//! Tests of the isobaric peak-filter setting.

use super::isobaric_peak_filter;

#[test]
fn high_res_turns_the_isobaric_filter_off_unless_set() {
    assert_eq!(isobaric_peak_filter(None, false), None);
    assert_eq!(isobaric_peak_filter(None, true), Some((0.0, 0)));
    assert_eq!(
        isobaric_peak_filter(Some((100.0, 20)), true),
        Some((100.0, 20))
    );
    assert_eq!(
        isobaric_peak_filter(Some((50.0, 10)), false),
        Some((50.0, 10))
    );
}
