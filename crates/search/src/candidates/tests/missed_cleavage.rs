use super::missed_cleavages_between;

#[test]
fn counts_only_strictly_internal_sites() {
    let cp = [0u32, 5, 10, 15, 20];
    assert_eq!(missed_cleavages_between(&cp, 0, 20), 3);
    assert_eq!(missed_cleavages_between(&cp, 5, 10), 0);
    assert_eq!(missed_cleavages_between(&cp, 0, 10), 1);
    assert_eq!(missed_cleavages_between(&cp, 3, 17), 3);
    assert_eq!(missed_cleavages_between(&cp, 10, 11), 0);
}

#[test]
fn matches_strict_offset_semantics() {
    let cp = [2u32, 9, 14, 30];
    assert_eq!(missed_cleavages_between(&cp, cp[0], cp[1]), 0);
    assert_eq!(missed_cleavages_between(&cp, cp[0], cp[2]), 1);
    assert_eq!(missed_cleavages_between(&cp, cp[0], cp[3]), 2);
}
