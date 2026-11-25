//! Tests of the crosslink precursor snap.

use super::outputs::xl_precursor_snap;

#[test]
fn the_precursor_snap_is_decided_per_file_for_ms2_crosslink_searches() {
    // An MGF next to an mzML: each file decides alone, whatever the input order.
    let files = [(false, false), (true, false), (false, true)];
    let got: Vec<bool> = files
        .iter()
        .map(|&(m, r)| xl_precursor_snap(true, 2, m, r))
        .collect();
    assert_eq!(got, [false, true, true]);
    let reversed: Vec<bool> = files
        .iter()
        .rev()
        .map(|&(m, r)| xl_precursor_snap(true, 2, m, r))
        .collect();
    assert_eq!(reversed, [true, true, false]);
    // An MS3 search keeps its own reader; no crosslink search, no snap.
    assert!(!xl_precursor_snap(true, 3, true, false));
    assert!(!xl_precursor_snap(false, 2, true, false));
}
