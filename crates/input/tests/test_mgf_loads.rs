//! Load `astral-speed/test-fixtures/test.mgf` (small fixture)
//! and assert basic invariants.

use std::io::BufReader;
use std::path::PathBuf;

use input::MgfReader;

fn fixture_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join("test-fixtures/test.mgf.gz")
        .canonicalize()
        .expect("canonicalize test.mgf path")
}

#[test]
fn test_mgf_parses_completely() {
    let path = fixture_path();
    let file = input::open_buf_maybe_gz(&path)
        .unwrap_or_else(|e| panic!("open {path:?}: {e}"));
    let reader = MgfReader::new(BufReader::new(file));
    let mut count = 0;
    for result in reader {
        let s = result.unwrap_or_else(|e| panic!("parse error: {e}"));
        assert!(!s.peaks.is_empty(), "spectrum {} has no peaks", count);
        count += 1;
    }
    assert!(count > 0, "test.mgf produced 0 spectra");
}

#[test]
fn test_mgf_first_spectrum_has_expected_shape() {
    let path = fixture_path();
    let file = input::open_buf_maybe_gz(&path).unwrap();
    let reader = MgfReader::new(BufReader::new(file));
    let first = reader.into_iter().next().unwrap().unwrap();
    assert!(!first.title.is_empty(), "first spectrum has empty title");
    assert!(first.precursor_mz > 0.0, "first spectrum precursor_mz <= 0");
    assert!(first.peaks.len() >= 5, "first spectrum has < 5 peaks");
    let mzs: Vec<_> = first.peaks.iter().map(|p| p.0).collect();
    let mut sorted = mzs.clone();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    assert_eq!(mzs, sorted, "peaks not sorted ascending");
}
