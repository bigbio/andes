//! Load `astral-speed/src/test/resources/iprg-2013/F13.mgf` (1,406
//! spectra) and assert count + wall-time budget.

use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;
use std::time::Instant;

use input::MgfReader;

fn fixture_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .join("src/test/resources/iprg-2013/F13.mgf")
        .canonicalize()
        .expect("canonicalize F13.mgf path")
}

#[test]
fn f13_mgf_parses_1406_spectra() {
    let path = fixture_path();
    let file = File::open(&path).unwrap_or_else(|e| panic!("open {path:?}: {e}"));
    let reader = MgfReader::new(BufReader::new(file));

    let start = Instant::now();
    let count = reader.into_iter().filter_map(|r| r.ok()).count();
    let elapsed = start.elapsed();

    assert_eq!(count, 1406, "expected 1406 spectra, got {count}");
    assert!(
        elapsed.as_secs_f32() < 3.0,
        "F13.mgf parse took {:.2}s, target < 3s", elapsed.as_secs_f32()
    );
}
