//! Phase 2 exit gate: load every bundled `.param` file and assert
//! structural invariants. Path is resolved via `CARGO_MANIFEST_DIR`
//! (`crates/engine/`) walked up to `astral-speed/`, then into
//! `src/main/resources/ionstat/`.

use std::fs;
use std::path::PathBuf;

use scoring::Param;

fn ionstat_dir() -> PathBuf {
    // CARGO_MANIFEST_DIR = astral-speed/rust/crates/engine
    // ../../../  → astral-speed/
    // src/main/resources/ionstat/
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .join("src/main/resources/ionstat")
        .canonicalize()
        .expect("canonicalize ionstat path")
}

fn collect_param_files() -> Vec<PathBuf> {
    let dir = ionstat_dir();
    let mut files: Vec<PathBuf> = fs::read_dir(&dir)
        .unwrap_or_else(|e| panic!("read_dir {dir:?}: {e}"))
        .filter_map(|entry| entry.ok().map(|e| e.path()))
        .filter(|p| p.extension().is_some_and(|ext| ext == "param"))
        .collect();
    files.sort();
    files
}

#[test]
fn all_39_bundled_param_files_load() {
    let files = collect_param_files();
    assert_eq!(
        files.len(), 39,
        "expected 39 .param files in {:?}, found {}",
        ionstat_dir(), files.len()
    );

    let mut failures = Vec::new();
    for path in &files {
        let bytes = fs::read(path).unwrap_or_else(|e| panic!("read {path:?}: {e}"));
        match Param::load_from_bytes(&bytes) {
            Ok(param) => {
                if param.version <= 0 {
                    failures.push(format!("{path:?}: bad version {}", param.version));
                }
                if param.partitions.is_empty() {
                    failures.push(format!("{path:?}: no partitions"));
                }
                if param.charge_hist.is_empty() {
                    failures.push(format!("{path:?}: empty charge_hist"));
                }
                if param.max_rank < 0 {
                    failures.push(format!("{path:?}: negative max_rank {}", param.max_rank));
                }
            }
            Err(e) => {
                failures.push(format!("{path:?}: load failed: {e}"));
            }
        }
    }

    if !failures.is_empty() {
        panic!("{} of {} .param files failed to load:\n{}",
            failures.len(), files.len(), failures.join("\n"));
    }
}

#[test]
fn each_param_round_trips_validation_marker() {
    let files = collect_param_files();
    for path in &files {
        let bytes = fs::read(path).unwrap();
        let result = Param::load_from_bytes(&bytes);
        assert!(result.is_ok(), "{path:?}: {:?}", result.err());
    }
}
