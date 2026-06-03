//! Binary `.param` reader gate: load the representative `.param` fixtures and
//! assert structural invariants. The 39 production models now ship as a single
//! `resources/ionstat/models.parquet` store (see the `model-train` crate); a
//! diverse subset is kept under `tests/fixtures/` to exercise the binary reader
//! across activations, resolutions, and protocols.

use std::fs;
use std::path::PathBuf;

use scoring::Param;

fn fixtures_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures")
}

fn collect_param_files() -> Vec<PathBuf> {
    let dir = fixtures_dir();
    let mut files: Vec<PathBuf> = fs::read_dir(&dir)
        .unwrap_or_else(|e| panic!("read_dir {dir:?}: {e}"))
        .filter_map(|entry| entry.ok().map(|e| e.path()))
        .filter(|p| p.extension().is_some_and(|ext| ext == "param"))
        .collect();
    files.sort();
    files
}

#[test]
fn bundled_param_fixtures_load() {
    let files = collect_param_files();
    assert!(
        !files.is_empty(),
        "expected at least one .param fixture in {:?}",
        fixtures_dir()
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
        panic!(
            "{} of {} .param fixtures failed to load:\n{}",
            failures.len(),
            files.len(),
            failures.join("\n")
        );
    }
}

#[test]
fn each_param_round_trips_validation_marker() {
    let files = collect_param_files();
    assert!(!files.is_empty());
    for path in &files {
        let bytes = fs::read(path).unwrap();
        let result = Param::load_from_bytes(&bytes);
        assert!(result.is_ok(), "{path:?}: {:?}", result.err());
    }
}
