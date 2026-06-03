//! End-to-end test for `msgf-rust train` and the resulting model store.
//!
//! Verifies:
//! 1. `msgf-rust train` exits 0 and writes a Parquet store.
//! 2. The store contains the trained model ID.
//! 3. A subsequent `msgf-rust --spectrum ... --model-store ... --model ...`
//!    search using that model exits 0 and produces a non-empty PIN file.

use std::path::PathBuf;
use std::process::Command;

use model_train::ModelStore;

/// Resolve a path relative to the workspace root (two levels above the crate
/// manifest dir: crates/msgf-rust → workspace root).
fn fixture(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(rel)
        .canonicalize()
        .unwrap_or_else(|e| panic!("canonicalize {rel}: {e}"))
}

#[test]
fn train_writes_model_and_search_uses_it() {
    let dir = tempfile::tempdir().expect("tempdir");
    let store_path = dir.path().join("m.parquet");
    let pin_path = dir.path().join("out.pin");

    let bsa_mgf = fixture("test-fixtures/test.mgf");
    let bsa_fasta = fixture("test-fixtures/BSA.fasta");

    // ── Step 1: run `msgf-rust train` ────────────────────────────────────────
    let train_status = Command::new(env!("CARGO_BIN_EXE_msgf-rust"))
        .arg("train")
        .arg("--spectra")
        .arg(&bsa_mgf)
        .arg("--database")
        .arg(&bsa_fasta)
        .arg("--out-store")
        .arg(&store_path)
        // Use a lenient FDR so the small BSA fixture yields confident labels.
        .arg("--train-fdr")
        .arg("0.5")
        .arg("--model-id")
        .arg("bsa_test")
        .status()
        .expect("run msgf-rust train");

    assert!(
        train_status.success(),
        "msgf-rust train should exit 0, got: {train_status}"
    );

    // ── Step 2: verify the store file exists and contains the model ──────────
    assert!(store_path.exists(), "model store should be written");

    let store = ModelStore::open(&store_path).expect("open trained model store");
    let ids = store.model_ids();
    assert!(
        ids.contains(&"bsa_test".to_string()),
        "store should contain model 'bsa_test'; found: {ids:?}"
    );

    // Verify the param loads without error.
    let param = store.load_param("bsa_test").expect("load bsa_test param");
    assert!(
        !param.partitions.is_empty(),
        "trained param should have at least one partition"
    );

    // ── Step 3: run search using the trained model ────────────────────────────
    let search_status = Command::new(env!("CARGO_BIN_EXE_msgf-rust"))
        .arg("--spectrum")
        .arg(&bsa_mgf)
        .arg("--database")
        .arg(&bsa_fasta)
        .arg("--output-pin")
        .arg(&pin_path)
        .arg("--model-store")
        .arg(&store_path)
        .arg("--model")
        .arg("bsa_test")
        .status()
        .expect("run search with trained model");

    assert!(
        search_status.success(),
        "search with trained model should exit 0, got: {search_status}"
    );
    assert!(pin_path.exists(), "PIN file should be written");

    // Verify PIN has a header + at least one data row.
    let pin_content = std::fs::read_to_string(&pin_path).expect("read PIN");
    let line_count = pin_content.lines().count();
    assert!(
        line_count > 1,
        "PIN should have header + at least 1 data row; got {line_count} line(s)"
    );
    let header = pin_content.lines().next().unwrap_or("");
    assert!(
        header.starts_with("SpecId\tLabel\tScanNr"),
        "unexpected PIN header: {header}"
    );
}
