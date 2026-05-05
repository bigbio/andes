//! End-to-end smoke test: invoke msgf-rust on BSA + test.mgf and verify
//! the PIN and TSV outputs exist with sensible content.

use std::path::PathBuf;
use std::process::Command;

/// Resolve a path relative to the workspace root (three levels above the
/// cli crate's manifest directory: cli → crates → rust → astral-speed).
fn fixture(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .join(rel)
        .canonicalize()
        .unwrap_or_else(|e| panic!("canonicalize {rel}: {e}"))
}

#[test]
fn cli_runs_end_to_end_on_bsa_test_mgf() {
    let dir = tempfile::tempdir().expect("tempdir");
    let pin_path = dir.path().join("rust.pin");
    let tsv_path = dir.path().join("rust.tsv");

    let status = Command::new(env!("CARGO_BIN_EXE_msgf-rust"))
        .arg("--spectrum")
        .arg(fixture("src/test/resources/test.mgf"))
        .arg("--database")
        .arg(fixture("src/test/resources/BSA.fasta"))
        .arg("--output-pin")
        .arg(&pin_path)
        .arg("--output-tsv")
        .arg(&tsv_path)
        .arg("--decoy-prefix")
        .arg("XXX_")
        .status()
        .expect("run msgf-rust");

    assert!(status.success(), "msgf-rust exit code: {status}");
    assert!(pin_path.exists(), "PIN output not written");
    assert!(tsv_path.exists(), "TSV output not written");

    // Validate PIN header and content.
    let pin_content = std::fs::read_to_string(&pin_path).unwrap();
    assert!(
        pin_content.lines().count() > 1,
        "PIN should have header + at least 1 row"
    );
    let pin_header = pin_content.lines().next().unwrap();
    assert!(
        pin_header.starts_with("SpecId\tLabel\tScanNr"),
        "unexpected PIN header: {pin_header}"
    );

    // Validate TSV header and content.
    let tsv_content = std::fs::read_to_string(&tsv_path).unwrap();
    assert!(
        tsv_content.lines().count() > 1,
        "TSV should have header + at least 1 row"
    );
    let tsv_header = tsv_content.lines().next().unwrap();
    assert!(
        tsv_header.starts_with("#SpecFile\tSpecID\tScanNum"),
        "unexpected TSV header: {tsv_header}"
    );
}
