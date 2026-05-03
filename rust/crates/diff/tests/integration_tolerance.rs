//! Tolerance-comparison tests.

use std::fs;
use std::process::Command;

fn diff_bin() -> std::path::PathBuf {
    let manifest = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    manifest
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("target/debug/msgf-diff")
}

const HEADER: &str = "SpecID\tScanNum\tCharge\tPeptide\tSpecEValue\n";

#[test]
fn tolerance_within_threshold_passes() {
    let dir = tempfile::tempdir().unwrap();
    let a = dir.path().join("a.pin");
    let b = dir.path().join("b.pin");
    fs::write(&a, format!("{HEADER}s1\t1\t2\tPEPTIDE\t1.0e-5\n")).unwrap();
    fs::write(&b, format!("{HEADER}s1\t1\t2\tPEPTIDE\t1.0005e-5\n")).unwrap();

    let status = Command::new(diff_bin())
        .args([
            "compare",
            a.to_str().unwrap(),
            b.to_str().unwrap(),
            "--tolerance",
            "SpecEValue:1e-3",
        ])
        .status()
        .expect("spawn msgf-diff");
    assert_eq!(status.code(), Some(0));
}

#[test]
fn tolerance_exceeded_fails() {
    let dir = tempfile::tempdir().unwrap();
    let a = dir.path().join("a.pin");
    let b = dir.path().join("b.pin");
    fs::write(&a, format!("{HEADER}s1\t1\t2\tPEPTIDE\t1.0e-5\n")).unwrap();
    fs::write(&b, format!("{HEADER}s1\t1\t2\tPEPTIDE\t1.05e-5\n")).unwrap();

    let output = Command::new(diff_bin())
        .args([
            "compare",
            a.to_str().unwrap(),
            b.to_str().unwrap(),
            "--tolerance",
            "SpecEValue:1e-3",
        ])
        .output()
        .expect("spawn msgf-diff");
    assert_eq!(output.status.code(), Some(1));
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("SpecEValue"),
        "stderr should name field: {stderr}"
    );
}

#[test]
fn missing_join_key_row_fails_with_clear_message() {
    let dir = tempfile::tempdir().unwrap();
    let a = dir.path().join("a.pin");
    let b = dir.path().join("b.pin");
    fs::write(&a, format!("{HEADER}s1\t1\t2\tPEPTIDE\t1.0e-5\n")).unwrap();
    fs::write(&b, format!("{HEADER}s2\t1\t2\tPEPTIDE\t1.0e-5\n")).unwrap();

    let output = Command::new(diff_bin())
        .args([
            "compare",
            a.to_str().unwrap(),
            b.to_str().unwrap(),
            "--tolerance",
            "SpecEValue:1e-3",
        ])
        .output()
        .expect("spawn msgf-diff");
    assert_eq!(output.status.code(), Some(1));
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("only in"),
        "stderr should report missing rows: {stderr}"
    );
}
