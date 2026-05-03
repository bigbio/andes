//! Schema-level diff tests.

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

#[test]
fn different_schema_exits_three_with_diff_message() {
    let dir = tempfile::tempdir().unwrap();
    let a = dir.path().join("a.pin");
    let b = dir.path().join("b.pin");
    fs::write(&a, "SpecID\tCharge\tPeptide\tSpecEValue\n").unwrap();
    fs::write(&b, "SpecID\tCharge\tPeptide\tQValue\n").unwrap();

    let output = Command::new(diff_bin())
        .args(["compare", a.to_str().unwrap(), b.to_str().unwrap()])
        .output()
        .expect("spawn msgf-diff");
    assert_eq!(output.status.code(), Some(3));
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("SpecEValue"),
        "stderr should name the missing column: {stderr}"
    );
    assert!(
        stderr.contains("QValue"),
        "stderr should name the new column: {stderr}"
    );
}

#[test]
fn empty_first_line_is_invalid_pin() {
    let dir = tempfile::tempdir().unwrap();
    let a = dir.path().join("a.pin");
    let b = dir.path().join("b.pin");
    fs::write(&a, "").unwrap();
    fs::write(&b, "SpecID\n").unwrap();

    let status = Command::new(diff_bin())
        .args(["compare", a.to_str().unwrap(), b.to_str().unwrap()])
        .status()
        .expect("spawn msgf-diff");
    assert_eq!(status.code(), Some(2)); // 2 = input error
}
