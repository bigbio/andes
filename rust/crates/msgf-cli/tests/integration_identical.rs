//! Integration tests for msgf-diff's file-level identical/different check.

use std::fs;
use std::process::Command;

/// Locate the built msgf-diff binary in target/debug.
fn diff_bin() -> std::path::PathBuf {
    let manifest = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    manifest
        .parent()
        .unwrap() // crates/
        .parent()
        .unwrap() // workspace root
        .join("target/debug/msgf-diff")
}

#[test]
fn identical_files_exit_zero() {
    let dir = tempfile::tempdir().unwrap();
    let a = dir.path().join("a.pin");
    let b = dir.path().join("b.pin");
    fs::write(&a, "header\nrow\n").unwrap();
    fs::write(&b, "header\nrow\n").unwrap();

    let status = Command::new(diff_bin())
        .args(["compare", a.to_str().unwrap(), b.to_str().unwrap()])
        .status()
        .expect("spawn msgf-diff");
    assert_eq!(status.code(), Some(0));
}

#[test]
fn different_files_exit_one() {
    let dir = tempfile::tempdir().unwrap();
    let a = dir.path().join("a.pin");
    let b = dir.path().join("b.pin");
    fs::write(&a, "header\nrow\n").unwrap();
    fs::write(&b, "header\nDIFFERENT\n").unwrap();

    let status = Command::new(diff_bin())
        .args(["compare", a.to_str().unwrap(), b.to_str().unwrap()])
        .status()
        .expect("spawn msgf-diff");
    assert_eq!(status.code(), Some(1));
}
