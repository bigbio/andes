//! The acquisition context on metadata-less input: an MGF carries no analyzer,
//! so the fragment resolution is read from its peaks, and the selected model
//! and the feature tolerances follow that class.

use std::path::PathBuf;
use std::process::Command;

fn fixture(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(rel)
        .canonicalize()
        .unwrap_or_else(|e| panic!("canonicalize {rel}: {e}"))
}

/// Run a search and return its stderr.
fn run(spectrum: &str, database: &str, extra: &[&str]) -> String {
    let dir = tempfile::tempdir().expect("tempdir");
    let out = Command::new(env!("CARGO_BIN_EXE_andes"))
        .arg("--spectrum")
        .arg(fixture(spectrum))
        .arg("--database")
        .arg(fixture(database))
        .arg("--output-pin")
        .arg(dir.path().join("out.pin"))
        .args(extra)
        .output()
        .expect("run andes");
    let stderr = String::from_utf8_lossy(&out.stderr).into_owned();
    assert!(out.status.success(), "andes failed:\n{stderr}");
    stderr
}

fn banner_line<'a>(stderr: &'a str, key: &str) -> &'a str {
    stderr
        .lines()
        .find(|l| l.trim_start().starts_with(key))
        .unwrap_or_else(|| panic!("no `{key}` line in:\n{stderr}"))
}

/// Orbitrap HCD spectra exported to MGF without any activation or analyzer
/// line: the peaks are read as high resolution and the high-res HCD model is
/// selected, instead of the ion-trap default.
#[test]
fn metadataless_orbitrap_mgf_is_read_as_high_resolution() {
    let err = run(
        "test-fixtures/glyco_fixture.mgf.gz",
        "test-fixtures/glyco_fixture.fasta",
        &[],
    );
    let acq = banner_line(&err, "acquisition");
    assert!(
        acq.contains("high-resolution MS2 [inferred from peaks]"),
        "{acq}"
    );
    let model = banner_line(&err, "model ");
    assert!(model.contains("hcd_qexactive_tryp"), "{model}");
    assert!(
        !err.contains("assuming CID / low-res / 0.5 Da"),
        "the low-res fallback warning fired:\n{err}"
    );
}

/// An explicit `--fragment-tol-da` names the low-resolution class and wins
/// over the peak inference.
#[test]
fn fragment_tol_flag_wins_over_the_peak_inference() {
    let err = run(
        "test-fixtures/glyco_fixture.mgf.gz",
        "test-fixtures/glyco_fixture.fasta",
        &["--fragment-tol-da", "0.5"],
    );
    let acq = banner_line(&err, "acquisition");
    assert!(
        acq.contains("low-resolution MS2 [from --fragment-tol-*]"),
        "{acq}"
    );
    assert!(banner_line(&err, "model ").contains("cid_lowres_tryp"));
}

/// Ion-trap spectra stay on the low-resolution model.
#[test]
fn metadataless_ion_trap_mgf_stays_low_resolution() {
    let err = run(
        "test-fixtures/test.mgf.gz",
        "test-fixtures/BSA.fasta",
        &["--max-spectra", "200"],
    );
    let acq = banner_line(&err, "acquisition");
    assert!(acq.contains("low-resolution MS2"), "{acq}");
    assert!(banner_line(&err, "model ").contains("cid_lowres_tryp"));
}

/// Without an override a high-resolution run measures its fragment error
/// before scoring. This fixture has too few confident PSMs, so the PIN
/// features keep the 20 ppm class default and say so.
#[test]
fn high_res_run_measures_or_reports_its_fragment_error() {
    let err = run(
        "test-fixtures/orbitrap_lumos_120.mzML.gz",
        "test-fixtures/BSA.fasta",
        &[],
    );
    let line = err
        .lines()
        .find(|l| l.starts_with("Fragment error:"))
        .unwrap_or_else(|| panic!("no measurement line:\n{err}"));
    assert!(
        line.contains("not measured") && line.contains("Ppm(20.0)"),
        "{line}"
    );
}
