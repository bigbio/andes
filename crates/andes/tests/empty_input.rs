//! An input file with no MS2 spectra is a clear error naming the file, not a
//! search of zero spectra on the default acquisition parameters.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const EMPTY_MZML: &str = r#"<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml" version="1.1.0">
  <run id="empty">
    <spectrumList count="0" defaultDataProcessingRef="dp"/>
  </run>
</mzML>
"#;

fn bsa_fasta() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../test-fixtures/BSA.fasta")
        .canonicalize()
        .expect("BSA fixture")
}

fn search(spectrum: &Path, dir: &Path) -> Output {
    Command::new(env!("CARGO_BIN_EXE_andes"))
        .arg("--spectrum")
        .arg(spectrum)
        .arg("--database")
        .arg(bsa_fasta())
        .arg("--output-pin")
        .arg(dir.join("out.pin"))
        .output()
        .expect("run andes")
}

fn assert_clear_error(out: &Output, file: &Path) {
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(!out.status.success(), "empty input must fail:\n{stderr}");
    let name = file.file_name().unwrap().to_str().unwrap();
    assert!(
        stderr.contains("no MS2 spectra could be read") && stderr.contains(name),
        "error should name the file and the missing spectra:\n{stderr}"
    );
}

#[test]
fn empty_mgf_is_an_error_naming_the_file() {
    let dir = tempfile::tempdir().expect("tempdir");
    let mgf = dir.path().join("empty_run.mgf");
    std::fs::write(&mgf, "").unwrap();
    assert_clear_error(&search(&mgf, dir.path()), &mgf);
}

#[test]
fn mzml_without_spectra_is_an_error_naming_the_file() {
    let dir = tempfile::tempdir().expect("tempdir");
    let mzml = dir.path().join("empty_run.mzML");
    std::fs::write(&mzml, EMPTY_MZML).unwrap();
    assert_clear_error(&search(&mzml, dir.path()), &mzml);
}

/// One empty file among several is not hidden by the spectra of the others.
#[test]
fn one_empty_file_among_several_is_an_error() {
    let dir = tempfile::tempdir().expect("tempdir");
    let good = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../test-fixtures/test.mgf.gz")
        .canonicalize()
        .expect("MGF fixture");
    let empty = dir.path().join("empty_run.mgf");
    std::fs::write(&empty, "").unwrap();
    let out = Command::new(env!("CARGO_BIN_EXE_andes"))
        .arg("--spectrum")
        .arg(&good)
        .arg("--spectrum")
        .arg(&empty)
        .arg("--database")
        .arg(bsa_fasta())
        .arg("--output-pin")
        .arg(dir.path().join("out.pin"))
        .output()
        .expect("run andes");
    assert_clear_error(&out, &empty);
}
