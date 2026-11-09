//! Input checks that run before a search: a database that already holds
//! decoys, spectrum files of different acquisition classes in one search, and
//! an mzML that names its analyzer but no activation.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn fixture(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(rel)
        .canonicalize()
        .unwrap_or_else(|e| panic!("canonicalize {rel}: {e}"))
}

fn search(dir: &Path, spectra: &[&Path], database: &Path, extra: &[&str]) -> Output {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_andes"));
    for s in spectra {
        cmd.arg("--spectrum").arg(s);
    }
    cmd.arg("--database")
        .arg(database)
        .arg("--output-pin")
        .arg(dir.join("out.pin"))
        .args(extra)
        .output()
        .expect("run andes")
}

fn stderr_of(out: &Output) -> String {
    String::from_utf8_lossy(&out.stderr).into_owned()
}

/// The BSA fixture followed by its reversed sequence under `decoy_accession`.
fn bsa_with_reversed_copy(dir: &Path, decoy_accession: &str) -> PathBuf {
    let text = std::fs::read_to_string(fixture("test-fixtures/BSA.fasta")).expect("read BSA");
    let sequence: String = text.lines().skip(1).map(str::trim).collect();
    let reversed: String = sequence.chars().rev().collect();
    let path = dir.join("with_decoys.fasta");
    std::fs::write(&path, format!("{text}\n>{decoy_accession}\n{reversed}\n")).expect("write");
    path
}

/// The high-resolution Orbitrap fixture, inflated, with every line holding
/// `drop` removed and every `from` replaced by `to`.
fn edited_orbitrap_mzml(
    dir: &Path,
    name: &str,
    drop: Option<&str>,
    swap: Option<(&str, &str)>,
) -> PathBuf {
    use std::io::Read;
    let gz = std::fs::File::open(fixture("test-fixtures/orbitrap_lumos_120.mzML.gz")).expect("gz");
    let mut xml = String::new();
    flate2::read::GzDecoder::new(gz)
        .read_to_string(&mut xml)
        .expect("inflate");
    let mut edited = String::with_capacity(xml.len());
    let mut changed = false;
    for line in xml.lines() {
        if drop.is_some_and(|d| line.contains(d)) {
            changed = true;
            continue;
        }
        match swap {
            Some((from, to)) if line.contains(from) => {
                changed = true;
                edited.push_str(&line.replace(from, to));
            }
            _ => edited.push_str(line),
        }
        edited.push('\n');
    }
    assert!(changed, "the edit changed nothing");
    let path = dir.join(name);
    std::fs::write(&path, edited).expect("write mzML");
    path
}

#[test]
fn database_with_named_decoys_stops_the_search() {
    let dir = tempfile::tempdir().expect("tempdir");
    let fasta = bsa_with_reversed_copy(dir.path(), "rev_sp|P02769|ALBU_BOVIN");
    let spectra = fixture("test-fixtures/orbitrap_lumos_120.mzML.gz");
    let out = search(dir.path(), &[&spectra], &fasta, &[]);
    let err = stderr_of(&out);
    assert!(!out.status.success(), "the search ran:\n{err}");
    assert!(
        err.contains("the database already contains decoys"),
        "{err}"
    );
    assert!(err.contains("rev_sp|P02769|ALBU_BOVIN"), "{err}");
}

/// No decoy word in the accession: the reversed sequence is the evidence.
#[test]
fn database_with_reversed_records_stops_the_search() {
    let dir = tempfile::tempdir().expect("tempdir");
    let fasta = bsa_with_reversed_copy(dir.path(), "zz9|P02769|ALBU_BOVIN");
    let spectra = fixture("test-fixtures/orbitrap_lumos_120.mzML.gz");
    let out = search(dir.path(), &[&spectra], &fasta, &[]);
    let err = stderr_of(&out);
    assert!(!out.status.success(), "the search ran:\n{err}");
    assert!(err.contains("exact reverse of another record"), "{err}");
}

/// The fix the message names: declare the decoys and generate none.
#[test]
fn declared_decoys_without_generation_are_accepted() {
    let dir = tempfile::tempdir().expect("tempdir");
    let fasta = bsa_with_reversed_copy(dir.path(), "rev_sp|P02769|ALBU_BOVIN");
    let spectra = fixture("test-fixtures/orbitrap_lumos_120.mzML.gz");
    let out = search(
        dir.path(),
        &[&spectra],
        &fasta,
        &["--decoy-prefix", "rev", "--decoy-strategy", "none"],
    );
    let err = stderr_of(&out);
    assert!(out.status.success(), "andes failed:\n{err}");
    assert!(!err.contains("contain decoys"), "{err}");
}

#[test]
fn target_only_database_raises_no_decoy_message() {
    let dir = tempfile::tempdir().expect("tempdir");
    let spectra = fixture("test-fixtures/orbitrap_lumos_120.mzML.gz");
    let out = search(
        dir.path(),
        &[&spectra],
        &fixture("test-fixtures/ecoli.fasta"),
        &["--max-spectra", "20"],
    );
    let err = stderr_of(&out);
    assert!(out.status.success(), "andes failed:\n{err}");
    assert!(!err.contains("contain decoys"), "{err}");
}

#[test]
fn inputs_of_different_resolution_classes_stop_the_search() {
    let dir = tempfile::tempdir().expect("tempdir");
    let high = fixture("test-fixtures/orbitrap_lumos_120.mzML.gz");
    let low = edited_orbitrap_mzml(
        dir.path(),
        "ion_trap.mzML",
        None,
        Some((
            r#"accession="MS:1000484" value="" name="orbitrap""#,
            r#"accession="MS:1000264" value="" name="ion trap""#,
        )),
    );
    let out = search(
        dir.path(),
        &[&high, &low],
        &fixture("test-fixtures/BSA.fasta"),
        &[],
    );
    let err = stderr_of(&out);
    assert!(!out.status.success(), "the search ran:\n{err}");
    assert!(err.contains("differ in fragment resolution"), "{err}");
    assert!(err.contains("ion_trap.mzML"), "{err}");
}

#[test]
fn inputs_of_the_same_class_are_searched_together() {
    let dir = tempfile::tempdir().expect("tempdir");
    let high = fixture("test-fixtures/orbitrap_lumos_120.mzML.gz");
    let out = search(
        dir.path(),
        &[&high, &high],
        &fixture("test-fixtures/BSA.fasta"),
        &[],
    );
    let err = stderr_of(&out);
    assert!(out.status.success(), "andes failed:\n{err}");
    assert!(!err.contains("differs from the first input"), "{err}");
}

/// Without an activation term the analyzer is still read: the file stays on
/// the high-resolution model instead of falling back to peak inference or the
/// ion-trap default.
#[test]
fn analyzer_is_read_from_an_mzml_without_activation() {
    let dir = tempfile::tempdir().expect("tempdir");
    let mzml = edited_orbitrap_mzml(
        dir.path(),
        "no_activation.mzML",
        Some(r#"accession="MS:1000422""#),
        None,
    );
    let out = search(
        dir.path(),
        &[&mzml],
        &fixture("test-fixtures/BSA.fasta"),
        &[],
    );
    let err = stderr_of(&out);
    assert!(out.status.success(), "andes failed:\n{err}");
    assert!(err.contains("no activation recorded"), "{err}");
    let acquisition = err
        .lines()
        .find(|l| l.trim_start().starts_with("acquisition"))
        .unwrap_or_else(|| panic!("no acquisition line:\n{err}"));
    assert!(
        acquisition.contains("high-resolution MS2") && !acquisition.contains("inferred from peaks"),
        "{acquisition}"
    );
    let model = err
        .lines()
        .find(|l| l.trim_start().starts_with("model "))
        .unwrap_or_else(|| panic!("no model line:\n{err}"));
    assert!(model.contains("hcd_"), "{model}");
}
