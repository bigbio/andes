//! `andes database` end to end: a search of the database it writes gives the same PIN and TSV,
//! byte for byte, as the search that generates its decoys with the same decoy key, and a search
//! with other decoy settings or another digest is refused with the setting named.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("workspace root")
}

fn andes(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_andes"))
        .args(args)
        .output()
        .expect("run andes")
}

fn ok(out: Output, what: &str) -> String {
    let stderr = String::from_utf8_lossy(&out.stderr).into_owned();
    assert!(out.status.success(), "{what} failed: {stderr}");
    stderr
}

fn s(p: &Path) -> &str {
    p.to_str().expect("utf-8 path")
}

fn build(dir: &Path, fasta: &Path, extra: &[&str]) -> PathBuf {
    let td = dir.join("td.fasta");
    let mut args = vec!["database", "--fasta", s(fasta), "--output", s(&td), "--threads", "2"];
    args.extend_from_slice(extra);
    ok(andes(&args), "andes database");
    td
}

/// Searches `database` into `dir` and returns the PIN, the TSV and the stderr.
fn search(dir: &Path, spectrum: &Path, database: &Path, flags: &[&str]) -> (Vec<u8>, Vec<u8>, String) {
    std::fs::create_dir_all(dir).unwrap();
    let (pin, tsv) = (dir.join("run.pin"), dir.join("run.tsv"));
    let mut args = vec![
        "--spectrum",
        s(spectrum),
        "--database",
        s(database),
        "--output-pin",
        s(&pin),
        "--output-tsv",
        s(&tsv),
        "--threads",
        "1",
    ];
    args.extend_from_slice(flags);
    let stderr = ok(andes(&args), "search");
    (
        std::fs::read(&pin).unwrap(),
        std::fs::read(&tsv).unwrap(),
        stderr,
    )
}

fn decoy_key(stderr: &str) -> String {
    stderr
        .lines()
        .find_map(|l| l.split("decoy key ").nth(1))
        .and_then(|rest| rest.split_whitespace().next())
        .unwrap_or_else(|| panic!("no decoy key in the log:\n{stderr}"))
        .trim_end_matches(';')
        .to_string()
}

fn assert_same_search(spectrum: &str, fasta: &str, flags: &[&str], database_flags: &[&str]) {
    let root = root();
    let (spectrum, fasta) = (root.join(spectrum), root.join(fasta));
    let dir = tempfile::tempdir().unwrap();
    let td = build(dir.path(), &fasta, database_flags);
    assert!(dir.path().join("td.database.parquet").is_file());
    assert!(dir.path().join("td.overlap.parquet").is_file());
    let (pin_a, tsv_a, log_a) = search(&dir.path().join("own"), &spectrum, &fasta, flags);
    let (pin_b, tsv_b, log_b) = search(&dir.path().join("prebuilt"), &spectrum, &td, flags);
    assert!(log_b.contains("written by andes database"), "{log_b}");
    assert!(!log_b.contains("decoy check:"), "the prebuilt search generated decoys");
    assert_eq!(decoy_key(&log_a), decoy_key(&log_b), "decoy keys differ");
    assert!(!pin_a.is_empty());
    assert!(pin_a == pin_b, "PIN differs between the generated and the prebuilt decoys");
    assert!(tsv_a == tsv_b, "TSV differs between the generated and the prebuilt decoys");
}

#[test]
fn prebuilt_decoys_search_byte_identical_bsa() {
    assert_same_search(
        "test-fixtures/test.mgf.gz",
        "test-fixtures/BSA.fasta",
        &["--predict", "none", "--max-missed-cleavages", "1"],
        &["--max-missed-cleavages", "1"],
    );
}

#[test]
fn prebuilt_decoys_search_byte_identical_ecoli_refine() {
    assert_same_search(
        "test-fixtures/astral_pxd070049_slice.mzML.gz",
        "test-fixtures/ecoli.fasta",
        &[
            "--precursor-cal",
            "off",
            "--fragment-index",
            "off",
            "--refine",
            "--psm-fdr",
            "1.0",
            "--predict",
            "none",
            "--max-missed-cleavages",
            "1",
            "--decoy-strategy",
            "reverse",
        ],
        &["--max-missed-cleavages", "1", "--decoy-strategy", "reverse"],
    );
}

#[test]
fn a_search_with_other_settings_is_refused_naming_them() {
    let root = root();
    let (spectrum, fasta) = (
        root.join("test-fixtures/test.mgf.gz"),
        root.join("test-fixtures/BSA.fasta"),
    );
    let dir = tempfile::tempdir().unwrap();
    let td = build(
        dir.path(),
        &fasta,
        &["--decoy-strategy", "reverse", "--check-digest", "lysc/fully/6-50/2"],
    );
    let pin = dir.path().join("x.pin");
    let refused = |flags: &[&str]| {
        let mut args = vec![
            "--spectrum",
            s(&spectrum),
            "--database",
            s(&td),
            "--output-pin",
            s(&pin),
            "--predict",
            "none",
        ];
        args.extend_from_slice(flags);
        let out = andes(&args);
        assert!(!out.status.success(), "accepted {flags:?}");
        String::from_utf8_lossy(&out.stderr).into_owned()
    };
    let e = refused(&["--decoy-strategy", "shuffle"]);
    assert!(e.contains("--decoy-strategy: database reverse, search shuffle"), "{e}");
    let e = refused(&["--decoy-strategy", "reverse", "--max-missed-cleavages", "3"]);
    assert!(e.contains("--max-missed-cleavages: database 2, search 3"), "{e}");
    let out = andes(&[
        "--spectrum",
        s(&spectrum),
        "--database",
        s(&td),
        "--output-pin",
        s(&pin),
        "--predict",
        "none",
        "--decoy-strategy",
        "reverse",
        "--enzyme",
        "lysc",
    ]);
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(out.status.success(), "the checked Lys-C digest is refused: {stderr}");

    std::fs::write(&td, ">P1\nPEPTIDEK\n").unwrap();
    let e = refused(&["--decoy-strategy", "reverse"]);
    assert!(e.contains("changed after andes database wrote it"), "{e}");
}

#[test]
fn database_refuses_a_fasta_that_holds_decoys() {
    let dir = tempfile::tempdir().unwrap();
    let fasta = dir.path().join("in.fasta");
    std::fs::write(&fasta, ">P1\nPEPTIDEKAAAK\n>XXX_P1\nEDITPEPKAAAK\n").unwrap();
    let td = dir.path().join("td.fasta");
    let out = andes(&["database", "--fasta", s(&fasta), "--output", s(&td)]);
    assert!(!out.status.success());
    let e = String::from_utf8_lossy(&out.stderr);
    assert!(e.contains("already carry the decoy prefix"), "{e}");
}
