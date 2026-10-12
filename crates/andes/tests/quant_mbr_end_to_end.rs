//! End-to-end tests for match-between-runs: the `andes quant
//! library|extract|finalize` stages and the in-process `andes search --lfq
//! --mbr`, on the three synthetic runs of `scripts/make_quant_fixtures.py
//! --mbr`.
//!
//! Rescoring uses a stand-in for the Percolator binary (a shell script passed
//! as `--percolator-bin`) that copies each PIN row's `RawScore` into the
//! `score` column and writes targets and decoys to the result files
//! Percolator would write. The stages under test only read those files.
#![cfg(unix)]

use std::collections::HashMap;
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

/// Designed quantity of HLVDEPQNLIK/2 (peptide 4): mono apex 8e5 times the
/// summed first four isotope fractions over the mono's (2.06764, pyOpenMS
/// coarse isotope pattern), times sigma (4 s) times sqrt(2 pi).
const PEPTIDE4_QUANTITY: f64 = 16_584_995.85;

const RUNS: [&str; 3] = ["mbr_run1", "mbr_run2", "mbr_run3"];

const FAKE_PERCOLATOR: &str = r#"#!/bin/sh
# Stand-in for percolator: score = the PIN's RawScore; targets and decoys go
# to the files named by --results-psms and --decoy-results-psms.
while [ $# -gt 1 ]; do
  case "$1" in
    --results-psms) T="$2"; shift 2 ;;
    --decoy-results-psms) D="$2"; shift 2 ;;
    --weights) W="$2"; shift 2 ;;
    *) shift ;;
  esac
done
PIN="$1"
awk -F'\t' -v T="$T" -v D="$D" '
NR == 1 {
  for (i = 1; i <= NF; i++) h[$i] = i
  hdr = "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds"
  print hdr > T
  print hdr > D
  next
}
{
  out = ($h["Label"] == 1) ? T : D
  print $h["SpecId"] "\t" $h["RawScore"] "\t0.001\t0.001\t" $h["Peptide"] "\t" $h["Proteins"] > out
}' "$PIN"
: > "$W"
"#;

fn fixture(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(rel)
        .canonicalize()
        .unwrap_or_else(|e| panic!("canonicalize {rel}: {e}"))
}

fn andes(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_andes"))
        .args(args)
        .output()
        .expect("run andes")
}

fn run_andes(args: &[&str]) {
    let out = andes(args);
    assert!(
        out.status.success(),
        "andes {args:?} failed: {}\n{}",
        out.status,
        String::from_utf8_lossy(&out.stderr)
    );
}

fn fake_percolator(dir: &Path) -> PathBuf {
    let path = dir.join("percolator");
    std::fs::write(&path, FAKE_PERCOLATOR).unwrap();
    std::fs::set_permissions(&path, std::fs::Permissions::from_mode(0o755)).unwrap();
    path
}

/// `--spectrum` for the three runs, the database and two threads.
fn search_args() -> Vec<String> {
    let mut v = Vec::new();
    for r in RUNS {
        v.push("--spectrum".to_string());
        v.push(
            fixture(&format!("test-fixtures/{r}.mzML.gz"))
                .display()
                .to_string(),
        );
    }
    v.extend([
        "--database".to_string(),
        fixture("test-fixtures/BSA.fasta").display().to_string(),
        "--threads".to_string(),
        "2".to_string(),
    ]);
    v
}

fn s(p: &Path) -> &str {
    p.to_str().unwrap()
}

/// mbr.tsv rows by `(precursor_key, run)`.
type MbrRows = HashMap<(String, String), Vec<String>>;

/// The header and rows of an mbr.tsv.
fn read_mbr(path: &Path) -> (Vec<String>, MbrRows) {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path:?}: {e}"));
    let mut lines = text.lines();
    let header: Vec<String> = lines
        .next()
        .unwrap()
        .split('\t')
        .map(String::from)
        .collect();
    let rows = lines
        .map(|l| {
            let f: Vec<String> = l.split('\t').map(String::from).collect();
            ((f[0].clone(), f[1].clone()), f)
        })
        .collect();
    (header, rows)
}

/// The search for the stage path: per-run artifacts and the Percolator
/// result files next to `<dir>/s.pin`.
fn staged_search(dir: &Path, perc: &Path) -> PathBuf {
    std::fs::create_dir_all(dir).unwrap();
    let q = dir.join("q");
    let pin = dir.join("s.pin");
    let mut args = search_args();
    args.extend(
        [
            "--output-pin",
            s(&pin),
            "--quant-dir",
            s(&q),
            "--rescore",
            "--percolator-bin",
            s(perc),
        ]
        .map(String::from),
    );
    run_andes(&args.iter().map(String::as_str).collect::<Vec<_>>());
    q
}

#[test]
fn stages_equal_all_in_one_and_transfer_the_untriggered_peptide() {
    let tmp = tempfile::tempdir().unwrap();
    let perc = fake_percolator(tmp.path());

    // Stage path: search, library, one extract per run, finalize.
    let staged = tmp.path().join("staged");
    let q = staged_search(&staged, &perc);
    let targets = staged.join("s.percolator.target.psms.txt");
    let decoys = staged.join("s.percolator.decoy.psms.txt");
    run_andes(&[
        "quant",
        "library",
        "--quant-dir",
        s(&q),
        "--psms",
        s(&targets),
        "--decoy-psms",
        s(&decoys),
    ]);
    assert!(q.join("library.parquet").is_file());
    assert!(q.join("alignment.parquet").is_file());
    // Extract in reverse order: the run index comes from the library.
    for r in RUNS.iter().rev() {
        run_andes(&["quant", "extract", "--quant-dir", s(&q), "--run", r]);
        assert!(q.join(format!("{r}.extract.parquet")).is_file());
    }
    let base = staged.join("out");
    let qpx = staged.join("bundle");
    run_andes(&[
        "quant",
        "finalize",
        "--quant-dir",
        s(&q),
        "-o",
        s(&base),
        "--output-parquet",
        s(&qpx),
    ]);
    assert!(qpx.join("quantms.feature.parquet").is_file());

    // All-in-one path, with a temporary quant directory.
    let aio = tmp.path().join("aio");
    std::fs::create_dir_all(&aio).unwrap();
    let pin = aio.join("a.pin");
    let mut args = search_args();
    args.extend(
        [
            "--output-pin",
            s(&pin),
            "--lfq",
            "--mbr",
            "--rescore",
            "--percolator-bin",
            s(&perc),
        ]
        .map(String::from),
    );
    let out = andes(&args.iter().map(String::as_str).collect::<Vec<_>>());
    let log = String::from_utf8_lossy(&out.stderr);
    assert!(out.status.success(), "all-in-one failed:\n{log}");
    assert!(log.contains("[mbr] library:"), "{log}");
    assert!(log.contains("[mbr] extract mbr_run3:"), "{log}");

    let staged_tsv = std::fs::read(staged.join("out.mbr.tsv")).unwrap();
    let aio_tsv = std::fs::read(aio.join("a.mbr.tsv")).unwrap();
    assert!(
        staged_tsv == aio_tsv,
        "mbr.tsv differs between the stage path and the all-in-one run"
    );
    assert_eq!(
        std::fs::read(staged.join("out.mbr_validation.tsv")).unwrap(),
        std::fs::read(aio.join("a.mbr_validation.tsv")).unwrap()
    );

    let (header, rows) = read_mbr(&aio.join("a.mbr.tsv"));
    let col = |name: &str| header.iter().position(|h| h == name).unwrap();
    let (ev, quantity) = (col("evidence_type"), col("quantity"));
    let row = |key: &str, run: &str| rows.get(&(key.to_string(), run.to_string()));

    // Peptide 4 was never triggered in run 3: it is transferred there, with
    // the designed quantity.
    let p4 = row("HLVDEPQNLIK/2", "mbr_run3").expect("peptide 4 in run 3");
    assert_eq!(p4[ev], "transfer_ms1", "{p4:?}");
    let qv: f64 = p4[quantity].parse().unwrap();
    assert!(
        (qv / PEPTIDE4_QUANTITY - 1.0).abs() < 0.10,
        "peptide 4 quantity {qv} vs designed {PEPTIDE4_QUANTITY}"
    );
    for run in ["mbr_run1", "mbr_run2"] {
        assert_eq!(row("HLVDEPQNLIK/2", run).unwrap()[ev], "identified");
    }
    // Peptide 6 has no MS1 signal anywhere and is identified in run 1 only.
    for key in ["DAFLGSFLYEYSR/2", "DAFLGSFLYEYSR/3"] {
        if let Some(r1) = row(key, "mbr_run1") {
            assert_eq!(r1[ev], "identified");
        }
        for run in ["mbr_run2", "mbr_run3"] {
            if let Some(r) = row(key, run) {
                assert_ne!(r[ev], "transfer_ms1", "{key} in {run}: {r:?}");
            }
        }
    }
    // Most untriggered precursors with signal are transferred.
    let transfers = rows.values().filter(|r| r[ev] == "transfer_ms1").count();
    assert!(transfers >= 50, "only {transfers} transfers");
}

#[test]
fn library_stage_errors_name_what_is_missing() {
    let tmp = tempfile::tempdir().unwrap();
    let perc = fake_percolator(tmp.path());
    let q = staged_search(tmp.path(), &perc);
    let targets = tmp.path().join("s.percolator.target.psms.txt");
    let decoys = tmp.path().join("s.percolator.decoy.psms.txt");

    // A run without artifacts is named.
    let out = andes(&[
        "quant",
        "library",
        "--quant-dir",
        s(&q),
        "--runs",
        "mbr_run1",
        "absent_run",
        "--psms",
        s(&targets),
        "--decoy-psms",
        s(&decoys),
    ]);
    assert!(!out.status.success());
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(err.contains("absent_run.quant.parquet"), "{err}");
    assert!(err.contains("absent_run.cache"), "{err}");

    // A Percolator file without the score column.
    let text = std::fs::read_to_string(&targets).unwrap();
    let no_score: String = text
        .lines()
        .map(|l| {
            let mut f: Vec<&str> = l.split('\t').collect();
            f.remove(1);
            f.join("\t") + "\n"
        })
        .collect();
    let stripped = tmp.path().join("no_score.tsv");
    std::fs::write(&stripped, no_score).unwrap();
    let out = andes(&[
        "quant",
        "library",
        "--quant-dir",
        s(&q),
        "--psms",
        s(&stripped),
        "--decoy-psms",
        s(&decoys),
    ]);
    assert!(!out.status.success());
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(
        err.contains("no_score.tsv") && err.contains("score"),
        "{err}"
    );

    // Extract before the library exists.
    let out = andes(&[
        "quant",
        "extract",
        "--quant-dir",
        s(&q),
        "--run",
        "mbr_run1",
    ]);
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("library.parquet"));
}

#[test]
fn mbr_needs_two_inputs_and_rescoring() {
    let tmp = tempfile::tempdir().unwrap();
    let pin = tmp.path().join("x.pin");
    let one = fixture("test-fixtures/mbr_run1.mzML.gz");
    let db = fixture("test-fixtures/BSA.fasta");
    let base = [
        "--spectrum",
        s(&one),
        "--database",
        s(&db),
        "--output-pin",
        s(&pin),
        "--lfq",
        "--mbr",
    ];
    let out = andes(&[&base[..], &["--rescore-native"]].concat());
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("two or more"));

    let two = fixture("test-fixtures/mbr_run2.mzML.gz");
    let out = andes(&[&base[..], &["--spectrum", s(&two)]].concat());
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("--rescore"));
}
