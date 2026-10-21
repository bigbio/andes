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
use std::fs::File;
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use arrow::array::{Array, Int16Array, StringArray};
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

/// Designed quantity of HLVDEPQNLIK/2 (peptide 4): mono apex 8e5 times the
/// summed first four isotope fractions over the mono's (2.06764, a coarse
/// isotope pattern), times sigma (4 s) times sqrt(2 pi).
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
    fake_percolator_with_pep(dir, "0.001")
}

/// [`FAKE_PERCOLATOR`] writing `pep` as every PSM's PEP.
fn fake_percolator_with_pep(dir: &Path, pep: &str) -> PathBuf {
    let path = dir.join(format!("percolator_pep{pep}"));
    let script = FAKE_PERCOLATOR.replace("\\t0.001\\t0.001\\t", &format!("\\t0.001\\t{pep}\\t"));
    assert!(script.contains(&format!("\\t0.001\\t{pep}\\t")));
    std::fs::write(&path, script).unwrap();
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

/// mbr.parquet rows by `(precursor_key, run)`, as text cells (a null is
/// `NA`).
type MbrRows = HashMap<(String, String), Vec<String>>;

/// The header and rows of an mbr.parquet.
fn read_mbr(path: &Path) -> (Vec<String>, MbrRows) {
    let header = output::mbr_out::MBR_COLUMNS.map(String::from).to_vec();
    let num = |v: Option<f64>| v.map_or_else(|| "NA".to_string(), |v| v.to_string());
    let rows = output::mbr_out::read_mbr_parquet(path)
        .unwrap_or_else(|e| panic!("{e}"))
        .into_iter()
        .map(|r| {
            let f = vec![
                r.precursor_key.clone(),
                r.run.clone(),
                r.evidence_type,
                num(r.quantity),
                num(r.quantity_se),
                num(r.quantity_normalized),
                num(r.transfer_q),
                num(r.quant_q),
                num(r.apex_rt_s),
                num(Some(r.library_q)),
            ];
            ((r.precursor_key, r.run), f)
        })
        .collect();
    (header, rows)
}

/// Whether two runs wrote the same MBR tables (`<a>`/`<b>` are output bases).
fn same_mbr_tables(a: &Path, b: &Path) -> bool {
    let with = |p: &Path, s: &str| PathBuf::from(format!("{}{s}", p.display()));
    let rows = |p: &Path| output::mbr_out::read_mbr_parquet(&with(p, ".mbr.parquet")).unwrap();
    let val = |p: &Path| {
        output::mbr_out::read_mbr_validation(&with(p, ".mbr_validation.parquet")).unwrap()
    };
    rows(a) == rows(b) && val(a) == val(b)
}

/// The column-data bytes of a Parquet file: everything between the leading
/// magic and the footer (the footer holds a per-write uuid and date).
fn parquet_column_data(path: &Path) -> Vec<u8> {
    let b = std::fs::read(path).unwrap_or_else(|e| panic!("{path:?}: {e}"));
    let n = b.len();
    let footer = u32::from_le_bytes(b[n - 8..n - 4].try_into().unwrap()) as usize;
    b[4..n - 8 - footer].to_vec()
}

/// `(run_file_name, peptidoform, charge, evidence_type)` of every feature
/// parquet row.
fn feature_rows(path: &Path) -> Vec<(String, String, String, String)> {
    let builder = ParquetRecordBatchReaderBuilder::try_new(File::open(path).unwrap()).unwrap();
    let mut out = Vec::new();
    for b in builder.build().unwrap() {
        let b = b.unwrap();
        let col = |name: &str| b.column_by_name(name).unwrap_or_else(|| panic!("{name}"));
        let text = |name: &str| -> Vec<String> {
            let s = col(name).as_any().downcast_ref::<StringArray>().unwrap();
            (0..s.len()).map(|i| s.value(i).to_string()).collect()
        };
        let z = col("charge").as_any().downcast_ref::<Int16Array>().unwrap();
        let z: Vec<String> = (0..z.len()).map(|i| z.value(i).to_string()).collect();
        let (run, pf, ev) = (
            text("run_file_name"),
            text("peptidoform"),
            text("evidence_type"),
        );
        for i in 0..b.num_rows() {
            out.push((run[i].clone(), pf[i].clone(), z[i].clone(), ev[i].clone()));
        }
    }
    out
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

/// `--mbr-fdr` of the transfer test. The fixture's transfer folds hold about
/// 30 targets each, so with the K-decoy estimate's `+1` correction no transfer
/// q-value can fall below about `1 / (2 * 30)`.
const FIXTURE_MBR_FDR: &str = "0.05";

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
        "--mbr-fdr",
        FIXTURE_MBR_FDR,
        "--quant-feature-fdr",
        "0.05",
    ]);
    assert!(qpx.join("quantms.feature.parquet").is_file());

    // All-in-one path, with a temporary quant directory.
    let aio = tmp.path().join("aio");
    std::fs::create_dir_all(&aio).unwrap();
    let pin = aio.join("a.pin");
    let aio_qpx = aio.join("bundle");
    let mut args = search_args();
    args.extend(
        [
            "--output-pin",
            s(&pin),
            "--output-parquet",
            s(&aio_qpx),
            "--lfq",
            "--mbr",
            "--rescore",
            "--percolator-bin",
            s(&perc),
            "--mbr-fdr",
            FIXTURE_MBR_FDR,
            "--quant-feature-fdr",
            "0.05",
        ]
        .map(String::from),
    );
    let out = andes(&args.iter().map(String::as_str).collect::<Vec<_>>());
    let log = String::from_utf8_lossy(&out.stderr);
    assert!(out.status.success(), "all-in-one failed:\n{log}");
    assert!(log.contains("[mbr] library:"), "{log}");
    assert!(log.contains("[mbr] extract mbr_run3:"), "{log}");

    assert!(
        same_mbr_tables(&staged.join("out"), &aio.join("a")),
        "mbr.parquet or mbr_validation.parquet differs between the stage path and the all-in-one run"
    );
    // The feature parquets hold the same column data.
    let staged_parquet = qpx.join("quantms.feature.parquet");
    assert!(
        parquet_column_data(&staged_parquet)
            == parquet_column_data(&aio_qpx.join("quantms.feature.parquet")),
        "quantms.feature.parquet differs between the stage path and the all-in-one run"
    );

    let (header, rows) = read_mbr(&aio.join("a.mbr.parquet"));
    let col = |name: &str| header.iter().position(|h| h == name).unwrap();
    let (ev, quantity) = (col("evidence_type"), col("quantity"));
    let row = |key: &str, run: &str| rows.get(&(key.to_string(), run.to_string()));

    // Peptide 4 was never triggered in run 3: it is transferred there, with
    // the designed quantity.
    let p4 = row("HLVDEPQNLIK/2", "mbr_run3").expect("peptide 4 in run 3");
    assert!(p4[ev].starts_with("transfer_ms"), "{p4:?}");
    let qv: f64 = p4[quantity].parse().unwrap();
    assert!(
        (qv / PEPTIDE4_QUANTITY - 1.0).abs() < 0.10,
        "peptide 4 quantity {qv} vs designed {PEPTIDE4_QUANTITY}"
    );
    for run in ["mbr_run1", "mbr_run2"] {
        assert_eq!(row("HLVDEPQNLIK/2", run).unwrap()[ev], "identified");
    }
    // Peptide 6 has no MS1 signal anywhere and is identified in run 1 only:
    // its identified rows fail the identified quality q-value.
    for key in ["DAFLGSFLYEYSR/2", "DAFLGSFLYEYSR/3"] {
        let r1 = row(key, "mbr_run1").unwrap_or_else(|| panic!("{key} in run 1"));
        assert_eq!(r1[ev], "identified_below_threshold", "{key}: {r1:?}");
        for run in ["mbr_run2", "mbr_run3"] {
            let r = row(key, run).unwrap_or_else(|| panic!("{key} in {run}"));
            assert!(!r[ev].starts_with("transfer_ms"), "{key} in {run}: {r:?}");
        }
    }
    // Most untriggered precursors with signal are transferred.
    let transfers = rows
        .values()
        .filter(|r| r[ev].starts_with("transfer_ms"))
        .count();
    eprintln!(
        "transfers at --mbr-fdr {FIXTURE_MBR_FDR}: {transfers} of {} rows",
        rows.len()
    );
    assert!(transfers >= 50, "only {transfers} transfers");

    // The feature parquet holds every mbr.parquet row except the below-threshold
    // transfers and identified rows.
    let below: Vec<&Vec<String>> = rows
        .values()
        .filter(|r| r[ev] == "below_threshold" || r[ev] == "identified_below_threshold")
        .collect();
    assert!(
        !below.is_empty(),
        "the fixture has below-threshold transfers"
    );
    let features = feature_rows(&staged_parquet);
    assert_eq!(features.len(), rows.len() - below.len());
    assert!(features
        .iter()
        .all(|f| f.3 != "below_threshold" && f.3 != "identified_below_threshold"));
    for r in below {
        let (pf, z) = r[0].rsplit_once('/').unwrap();
        assert!(
            !features
                .iter()
                .any(|f| f.0 == r[1] && f.1 == pf && f.2 == z),
            "below-threshold {r:?} is a feature row"
        );
    }
}

#[test]
fn duplicate_run_stems_are_rejected_before_the_search() {
    let tmp = tempfile::tempdir().unwrap();
    let src = fixture("test-fixtures/mbr_run1.mzML.gz");
    let (a, b) = (tmp.path().join("a"), tmp.path().join("b"));
    for d in [&a, &b] {
        std::fs::create_dir_all(d).unwrap();
        std::fs::copy(&src, d.join("mbr_run1.mzML.gz")).unwrap();
    }
    let (ra, rb) = (a.join("mbr_run1.mzML.gz"), b.join("mbr_run1.mzML.gz"));
    let db = fixture("test-fixtures/BSA.fasta");
    let out_dir = tmp.path().join("out");
    let pin = out_dir.join("x.pin");
    let q = out_dir.join("q");
    let common = [
        "--spectrum",
        s(&ra),
        "--spectrum",
        s(&rb),
        "--database",
        s(&db),
        "--output-pin",
        s(&pin),
    ];
    for extra in [
        &["--quant-dir", s(&q)][..],
        &["--lfq", "--mbr", "--rescore-native"][..],
    ] {
        let out = andes(&[&common[..], extra].concat());
        assert!(!out.status.success(), "{extra:?}");
        let err = String::from_utf8_lossy(&out.stderr);
        assert!(
            err.contains(s(&ra)) && err.contains(s(&rb)) && err.contains("mbr_run1"),
            "{extra:?}: {err}"
        );
        assert!(!err.contains("Loaded"), "searched before failing: {err}");
        assert!(!out_dir.exists(), "{extra:?} wrote {out_dir:?}");
    }
}

/// Splits a merged Percolator result of a multi-input search (SpecIds
/// `<file>/<id>`) into one file per run with the run's own SpecIds `<id>`.
fn split_by_run(merged: &Path, out_dir: &Path, kind: &str) -> Vec<String> {
    let text = std::fs::read_to_string(merged).unwrap();
    let mut lines = text.lines();
    let header = lines.next().unwrap();
    let mut per_run: HashMap<&str, String> =
        RUNS.iter().map(|r| (*r, format!("{header}\n"))).collect();
    for l in lines {
        // The prefix is the input's file name without its last extension; a
        // separate search of the same input writes the same SpecId.
        let (prefix, _) = l.split_once('/').expect("SpecId has a file prefix");
        let stem = prefix.split('.').next().unwrap();
        per_run
            .get_mut(stem)
            .unwrap_or_else(|| panic!("SpecId prefix {prefix:?}"))
            .push_str(&format!("{l}\n"));
    }
    RUNS.iter()
        .map(|r| {
            let path = out_dir.join(format!("{r}.{kind}.tsv"));
            std::fs::write(&path, &per_run[r]).unwrap();
            format!("{r}={}", path.display())
        })
        .collect()
}

#[test]
fn separate_searches_join_on_run_and_spec_id() {
    let tmp = tempfile::tempdir().unwrap();
    let perc = fake_percolator(tmp.path());

    // All-in-one over the three runs; its Percolator files give the scores.
    let aio = tmp.path().join("aio");
    std::fs::create_dir_all(&aio).unwrap();
    let mut args = search_args();
    args.extend(
        [
            "--output-pin",
            s(&aio.join("a.pin")),
            "--lfq",
            "--lfq-tol",
            "8ppm",
            "--mbr",
            "--rescore",
            "--percolator-bin",
            s(&perc),
        ]
        .map(String::from),
    );
    run_andes(&args.iter().map(String::as_str).collect::<Vec<_>>());

    // Each run searched on its own into one quant directory, with its own
    // Percolator result. The runs share scan titles; the file stem in front of
    // every SpecId keeps them apart.
    let sep = tmp.path().join("sep");
    let q = sep.join("q");
    std::fs::create_dir_all(&sep).unwrap();
    for r in RUNS {
        let spectrum = fixture(&format!("test-fixtures/{r}.mzML.gz"));
        let db = fixture("test-fixtures/BSA.fasta");
        run_andes(&[
            "--spectrum",
            s(&spectrum),
            "--database",
            s(&db),
            "--threads",
            "2",
            "--output-pin",
            s(&sep.join(format!("{r}.pin"))),
            "--quant-dir",
            s(&q),
            "--rescore",
            "--percolator-bin",
            s(&perc),
        ]);
    }
    let own = |kind: &str| -> Vec<String> {
        RUNS.iter()
            .map(|r| {
                format!(
                    "{r}={}",
                    sep.join(format!("{r}.percolator.{kind}.psms.txt"))
                        .display()
                )
            })
            .collect()
    };

    // One merged pair over every run joins: what one Percolator run over the
    // concatenated PINs of separate searches gives.
    let merged = |kind: &str| -> PathBuf {
        let path = sep.join(format!("merged.{kind}.psms.txt"));
        let mut text = String::new();
        for (i, r) in RUNS.iter().enumerate() {
            let part = std::fs::read_to_string(sep.join(format!("{r}.percolator.{kind}.psms.txt")))
                .unwrap();
            let mut lines = part.lines();
            let header = lines.next().unwrap();
            if i == 0 {
                text.push_str(header);
                text.push('\n');
            }
            for l in lines {
                text.push_str(l);
                text.push('\n');
            }
        }
        std::fs::write(&path, text).unwrap();
        path
    };
    let (merged_t, merged_d) = (merged("target"), merged("decoy"));
    run_andes(&[
        "quant",
        "library",
        "--quant-dir",
        s(&q),
        "--psms",
        s(&merged_t),
        "--decoy-psms",
        s(&merged_d),
    ]);
    assert!(q.join("library.parquet").is_file());

    // Swapped per-run files hold the other run's SpecIds, so nothing joins.
    let mut swapped: Vec<String> = ["quant", "library", "--quant-dir", s(&q), "--psms"]
        .map(String::from)
        .to_vec();
    let mut t = own("target");
    let a = t[0].split_once('=').unwrap().1.to_string();
    let b = t[1].split_once('=').unwrap().1.to_string();
    t[0] = format!("mbr_run1={b}");
    t[1] = format!("mbr_run2={a}");
    swapped.extend(t);
    swapped.push("--decoy-psms".into());
    swapped.extend(own("decoy"));
    let out = andes(&swapped.iter().map(String::as_str).collect::<Vec<_>>());
    assert!(!out.status.success());
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(
        err.contains("run mbr_run1")
            && err.contains("mbr_run2.percolator.target.psms.txt")
            && err.contains("has a SpecId of this run"),
        "{err}"
    );

    // Per-run results from the separate Percolator runs build a library.
    let mut lib_args: Vec<String> = ["quant", "library", "--quant-dir", s(&q), "--psms"]
        .map(String::from)
        .to_vec();
    lib_args.extend(own("target"));
    lib_args.push("--decoy-psms".into());
    lib_args.extend(own("decoy"));
    run_andes(&lib_args.iter().map(String::as_str).collect::<Vec<_>>());
    assert!(q.join("library.parquet").is_file());

    // The same per-PSM scores as the all-in-one run give its mbr.parquet exactly.
    let mut lib_args: Vec<String> = ["quant", "library", "--quant-dir", s(&q), "--psms"]
        .map(String::from)
        .to_vec();
    lib_args.extend(split_by_run(
        &aio.join("a.percolator.target.psms.txt"),
        &sep,
        "target",
    ));
    lib_args.push("--decoy-psms".into());
    lib_args.extend(split_by_run(
        &aio.join("a.percolator.decoy.psms.txt"),
        &sep,
        "decoy",
    ));
    run_andes(&lib_args.iter().map(String::as_str).collect::<Vec<_>>());
    for r in RUNS {
        run_andes(&[
            "quant",
            "extract",
            "--quant-dir",
            s(&q),
            "--run",
            r,
            "--lfq-tol",
            "8ppm",
        ]);
    }
    let base = sep.join("out");
    run_andes(&["quant", "finalize", "--quant-dir", s(&q), "-o", s(&base)]);
    assert!(
        read_mbr(&sep.join("out.mbr.parquet")).1 == read_mbr(&aio.join("a.mbr.parquet")).1,
        "mbr.parquet of the separate searches differs from the all-in-one run"
    );
    // The default `--quant-feature-fdr off` keeps identified rows whose MS1
    // feature fails its decoys (peptide 6 has no MS1 signal) as `identified`.
    // With no signal its quantity is missing (NA), not 0.
    let (header, rows) = read_mbr(&sep.join("out.mbr.parquet"));
    let ev = header.iter().position(|h| h == "evidence_type").unwrap();
    let qty = header.iter().position(|h| h == "quantity").unwrap();
    for key in ["DAFLGSFLYEYSR/2", "DAFLGSFLYEYSR/3"] {
        let r1 = rows
            .get(&(key.to_string(), "mbr_run1".to_string()))
            .unwrap_or_else(|| panic!("{key} in run 1"));
        assert_eq!(r1[ev], "identified", "{key}: {r1:?}");
        assert_ne!(r1[qty], "0", "{key}: {r1:?}");
    }
    assert!(rows
        .values()
        .all(|r| r[qty] == "NA" || r[qty].parse::<f64>().unwrap() > 0.0));
    assert!(rows.values().all(|r| r[ev] != "identified_below_threshold"));
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
    assert!(err.contains("absent_run.ms1.parquet"), "{err}");

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

    // Finalize checks every extract against the library and the tolerance.
    let library = |fdr: &str| {
        run_andes(&[
            "quant",
            "library",
            "--quant-dir",
            s(&q),
            "--psms",
            s(&targets),
            "--decoy-psms",
            s(&decoys),
            "--mbr-library-fdr",
            fdr,
        ])
    };
    let extract = |run: &str, tol: &str| {
        run_andes(&[
            "quant",
            "extract",
            "--quant-dir",
            s(&q),
            "--run",
            run,
            "--lfq-tol",
            tol,
        ])
    };
    let base = tmp.path().join("out");
    let finalize = || andes(&["quant", "finalize", "--quant-dir", s(&q), "-o", s(&base)]);
    library("0.01");
    for r in RUNS {
        extract(r, "10ppm");
    }
    extract("mbr_run2", "5ppm");
    let out = finalize();
    assert!(!out.status.success());
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(
        err.contains("mbr_run2.extract.parquet") && err.contains("--lfq-tol"),
        "{err}"
    );
    extract("mbr_run2", "10ppm");
    assert!(finalize().status.success());
    library("0.5");
    let out = finalize();
    assert!(!out.status.success());
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(
        err.contains("mbr_run1.extract.parquet") && err.contains("stale"),
        "{err}"
    );
}

/// `mbr.parquet` rows whose evidence type is `identified`.
fn identified_rows(path: &Path) -> usize {
    let (header, rows) = read_mbr(path);
    let ev = header.iter().position(|h| h == "evidence_type").unwrap();
    rows.values().filter(|r| r[ev] == "identified").count()
}

#[test]
fn mbr_donors_honour_the_pep_cap() {
    // Every PSM at q = 0.001 with PEP = 0.5: `--pep 0.01` leaves no donor.
    let tmp = tempfile::tempdir().unwrap();
    let perc = fake_percolator_with_pep(tmp.path(), "0.5");

    // All-in-one.
    let aio = tmp.path().join("aio");
    std::fs::create_dir_all(&aio).unwrap();
    let mut args = search_args();
    args.extend(
        [
            "--output-pin",
            s(&aio.join("a.pin")),
            "--lfq",
            "--mbr",
            "--rescore",
            "--percolator-bin",
            s(&perc),
            "--pep",
            "0.01",
        ]
        .map(String::from),
    );
    let out = andes(&args.iter().map(String::as_str).collect::<Vec<_>>());
    let log = String::from_utf8_lossy(&out.stderr);
    assert!(out.status.success(), "all-in-one failed:\n{log}");
    let n = identified_rows(&aio.join("a.mbr.parquet"));
    eprintln!("all-in-one --pep 0.01: {n} identified rows");
    assert_eq!(n, 0, "identified rows from PSMs above the PEP cap");

    // Stages: `quant library --pep`.
    let staged = tmp.path().join("staged");
    let q = staged_search(&staged, &perc);
    let library = |pep: Option<&str>| {
        let mut a = vec!["quant", "library", "--quant-dir", s(&q)];
        let (t, d) = (
            staged.join("s.percolator.target.psms.txt"),
            staged.join("s.percolator.decoy.psms.txt"),
        );
        let (t, d) = (t.display().to_string(), d.display().to_string());
        let mut owned = vec!["--psms".to_string(), t, "--decoy-psms".to_string(), d];
        if let Some(p) = pep {
            owned.extend(["--pep".to_string(), p.to_string()]);
        }
        a.extend(owned.iter().map(String::as_str));
        run_andes(&a);
        for r in RUNS {
            run_andes(&["quant", "extract", "--quant-dir", s(&q), "--run", r]);
        }
        let base = staged.join("out");
        run_andes(&["quant", "finalize", "--quant-dir", s(&q), "-o", s(&base)]);
        identified_rows(&staged.join("out.mbr.parquet"))
    };
    let uncapped = library(None);
    let capped = library(Some("0.01"));
    eprintln!("stages: {uncapped} identified rows without --pep, {capped} with --pep 0.01");
    assert!(uncapped > 0);
    assert_eq!(capped, 0, "identified rows from PSMs above the PEP cap");
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
