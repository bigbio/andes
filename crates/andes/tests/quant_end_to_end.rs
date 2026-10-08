//! End-to-end tests for `--tmt` and `--lfq` on the synthetic fixtures built by
//! `scripts/make_quant_fixtures.py` (six BSA peptides with known intensities).
//!
//! The assertions hold the designed numbers of that script: reporter channels
//! 6..10 at twice channels 1..5 in every MS2, a 1..10 ladder in the MS3, a
//! precursor purity of 2/3 for the second peptide, and label-free areas
//! proportional to the designed apex intensities.

use std::collections::HashMap;
use std::path::PathBuf;
use std::process::Command;

use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

fn fixture(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(rel)
        .canonicalize()
        .unwrap_or_else(|e| panic!("canonicalize {rel}: {e}"))
}

fn read_tsv(path: &PathBuf) -> (Vec<String>, Vec<Vec<String>>) {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("read {path:?}: {e}"));
    let mut lines = text.lines();
    let header: Vec<String> = lines
        .next()
        .expect("header")
        .split('\t')
        .map(String::from)
        .collect();
    let rows = lines
        .map(|l| l.split('\t').map(String::from).collect::<Vec<_>>())
        .collect();
    (header, rows)
}

/// Peptidoform without its `[...]` modification tags.
fn bare(peptidoform: &str) -> String {
    let mut out = String::new();
    let mut depth = 0usize;
    for c in peptidoform.chars() {
        match c {
            '[' => depth += 1,
            ']' => depth = depth.saturating_sub(1),
            _ if depth == 0 => out.push(c),
            _ => {}
        }
    }
    out
}

fn col(header: &[String], name: &str) -> usize {
    header
        .iter()
        .position(|h| h == name)
        .unwrap_or_else(|| panic!("column {name} missing in {header:?}"))
}

fn run_andes(args: &[&str]) {
    let status = Command::new(env!("CARGO_BIN_EXE_andes"))
        .args(args)
        .status()
        .expect("run andes");
    assert!(status.success(), "andes failed: {status}");
}

fn parquet_rows(path: &PathBuf) -> (usize, Vec<String>) {
    let file = std::fs::File::open(path).unwrap_or_else(|e| panic!("open {path:?}: {e}"));
    let builder = ParquetRecordBatchReaderBuilder::try_new(file).expect("parquet");
    let names: Vec<String> = builder
        .schema()
        .fields()
        .iter()
        .map(|f| f.name().clone())
        .collect();
    assert_eq!(
        builder
            .schema()
            .metadata()
            .get("file_type")
            .map(String::as_str),
        Some("feature_file")
    );
    let mut n = 0;
    for b in builder.build().unwrap() {
        n += b.unwrap().num_rows();
    }
    (n, names)
}

#[test]
fn tmt10_ms2_reporters_purity_and_correction() {
    let dir = tempfile::tempdir().unwrap();
    let pin = dir.path().join("tmt.pin");
    let qpx = dir.path().join("tmt.idparquet");
    run_andes(&[
        "--spectrum",
        fixture("test-fixtures/tmt10_synthetic.mzML.gz")
            .to_str()
            .unwrap(),
        "--database",
        fixture("test-fixtures/BSA.fasta").to_str().unwrap(),
        "--mods",
        fixture("test-fixtures/tmt10_mods.txt").to_str().unwrap(),
        "--output-pin",
        pin.to_str().unwrap(),
        "--output-parquet",
        qpx.to_str().unwrap(),
        "--tmt",
        "tmt10",
        "--tmt-correction",
        fixture("test-fixtures/tmt10_correction.txt")
            .to_str()
            .unwrap(),
        "--threads",
        "2",
    ]);
    let (header, rows) = read_tsv(&dir.path().join("tmt.tmt.tsv"));
    assert!(
        rows.len() >= 4,
        "expected most of the 6 peptides quantified, got {}",
        rows.len()
    );
    let c126 = col(&header, "tmt_126");
    let c129c = col(&header, "tmt_129C");
    let raw126 = col(&header, "raw_tmt_126");
    let raw129c = col(&header, "raw_tmt_129C");
    let purity = col(&header, "purity");
    let peptide = col(&header, "peptide");
    let level = col(&header, "quant_ms_level");
    let mut saw_interfered = false;
    for r in &rows {
        assert_eq!(r[level], "2");
        let raw_a: f64 = r[raw126].parse().unwrap();
        let raw_b: f64 = r[raw129c].parse().unwrap();
        // Designed pattern: channels 6..10 at twice channels 1..5.
        assert!(
            (raw_b / raw_a - 2.0).abs() < 0.02,
            "{}: {raw_a} vs {raw_b}",
            r[peptide]
        );
        // The correction undoes the 5 % +1 spill the fixture does NOT contain:
        // corrected 126 (which receives nothing) is raw/0.95, while 129C is
        // inflated less than that because a share of its signal is attributed
        // to the 128C spill.
        let cor_a: f64 = r[c126].parse().unwrap();
        let cor_b: f64 = r[c129c].parse().unwrap();
        assert!(
            (cor_a / raw_a - 1.0 / 0.95).abs() < 0.01,
            "{}: corrected 126 {cor_a} raw {raw_a}",
            r[peptide]
        );
        assert!(
            cor_b > raw_b && cor_b / raw_b < 1.0 / 0.95 - 0.005,
            "{}: corrected 129C {cor_b} raw {raw_b}",
            r[peptide]
        );
        let p: f64 = r[purity].parse().unwrap();
        if bare(&r[peptide]) == "YLYEIAR" {
            // Envelope inside the ±0.7 Th window (M, M+1) ≈ 1.5 × mono against an
            // interferer of 0.5 × mono → purity ≈ 0.75 (plus a little background).
            assert!((0.70..0.80).contains(&p), "interfered purity {p}");
            saw_interfered = true;
        } else {
            assert!(p > 0.97, "{}: purity {p}", r[peptide]);
        }
    }
    assert!(
        saw_interfered,
        "YLYEIAR (the interfered peptide) should be quantified"
    );
    let (n, names) = parquet_rows(&qpx.join("quantms.feature.parquet"));
    assert_eq!(n, rows.len());
    assert!(names.contains(&"intensities".to_string()));
}

#[test]
fn tmt10_ms3_reporters_follow_the_sps_scan() {
    let dir = tempfile::tempdir().unwrap();
    let pin = dir.path().join("tmt3.pin");
    run_andes(&[
        "--spectrum",
        fixture("test-fixtures/tmt10_synthetic.mzML.gz")
            .to_str()
            .unwrap(),
        "--database",
        fixture("test-fixtures/BSA.fasta").to_str().unwrap(),
        "--mods",
        fixture("test-fixtures/tmt10_mods.txt").to_str().unwrap(),
        "--output-pin",
        pin.to_str().unwrap(),
        "--tmt",
        "tmt10",
        "--tmt-level",
        "3",
        "--threads",
        "2",
    ]);
    let (header, rows) = read_tsv(&dir.path().join("tmt3.tmt.tsv"));
    assert!(rows.len() >= 4, "got {} rows", rows.len());
    let c126 = col(&header, "tmt_126");
    let c131 = col(&header, "tmt_131");
    let level = col(&header, "quant_ms_level");
    let quant_scan = col(&header, "quant_scan");
    let spec_id = col(&header, "spec_id");
    assert!(
        !header.iter().any(|h| h.starts_with("raw_")),
        "no correction → no raw block"
    );
    for r in &rows {
        assert_eq!(r[level], "3");
        assert_ne!(
            r[quant_scan], r[spec_id],
            "reporters come from the MS3, not the MS2"
        );
        let a: f64 = r[c126].parse().unwrap();
        let b: f64 = r[c131].parse().unwrap();
        // MS3 ladder: channel 10 is ten times channel 1.
        assert!((b / a - 10.0).abs() < 0.05, "{a} vs {b}");
    }
}

#[test]
fn lfq_areas_track_the_designed_intensities() {
    let dir = tempfile::tempdir().unwrap();
    let pin = dir.path().join("lfq.pin");
    let qpx = dir.path().join("lfq.idparquet");
    run_andes(&[
        "--spectrum",
        fixture("test-fixtures/lfq_synthetic.mzML.gz")
            .to_str()
            .unwrap(),
        "--database",
        fixture("test-fixtures/BSA.fasta").to_str().unwrap(),
        "--output-pin",
        pin.to_str().unwrap(),
        "--output-parquet",
        qpx.to_str().unwrap(),
        "--lfq",
        // Six targets: the conservative (decoys + 1) / targets estimate cannot
        // go below 1/6, so the wide table needs a looser cut on this fixture.
        "--lfq-feature-fdr",
        "0.2",
        "--threads",
        "2",
    ]);
    let (header, rows) = read_tsv(&dir.path().join("lfq.lfq_features.tsv"));
    assert!(
        rows.len() >= 4,
        "expected most of the 6 peptides, got {}",
        rows.len()
    );
    let peptide = col(&header, "peptide");
    let area = col(&header, "area");
    let cosine = col(&header, "cosine");
    let apex = col(&header, "rt_apex");
    let q = col(&header, "feature_q_value");
    let n_iso = col(&header, "n_isotopes");
    // Designed apex intensities (monoisotopic) and apex RTs.
    let design: HashMap<&str, (f64, f64)> = HashMap::from([
        ("LVNELTEFAK", (1.0e6, 40.0)),
        ("YLYEIAR", (5.0e5, 70.0)),
        ("AEFVEVTK", (2.0e6, 100.0)),
        ("HLVDEPQNLIK", (8.0e5, 130.0)),
        ("LGEYGFQNALIVR", (3.0e5, 160.0)),
        ("DAFLGSFLYEYSR", (1.5e6, 190.0)),
    ]);
    let iso_areas = col(&header, "isotope_areas");
    // The mono trace is a Gaussian of the designed height and sigma 4 s, so its
    // area is height · sigma · sqrt(2π) ≈ 10.03 · height for every peptide.
    let gaussian_area_per_height = 4.0 * (2.0 * std::f64::consts::PI).sqrt();
    for r in &rows {
        let seq_owned = bare(&r[peptide]);
        let seq = seq_owned.as_str();
        let (apex_int, apex_rt) = design[seq];
        let a: f64 = r[area].parse().unwrap();
        let mono_area: f64 = r[iso_areas].split(';').next().unwrap().parse().unwrap();
        let rt: f64 = r[apex].parse().unwrap();
        let cos: f64 = r[cosine].parse().unwrap();
        assert!((rt - apex_rt).abs() <= 1.5, "{seq}: apex {rt} vs {apex_rt}");
        assert!(cos > 0.95, "{seq}: cosine {cos}");
        assert!(
            r[n_iso].parse::<u32>().unwrap() >= 3,
            "{seq}: isotopes {}",
            r[n_iso]
        );
        let qv: f64 = r[q].parse().unwrap_or(1.0);
        assert!(qv <= 0.2, "{seq}: feature q {qv}");
        assert!(
            (mono_area / apex_int / gaussian_area_per_height - 1.0).abs() < 0.05,
            "{seq}: mono area {mono_area} for height {apex_int}"
        );
        assert!(
            a > mono_area,
            "{seq}: summed area {a} must exceed the mono area {mono_area}"
        );
    }
    // Wide table: one column for the run, every confident row has a value.
    let (wide_header, wide_rows) = read_tsv(&dir.path().join("lfq.lfq.tsv"));
    assert_eq!(
        wide_header.last().map(String::as_str),
        Some("lfq_synthetic")
    );
    assert_eq!(wide_rows.len(), rows.len());
    assert!(wide_rows.iter().all(|r| !r.last().unwrap().is_empty()));
    let (n, _) = parquet_rows(&qpx.join("quantms.feature.parquet"));
    assert_eq!(n, rows.len());
}

#[test]
fn lfq_rejects_mgf_input() {
    let dir = tempfile::tempdir().unwrap();
    let pin = dir.path().join("x.pin");
    let status = Command::new(env!("CARGO_BIN_EXE_andes"))
        .args([
            "--spectrum",
            fixture("test-fixtures/test.mgf.gz").to_str().unwrap(),
            "--database",
            fixture("test-fixtures/BSA.fasta").to_str().unwrap(),
            "--output-pin",
            pin.to_str().unwrap(),
            "--lfq",
        ])
        .status()
        .expect("run andes");
    assert!(!status.success(), "--lfq on MGF must fail");
}
