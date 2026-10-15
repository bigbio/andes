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
fn tmt35_reads_the_deuterated_kit_in_openms_order() {
    // The TMT10 fixture's reporters sit on the 35-plex channels of the same
    // name (TMT10's 131 is 131N); the deuterated channels, 2.9 mDa above
    // their twins, read nothing at the 10 ppm default.
    let dir = tempfile::tempdir().unwrap();
    let pin = dir.path().join("tmt35.pin");
    let qpx = dir.path().join("tmt35.idparquet");
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
        "tmt35",
        "--threads",
        "2",
    ]);
    let (header, rows) = read_tsv(&dir.path().join("tmt35.tmt.tsv"));
    assert!(rows.len() >= 4, "got {} rows", rows.len());
    let first = col(&header, "tmt_126");
    let channels: Vec<&str> = header[first..].iter().map(String::as_str).collect();
    assert_eq!(channels.len(), 35, "{channels:?}");
    assert_eq!(
        &channels[..8],
        &[
            "tmt_126",
            "tmt_127N",
            "tmt_127C",
            "tmt_127D",
            "tmt_128N",
            "tmt_128C",
            "tmt_128ND",
            "tmt_128CD"
        ]
    );
    assert_eq!(channels[34], "tmt_135CD");
    let c126 = col(&header, "tmt_126");
    let c129c = col(&header, "tmt_129C");
    let c131n = col(&header, "tmt_131N");
    let deuterated: Vec<usize> = header
        .iter()
        .enumerate()
        .filter(|(_, h)| h.starts_with("tmt_") && h.ends_with('D'))
        .map(|(i, _)| i)
        .collect();
    assert_eq!(deuterated.len(), 17);
    for r in &rows {
        let a: f64 = r[c126].parse().unwrap();
        let b: f64 = r[c129c].parse().unwrap();
        assert!(a > 0.0 && (b / a - 2.0).abs() < 0.02, "{a} vs {b}");
        assert!(r[c131n].parse::<f64>().unwrap() > 0.0);
        for &d in &deuterated {
            assert_eq!(r[d].parse::<f64>().unwrap(), 0.0, "{}", header[d]);
        }
    }
    let (n, _) = parquet_rows(&qpx.join("quantms.feature.parquet"));
    assert_eq!(n, rows.len());
}

#[test]
fn tmt35_refuses_a_correction_sheet_and_a_20ppm_window() {
    let dir = tempfile::tempdir().unwrap();
    let pin = dir.path().join("tmt35.pin");
    let refused = |extra: &[&str]| {
        let spectrum = fixture("test-fixtures/tmt10_synthetic.mzML.gz");
        let database = fixture("test-fixtures/BSA.fasta");
        let mods = fixture("test-fixtures/tmt10_mods.txt");
        let mut args = vec![
            "--spectrum",
            spectrum.to_str().unwrap(),
            "--database",
            database.to_str().unwrap(),
            "--mods",
            mods.to_str().unwrap(),
            "--output-pin",
            pin.to_str().unwrap(),
            "--tmt",
            "tmt35",
            "--threads",
            "2",
        ];
        args.extend_from_slice(extra);
        let out = Command::new(env!("CARGO_BIN_EXE_andes"))
            .args(&args)
            .output()
            .expect("run andes");
        assert!(!out.status.success(), "andes accepted {extra:?}");
        String::from_utf8_lossy(&out.stderr).into_owned()
    };
    let sheet = dir.path().join("lot.txt");
    std::fs::write(&sheet, "NA/NA/NA/NA/NA/NA/NA/NA\n".repeat(35)).unwrap();
    let stderr = refused(&["--tmt-correction", sheet.to_str().unwrap()]);
    assert!(stderr.contains("not supported for tmt35plex"), "{stderr}");
    let stderr = refused(&["--tmt-tol", "20ppm"]);
    assert!(stderr.contains("e.g. 10ppm"), "{stderr}");
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

/// The LFQ fixture search; `extra` is appended.
fn run_lfq_fixture(pin: &std::path::Path, qpx: &std::path::Path, extra: &[&str]) {
    let spectrum = fixture("test-fixtures/lfq_synthetic.mzML.gz");
    let database = fixture("test-fixtures/BSA.fasta");
    let mut args = vec![
        "--spectrum",
        spectrum.to_str().unwrap(),
        "--database",
        database.to_str().unwrap(),
        "--output-pin",
        pin.to_str().unwrap(),
        "--output-parquet",
        qpx.to_str().unwrap(),
        "--lfq",
        "--lfq-feature-fdr",
        "0.2",
        "--threads",
        "2",
    ];
    args.extend_from_slice(extra);
    run_andes(&args);
}

/// SpecIds of a PIN file.
fn pin_spec_ids(pin: &std::path::Path) -> std::collections::HashSet<String> {
    std::fs::read_to_string(pin)
        .unwrap()
        .lines()
        .skip(1)
        .filter_map(|l| l.split('\t').next().map(String::from))
        .collect()
}

#[test]
fn quant_dir_writes_ms1_and_psm_tables() {
    let out = tempfile::tempdir().unwrap();
    let qdir = out.path().join("q");
    let pin = out.path().join("lfq.pin");
    let qpx = out.path().join("lfq.idparquet");
    run_lfq_fixture(&pin, &qpx, &["--quant-dir", qdir.to_str().unwrap()]);
    let stem = "lfq_synthetic";
    let ms1 = output::ms1_table::read_ms1_table(&qdir.join(format!("{stem}.ms1.parquet"))).unwrap();
    assert!(!ms1.is_empty());
    assert!(!qdir.join(format!("{stem}.ms1.parquet.partial")).exists());
    let (run, rows) =
        output::run_quant_table::read_run_quant_table(&qdir.join(format!("{stem}.quant.parquet")))
            .unwrap();
    assert_eq!(run, stem);
    assert!(
        rows.iter().filter(|r| !r.is_decoy).count() >= 6,
        "the six designed BSA peptides"
    );
    let pin_ids = pin_spec_ids(&pin);
    for r in &rows {
        assert!(pin_ids.contains(&r.spec_id), "{} not in PIN", r.spec_id);
        assert!((r.envelope.iter().sum::<f32>() - 1.0).abs() < 1e-3);
        assert_eq!(r.precursor_key, format!("{}/{}", r.peptidoform, r.charge));
        assert!(r.rt_s > 0.0);
    }
}

/// The real-data slice has rank-1 decoys: `--quant-dir` alone (no `--lfq`)
/// writes them, unfiltered, with PIN SpecIds.
#[test]
fn quant_dir_table_contains_decoys() {
    let out = tempfile::tempdir().unwrap();
    let qdir = out.path().join("q");
    let pin = out.path().join("slice.pin");
    run_andes(&[
        "--spectrum",
        fixture("test-fixtures/astral_pxd070049_slice.mzML.gz")
            .to_str()
            .unwrap(),
        "--database",
        fixture("test-fixtures/ecoli.fasta").to_str().unwrap(),
        "--output-pin",
        pin.to_str().unwrap(),
        "--quant-dir",
        qdir.to_str().unwrap(),
        "--output-parquet",
        out.path().join("slice.idparquet").to_str().unwrap(),
        "--threads",
        "2",
    ]);
    // --quant-dir alone quantifies nothing: no feature parquet.
    assert!(out.path().join("slice.idparquet").is_dir());
    assert!(!out
        .path()
        .join("slice.idparquet/quantms.feature.parquet")
        .exists());
    let stem = "astral_pxd070049_slice";
    assert!(qdir.join(format!("{stem}.ms1.parquet")).exists());
    let (run, rows) =
        output::run_quant_table::read_run_quant_table(&qdir.join(format!("{stem}.quant.parquet")))
            .unwrap();
    assert_eq!(run, stem);
    assert!(
        rows.iter().any(|r| r.is_decoy),
        "decoy PSMs must be written"
    );
    assert!(rows.iter().any(|r| !r.is_decoy));
    let pin_ids = pin_spec_ids(&pin);
    for r in &rows {
        assert!(pin_ids.contains(&r.spec_id), "{} not in PIN", r.spec_id);
        assert!((r.envelope.iter().sum::<f32>() - 1.0).abs() < 1e-3);
    }

    // The MS2 table: every PSM's scan, reduced, with its isolation window;
    // the best target PSM's scan gives a library spectrum (the slice is
    // searched against one of its three proteomes, so most lower-ranked
    // targets are random matches with few fragments).
    let ms2 = output::ms2_table::read_ms2_table(&qdir.join(format!("{stem}.ms2.parquet"))).unwrap();
    assert!(!qdir.join(format!("{stem}.ms2.parquet.partial")).exists());
    assert_eq!(ms2.run, stem);
    assert!(ms2.scans.windows(2).all(|w| w[0].rt <= w[1].rt));
    assert!(ms2
        .scans
        .iter()
        .all(|s| s.peaks.len() <= ms2.top_peaks && s.peaks.windows(2).all(|w| w[0].0 <= w[1].0)));
    assert!(
        ms2.scans.iter().all(|s| s.isolation.is_some()),
        "isolation windows recorded"
    );
    let by_key: std::collections::HashMap<&str, &quant::ms2::Ms2Scan> =
        ms2.scans.iter().map(|s| (s.spec_key.as_str(), s)).collect();
    let mods =
        output::run_quant_table::read_run_quant_mods(&qdir.join(format!("{stem}.quant.parquet")))
            .unwrap();
    let mut targets: Vec<&output::run_quant_table::RunPsmRow> =
        rows.iter().filter(|r| !r.is_decoy).collect();
    targets.sort_by(|a, b| b.search_score.total_cmp(&a.search_score));
    let best = &targets[..targets.len().min(20)];
    let mut with_spectrum = 0;
    for r in best {
        // `<spec_key>_<rank>` or `<spec_key>_<rank>_<row>`.
        let one = r.spec_id.rsplit_once('_').unwrap().0;
        let scan = by_key
            .get(one)
            .or_else(|| by_key.get(one.rsplit_once('_')?.0))
            .unwrap_or_else(|| panic!("{} has no scan", r.spec_id));
        assert!((scan.rt - r.rt_s).abs() < 1e-9);
        let (lo, hi) = scan.window();
        assert!(
            (lo - 1.5..hi + 1.5).contains(&r.mz),
            "{} outside {lo}-{hi}",
            r.mz
        );
        let masses = quant::ms2::residue_masses(&r.peptidoform, &mods).unwrap();
        let ladder = quant::ms2::fragment_ladder(&masses, r.charge);
        with_spectrum += usize::from(
            quant::ms2::library_spectrum(&ladder, &scan.peaks, ms2.fragment_tol).is_some(),
        );
    }
    assert!(with_spectrum >= 1, "{with_spectrum} of {}", best.len());
    let first = best[0];
    let one = first.spec_id.rsplit_once('_').unwrap().0;
    let scan = by_key
        .get(one)
        .or_else(|| by_key.get(one.rsplit_once('_')?.0))
        .unwrap();
    let masses = quant::ms2::residue_masses(&first.peptidoform, &mods).unwrap();
    let ladder = quant::ms2::fragment_ladder(&masses, first.charge);
    let lib = quant::ms2::library_spectrum(&ladder, &scan.peaks, ms2.fragment_tol)
        .expect("the best target PSM's scan gives a library spectrum");
    assert!(quant::ms2::score_scan(&lib, &scan.peaks, ms2.fragment_tol).angle > 0.99);
}

#[test]
fn outputs_unchanged_without_quant_dir() {
    let plain = tempfile::tempdir().unwrap();
    let with = tempfile::tempdir().unwrap();
    for (dir, quant_dir) in [(&plain, false), (&with, true)] {
        let pin = dir.path().join("lfq.pin");
        let qpx = dir.path().join("lfq.idparquet");
        let qd = dir.path().join("q");
        if quant_dir {
            run_lfq_fixture(&pin, &qpx, &["--quant-dir", qd.to_str().unwrap()]);
        } else {
            run_lfq_fixture(&pin, &qpx, &[]);
        }
    }
    for name in ["lfq.lfq.tsv", "lfq.lfq_features.tsv", "lfq.pin"] {
        let a = std::fs::read(plain.path().join(name)).unwrap();
        let b = std::fs::read(with.path().join(name)).unwrap();
        assert!(!a.is_empty(), "{name} is empty");
        assert_eq!(a, b, "{name} differs with --quant-dir");
    }
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
