//! End-to-end `--label`: two synthetic SILAC duplex runs
//! (`scripts/make_quant_fixtures.py --silac`), searched with the heavy
//! Lys8/Arg10 labels, rescored natively and quantified in both channels with
//! match between channels and match between runs.
//!
//! The fixture's peptides have H/L ratios of 1, 0.5, 2 or 0.25 (of their
//! monoisotopes; the molar ratio andes reports is 0.10-0.33 log2 lower,
//! since a heavy peptide's natural envelope is narrower), every tenth has no
//! heavy channel, and each channel is triggered for MS2 with probability 0.6,
//! so many groups are identified in one channel or one run only. The labels
//! are fully enriched, so the test sets `--label-heavy-enrichment 1`.

use std::collections::HashMap;
use std::fs::File;
use std::path::{Path, PathBuf};
use std::process::Command;

use arrow::array::{Array, BooleanArray, Float64Array, StringArray, UInt8Array};
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

fn fixture(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(rel)
        .canonicalize()
        .unwrap_or_else(|e| panic!("canonicalize {rel}: {e}"))
}

struct Row {
    run: String,
    peptide: String,
    charge: u8,
    channel: u8,
    status: String,
    transfer: bool,
    log2_ratio: Option<f64>,
}

fn read_rows(path: &Path) -> Vec<Row> {
    let file = File::open(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let reader = ParquetRecordBatchReaderBuilder::try_new(file)
        .unwrap()
        .build()
        .unwrap();
    let mut out = Vec::new();
    for batch in reader {
        let b = batch.unwrap();
        let col = |n: &str| b.column(b.schema().index_of(n).unwrap()).clone();
        let s = |n: &str| {
            col(n)
                .as_any()
                .downcast_ref::<StringArray>()
                .unwrap()
                .clone()
        };
        let run = s("filename");
        let peptide = s("peptide");
        let status = s("status");
        let charge = col("charge");
        let charge = charge.as_any().downcast_ref::<UInt8Array>().unwrap();
        let channel = col("channel");
        let channel = channel.as_any().downcast_ref::<UInt8Array>().unwrap();
        let transfer = col("transfer");
        let transfer = transfer.as_any().downcast_ref::<BooleanArray>().unwrap();
        let ratio = col("log2_ratio");
        let ratio = ratio.as_any().downcast_ref::<Float64Array>().unwrap();
        for i in 0..b.num_rows() {
            out.push(Row {
                run: run.value(i).to_string(),
                peptide: peptide.value(i).to_string(),
                charge: charge.value(i),
                channel: channel.value(i),
                status: status.value(i).to_string(),
                transfer: transfer.value(i),
                log2_ratio: (!ratio.is_null(i)).then(|| ratio.value(i)),
            });
        }
    }
    out
}

#[test]
fn silac_duplex_quantifies_matches_and_transfers() {
    let tmp = tempfile::tempdir().unwrap();
    let base = tmp.path().join("s");
    let pin = tmp.path().join("s.pin");
    let qpx = tmp.path().join("qpx");
    let out = Command::new(env!("CARGO_BIN_EXE_andes"))
        .args([
            "--spectrum",
            fixture("test-fixtures/silac_run1.mzML.gz")
                .to_str()
                .unwrap(),
            "--spectrum",
            fixture("test-fixtures/silac_run2.mzML.gz")
                .to_str()
                .unwrap(),
            "--database",
            fixture("test-fixtures/BSA.fasta").to_str().unwrap(),
            "--threads",
            "2",
            "--output-pin",
            pin.to_str().unwrap(),
            "--output-parquet",
            qpx.to_str().unwrap(),
            "--label",
            "silac2",
            "--label-heavy-enrichment",
            "1",
            "--mbr",
            "--rescore-native",
        ])
        .output()
        .expect("run andes");
    let log = String::from_utf8_lossy(&out.stderr).to_string();
    assert!(out.status.success(), "andes failed:\n{log}");
    let _ = base;
    let table = tmp.path().join("s.label.parquet");
    let rows = read_rows(&table);
    assert!(!rows.is_empty(), "{log}");

    let count = |ch: u8, st: &str| {
        rows.iter()
            .filter(|r| r.channel == ch && r.status == st)
            .count()
    };
    // Match between channels on both sides, transfers between runs, and the
    // heavy-absent peptides below detection.
    for ch in [1, 2] {
        assert!(count(ch, "identified") >= 150, "channel {ch}:\n{log}");
        assert!(count(ch, "match_channel") >= 50, "channel {ch}:\n{log}");
        assert!(count(ch, "transfer") >= 20, "channel {ch}:\n{log}");
    }
    assert!(count(2, "below_detection") >= 15, "{log}");
    assert_eq!(count(1, "below_detection"), 0, "{log}");
    // A heavy channel below detection has its light channel measured.
    let status_of: HashMap<(String, String, u8, u8), &str> = rows
        .iter()
        .map(|r| {
            (
                (r.run.clone(), r.peptide.clone(), r.charge, r.channel),
                r.status.as_str(),
            )
        })
        .collect();
    for r in rows
        .iter()
        .filter(|r| r.channel == 2 && r.status == "below_detection")
    {
        let light = status_of[&(r.run.clone(), r.peptide.clone(), r.charge, 1)];
        assert!(
            ["identified", "match_channel", "transfer"].contains(&light),
            "{} light {light}",
            r.peptide
        );
    }
    // Every ratio sits near a design level: the molar ratio is 0.10 to 0.33
    // log2 below the monoisotope design for these peptides (a heavy
    // peptide's envelope is narrower, so its monoisotope holds more of it).
    let ratios: Vec<&Row> = rows
        .iter()
        .filter(|r| r.channel == 2 && r.log2_ratio.is_some())
        .collect();
    assert!(ratios.len() >= 200, "{}", ratios.len());
    let off: Vec<f64> = ratios
        .iter()
        .map(|r| {
            let v = r.log2_ratio.unwrap();
            [0.0, -1.0, 1.0, -2.0]
                .iter()
                .map(|l| v - l)
                .min_by(|a, b| a.abs().total_cmp(&b.abs()))
                .unwrap()
        })
        .collect();
    let within = off.iter().filter(|d| (-0.36..=-0.07).contains(*d)).count();
    assert!(
        within as f64 >= 0.97 * off.len() as f64,
        "{within} of {} ratios near a design level",
        off.len()
    );
    // The same group gives the same ratio in both runs, whether identified,
    // matched between channels or transferred.
    let mut by_group: HashMap<(String, u8), Vec<(f64, bool)>> = HashMap::new();
    for r in &ratios {
        by_group
            .entry((r.peptide.clone(), r.charge))
            .or_default()
            .push((r.log2_ratio.unwrap(), r.transfer));
    }
    let mut pairs = 0;
    let mut with_transfer = 0;
    for v in by_group.values().filter(|v| v.len() == 2) {
        pairs += 1;
        with_transfer += usize::from(v.iter().any(|x| x.1));
        assert!((v[0].0 - v[1].0).abs() < 0.05, "{v:?}");
    }
    assert!(
        pairs >= 80 && with_transfer >= 15,
        "{pairs} {with_transfer}"
    );

    // The QPX features carry the SILAC channel labels.
    let feature = qpx.join("quantms.feature.parquet");
    let file = File::open(&feature).unwrap();
    let reader = ParquetRecordBatchReaderBuilder::try_new(file)
        .unwrap()
        .build()
        .unwrap();
    let n: usize = reader.map(|b| b.unwrap().num_rows()).sum();
    assert!(n >= 200, "{n} feature rows");
    assert!(log.contains("quant: label: match between runs:"), "{log}");
}

/// Dimethyl triplex (`make_quant_fixtures.py --dimethyl`): every amine
/// labelled, L:M:H of 1:1:1, 1:0.5:0.25 or 0.25:1:2 (of the fully labelled
/// monoisotopes), the deuterated channels eluting 1.5 s and 3 s early and
/// carrying their lower isotopologues at 99 % 2H, the medium overlapping the
/// light M+4.
#[test]
fn dimethyl_triplex_learns_shifts_and_unmixes_channels() {
    let tmp = tempfile::tempdir().unwrap();
    let pin = tmp.path().join("d.pin");
    let out = Command::new(env!("CARGO_BIN_EXE_andes"))
        .args([
            "--spectrum",
            fixture("test-fixtures/dimethyl_run1.mzML.gz")
                .to_str()
                .unwrap(),
            "--database",
            fixture("test-fixtures/BSA.fasta").to_str().unwrap(),
            "--threads",
            "2",
            "--output-pin",
            pin.to_str().unwrap(),
            "--label",
            "dimethyl3",
            "--label-heavy-enrichment",
            "1",
            "--rescore-native",
        ])
        .output()
        .expect("run andes");
    let log = String::from_utf8_lossy(&out.stderr).to_string();
    assert!(out.status.success(), "andes failed:\n{log}");

    // The learned RT offsets of the deuterated channels.
    let line = log
        .lines()
        .find(|l| l.contains("channel RT offsets"))
        .unwrap_or_else(|| panic!("{log}"));
    let offsets: Vec<f64> = line
        .split("channel RT offsets")
        .nth(1)
        .unwrap()
        .split_whitespace()
        .filter_map(|t| t.parse().ok())
        .collect();
    assert_eq!(offsets.len(), 3, "{line}");
    assert!((offsets[1] + 1.5).abs() < 0.2, "{line}");
    assert!((offsets[2] + 3.0).abs() < 0.2, "{line}");

    let rows = read_rows(&tmp.path().join("d.label.parquet"));
    let count = |ch: u8, st: &str| {
        rows.iter()
            .filter(|r| r.channel == ch && r.status == st)
            .count()
    };
    for ch in [1, 2, 3] {
        assert!(count(ch, "match_channel") >= 40, "channel {ch}:\n{log}");
        assert_eq!(count(ch, "below_detection"), 0, "channel {ch}");
    }

    // Monoisotope ratios from the molar ones: the fully labelled
    // monoisotope's share of a channel is its natural π0 (the formula the fit
    // uses, carbamidomethyl C as in the search) times 0.99 per 2H atom.
    let mut spec = quant::label::LabelSpec::parse("dimethyl3").unwrap();
    spec.enrichment_heavy = 1.0;
    spec.enrichment_2h_chemical = 1.0;
    let sites = spec.sites();
    let mono_share = |seq: &str, c: usize| {
        let f = quant::formula::Formula::peptide(
            seq.bytes()
                .map(|r| (r, if r == b'C' { 57.02146 } else { 0.0 })),
        );
        let counts: Vec<u32> = sites
            .iter()
            .map(|s| match s {
                quant::label::Site::NTerm => 1,
                quant::label::Site::Residue(r) => seq.bytes().filter(|b| b == r).count() as u32,
            })
            .collect();
        let d = spec.channel_atoms(c, &sites, &counts).deuterium;
        spec.channel_envelope(c, &sites, &counts, &f).at(0) * 0.99f64.powi(d as i32)
    };
    let mut groups: HashMap<(String, u8), [Option<f64>; 3]> = HashMap::new();
    let seq_of: HashMap<(String, u8), String> = rows
        .iter()
        .map(|r| ((r.peptide.clone(), r.charge), r.peptide.clone()))
        .collect();
    for r in rows.iter().filter(|r| r.channel > 1) {
        if let Some(v) = r.log2_ratio {
            groups.entry((r.peptide.clone(), r.charge)).or_default()[r.channel as usize - 1] =
                Some(v);
        }
    }
    let designs: [[f64; 2]; 3] = [[0.0, 0.0], [-1.0, -2.0], [2.0, 3.0]];
    let (mut good, mut total) = (0usize, 0usize);
    for (key, v) in &groups {
        let (Some(m), Some(h)) = (v[1], v[2]) else {
            continue;
        };
        // The bare sequence: the peptide identity without its brackets.
        let seq: String = seq_of[key]
            .split('[')
            .enumerate()
            .map(|(i, p)| {
                if i == 0 {
                    p
                } else {
                    p.split_once(']').map_or(p, |x| x.1)
                }
            })
            .collect();
        let l0 = mono_share(&seq, 0);
        let mono = [
            m - (l0 / mono_share(&seq, 1)).log2(),
            h - (l0 / mono_share(&seq, 2)).log2(),
        ];
        let dev = designs
            .iter()
            .map(|d| (mono[0] - d[0]).abs().max((mono[1] - d[1]).abs()))
            .fold(f64::INFINITY, f64::min);
        total += 1;
        good += usize::from(dev < 0.03);
    }
    assert!(total >= 120, "{total} groups with both ratios");
    assert!(
        good as f64 >= 0.9 * total as f64,
        "{good} of {total} groups at a design ratio"
    );
}
