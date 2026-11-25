//! Goldens for the default DDA `--rescore` path: the rescored target and decoy PSMs, the protein
//! groups and the QPX PSM Parquet of one E. coli search. Each golden is (rows, digest of every
//! non-float column), so scores may move in the last digits while identities stay pinned.

use std::path::{Path, PathBuf};
use std::process::Command;

use arrow::array::{Array, RecordBatch};
use arrow::datatypes::DataType;
use arrow::util::display::array_value_to_string;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

type PathGolden = (usize, u64);

const RESCORE_TARGETS: PathGolden = (1636, 0x3b8124c7918e16dd);
const RESCORE_DECOYS: PathGolden = (1712, 0x63b797d633b82619);
const PROTEIN_GROUPS: PathGolden = (1264, 0xf23641f48de4283e);
const QPX_PSMS: PathGolden = (3348, 0x00a6175aafd4b9cc);

fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("canonicalize workspace root")
}

fn is_float_field(v: &str) -> bool {
    v.parse::<f64>().is_ok() && v.parse::<i64>().is_err()
}

fn float_columns(rows: &[String]) -> Vec<bool> {
    let mut float = Vec::new();
    for row in rows {
        for (i, v) in row.split('\t').enumerate() {
            if i >= float.len() {
                float.resize(i + 1, false);
            }
            float[i] |= is_float_field(v);
        }
    }
    float
}

fn identity_digest(header: &str, rows: &[String]) -> (usize, u64) {
    let float = float_columns(rows);
    let mut projected: Vec<String> = rows
        .iter()
        .map(|row| {
            row.split('\t')
                .enumerate()
                .filter(|(i, _)| !float[*i])
                .map(|(_, v)| v)
                .collect::<Vec<&str>>()
                .join("\t")
        })
        .collect();
    projected.sort();
    let mut h = model::hash::FNV1A64_OFFSET;
    let mut feed = |bytes: &[u8]| h = model::hash::fnv1a64_extend(h, bytes);
    feed(header.as_bytes());
    feed(b"\n");
    for row in &projected {
        feed(row.as_bytes());
        feed(b"\n");
    }
    (rows.len(), h)
}

fn text_digest(path: &Path) -> PathGolden {
    let text = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    let mut lines = text.lines();
    let header = lines.next().expect("header line");
    let rows: Vec<String> = lines.map(str::to_string).collect();
    identity_digest(header, &rows)
}

fn batch_rows(batch: &RecordBatch, keep: &[usize]) -> Vec<(String, String)> {
    (0..batch.num_rows())
        .map(|r| {
            let mut key = String::new();
            let mut row = String::new();
            for (n, &c) in keep.iter().enumerate() {
                let col = batch.column(c);
                let v = if col.is_null(r) {
                    String::new()
                } else {
                    array_value_to_string(col, r).expect("format cell")
                };
                if n == 0 {
                    key = v.clone();
                }
                if n > 0 {
                    row.push('\t');
                }
                row.push_str(&v);
            }
            (key, row)
        })
        .collect()
}

fn parquet_digest(path: &Path) -> PathGolden {
    let file = std::fs::File::open(path).unwrap_or_else(|e| panic!("open {}: {e}", path.display()));
    let builder = ParquetRecordBatchReaderBuilder::try_new(file).expect("parquet reader");
    let schema = builder.schema().clone();
    let keep: Vec<usize> = schema
        .fields()
        .iter()
        .enumerate()
        .filter(|(_, f)| {
            matches!(
                f.data_type(),
                DataType::Utf8
                    | DataType::LargeUtf8
                    | DataType::Boolean
                    | DataType::Int8
                    | DataType::Int16
                    | DataType::Int32
                    | DataType::Int64
                    | DataType::UInt8
                    | DataType::UInt16
                    | DataType::UInt32
                    | DataType::UInt64
            )
        })
        .map(|(i, _)| i)
        .collect();
    assert!(
        matches!(
            schema.field(keep[0]).data_type(),
            DataType::Utf8 | DataType::LargeUtf8
        ),
        "the first digested column must be the first string column"
    );
    let header: Vec<&str> = keep.iter().map(|&i| schema.field(i).name().as_str()).collect();
    let mut rows: Vec<(String, String)> = Vec::new();
    for batch in builder.build().expect("parquet batches") {
        rows.extend(batch_rows(&batch.expect("read batch"), &keep));
    }
    rows.sort();
    let rows: Vec<String> = rows.into_iter().map(|(_, row)| row).collect();
    identity_digest(&header.join("\t"), &rows)
}

fn run_rescore(threads: usize, dir: &Path) -> [PathGolden; 4] {
    let qpx = dir.join("qpx");
    let output = Command::new(env!("CARGO_BIN_EXE_andes"))
        .current_dir(workspace_root())
        .args([
            "--spectrum",
            "test-fixtures/test.mgf.gz",
            "--database",
            "test-fixtures/ecoli.fasta",
            "--fragmentation",
            "HCD",
            "--fragment-tol-ppm",
            "20",
            "--predict",
            "none",
            "--max-missed-cleavages",
            "1",
            "--rescore",
        ])
        .arg("--output-pin")
        .arg(dir.join("run.pin"))
        .arg("--output-parquet")
        .arg(&qpx)
        .args(["--threads", &threads.to_string()])
        .output()
        .expect("spawn andes");
    assert!(
        output.status.success(),
        "andes --rescore --threads {threads} failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    [
        text_digest(&dir.join("run.percolator.target.psms.txt")),
        text_digest(&dir.join("run.percolator.decoy.psms.txt")),
        text_digest(&dir.join("run.protein_groups.tsv")),
        parquet_digest(&qpx.join("psms.parquet")),
    ]
}

fn show(found: &[PathGolden; 4]) -> String {
    let names = ["RESCORE_TARGETS", "RESCORE_DECOYS", "PROTEIN_GROUPS", "QPX_PSMS"];
    names
        .iter()
        .zip(found)
        .map(|(n, (rows, digest))| format!("const {n}: PathGolden = ({rows}, {digest:#x});"))
        .collect::<Vec<_>>()
        .join("\n")
}

#[test]
#[ignore = "runs a full E. coli search with rescoring"]
fn ecoli_rescore_path_matches_golden() {
    let dir = tempfile::tempdir().expect("tempdir");
    let found = run_rescore(1, dir.path());
    let golden = [RESCORE_TARGETS, RESCORE_DECOYS, PROTEIN_GROUPS, QPX_PSMS];
    assert_eq!(
        found,
        golden,
        "(rows, digest) differ from the goldens; found:\n{}",
        show(&found)
    );
}

#[test]
#[ignore = "runs a full E. coli search with rescoring"]
fn ecoli_rescore_path_is_thread_stable() {
    let one = tempfile::tempdir().expect("tempdir");
    let four = tempfile::tempdir().expect("tempdir");
    let single = run_rescore(1, one.path());
    let multi = run_rescore(4, four.path());
    assert_eq!(
        single,
        multi,
        "--threads 4 differs from --threads 1; 1 thread:\n{}\n4 threads:\n{}",
        show(&single),
        show(&multi)
    );
}

#[test]
#[ignore = "runs a full E. coli search with rescoring"]
fn ecoli_rescore_path_reports_fdr_levels() {
    let dir = tempfile::tempdir().expect("tempdir");
    run_rescore(1, dir.path());
    let log = dir.path().join("run.statistics.log");
    let text = std::fs::read_to_string(&log).expect("statistics log exists");
    assert!(!text.is_empty(), "run.statistics.log is empty");
    for title in ["Inputs", "Parameters", "Stages", "Counts", "Software", "Methods", "References"]
    {
        let header = format!("──────── {title} ────────");
        assert!(text.contains(&header), "no {header} in run.statistics.log:\n{text}");
    }
    let fasta = workspace_root().join("test-fixtures/ecoli.fasta");
    let sha = search::decoy::file_sha256(&fasta).expect("hash the fixture FASTA");
    let database = text
        .lines()
        .skip_while(|l| !l.trim_start().starts_with("database"))
        .nth(1)
        .unwrap_or_default();
    assert!(database.contains(&format!("sha256 {sha}")), "database checksum line: {database}");
    assert!(text.contains("Percolator ["), "no Percolator citation in:\n{text}");
    let parquet = std::fs::File::open(dir.path().join("run.statistics.parquet"))
        .expect("run.statistics.parquet exists");
    let rows: usize = ParquetRecordBatchReaderBuilder::try_new(parquet)
        .expect("read run.statistics.parquet")
        .build()
        .expect("parquet batches")
        .map(|b| b.expect("read batch").num_rows())
        .sum();
    assert!(rows > 100, "run.statistics.parquet has {rows} rows; every parameter is a row");
}
