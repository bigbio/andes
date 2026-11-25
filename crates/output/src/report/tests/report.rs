use super::*;
use arrow::array::Array;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

fn param(id: &str, value: &str, source: &str) -> ParamRow {
    ParamRow {
        id: id.to_string(),
        value: value.to_string(),
        source: source.to_string(),
        number: None,
    }
}

fn sample() -> RunReport {
    RunReport {
        andes_version: "0.4.0".to_string(),
        command_line: "andes --spectrum run.mzML --database db.fasta".to_string(),
        inputs: vec![
            InputFile {
                role: "spectrum".to_string(),
                path: "run.mzML".to_string(),
                bytes: Some(1234),
                sha256: Ok("ab".repeat(32)),
            },
            InputFile {
                role: "database".to_string(),
                path: "gone.fasta".to_string(),
                bytes: None,
                sha256: Err("No such file".to_string()),
            },
        ],
        params: vec![
            param("max_missed_cleavages", "2", "default"),
            param("precursor_tol", "20ppm", "command line"),
        ],
        stages: vec![StageRow {
            name: "search".to_string(),
            seconds: 1.5,
            peak_rss_bytes: Some(1 << 30),
        }],
        counts: vec![
            CountRow {
                what: "spectra".to_string(),
                fdr: None,
                n: 3348,
            },
            CountRow {
                what: "PSMs".to_string(),
                fdr: Some(0.01),
                n: 1636,
            },
        ],
        models: vec![("scoring".to_string(), "hcd_highres_tryp".to_string())],
        decoys: "shuffled sequences (seed 42)".to_string(),
        decoy_key: "0123".to_string(),
        methods: "Spectra were searched with andes 0.4.0 [1].".to_string(),
        references: vec![Reference {
            key: "kim2014_msgfplus".to_string(),
            text: "Kim S, Pevzner PA (2014). MS-GF+. Nature Communications. doi:10.1038/ncomms6277"
                .to_string(),
        }],
        parquet_name: "run.statistics.parquet".to_string(),
    }
}

#[test]
fn render_lists_only_non_default_parameters() {
    let text = sample().render();
    assert!(text.contains("precursor_tol  20ppm  (command line)"), "{text}");
    assert!(!text.contains("max_missed_cleavages"), "{text}");
    assert!(
        text.contains("(1 parameters at their default; all are in run.statistics.parquet)"),
        "{text}"
    );
}

#[test]
fn render_shows_a_failed_checksum() {
    let text = sample().render();
    assert!(text.contains("checksum unavailable: No such file"), "{text}");
    assert!(text.contains(&format!("sha256 {}", "ab".repeat(32))), "{text}");
}

#[test]
fn render_sections_are_in_order() {
    let text = sample().render();
    let at: Vec<usize> = [
        "Inputs",
        "Parameters",
        "Stages",
        "Counts",
        "Software",
        "Methods",
        "References",
    ]
    .iter()
    .map(|t| {
        text.find(&format!("──────── {t} ────────"))
            .unwrap_or_else(|| panic!("no {t} header in\n{text}"))
    })
    .collect();
    assert!(at.windows(2).all(|w| w[0] < w[1]), "{at:?}");
}

#[test]
fn render_wraps_long_lines() {
    let mut r = sample();
    r.methods = "word ".repeat(60);
    let text = r.render();
    assert!(text.lines().all(|l| l.chars().count() <= LINE_WIDTH_CHARS), "{text}");
}

#[test]
fn parquet_round_trips_every_row() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("run.statistics.parquet");
    let r = sample();
    r.write_parquet(&path).unwrap();
    assert!(!crate::atomic::partial(&path).exists());
    let file = std::fs::File::open(&path).unwrap();
    let builder = ParquetRecordBatchReaderBuilder::try_new(file).unwrap();
    let version = builder
        .metadata()
        .file_metadata()
        .key_value_metadata()
        .and_then(|kv| kv.iter().find(|k| k.key == "andes.run_report.schema").cloned())
        .and_then(|k| k.value);
    assert_eq!(version.as_deref(), Some("1.0"));
    let batches: Vec<RecordBatch> = builder.build().unwrap().map(|b| b.unwrap()).collect();
    let rows: usize = batches.iter().map(RecordBatch::num_rows).sum();
    let expected = r.inputs.len()
        + r.params.len()
        + 2 * r.stages.len()
        + r.counts.len()
        + r.models.len()
        + 4
        + 1
        + r.references.len();
    assert_eq!(rows, expected);
    let b = &batches[0];
    let section = b.column(0).as_any().downcast_ref::<StringArray>().unwrap();
    let number = b.column(4).as_any().downcast_ref::<Float64Array>().unwrap();
    let fdr = b.column(5).as_any().downcast_ref::<Float64Array>().unwrap();
    let psms = (0..b.num_rows())
        .find(|&i| section.value(i) == "count" && !fdr.is_null(i))
        .unwrap();
    assert_eq!((number.value(psms), fdr.value(psms)), (1636.0, 0.01));
}
