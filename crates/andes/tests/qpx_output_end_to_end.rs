//! End-to-end test for the QPX `--output-parquet` `.idparquet/` writer.
//!
//! Runs a default search over the BSA fixture with `--output-parquet`, then
//! opens `psms.parquet` with the `arrow`/`parquet` reader and asserts the
//! column names + Arrow types match the OpenMS QPX 1.0 PSM schema and that at
//! least one row was written. Also checks the QPX schema metadata is present and
//! that the sibling `search_params.parquet` was produced.

use std::collections::HashMap;
use std::path::PathBuf;
use std::process::Command;

use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

/// Resolve a path relative to the workspace root (crates/andes → workspace root).
fn fixture(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(rel)
        .canonicalize()
        .unwrap_or_else(|e| panic!("canonicalize {rel}: {e}"))
}

#[test]
fn search_writes_qpx_idparquet_bundle() {
    let dir = tempfile::tempdir().expect("tempdir");
    let pin_path = dir.path().join("out.pin");
    let qpx_dir = dir.path().join("out.idparquet");

    let bsa_mgf = fixture("test-fixtures/test.mgf.gz");
    let bsa_fasta = fixture("test-fixtures/BSA.fasta");

    // ── Run a default search asking for the QPX bundle ───────────────────────
    let status = Command::new(env!("CARGO_BIN_EXE_andes"))
        .arg("--spectrum")
        .arg(&bsa_mgf)
        .arg("--database")
        .arg(&bsa_fasta)
        .arg("--output-pin")
        .arg(&pin_path)
        .arg("--output-parquet")
        .arg(&qpx_dir)
        .status()
        .expect("run andes search");
    assert!(status.success(), "search should exit 0, got: {status}");

    // ── The bundle directory + its three members must exist ──────────────────
    assert!(qpx_dir.is_dir(), "QPX bundle dir should be created");
    let psms_path = qpx_dir.join("psms.parquet");
    let sp_path = qpx_dir.join("search_params.parquet");
    let prot_path = qpx_dir.join("proteins.parquet");
    assert!(psms_path.exists(), "psms.parquet should be written");
    assert!(sp_path.exists(), "search_params.parquet should be written");
    assert!(prot_path.exists(), "proteins.parquet should be written");

    // ── Open psms.parquet and assert schema column names + Arrow types ───────
    let file = std::fs::File::open(&psms_path).expect("open psms.parquet");
    let builder = ParquetRecordBatchReaderBuilder::try_new(file).expect("parquet reader builder");
    let schema = builder.schema().clone();

    // QPX schema metadata must be present on the file.
    let md = schema.metadata();
    assert_eq!(md.get("qpx_version").map(String::as_str), Some("1.0"));
    assert_eq!(md.get("creator").map(String::as_str), Some("andes"));
    assert_eq!(md.get("software_provider").map(String::as_str), Some("andes"));
    assert_eq!(md.get("file_type").map(String::as_str), Some("psms"));
    assert!(md.contains_key("creation_date"), "creation_date metadata missing");
    assert!(md.contains_key("uuid"), "uuid metadata missing");

    // Expected (column name, Arrow DataType string) for the scalar columns +
    // the list/struct columns (checked by name + top-level DataType variant).
    let by_name: HashMap<&str, &arrow::datatypes::DataType> = schema
        .fields()
        .iter()
        .map(|f| (f.name().as_str(), f.data_type()))
        .collect();

    use arrow::datatypes::DataType;
    let scalar_expected: &[(&str, DataType)] = &[
        ("sequence", DataType::Utf8),
        ("peptidoform", DataType::Utf8),
        ("precursor_charge", DataType::Int32),
        ("posterior_error_probability", DataType::Float64),
        ("is_decoy", DataType::Boolean),
        ("calculated_mz", DataType::Float64),
        ("observed_mz", DataType::Float64),
        ("predicted_rt", DataType::Float64),
        ("reference_file_name", DataType::Utf8),
        ("cv_params", DataType::Utf8),
        ("scan", DataType::Int32),
        ("rt", DataType::Float64),
        ("ion_mobility", DataType::Float64),
        ("spectrum_reference", DataType::Utf8),
        ("score", DataType::Float64),
        ("score_type", DataType::Utf8),
        ("higher_score_better", DataType::Boolean),
        ("hit_index", DataType::Int32),
        ("peptide_identification_index", DataType::Int32),
        ("run_identifier", DataType::Utf8),
    ];
    for (name, dt) in scalar_expected {
        let actual = by_name
            .get(name)
            .unwrap_or_else(|| panic!("psms.parquet missing column `{name}`"));
        assert_eq!(*actual, dt, "column `{name}` has wrong Arrow type");
    }

    // List/struct columns: present + the top-level type is a List.
    for name in [
        "modifications",
        "additional_scores",
        "protein_accessions",
        "psm_metavalues",
        "spectrum_metavalues",
        "mz_array",
        "intensity_array",
        "charge_array",
        "ion_type_array",
    ] {
        let dt = by_name
            .get(name)
            .unwrap_or_else(|| panic!("psms.parquet missing list column `{name}`"));
        let DataType::List(elem) = dt else {
            panic!("column `{name}` should be a List, got {dt:?}");
        };
        // QPX/Parquet LIST convention: the list element field must be named
        // `element` (arrow-rs defaults to `item`; regression guard for the
        // OpenMS interop fix).
        assert_eq!(
            elem.name(),
            "element",
            "column `{name}` list element field must be named `element`, got `{}`",
            elem.name()
        );
    }

    // Column count must equal the full QPX PSM schema (29 columns).
    assert_eq!(schema.fields().len(), 29, "psms.parquet column count");

    // ── At least one data row ────────────────────────────────────────────────
    let mut reader = builder.build().expect("build parquet reader");
    let total: usize = std::iter::from_fn(|| reader.next())
        .map(|b| b.expect("read batch").num_rows())
        .sum();
    assert!(total >= 1, "psms.parquet should have >=1 row, got {total}");
}
