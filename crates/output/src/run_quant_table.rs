//! Per-run PSM table `<stem>.quant.parquet`: every rank-1 PSM of one run
//! (targets and decoys) with the fields the multi-run quantification stage
//! needs. The file carries a schema version in its key-value metadata; a
//! reader rejects a different major version.

use std::fs::File;
use std::path::Path;
use std::sync::Arc;

use arrow::array::{
    Array, ArrayRef, BooleanArray, BooleanBuilder, Float32Array, Float32Builder, Float64Array,
    Float64Builder, ListArray, ListBuilder, StringArray, StringBuilder, UInt8Array, UInt8Builder,
};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::ArrowWriter;
use parquet::basic::{Compression, ZstdLevel};
use parquet::file::properties::WriterProperties;
use parquet::format::KeyValue;

/// Major version of the table layout. A reader rejects any other major.
pub const RUN_QUANT_SCHEMA_MAJOR: u32 = 1;

const SCHEMA_KEY: &str = "andes.run_quant.schema";
const RUN_KEY: &str = "andes.run_quant.run";
const SCHEMA_VERSION: &str = "1.0";
/// Rows per Parquet row group / record batch.
const BATCH_ROWS: usize = 65_536;

/// One rank-1 PSM of a run, target or decoy.
#[derive(Debug, Clone, PartialEq)]
pub struct RunPsmRow {
    /// PIN SpecId, joins the Percolator PSMId.
    pub spec_id: String,
    /// MS2 retention time, seconds.
    pub rt_s: f64,
    /// `"<peptidoform>/<charge>"`.
    pub precursor_key: String,
    pub peptidoform: String,
    pub charge: u8,
    /// Theoretical monoisotopic m/z.
    pub mz: f64,
    pub is_decoy: bool,
    /// andes RawScore of the rank-1 match.
    pub search_score: f32,
    pub faims_cv: Option<f32>,
    /// Theoretical isotope fractions (sum 1).
    pub envelope: Vec<f32>,
    /// `;`-joined accessions.
    pub proteins: String,
}

fn schema() -> Schema {
    Schema::new(vec![
        Field::new("spec_id", DataType::Utf8, false),
        Field::new("rt_s", DataType::Float64, false),
        Field::new("precursor_key", DataType::Utf8, false),
        Field::new("peptidoform", DataType::Utf8, false),
        Field::new("charge", DataType::UInt8, false),
        Field::new("mz", DataType::Float64, false),
        Field::new("is_decoy", DataType::Boolean, false),
        Field::new("search_score", DataType::Float32, false),
        Field::new("faims_cv", DataType::Float32, true),
        Field::new(
            "envelope",
            DataType::List(Arc::new(Field::new("item", DataType::Float32, true))),
            false,
        ),
        Field::new("proteins", DataType::Utf8, false),
    ])
}

fn batch(schema: &Arc<Schema>, rows: &[RunPsmRow]) -> Result<RecordBatch, String> {
    let mut spec_id = StringBuilder::new();
    let mut rt_s = Float64Builder::new();
    let mut key = StringBuilder::new();
    let mut pepf = StringBuilder::new();
    let mut charge = UInt8Builder::new();
    let mut mz = Float64Builder::new();
    let mut decoy = BooleanBuilder::new();
    let mut score = Float32Builder::new();
    let mut cv = Float32Builder::new();
    let mut env = ListBuilder::new(Float32Builder::new());
    let mut prot = StringBuilder::new();
    for r in rows {
        spec_id.append_value(&r.spec_id);
        rt_s.append_value(r.rt_s);
        key.append_value(&r.precursor_key);
        pepf.append_value(&r.peptidoform);
        charge.append_value(r.charge);
        mz.append_value(r.mz);
        decoy.append_value(r.is_decoy);
        score.append_value(r.search_score);
        cv.append_option(r.faims_cv);
        env.values().append_slice(&r.envelope);
        env.append(true);
        prot.append_value(&r.proteins);
    }
    let cols: Vec<ArrayRef> = vec![
        Arc::new(spec_id.finish()),
        Arc::new(rt_s.finish()),
        Arc::new(key.finish()),
        Arc::new(pepf.finish()),
        Arc::new(charge.finish()),
        Arc::new(mz.finish()),
        Arc::new(decoy.finish()),
        Arc::new(score.finish()),
        Arc::new(cv.finish()),
        Arc::new(env.finish()),
        Arc::new(prot.finish()),
    ];
    RecordBatch::try_new(schema.clone(), cols).map_err(|e| e.to_string())
}

fn write_versioned(
    path: &Path,
    run_stem: &str,
    rows: &[RunPsmRow],
    version: &str,
) -> Result<(), String> {
    let err = |e: &dyn std::fmt::Display| format!("{}: {e}", path.display());
    let schema = Arc::new(schema());
    let props = WriterProperties::builder()
        .set_compression(Compression::ZSTD(ZstdLevel::default()))
        .set_max_row_group_size(BATCH_ROWS)
        .set_key_value_metadata(Some(vec![
            KeyValue::new(SCHEMA_KEY.to_string(), version.to_string()),
            KeyValue::new(RUN_KEY.to_string(), run_stem.to_string()),
        ]))
        .build();
    let file = File::create(path).map_err(|e| err(&e))?;
    let mut writer =
        ArrowWriter::try_new(file, schema.clone(), Some(props)).map_err(|e| err(&e))?;
    for chunk in rows.chunks(BATCH_ROWS) {
        let b = batch(&schema, chunk).map_err(|e| err(&e))?;
        writer.write(&b).map_err(|e| err(&e))?;
    }
    writer.close().map_err(|e| err(&e))?;
    Ok(())
}

/// Write the per-run PSM table for run `run_stem`.
pub fn write_run_quant_table(
    path: &Path,
    run_stem: &str,
    rows: &[RunPsmRow],
) -> Result<(), String> {
    write_versioned(path, run_stem, rows, SCHEMA_VERSION)
}

#[cfg(test)]
pub(crate) fn write_with_schema_version(
    path: &Path,
    run_stem: &str,
    rows: &[RunPsmRow],
    version: &str,
) -> Result<(), String> {
    write_versioned(path, run_stem, rows, version)
}

fn column<'a, T: 'static>(b: &'a RecordBatch, i: usize, path: &Path) -> Result<&'a T, String> {
    b.column(i).as_any().downcast_ref::<T>().ok_or_else(|| {
        format!(
            "{}: column {} has unexpected type {}",
            path.display(),
            b.schema().field(i).name(),
            b.column(i).data_type()
        )
    })
}

/// Read a per-run PSM table; returns `(run_stem, rows)`. Every error names
/// `path`.
pub fn read_run_quant_table(path: &Path) -> Result<(String, Vec<RunPsmRow>), String> {
    let err = |e: &dyn std::fmt::Display| format!("{}: {e}", path.display());
    let file = File::open(path).map_err(|e| err(&e))?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file).map_err(|e| err(&e))?;
    let meta = builder
        .metadata()
        .file_metadata()
        .key_value_metadata()
        .cloned()
        .unwrap_or_default();
    let get = |key: &str| {
        meta.iter()
            .find(|kv| kv.key == key)
            .and_then(|kv| kv.value.clone())
    };
    let version =
        get(SCHEMA_KEY).ok_or_else(|| err(&format!("missing metadata key {SCHEMA_KEY}")))?;
    let major: u32 = version
        .split('.')
        .next()
        .and_then(|m| m.parse().ok())
        .ok_or_else(|| err(&format!("unparseable schema version {version:?}")))?;
    if major != RUN_QUANT_SCHEMA_MAJOR {
        return Err(err(&format!(
            "unsupported run-quant table schema {version} (this build reads major {RUN_QUANT_SCHEMA_MAJOR})"
        )));
    }
    let run = get(RUN_KEY).ok_or_else(|| err(&format!("missing metadata key {RUN_KEY}")))?;

    let reader = builder.build().map_err(|e| err(&e))?;
    let mut rows = Vec::new();
    for b in reader {
        let b = b.map_err(|e| err(&e))?;
        if b.num_columns() != 11 {
            return Err(err(&format!(
                "expected 11 columns, found {}",
                b.num_columns()
            )));
        }
        let spec_id = column::<StringArray>(&b, 0, path)?;
        let rt_s = column::<Float64Array>(&b, 1, path)?;
        let key = column::<StringArray>(&b, 2, path)?;
        let pepf = column::<StringArray>(&b, 3, path)?;
        let charge = column::<UInt8Array>(&b, 4, path)?;
        let mz = column::<Float64Array>(&b, 5, path)?;
        let decoy = column::<BooleanArray>(&b, 6, path)?;
        let score = column::<Float32Array>(&b, 7, path)?;
        let cv = column::<Float32Array>(&b, 8, path)?;
        let env = column::<ListArray>(&b, 9, path)?;
        let prot = column::<StringArray>(&b, 10, path)?;
        for i in 0..b.num_rows() {
            let e = env.value(i);
            let e = e
                .as_any()
                .downcast_ref::<Float32Array>()
                .ok_or_else(|| err(&"envelope items are not Float32"))?;
            rows.push(RunPsmRow {
                spec_id: spec_id.value(i).to_string(),
                rt_s: rt_s.value(i),
                precursor_key: key.value(i).to_string(),
                peptidoform: pepf.value(i).to_string(),
                charge: charge.value(i),
                mz: mz.value(i),
                is_decoy: decoy.value(i),
                search_score: score.value(i),
                faims_cv: (!cv.is_null(i)).then(|| cv.value(i)),
                envelope: e.values().to_vec(),
                proteins: prot.value(i).to_string(),
            });
        }
    }
    Ok((run, rows))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn row(i: usize) -> RunPsmRow {
        RunPsmRow {
            spec_id: format!("a/scan={i}_2_1"),
            rt_s: 1234.5 + i as f64,
            precursor_key: "PEPTIDEK/2".into(),
            peptidoform: "PEPTIDEK".into(),
            charge: 2,
            mz: 465.7380,
            is_decoy: i % 2 == 1,
            search_score: 41.0,
            faims_cv: if i % 2 == 0 { Some(-45.0) } else { None },
            envelope: vec![0.6, 0.3, 0.08, 0.02],
            proteins: "P1;P2".into(),
        }
    }

    #[test]
    fn round_trip_preserves_rows_and_run() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("a.quant.parquet");
        let rows = vec![RunPsmRow {
            spec_id: "a/scan=7_2_1".into(),
            rt_s: 1234.5,
            precursor_key: "PEPTIDEK/2".into(),
            peptidoform: "PEPTIDEK".into(),
            charge: 2,
            mz: 465.7380,
            is_decoy: false,
            search_score: 41.0,
            faims_cv: Some(-45.0),
            envelope: vec![0.6, 0.3, 0.08, 0.02],
            proteins: "P1;P2".into(),
        }];
        write_run_quant_table(&p, "a", &rows).unwrap();
        let (stem, back) = read_run_quant_table(&p).unwrap();
        assert_eq!(stem, "a");
        assert_eq!(back, rows);
    }

    #[test]
    fn round_trip_mixed_nulls_decoys_and_empty() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("m.quant.parquet");
        let rows: Vec<RunPsmRow> = (0..5).map(row).collect();
        write_run_quant_table(&p, "m", &rows).unwrap();
        assert_eq!(read_run_quant_table(&p).unwrap(), ("m".to_string(), rows));

        let e = dir.path().join("e.quant.parquet");
        write_run_quant_table(&e, "e", &[]).unwrap();
        assert_eq!(read_run_quant_table(&e).unwrap(), ("e".to_string(), vec![]));
    }

    #[test]
    fn quant_table_rejects_unknown_major_version() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("b.quant.parquet");
        write_with_schema_version(&p, "b", &[], "2.0").unwrap();
        let e = read_run_quant_table(&p).unwrap_err();
        assert!(
            e.contains("b.quant.parquet") && e.contains("schema 2.0"),
            "{e}"
        );
    }

    #[test]
    fn missing_and_non_parquet_files_name_the_path() {
        let dir = tempfile::tempdir().unwrap();
        let missing = dir.path().join("nope.quant.parquet");
        let e = read_run_quant_table(&missing).unwrap_err();
        assert!(e.contains("nope.quant.parquet"), "{e}");

        let junk = dir.path().join("junk.quant.parquet");
        std::fs::write(&junk, b"not parquet").unwrap();
        let e = read_run_quant_table(&junk).unwrap_err();
        assert!(e.contains("junk.quant.parquet"), "{e}");
    }
}
