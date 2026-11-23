//! Per-run DIA precursor quantities of a study: one row per (run, precursor) reported at the
//! study's q threshold, written by `andes quant finalize` on a DIA study as
//! `<base>.dia_quant.parquet`.
//!
//! Columns: `filename` (run stem), `peptide` (as in the precursor PIN), `charge`, `proteins`,
//! `evidence_type` (`identified` for the run's own nominees, `transfer` for candidates
//! transferred from another run), `q_value` (the reported precursor q), `quantity` (the default
//! quantity), `quantity_frag`, `quantity_ms1` and `quantity_normalized` (`quantity` scaled by
//! the run's median-level factor). Missing quantities are null, never 0.

use std::fs::File;
use std::path::Path;
use std::sync::Arc;

use arrow::array::{Array, ArrayRef, Float64Array, StringArray, UInt8Array};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

use crate::versioned::{write_table, TableVersion, BATCH_ROWS};

pub const DIA_QUANT_SCHEMA_MAJOR: u32 = 1;
const SCHEMA_VERSION: &str = "1.0";
const TABLE: TableVersion = TableVersion {
    key: "andes.dia_quant.schema",
    version: SCHEMA_VERSION,
    major: DIA_QUANT_SCHEMA_MAJOR,
    what: "DIA quantity table",
};

/// `evidence_type` of a run's own nominee.
pub const EVIDENCE_IDENTIFIED: &str = "identified";
/// `evidence_type` of a candidate transferred from another run.
pub const EVIDENCE_TRANSFER: &str = "transfer";

/// One (run, precursor) row.
#[derive(Debug, Clone, PartialEq)]
pub struct DiaQuantRow {
    pub filename: String,
    pub peptide: String,
    pub charge: u8,
    pub proteins: String,
    pub transfer: bool,
    pub q_value: f64,
    pub quantity: Option<f64>,
    pub quantity_frag: Option<f64>,
    pub quantity_ms1: Option<f64>,
    pub quantity_normalized: Option<f64>,
}

fn schema() -> Arc<Schema> {
    Arc::new(Schema::new(vec![
        Field::new("filename", DataType::Utf8, false),
        Field::new("peptide", DataType::Utf8, false),
        Field::new("charge", DataType::UInt8, false),
        Field::new("proteins", DataType::Utf8, false),
        Field::new("evidence_type", DataType::Utf8, false),
        Field::new("q_value", DataType::Float64, false),
        Field::new("quantity", DataType::Float64, true),
        Field::new("quantity_frag", DataType::Float64, true),
        Field::new("quantity_ms1", DataType::Float64, true),
        Field::new("quantity_normalized", DataType::Float64, true),
    ]))
}

/// Writes the table; `meta` is extra key-value metadata (e.g. the default quantity).
pub fn write_dia_quant(
    path: &Path,
    rows: &[DiaQuantRow],
    meta: &[(&str, String)],
) -> Result<(), String> {
    let sch = schema();
    let n = rows.len();
    let batches = (0..n.div_ceil(BATCH_ROWS).max(1)).map(|b| {
        let r = &rows[b * BATCH_ROWS..((b + 1) * BATCH_ROWS).min(n)];
        let f = |g: &dyn Fn(&DiaQuantRow) -> Option<f64>| -> ArrayRef {
            Arc::new(Float64Array::from(r.iter().map(g).collect::<Vec<_>>()))
        };
        let cols: Vec<ArrayRef> = vec![
            Arc::new(StringArray::from_iter_values(r.iter().map(|x| &x.filename))),
            Arc::new(StringArray::from_iter_values(r.iter().map(|x| &x.peptide))),
            Arc::new(UInt8Array::from(r.iter().map(|x| x.charge).collect::<Vec<_>>())),
            Arc::new(StringArray::from_iter_values(r.iter().map(|x| &x.proteins))),
            Arc::new(StringArray::from_iter_values(r.iter().map(|x| {
                if x.transfer {
                    EVIDENCE_TRANSFER
                } else {
                    EVIDENCE_IDENTIFIED
                }
            }))),
            Arc::new(Float64Array::from(r.iter().map(|x| x.q_value).collect::<Vec<_>>())),
            f(&|x| x.quantity),
            f(&|x| x.quantity_frag),
            f(&|x| x.quantity_ms1),
            f(&|x| x.quantity_normalized),
        ];
        RecordBatch::try_new(sch.clone(), cols).map_err(|e| e.to_string())
    });
    write_table(path, &TABLE, SCHEMA_VERSION, sch.clone(), batches, meta)
}

/// Reads a DIA quantity table written by [`write_dia_quant`], checking its schema version.
pub fn read_dia_quant(path: &Path) -> Result<Vec<DiaQuantRow>, String> {
    let err = |e: &dyn std::fmt::Display| format!("{}: {e}", path.display());
    let f = File::open(path).map_err(|e| err(&e))?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(f).map_err(|e| err(&e))?;
    let version = builder
        .metadata()
        .file_metadata()
        .key_value_metadata()
        .and_then(|kv| kv.iter().find(|k| k.key == TABLE.key))
        .and_then(|k| k.value.clone())
        .ok_or_else(|| format!("{}: not a {}", path.display(), TABLE.what))?;
    let major: u32 = version.split('.').next().and_then(|m| m.parse().ok()).unwrap_or(0);
    if major != DIA_QUANT_SCHEMA_MAJOR {
        return Err(format!(
            "{}: {} schema {version}, this andes reads {DIA_QUANT_SCHEMA_MAJOR}.x",
            path.display(),
            TABLE.what
        ));
    }
    let reader = builder.with_batch_size(BATCH_ROWS).build().map_err(|e| err(&e))?;
    let mut out = Vec::new();
    for b in reader {
        let b = b.map_err(|e| err(&e))?;
        let bad = |n: &str| format!("{}: column {n} is missing or has the wrong type", path.display());
        let get = |n: &str| b.column_by_name(n).map(|c| c.as_any());
        let s = |n: &str| get(n).and_then(|c| c.downcast_ref::<StringArray>()).ok_or_else(|| bad(n));
        let d = |n: &str| get(n).and_then(|c| c.downcast_ref::<Float64Array>()).ok_or_else(|| bad(n));
        let (file, pep, prot, ev) = (s("filename")?, s("peptide")?, s("proteins")?, s("evidence_type")?);
        let charge = get("charge")
            .and_then(|c| c.downcast_ref::<UInt8Array>())
            .ok_or_else(|| bad("charge"))?;
        let q = d("q_value")?;
        let (qd, qf, qm, qn) = (d("quantity")?, d("quantity_frag")?, d("quantity_ms1")?, d("quantity_normalized")?);
        let opt = |a: &Float64Array, i: usize| a.is_valid(i).then(|| a.value(i));
        for i in 0..b.num_rows() {
            out.push(DiaQuantRow {
                filename: file.value(i).to_string(),
                peptide: pep.value(i).to_string(),
                charge: charge.value(i),
                proteins: prot.value(i).to_string(),
                transfer: ev.value(i) == EVIDENCE_TRANSFER,
                q_value: q.value(i),
                quantity: opt(qd, i),
                quantity_frag: opt(qf, i),
                quantity_ms1: opt(qm, i),
                quantity_normalized: opt(qn, i),
            });
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn schema_version_metadata() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("study.dia_quant.parquet");
        let rows = vec![
            DiaQuantRow {
                filename: "run_A".into(),
                peptide: "_.C+57.02146PEPTIDEK.-".into(),
                charge: 2,
                proteins: "P1;P2".into(),
                transfer: false,
                q_value: 0.001,
                quantity: Some(2.0e6),
                quantity_frag: Some(2.0e6),
                quantity_ms1: None,
                quantity_normalized: Some(1.8e6),
            },
            DiaQuantRow {
                filename: "run_B".into(),
                peptide: "_.PEPTIDEK.-".into(),
                charge: 3,
                proteins: "P1".into(),
                transfer: true,
                q_value: 0.008,
                quantity: None,
                quantity_frag: None,
                quantity_ms1: Some(4.0e7),
                quantity_normalized: None,
            },
        ];
        write_dia_quant(&p, &rows, &[("andes.dia_quant.default", "quantity_frag".to_string())]).unwrap();
        let f = File::open(&p).unwrap();
        let b = ParquetRecordBatchReaderBuilder::try_new(f).unwrap();
        let kv = b.metadata().file_metadata().key_value_metadata().unwrap().clone();
        let get = |k: &str| kv.iter().find(|x| x.key == k).and_then(|x| x.value.clone());
        assert_eq!(get("andes.dia_quant.schema").as_deref(), Some("1.0"));
        assert_eq!(get("andes.dia_quant.default").as_deref(), Some("quantity_frag"));
        let names: Vec<String> = b.schema().fields().iter().map(|f| f.name().clone()).collect();
        assert_eq!(
            names,
            [
                "filename", "peptide", "charge", "proteins", "evidence_type", "q_value", "quantity",
                "quantity_frag", "quantity_ms1", "quantity_normalized"
            ]
        );
        assert_eq!(read_dia_quant(&p).unwrap(), rows);
        let other = dir.path().join("pass2.parquet");
        crate::dia::pass2::write_pass2(&other, &Default::default()).unwrap();
        assert!(read_dia_quant(&other).unwrap_err().contains("not a DIA quantity table"));
    }
}
