//! Per-run DIA nominee table: the pass-2 inputs of every nominated precursor of a run (targets
//! and decoys), so that a study step can score one run's nominees in another run. Written by
//! `andes --dia --quant-dir` and read by `andes quant extract` on a DIA study.
//!
//! One row per nominee: its precursor (peptide string, charge, decoy, m/z), the predicted RT
//! (calibrated to the run), the raw andes-rt prediction (NaN without one) and whether the PIN
//! row had no prediction, the RT of the scan that nominated it, its pool `dia_llr`, the ion mobility of that scan's
//! isolation window (NaN without one), and its top predicted fragments (m/z and predicted
//! intensity, most intense first) as two list columns.

use std::fs::File;
use std::path::Path;
use std::sync::Arc;

use arrow::array::{
    Array, ArrayRef, BooleanArray, Float32Array, Float32Builder, Float64Array, Float64Builder,
    ListArray, ListBuilder, StringArray, UInt8Array,
};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::ProjectionMask;

use crate::versioned::{write_table, TableVersion, BATCH_ROWS};

pub const DIA_NOMINEE_SCHEMA_MAJOR: u32 = 1;
const SCHEMA_VERSION: &str = "1.2";
const TABLE: TableVersion = TableVersion {
    key: "andes.dia_nominees.schema",
    version: SCHEMA_VERSION,
    major: DIA_NOMINEE_SCHEMA_MAJOR,
    what: "DIA nominee table",
};

/// One nominee as written; borrowed from the caller's rows.
pub struct NomineeRow<'a> {
    pub spec_id: &'a str,
    pub decoy: bool,
    pub charge: u8,
    pub peptide: &'a str,
    pub proteins: &'a str,
    pub seq: &'a str,
    pub prec_mz: f64,
    pub pred_rt: f64,
    pub raw_pred_rt: f64,
    pub rt_missing: bool,
    pub apex_rt: f64,
    pub dia_llr: f64,
    /// Centre of the nominating scan's mobility window (scan index), NaN without one.
    pub mobility: f64,
    pub ions: &'a [(f64, f32)],
}

/// Column-wise nominees as read; the ions of row `i` are `ion_mz[ion_off[i]..ion_off[i + 1]]`
/// (and the same range of `ion_int`).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct NomineeTable {
    pub spec_id: Vec<String>,
    pub decoy: Vec<bool>,
    pub charge: Vec<u8>,
    pub peptide: Vec<String>,
    pub proteins: Vec<String>,
    pub seq: Vec<String>,
    pub prec_mz: Vec<f64>,
    pub pred_rt: Vec<f64>,
    pub raw_pred_rt: Vec<f64>,
    pub rt_missing: Vec<bool>,
    pub apex_rt: Vec<f64>,
    pub dia_llr: Vec<f64>,
    pub mobility: Vec<f64>,
    pub ion_off: Vec<usize>,
    pub ion_mz: Vec<f64>,
    pub ion_int: Vec<f32>,
}

impl NomineeTable {
    pub fn len(&self) -> usize {
        self.decoy.len()
    }

    pub fn is_empty(&self) -> bool {
        self.decoy.is_empty()
    }

    /// The ions of row `i`, (m/z, predicted intensity).
    pub fn ions(&self, i: usize) -> impl Iterator<Item = (f64, f32)> + '_ {
        let r = self.ion_off[i]..self.ion_off[i + 1];
        self.ion_mz[r.clone()]
            .iter()
            .copied()
            .zip(self.ion_int[r].iter().copied())
    }
}

/// The precursor keys of a nominee table with each row's `dia_llr`, for study-level joins.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct NomineeKeys {
    pub decoy: Vec<bool>,
    pub charge: Vec<u8>,
    pub peptide: Vec<String>,
    pub dia_llr: Vec<f64>,
}

fn list(item: DataType) -> DataType {
    DataType::List(Arc::new(Field::new("item", item, true)))
}

fn schema() -> Arc<Schema> {
    Arc::new(Schema::new(vec![
        Field::new("spec_id", DataType::Utf8, false),
        Field::new("decoy", DataType::Boolean, false),
        Field::new("charge", DataType::UInt8, false),
        Field::new("peptide", DataType::Utf8, false),
        Field::new("proteins", DataType::Utf8, false),
        Field::new("seq", DataType::Utf8, false),
        Field::new("prec_mz", DataType::Float64, false),
        Field::new("pred_rt", DataType::Float64, false),
        Field::new("raw_pred_rt", DataType::Float64, false),
        Field::new("rt_missing", DataType::Boolean, false),
        Field::new("apex_rt", DataType::Float64, false),
        Field::new("dia_llr", DataType::Float64, false),
        Field::new("mobility", DataType::Float64, false),
        Field::new("ion_mz", list(DataType::Float64), false),
        Field::new("ion_int", list(DataType::Float32), false),
    ]))
}

/// Writes `n` nominees, row `i` given by `row(i)`, with extra key-value metadata.
pub fn write_nominees<'a>(
    path: &Path,
    n: usize,
    row: impl Fn(usize) -> NomineeRow<'a>,
    meta: &[(&str, String)],
) -> Result<(), String> {
    let sch = schema();
    let batches = (0..n.div_ceil(BATCH_ROWS).max(1)).map(|b| {
        let r = b * BATCH_ROWS..((b + 1) * BATCH_ROWS).min(n);
        let rows: Vec<NomineeRow> = r.map(&row).collect();
        let mut mz = ListBuilder::new(Float64Builder::new());
        let mut int = ListBuilder::new(Float32Builder::new());
        for x in &rows {
            for &(m, p) in x.ions {
                mz.values().append_value(m);
                int.values().append_value(p);
            }
            mz.append(true);
            int.append(true);
        }
        let cols: Vec<ArrayRef> = vec![
            Arc::new(StringArray::from_iter_values(rows.iter().map(|x| x.spec_id))),
            Arc::new(BooleanArray::from_iter(rows.iter().map(|x| Some(x.decoy)))),
            Arc::new(UInt8Array::from_iter_values(rows.iter().map(|x| x.charge))),
            Arc::new(StringArray::from_iter_values(rows.iter().map(|x| x.peptide))),
            Arc::new(StringArray::from_iter_values(rows.iter().map(|x| x.proteins))),
            Arc::new(StringArray::from_iter_values(rows.iter().map(|x| x.seq))),
            Arc::new(Float64Array::from_iter_values(rows.iter().map(|x| x.prec_mz))),
            Arc::new(Float64Array::from_iter_values(rows.iter().map(|x| x.pred_rt))),
            Arc::new(Float64Array::from_iter_values(rows.iter().map(|x| x.raw_pred_rt))),
            Arc::new(BooleanArray::from_iter(rows.iter().map(|x| Some(x.rt_missing)))),
            Arc::new(Float64Array::from_iter_values(rows.iter().map(|x| x.apex_rt))),
            Arc::new(Float64Array::from_iter_values(rows.iter().map(|x| x.dia_llr))),
            Arc::new(Float64Array::from_iter_values(rows.iter().map(|x| x.mobility))),
            Arc::new(mz.finish()),
            Arc::new(int.finish()),
        ];
        RecordBatch::try_new(sch.clone(), cols).map_err(|e| e.to_string())
    });
    write_table(path, &TABLE, SCHEMA_VERSION, sch.clone(), batches, meta)
}

fn open(path: &Path) -> Result<ParquetRecordBatchReaderBuilder<File>, String> {
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
    let major: u32 = version
        .split('.')
        .next()
        .and_then(|m| m.parse().ok())
        .unwrap_or(0);
    if major != DIA_NOMINEE_SCHEMA_MAJOR {
        return Err(format!(
            "{}: {} schema {version}, this andes reads {DIA_NOMINEE_SCHEMA_MAJOR}.x",
            path.display(),
            TABLE.what
        ));
    }
    Ok(builder)
}

/// The key-value metadata of a nominee table.
pub fn read_nominee_meta(path: &Path) -> Result<Vec<(String, String)>, String> {
    let builder = open(path)?;
    Ok(builder
        .metadata()
        .file_metadata()
        .key_value_metadata()
        .map(|kv| {
            kv.iter()
                .map(|k| (k.key.clone(), k.value.clone().unwrap_or_default()))
                .collect()
        })
        .unwrap_or_default())
}

fn projected(
    path: &Path,
    want: &[&str],
) -> Result<parquet::arrow::arrow_reader::ParquetRecordBatchReader, String> {
    let err = |e: &dyn std::fmt::Display| format!("{}: {e}", path.display());
    let builder = open(path)?;
    let fields = builder.schema().fields();
    let mut keep = Vec::with_capacity(want.len());
    for w in want {
        keep.push(
            fields
                .iter()
                .position(|f| f.name() == w)
                .ok_or_else(|| format!("{}: no column {w}", path.display()))?,
        );
    }
    let mask = ProjectionMask::roots(builder.parquet_schema(), keep.iter().copied());
    builder
        .with_projection(mask)
        .with_batch_size(BATCH_ROWS)
        .build()
        .map_err(|e| err(&e))
}

fn column<'b, A: 'static>(b: &'b RecordBatch, n: &str, path: &Path) -> Result<&'b A, String> {
    b.column_by_name(n)
        .and_then(|c| c.as_any().downcast_ref::<A>())
        .ok_or_else(|| format!("{}: column {n} is missing or has the wrong type", path.display()))
}

/// Reads only `decoy`, `charge`, `peptide` and `dia_llr` of a nominee table.
pub fn read_nominee_keys(path: &Path) -> Result<NomineeKeys, String> {
    let reader = projected(path, &["decoy", "charge", "peptide", "dia_llr"])?;
    let mut t = NomineeKeys::default();
    for b in reader {
        let b = b.map_err(|e| format!("{}: {e}", path.display()))?;
        let decoy = column::<BooleanArray>(&b, "decoy", path)?;
        let charge = column::<UInt8Array>(&b, "charge", path)?;
        let pep = column::<StringArray>(&b, "peptide", path)?;
        let llr = column::<Float64Array>(&b, "dia_llr", path)?;
        for i in 0..b.num_rows() {
            t.decoy.push(decoy.value(i));
            t.charge.push(charge.value(i));
            t.peptide.push(pep.value(i).to_string());
        }
        t.dia_llr.extend(llr.values().iter().copied());
    }
    Ok(t)
}

/// Reads a nominee table; with `keep` (sorted row indices), only those rows.
pub fn read_nominees(path: &Path, keep: Option<&[u32]>) -> Result<NomineeTable, String> {
    let names: Vec<String> = schema().fields().iter().map(|f| f.name().clone()).collect();
    let want: Vec<&str> = names.iter().map(String::as_str).collect();
    let reader = projected(path, &want)?;
    let mut t = NomineeTable {
        ion_off: vec![0],
        ..NomineeTable::default()
    };
    let mut row0 = 0usize;
    let mut k = 0usize;
    for b in reader {
        let b = b.map_err(|e| format!("{}: {e}", path.display()))?;
        let n = b.num_rows();
        let sid = column::<StringArray>(&b, "spec_id", path)?;
        let decoy = column::<BooleanArray>(&b, "decoy", path)?;
        let charge = column::<UInt8Array>(&b, "charge", path)?;
        let pep = column::<StringArray>(&b, "peptide", path)?;
        let prot = column::<StringArray>(&b, "proteins", path)?;
        let seq = column::<StringArray>(&b, "seq", path)?;
        let f64col = |c: &str| column::<Float64Array>(&b, c, path);
        let (mz, pred, apex, llr) = (f64col("prec_mz")?, f64col("pred_rt")?, f64col("apex_rt")?, f64col("dia_llr")?);
        let mobility = f64col("mobility")?;
        let raw = f64col("raw_pred_rt")?;
        let rt_missing = column::<BooleanArray>(&b, "rt_missing", path)?;
        let ion_mz = column::<ListArray>(&b, "ion_mz", path)?;
        let ion_int = column::<ListArray>(&b, "ion_int", path)?;
        let mut take = |i: usize| -> Result<(), String> {
            t.spec_id.push(sid.value(i).to_string());
            t.decoy.push(decoy.value(i));
            t.charge.push(charge.value(i));
            t.peptide.push(pep.value(i).to_string());
            t.proteins.push(prot.value(i).to_string());
            t.seq.push(seq.value(i).to_string());
            t.prec_mz.push(mz.value(i));
            t.pred_rt.push(pred.value(i));
            t.raw_pred_rt.push(raw.value(i));
            t.rt_missing.push(rt_missing.value(i));
            t.apex_rt.push(apex.value(i));
            t.dia_llr.push(llr.value(i));
            t.mobility.push(mobility.value(i));
            let m = ion_mz.value(i);
            let p = ion_int.value(i);
            let m = m
                .as_any()
                .downcast_ref::<Float64Array>()
                .ok_or_else(|| format!("{}: ion_mz items are not f64", path.display()))?;
            let p = p
                .as_any()
                .downcast_ref::<Float32Array>()
                .ok_or_else(|| format!("{}: ion_int items are not f32", path.display()))?;
            if m.len() != p.len() || m.null_count() + p.null_count() > 0 {
                return Err(format!("{}: row {} has unpaired ions", path.display(), row0 + i));
            }
            t.ion_mz.extend(m.values().iter().copied());
            t.ion_int.extend(p.values().iter().copied());
            t.ion_off.push(t.ion_mz.len());
            Ok(())
        };
        match keep {
            None => {
                for i in 0..n {
                    take(i)?;
                }
            }
            Some(rows) => {
                while k < rows.len() && (rows[k] as usize) < row0 + n {
                    take(rows[k] as usize - row0)?;
                    k += 1;
                }
            }
        }
        row0 += n;
    }
    if let Some(rows) = keep {
        if k != rows.len() {
            return Err(format!(
                "{}: {} of {} wanted rows found ({row0} rows)",
                path.display(),
                k,
                rows.len()
            ));
        }
    }
    Ok(t)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nominee_table_round_trips_and_selects_rows() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("run.parquet");
        let ions: Vec<Vec<(f64, f32)>> = vec![
            vec![(300.1, 1.0), (400.2, 0.5)],
            vec![],
            vec![(500.3, 0.25)],
        ];
        let peps = ["-.PEPTIDEK.-", "-.EPTPIDEK.-", "-.PEPTIDER.-"];
        write_nominees(
            &p,
            3,
            |i| NomineeRow {
                spec_id: ["a", "b", "c"][i],
                decoy: i == 1,
                charge: 2 + i as u8,
                peptide: peps[i],
                proteins: "P1",
                seq: &peps[i][2..10],
                prec_mz: 400.0 + i as f64,
                pred_rt: 10.0 + i as f64,
                raw_pred_rt: if i == 2 { f64::NAN } else { 0.3 + i as f64 },
                rt_missing: i == 2,
                apex_rt: 11.0 + i as f64,
                dia_llr: 20.0 - 10.0 * i as f64,
                mobility: if i == 0 { 702.5 } else { f64::NAN },
                ions: &ions[i],
            },
            &[("andes.test", "x".into())],
        )
        .unwrap();
        let t = read_nominees(&p, None).unwrap();
        assert_eq!(t.len(), 3);
        assert_eq!(t.decoy, vec![false, true, false]);
        assert_eq!(t.charge, vec![2, 3, 4]);
        assert_eq!(t.mobility[0], 702.5);
        assert_eq!(t.raw_pred_rt[..2], [0.3, 1.3]);
        assert!(t.raw_pred_rt[2].is_nan());
        assert_eq!(t.rt_missing, vec![false, false, true]);
        assert!(t.mobility[1].is_nan() && t.mobility[2].is_nan());
        assert_eq!(t.ions(0).collect::<Vec<_>>(), ions[0]);
        assert_eq!(t.ions(1).count(), 0);
        assert_eq!(t.ions(2).collect::<Vec<_>>(), ions[2]);
        let s = read_nominees(&p, Some(&[0, 2])).unwrap();
        assert_eq!(s.spec_id, vec!["a".to_string(), "c".to_string()]);
        assert_eq!(s.ions(1).collect::<Vec<_>>(), ions[2]);
        assert!(read_nominees(&p, Some(&[5])).is_err());
        let k = read_nominee_keys(&p).unwrap();
        assert_eq!(k.dia_llr, vec![20.0, 10.0, 0.0]);
        assert_eq!(k.peptide[2], "-.PEPTIDER.-");
        assert!(read_nominee_meta(&p)
            .unwrap()
            .iter()
            .any(|(k, v)| k == "andes.test" && v == "x"));
    }
}
