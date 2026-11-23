//! Per-run DIA pass-2 table: one row per nominated precursor (targets and decoys) with its
//! pass-2 evidence columns and its quantities. Written by `andes --dia` and read by the
//! study-level steps.
//!
//! Schema 1.1 adds the quantity columns [`QUANT_COLUMNS`] (nullable; a 1.0 table reads as
//! having no quantities). They are kept apart from the evidence features, so they never enter
//! the precursor PIN or the learned score.

use std::fs::File;
use std::path::Path;
use std::sync::Arc;

use arrow::array::{
    Array, ArrayRef, BooleanArray, Float32Array, Float64Array, StringArray, UInt32Array, UInt8Array,
};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::ProjectionMask;

use crate::versioned::{write_table, TableVersion, BATCH_ROWS};

pub const DIA_PASS2_SCHEMA_MAJOR: u32 = 1;
const SCHEMA_VERSION: &str = "1.1";
const TABLE: TableVersion = TableVersion {
    key: "andes.dia_pass2.schema",
    version: SCHEMA_VERSION,
    major: DIA_PASS2_SCHEMA_MAJOR,
    what: "DIA pass-2 table",
};
const META: [&str; 8] = [
    "spec_id", "run", "decoy", "charge", "peptide", "proteins", "seq", "prec_mz",
];
/// Quantity columns of the table: the fragment quantity, the MS1 quantity and the number of
/// fragments the fragment quantity sums.
pub const QUANT_COLUMNS: [&str; 3] = ["quant_frag", "quant_ms1", "quant_n_frag"];

/// Column-wise pass-2 rows; `features` hold the evidence columns in a fixed order. The
/// quantity vectors are either empty (no quantities: written as nulls) or one value per row.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Pass2Table {
    pub spec_id: Vec<String>,
    pub run: Vec<u32>,
    pub decoy: Vec<bool>,
    pub charge: Vec<u8>,
    pub peptide: Vec<String>,
    pub proteins: Vec<String>,
    pub seq: Vec<String>,
    pub prec_mz: Vec<f64>,
    pub features: Vec<(String, Vec<f32>)>,
    /// Fragment quantity of the row in its run; `None` without a contained elution peak.
    pub quant_frag: Vec<Option<f64>>,
    /// MS1 quantity (monoisotopic and M+1 traces over the same window); `None` without signal.
    pub quant_ms1: Vec<Option<f64>>,
    /// Fragments summed into `quant_frag` (0 when it is `None`).
    pub quant_n_frag: Vec<u8>,
}

impl Pass2Table {
    pub fn len(&self) -> usize {
        self.run.len()
    }

    pub fn is_empty(&self) -> bool {
        self.run.is_empty()
    }

    pub fn feature(&self, name: &str) -> Option<&[f32]> {
        self.features
            .iter()
            .find(|(n, _)| n == name)
            .map(|(_, v)| v.as_slice())
    }
}

fn schema(features: &[(String, Vec<f32>)]) -> Arc<Schema> {
    let mut fields = vec![
        Field::new("spec_id", DataType::Utf8, false),
        Field::new("run", DataType::UInt32, false),
        Field::new("decoy", DataType::Boolean, false),
        Field::new("charge", DataType::UInt8, false),
        Field::new("peptide", DataType::Utf8, false),
        Field::new("proteins", DataType::Utf8, false),
        Field::new("seq", DataType::Utf8, false),
        Field::new("prec_mz", DataType::Float64, false),
    ];
    fields.extend(
        features
            .iter()
            .map(|(n, _)| Field::new(n.as_str(), DataType::Float32, false)),
    );
    fields.extend([
        Field::new(QUANT_COLUMNS[0], DataType::Float64, true),
        Field::new(QUANT_COLUMNS[1], DataType::Float64, true),
        Field::new(QUANT_COLUMNS[2], DataType::UInt8, false),
    ]);
    Arc::new(Schema::new(fields))
}

pub fn write_pass2(path: &Path, t: &Pass2Table) -> Result<(), String> {
    write_pass2_with(path, t, &[])
}

/// Writes a pass-2 table with extra key-value metadata (e.g. the rescoring settings of a run).
pub fn write_pass2_with(path: &Path, t: &Pass2Table, meta: &[(&str, String)]) -> Result<(), String> {
    let sch = schema(&t.features);
    let n = t.len();
    for (name, len) in [
        (QUANT_COLUMNS[0], t.quant_frag.len()),
        (QUANT_COLUMNS[1], t.quant_ms1.len()),
        (QUANT_COLUMNS[2], t.quant_n_frag.len()),
    ] {
        if len != 0 && len != n {
            return Err(format!(
                "{}: {name} has {len} values for {n} rows",
                path.display()
            ));
        }
    }
    let opt = |v: &[Option<f64>], r: std::ops::Range<usize>| -> ArrayRef {
        if v.is_empty() {
            Arc::new(Float64Array::from(vec![None; r.len()]))
        } else {
            Arc::new(Float64Array::from(v[r].to_vec()))
        }
    };
    let batches = (0..n.div_ceil(BATCH_ROWS).max(1)).map(|b| {
        let r = b * BATCH_ROWS..((b + 1) * BATCH_ROWS).min(n);
        let mut cols: Vec<ArrayRef> = vec![
            Arc::new(StringArray::from_iter_values(t.spec_id[r.clone()].iter())),
            Arc::new(UInt32Array::from(t.run[r.clone()].to_vec())),
            Arc::new(BooleanArray::from(t.decoy[r.clone()].to_vec())),
            Arc::new(UInt8Array::from(t.charge[r.clone()].to_vec())),
            Arc::new(StringArray::from_iter_values(t.peptide[r.clone()].iter())),
            Arc::new(StringArray::from_iter_values(t.proteins[r.clone()].iter())),
            Arc::new(StringArray::from_iter_values(t.seq[r.clone()].iter())),
            Arc::new(Float64Array::from(t.prec_mz[r.clone()].to_vec())),
        ];
        cols.extend(
            t.features
                .iter()
                .map(|(_, v)| Arc::new(Float32Array::from(v[r.clone()].to_vec())) as ArrayRef),
        );
        cols.push(opt(&t.quant_frag, r.clone()));
        cols.push(opt(&t.quant_ms1, r.clone()));
        cols.push(Arc::new(UInt8Array::from(if t.quant_n_frag.is_empty() {
            vec![0; r.len()]
        } else {
            t.quant_n_frag[r.clone()].to_vec()
        })));
        RecordBatch::try_new(sch.clone(), cols).map_err(|e| e.to_string())
    });
    write_table(path, &TABLE, SCHEMA_VERSION, sch.clone(), batches, meta)
}

/// Opens a pass-2 table and checks its schema version.
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
    if major != DIA_PASS2_SCHEMA_MAJOR {
        return Err(format!(
            "{}: {} schema {version}, this andes reads {DIA_PASS2_SCHEMA_MAJOR}.x",
            path.display(),
            TABLE.what
        ));
    }
    Ok(builder)
}

/// The key-value metadata of a pass-2 table.
pub fn read_pass2_meta(path: &Path) -> Result<Vec<(String, String)>, String> {
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

/// The precursor keys of a pass-2 table with two of its features, for study-level joins.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Pass2Keys {
    pub decoy: Vec<bool>,
    pub charge: Vec<u8>,
    pub peptide: Vec<String>,
    pub features: [Vec<f32>; 2],
}

/// Reads only `decoy`, `charge`, `peptide` and the two named features of a pass-2 table.
pub fn read_pass2_keys(path: &Path, features: [&str; 2]) -> Result<Pass2Keys, String> {
    let err = |e: &dyn std::fmt::Display| format!("{}: {e}", path.display());
    let builder = open(path)?;
    let want = ["decoy", "charge", "peptide", features[0], features[1]];
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
    let reader = builder
        .with_projection(mask)
        .with_batch_size(BATCH_ROWS)
        .build()
        .map_err(|e| err(&e))?;
    let mut t = Pass2Keys::default();
    for b in reader {
        let b = b.map_err(|e| err(&e))?;
        let bad = |n: &str| format!("{}: column {n} is missing or has the wrong type", path.display());
        let get = |n: &str| b.column_by_name(n).map(|c| c.as_any());
        let decoy = get("decoy").and_then(|c| c.downcast_ref::<BooleanArray>()).ok_or_else(|| bad("decoy"))?;
        let charge = get("charge").and_then(|c| c.downcast_ref::<UInt8Array>()).ok_or_else(|| bad("charge"))?;
        let pep = get("peptide").and_then(|c| c.downcast_ref::<StringArray>()).ok_or_else(|| bad("peptide"))?;
        for i in 0..b.num_rows() {
            t.decoy.push(decoy.value(i));
            t.charge.push(charge.value(i));
            t.peptide.push(pep.value(i).to_string());
        }
        for (k, n) in features.iter().enumerate() {
            let a = get(n).and_then(|c| c.downcast_ref::<Float32Array>()).ok_or_else(|| bad(n))?;
            t.features[k].extend(a.values().iter().copied());
        }
    }
    Ok(t)
}

/// Reads a pass-2 table; with `only`, the metadata columns and just the named features.
pub fn read_pass2(path: &Path, only: Option<&[&str]>) -> Result<Pass2Table, String> {
    let err = |e: &dyn std::fmt::Display| format!("{}: {e}", path.display());
    let builder = open(path)?;
    let names: Vec<String> = builder
        .schema()
        .fields()
        .iter()
        .map(|f| f.name().clone())
        .collect();
    let keep: Vec<usize> = names
        .iter()
        .enumerate()
        .filter(|(_, n)| {
            META.contains(&n.as_str())
                || QUANT_COLUMNS.contains(&n.as_str())
                || only.is_none_or(|o| o.contains(&n.as_str()))
        })
        .map(|(i, _)| i)
        .collect();
    let mask = ProjectionMask::roots(builder.parquet_schema(), keep.iter().copied());
    let reader = builder
        .with_projection(mask)
        .with_batch_size(BATCH_ROWS)
        .build()
        .map_err(|e| err(&e))?;
    let mut t = Pass2Table::default();
    let feat_names: Vec<String> = keep
        .iter()
        .map(|&i| names[i].clone())
        .filter(|n| !META.contains(&n.as_str()) && !QUANT_COLUMNS.contains(&n.as_str()))
        .collect();
    t.features = feat_names.iter().map(|n| (n.clone(), Vec::new())).collect();
    for b in reader {
        let b = b.map_err(|e| err(&e))?;
        let col = |n: &str| -> Result<&ArrayRef, String> {
            b.column_by_name(n)
                .ok_or_else(|| format!("{}: no column {n}", path.display()))
        };
        let down = |n: &str| col(n).map(|c| c.as_any());
        let bad = |n: &str| format!("{}: column {n} has the wrong type", path.display());
        let sid = down("spec_id")?.downcast_ref::<StringArray>().ok_or_else(|| bad("spec_id"))?;
        let run = down("run")?.downcast_ref::<UInt32Array>().ok_or_else(|| bad("run"))?;
        let decoy = down("decoy")?.downcast_ref::<BooleanArray>().ok_or_else(|| bad("decoy"))?;
        let charge = down("charge")?.downcast_ref::<UInt8Array>().ok_or_else(|| bad("charge"))?;
        let pep = down("peptide")?.downcast_ref::<StringArray>().ok_or_else(|| bad("peptide"))?;
        let prot = down("proteins")?.downcast_ref::<StringArray>().ok_or_else(|| bad("proteins"))?;
        let seq = down("seq")?.downcast_ref::<StringArray>().ok_or_else(|| bad("seq"))?;
        let mz = down("prec_mz")?.downcast_ref::<Float64Array>().ok_or_else(|| bad("prec_mz"))?;
        for i in 0..b.num_rows() {
            t.spec_id.push(sid.value(i).to_string());
            t.run.push(run.value(i));
            t.decoy.push(decoy.value(i));
            t.charge.push(charge.value(i));
            t.peptide.push(pep.value(i).to_string());
            t.proteins.push(prot.value(i).to_string());
            t.seq.push(seq.value(i).to_string());
            t.prec_mz.push(mz.value(i));
        }
        for (n, v) in t.features.iter_mut() {
            let a = down(n)?.downcast_ref::<Float32Array>().ok_or_else(|| bad(n))?;
            v.extend(a.values().iter().copied());
            debug_assert_eq!(a.null_count(), 0);
        }
        for (n, v) in [
            (QUANT_COLUMNS[0], &mut t.quant_frag),
            (QUANT_COLUMNS[1], &mut t.quant_ms1),
        ] {
            match b.column_by_name(n) {
                Some(c) => {
                    let a = c.as_any().downcast_ref::<Float64Array>().ok_or_else(|| bad(n))?;
                    v.extend((0..a.len()).map(|i| a.is_valid(i).then(|| a.value(i))));
                }
                None => v.extend(std::iter::repeat_n(None, b.num_rows())),
            }
        }
        match b.column_by_name(QUANT_COLUMNS[2]) {
            Some(c) => {
                let a = c
                    .as_any()
                    .downcast_ref::<UInt8Array>()
                    .ok_or_else(|| bad(QUANT_COLUMNS[2]))?;
                t.quant_n_frag.extend(a.values().iter().copied());
            }
            None => t.quant_n_frag.extend(std::iter::repeat_n(0, b.num_rows())),
        }
    }
    Ok(t)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn table() -> Pass2Table {
        Pass2Table {
            spec_id: vec!["a_1_1".into(), "a_1_2".into(), "b_7_1".into()],
            run: vec![0, 0, 1],
            decoy: vec![false, true, false],
            charge: vec![2, 3, 2],
            peptide: vec!["-.PEPTIDEK.-".into(), "-.EPTPIDEK.-".into(), "-.PEPTIDEK.-".into()],
            proteins: vec!["P1".into(), "XXX_P1".into(), "P1".into()],
            seq: vec!["PEPTIDEK".into(), "EPTPIDEK".into(), "PEPTIDEK".into()],
            prec_mz: vec![464.7, 310.1, 464.7],
            features: vec![
                ("p2_llr".into(), vec![12.5, -3.0, 9.0]),
                ("p2_rtstar".into(), vec![10.1, 11.0, 10.2]),
            ],
            quant_frag: vec![Some(1.5e6), None, Some(2.0e5)],
            quant_ms1: vec![Some(3.0e7), Some(1.0e4), None],
            quant_n_frag: vec![6, 0, 3],
        }
    }

    #[test]
    fn pass2_table_round_trips() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("run.dia_pass2.parquet");
        let t = table();
        write_pass2(&p, &t).unwrap();
        assert_eq!(read_pass2(&p, None).unwrap(), t);
    }

    #[test]
    fn pass2_keys_and_metadata_round_trip() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("run.parquet");
        write_pass2_with(&p, &table(), &[("andes.test", "a b".to_string())]).unwrap();
        let k = read_pass2_keys(&p, ["p2_llr", "p2_rtstar"]).unwrap();
        assert_eq!(k.peptide, table().peptide);
        assert_eq!(k.decoy, vec![false, true, false]);
        assert_eq!(k.charge, vec![2, 3, 2]);
        assert_eq!(k.features[0], vec![12.5, -3.0, 9.0]);
        assert_eq!(k.features[1], vec![10.1, 11.0, 10.2]);
        let meta = read_pass2_meta(&p).unwrap();
        assert!(meta.iter().any(|(k, v)| k == "andes.test" && v == "a b"));
        assert!(read_pass2_keys(&p, ["p2_llr", "p2_none"]).is_err());
    }

    #[test]
    fn pass2_projection_keeps_named_features() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("run.dia_pass2.parquet");
        write_pass2(&p, &table()).unwrap();
        let r = read_pass2(&p, Some(&["p2_rtstar"])).unwrap();
        assert_eq!(r.features.len(), 1);
        assert_eq!(r.feature("p2_rtstar").unwrap(), &[10.1, 11.0, 10.2]);
        assert_eq!(r.len(), 3);
        assert_eq!(r.quant_frag, table().quant_frag, "quantities are read with any projection");
    }

    #[test]
    fn quantities_stay_out_of_the_features_and_default_to_missing() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("run.dia_pass2.parquet");
        let mut t = table();
        t.quant_frag.clear();
        t.quant_ms1.clear();
        t.quant_n_frag.clear();
        write_pass2(&p, &t).unwrap();
        let r = read_pass2(&p, None).unwrap();
        assert_eq!(r.features, t.features);
        assert_eq!(r.quant_frag, vec![None; 3]);
        assert_eq!(r.quant_ms1, vec![None; 3]);
        assert_eq!(r.quant_n_frag, vec![0; 3]);
        t.quant_frag = vec![Some(1.0)];
        assert!(write_pass2(&p, &t).is_err(), "a quantity column of the wrong length");
    }
}
