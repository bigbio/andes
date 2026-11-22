//! Checks of the default protein quantities of the QPX pg view against the feature view they
//! are rolled up from, shared by the finalize end-to-end tests.

#![allow(dead_code)]

use std::fs::File;
use std::path::Path;

use arrow::array::{Array, BooleanArray, Float32Array, Float64Array, ListArray, StringArray, StructArray};
use arrow::record_batch::RecordBatch;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

fn batches(path: &Path) -> Vec<RecordBatch> {
    let f = File::open(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    ParquetRecordBatchReaderBuilder::try_new(f)
        .unwrap()
        .build()
        .unwrap()
        .map(|b| b.unwrap())
        .collect()
}

fn col<'a, T: 'static>(b: &'a RecordBatch, name: &str) -> &'a T {
    b.column_by_name(name)
        .unwrap_or_else(|| panic!("no column {name}"))
        .as_any()
        .downcast_ref::<T>()
        .unwrap_or_else(|| panic!("column {name}: unexpected type"))
}

/// Totals of the two views: (sum of the pg rows' `intensity`, sum of the positive intensities
/// of the feature rows of unique peptides, quantified pg rows, labels seen on them). With
/// `gated`, feature rows whose `quant_weight` score is 0 are left out, as the roll-up does.
pub fn totals(qpx: &Path, gated: bool) -> (f64, f64, usize, Vec<String>) {
    let mut pg_total = 0.0;
    let mut quantified = 0;
    let mut labels: Vec<String> = Vec::new();
    for b in batches(&qpx.join("quantms.pg.parquet")) {
        let label = col::<StringArray>(&b, "label");
        let intensity = col::<Float32Array>(&b, "intensity");
        for i in 0..b.num_rows() {
            if label.is_null(i) {
                continue;
            }
            quantified += 1;
            if !labels.iter().any(|l| l == label.value(i)) {
                labels.push(label.value(i).to_string());
            }
            if !intensity.is_null(i) {
                pg_total += f64::from(intensity.value(i));
            }
        }
    }
    let mut feature_total = 0.0;
    for b in batches(&qpx.join("quantms.feature.parquet")) {
        let unique = col::<BooleanArray>(&b, "unique");
        let ints = col::<ListArray>(&b, "intensities");
        let scores = col::<ListArray>(&b, "additional_scores");
        for i in 0..b.num_rows() {
            if unique.is_null(i) || !unique.value(i) {
                continue;
            }
            if gated {
                let s = scores.value(i);
                let s = s.as_any().downcast_ref::<StructArray>().unwrap();
                let names = s.column(0).as_any().downcast_ref::<StringArray>().unwrap();
                let values = s.column(1).as_any().downcast_ref::<Float64Array>().unwrap();
                let w = (0..s.len())
                    .find(|&k| names.value(k) == "quant_weight")
                    .map(|k| values.value(k));
                if w.is_some_and(|w| w <= 0.0) {
                    continue;
                }
            }
            let e = ints.value(i);
            let e = e.as_any().downcast_ref::<StructArray>().unwrap();
            let v = e.column(1).as_any().downcast_ref::<Float32Array>().unwrap();
            for k in 0..v.len() {
                let x = f64::from(v.value(k));
                if x.is_finite() && x > 0.0 {
                    feature_total += x;
                }
            }
        }
    }
    (pg_total, feature_total, quantified, labels)
}
