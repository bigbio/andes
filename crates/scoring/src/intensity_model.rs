//! Context-conditioned relative fragment intensity model (strong-score numerator).
//!
//! Loaded from `intensity_model.parquet` produced by `msnet_intensity_agg.py` /
//! `andes train-intensity`. Lookup backs off sparse keys: drop `nce_bin`, then
//! flanking residues.

use std::fs::File;
use std::path::Path;

use arrow::array::{Array, Float64Array, Int32Array, Int64Array, StringArray};
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use rustc_hash::FxHashMap;
use thiserror::Error;

/// Minimum observations before trusting a fine-grained key (matches T1 backoff).
pub const DEFAULT_MIN_COUNT: u64 = 30;

/// Ion type for intensity lookup (plain b/y only in Phase T).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum IntensityIonType {
    B,
    Y,
}

impl IntensityIonType {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::B => "b",
            Self::Y => "y",
        }
    }

    pub fn parse(s: &str) -> Option<Self> {
        match s {
            "b" => Some(Self::B),
            "y" => Some(Self::Y),
            _ => None,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct IntensityKey {
    ion_type: IntensityIonType,
    flank_n: u8,
    flank_c: u8,
    pos_bin: i32,
    charge: i32,
    nce_bin: String,
}

#[derive(Debug, Clone, Copy)]
struct Stats {
    count: u64,
    mean_log_rel: f64,
    var_log_rel: f64,
}

impl Stats {
    fn spread(self) -> f64 {
        self.var_log_rel.max(0.0).sqrt()
    }
}

#[derive(Debug, Error)]
pub enum IntensityModelError {
    #[error("parquet read: {0}")]
    Parquet(String),
    #[error("missing column {0}")]
    MissingColumn(&'static str),
    #[error("empty intensity model")]
    Empty,
}

/// Relative-intensity context model with sparse-key backoff.
#[derive(Debug, Clone)]
pub struct IntensityModel {
    entries: FxHashMap<IntensityKey, Stats>,
    min_count: u64,
    /// Global fallback when all backoff keys are sparse.
    global: Stats,
}

impl IntensityModel {
    pub fn min_count(&self) -> u64 {
        self.min_count
    }

    /// Load a finalized `intensity_model.parquet` (mean/var columns).
    pub fn load(path: &Path) -> Result<Self, IntensityModelError> {
        Self::load_with_min_count(path, DEFAULT_MIN_COUNT)
    }

    pub fn load_with_min_count(path: &Path, min_count: u64) -> Result<Self, IntensityModelError> {
        let file = File::open(path)
            .map_err(|e| IntensityModelError::Parquet(format!("open {}: {e}", path.display())))?;
        let builder = ParquetRecordBatchReaderBuilder::try_new(file)
            .map_err(|e| IntensityModelError::Parquet(format!("reader: {e}")))?;
        let mut entries = FxHashMap::default();
        let mut total_count: u64 = 0;
        let mut total_mean_sum = 0.0;

        for batch_result in builder.build().map_err(|e| IntensityModelError::Parquet(e.to_string()))? {
            let batch = batch_result.map_err(|e| IntensityModelError::Parquet(e.to_string()))?;
            let ion_col = str_col(&batch, "ion_type")?;
            let flank_n_col = str_col(&batch, "flank_n")?;
            let flank_c_col = str_col(&batch, "flank_c")?;
            let pos_col = batch
                .column_by_name("pos_bin")
                .and_then(|c| c.as_any().downcast_ref::<Int32Array>())
                .ok_or(IntensityModelError::MissingColumn("pos_bin"))?;
            let charge_col = batch
                .column_by_name("charge")
                .and_then(|c| c.as_any().downcast_ref::<Int32Array>())
                .ok_or(IntensityModelError::MissingColumn("charge"))?;
            let nce_col = str_col(&batch, "nce_bin")?;
            let count_col = batch
                .column_by_name("count")
                .and_then(|c| c.as_any().downcast_ref::<Int64Array>())
                .ok_or(IntensityModelError::MissingColumn("count"))?;
            let mean_col = batch
                .column_by_name("mean_log_rel")
                .and_then(|c| c.as_any().downcast_ref::<Float64Array>())
                .ok_or(IntensityModelError::MissingColumn("mean_log_rel"))?;
            let var_col = batch
                .column_by_name("var_log_rel")
                .and_then(|c| c.as_any().downcast_ref::<Float64Array>())
                .ok_or(IntensityModelError::MissingColumn("var_log_rel"))?;

            for i in 0..batch.num_rows() {
                if ion_col.is_null(i) || flank_n_col.is_null(i) || flank_c_col.is_null(i) {
                    continue;
                }
                // Skip rows with null numeric/position fields: arrow `.value(i)`
                // on a null returns a default (0), which would silently corrupt
                // the key/stats. A well-formed `train-intensity` parquet has no
                // nulls here, but be defensive against externally-produced files.
                if pos_col.is_null(i)
                    || charge_col.is_null(i)
                    || count_col.is_null(i)
                    || mean_col.is_null(i)
                    || var_col.is_null(i)
                {
                    continue;
                }
                let ion_s = ion_col.value(i);
                let Some(ion_type) = IntensityIonType::parse(ion_s) else {
                    continue;
                };
                let flank_n = flank_n_col.value(i).as_bytes().first().copied().unwrap_or(b'X');
                let flank_c = flank_c_col.value(i).as_bytes().first().copied().unwrap_or(b'X');
                let count = count_col.value(i) as u64;
                let mean = mean_col.value(i);
                let var = var_col.value(i);
                let key = IntensityKey {
                    ion_type,
                    flank_n,
                    flank_c,
                    pos_bin: pos_col.value(i),
                    charge: charge_col.value(i),
                    nce_bin: if nce_col.is_null(i) {
                        "unknown".to_string()
                    } else {
                        nce_col.value(i).to_string()
                    },
                };
                entries.insert(
                    key,
                    Stats {
                        count,
                        mean_log_rel: mean,
                        var_log_rel: var,
                    },
                );
                total_count = total_count.saturating_add(count);
                total_mean_sum += mean * count as f64;
            }
        }

        if entries.is_empty() {
            return Err(IntensityModelError::Empty);
        }

        let global_mean = if total_count > 0 {
            total_mean_sum / total_count as f64
        } else {
            -2.0
        };
        let global = Stats {
            count: total_count,
            mean_log_rel: global_mean,
            var_log_rel: 1.0,
        };

        Ok(Self {
            entries,
            min_count,
            global,
        })
    }

    /// Predict log relative intensity (mean, spread) with sparse-key backoff.
    pub fn predict_log_rel(
        &self,
        ion_type: IntensityIonType,
        flank_n: u8,
        flank_c: u8,
        pos_bin: i32,
        charge: i32,
        nce_bin: &str,
    ) -> (f64, f64) {
        let candidates = backoff_keys(ion_type, flank_n, flank_c, pos_bin, charge, nce_bin);
        for key in &candidates {
            if let Some(stats) = self.entries.get(key) {
                if stats.count >= self.min_count {
                    return (stats.mean_log_rel, stats.spread());
                }
            }
        }
        (self.global.mean_log_rel, self.global.spread())
    }

    /// Number of stored context keys (for tests / diagnostics).
    pub fn key_count(&self) -> usize {
        self.entries.len()
    }
}

fn str_col<'a>(
    batch: &'a arrow::record_batch::RecordBatch,
    name: &'static str,
) -> Result<&'a StringArray, IntensityModelError> {
    batch
        .column_by_name(name)
        .and_then(|c| c.as_any().downcast_ref::<StringArray>())
        .ok_or(IntensityModelError::MissingColumn(name))
}

fn backoff_keys(
    ion_type: IntensityIonType,
    flank_n: u8,
    flank_c: u8,
    pos_bin: i32,
    charge: i32,
    nce_bin: &str,
) -> Vec<IntensityKey> {
    let any_nce = IntensityKey {
        ion_type,
        flank_n,
        flank_c,
        pos_bin,
        charge,
        nce_bin: "__any__".to_string(),
    };
    let any_flank = IntensityKey {
        ion_type,
        flank_n: b'*',
        flank_c: b'*',
        pos_bin,
        charge,
        nce_bin: nce_bin.to_string(),
    };
    let any_both = IntensityKey {
        ion_type,
        flank_n: b'*',
        flank_c: b'*',
        pos_bin,
        charge,
        nce_bin: "__any__".to_string(),
    };
    vec![
        IntensityKey {
            ion_type,
            flank_n,
            flank_c,
            pos_bin,
            charge,
            nce_bin: nce_bin.to_string(),
        },
        any_nce,
        any_flank,
        any_both,
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use arrow::array::{Float64Array, Int32Array, Int64Array, StringArray};
    use arrow::datatypes::{DataType, Field, Schema};
    use arrow::record_batch::RecordBatch;
    use parquet::arrow::ArrowWriter;
    use tempfile::NamedTempFile;

    /// (ion_type, flank_n, flank_c, pos_bin, charge, nce_bin, count, mean, var)
    type FixtureRow<'a> = (&'a str, &'a str, &'a str, i32, i32, &'a str, i64, f64, f64);

    fn write_fixture(path: &Path, rows: &[FixtureRow]) {
        let schema = Schema::new(vec![
            Field::new("ion_type", DataType::Utf8, false),
            Field::new("flank_n", DataType::Utf8, false),
            Field::new("flank_c", DataType::Utf8, false),
            Field::new("pos_bin", DataType::Int32, false),
            Field::new("charge", DataType::Int32, false),
            Field::new("nce_bin", DataType::Utf8, false),
            Field::new("count", DataType::Int64, false),
            Field::new("mean_log_rel", DataType::Float64, false),
            Field::new("var_log_rel", DataType::Float64, false),
        ]);
        let ion: Vec<_> = rows.iter().map(|r| r.0).collect();
        let fn_: Vec<_> = rows.iter().map(|r| r.1).collect();
        let fc: Vec<_> = rows.iter().map(|r| r.2).collect();
        let pb: Vec<_> = rows.iter().map(|r| r.3).collect();
        let ch: Vec<_> = rows.iter().map(|r| r.4).collect();
        let nce: Vec<_> = rows.iter().map(|r| r.5).collect();
        let cnt: Vec<_> = rows.iter().map(|r| r.6).collect();
        let mean: Vec<_> = rows.iter().map(|r| r.7).collect();
        let var: Vec<_> = rows.iter().map(|r| r.8).collect();
        let batch = RecordBatch::try_new(
            std::sync::Arc::new(schema),
            vec![
                std::sync::Arc::new(StringArray::from(ion)),
                std::sync::Arc::new(StringArray::from(fn_)),
                std::sync::Arc::new(StringArray::from(fc)),
                std::sync::Arc::new(Int32Array::from(pb)),
                std::sync::Arc::new(Int32Array::from(ch)),
                std::sync::Arc::new(StringArray::from(nce)),
                std::sync::Arc::new(Int64Array::from(cnt)),
                std::sync::Arc::new(Float64Array::from(mean)),
                std::sync::Arc::new(Float64Array::from(var)),
            ],
        )
        .unwrap();
        let file = File::create(path).unwrap();
        let mut writer = ArrowWriter::try_new(file, batch.schema(), None).unwrap();
        writer.write(&batch).unwrap();
        writer.close().unwrap();
    }

    #[test]
    fn round_trip_and_known_key_lookup() {
        let tmp = NamedTempFile::new().unwrap();
        write_fixture(
            tmp.path(),
            &[
                ("y", "K", "R", 5, 2, "25", 100, -0.5, 0.2),
                ("b", "A", "L", 1, 2, "25", 5, -3.0, 0.5),
                ("y", "*", "*", 5, 2, "__any__", 500, -1.0, 0.4),
            ],
        );
        let model = IntensityModel::load(tmp.path()).unwrap();
        assert_eq!(model.key_count(), 3);
        let (mean, spread) = model.predict_log_rel(IntensityIonType::Y, b'K', b'R', 5, 2, "25");
        assert!((mean - (-0.5)).abs() < 1e-9);
        assert!((spread - 0.2_f64.sqrt()).abs() < 1e-9);
        // Sparse b-ion (count=5) backs off past the fine key to the global mean.
        let (mean_b, _) =
            model.predict_log_rel(IntensityIonType::B, b'A', b'L', 1, 2, "25");
        let global = (-0.5 * 100.0 - 3.0 * 5.0 - 1.0 * 500.0) / 605.0;
        assert!((mean_b - global).abs() < 1e-9);
        assert!((mean_b - (-3.0)).abs() > 0.1);
    }

    #[test]
    fn backoff_drops_nce_when_fine_key_sparse() {
        let tmp = NamedTempFile::new().unwrap();
        write_fixture(
            tmp.path(),
            &[
                ("y", "K", "R", 5, 2, "unknown", 2, -0.1, 0.1),
                ("y", "K", "R", 5, 2, "__any__", 80, -0.8, 0.3),
            ],
        );
        let model = IntensityModel::load(tmp.path()).unwrap();
        let (mean, _) = model.predict_log_rel(IntensityIonType::Y, b'K', b'R', 5, 2, "25");
        assert!((mean - (-0.8)).abs() < 1e-9);
    }
}
