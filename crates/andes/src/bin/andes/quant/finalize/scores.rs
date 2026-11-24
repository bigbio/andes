//! Temporary Parquet files for transfer scores and the score cache reused between passes.

use std::cell::RefCell;
use std::fs::File;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

use arrow::array::{Array, ArrayRef, Float32Array};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::arrow_reader::{ParquetRecordBatchReader, ParquetRecordBatchReaderBuilder};
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;
use parquet::format::KeyValue;
use rayon::prelude::*;

use super::index::{EntryIndex, Table};
use super::SCORE_GROUP_ROWS;
use crate::quant::mbr::ExtractRow;
use crate::transfers::{FoldScores, QValues};

const SCORES_KEY: &str = "andes.finalize_scores.schema";
pub(super) const ROWS_KEY: &str = "andes.finalize_rows.schema";
const TEMP_VERSION: &str = "1.0";
const TEMP_MAJOR: &str = "1";

pub(super) fn temp_writer(
    path: &Path,
    key: &str,
    schema: Arc<Schema>,
    compression: Compression,
    group_rows: usize,
) -> Result<ArrowWriter<File>, String> {
    let err = |e: &dyn std::fmt::Display| format!("{}: {e}", path.display());
    let props = WriterProperties::builder()
        .set_compression(compression)
        .set_max_row_group_size(group_rows)
        .set_key_value_metadata(Some(vec![KeyValue::new(
            key.to_string(),
            TEMP_VERSION.to_string(),
        )]))
        .build();
    let file = File::create(path).map_err(|e| err(&e))?;
    ArrowWriter::try_new(file, schema, Some(props)).map_err(|e| err(&e))
}

pub(super) fn temp_reader(
    path: &Path,
    key: &str,
    batch_rows: usize,
) -> Result<ParquetRecordBatchReader, String> {
    let err = |e: &dyn std::fmt::Display| format!("{}: {e}", path.display());
    let file = File::open(path).map_err(|e| err(&e))?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file).map_err(|e| err(&e))?;
    let version = builder
        .metadata()
        .file_metadata()
        .key_value_metadata()
        .and_then(|kv| kv.iter().find(|k| k.key == key))
        .and_then(|k| k.value.clone());
    if version.as_deref().and_then(|v| v.split('.').next()) != Some(TEMP_MAJOR) {
        return Err(err(&format!(
            "not a temporary table of schema {key} {TEMP_MAJOR}.x"
        )));
    }
    builder
        .with_batch_size(batch_rows)
        .build()
        .map_err(|e| err(&e))
}

pub(super) fn column<'a, T: 'static>(
    b: &'a RecordBatch,
    i: usize,
    path: &Path,
) -> Result<&'a T, String> {
    b.column(i)
        .as_any()
        .downcast_ref::<T>()
        .ok_or_else(|| format!("{}: column {i} has an unexpected type", path.display()))
}

pub(super) fn write_scores(path: &Path, scores: &[f32]) -> Result<(), String> {
    let err = |e: &dyn std::fmt::Display| format!("{}: {e}", path.display());
    let schema = Arc::new(Schema::new(vec![Field::new(
        "score",
        DataType::Float32,
        false,
    )]));
    let mut w = temp_writer(
        path,
        SCORES_KEY,
        schema.clone(),
        Compression::UNCOMPRESSED,
        SCORE_GROUP_ROWS,
    )?;
    for c in scores.chunks(SCORE_GROUP_ROWS) {
        let col: ArrayRef = Arc::new(Float32Array::from(c.to_vec()));
        let b = RecordBatch::try_new(schema.clone(), vec![col]).map_err(|e| err(&e))?;
        w.write(&b).map_err(|e| err(&e))?;
    }
    w.close().map_err(|e| err(&e))?;
    Ok(())
}

pub(super) fn is_cached(idx: &EntryIndex<'_>, r: &ExtractRow) -> bool {
    matches!(idx.table(r), Some((Table::Transfer | Table::Identified, _)))
}

pub(super) const NO_SCORE: f32 = f32::INFINITY;

pub(super) fn cached_scores_path(dir: &Path, chunk: usize) -> PathBuf {
    dir.join(format!("transfer-scores-{chunk:06}.parquet"))
}

pub(super) fn read_cached_scores(path: &Path) -> Result<Vec<f32>, String> {
    let error = RefCell::new(None);
    let scores: Vec<f32> = ScoreFile::new(path, &error).collect();
    match error.into_inner() {
        Some(e) => Err(e),
        None => Ok(scores),
    }
}

pub(super) struct ScoreFile<'a> {
    pub(super) path: PathBuf,
    pub(super) reader: Option<ParquetRecordBatchReader>,
    batch: Option<Float32Array>,
    pos: usize,
    error: &'a RefCell<Option<String>>,
}

impl<'a> ScoreFile<'a> {
    pub(super) fn new(path: &Path, error: &'a RefCell<Option<String>>) -> Self {
        let reader = match temp_reader(path, SCORES_KEY, SCORE_GROUP_ROWS) {
            Ok(r) => Some(r),
            Err(e) => {
                error.borrow_mut().get_or_insert(e);
                None
            }
        };
        Self {
            path: path.to_path_buf(),
            reader,
            batch: None,
            pos: 0,
            error,
        }
    }

    fn fail(&mut self, e: String) -> Option<f32> {
        self.error.borrow_mut().get_or_insert(e);
        self.reader = None;
        None
    }
}

impl Iterator for ScoreFile<'_> {
    type Item = f32;

    fn next(&mut self) -> Option<f32> {
        loop {
            if let Some(b) = &self.batch {
                if self.pos < b.len() {
                    self.pos += 1;
                    return Some(b.value(self.pos - 1));
                }
            }
            let next = self.reader.as_mut()?.next();
            match next {
                None => {
                    self.reader = None;
                    return None;
                }
                Some(Err(e)) => return self.fail(format!("{}: {e}", self.path.display())),
                Some(Ok(b)) => match column::<Float32Array>(&b, 0, &self.path) {
                    Ok(a) => {
                        self.batch = Some(a.clone());
                        self.pos = 0;
                    }
                    Err(e) => return self.fail(e),
                },
            }
        }
    }
}

pub(super) struct ScoreStore<'d> {
    pub(super) name: &'static str,
    pub(super) dir: &'d Path,
    cap: usize,
    pub(super) scores: Mutex<FoldScores>,
    spills: Mutex<Vec<Vec<[PathBuf; 4]>>>,
    pub(super) next: AtomicUsize,
}

impl<'d> ScoreStore<'d> {
    pub(super) fn new(
        name: &'static str,
        dir: &'d Path,
        cap: usize,
        seed: u64,
        strata: usize,
    ) -> Self {
        Self {
            name,
            dir,
            cap,
            scores: Mutex::new(FoldScores::with_strata(seed, strata)),
            spills: Mutex::new(Vec::new()),
            next: AtomicUsize::new(0),
        }
    }

    pub(super) fn add(&self, s: FoldScores) -> Result<(), String> {
        let held = {
            let mut all = self.scores.lock().expect("score store lock");
            all.merge(s);
            (all.held() > self.cap).then(|| all.take_held())
        };
        let Some(held) = held else { return Ok(()) };
        let k = self.next.fetch_add(1, Ordering::Relaxed);
        let mut paths = Vec::with_capacity(held.len());
        for (fold, by_family) in held.into_iter().enumerate() {
            let mut fold_paths: [PathBuf; 4] = Default::default();
            for (f, mut v) in by_family.into_iter().enumerate() {
                v.par_sort_unstable_by(f32::total_cmp);
                let p = self
                    .dir
                    .join(format!("scores-{}-{k:04}-{fold}-{f}.parquet", self.name));
                write_scores(&p, &v)?;
                fold_paths[f] = p;
            }
            paths.push(fold_paths);
        }
        self.spills.lock().expect("score store lock").push(paths);
        Ok(())
    }

    pub(super) fn finish(
        self,
        failed_folds: &[String],
        thresholds: &[f64],
    ) -> Result<(QValues, usize), String> {
        let spills = self.spills.into_inner().expect("score store lock");
        let scores = self.scores.into_inner().expect("score store lock");
        let error = RefCell::new(None);
        let q = scores.finish_with(failed_folds, Some(thresholds), |fold, f| {
            spills
                .iter()
                .map(|s| {
                    Box::new(ScoreFile::new(&s[fold][f], &error))
                        as Box<dyn Iterator<Item = f32> + '_>
                })
                .collect()
        });
        if let Some(e) = error.into_inner() {
            return Err(e);
        }
        for p in spills.iter().flatten().flatten() {
            let _ = std::fs::remove_file(p);
        }
        Ok((q, spills.len()))
    }
}
