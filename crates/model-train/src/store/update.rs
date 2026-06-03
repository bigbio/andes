//! Incremental update operations for the Parquet model store.
//!
//! These functions allow adding, removing, or reweighting per-source
//! [`CountStats`] and then re-estimating a candidate [`Param`] via exact
//! count arithmetic — WITHOUT modifying the on-disk store until the caller
//! calls [`commit_update`].
//!
//! # Exact add / remove semantics
//!
//! Each source stores its **unscaled** [`CountStats`] in the Parquet file.
//! When building the combined statistics, each source's counts are scaled by
//! its `weight` before summing:
//!
//! ```text
//! combined = Σ source_i.stats.scaled(source_i.ledger.weight)
//! ```
//!
//! This means `update_add` followed by `update_remove` of the same source is
//! exact: `scaled(w) + scaled(w) - scaled(w) = scaled(w)`, and the re-sum
//! over only the original sources reproduces the original combined counts
//! **to integer precision** (no floating-point drift: `scaled` rounds to `u64`).
//! The add-then-remove restore test validates this end-to-end.
//!
//! # Decay without dates
//!
//! [`update_decay`] reads each source's `date` field and computes age decay
//! `weight *= exp(-ln(2) * age_days / half_life_days)`. When a source's `date`
//! field is empty (as written by `train` when `--date` is omitted), its age
//! cannot be computed; the function **leaves that source's weight unchanged**
//! and emits an `eprintln!` warning noting the source ID. This is documented
//! as the decay-without-dates contract.

use std::path::Path;

use crate::counts::CountStats;
use crate::estimate::{Estimator, EstimatorConfig};
use crate::store::read::ModelStore;
use crate::store::write::{write_model_with_sources, SourceLedger};
use crate::TrainError;
use scoring_crate::param_model::Param;

// ---------------------------------------------------------------------------
// Type aliases (reduce type-complexity lint noise)
// ---------------------------------------------------------------------------

/// One entry in the all-models list: `(model_id, Param, sources)`.
type ModelEntry = (String, Param, Vec<(SourceLedger, CountStats)>);

/// A slice of borrowed model entries for writing.
type ModelSlice<'a> = &'a [(&'a str, &'a Param, &'a [(SourceLedger, CountStats)])];

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Append a new source to `model_id`, re-sum all per-source stats (applying
/// each source's `weight` via `stats.scaled(weight)`), and re-estimate a
/// candidate [`Param`].
///
/// **Does NOT persist** the candidate.  Call [`commit_update`] to write it.
pub fn update_add(
    path: &Path,
    model_id: &str,
    ledger: SourceLedger,
    stats: CountStats,
    est_cfg: EstimatorConfig,
) -> Result<(Param, Vec<(SourceLedger, CountStats)>), TrainError> {
    let store = ModelStore::open(path)?;
    let template = store.load_param(model_id)?;
    let existing_ledgers = store.load_sources(model_id)?;

    // Load per-source stats for all existing sources.
    let mut sources: Vec<(SourceLedger, CountStats)> = Vec::new();
    for l in existing_ledgers {
        let s = store.load_source_stats(model_id, &l.source_id)?;
        sources.push((l, s));
    }

    // Append the new source.
    sources.push((ledger, stats));

    let candidate = estimate_from_sources(&sources, &template, est_cfg)?;
    Ok((candidate, sources))
}

/// Remove the source with `source_id` from `model_id`, re-sum remaining
/// sources, and re-estimate a candidate [`Param`].
///
/// Returns [`TrainError::NoModel`] if `source_id` is not present.
///
/// **Does NOT persist** the candidate.  Call [`commit_update`] to write it.
pub fn update_remove(
    path: &Path,
    model_id: &str,
    source_id: &str,
    est_cfg: EstimatorConfig,
) -> Result<(Param, Vec<(SourceLedger, CountStats)>), TrainError> {
    let store = ModelStore::open(path)?;
    let template = store.load_param(model_id)?;
    let existing_ledgers = store.load_sources(model_id)?;

    // Check the source exists.
    if !existing_ledgers.iter().any(|l| l.source_id == source_id) {
        return Err(TrainError::NoModel(format!(
            "source '{source_id}' not found in model '{model_id}'"
        )));
    }

    // Load remaining sources (excluding the one being removed).
    let mut sources: Vec<(SourceLedger, CountStats)> = Vec::new();
    for l in existing_ledgers {
        if l.source_id == source_id {
            continue;
        }
        let s = store.load_source_stats(model_id, &l.source_id)?;
        sources.push((l, s));
    }

    let candidate = estimate_from_sources(&sources, &template, est_cfg)?;
    Ok((candidate, sources))
}

/// Update the `weight` field of `source_id` in `model_id`, re-sum all sources
/// with the new weight, and re-estimate a candidate [`Param`].
///
/// Returns [`TrainError::NoModel`] if `source_id` is not present.
///
/// **Does NOT persist** the candidate.  Call [`commit_update`] to write it.
pub fn update_reweight(
    path: &Path,
    model_id: &str,
    source_id: &str,
    weight: f32,
    est_cfg: EstimatorConfig,
) -> Result<(Param, Vec<(SourceLedger, CountStats)>), TrainError> {
    let store = ModelStore::open(path)?;
    let template = store.load_param(model_id)?;
    let existing_ledgers = store.load_sources(model_id)?;

    if !existing_ledgers.iter().any(|l| l.source_id == source_id) {
        return Err(TrainError::NoModel(format!(
            "source '{source_id}' not found in model '{model_id}'"
        )));
    }

    let mut sources: Vec<(SourceLedger, CountStats)> = Vec::new();
    for mut l in existing_ledgers {
        if l.source_id == source_id {
            l.weight = weight;
        }
        let s = store.load_source_stats(model_id, &l.source_id)?;
        sources.push((l, s));
    }

    let candidate = estimate_from_sources(&sources, &template, est_cfg)?;
    Ok((candidate, sources))
}

/// Apply exponential age-decay to each source's weight and re-estimate a
/// candidate [`Param`].
///
/// For each source whose `date` field is a valid `YYYY-MM-DD` string, the new
/// weight is:
/// ```text
/// w' = w * exp(-ln(2) * age_days / half_life_days)
/// ```
///
/// Sources with an **empty or unparseable `date`** are left with their current
/// weight unchanged, and a warning is printed to stderr.
///
/// **Does NOT persist** the candidate.  Call [`commit_update`] to write it.
pub fn update_decay(
    path: &Path,
    model_id: &str,
    half_life_days: f32,
    est_cfg: EstimatorConfig,
) -> Result<(Param, Vec<(SourceLedger, CountStats)>), TrainError> {
    // Guard against a non-positive half-life: `exp(-ln2 * age / 0)` underflows to
    // 0, silently zeroing every source weight and collapsing the model to a flat
    // distribution. Reject it (and NaN) instead of corrupting the store.
    if half_life_days <= 0.0 || half_life_days.is_nan() {
        return Err(TrainError::Other(format!(
            "decay half-life must be > 0 days, got {half_life_days}"
        )));
    }
    let today = today_naive();

    let store = ModelStore::open(path)?;
    let template = store.load_param(model_id)?;
    let existing_ledgers = store.load_sources(model_id)?;

    let mut sources: Vec<(SourceLedger, CountStats)> = Vec::new();
    for mut l in existing_ledgers {
        let s = store.load_source_stats(model_id, &l.source_id)?;

        if l.date.is_empty() {
            eprintln!(
                "update_decay: source '{}' has no date — weight unchanged",
                l.source_id
            );
        } else {
            match parse_date(&l.date) {
                Some(date_days) => {
                    let age_days = today.saturating_sub(date_days) as f32;
                    let decay = (-std::f32::consts::LN_2 * age_days / half_life_days).exp();
                    l.weight *= decay;
                }
                None => {
                    eprintln!(
                        "update_decay: source '{}' has unparseable date '{}' — weight unchanged",
                        l.source_id, l.date
                    );
                }
            }
        }

        sources.push((l, s));
    }

    let candidate = estimate_from_sources(&sources, &template, est_cfg)?;
    Ok((candidate, sources))
}

/// Persist a candidate [`Param`] and its updated source list to the store at
/// `path` for `model_id`.
///
/// All other models in the store are read, kept unchanged, and rewritten
/// together with the new state for `model_id` (read-all-others + rewrite-all).
pub fn commit_update(
    path: &Path,
    model_id: &str,
    new_param: &Param,
    new_sources: &[(SourceLedger, CountStats)],
) -> Result<(), TrainError> {
    // Read all existing models into memory.
    let mut other_models: Vec<ModelEntry> = Vec::new();

    if path.exists() {
        let store = ModelStore::open(path)?;
        for id in store.model_ids() {
            if id == model_id {
                continue; // will be replaced
            }
            let p = store.load_param(&id)?;
            let ledgers = store.load_sources(&id)?;
            let mut src = Vec::new();
            for l in ledgers {
                let s = store.load_source_stats(&id, &l.source_id)?;
                src.push((l, s));
            }
            other_models.push((id, p, src));
        }
    }

    // Write to a temp file, then rename atomically.
    let tmp = path.with_extension("parquet.tmp");

    // Write the updated model first, then the other models.
    // We reuse write_model_with_sources for each model sequentially by
    // writing into the same file via a multi-model approach.
    //
    // Current write_model_with_sources only supports ONE model per call.
    // We combine all models into a single file by writing each sequentially
    // into a temp file and then combining.  The simplest approach: write
    // all (model_id, param, sources) tuples into the same file.
    //
    // Since write_model_with_sources creates a new file each time, we write
    // each model to a separate temp file and then merge.  Actually, looking at
    // the write path: the combined schema supports MULTIPLE models in one file
    // (each row is keyed by model_id). We can call write_models for the plain
    // params and then write source rows for each model.
    //
    // The simplest correct path: write ALL models (including the updated one)
    // using a single write_multi_model_with_sources helper.
    write_all_models_with_sources(
        &tmp,
        std::iter::once((model_id.to_string(), new_param, new_sources.to_vec()))
            .chain(other_models.iter().map(|(id, p, s)| (id.clone(), p, s.clone())))
            .collect::<Vec<_>>()
            .iter()
            .map(|(id, p, s)| (id.as_str(), *p, s.as_slice()))
            .collect::<Vec<_>>()
            .as_slice(),
    )?;

    // Atomic rename.
    std::fs::rename(&tmp, path)?;
    Ok(())
}

// ---------------------------------------------------------------------------
// Private helpers
// ---------------------------------------------------------------------------

/// Re-sum all per-source stats (each scaled by its `weight`) and estimate a
/// new [`Param`] from the combined counts.
fn estimate_from_sources(
    sources: &[(SourceLedger, CountStats)],
    template: &Param,
    est_cfg: EstimatorConfig,
) -> Result<Param, TrainError> {
    let mut combined = CountStats::new();
    for (ledger, stats) in sources {
        combined.add(&stats.scaled(ledger.weight));
    }
    let estimator = Estimator::new(est_cfg);
    Ok(estimator.estimate(&combined, template))
}

/// Write multiple models (each with their per-source stats) to a single
/// Parquet file.
///
/// This is the multi-model equivalent of [`write_model_with_sources`], which
/// only supports a single model.  We implement it by serialising each model
/// into its own temp file and then combining the Arrow batches.
///
/// The combined schema supports multiple models (all rows are keyed by
/// `model_id`).  We therefore write each model to a temporary Parquet file
/// and then concatenate via the Arrow reader.
///
/// # Correctness
///
/// The schema is identical for all temp files (same `combined_schema()`), so
/// concatenation via the Arrow batch API is correct and lossless.
///
/// This function is `pub` so the binary crate can call it directly for the
/// initial-training path which needs to record sources.
pub fn write_all_models_with_sources_pub(
    path: &Path,
    models: ModelSlice<'_>,
) -> Result<(), TrainError> {
    write_all_models_with_sources(path, models)
}

fn write_all_models_with_sources(
    path: &Path,
    models: ModelSlice<'_>,
) -> Result<(), TrainError> {
    use crate::store::schema::combined_schema;
    use parquet::arrow::ArrowWriter;
    use parquet::file::properties::WriterProperties;
    use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

    let schema = combined_schema();
    let props = WriterProperties::builder().build();
    let file = std::fs::File::create(path)?;
    let mut writer = ArrowWriter::try_new(file, schema.clone(), Some(props))
        .map_err(|e| TrainError::Parquet(e.to_string()))?;

    for &(model_id, param, sources) in models {
        // Write each model to a temp buffer and read back the batches.
        let tmp_path = path.with_extension(format!("tmp_{model_id}.parquet"));
        write_model_with_sources(&tmp_path, model_id, param, sources)?;

        let tmp_file = std::fs::File::open(&tmp_path)?;
        let builder = ParquetRecordBatchReaderBuilder::try_new(tmp_file)
            .map_err(|e| TrainError::Parquet(e.to_string()))?;
        let reader = builder.build().map_err(|e| TrainError::Parquet(e.to_string()))?;
        for batch in reader {
            let batch = batch.map_err(|e| TrainError::Parquet(e.to_string()))?;
            if batch.num_rows() > 0 {
                writer.write(&batch).map_err(|e| TrainError::Parquet(e.to_string()))?;
            }
        }
        let _ = std::fs::remove_file(&tmp_path); // best-effort cleanup
    }

    writer.close().map_err(|e| TrainError::Parquet(e.to_string()))?;
    Ok(())
}

// ---------------------------------------------------------------------------
// Date helpers (no-std-time-free: we use std::time for "today")
// ---------------------------------------------------------------------------

/// Return today as "days since 1970-01-01" (Julian day count).
/// Uses `std::time::SystemTime`. Falls back to 0 on error.
fn today_naive() -> u32 {
    match std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH) {
        Ok(d) => (d.as_secs() / 86400) as u32,
        Err(_) => 0,
    }
}

/// Parse a `YYYY-MM-DD` date string into days since 1970-01-01.
/// Returns `None` if the string is empty or malformed.
fn parse_date(date: &str) -> Option<u32> {
    if date.len() != 10 {
        return None;
    }
    let y: i32 = date[..4].parse().ok()?;
    let m: u32 = date[5..7].parse().ok()?;
    let d: u32 = date[8..10].parse().ok()?;

    if !(1..=12).contains(&m) || !(1..=31).contains(&d) {
        return None;
    }

    // Convert to Julian Day Number and subtract offset for 1970-01-01.
    // Using the Gregorian formula (valid for 1970+):
    let jdn = julian_day_number(y, m, d);
    let epoch_jdn = julian_day_number(1970, 1, 1);
    if jdn >= epoch_jdn {
        Some((jdn - epoch_jdn) as u32)
    } else {
        Some(0)
    }
}

fn julian_day_number(y: i32, m: u32, d: u32) -> i64 {
    let a = (14 - m as i32) / 12;
    let y2 = y + 4800 - a;
    let m2 = m as i32 + 12 * a - 3;
    d as i64 + (153 * m2 + 2) as i64 / 5
        + 365 * y2 as i64
        + y2 as i64 / 4
        - y2 as i64 / 100
        + y2 as i64 / 400
        - 32045
}
