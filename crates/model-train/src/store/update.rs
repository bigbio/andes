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
//! its effective weight before summing:
//!
//! ```text
//! combined = Σ source_i.stats.scaled(effective_weight(source_i.ledger))
//! ```
//!
//! This means `update_add` followed by `update_remove` of the same source is
//! exact: `scaled(w) + scaled(w) - scaled(w) = scaled(w)`, and the re-sum
//! over only the original sources reproduces the original combined counts
//! **to integer precision** (no floating-point drift: `scaled` rounds to `u64`).
//! The add-then-remove restore test validates this end-to-end.
//!
//! # Decay
//!
//! A ledger's `weight` is the undecayed base weight. [`update_decay`] records a
//! half-life on every source (`decay_half_life_days`); every re-estimation then
//! uses
//!
//! ```text
//! effective_weight = weight * exp(-ln(2) * age_days / half_life_days)
//! ```
//!
//! with `age_days` measured from the source's `date` to the day of estimation.
//! Decay is therefore a function of age, not of how many times `--decay` ran:
//! two decays on the same day give the same model as one, and a later
//! `--add` / `--reweight` / `--remove-source` keeps the decay applied.
//!
//! A source whose `date` is empty (as written by `train` when `--date` is
//! omitted) or unparseable has no computable age; its weight is used
//! undecayed and [`update_decay`] prints a warning naming the source.

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
// Weight validation (finding 3.5)
// ---------------------------------------------------------------------------

/// Upper bound on a source weight. A weight is a count multiplier; a value this
/// large would either be a units mistake or risk saturating the `u64` counts in
/// [`CountStats::scaled`]. Far above any legitimate re-weighting.
const MAX_SOURCE_WEIGHT: f32 = 1.0e6;

/// Validate a source weight before it scales count statistics.
///
/// An unchecked weight silently corrupts the combined model: `NaN`/negative
/// rounds counts to `0` (erasing the source), and a huge value saturates the
/// `u64` counts in [`CountStats::scaled`]. Require a finite, non-negative,
/// bounded weight at every ingestion point (the `--add`/`--reweight` API and
/// the Parquet-loaded ledgers funneled through `estimate_from_sources`).
fn validate_weight(source_id: &str, weight: f32) -> Result<(), TrainError> {
    if !weight.is_finite() || !(0.0..=MAX_SOURCE_WEIGHT).contains(&weight) {
        return Err(TrainError::Other(format!(
            "source '{source_id}' has invalid weight {weight}: must be finite, \
             >= 0.0, and <= {MAX_SOURCE_WEIGHT}"
        )));
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Reject a store path the update and attach writers cannot rewrite in place.
///
/// These writers replace a single Parquet file atomically. A partitioned
/// (directory) store such as the bundled `resources/models/` is laid out per
/// protocol by `scripts/split_store_by_protocol.py`, so it is read-only here:
/// copy the model into a single-file store, update that, and re-split.
pub fn ensure_single_file_store(path: &Path) -> Result<(), TrainError> {
    if path.is_dir() {
        return Err(TrainError::Other(format!(
            "{} is a partitioned (directory) model store, which cannot be rewritten \
             in place; write the model into a single-file store (e.g. models.parquet) \
             and re-split it with scripts/split_store_by_protocol.py",
            path.display()
        )));
    }
    Ok(())
}

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
    // Validate the incoming weight up front (finding 3.5) so a NaN/negative/huge
    // value fails loudly at the API rather than silently erasing or saturating
    // counts.
    validate_weight(&ledger.source_id, ledger.weight)?;

    ensure_single_file_store(path)?;
    let store = ModelStore::open(path)?;
    let template = store.load_param(model_id)?;
    let existing_ledgers = store.load_sources(model_id)?;

    // Reject a duplicate source_id: appending it would list the source twice and
    // `estimate_from_sources` would double-count its weighted stats. Callers that
    // want to change an existing source should use --reweight / --remove-source.
    if existing_ledgers
        .iter()
        .any(|l| l.source_id == ledger.source_id)
    {
        return Err(TrainError::Other(format!(
            "source '{}' already exists in model '{model_id}'; \
             use --reweight to change its weight or --remove-source first",
            ledger.source_id
        )));
    }

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
    ensure_single_file_store(path)?;
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
    // Validate the new weight up front (finding 3.5).
    validate_weight(source_id, weight)?;

    ensure_single_file_store(path)?;
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

/// Set an exponential age-decay half-life on every source of `model_id` and
/// re-estimate a candidate [`Param`].
///
/// The ledger weights stay the undecayed base weights; the decay is applied by
/// [`estimate_from_sources`] from each source's `date` (see the module docs),
/// so repeating the call with the same half-life on the same day is a no-op.
/// Sources with an **empty or unparseable `date`** stay undecayed, and a
/// warning is printed to stderr.
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
    // distribution. Reject it (and NaN/inf) instead of corrupting the store.
    if !(half_life_days.is_finite() && half_life_days > 0.0) {
        return Err(TrainError::Other(format!(
            "decay half-life must be a finite number of days > 0, got {half_life_days}"
        )));
    }

    ensure_single_file_store(path)?;
    let store = ModelStore::open(path)?;
    let template = store.load_param(model_id)?;
    let existing_ledgers = store.load_sources(model_id)?;

    let mut sources: Vec<(SourceLedger, CountStats)> = Vec::new();
    for mut l in existing_ledgers {
        let s = store.load_source_stats(model_id, &l.source_id)?;
        if l.date.is_empty() {
            eprintln!(
                "update_decay: source '{}' has no date — weight not decayed",
                l.source_id
            );
        } else if parse_date(&l.date).is_none() {
            eprintln!(
                "update_decay: source '{}' has unparseable date '{}' — weight not decayed",
                l.source_id, l.date
            );
        }
        l.decay_half_life_days = Some(half_life_days);
        sources.push((l, s));
    }

    let candidate = estimate_from_sources(&sources, &template, est_cfg)?;
    Ok((candidate, sources))
}

/// Weight a source contributes at estimation time on day `today` (days since
/// 1970-01-01): the base `weight`, decayed by age when the ledger carries a
/// half-life and a parseable `date`.
fn effective_weight(ledger: &SourceLedger, today: u32) -> f32 {
    match (ledger.decay_half_life_days, parse_date(&ledger.date)) {
        (Some(half_life), Some(date_days)) => {
            let age_days = today.saturating_sub(date_days) as f32;
            ledger.weight * (-std::f32::consts::LN_2 * age_days / half_life).exp()
        }
        _ => ledger.weight,
    }
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
    ensure_single_file_store(path)?;
    // Read all existing models into memory, preserving their GBDT blobs.
    let mut other_models: Vec<ModelEntry> = Vec::new();
    let mut other_blobs: Vec<Option<Vec<u8>>> = Vec::new();

    if path.exists() {
        let store = ModelStore::open(path)?;
        for id in store.model_ids() {
            if id == model_id {
                continue; // will be replaced
            }
            let p = store.load_param(&id)?;
            // Preserve the GBDT blob (if any) for the existing model.
            let blob = p.gbdt_peak_model.as_ref().map(|m| m.to_bytes());
            let ledgers = store.load_sources(&id)?;
            let mut src = Vec::new();
            for l in ledgers {
                let s = store.load_source_stats(&id, &l.source_id)?;
                src.push((l, s));
            }
            other_models.push((id, p, src));
            other_blobs.push(blob);
        }
    }

    // Write to a temp file, then rename atomically.
    let tmp = path.with_extension("parquet.tmp");

    // Write ALL models (including the updated one) using a single
    // write_all_models_with_sources_inner call. The updated model inherits
    // any GBDT blob from new_param.gbdt_peak_model.
    let all_entries: Vec<_> =
        std::iter::once((model_id.to_string(), new_param, new_sources.to_vec()))
            .chain(
                other_models
                    .iter()
                    .map(|(id, p, s)| (id.clone(), p, s.clone())),
            )
            .collect();
    let new_blob = new_param.gbdt_peak_model.as_ref().map(|m| m.to_bytes());
    let all_blobs: Vec<Option<Vec<u8>>> = std::iter::once(new_blob).chain(other_blobs).collect();
    let entries_slice: Vec<_> = all_entries
        .iter()
        .map(|(id, p, s)| (id.as_str(), *p, s.as_slice()))
        .collect();
    write_all_models_with_sources_inner(&tmp, &entries_slice, &all_blobs)?;

    // Commit the temp file. `rename` replaces atomically on Unix, but on Windows
    // it fails when the destination already exists — so fall back to
    // remove-then-rename there (only when the rename actually failed and the
    // destination is present), keeping the Unix path atomic.
    match std::fs::rename(&tmp, path) {
        Ok(()) => {}
        Err(_) if path.exists() => {
            std::fs::remove_file(path)?;
            std::fs::rename(&tmp, path)?;
        }
        Err(e) => return Err(e.into()),
    }
    Ok(())
}

/// Attach separately trained components (e.g. a fragment-intensity or rich-ion
/// GBDT) to `model_id` in the single-file store at `path`, keeping everything
/// else about the model.
///
/// When the store already holds `model_id`, that model — its rank tables,
/// its other GBDT blobs and its source ledger with per-source stats — is the
/// base. Otherwise `fallback` (typically the seed model the components were
/// trained against) is the base and the model starts with no sources. `attach`
/// edits the base, and the result is committed with [`commit_update`], so all
/// other models in the store are preserved. Returns whether `model_id`
/// already existed.
pub fn attach_to_model(
    path: &Path,
    model_id: &str,
    fallback: Param,
    attach: impl FnOnce(&mut Param),
) -> Result<bool, TrainError> {
    ensure_single_file_store(path)?;
    let existing = if path.exists() {
        let store = ModelStore::open(path)?;
        if store.model_ids().iter().any(|id| id == model_id) {
            let param = store.load_param(model_id)?;
            let mut sources = Vec::new();
            for l in store.load_sources(model_id)? {
                let s = store.load_source_stats(model_id, &l.source_id)?;
                sources.push((l, s));
            }
            Some((param, sources))
        } else {
            None
        }
    } else {
        None
    };
    let existed = existing.is_some();
    let (mut param, sources) = existing.unwrap_or((fallback, Vec::new()));
    attach(&mut param);
    commit_update(path, model_id, &param, &sources)?;
    Ok(existed)
}

// ---------------------------------------------------------------------------
// Private helpers
// ---------------------------------------------------------------------------

/// Re-sum all per-source stats (each scaled by its effective weight, see
/// [`effective_weight`]) and estimate a
/// new [`Param`] from the combined counts.
///
/// The template's GBDT models (peak, fragment-intensity, rich-ion) are carried
/// onto the candidate unchanged. Their features are spectrum-derived (intensity
/// ranks, peak geometry, fragment m/z errors) and never read the rank,
/// error or existence tables re-estimated here, so they stay valid for the
/// candidate. Retraining them on new data is a separate step
/// (`train-intensity-gbdt`, `train-rich-ion-llr`, `train-msnet --gbdt`).
fn estimate_from_sources(
    sources: &[(SourceLedger, CountStats)],
    template: &Param,
    est_cfg: EstimatorConfig,
) -> Result<Param, TrainError> {
    let today = today_naive();
    let mut combined = CountStats::new();
    for (ledger, stats) in sources {
        // Validate every weight before scaling (finding 3.5) — this is the
        // chokepoint that also covers weights loaded from Parquet ledgers, not
        // just the API entry points.
        validate_weight(&ledger.source_id, ledger.weight)?;
        if let Some(h) = ledger.decay_half_life_days {
            if !(h.is_finite() && h > 0.0) {
                return Err(TrainError::Other(format!(
                    "source '{}' has invalid decay half-life {h}",
                    ledger.source_id
                )));
            }
        }
        combined.add(&stats.scaled(effective_weight(ledger, today)));
    }
    let estimator = Estimator::new(est_cfg);
    let mut candidate = estimator.estimate(&combined, template);
    candidate.gbdt_peak_model = template.gbdt_peak_model.clone();
    candidate.frag_intensity_model = template.frag_intensity_model.clone();
    candidate.rich_ion_model = template.rich_ion_model.clone();
    Ok(candidate)
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
    let blobs = vec![None; models.len()];
    write_all_models_with_sources_inner(path, models, &blobs)
}

/// Like [`write_all_models_with_sources_pub`] but each model may carry a GBDT
/// blob (parallel to `models`). `None` ⇒ no blob (rank-core-only).
///
/// `blobs` must have the same length as `models` (enforced by a `debug_assert`).
/// Each entry corresponds to the model at the same index in `models`.
pub fn write_all_models_with_sources_and_gbdt_pub(
    path: &Path,
    models: ModelSlice<'_>,
    blobs: &[Option<Vec<u8>>],
) -> Result<(), TrainError> {
    write_all_models_with_sources_inner(path, models, blobs)
}

fn write_all_models_with_sources_inner(
    path: &Path,
    models: ModelSlice<'_>,
    blobs: &[Option<Vec<u8>>],
) -> Result<(), TrainError> {
    debug_assert_eq!(
        blobs.len(),
        models.len(),
        "blobs slice must be 1:1 with models"
    );
    use crate::store::schema::combined_schema;
    use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
    use parquet::arrow::ArrowWriter;
    use parquet::file::properties::WriterProperties;

    let schema = combined_schema();
    let props = WriterProperties::builder().build();
    let file = std::fs::File::create(path)?;
    let mut writer = ArrowWriter::try_new(file, schema.clone(), Some(props))
        .map_err(|e| TrainError::Parquet(e.to_string()))?;

    for (i, &(model_id, param, sources)) in models.iter().enumerate() {
        // Write each model to a temp buffer and read back the batches.
        let tmp_path = path.with_extension(format!("tmp_{model_id}.parquet"));
        let blob = blobs.get(i).cloned().flatten();
        write_model_with_sources(&tmp_path, model_id, param, sources, blob)?;

        let tmp_file = std::fs::File::open(&tmp_path)?;
        let builder = ParquetRecordBatchReaderBuilder::try_new(tmp_file)
            .map_err(|e| TrainError::Parquet(e.to_string()))?;
        let reader = builder
            .build()
            .map_err(|e| TrainError::Parquet(e.to_string()))?;
        for batch in reader {
            let batch = batch.map_err(|e| TrainError::Parquet(e.to_string()))?;
            if batch.num_rows() > 0 {
                writer
                    .write(&batch)
                    .map_err(|e| TrainError::Parquet(e.to_string()))?;
            }
        }
        let _ = std::fs::remove_file(&tmp_path); // best-effort cleanup
    }

    writer
        .close()
        .map_err(|e| TrainError::Parquet(e.to_string()))?;
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
    // Pre-1970 dates are out of the supported range: treat them as invalid
    // (→ weight unchanged) rather than as 1970-01-01 (which would apply an
    // extreme decay).
    if jdn >= epoch_jdn {
        Some((jdn - epoch_jdn) as u32)
    } else {
        None
    }
}

fn julian_day_number(y: i32, m: u32, d: u32) -> i64 {
    let a = (14 - m as i32) / 12;
    let y2 = y + 4800 - a;
    let m2 = m as i32 + 12 * a - 3;
    d as i64 + (153 * m2 + 2) as i64 / 5 + 365 * y2 as i64 + y2 as i64 / 4 - y2 as i64 / 100
        + y2 as i64 / 400
        - 32045
}

#[cfg(test)]
mod weight_validation_tests {
    use super::{validate_weight, MAX_SOURCE_WEIGHT};

    #[test]
    fn accepts_valid_weights() {
        assert!(validate_weight("s", 0.0).is_ok());
        assert!(validate_weight("s", 1.0).is_ok());
        assert!(validate_weight("s", 0.5).is_ok());
        assert!(validate_weight("s", MAX_SOURCE_WEIGHT).is_ok());
    }

    #[test]
    fn rejects_invalid_weights() {
        assert!(validate_weight("s", f32::NAN).is_err(), "NaN");
        assert!(validate_weight("s", f32::INFINITY).is_err(), "inf");
        assert!(validate_weight("s", f32::NEG_INFINITY).is_err(), "-inf");
        assert!(validate_weight("s", -1.0).is_err(), "negative");
        assert!(
            validate_weight("s", MAX_SOURCE_WEIGHT * 2.0).is_err(),
            "too large"
        );
    }
}
