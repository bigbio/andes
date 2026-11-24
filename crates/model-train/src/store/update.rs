use std::path::Path;

use crate::counts::CountStats;
use crate::estimate::{Estimator, EstimatorConfig};
use crate::store::read::ModelStore;
use crate::store::write::{write_model_with_sources, SourceLedger};
use crate::TrainError;
use scoring_crate::params::Param;

type ModelEntry = (String, Param, Vec<(SourceLedger, CountStats)>);

type ModelSlice<'a> = &'a [(&'a str, &'a Param, &'a [(SourceLedger, CountStats)])];

const MAX_SOURCE_WEIGHT: f32 = 1.0e6;

fn validate_weight(source_id: &str, weight: f32) -> Result<(), TrainError> {
    if !weight.is_finite() || !(0.0..=MAX_SOURCE_WEIGHT).contains(&weight) {
        return Err(TrainError::Other(format!(
            "source '{source_id}' has invalid weight {weight}: must be finite, \
             >= 0.0, and <= {MAX_SOURCE_WEIGHT}"
        )));
    }
    Ok(())
}

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

pub fn update_add(
    path: &Path,
    model_id: &str,
    ledger: SourceLedger,
    stats: CountStats,
    est_cfg: EstimatorConfig,
) -> Result<(Param, Vec<(SourceLedger, CountStats)>), TrainError> {
    validate_weight(&ledger.source_id, ledger.weight)?;

    ensure_single_file_store(path)?;
    let store = ModelStore::open(path)?;
    let template = store.load_param(model_id)?;
    let existing_ledgers = store.load_sources(model_id)?;

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

    let mut sources: Vec<(SourceLedger, CountStats)> = Vec::new();
    for l in existing_ledgers {
        let s = store.load_source_stats(model_id, &l.source_id)?;
        sources.push((l, s));
    }

    sources.push((ledger, stats));

    let candidate = estimate_from_sources(&sources, &template, est_cfg)?;
    Ok((candidate, sources))
}

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

    if !existing_ledgers.iter().any(|l| l.source_id == source_id) {
        return Err(TrainError::NoModel(format!(
            "source '{source_id}' not found in model '{model_id}'"
        )));
    }

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

pub fn update_reweight(
    path: &Path,
    model_id: &str,
    source_id: &str,
    weight: f32,
    est_cfg: EstimatorConfig,
) -> Result<(Param, Vec<(SourceLedger, CountStats)>), TrainError> {
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

pub fn update_decay(
    path: &Path,
    model_id: &str,
    half_life_days: f32,
    est_cfg: EstimatorConfig,
) -> Result<(Param, Vec<(SourceLedger, CountStats)>), TrainError> {
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
            log::warn!(
                target: "update-decay",
                "source '{}' has no date — weight not decayed",
                l.source_id
            );
        } else if parse_date(&l.date).is_none() {
            log::warn!(
                target: "update-decay",
                "source '{}' has unparseable date '{}' — weight not decayed",
                l.source_id, l.date
            );
        }
        l.decay_half_life_days = Some(half_life_days);
        sources.push((l, s));
    }

    let candidate = estimate_from_sources(&sources, &template, est_cfg)?;
    Ok((candidate, sources))
}

fn effective_weight(ledger: &SourceLedger, today: u32) -> f32 {
    match (ledger.decay_half_life_days, parse_date(&ledger.date)) {
        (Some(half_life), Some(date_days)) => {
            let age_days = today.saturating_sub(date_days) as f32;
            ledger.weight * (-std::f32::consts::LN_2 * age_days / half_life).exp()
        }
        _ => ledger.weight,
    }
}

pub fn commit_update(
    path: &Path,
    model_id: &str,
    new_param: &Param,
    new_sources: &[(SourceLedger, CountStats)],
) -> Result<(), TrainError> {
    ensure_single_file_store(path)?;
    let mut other_models: Vec<ModelEntry> = Vec::new();
    let mut other_blobs: Vec<Option<Vec<u8>>> = Vec::new();

    if path.exists() {
        let store = ModelStore::open(path)?;
        for id in store.model_ids() {
            if id == model_id {
                continue;
            }
            let p = store.load_param(&id)?;
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

    let tmp = path.with_extension("parquet.tmp");

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

fn estimate_from_sources(
    sources: &[(SourceLedger, CountStats)],
    template: &Param,
    est_cfg: EstimatorConfig,
) -> Result<Param, TrainError> {
    let today = today_naive();
    let mut combined = CountStats::new();
    for (ledger, stats) in sources {
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

pub fn write_all_models_with_sources_pub(
    path: &Path,
    models: ModelSlice<'_>,
) -> Result<(), TrainError> {
    let blobs = vec![None; models.len()];
    write_all_models_with_sources_inner(path, models, &blobs)
}

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
        let _ = std::fs::remove_file(&tmp_path);
    }

    writer
        .close()
        .map_err(|e| TrainError::Parquet(e.to_string()))?;
    Ok(())
}

fn today_naive() -> u32 {
    match std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH) {
        Ok(d) => (d.as_secs() / 86400) as u32,
        Err(_) => 0,
    }
}

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

    let jdn = julian_day_number(y, m, d);
    let epoch_jdn = julian_day_number(1970, 1, 1);
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
