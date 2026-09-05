//! Intensity-model, fragment-intensity GBDT and rich-ion LLR training.

use std::path::{Path, PathBuf};

use crate::model_select::{load_seed_param, ModelEntryOwned};
use crate::train::{read_msnet_parquet, MsnetPsm};
use clap::{Args, ValueEnum};
use model_train::{store::write_all_models_with_sources_and_gbdt_pub, ModelStore};
use scoring_crate::RankScorer;

/// Training arguments for `andes train-intensity`.
///
/// Merges one or more partial intensity aggregation parquets (from
/// `msnet_intensity_agg.py`) into a finalized `intensity_model.parquet` with
/// `mean_log_rel` / `var_log_rel` columns for runtime lookup.
#[derive(Args, Debug)]
pub(crate) struct TrainIntensityArgs {
    /// Input partial or finalized intensity parquets. Repeatable; stats merge
    /// across all inputs.
    #[arg(long = "in", required = true)]
    pub(crate) inputs: Vec<PathBuf>,

    /// Output path for the finalized intensity model parquet.
    #[arg(long = "out", required = true)]
    pub(crate) out: PathBuf,
}

/// GBDT training mode for `andes train`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum)]
pub(crate) enum GbdtMode {
    /// Train and embed a GBDT peak model (default).
    #[default]
    #[clap(name = "on")]
    On,
    /// Skip GBDT; write rank-core only (byte-identical to pre-GBDT path).
    #[clap(name = "off")]
    Off,
}

/// Training arguments for `andes train-intensity-gbdt`.
///
/// Reads flat training parquets (same schema as `train-from-msnet`) and fits a
/// GBDT fragment-intensity regressor (`v3 frag model`).  The trained model is
/// written into `--out-store` alongside any existing models under `--model-id`.
#[derive(Args, Debug)]
pub(crate) struct TrainIntensityGbdtArgs {
    /// Input flat training parquet(s). Repeatable; data accumulate across all
    /// inputs into a single frag-intensity model.
    #[arg(long = "in", required = true)]
    pub(crate) inputs: Vec<PathBuf>,

    /// Path to the Parquet model store to write (created if absent; existing
    /// models are preserved and re-written alongside the new one). REQUIRED.
    #[arg(long = "out-store", required = true)]
    pub(crate) out_store: PathBuf,

    /// Model ID written into the store. Default: `default`.
    #[arg(long = "model-id", default_value = "default")]
    pub(crate) model_id: String,

    /// Seed model: slug from the bundled store (e.g. `hcd_qexactive_tryp`).
    /// Supplies structural hyperparameters
    /// (fragment tolerance, charge range) used when building the frag dataset.
    #[arg(long = "seed-model", default_value = "hcd_qexactive_tryp")]
    pub(crate) seed_model: String,

    /// Number of worker threads (Rayon). Default: 8.
    #[arg(long, default_value_t = 8usize)]
    pub(crate) threads: usize,

    /// Opt-in fallback (finding 3.6): when set, a failed GBDT quality gate
    /// (too few rows / no held-out signal / empty ensemble) is downgraded from a
    /// hard error to a warning and the degenerate model is still written. Default
    /// off — gate failures abort with a non-zero exit. Intended for small
    /// synthetic fixtures / benchmarking only.
    #[arg(long = "allow-degenerate-model", hide = true, default_value_t = false)]
    pub(crate) allow_degenerate_model: bool,
}

/// Training arguments for `andes train-rich-ion-llr`.
///
/// Reads flat training parquets (same schema as `train-from-msnet`) and fits a
/// GBDT rich-ion LLR classifier (logistic; decoy-aware).  The trained model is
/// written into `--out-store` alongside any existing models under `--model-id`.
#[derive(Args, Debug)]
pub(crate) struct TrainRichIonLlrArgs {
    /// Input flat training parquet(s). Repeatable; data accumulate across all
    /// inputs into a single rich-ion model.
    #[arg(long = "in", required = true)]
    pub(crate) inputs: Vec<PathBuf>,

    /// Path to the Parquet model store to write (created if absent; existing
    /// models are preserved and re-written alongside the new one). REQUIRED.
    #[arg(long = "out-store", required = true)]
    pub(crate) out_store: PathBuf,

    /// Model ID written into the store. Default: `default`.
    #[arg(long = "model-id", default_value = "default")]
    pub(crate) model_id: String,

    /// Seed model: slug from the bundled store (e.g. `hcd_qexactive_tryp`).
    /// Supplies structural hyperparameters
    /// (fragment tolerance, charge range) used when building the ion dataset.
    #[arg(long = "seed-model", default_value = "hcd_qexactive_tryp")]
    pub(crate) seed_model: String,

    /// Number of worker threads (Rayon). Default: 8.
    #[arg(long, default_value_t = 8usize)]
    pub(crate) threads: usize,

    /// Opt-in fallback (finding 3.6): downgrade a failed GBDT quality gate to a
    /// warning and write the degenerate model anyway. Default off.
    #[arg(long = "allow-degenerate-model", hide = true, default_value_t = false)]
    pub(crate) allow_degenerate_model: bool,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub(crate) struct IntensityAggKey {
    pub(crate) ion_type: String,
    pub(crate) flank_n: String,
    pub(crate) flank_c: String,
    pub(crate) pos_bin: i32,
    pub(crate) charge: i32,
    pub(crate) nce_bin: String,
}

#[derive(Debug, Clone, Default)]
pub(crate) struct IntensityAggStats {
    pub(crate) count: i64,
    pub(crate) sum_log_rel: f64,
    pub(crate) sum_log_rel_sq: f64,
}

pub(crate) fn read_intensity_partial(
    path: &Path,
) -> Result<Vec<(IntensityAggKey, IntensityAggStats)>, Box<dyn std::error::Error>> {
    use arrow::array::{Array, Float64Array, Int32Array, Int64Array, StringArray};
    use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

    let file = std::fs::File::open(path).map_err(|e| format!("open {}: {e}", path.display()))?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file)
        .map_err(|e| format!("parquet reader for {}: {e}", path.display()))?;
    let has_mean = builder.schema().field_with_name("mean_log_rel").is_ok();
    let mut rows = Vec::new();

    for batch_result in builder.build().map_err(|e| format!("build reader: {e}"))? {
        let batch = batch_result?;
        let ion_col = batch
            .column_by_name("ion_type")
            .and_then(|c| c.as_any().downcast_ref::<StringArray>())
            .ok_or("missing ion_type")?;
        let flank_n_col = batch
            .column_by_name("flank_n")
            .and_then(|c| c.as_any().downcast_ref::<StringArray>())
            .ok_or("missing flank_n")?;
        let flank_c_col = batch
            .column_by_name("flank_c")
            .and_then(|c| c.as_any().downcast_ref::<StringArray>())
            .ok_or("missing flank_c")?;
        let pos_col = batch
            .column_by_name("pos_bin")
            .and_then(|c| c.as_any().downcast_ref::<Int32Array>())
            .ok_or("missing pos_bin")?;
        let charge_col = batch
            .column_by_name("charge")
            .and_then(|c| c.as_any().downcast_ref::<Int32Array>())
            .ok_or("missing charge")?;
        let nce_col = batch
            .column_by_name("nce_bin")
            .and_then(|c| c.as_any().downcast_ref::<StringArray>())
            .ok_or("missing nce_bin")?;
        let count_col = batch
            .column_by_name("count")
            .and_then(|c| c.as_any().downcast_ref::<Int64Array>())
            .ok_or("missing count")?;

        let sum_col = if has_mean {
            None
        } else {
            Some(
                batch
                    .column_by_name("sum_log_rel")
                    .and_then(|c| c.as_any().downcast_ref::<Float64Array>())
                    .ok_or("missing sum_log_rel")?,
            )
        };
        let sum_sq_col = if has_mean {
            None
        } else {
            Some(
                batch
                    .column_by_name("sum_log_rel_sq")
                    .and_then(|c| c.as_any().downcast_ref::<Float64Array>())
                    .ok_or("missing sum_log_rel_sq")?,
            )
        };
        let mean_col = if has_mean {
            Some(
                batch
                    .column_by_name("mean_log_rel")
                    .and_then(|c| c.as_any().downcast_ref::<Float64Array>())
                    .ok_or("missing mean_log_rel")?,
            )
        } else {
            None
        };
        let var_col = if has_mean {
            Some(
                batch
                    .column_by_name("var_log_rel")
                    .and_then(|c| c.as_any().downcast_ref::<Float64Array>())
                    .ok_or("missing var_log_rel")?,
            )
        } else {
            None
        };

        for i in 0..batch.num_rows() {
            let count = count_col.value(i);
            let (sum, sum_sq) = if let (Some(sum_c), Some(sq_c)) = (sum_col, sum_sq_col) {
                (sum_c.value(i), sq_c.value(i))
            } else {
                let mean = mean_col
                    .ok_or("intensity partial: mean_log_rel column missing while in mean/var path")?
                    .value(i);
                let var = var_col
                    .ok_or("intensity partial: var_log_rel column missing while in mean/var path")?
                    .value(i);
                (mean * count as f64, (var + mean * mean) * count as f64)
            };
            let key = IntensityAggKey {
                ion_type: ion_col.value(i).to_string(),
                flank_n: flank_n_col.value(i).to_string(),
                flank_c: flank_c_col.value(i).to_string(),
                pos_bin: pos_col.value(i),
                charge: charge_col.value(i),
                nce_bin: nce_col.value(i).to_string(),
            };
            rows.push((
                key,
                IntensityAggStats {
                    count,
                    sum_log_rel: sum,
                    sum_log_rel_sq: sum_sq,
                },
            ));
        }
    }
    Ok(rows)
}

/// Finalize one aggregation cell into `(mean_log_rel, var_log_rel)`.
/// Returns `None` for `count <= 0` so empty cells (e.g. from a partial
/// aggregation parquet) are dropped instead of writing NaN. Variance is
/// clamped at 0 to absorb floating-point round-off in `E[x²] − E[x]²`.
pub(crate) fn finalize_intensity_stats(
    sum_log_rel: f64,
    sum_log_rel_sq: f64,
    count: i64,
) -> Option<(f64, f64)> {
    if count <= 0 {
        return None;
    }
    let n = count as f64;
    let mean = sum_log_rel / n;
    let var = (sum_log_rel_sq / n - mean * mean).max(0.0);
    Some((mean, var))
}

/// Sentinel flank residue for flank-marginalized backoff rows (matches the
/// `b'*'` key built by `IntensityModel::predict_log_rel`'s backoff).
pub(crate) const ANY_FLANK: &str = "*";
/// Sentinel nce bin for nce-marginalized backoff rows (matches `"__any__"`).
const ANY_NCE: &str = "__any__";

/// Emit backoff marginal rows so `IntensityModel::predict_log_rel`'s documented
/// sparse-key backoff (drop `nce_bin`, then flank residues) actually finds rows.
///
/// The aggregator only writes exact `(ion, flank_n, flank_c, pos_bin, charge,
/// nce_bin)` cells with real flanks and a real (numeric/`"unknown"`) nce bin.
/// The backoff probes `nce_bin="__any__"` and `flank=b'*'` keys that the
/// aggregator never produces, so without these marginals every lookup whose
/// exact key is absent falls through to the single global mean — making the
/// trained per-context intensities unused. In particular the inference path
/// (`compute_psm_features`) passes `nce_bin="unknown"`, which matches no trained
/// numeric bin, so *every* inference lookup hit the global mean and the
/// `IntensitySignal` cosine became model-value-independent.
///
/// We accumulate each real cell into three marginals: drop nce, drop flanks,
/// and drop both. Marginal keys are disjoint from real keys (real flanks are
/// residues, real nce bins are numeric or `"unknown"`), so snapshotting the
/// real cells first avoids double-counting.
pub(crate) fn add_backoff_marginals(
    merged: &mut rustc_hash::FxHashMap<IntensityAggKey, IntensityAggStats>,
) {
    // Drop any marginal rows already present (a merged map may come from partials
    // that carried them) and rebuild them from real cells only, so calling this
    // twice does not double the marginal counts.
    merged.retain(|k, _| k.nce_bin != ANY_NCE && k.flank_n != ANY_FLANK && k.flank_c != ANY_FLANK);
    let base: Vec<(IntensityAggKey, IntensityAggStats)> =
        merged.iter().map(|(k, s)| (k.clone(), s.clone())).collect();

    fn bump(
        merged: &mut rustc_hash::FxHashMap<IntensityAggKey, IntensityAggStats>,
        key: IntensityAggKey,
        s: &IntensityAggStats,
    ) {
        let slot = merged.entry(key).or_default();
        slot.count += s.count;
        slot.sum_log_rel += s.sum_log_rel;
        slot.sum_log_rel_sq += s.sum_log_rel_sq;
    }

    for (k, s) in &base {
        // Drop nce: keep flank/pos/charge structure, marginalize over collision energy.
        bump(
            merged,
            IntensityAggKey {
                nce_bin: ANY_NCE.to_string(),
                ..k.clone()
            },
            s,
        );
        // Drop flanks: keep pos/charge/nce, marginalize over the residue pair.
        bump(
            merged,
            IntensityAggKey {
                flank_n: ANY_FLANK.to_string(),
                flank_c: ANY_FLANK.to_string(),
                ..k.clone()
            },
            s,
        );
        // Drop both flanks and nce.
        bump(
            merged,
            IntensityAggKey {
                flank_n: ANY_FLANK.to_string(),
                flank_c: ANY_FLANK.to_string(),
                nce_bin: ANY_NCE.to_string(),
                ..k.clone()
            },
            s,
        );
    }
}

pub(crate) fn write_intensity_model(
    path: &Path,
    merged: &rustc_hash::FxHashMap<IntensityAggKey, IntensityAggStats>,
) -> Result<(), Box<dyn std::error::Error>> {
    use arrow::array::{Float64Array, Int32Array, Int64Array, StringArray};
    use arrow::datatypes::{DataType, Field, Schema};
    use arrow::record_batch::RecordBatch;
    use parquet::arrow::ArrowWriter;

    let mut keys: Vec<_> = merged.iter().collect();
    // Drop empty cells so they never write NaN mean/var (a partial aggregation
    // parquet can carry count==0 rows). finalize_intensity_stats gates on this.
    keys.retain(|(_, s)| s.count > 0);
    keys.sort_by(|a, b| {
        (
            &a.0.ion_type,
            &a.0.flank_n,
            &a.0.flank_c,
            a.0.pos_bin,
            a.0.charge,
            &a.0.nce_bin,
        )
            .cmp(&(
                &b.0.ion_type,
                &b.0.flank_n,
                &b.0.flank_c,
                b.0.pos_bin,
                b.0.charge,
                &b.0.nce_bin,
            ))
    });

    let ion: Vec<_> = keys.iter().map(|(k, _)| k.ion_type.as_str()).collect();
    let flank_n: Vec<_> = keys.iter().map(|(k, _)| k.flank_n.as_str()).collect();
    let flank_c: Vec<_> = keys.iter().map(|(k, _)| k.flank_c.as_str()).collect();
    let pos_bin: Vec<_> = keys.iter().map(|(k, _)| k.pos_bin).collect();
    let charge: Vec<_> = keys.iter().map(|(k, _)| k.charge).collect();
    let nce: Vec<_> = keys.iter().map(|(k, _)| k.nce_bin.as_str()).collect();
    let count: Vec<_> = keys.iter().map(|(_, s)| s.count).collect();
    // Safe to unwrap: zero-count keys were retained out above.
    let mean: Vec<_> = keys
        .iter()
        .map(|(_, s)| {
            finalize_intensity_stats(s.sum_log_rel, s.sum_log_rel_sq, s.count)
                .unwrap()
                .0
        })
        .collect();
    let var: Vec<_> = keys
        .iter()
        .map(|(_, s)| {
            finalize_intensity_stats(s.sum_log_rel, s.sum_log_rel_sq, s.count)
                .unwrap()
                .1
        })
        .collect();

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
    let batch = RecordBatch::try_new(
        std::sync::Arc::new(schema.clone()),
        vec![
            std::sync::Arc::new(StringArray::from(ion)),
            std::sync::Arc::new(StringArray::from(flank_n)),
            std::sync::Arc::new(StringArray::from(flank_c)),
            std::sync::Arc::new(Int32Array::from(pos_bin)),
            std::sync::Arc::new(Int32Array::from(charge)),
            std::sync::Arc::new(StringArray::from(nce)),
            std::sync::Arc::new(Int64Array::from(count)),
            std::sync::Arc::new(Float64Array::from(mean)),
            std::sync::Arc::new(Float64Array::from(var)),
        ],
    )?;
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let file = std::fs::File::create(path)?;
    let mut writer = ArrowWriter::try_new(file, std::sync::Arc::new(schema), None)?;
    writer.write(&batch)?;
    writer.close()?;
    Ok(())
}

pub(crate) fn sanity_check_intensity_model(
    model: &scoring_crate::IntensityModel,
) -> Result<(), Box<dyn std::error::Error>> {
    use scoring_crate::IntensityIonType;

    // y after K/R should be brighter than b1 at N-terminus (when keys exist).
    let (y_kr, _) = model.predict_log_rel(IntensityIonType::Y, b'K', b'R', 5, 2, "25");
    let (b1, _) = model.predict_log_rel(IntensityIonType::B, b'A', b'L', 1, 2, "25");
    if y_kr <= b1 {
        eprintln!(
            "train-intensity: warning: y(K|R) mean {y_kr:.3} not above b1 mean {b1:.3} (sparse training data?)"
        );
    } else {
        eprintln!("train-intensity: sanity OK: y(K|R)={y_kr:.3} > b1={b1:.3}");
    }
    Ok(())
}

/// `andes train-intensity`: merge partial intensity aggregation parquets.
pub(crate) fn run_train_intensity(
    args: TrainIntensityArgs,
) -> Result<(), Box<dyn std::error::Error>> {
    use rustc_hash::FxHashMap;
    use scoring_crate::IntensityModel;

    let t0 = std::time::Instant::now();
    let mut merged: FxHashMap<IntensityAggKey, IntensityAggStats> = FxHashMap::default();
    let mut rows_read = 0usize;

    for input in &args.inputs {
        eprintln!("train-intensity: reading {} ...", input.display());
        let part = read_intensity_partial(input)?;
        rows_read += part.len();
        for (key, stats) in part {
            let slot = merged.entry(key).or_default();
            slot.count += stats.count;
            slot.sum_log_rel += stats.sum_log_rel;
            slot.sum_log_rel_sq += stats.sum_log_rel_sq;
        }
        eprintln!("train-intensity:   {} key rows", rows_read);
    }
    if merged.is_empty() {
        return Err("no intensity key rows read from any --in parquet".into());
    }

    let exact_keys = merged.len();
    add_backoff_marginals(&mut merged);
    eprintln!(
        "train-intensity: {} exact keys + {} backoff marginals = {} total",
        exact_keys,
        merged.len() - exact_keys,
        merged.len()
    );

    write_intensity_model(&args.out, &merged)?;
    eprintln!(
        "train-intensity: wrote {} keys -> {}",
        merged.len(),
        args.out.display()
    );

    let model = IntensityModel::load(&args.out)?;
    sanity_check_intensity_model(&model)?;
    eprintln!(
        "train-intensity: done in {:.1}s",
        t0.elapsed().as_secs_f64()
    );
    Ok(())
}

/// `andes train-intensity-gbdt`: fit a v3 GBDT fragment-intensity regressor
/// from externally-labeled PSM parquets and embed it in a Parquet model store.
///
/// The function reuses `read_msnet_parquet` / `load_seed_param` / `RankScorer`
/// from the `run_train` path and delegates the store write to
/// `write_all_models_with_sources_and_gbdt_pub`, preserving all other models.
pub(crate) fn run_train_intensity_gbdt(
    args: TrainIntensityGbdtArgs,
) -> Result<(), Box<dyn std::error::Error>> {
    use model_train::gbdt::dataset::PsmRow;
    use model_train::gbdt::frag_dataset::build_frag_dataset;
    use model_train::gbdt::train::{train_gbdt_regression, TrainParams};
    use std::sync::Arc;

    let n_files = args.inputs.len();
    let model_id = args.model_id.clone();
    let seed = args.seed_model.clone();
    let t = args.threads;
    eprintln!("train-intensity-gbdt: in={n_files} model_id={model_id} seed={seed} threads={t}");

    let t0 = std::time::Instant::now();

    // ── 1. Configure Rayon thread pool ────────────────────────────────────────
    static POOL_INIT_FRAG: std::sync::OnceLock<()> = std::sync::OnceLock::new();
    POOL_INIT_FRAG.get_or_init(|| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .expect("build_global");
    });

    // ── 2. Read all input parquets ────────────────────────────────────────────
    let mut psms: Vec<MsnetPsm> = Vec::new();
    for input in &args.inputs {
        eprintln!("train-intensity-gbdt: reading {} ...", input.display());
        let part = read_msnet_parquet(input)?;
        eprintln!("train-intensity-gbdt:   {} PSM rows", part.len());
        psms.extend(part);
    }
    if psms.is_empty() {
        return Err("no PSM rows read from any --in parquet".into());
    }
    eprintln!(
        "train-intensity-gbdt: {} total PSM rows across {} file(s)",
        psms.len(),
        args.inputs.len()
    );

    // ── 3. Load seed Param and build the scorer ───────────────────────────────
    let (seed_model_id, seed_param) = load_seed_param(&Some(args.seed_model.clone()))?;
    eprintln!("train-intensity-gbdt: seed model = {seed_model_id}");
    let seed_scorer = RankScorer::new(&seed_param);

    // ── 4. Build frag-intensity regression dataset ────────────────────────────
    eprintln!(
        "train-intensity-gbdt: building frag dataset from {} PSMs ...",
        psms.len()
    );
    let rows: Vec<PsmRow<'_>> = psms
        .iter()
        .map(|psm| PsmRow {
            spectrum: &psm.spectrum,
            peptide: &psm.peptide,
            charge: psm.charge,
        })
        .collect();
    let ds = build_frag_dataset(&rows, &seed_scorer);
    eprintln!(
        "train-intensity-gbdt: dataset: {} ion rows, {} features",
        ds.y.len(),
        ds.n_features,
    );

    // ── 5. Train the GBDT regressor ───────────────────────────────────────────
    // Hard-error if the trainer's quality gate fails (finding 3.6): a training
    // subcommand must not silently emit a non-deployable model, unless the
    // operator opts into the degenerate fallback.
    let train_params = TrainParams {
        allow_degenerate: args.allow_degenerate_model,
        ..TrainParams::default()
    };
    let trained_frag = train_gbdt_regression(&ds, &train_params, 42)?;
    eprintln!(
        "train-intensity-gbdt: trained frag model: {} trees",
        trained_frag.trees.len()
    );

    // ── 6. Embed the trained model in the seed Param ──────────────────────────
    // The seed Param supplies selection columns (activation/instrument/enzyme/
    // protocol) so model routing works; we stamp the frag model onto it.
    let mut out_param = seed_param;
    out_param.frag_intensity_model = Some(Arc::new(trained_frag));

    // ── 7. Write to store, preserving existing models ─────────────────────────
    let store_path = &args.out_store;

    {
        let mut existing_other: Vec<ModelEntryOwned> = Vec::new();
        let mut existing_blobs: Vec<Option<Vec<u8>>> = Vec::new();
        if store_path.exists() {
            let store = ModelStore::open(store_path)
                .map_err(|e| format!("opening existing store {}: {e}", store_path.display()))?;
            for id in store.model_ids() {
                if id == args.model_id {
                    eprintln!("train-intensity-gbdt: overwriting existing model '{id}' in store");
                    continue;
                }
                let p = store
                    .load_param(&id)
                    .map_err(|e| format!("reading model '{id}': {e}"))?;
                let blob = p.gbdt_peak_model.as_ref().map(|m| m.to_bytes());
                let src_ledgers = store.load_sources(&id).unwrap_or_default();
                let mut src = Vec::new();
                for l in src_ledgers {
                    if let Ok(s) = store.load_source_stats(&id, &l.source_id) {
                        src.push((l, s));
                    }
                }
                existing_other.push((id, p, src));
                existing_blobs.push(blob);
            }
        }

        // New model has no rank-core sources; sources slice is empty.
        let mut all_entries: Vec<ModelEntryOwned> = Vec::new();
        all_entries.push((args.model_id.clone(), out_param, vec![]));
        for (id, p, src) in existing_other {
            all_entries.push((id, p, src));
        }

        // New model carries no separate GBDT peak-model blob (the frag-intensity
        // model is embedded directly on Param.frag_intensity_model, not in the
        // gbdt_model_bytes column).  Existing models' blobs are preserved.
        let mut all_blobs: Vec<Option<Vec<u8>>> = vec![None];
        all_blobs.extend(existing_blobs);

        write_all_models_with_sources_and_gbdt_pub(
            store_path,
            &all_entries
                .iter()
                .map(|(id, p, s)| (id.as_str(), p, s.as_slice()))
                .collect::<Vec<_>>(),
            &all_blobs,
        )
        .map_err(|e| format!("writing model store {}: {e}", store_path.display()))?;
    }

    eprintln!(
        "train-intensity-gbdt: wrote model '{model_id}' to {} [{:.2}s]",
        store_path.display(),
        t0.elapsed().as_secs_f64(),
    );

    Ok(())
}

/// `andes train-rich-ion-llr`: fit a GBDT rich-ion LLR classifier (logistic;
/// decoy-aware) from externally-labeled PSM parquets and embed it in a Parquet
/// model store.
///
/// The function reuses `read_msnet_parquet` / `load_seed_param` / `RankScorer`
/// from the `run_train` path and delegates the store write to
/// `write_all_models_with_sources_and_gbdt_pub`, preserving all other models.
pub(crate) fn run_train_rich_ion_llr(
    args: TrainRichIonLlrArgs,
) -> Result<(), Box<dyn std::error::Error>> {
    use model_train::gbdt::dataset::PsmRow;
    use model_train::gbdt::ion_dataset::build_ion_dataset;
    use model_train::gbdt::train::{train_gbdt, TrainParams};
    use std::sync::Arc;

    let n_files = args.inputs.len();
    let model_id = args.model_id.clone();
    let seed = args.seed_model.clone();
    let t = args.threads;
    eprintln!("train-rich-ion-llr: in={n_files} model_id={model_id} seed={seed} threads={t}");

    let t0 = std::time::Instant::now();

    // ── 1. Configure Rayon thread pool ────────────────────────────────────────
    static POOL_INIT_RICH_ION: std::sync::OnceLock<()> = std::sync::OnceLock::new();
    POOL_INIT_RICH_ION.get_or_init(|| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .expect("build_global");
    });

    // ── 2. Read all input parquets ────────────────────────────────────────────
    let mut psms: Vec<MsnetPsm> = Vec::new();
    for input in &args.inputs {
        eprintln!("train-rich-ion-llr: reading {} ...", input.display());
        let part = read_msnet_parquet(input)?;
        eprintln!("train-rich-ion-llr:   {} PSM rows", part.len());
        psms.extend(part);
    }
    if psms.is_empty() {
        return Err("no PSM rows read from any --in parquet".into());
    }
    eprintln!(
        "train-rich-ion-llr: {} total PSM rows across {} file(s)",
        psms.len(),
        args.inputs.len()
    );

    // ── 3. Load seed Param and build the scorer ───────────────────────────────
    let (seed_model_id, seed_param) = load_seed_param(&Some(args.seed_model.clone()))?;
    eprintln!("train-rich-ion-llr: seed model = {seed_model_id}");
    let seed_scorer = RankScorer::new(&seed_param);

    // ── 4. Build rich-ion classification dataset ──────────────────────────────
    eprintln!(
        "train-rich-ion-llr: building ion dataset from {} PSMs ...",
        psms.len()
    );
    let rows: Vec<PsmRow<'_>> = psms
        .iter()
        .map(|psm| PsmRow {
            spectrum: &psm.spectrum,
            peptide: &psm.peptide,
            charge: psm.charge,
        })
        .collect();
    let ds = build_ion_dataset(&rows, &seed_scorer);
    eprintln!(
        "train-rich-ion-llr: dataset: {} ion rows, {} features",
        ds.y.len(),
        ds.n_features,
    );

    // ── 5. Train the GBDT classifier (logits held-out AUC) ────────────────────
    // Hard-error on quality-gate failure (finding 3.6) unless opted out.
    let train_params = TrainParams {
        allow_degenerate: args.allow_degenerate_model,
        ..TrainParams::default()
    };
    let trained_rich_ion = train_gbdt(&ds, &train_params, 42)?;
    eprintln!(
        "train-rich-ion-llr: trained rich-ion model: {} trees",
        trained_rich_ion.trees.len()
    );

    // ── 6. Embed the trained model in the seed Param ──────────────────────────
    // The seed Param supplies selection columns (activation/instrument/enzyme/
    // protocol) so model routing works; we stamp the rich-ion model onto it.
    let mut out_param = seed_param;
    out_param.rich_ion_model = Some(Arc::new(trained_rich_ion));

    // ── 7. Write to store, preserving existing models ─────────────────────────
    let store_path = &args.out_store;

    {
        let mut existing_other: Vec<ModelEntryOwned> = Vec::new();
        let mut existing_blobs: Vec<Option<Vec<u8>>> = Vec::new();
        if store_path.exists() {
            let store = ModelStore::open(store_path)
                .map_err(|e| format!("opening existing store {}: {e}", store_path.display()))?;
            for id in store.model_ids() {
                if id == args.model_id {
                    eprintln!("train-rich-ion-llr: overwriting existing model '{id}' in store");
                    continue;
                }
                let p = store
                    .load_param(&id)
                    .map_err(|e| format!("reading model '{id}': {e}"))?;
                let blob = p.gbdt_peak_model.as_ref().map(|m| m.to_bytes());
                let src_ledgers = store.load_sources(&id).unwrap_or_default();
                let mut src = Vec::new();
                for l in src_ledgers {
                    if let Ok(s) = store.load_source_stats(&id, &l.source_id) {
                        src.push((l, s));
                    }
                }
                existing_other.push((id, p, src));
                existing_blobs.push(blob);
            }
        }

        // New model has no rank-core sources; sources slice is empty.
        let mut all_entries: Vec<ModelEntryOwned> = Vec::new();
        all_entries.push((args.model_id.clone(), out_param, vec![]));
        for (id, p, src) in existing_other {
            all_entries.push((id, p, src));
        }

        // New model carries no separate GBDT peak-model blob (the rich-ion
        // model is embedded directly on Param.rich_ion_model, not in the
        // gbdt_model_bytes column).  Existing models' blobs are preserved.
        let mut all_blobs: Vec<Option<Vec<u8>>> = vec![None];
        all_blobs.extend(existing_blobs);

        write_all_models_with_sources_and_gbdt_pub(
            store_path,
            &all_entries
                .iter()
                .map(|(id, p, s)| (id.as_str(), p, s.as_slice()))
                .collect::<Vec<_>>(),
            &all_blobs,
        )
        .map_err(|e| format!("writing model store {}: {e}", store_path.display()))?;
    }

    eprintln!(
        "train-rich-ion-llr: wrote model '{model_id}' to {} [{:.2}s]",
        store_path.display(),
        t0.elapsed().as_secs_f64(),
    );

    Ok(())
}

#[cfg(test)]
mod intensity_marginal_tests {
    use super::*;
    use rustc_hash::FxHashMap;

    fn cell(
        ion: &str,
        fn_: &str,
        fc: &str,
        nce: &str,
        count: i64,
        sum: f64,
    ) -> (IntensityAggKey, IntensityAggStats) {
        (
            IntensityAggKey {
                ion_type: ion.to_string(),
                flank_n: fn_.to_string(),
                flank_c: fc.to_string(),
                pos_bin: 5,
                charge: 2,
                nce_bin: nce.to_string(),
            },
            IntensityAggStats {
                count,
                sum_log_rel: sum,
                sum_log_rel_sq: sum * sum / count as f64,
            },
        )
    }

    /// The finalizer must emit the backoff marginals that `predict_log_rel`
    /// probes, or the trained per-context intensities go unused at inference
    /// (which passes `nce_bin="unknown"`).
    #[test]
    fn backoff_marginals_are_emitted_with_summed_stats() {
        let mut merged: FxHashMap<IntensityAggKey, IntensityAggStats> = FxHashMap::default();
        // Two real cells: same flank/pos/charge, different numeric nce bins.
        for (k, s) in [
            cell("y", "K", "R", "30", 100, -20.0),
            cell("y", "K", "R", "40", 50, -15.0),
        ] {
            merged.insert(k, s);
        }
        let exact = merged.len();
        add_backoff_marginals(&mut merged);

        // nce-marginal (`__any__`) keeps flanks, sums across both nce bins.
        let any_nce = merged
            .get(&cell("y", "K", "R", ANY_NCE, 0, 0.0).0)
            .expect("nce-marginal must exist");
        assert_eq!(any_nce.count, 150, "summed over the two nce bins");
        assert!((any_nce.sum_log_rel - (-35.0)).abs() < 1e-9);

        // flank-marginal (`*`/`*`) keeps each nce bin separately.
        assert!(merged.contains_key(&cell("y", ANY_FLANK, ANY_FLANK, "30", 0, 0.0).0));
        // both-marginal exists too.
        assert!(merged.contains_key(&cell("y", ANY_FLANK, ANY_FLANK, ANY_NCE, 0, 0.0).0));

        // Real cells are untouched (no double counting).
        assert_eq!(
            merged
                .get(&cell("y", "K", "R", "30", 0, 0.0).0)
                .unwrap()
                .count,
            100
        );
        assert!(merged.len() > exact, "marginals added new keys");
    }
}
