//! Fragment-intensity GBDT and rich-ion LLR training.

use std::path::PathBuf;

use crate::model_select::load_seed_param;
use crate::train::{read_msnet_parquet, MsnetPsm};
use clap::{Args, ValueEnum};
use model_train::store::{attach_to_model, ensure_single_file_store};
use scoring_crate::RankScorer;

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

    /// Single-file Parquet model store to write (created if absent; other
    /// models are preserved). When `--model-id` already exists there, only its
    /// fragment-intensity GBDT is replaced; its rank tables, other GBDTs and source
    /// ledger are kept. REQUIRED.
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

    /// Single-file Parquet model store to write (created if absent; other
    /// models are preserved). When `--model-id` already exists there, only its
    /// rich-ion GBDT is replaced; its rank tables, other GBDTs and source
    /// ledger are kept. REQUIRED.
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

/// `andes train-intensity-gbdt`: fit a v3 GBDT fragment-intensity regressor
/// from externally-labeled PSM parquets and embed it in a Parquet model store.
///
/// The function reuses `read_msnet_parquet` / `load_seed_param` / `RankScorer`
/// from the `run_train` path and delegates the store write to
/// `attach_to_model`, preserving the target model's other content and all
/// other models.
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
    ensure_single_file_store(&args.out_store)?;

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

    // ── 6. Attach the trained model to --model-id in --out-store ─────────────
    // An existing model keeps its rank tables, other GBDTs and source ledger;
    // a new one starts from the seed Param, whose selection columns
    // (activation/instrument/enzyme/protocol) route it.
    let store_path = &args.out_store;
    let existed = attach_to_model(store_path, &model_id, seed_param, |p| {
        p.frag_intensity_model = Some(Arc::new(trained_frag));
    })
    .map_err(|e| format!("writing model store {}: {e}", store_path.display()))?;
    if existed {
        eprintln!("train-intensity-gbdt: replaced the frag model of existing '{model_id}'");
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
/// `attach_to_model`, preserving the target model's other content and all
/// other models.
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
    ensure_single_file_store(&args.out_store)?;

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

    // ── 6. Attach the trained model to --model-id in --out-store ─────────────
    // An existing model keeps its rank tables, other GBDTs and source ledger;
    // a new one starts from the seed Param, whose selection columns
    // (activation/instrument/enzyme/protocol) route it.
    let store_path = &args.out_store;
    let existed = attach_to_model(store_path, &model_id, seed_param, |p| {
        p.rich_ion_model = Some(Arc::new(trained_rich_ion));
    })
    .map_err(|e| format!("writing model store {}: {e}", store_path.display()))?;
    if existed {
        eprintln!("train-rich-ion-llr: replaced the rich-ion model of existing '{model_id}'");
    }

    eprintln!(
        "train-rich-ion-llr: wrote model '{model_id}' to {} [{:.2}s]",
        store_path.display(),
        t0.elapsed().as_secs_f64(),
    );

    Ok(())
}
