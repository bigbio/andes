//! Arguments of `andes train` and `andes train-from-search`.

use std::path::PathBuf;

use clap::Args;

use crate::cli::parse_unit_fraction;
use crate::train::intensity::GbdtMode;

/// Training arguments for `andes train-from-search`.
#[derive(Args, Debug)]
pub(crate) struct TrainFromSearchArgs {
    /// Reuse the seed model's geometry instead of deriving one from the corpus.
    /// Deriving (the default) fits segments and mass tiers to the training data.
    #[arg(long = "seed-geometry", default_value_t = false)]
    pub(crate) seed_geometry: bool,
    /// Input spectrum file (training data). Same format dispatch as for search:
    /// `.mzML`/`.mzml` → mzML reader; anything else → MGF reader.
    ///
    /// Repeat the flag to train one model from several runs. Each run is
    /// searched and labelled on its own (its own target-decoy q-values), keeps
    /// only its labelled spectra in memory, and becomes its own source in the
    /// store ledger (source id = file stem), so it can later be removed or
    /// reweighted with `--update`. `--labels` and `--update --add` take a
    /// single file.
    ///
    /// Required for initial training.  In `--update` mode with `--remove-source`
    /// or `--reweight` / `--decay`, `--spectra` is only required when
    /// `--validate` is also given (to run the acceptance gate).
    #[arg(long)]
    pub(crate) spectra: Vec<PathBuf>,

    /// Input FASTA target database (decoys are generated automatically).
    ///
    /// Required for initial training and for `--update --add`.
    /// In `--update` mode without `--add`, only required when `--validate` is
    /// given.
    #[arg(long)]
    pub(crate) database: Option<PathBuf>,

    /// EXTERNAL LABELS (SP-B glyco training): a TSV with columns `scan`,
    /// `peptide`, `charge`. When given, labels come from this file instead of a
    /// seed search (`bootstrap_labels`), and `--database` is not required. The
    /// `peptide` is the BARE backbone sequence (the glycan is stripped — the
    /// rank model scores backbone b/y ions, which in HCD come from the
    /// deglycosylated backbone); Cam-C and any `--fixed-mods` / `--variable-mods` are applied via the
    /// AminoAcidSet. Rows whose scan is absent from `--spectra` or whose peptide
    /// fails to parse are skipped. Use to train a glyco-regime rank model from
    /// MSFragger/pGlyco backbone IDs.
    ///
    /// With several `--spectra`, give one `--labels` per run, in the same order.
    #[arg(long = "labels")]
    pub(crate) labels: Vec<PathBuf>,

    /// Seed model: slug from the bundled store (e.g. `hcd_qexactive_tryp`).
    /// When omitted, the bundled `hcd_qexactive_tryp` model is used as the seed.
    #[arg(long = "seed-model")]
    pub(crate) seed_model: Option<String>,

    /// Target-decoy q-value threshold for accepting PSMs as confident training
    /// labels. Use a lenient value (e.g. 0.1 or 0.5) for small fixtures.
    #[arg(long = "train-fdr", default_value = "0.01", value_parser = parse_unit_fraction)]
    pub(crate) train_fdr: f64,

    /// Instrument of the trained model (`LowRes`, `HighRes`, `QExactive`,
    /// `OrbitrapAstral`, `TOF`, `TimsTOF`). It is part of the model's selection
    /// key, so a search on data from this instrument picks the model. Default:
    /// the seed model's instrument (`QExactive` for the default seed).
    #[arg(long)]
    pub(crate) instrument: Option<String>,

    /// Activation of the trained model (`CID`, `HCD`, `ETD`, `EThcD`, `UVPD`, `PQD`),
    /// part of its selection key. It also sets the ion series of the derived geometry
    /// (EThcD: b, y, c and z•), and only spectra recorded with this activation are
    /// searched and trained on (spectra without a recorded activation are kept).
    /// Default: the seed model's activation, all spectra.
    #[arg(long, conflicts_with = "seed_geometry")]
    pub(crate) activation: Option<String>,

    /// Experiment-class / protocol tag (e.g. `Automatic`, `TMT`). Default: `Automatic`.
    #[arg(long, default_value = "Automatic")]
    pub(crate) protocol: String,

    /// Path to the Parquet model store to write (created if absent, appended
    /// otherwise). REQUIRED.
    #[arg(long = "out-store")]
    pub(crate) out_store: PathBuf,

    /// Directory of the candidate-index cache of the seed searches, created if
    /// missing. Default: the directory of `--out-store`. A directory that
    /// cannot be written is an error.
    #[arg(long = "index-cache-dir", value_name = "DIR")]
    pub(crate) index_cache_dir: Option<PathBuf>,

    /// Model ID written into the store. Default: `trained_<instrument>_<protocol>`.
    #[arg(long = "model-id")]
    pub(crate) model_id: Option<String>,

    /// Fixed modifications, a list of Name (sites) as for the search's
    /// `--fixed-mods` [default: Carbamidomethyl (C)].
    #[arg(long, value_name = "LIST", value_parser = crate::cli::parse_fixed_mods)]
    pub(crate) fixed_mods: Option<crate::cli::ModList>,

    /// Variable modifications, a list of Name (sites) as for the search's
    /// `--variable-mods` [default: Oxidation (M)].
    #[arg(long, value_name = "LIST", value_parser = crate::cli::parse_variable_mods)]
    pub(crate) variable_mods: Option<crate::cli::ModList>,

    /// Number of worker threads. Defaults to logical CPU count.
    #[arg(long, default_value_t = num_cpus::get())]
    pub(crate) threads: usize,

    /// ISO 8601 date string (e.g. `2026-01-01`) recorded in the source ledger.
    /// When omitted, the current date is used for initial training; empty string
    /// is stored when `--date ""` is explicitly passed.
    #[arg(long)]
    pub(crate) date: Option<String>,

    /// Fragment match tolerance in ppm for initial training: the window in
    /// which training matches fragments, stored as the trained model's `mme`.
    /// The seed search that labels the PSMs keeps the seed's `mme`. When
    /// neither this nor `--fragment-tol-da` is given, the seed's `mme` is
    /// kept. Training reconstructs fragment m/z from nominal masses, so a ppm
    /// window this narrow marks most fragments as missing; keep the seed `mme`
    /// unless the comparison is the point.
    #[arg(long = "fragment-tol-ppm", conflicts_with = "fragment_tol_da")]
    pub(crate) fragment_tol_ppm: Option<f64>,

    /// Fragment match tolerance in Da for initial training. Overwrites the seed
    /// model's `mme` (see `--fragment-tol-ppm`).
    #[arg(long = "fragment-tol-da")]
    pub(crate) fragment_tol_da: Option<f64>,

    /// Switch to incremental update mode for this model ID.
    /// When set, one of `--add`, `--remove-source`, `--reweight`, or `--decay`
    /// must be provided.
    #[arg(long = "update", value_name = "MODEL_ID")]
    pub(crate) update_model: Option<String>,

    /// (Update mode) Add a new source from `--spectra`.
    /// Requires `--source-id` and `--database`.
    #[arg(long, requires = "update_model")]
    pub(crate) add: bool,

    /// (Update mode) Source identifier for the new source being added
    /// (used with `--add`).
    #[arg(long = "source-id", requires = "add", value_name = "ID")]
    pub(crate) source_id: Option<String>,

    /// (Update mode) Remove the source with this ID from the model.
    #[arg(long = "remove-source", requires = "update_model", value_name = "ID")]
    pub(crate) remove_source: Option<String>,

    /// (Update mode) Set a source's weight.  Format: `<source-id>=<weight>`,
    /// e.g. `--reweight s0=0.5`.
    #[arg(long = "reweight", requires = "update_model", value_name = "ID=W")]
    pub(crate) reweight: Option<String>,

    /// (Update mode) Apply exponential age-decay to all sources with this
    /// half-life in days.
    #[arg(long = "decay", requires = "update_model", value_name = "DAYS")]
    pub(crate) decay: Option<f32>,

    /// (Update mode) Held-out validation spectra for the acceptance gate.
    /// When omitted the gate is skipped (a warning is printed).
    #[arg(long = "validate", requires = "update_model")]
    pub(crate) validate: Option<PathBuf>,

    /// (Update mode) Commit the update even if the acceptance gate fails.
    #[arg(long, requires = "update_model")]
    pub(crate) force: bool,
}

/// Training arguments for `andes train`.
///
/// Trains a scoring model directly from externally-labeled, high-confidence
/// PSMs supplied as a "flat training parquet" (one row per PSM, each carrying
/// the spectrum peaks + identified peptide + resolved mod mass-deltas). This
/// bypasses the bootstrap-search label step entirely: every input row is a
/// label. The seed model supplies only structural hyperparameters (`mme`,
/// deconvolution, segments, frag/precursor offset tables, max_rank); all
/// learned distributions come from the input data.
#[derive(Args, Debug)]
pub(crate) struct TrainArgs {
    /// Reuse the seed model's geometry instead of deriving one from the corpus.
    /// Deriving (the default) fits segments and mass tiers to the training data.
    #[arg(long = "seed-geometry", default_value_t = false)]
    pub(crate) seed_geometry: bool,
    /// Input flat training parquet(s). Repeatable; stats accumulate across all
    /// inputs into a single model.
    #[arg(long = "in", required = true)]
    pub(crate) inputs: Vec<PathBuf>,

    /// Path to the Parquet model store to write (created if absent; existing
    /// models are preserved and re-written alongside the new one). REQUIRED.
    #[arg(long = "out-store")]
    pub(crate) out_store: PathBuf,

    /// Model ID written into the store. Default: `default`.
    #[arg(long = "model-id", default_value = "default")]
    pub(crate) model_id: String,

    /// Seed model: slug from the bundled store (e.g. `hcd_qexactive_tryp`).
    /// Supplies structural hyperparameters only.
    #[arg(long = "seed-model", default_value = "hcd_qexactive_tryp")]
    pub(crate) seed_model: String,

    /// Override the trained model's activation in the store `data_type`
    /// (e.g. `CID`, `HCD`, `ETD`, `UVPD`). Defaults to the seed's value.
    /// Together with `--instrument/--enzyme/--protocol` this lets a new slug
    /// carry the correct selection columns even when seeded from a related model.
    #[arg(long = "activation")]
    pub(crate) activation: Option<String>,

    /// Override the trained model's instrument in the store `data_type`
    /// (e.g. `LowRes`, `HighRes`, `QExactive`, `TOF`). Defaults to the seed's value.
    #[arg(long = "instrument")]
    pub(crate) instrument: Option<String>,

    /// Override the trained model's enzyme in the store `data_type`
    /// (e.g. `Trypsin`, `LysC`, `LysN`). Defaults to the seed's value.
    #[arg(long = "enzyme")]
    pub(crate) enzyme: Option<String>,

    /// Override the trained model's protocol in the store `data_type`
    /// (e.g. `TMT`, `iTRAQ`, `Phosphorylation`, `Automatic`). Drives
    /// `experiment_class` model selection. Defaults to the seed's value.
    #[arg(long = "protocol")]
    pub(crate) protocol: Option<String>,

    /// Fragment match tolerance in ppm. Overwrites the seed model's `mme`
    /// before training. Mutually exclusive with `--fragment-tol-da`. When
    /// neither is given, the seed model's `mme` is kept.
    #[arg(long = "fragment-tol-ppm", conflicts_with = "fragment_tol_da")]
    pub(crate) fragment_tol_ppm: Option<f64>,

    /// Fragment match tolerance in Da. Overwrites the seed model's `mme`
    /// before training. Mutually exclusive with `--fragment-tol-ppm`.
    #[arg(long = "fragment-tol-da")]
    pub(crate) fragment_tol_da: Option<f64>,

    /// Number of worker threads. Defaults to logical CPU count.
    #[arg(long, default_value_t = num_cpus::get())]
    pub(crate) threads: usize,

    /// Laplace pseudo-count for rank/error tables (lower = sharper; default 1.0).
    #[arg(long = "train-pseudo", default_value_t = 1.0)]
    pub(crate) train_pseudo: f32,

    /// Laplace pseudo-count for the NOISE rank distribution (lower = sharper).
    /// Noise is abundant and concentrated, so it needs far less smoothing than
    /// signal ions; the signal `--train-pseudo` over-flattens it. Default 0.05.
    #[arg(long = "train-noise-pseudo", default_value_t = 0.05)]
    pub(crate) train_noise_pseudo: f32,

    /// Partition backoff prior weight (lower = less smoothing toward parent; default 20).
    #[arg(long = "train-backoff-weight", default_value_t = 20.0)]
    pub(crate) train_backoff_weight: f32,

    /// Minimum partition count before backoff blending (default 50).
    #[arg(long = "train-min-count", default_value_t = 50)]
    pub(crate) train_min_count: u64,

    /// Optional path to an independent prior model store. Sparse partitions in
    /// the trained model shrink toward the matching prior model instead of the
    /// corpus-internal pool. Must be own-data (NOT a bundled seed model) to stay
    /// relicense-safe.
    #[arg(long)]
    pub(crate) prior_model_store: Option<PathBuf>,

    /// Model id to load from `--prior-model-store` (defaults to the trained
    /// model id when omitted).
    #[arg(long)]
    pub(crate) prior_model: Option<String>,

    /// Apply widening rank-window smoothing to signal rank distributions
    /// (Kim et al., Nat Commun 5:5277, 2014).
    #[arg(long)]
    pub(crate) rank_smoothing: bool,

    /// Source identifier for the source ledger. Defaults to "msnet".
    #[arg(long, default_value = "msnet")]
    pub(crate) source: String,

    /// Whether to also train and embed a GBDT peak model. `on` (default) trains
    /// GBDT and writes the blob; `off` writes rank-core only (byte-identical to
    /// the pre-GBDT path).
    #[arg(long, default_value = "on")]
    pub(crate) gbdt: GbdtMode,

    /// Opt-in fallback (finding 3.6): downgrade a failed GBDT quality gate to a
    /// warning and embed the degenerate model anyway. Default off.
    #[arg(long = "allow-degenerate-model", hide = true, default_value_t = false)]
    pub(crate) allow_degenerate_model: bool,
}
