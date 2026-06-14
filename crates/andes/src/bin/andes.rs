//! andes: end-to-end peptide-spectrum database search.
//!
//! Loads an MGF or mzML spectrum file and a FASTA target database, runs a
//! tryptic database search and writes output
//! in Percolator `.pin` format (and optionally `.tsv` format).
//!
//! Format dispatch by `--spectrum` extension: `.mzML`/`.mzml` → `MzMLReader`;
//! `.d` → `TimsTofReader` (native Bruker timsTOF, only under `--features
//! timstof`); otherwise `MgfReader` (default / backwards-compatible).

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::mpsc::{sync_channel, SyncSender};
use std::sync::Arc;
use std::thread;

use clap::{Args, Parser, Subcommand, ValueEnum};
use model::{
    activation::ActivationMethod, AminoAcidSetBuilder, InstrumentType, ModLocation, Modification,
    PrecursorTolerance, ResidueSpec, Spectrum, Tolerance,
};
use model_train::{
    ModelStore,
    accumulate::{merge, StatsAccumulator},
    counts::CountStats,
    estimate::{Estimator, EstimatorConfig},
    gate::evaluate_candidate,
    labeled::bootstrap_labels,
    select::{select, SelectionKey},
    protocol_to_experiment_class as store_protocol_to_experiment_class,
    store::{
        SourceLedger,
        update_add, update_remove, update_reweight, update_decay, commit_update,
        write_all_models_with_sources_pub,
    },
};
use scoring_crate::{Param, RankScorer};
use search::{
    apply_shift_for_mode, apply_tightened_precursor_tolerance, build_spec_keys,
    learn_calibration_stats, CalibrationStats,
    PreparedSearch, PrecursorCalMode, SearchIndex, SearchParams, SpecKey, TopNQueue,
};
use search::precursor_cal::{constants as cal_constants, sample_every_nth};
use input::{detect_instrument_type, FastaReader, MgfReader, Ms1Link, MzMLReader};

// Type alias to reduce clippy type_complexity warnings in the train path.
type ModelEntryOwned = (String, Param, Vec<(SourceLedger, CountStats)>);

/// Fragmentation method. `Auto` detects from the mzML's activation block and
/// falls back to the bundled `HCD_QExactive_Tryp.param` when nothing is detected.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum Fragmentation {
    #[clap(name = "auto")] Auto,
    #[clap(name = "CID")]  Cid,
    #[clap(name = "ETD")]  Etd,
    #[clap(name = "HCD")]  Hcd,
    #[clap(name = "UVPD")] Uvpd,
}

/// Search protocol: sample labeling or enrichment strategy applied during the experiment.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum Protocol {
    #[clap(name = "auto")]          Auto,
    #[clap(name = "phospho")]       Phospho,
    #[clap(name = "iTRAQ")]         Itraq,
    #[clap(name = "iTRAQ-phospho")] ItraqPhospho,
    #[clap(name = "TMT")]           Tmt,
    #[clap(name = "standard")]      Standard,
}

/// Enzymatic-cleavage enforcement at peptide span boundaries:
/// 2=fully, 1=semi, 0=non-specific.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum EnzymeSpecificity {
    #[clap(name = "non-specific")] NonSpecific,
    #[clap(name = "semi")]         Semi,
    #[clap(name = "fully")]        Fully,
}

/// Primary ranking mode: inherited RawScore (`rank`) or fused strong score (`strong`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum)]
enum ScoreFlag {
    #[default]
    Rank,
    Strong,
}

/// Search arguments (shared by the default search path and exposed as a
/// flat arg group so that `andes --spectrum X --database Y --output-pin Z`
/// keeps working unchanged).
///
/// Note: `spectrum`, `database`, and `output_pin` are declared `Option<PathBuf>`
/// at the clap level so that they are not required when a subcommand (e.g.
/// `train`) is given.  When no subcommand is present, `run()` validates them
/// manually and returns an early error if they are missing.
#[derive(Args, Debug)]
struct SearchArgs {
    /// Input spectrum file(s). Repeat `--spectrum` for multiple inputs (one PIN).
    /// Format is auto-detected per file by extension.
    #[arg(long)]
    spectrum: Vec<PathBuf>,

    /// Input FASTA database (target sequences only; decoys are generated automatically).
    #[arg(long)]
    database: Option<PathBuf>,

    /// Output Percolator PIN file path.
    #[arg(long)]
    output_pin: Option<PathBuf>,

    /// Output TSV file path (optional).
    #[arg(long)]
    output_tsv: Option<PathBuf>,

    /// Decoy prefix used when generating reversed decoy sequences.
    #[arg(long, default_value = "XXX_")]
    decoy_prefix: String,

    /// Minimum isotope-error offset to try.
    #[arg(long, default_value = "-1")]
    isotope_error_min: i8,

    /// Maximum isotope-error offset to try.
    #[arg(long, default_value = "2")]
    isotope_error_max: i8,

    /// Precursor-mass calibration: `off`, `auto`, or `on`. `auto`/`on` learn a
    /// systematic ppm shift from confident PSMs in a pre-pass and tighten the
    /// precursor tolerance for the main search; `auto` skips the correction when
    /// the sample is too small to be reliable.
    #[arg(long = "precursor-cal", default_value = "auto", value_parser = parse_precursor_cal)]
    precursor_cal: PrecursorCalMode,

    /// Precursor mass tolerance in ppm.
    #[arg(long, default_value = "20.0")]
    precursor_tol_ppm: f64,

    /// Precursor tolerance in Da (overrides --precursor-tol-ppm; for low-res
    /// precursor selection). The asymmetric ppm flags below take precedence.
    #[arg(long = "precursor-tol-da")]
    precursor_tol_da: Option<f64>,

    /// Asymmetric precursor tolerance, left (lower) window in ppm. Requires
    /// --precursor-tol-right-ppm; together they override the symmetric forms
    /// (for a known systematic precursor offset).
    #[arg(long = "precursor-tol-left-ppm")]
    precursor_tol_left_ppm: Option<f64>,

    /// Asymmetric precursor tolerance, right (upper) window in ppm. Requires
    /// --precursor-tol-left-ppm.
    #[arg(long = "precursor-tol-right-ppm")]
    precursor_tol_right_ppm: Option<f64>,

    /// Minimum precursor charge to try when not specified in the spectrum.
    #[arg(long, default_value = "2")]
    charge_min: u8,

    /// Maximum precursor charge to try when not specified in the spectrum.
    #[arg(long, default_value = "5")]
    charge_max: u8,

    /// Maximum number of PSMs to retain per spectrum.
    #[arg(long, default_value = "10")]
    top_n: u32,

    /// Number of Tolerable Termini (enzymatic-cleavage enforcement at span
    /// boundaries). `fully`: both termini must be cleavage sites (strict,
    /// equivalent to legacy `-ntt 2`). `semi`: at least one terminus must be a
    /// cleavage site (legacy `-ntt 1`). `non-specific`: neither terminus needs
    /// to be a cleavage site (legacy `-ntt 0`). Legacy numeric 0/1/2 still accepted.
    #[arg(long = "enzyme-specificity", alias = "ntt",
          default_value = "fully", value_parser = parse_enzyme_specificity)]
    enzyme_specificity: EnzymeSpecificity,

    /// Proteolytic enzyme for in-silico digestion. Named values: trypsin
    /// (default), chymotrypsin, lysc, aspn, gluc, lysn, argc, alphalp,
    /// nocleavage, nonspecific. A wrong enzyme yields ~no PSMs (fails loud,
    /// not silent). Previously hardcoded to trypsin with no override.
    #[arg(long, default_value = "trypsin")]
    enzyme: String,

    /// Maximum number of missed cleavages per peptide.
    #[arg(long, default_value = "1")]
    max_missed_cleavages: u32,

    /// Minimum number of peaks an MS2 spectrum must have to be scored; spectra
    /// with fewer peaks are skipped.
    #[arg(long, default_value = "10")]
    min_peaks: u32,

    /// Minimum peptide length, in residues.
    #[arg(long, default_value = "6")]
    min_length: u32,

    /// Maximum peptide length, in residues. (50 matches MSFragger/Comet defaults;
    /// 40 dropped long tryptic peptides.)
    #[arg(long, default_value = "50")]
    max_length: u32,

    /// Maximum number of variable modifications per peptide. A `NumMods=N` line
    /// in a --mods file overrides this.
    #[arg(long = "max-mods", default_value = "3")]
    max_mods: u32,

    /// Path to the .param scoring model file.
    ///
    /// If not supplied, a scoring model is selected from the bundled
    /// `models.parquet` store. For mzML/.raw/.d the activation method and
    /// analyzer resolution are auto-detected from metadata; for MGF the
    /// `--fragmentation` and `--fragment-tol-ppm/-da` flags drive selection
    /// (default: CID / low-res). When running the binary outside the source
    /// tree the bundled store may not exist; supply --param-file explicitly
    /// in that case.
    #[arg(long)]
    param_file: Option<PathBuf>,

    /// Path to a mods.txt file describing fixed and variable modifications.
    /// Format: each non-comment line is
    /// `<mass>,<aa>,<fix|opt>,<location>,<name>`, where:
    ///   - `<mass>` is a numeric monoisotopic mass delta (Da). Composition
    ///     strings (e.g. `C2H3N1O1`) are **not** yet supported.
    ///   - `<aa>` is a single uppercase letter or `*` (wildcard).
    ///   - `<location>` is one of `any|N-term|C-term|Prot-N-term|Prot-C-term`.
    ///
    /// A single `NumMods=N` line sets the max variable mods per peptide.
    /// Inline `#`-comments are stripped. Blank lines and full-line `#`-comments
    /// are ignored. When omitted, the binary uses its built-in defaults
    /// (Carbamidomethyl-C fixed, Oxidation-M variable). The deprecated
    /// `--mod` form (singular) is still accepted as a hidden alias.
    #[arg(long = "mods", alias = "mod", value_name = "MODFILE")]
    mods: Option<PathBuf>,

    /// Fragmentation/activation method for MGF input only. mzML/.raw/.d
    /// auto-detect this. Named values: auto, CID, ETD, HCD, UVPD.
    /// Legacy numeric CLI indices: 0=auto, 1=CID, 2=ETD, 3=HCD, 4=UVPD.
    #[arg(long, hide = true, default_value = "auto", value_parser = parse_fragmentation)]
    fragmentation: Fragmentation,

    /// Search protocol. Named values: auto, phospho, iTRAQ, iTRAQ-phospho, TMT, standard.
    /// Legacy numeric CLI indices: 0=auto, 1=phospho, 2=iTRAQ, 3=iTRAQ-phospho, 4=TMT, 5=standard.
    #[arg(long, default_value = "auto", value_parser = parse_protocol)]
    protocol: Protocol,

    /// Fragment-matching tolerance in ppm for **MGF input only** (high-resolution
    /// MS/MS). Has no effect on mzML/.raw/.d (analyzer auto-detected). Mutually
    /// exclusive with `--fragment-tol-da`.
    #[arg(long = "fragment-tol-ppm", hide = true, conflicts_with = "fragment_tol_da")]
    fragment_tol_ppm: Option<f64>,

    /// Fragment-matching tolerance in Da for **MGF input only** (low-resolution
    /// ion-trap MS/MS). Has no effect on mzML/.raw/.d. Mutually exclusive with
    /// `--fragment-tol-ppm`.
    #[arg(long = "fragment-tol-da", hide = true, conflicts_with = "fragment_tol_ppm")]
    fragment_tol_da: Option<f64>,

    /// Number of worker threads for the search loop. Defaults to logical CPU count.
    #[arg(long, default_value_t = num_cpus::get())]
    threads: usize,

    /// Debug/benchmark cap: process only the first N spectra (0 = no cap).
    #[arg(long, default_value = "0")]
    max_spectra: usize,

    /// MS level to search. Defaults to MS2 (identification); MS1 and any higher
    /// levels (e.g. TMT SPS-MS3 reporter-quant scans) are filtered out at load
    /// time so they never enter the search loop. Override only if you explicitly
    /// want a different level. Applies to mzML and Thermo `.raw`; MGF files do
    /// not encode MS level and are always treated as MS2. The chimeric cascade
    /// always searches MS2 (it pairs MS2 with its preceding MS1).
    #[arg(long, default_value = "2")]
    ms_level: u8,

    /// Enable the two-pass chimeric cascade for co-isolated (co-fragmented)
    /// peptides. Pass 1 is the normal top-1 search; Pass 2 detects co-isolated
    /// precursors in each scan's MS1 isolation window and runs a targeted search
    /// for the second peptide on the residual spectrum, emitting it as an extra
    /// PSM. Requires mzML (MS1 scans); has no effect on MGF input.
    #[arg(long, default_value = "false")]
    chimeric: bool,

    /// Chimeric mode: fallback isolation half-width in Da when the mzML omits the
    /// per-scan isolation-window offsets.
    #[arg(long, default_value = "1.5")]
    isolation_halfwidth: f64,

    /// Path to a Parquet model store to use instead of the bundled
    /// `resources/ionstat/models.parquet`. When set, model selection reads from
    /// this store; when unset, the bundled store is used.
    #[arg(long = "model-store")]
    model_store: Option<PathBuf>,

    /// Exact model ID to load from the model store (bundled or `--model-store`).
    /// When set, skips automatic selection (metadata detection / `--fragmentation`
    /// / `--protocol`) and loads this ID directly. Useful after `andes train`
    /// to search with the freshly-trained model.
    #[arg(long = "model")]
    model_id_override: Option<String>,

    /// Path to a trained intensity model parquet (`andes train-intensity` output).
    /// Populates the additive `IntensitySignal` PIN column; ranking stays on RawScore
    /// until `--score strong` is enabled in a later phase. When unset, the column is 0.0.
    #[arg(long = "intensity-model")]
    intensity_model: Option<PathBuf>,

    /// Ranking / PIN RawScore source: `rank` (default, byte-identical) or `strong`
    /// (fused intensity + competition score from S1–S3).
    #[arg(long = "score", default_value = "rank")]
    score: ScoreFlag,
}

/// Training arguments for `andes train`.
#[derive(Args, Debug)]
struct TrainArgs {
    /// Input spectrum file (training data). Same format dispatch as for search:
    /// `.mzML`/`.mzml` → mzML reader; anything else → MGF reader.
    ///
    /// Required for initial training.  In `--update` mode with `--remove-source`
    /// or `--reweight` / `--decay`, `--spectra` is only required when
    /// `--validate` is also given (to run the acceptance gate).
    #[arg(long)]
    spectra: Option<PathBuf>,

    /// Input FASTA target database (decoys are generated automatically).
    ///
    /// Required for initial training and for `--update --add`.
    /// In `--update` mode without `--add`, only required when `--validate` is
    /// given.
    #[arg(long)]
    database: Option<PathBuf>,

    /// Seed model: slug from the bundled store (e.g. `hcd_qexactive_tryp`) or
    /// a path to a binary `.param` file. When omitted, the bundled
    /// `hcd_qexactive_tryp` model is used as the seed.
    #[arg(long = "seed-model")]
    seed_model: Option<String>,

    /// Target-decoy q-value threshold for accepting PSMs as confident training
    /// labels. Use a lenient value (e.g. 0.1 or 0.5) for small fixtures.
    #[arg(long = "train-fdr", default_value = "0.01")]
    train_fdr: f64,

    /// Instrument tag to embed in the trained model's metadata. Default: `QExactive`.
    #[arg(long, default_value = "QExactive")]
    instrument: String,

    /// Experiment-class / protocol tag (e.g. `Automatic`, `TMT`). Default: `Automatic`.
    #[arg(long, default_value = "Automatic")]
    protocol: String,

    /// Path to the Parquet model store to write (created if absent, appended
    /// otherwise). REQUIRED.
    #[arg(long = "out-store")]
    out_store: PathBuf,

    /// Model ID written into the store. Default: `trained_<instrument>_<protocol>`.
    #[arg(long = "model-id")]
    model_id: Option<String>,

    /// Path to a mods.txt file (same format as `--mods` for search). When
    /// omitted, uses built-in defaults (Carbamidomethyl-C fixed, Oxidation-M
    /// variable).
    #[arg(long)]
    mods: Option<PathBuf>,

    /// Number of worker threads. Defaults to logical CPU count.
    #[arg(long, default_value_t = num_cpus::get())]
    threads: usize,

    /// ISO 8601 date string (e.g. `2026-01-01`) recorded in the source ledger.
    /// When omitted, the current date is used for initial training; empty string
    /// is stored when `--date ""` is explicitly passed.
    #[arg(long)]
    date: Option<String>,

    // ── Update mode ──────────────────────────────────────────────────────────

    /// Switch to incremental update mode for this model ID.
    /// When set, one of `--add`, `--remove-source`, `--reweight`, or `--decay`
    /// must be provided.
    #[arg(long = "update", value_name = "MODEL_ID")]
    update_model: Option<String>,

    /// (Update mode) Add a new source from `--spectra`.
    /// Requires `--source-id` and `--database`.
    #[arg(long, requires = "update_model")]
    add: bool,

    /// (Update mode) Source identifier for the new source being added
    /// (used with `--add`).
    #[arg(long = "source-id", requires = "add", value_name = "ID")]
    source_id: Option<String>,

    /// (Update mode) Remove the source with this ID from the model.
    #[arg(long = "remove-source", requires = "update_model", value_name = "ID")]
    remove_source: Option<String>,

    /// (Update mode) Set a source's weight.  Format: `<source-id>=<weight>`,
    /// e.g. `--reweight s0=0.5`.
    #[arg(long = "reweight", requires = "update_model", value_name = "ID=W")]
    reweight: Option<String>,

    /// (Update mode) Apply exponential age-decay to all sources with this
    /// half-life in days.
    #[arg(long = "decay", requires = "update_model", value_name = "DAYS")]
    decay: Option<f32>,

    /// (Update mode) Held-out validation spectra for the acceptance gate.
    /// When omitted the gate is skipped (a warning is printed).
    #[arg(long = "validate", requires = "update_model")]
    validate: Option<PathBuf>,

    /// (Update mode) Commit the update even if the acceptance gate fails.
    #[arg(long, requires = "update_model")]
    force: bool,
}

/// Training arguments for `andes train-from-msnet`.
///
/// Trains a scoring model directly from externally-labeled, high-confidence
/// PSMs supplied as a "flat training parquet" (one row per PSM, each carrying
/// the spectrum peaks + identified peptide + resolved mod mass-deltas). This
/// bypasses the bootstrap-search label step entirely: every input row is a
/// label. The seed model supplies only structural hyperparameters (`mme`,
/// deconvolution, segments, frag/precursor offset tables, max_rank); all
/// learned distributions come from the input data.
#[derive(Args, Debug)]
struct TrainFromMsnetArgs {
    /// Input flat training parquet(s). Repeatable; stats accumulate across all
    /// inputs into a single model.
    #[arg(long = "in", required = true)]
    inputs: Vec<PathBuf>,

    /// Path to the Parquet model store to write (created if absent; existing
    /// models are preserved and re-written alongside the new one). REQUIRED.
    #[arg(long = "out-store")]
    out_store: PathBuf,

    /// Model ID written into the store. Default: `default`.
    #[arg(long = "model-id", default_value = "default")]
    model_id: String,

    /// Seed model: slug from the bundled store (e.g. `hcd_qexactive_tryp`) or a
    /// path to a binary `.param` file. Supplies structural hyperparameters only.
    #[arg(long = "seed-model", default_value = "hcd_qexactive_tryp")]
    seed_model: String,

    /// Override the trained model's activation in the store `data_type`
    /// (e.g. `CID`, `HCD`, `ETD`, `UVPD`). Defaults to the seed's value.
    /// Together with `--instrument/--enzyme/--protocol` this lets a new slug
    /// carry the correct selection columns even when seeded from a related model.
    #[arg(long = "activation")]
    activation: Option<String>,

    /// Override the trained model's instrument in the store `data_type`
    /// (e.g. `LowRes`, `HighRes`, `QExactive`, `TOF`). Defaults to the seed's value.
    #[arg(long = "instrument")]
    instrument: Option<String>,

    /// Override the trained model's enzyme in the store `data_type`
    /// (e.g. `Trypsin`, `LysC`, `LysN`). Defaults to the seed's value.
    #[arg(long = "enzyme")]
    enzyme: Option<String>,

    /// Override the trained model's protocol in the store `data_type`
    /// (e.g. `TMT`, `iTRAQ`, `Phosphorylation`, `Automatic`). Drives
    /// `experiment_class` model selection. Defaults to the seed's value.
    #[arg(long = "protocol")]
    protocol: Option<String>,

    /// Fragment match tolerance in ppm. Overwrites the seed model's `mme`
    /// before training. Mutually exclusive with `--fragment-tol-da`. When
    /// neither is given, the seed model's `mme` is kept.
    #[arg(long = "fragment-tol-ppm", conflicts_with = "fragment_tol_da")]
    fragment_tol_ppm: Option<f64>,

    /// Fragment match tolerance in Da. Overwrites the seed model's `mme`
    /// before training. Mutually exclusive with `--fragment-tol-ppm`.
    #[arg(long = "fragment-tol-da")]
    fragment_tol_da: Option<f64>,

    /// Number of worker threads. Defaults to logical CPU count.
    #[arg(long, default_value_t = num_cpus::get())]
    threads: usize,

    /// Laplace pseudo-count for rank/error tables (lower = sharper; default 1.0).
    #[arg(long = "train-pseudo", default_value_t = 1.0)]
    train_pseudo: f32,

    /// Laplace pseudo-count for the NOISE rank distribution (lower = sharper).
    /// Noise is abundant and concentrated, so it needs far less smoothing than
    /// signal ions; the signal `--train-pseudo` over-flattens it. Default 0.05.
    #[arg(long = "train-noise-pseudo", default_value_t = 0.05)]
    train_noise_pseudo: f32,

    /// Partition backoff prior weight (lower = less smoothing toward parent; default 20).
    #[arg(long = "train-backoff-weight", default_value_t = 20.0)]
    train_backoff_weight: f32,

    /// Minimum partition count before backoff blending (default 50).
    #[arg(long = "train-min-count", default_value_t = 50)]
    train_min_count: u64,

    /// Optional path to an independent prior model store. Sparse partitions in
    /// the trained model shrink toward the matching prior model instead of the
    /// corpus-internal pool. Must be own-data (NOT a bundled seed model) to stay
    /// relicense-safe.
    #[arg(long)]
    prior_model_store: Option<PathBuf>,

    /// Model id to load from `--prior-model-store` (defaults to the trained
    /// model id when omitted).
    #[arg(long)]
    prior_model: Option<String>,

    /// Apply widening rank-window smoothing to signal rank distributions
    /// (Kim et al., Nat Commun 5:5277, 2014).
    #[arg(long)]
    rank_smoothing: bool,
}

/// Training arguments for `andes train-intensity`.
///
/// Merges one or more partial intensity aggregation parquets (from
/// `msnet_intensity_agg.py`) into a finalized `intensity_model.parquet` with
/// `mean_log_rel` / `var_log_rel` columns for runtime lookup.
#[derive(Args, Debug)]
struct TrainIntensityArgs {
    /// Input partial or finalized intensity parquets. Repeatable; stats merge
    /// across all inputs.
    #[arg(long = "in", required = true)]
    inputs: Vec<PathBuf>,

    /// Output path for the finalized intensity model parquet.
    #[arg(long = "out", required = true)]
    out: PathBuf,
}

/// Available subcommands.
#[derive(Subcommand, Debug)]
enum Command {
    /// Train a scoring model from spectra and a FASTA database, writing the
    /// result to a Parquet model store.
    ///
    /// Boxed to keep the `Command` enum compact (clippy `large_enum_variant`).
    Train(Box<TrainArgs>),

    /// Train a scoring model directly from externally-labeled, high-confidence
    /// PSMs supplied as flat training parquet(s), bypassing the bootstrap
    /// search. Used for the Phase-3 "own models" path.
    ///
    /// Boxed to keep the `Command` enum compact (clippy `large_enum_variant`):
    /// the largest variant (`Train`) dominates the size otherwise.
    #[command(name = "train-from-msnet")]
    TrainFromMsnet(Box<TrainFromMsnetArgs>),

    /// Merge MSNet intensity aggregation parquets into a finalized intensity
    /// model for the strong-score numerator.
    #[command(name = "train-intensity")]
    TrainIntensity(Box<TrainIntensityArgs>),
}

/// Top-level CLI.  When no subcommand is given, the flattened `SearchArgs`
/// drive the existing search path (byte-identical to the pre-subcommand
/// behaviour).
#[derive(Parser, Debug)]
#[command(
    name = "andes",
    about = "andes: database search of MGF/mzML spectra against FASTA",
    allow_hyphen_values = true,
)]
struct TopCli {
    #[command(subcommand)]
    command: Option<Command>,

    #[command(flatten)]
    search: SearchArgs,
}

// Backward-compatibility alias: keep old name pointing to the search args.
type Cli = SearchArgs;

fn main() -> ExitCode {
    #[cfg(feature = "thermo")]
    configure_bundled_dotnet();
    let top = TopCli::parse();
    let result = match top.command {
        Some(Command::Train(args)) => run_train(*args),
        Some(Command::TrainFromMsnet(args)) => run_train_from_msnet(*args),
        Some(Command::TrainIntensity(args)) => run_train_intensity(*args),
        None => {
            // Validate required search args that are Option<> at the clap level.
            let search = top.search;
            if search.spectrum.is_empty() {
                eprintln!("error: --spectrum is required for search (or use `andes train`)");
                return ExitCode::from(2);
            }
            let database = match search.database {
                Some(p) => p,
                None => {
                    eprintln!("error: --database is required for search");
                    return ExitCode::from(2);
                }
            };
            let output_pin = match search.output_pin {
                Some(p) => p,
                None => {
                    eprintln!("error: --output-pin is required for search");
                    return ExitCode::from(2);
                }
            };
            run(Cli {
                database: Some(database),
                output_pin: Some(output_pin),
                ..search
            })
        }
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("andes: {e}");
            ExitCode::from(1)
        }
    }
}

/// Make Thermo `.raw` reading work with zero setup when the runtime is bundled.
///
/// If a .NET runtime ships next to the executable (`<exe_dir>/dotnet`, as the
/// release archives do), point `DOTNET_ROOT` at it so opening a `.raw` "just
/// works". An existing `DOTNET_ROOT` or a system-wide .NET install is left
/// untouched (it takes precedence). No effect on mzML/MGF, which never load .NET.
#[cfg(feature = "thermo")]
fn configure_bundled_dotnet() {
    if std::env::var_os("DOTNET_ROOT").is_some() {
        return;
    }
    if let Ok(exe) = std::env::current_exe() {
        if let Some(dir) = exe.parent() {
            let bundled = dir.join("dotnet");
            if bundled.join("shared").join("Microsoft.NETCore.App").is_dir() {
                std::env::set_var("DOTNET_ROOT", &bundled);
            }
        }
    }
}

/// Print VmRSS for the current process when `ANDES_RSS_PROBE=1`. No-op
/// otherwise and a no-op on non-Linux platforms regardless of the env var.
///
/// We gate behind an env var so production runs stay quiet; flip the var on
/// when debugging memory regressions.
fn log_rss(tag: &str) {
    let probe_set = std::env::var_os("ANDES_RSS_PROBE").is_some()
        || std::env::var_os("MSGF_RSS_PROBE").is_some()
        || std::env::var_os("MSGFRUST_RSS_PROBE").is_some();
    if !probe_set {
        return;
    }
    #[cfg(target_os = "linux")]
    {
        if let Ok(s) = std::fs::read_to_string("/proc/self/status") {
            for line in s.lines() {
                if line.starts_with("VmRSS:") {
                    eprintln!(
                        "[RSS {tag}] {}",
                        line.trim_start_matches("VmRSS:").trim()
                    );
                    return;
                }
            }
        }
    }
    #[cfg(not(target_os = "linux"))]
    {
        let _ = tag;
    }
}

/// Statistics returned by the parser-thread helper.
#[derive(Debug, Default)]
struct ParseStats {
    error_count: usize,
    first_errors: Vec<String>,
}

fn input_format_flags(path: &Path) -> (bool, bool, bool, bool) {
    let ext = path
        .extension()
        .and_then(|e| e.to_str())
        .map(|s| s.to_lowercase());
    let is_mzml = matches!(ext.as_deref(), Some("mzml"));
    let is_raw = matches!(ext.as_deref(), Some("raw"));
    let is_d = matches!(ext.as_deref(), Some("d"));
    let is_mgf = !is_mzml && !is_raw && !is_d;
    (is_mzml, is_raw, is_d, is_mgf)
}

/// Prefix spectrum titles so pooled multi-file PIN SpecIds stay unique.
/// Decide the per-file SpecId/title prefix. Returns `None` for a single-file
/// search so its PIN output stays byte-identical to the pre-multi-file path;
/// returns `Some("<stem>/")` only when disambiguating across multiple inputs.
fn title_prefix_for(num_files: usize, file_stem: &str) -> Option<String> {
    (num_files > 1).then(|| format!("{file_stem}/"))
}

fn prefix_spectrum_titles(chunk: &mut [Spectrum], prefix: &str) {
    for spec in chunk.iter_mut() {
        if spec.title.is_empty() {
            spec.title = format!("{prefix}scan={}", spec.scan.unwrap_or(0));
        } else {
            spec.title = format!("{prefix}{}", spec.title);
        }
    }
}

fn merge_parse_stats(acc: &mut ParseStats, part: ParseStats) {
    acc.error_count += part.error_count;
    for e in part.first_errors {
        if acc.first_errors.len() < 10 {
            acc.first_errors.push(e);
        }
    }
}

/// Producer helper: drains `reader` into fixed-size chunks of `Spectrum`
/// and sends them through `tx`. Stops at `bench_cap` total spectra (or
/// `usize::MAX` for unbounded). Parse errors are counted and the first few
/// captured for downstream reporting; the channel is closed when the
/// reader is exhausted or the consumer hangs up.
///
/// Generic over the reader's error type so the same helper serves both
/// MGF and mzML.
///
/// Runs on a dedicated thread so chunk N+1 is PARSED while chunk N is SCORED.
/// Channel capacity is 2 (one in-flight + one queued) so the producer stays at
/// most one chunk ahead.
fn send_chunks<R, E>(
    reader: R,
    chunk_size: usize,
    bench_cap: usize,
    tx: SyncSender<Vec<Spectrum>>,
) -> ParseStats
where
    R: Iterator<Item = Result<Spectrum, E>>,
    E: std::fmt::Display,
{
    let mut stats = ParseStats::default();
    let mut chunk: Vec<Spectrum> = Vec::with_capacity(chunk_size);
    let mut total = 0usize;
    for result in reader {
        if total >= bench_cap {
            break;
        }
        match result {
            Ok(s) => {
                chunk.push(s);
                total += 1;
                if chunk.len() >= chunk_size {
                    // If the consumer hung up, stop. Sender is moved into the
                    // function, so dropping returns `Err(SendError(chunk))`.
                    let payload = std::mem::replace(&mut chunk, Vec::with_capacity(chunk_size));
                    if tx.send(payload).is_err() {
                        return stats;
                    }
                }
            }
            Err(e) => {
                stats.error_count += 1;
                if stats.first_errors.len() < 3 {
                    stats.first_errors.push(format!("{e}"));
                }
            }
        }
    }
    if !chunk.is_empty() {
        let _ = tx.send(chunk);
    }
    stats
}

/// Lightweight metadata collected in one linear file scan for precursorCal.
#[derive(Debug, Clone)]
struct SpectrumMeta {
    precursor_mz: f64,
    precursor_charge: Option<i32>,
    num_peaks: usize,
}

fn scan_spectrum_metadata(
    path: &Path,
    is_mzml: bool,
    ms_level: u32,
    bench_cap: usize,
) -> Result<Vec<SpectrumMeta>, Box<dyn std::error::Error>> {
    let mut out = Vec::new();
    if is_mzml {
        let f = File::open(path)?;
        let reader = MzMLReader::new(BufReader::new(f)).with_ms_level_range(ms_level, ms_level);
        for result in reader {
            if out.len() >= bench_cap {
                break;
            }
            let spec = result.map_err(|e| format!("mzML parse: {e}"))?;
            out.push(SpectrumMeta {
                precursor_mz: spec.precursor_mz,
                precursor_charge: spec.precursor_charge,
                num_peaks: spec.peaks.len(),
            });
        }
    } else {
        let f = File::open(path)?;
        let reader = MgfReader::new(BufReader::new(f));
        for result in reader {
            if out.len() >= bench_cap {
                break;
            }
            let spec = result.map_err(|e| format!("MGF parse: {e}"))?;
            out.push(SpectrumMeta {
                precursor_mz: spec.precursor_mz,
                precursor_charge: spec.precursor_charge,
                num_peaks: spec.peaks.len(),
            });
        }
    }
    Ok(out)
}

fn build_spec_keys_from_metadata(
    meta: &[SpectrumMeta],
    charge_range: std::ops::RangeInclusive<u8>,
    min_peaks: u32,
) -> Vec<SpecKey> {
    let spectra: Vec<Spectrum> = meta
        .iter()
        .map(|m| Spectrum {
            title: String::new(),
            precursor_mz: m.precursor_mz,
            precursor_intensity: None,
            precursor_charge: m.precursor_charge,
            rt_seconds: None,
            scan: None,
            peaks: vec![(0.0, 0.0); m.num_peaks],
            activation_method: None,
            isolation_lower_offset: None,
            isolation_upper_offset: None,
        })
        .collect();
    build_spec_keys(&spectra, &charge_range, min_peaks)
}

fn load_spectra_by_index(
    path: &Path,
    is_mzml: bool,
    ms_level: u32,
    indices: &HashSet<usize>,
    bench_cap: usize,
) -> Result<HashMap<usize, Spectrum>, Box<dyn std::error::Error>> {
    let mut loaded = HashMap::new();
    if indices.is_empty() {
        return Ok(loaded);
    }
    if is_mzml {
        let f = File::open(path)?;
        let reader = MzMLReader::new(BufReader::new(f)).with_ms_level_range(ms_level, ms_level);
        for (idx, result) in reader.enumerate() {
            if idx >= bench_cap {
                break;
            }
            if !indices.contains(&idx) {
                continue;
            }
            let spec = result.map_err(|e| format!("mzML parse: {e}"))?;
            loaded.insert(idx, spec);
            if loaded.len() == indices.len() {
                break;
            }
        }
    } else {
        let f = File::open(path)?;
        let reader = MgfReader::new(BufReader::new(f));
        for (idx, result) in reader.enumerate() {
            if idx >= bench_cap {
                break;
            }
            if !indices.contains(&idx) {
                continue;
            }
            let spec = result.map_err(|e| format!("MGF parse: {e}"))?;
            loaded.insert(idx, spec);
            if loaded.len() == indices.len() {
                break;
            }
        }
    }
    Ok(loaded)
}

/// Auto-detect an isobaric label (TMT/iTRAQ) by sampling the first `SAMPLE_N`
/// MS2 spectra and inspecting their reporter-ion region. Used only when
/// `--protocol auto` is left at its default, to engage the isobaric windowed
/// peak filter with zero config.
///
/// Returns `None` for `.raw`/`.d` (the sampling reader here is mzML/MGF only —
/// the protocol then stays as-is, byte-identical) and for label-free data, so
/// non-isobaric runs are unchanged. The mzML benchmark datasets (Astral, UPS1,
/// TMT) all flow through the mzML branch.
fn detect_isobaric_sampled(
    path: &Path,
    is_mzml: bool,
    is_mgf: bool,
    ms_level: u32,
    high_res: bool,
) -> Option<input::IsobaricLabel> {
    const SAMPLE_N: usize = 1000;
    if !(is_mzml || is_mgf) {
        return None;
    }
    let indices: HashSet<usize> = (0..SAMPLE_N).collect();
    let loaded = load_spectra_by_index(path, is_mzml, ms_level, &indices, usize::MAX).ok()?;
    let sample: Vec<Spectrum> = loaded.into_values().collect();
    input::detect_isobaric(&sample, high_res)
}

fn tolerance_ppm_display(t: Tolerance) -> Option<f64> {
    match t {
        Tolerance::Ppm(v) => Some(v),
        Tolerance::Da(_) => None,
    }
}

fn run_precursor_calibration(
    spectrum_path: &Path,
    is_mzml: bool,
    ms_level: u32,
    bench_cap: usize,
    params: &SearchParams,
    prepared: &PreparedSearch<'_>,
) -> Result<CalibrationStats, Box<dyn std::error::Error>> {
    if params.precursor_cal_mode == PrecursorCalMode::Off {
        return Ok(CalibrationStats::default());
    }

    let t_cal = std::time::Instant::now();
    let meta = scan_spectrum_metadata(spectrum_path, is_mzml, ms_level, bench_cap)?;
    let spec_keys = build_spec_keys_from_metadata(&meta, params.charge_range.clone(), params.min_peaks);

    if spec_keys.len() < cal_constants::MIN_SPECKEYS_FOR_PREPASS {
        eprintln!(
            "Precursor mass calibration skipped ({} SpecKeys < {} threshold; elapsed: {:.2}s)",
            spec_keys.len(),
            cal_constants::MIN_SPECKEYS_FOR_PREPASS,
            t_cal.elapsed().as_secs_f64()
        );
        return Ok(CalibrationStats::default());
    }

    let sampled = sample_every_nth(
        &spec_keys,
        cal_constants::SAMPLING_STRIDE,
        cal_constants::MAX_SAMPLED,
    );
    let needed: HashSet<usize> = sampled.iter().map(|k| k.spectrum_idx).collect();
    let originals = load_spectra_by_index(spectrum_path, is_mzml, ms_level, &needed, bench_cap)?;

    let stats = learn_calibration_stats(&spec_keys, &originals, prepared, params);

    if stats.has_reliable_stats() {
        eprintln!(
            "Precursor mass shift learned: {:.3} ppm from {} confident PSMs (robust sigma {:.3} ppm; elapsed: {:.2}s)",
            stats.shift_ppm,
            stats.confident_psm_count,
            stats.robust_sigma_ppm,
            t_cal.elapsed().as_secs_f64()
        );
    } else {
        eprintln!(
            "Precursor mass calibration skipped (insufficient confident PSMs: {} with PSMs, {} below RawScore floor, {} failed |residual|>50ppm; elapsed: {:.2}s)",
            stats.queues_with_psm,
            stats.rejected_low_score,
            stats.rejected_residual,
            t_cal.elapsed().as_secs_f64()
        );
    }
    Ok(stats)
}

fn run(cli: Cli) -> Result<(), Box<dyn std::error::Error>> {
    // These three were validated as Some(..) by main() before calling run().
    if cli.spectrum.is_empty() {
        return Err("no --spectrum inputs".into());
    }
    // Parse the digestion enzyme once — drives BOTH model selection
    // (build_selection_key) and digestion (params.enzyme). Previously trypsin
    // was hardcoded in both places.
    let search_enzyme = model::enzyme::Enzyme::from_name(&cli.enzyme)
        .ok_or_else(|| format!(
            "unknown --enzyme '{}' (expected trypsin/chymotrypsin/lysc/aspn/gluc/lysn/argc/alphalp/nocleavage/nonspecific)",
            cli.enzyme
        ))?;
    let spectrum_paths = &cli.spectrum;
    let spectrum_path: PathBuf = spectrum_paths[0].clone();
    let database_path: PathBuf = cli.database.expect("database validated in main");
    let output_pin_path: PathBuf = cli.output_pin.expect("output_pin validated in main");
    if spectrum_paths.len() > 1 {
        eprintln!(
            "Multi-spectrum search: {} inputs → one PIN",
            spectrum_paths.len()
        );
    }

    log_rss("startup");
    let t_total = std::time::Instant::now();
    let t_phase = std::time::Instant::now();
    // ── 1. Load FASTA target database ────────────────────────────────────────
    let target_db =
        FastaReader::load_all(BufReader::new(File::open(&database_path)?))?;
    eprintln!(
        "Loaded {} target proteins from {} [PHASE fasta_load: {:.2}s]",
        target_db.proteins.len(),
        database_path.display(),
        t_phase.elapsed().as_secs_f64()
    );
    log_rss("after_fasta_load");

    // ── 2. Build SearchIndex (target + reversed decoys) ───────────────────────
    let t_phase = std::time::Instant::now();
    let idx = SearchIndex::from_target_db(&target_db, &cli.decoy_prefix);
    eprintln!("[PHASE search_index_build: {:.2}s]", t_phase.elapsed().as_secs_f64());
    log_rss("after_search_index_build");

    // ── 3. Build AminoAcidSet ────────────────────────────────────────────────
    //
    // If --mod is given, parse the mods.txt file. Otherwise
    // fall back to andes's historical defaults (CAM fixed on C,
    // Oxidation variable on M) so existing tests keep their behaviour.
    //
    // `num_mods_from_file` is populated only when --mod is given and the
    // file contains a `NumMods=N` line; it overrides the default
    // `max_variable_mods_per_peptide` (3) below.
    let (mut aa, num_mods_from_file) = match &cli.mods {
        Some(path) => {
            let n = AminoAcidSetBuilder::parse_num_mods_from_file(path)
                .map_err(|e| format!("parsing NumMods= from {}: {e}", path.display()))?;
            let set = AminoAcidSetBuilder::new_standard()
                .add_mods_from_file(path)
                .map_err(|e| format!("loading mods from {}: {e}", path.display()))?
                .build()
                .map_err(|e| format!("building amino-acid set from {}: {e}", path.display()))?;
            eprintln!(
                "Loaded modifications from {} (NumMods={})",
                path.display(),
                n.map(|v| v.to_string()).unwrap_or_else(|| "default".into()),
            );
            (set, n)
        }
        // No --mods: andes defaults (CAM-C fixed, Ox-M variable). The isobaric
        // tag (TMT/iTRAQ) is injected later, after protocol detection (C1).
        None => (default_aa_set_with_tag(None)?, None),
    };

    // ── 4. Load Param scoring model ───────────────────────────────────────────
    //
    // `--param-file` wins outright. Otherwise the model is selected from the
    // Parquet store: for mzML/.raw/.d the activation+analyzer are auto-detected
    // from metadata; for MGF (metadata-less) the `--fragmentation` /
    // `--fragment-tol-*` flags drive `resolve_metadataless_selection`.
    let spectrum_ext = spectrum_path
        .extension()
        .and_then(|e| e.to_str())
        .map(|s| s.to_lowercase());
    let is_mzml = matches!(spectrum_ext.as_deref(), Some("mzml"));
    // Native Thermo `.raw` (feature-gated; needs the .NET 8 runtime).
    let is_raw = matches!(spectrum_ext.as_deref(), Some("raw"));
    // Native Bruker timsTOF `.d` (feature-gated; pure Rust, no vendor runtime).
    // A `.d` is a directory, but the path still carries the `.d` extension.
    let is_d = matches!(spectrum_ext.as_deref(), Some("d"));
    // Anything that is neither mzML, `.raw`, nor `.d` is treated as MGF (default).
    let is_mgf = !is_mzml && !is_raw && !is_d;

    // Detect (activation, instrument) from the input for auto-routing.
    // mzML peeks the file; Thermo `.raw` reads vendor metadata; Bruker `.d`
    // is always CID/TimsTOF (DDA-PASEF). Detection runs for every metadata-
    // bearing format and always wins over the MGF-only `--fragmentation` /
    // `--fragment-tol-*` flags (which carry no metadata of their own).
    let auto_route_eligible = is_mzml || is_raw || is_d;
    let detected_activation_instrument: Option<(ActivationMethod, Option<InstrumentType>)> =
        if !auto_route_eligible {
            None
        } else if is_mzml {
            detect_dominant_activation(&spectrum_path)
                .map(|m| (m, detect_instrument_type_for_path(&spectrum_path)))
        } else if is_raw {
            #[cfg(feature = "thermo")]
            {
                input::thermo::detect_activation_instrument(&spectrum_path, 64)
            }
            #[cfg(not(feature = "thermo"))]
            {
                None
            }
        } else {
            // is_d — timsTOF DDA-PASEF: CID fragmentation on a TOF analyzer.
            Some((ActivationMethod::CID, Some(InstrumentType::TimsTOF)))
        };
    // Pre-compute before the routing match consumes `detected_activation_instrument`.
    let instrument_was_detected = detected_activation_instrument
        .map(|(_, inst)| inst.is_some())
        .unwrap_or(false);

    let t_phase = std::time::Instant::now();
    let mut param = if let Some(ref override_path) = cli.param_file {
        // ── Override path: load binary .param directly (unchanged behaviour). ──
        eprintln!("Param file (override): {}", override_path.display());
        Param::load_from_file(override_path)
            .map_err(|e| format!("loading param file {}: {e}", override_path.display()))?
    } else {
        // ── Resolve (activation, instrument) for the Parquet model store. ─────
        //
        // Metadata-first precedence: a fully detected (activation, instrument)
        // wins outright. When only the activation method is detected (analyzer
        // unknown), or nothing is detected (MGF / metadata-less mzML/.raw), the
        // metadata-less resolver folds in the MGF-only `--fragmentation` and
        // `--fragment-tol-*` flags (decision E default: CID / low-res).
        let (activation, instrument_opt): (ActivationMethod, Option<InstrumentType>) =
            match detected_activation_instrument {
                Some((method, Some(inst))) => {
                    eprintln!(
                        "Param resolver: auto-detected activation = {} (instrument = {}) from {}",
                        method.name(), inst.name(), spectrum_path.display()
                    );
                    (method, Some(inst))
                }
                Some((method, None)) => resolve_metadataless_selection(
                    Some(method), cli.fragmentation, cli.fragment_tol_ppm, cli.fragment_tol_da,
                ),
                None => resolve_metadataless_selection(
                    None, cli.fragmentation, cli.fragment_tol_ppm, cli.fragment_tol_da,
                ),
            };

        let (model_id, p) = load_param_from_store(
            activation,
            instrument_opt,
            cli.protocol,
            search_enzyme,
            cli.model_store.as_deref(),
            cli.model_id_override.as_deref(),
        )?;
        eprintln!("Param model: {model_id} (from store)");
        p
    };
    // Stamp the requested isobaric protocol onto the loaded model so the dense-
    // spectrum windowed peak filter (ScoredSpectrum) engages on TMT/iTRAQ
    // searches even when model selection fell back to a non-isobaric table
    // (there is no bundled CID-TMT model, so `--protocol TMT` resolves to
    // `cid_lowres_tryp`, whose stored protocol is Standard).
    // An explicit `--protocol` wins outright. When left at `auto` (the default),
    // auto-detect TMT/iTRAQ from MS2 reporter ions (mzML/MGF) so the dense-peak
    // windowed filter engages with zero config — the same path `--protocol TMT`
    // takes today. Detection returns None for label-free data, so non-isobaric
    // runs stay byte-identical.
    match cli.protocol {
        Protocol::Tmt => param.data_type.protocol = model::protocol::Protocol::TMT,
        Protocol::Itraq => param.data_type.protocol = model::protocol::Protocol::ITRAQ,
        Protocol::ItraqPhospho => param.data_type.protocol = model::protocol::Protocol::ITRAQPhospho,
        Protocol::Auto => {
            let high_res = param.data_type.instrument.is_high_resolution();
            match detect_isobaric_sampled(&spectrum_path, is_mzml, is_mgf, cli.ms_level as u32, high_res) {
                Some(input::IsobaricLabel::Tmt) => {
                    eprintln!("Protocol resolver: auto-detected TMT reporter ions → engaging isobaric windowed peak filter");
                    param.data_type.protocol = model::protocol::Protocol::TMT;
                }
                Some(input::IsobaricLabel::Itraq) => {
                    eprintln!("Protocol resolver: auto-detected iTRAQ reporter ions → engaging isobaric windowed peak filter");
                    param.data_type.protocol = model::protocol::Protocol::ITRAQ;
                }
                None => {}
            }
        }
        _ => {}
    }
    // C1: parameter-free path only (no explicit --mods). When the protocol
    // resolves to TMT/iTRAQ, inject the tag as a fixed mod on K + peptide
    // N-term so labeled peptides match their precursor mass — otherwise the
    // reporter filter engages but every labeled candidate is +tag Da off and
    // misses. With explicit --mods the user owns the mod set (they may already
    // supply the tag), so those runs stay byte-identical.
    if cli.mods.is_none() {
        let tag = match param.data_type.protocol {
            model::protocol::Protocol::TMT => Some(("TMT6plex", 229.162932_f64)),
            model::protocol::Protocol::ITRAQ | model::protocol::Protocol::ITRAQPhospho => {
                Some(("iTRAQ4plex", 144.102063_f64))
            }
            _ => None,
        };
        if let Some((name, mass)) = tag {
            aa = default_aa_set_with_tag(Some((name, mass)))?;
            eprintln!(
                "Protocol resolver: injected {name} fixed mod (+{mass:.4} on K + peptide N-term) \
                 into the candidate set (no --mods given)"
            );
        }
    }
    let mut scorer = RankScorer::new(&param);
    // Fragment-tol override applies to metadata-less (MGF) input only. For
    // mzML/.raw/.d the analyzer is auto-detected, so the override is ignored.
    let frag_tol_override = cli_fragment_tol_override(cli.fragment_tol_ppm, cli.fragment_tol_da);
    if frag_tol_override.is_some() {
        if instrument_was_detected {
            eprintln!("WARN: --fragment-tol-* ignored — instrument auto-detected from metadata (use --fragment-tol-ppm/-da with MGF input only)");
        } else {
            scorer.set_fragment_tol_override(frag_tol_override);
        }
    }
    eprintln!("[PHASE param_and_scorer: {:.2}s]", t_phase.elapsed().as_secs_f64());

    // ── 5. Build SearchParams ─────────────────────────────────────────────────
    let mut params = SearchParams::default_tryptic(aa);
    params.precursor_tolerance =
        match (cli.precursor_tol_left_ppm, cli.precursor_tol_right_ppm, cli.precursor_tol_da) {
            (Some(l), Some(r), _) => {
                PrecursorTolerance::asymmetric(Tolerance::Ppm(l), Tolerance::Ppm(r))
            }
            (Some(_), None, _) | (None, Some(_), _) => {
                return Err("--precursor-tol-left-ppm and --precursor-tol-right-ppm must be given together".into());
            }
            (None, None, Some(da)) => PrecursorTolerance::symmetric(Tolerance::Da(da)),
            (None, None, None) => PrecursorTolerance::symmetric(Tolerance::Ppm(cli.precursor_tol_ppm)),
        };
    params.charge_range = cli.charge_min..=cli.charge_max;
    if cli.charge_min > cli.charge_max {
        return Err(format!(
            "invalid charge range: --charge-min {} > --charge-max {}",
            cli.charge_min, cli.charge_max
        ).into());
    }
    params.isotope_error_range = cli.isotope_error_min..=cli.isotope_error_max;
    if cli.isotope_error_min > cli.isotope_error_max {
        return Err(format!(
            "invalid isotope error range: --isotope-error-min {} > --isotope-error-max {}",
            cli.isotope_error_min, cli.isotope_error_max
        ).into());
    }
    // Pass 2 co-isolation requires MS1 scans, captured by the mzML and Thermo
    // `.raw` readers. MGF (no MS1) and the Bruker `.d` reader (DDA MS2 only;
    // chimeric on `.d` is out of scope) make `--chimeric` inert, so keep
    // `params.chimeric` FALSE to turn the ENTIRE chimeric path off (Pass 2, PIN
    // column/SpecId gates, top-N forcing) — the run is then identical to a normal search.
    let chimeric_active = cli.chimeric && (is_mzml || is_raw);
    if cli.chimeric && !(is_mzml || is_raw) {
        eprintln!(
            "WARN: --chimeric requires MS1 data (mzML or Thermo .raw); the input is {} \
             so the co-isolation cascade is disabled and the search runs normally.",
            if is_d { "Bruker .d (DDA MS2 only)," } else { "MGF," }
        );
    }
    // The cascade pairs MS2 with its preceding MS1 — it is MS2-only by
    // construction. Ignore a non-2 `--ms-level` under `--chimeric` so MS3+
    // (e.g. TMT SPS-MS3) can never enter the search on any input format.
    if chimeric_active && cli.ms_level != 2 {
        eprintln!(
            "WARN: --ms-level={} is ignored under --chimeric; the cascade always searches MS2.",
            cli.ms_level
        );
    }
    params.chimeric = chimeric_active;
    params.chimeric_isolation_halfwidth_da = cli.isolation_halfwidth;
    // FORCE top-1 under the cascade: Pass 1 emits only the best primary per scan;
    // secondaries come from Pass 2. The default top_n (10) would make Pass 1 emit
    // blind multi-emission per scan = inflated FDR.
    params.top_n_psms_per_spectrum = if chimeric_active { 1 } else { cli.top_n };
    params.enzyme = search_enzyme;
    params.num_tolerable_termini = match cli.enzyme_specificity {
        EnzymeSpecificity::Fully => 2,
        EnzymeSpecificity::Semi => 1,
        EnzymeSpecificity::NonSpecific => 0,
    };
    params.max_missed_cleavages = cli.max_missed_cleavages;
    params.min_peaks = cli.min_peaks;
    params.min_length = cli.min_length;
    params.max_length = cli.max_length;
    params.max_variable_mods_per_peptide = cli.max_mods;
    if let Some(n) = num_mods_from_file {
        params.max_variable_mods_per_peptide = n; // NumMods= in --mods overrides --max-mods
    }
    params.precursor_cal_mode = cli.precursor_cal;
    params.precursor_mass_shift_ppm = 0.0;
    params.score_mode = match cli.score {
        ScoreFlag::Rank => search::ScoreMode::Rank,
        ScoreFlag::Strong => search::ScoreMode::Strong,
    };
    if params.score_mode == search::ScoreMode::Strong {
        eprintln!("score mode: strong (ranking + PIN RawScore use StrongScore)");
    }

    // ── Resolved-parameter banner (reanalysis auditability) ───────────────────
    // One consolidated record of every resolved search parameter, so a
    // (zero-config) run is fully reproducible/auditable from its log. Values that
    // were auto-detected from the data/store are tagged [detected].
    eprintln!("──────── andes resolved parameters ────────");
    eprintln!("  spectra        : {}", spectrum_path.display());
    eprintln!("  model          : (see 'Param model:' line above) [detected]");
    eprintln!("  activation     : {:?} [detected]", param.data_type.activation);
    eprintln!("  instrument     : {:?} [detected]", param.data_type.instrument);
    eprintln!("  protocol       : {:?}", param.data_type.protocol);
    eprintln!("  enzyme         : {} ({:?} termini, <={} missed cleavages)",
              search_enzyme.name(), cli.enzyme_specificity, params.max_missed_cleavages);
    eprintln!("  mods           : {}",
              if cli.mods.is_some() { "from --mods file" }
              else { "defaults (Cam-C fixed, Ox-M variable) + isobaric tag if detected" });
    eprintln!("  max var-mods   : {} per peptide", params.max_variable_mods_per_peptide);
    eprintln!("  peptide length : {}-{}", params.min_length, params.max_length);
    eprintln!("  precursor tol  : {:?} (calibration: {:?})", params.precursor_tolerance, params.precursor_cal_mode);
    eprintln!("  charge range   : {}-{}", params.charge_range.start(), params.charge_range.end());
    eprintln!("  isotope errors : {}..={}", params.isotope_error_range.start(), params.isotope_error_range.end());
    eprintln!("  decoy prefix   : {}   chimeric: {}", cli.decoy_prefix, params.chimeric);
    eprintln!("───────────────────────────────────────────");

    // ── 6+7. Stream-load + chunked search ─────────────────────────────────
    //
    // Spectra are parsed and scored in chunks of CHUNK_SIZE. Each chunk's
    // peak data lives in RAM only for the time it takes to score the chunk,
    // then is dropped before the next chunk is read. The Vec<Spectrum> that
    // survives into the PIN/TSV writers retains scan/title/precursor_mz/scan
    // (the only fields the writers read) but has empty peaks.
    //
    // This bounds peak-data memory to ~CHUNK_SIZE × per-spectrum peak size
    // regardless of dataset size — fixes the Astral-scale OOM where loading
    // all 123k spectra at once pushed RSS to 28 GB on a 31 GB VM.
    const CHUNK_SIZE: usize = 5000;

    let t_phase = std::time::Instant::now();

    // Configure the global Rayon worker pool BEFORE we build PreparedSearch
    // or run any chunks. `build_global()` panics if called twice; guard with
    // `OnceLock` so repeated CLI invocations within a single test process
    // don't blow up.
    static POOL_INIT: std::sync::OnceLock<()> = std::sync::OnceLock::new();
    POOL_INIT.get_or_init(|| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(cli.threads)
            .build_global()
            .expect("build_global");
    });
    eprintln!("Using {} worker threads", cli.threads);

    // Fragment tolerance of 0.5 Da matches the gf_bsa_parity integration test
    // (and the canonical HCD default).
    let fragment_tol_da = 0.5_f64;

    let bench_cap = if cli.max_spectra > 0 {
        cli.max_spectra
    } else {
        usize::MAX
    };
    let ms_level_u32 = cli.ms_level as u32;
    // Native `.raw`/`.d` readers search MS2 (identification) scans only. A non-2
    // `--ms-level` would otherwise make the Thermo iterator emit MS3 scans (which
    // carry a precursor and would be searched — e.g. TMT SPS-MS3 reporter scans)
    // or MS1 (no precursor → an empty run). Force MS2 + warn. `--ms-level` still
    // applies to mzML.
    if (is_raw || is_d) && cli.ms_level != 2 {
        eprintln!(
            "WARN: --ms-level={} is ignored for native .raw/.d input; these formats \
             search MS2 (identification) scans only.",
            cli.ms_level
        );
    }

    // The precursor-calibration pre-pass currently reads only mzML/MGF. For a
    // Thermo `.raw` or Bruker `.d` it would be misread as MGF, so skip
    // calibration and warn (native-format calibration support is a follow-up).
    if (is_raw || is_d) && params.precursor_cal_mode != PrecursorCalMode::Off {
        let fmt = if is_raw { "Thermo .raw" } else { "Bruker .d" };
        eprintln!(
            "WARN: --precursor-cal is not yet supported for {fmt} input; \
             proceeding without calibration."
        );
        params.precursor_cal_mode = PrecursorCalMode::Off;
    }

    // Calibration pre-pass. Candidate enumeration is precursor-tolerance
    // independent, so keep the cal pass's `PreparedParts` and reuse them for the
    // main pass instead of re-enumerating all 16.8M candidates (~15s saved on
    // Astral). `into_parts()` runs BEFORE tightening so the owned parts outlive
    // the `params` borrow the cal `PreparedSearch` held.
    let reuse_parts = if params.precursor_cal_mode != PrecursorCalMode::Off {
        let cal_prepared = PreparedSearch::prepare(
            &idx,
            &params,
            &scorer,
            fragment_tol_da,
            &cli.decoy_prefix,
        );
        let cal_stats = run_precursor_calibration(
            &spectrum_path,
            is_mzml,
            ms_level_u32,
            bench_cap,
            &params,
            &cal_prepared,
        )?;
        let parts = cal_prepared.into_parts();
        params.precursor_mass_shift_ppm = apply_shift_for_mode(params.precursor_cal_mode, cal_stats);
        let tol_before = params.precursor_tolerance;
        apply_tightened_precursor_tolerance(&mut params, cal_stats);
        if cal_stats.has_reliable_stats() {
            let left_before = tolerance_ppm_display(tol_before.left);
            let right_before = tolerance_ppm_display(tol_before.right);
            let left_after = tolerance_ppm_display(params.precursor_tolerance.left);
            let right_after = tolerance_ppm_display(params.precursor_tolerance.right);
            if left_after.is_some()
                && right_after.is_some()
                && (left_after != left_before || right_after != right_before)
            {
                eprintln!(
                    "Tightened precursor tolerance for main pass: left {:.3} ppm -> {:.3} ppm, right {:.3} ppm -> {:.3} ppm",
                    left_before.unwrap_or(0.0),
                    left_after.unwrap_or(0.0),
                    right_before.unwrap_or(0.0),
                    right_after.unwrap_or(0.0),
                );
            }
        }
        Some(parts)
    } else {
        None
    };

    let intensity_model: Option<Arc<scoring_crate::IntensityModel>> = cli
        .intensity_model
        .as_ref()
        .map(|path| {
            eprintln!("loading intensity model from {} ...", path.display());
            scoring_crate::IntensityModel::load(path)
                .map(Arc::new)
                .map_err(|e| format!("intensity model {}: {e}", path.display()))
        })
        .transpose()?;

    let prepared = match reuse_parts {
        Some(parts) => {
            PreparedSearch::from_parts(&idx, &params, &scorer, fragment_tol_da, parts)
        }
        None => PreparedSearch::prepare(&idx, &params, &scorer, fragment_tol_da, &cli.decoy_prefix),
    }
    .with_intensity_model(intensity_model);
    log_rss("after_prepared_search");
    eprintln!(
        "PreparedSearch: {} candidates, {} mass buckets",
        prepared.candidates.len(),
        prepared.bucket_index.len(),
    );

    let bench_mode = cli.max_spectra > 0;

    let mut all_spectra: Vec<Spectrum> = Vec::new();
    let mut all_queues: Vec<TopNQueue> = Vec::new();

    let t_search_start = std::time::Instant::now();

    // Pipeline mzML/MGF parsing with Rayon scoring via a bounded sync_channel.
    // The parser runs on a dedicated thread and pushes
    // CHUNK_SIZE-sized `Vec<Spectrum>` payloads through the channel; the main
    // thread (this one) drains the channel and calls `prepared.run_chunk` on
    // each chunk (which is itself Rayon-parallel internally). With capacity 2
    // the parser stays at most one chunk ahead of the scorer, overlapping
    // parse-of-chunk-(N+1) with score-of-chunk-N — so parse time (Astral
    // ~2-3s per chunk) overlaps scoring instead of running serially.
    // MGF carries no MS level (always treated as MS2). (Native `.raw`/`.d` are
    // warned separately above.)
    let mzml_warn_ms_level_emitted = if is_mgf && cli.ms_level != 2 {
        eprintln!(
            "WARN: --ms-level={} requested for an MGF input; MGF carries no MS \
             level (always treated as MS2). The flag has no effect on this input.",
            cli.ms_level
        );
        true
    } else {
        false
    };
    let _ = mzml_warn_ms_level_emitted; // silenced — unused for now.

    // Under `--chimeric` (mzML or Thermo `.raw`), stream MS2 in CHUNK_SIZE batches,
    // each paired with a bounded per-chunk MS1 link (`read_with_ms1_chunked`). Pass 1 + Pass 2
    // run per chunk on this thread and peaks are dropped immediately, so RSS stays
    // bounded to ~CHUNK_SIZE spectra (NOT the whole file). The parser runs on a
    // dedicated thread so chunk N+1 parses while chunk N scores. The streaming
    // pipeline in the `else` branch handles the (default) non-chimeric path.
    let chimeric_input = chimeric_active;
    let mut parse_stats = ParseStats::default();

    for (file_idx, input_path) in spectrum_paths.iter().enumerate() {
        if bench_mode && all_spectra.len() >= bench_cap {
            break;
        }
        let remaining_cap = if bench_mode {
            bench_cap.saturating_sub(all_spectra.len())
        } else {
            usize::MAX
        };

        let (file_is_mzml, file_is_raw, file_is_d, _file_is_mgf) =
            input_format_flags(input_path);
        let file_stem = input_path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("spectrum");
        let title_prefix = title_prefix_for(spectrum_paths.len(), file_stem);

        if spectrum_paths.len() > 1 {
            eprintln!(
                "=== spectrum [{}/{}] {} ===",
                file_idx + 1,
                spectrum_paths.len(),
                input_path.display()
            );
        }

        if chimeric_input && !(file_is_mzml || file_is_raw) {
            return Err(format!(
                "--chimeric only supports mzML/.raw inputs, got {}",
                input_path.display()
            )
            .into());
        }

        let file_stats = if chimeric_input {
            let (tx, rx) = sync_channel::<(Vec<Spectrum>, Ms1Link)>(2);
            let spectrum_path = input_path.clone();
            let cap = remaining_cap;
            // The cascade is MS2-only by construction (MS2 paired with its preceding
            // MS1); hardcode MS2 so `--ms-level 3` can never widen the mzML reader's
            // range to admit MS3 (the .raw chunked reader is already MS2-only).
            let mslevel = 2;
            let parser_handle = thread::spawn(move || -> Result<(usize, Vec<String>), String> {
                if file_is_mzml {
                    let f = File::open(&spectrum_path).map_err(|e| format!("open mzML: {e}"))?;
                    let reader = MzMLReader::new(BufReader::new(f))
                        .with_ms_level_range(mslevel, mslevel)
                        .with_ms1_capture(true);
                    let (errc, errs) =
                        reader.read_with_ms1_chunked(CHUNK_SIZE, cap, |chunk, link| {
                            let _ = tx.send((chunk, link));
                        });
                    Ok((errc, errs))
                } else {
                    #[cfg(feature = "thermo")]
                    {
                        let reader = input::ThermoRawReader::open(&spectrum_path)
                            .map_err(|e| format!("open Thermo .raw: {e}"))?;
                        let (errc, errs) =
                            reader.read_with_ms1_chunked(CHUNK_SIZE, cap, |chunk, link| {
                                let _ = tx.send((chunk, link));
                            });
                        Ok((errc, errs))
                    }
                    #[cfg(not(feature = "thermo"))]
                    {
                        Err("this andes build has no Thermo .raw support; \
                             rebuild with `--features thermo`."
                            .to_string())
                    }
                }
            });

            let mut file_offset = 0usize;
            let mut ms1_linked = 0usize;
            for (mut chunk_spectra, chunk_link) in rx {
                if let Some(prefix) = &title_prefix {
                    prefix_spectrum_titles(&mut chunk_spectra, prefix);
                }
                let offset = all_spectra.len();
                let mut queues = prepared.run_chunk(&chunk_spectra, offset);
                search::match_engine::run_pass2_coisolation(
                    &prepared,
                    &chunk_spectra,
                    &mut queues,
                    &params,
                    &chunk_link,
                    offset,
                );
                file_offset += chunk_spectra.len();
                ms1_linked += chunk_link.ms1_peaks.len();
                all_queues.extend(queues);
                for mut spec in chunk_spectra.into_iter() {
                    spec.peaks = Vec::new();
                    all_spectra.push(spec);
                }
            }
            let (err_count, first_errors) = parser_handle
                .join()
                .map_err(|_| "chimeric parser thread panicked".to_string())??;
            eprintln!(
                "chimeric mode: streamed {} MS2 spectra ({} MS1 scans linked) from {}",
                file_offset, ms1_linked, input_path.display()
            );
            log_rss("after_chimeric_stream_search");
            ParseStats {
                error_count: err_count,
                first_errors,
            }
        } else {
            let (tx, rx) = sync_channel::<Vec<Spectrum>>(2);
            let spectrum_path = input_path.clone();
            let parser_handle = thread::spawn(
                move || -> Result<ParseStats, Box<dyn std::error::Error + Send + Sync>> {
                    if file_is_mzml {
                        let f = File::open(&spectrum_path)
                            .map_err(|e| format!("open mzML: {e}"))?;
                        let reader = MzMLReader::new(BufReader::new(f))
                            .with_ms_level_range(ms_level_u32, ms_level_u32);
                        Ok(send_chunks(reader, CHUNK_SIZE, remaining_cap, tx))
                    } else if file_is_raw {
                        #[cfg(feature = "thermo")]
                        {
                            let reader = input::ThermoRawReader::open(&spectrum_path)
                                .map_err(|e| format!("open Thermo .raw: {e}"))?
                                .with_ms_level(Some(2));
                            Ok(send_chunks(reader, CHUNK_SIZE, remaining_cap, tx))
                        }
                        #[cfg(not(feature = "thermo"))]
                        {
                            Err("this andes build has no Thermo .raw support; \
                                 rebuild with `--features thermo` (and run with the \
                                 .NET 8 runtime installed). mzML/MGF inputs work without it."
                                .into())
                        }
                    } else if file_is_d {
                        #[cfg(feature = "timstof")]
                        {
                            let reader = input::TimsTofReader::open(&spectrum_path)
                                .map_err(|e| format!("open Bruker .d: {e}"))?;
                            Ok(send_chunks(reader, CHUNK_SIZE, remaining_cap, tx))
                        }
                        #[cfg(not(feature = "timstof"))]
                        {
                            Err("this andes build has no Bruker .d (timsTOF) support; \
                                 rebuild with `--features timstof`. mzML/MGF inputs work \
                                 without it."
                                .into())
                        }
                    } else {
                        let f = File::open(&spectrum_path)
                            .map_err(|e| format!("open MGF: {e}"))?;
                        let reader = MgfReader::new(BufReader::new(f));
                        Ok(send_chunks(reader, CHUNK_SIZE, remaining_cap, tx))
                    }
                },
            );

            log_rss("after_parser_thread_spawn");

            for mut chunk in rx {
                if chunk.is_empty() {
                    continue;
                }
                if let Some(prefix) = &title_prefix {
                    prefix_spectrum_titles(&mut chunk, prefix);
                }
                let offset = all_spectra.len();
                let queues = prepared.run_chunk(&chunk, offset);
                all_queues.extend(queues);
                for mut spec in chunk.into_iter() {
                    spec.peaks = Vec::new();
                    all_spectra.push(spec);
                }
                log_rss(&format!("after_chunk_{:06}_specs", all_spectra.len()));
            }

            match parser_handle.join() {
                Ok(Ok(stats)) => stats,
                Ok(Err(e)) => return Err(format!("parser thread error: {e}").into()),
                Err(_) => return Err("parser thread panicked".into()),
            }
        };
        merge_parse_stats(&mut parse_stats, file_stats);
    }

    if parse_stats.error_count > 0 {
        eprintln!(
            "WARN: {} spectra failed to parse{}",
            parse_stats.error_count,
            if !parse_stats.first_errors.is_empty() {
                format!(" (first {}):", parse_stats.first_errors.len())
            } else {
                String::new()
            }
        );
        for e in &parse_stats.first_errors {
            eprintln!("  - {e}");
        }
    }

    if is_mzml || is_raw {
        eprintln!(
            "MS-level filter: {} (only MS{} spectra entered the search)",
            cli.ms_level, cli.ms_level
        );
    }

    if all_spectra.is_empty() {
        let paths: Vec<String> = spectrum_paths
            .iter()
            .map(|p| p.display().to_string())
            .collect();
        return Err(format!("no spectra parsed from {}", paths.join(", ")).into());
    }

    log_rss("after_all_spectra");
    let search_elapsed = t_search_start.elapsed();
    eprintln!(
        "Loaded+scored {} spectra from {} in chunks of {} [PHASE stream_search: {:.2}s]",
        all_spectra.len(),
        spectrum_path.display(),
        CHUNK_SIZE,
        t_phase.elapsed().as_secs_f64()
    );
    if bench_mode {
        eprintln!("Bench mode: capped at {} spectra", cli.max_spectra);
    }

    // Downstream code uses these names.
    let spectra = all_spectra;
    let queues = all_queues;

    let non_empty = queues.iter().filter(|q| !q.is_empty()).count();
    eprintln!(
        "Search complete: {non_empty} / {} spectra have PSMs (match_spectra wall: {:.2}s)",
        spectra.len(),
        search_elapsed.as_secs_f64()
    );

    // ── 8. Write PIN ─────────────────────────────────────────────────────────
    // Bench mode still writes PIN (so we can diff against the reference
    // fixture) but skips TSV.
    let t_phase = std::time::Instant::now();
    output::write_pin(&output_pin_path, &spectra, &queues, &prepared.candidates, &params, &idx)?;
    eprintln!(
        "Wrote PIN: {} [PHASE pin_write: {:.2}s] [PHASE TOTAL: {:.2}s]",
        output_pin_path.display(),
        t_phase.elapsed().as_secs_f64(),
        t_total.elapsed().as_secs_f64()
    );
    log_rss("after_pin_write");

    if bench_mode {
        eprintln!("Bench mode: skipping TSV write.");
        return Ok(());
    }

    // ── 9. Write TSV (optional) ───────────────────────────────────────────────
    if let Some(ref tsv_path) = cli.output_tsv {
        let spec_file_name = if spectrum_paths.len() == 1 {
            spectrum_path
                .file_name()
                .map(|n| n.to_string_lossy().into_owned())
                .unwrap_or_else(|| spectrum_path.display().to_string())
        } else {
            spectrum_paths
                .iter()
                .filter_map(|p| p.file_name().map(|n| n.to_string_lossy().into_owned()))
                .collect::<Vec<_>>()
                .join("+")
        };
        output::write_tsv(tsv_path, &spectra, &queues, &prepared.candidates, &params, &idx, &spec_file_name, is_mgf)?;
        eprintln!("Wrote TSV: {}", tsv_path.display());
    }

    Ok(())
}

// ── Training pipeline ─────────────────────────────────────────────────────────

/// Load all MS2 spectra from a path using the same format-dispatch logic as
/// the search path (mzML by extension, otherwise MGF).
fn load_spectra_for_train(
    path: &Path,
) -> Result<Vec<Spectrum>, Box<dyn std::error::Error>> {
    let ext_lower = path
        .extension()
        .and_then(|e| e.to_str())
        .map(|s| s.to_lowercase());
    let mut spectra = Vec::new();
    match ext_lower.as_deref() {
        Some("mzml") => {
            let f = File::open(path)?;
            let reader = MzMLReader::new(BufReader::new(f)).with_ms_level_range(2, 2);
            for item in reader {
                match item {
                    Ok(s) => spectra.push(s),
                    Err(e) => eprintln!("WARN: mzML parse: {e}"),
                }
            }
        }
        // Native Thermo `.raw` — MS2 only, same Spectrum stream as the search
        // path. Requires building with `--features thermo`.
        Some("raw") => {
            #[cfg(feature = "thermo")]
            {
                let reader = input::ThermoRawReader::open(path)
                    .map_err(|e| format!("open Thermo .raw {}: {e}", path.display()))?
                    .with_ms_level(Some(2));
                for item in reader {
                    match item {
                        Ok(s) => spectra.push(s),
                        Err(e) => eprintln!("WARN: .raw parse: {e}"),
                    }
                }
            }
            #[cfg(not(feature = "thermo"))]
            {
                return Err(format!(
                    "native Thermo `.raw` training input requires building with \
                     `--features thermo` (and the .NET 8 runtime): {}",
                    path.display()
                )
                .into());
            }
        }
        // Native Bruker timsTOF `.d` (a directory). Requires `--features timstof`.
        Some("d") => {
            #[cfg(feature = "timstof")]
            {
                let reader = input::TimsTofReader::open(path)
                    .map_err(|e| format!("open Bruker .d {}: {e}", path.display()))?;
                for item in reader {
                    match item {
                        Ok(s) => spectra.push(s),
                        Err(e) => eprintln!("WARN: .d parse: {e}"),
                    }
                }
            }
            #[cfg(not(feature = "timstof"))]
            {
                return Err(format!(
                    "native Bruker `.d` training input requires building with \
                     `--features timstof`: {}",
                    path.display()
                )
                .into());
            }
        }
        // MGF (default / backwards-compatible).
        _ => {
            let f = File::open(path)?;
            let reader = MgfReader::new(BufReader::new(f));
            for item in reader {
                match item {
                    Ok(s) => spectra.push(s),
                    Err(e) => eprintln!("WARN: MGF parse: {e}"),
                }
            }
        }
    }
    Ok(spectra)
}

/// Build the `SearchParams` used by every training mode (initial bootstrap,
/// `--add`, and the acceptance gate) in one place, so they stay consistent with
/// each other and with the production search:
///   - charge span `2..=5` (the search binary's default; `default_tryptic` alone
///     is the narrow `2..=3`, which drops z=4/5 labels common on Astral/timsTOF),
///   - the `NumMods=` variable-mod limit from the `--mods` file when present
///     (the search path applies it too).
fn build_train_search_params(
    mods: &Option<PathBuf>,
) -> Result<SearchParams, Box<dyn std::error::Error>> {
    let aa = build_aa_set(mods)?;
    let mut params = SearchParams::default_tryptic(aa);
    params.charge_range = 2..=5;
    if let Some(path) = mods {
        if let Some(n) = AminoAcidSetBuilder::parse_num_mods_from_file(path)
            .map_err(|e| format!("parsing NumMods= from {}: {e}", path.display()))?
        {
            params.max_variable_mods_per_peptide = n;
        }
    }
    Ok(params)
}

/// Run the full training pipeline and write a model to a Parquet store.
///
/// When `args.update_model` is set, runs in incremental update mode (Part D).
/// Otherwise runs the standard initial-training pipeline (Part A).
fn run_train(args: TrainArgs) -> Result<(), Box<dyn std::error::Error>> {
    let t0 = std::time::Instant::now();

    // ── 1. Configure Rayon thread pool ────────────────────────────────────────
    static POOL_INIT: std::sync::OnceLock<()> = std::sync::OnceLock::new();
    POOL_INIT.get_or_init(|| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .expect("build_global");
    });

    if let Some(ref update_model_id) = args.update_model.clone() {
        return run_train_update(args, update_model_id, t0);
    }

    // ── Standard training path ────────────────────────────────────────────────

    // ── 2. Load spectra ───────────────────────────────────────────────────────
    let spectra_path = args.spectra.clone().ok_or("--spectra is required for initial training")?;
    eprintln!("train: loading spectra from {} ...", spectra_path.display());
    let spectra = load_spectra_for_train(&spectra_path)?;
    eprintln!("train: loaded {} spectra", spectra.len());

    let database = args.database.clone().ok_or("--database is required for initial training")?;

    // ── 3. Load seed Param + RankScorer ──────────────────────────────────────
    let (seed_model_id, seed_param): (String, Param) = load_seed_param(&args.seed_model)?;
    eprintln!("train: seed model = {seed_model_id}");
    let seed_scorer = RankScorer::new(&seed_param);

    // ── 4-5. Build search params (charge span + NumMods) + bootstrap labels ───
    let search_params = build_train_search_params(&args.mods)?;
    eprintln!(
        "train: running seed search (train-fdr = {}) ...",
        args.train_fdr
    );
    let labels = bootstrap_labels(
        &spectra,
        &database,
        &seed_scorer,
        &search_params,
        args.train_fdr,
    )
    .map_err(|e| format!("bootstrap_labels: {e}"))?;
    eprintln!("train: {} confident labels at q <= {}", labels.len(), args.train_fdr);

    if labels.is_empty() {
        return Err(format!(
            "no confident labels found at train-fdr={} — try a higher --train-fdr",
            args.train_fdr
        )
        .into());
    }

    // ── 6. Accumulate stats ───────────────────────────────────────────────────
    eprintln!("train: accumulating ion-match statistics ...");
    let accumulator = StatsAccumulator::new(&seed_scorer);
    let mut stats = CountStats::new();
    for label in &labels {
        let spec = &spectra[label.spectrum_index];
        accumulator.accumulate(&mut stats, spec, &label.peptide, label.charge);
    }
    let stats = merge(vec![stats]);

    // ── 7. Estimate model ─────────────────────────────────────────────────────
    eprintln!("train: estimating model parameters ...");
    let cfg = EstimatorConfig::default();
    let estimator = Estimator::new(cfg);
    let trained_param = estimator.estimate(&stats, &seed_param);
    let n_partitions = trained_param.partitions.len();
    eprintln!("train: trained model has {} partitions", n_partitions);

    // ── 8. Determine model ID ─────────────────────────────────────────────────
    let model_id = args
        .model_id
        .clone()
        .unwrap_or_else(|| format!("trained_{}_{}", args.instrument, args.protocol));
    eprintln!("train: model ID = {model_id}");

    // ── 9. Build source ledger (Part A) ───────────────────────────────────────
    // Determine date: use args.date if provided, else today's date in ISO 8601.
    let date_str = args.date.clone().unwrap_or_else(|| {
        // Format today as YYYY-MM-DD using std::time.
        format_today_iso8601()
    });
    let spectra_filename = spectra_path
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_else(|| spectra_path.display().to_string());
    let source_id = format!("bootstrap_{model_id}");
    let ledger = SourceLedger {
        source_id: source_id.clone(),
        dataset: spectra_filename,
        n_psms: labels.len() as i64,
        date: date_str,
        weight: 1.0,
        train_fdr: args.train_fdr as f32,
        instrument: args.instrument.clone(),
        experiment_class: args.protocol.clone(),
    };

    // ── 10. Write to store with source tracking ────────────────────────────────
    // Read existing OTHER models from the store (preserve them on append).
    let store_path = &args.out_store;
    let mut existing_other: Vec<ModelEntryOwned> = Vec::new();
    if store_path.exists() {
        let store = ModelStore::open(store_path)
            .map_err(|e| format!("opening existing store {}: {e}", store_path.display()))?;
        for id in store.model_ids() {
            if id == model_id {
                eprintln!("train: overwriting existing model '{id}' in store");
                continue;
            }
            let p = store.load_param(&id)
                .map_err(|e| format!("reading model '{id}': {e}"))?;
            let src_ledgers = store.load_sources(&id)
                .unwrap_or_default();
            let mut src = Vec::new();
            for l in src_ledgers {
                if let Ok(s) = store.load_source_stats(&id, &l.source_id) {
                    src.push((l, s));
                }
            }
            existing_other.push((id, p, src));
        }
    }

    // Combine all models: write the trained model + existing others together.
    let mut all_entries: Vec<ModelEntryOwned> = Vec::new();
    all_entries.push((model_id.clone(), trained_param.clone(), vec![(ledger, stats)]));
    for (id, p, src) in existing_other {
        all_entries.push((id, p, src));
    }

    write_all_models_with_sources_pub(
        store_path,
        &all_entries.iter()
            .map(|(id, p, s)| (id.as_str(), p, s.as_slice()))
            .collect::<Vec<_>>(),
    )
    .map_err(|e| format!("writing model store {}: {e}", store_path.display()))?;

    eprintln!(
        "train: wrote model '{model_id}' to {} (source '{source_id}') [{:.2}s]",
        store_path.display(),
        t0.elapsed().as_secs_f64(),
    );
    eprintln!(
        "train: summary — labels={}, partitions={}, store={}",
        labels.len(),
        n_partitions,
        store_path.display(),
    );

    Ok(())
}

// ════════════════════════════════════════════════════════════════════════════
// train-from-msnet: supervised training from externally-labeled PSM parquets
// ════════════════════════════════════════════════════════════════════════════

/// One confident, externally-labeled PSM read from a flat training parquet.
///
/// This is the in-memory form of one parquet row (see the column contract in
/// the `train-from-msnet` CLI docs). The peaks are stored already sorted
/// ascending by m/z (the parquet stores acquisition order; the reader sorts).
struct MsnetPsm {
    spectrum: Spectrum,
    peptide: model::peptide::Peptide,
    charge: u8,
}

/// Build a [`model::peptide::Peptide`] from a bare uppercase sequence plus
/// resolved modification mass-deltas.
///
/// Modifications are applied two ways, matching how the scoring code computes
/// peptide mass (`Peptide::new` sums `aa.mass + mod_.mass_delta` per residue):
///
/// - **Residue mods** (`res_mod_pos` / `res_mod_delta`, 1-based positions):
///   attached to the residue at that position via `AminoAcid::with_mod` with a
///   `ModLocation::Anywhere`, `ResidueSpec::Specific(residue)` modification.
/// - **Terminal mods** (`nterm_delta` / `cterm_delta`): folded onto the first /
///   last residue's mass-delta. Because `Peptide::new` only sums per-residue
///   `mod_.mass_delta`, a terminal delta must be carried by a residue's `mod_`
///   to be counted. If a residue already carries a residue mod, the terminal
///   delta is *added* to that residue's existing delta (a single combined
///   `Modification`); otherwise a fresh terminal `Modification` is attached.
///   This keeps `peptide.mass()` correct regardless of overlap.
///
/// Returns an error if `seq` contains a non-standard residue or the mod arrays
/// are misaligned.
fn build_msnet_peptide(
    seq: &str,
    res_mod_pos: &[i32],
    res_mod_delta: &[f64],
    nterm_delta: f64,
    cterm_delta: f64,
) -> Result<model::peptide::Peptide, Box<dyn std::error::Error>> {
    use std::sync::Arc;

    if res_mod_pos.len() != res_mod_delta.len() {
        return Err(format!(
            "res_mod_pos ({}) and res_mod_delta ({}) length mismatch",
            res_mod_pos.len(),
            res_mod_delta.len()
        )
        .into());
    }
    let bytes = seq.as_bytes();
    if bytes.is_empty() {
        return Err("empty peptide sequence".into());
    }
    let n = bytes.len();

    // Accumulate the total mod delta to apply to each residue (1-based -> 0-based).
    // Residue mods first, then terminal deltas folded onto the end residues.
    let mut residue_delta = vec![0.0f64; n];
    let mut residue_modded = vec![false; n];
    for (&pos1, &delta) in res_mod_pos.iter().zip(res_mod_delta.iter()) {
        if pos1 < 1 || (pos1 as usize) > n {
            return Err(format!(
                "res_mod_pos {pos1} out of range for sequence of length {n}"
            )
            .into());
        }
        let idx = (pos1 - 1) as usize;
        residue_delta[idx] += delta;
        residue_modded[idx] = true;
    }
    if nterm_delta != 0.0 {
        residue_delta[0] += nterm_delta;
        residue_modded[0] = true;
    }
    if cterm_delta != 0.0 {
        residue_delta[n - 1] += cterm_delta;
        residue_modded[n - 1] = true;
    }

    let mut residues = Vec::with_capacity(n);
    for (i, &r) in bytes.iter().enumerate() {
        let aa = model::AminoAcid::standard(r)
            .ok_or_else(|| format!("non-standard residue {:?} at position {}", r as char, i + 1))?;
        if residue_modded[i] {
            let m = Modification {
                name: "msnet".to_string(),
                mass_delta: residue_delta[i],
                residue: ResidueSpec::Specific(r),
                location: ModLocation::Anywhere,
                fixed: false,
                accession: None,
                neutral_losses: Vec::new(),
                loss_class: 0,
            };
            residues.push(aa.with_mod(Arc::new(m)));
        } else {
            residues.push(aa);
        }
    }

    // Flanking residues per the spec: pre=`_`, post=`-`.
    Ok(model::peptide::Peptide::new(residues, b'_', b'-'))
}

/// Read one flat training parquet into a vector of [`MsnetPsm`].
///
/// Reads via the workspace `parquet`/`arrow` crates in record-batch chunks.
/// List columns (`res_mod_pos`, `res_mod_delta`, `mz`, `intensity`) are
/// decoded per-row from their `ListArray` offsets. Peaks are sorted ascending
/// by m/z (the parquet stores acquisition order).
fn read_msnet_parquet(path: &Path) -> Result<Vec<MsnetPsm>, Box<dyn std::error::Error>> {
    use arrow::array::{Array, Float32Array, Float64Array, Int32Array, ListArray, StringArray};
    use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

    let file = File::open(path)
        .map_err(|e| format!("opening {}: {e}", path.display()))?;
    let reader = ParquetRecordBatchReaderBuilder::try_new(file)
        .map_err(|e| format!("parquet reader for {}: {e}", path.display()))?
        .build()
        .map_err(|e| format!("building parquet reader for {}: {e}", path.display()))?;

    let mut out = Vec::new();
    for batch in reader {
        let batch = batch.map_err(|e| format!("reading batch from {}: {e}", path.display()))?;

        let col = |name: &str| -> Result<&std::sync::Arc<dyn Array>, Box<dyn std::error::Error>> {
            batch
                .column_by_name(name)
                .ok_or_else(|| format!("missing column '{name}' in {}", path.display()).into())
        };

        let seq = col("seq")?.as_any().downcast_ref::<StringArray>()
            .ok_or("column 'seq' is not a STRING column")?;
        let charge = col("charge")?.as_any().downcast_ref::<Int32Array>()
            .ok_or("column 'charge' is not an INT32 column")?;
        let prec_mz = col("prec_mz")?.as_any().downcast_ref::<Float64Array>()
            .ok_or("column 'prec_mz' is not a DOUBLE column")?;
        let res_mod_pos = col("res_mod_pos")?.as_any().downcast_ref::<ListArray>()
            .ok_or("column 'res_mod_pos' is not a LIST column")?;
        let res_mod_delta = col("res_mod_delta")?.as_any().downcast_ref::<ListArray>()
            .ok_or("column 'res_mod_delta' is not a LIST column")?;
        let nterm = col("nterm_delta")?.as_any().downcast_ref::<Float64Array>()
            .ok_or("column 'nterm_delta' is not a DOUBLE column")?;
        let cterm = col("cterm_delta")?.as_any().downcast_ref::<Float64Array>()
            .ok_or("column 'cterm_delta' is not a DOUBLE column")?;
        let mz = col("mz")?.as_any().downcast_ref::<ListArray>()
            .ok_or("column 'mz' is not a LIST column")?;
        let intensity = col("intensity")?.as_any().downcast_ref::<ListArray>()
            .ok_or("column 'intensity' is not a LIST column")?;

        // Helper to pull a Vec<i32> out of one ListArray row.
        let list_i32 = |list: &ListArray, i: usize| -> Result<Vec<i32>, Box<dyn std::error::Error>> {
            if list.is_null(i) {
                return Ok(Vec::new());
            }
            let v = list.value(i);
            let a = v.as_any().downcast_ref::<Int32Array>()
                .ok_or("list element is not INT32")?;
            Ok((0..a.len()).map(|j| a.value(j)).collect())
        };
        let list_f64 = |list: &ListArray, i: usize| -> Result<Vec<f64>, Box<dyn std::error::Error>> {
            if list.is_null(i) {
                return Ok(Vec::new());
            }
            let v = list.value(i);
            let a = v.as_any().downcast_ref::<Float64Array>()
                .ok_or("list element is not DOUBLE")?;
            Ok((0..a.len()).map(|j| a.value(j)).collect())
        };
        let list_f32 = |list: &ListArray, i: usize| -> Result<Vec<f32>, Box<dyn std::error::Error>> {
            if list.is_null(i) {
                return Ok(Vec::new());
            }
            let v = list.value(i);
            let a = v.as_any().downcast_ref::<Float32Array>()
                .ok_or("list element is not FLOAT")?;
            Ok((0..a.len()).map(|j| a.value(j)).collect())
        };

        for i in 0..batch.num_rows() {
            let seq_s = seq.value(i);
            let ch = charge.value(i);
            if !(1..=255).contains(&ch) {
                return Err(format!("invalid charge {ch} at row {i} of {}", path.display()).into());
            }
            let charge_u8 = ch as u8;

            let positions = list_i32(res_mod_pos, i)?;
            let deltas = list_f64(res_mod_delta, i)?;
            let peptide = build_msnet_peptide(
                seq_s,
                &positions,
                &deltas,
                nterm.value(i),
                cterm.value(i),
            )?;

            let mzs = list_f32(mz, i)?;
            let ints = list_f32(intensity, i)?;
            if mzs.len() != ints.len() {
                return Err(format!(
                    "mz ({}) and intensity ({}) length mismatch at row {i} of {}",
                    mzs.len(),
                    ints.len(),
                    path.display()
                )
                .into());
            }
            let mut peaks: Vec<(f64, f32)> =
                mzs.iter().zip(ints.iter()).map(|(&m, &it)| (m as f64, it)).collect();
            // Input is acquisition order; the scoring path requires ascending m/z.
            peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

            let spectrum = Spectrum {
                title: format!("row{i}"),
                precursor_mz: prec_mz.value(i),
                precursor_intensity: None,
                precursor_charge: Some(charge_u8 as i32),
                rt_seconds: None,
                scan: None,
                peaks,
                activation_method: None,
                isolation_lower_offset: None,
                isolation_upper_offset: None,
            };

            out.push(MsnetPsm { spectrum, peptide, charge: charge_u8 });
        }
    }
    Ok(out)
}

// train-intensity: merge partial intensity stats into a finalized model parquet
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct IntensityAggKey {
    ion_type: String,
    flank_n: String,
    flank_c: String,
    pos_bin: i32,
    charge: i32,
    nce_bin: String,
}

#[derive(Debug, Clone, Default)]
struct IntensityAggStats {
    count: i64,
    sum_log_rel: f64,
    sum_log_rel_sq: f64,
}

fn read_intensity_partial(path: &Path) -> Result<Vec<(IntensityAggKey, IntensityAggStats)>, Box<dyn std::error::Error>> {
    use arrow::array::{Array, Float64Array, Int32Array, Int64Array, StringArray};
    use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

    let file = std::fs::File::open(path)
        .map_err(|e| format!("open {}: {e}", path.display()))?;
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
                let mean = mean_col.unwrap().value(i);
                let var = var_col.unwrap().value(i);
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
fn finalize_intensity_stats(sum_log_rel: f64, sum_log_rel_sq: f64, count: i64) -> Option<(f64, f64)> {
    if count <= 0 {
        return None;
    }
    let n = count as f64;
    let mean = sum_log_rel / n;
    let var = (sum_log_rel_sq / n - mean * mean).max(0.0);
    Some((mean, var))
}

fn write_intensity_model(
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
        (&a.0.ion_type, &a.0.flank_n, &a.0.flank_c, a.0.pos_bin, a.0.charge, &a.0.nce_bin)
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
        .map(|(_, s)| finalize_intensity_stats(s.sum_log_rel, s.sum_log_rel_sq, s.count).unwrap().0)
        .collect();
    let var: Vec<_> = keys
        .iter()
        .map(|(_, s)| finalize_intensity_stats(s.sum_log_rel, s.sum_log_rel_sq, s.count).unwrap().1)
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

fn sanity_check_intensity_model(model: &scoring_crate::IntensityModel) -> Result<(), Box<dyn std::error::Error>> {
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
fn run_train_intensity(args: TrainIntensityArgs) -> Result<(), Box<dyn std::error::Error>> {
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

    write_intensity_model(&args.out, &merged)?;
    eprintln!(
        "train-intensity: wrote {} keys -> {}",
        merged.len(),
        args.out.display()
    );

    let model = IntensityModel::load(&args.out)?;
    sanity_check_intensity_model(&model)?;
    eprintln!("train-intensity: done in {:.1}s", t0.elapsed().as_secs_f64());
    Ok(())
}

/// `andes train-from-msnet`: train a scoring model directly from
/// externally-labeled PSM parquets, reusing the existing
/// accumulate → estimate → store machinery but bypassing the bootstrap search.
fn run_train_from_msnet(
    args: TrainFromMsnetArgs,
) -> Result<(), Box<dyn std::error::Error>> {
    use rayon::prelude::*;

    let t0 = std::time::Instant::now();

    // ── 1. Configure Rayon thread pool ────────────────────────────────────────
    static POOL_INIT: std::sync::OnceLock<()> = std::sync::OnceLock::new();
    POOL_INIT.get_or_init(|| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .expect("build_global");
    });

    // ── 2. Read all input parquets ────────────────────────────────────────────
    let mut psms: Vec<MsnetPsm> = Vec::new();
    let mut rows_read = 0usize;
    for input in &args.inputs {
        eprintln!("train-from-msnet: reading {} ...", input.display());
        let part = read_msnet_parquet(input)?;
        rows_read += part.len();
        eprintln!("train-from-msnet:   {} PSM rows", part.len());
        psms.extend(part);
    }
    if psms.is_empty() {
        return Err("no PSM rows read from any --in parquet".into());
    }
    eprintln!("train-from-msnet: {rows_read} total PSM rows across {} file(s)", args.inputs.len());

    // ── 3. Load seed Param and apply the fragment-tolerance override ──────────
    let (seed_model_id, mut seed_param): (String, Param) =
        load_seed_param(&Some(args.seed_model.clone()))?;
    eprintln!("train-from-msnet: seed model = {seed_model_id}");
    if let Some(ppm) = args.fragment_tol_ppm {
        seed_param.mme = Tolerance::Ppm(ppm);
        eprintln!("train-from-msnet: fragment tolerance overridden to {ppm} ppm");
    } else if let Some(da) = args.fragment_tol_da {
        seed_param.mme = Tolerance::Da(da);
        eprintln!("train-from-msnet: fragment tolerance overridden to {da} Da");
    } else {
        eprintln!("train-from-msnet: using seed fragment tolerance {:?}", seed_param.mme);
    }

    // Build the scorer AFTER the tolerance override so accumulation uses it.
    let seed_scorer = RankScorer::new(&seed_param);

    // ── 4. Accumulate ion-match statistics (parallel; per-worker CountStats) ──
    eprintln!("train-from-msnet: accumulating ion-match statistics ...");
    let stats = psms
        .par_iter()
        .fold(
            CountStats::new,
            |mut acc, psm| {
                let accumulator = StatsAccumulator::new(&seed_scorer);
                accumulator.accumulate(&mut acc, &psm.spectrum, &psm.peptide, psm.charge);
                acc
            },
        )
        .collect::<Vec<_>>();
    let stats = merge(stats);
    eprintln!("train-from-msnet: accumulated {} PSMs", psms.len());

    // ── 5. Estimate the model (replaces all learned tables in the seed) ───────
    eprintln!("train-from-msnet: estimating model parameters ...");
    let cfg = EstimatorConfig {
        pseudo: args.train_pseudo,
        noise_pseudo: args.train_noise_pseudo,
        min_count: args.train_min_count,
        backoff_weight: args.train_backoff_weight,
        error_scaling_factor_override: None,
        rank_smoothing: args.rank_smoothing,
    };
    eprintln!(
        "train-from-msnet: estimator pseudo={} noise_pseudo={} backoff_weight={} min_count={}",
        cfg.pseudo, cfg.noise_pseudo, cfg.backoff_weight, cfg.min_count
    );
    let estimator = Estimator::new(cfg);

    // Optional independent prior: sparse partitions shrink toward this model
    // (Level 0 of the backoff hierarchy) instead of the corpus-internal pool.
    // When `--prior-model` is omitted, default to the trained model id. The
    // selection columns passed to `load_param_from_store` are inert here because
    // `model_id_override` is `Some` (it loads that exact id).
    let prior_param: Option<Param> = match &args.prior_model_store {
        Some(store_path) => {
            let prior_id = args.prior_model.clone().unwrap_or_else(|| args.model_id.clone());
            // `load_param_from_store`'s activation/instrument/protocol are only
            // consulted for automatic selection; passing an explicit
            // `model_id_override` makes them inert, so `Protocol::Auto` (the CLI
            // enum the signature expects) is a harmless placeholder.
            let (_pid, p) = load_param_from_store(
                seed_param.data_type.activation,
                Some(seed_param.data_type.instrument),
                Protocol::Auto,
                model::enzyme::Enzyme::Trypsin, // inert (model_id_override below makes selection columns unused)
                Some(store_path.as_path()),
                Some(&prior_id),
            )
            .map_err(|e| format!("loading --prior-model '{prior_id}': {e}"))?;
            eprintln!("train-from-msnet: prior model = {prior_id} (from {})", store_path.display());
            Some(p)
        }
        None => None,
    };

    let mut trained_param =
        estimator.estimate_with_prior(&stats, &seed_param, prior_param.as_ref());
    let n_partitions = trained_param.partitions.len();
    eprintln!("train-from-msnet: trained model has {n_partitions} partitions");

    // ── 5b. Override the selection-relevant data_type from flags ──────────────
    // The trained model inherits the seed's data_type; minting a NEW slug whose
    // (activation, instrument, enzyme, protocol) differs from the seed requires
    // overriding those columns explicitly, otherwise model selection (which keys
    // on these columns, not the model_id string) would never route to it.
    if let Some(act) = &args.activation {
        trained_param.data_type.activation = ActivationMethod::from_name(act)
            .ok_or_else(|| format!("unknown --activation '{act}' (expected CID/HCD/ETD/UVPD/PQD)"))?;
    }
    if let Some(inst) = &args.instrument {
        trained_param.data_type.instrument = InstrumentType::from_name(inst)
            .ok_or_else(|| format!("unknown --instrument '{inst}' (expected LowRes/HighRes/QExactive/TOF)"))?;
    }
    if let Some(enz) = &args.enzyme {
        trained_param.data_type.enzyme = Some(
            model::enzyme::Enzyme::from_name(enz)
                .ok_or_else(|| format!("unknown --enzyme '{enz}' (e.g. Trypsin/LysC/LysN/AspN/GluC/ArgC)"))?,
        );
    }
    if let Some(prot) = &args.protocol {
        trained_param.data_type.protocol = model::protocol::Protocol::from_name(prot)
            .ok_or_else(|| format!("unknown --protocol '{prot}' (expected Automatic/TMT/iTRAQ/iTRAQPhospho/Phosphorylation/Standard)"))?;
    }
    eprintln!(
        "train-from-msnet: model data_type = {:?}/{:?}/{:?}/{:?}",
        trained_param.data_type.activation,
        trained_param.data_type.instrument,
        trained_param.data_type.enzyme,
        trained_param.data_type.protocol,
    );

    // ── 6. Build the source ledger (sentinel train_fdr; pre-labeled input) ────
    let dataset = args
        .inputs
        .first()
        .and_then(|p| p.file_name())
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_else(|| "msnet".to_string());
    let ledger = SourceLedger {
        source_id: "msnet".to_string(),
        dataset,
        n_psms: psms.len() as i64,
        date: format_today_iso8601(),
        weight: 1.0,
        train_fdr: 1.0, // sentinel: input is pre-labeled, no q-value filtering here
        instrument: trained_param.data_type.instrument.name().to_string(),
        experiment_class: trained_param.data_type.protocol.name().to_string(),
    };

    // ── 7. Write to store, preserving any other existing models ───────────────
    let store_path = &args.out_store;
    let model_id = args.model_id.clone();
    let mut existing_other: Vec<ModelEntryOwned> = Vec::new();
    if store_path.exists() {
        let store = ModelStore::open(store_path)
            .map_err(|e| format!("opening existing store {}: {e}", store_path.display()))?;
        for id in store.model_ids() {
            if id == model_id {
                eprintln!("train-from-msnet: overwriting existing model '{id}' in store");
                continue;
            }
            let p = store.load_param(&id)
                .map_err(|e| format!("reading model '{id}': {e}"))?;
            let src_ledgers = store.load_sources(&id).unwrap_or_default();
            let mut src = Vec::new();
            for l in src_ledgers {
                if let Ok(s) = store.load_source_stats(&id, &l.source_id) {
                    src.push((l, s));
                }
            }
            existing_other.push((id, p, src));
        }
    }

    let mut all_entries: Vec<ModelEntryOwned> = Vec::new();
    all_entries.push((model_id.clone(), trained_param, vec![(ledger, stats)]));
    for (id, p, src) in existing_other {
        all_entries.push((id, p, src));
    }

    write_all_models_with_sources_pub(
        store_path,
        &all_entries.iter()
            .map(|(id, p, s)| (id.as_str(), p, s.as_slice()))
            .collect::<Vec<_>>(),
    )
    .map_err(|e| format!("writing model store {}: {e}", store_path.display()))?;

    eprintln!(
        "train-from-msnet: wrote model '{model_id}' to {} (source 'msnet', {} PSMs, {n_partitions} partitions) [{:.2}s]",
        store_path.display(),
        psms.len(),
        t0.elapsed().as_secs_f64(),
    );

    Ok(())
}

/// Incremental update mode (Part D): `--update <MODEL_ID>` plus one of
/// `--add`, `--remove-source`, `--reweight`, `--decay`.
fn run_train_update(
    args: TrainArgs,
    model_id: &str,
    t0: std::time::Instant,
) -> Result<(), Box<dyn std::error::Error>> {
    let store_path = &args.out_store;
    let cfg = EstimatorConfig::default();

    // ── Dispatch to the right update operation ────────────────────────────────
    let (candidate, new_sources) = if args.add {
        // --add mode: search spectra, accumulate stats, call update_add.
        let spectra_path = args.spectra.clone()
            .ok_or("--spectra is required with --add")?;
        let database = args.database.clone()
            .ok_or("--database is required with --add")?;
        let source_id = args.source_id.clone()
            .ok_or("--source-id is required with --add")?;

        eprintln!("train update: loading spectra from {} ...", spectra_path.display());
        let spectra = load_spectra_for_train(&spectra_path)?;
        eprintln!("train update: loaded {} spectra", spectra.len());

        // Load the current stored model as the seed.
        let store = ModelStore::open(store_path)
            .map_err(|e| format!("opening store {}: {e}", store_path.display()))?;
        let current_param = store.load_param(model_id)
            .map_err(|e| format!("loading model '{model_id}': {e}"))?;
        let current_scorer = RankScorer::new(&current_param);

        let search_params = build_train_search_params(&args.mods)?;

        eprintln!("train update: running seed search (train-fdr={}) ...", args.train_fdr);
        let labels = bootstrap_labels(
            &spectra,
            &database,
            &current_scorer,
            &search_params,
            args.train_fdr,
        )
        .map_err(|e| format!("bootstrap_labels: {e}"))?;
        eprintln!("train update: {} confident labels", labels.len());

        if labels.is_empty() {
            return Err(format!(
                "no confident labels at train-fdr={} — try a higher --train-fdr",
                args.train_fdr
            ).into());
        }

        let accumulator = StatsAccumulator::new(&current_scorer);
        let mut stats = CountStats::new();
        for label in &labels {
            accumulator.accumulate(&mut stats, &spectra[label.spectrum_index], &label.peptide, label.charge);
        }
        let stats = merge(vec![stats]);

        let spectra_filename = spectra_path
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_else(|| spectra_path.display().to_string());
        let date_str = args.date.clone().unwrap_or_else(format_today_iso8601);
        let ledger = SourceLedger {
            source_id: source_id.clone(),
            dataset: spectra_filename,
            n_psms: labels.len() as i64,
            date: date_str,
            weight: 1.0,
            train_fdr: args.train_fdr as f32,
            instrument: args.instrument.clone(),
            experiment_class: args.protocol.clone(),
        };

        update_add(store_path, model_id, ledger, stats, cfg)
            .map_err(|e| format!("update_add: {e}"))?

    } else if let Some(ref sid) = args.remove_source.clone() {
        update_remove(store_path, model_id, sid, cfg)
            .map_err(|e| format!("update_remove: {e}"))?

    } else if let Some(ref spec) = args.reweight.clone() {
        // Parse "source-id=weight"
        let (sid, weight) = parse_reweight_spec(spec)?;
        update_reweight(store_path, model_id, &sid, weight, cfg)
            .map_err(|e| format!("update_reweight: {e}"))?

    } else if let Some(half_life) = args.decay {
        update_decay(store_path, model_id, half_life, cfg)
            .map_err(|e| format!("update_decay: {e}"))?

    } else {
        return Err(
            "update mode requires one of: --add, --remove-source, --reweight, --decay".into()
        );
    };

    // ── Acceptance gate (Part D) ──────────────────────────────────────────────
    let commit = if let Some(ref validate_path) = args.validate.clone() {
        let database = args.database.clone()
            .ok_or("--database is required with --validate")?;

        eprintln!("train update: running acceptance gate on {} ...", validate_path.display());
        let val_spectra = load_spectra_for_train(validate_path)?;

        let store = ModelStore::open(store_path)
            .map_err(|e| format!("opening store for gate: {e}"))?;
        let current_param = store.load_param(model_id)
            .map_err(|e| format!("loading current model for gate: {e}"))?;
        let current_scorer = RankScorer::new(&current_param);
        let candidate_scorer = RankScorer::new(&candidate);

        let search_params = build_train_search_params(&args.mods)?;

        let delta = evaluate_candidate(
            &val_spectra,
            &database,
            &current_scorer,
            &candidate_scorer,
            &search_params,
            args.train_fdr,
        )
        .map_err(|e| format!("evaluate_candidate: {e}"))?;

        eprintln!(
            "train update: gate — current={} PSMs, candidate={} PSMs at FDR={}",
            delta.current_count, delta.candidate_count, args.train_fdr
        );

        if delta.is_accepted() {
            eprintln!("train update: ACCEPTED (candidate >= current)");
            true
        } else {
            eprintln!("train update: REJECTED (candidate < current)");
            if args.force {
                eprintln!("train update: --force set, committing anyway");
                true
            } else {
                eprintln!("train update: skipping commit (use --force to override)");
                false
            }
        }
    } else {
        eprintln!("train update: no --validate dataset; skipping acceptance gate");
        if args.force {
            eprintln!("train update: --force set, committing unconditionally");
        }
        // Without --validate, commit unless user explicitly uses --force to control.
        // Default: commit (no gate run = no evidence of regression).
        true
    };

    if commit {
        commit_update(store_path, model_id, &candidate, &new_sources)
            .map_err(|e| format!("commit_update: {e}"))?;
        eprintln!(
            "train update: committed model '{model_id}' to {} [{:.2}s]",
            store_path.display(),
            t0.elapsed().as_secs_f64(),
        );
    }

    Ok(())
}

/// Parse `"source-id=weight"` from a `--reweight` argument.
fn parse_reweight_spec(spec: &str) -> Result<(String, f32), Box<dyn std::error::Error>> {
    let pos = spec.rfind('=').ok_or_else(|| {
        format!("--reweight value must be <source-id>=<weight>, got '{spec}'")
    })?;
    let sid = spec[..pos].to_string();
    let weight: f32 = spec[pos + 1..].parse()
        .map_err(|e| format!("invalid weight in --reweight '{spec}': {e}"))?;
    Ok((sid, weight))
}

/// Load the seed Param from the optional seed model specifier.
fn load_seed_param(seed_model: &Option<String>) -> Result<(String, Param), Box<dyn std::error::Error>> {
    match seed_model {
        None => {
            let store_path = bundled_store_path();
            let store = ModelStore::open(&store_path)
                .map_err(|e| format!("opening bundled store: {e}"))?;
            let p = store.load_param("hcd_qexactive_tryp")
                .map_err(|e| format!("loading seed model: {e}"))?;
            Ok(("hcd_qexactive_tryp".to_string(), p))
        }
        Some(seed) => {
            let as_path = Path::new(seed);
            if as_path.is_file() {
                let p = Param::load_from_file(as_path)
                    .map_err(|e| format!("loading seed param file {}: {e}", as_path.display()))?;
                Ok((seed.clone(), p))
            } else {
                let store_path = bundled_store_path();
                let store = ModelStore::open(&store_path)
                    .map_err(|e| format!("opening bundled store: {e}"))?;
                let p = store.load_param(seed)
                    .map_err(|e| format!("loading seed model '{seed}': {e}"))?;
                Ok((seed.clone(), p))
            }
        }
    }
}

/// Build an `AminoAcidSet` from an optional mods file, defaulting to
/// Carbamidomethyl-C fixed + Oxidation-M variable.
fn build_aa_set(
    mods: &Option<PathBuf>,
) -> Result<model::AminoAcidSet, Box<dyn std::error::Error>> {
    match mods {
        Some(path) => {
            let set = AminoAcidSetBuilder::new_standard()
                .add_mods_from_file(path)
                .map_err(|e| format!("loading mods from {}: {e}", path.display()))?
                .build()
                .map_err(|e| format!("building amino-acid set: {e}"))?;
            Ok(set)
        }
        None => {
            let cam = Modification {
                name: "Carbamidomethyl".into(),
                mass_delta: 57.02146,
                residue: ResidueSpec::Specific(b'C'),
                location: ModLocation::Anywhere,
                fixed: true,
                accession: None,
                neutral_losses: Vec::new(),
                loss_class: 0,
            };
            let ox = Modification {
                name: "Oxidation".into(),
                mass_delta: 15.99491,
                residue: ResidueSpec::Specific(b'M'),
                location: ModLocation::Anywhere,
                fixed: false,
                accession: None,
                neutral_losses: Vec::new(),
                loss_class: 0,
            };
            Ok(AminoAcidSetBuilder::new_standard()
                .add_fixed_mod(cam)
                .add_variable_mod(ox)
                .build()?)
        }
    }
}

/// Build the default `AminoAcidSet` (CAM-C fixed + Ox-M variable), optionally
/// with an isobaric tag (TMT/iTRAQ) as a FIXED mod on K + peptide N-term.
///
/// Used by the no-`--mods` (parameter-free) path: when the protocol resolves to
/// TMT/iTRAQ the tag MUST be in the candidate set, or every labeled peptide is
/// `+tag·(nK+1)` Da off at the precursor and silently misses (C1).
fn default_aa_set_with_tag(
    tag: Option<(&str, f64)>,
) -> Result<model::AminoAcidSet, Box<dyn std::error::Error>> {
    let cam = Modification {
        name: "Carbamidomethyl".into(),
        mass_delta: 57.02146,
        residue: ResidueSpec::Specific(b'C'),
        location: ModLocation::Anywhere,
        fixed: true,
        accession: None,
        neutral_losses: Vec::new(),
        loss_class: 0,
    };
    let ox = Modification {
        name: "Oxidation".into(),
        mass_delta: 15.99491,
        residue: ResidueSpec::Specific(b'M'),
        location: ModLocation::Anywhere,
        fixed: false,
        accession: None,
        neutral_losses: Vec::new(),
        loss_class: 0,
    };
    let mut b = AminoAcidSetBuilder::new_standard()
        .add_fixed_mod(cam)
        .add_variable_mod(ox);
    if let Some((name, mass)) = tag {
        let tag_k = Modification {
            name: name.into(),
            mass_delta: mass,
            residue: ResidueSpec::Specific(b'K'),
            location: ModLocation::Anywhere,
            fixed: true,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
        };
        let tag_nterm = Modification {
            name: name.into(),
            mass_delta: mass,
            residue: ResidueSpec::Wildcard,
            location: ModLocation::NTerm,
            fixed: true,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
        };
        b = b.add_fixed_mod(tag_k).add_fixed_mod(tag_nterm);
    }
    Ok(b.build()?)
}

/// Format today's date as `YYYY-MM-DD` using `std::time::SystemTime`.
fn format_today_iso8601() -> String {
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    // Simple Gregorian calendar conversion from Unix timestamp (days since epoch).
    let days = secs / 86400;
    unix_days_to_iso8601(days)
}

fn unix_days_to_iso8601(days: u64) -> String {
    // Algorithm: Gregorian calendar from Julian Day Number.
    // JDN for 1970-01-01 = 2440588.
    let jdn = days as i64 + 2_440_588;
    let a = jdn + 32044;
    let b = (4 * a + 3) / 146097;
    let c = a - (146097 * b) / 4;
    let d = (4 * c + 3) / 1461;
    let e = c - (1461 * d) / 4;
    let m = (5 * e + 2) / 153;
    let day = e - (153 * m + 2) / 5 + 1;
    let month = m + 3 - 12 * (m / 10);
    let year = 100 * b + d - 4800 + m / 10;
    format!("{:04}-{:02}-{:02}", year, month, day)
}

/// Convert the CLI `Fragmentation` enum to `Option<ActivationMethod>`.
///
/// `Fragmentation::Auto` returns `None` (no activation explicitly requested);
/// every concrete variant maps to its `ActivationMethod`. Used by
/// [`resolve_metadataless_selection`] so that an unset `--fragmentation`
/// defers to detection or the class-consistent default.
fn cli_fragmentation_to_activation_opt(f: Fragmentation) -> Option<ActivationMethod> {
    match f {
        Fragmentation::Auto => None,
        Fragmentation::Cid  => Some(ActivationMethod::CID),
        Fragmentation::Etd  => Some(ActivationMethod::ETD),
        Fragmentation::Hcd  => Some(ActivationMethod::HCD),
        Fragmentation::Uvpd => Some(ActivationMethod::UVPD),
    }
}

/// Resolve the CLI fragment-tolerance override (MGF only) into a `Tolerance`.
/// `--fragment-tol-ppm` ⇒ `Ppm`; `--fragment-tol-da` ⇒ `Da`; none ⇒ `None`.
fn cli_fragment_tol_override(
    fragment_tol_ppm: Option<f64>,
    fragment_tol_da: Option<f64>,
) -> Option<model::tolerance::Tolerance> {
    use model::tolerance::Tolerance;
    fragment_tol_ppm
        .map(Tolerance::Ppm)
        .or_else(|| fragment_tol_da.map(Tolerance::Da))
}

/// Resolve (activation, instrument) for model selection on metadata-less input
/// (MGF, or mzML/.raw with no analyzer metadata). Resolution class comes from
/// the `--fragment-tol-*` unit; activation from detected method, else
/// `--fragmentation`, else the class-consistent default. When nothing
/// disambiguates, decision E: CID / LowRes (→ `cid_lowres_tryp`) + a warning.
fn resolve_metadataless_selection(
    detected_activation: Option<ActivationMethod>,
    fragmentation: Fragmentation,
    fragment_tol_ppm: Option<f64>,
    fragment_tol_da: Option<f64>,
) -> (ActivationMethod, Option<InstrumentType>) {
    let instrument: Option<InstrumentType> = if fragment_tol_ppm.is_some() {
        Some(InstrumentType::QExactive)
    } else if fragment_tol_da.is_some() {
        Some(InstrumentType::LowRes)
    } else {
        None
    };
    let explicit = cli_fragmentation_to_activation_opt(fragmentation);
    // Class-consistent default when neither detection nor `--fragmentation`
    // names an activation: high-res classes imply HCD, otherwise CID.
    let class_default = match instrument {
        Some(InstrumentType::QExactive)
        | Some(InstrumentType::HighRes)
        | Some(InstrumentType::TOF) => ActivationMethod::HCD,
        _ => ActivationMethod::CID,
    };
    let activation = detected_activation.or(explicit).unwrap_or(class_default);
    if detected_activation.is_none() && explicit.is_none() && instrument.is_none() {
        eprintln!(
            "WARN: MGF input with no --fragmentation/--fragment-tol; assuming \
             CID / low-res / 0.5 Da. Pass --fragmentation and --fragment-tol-ppm/-da \
             to override."
        );
    }
    (activation, instrument)
}

/// Peek the spectrum file and return the dominant
/// `ActivationMethod` across the first several MS2 spectra.
///
/// Reads up to `MAX_PEEK` spectra (early-exit) and tallies a histogram of
/// activation methods. Returns the most-common method, or `None` when no
/// spectra carry an activation cvParam (older mzMLs, MGF, etc.).
///
/// Currently only mzML files (`.mzml` / `.mzML` extension) carry an
/// `<activation>` block. For anything else (MGF, unknown extension) we
/// return `None` and the caller falls back to the historical default.
///
/// When multiple activation methods are present, prints a single
/// `eprintln!` warning naming the runner-up and its count.
fn detect_dominant_activation(spectrum_path: &std::path::Path) -> Option<ActivationMethod> {
    // Only mzML carries `<activation>`. Other formats: caller falls back.
    let ext_lower = spectrum_path
        .extension()
        .and_then(|s| s.to_str())
        .map(|s| s.to_ascii_lowercase());
    if ext_lower.as_deref() != Some("mzml") {
        return None;
    }

    const MAX_PEEK: usize = 64;

    let file = File::open(spectrum_path).ok()?;
    let reader = MzMLReader::new(BufReader::new(file));

    // Tally counts keyed by ActivationMethod variant.
    let mut counts: std::collections::HashMap<ActivationMethod, usize> =
        std::collections::HashMap::new();
    for (seen, item) in reader.enumerate() {
        if seen >= MAX_PEEK {
            break;
        }
        if let Ok(spec) = item {
            if let Some(m) = spec.activation_method {
                *counts.entry(m).or_insert(0) += 1;
            }
        }
    }

    if counts.is_empty() {
        return None;
    }

    // Find the dominant method. Ties are broken by ActivationMethod's
    // declaration order via match below, which is stable.
    let dominant = counts
        .iter()
        .max_by_key(|(_, &n)| n)
        .map(|(&m, _)| m)?;

    // Warn on mixed activation. The dominant method still wins; this is
    // purely informational so the user can spot heterogeneous mzMLs.
    if counts.len() > 1 {
        let mut other_pairs: Vec<(ActivationMethod, usize)> = counts
            .iter()
            .filter(|(&m, _)| m != dominant)
            .map(|(&m, &n)| (m, n))
            .collect();
        other_pairs.sort_by(|a, b| b.1.cmp(&a.1));
        let total: usize = counts.values().sum();
        let dominant_count = counts[&dominant];
        eprintln!(
            "Param resolver: mixed activation methods in input ({} different methods \
             across {} peeked MS2 spectra). Using dominant = {} ({}/{}); other methods \
             present: {}",
            counts.len(),
            total,
            dominant.name(),
            dominant_count,
            total,
            other_pairs
                .iter()
                .map(|(m, n)| format!("{}={}", m.name(), n))
                .collect::<Vec<_>>()
                .join(", "),
        );
    }

    Some(dominant)
}

/// Helper to call `input::detect_instrument_type` on an mzML path.
///
/// Mirrors the structure of `detect_dominant_activation` so the two
/// detection passes look symmetric at the call site. Returns `None` for
/// non-mzML inputs or when the mzML has no recoverable instrument metadata.
fn detect_instrument_type_for_path(spectrum_path: &std::path::Path) -> Option<InstrumentType> {
    let ext_lower = spectrum_path
        .extension()
        .and_then(|s| s.to_str())
        .map(|s| s.to_ascii_lowercase());
    if ext_lower.as_deref() != Some("mzml") {
        return None;
    }

    let file = File::open(spectrum_path).ok()?;
    detect_instrument_type(BufReader::new(file))
}

/// Resolve the path to the bundled `models.parquet` store.
///
/// A packaged release ships `resources/` next to the binary, so prefer
/// `<exe_dir>/resources/ionstat/models.parquet` when it exists — that makes an
/// installed binary self-contained regardless of where it runs. Fall back to the
/// compile-time source tree (`CARGO_MANIFEST_DIR`) for `cargo run` / tests.
fn bundled_store_path() -> PathBuf {
    if let Ok(exe) = std::env::current_exe() {
        if let Some(dir) = exe.parent() {
            let next_to_binary = dir.join("resources/ionstat/models.parquet");
            if next_to_binary.exists() {
                return next_to_binary;
            }
        }
    }
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../resources/ionstat/models.parquet")
}

/// Build a [`SelectionKey`] from `(activation, instrument, protocol)` applying
/// all old-ladder normalizations. This is the new entry point used by the
/// search binary, replacing the former filename-based resolution ladder.
///
/// `activation`: the detected or explicitly set `ActivationMethod`.
/// `instrument`: the detected or explicitly set `InstrumentType` (None = undetected → LowRes).
/// `protocol`:   the CLI `Protocol` value.
fn build_selection_key(
    activation: ActivationMethod,
    instrument: Option<InstrumentType>,
    protocol: Protocol,
    enzyme: model::enzyme::Enzyme,
) -> SelectionKey {
    use std::collections::BTreeSet;

    // 1. PQD → CID for model routing.
    let act_str: &str = match activation {
        ActivationMethod::PQD => "CID",
        other                 => other.name(),
    };
    // 2. Apply family fallback (OrbitrapAstral → QExactive, TimsTOF → TOF).
    let inst_after_family: &str = match instrument {
        Some(i) => i.family_fallback().name(),
        None    => "LowRes",
    };

    // 3. Apply old-ladder (activation, instrument) normalization.
    //    Because `normalize_for_store` returns `&'static str` only for the
    //    normalizing arms (to avoid lifetime issues), we handle identity
    //    inline here.
    let (final_act, final_inst, drop_protocol): (&str, &str, bool) =
        match (act_str, inst_after_family) {
            ("HCD", "LowRes")    => ("HCD", "QExactive", false),
            ("HCD", "TOF")       => ("CID", "TOF",       true),
            ("CID", "QExactive") => ("CID", "LowRes",    true),
            ("ETD", i) if !matches!(i, "LowRes" | "HighRes") => ("ETD", "LowRes", true),
            ("UVPD", i) if i != "QExactive" => ("CID", "LowRes", true),
            _ => (act_str, inst_after_family, false),
        };

    // 4. Build experiment_class from protocol (unless the final fallback dropped it).
    //    Protocol → experiment_class mapping matches the parquet's `protocol` column.
    let protocol_for_store: &str = match protocol {
        Protocol::Auto | Protocol::Standard => "Automatic",
        Protocol::Tmt          => "TMT",
        Protocol::Phospho      => "Phosphorylation",
        Protocol::Itraq        => "iTRAQ",
        Protocol::ItraqPhospho => "iTRAQPhospho",
    };
    let experiment_class: BTreeSet<String> = if drop_protocol {
        BTreeSet::new()
    } else {
        store_protocol_to_experiment_class(protocol_for_store)
    };

    SelectionKey {
        activation: final_act.to_string(),
        instrument: final_inst.to_string(),
        // Parquet stores the enzyme as its `Enzyme::name()` ("Trypsin", "LysC", ...).
        enzyme: enzyme.name().to_string(),
        experiment_class,
    }
}

/// Load the scoring [`Param`] from the bundled Parquet store for the given
/// `(activation, instrument, protocol)` combination.
///
/// This is the new model-resolution path, replacing the former
/// filename-based resolution ladder. The `model_id`
/// selected from the store will be identical to the lowercased filename
/// stem of the old `.param` path (guaranteed by the equivalence gate test
/// `store_selection_matches_old_ladder_for_all_combos`).
///
/// `custom_store_path`: when `Some`, use that Parquet file instead of the
/// bundled `resources/ionstat/models.parquet` (honours `--model-store`).
///
/// `model_id_override`: when `Some`, skip automatic selection and load this
/// exact model ID (honours `--model`).
fn load_param_from_store(
    activation: ActivationMethod,
    instrument: Option<InstrumentType>,
    protocol: Protocol,
    enzyme: model::enzyme::Enzyme,
    custom_store_path: Option<&Path>,
    model_id_override: Option<&str>,
) -> Result<(String, Param), Box<dyn std::error::Error>> {
    let store_path = custom_store_path
        .map(|p| p.to_owned())
        .unwrap_or_else(bundled_store_path);
    let store = ModelStore::open(&store_path)
        .map_err(|e| format!("opening model store {}: {e}", store_path.display()))?;

    let model_id: String = if let Some(id) = model_id_override {
        id.to_string()
    } else {
        let entries = store.selection_entries();
        let key = build_selection_key(activation, instrument, protocol, enzyme);

        // Forward-compat: `build_selection_key` collapses instruments with a real
        // family fallback (OrbitrapAstral → QExactive, TimsTOF → TOF) so the
        // bundled models resolve correctly. But that also hides a model trained
        // for the EXACT instrument (e.g. a user-trained OrbitrapAstral model).
        // Try the exact detected instrument FIRST; only when no such model exists
        // (the bundled case) do we fall through to the normalized family ladder —
        // so bundled selection (and the equivalence gate) is unchanged.
        let exact_id: Option<String> = match instrument {
            Some(i) if i.family_fallback().name() != i.name() => {
                let raw_key = SelectionKey {
                    instrument: i.name().to_string(),
                    ..key.clone()
                };
                select(&entries, &raw_key, |s| s.to_string(), None).map(|s| s.to_string())
            }
            _ => None,
        };

        exact_id.unwrap_or_else(|| {
            // `build_selection_key` already applies family fallback + all
            // normalizations, so the family_fn here is the identity. (L6) Pass
            // `None` for the generic so a true no-match returns `None`, letting us
            // WARN that the chosen model is a last-resort fallback rather than
            // silently emitting `hcd_qexactive_tryp` for mis-detected data.
            match select(&entries, &key, |i| i.to_string(), None) {
                Some(id) => id.to_string(),
                None => {
                    eprintln!(
                        "WARN: no model matched (activation={}, instrument={}, enzyme={}, class={:?}) \
                         — falling back to the generic 'hcd_qexactive_tryp'; scores may be \
                         mis-calibrated for this data. Pin a model with --model if this is wrong.",
                        key.activation, key.instrument, key.enzyme, key.experiment_class
                    );
                    "hcd_qexactive_tryp".to_string()
                }
            }
        })
    };

    let param = store.load_param(&model_id)
        .map_err(|e| format!("loading model '{model_id}' from store: {e}"))?;

    Ok((model_id, param))
}

/// Parse `--fragmentation` value. Accepts named (case-insensitive: auto, CID,
/// ETD, HCD, UVPD) or legacy numeric (0=Auto, 1=CID, 2=ETD, 3=HCD, 4=UVPD).
fn parse_fragmentation(s: &str) -> Result<Fragmentation, String> {
    if let Ok(v) = <Fragmentation as ValueEnum>::from_str(s, true) { return Ok(v); }
    match s.parse::<u8>() {
        Ok(0) => Ok(Fragmentation::Auto),
        Ok(1) => Ok(Fragmentation::Cid),
        Ok(2) => Ok(Fragmentation::Etd),
        Ok(3) => Ok(Fragmentation::Hcd),
        Ok(4) => Ok(Fragmentation::Uvpd),
        _ => Err(format!(
            "invalid fragmentation `{s}`: expected auto|CID|ETD|HCD|UVPD \
             (or legacy 0..=4)"
        )),
    }
}

/// Parse `--protocol` value. Accepts named or legacy numeric
/// (0=Auto, 1=Phospho, 2=iTRAQ, 3=iTRAQ-phospho, 4=TMT, 5=Standard).
fn parse_protocol(s: &str) -> Result<Protocol, String> {
    if let Ok(v) = <Protocol as ValueEnum>::from_str(s, true) { return Ok(v); }
    match s.parse::<u8>() {
        Ok(0) => Ok(Protocol::Auto),
        Ok(1) => Ok(Protocol::Phospho),
        Ok(2) => Ok(Protocol::Itraq),
        Ok(3) => Ok(Protocol::ItraqPhospho),
        Ok(4) => Ok(Protocol::Tmt),
        Ok(5) => Ok(Protocol::Standard),
        _ => Err(format!(
            "invalid --protocol `{s}`: valid range is 0..=5 \
             (0=Automatic, 1=Phosphorylation, 2=iTRAQ, 3=iTRAQPhospho, \
              4=TMT, 5=Standard) or named auto|phospho|iTRAQ|iTRAQ-phospho|TMT|standard"
        )),
    }
}

/// Parse `--enzyme-specificity` (`--ntt`) value. Accepts named
/// (non-specific, semi, fully) or legacy numeric (0=non-specific,
/// 1=semi, 2=fully).
fn parse_precursor_cal(s: &str) -> Result<PrecursorCalMode, String> {
    match s.to_ascii_lowercase().as_str() {
        "auto" => Ok(PrecursorCalMode::Auto),
        "on" => Ok(PrecursorCalMode::On),
        "off" => Ok(PrecursorCalMode::Off),
        _ => Err(format!(
            "invalid precursor-cal `{s}`: expected auto|on|off (Java -precursorCal)"
        )),
    }
}

fn parse_enzyme_specificity(s: &str) -> Result<EnzymeSpecificity, String> {
    if let Ok(v) = <EnzymeSpecificity as ValueEnum>::from_str(s, true) { return Ok(v); }
    match s.parse::<u8>() {
        Ok(0) => Ok(EnzymeSpecificity::NonSpecific),
        Ok(1) => Ok(EnzymeSpecificity::Semi),
        Ok(2) => Ok(EnzymeSpecificity::Fully),
        _ => Err(format!(
            "invalid enzyme specificity `{s}`: expected non-specific|semi|fully \
             (or legacy 0..=2)"
        )),
    }
}

#[cfg(test)]
mod param_resolver_tests {
    use super::*;

    // ── Tests of the resolve_bundled_param / resolve_bundled_param_for_activation
    //    ladder were removed in the model-store migration: those functions are now
    //    dead code (#[cfg_attr(not(test), allow(dead_code))]) and the bundled .param
    //    files no longer ship on disk (they live in resources/ionstat/models.parquet).
    //    The store_selection_equivalence integration test covers the same
    //    correctness invariant without requiring physical .param files.

    #[test]
    fn parse_fragmentation_rejects_out_of_range_numeric() {
        let err = parse_fragmentation("99").unwrap_err();
        assert!(err.contains("0..=4"), "error message should mention range, got: {err}");
    }

    #[test]
    fn parse_protocol_rejects_out_of_range_numeric() {
        let err = parse_protocol("99").unwrap_err();
        assert!(err.contains("0..=5"), "got: {err}");
    }

    #[test]
    fn parse_precursor_cal_accepts_named_modes() {
        assert_eq!(parse_precursor_cal("auto").unwrap(), PrecursorCalMode::Auto);
        assert_eq!(parse_precursor_cal("OFF").unwrap(), PrecursorCalMode::Off);
        assert_eq!(parse_precursor_cal("on").unwrap(), PrecursorCalMode::On);
        assert!(parse_precursor_cal("bogus").is_err());
    }

    #[test]
    fn finalize_intensity_stats_drops_zero_count() {
        // A count==0 key (e.g. from a partial aggregation parquet) must not
        // produce NaN mean/var in the finalized model — it carries no signal.
        assert_eq!(finalize_intensity_stats(0.0, 0.0, 0), None);
    }

    #[test]
    fn finalize_intensity_stats_computes_mean_and_clamped_var() {
        // sum=6, sum_sq=14, count=2 -> mean=3, var=max(0, 14/2 - 9)= max(0,-2)=0.
        let (mean, var) = finalize_intensity_stats(6.0, 14.0, 2).unwrap();
        assert_eq!(mean, 3.0);
        assert_eq!(var, 0.0);
    }

    #[test]
    fn single_file_search_does_not_prefix_titles() {
        // Regression (default-path PIN parity): with one input file the
        // SpecId/title must stay unprefixed so PIN output is byte-identical to
        // the single-file behavior that predates multi-`--spectrum` support.
        assert_eq!(title_prefix_for(1, "myfile"), None);
    }

    #[test]
    fn multi_file_search_prefixes_titles_with_file_stem() {
        assert_eq!(title_prefix_for(2, "myfile").as_deref(), Some("myfile/"));
    }
}

#[cfg(test)]
mod train_from_msnet_tests {
    use super::*;
    use model::mass::H2O;

    /// Reference: unmodified peptide mass = sum of residue masses + H2O.
    fn unmod_mass(seq: &[u8]) -> f64 {
        let rsum: f64 = seq
            .iter()
            .map(|&r| model::AminoAcid::standard(r).unwrap().mass)
            .sum();
        rsum + H2O
    }

    #[test]
    fn unmodified_peptide_mass_matches_reference() {
        let p = build_msnet_peptide("PEPTIDEK", &[], &[], 0.0, 0.0).unwrap();
        let expected = unmod_mass(b"PEPTIDEK");
        assert_eq!(p.mass().to_bits(), expected.to_bits());
        assert!(p.residues.iter().all(|aa| !aa.is_modified()));
    }

    /// Oxidation on residue 1 of "MPEPTIDE" must add exactly +15.994915 Da to
    /// the unmodified mass — verifying mods are actually applied, not dropped.
    #[test]
    fn residue_mod_mass_is_correct() {
        const OX: f64 = 15.994915;
        let p = build_msnet_peptide("MPEPTIDE", &[1], &[OX], 0.0, 0.0).unwrap();
        let expected = unmod_mass(b"MPEPTIDE") + OX;
        // Exact: Peptide::new sums aa.mass + mod_.mass_delta, so the delta is
        // added with the same arithmetic as the reference.
        assert_eq!(p.mass().to_bits(), expected.to_bits());
        assert!(p.residues[0].is_modified(), "residue 1 should carry the mod");
        assert_eq!(p.residues[0].mod_.as_ref().unwrap().mass_delta, OX);
    }

    /// N-terminal Acetyl folds onto the first residue's mass-delta.
    #[test]
    fn nterm_mod_mass_is_correct() {
        const ACETYL: f64 = 42.010565;
        let p = build_msnet_peptide("PEPTIDEK", &[], &[], ACETYL, 0.0).unwrap();
        let expected = unmod_mass(b"PEPTIDEK") + ACETYL;
        assert_eq!(p.mass().to_bits(), expected.to_bits());
        assert!(p.residues[0].is_modified());
    }

    /// A residue mod and a terminal mod on the SAME residue must sum (not
    /// clobber) so the total mass stays correct.
    #[test]
    fn overlapping_residue_and_nterm_mods_sum() {
        const OX: f64 = 15.994915;
        const ACETYL: f64 = 42.010565;
        // Oxidation on residue 1 (M) AND N-term acetyl on the same first residue.
        let p = build_msnet_peptide("MPEPTIDE", &[1], &[OX], ACETYL, 0.0).unwrap();
        let expected = unmod_mass(b"MPEPTIDE") + OX + ACETYL;
        assert_eq!(p.mass().to_bits(), expected.to_bits());
        assert_eq!(p.residues[0].mod_.as_ref().unwrap().mass_delta, OX + ACETYL);
    }

    #[test]
    fn mismatched_mod_arrays_error() {
        let r = build_msnet_peptide("PEPTIDE", &[1, 2], &[1.0], 0.0, 0.0);
        assert!(r.is_err());
    }

    #[test]
    fn nonstandard_residue_errors() {
        let r = build_msnet_peptide("PEPTBDE", &[], &[], 0.0, 0.0);
        assert!(r.is_err());
    }
}
