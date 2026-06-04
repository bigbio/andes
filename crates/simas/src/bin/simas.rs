//! simas: end-to-end peptide-spectrum database search.
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

/// Instrument class. Drives the `LowRes`/`HighRes`/`TOF`/`QExactive`
/// classification used to pick the bundled param file.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum Instrument {
    #[clap(name = "low-res")]   LowRes,
    #[clap(name = "high-res")]  HighRes,
    #[clap(name = "TOF")]       Tof,
    #[clap(name = "QExactive")] QExactive,
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

/// Search arguments (shared by the default search path and exposed as a
/// flat arg group so that `simas --spectrum X --database Y --output-pin Z`
/// keeps working unchanged).
///
/// Note: `spectrum`, `database`, and `output_pin` are declared `Option<PathBuf>`
/// at the clap level so that they are not required when a subcommand (e.g.
/// `train`) is given.  When no subcommand is present, `run()` validates them
/// manually and returns an early error if they are missing.
#[derive(Args, Debug)]
struct SearchArgs {
    /// Input spectrum file. Format is auto-detected by extension:
    /// `.mzML`/`.mzml` is read as mzML, anything else as MGF.
    #[arg(long)]
    spectrum: Option<PathBuf>,

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
    #[arg(long = "precursor-cal", default_value = "off", value_parser = parse_precursor_cal)]
    precursor_cal: PrecursorCalMode,

    /// Precursor mass tolerance in ppm.
    #[arg(long, default_value = "20.0")]
    precursor_tol_ppm: f64,

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
    /// equivalent to Java -ntt 2). `semi`: at least one terminus must be a
    /// cleavage site (Java -ntt 1). `non-specific`: neither terminus needs
    /// to be a cleavage site (Java -ntt 0). Legacy numeric 0/1/2 still accepted.
    #[arg(long = "enzyme-specificity", alias = "ntt",
          default_value = "fully", value_parser = parse_enzyme_specificity)]
    enzyme_specificity: EnzymeSpecificity,

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

    /// Maximum peptide length, in residues.
    #[arg(long, default_value = "40")]
    max_length: u32,

    /// Path to the .param scoring model file.
    ///
    /// If not supplied, a bundled file under
    /// `resources/ionstat/` is selected from
    /// `(--fragmentation, --instrument, --protocol)` (default
    /// `HCD_QExactive_Tryp.param`). When running the binary outside the source
    /// tree this path may not exist; supply --param-file explicitly in that
    /// case.
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

    /// Fragmentation method. Named values: auto, CID, ETD, HCD, UVPD.
    /// Legacy numeric (Java MS-GF+ `-m`): 0=auto, 1=CID, 2=ETD, 3=HCD, 4=UVPD.
    #[arg(long, default_value = "auto", value_parser = parse_fragmentation)]
    fragmentation: Fragmentation,

    /// Instrument class. Named values: low-res, high-res, TOF, QExactive.
    /// Legacy numeric (Java MS-GF+ `-inst`): 0=low-res, 1=high-res, 2=TOF, 3=QExactive.
    #[arg(long, default_value = "low-res", value_parser = parse_instrument)]
    instrument: Instrument,

    /// Search protocol. Named values: auto, phospho, iTRAQ, iTRAQ-phospho, TMT, standard.
    /// Legacy numeric (Java MS-GF+ `-protocol`): 0=auto, 1=phospho, 2=iTRAQ, 3=iTRAQ-phospho, 4=TMT, 5=standard.
    #[arg(long, default_value = "auto", value_parser = parse_protocol)]
    protocol: Protocol,

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
    /// When set, skips automatic selection by `(--fragmentation, --instrument,
    /// --protocol)` and loads this ID directly. Useful after `simas train`
    /// to search with the freshly-trained model.
    #[arg(long = "model")]
    model_id_override: Option<String>,
}

/// Training arguments for `simas train`.
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

/// Available subcommands.
#[derive(Subcommand, Debug)]
enum Command {
    /// Train a scoring model from spectra and a FASTA database, writing the
    /// result to a Parquet model store.
    Train(TrainArgs),
}

/// Top-level CLI.  When no subcommand is given, the flattened `SearchArgs`
/// drive the existing search path (byte-identical to the pre-subcommand
/// behaviour).
#[derive(Parser, Debug)]
#[command(
    name = "simas",
    about = "simas: database search of MGF/mzML spectra against FASTA",
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
        Some(Command::Train(args)) => run_train(args),
        None => {
            // Validate required search args that are Option<> at the clap level.
            let search = top.search;
            let spectrum = match search.spectrum {
                Some(p) => p,
                None => {
                    eprintln!("error: --spectrum is required for search (or use `simas train`)");
                    return ExitCode::from(2);
                }
            };
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
            // Reconstruct a Cli (= SearchArgs) with the validated paths.
            run(Cli {
                spectrum: Some(spectrum),
                database: Some(database),
                output_pin: Some(output_pin),
                ..search
            })
        }
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("simas: {e}");
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

/// Print VmRSS for the current process under MSGF_RSS_PROBE=1. No-op
/// otherwise and a no-op on non-Linux platforms regardless of the env var.
/// (Legacy name MSGFRUST_RSS_PROBE is accepted with a deprecation warning.)
///
/// We gate behind an env var so production runs stay quiet; flip the var on
/// when debugging memory regressions.
fn log_rss(tag: &str) {
    // Accept both new and legacy env var names. Legacy emits the
    // deprecation warning once per process (sync::Once guard).
    let new_set = std::env::var_os("MSGF_RSS_PROBE").is_some();
    let legacy_set = std::env::var_os("MSGFRUST_RSS_PROBE").is_some();
    if legacy_set && !new_set {
        static LEGACY_WARN_ONCE: std::sync::Once = std::sync::Once::new();
        LEGACY_WARN_ONCE.call_once(|| {
            eprintln!(
                "WARN: MSGFRUST_RSS_PROBE is deprecated; use MSGF_RSS_PROBE \
                 (legacy name accepted in this release, will be removed next)"
            );
        });
    }
    if !new_set && !legacy_set {
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
            "Precursor mass calibration skipped (insufficient confident PSMs: {} with PSMs, {} failed SpecE, {} failed |residual|>50ppm; elapsed: {:.2}s)",
            stats.queues_with_psm,
            stats.rejected_spec_e,
            stats.rejected_residual,
            t_cal.elapsed().as_secs_f64()
        );
    }
    Ok(stats)
}

fn run(cli: Cli) -> Result<(), Box<dyn std::error::Error>> {
    // These three were validated as Some(..) by main() before calling run().
    let spectrum_path: PathBuf = cli.spectrum.expect("spectrum validated in main");
    let database_path: PathBuf = cli.database.expect("database validated in main");
    let output_pin_path: PathBuf = cli.output_pin.expect("output_pin validated in main");

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
    // If --mod is given, parse the Java-format mods.txt file. Otherwise
    // fall back to simas's historical defaults (CAM fixed on C,
    // Oxidation variable on M) so existing tests keep their behaviour.
    //
    // `num_mods_from_file` is populated only when --mod is given and the
    // file contains a `NumMods=N` line; it overrides the default
    // `max_variable_mods_per_peptide` (3) below.
    let (aa, num_mods_from_file) = match &cli.mods {
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
        None => {
            let cam = Modification {
                name: "Carbamidomethyl".into(),
                mass_delta: 57.02146,
                residue: ResidueSpec::Specific(b'C'),
                location: ModLocation::Anywhere,
                fixed: true,
                accession: None,
            };
            let ox = Modification {
                name: "Oxidation".into(),
                mass_delta: 15.99491,
                residue: ResidueSpec::Specific(b'M'),
                location: ModLocation::Anywhere,
                fixed: false,
                accession: None,
            };
            let set = AminoAcidSetBuilder::new_standard()
                .add_fixed_mod(cam)
                .add_variable_mod(ox)
                .build()?;
            (set, None)
        }
    };

    // ── 4. Load Param scoring model ───────────────────────────────────────────
    //
    // `--param-file` wins outright. Otherwise, for mzML with `--fragmentation auto`,
    // peek the file's dominant activation method and pick the bundled `.param`.
    // MGF and explicit fragmentation/instrument flags use `resolve_bundled_param`.
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
    // is always CID/TimsTOF (DDA-PASEF). Detection only runs when
    // `--fragmentation auto` is set (otherwise the CLI flags override).
    let auto_route_eligible = cli.fragmentation == Fragmentation::Auto && (is_mzml || is_raw || is_d);
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

    let t_phase = std::time::Instant::now();
    let param = if let Some(ref override_path) = cli.param_file {
        // ── Override path: load binary .param directly (unchanged behaviour). ──
        eprintln!("Param file (override): {}", override_path.display());
        Param::load_from_file(override_path)
            .map_err(|e| format!("loading param file {}: {e}", override_path.display()))?
    } else {
        // ── Auto / explicit flags: resolve from the Parquet model store. ──────
        //
        // For `--fragmentation auto` with a detectable input the detected
        // (activation, instrument) is used; for explicit flags or MGF (no
        // detection) the CLI enum values are converted to ActivationMethod /
        // InstrumentType for the store lookup.
        let (activation, instrument_opt): (ActivationMethod, Option<InstrumentType>) =
            if auto_route_eligible {
                match detected_activation_instrument {
                    Some((method, inst)) => {
                        eprintln!(
                            "Param resolver: auto-detected dominant activation \
                             method = {} (instrument = {}) from {}",
                            method.name(),
                            inst.map(|i| i.name()).unwrap_or("unknown/default"),
                            spectrum_path.display()
                        );
                        (method, inst)
                    }
                    None => {
                        // No detectable activation — fall back to CLI flags.
                        // For the all-defaults case (Auto+LowRes+Auto) this
                        // returns HCD/QExactive to match the historical default.
                        cli_flags_to_activation_instrument(
                            cli.fragmentation, cli.instrument, cli.protocol,
                        )
                    }
                }
            } else {
                // Explicit `--fragmentation` / `--instrument` flags (or MGF
                // where auto-detection is not eligible).
                cli_flags_to_activation_instrument(
                    cli.fragmentation, cli.instrument, cli.protocol,
                )
            };

        let (model_id, p) = load_param_from_store(
            activation,
            instrument_opt,
            cli.protocol,
            cli.model_store.as_deref(),
            cli.model_id_override.as_deref(),
        )?;
        eprintln!("Param model: {model_id} (from store)");
        p
    };
    let scorer = RankScorer::new(&param);
    eprintln!("[PHASE param_and_scorer: {:.2}s]", t_phase.elapsed().as_secs_f64());

    // ── 5. Build SearchParams ─────────────────────────────────────────────────
    let mut params = SearchParams::default_tryptic(aa);
    params.precursor_tolerance =
        PrecursorTolerance::symmetric(Tolerance::Ppm(cli.precursor_tol_ppm));
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
    params.num_tolerable_termini = match cli.enzyme_specificity {
        EnzymeSpecificity::Fully => 2,
        EnzymeSpecificity::Semi => 1,
        EnzymeSpecificity::NonSpecific => 0,
    };
    params.max_missed_cleavages = cli.max_missed_cleavages;
    params.min_peaks = cli.min_peaks;
    params.min_length = cli.min_length;
    params.max_length = cli.max_length;
    if let Some(n) = num_mods_from_file {
        params.max_variable_mods_per_peptide = n;
    }
    params.precursor_cal_mode = cli.precursor_cal;
    params.precursor_mass_shift_ppm = 0.0;

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

    let prepared = match reuse_parts {
        Some(parts) => {
            PreparedSearch::from_parts(&idx, &params, &scorer, fragment_tol_da, parts)
        }
        None => PreparedSearch::prepare(&idx, &params, &scorer, fragment_tol_da, &cli.decoy_prefix),
    };
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
    let parse_stats = if chimeric_input {
        let (tx, rx) = sync_channel::<(Vec<Spectrum>, Ms1Link)>(2);
        let spectrum_path = spectrum_path.clone();
        let cap = bench_cap;
        // The cascade is MS2-only by construction (MS2 paired with its preceding
        // MS1); hardcode MS2 so `--ms-level 3` can never widen the mzML reader's
        // range to admit MS3 (the .raw chunked reader is already MS2-only).
        let mslevel = 2;
        let parser_handle = thread::spawn(move || -> Result<(usize, Vec<String>), String> {
            // The inner closures borrow `tx` (SyncSender::send takes &self), so
            // both reader branches can use it without moving it twice.
            if is_mzml {
                let f = File::open(&spectrum_path).map_err(|e| format!("open mzML: {e}"))?;
                let reader = MzMLReader::new(BufReader::new(f))
                    .with_ms_level_range(mslevel, mslevel)
                    .with_ms1_capture(true);
                let (errc, errs) = reader.read_with_ms1_chunked(CHUNK_SIZE, cap, |chunk, link| {
                    // If the consumer hung up, sending fails; nothing more to do.
                    let _ = tx.send((chunk, link));
                });
                Ok((errc, errs))
            } else {
                // is_raw — same MS2 + bounded-MS1 chunk stream from the .raw.
                #[cfg(feature = "thermo")]
                {
                    let reader = input::ThermoRawReader::open(&spectrum_path)
                        .map_err(|e| format!("open Thermo .raw: {e}"))?;
                    let (errc, errs) = reader.read_with_ms1_chunked(CHUNK_SIZE, cap, |chunk, link| {
                        let _ = tx.send((chunk, link));
                    });
                    Ok((errc, errs))
                }
                #[cfg(not(feature = "thermo"))]
                {
                    Err("this simas build has no Thermo .raw support; \
                         rebuild with `--features thermo`."
                        .to_string())
                }
            }
        });

        let mut offset = 0usize;
        let mut ms1_linked = 0usize;
        for (chunk_spectra, chunk_link) in rx {
            // Pass 1 (offset → global spectrum_idx for PIN), then Pass 2 on the same
            // chunk (BEFORE peaks are dropped — the residual needs them). The chunk's
            // own `Ms1Link` and `offset` keep chunk-local indexing aligned.
            let mut queues = prepared.run_chunk(&chunk_spectra, offset);
            search::match_engine::run_pass2_coisolation(
                &prepared,
                &chunk_spectra,
                &mut queues,
                &params,
                &chunk_link,
                offset,
            );
            offset += chunk_spectra.len();
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
            "chimeric mode: streamed {} MS2 spectra ({} MS1 scans linked) in chunks of {}",
            offset, ms1_linked, CHUNK_SIZE
        );
        log_rss("after_chimeric_stream_search");
        ParseStats { error_count: err_count, first_errors }
    } else {
        // ── Non-chimeric streaming pipeline (UNCHANGED) ───────────────────
        let (tx, rx) = sync_channel::<Vec<Spectrum>>(2);

        // Spawn the parser thread. It owns the reader (paths + flags moved in).
        // The thread returns ParseStats with the error count + sample messages.
        let spectrum_path = spectrum_path.clone();
        let parser_handle = thread::spawn(move || -> Result<ParseStats, Box<dyn std::error::Error + Send + Sync>> {
            if is_mzml {
                let f = File::open(&spectrum_path)
                    .map_err(|e| format!("open mzML: {e}"))?;
                let reader = MzMLReader::new(BufReader::new(f))
                    .with_ms_level_range(ms_level_u32, ms_level_u32);
                Ok(send_chunks(reader, CHUNK_SIZE, bench_cap, tx))
            } else if is_raw {
                // Native Thermo .raw (feature-gated). The reader spawns no
                // thread of its own; it yields the same Spectrum stream as the
                // mzML/MGF readers. MS2-only: identification scans, never MS3
                // (e.g. TMT SPS-MS3 reporter scans) regardless of `--ms-level`.
                #[cfg(feature = "thermo")]
                {
                    let reader = input::ThermoRawReader::open(&spectrum_path)
                        .map_err(|e| format!("open Thermo .raw: {e}"))?
                        .with_ms_level(Some(2));
                    Ok(send_chunks(reader, CHUNK_SIZE, bench_cap, tx))
                }
                #[cfg(not(feature = "thermo"))]
                {
                    Err("this simas build has no Thermo .raw support; \
                         rebuild with `--features thermo` (and run with the \
                         .NET 8 runtime installed). mzML/MGF inputs work without it."
                        .into())
                }
            } else if is_d {
                // Native Bruker timsTOF `.d` (feature-gated). The reader opens
                // the `.d` directory and yields the same MS2 `Spectrum` stream
                // as the mzML/MGF readers; it spawns no thread of its own.
                #[cfg(feature = "timstof")]
                {
                    let reader = input::TimsTofReader::open(&spectrum_path)
                        .map_err(|e| format!("open Bruker .d: {e}"))?;
                    Ok(send_chunks(reader, CHUNK_SIZE, bench_cap, tx))
                }
                #[cfg(not(feature = "timstof"))]
                {
                    Err("this simas build has no Bruker .d (timsTOF) support; \
                         rebuild with `--features timstof`. mzML/MGF inputs work \
                         without it."
                        .into())
                }
            } else {
                let f = File::open(&spectrum_path)
                    .map_err(|e| format!("open MGF: {e}"))?;
                let reader = MgfReader::new(BufReader::new(f));
                Ok(send_chunks(reader, CHUNK_SIZE, bench_cap, tx))
            }
        });

        log_rss("after_parser_thread_spawn");

        // Consumer loop: drain chunks from the channel as they arrive. Each
        // received chunk is processed via `prepared.run_chunk` (Rayon-parallel)
        // synchronously on this thread; while the inner Rayon runs, the parser
        // thread is filling the next chunk concurrently.
        for chunk in rx {
            if chunk.is_empty() {
                continue;
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

        // Reap the parser thread for its stats. join() should never block here
        // (channel close has already fired on parser exit).
        match parser_handle.join() {
            Ok(Ok(stats)) => stats,
            Ok(Err(e)) => return Err(format!("parser thread error: {e}").into()),
            Err(_) => return Err("parser thread panicked".into()),
        }
    };

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
        return Err(format!(
            "no spectra parsed from {}",
            spectrum_path.display()
        )
        .into());
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
        let spec_file_name = spectrum_path
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_else(|| spectrum_path.display().to_string());
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
            };
            let ox = Modification {
                name: "Oxidation".into(),
                mass_delta: 15.99491,
                residue: ResidueSpec::Specific(b'M'),
                location: ModLocation::Anywhere,
                fixed: false,
                accession: None,
            };
            Ok(AminoAcidSetBuilder::new_standard()
                .add_fixed_mod(cam)
                .add_variable_mod(ox)
                .build()?)
        }
    }
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

/// Convert the CLI `Fragmentation` enum to `ActivationMethod`.
///
/// `Fragmentation::Auto` is treated as `ActivationMethod::CID` here because
/// the old ladder normalised `Auto → CID` when explicit flags were used
/// (without mzML peek). When `--fragmentation auto` is combined with mzML/raw
/// input the detection path is taken instead of this function.
fn cli_fragmentation_to_activation(f: Fragmentation) -> ActivationMethod {
    match f {
        Fragmentation::Auto => ActivationMethod::CID,
        Fragmentation::Cid  => ActivationMethod::CID,
        Fragmentation::Etd  => ActivationMethod::ETD,
        Fragmentation::Hcd  => ActivationMethod::HCD,
        Fragmentation::Uvpd => ActivationMethod::UVPD,
    }
}

/// Convert the CLI `Instrument` enum to `InstrumentType`.
fn cli_instrument_to_instrument_type(i: Instrument) -> InstrumentType {
    match i {
        Instrument::LowRes    => InstrumentType::LowRes,
        Instrument::HighRes   => InstrumentType::HighRes,
        Instrument::Tof       => InstrumentType::TOF,
        Instrument::QExactive => InstrumentType::QExactive,
    }
}

/// Resolve `(Fragmentation, Instrument, Protocol)` from CLI flags to
/// `(ActivationMethod, InstrumentType, Protocol)` for store lookup.
///
/// Handles the historical all-defaults short-circuit: when the user omits
/// all scoring-model flags (`--fragmentation auto`, `--instrument low-res`,
/// `--protocol auto`) the old ladder returned `HCD_QExactive_Tryp.param`.
/// We replicate this by returning `(HCD, QExactive, Auto)` for that case
/// instead of `(CID, LowRes, Auto)` (which would resolve to `cid_lowres_tryp`).
fn cli_flags_to_activation_instrument(
    fragmentation: Fragmentation,
    instrument: Instrument,
    protocol: Protocol,
) -> (ActivationMethod, Option<InstrumentType>) {
    // Historical all-defaults short-circuit (mirrors resolve_bundled_param step 0).
    if fragmentation == Fragmentation::Auto
        && instrument == Instrument::LowRes
        && protocol == Protocol::Auto
    {
        return (ActivationMethod::HCD, Some(InstrumentType::QExactive));
    }
    (
        cli_fragmentation_to_activation(fragmentation),
        Some(cli_instrument_to_instrument_type(instrument)),
    )
}

/// Translate `(--fragmentation, --instrument, --protocol)` into a bundled
/// `.param` filename and resolve it under
/// `resources/ionstat/` relative to the cargo manifest dir.
///
/// CLI indices match Java's:
/// - fragmentation: 0=Auto/CID, 1=CID, 2=ETD, 3=HCD, 4=UVPD
/// - instrument:    0=LowRes,   1=HighRes, 2=TOF, 3=QExactive
/// - protocol:      0=Automatic,1=Phosphorylation, 2=iTRAQ,
///   3=iTRAQPhospho, 4=TMT, 5=Standard
///
/// When all three are `None`, the historical default
/// `HCD_QExactive_Tryp.param` is returned (preserving existing tests'
/// behaviour). Only Tryp is supported as the enzyme component for now;
/// other enzymes require the user to pass `--param-file` directly.
///
/// Walks Java's `NewScorerFactory.get(...)` fallback ladder: try the exact
/// `{frag}_{inst}_Tryp{protocol}.param` first; if that doesn't resolve, drop
/// the protocol suffix; if that also doesn't resolve, use the final
/// `(frag, inst)`-keyed ladder. Returns an error only if even the
/// last-resort `CID_LowRes_Tryp.param` is missing from the bundled
/// resources (a packaging defect, not a CLI input error).
///
/// Reference implementation of the historical filename-based resolution ladder.
/// The search path now goes through [`load_param_from_store`]; the store-selection
/// equivalence test validates the store-based selection against an independent
/// copy of this logic.
#[allow(dead_code)]
fn resolve_bundled_param(
    fragmentation: Fragmentation,
    instrument:    Instrument,
    protocol:      Protocol,
) -> Result<PathBuf, String> {
    // Step 0: default-to-bundled short-circuit. When the caller passes all
    // defaults (Fragmentation::Auto, Instrument::LowRes, Protocol::Auto),
    // fall back to the bundled HCD_QExactive_Tryp.param — the behavior of
    // omitting all three flags.
    if fragmentation == Fragmentation::Auto
        && instrument == Instrument::LowRes
        && protocol == Protocol::Auto {
        return canonicalize_bundled("HCD_QExactive_Tryp.param");
    }

    // Step 1: Normalize. Java's normalization rules mirrored here:
    //   - Auto fragmentation → CID (Java's "null/PQD → CID")
    //   - HCD with low-res inst → upgrade to QExactive (Java's HCD-upgrade rule)
    let frag = match fragmentation {
        Fragmentation::Auto => "CID",
        Fragmentation::Cid  => "CID",
        Fragmentation::Etd  => "ETD",
        Fragmentation::Hcd  => "HCD",
        Fragmentation::Uvpd => "UVPD",
    };
    let mut inst = match instrument {
        Instrument::LowRes    => "LowRes",
        Instrument::HighRes   => "HighRes",
        Instrument::Tof       => "TOF",
        Instrument::QExactive => "QExactive",
    };
    // HCD-upgrade rule: HCD with low-res inst → upgrade to QExactive.
    if frag == "HCD" && inst == "LowRes" {
        inst = "QExactive";
    }

    let prot_suffix: &str = match protocol {
        Protocol::Auto         => "",          // empty: no protocol suffix
        Protocol::Phospho      => "_Phosphorylation",
        Protocol::Itraq        => "_iTRAQ",
        Protocol::ItraqPhospho => "_iTRAQPhospho",
        Protocol::Tmt          => "_TMT",
        Protocol::Standard     => "",          // standard = no suffix
    };

    // Step 1: Try the exact requested combination first.
    //   `{frag}_{inst}_Tryp{prot_suffix}.param`
    let exact = format!("{frag}_{inst}_Tryp{prot_suffix}.param");
    if let Ok(path) = canonicalize_bundled(&exact) {
        return Ok(path);
    }

    // Step 2: Drop protocol — try `{frag}_{inst}_Tryp.param`.
    // This mirrors Java parity: `return get(method, instType, enzyme)` fallback
    // (drop protocol suffix when exact match is missing). For (CID, HighRes, Tryp, TMT) this
    // lands on `CID_HighRes_Tryp.param`, which IS what Java would pick when
    // the protocol-specific file is missing.
    if !prot_suffix.is_empty() {
        let no_protocol = format!("{frag}_{inst}_Tryp.param");
        if let Ok(path) = canonicalize_bundled(&no_protocol) {
            eprintln!(
                "Param resolver: `{exact}` not bundled; falling back to `{no_protocol}` \
                 (Java NewScorerFactory drops protocol suffix when exact match missing)",
            );
            return Ok(path);
        }
    }

    // Step 3: Alternate enzyme — Java tries Trypsin (for C-term enzymes) or
    // LysN (for N-term enzymes). We always use Tryp here, so this step is
    // a no-op for now. If/when N-term enzyme support lands, replicate this.

    // Step 4: Final fallback ladder (Java parity for scorer factory fallback).
    //   - HCD + (TOF|HighRes) + C-term → CID_TOF_Tryp
    //   - ETD + C-term                  → ETD_LowRes_Tryp
    //   - Non-electron + N-term         → CID_LowRes_LysN  (skipped; N-term TBD)
    //   - default                        → CID_LowRes_Tryp
    //
    // For our currently-supported (frag, inst) combos:
    let final_fallback = match (frag, inst) {
        ("HCD", "TOF") | ("HCD", "HighRes") => "CID_TOF_Tryp.param",
        ("ETD", _) => "ETD_LowRes_Tryp.param",
        _ => "CID_LowRes_Tryp.param",
    };
    eprintln!(
        "Param resolver: `{exact}` not bundled and protocol-less drop also missing; \
         using final fallback `{final_fallback}` (Java NewScorerFactory final ladder)",
    );
    canonicalize_bundled(final_fallback)
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

/// Resolve a bundled `.param` file for the given activation method.
///
/// This is the auto-detect path: we already know the activation, and we
/// pick the bundled instrument+enzyme pair that best matches the dataset.
/// Mirrors Java parity for per-spectrum param dispatch when the user passes `-m 0`
/// (ASWRITTEN), but applied at file-wide granularity here.
///
/// The `detected_instrument` argument is the instrument type detected by
/// scanning the mzML's `<instrumentConfiguration>` blocks (see
/// `input::detect_instrument_type`). `None` means we couldn't detect it
/// (older mzML, MGF, etc.) — in that case we mirror Java's
/// `NewScorerFactory.get` default of `LOW_RESOLUTION_LTQ`.
///
/// Mapping (Tryp / no-protocol unless protocol overrides):
///   - CID  → frag=1, inst=detected (LowRes when none).
///     LowRes for LTQ Velos / ion-trap data; HighRes / QExactive
///     for Orbitrap data. Matches Java's default + the user-supplied
///     `-inst` path.
///   - HCD  → frag=3, inst=detected. `resolve_bundled_param`'s Java-mirror
///     normalization upgrades HCD with non-(HighRes|QExactive) to
///     QExactive, so HCD on LTQ data still routes to a QExactive
///     model (Java does the same).
///   - ETD  → frag=2, inst=detected.
///   - PQD  → CID (Java collapses PQD → CID in `NewScorerFactory.get`).
///   - UVPD → frag=4, inst=QExactive (only QExactive variant exists bundled).
///
/// Reference implementation of the historical activation-aware resolution ladder
/// (kept alongside [`resolve_bundled_param`]); the search path now uses
/// [`load_param_from_store`].
#[allow(dead_code)]
fn resolve_bundled_param_for_activation(
    method:               ActivationMethod,
    detected_instrument:  Option<InstrumentType>,
    protocol:             Protocol,
) -> Result<PathBuf, String> {
    let frag = match method {
        ActivationMethod::CID  => Fragmentation::Cid,
        ActivationMethod::ETD  => Fragmentation::Etd,
        ActivationMethod::HCD  => Fragmentation::Hcd,
        ActivationMethod::UVPD => Fragmentation::Uvpd,
        // PQD → CID (Java's NewScorerFactory rule: "PQD or null → CID").
        ActivationMethod::PQD  => Fragmentation::Cid,
    };
    // New instrument classes fall back to their family for model resolution
    // (no Astral-specific or TimsTOF-specific .param file bundled yet).
    let inst = match detected_instrument.map(|i| i.family_fallback()) {
        Some(InstrumentType::LowRes)    => Instrument::LowRes,
        Some(InstrumentType::HighRes)   => Instrument::HighRes,
        Some(InstrumentType::TOF)       => Instrument::Tof,
        Some(InstrumentType::QExactive) => Instrument::QExactive,
        // OrbitrapAstral → QExactive and TimsTOF → TOF via family_fallback above;
        // these arms are unreachable after fallback but keep the match exhaustive.
        Some(InstrumentType::OrbitrapAstral) => Instrument::QExactive,
        Some(InstrumentType::TimsTOF)        => Instrument::Tof,
        None                                 => Instrument::LowRes,
    };
    resolve_bundled_param(frag, inst, protocol)
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

/// Resolve a bundled `.param` filename under
/// `resources/ionstat/` relative to the crate's cargo manifest
/// dir (set at compile time). Returns a helpful error if the file does
/// not exist.
#[allow(dead_code)]
fn canonicalize_bundled(filename: &str) -> Result<PathBuf, String> {
    let candidate = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join("resources/ionstat")
        .join(filename);
    candidate.canonicalize().map_err(|e| format!(
        "bundled param file not found at `{}`: {e}\n\
         Hint: not every (fragmentation, instrument, protocol) combination \
         has a bundled .param file. Supply --param-file <PATH> to specify \
         the scoring model explicitly, or list available files under \
         `resources/ionstat/`.",
        candidate.display()
    ))
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
/// search binary in place of `resolve_bundled_param*`.
///
/// `activation`: the detected or explicitly set `ActivationMethod`.
/// `instrument`: the detected or explicitly set `InstrumentType` (None = undetected → LowRes).
/// `protocol`:   the CLI `Protocol` value.
fn build_selection_key(
    activation: ActivationMethod,
    instrument: Option<InstrumentType>,
    protocol: Protocol,
) -> SelectionKey {
    use std::collections::BTreeSet;

    // 1. PQD → CID (Java's NewScorerFactory rule).
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
        // Parquet stores enzyme as "Trypsin" for the tryptic models.
        enzyme: "Trypsin".to_string(),
        experiment_class,
    }
}

/// Load the scoring [`Param`] from the bundled Parquet store for the given
/// `(activation, instrument, protocol)` combination.
///
/// This is the new model-resolution path that replaces
/// `Param::load_from_file(resolve_bundled_param*(...))`. The `model_id`
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
        let key = build_selection_key(activation, instrument, protocol);

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
            select(
                &entries,
                &key,
                // `build_selection_key` already applies family fallback + all
                // normalizations, so the family_fn here is the identity.
                |i| i.to_string(),
                Some("hcd_qexactive_tryp"),
            )
            .unwrap_or("hcd_qexactive_tryp")
            .to_string()
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

/// Parse `--instrument` value. Accepts named (low-res, high-res, TOF,
/// QExactive) or legacy numeric (0=LowRes, 1=HighRes, 2=TOF, 3=QExactive).
fn parse_instrument(s: &str) -> Result<Instrument, String> {
    if let Ok(v) = <Instrument as ValueEnum>::from_str(s, true) { return Ok(v); }
    match s.parse::<u8>() {
        Ok(0) => Ok(Instrument::LowRes),
        Ok(1) => Ok(Instrument::HighRes),
        Ok(2) => Ok(Instrument::Tof),
        Ok(3) => Ok(Instrument::QExactive),
        _ => Err(format!(
            "invalid instrument `{s}`: expected low-res|high-res|TOF|QExactive \
             (or legacy 0..=3)"
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
    fn parse_instrument_rejects_out_of_range_numeric() {
        let err = parse_instrument("99").unwrap_err();
        assert!(err.contains("0..=3"), "got: {err}");
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
}
