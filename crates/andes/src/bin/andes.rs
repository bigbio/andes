//! andes: end-to-end peptide-spectrum database search.
//!
//! Loads an MGF or mzML spectrum file and a FASTA target database, runs a
//! tryptic database search and writes output
//! in Percolator `.pin` format (and optionally `.tsv` format).
//!
//! Format dispatch by `--spectrum` extension: `.mzML`/`.mzml` → `MzMLReader`;
//! `.d` → `TimsTofReader` (native Bruker timsTOF, only under `--features
//! timstof`); otherwise `MgfReader` (default reader).

use std::process::ExitCode;

use clap::{Parser, Subcommand};

#[path = "andes/cli.rs"]
mod cli;
#[path = "../config.rs"]
mod config;
#[path = "andes/glyco_run.rs"]
mod glyco_run;
#[path = "andes/model_select.rs"]
mod model_select;
#[path = "../rescore.rs"]
mod rescore;
#[path = "andes/search.rs"]
mod search;
#[path = "andes/spectra.rs"]
mod spectra;
#[path = "andes/train.rs"]
mod train;
#[path = "andes/train_intensity.rs"]
mod train_intensity;

// Re-imported so `config.rs` can keep its `crate::…` paths unchanged.
use crate::cli::{
    parse_charge_range, parse_enzyme_specificity, parse_fragmentation, parse_isotope_error_range,
    parse_precursor_cal, parse_precursor_tol, parse_protocol, CandidateIndexFlag, Cli, ScoreFlag,
    SearchArgs,
};
use crate::search::run;
use crate::train::{
    run_rescore_pin, run_train, run_train_from_search, RescorePinArgs, TrainArgs,
    TrainFromSearchArgs,
};
use crate::train_intensity::{
    run_train_intensity, run_train_intensity_gbdt, run_train_rich_ion_llr, TrainIntensityArgs,
    TrainIntensityGbdtArgs, TrainRichIonLlrArgs,
};

/// Available system memory in bytes, from Linux `/proc/meminfo` (`MemAvailable`).
/// Returns `None` when it cannot be determined (non-Linux, restricted sandbox);
/// callers should then keep the conservative RAM default.
fn available_memory_bytes() -> Option<u64> {
    let meminfo = std::fs::read_to_string("/proc/meminfo").ok()?;
    for line in meminfo.lines() {
        if let Some(rest) = line.strip_prefix("MemAvailable:") {
            // e.g. "MemAvailable:   12345678 kB"
            let kb: u64 = rest.split_whitespace().next()?.parse().ok()?;
            return Some(kb.saturating_mul(1024));
        }
    }
    None
}

/// Emit a one-line search-progress update after each scored chunk: cumulative
/// spectra scored, throughput, and elapsed time. The streaming search has no
/// upfront total for mzML/MGF, so this reports a running count + rate rather than
/// a percentage (use it to see the search is live and how fast it is going).
fn report_search_progress(scored: usize, start: std::time::Instant) {
    let secs = start.elapsed().as_secs_f64();
    let rate = if secs > 0.0 {
        scored as f64 / secs
    } else {
        0.0
    };
    eprintln!("[search] {scored} spectra scored (~{rate:.0}/s, {secs:.0}s elapsed)");
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Rescore an existing PIN with the built-in rescorer and write Percolator-shaped
    /// `.psms` / `.dpsms` files.
    ///
    /// Exists so two rescorers can be compared on IDENTICAL input: without it the only
    /// way to try the native rescorer was to re-run the whole search, which changes the
    /// PIN as well as the rescoring and confounds the comparison.
    #[command(name = "rescore-pin", hide = true)]
    RescorePin(RescorePinArgs),

    /// Train a scoring model directly from externally-labeled, high-confidence
    /// PSMs supplied as flat training parquet(s), bypassing the bootstrap
    /// search. This is the primary training path for the Phase-3 "own models".
    ///
    /// Boxed to keep the `Command` enum compact (clippy `large_enum_variant`).
    #[command(hide = true)]
    Train(Box<TrainArgs>),

    /// Train a scoring model from spectra and a FASTA database, writing the
    /// result to a Parquet model store.
    ///
    /// Boxed to keep the `Command` enum compact.
    #[command(name = "train-from-search", hide = true)]
    TrainFromSearch(Box<TrainFromSearchArgs>),

    /// Merge MSNet intensity aggregation parquets into a finalized intensity
    /// model for the strong-score numerator.
    #[command(name = "train-intensity", hide = true)]
    TrainIntensity(Box<TrainIntensityArgs>),

    /// Train a v3 GBDT fragment-intensity regressor from flat training
    /// parquets and embed it in a Parquet model store alongside existing
    /// rank-core models.
    ///
    /// Boxed to keep the `Command` enum compact.
    #[command(name = "train-intensity-gbdt", hide = true)]
    TrainIntensityGbdt(Box<TrainIntensityGbdtArgs>),

    /// Train a GBDT rich-ion LLR classifier (logistic; decoy-aware) from flat
    /// training parquets and embed it in a Parquet model store alongside
    /// existing rank-core models.
    ///
    /// Boxed to keep the `Command` enum compact.
    #[command(name = "train-rich-ion-llr", hide = true)]
    TrainRichIonLlr(Box<TrainRichIonLlrArgs>),
}

/// Top-level CLI.  When no subcommand is given, the flattened `SearchArgs`
/// drive the existing search path (byte-identical to the pre-subcommand
/// behaviour).
#[derive(Parser, Debug)]
#[command(
    name = "andes",
    // `--version` / `-V` from the crate version. Every CLI is expected to answer it,
    // and a search result is only reproducible if you can record which build produced
    // it — so this is provenance, not just convention.
    version = env!("CARGO_PKG_VERSION"),
    about = "andes: database search of MGF/mzML spectra against FASTA",
    allow_hyphen_values = true,
)]
struct TopCli {
    #[command(subcommand)]
    command: Option<Command>,

    #[command(flatten)]
    search: SearchArgs,
}

fn main() -> ExitCode {
    #[cfg(feature = "thermo")]
    configure_bundled_dotnet();
    // Parse via get_matches so we can query each flag's ValueSource (for the
    // --config merge: an explicit CLI flag must override the YAML value).
    let matches = <TopCli as clap::CommandFactory>::command().get_matches();
    let mut top =
        <TopCli as clap::FromArgMatches>::from_arg_matches(&matches).unwrap_or_else(|e| e.exit());
    // Record whether the user typed --max-missed-cleavages. `--glyco` raises the
    // floor to 3, which costs real memory (measured +4.4 GB on a 20k-protein human
    // FASTA), so an explicit lower value must be allowed to win. The default (1) is
    // itself a plausible explicit value, so comparing against it cannot distinguish
    // the two cases — only the ValueSource can.
    let _ = EXPLICIT_MISSED_CLEAVAGES.set(matches!(
        matches.value_source("max_missed_cleavages"),
        Some(clap::parser::ValueSource::CommandLine)
    ));
    let result = match top.command.take() {
        Some(Command::RescorePin(args)) => run_rescore_pin(args),
        Some(Command::Train(args)) => run_train(*args),
        Some(Command::TrainFromSearch(args)) => run_train_from_search(*args),
        Some(Command::TrainIntensity(args)) => run_train_intensity(*args),
        Some(Command::TrainIntensityGbdt(args)) => run_train_intensity_gbdt(*args),
        Some(Command::TrainRichIonLlr(args)) => run_train_rich_ion_llr(*args),
        None => {
            // --config: fill any parameter the user did not type on the CLI.
            if let Some(cfg_path) = top.search.config.clone() {
                if let Err(e) = config::RunConfig::load(&cfg_path)
                    .and_then(|c| config::apply(c, &mut top.search, &matches))
                {
                    eprintln!("error: {e}");
                    return ExitCode::from(2);
                }
            }
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
            // --output-pin is required UNLESS --rescore is set: rescore can route
            // the search through a temporary PIN (deleted afterwards when
            // --keep-pin false), so a PIN path is not mandatory in that mode.
            if search.output_pin.is_none() && !search.rescore && !search.rescore_native {
                eprintln!("error: --output-pin is required for search (or use --rescore)");
                return ExitCode::from(2);
            }
            // --fdr / --pep are thresholds APPLIED BY a rescoring run; they do not
            // start one on their own. Warn (don't silently launch a backend) if set
            // without an explicit --rescore / --rescore-native.
            if (search.fdr.is_some() || search.pep.is_some())
                && !search.rescore
                && !search.rescore_native
            {
                eprintln!(
                    "warning: --fdr/--pep set without --rescore or --rescore-native; \
                     they are ignored (add --rescore to run Percolator, or feed the PIN \
                     to Percolator yourself)."
                );
            }
            run(Cli {
                database: Some(database),
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
            if bundled
                .join("shared")
                .join("Microsoft.NETCore.App")
                .is_dir()
            {
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
    let probe_set = RSS_PROBE.get().copied().unwrap_or(false);
    if !probe_set {
        return;
    }
    #[cfg(target_os = "linux")]
    {
        if let Ok(s) = std::fs::read_to_string("/proc/self/status") {
            for line in s.lines() {
                if line.starts_with("VmRSS:") {
                    eprintln!("[RSS {tag}] {}", line.trim_start_matches("VmRSS:").trim());
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

/// True when `flag` appears literally in the process arguments. Used to tell a
/// user-supplied value from a clap default, where the default itself is
/// context-dependent (see `--gbdt-max-trees`, which is exempt under `--glyco`).
fn arg_present(flag: &str) -> bool {
    std::env::args().any(|a| a == flag || a.starts_with(&format!("{flag}=")))
}

/// Diagnostic RSS logging, installed from `--rss-probe`.
static RSS_PROBE: std::sync::OnceLock<bool> = std::sync::OnceLock::new();

/// Set in `main` from clap's `ValueSource`: did the user type
/// `--max-missed-cleavages` themselves? See the glyco floor below.
static EXPLICIT_MISSED_CLEAVAGES: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
