use std::process::ExitCode;

use clap::{Parser, Subcommand};

#[path = "andes/acquisition.rs"]
mod acquisition;
#[path = "andes/cli.rs"]
mod cli;
#[path = "andes/cli_parse.rs"]
mod cli_parse;
#[path = "../cohort.rs"]
mod cohort;
#[path = "andes/cohort_cmd.rs"]
mod cohort_cmd;
#[path = "../config.rs"]
mod config;
#[path = "../config_mods.rs"]
mod config_mods;
#[path = "andes/database_cmd.rs"]
mod database_cmd;
#[path = "andes/dia.rs"]
mod dia;
#[path = "andes/dia_joint.rs"]
mod dia_joint;
#[path = "andes/dia_llr.rs"]
mod dia_llr;
#[path = "andes/dia_pass2.rs"]
mod dia_pass2;
#[path = "andes/dia_score.rs"]
mod dia_score;
#[path = "andes/dia_xfer.rs"]
mod dia_xfer;
#[path = "andes/dia_xrun.rs"]
mod dia_xrun;
#[path = "andes/glyco_match.rs"]
mod glyco_match;
#[path = "andes/glyco_quant.rs"]
mod glyco_quant;
#[path = "andes/glyco_run.rs"]
mod glyco_run;
#[path = "andes/inference_cmd.rs"]
mod inference_cmd;
#[path = "andes/logging.rs"]
mod logging;
#[path = "andes/memlimit.rs"]
mod memlimit;
#[path = "andes/model_select.rs"]
mod model_select;
#[path = "andes/mono.rs"]
mod mono;
#[path = "andes/predicted.rs"]
mod predicted;
#[path = "andes/quant_cmd.rs"]
mod quant_cmd;
#[path = "andes/quant_finalize.rs"]
mod quant_finalize;
#[path = "andes/quant_mbr.rs"]
mod quant_mbr;
#[path = "andes/quant_run.rs"]
mod quant_run;
#[path = "../sdrf.rs"]
mod sdrf;
#[path = "andes/optimize.rs"]
mod optimize;
#[path = "../param_sources.rs"]
mod param_sources;
#[path = "../params.rs"]
mod params;
#[path = "../params_render.rs"]
mod params_render;
#[path = "andes/search.rs"]
mod search;
#[path = "andes/spectra.rs"]
mod spectra;
#[path = "andes/tmt_mbr.rs"]
mod tmt_mbr;
#[path = "andes/tmt_study.rs"]
mod tmt_study;
#[path = "andes/train.rs"]
mod train;
#[path = "andes/train_intensity.rs"]
mod train_intensity;
#[path = "../transfer_fdr.rs"]
mod transfer_fdr;
#[path = "andes/xl_ambiguity.rs"]
mod xl_ambiguity;
#[path = "andes/xl_fdr.rs"]
mod xl_fdr;
#[path = "andes/xl_run.rs"]
mod xl_run;

use crate::cli::{
    parse_charge_range, parse_isotope_error_range, parse_open_search_delta, parse_precursor_cal,
    parse_precursor_tol,
    parse_uniprot_evidence, Cli, ScoreFlag, SearchArgs,
};
use crate::search::run;
use crate::train::{run_train, run_train_from_search, TrainArgs, TrainFromSearchArgs};
use crate::train_intensity::{
    run_train_intensity_gbdt, run_train_rich_ion_llr, TrainIntensityGbdtArgs, TrainRichIonLlrArgs,
};

/// Time between two progress lines at the default level; every chunk's line is a debug
/// line.
const PROGRESS_INTERVAL: std::time::Duration = std::time::Duration::from_secs(10);

/// When the last default-level progress line finished writing; `None` until the first one.
static LAST_PROGRESS: std::sync::Mutex<Option<std::time::Instant>> = std::sync::Mutex::new(None);
static PROCESS_START: std::sync::OnceLock<std::time::Instant> = std::sync::OnceLock::new();

fn report_search_progress(scored: usize, start: std::time::Instant) {
    let secs = start.elapsed().as_secs_f64();
    let rate = if secs > 0.0 {
        scored as f64 / secs
    } else {
        0.0
    };
    let mut last = LAST_PROGRESS.lock().unwrap_or_else(|e| e.into_inner());
    let since = last.unwrap_or(*PROCESS_START.get_or_init(std::time::Instant::now));
    if since.elapsed() >= PROGRESS_INTERVAL {
        log::info!(target: "search", "{scored} spectra scored (~{rate:.0}/s)");
        *last = Some(std::time::Instant::now());
    } else {
        log::debug!(target: "search", "{scored} spectra scored (~{rate:.0}/s)");
    }
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Train a scoring model from externally labelled, high-confidence PSMs
    /// given as flat training parquets, without a bootstrap search.
    #[command(hide = true)]
    Train(Box<TrainArgs>),

    /// Train a scoring model from spectra and a FASTA database, writing the
    /// result to a Parquet model store.
    #[command(name = "train-from-search", hide = true)]
    TrainFromSearch(Box<TrainFromSearchArgs>),

    /// Train a v3 GBDT fragment-intensity regressor from flat training
    /// parquets and embed it in a Parquet model store alongside existing
    /// rank-core models.
    #[command(name = "train-intensity-gbdt", hide = true)]
    TrainIntensityGbdt(Box<TrainIntensityGbdtArgs>),

    /// Train a GBDT rich-ion LLR classifier (logistic; decoy-aware) from flat
    /// training parquets and embed it in a Parquet model store alongside
    /// existing rank-core models.
    #[command(name = "train-rich-ion-llr", hide = true)]
    TrainRichIonLlr(Box<TrainRichIonLlrArgs>),

    /// Match between runs in stages over a --quant-dir: library, extract (per run), finalize
    ///
    /// The same result as a search with --lfq --mbr on the same inputs.
    #[command(subcommand)]
    Quant(cli::QuantCommand),

    /// Analyse the SDRFs of several projects together
    ///
    /// `plan` reads the `cohort` section of a config file and the SDRFs it lists, checks
    /// them, and writes every run's settings and the groups that share one FDR as JSON.
    #[command(subcommand)]
    Cohort(cli::CohortCommand),

    /// Target + decoy database of a study or cohort, with the decoy-target overlap report
    ///
    /// Writes the FASTA of the targets and their decoys, a settings table and an overlap
    /// report. Every search given the FASTA as --database searches these decoys, and stops
    /// when its decoy settings or digest are not the database's.
    Database(Box<cli::DatabaseCmdArgs>),

    /// Protein groups with protein-level q-values over the rescored PSMs of one or more runs
    ///
    /// Writes the QPX protein-group view.
    Inference(Box<cli::InferenceArgs>),

    #[command(
        hide = true,
        about = "Run the native rescoring engine on a Percolator command line",
        long_about = "Run the native rescoring engine on a Percolator command line\n\n\
The engine implements the Percolator algorithm (Käll et al., Nat. Methods 2007) and tracks \
upstream Percolator master (commit a9783ee8); it takes Percolator's arguments and writes \
Percolator's result files."
    )]
    Rescore(Box<cli::RescoreArgs>),
}

const USAGE: &str = "andes --spectrum <FILE>... --database <FILE> --output-pin <FILE> [OPTIONS]
       andes <COMMAND> [OPTIONS]";

const AFTER_HELP: &str = "Example:
  andes --spectrum run.mzML --database human.fasta --output-pin run.pin --rescore

`andes --help` shows every option.";

const AFTER_LONG_HELP: &str = "Example:
  andes --spectrum run.mzML --database human.fasta --output-pin run.pin --rescore

The full reference is DOCS.md.";

#[derive(Parser, Debug)]
#[command(
    name = "andes",
    version = env!("CARGO_PKG_VERSION"),
    about = "andes: peptide identification and quantification from MS/MS spectra \
             (database search, rescoring, protein groups, TMT, label-free, SILAC)",
    override_usage = USAGE,
    after_help = AFTER_HELP,
    after_long_help = AFTER_LONG_HELP,
    allow_hyphen_values = true,
)]
struct TopCli {
    #[command(subcommand)]
    command: Option<Command>,

    /// Print only warnings and errors to stderr
    #[arg(long, global = true)]
    quiet: bool,

    /// Also print debug lines to stderr
    #[arg(long, global = true)]
    verbose: bool,

    #[command(flatten)]
    search: SearchArgs,
}

fn main() -> ExitCode {
    PROCESS_START.get_or_init(std::time::Instant::now);
    #[cfg(feature = "thermo")]
    configure_bundled_dotnet();
    let matches = <TopCli as clap::CommandFactory>::command().get_matches();
    let mut top =
        <TopCli as clap::FromArgMatches>::from_arg_matches(&matches).unwrap_or_else(|e| e.exit());
    logging::init(logging::level_from_flags(top.quiet, top.verbose));
    let result = match top.command.take() {
        Some(Command::Train(args)) => run_train(*args),
        Some(Command::TrainFromSearch(args)) => run_train_from_search(*args),
        Some(Command::TrainIntensityGbdt(args)) => run_train_intensity_gbdt(*args),
        Some(Command::TrainRichIonLlr(args)) => run_train_rich_ion_llr(*args),
        Some(Command::Quant(mut cmd)) => {
            if let cli::QuantCommand::Finalize(a) = &mut cmd {
                let sub = matches
                    .subcommand_matches("quant")
                    .and_then(|q| q.subcommand_matches("finalize"));
                if let Some(Err(e)) = sub.map(|m| config::apply_protein_config(&mut a.proteins, m))
                {
                    log::error!("{e}");
                    return ExitCode::from(2);
                }
            }
            quant_cmd::run_quant_command(cmd)
        }
        Some(Command::Inference(mut args)) => {
            if let Some(Err(e)) = matches
                .subcommand_matches("inference")
                .map(|m| config::apply_protein_config(&mut args.proteins, m))
            {
                log::error!("{e}");
                return ExitCode::from(2);
            }
            inference_cmd::run_inference(*args)
        }
        Some(Command::Rescore(args)) => run_rescore(*args),
        Some(Command::Database(args)) => database_cmd::run(*args),
        Some(Command::Cohort(cmd)) => cohort_cmd::run(cmd),
        None => {
            let mut cfg_missed_cleavages = false;
            let mut sources = param_sources::from_command_line(&matches);
            if let Some(cfg_path) = top.search.config.clone() {
                if let Err(e) = config::RunConfig::load(&cfg_path).and_then(|c| {
                    cfg_missed_cleavages = c.search.max_missed_cleavages.is_some();
                    config::apply(c, &mut top.search, &matches, &mut sources)
                }) {
                    log::error!("{e}");
                    return ExitCode::from(2);
                }
            }
            if let Some(sdrf_path) = top.search.sdrf.clone() {
                match sdrf::read(&sdrf_path)
                    .and_then(|s| sdrf::apply(&s, &mut top.search, &matches, &mut sources))
                {
                    Ok(lines) => lines.iter().for_each(|l| match l.strip_prefix("WARN: ") {
                        Some(w) => log::warn!(target: "sdrf", "{w}"),
                        None => log::info!(target: "sdrf", "{l}"),
                    }),
                    Err(e) => {
                        log::error!("{e}");
                        return ExitCode::from(2);
                    }
                }
            }
            if let Err(e) = config::check_resolved(&top.search) {
                log::error!("{e}");
                return ExitCode::from(2);
            }
            if let Some(dest) = top.search.print_config.as_deref() {
                return print_config(dest, &sources);
            }
            let _ = EXPLICIT_MISSED_CLEAVAGES.set(missed_cleavages_explicit(
                matches.value_source("max_missed_cleavages"),
                cfg_missed_cleavages,
            ));
            let search = top.search;
            if search.spectrum.is_empty() {
                log::error!("--spectrum is required for search (or use `andes train`)");
                return ExitCode::from(2);
            }
            let database = match search.database {
                Some(p) => p,
                None => {
                    log::error!("--database is required for search");
                    return ExitCode::from(2);
                }
            };
            if search.output_pin.is_none() && !search.rescore {
                log::error!("--output-pin is required for search (or use --rescore)");
                return ExitCode::from(2);
            }
            if (search.psm_fdr.is_some() || search.pep.is_some()) && !search.rescore {
                log::warn!(
                    "--psm-fdr/--pep set without --rescore; they are ignored (add --rescore \
                     to rescore with Percolator, or feed the PIN to Percolator yourself)."
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
            log::error!("{e}");
            ExitCode::from(1)
        }
    }
}

fn print_config(dest: &std::path::Path, sources: &param_sources::ParamSources) -> ExitCode {
    let text = match params_render::effective_config(sources) {
        Ok(text) => text,
        Err(e) => {
            log::error!("{e}");
            return ExitCode::from(2);
        }
    };
    let cli_only = params_render::command_line_only(sources);
    if !cli_only.is_empty() {
        log::warn!(
            "--print-config: {} has no config key; the printed config lists it as a comment, \
             pass it again on the command line to reproduce the run",
            cli_only.join(" ")
        );
    }
    let written = if dest == std::path::Path::new("-") {
        use std::io::Write;
        std::io::stdout().write_all(text.as_bytes())
    } else {
        std::fs::write(dest, &text)
    };
    match written {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            log::error!("cannot write the config to {}: {e}", dest.display());
            ExitCode::from(1)
        }
    }
}

fn run_rescore(args: cli::RescoreArgs) -> Result<(), Box<dyn std::error::Error>> {
    cli::ensure_scratch_dir(&args.temp_dir)?;
    let limit = memlimit::rescore_memory(args.rescore_memory_gb, args.temp_dir);
    output::run_native_percolator(&args.args, limit)?;
    Ok(())
}

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

fn arg_present(flag: &str) -> bool {
    std::env::args().any(|a| a == flag || a.starts_with(&format!("{flag}=")))
}

static EXPLICIT_MISSED_CLEAVAGES: std::sync::OnceLock<bool> = std::sync::OnceLock::new();

fn missed_cleavages_explicit(
    cli_source: Option<clap::parser::ValueSource>,
    in_config: bool,
) -> bool {
    in_config || matches!(cli_source, Some(clap::parser::ValueSource::CommandLine))
}

#[cfg(test)]
mod missed_cleavages_explicit_tests {
    use super::missed_cleavages_explicit;
    use clap::parser::ValueSource;

    #[test]
    fn config_value_is_explicit() {
        assert!(missed_cleavages_explicit(
            Some(ValueSource::DefaultValue),
            true
        ));
        assert!(missed_cleavages_explicit(
            Some(ValueSource::CommandLine),
            false
        ));
        assert!(!missed_cleavages_explicit(
            Some(ValueSource::DefaultValue),
            false
        ));
        assert!(!missed_cleavages_explicit(None, false));
    }
}

#[cfg(test)]
#[path = "../tests/params_table.rs"]
mod params_table;

#[cfg(test)]
#[path = "../tests/params_sources.rs"]
mod params_sources;

#[cfg(test)]
#[path = "../tests/params_render.rs"]
mod params_render_tests;
