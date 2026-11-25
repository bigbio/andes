//! The start and the input checks of a search: the thread pool and the cross-link flags, the MGF
//! MS-level warning, and the parse report after the spectra are read.

use std::path::PathBuf;

use model::Spectrum;

use super::files::{search_file, FileInputs, FileOutputs};
use super::{check_flag_combinations, resolve};
use crate::cli::Cli;
use crate::cli::Fragmentation;
use crate::spectra::{learn_file_calibrations, ms3_carries_reporters, ParseStats};
use model::activation::ActivationMethod;
use scoring_crate::RankScorer;
use search::{PrecursorCalMode, PreparedSearch};

/// Checks the inputs, builds the thread pool, and applies the cross-link flags (the default
/// variable mods, `--xl-linear-only`), then the flag combination checks.
pub(super) fn prepare_cli(cli: Cli) -> Result<Cli, Box<dyn std::error::Error>> {
    if cli.spectrum.is_empty() {
        return Err("no --spectrum inputs".into());
    }
    static POOL_INIT: std::sync::OnceLock<()> = std::sync::OnceLock::new();
    POOL_INIT.get_or_init(|| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(cli.threads)
            .build_global()
            .expect("build_global");
    });
    let mut cli = cli;
    if let Some(x) = &cli.xl {
        if let Some(v) = crate::xl::xl_default_variable_mods(cli.db.variable_mods.as_ref(), &x.0) {
            cli.db.variable_mods = Some(v);
        }
        if cli.xl_linear_only {
            log::info!(
                target: "xl",
                "--xl-linear-only: linear and monolinked peptides only; no cross-link search"
            );
            cli.xl = None;
        }
    }
    check_flag_combinations(&cli)?;
    Ok(cli)
}

/// The input files of the run report, as `(role, path)`: each spectrum file as the user gave it,
/// the database, and the SDRF and config file when given.
pub(super) fn report_inputs(
    cli: &Cli,
    spectra: &[PathBuf],
    database: &std::path::Path,
) -> Vec<(String, PathBuf)> {
    let role = |r: &str, p: &std::path::Path| (r.to_string(), p.to_path_buf());
    spectra
        .iter()
        .map(|p| role("spectrum", p))
        .chain(std::iter::once(role("database", database)))
        .chain(cli.sdrf.iter().map(|p| role("sdrf", p)))
        .chain(cli.config.iter().map(|p| role("config", p)))
        .collect()
}

/// Warns that `--ms-level` has no effect on an MGF input, which carries no MS level.
pub(super) fn warn_mgf_ms_level(cli: &Cli, is_mgf: bool) {
    if is_mgf && cli.ms_level != 2 {
        log::warn!(
            "--ms-level={} requested for an MGF input; MGF carries no MS \
             level (always treated as MS2). The flag has no effect on this input.",
            cli.ms_level
        );
    }
}

/// Reports the spectra that failed to parse and the MS-level filter, and fails when no spectrum
/// was read.
pub(super) fn check_parsed(
    cli: &Cli,
    parse_stats: &ParseStats,
    ms_level_filtered: bool,
    all_spectra: &[Spectrum],
    spectrum_paths: &[PathBuf],
) -> Result<(), Box<dyn std::error::Error>> {
    if parse_stats.error_count > 0 {
        log::warn!(
            target: "input",
            "{} spectra failed to parse{}",
            parse_stats.error_count,
            if !parse_stats.first_errors.is_empty() {
                format!(" (first {}):", parse_stats.first_errors.len())
            } else {
                String::new()
            }
        );
        for e in &parse_stats.first_errors {
            log::warn!(target: "input", "- {e}");
        }
    }

    if ms_level_filtered {
        log::info!(
            target: "input",
            "MS-level filter {} (only MS{} spectra entered the search)",
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
    Ok(())
}

/// What routing the per-scan scorers by activation reads.
pub(super) struct RouteInputs<'a> {
    pub(super) cli: &'a Cli,
    pub(super) standalone: bool,
    pub(super) instrument_was_detected: bool,
    pub(super) resolved_model_id: &'a String,
    pub(super) param: &'a scoring_crate::Param,
    pub(super) acquisition: &'a crate::acquisition::Acquisition,
    pub(super) params: &'a search::SearchParams,
    pub(super) gbdt_k: usize,
    pub(super) routed_instrument: Option<model::instrument::InstrumentType>,
    pub(super) route_protocol: crate::cli::Protocol,
    pub(super) search_enzyme: model::enzyme::Enzyme,
    pub(super) frag_tol_override: Option<model::tolerance::Tolerance>,
    pub(super) routed_activation: ActivationMethod,
}

/// With an auto-detected instrument and a resolved model, loads a scorer per activation on
/// demand, so scans of another activation are scored with their own model.
pub(super) fn route_per_scan_scorers<'p>(
    prepared: search::PreparedSearch<'p>,
    inputs: RouteInputs<'_>,
) -> search::PreparedSearch<'p> {
    let RouteInputs {
        cli,
        standalone,
        instrument_was_detected,
        resolved_model_id,
        param,
        acquisition,
        params,
        gbdt_k,
        routed_instrument,
        route_protocol,
        search_enzyme,
        frag_tol_override,
        routed_activation,
    } = inputs;
    if cli.model_id_override.is_none()
        && cli.fragmentation == Fragmentation::Auto
        && !standalone
        && instrument_was_detected
        && !resolved_model_id.is_empty()
    {
        let store_path = cli
            .model_store
            .clone()
            .unwrap_or_else(crate::models::bundled_store_path);
        let run_model_id = resolved_model_id.clone();
        let run_analyser =
            crate::models::analyser_of(param.data_type.instrument.name()).to_string();
        let run_protocol = param.data_type.protocol;
        let high_res = acquisition.high_res;
        let index_on = params.fragment_index_top_k > 0;
        let gbdt_max_trees = gbdt_k;
        let no_search_gbdt = cli.no_search_gbdt;
        let no_search_rich_ion = cli.no_search_rich_ion;
        let loader: search::ScorerLoader = Box::new(move |activation: ActivationMethod| {
            let route = match crate::models::load_activation_route(
                activation,
                routed_instrument,
                route_protocol,
                search_enzyme,
                &store_path,
                &run_model_id,
                &run_analyser,
            ) {
                Ok(r) => r?,
                Err(e) => {
                    log::warn!(
                        target: "routing",
                        "{e}; {} scans use the run's model",
                        activation.name()
                    );
                    return None;
                }
            };
            let (model_id, mut p) = route;
            if index_on && !search::fragments::index::covers_activation(p.data_type.activation) {
                return None;
            }
            if gbdt_max_trees > 0 {
                for slot in [&mut p.frag_intensity_model, &mut p.rich_ion_model] {
                    resolve::truncate_gbdt(slot, gbdt_max_trees);
                }
            }
            if no_search_gbdt {
                p.frag_intensity_model = None;
                p.rich_ion_model = None;
            }
            if no_search_rich_ion {
                p.rich_ion_model = None;
            }
            p.data_type.protocol = run_protocol;
            let mut sc = RankScorer::new(&p);
            sc.set_acquisition_high_res(Some(high_res));
            if let Some(t) = frag_tol_override {
                sc.set_fragment_tol_override(Some(t));
            }
            log::info!(
                target: "routing",
                "{} scans -> model {model_id} (the run is routed {})",
                activation.name(),
                routed_activation.name()
            );
            Some(sc)
        });
        return prepared.with_per_scan_scorers(routed_activation, loader);
    }
    prepared
}

/// The MS level the TMT reporters are read from: `--tmt-level`, else MS3 when the MS2 is
/// low-resolution for the plex and the MS3 scans carry reporters (SPS-MS3), else MS2.
pub(super) fn tmt_level_for(
    cli: &Cli,
    high_res: bool,
    spectrum_path: &std::path::Path,
    is_mzml: bool,
    is_raw: bool,
) -> u8 {
    match (cli.tmt_level, cli.tmt) {
        (Some(level), _) => level,
        (None, Some(plex))
            if !high_res
                && plex
                    .check_tolerance(model::tolerance::Tolerance::Da(0.3))
                    .is_err()
                && ms3_carries_reporters(spectrum_path, is_mzml, is_raw) =>
        {
            log::info!(
                target: "protocol",
                "--tmt {}: the MS2 is low-resolution and the MS3 scans carry \
                 reporters (SPS-MS3); reading reporters from the MS3 (--tmt-level 3)",
                plex.name()
            );
            3
        }
        _ => 2,
    }
}

/// Whether the search captures MS3 scans for a cleavable cross-linker.
pub(super) fn xl_ms3_capture(cli: &Cli) -> bool {
    cli.xl
        .as_ref()
        .is_some_and(|x| x.0.iter().any(|l| l.cleavage.is_some()))
}

/// The precursor calibration prepass: a calibration per file (mass shift, tolerance) learned on a
/// pass over each file, with the first file's applied to `params`; empty when calibration is off.
pub(super) fn calibration_prepass(
    cli: &Cli,
    (idx, scorer, fragment_tol_da, index_path): (
        &search::SearchIndex,
        &RankScorer,
        f64,
        &std::path::Path,
    ),
    params: &mut search::SearchParams,
    (spectrum_paths, ms_level_u32, bench_cap): (&[PathBuf], u32, usize),
) -> Result<Vec<(f64, model::tolerance::PrecursorTolerance)>, Box<dyn std::error::Error>> {
    let mut file_cal: Vec<(f64, model::tolerance::PrecursorTolerance)> = Vec::new();
    if params.precursor_cal_mode != PrecursorCalMode::Off {
        let t_cal = std::time::Instant::now();
        let mut cal_prepared = PreparedSearch::prepare_mmap(
            idx,
            &*params,
            scorer,
            fragment_tol_da,
            &cli.db.decoy_prefix,
            index_path,
        )
        .map_err(|e| format!("build out-of-core candidate index: {e}"))?;
        file_cal = learn_file_calibrations(
            spectrum_paths,
            ms_level_u32,
            bench_cap,
            &*params,
            &mut cal_prepared,
        )?;
        drop(cal_prepared);
        log::info!(
            target: "calibration_prepass",
            "done in {:.2} s (includes opening or building the candidate index)",
            t_cal.elapsed().as_secs_f64()
        );
        (params.precursor_mass_shift_ppm, params.precursor_tolerance) = file_cal[0];
    }
    Ok(file_cal)
}

/// Opens `--precursor-mono-dump` when the monoisotopic correction is active (and warns that it is
/// ignored otherwise).
pub(super) fn open_mono_dump(
    cli: &Cli,
    mono_active: bool,
    mono_params: &search::precursor::monoisotopic::MonoParams,
) -> Result<Option<crate::monoisotopic::MonoDump>, Box<dyn std::error::Error>> {
    Ok(match (&cli.precursor_mono_dump, mono_active) {
        (Some(p), true) => Some(
            crate::monoisotopic::MonoDump::create(p, mono_params.max_shift)
                .map_err(|e| format!("create --precursor-mono-dump {}: {e}", p.display()))?,
        ),
        (Some(p), false) => {
            log::warn!(
                "--precursor-mono-dump {} ignored: the correction is not active.",
                p.display()
            );
            None
        }
        (None, _) => None,
    })
}

/// The spectrum cap of a `--max-spectra` bench run (none otherwise); warns that `--ms-level` is
/// ignored for native .raw/.d inputs.
pub(super) fn bench_cap(cli: &Cli, native_input: bool) -> usize {
    let bench_cap = if cli.max_spectra > 0 {
        cli.max_spectra
    } else {
        usize::MAX
    };
    if native_input && cli.ms_level != 2 {
        log::warn!(
            "--ms-level={} is ignored for native .raw/.d input; these formats \
             search MS2 (identification) scans only.",
            cli.ms_level
        );
    }
    bench_cap
}

/// The search parameters of each file: the run's, with the file's calibrated mass shift, and its
/// tolerance when the run kept the first file's.
pub(super) fn file_params(
    params: &search::SearchParams,
    file_cal: &[(f64, model::tolerance::PrecursorTolerance)],
) -> Vec<search::SearchParams> {
    file_cal
        .iter()
        .map(|&(shift, tolerance)| {
            let mut p = params.clone();
            p.precursor_mass_shift_ppm = shift;
            if params.precursor_tolerance == file_cal[0].1 {
                p.precursor_tolerance = tolerance;
            }
            p
        })
        .collect()
}

/// Ends opening the candidate index: the DIA RT form model, and the index logs.
pub(super) fn finish_index_open(
    cli: &Cli,
    prepared: &mut PreparedSearch<'_>,
    t_index: std::time::Instant,
) -> Result<(), Box<dyn std::error::Error>> {
    if cli.dia {
        prepared.set_form_rt_model(crate::dia::rt::andes_rt_form_model()?);
    }
    log::info!(
        target: "candidate_index_open",
        "done in {:.2} s",
        t_index.elapsed().as_secs_f64()
    );
    log::info!(
        target: "candidate-index",
        "out-of-core mmap (base peptides resolved lazily per spectrum)"
    );
    Ok(())
}

/// The predictors of the run (none for glyco or cross-link searches) and whether the spectra's
/// peaks are kept after the search for them, refine, DIA or IPF.
pub(super) fn predictors(
    cli: &Cli,
    standalone: bool,
    acquisition: &crate::acquisition::Acquisition,
    instrument: model::instrument::InstrumentType,
) -> Result<(Option<crate::dia::rt::Predictors>, bool), Box<dyn std::error::Error>> {
    let predictors = if standalone {
        None
    } else {
        crate::dia::rt::Predictors::from_cli(cli, acquisition, instrument)?
    };
    if let Some(p) = &predictors {
        log::info!(target: "predict", "{}", p.describe());
    }
    let keep_peaks = cli.refine
        || cli.dia
        || cli.localize.is_some()
        || predictors.as_ref().is_some_and(|p| p.needs_peaks());
    Ok((predictors, keep_peaks))
}

/// Searches each spectrum file in turn, stopping at the bench cap.
pub(super) fn search_files<'a>(
    file_inputs: &FileInputs<'a>,
    outputs: &mut FileOutputs,
    prepared: &mut PreparedSearch<'a>,
    measurer: &mut super::calibrate::FragmentErrorMeasurer,
    quant: &mut Option<crate::quant::run::QuantCollector>,
    (bench_mode, bench_cap): (bool, usize),
) -> Result<(), Box<dyn std::error::Error>> {
    for (file_idx, input_path) in file_inputs.spectrum_paths.iter().enumerate() {
        if bench_mode && outputs.all_spectra.len() >= bench_cap {
            break;
        }
        search_file(
            file_inputs,
            outputs,
            prepared,
            measurer,
            quant,
            file_idx,
            input_path,
        )?;
    }
    Ok(())
}

/// Each file's precursor calibration as the PIN records it: where its spectra end, its mass shift
/// and tolerance.
pub(super) fn file_calibrations(
    file_params: &[search::SearchParams],
    run_ends: &[usize],
) -> Vec<search::FileCalibration> {
    file_params
        .iter()
        .zip(run_ends)
        .map(|(p, &end)| search::FileCalibration {
            end,
            shift_ppm: p.precursor_mass_shift_ppm,
            tolerance: p.precursor_tolerance,
        })
        .collect()
}
