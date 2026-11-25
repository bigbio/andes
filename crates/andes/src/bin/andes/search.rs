//! `andes search`: resolve the parameters, search every input file, rescore, and write the PIN,
//! the identifications and the quantification tables. `run` holds the stages in order; each
//! stage's work is in a child module.

mod calibrate;
mod collect;
mod files;
mod finish;
mod outputs;
mod rescored;
mod rescoring;
mod resolve;
mod setup;
mod stream;

use std::path::{Path, PathBuf};

use model::activation::ActivationMethod;
use scoring_crate::RankScorer;
use search::candidates::records::index_cache_path;
use search::precursor::monoisotopic::MonoCorrection;
use search::PreparedSearch;

use crate::cli::Cli;
use crate::glyco::run_glyco;
use crate::models::search_aa_set;
use crate::monoisotopic::{MonoDump, MonoStats};
use calibrate::{finish_mono_correction, optimize_on_sample, FragmentErrorMeasurer};
use collect::quant_collector;
use files::{FileInputs, FileOutputs};
use finish::{post_search, PostSearch, PostSearchOwned};
pub(crate) use outputs::load_claim_database;
use outputs::{
    infer_and_write_groups, write_filtered_tsv, write_open_offsets, write_qpx_bundle, ClaimDatabase,
};
use rescored::{refine_rescored, RescoredInputs};
use resolve::{
    check_flag_combinations, detect_inputs, load_search_database, log_resolved_parameters,
    prepare_run, resolve_scoring_model, resolve_search_params, ResolvedSearch, RunInputs, RunPaths,
    ScoringModel, SearchDatabase,
};

type Ms1Peaks = Vec<(f64, f32)>;

fn scratch_tempdir(prefix: &str, scratch_dir: &Path) -> Result<tempfile::TempDir, String> {
    crate::cli::ensure_scratch_dir(scratch_dir)?;
    tempfile::Builder::new()
        .prefix(prefix)
        .tempdir_in(scratch_dir)
        .map_err(|e| crate::cli::scratch_dir_error(scratch_dir, &e))
}

pub(crate) fn run(
    cli: Cli,
    sources: crate::params::sources::ParamSources,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut clock = crate::report::StageClock::new();
    clock.stage("inputs");
    let mut cli = setup::prepare_cli(cli)?;
    let RunPaths {
        scratch_dir,
        original_spectrum_paths,
        output_pin_path,
        report_base,
        temp_dirs: _temp_dirs,
    } = prepare_run(&mut cli)?;
    let standalone = cli.glyco || cli.xl.is_some();
    let database_path: PathBuf = cli.database.clone().expect("database validated in main");
    let hasher = crate::report::InputHasher::start(setup::report_inputs(
        &cli,
        &original_spectrum_paths,
        &database_path,
    ));

    let t_total = std::time::Instant::now();
    clock.stage("database");
    let db = load_search_database(&cli, &database_path)?;
    clock.stage("model");

    let mut aa = search_aa_set(&cli.db, None)?;
    let inputs = detect_inputs(&cli, &aa)?;
    let model = resolve_scoring_model(&mut cli, &inputs, &mut aa, db.search_enzyme)?;
    let spectrum_paths = &cli.spectrum;
    let spectrum_path: PathBuf = spectrum_paths[0].clone();
    let ResolvedSearch {
        mut params,
        chimeric_active,
        mono_active,
        mono_params,
        restricted_site_masses,
        fixed_terminal_mods,
    } = resolve_search_params(&cli, aa, &inputs, &model, &db, standalone)?;
    log_resolved_parameters(&cli, &inputs, &model, &db, &params);
    let SearchDatabase {
        target_db,
        idx,
        search_enzyme,
        decoy_strategy,
        decoy_key,
        ..
    } = db;

    let t_phase = std::time::Instant::now();

    log::info!(target: "run", "Using {} worker threads", cli.threads);

    let fragment_tol_da = 0.5_f64;

    let bench_cap = setup::bench_cap(&cli, inputs.is_raw || inputs.is_d);
    let ms_level_u32 = cli.ms_level as u32;

    let mut trace = run_trace(
        (clock, hasher, sources),
        &model.resolved_model_id,
        decoy_strategy,
        decoy_key,
    );
    trace.clock.stage("calibration");
    let cache_dir = crate::cli::dir_or_output_dir(cli.index_cache_dir.as_deref(), &report_base);
    let index_path = index_cache_path(&idx, &params, &cache_dir);
    let file_cal = calibrate_params(
        &cli,
        (&inputs, &spectrum_path, &scratch_dir),
        (&idx, &model.scorer, fragment_tol_da, &index_path),
        &mut params,
        (spectrum_paths, ms_level_u32, bench_cap),
    )?;

    let file_params = setup::file_params(&params, &file_cal);
    let bench_mode = cli.max_spectra > 0;
    let mut quant = quant_for_run(&cli, &inputs, &model.scorer, &spectrum_path)?;
    let xl_ms3_capture = setup::xl_ms3_capture(&cli);
    let ms1_input = chimeric_active || mono_active;
    let prefetched = stream::prefetch_first_file(
        &cli,
        (&params, standalone, ms1_input, xl_ms3_capture),
        (bench_mode, bench_cap),
        quant.as_ref(),
    );
    trace.clock.stage("index");
    let mut prepared = prepare_index(
        (&idx, &params, &model.scorer, fragment_tol_da, &index_path),
        route_inputs(&cli, standalone, &inputs, &model, &params, search_enzyme),
    )?;

    let (predictors, keep_peaks) = setup::predictors(
        &cli,
        standalone,
        &inputs.acquisition,
        model.param.data_type.instrument,
    )?;
    let t_search_start = std::time::Instant::now();

    setup::warn_mgf_ms_level(&cli, inputs.is_mgf);

    let mut measurer = FragmentErrorMeasurer {
        enabled: inputs.acquisition.high_res && model.frag_tol_override.is_none() && !standalone,
        done: false,
        per_file: Vec::with_capacity(spectrum_paths.len()),
    };
    let file_inputs = FileInputs {
        cli: &cli,
        params: &params,
        file_params: &file_params,
        spectrum_paths,
        bench_cap,
        bench_mode,
        chimeric_active,
        mono_active,
        mono_params: &mono_params,
        ms1_input,
        ms_level_u32,
        keep_peaks,
        standalone,
        t_search_start,
        xl_ms3_capture,
    };
    let mut outputs = FileOutputs::new(
        spectrum_paths.len(),
        setup::open_mono_dump(&cli, mono_active, &mono_params)?,
    );
    outputs.prefetched = prefetched;
    trace.clock.stage("search");
    let FileOutputs {
        all_spectra,
        all_queues,
        run_ends,
        dia_ms1,
        xl_ms3,
        xl_snapped,
        xl_snap_files,
        mono_table,
        mono_stats,
        mono_dump,
        ..
    } = search_all_files(
        &file_inputs,
        outputs,
        (&mut prepared, &mut measurer, &mut quant),
        inputs.is_mzml || inputs.is_raw,
    )?;
    prepared.params = &params;
    let file_calibrations = setup::file_calibrations(&file_params, &run_ends);

    trace.clock.stage(if cli.refine { "refine" } else { "post-search" });
    let Some((spectra, mut queues, refine_output)) = post_search(
        &PostSearch {
            t_phase,
            xl_snapped,
            cli: &cli,
            acquisition: &inputs.acquisition,
            decoy_strategy,
            fragment_tol_da,
            mono_active,
            mono_params: &mono_params,
            routed_activation: model.routed_activation,
            spectrum_paths,
            standalone,
            t_search_start,
            xl_snap_files,
            bench_mode,
            file_calibrations: &file_calibrations,
            idx: &idx,
            output_pin_path: &output_pin_path,
            params: &params,
            quant: &quant,
            report_base: &report_base,
            run_ends: &run_ends,
            scorer: &model.scorer,
            spectrum_path: &spectrum_path,
            t_total,
        },
        PostSearchOwned {
            all_spectra,
            all_queues,
            mono_table,
            mono_dump,
            mono_stats,
            xl_ms3,
        },
        &mut prepared,
    )?
    else {
        return Ok(());
    };
    let t_phase = std::time::Instant::now();
    let (refine_merged, ipf_scoring, pin_candidates, merged_index) = pin_tables(
        &cli,
        refine_output,
        &mut queues,
        prepared,
        (&idx, &model.scorer),
    );
    let pin_index = merged_index.unwrap_or(idx);

    finish::write_pin_and_rescore(finish::PinAndRescore {
        cli,
        t_total,
        spectrum_path,
        bench_mode,
        quant,
        predictors,
        run_ends,
        dia_ms1,
        measurer,
        file_calibrations,
        spectra,
        queues,
        target_db,
        scratch_dir,
        original_spectrum_paths,
        output_pin_path,
        report_base,
        params,
        restricted_site_masses,
        fixed_terminal_mods,
        scorer: &model.scorer,
        is_mgf: inputs.is_mgf,
        t_phase,
        refine_merged,
        ipf_scoring,
        pin_candidates,
        pin_index,
        trace,
    })
}

/// From `pin_tables`: whether refinement candidates were merged, the IPF scoring, the PIN's
/// candidates, and the merged index when refinement ran.
type PinTables<'s> = (
    bool,
    Option<search::engine::IpfScoring<'s>>,
    Vec<search::candidates::Candidate>,
    Option<search::SearchIndex>,
);

/// The run report's trace: the stage clock, the input hasher and the parameter sources, with the
/// scoring model and the decoys of the search.
fn run_trace(
    (clock, hasher, sources): (
        crate::report::StageClock,
        crate::report::InputHasher,
        crate::params::sources::ParamSources,
    ),
    resolved_model_id: &str,
    decoy_strategy: search::decoy::DecoyStrategy,
    decoy_key: Option<String>,
) -> crate::report::build::Trace {
    crate::report::build::Trace {
        clock,
        hasher,
        sources,
        command_line: crate::report::command_line(),
        models: vec![("scoring".to_string(), resolved_model_id.to_string())],
        decoy_strategy,
        decoy_key,
    }
}

/// Calibrates the precursor tolerance per file (the prepass), then, with `--optimize`, tunes the
/// search parameters on a sample of the first file; returns the per-file calibration.
fn calibrate_params(
    cli: &Cli,
    (inputs, spectrum_path, scratch_dir): (&RunInputs, &Path, &Path),
    index: (&search::SearchIndex, &RankScorer, f64, &Path),
    params: &mut search::SearchParams,
    files: (&[PathBuf], u32, usize),
) -> Result<Vec<(f64, model::tolerance::PrecursorTolerance)>, Box<dyn std::error::Error>> {
    let file_cal = setup::calibration_prepass(cli, index, params, files)?;
    if cli.optimize {
        optimize_on_sample(
            cli,
            spectrum_path,
            (inputs.is_mzml, inputs.is_raw, inputs.is_mgf),
            scratch_dir,
            index,
            params,
        )?;
    }
    Ok(file_cal)
}

/// The run's quantification collector (label-free, TMT at the level the input supports, or
/// none), from the flags and the detected input.
fn quant_for_run(
    cli: &Cli,
    inputs: &RunInputs,
    scorer: &RankScorer,
    spectrum_path: &Path,
) -> Result<Option<crate::quant::run::QuantCollector>, Box<dyn std::error::Error>> {
    let high_res = inputs.acquisition.high_res;
    let tmt_level =
        setup::tmt_level_for(cli, high_res, spectrum_path, inputs.is_mzml, inputs.is_raw);
    quant_collector(
        cli,
        high_res,
        tmt_level,
        scorer,
        inputs.is_mgf,
        inputs.is_d,
        &cli.spectrum,
    )
}

/// What `setup::route_per_scan_scorers` reads: the flags, the detected input, the scoring model
/// and the search parameters.
fn route_inputs<'r>(
    cli: &'r Cli,
    standalone: bool,
    inputs: &'r RunInputs,
    model: &'r ScoringModel,
    params: &'r search::SearchParams,
    search_enzyme: model::enzyme::Enzyme,
) -> setup::RouteInputs<'r> {
    setup::RouteInputs {
        cli,
        standalone,
        instrument_was_detected: inputs.instrument_was_detected,
        resolved_model_id: &model.resolved_model_id,
        param: &model.param,
        acquisition: &inputs.acquisition,
        params,
        gbdt_k: model.gbdt_k,
        routed_instrument: model.routed_instrument,
        route_protocol: inputs.route_protocol,
        search_enzyme,
        frag_tol_override: model.frag_tol_override,
        routed_activation: model.routed_activation,
    }
}

/// Opens (or builds) the out-of-core candidate index, routes a scorer per scan, and finishes the
/// opening (the DIA retention-time model and the index log line).
fn prepare_index<'a>(
    (idx, params, scorer, fragment_tol_da, index_path): (
        &'a search::SearchIndex,
        &'a search::SearchParams,
        &'a RankScorer,
        f64,
        &Path,
    ),
    route: setup::RouteInputs<'_>,
) -> Result<PreparedSearch<'a>, Box<dyn std::error::Error>> {
    let cli = route.cli;
    let t_index = std::time::Instant::now();
    let prepared = PreparedSearch::prepare_mmap(
        idx,
        params,
        scorer,
        fragment_tol_da,
        &cli.db.decoy_prefix,
        index_path,
    )
    .map_err(|e| format!("build out-of-core candidate index: {e}"))?;
    let mut prepared = setup::route_per_scan_scorers(prepared, route);
    setup::finish_index_open(cli, &mut prepared, t_index)?;
    Ok(prepared)
}

/// Searches every input file into `outputs`; the run ends then cover every file, and the parse
/// is checked (`ms_level_filtered`: the reader filters by MS level).
fn search_all_files<'a>(
    file_inputs: &FileInputs<'a>,
    mut outputs: FileOutputs,
    (prepared, measurer, quant): (
        &mut PreparedSearch<'a>,
        &mut FragmentErrorMeasurer,
        &mut Option<crate::quant::run::QuantCollector>,
    ),
    ms_level_filtered: bool,
) -> Result<FileOutputs, Box<dyn std::error::Error>> {
    let spectrum_paths = file_inputs.spectrum_paths;
    setup::search_files(
        file_inputs,
        &mut outputs,
        prepared,
        measurer,
        quant,
        (file_inputs.bench_mode, file_inputs.bench_cap),
    )?;
    outputs
        .run_ends
        .resize(spectrum_paths.len(), outputs.all_spectra.len());
    setup::check_parsed(
        file_inputs.cli,
        &outputs.parse_stats,
        ms_level_filtered,
        &outputs.all_spectra,
        spectrum_paths,
    )?;
    Ok(outputs)
}

/// The PIN's inputs: whether refinement candidates were merged, the IPF scoring (with
/// `--localize`), and pass 1's candidates with the refinement's merged in when it ran (then also
/// the merged index). The prepared search is released here.
fn pin_tables<'s>(
    cli: &Cli,
    refine_output: Option<search::refinement::RefinementOutput>,
    queues: &mut [search::TopNQueue],
    mut prepared: PreparedSearch<'_>,
    (idx, scorer): (&search::SearchIndex, &'s RankScorer),
) -> PinTables<'s> {
    let refine_merged = refine_output
        .as_ref()
        .is_some_and(|o| !o.candidates.is_empty());
    let ipf_scoring = cli.localize.is_some().then(|| prepared.ipf_scoring(scorer));
    match refine_output {
        Some(out) => {
            let merged = search::refinement::merge_into_pass1(
                queues,
                std::mem::take(&mut prepared.candidates),
                idx,
                out,
            );
            drop(prepared);
            (
                refine_merged,
                ipf_scoring,
                merged.candidates,
                Some(merged.index),
            )
        }
        None => {
            let candidates = std::mem::take(&mut prepared.candidates);
            drop(prepared);
            (refine_merged, ipf_scoring, candidates, None)
        }
    }
}

/// Spectra per chunk handed from the reader to the search threads.
const CHUNK_SIZE: usize = 5000;

/// The knockoff rows and the IPF table for the QPX output, from `refine_rescored`.
type RescoredOutputs = (
    std::collections::HashSet<String>,
    Option<std::collections::HashMap<String, output::ipf::IpfQpx>>,
);

/// The isotope-error window of the glyco search and the spectra it applies to (those with a verified
/// monoisotope).
type GlycoIsotopeOverride = Option<(std::ops::RangeInclusive<i8>, Vec<bool>)>;

#[cfg(test)]
#[path = "search/tests/snap.rs"]
mod xl_snap_tests;

#[cfg(test)]
#[path = "search/tests/filtered.rs"]
mod filtered_tsv_tests;

#[cfg(test)]
#[path = "search/tests/rtcal.rs"]
mod dia_rt_cal_tests;
