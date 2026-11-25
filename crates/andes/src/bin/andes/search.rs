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
    let ScoringModel {
        param,
        scorer,
        routed_activation,
        routed_instrument,
        resolved_model_id,
        gbdt_k,
        frag_tol_override,
    } = model;
    let RunInputs {
        is_mzml,
        is_raw,
        is_d,
        is_mgf,
        instrument_was_detected,
        acquisition,
        route_protocol,
        ..
    } = inputs;

    let t_phase = std::time::Instant::now();

    log::info!(target: "run", "Using {} worker threads", cli.threads);

    let fragment_tol_da = 0.5_f64;

    let bench_cap = setup::bench_cap(&cli, is_raw || is_d);
    let ms_level_u32 = cli.ms_level as u32;

    let mut trace = crate::report::build::Trace {
        clock,
        hasher,
        sources,
        command_line: crate::report::command_line(),
        models: vec![("scoring".to_string(), resolved_model_id.clone())],
        decoy_strategy,
        decoy_key,
    };
    trace.clock.stage("calibration");
    let cache_dir = crate::cli::dir_or_output_dir(cli.index_cache_dir.as_deref(), &report_base);
    let index_path = index_cache_path(&idx, &params, &cache_dir);
    let file_cal = setup::calibration_prepass(
        &cli,
        (&idx, &scorer, fragment_tol_da, &index_path),
        &mut params,
        (spectrum_paths, ms_level_u32, bench_cap),
    )?;

    if cli.optimize {
        optimize_on_sample(
            &cli,
            &spectrum_path,
            (is_mzml, is_raw, is_mgf),
            &scratch_dir,
            (&idx, &scorer, fragment_tol_da, &index_path),
            &mut params,
        )?;
    }

    let file_params = setup::file_params(&params, &file_cal);
    trace.clock.stage("index");
    let t_index = std::time::Instant::now();
    let mut prepared = PreparedSearch::prepare_mmap(
        &idx,
        &params,
        &scorer,
        fragment_tol_da,
        &cli.db.decoy_prefix,
        &index_path,
    )
    .map_err(|e| format!("build out-of-core candidate index: {e}"))?;
    prepared = setup::route_per_scan_scorers(
        prepared,
        setup::RouteInputs {
            cli: &cli,
            standalone,
            instrument_was_detected,
            resolved_model_id: &resolved_model_id,
            param: &param,
            acquisition: &acquisition,
            params: &params,
            gbdt_k,
            routed_instrument,
            route_protocol,
            search_enzyme,
            frag_tol_override,
            routed_activation,
        },
    );
    setup::finish_index_open(&cli, &mut prepared, t_index)?;

    let bench_mode = cli.max_spectra > 0;

    let high_res = acquisition.high_res;
    let tmt_level = setup::tmt_level_for(&cli, high_res, &spectrum_path, is_mzml, is_raw);
    let mut quant = quant_collector(
        &cli,
        high_res,
        tmt_level,
        &scorer,
        is_mgf,
        is_d,
        spectrum_paths,
    )?;
    let (predictors, keep_peaks) =
        setup::predictors(&cli, standalone, &acquisition, param.data_type.instrument)?;
    let xl_ms3_capture = cli
        .xl
        .as_ref()
        .is_some_and(|x| x.0.iter().any(|l| l.cleavage.is_some()));

    let t_search_start = std::time::Instant::now();

    setup::warn_mgf_ms_level(&cli, is_mgf);

    let ms1_input = chimeric_active || mono_active;

    let mut measurer = FragmentErrorMeasurer {
        enabled: acquisition.high_res && frag_tol_override.is_none() && !standalone,
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
    trace.clock.stage("search");
    setup::search_files(
        &file_inputs,
        &mut outputs,
        &mut prepared,
        &mut measurer,
        &mut quant,
        (bench_mode, bench_cap),
    )?;
    let FileOutputs {
        all_spectra,
        all_queues,
        mut run_ends,
        dia_ms1,
        xl_ms3,
        xl_snapped,
        xl_snap_files,
        parse_stats,
        mono_table,
        mono_stats,
        mono_dump,
    } = outputs;
    run_ends.resize(spectrum_paths.len(), all_spectra.len());
    prepared.params = &params;
    let file_calibrations = setup::file_calibrations(&file_params, &run_ends);

    setup::check_parsed(
        &cli,
        &parse_stats,
        is_mzml || is_raw,
        &all_spectra,
        spectrum_paths,
    )?;

    trace.clock.stage(if cli.refine { "refine" } else { "post-search" });
    let Some((spectra, mut queues, refine_output)) = post_search(
        &PostSearch {
            t_phase,
            xl_snapped,
            cli: &cli,
            acquisition: &acquisition,
            decoy_strategy,
            fragment_tol_da,
            mono_active,
            mono_params: &mono_params,
            routed_activation,
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
            scorer: &scorer,
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
    let refine_merged = refine_output
        .as_ref()
        .is_some_and(|o| !o.candidates.is_empty());
    let ipf_scoring = cli
        .localize
        .is_some()
        .then(|| prepared.ipf_scoring(&scorer));
    let pin_candidates;
    let pin_index;
    if let Some(out) = refine_output {
        let merged = search::refinement::merge_into_pass1(
            &mut queues,
            std::mem::take(&mut prepared.candidates),
            &idx,
            out,
        );
        pin_candidates = merged.candidates;
        drop(prepared);
        pin_index = merged.index;
    } else {
        pin_candidates = std::mem::take(&mut prepared.candidates);
        drop(prepared);
        pin_index = idx;
    }

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
        scorer: &scorer,
        is_mgf,
        t_phase,
        refine_merged,
        ipf_scoring,
        pin_candidates,
        pin_index,
        trace,
    })
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
