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

use std::path::{Path, PathBuf};

use model::activation::ActivationMethod;
use model::Spectrum;
use scoring_crate::RankScorer;
use search::candidates::records::index_cache_path;
use search::precursor::monoisotopic::MonoCorrection;
use search::{FileCalibration, PrecursorCalMode, PreparedSearch, TopNQueue};

use crate::cli::{Cli, Fragmentation};
use crate::glyco::run_glyco;
use crate::models::search_aa_set;
use crate::monoisotopic::{MonoDump, MonoStats};
use crate::spectra::{learn_file_calibrations, ms3_carries_reporters, ParseStats};
use calibrate::{finish_mono_correction, optimize_on_sample, FragmentErrorMeasurer};
use collect::quant_collector;
use files::{search_file, FileInputs, FileOutputs};
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

pub(crate) fn run(cli: Cli) -> Result<(), Box<dyn std::error::Error>> {
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
    let RunPaths {
        scratch_dir,
        original_spectrum_paths,
        output_pin_path,
        report_base,
        temp_dirs: _temp_dirs,
    } = prepare_run(&mut cli)?;
    let standalone = cli.glyco || cli.xl.is_some();
    let database_path: PathBuf = cli.database.clone().expect("database validated in main");

    let t_total = std::time::Instant::now();
    let db = load_search_database(&cli, &database_path)?;

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

    let bench_cap = if cli.max_spectra > 0 {
        cli.max_spectra
    } else {
        usize::MAX
    };
    let ms_level_u32 = cli.ms_level as u32;
    if (is_raw || is_d) && cli.ms_level != 2 {
        log::warn!(
            "--ms-level={} is ignored for native .raw/.d input; these formats \
             search MS2 (identification) scans only.",
            cli.ms_level
        );
    }

    let cache_dir = crate::cli::dir_or_output_dir(cli.index_cache_dir.as_deref(), &report_base);
    let index_path = index_cache_path(&idx, &params, &cache_dir);
    let mut file_cal: Vec<(f64, model::tolerance::PrecursorTolerance)> = Vec::new();
    if params.precursor_cal_mode != PrecursorCalMode::Off {
        let t_cal = std::time::Instant::now();
        let mut cal_prepared = PreparedSearch::prepare_mmap(
            &idx,
            &params,
            &scorer,
            fragment_tol_da,
            &cli.db.decoy_prefix,
            &index_path,
        )
        .map_err(|e| format!("build out-of-core candidate index: {e}"))?;
        file_cal = learn_file_calibrations(
            spectrum_paths,
            ms_level_u32,
            bench_cap,
            &params,
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

    let file_params: Vec<search::SearchParams> = file_cal
        .iter()
        .map(|&(shift, tolerance)| {
            let mut p = params.clone();
            p.precursor_mass_shift_ppm = shift;
            if params.precursor_tolerance == file_cal[0].1 {
                p.precursor_tolerance = tolerance;
            }
            p
        })
        .collect();
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
        prepared = prepared.with_per_scan_scorers(routed_activation, loader);
    }
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

    let bench_mode = cli.max_spectra > 0;

    let high_res = acquisition.high_res;
    let tmt_level = match (cli.tmt_level, cli.tmt) {
        (Some(level), _) => level,
        (None, Some(plex))
            if !high_res
                && plex
                    .check_tolerance(model::tolerance::Tolerance::Da(0.3))
                    .is_err()
                && ms3_carries_reporters(&spectrum_path, is_mzml, is_raw) =>
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
    };
    let mut quant = quant_collector(
        &cli,
        high_res,
        tmt_level,
        &scorer,
        is_mgf,
        is_d,
        spectrum_paths,
    )?;
    let predictors = if standalone {
        None
    } else {
        crate::dia::rt::Predictors::from_cli(&cli, &acquisition, param.data_type.instrument)?
    };
    if let Some(p) = &predictors {
        log::info!(target: "predict", "{}", p.describe());
    }
    let keep_peaks = cli.refine
        || cli.dia
        || cli.localize.is_some()
        || predictors.as_ref().is_some_and(|p| p.needs_peaks());
    let all_spectra: Vec<Spectrum> = Vec::new();
    let all_queues: Vec<TopNQueue> = Vec::new();
    let run_ends: Vec<usize> = Vec::with_capacity(spectrum_paths.len());
    let dia_ms1: Vec<Vec<model::scan::Ms1Scan>> = Vec::new();
    let xl_ms3_capture = cli
        .xl
        .as_ref()
        .is_some_and(|x| x.0.iter().any(|l| l.cleavage.is_some()));
    let xl_ms3: Vec<Vec<model::scan::ProductScan>> = Vec::new();

    let t_search_start = std::time::Instant::now();

    let mzml_warn_ms_level_emitted = if is_mgf && cli.ms_level != 2 {
        log::warn!(
            "--ms-level={} requested for an MGF input; MGF carries no MS \
             level (always treated as MS2). The flag has no effect on this input.",
            cli.ms_level
        );
        true
    } else {
        false
    };
    let _ = mzml_warn_ms_level_emitted;

    let xl_snapped = 0usize;
    let xl_snap_files = 0usize;
    let ms1_input = chimeric_active || mono_active;
    let parse_stats = ParseStats::default();
    let mono_table: Vec<Option<MonoCorrection>> = Vec::new();
    let mono_stats = MonoStats::default();
    let mono_dump: Option<MonoDump> = match (&cli.precursor_mono_dump, mono_active) {
        (Some(p), true) => Some(
            MonoDump::create(p, mono_params.max_shift)
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
    };

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
    let mut outputs = FileOutputs {
        all_spectra,
        all_queues,
        run_ends,
        dia_ms1,
        xl_ms3,
        xl_snapped,
        xl_snap_files,
        parse_stats,
        mono_table,
        mono_stats,
        mono_dump,
    };
    for (file_idx, input_path) in spectrum_paths.iter().enumerate() {
        if bench_mode && outputs.all_spectra.len() >= bench_cap {
            break;
        }
        search_file(
            &file_inputs,
            &mut outputs,
            &mut prepared,
            &mut measurer,
            &mut quant,
            file_idx,
            input_path,
        )?;
    }
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
    let file_calibrations: Vec<FileCalibration> = file_params
        .iter()
        .zip(&run_ends)
        .map(|(p, &end)| FileCalibration {
            end,
            shift_ppm: p.precursor_mass_shift_ppm,
            tolerance: p.precursor_tolerance,
        })
        .collect();

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

    if is_mzml || is_raw {
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
