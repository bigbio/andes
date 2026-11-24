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
        mut mono_table,
        mono_stats,
        mut mono_dump,
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

    let search_elapsed = t_search_start.elapsed();
    log::info!(
        target: "stream_search",
        "Loaded and scored {} spectra from {} in chunks of {}, done in {:.2} s",
        all_spectra.len(),
        spectrum_path.display(),
        CHUNK_SIZE,
        t_phase.elapsed().as_secs_f64()
    );
    log::info!(
        target: "chunk_candidate_prep",
        "done in {:.2} s (part of stream_search: fragment-index build or record expansion, summed over chunks)",
        prepared.chunk_prep_seconds()
    );
    for (activation, own_model, scans) in prepared.per_scan_routes() {
        log::info!(
            target: "routing",
            "{scans} {} scan scorings with {} (run routed {})",
            activation.name(),
            if own_model {
                "the model of their own activation"
            } else {
                "the run's model"
            },
            routed_activation.name()
        );
    }
    if bench_mode {
        log::info!(target: "bench", "capped at {} spectra", cli.max_spectra);
    }

    prepared.sync_materialized_candidates();

    let mono_columns = mono_active;
    let (mono_active, glyco_isotope_override) = finish_mono_correction(
        &cli,
        &params,
        &mono_params,
        mono_active,
        &mono_stats,
        &mut mono_dump,
        &mut mono_table,
    )?;
    if mono_active {
        debug_assert_eq!(mono_table.len(), all_spectra.len());
    }

    let mut spectra = all_spectra;
    let mut queues = all_queues;

    if standalone {
        log::info!(
            "Standard peptide search skipped in --{} mode (PSMs unused; the driver \
             derives its own candidates); {} spectra loaded ({:.2}s)",
            if cli.glyco { "glyco" } else { "xl" },
            spectra.len(),
            search_elapsed.as_secs_f64()
        );
    } else {
        let non_empty = queues.iter().filter(|q| !q.is_empty()).count();
        log::info!(
            "Search complete: {non_empty} / {} spectra have PSMs (search wall: {:.2}s)",
            spectra.len(),
            search_elapsed.as_secs_f64()
        );
    }
    if cli.dia {
        crate::memory::release_free_memory();
    }
    crate::memory::mark_stage("search");
    if cli.open_search {
        let (attempts, wins) = prepared.open_search_stats();
        log::info!(
            target: "open-search",
            "{attempts} modification-aware re-scores ({wins} improved the node \
             score)"
        );
    }

    if cli.xl.is_some() {
        if xl_snap_files > 0 {
            log::info!(
                target: "xl",
                "{xl_snapped} precursor m/z values moved to the dominant MS1 peak within {} ppm",
                search::crosslinks::PRECURSOR_SNAP_PPM
            );
        }
        crate::xl::run_xl(
            &cli,
            &spectra,
            &xl_ms3,
            &prepared,
            &report_base,
            spectrum_paths,
            &run_ends,
            t_total,
        )?;
        return Ok(());
    }

    if cli.glyco {
        let unfitted: Vec<Option<MonoCorrection>>;
        let mono_for_pin = if mono_active {
            Some(mono_table.as_slice())
        } else if mono_columns {
            unfitted = vec![None; spectra.len()];
            Some(unfitted.as_slice())
        } else {
            None
        };
        run_glyco(
            &cli,
            &spectra,
            &prepared,
            &params,
            &idx,
            &output_pin_path,
            spectrum_paths,
            acquisition.high_res,
            mono_for_pin,
            glyco_isotope_override,
            quant.as_ref(),
            t_total,
        )?;
        return Ok(());
    }

    let refine_output = if cli.refine {
        let base_cfg = match &cli.refine_config {
            Some(p) => search::RefineConfig::from_yaml_str(&std::fs::read_to_string(p)?)
                .map_err(|e| format!("parsing --refine-config {}: {e}", p.display()))?,
            None => search::RefineConfig::default_tier(),
        };
        let cfg = base_cfg;

        let high_res = acquisition.high_res;

        let high_res_only = cfg.high_res_only;
        if high_res_only && !high_res {
            log::warn!("refine is high-res-only and the data is low-res; skipping refinement.");
            None
        } else {
            let target_db = idx.target_db();
            let t_refine = std::time::Instant::now();
            let out = search::refinement::run_refinement(
                &queues,
                &spectra,
                &prepared.candidates,
                &target_db,
                &params,
                &scorer,
                &cfg,
                0.01,
                high_res,
                fragment_tol_da,
                &cli.db.decoy_prefix,
                decoy_strategy,
                cli.decoy_seed,
                &file_calibrations,
            );
            log::info!(
                target: "refinement",
                "done in {:.2} s",
                t_refine.elapsed().as_secs_f64()
            );
            out
        }
    } else {
        None
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

    output::populate_rt_features_per_run(&spectra, &mut queues, &pin_candidates, &run_ends);

    output::populate_provenance_features(&mut queues, &pin_candidates, &pin_index);

    output::write_pin(
        &output_pin_path,
        &spectra,
        &queues,
        &pin_candidates,
        &params,
        &pin_index,
    )?;
    log::info!(
        target: "output",
        "Wrote PIN: {} in {:.2} s",
        output_pin_path.display(),
        t_phase.elapsed().as_secs_f64()
    );
    crate::memory::mark_stage("PIN written");
    let run_stats = output::RunStatistics::compute(
        &queues,
        &pin_candidates,
        &params,
        &scorer.feature_match_tolerance(),
        &scorer.pin_match_tolerance(),
        &scorer.node_score_window(),
    );
    let has_site_claims = pin_index.iter_target_proteins().any(|p| {
        p.annotation
            .as_deref()
            .is_some_and(|a| !a.site_claims.is_empty())
    });
    let fact_rows = (cli.dia
        && !has_site_claims
        && ipf_scoring.is_none()
        && quant.is_none()
        && cli.output_parquet.is_none()
        && cli.output_tsv.is_none())
    .then(|| output::fact_rows(&spectra, &queues, &pin_candidates));
    if fact_rows.is_some() {
        queues = Vec::new();
        crate::memory::release_free_memory();
    }
    log::info!(
        target: "total",
        "done in {:.2} s",
        t_total.elapsed().as_secs_f64()
    );
    let pooled = cli.dia && params.top_n_psms_per_spectrum > cli.top_n;
    let n_dia_ms1: usize = dia_ms1.iter().map(Vec::len).sum();
    let dia_ms1 = if cli.dia {
        crate::dia::ms1_indexes(dia_ms1)
    } else {
        Vec::new()
    };
    let mut dia_pass2: Option<output::dia::pass2::Pass2Table> = None;
    if let Some(p) = &predictors {
        let kept = p.add_to_pin(&output_pin_path, &spectra, &run_ends, pooled)?;
        crate::memory::mark_stage("predictions");
        if let Some((predicted, table)) = kept.filter(|_| pooled) {
            let pool = crate::dia::llr::add_to_pin(
                &output_pin_path,
                &spectra,
                &run_ends,
                &dia_ms1,
                table,
                &predicted,
                p.tolerance(),
            )?;
            crate::memory::mark_stage("pool scores");
            let nominees_out = cli
                .quant_dir
                .as_deref()
                .filter(|_| run_ends.len() == 1)
                .map(|d| {
                    crate::dia::transfer::nominee_path(
                        d,
                        &crate::quant::run::run_stem(&spectrum_path),
                    )
                });
            dia_pass2 = Some(crate::dia::pass2::run(
                &output_pin_path,
                &spectra,
                &run_ends,
                &dia_ms1,
                pool,
                predicted,
                nominees_out
                    .as_deref()
                    .map(|path| crate::dia::transfer::NomineesOut {
                        path,
                        spectrum: spectrum_path.as_path(),
                        mobility_slices: cli.dia_mobility_slices,
                    }),
            )?);
            crate::memory::mark_stage("pass 2");
        }
    }
    if cli.dia {
        crate::memory::release_free_memory();
    }
    crate::memory::mark_stage("pool released");
    let mut rerank_dropped = std::collections::HashSet::new();
    if pooled {
        rerank_dropped = crate::dia::rerank_pin(
            &output_pin_path,
            cli.top_n as usize,
            cli.output_parquet.is_some() || quant.is_some(),
        )?;
        crate::memory::release_free_memory();
        crate::memory::mark_stage("per-scan cut");
    }
    if cli.dia {
        crate::dia::add_to_pin(&output_pin_path, &spectra, &run_ends, &dia_ms1, n_dia_ms1)?;
        crate::memory::mark_stage("co-elution");
    }
    drop(dia_ms1);
    if cli.dia {
        crate::memory::release_free_memory();
    }
    if cli.refine {
        if refine_merged {
            log::info!("Refinement PSMs merged into unified PIN (Pass-1 ⊕ Pass-2).");
        } else {
            log::info!("Refinement produced no Pass-2 PSMs; unified PIN contains Pass-1 only.");
        }
    }

    let fdr_threshold = cli.psm_fdr();
    let proteoform_fdr = cli.proteoform_fdr();
    let use_native = !cli.external_percolator();
    let use_percolator = cli.rescore;
    if !use_percolator && (cli.psm_fdr.is_some() || cli.pep.is_some()) {
        let which = match (cli.psm_fdr.is_some(), cli.pep.is_some()) {
            (true, true) => "--psm-fdr and --pep",
            (true, false) => "--psm-fdr",
            _ => "--pep",
        };
        log::warn!(
            "{which} set without --rescore, so no rescoring runs and no filtered .tsv \
             is written. andes does not compute FDR on its own; the PIN is written for \
             Percolator. Add --rescore to rescore with the native Percolator engine."
        );
    }
    let mbr_decoys: std::collections::HashMap<String, output::PercolatorPsm> =
        std::collections::HashMap::new();
    let rescoring = rescoring::rescore_pin(
        &rescoring::RescoreInputs {
            cli: &cli,
            output_pin_path: &output_pin_path,
            spectrum_path: &spectrum_path,
            scratch_dir: &scratch_dir,
            use_native,
            use_percolator,
            free_peaks: fact_rows.is_some(),
            has_site_claims,
        },
        &mut spectra,
        &mut dia_pass2,
    )?;
    let rescoring::Rescoring {
        rescore_map,
        mbr_source,
        rescore_decoys,
        rescore_decoy_psms,
        row_scores,
        dia_replaced,
        dia_precursor_ids,
    } = rescoring;
    finish::after_rescoring(finish::AfterRescoring {
        cli,
        spectrum_path,
        bench_mode,
        quant,
        run_ends,
        measurer,
        file_calibrations,
        spectra,
        queues,
        target_db,
        ipf_scoring,
        run_stats,
        fact_rows,
        dia_pass2,
        rerank_dropped,
        fdr_threshold,
        proteoform_fdr,
        use_native,
        mbr_decoys,
        original_spectrum_paths,
        output_pin_path,
        report_base,
        params,
        restricted_site_masses,
        fixed_terminal_mods,
        scorer: &scorer,
        is_mgf,
        rescore_map,
        mbr_source,
        rescore_decoys,
        rescore_decoy_psms,
        row_scores,
        dia_replaced,
        dia_precursor_ids,
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
