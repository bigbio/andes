use std::fs::File;
use std::path::{Path, PathBuf};
use std::sync::mpsc::sync_channel;
use std::thread;

use crate::cli::{Cli, Fragmentation};
use crate::glyco_run::run_glyco;
use crate::memlimit::available_memory_budget;
use crate::model_select::search_aa_set;
use crate::mono::{correct_chunk, MonoDump, MonoStats};
use crate::quant_run::{QuantCollector, QuantSettings, RunScans};
use crate::report_search_progress;
use crate::spectra::{
    add_mzml_skipped, input_format_flags, merge_parse_stats, ms3_carries_reporters,
    prefix_spectrum_titles, require_spectra, run_precursor_calibration, sample_spectra,
    send_chunks, send_linked_chunks, title_prefix_for, tolerance_ppm_display, ParseStats,
};
use input::{MgfReader, Ms1Link, MzMLReader};
use model::{activation::ActivationMethod, Spectrum};
use scoring_crate::RankScorer;
use search::candidate_index::index_cache_path;
use search::fragment_index::{index_chunk_len, index_order_mass, INDEX_CHUNK_SIZE};
use search::precursor_mono::MonoCorrection;
use search::{
    apply_shift_for_mode, apply_tightened_precursor_tolerance, PrecursorCalMode, PreparedSearch,
    TopNQueue,
};

#[path = "search/resolve.rs"]
mod resolve;
use resolve::{
    check_flag_combinations, detect_inputs, load_search_database, log_resolved_parameters,
    prepare_run, resolve_scoring_model, resolve_search_params, ResolvedSearch, RunInputs, RunPaths,
    ScoringModel, SearchDatabase,
};

type Ms1Peaks = Vec<(f64, f32)>;

fn index_chunk_span_da(flag: Option<f64>) -> f64 {
    const GIB: f64 = 1024.0 * 1024.0 * 1024.0;
    match flag {
        Some(d) => d.max(1.0),
        None => {
            let budget = available_memory_budget()
                .map(|b| b.bytes as f64)
                .unwrap_or(32.0 * GIB);
            ((budget * 0.6 - 4.0 * GIB) / (0.13 * GIB)).clamp(10.0, 150.0)
        }
    }
}

fn is_gzip_path(path: &Path) -> bool {
    path.extension()
        .and_then(|e| e.to_str())
        .is_some_and(|e| e.eq_ignore_ascii_case("gz"))
}

#[allow(clippy::too_many_arguments)]
fn score_mass_ordered(
    prepared: &PreparedSearch,
    mut pending: Vec<Spectrum>,
    span_da: f64,
    all_spectra: &mut Vec<Spectrum>,
    all_queues: &mut Vec<TopNQueue>,
    t_search_start: std::time::Instant,
    keep_peaks: bool,
    mut quant: Option<&mut QuantCollector>,
    measurer: &mut FragmentErrorMeasurer,
) {
    if prepared.params.dia {
        prepared.set_dia_rt_prior(None);
    }
    measurer.measure(prepared, pending.iter(), Some(span_da));
    if prepared.params.dia {
        calibrate_dia_rt(prepared, &pending, span_da);
    }
    pending.sort_by(|a, b| index_order_mass(a).total_cmp(&index_order_mass(b)));
    log::info!(
        target: "fragment-index",
        "{} spectra sorted by precursor mass, chunks of <= {} spectra and <= {:.0} Da",
        pending.len(),
        INDEX_CHUNK_SIZE,
        span_da
    );
    let mut rest = pending;
    while !rest.is_empty() {
        let take = index_chunk_len(&rest, index_order_mass, span_da);
        let chunk: Vec<Spectrum> = rest.drain(..take).collect();
        let offset = all_spectra.len();
        let queues = prepared.run_chunk(&chunk, offset);
        all_queues.extend(queues);
        for mut spec in chunk.into_iter() {
            if let Some(q) = quant.as_deref_mut() {
                q.on_spectrum(&spec);
            }
            if !keep_peaks {
                spec.peaks = Vec::new();
            }
            all_spectra.push(spec);
        }
        report_search_progress(all_spectra.len(), t_search_start);
    }
    if prepared.params.dia {
        prepared.set_dia_rt_prior(None);
    }
}

/// Share of a DIA run's spectra searched first, with the plain vote rank, to fit the RT prior.
const DIA_RT_CAL_FRACTION: f64 = 0.02;
/// Fewest calibration spectra (all of them in a smaller run).
const DIA_RT_CAL_MIN_SPECTRA: usize = 2_000;
/// Target-decoy q-value of the calibration winners that become RT anchors.
const DIA_RT_CAL_Q: f64 = 0.01;
const SECS_PER_MIN: f64 = 60.0;

/// Fits the run's RT prior for the DIA vote cut and installs it on `prepared`: a deterministic
/// subsample of the run is searched with the plain vote rank, the top match of each sample
/// spectrum that is a target at q <= `DIA_RT_CAL_Q` becomes an anchor (andes-rt raw prediction
/// of its peptide, scan minutes), and [`search::dia_rt_rank::DiaRtPrior::fit`] maps raw
/// predictions to minutes. With too few anchors no prior is installed and the cut keeps the
/// vote rank.
fn calibrate_dia_rt(prepared: &PreparedSearch, spectra: &[Spectrum], span_da: f64) {
    let Some(model) = prepared.form_rt_model() else {
        return;
    };
    if spectra.is_empty() {
        return;
    }
    let t = std::time::Instant::now();
    let n = ((spectra.len() as f64 * DIA_RT_CAL_FRACTION).round() as usize)
        .max(DIA_RT_CAL_MIN_SPECTRA)
        .min(spectra.len());
    let mut order: Vec<(u64, usize)> = (0..spectra.len())
        .map(|i| (splitmix64(i as u64), i))
        .collect();
    if n < order.len() {
        order.select_nth_unstable(n - 1);
        order.truncate(n);
    }
    let mut sample: Vec<Spectrum> = order.iter().map(|&(_, i)| spectra[i].clone()).collect();
    sample.sort_by(|a, b| index_order_mass(a).total_cmp(&index_order_mass(b)));
    let mut queues: Vec<TopNQueue> = Vec::with_capacity(sample.len());
    let mut start = 0;
    while start < sample.len() {
        let take = index_chunk_len(&sample[start..], index_order_mass, span_da);
        queues.extend(prepared.run_chunk(&sample[start..start + take], start));
        start += take;
    }
    let mut tops: Vec<(usize, u32)> = Vec::new();
    let mut labels: Vec<search::tdc::ScoredLabel> = Vec::new();
    for (i, q) in queues.iter().enumerate() {
        if let Some(p) = q.peek_top() {
            let idx = p.primary_candidate_idx();
            tops.push((i, idx));
            labels.push(search::tdc::ScoredLabel {
                score: p.score,
                is_decoy: prepared.candidate_is_decoy(idx),
            });
        }
    }
    let q = search::tdc::qvalues(&labels);
    let anchors: Vec<(f64, f64)> = tops
        .iter()
        .zip(&labels)
        .zip(&q)
        .filter(|((_, l), q)| !l.is_decoy && **q <= DIA_RT_CAL_Q)
        .filter_map(|(((si, idx), _), _)| {
            let rt = sample[*si].rt_seconds.filter(|t| t.is_finite())? / SECS_PER_MIN;
            Some((f64::from(model(&prepared.candidate_peptide(*idx))), rt))
        })
        .collect();
    match search::dia_rt_rank::DiaRtPrior::fit(&anchors) {
        Some(prior) => {
            log::info!(
                target: "dia-rt-rank",
                "RT prior from {} of {} spectra: {} anchors at q <= {DIA_RT_CAL_Q} ({} within \
                 the trimming band), {} knots, sigma {:.3} min; the vote cut ranks by votes \
                 and RT [{:.1}s]",
                sample.len(),
                spectra.len(),
                prior.anchors(),
                prior.kept(),
                prior.knots(),
                prior.sigma_min(),
                t.elapsed().as_secs_f64()
            );
            prepared.set_dia_rt_prior(Some(prior));
        }
        None => log::info!(
            target: "dia-rt-rank",
            "no RT prior: {} anchors at q <= {DIA_RT_CAL_Q} from {} sample spectra (at least \
             {} needed); the vote cut ranks by votes [{:.1}s]",
            anchors.len(),
            sample.len(),
            search::dia_rt_rank::MIN_ANCHORS,
            t.elapsed().as_secs_f64()
        ),
    }
}

fn splitmix64(x: u64) -> u64 {
    let mut z = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

fn recalibrate_dia_ms1(prepared: &mut PreparedSearch, pending: &[Spectrum], span_da: f64) {
    const SAMPLE: usize = 20_000;
    const WINDOW_PPM: f64 = 25.0;
    let t = std::time::Instant::now();
    let mut refs: Vec<&Spectrum> = pending.iter().filter(|s| s.peaks.len() >= 10).collect();
    refs.sort_by(|a, b| index_order_mass(a).total_cmp(&index_order_mass(b)));
    let n = SAMPLE.min(refs.len());
    let lo = (refs.len() / 2).saturating_sub(n / 2);
    let sample: Vec<Spectrum> = refs[lo..lo + n].iter().map(|s| (*s).clone()).collect();
    drop(refs);
    let mut wide = prepared.params.clone();
    wide.dia_ms1_tol_ppm = WINDOW_PPM;
    let mut queues: Vec<TopNQueue> = Vec::with_capacity(sample.len());
    let mut start = 0;
    while start < sample.len() {
        let take = index_chunk_len(&sample[start..], index_order_mass, span_da);
        queues.extend(prepared.run_chunk_with_params(&sample[start..start + take], start, &wide));
        start += take;
    }
    match search::fragment_error::dia_ms1_offset_ppm(prepared, &sample, &queues, WINDOW_PPM) {
        Some((ppm, n)) => {
            prepared.shift_dia_ms1_ppm(ppm);
            log::info!(
                target: "dia",
                "MS1 m/z offset {ppm:+.2} ppm from {n} confident PSMs; survey peaks recentred ({:.1}s)",
                t.elapsed().as_secs_f64()
            );
        }
        None => log::info!(
            target: "dia",
            "MS1 m/z offset not measured (too few confident PSMs) ({:.1}s)",
            t.elapsed().as_secs_f64()
        ),
    }
}

struct FragmentErrorMeasurer {
    enabled: bool,
    done: bool,
    per_file: Vec<Option<scoring_crate::scoring::FragmentErrorModel>>,
}

impl FragmentErrorMeasurer {
    fn start_file(&mut self, prepared: &PreparedSearch) {
        self.done = false;
        self.per_file.push(None);
        if self.enabled {
            prepared.set_fragment_error(None);
        }
    }

    fn measure<'s>(
        &mut self,
        prepared: &PreparedSearch,
        spectra: impl Iterator<Item = &'s Spectrum>,
        index_span_da: Option<f64>,
    ) {
        if !self.enabled || self.done {
            return;
        }
        self.done = true;
        let t = std::time::Instant::now();
        let mut refs: Vec<&Spectrum> = spectra.filter(|s| s.peaks.len() >= 10).collect();
        refs.sort_by(|a, b| index_order_mass(a).total_cmp(&index_order_mass(b)));
        let n = search::fragment_error::SAMPLE_SPECTRA.min(refs.len());
        let lo = (refs.len() / 2).saturating_sub(n / 2);
        let sample: Vec<Spectrum> = refs[lo..lo + n].iter().map(|s| (*s).clone()).collect();
        drop(refs);
        let mut queues: Vec<TopNQueue> = Vec::with_capacity(sample.len());
        match index_span_da {
            Some(span_da) => {
                let mut start = 0;
                while start < sample.len() {
                    let take = index_chunk_len(&sample[start..], index_order_mass, span_da);
                    queues.extend(prepared.run_chunk(&sample[start..start + take], start));
                    start += take;
                }
            }
            None => queues.extend(prepared.run_chunk(&sample, 0)),
        }
        match search::fragment_error::measure_on_queues(prepared, &sample, &queues) {
            Ok(m) => {
                prepared.set_fragment_error(Some(m));
                if let Some(last) = self.per_file.last_mut() {
                    *last = Some(m);
                }
                log::info!(
                    target: "fragment-error",
                    "offset {:+.2} ppm, spread {:.2} ppm from {} ions of {} PSMs \
                     at 1% q in {} sampled spectra -> fragment peaks recentred by {:+.2} ppm, \
                     PIN feature tolerance {:.1} ppm, ppm-Gaussian width {:.2} ppm, \
                     node-score window {:?} ({:.1}s)",
                    m.offset_ppm,
                    m.sigma_ppm,
                    m.ions,
                    m.psms,
                    sample.len(),
                    prepared.scorer.fragment_recentre_ppm(),
                    m.feature_tolerance_ppm(),
                    m.gaussian_sigma_ppm(),
                    prepared.scorer.node_score_window(),
                    t.elapsed().as_secs_f64()
                );
            }
            Err(why) => log::info!(
                target: "fragment-error",
                "not measured ({why}); the PIN features keep {:?} ({:.1}s)",
                prepared.scorer.feature_match_tolerance(),
                t.elapsed().as_secs_f64()
            ),
        }
    }
}

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
        if let Some(v) =
            crate::xl_run::xl_default_variable_mods(cli.db.variable_mods.as_ref(), &x.0)
        {
            cli.db.variable_mods = Some(v);
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

    const CHUNK_SIZE: usize = 5000;

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

    if is_raw && params.precursor_cal_mode != PrecursorCalMode::Off {
        log::warn!(
            "--precursor-cal is not yet supported for Thermo .raw input; \
             proceeding without calibration."
        );
        params.precursor_cal_mode = PrecursorCalMode::Off;
    }

    let cache_dir = crate::cli::dir_or_output_dir(cli.index_cache_dir.as_deref(), &report_base);
    let index_path = index_cache_path(&idx, &params, &cache_dir);
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
        let cal_stats = run_precursor_calibration(
            &spectrum_path,
            is_mzml,
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
        params.precursor_mass_shift_ppm =
            apply_shift_for_mode(params.precursor_cal_mode, cal_stats);
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
                log::info!(
                    target: "calibration",
                    "Tightened precursor tolerance for main pass: left {:.3} ppm -> {:.3} ppm, right {:.3} ppm -> {:.3} ppm",
                    left_before.unwrap_or(0.0),
                    left_after.unwrap_or(0.0),
                    right_before.unwrap_or(0.0),
                    right_after.unwrap_or(0.0),
                );
            }
        }
    }

    if cli.optimize {
        if cli.glyco {
            return Err("--optimize is not supported with --glyco".into());
        }
        let sample = crate::spectra::sample_spectra_strided(
            &spectrum_path,
            is_mzml,
            cli.ms_level as u32,
            cli.optimize_sample,
        )
        .or_else(|| {
            sample_spectra(&spectrum_path, is_mzml, is_mgf, is_raw, cli.ms_level as u32)
                .map(|v| v.into_iter().take(cli.optimize_sample).collect())
        });
        match sample {
            Some(sample) => {
                let optimize_dir = scratch_tempdir("andes-optimize-", &scratch_dir)?;
                let rescorer = crate::optimize::SampleRescorer {
                    backend: output::PercolatorBackend::Native {
                        memory: output::MemLimit {
                            bytes: cli.rescore_memory_gb.saturating_mul(1 << 30),
                            spill_dir: optimize_dir.path().to_path_buf(),
                        },
                    },
                    args: output::percolator_args_with_threads(&[], cli.threads),
                    dir: optimize_dir.path().to_path_buf(),
                };
                crate::optimize::optimize_parameters(
                    &idx,
                    &scorer,
                    fragment_tol_da,
                    &cli.db.decoy_prefix,
                    &index_path,
                    &mut params,
                    &sample,
                    &rescorer,
                );
                let optimize_path = optimize_dir.path().to_path_buf();
                if let Err(e) = optimize_dir.close() {
                    log::warn!(
                        "optimize: could not remove {}: {e}",
                        optimize_path.display()
                    );
                }
            }
            None => {
                log::warn!("--optimize could not read a sample; keeping the resolved parameters")
            }
        }
    }

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
            .unwrap_or_else(crate::model_select::bundled_store_path);
        let run_model_id = resolved_model_id.clone();
        let run_analyser =
            crate::model_select::analyser_of(param.data_type.instrument.name()).to_string();
        let run_protocol = param.data_type.protocol;
        let high_res = acquisition.high_res;
        let index_on = params.fragment_index_top_k > 0;
        let gbdt_max_trees = gbdt_k;
        let loader: search::ScorerLoader = Box::new(move |activation: ActivationMethod| {
            let route = match crate::model_select::load_activation_route(
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
            if index_on && !search::fragment_index::covers_activation(p.data_type.activation) {
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
        prepared.set_form_rt_model(crate::predicted::andes_rt_form_model()?);
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
    if cli.glyco && cli.mbr && !cli.lfq {
        return Err(
            "--glyco --mbr needs --lfq: match-between-runs transfers the identified intact \
             glycopeptides into the other runs"
                .into(),
        );
    }
    let mut quant: Option<QuantCollector> = QuantSettings::from_cli(&cli, high_res, tmt_level)?
        .map(|mut settings| {
            settings.fragment_tol = Some(scorer.feature_match_tolerance());
            QuantCollector::new(settings)
        });
    if let Some(q) = &quant {
        let s = &q.settings;
        if is_mgf && (s.needs_ms1_series() || s.needs_ms3()) {
            return Err(
                "--lfq, --quant-dir and --tmt-level 3 need MS1/MS3 scans; MGF input carries MS2 only".into(),
            );
        }
        if is_d && s.needs_ms3() {
            return Err(
                "--tmt-level 3 reads reporter ions from MS3 scans, and a Bruker .d has none \
                 (a TDF records MS1 and PASEF MS2 frames only); quantify the reporters from \
                 the MS2 (--tmt-level 2)"
                    .into(),
            );
        }
        if cli.chimeric && s.needs_ms3() {
            return Err(
                "--tmt-level 3 cannot be combined with --chimeric: the cascade reads MS2 only"
                    .into(),
            );
        }
        if cli.glyco {
            let mut unsupported: Vec<&str> = Vec::new();
            if s.plex.is_some() {
                unsupported.push("--tmt");
            }
            if cli.quant_dir.is_some() {
                unsupported.push("--quant-dir");
            }
            if !unsupported.is_empty() {
                return Err(format!(
                    "--glyco does not support: {} (label-free quantification of intact \
                     glycopeptides is `--glyco --lfq`; the isobaric and match-between-runs \
                     stages are not implemented for glyco mode)",
                    unsupported.join(", ")
                )
                .into());
            }
        }
        if !cli.glyco && !cli.rescore && !cli.rescore_native && (s.plex.is_some() || s.lfq) {
            log::warn!(
                "--tmt/--lfq without --rescore or --rescore-native: every rank-1 target \
                 PSM is quantified; filter the tables by the q-values you compute downstream, \
                 or add --rescore so only PSMs at q<=--psm-fdr are quantified."
            );
        }
        if !cli.glyco {
            log::info!(target: "quant", "{}", s.describe());
        }
    }
    let predictors = if standalone {
        None
    } else {
        crate::predicted::Predictors::from_cli(&cli, &acquisition, param.data_type.instrument)?
    };
    if let Some(p) = &predictors {
        log::info!(target: "predict", "{}", p.describe());
    }
    let keep_peaks = cli.refine
        || cli.dia
        || cli.localize.is_some()
        || predictors.as_ref().is_some_and(|p| p.needs_peaks());
    let mut all_spectra: Vec<Spectrum> = Vec::new();
    let mut all_queues: Vec<TopNQueue> = Vec::new();
    let mut run_ends: Vec<usize> = Vec::with_capacity(spectrum_paths.len());
    let mut dia_ms1: Vec<Vec<model::scan::Ms1Scan>> = Vec::new();

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

    let mut xl_snapped = 0usize;
    let mut xl_snap_files = 0usize;
    let ms1_input = chimeric_active || mono_active;
    let mut parse_stats = ParseStats::default();
    let mut mono_table: Vec<Option<MonoCorrection>> = Vec::new();
    let mut mono_stats = MonoStats::default();
    let mut mono_dump: Option<MonoDump> = match (&cli.precursor_mono_dump, mono_active) {
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
    for (file_idx, input_path) in spectrum_paths.iter().enumerate() {
        if bench_mode && all_spectra.len() >= bench_cap {
            break;
        }
        measurer.start_file(&prepared);
        let remaining_cap = if bench_mode {
            bench_cap.saturating_sub(all_spectra.len())
        } else {
            usize::MAX
        };

        let (file_is_mzml, file_is_raw, file_is_d, _file_is_mgf) = input_format_flags(input_path);
        let file_stem = input_path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("spectrum");
        let title_prefix = title_prefix_for(file_stem);

        if spectrum_paths.len() > 1 {
            log::info!(
                target: "input",
                "spectrum file {}/{}: {}",
                file_idx + 1,
                spectrum_paths.len(),
                input_path.display()
            );
        }

        let file_ms1 = file_is_mzml || file_is_raw || file_is_d;
        let file_snap =
            xl_precursor_snap(cli.xl.is_some(), cli.ms_level, file_is_mzml, file_is_raw);
        xl_snap_files += usize::from(file_snap);
        // This file reads its MS1 scans with its MS2 scans.
        let file_linked = (ms1_input && file_ms1) || file_snap;
        if (chimeric_active || mono_active) && !file_ms1 {
            return Err(format!(
                "{} only supports mzML/.raw/.d inputs, got {}",
                if chimeric_active {
                    "--chimeric"
                } else {
                    "--precursor-mono"
                },
                input_path.display()
            )
            .into());
        }

        let spectra_before = all_spectra.len();
        let mut run_scans = RunScans::default();
        let file_stats = if file_linked {
            let (tx, rx) = sync_channel::<(Vec<Spectrum>, Ms1Link)>(2);
            let spectrum_path = input_path.clone();
            let cap = remaining_cap;
            let mslevel = 2;
            let capture_ms1 = quant.as_ref().is_some_and(|q| q.settings.needs_ms1());
            let capture_context = quant.is_some();
            let capture_noise = quant.as_ref().is_some_and(|q| q.settings.plex.is_some());
            let index_ordered = chimeric_active && params.fragment_index_top_k > 0 && !standalone;
            let parallel_read = index_ordered
                && file_is_mzml
                && capture_ms1
                && cap == usize::MAX
                && !is_gzip_path(input_path);
            let threads = cli.threads;
            type Ms1Parsed = (usize, Vec<String>, RunScans);
            let parser_handle = thread::spawn(move || -> Result<Ms1Parsed, String> {
                let linked = if parallel_read {
                    input::mzml::read_ms1_linked_parallel(
                        &spectrum_path,
                        threads,
                        capture_context,
                        capture_noise,
                    )
                    .map_err(|e| format!("read mzML: {e}"))?
                } else {
                    None
                };
                if let Some(read) = linked {
                    send_linked_chunks(
                        read.spectra,
                        &read.ms2_to_ms1,
                        &read.run_ms1,
                        CHUNK_SIZE,
                        &tx,
                    );
                    let scans = RunScans {
                        ms1: read.run_ms1,
                        ms3: Vec::new(),
                        contexts: read.contexts,
                    };
                    Ok((0, Vec::new(), scans))
                } else if file_is_mzml {
                    let mut reader = MzMLReader::new(
                        input::open_buf_maybe_gz(&spectrum_path)
                            .map_err(|e| format!("open mzML: {e}"))?,
                    )
                    .with_ms_level_range(mslevel, mslevel)
                    .with_ms1_capture(true)
                    .with_run_ms1_capture(capture_ms1)
                    .with_scan_context_capture(capture_context)
                    .with_noise_capture(capture_noise);
                    let (errc, errs) =
                        reader.read_with_ms1_chunked(CHUNK_SIZE, cap, |chunk, link| {
                            let _ = tx.send((chunk, link));
                        });
                    let scans = RunScans {
                        ms1: reader.take_run_ms1(),
                        ms3: Vec::new(),
                        contexts: reader.take_scan_contexts(),
                    };
                    Ok((errc, errs, scans))
                } else if file_is_d {
                    #[cfg(feature = "timstof")]
                    {
                        let mut reader = input::TimsTofReader::open(&spectrum_path)
                            .map_err(|e| format!("open Bruker .d: {e}"))?
                            .with_dia_mobility_slices(cli.dia_mobility_slices)
                            .with_run_ms1_capture(capture_ms1)
                            .with_ms1_mobility(true)
                            .with_scan_context_capture(capture_context);
                        let (errc, errs) =
                            reader.read_with_ms1_chunked(CHUNK_SIZE, cap, |chunk, link| {
                                let _ = tx.send((chunk, link));
                            });
                        let scans = RunScans {
                            ms1: reader.take_run_ms1(),
                            ms3: Vec::new(),
                            contexts: reader.take_scan_contexts(),
                        };
                        Ok((errc, errs, scans))
                    }
                    #[cfg(not(feature = "timstof"))]
                    {
                        Err("this andes build has no Bruker .d (timsTOF) support; \
                             rebuild with `--features timstof`."
                            .to_string())
                    }
                } else {
                    #[cfg(feature = "thermo")]
                    {
                        let mut reader = input::ThermoRawReader::open(&spectrum_path)
                            .map_err(|e| format!("open Thermo .raw: {e}"))?
                            .with_run_ms1_capture(capture_ms1)
                            .with_scan_context_capture(capture_context)
                            .with_noise_capture(capture_noise);
                        let (errc, errs) =
                            reader.read_with_ms1_chunked(CHUNK_SIZE, cap, |chunk, link| {
                                let _ = tx.send((chunk, link));
                            });
                        let scans = RunScans {
                            ms1: reader.take_run_ms1(),
                            ms3: Vec::new(),
                            contexts: reader.take_scan_contexts(),
                        };
                        Ok((errc, errs, scans))
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
            let mut pending: Vec<(Spectrum, Option<Ms1Peaks>)> = Vec::new();
            for (mut chunk_spectra, chunk_link) in rx {
                if let Some(prefix) = &title_prefix {
                    prefix_spectrum_titles(&mut chunk_spectra, prefix);
                }
                if file_snap {
                    for (i, spec) in chunk_spectra.iter_mut().enumerate() {
                        let ms1 = chunk_link
                            .ms2_to_ms1
                            .get(i)
                            .copied()
                            .flatten()
                            .and_then(|m| chunk_link.ms1_peaks.get(m));
                        if let Some(mz) = ms1.and_then(|p| {
                            search::precursor_mono::snap_precursor_mz(
                                p,
                                spec.precursor_mz,
                                search::xl_search::PRECURSOR_SNAP_PPM,
                                search::xl_search::PRECURSOR_SNAP_RATIO,
                            )
                        }) {
                            spec.precursor_mz = mz;
                            xl_snapped += 1;
                        }
                    }
                }
                if mono_active {
                    correct_chunk(
                        &mut chunk_spectra,
                        &chunk_link,
                        &mono_params,
                        &mut mono_table,
                        &mut mono_stats,
                        mono_dump.as_mut(),
                    )
                    .map_err(|e| format!("write --precursor-mono-dump: {e}"))?;
                }
                if index_ordered {
                    file_offset += chunk_spectra.len();
                    ms1_linked += chunk_link.ms1_peaks.len();
                    for (i, spec) in chunk_spectra.into_iter().enumerate() {
                        let ms1 = chunk_link
                            .ms2_to_ms1
                            .get(i)
                            .copied()
                            .flatten()
                            .and_then(|m| chunk_link.ms1_peaks.get(m))
                            .map(|peaks| {
                                let (lo, hi) =
                                    search::match_engine::coisolation_ms1_bounds(&spec, &params);
                                let a = peaks.partition_point(|p| p.0 < lo);
                                let b = peaks.partition_point(|p| p.0 <= hi);
                                peaks[a..b.max(a)].to_vec()
                            });
                        pending.push((spec, ms1));
                    }
                    continue;
                }
                let offset = all_spectra.len();
                measurer.measure(&prepared, chunk_spectra.iter(), None);
                if chimeric_active {
                    let mut queues = prepared.run_chunk(&chunk_spectra, offset);
                    search::match_engine::run_pass2_coisolation(
                        &prepared,
                        &chunk_spectra,
                        &mut queues,
                        &params,
                        &chunk_link,
                        offset,
                    );
                    all_queues.extend(queues);
                } else if !standalone {
                    let queues = prepared.run_chunk(&chunk_spectra, offset);
                    all_queues.extend(queues);
                }
                file_offset += chunk_spectra.len();
                ms1_linked += chunk_link.ms1_peaks.len();
                for mut spec in chunk_spectra.into_iter() {
                    if let Some(q) = quant.as_mut() {
                        q.on_spectrum(&spec);
                    }
                    if !keep_peaks && !standalone {
                        spec.peaks = Vec::new();
                    }
                    all_spectra.push(spec);
                }
                report_search_progress(all_spectra.len(), t_search_start);
            }
            if index_ordered {
                let span_da = index_chunk_span_da(cli.fragment_index_slice_da);
                measurer.measure(&prepared, pending.iter().map(|e| &e.0), Some(span_da));
                pending.sort_by(|a, b| index_order_mass(&a.0).total_cmp(&index_order_mass(&b.0)));
                log::info!(
                    target: "fragment-index",
                    "{} spectra sorted by precursor mass, chunks of <= {} spectra and <= {:.0} Da",
                    pending.len(),
                    INDEX_CHUNK_SIZE,
                    span_da
                );
                let base = all_spectra.len();
                let mut scored: Vec<Spectrum> = Vec::with_capacity(pending.len());
                let mut queues: Vec<search::psm::TopNQueue> = Vec::with_capacity(pending.len());
                let mut link = Ms1Link::default();
                let mut rest = pending;
                while !rest.is_empty() {
                    let take = index_chunk_len(&rest, |e| index_order_mass(&e.0), span_da);
                    let mut chunk: Vec<Spectrum> = Vec::with_capacity(take);
                    for (spec, ms1) in rest.drain(..take) {
                        chunk.push(spec);
                        let ms1_idx = ms1.map(|peaks| {
                            link.ms1_peaks.push(peaks);
                            link.ms1_peaks.len() - 1
                        });
                        link.ms2_to_ms1.push(ms1_idx);
                    }
                    queues.extend(prepared.run_chunk(&chunk, base + scored.len()));
                    scored.extend(chunk);
                    report_search_progress(base + scored.len(), t_search_start);
                }
                search::match_engine::run_pass2_coisolation(
                    &prepared,
                    &scored,
                    &mut queues,
                    &params,
                    &link,
                    base,
                );
                all_queues.extend(queues);
                for mut spec in scored.into_iter() {
                    if let Some(q) = quant.as_mut() {
                        q.on_spectrum(&spec);
                    }
                    if !keep_peaks {
                        spec.peaks = Vec::new();
                    }
                    all_spectra.push(spec);
                }
            }
            let (err_count, first_errors, scans) = parser_handle
                .join()
                .map_err(|_| "MS1-linked parser thread panicked".to_string())??;
            run_scans = scans;
            log::info!(
                target: if chimeric_active {
                    "chimeric"
                } else if mono_active {
                    "precursor-mono"
                } else {
                    "xl"
                },
                "streamed {} MS2 spectra ({} MS1 scans linked) from {}",
                file_offset,
                ms1_linked,
                input_path.display()
            );
            ParseStats {
                error_count: err_count,
                first_errors,
            }
        } else if params.fragment_index_top_k > 0
            && !standalone
            && !cli.dia
            && file_is_mzml
            && !is_gzip_path(input_path)
            && !quant
                .as_ref()
                .is_some_and(|q| q.settings.needs_ms1() || q.settings.needs_ms3())
        {
            let (mut spectra, errors) =
                input::read_spectra_parallel(input_path, ms_level_u32, ms_level_u32, cli.threads)
                    .map_err(|e| format!("read mzML {}: {e}", input_path.display()))?;
            spectra.truncate(remaining_cap);
            if let Some(prefix) = &title_prefix {
                prefix_spectrum_titles(&mut spectra, prefix);
            }
            score_mass_ordered(
                &prepared,
                spectra,
                index_chunk_span_da(cli.fragment_index_slice_da),
                &mut all_spectra,
                &mut all_queues,
                t_search_start,
                keep_peaks,
                quant.as_mut(),
                &mut measurer,
            );
            ParseStats {
                error_count: errors.len(),
                first_errors: errors.into_iter().take(3).collect(),
            }
        } else {
            let (tx, rx) = sync_channel::<Vec<Spectrum>>(2);
            let spectrum_path = input_path.clone();
            let capture_ms1 = cli.dia || quant.as_ref().is_some_and(|q| q.settings.needs_ms1());
            let capture_ms3 = quant.as_ref().is_some_and(|q| q.settings.needs_ms3());
            let capture_context = quant.is_some();
            let capture_noise = quant.as_ref().is_some_and(|q| q.settings.plex.is_some());
            let parser_handle = thread::spawn(
                move || -> Result<(ParseStats, RunScans), Box<dyn std::error::Error + Send + Sync>> {
                    if file_is_mzml {
                        let mut reader = MzMLReader::new(
                            input::open_buf_maybe_gz(&spectrum_path)
                                .map_err(|e| format!("open mzML: {e}"))?,
                        )
                        .with_ms_level_range(ms_level_u32, ms_level_u32)
                        .with_run_ms1_capture(capture_ms1)
                        .with_product_scan_capture(capture_ms3)
                        .with_scan_context_capture(capture_context)
                        .with_noise_capture(capture_noise);
                        let mut stats =
                            send_chunks(reader.by_ref(), CHUNK_SIZE, remaining_cap, tx);
                        add_mzml_skipped(&mut stats, &reader);
                        let scans = RunScans {
                            ms1: reader.take_run_ms1(),
                            ms3: reader.take_product_scans(),
                            contexts: reader.take_scan_contexts(),
                        };
                        Ok((stats, scans))
                    } else if file_is_raw {
                        #[cfg(feature = "thermo")]
                        {
                            let mut reader = input::ThermoRawReader::open(&spectrum_path)
                                .map_err(|e| format!("open Thermo .raw: {e}"))?
                                .with_ms_level(Some(2))
                                .with_run_ms1_capture(capture_ms1)
                                .with_product_scan_capture(capture_ms3)
                                .with_scan_context_capture(capture_context)
                                .with_noise_capture(capture_noise);
                            let stats =
                                send_chunks(reader.by_ref(), CHUNK_SIZE, remaining_cap, tx);
                            let scans = RunScans {
                                ms1: reader.take_run_ms1(),
                                ms3: reader.take_product_scans(),
                                contexts: reader.take_scan_contexts(),
                            };
                            Ok((stats, scans))
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
                            let mut reader = input::TimsTofReader::open(&spectrum_path)
                                .map_err(|e| format!("open Bruker .d: {e}"))?
                                .with_dia_mobility_slices(cli.dia_mobility_slices)
                                .with_run_ms1_capture(capture_ms1)
                                .with_ms1_mobility(true)
                                .with_scan_context_capture(capture_context);
                            let stats =
                                send_chunks(reader.by_ref(), CHUNK_SIZE, remaining_cap, tx);
                            let scans = RunScans {
                                ms1: reader.take_run_ms1(),
                                ms3: Vec::new(),
                                contexts: reader.take_scan_contexts(),
                            };
                            Ok((stats, scans))
                        }
                        #[cfg(not(feature = "timstof"))]
                        {
                            Err("this andes build has no Bruker .d (timsTOF) support; \
                                 rebuild with `--features timstof`. mzML/MGF inputs work \
                                 without it."
                                .into())
                        }
                    } else {
                        let reader = MgfReader::new(
                            input::open_buf_maybe_gz(&spectrum_path)
                                .map_err(|e| format!("open MGF: {e}"))?,
                        );
                        Ok((
                            send_chunks(reader, CHUNK_SIZE, remaining_cap, tx),
                            RunScans::default(),
                        ))
                    }
                },
            );

            let mass_ordered = (params.fragment_index_top_k > 0 || cli.dia) && !standalone;
            let mut pending: Vec<Spectrum> = Vec::new();
            for mut chunk in rx {
                if chunk.is_empty() {
                    continue;
                }
                if let Some(prefix) = &title_prefix {
                    prefix_spectrum_titles(&mut chunk, prefix);
                }
                if mass_ordered {
                    pending.extend(chunk);
                    continue;
                }
                if !standalone {
                    measurer.measure(&prepared, chunk.iter(), None);
                    let offset = all_spectra.len();
                    let queues = prepared.run_chunk(&chunk, offset);
                    all_queues.extend(queues);
                }
                for mut spec in chunk.into_iter() {
                    if let Some(q) = quant.as_mut() {
                        q.on_spectrum(&spec);
                    }
                    if !keep_peaks && !standalone {
                        spec.peaks = Vec::new();
                    }
                    all_spectra.push(spec);
                }
                report_search_progress(all_spectra.len(), t_search_start);
            }
            let stats = match parser_handle.join() {
                Ok(Ok((stats, scans))) => {
                    run_scans = scans;
                    stats
                }
                Ok(Err(e)) => return Err(format!("parser thread error: {e}").into()),
                Err(_) => return Err("parser thread panicked".into()),
            };
            if cli.dia && mass_ordered {
                if run_scans.ms1.is_empty() {
                    return Err(format!(
                        "--dia needs MS1 survey scans, and {} has none",
                        input_path.display()
                    )
                    .into());
                }
                for s in &mut pending {
                    s.precursor_charge = None;
                }
                prepared.set_dia_ms1(if quant.is_some() {
                    run_scans.ms1.clone()
                } else {
                    std::mem::take(&mut run_scans.ms1)
                });
                if file_is_d {
                    recalibrate_dia_ms1(
                        &mut prepared,
                        &pending,
                        index_chunk_span_da(cli.fragment_index_slice_da),
                    );
                }
            }
            if mass_ordered {
                score_mass_ordered(
                    &prepared,
                    pending,
                    index_chunk_span_da(cli.fragment_index_slice_da),
                    &mut all_spectra,
                    &mut all_queues,
                    t_search_start,
                    keep_peaks,
                    quant.as_mut(),
                    &mut measurer,
                );
            }
            stats
        };
        require_spectra(
            input_path,
            if file_linked { 2 } else { cli.ms_level },
            all_spectra.len() - spectra_before,
        )
        .map_err(|e| crate::spectra::with_first_read_error(e, &file_stats.first_errors))?;
        merge_parse_stats(&mut parse_stats, file_stats);
        if cli.dia {
            let gated = prepared.take_dia_ms1();
            dia_ms1.push(if !gated.is_empty() {
                gated
            } else if quant.is_some() {
                run_scans.ms1.clone()
            } else {
                std::mem::take(&mut run_scans.ms1)
            });
        }
        if let Some(q) = quant.as_mut() {
            q.finish_file(
                input_path,
                title_prefix.as_deref(),
                run_scans,
                all_spectra.len(),
            )?;
        }
        run_ends.push(all_spectra.len());
    }
    run_ends.resize(spectrum_paths.len(), all_spectra.len());

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
    let mono_active = mono_active && {
        if mono_stats.fitted == 0 {
            log::warn!(
                "--precursor-mono auto: no MS2 could be fitted ({} MS2 seen, none with \
                 a linked MS1 and a known charge); precursors are left as recorded and the \
                 glyco PIN Mono columns are 0.",
                mono_stats.seen
            );
            if let Some(d) = mono_dump.take() {
                d.finish()
                    .map_err(|e| format!("flush --precursor-mono-dump: {e}"))?;
            }
            mono_table.clear();
        }
        mono_stats.fitted > 0
    };
    if mono_active {
        let hist: Vec<String> = mono_stats
            .by_shift
            .iter()
            .enumerate()
            .skip(1)
            .filter(|(_, &n)| n > 0)
            .map(|(k, n)| format!("{k}:{n}"))
            .collect();
        log::info!(
            target: "precursor-mono",
            "{} of {} MS2 fitted (linked MS1 + known charge), {} shifted \
             down [shift:count {}]",
            mono_stats.fitted,
            mono_stats.seen,
            mono_stats.shifted(),
            if hist.is_empty() {
                "-".to_string()
            } else {
                hist.join(" ")
            }
        );
        if let Some(d) = mono_dump.take() {
            d.finish()
                .map_err(|e| format!("flush --precursor-mono-dump: {e}"))?;
        }
        debug_assert_eq!(mono_table.len(), all_spectra.len());
    }
    let (window, narrowed) = crate::mono::coupled_isotope_window(
        mono_active,
        cli.isotope_error,
        params.isotope_error_range.clone(),
    );
    let glyco_isotope_override = if narrowed {
        let mask: Vec<bool> = mono_table
            .iter()
            .map(|c| {
                c.as_ref()
                    .is_some_and(|c| c.verified_within(mono_params.min_fit, *window.end()))
            })
            .collect();
        log::info!(
            target: "precursor-mono",
            "isotope-error window {}..={} -> {}..={} for the {} of {} MS2 \
             with a verified monoisotope (the others keep {}..={}; a verified monoisotope \
             needs no +2 step, which only admits the Hex+Fuc/NeuGc composition degeneracy — \
             pass --isotope-error to override)",
            params.isotope_error_range.start(),
            params.isotope_error_range.end(),
            window.start(),
            window.end(),
            mask.iter().filter(|&&v| v).count(),
            mono_stats.seen,
            params.isotope_error_range.start(),
            params.isotope_error_range.end()
        );
        Some((window, mask))
    } else {
        None
    };

    let spectra = all_spectra;
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
                search::xl_search::PRECURSOR_SNAP_PPM
            );
        }
        crate::xl_run::run_xl(
            &cli,
            &spectra,
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
                1,
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
    log::info!(
        target: "total",
        "done in {:.2} s",
        t_total.elapsed().as_secs_f64()
    );
    let pooled = cli.dia && params.top_n_psms_per_spectrum > cli.top_n;
    let mut dia_pass2: Option<output::dia_pass2_table::Pass2Table> = None;
    if let Some(p) = &predictors {
        let kept = p.add_to_pin(&output_pin_path, &spectra, &run_ends, pooled)?;
        if let Some(predicted) = kept.filter(|_| pooled) {
            let pool = crate::dia_llr::add_to_pin(
                &output_pin_path,
                &spectra,
                &run_ends,
                &dia_ms1,
                &predicted,
                p.tolerance(),
            )?;
            dia_pass2 = Some(crate::dia_pass2::run(
                &output_pin_path,
                &spectra,
                &run_ends,
                &dia_ms1,
                &pool,
                &predicted,
            )?);
        }
    }
    let mut rerank_dropped = std::collections::HashSet::new();
    if pooled {
        rerank_dropped = crate::dia::rerank_pin(&output_pin_path, cli.top_n as usize)?;
    }
    if cli.dia {
        crate::dia::add_to_pin(&output_pin_path, &spectra, &run_ends, dia_ms1)?;
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
    let use_native = cli.rescore_native;
    let use_percolator = cli.rescore || cli.rescore_native;
    if !use_percolator && (cli.psm_fdr.is_some() || cli.pep.is_some()) {
        let which = match (cli.psm_fdr.is_some(), cli.pep.is_some()) {
            (true, true) => "--psm-fdr and --pep",
            (true, false) => "--psm-fdr",
            _ => "--pep",
        };
        log::warn!(
            "{which} set without --rescore or --rescore-native, so no rescoring \
             runs and no filtered .tsv is written. andes does not compute FDR on its \
             own; the PIN is written for Percolator. Add --rescore to run Percolator \
             in-process."
        );
    }
    let mut mbr_decoys: std::collections::HashMap<String, output::PercolatorPsm> =
        std::collections::HashMap::new();
    let mut mbr_source = String::new();
    let mut rescore_decoys: Option<std::collections::HashSet<String>> = None;
    let mut rescore_decoy_psms: std::collections::HashMap<String, output::PercolatorPsm> =
        std::collections::HashMap::new();
    let has_site_claims = pin_index.iter_target_proteins().any(|p| {
        p.annotation
            .as_deref()
            .is_some_and(|a| !a.site_claims.is_empty())
    });
    let mut row_scores: Option<output::RowScores> = None;
    let mut dia_replaced: std::collections::HashSet<String> = std::collections::HashSet::new();
    let mut dia_precursor_ids: std::collections::HashSet<String> = std::collections::HashSet::new();
    let mut rescore_map: Option<std::collections::HashMap<String, output::PercolatorPsm>> =
        if use_percolator {
            let backend = if use_native {
                crate::cli::ensure_scratch_dir(&scratch_dir)?;
                output::PercolatorBackend::Native {
                    memory: output::MemLimit {
                        bytes: cli.rescore_memory_gb.saturating_mul(1 << 30),
                        spill_dir: scratch_dir.clone(),
                    },
                }
            } else {
                output::resolve_backend(
                    cli.percolator_bin.as_deref(),
                    cli.percolator_docker,
                    &cli.percolator_image,
                )
                .map_err(|e| format!("{e}"))?
            };
            log::info!(target: "rescore", "Percolator backend = {}", backend.describe());
            let mut user_args: Vec<String> = cli
                .percolator_args
                .split_whitespace()
                .map(|s| s.to_string())
                .collect();
            if cli.dia {
                crate::dia::default_train_fdr(&mut user_args);
            }
            let extra = output::percolator_args_with_threads(&user_args, cli.threads);
            let t_resc = std::time::Instant::now();
            let (mut map, mut decoys) =
                output::run_percolator_with_decoys(&backend, &output_pin_path, &extra)
                    .map_err(|e| format!("{e}"))?;
            if cli.dia
                && crate::dia::refine_rt(&output_pin_path, &spectra, &|id: &str| {
                    map.get(id).is_some_and(|p| p.q_value <= 0.01)
                })?
            {
                (map, decoys) =
                    output::run_percolator_with_decoys(&backend, &output_pin_path, &extra)
                        .map_err(|e| format!("{e}"))?;
            }
            let path = output::decoy_results_path(&output_pin_path);
            let text =
                std::fs::read_to_string(&path).map_err(|e| format!("{}: {e}", path.display()))?;
            rescore_decoy_psms =
                output::parse_psm_results(&text).map_err(|e| format!("{}: {e}", path.display()))?;
            let mut pass1_q = 0;
            let mut pass1_best = Vec::new();
            if cli.dia {
                let ((t, d, n), best) = crate::dia::precursor_q(
                    &output_pin_path,
                    &spectra,
                    &mut map,
                    &rescore_decoy_psms,
                )?;
                pass1_q = n;
                pass1_best = best;
                log::info!(
                    target: "rescore",
                    "Precursor-level q (DIA): {t} target and {d} decoy precursors, {n} at q <= 0.01"
                );
            }
            let mut results_pin = output_pin_path.clone();
            if let Some(t) = dia_pass2.take() {
                let seed = crate::dia_score::percolator_seed(&user_args);
                let rescore = |pin: &Path| -> Result<crate::dia_score::Rescored, String> {
                    let (targets, decoy_ids) =
                        output::run_percolator_with_decoys(&backend, pin, &extra)
                            .map_err(|e| format!("{e}"))?;
                    let path = output::decoy_results_path(pin);
                    let text = std::fs::read_to_string(&path)
                        .map_err(|e| format!("{}: {e}", path.display()))?;
                    let decoys = output::parse_psm_results(&text)
                        .map_err(|e| format!("{}: {e}", path.display()))?;
                    Ok(crate::dia_score::Rescored {
                        targets,
                        decoy_ids,
                        decoys,
                    })
                };
                match crate::dia_score::rescore_precursors(
                    &output_pin_path,
                    &t,
                    seed,
                    &pass1_best,
                    rescore,
                ) {
                    Ok(r) => {
                        log::info!(
                            target: "dia",
                            "precursor PIN used: {} target precursors at q <= 0.01 (pass 1: \
                             {pass1_q}), {} degenerate learned-score folds",
                            r.n_q,
                            r.ids.degenerate_folds
                        );
                        map = r.rescored.targets;
                        decoys = r.rescored.decoy_ids;
                        rescore_decoy_psms = r.rescored.decoys;
                        dia_replaced = r.replaced;
                        dia_precursor_ids = r.ids.targets.into_iter().chain(r.ids.decoys).collect();
                        results_pin = r.pin;
                    }
                    Err(why) => log::info!(
                        target: "dia",
                        "pass-1 results kept, the precursor PIN is not used: {why}"
                    ),
                }
            }
            rescore_decoys = Some(decoys);
            log::info!(
                target: "rescore",
                "Percolator produced {} target PSMs [{:.2}s]",
                map.len(),
                t_resc.elapsed().as_secs_f64()
            );
            if has_site_claims {
                row_scores = Some(
                    output::score_all_rows(&backend, &results_pin, &extra).map_err(|e| {
                        format!("scoring every PIN row for the site posteriors: {e}")
                    })?,
                );
            }
            if cli.mbr {
                mbr_source = format!("the Percolator results of {}", results_pin.display());
            }
            Some(map)
        } else {
            None
        };

    let mut knockoff_rows: std::collections::HashSet<String> = std::collections::HashSet::new();
    let mut ipf_qpx: Option<std::collections::HashMap<String, output::ipf::IpfQpx>> = None;
    if let Some(ref mut map) = rescore_map {
        if pin_index.has_alternative_entries {
            let shared = output::share_isobaric_references(
                map,
                &pin_index,
                fdr_threshold.max(proteoform_fdr),
            );
            if shared > 0 {
                log::info!(
                    target: "rescore",
                    "{shared} variant-class PSMs share their peptide with a reference \
                     protein once I and L are taken as equal; they now name it"
                );
            }
        }
        let counts = output::apply_class_fdr(map, &rescore_decoy_psms, &pin_index);
        if counts.len() > 1 {
            let line: Vec<String> = counts
                .iter()
                .map(|c| {
                    format!(
                        "{} {}/{} {} pi0 {:.2}",
                        c.class,
                        c.targets,
                        c.decoys,
                        c.method.name(),
                        c.pi0
                    )
                })
                .collect();
            log::info!(
                target: "rescore",
                "FDR within each class (targets/decoys, method, pi0): {}",
                line.join(", ")
            );
        }
        let outcome =
            output::apply_knockoff_fdr(map, &pin_index, proteoform_fdr, cli.db.knockoffs > 0);
        knockoff_rows = outcome.knockoff_rows;
        let untested: usize = outcome
            .counts
            .iter()
            .map(|c| c.real_before - c.tested_before)
            .sum();
        if untested > 0 {
            log::warn!(
                "{untested} variant PSMs at q <= {proteoform_fdr} have no knockoff \
                 (the enzyme or the variant kind leaves no position for one) and are not \
                 accepted; --knockoffs 0 leaves variants to the class FDR alone"
            );
        }
        if !outcome.counts.is_empty() {
            let line: Vec<String> = outcome
                .counts
                .iter()
                .map(|c| {
                    format!(
                        "{} {} -> {} ({} tested; knockoffs {:.1}, false share {:.2})",
                        c.stratum,
                        c.real_before,
                        c.real_after,
                        c.tested_before,
                        c.knockoff_weighted,
                        c.false_share()
                    )
                })
                .collect();
            log::info!(
                target: "rescore",
                "variant PSMs at q<={proteoform_fdr} per substitution class, before -> \
                 after knockoff q-values: {}",
                line.join(", ")
            );
        }
        let posteriors = row_scores.as_ref().map(|(target_scores, decoy_scores)| {
            let rows: std::collections::HashMap<String, output::site_identity::RowEvidence> =
                queues
                    .iter()
                    .enumerate()
                    .filter(|(_, q)| q.iter_psms().any(|p| p.site_evidence.is_some()))
                    .flat_map(|(i, q)| {
                        output::emitted_row_evidence(&spectra[i], q, &pin_candidates)
                            .into_iter()
                            .map(move |(id, isotope_offset, evidence)| {
                                (
                                    id,
                                    output::site_identity::RowEvidence {
                                        spectrum: i as u32,
                                        isotope_offset,
                                        evidence,
                                    },
                                )
                            })
                    })
                    .collect();
            output::site_identity::row_posteriors(target_scores, decoy_scores, &rows)
        });
        let sites = output::site_claims(
            map,
            &pin_index,
            fdr_threshold,
            proteoform_fdr,
            &restricted_site_masses,
            posteriors.as_ref(),
            &fixed_terminal_mods,
        );
        knockoff_rows.extend(sites.knockoff_rows.iter().cloned());
        if !sites.counts.is_empty() {
            let line: Vec<String> = sites
                .counts
                .iter()
                .map(|c| {
                    format!(
                        "{} {} -> {} ({} tested; knockoffs {:.1})",
                        c.stratum, c.claims, c.accepted, c.tested, c.knockoff_weighted
                    )
                })
                .collect();
            log::info!(
                target: "rescore",
                "site-claim PSMs at q<={fdr_threshold} per modification and residue, \
                 before -> after site q<={proteoform_fdr}: {}",
                line.join(", ")
            );
            let stem = report_base
                .file_stem()
                .map(|s| s.to_string_lossy().into_owned())
                .unwrap_or_else(|| "andes".to_string());
            let dir = report_base.parent().unwrap_or_else(|| Path::new("."));
            let path = dir.join(format!("{stem}.sites.tsv"));
            match output::write_sites_tsv(&path, &sites.rows) {
                Ok(()) => log::info!(
                    target: "output",
                    "Wrote site claims: {} ({} sites)",
                    path.display(),
                    sites.rows.len()
                ),
                Err(e) => log::warn!(target: "output", "could not write {}: {e}", path.display()),
            }
        }
        if let Some(scoring) = ipf_scoring.as_ref() {
            use rayon::prelude::*;
            let t_ipf = std::time::Instant::now();
            let map_ref: &std::collections::HashMap<String, output::percolator::PercolatorPsm> =
                map;
            let cutoffs = output::Cutoffs::new(map_ref, &pin_index, fdr_threshold, proteoform_fdr);
            let searched_with = scoring.fragment_error();
            let mut evidence: std::collections::HashMap<String, search::ipf::IpfEvidence> =
                std::collections::HashMap::new();
            let mut lo = 0usize;
            for (f, &hi) in run_ends.iter().enumerate() {
                scoring.set_fragment_error(measurer.per_file.get(f).copied().flatten());
                let part: Vec<(String, search::ipf::IpfEvidence)> = queues[lo..hi]
                    .par_iter()
                    .enumerate()
                    .map(|(j, q)| (lo + j, q))
                    .flat_map_iter(|(i, q)| {
                        let rows: Vec<(String, u32, u8)> =
                            output::emitted_row_candidates(&spectra[i], q, &pin_candidates)
                                .into_iter()
                                .filter(|(id, c, _)| {
                                    !pin_candidates[*c as usize].is_decoy
                                        && !knockoff_rows.contains(id)
                                        && map_ref
                                            .get(id)
                                            .is_some_and(|p| p.q_value <= cutoffs.of(id))
                                })
                                .collect();
                        let peptides: Vec<(&model::peptide::Peptide, u8)> = rows
                            .iter()
                            .map(|(_, c, z)| (&pin_candidates[*c as usize].peptide, *z))
                            .collect();
                        let found = if peptides.is_empty() {
                            Vec::new()
                        } else {
                            scoring.evidence(&spectra[i], &peptides)
                        };
                        rows.into_iter()
                            .zip(found)
                            .filter_map(|((id, _, _), e)| e.map(|e| (id, e)))
                            .collect::<Vec<_>>()
                    })
                    .collect();
                evidence.extend(part);
                lo = hi;
            }
            scoring.set_fragment_error(searched_with);
            let (rows, model) = output::ipf::localize_psms(map_ref, &evidence, cli.ipf_debug);
            let stem = report_base
                .file_stem()
                .map(|s| s.to_string_lossy().into_owned())
                .unwrap_or_else(|| "andes".to_string());
            let dir = report_base.parent().unwrap_or_else(|| Path::new("."));
            let path = dir.join(format!("{stem}.ipf.tsv"));
            output::ipf::write_ipf_tsv(&path, &rows)
                .map_err(|e| format!("could not write {}: {e}", path.display()))?;
            if rows.is_empty() {
                log::warn!(
                    "--localize ipf: no identified target PSM has a variable \
                     modification with more than one possible site; {} is empty",
                    path.display()
                );
            }
            let truncated = rows.iter().filter(|r| r.truncated).count();
            log::info!(
                target: "output",
                "Wrote IPF localization: {} ({} PSMs, {} with placements limited to the \
                 {} nearest the search's; ion PEP by intensity quartile {:?}, unmatched \
                 {:.3}) [{:.2}s]",
                path.display(),
                rows.len(),
                truncated,
                search::ipf::MAX_PLACEMENTS,
                model.matched.map(|p| (p * 1000.0).round() / 1000.0),
                model.unmatched,
                t_ipf.elapsed().as_secs_f64()
            );
            if cli.output_parquet.is_some() {
                ipf_qpx = Some(output::ipf::qpx_map(&rows));
            }
        }
    }
    let protein_groups = rescore_map.as_ref().map(|map| {
        let cutoffs = output::Cutoffs::new(map, &pin_index, fdr_threshold, proteoform_fdr);
        let inf = output::infer_protein_groups(map, &rescore_decoy_psms, &pin_index, |id| {
            map[id].q_value <= cutoffs.of(id) && !knockoff_rows.contains(id)
        });
        let stem = report_base
            .file_stem()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| "andes".to_string());
        let dir = report_base.parent().unwrap_or_else(|| Path::new("."));
        let path = dir.join(format!("{stem}.protein_groups.tsv"));
        let at = |q: f64| inf.groups.iter().filter(|g| g.q_value <= q).count();
        match output::write_protein_groups_tsv(&path, &inf) {
            Ok(()) => log::info!(
                target: "output",
                "Wrote protein groups: {} ({} groups; {} at protein q<=0.01, {} at q<=0.05)",
                path.display(),
                inf.groups.len(),
                at(0.01),
                at(0.05)
            ),
            Err(e) => log::warn!(target: "output", "could not write {}: {e}", path.display()),
        }
        inf
    });
    if cli.mbr {
        mbr_decoys = std::mem::take(&mut rescore_decoy_psms);
    }
    let cutoffs = rescore_map
        .as_ref()
        .map(|m| output::Cutoffs::new(m, &pin_index, fdr_threshold, proteoform_fdr));

    let mut identified: Option<(Vec<&output::PercolatorPsm>, String)> = None;
    let mut open_search_summary: Option<String> = None;
    if let (Some(map), Some(cutoffs)) = (rescore_map.as_ref(), cutoffs.as_ref()) {
        let fdr = fdr_threshold;
        let stem = report_base
            .file_stem()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| "andes".to_string());
        let dir = report_base.parent().unwrap_or_else(|| Path::new("."));
        let q_tag = format!("{fdr}").replace('.', "p");
        let filtered_path = dir.join(format!("{stem}.q{q_tag}.tsv"));
        let pep_cap = cli.pep;
        let mut accepted: Vec<&output::PercolatorPsm> = map
            .iter()
            .filter(|(id, p)| {
                p.q_value <= cutoffs.of(id)
                    && pep_cap.is_none_or(|t| p.pep <= t)
                    && !knockoff_rows.contains(*id)
                    && cli.protein_fdr.is_none_or(|t| {
                        protein_groups
                            .as_ref()
                            .and_then(|inf| inf.group_of(p))
                            .is_some_and(|g| g.q_value <= t)
                    })
            })
            .map(|(_, p)| p)
            .collect();
        accepted.sort_by(|a, b| a.report_cmp(b));
        let mut pep_note = pep_cap
            .map(|t| format!(" and PEP<={t}"))
            .unwrap_or_default();
        if pin_index.has_alternative_entries && proteoform_fdr != fdr {
            pep_note.push_str(&format!(", proteoform claims at q<={proteoform_fdr}"));
        }
        if let Some(t) = cli.protein_fdr {
            pep_note.push_str(&format!(", protein groups at q<={t}"));
        }
        match write_filtered_tsv(&filtered_path, &accepted, protein_groups.as_ref()) {
            Ok(()) => log::info!(
                target: "output",
                "Wrote filtered TSV: {} ({} PSMs at q<={}{})",
                filtered_path.display(),
                accepted.len(),
                fdr,
                pep_note
            ),
            Err(e) => {
                log::warn!(target: "output", "could not write {}: {e}", filtered_path.display())
            }
        }
        if cli.open_search {
            open_search_summary =
                write_open_offsets(&output_pin_path, dir, &stem, &accepted, &rescore_decoy_psms);
        }
        identified = Some((accepted, pep_note));
    }

    let run_stats = output::RunStatistics::compute(
        &queues,
        &pin_candidates,
        &params,
        &scorer.feature_match_tolerance(),
        &scorer.pin_match_tolerance(),
        &scorer.node_score_window(),
    );
    let mut id_summary = identified.as_ref().map(|(accepted, cut_note)| {
        let wanted: std::collections::HashSet<&str> =
            accepted.iter().map(|p| p.psm_id.as_str()).collect();
        let facts = output::psm_facts(
            &spectra,
            &queues,
            &pin_candidates,
            &params,
            &run_ends,
            &wanted,
        );
        let linked: Vec<output::PsmFacts> = accepted
            .iter()
            .filter_map(|p| facts.get(p.psm_id.as_str()).cloned())
            .collect();
        let decoys = if cli.mbr {
            &mbr_decoys
        } else {
            &rescore_decoy_psms
        };
        let decoys_at_cut = output::Cutoffs::new(decoys, &pin_index, fdr_threshold, proteoform_fdr)
            .count_passing(decoys, cli.pep);
        let run_names: Vec<String> = spectrum_paths
            .iter()
            .map(|p| crate::quant_run::run_stem(p))
            .collect();
        let fixed = output::FixedModSpec::declared(&params);
        output::IdentificationSummary::compute(
            &output::SpectraTally::of(&spectra, &run_ends),
            &run_names,
            &linked,
            accepted.len() - linked.len(),
            &output::SummaryContext {
                rescorer: if use_native {
                    "Percolator (native engine)"
                } else {
                    "Percolator"
                },
                cut: &format!("q ≤ {fdr_threshold}{cut_note}"),
                decoys_at_cut,
                protein_groups: protein_groups.as_ref().map(|inf| {
                    (
                        fdr_threshold,
                        inf.groups
                            .iter()
                            .filter(|g| g.q_value <= fdr_threshold)
                            .count(),
                    )
                }),
                fixed: &fixed,
                calibration_shift_ppm: params.precursor_mass_shift_ppm,
            },
        )
    });
    if let Some(summary) = id_summary.as_mut() {
        summary.open_search = open_search_summary.take();
    }
    log::info!(target: "summary", "{}", run_stats.render_with(id_summary.as_ref()));
    let stats_path = {
        let stem = report_base
            .file_stem()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| "andes".to_string());
        let dir = report_base.parent().unwrap_or_else(|| Path::new("."));
        dir.join(format!("{stem}.statistics.log"))
    };
    match output::write_statistics_log(&run_stats, id_summary.as_ref(), &stats_path) {
        Ok(()) => log::info!(target: "output", "Wrote statistics: {}", stats_path.display()),
        Err(e) => log::warn!(target: "output", "could not write {}: {e}", stats_path.display()),
    }

    knockoff_rows.extend(rerank_dropped);
    if !dia_precursor_ids.is_empty() {
        knockoff_rows.extend(dia_replaced);
        knockoff_rows.retain(|id| !dia_precursor_ids.contains(id));
    }
    if let Some(ref parquet_dir) = cli.output_parquet {
        let run_id = spectrum_path
            .file_stem()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| "andes_run".to_string());
        let primary_paths: Vec<String> = original_spectrum_paths
            .iter()
            .map(|p| p.display().to_string())
            .collect();
        match output::write_qpx(
            parquet_dir,
            &spectra,
            &queues,
            &pin_candidates,
            &params,
            &pin_index,
            &scorer.feature_match_tolerance(),
            &run_id,
            &primary_paths,
            &run_ends,
            rescore_map.as_ref(),
            rescore_decoys.as_ref(),
            &knockoff_rows,
            ipf_qpx.as_ref(),
        ) {
            Ok(()) => log::info!(target: "output", "Wrote QPX bundle: {}", parquet_dir.display()),
            Err(e) => {
                return Err(format!(
                    "failed to write requested --output-parquet {}: {e}",
                    parquet_dir.display()
                )
                .into());
            }
        }
    }

    let omitted_psms = match &quant {
        Some(_) => output::omitted_psm_rows(
            &spectra,
            &queues,
            &pin_candidates,
            rescore_map.as_ref(),
            rescore_decoys.as_ref(),
            &knockoff_rows,
        )
        .map_err(|e| format!("quant: {e}"))?,
        None => std::collections::HashSet::new(),
    };
    if let Some(ref mut map) = rescore_map {
        map.retain(|id, _| !knockoff_rows.contains(id));
    }

    let qpx_groups = protein_groups
        .as_ref()
        .map(|inf| crate::quant_cmd::ProteinGroups {
            inf,
            max_q: cli.protein_fdr.unwrap_or(0.01),
        });
    let tmt_mbr_input = (cli.mbr && cli.tmt.is_some()).then(|| crate::quant_run::TmtMbrInput {
        decoys: &mbr_decoys,
        source: &mbr_source,
        library_fdr: cli.mbr_donor_fdr(),
        accept_q: cli.quant_transfer_fdr(),
        wrong_log2: cli.mbr_wrong_log2,
    });
    let quant_summary = match &quant {
        Some(q) => Some(crate::quant_run::run_quant(
            q,
            &spectra,
            &queues,
            &pin_candidates,
            &pin_index,
            rescore_map.as_ref(),
            &omitted_psms,
            cutoffs.as_ref(),
            &report_base,
            cli.output_parquet.as_deref(),
            qpx_groups.as_ref(),
            tmt_mbr_input.as_ref(),
            (cli.mbr && cli.label.is_some()).then_some(&mbr_decoys),
        )?),
        None => None,
    };

    if let (Some(map), Some(cutoffs), true) = (
        rescore_map.as_ref(),
        cutoffs.as_ref(),
        pin_index.has_alternative_entries,
    ) {
        let mut rows = output::proteoform_table(map, &pin_index, cutoffs, quant_summary.as_ref());
        if let (Some(t), Some(inf)) = (cli.protein_fdr, protein_groups.as_ref()) {
            let fq = inf.family_q();
            rows.retain(|r| fq.get(r.family.as_str()).is_some_and(|&q| q <= t));
        }
        let stem = report_base
            .file_stem()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| "andes".to_string());
        let dir = report_base.parent().unwrap_or_else(|| Path::new("."));
        let path = dir.join(format!("{stem}.proteoforms.tsv"));
        let columns = quant_summary
            .as_ref()
            .map(|q| q.columns.clone())
            .unwrap_or_default();
        match output::write_proteoform_tsv(&path, &rows, &columns) {
            Ok(()) => log::info!(
                target: "output",
                "Wrote proteoform table: {} ({} family members, {} families)",
                path.display(),
                rows.len(),
                rows.iter()
                    .map(|r| r.family.as_str())
                    .collect::<std::collections::BTreeSet<_>>()
                    .len()
            ),
            Err(e) => log::warn!(target: "output", "could not write {}: {e}", path.display()),
        }
    }

    let mbr_job = if cli.mbr && cli.tmt.is_none() && cli.label.is_none() {
        let (Some(dir), Some(map)) = (cli.quant_dir.clone(), rescore_map.as_ref()) else {
            return Err("--mbr needs --quant-dir and rescoring".into());
        };
        let stems: Vec<String> = spectrum_paths
            .iter()
            .map(|p| crate::quant_run::run_stem(p))
            .collect();
        let scores = crate::quant_cmd::rescore_scores(map, &mbr_decoys, &mbr_source)?;
        let ep = crate::quant_cmd::evidence_params(cli.lfq_tol)?;
        let (donor_q, donor_pep) = quant
            .as_ref()
            .map(|q| (q.settings.psm_fdr, q.settings.pep_cap))
            .ok_or("--mbr needs --lfq")?;
        Some((dir, stems, scores, ep, donor_q, donor_pep))
    } else {
        None
    };
    drop(mbr_decoys);

    if bench_mode {
        log::info!(target: "bench", "skipping TSV write.");
    } else if let Some(ref tsv_path) = cli.output_tsv {
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
        output::write_tsv(
            tsv_path,
            &spectra,
            &queues,
            &pin_candidates,
            &params,
            &pin_index,
            &spec_file_name,
            is_mgf,
        )?;
        log::info!(target: "output", "Wrote TSV: {}", tsv_path.display());
    }

    let Some((dir, stems, scores, ep, donor_q, donor_pep)) = mbr_job else {
        release_at_exit(pin_candidates);
        return Ok(());
    };
    drop(spectra);
    drop(queues);
    drop(pin_candidates);
    drop(pin_index);
    drop(quant);
    drop(rescore_map);
    drop(target_db);
    crate::quant_cmd::run_in_process(
        &dir,
        &stems,
        scores,
        &ep,
        cli.mbr_donor_fdr(),
        donor_q,
        donor_pep,
        cli.quant_transfer_fdr(),
        &report_base.with_extension(""),
        cli.output_parquet.as_deref(),
        qpx_groups.as_ref(),
    )?;
    Ok(())
}

pub(crate) struct ClaimDatabase {
    pub(crate) loaded: input::LoadedDatabase,
    pub(crate) annotations: Option<input::AnnotationLoadStats>,
    pub(crate) site_knockoffs: input::SiteKnockoffStats,
}

pub(crate) fn load_claim_database(
    path: &Path,
    annotations: &[PathBuf],
    opts: &input::UniProtOptions,
) -> Result<ClaimDatabase, Box<dyn std::error::Error>> {
    let mut loaded = input::load_protein_database(path, opts)?;
    let annotations = if annotations.is_empty() {
        None
    } else {
        Some(input::load_annotations(
            annotations,
            &mut loaded.db,
            opts.knockoffs,
        )?)
    };
    let site_knockoffs = opts
        .knockoffs
        .map(|ko| input::add_site_knockoffs(&mut loaded.db, ko))
        .unwrap_or_default();
    Ok(ClaimDatabase {
        loaded,
        annotations,
        site_knockoffs,
    })
}

fn release_at_exit(candidates: Vec<search::candidate_gen::Candidate>) {
    std::mem::forget(candidates);
}

fn write_open_offsets(
    pin: &Path,
    dir: &Path,
    stem: &str,
    accepted: &[&output::PercolatorPsm],
    decoys: &std::collections::HashMap<String, output::PercolatorPsm>,
) -> Option<String> {
    let report = match output::open_offsets::report_from_results(
        pin,
        accepted,
        decoys,
        model::unimod::bundled(),
    ) {
        Ok(r) => r,
        Err(e) => {
            log::warn!("open-search offset report not written: {e}");
            return None;
        }
    };
    let peaks_path = dir.join(format!("{stem}.open_offsets.parquet"));
    let psms_path = dir.join(format!("{stem}.open_psms.parquet"));
    let written = output::open_offsets::write_peaks_table(&peaks_path, &report)
        .and_then(|()| output::open_offsets::write_psms_table(&psms_path, &report));
    if let Err(e) = written {
        let mut left = Vec::new();
        for stale in [&peaks_path, &psms_path] {
            if let Err(r) = std::fs::remove_file(stale) {
                if r.kind() != std::io::ErrorKind::NotFound {
                    log::warn!(target: "output", "could not remove {}: {r}", stale.display());
                    left.push(stale.display().to_string());
                }
            }
        }
        if left.is_empty() {
            log::warn!(
                "open-search report not written ({e}); neither {} nor {} is left",
                peaks_path.display(),
                psms_path.display()
            );
        } else {
            log::warn!(
                "open-search report not written ({e}); incomplete file left: {}",
                left.join(", ")
            );
        }
        return None;
    }
    log::info!(
        target: "output",
        "Wrote open-search report: {} ({} offset peaks) and {} ({} PSMs)",
        peaks_path.display(),
        report.peaks.len(),
        psms_path.display(),
        report.assignments.len()
    );
    Some(output::open_offsets::summary_block(&report))
}

pub(crate) fn write_filtered_tsv(
    path: &Path,
    psms: &[&output::PercolatorPsm],
    groups: Option<&output::Inference>,
) -> std::io::Result<()> {
    use std::io::Write;
    output::write_atomic(path, |path| {
        let mut w = std::io::BufWriter::new(File::create(path)?);
        writeln!(
            w,
            "PSMId\tq-value\tposterior_error_probability\tpeptide\tproteinIds\tprotein_group\tprotein_q"
        )?;
        for p in psms {
            let (group, q) = groups
                .and_then(|g| g.group_of(p))
                .map_or(("NA".to_string(), "NA".to_string()), |g| {
                    (g.id.to_string(), g.q_value.to_string())
                });
            writeln!(
                w,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                p.psm_id, p.q_value, p.pep, p.peptide, p.proteins, group, q
            )?;
        }
        w.flush()
    })
}

/// Whether a crosslink search moves this file's precursor m/z values to the dominant MS1 peak
/// (`search::xl_search::PRECURSOR_SNAP_PPM`): an MS2 search of a file that carries its MS1
/// scans (mzML, Thermo `.raw`). Decided per file, so the order of the inputs does not matter.
pub(crate) fn xl_precursor_snap(
    xl: bool,
    ms_level: u8,
    file_is_mzml: bool,
    file_is_raw: bool,
) -> bool {
    xl && ms_level == 2 && (file_is_mzml || file_is_raw)
}

#[cfg(test)]
mod xl_snap_tests {
    use super::xl_precursor_snap;

    #[test]
    fn the_precursor_snap_is_decided_per_file_for_ms2_crosslink_searches() {
        // An MGF next to an mzML: each file decides alone, whatever the input order.
        let files = [(false, false), (true, false), (false, true)];
        let got: Vec<bool> = files
            .iter()
            .map(|&(m, r)| xl_precursor_snap(true, 2, m, r))
            .collect();
        assert_eq!(got, [false, true, true]);
        let reversed: Vec<bool> = files
            .iter()
            .rev()
            .map(|&(m, r)| xl_precursor_snap(true, 2, m, r))
            .collect();
        assert_eq!(reversed, [true, true, false]);
        // An MS3 search keeps its own reader; no crosslink search, no snap.
        assert!(!xl_precursor_snap(true, 3, true, false));
        assert!(!xl_precursor_snap(false, 2, true, false));
    }
}

#[cfg(test)]
mod filtered_tsv_tests {
    use super::*;

    #[test]
    fn filtered_tsv_is_written_through_a_partial_file() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("run.q0p01.tsv");
        let psm = output::PercolatorPsm {
            pooled_q_value: None,
            psm_id: "run/scan=1_1_1".into(),
            q_value: 0.001,
            pep: 0.01,
            score: Some(2.0),
            peptide: "K.PEPTIDEK.A".into(),
            proteins: "P1".into(),
        };
        write_filtered_tsv(&path, &[&psm], None).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        assert!(text.starts_with("PSMId\tq-value\t"));
        assert!(text.contains("run/scan=1_1_1\t0.001\t0.01\tK.PEPTIDEK.A\tP1\tNA\tNA"));
        assert!(!output::ms1_table::partial_path(&path).exists());
    }
}
