use std::fs::File;
use std::path::{Path, PathBuf};
use std::sync::mpsc::sync_channel;
use std::sync::Arc;
use std::thread;

use crate::acquisition::Acquisition;
use crate::cli::{
    Cli, EnzymeSpecificity, EthcdActivationFlag, FragmentIndexFlag, Fragmentation,
    PrecursorMonoFlag, Protocol, ScoreFlag,
};
use crate::glyco_run::run_glyco;
use crate::memlimit::available_memory_budget;
use crate::model_select::{
    cli_fragment_tol_override, cli_fragmentation_to_activation_opt, label_from_mods,
    load_param_from_store, parse_enzymes, resolve_metadataless_selection, routing_protocol,
    search_aa_set, warn_if_universal_protease_combo,
};
use crate::mono::{correct_chunk, MonoDump, MonoStats};
use crate::quant_run::{QuantCollector, QuantSettings, RunScans};
use crate::spectra::{
    add_mzml_skipped, compare_with_first_input, detect_acquisition_metadata,
    detect_isobaric_sampled, input_format_flags, merge_parse_stats, ms3_carries_reporters,
    prefix_spectrum_titles, require_spectra, run_precursor_calibration, sample_ms2_meta,
    sample_spectra, send_chunks, send_linked_chunks, title_prefix_for, tolerance_ppm_display,
    ParseStats,
};
use crate::{arg_present, report_search_progress, EXPLICIT_MISSED_CLEAVAGES};
use input::{MgfReader, Ms1Link, MzMLReader};
use model::{
    activation::ActivationMethod, AminoAcidSetBuilder, InstrumentType, PrecursorTolerance, Spectrum,
};
use scoring_crate::RankScorer;
use search::candidate_index::index_cache_path;
use search::fragment_index::{index_chunk_len, index_order_mass, INDEX_CHUNK_SIZE};
use search::precursor_mono::{MonoCorrection, MonoParams};
use search::{
    apply_shift_for_mode, apply_tightened_precursor_tolerance, PrecursorCalMode, PreparedSearch,
    SearchIndex, SearchParams, TopNQueue,
};

pub(crate) struct RetrievalChoice {
    pub use_index: bool,
    pub report: Option<String>,
    pub refused_because: Option<&'static str>,
}

pub(crate) fn retrieval_choice(
    fragment_tol: model::tolerance::Tolerance,
    excluded_because: Option<&'static str>,
    flag: FragmentIndexFlag,
) -> RetrievalChoice {
    use model::tolerance::Tolerance;
    let tol_text = match fragment_tol {
        Tolerance::Ppm(v) => format!("{v} ppm"),
        Tolerance::Da(v) => format!("{v} Da"),
    };
    if let Some(reason) = excluded_because {
        return RetrievalChoice {
            use_index: false,
            report: None,
            refused_because: (flag == FragmentIndexFlag::On).then_some(reason),
        };
    }
    match flag {
        FragmentIndexFlag::Off => RetrievalChoice {
            use_index: false,
            report: Some("out-of-core index; the fragment-ion index was switched off".to_string()),
            refused_because: None,
        },
        FragmentIndexFlag::Auto | FragmentIndexFlag::On => {
            let report = match search::fragment_index::vote_peaks_for(fragment_tol) {
                0 => format!("out-of-core index, high-resolution fragments at {tol_text}"),
                n => format!(
                    "out-of-core index, low-resolution fragments at {tol_text}; \
                     the {n} most intense peaks vote"
                ),
            };
            RetrievalChoice {
                use_index: true,
                report: Some(report),
                refused_because: None,
            }
        }
    }
}

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
    measurer.measure(prepared, pending.iter(), Some(span_da));
    pending.sort_by(|a, b| index_order_mass(a).total_cmp(&index_order_mass(b)));
    eprintln!(
        "fragment-index: {} spectra sorted by precursor mass, chunks of <= {} spectra and <= {:.0} Da",
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
                eprintln!(
                    "Fragment error: offset {:+.2} ppm, spread {:.2} ppm from {} ions of {} PSMs \
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
            Err(why) => eprintln!(
                "Fragment error: not measured ({why}); the PIN features keep {:?} ({:.1}s)",
                prepared.scorer.feature_match_tolerance(),
                t.elapsed().as_secs_f64()
            ),
        }
    }
}

fn isobaric_peak_filter(explicit: Option<(f64, usize)>, high_res: bool) -> Option<(f64, usize)> {
    explicit.or(high_res.then_some((0.0, 0)))
}

fn output_base(
    output_pin: Option<&Path>,
    output_tsv: Option<&Path>,
    output_parquet: Option<&Path>,
) -> PathBuf {
    output_pin
        .or(output_tsv)
        .or(output_parquet)
        .map(Path::to_path_buf)
        .unwrap_or_else(|| PathBuf::from("andes.pin"))
}

fn scratch_tempdir(prefix: &str, scratch_dir: &Path) -> Result<tempfile::TempDir, String> {
    crate::cli::ensure_scratch_dir(scratch_dir)?;
    tempfile::Builder::new()
        .prefix(prefix)
        .tempdir_in(scratch_dir)
        .map_err(|e| crate::cli::scratch_dir_error(scratch_dir, &e))
}

fn rescore_pin_destination(
    output_pin: Option<&Path>,
    output_tsv: Option<&Path>,
    output_parquet: Option<&Path>,
    keep_pin: bool,
    scratch_dir: &Path,
    tmp: &mut Option<tempfile::TempDir>,
) -> Result<(PathBuf, PathBuf), Box<dyn std::error::Error>> {
    let base = output_base(output_pin, output_tsv, output_parquet);
    if output_pin.is_some() {
        return Ok((base.clone(), base));
    }
    if keep_pin {
        return Ok((base.with_extension("pin"), base));
    }
    let dir = scratch_tempdir("andes-pin-", scratch_dir)?;
    let path = dir.path().join("andes.pin");
    *tmp = Some(dir);
    Ok((path, base))
}

#[cfg(feature = "timstof")]
fn check_distinct_d_archive_runs(inputs: &[PathBuf]) -> Result<(), String> {
    let mut seen: std::collections::HashMap<String, &PathBuf> = std::collections::HashMap::new();
    for p in inputs {
        let Some((_, run)) = input::d_archive_kind(p) else {
            continue;
        };
        if let Some(first) = seen.insert(run.clone(), p) {
            return Err(format!(
                "the archives {} and {} both hold run '{run}'; give each run one archive",
                first.display(),
                p.display()
            ));
        }
    }
    Ok(())
}

fn unpack_d_archives(
    inputs: &mut [PathBuf],
    scratch_dir: &Path,
) -> Result<Option<tempfile::TempDir>, Box<dyn std::error::Error>> {
    if !inputs.iter().any(|p| crate::spectra::is_d_archive(p)) {
        return Ok(None);
    }
    #[cfg(feature = "timstof")]
    {
        check_distinct_d_archive_runs(inputs)?;
        let dir = scratch_tempdir("andes-d-", scratch_dir)?;
        for p in inputs
            .iter_mut()
            .filter(|p| crate::spectra::is_d_archive(p))
        {
            let t = std::time::Instant::now();
            let run = input::unpack_d_archive(p, dir.path())?;
            eprintln!(
                "Bruker .d archive {}: unpacked to {} ({:.1}s)",
                p.display(),
                run.display(),
                t.elapsed().as_secs_f64()
            );
            *p = run;
        }
        Ok(Some(dir))
    }
    #[cfg(not(feature = "timstof"))]
    {
        let _ = scratch_dir;
        Err(
            "this andes build has no Bruker .d (timsTOF) support; rebuild with \
             `--features timstof` to read archived .d inputs"
                .into(),
        )
    }
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
        if let Some(v) = crate::xl_run::xl_default_variable_mods(cli.db.variable_mods.as_ref(), x) {
            cli.db.variable_mods = Some(v);
        }
    }
    if cli.top_n == 0 {
        return Err("--top-n must be at least 1".into());
    }
    if cli.localize.is_some() && !cli.rescore && !cli.rescore_native {
        return Err("--localize ipf needs --rescore or --rescore-native".into());
    }
    if cli.rescore_native && !cli.percolator_args.trim().is_empty() {
        return Err(
            "native rescoring takes no percolator_args; drop --rescore-native to run Percolator"
                .into(),
        );
    }
    let scratch_dir = crate::cli::dir_or_output_dir(
        cli.temp_dir.as_deref(),
        &output_base(
            cli.output_pin.as_deref(),
            cli.output_tsv.as_deref(),
            cli.output_parquet.as_deref(),
        ),
    );
    cli.temp_dir = Some(scratch_dir.clone());
    let original_spectrum_paths = cli.spectrum.clone();
    let standalone = cli.glyco || cli.xl.is_some();
    let _d_archives = unpack_d_archives(&mut cli.spectrum, &scratch_dir)?;
    let mut _mbr_tmp: Option<tempfile::TempDir> = None;
    if cli.mbr && !cli.glyco {
        check_mbr_args(&cli)?;
        crate::quant_run::check_unique_stems(&cli.spectrum)?;
        if cli.quant_dir.is_none() && cli.label.is_none() {
            let dir = scratch_tempdir("andes-quant-", &scratch_dir)?;
            cli.quant_dir = Some(dir.path().to_path_buf());
            _mbr_tmp = Some(dir);
        }
    }
    if cli.quant_dir.is_some() {
        crate::quant_run::check_unique_stems(&cli.spectrum)?;
    }
    QuantSettings::check_args(&cli)?;
    if let Some(warning) = deuterated_tag_warning(cli.tmt, cli.precursor_tol) {
        eprintln!("{warning}");
    }
    let (search_enzyme, extra_enzymes) = parse_enzymes(&cli.db.enzyme)?;
    warn_if_universal_protease_combo(search_enzyme, &extra_enzymes);
    let spectrum_paths = &cli.spectrum;
    let spectrum_path: PathBuf = spectrum_paths[0].clone();
    let database_path: PathBuf = cli.database.clone().expect("database validated in main");
    let mut _rescore_tmp: Option<tempfile::TempDir> = None;
    let (output_pin_path, report_base) = rescore_pin_destination(
        cli.output_pin.as_deref(),
        cli.output_tsv.as_deref(),
        cli.output_parquet.as_deref(),
        cli.keep_pin,
        &scratch_dir,
        &mut _rescore_tmp,
    )?;
    if spectrum_paths.len() > 1 {
        eprintln!(
            "Multi-spectrum search: {} inputs → one PIN",
            spectrum_paths.len()
        );
    }

    let t_total = std::time::Instant::now();
    let t_phase = std::time::Instant::now();
    let uniprot_opts = input::UniProtOptions {
        variants: cli.db.uniprot_variants,
        isoforms: cli.db.uniprot_isoforms,
        chains: cli.db.uniprot_chains,
        evidence: cli.db.uniprot_evidence,
        knockoffs: Some(input::KnockoffOptions {
            per_variant: cli.db.knockoffs,
            enzyme: search_enzyme,
        }),
    };
    let ClaimDatabase {
        loaded,
        annotations,
        site_knockoffs,
    } = load_claim_database(&database_path, &cli.db.annotations, &uniprot_opts)?;
    let target_db = loaded.db;
    if let Some(peff) = &loaded.peff {
        log_claims("PEFF", &peff.claims);
        if peff.sites_without_mass + peff.processing_not_chain > 0 {
            eprintln!(
                "PEFF: {} entries; left out {} sites with no known mass and {} processed spans \
                 that are not mature chains",
                peff.entries, peff.sites_without_mass, peff.processing_not_chain
            );
        }
    }
    if let Some(stats) = &annotations {
        log_claims("Annotations", stats);
    }
    if site_knockoffs.knockoffs > 0 {
        eprintln!(
            "Site claims: {} annotated sites, {} knockoff sites ({} sites with no position for \
             one)",
            site_knockoffs.sites, site_knockoffs.knockoffs, site_knockoffs.untested
        );
    }
    eprintln!(
        "Loaded {} target proteins from {} [PHASE fasta_load: {:.2}s]",
        target_db.proteins.len(),
        database_path.display(),
        t_phase.elapsed().as_secs_f64()
    );
    match &loaded.uniprot {
        Some(stats) => log_uniprot_load(stats, &uniprot_opts),
        None => {
            for (set, flag) in [
                (cli.db.uniprot_variants, "--uniprot-variants"),
                (cli.db.uniprot_isoforms, "--uniprot-isoforms"),
                (cli.db.uniprot_chains, "--uniprot-chains"),
            ] {
                if set {
                    eprintln!(
                        "WARN: {flag} needs a UniProt XML database (*.xml[.gz]); {} is FASTA, \
                         so it does nothing.",
                        database_path.display()
                    );
                }
            }
        }
    }

    let decoy_strategy = cli.decoy_strategy;
    let proteases = search::decoy::Proteases {
        primary: search_enzyme,
        extra: extra_enzymes.clone(),
        num_tolerable_termini: match cli.enzyme_specificity {
            EnzymeSpecificity::Fully => 2,
            EnzymeSpecificity::Semi => 1,
            EnzymeSpecificity::NonSpecific => 0,
        },
    };
    if decoy_strategy == search::decoy::DecoyStrategy::Shuffle
        && proteases.is_non_specific()
        && target_db.proteins.iter().any(|p| {
            matches!(
                p.entry_kind(),
                model::EntryKind::Isoform | model::EntryKind::Chain
            )
        })
    {
        return Err(
            "--decoy-strategy shuffle cannot search isoform or mature-chain entries \
                    (--uniprot-isoforms, --uniprot-chains) non-specifically: shuffled \
                    whole-sequence decoys do not keep the peptides these entries share with \
                    their canonical entry. Use --decoy-strategy reverse."
                .into(),
        );
    }
    if decoy_strategy != search::decoy::DecoyStrategy::None {
        if let Some(warning) = search::decoy::ForeignDecoys::scan(&target_db).verdict()? {
            eprintln!("{warning}");
        }
    }
    let t_phase = std::time::Instant::now();
    let idx = SearchIndex::from_target_db_for(
        &target_db,
        &cli.db.decoy_prefix,
        cli.db.decoy_suffix.as_deref(),
        decoy_strategy,
        cli.decoy_seed,
        &proteases,
    );
    eprintln!(
        "[PHASE search_index_build: {:.2}s]",
        t_phase.elapsed().as_secs_f64()
    );

    let mut aa = search_aa_set(&cli.db, None)?;

    let (is_mzml, is_raw, is_d, is_mgf) = input_format_flags(&spectrum_path);

    input::mzml::init_ethcd_as_etd(cli.ethcd_activation == EthcdActivationFlag::Etd);

    if is_raw {
        #[cfg(feature = "thermo")]
        input::ThermoRawReader::open(&spectrum_path).map_err(|e| e.to_string())?;
        #[cfg(not(feature = "thermo"))]
        return Err("this andes build has no Thermo .raw support; rebuild with \
                    `--features thermo` (and run with the .NET 8 runtime installed). \
                    mzML/MGF inputs work without it."
            .into());
    }

    let (detected_activation, detected_instrument) = detect_acquisition_metadata(&spectrum_path);
    let instrument_was_detected = detected_instrument.is_some();
    for other_path in &spectrum_paths[1..] {
        let other = detect_acquisition_metadata(other_path);
        if let Some(warning) = compare_with_first_input(
            &spectrum_path,
            (detected_activation, detected_instrument),
            other_path,
            other,
        )? {
            eprintln!("{warning}");
        }
    }

    let peak_evidence = if !instrument_was_detected
        && cli.fragment_tol_ppm.is_none()
        && cli.fragment_tol_da.is_none()
    {
        let ev = sample_spectra(&spectrum_path, is_mzml, is_mgf, is_raw, cli.ms_level as u32)
            .map(|sample| input::infer_fragment_resolution(&sample));
        if let Some(e) = &ev {
            eprintln!(
                "Acquisition: no analyzer metadata; fragment resolution from the peaks: {} -> {}",
                e.describe(),
                match e.high_res {
                    Some(true) => "high resolution (HCD high-res model)",
                    Some(false) => "low resolution",
                    None => "undecided",
                }
            );
        }
        ev
    } else {
        None
    };
    let acquisition = Acquisition::resolve(
        detected_activation,
        detected_instrument,
        cli.fragment_tol_ppm,
        cli.fragment_tol_da,
        peak_evidence,
        sample_ms2_meta(&spectrum_path, is_mzml, is_raw, is_d),
    );
    if let (Some(said), Some(read)) = (cli.sdrf_dissociation, acquisition.activation) {
        if said != read {
            eprintln!(
                "WARN: the SDRF dissociation method is {} but {} records {}; andes follows the file",
                said.name(),
                spectrum_path.display(),
                read.name()
            );
        }
    }

    let detect = || {
        detect_isobaric_sampled(
            &spectrum_path,
            is_mzml,
            is_mgf,
            is_raw,
            cli.ms_level as u32,
            acquisition.high_res,
        )
    };
    let label: Option<input::IsobaricLabel> = match (cli.tmt, cli.protocol) {
        (Some(plex), _) => Some(isobaric_label_of_plex(plex)),
        (None, Protocol::Auto | Protocol::Tmt | Protocol::Itraq | Protocol::ItraqPhospho) => {
            detect()
        }
        (None, _) => None,
    };
    let label = label.or_else(|| {
        if cli.protocol == Protocol::Auto {
            label_from_mods(&aa)
        } else {
            None
        }
    });
    let route_protocol = routing_protocol(cli.protocol, label);
    let t_phase = std::time::Instant::now();
    let (mut param, routed_activation, routed_instrument, resolved_model_id) = {
        let detected = (detected_activation, detected_instrument);
        let (activation, instrument_opt): (ActivationMethod, Option<InstrumentType>) =
            match detected {
                (Some(method), Some(inst)) => {
                    let chosen = match cli.fragmentation {
                        Fragmentation::Auto => method,
                        Fragmentation::Cid => ActivationMethod::CID,
                        Fragmentation::Etd => ActivationMethod::ETD,
                        Fragmentation::Hcd => ActivationMethod::HCD,
                        Fragmentation::Uvpd => ActivationMethod::UVPD,
                    };
                    if chosen != method {
                        eprintln!(
                            "Param resolver: auto-detected activation = {} (instrument = {}) from {}, \
                             OVERRIDDEN by --fragmentation {}",
                            method.name(), inst.name(), spectrum_path.display(), chosen.name()
                        );
                    } else {
                        eprintln!(
                            "Param resolver: auto-detected activation = {} (instrument = {}) from {}",
                            method.name(), inst.name(), spectrum_path.display()
                        );
                    }
                    (chosen, Some(inst))
                }
                (None, Some(inst)) => {
                    let explicit = cli_fragmentation_to_activation_opt(cli.fragmentation);
                    let chosen = explicit.unwrap_or(if inst.is_high_resolution() {
                        ActivationMethod::HCD
                    } else {
                        ActivationMethod::CID
                    });
                    eprintln!(
                        "Param resolver: instrument = {} from {}, no activation recorded; \
                         activation = {} [{}]",
                        inst.name(),
                        spectrum_path.display(),
                        chosen.name(),
                        if explicit.is_some() {
                            "--fragmentation"
                        } else {
                            "default of the analyzer class"
                        }
                    );
                    (chosen, Some(inst))
                }
                (method, None) => resolve_metadataless_selection(
                    method,
                    cli.fragmentation,
                    cli.fragment_tol_ppm,
                    cli.fragment_tol_da,
                    acquisition.inferred_high_res(),
                ),
            };

        let (model_id, p) = load_param_from_store(
            activation,
            instrument_opt,
            route_protocol,
            search_enzyme,
            cli.model_store.as_deref(),
            cli.model_id_override.as_deref(),
        )?;
        eprintln!("Param model: {model_id} (from store)");
        (p, activation, instrument_opt, model_id)
    };
    let gbdt_k = cli.gbdt_max_trees.unwrap_or(if cli.glyco {
        0
    } else {
        crate::cli::DEFAULT_GBDT_MAX_TREES
    });
    if gbdt_k > 0 {
        let k = gbdt_k;
        let truncate = |slot: &mut Option<Arc<scoring_crate::gbdt_eval::GbdtPeakModel>>,
                        name: &str| {
            if let Some(g) = slot.as_ref() {
                let before = g.trees.len();
                if k < before {
                    let mut t = (**g).clone();
                    t.trees.truncate(k);
                    *slot = Some(Arc::new(t));
                    eprintln!(
                        "--gbdt-max-trees {k}: {name} ensemble {before} -> {k} trees \
                         (changes predictions; not byte-identical)"
                    );
                }
            }
        };
        truncate(&mut param.frag_intensity_model, "frag-intensity");
        truncate(&mut param.rich_ion_model, "rich-ion");
    }
    match cli.protocol {
        Protocol::Tmt => param.data_type.protocol = model::protocol::Protocol::TMT,
        Protocol::Itraq => param.data_type.protocol = model::protocol::Protocol::ITRAQ,
        Protocol::ItraqPhospho => {
            param.data_type.protocol = model::protocol::Protocol::ITRAQPhospho
        }
        Protocol::Auto => match label {
            Some(l @ (input::IsobaricLabel::Tmt | input::IsobaricLabel::TmtPro)) => {
                eprintln!(
                    "Protocol resolver: {} {l:?} → isobaric protocol (windowed peak filter on low-res MS2)",
                    if cli.tmt.is_some() {
                        "--tmt"
                    } else {
                        "detected"
                    }
                );
                param.data_type.protocol = model::protocol::Protocol::TMT;
            }
            Some(l @ (input::IsobaricLabel::Itraq | input::IsobaricLabel::Itraq8)) => {
                eprintln!(
                    "Protocol resolver: {} {l:?} → isobaric protocol (windowed peak filter on low-res MS2)",
                    if cli.tmt.is_some() {
                        "--tmt"
                    } else {
                        "detected"
                    }
                );
                param.data_type.protocol = model::protocol::Protocol::ITRAQ;
            }
            None => {}
        },
        _ => {}
    }
    if cli.db.fixed_mods.is_none() {
        let tag = match (param.data_type.protocol, label) {
            (
                model::protocol::Protocol::TMT,
                Some(l @ (input::IsobaricLabel::Tmt | input::IsobaricLabel::TmtPro)),
            ) => Some(l.tag()),
            (model::protocol::Protocol::TMT, _) => Some(input::IsobaricLabel::Tmt.tag()),
            (
                model::protocol::Protocol::ITRAQ | model::protocol::Protocol::ITRAQPhospho,
                Some(l @ (input::IsobaricLabel::Itraq | input::IsobaricLabel::Itraq8)),
            ) => Some(l.tag()),
            (model::protocol::Protocol::ITRAQ | model::protocol::Protocol::ITRAQPhospho, _) => {
                Some(input::IsobaricLabel::Itraq.tag())
            }
            (_, Some(l)) if cli.tmt.is_some() => Some(l.tag()),
            _ => None,
        };
        if let Some((name, mass)) = tag {
            aa = search_aa_set(&cli.db, Some((name, mass)))?;
            eprintln!(
                "Protocol resolver: injected {name} fixed mod (+{mass:.4} on K + peptide N-term) \
                 into the candidate set (no --fixed-mods given)"
            );
        }
    }
    let mut scorer = RankScorer::new(&param);
    scorer.set_acquisition_high_res(Some(acquisition.high_res));
    let frag_tol_override = cli_fragment_tol_override(cli.fragment_tol_ppm, cli.fragment_tol_da);
    if let Some(t) = frag_tol_override {
        if instrument_was_detected {
            eprintln!(
                "--fragment-tol-* {t:?} overrides the feature tolerance of the detected {} analyzer",
                acquisition.analyzer.map_or("", |i| i.name())
            );
        }
        scorer.set_fragment_tol_override(frag_tol_override);
    }
    eprintln!(
        "[PHASE param_and_scorer: {:.2}s]",
        t_phase.elapsed().as_secs_f64()
    );

    report_annotated_site_mods(&aa, &target_db);
    warn_restricted_artefacts(&aa, cli.db.uniprot_variants);

    let restricted_site_masses: Vec<f64> = aa
        .distinct_mods()
        .iter()
        .filter_map(|m| m.annotated_site_mass)
        .collect();
    let label_spec = match &cli.label {
        Some(text) => Some(crate::quant_run::label_spec(&cli, text)?),
        None => None,
    };
    if let Some(spec) = &label_spec {
        aa = with_label_mods(&aa, spec)?;
    }
    let fixed_terminal_mods = aa.fixed_terminal_deltas();
    let mut params = SearchParams::default_tryptic(aa);
    params.precursor_tolerance = PrecursorTolerance::symmetric(cli.precursor_tol);
    let (charge_min, charge_max) = cli.charge;
    params.charge_range = charge_min..=charge_max;
    let iso_default = if cli.glyco { (0, 2) } else { (-1, 2) };
    let (iso_min, iso_max) = cli.isotope_error.unwrap_or(iso_default);
    params.isotope_error_range = iso_min..=iso_max;
    let has_ms1 = is_mzml || is_raw || is_d;
    let chimeric_active = cli.chimeric && has_ms1;
    if cli.chimeric && !has_ms1 {
        eprintln!(
            "WARN: --chimeric requires MS1 data (mzML, Thermo .raw or Bruker .d); the input \
             is MGF, so the co-isolation cascade is disabled and the search runs normally."
        );
    }
    if chimeric_active && cli.ms_level != 2 {
        eprintln!(
            "WARN: --ms-level={} is ignored under --chimeric; the cascade always searches MS2.",
            cli.ms_level
        );
    }
    params.chimeric = chimeric_active;
    if cli.glyco_index_sequon_only && !cli.glyco {
        return Err(
            "--glyco-index-sequon-only requires --glyco (the standard peptide \
                    search scores every candidate)"
                .into(),
        );
    }
    params.require_nxst_sequon = cli.glyco && cli.glyco_index_sequon_only;
    if params.require_nxst_sequon {
        eprintln!(
            "candidate index: N-X-S/T sequon-bearing peptides only (--glyco-index-sequon-only)"
        );
    }
    let mono_active = cli.precursor_mono == PrecursorMonoFlag::Auto
        && (cli.glyco || cli.open_search)
        && has_ms1;
    if cli.precursor_mono == PrecursorMonoFlag::Auto
        && (cli.glyco || cli.open_search)
        && !mono_active
    {
        eprintln!(
            "WARN: --precursor-mono auto needs an input with MS1 (mzML, Thermo .raw or \
             Bruker .d); \
             {} has none — precursors are left as recorded.",
            spectrum_path.display()
        );
    }
    let mono_params = MonoParams::default();
    if mono_active {
        eprintln!(
            "precursor-mono: auto (max shift {}, min fit {:.2}, min gain {:.2}, min SNR {:.1}, \
             {} ppm, back-off {})",
            mono_params.max_shift,
            mono_params.min_fit,
            mono_params.min_gain,
            mono_params.min_snr,
            mono_params.tol_ppm,
            mono_params.backoff
        );
    }
    params.chimeric_isolation_halfwidth_da = 1.5;
    params.chimeric_max_coisolated = cli.chimeric_max_coisolated;
    params.chimeric_max_kl = cli.chimeric_max_kl;
    params.top_n_psms_per_spectrum = cli.top_n;
    params.enzyme = search_enzyme;
    params.extra_enzymes = extra_enzymes.clone();
    params.num_tolerable_termini = proteases.num_tolerable_termini;
    let peak_filter = match cli.peak_filter.as_deref() {
        None => None,
        Some(spec) => {
            let (w, k) = spec.split_once(':').ok_or_else(|| {
                format!("--peak-filter expects WINDOW_DA:PEAKS (e.g. 100:20), got `{spec}`")
            })?;
            let w: f64 = w
                .trim()
                .parse()
                .map_err(|_| format!("--peak-filter window `{w}` is not a number"))?;
            let k: usize = k
                .trim()
                .parse()
                .map_err(|_| format!("--peak-filter count `{k}` is not an integer"))?;
            Some((w, k))
        }
    };
    let peak_filter = isobaric_peak_filter(peak_filter, acquisition.high_res);
    if acquisition.high_res
        && cli.peak_filter.is_none()
        && matches!(
            param.data_type.protocol,
            model::protocol::Protocol::TMT
                | model::protocol::Protocol::ITRAQ
                | model::protocol::Protocol::ITRAQPhospho
        )
    {
        eprintln!(
            "isobaric windowed peak filter: off (high-resolution MS2; --peak-filter WINDOW:K turns it on)"
        );
    }
    scoring_crate::scoring::init_scoring_settings(scoring_crate::scoring::ScoringSettings {
        peak_filter,
    });

    params.max_missed_cleavages = if cli.glyco {
        let explicit = *EXPLICIT_MISSED_CLEAVAGES.get().unwrap_or(&false);
        if explicit && cli.max_missed_cleavages < 3 {
            eprintln!(
                "note: --glyco normally raises missed cleavages to 3 (sequon-bearing \
                 tryptic peptides often need a third to reach a scoreable length); \
                 honouring your explicit --max-missed-cleavages {}. Expect fewer \
                 glycopeptide IDs in exchange for a smaller candidate index.",
                cli.max_missed_cleavages
            );
            cli.max_missed_cleavages
        } else {
            cli.max_missed_cleavages.max(3)
        }
    } else {
        cli.max_missed_cleavages
    };
    params.min_peaks = cli.min_peaks;
    params.min_length = cli.min_length;
    params.max_length = cli.max_length;
    params.max_variable_mods_per_peptide = cli.max_mods;
    if let Some(spec) = &label_spec {
        let enzymes: Vec<_> = std::iter::once(params.enzyme)
            .chain(params.extra_enzymes.iter().copied())
            .collect();
        let extra = crate::quant_run::label_run::max_label_sites(
            spec,
            &enzymes,
            params.num_tolerable_termini,
            params.max_missed_cleavages,
            &idx.db,
            params.max_length,
        );
        params.max_label_sites = extra;
        eprintln!(
            "--label: channels {}; label mods added to the search (up to {extra} labelled sites \
             per peptide, apart from the {} other variable mods)",
            spec.describe(),
            params.max_variable_mods_per_peptide
        );
    }
    params.precursor_cal_mode = cli.precursor_cal;
    if let Some(n) = cli.cal_min_spec_keys {
        params.cal_min_spec_keys = n;
    }
    params.precursor_mass_shift_ppm = 0.0;
    params.refine_select_psm_fdr = cli.psm_fdr();
    params.score_mode = match cli.score {
        ScoreFlag::Rank => search::ScoreMode::Rank,
        ScoreFlag::Strong => search::ScoreMode::Strong,
        ScoreFlag::Auto => {
            if acquisition.high_res {
                search::ScoreMode::Strong
            } else {
                search::ScoreMode::Rank
            }
        }
    };
    let retrieval_excluded: Option<&'static str> = if cli.glyco {
        Some("--glyco looks candidates up by backbone mass")
    } else if cli.xl.is_some() {
        Some("--xl looks candidates up by peptide mass")
    } else if !search::fragment_index::covers_activation(param.data_type.activation) {
        Some(
            "it holds only b/y ions, and this activation method fragments into \
             c/z or a/x ions",
        )
    } else {
        None
    };
    if let (false, Some(reason)) = (standalone, retrieval_excluded) {
        eprintln!(
            "[auto] {} fragmentation -> per-spectrum candidate enumeration (the fragment-ion \
             index is not used: {reason})",
            param.data_type.activation.name()
        );
    }
    let choice = retrieval_choice(
        scorer.feature_match_tolerance(),
        retrieval_excluded,
        cli.fragment_index,
    );
    params.fragment_index_top_k = if choice.use_index {
        search::search_params::FRAGMENT_INDEX_TOP_K
    } else {
        0
    };
    if cli.open_search {
        if cli.glyco {
            return Err("--open-search is not supported with --glyco".into());
        }
        if cli.chimeric {
            return Err("--open-search is not combined with --chimeric".into());
        }
        if cli.label.is_some() {
            return Err("--open-search is not combined with --label".into());
        }
        if cli.tmt.is_some() {
            return Err("--open-search is not combined with --tmt".into());
        }
        if cli.lfq || cli.mbr {
            return Err("--open-search is not combined with --lfq or --mbr".into());
        }
        if cli.refine {
            return Err("--open-search is not combined with --refine".into());
        }
        if cli.xl.is_some() {
            return Err("--open-search is not combined with --xl".into());
        }
        if cli.isotope_error.is_some() {
            return Err(
                "--open-search chooses the isotope offset from the fragments (range 0..=2); \
                 remove --isotope-error"
                    .into(),
            );
        }
        if !choice.use_index {
            let why = retrieval_excluded.unwrap_or("it was disabled with --fragment-index off");
            return Err(format!(
                "--open-search needs the fragment-ion index (it cannot enumerate the wide \
                 precursor window): {why}"
            )
            .into());
        }
        let (delta_min, delta_max) = cli.open_search_delta;
        let top_k = cli.open_search_top_k.max(1);
        let open = search::search_params::OpenSearch {
            delta_min,
            delta_max,
            top_k,
        };
        params.precursor_tolerance = open.precursor_tolerance();
        params.open_search = Some(open);
        params.isotope_error_range = 0..=2;
        params.fragment_index_top_k = top_k;
        params.precursor_cal_mode = PrecursorCalMode::Off;
        eprintln!(
            "open search: mass-offset window {delta_min:+}..{delta_max:+} Da, {top_k} candidates \
             per spectrum (fragment-ion shortlist), isotope range 0..=2 (reported offset is \
             isotope-corrected)"
        );
    }
    if let Some(refusal) = choice.refused_because {
        eprintln!(
            "WARN: the fragment-ion index cannot be used here ({refusal}); using enumeration."
        );
    }
    if let Some(reason) = &choice.report {
        if choice.use_index {
            eprintln!(
                "candidate retrieval: fragment-ion index ({reason}; {} candidates per spectrum, \
                 >= {} matched b/y ions)",
                params.fragment_index_top_k,
                search::search_params::FRAGMENT_INDEX_MIN_MATCHED
            );
        } else {
            eprintln!("candidate retrieval: per-spectrum enumeration ({reason})");
        }
    }

    if !cli.glyco && arg_present("--glyco-hcd-pair") {
        eprintln!(
            "WARN: --glyco-hcd-pair has no effect without --glyco (it only drives \
             paired-scan candidate generation in glyco mode); ignoring it."
        );
    }
    if cli.xl.is_some() {
        let mut unsupported: Vec<&str> = Vec::new();
        for (on, flag) in [
            (cli.output_tsv.is_some(), "--output-tsv"),
            (cli.output_parquet.is_some(), "--output-parquet"),
            (cli.rescore_native, "--rescore-native"),
            (cli.refine, "--refine"),
            (cli.psm_fdr.is_some(), "--psm-fdr"),
            (cli.pep.is_some(), "--pep"),
            (cli.lfq, "--lfq"),
            (cli.mbr, "--mbr"),
            (cli.quant_dir.is_some(), "--quant-dir"),
            (cli.chimeric, "--chimeric"),
            (cli.localize.is_some(), "--localize"),
            (cli.tmt.is_some(), "--tmt"),
            (cli.label.is_some(), "--label"),
        ] {
            if on {
                unsupported.push(flag);
            }
        }
        if !unsupported.is_empty() {
            return Err(format!(
                "--xl does not support: {} (crosslink mode writes its own .xl.pin and \
                 tables, scored by the pair score)",
                unsupported.join(", ")
            )
            .into());
        }
    }
    if cli.glyco {
        let mut unsupported: Vec<&str> = Vec::new();
        if cli.output_tsv.is_some() {
            unsupported.push("--output-tsv");
        }
        if cli.output_parquet.is_some() {
            unsupported.push("--output-parquet");
        }
        if cli.rescore {
            unsupported.push("--rescore");
        }
        if cli.rescore_native {
            unsupported.push("--rescore-native");
        }
        if cli.refine {
            unsupported.push("--refine");
        }
        if cli.psm_fdr.is_some() {
            unsupported.push("--psm-fdr");
        }
        if cli.pep.is_some() {
            unsupported.push("--pep");
        }
        if !unsupported.is_empty() {
            return Err(format!(
                "--glyco does not support: {} (glyco mode writes a standalone \
                 .glyco.pin and skips the standard PIN/rescore/TSV/Parquet/refine \
                 pipeline; run Percolator on the .glyco.pin separately)",
                unsupported.join(", ")
            )
            .into());
        }
        if !matches!(
            decoy_strategy,
            search::decoy::DecoyStrategy::None | search::decoy::DecoyStrategy::Sequon
        ) {
            eprintln!(
                "WARN: --glyco with --decoy-strategy {ds} GENERATES decoys ({ds} of the \
                 target proteome), which does NOT preserve N-X-S/T sequon density that \
                 glyco scoring gates on — the resulting .glyco.pin FDR is \
                 ANTI-CONSERVATIVE and its q-values should NOT be trusted. Use \
                 `--decoy-strategy sequon`, or supply an EXTERNAL target+decoy \
                 FASTA with `--decoy-strategy none --decoy-prefix <PREFIX>`.",
                ds = crate::cli::decoy_strategy_name(decoy_strategy)
            );
        }
        if cli.glyco_hcd_pair && spectrum_paths.len() > 1 {
            eprintln!(
                "WARN: --glyco-hcd-pair is disabled for {} --spectrum inputs: paired-scan \
                 generation pairs ETD<->HCD over the concatenated spectra and would cross \
                 file boundaries. Run one file at a time to use paired generation.",
                spectrum_paths.len()
            );
        }
    }
    if params.chimeric && cli.refine {
        eprintln!(
            "WARN: --refine + --chimeric do not fully stack in this release — chimeric \
             already explains many of the spectra refinement would otherwise rescue, so \
             refinement adds much less on top of chimeric than it does alone. See \
             docs/benchmarks/README.md for the measured combination."
        );
    }
    if params.score_mode == search::ScoreMode::Strong {
        eprintln!("score mode: strong (ranking + PIN RawScore use StrongScore)");
    }

    eprintln!("──────── andes resolved parameters ────────");
    eprintln!("  spectra        : {}", spectrum_path.display());
    let model_source = if cli.model_id_override.is_some() {
        "--model"
    } else {
        "selected from store"
    };
    eprintln!("  model          : {resolved_model_id} [{model_source}]");
    let activation_source = if cli.fragmentation != Fragmentation::Auto {
        "--fragmentation"
    } else if detected_activation.is_some() {
        "detected"
    } else {
        "default"
    };
    eprintln!(
        "  activation     : {} [{activation_source}] (model trained on {})",
        routed_activation.name(),
        param.data_type.activation.name()
    );
    let instrument_source = acquisition.resolution_source.tag();
    eprintln!(
        "  instrument     : {} [{instrument_source}] (model trained on {})",
        routed_instrument.map_or("LowRes", |i| i.name()),
        param.data_type.instrument.name()
    );
    eprintln!("  acquisition    : {}", acquisition.describe());
    eprintln!("  protocol       : {:?}", param.data_type.protocol);
    if extra_enzymes.is_empty() {
        eprintln!(
            "  enzyme         : {} ({:?} termini, <={} missed cleavages)",
            search_enzyme.name(),
            cli.enzyme_specificity,
            params.max_missed_cleavages
        );
    } else {
        let extras: Vec<&str> = extra_enzymes.iter().map(|e| e.name()).collect();
        eprintln!("  enzyme         : {} (primary) + {} (multi-protease union) ({:?} termini, <={} missed cleavages)",
                  search_enzyme.name(), extras.join(","), cli.enzyme_specificity, params.max_missed_cleavages);
    }
    let mod_names = |fixed: bool| {
        let names: Vec<String> = params
            .aa_set
            .distinct_mods()
            .iter()
            .filter(|m| m.fixed == fixed && m.parts.is_empty())
            .map(|m| model::mod_spec::describe(m))
            .collect();
        if names.is_empty() {
            "none".to_string()
        } else {
            names.join(", ")
        }
    };
    eprintln!("  fixed mods     : {}", mod_names(true));
    eprintln!("  variable mods  : {}", mod_names(false));
    eprintln!(
        "  max var-mods   : {} per peptide",
        params.max_variable_mods_per_peptide
    );
    eprintln!(
        "  peptide length : {}-{}",
        params.min_length, params.max_length
    );
    eprintln!(
        "  precursor tol  : {:?} (calibration: {:?})",
        params.precursor_tolerance, params.precursor_cal_mode
    );
    eprintln!(
        "  charge range   : {}-{}",
        params.charge_range.start(),
        params.charge_range.end()
    );
    eprintln!(
        "  isotope errors : {}..={}",
        params.isotope_error_range.start(),
        params.isotope_error_range.end()
    );
    eprintln!(
        "  decoy          : {:?} (prefix {})   chimeric: {}",
        decoy_strategy, cli.db.decoy_prefix, params.chimeric
    );
    eprintln!("───────────────────────────────────────────");

    const CHUNK_SIZE: usize = 5000;

    let t_phase = std::time::Instant::now();

    eprintln!("Using {} worker threads", cli.threads);

    let fragment_tol_da = 0.5_f64;

    let bench_cap = if cli.max_spectra > 0 {
        cli.max_spectra
    } else {
        usize::MAX
    };
    let ms_level_u32 = cli.ms_level as u32;
    if (is_raw || is_d) && cli.ms_level != 2 {
        eprintln!(
            "WARN: --ms-level={} is ignored for native .raw/.d input; these formats \
             search MS2 (identification) scans only.",
            cli.ms_level
        );
    }

    if is_raw && params.precursor_cal_mode != PrecursorCalMode::Off {
        eprintln!(
            "WARN: --precursor-cal is not yet supported for Thermo .raw input; \
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
        eprintln!(
            "[PHASE calibration_prepass: {:.2}s] (includes opening or building the candidate index)",
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
                eprintln!(
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
            sample_spectra(
                &spectrum_path,
                is_mzml,
                is_mgf,
                is_raw,
                cli.ms_level as u32,
            )
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
                drop(rescorer);
                let optimize_path = optimize_dir.path().to_path_buf();
                if let Err(e) = optimize_dir.close() {
                    eprintln!(
                        "WARN: optimize: could not remove {}: {e}",
                        optimize_path.display()
                    );
                }
            }
            None => eprintln!(
                "WARN: --optimize could not read a sample; keeping the resolved parameters"
            ),
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
                    eprintln!(
                        "WARN: per-scan model routing: {e}; {} scans use the run's model",
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
                    if let Some(g) = slot.as_ref() {
                        if gbdt_max_trees < g.trees.len() {
                            let mut t = (**g).clone();
                            t.trees.truncate(gbdt_max_trees);
                            *slot = Some(Arc::new(t));
                        }
                    }
                }
            }
            p.data_type.protocol = run_protocol;
            let mut sc = RankScorer::new(&p);
            sc.set_acquisition_high_res(Some(high_res));
            if let Some(t) = frag_tol_override {
                sc.set_fragment_tol_override(Some(t));
            }
            eprintln!(
                "Per-scan routing: {} scans -> model {model_id} (the run is routed {})",
                activation.name(),
                routed_activation.name()
            );
            Some(sc)
        });
        prepared = prepared.with_per_scan_scorers(routed_activation, loader);
    }
    eprintln!(
        "[PHASE candidate_index_open: {:.2}s]",
        t_index.elapsed().as_secs_f64()
    );
    eprintln!(
        "PreparedSearch: out-of-core candidate-index: mmap \
         (base peptides resolved lazily per spectrum)"
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
            eprintln!(
                "Protocol resolver: --tmt {}: the MS2 is low-resolution and the MS3 scans carry \
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
            eprintln!(
                "WARN: --tmt/--lfq without --rescore or --rescore-native: every rank-1 target \
                 PSM is quantified; filter the tables by the q-values you compute downstream, \
                 or add --rescore so only PSMs at q<=--psm-fdr are quantified."
            );
        }
        if !cli.glyco {
            eprintln!("Quantification: {}", s.describe());
        }
    }
    let predictors = if standalone {
        None
    } else {
        crate::predicted::Predictors::from_cli(&cli, &acquisition, param.data_type.instrument)?
    };
    if let Some(p) = &predictors {
        eprintln!("Predict: {}", p.describe());
    }
    let keep_peaks = cli.refine
        || cli.localize.is_some()
        || predictors.as_ref().is_some_and(|p| p.needs_peaks());
    let mut all_spectra: Vec<Spectrum> = Vec::new();
    let mut all_queues: Vec<TopNQueue> = Vec::new();
    let mut run_ends: Vec<usize> = Vec::with_capacity(spectrum_paths.len());

    let t_search_start = std::time::Instant::now();

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
    let _ = mzml_warn_ms_level_emitted;

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
            eprintln!(
                "WARN: --precursor-mono-dump {} ignored: the correction is not active.",
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
            eprintln!(
                "=== spectrum [{}/{}] {} ===",
                file_idx + 1,
                spectrum_paths.len(),
                input_path.display()
            );
        }

        if ms1_input && !(file_is_mzml || file_is_raw || file_is_d) {
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
        let file_stats = if ms1_input {
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
                eprintln!(
                    "fragment-index: {} spectra sorted by precursor mass, chunks of <= {} spectra and <= {:.0} Da",
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
            eprintln!(
                "{}: streamed {} MS2 spectra ({} MS1 scans linked) from {}",
                if chimeric_active {
                    "chimeric mode"
                } else {
                    "precursor-mono"
                },
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
            let capture_ms1 = quant.as_ref().is_some_and(|q| q.settings.needs_ms1());
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

            let mass_ordered = params.fragment_index_top_k > 0 && !standalone;
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

            match parser_handle.join() {
                Ok(Ok((stats, scans))) => {
                    run_scans = scans;
                    stats
                }
                Ok(Err(e)) => return Err(format!("parser thread error: {e}").into()),
                Err(_) => return Err("parser thread panicked".into()),
            }
        };
        require_spectra(
            input_path,
            if ms1_input { 2 } else { cli.ms_level },
            all_spectra.len() - spectra_before,
        )
        .map_err(|e| crate::spectra::with_first_read_error(e, &file_stats.first_errors))?;
        merge_parse_stats(&mut parse_stats, file_stats);
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

    let search_elapsed = t_search_start.elapsed();
    eprintln!(
        "Loaded+scored {} spectra from {} in chunks of {} [PHASE stream_search: {:.2}s]",
        all_spectra.len(),
        spectrum_path.display(),
        CHUNK_SIZE,
        t_phase.elapsed().as_secs_f64()
    );
    eprintln!(
        "[PHASE chunk_candidate_prep: {:.2}s] (part of stream_search: fragment-index build or record expansion, summed over chunks)",
        prepared.chunk_prep_seconds()
    );
    for (activation, own_model, scans) in prepared.per_scan_routes() {
        eprintln!(
            "Per-scan routing: {scans} {} scan scorings with {} (run routed {})",
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
        eprintln!("Bench mode: capped at {} spectra", cli.max_spectra);
    }

    prepared.sync_materialized_candidates();

    let mono_columns = mono_active;
    let mono_active = mono_active && {
        if mono_stats.fitted == 0 {
            eprintln!(
                "WARN: --precursor-mono auto: no MS2 could be fitted ({} MS2 seen, none with \
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
        eprintln!(
            "precursor-mono: {} of {} MS2 fitted (linked MS1 + known charge), {} shifted \
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
        eprintln!(
            "precursor-mono: isotope-error window {}..={} -> {}..={} for the {} of {} MS2 \
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
        eprintln!(
            "Standard peptide search skipped in --{} mode (PSMs unused; the driver \
             derives its own candidates); {} spectra loaded ({:.2}s)",
            if cli.glyco { "glyco" } else { "xl" },
            spectra.len(),
            search_elapsed.as_secs_f64()
        );
    } else {
        let non_empty = queues.iter().filter(|q| !q.is_empty()).count();
        eprintln!(
            "Search complete: {non_empty} / {} spectra have PSMs (search wall: {:.2}s)",
            spectra.len(),
            search_elapsed.as_secs_f64()
        );
    }
    if cli.open_search {
        let (attempts, wins) = prepared.open_search_stats();
        eprintln!(
            "open search: {attempts} modification-aware re-scores ({wins} improved the node \
             score)"
        );
    }

    if cli.xl.is_some() {
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
            eprintln!(
                "WARN: refine is high-res-only and the data is low-res; skipping refinement."
            );
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
            eprintln!(
                "[PHASE refinement: {:.2}s]",
                t_refine.elapsed().as_secs_f64()
            );
            out
        }
    } else {
        None
    };

    let t_phase = std::time::Instant::now();
    let refine_merged = refine_output.is_some();
    let ipf_scoring = cli.localize.is_some().then(|| prepared.ipf_scoring(&scorer));
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
    eprintln!(
        "Wrote PIN: {} [PHASE pin_write: {:.2}s] [PHASE TOTAL: {:.2}s]",
        output_pin_path.display(),
        t_phase.elapsed().as_secs_f64(),
        t_total.elapsed().as_secs_f64()
    );
    if let Some(p) = &predictors {
        p.add_to_pin(&output_pin_path, &spectra, &run_ends)?;
    }
    if cli.refine {
        if refine_merged {
            eprintln!("Refinement PSMs merged into unified PIN (Pass-1 ⊕ Pass-2).");
        } else {
            eprintln!("Refinement produced no Pass-2 PSMs; unified PIN contains Pass-1 only.");
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
        eprintln!(
            "WARN: {which} set without --rescore or --rescore-native, so no rescoring \
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
            eprintln!("Rescore: Percolator backend = {}", backend.describe());
            let user_args: Vec<String> = cli
                .percolator_args
                .split_whitespace()
                .map(|s| s.to_string())
                .collect();
            let extra = output::percolator_args_with_threads(&user_args, cli.threads);
            let t_resc = std::time::Instant::now();
            let (map, decoys) =
                output::run_percolator_with_decoys(&backend, &output_pin_path, &extra)
                    .map_err(|e| format!("{e}"))?;
            rescore_decoys = Some(decoys);
            eprintln!(
                "Rescore: Percolator produced {} target PSMs [{:.2}s]",
                map.len(),
                t_resc.elapsed().as_secs_f64()
            );
            let path = output::decoy_results_path(&output_pin_path);
            let text =
                std::fs::read_to_string(&path).map_err(|e| format!("{}: {e}", path.display()))?;
            rescore_decoy_psms =
                output::parse_psm_results(&text).map_err(|e| format!("{}: {e}", path.display()))?;
            if has_site_claims {
                row_scores = Some(
                    output::score_all_rows(&backend, &output_pin_path, &extra)
                        .map_err(|e| format!("scoring every PIN row for the site posteriors: {e}"))?,
                );
            }
            if cli.mbr {
                mbr_source = format!("the Percolator results of {}", output_pin_path.display());
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
                eprintln!(
                    "Rescore: {shared} variant-class PSMs share their peptide with a reference \
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
            eprintln!(
                "Rescore: FDR within each class (targets/decoys, method, pi0): {}",
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
            eprintln!(
                "Warning: {untested} variant PSMs at q <= {proteoform_fdr} have no knockoff \
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
            eprintln!(
                "Rescore: variant PSMs at q<={proteoform_fdr} per substitution class, before -> \
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
            eprintln!(
                "Rescore: site-claim PSMs at q<={fdr_threshold} per modification and residue, \
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
                Ok(()) => eprintln!(
                    "Wrote site claims: {} ({} sites)",
                    path.display(),
                    sites.rows.len()
                ),
                Err(e) => eprintln!("WARN: could not write {}: {e}", path.display()),
            }
        }
        if let Some(scoring) = ipf_scoring.as_ref() {
            use rayon::prelude::*;
            let t_ipf = std::time::Instant::now();
            let map_ref: &std::collections::HashMap<String, output::percolator::PercolatorPsm> = map;
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
                                    && map_ref.get(id).is_some_and(|p| p.q_value <= cutoffs.of(id))
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
                eprintln!(
                    "WARN: --localize ipf: no identified target PSM has a variable \
                     modification with more than one possible site; {} is empty",
                    path.display()
                );
            }
            let truncated = rows.iter().filter(|r| r.truncated).count();
            eprintln!(
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
            Ok(()) => eprintln!(
                "Wrote protein groups: {} ({} groups; {} at protein q<=0.01, {} at q<=0.05)",
                path.display(),
                inf.groups.len(),
                at(0.01),
                at(0.05)
            ),
            Err(e) => eprintln!("WARN: could not write {}: {e}", path.display()),
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
            Ok(()) => eprintln!(
                "Wrote filtered TSV: {} ({} PSMs at q<={}{})",
                filtered_path.display(),
                accepted.len(),
                fdr,
                pep_note
            ),
            Err(e) => eprintln!("WARN: could not write {}: {e}", filtered_path.display()),
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
    let id_summary = identified.as_ref().map(|(accepted, cut_note)| {
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
    eprint!("{}", run_stats.render_with(id_summary.as_ref()));
    let stats_path = {
        let stem = report_base
            .file_stem()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| "andes".to_string());
        let dir = report_base.parent().unwrap_or_else(|| Path::new("."));
        dir.join(format!("{stem}.statistics.log"))
    };
    match output::write_statistics_log(&run_stats, id_summary.as_ref(), &stats_path) {
        Ok(()) => eprintln!("Wrote statistics: {}", stats_path.display()),
        Err(e) => eprintln!("WARN: could not write {}: {e}", stats_path.display()),
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
            Ok(()) => eprintln!("Wrote QPX bundle: {}", parquet_dir.display()),
            Err(e) => {
                return Err(format!(
                    "failed to write requested --output-parquet {}: {e}",
                    parquet_dir.display()
                )
                .into());
            }
        }
    }

    let omitted_psms = match (&quant, rescore_map.as_ref()) {
        (Some(_), Some(map)) => output::omitted_psm_rows(
            &spectra,
            &queues,
            &pin_candidates,
            Some(map),
            rescore_decoys.as_ref(),
            &knockoff_rows,
        )
        .map_err(|e| format!("quant: {e}"))?,
        _ => std::collections::HashSet::new(),
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
            Ok(()) => eprintln!(
                "Wrote proteoform table: {} ({} family members, {} families)",
                path.display(),
                rows.len(),
                rows.iter()
                    .map(|r| r.family.as_str())
                    .collect::<std::collections::BTreeSet<_>>()
                    .len()
            ),
            Err(e) => eprintln!("WARN: could not write {}: {e}", path.display()),
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
        eprintln!("Bench mode: skipping TSV write.");
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
        eprintln!("Wrote TSV: {}", tsv_path.display());
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

fn log_claims(source: &str, stats: &input::AnnotationLoadStats) {
    let dropped: Vec<String> = stats
        .dropped
        .iter()
        .map(|(reason, n)| format!("{reason} ({n})"))
        .collect();
    let files = if stats.files > 0 {
        format!(" from {} files", stats.files)
    } else {
        String::new()
    };
    eprintln!(
        "{source}: {} rows{files}: {} variant, {} knockoff and {} chain entries, {} sites on \
         reference proteins, {} rows merged into an equal claim, {} into an equal database \
         entry{}{}",
        stats.rows,
        stats.variant_entries,
        stats.knockoff_entries,
        stats.chain_entries,
        stats.sites,
        stats.merged,
        stats.merged_into_database,
        if dropped.is_empty() {
            ""
        } else {
            "; left out: "
        },
        dropped.join(", "),
    );
    if stats.knockoff_entries > 0 {
        eprintln!(
            "{source}: {} single-residue variants with no position for a knockoff",
            stats.variants_without_knockoff
        );
    }
    if stats.knockoffs_dropped > 0 {
        eprintln!(
            "{source}: {} earlier knockoff entries dropped, within the knockoff distance of a \
             variant these files list",
            stats.knockoffs_dropped
        );
    }
}

fn log_uniprot_load(stats: &input::UniProtLoadStats, opts: &input::UniProtOptions) {
    eprintln!(
        "UniProt XML: {} entries, {} annotated modification sites, {} variant entries \
         (evidence >= {}; {} features below it, {} not fitting their sequence, \
         {} on isoforms not searched)",
        stats.entries,
        stats.mod_sites,
        stats.variant_proteins,
        opts.evidence.as_str(),
        stats.below_evidence,
        stats.unusable_features,
        stats.isoform_features,
    );
    if opts.variants && opts.knockoffs.is_some_and(|k| k.per_variant > 0) {
        eprintln!(
            "UniProt XML: {} knockoff variant entries ({} single-residue variants with no \
             position for a knockoff)",
            stats.knockoff_proteins, stats.variants_without_knockoff,
        );
    }
    if opts.isoforms || opts.chains {
        eprintln!(
            "UniProt XML: {} isoform entries ({} isoforms not rebuildable), {} mature-chain entries",
            stats.isoform_proteins, stats.unbuilt_isoforms, stats.chain_proteins,
        );
    }
    if !stats.unknown_modifications.is_empty() {
        let mut names: Vec<(&String, &usize)> = stats.unknown_modifications.iter().collect();
        names.sort_by(|a, b| b.1.cmp(a.1).then(a.0.cmp(b.0)));
        let total: usize = names.iter().map(|(_, n)| **n).sum();
        let examples: Vec<String> = names
            .iter()
            .take(5)
            .map(|(name, n)| format!("{name} ({n})"))
            .collect();
        eprintln!(
            "UniProt XML: {total} modification features have no single known mass and are not \
             used as sites, e.g. {}",
            examples.join(", ")
        );
    }
}

const ARTEFACT_MASSES: [f64; 4] = [0.984016, 15.994915, 43.005814, 14.015650];

fn warn_restricted_artefacts(aa: &model::AminoAcidSet, variants_on: bool) {
    if !variants_on {
        return;
    }
    for m in aa.distinct_mods() {
        if m.annotated_site_mass.is_some()
            && ARTEFACT_MASSES
                .iter()
                .any(|a| (m.mass_delta - a).abs() <= model::protein::SITE_MASS_TOLERANCE)
        {
            eprintln!(
                "WARN: {} ({:+.4} Da) is restricted to annotated sites while \
                 --uniprot-variants is on: variants of the same mass lose the \
                 competing explanation at unannotated residues.",
                m.name, m.mass_delta
            );
        }
    }
}

fn report_annotated_site_mods(aa: &model::AminoAcidSet, target_db: &model::ProteinDb) {
    let mut mods: Vec<(&str, f64)> = aa
        .distinct_mods()
        .iter()
        .filter_map(|m| m.annotated_site_mass.map(|mass| (m.name.as_str(), mass)))
        .collect();
    mods.sort_by(|a, b| a.0.cmp(b.0).then(a.1.total_cmp(&b.1)));
    mods.dedup_by(|b, a| a.0 == b.0 && a.1.to_bits() == b.1.to_bits());
    for (name, mass) in mods {
        let (mut sites, mut proteins) = (0usize, 0usize);
        let same = |m: f64| (m - mass).abs() <= model::protein::SITE_MASS_TOLERANCE;
        for p in &target_db.proteins {
            let knockoffs: Vec<u32> = p
                .annotation
                .as_deref()
                .map(|a| {
                    a.site_claims
                        .iter()
                        .filter(|c| c.knockoff_of.is_some() && same(c.mass_delta))
                        .map(|c| c.position)
                        .collect()
                })
                .unwrap_or_default();
            let n = p
                .mod_sites()
                .iter()
                .filter(|s| same(s.mass_delta) && !knockoffs.contains(&s.position))
                .count();
            if n > 0 {
                sites += n;
                proteins += 1;
            }
        }
        if sites == 0 {
            eprintln!(
                "WARN: {name} ({mass:+.4} Da) is restricted to annotated sites (restrict_to: annotated), \
                 but no target protein annotates a site of that mass, so it is never placed. \
                 Annotated sites come from a UniProt XML database (--database *.xml[.gz]); \
                 --uniprot-evidence sets which are used."
            );
        } else {
            eprintln!(
                "Annotated-site mod {name} ({mass:+.4} Da): {sites} sites on {proteins} target \
                 proteins"
            );
        }
    }
}

fn check_mbr_args(cli: &Cli) -> Result<(), String> {
    if cli.spectrum.len() < 2 {
        return Err(format!(
            "--mbr transfers identifications between runs and needs two or more --spectrum inputs (got {})",
            cli.spectrum.len()
        ));
    }
    if !cli.lfq && cli.tmt.is_none() && cli.label.is_none() {
        return Err(
            "--mbr needs --lfq (label-free), --tmt (reporter ions) or --label (labelled pairs)"
                .into(),
        );
    }
    if !cli.rescore && !cli.rescore_native {
        return Err(
            "--mbr needs --rescore (or --rescore-native): the library ranks precursors by the rescored PSM scores"
                .into(),
        );
    }
    if cli.tmt.is_some() {
        if cli.lfq {
            return Err(
                "--mbr transfers either label-free MS1 quantities (--lfq) or reporter ions \
                 (--tmt), not both; run them separately"
                    .into(),
            );
        }
        if cli.tmt_level.is_some_and(|level| level >= 3) {
            return Err(
                "--mbr with --tmt needs the reporter ions in the MS2 scans and does not apply \
                 to SPS-MS3 reporters (--tmt-level 3)"
                    .into(),
            );
        }
        return Ok(());
    }
    crate::quant_cmd::evidence_params(cli.lfq_tol)?;
    Ok(())
}

fn release_at_exit(candidates: Vec<search::candidate_gen::Candidate>) {
    std::mem::forget(candidates);
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

fn deuterated_tag_warning(
    plex: Option<quant::Plex>,
    precursor_tol: model::Tolerance,
) -> Option<String> {
    use model::Tolerance;
    let plex = plex.filter(|p| p.has_deuterated_channels())?;
    let shift = 3.0 * quant::isobaric::DEUTERATED_TAG_SHIFT_DA;
    if precursor_tol.as_da(1000.0) >= shift {
        return None;
    }
    let tol = match precursor_tol {
        Tolerance::Ppm(v) => format!("{v}ppm"),
        Tolerance::Da(v) => format!("{v}da"),
    };
    Some(format!(
        "WARN: --tmt {}: each deuterated TMTpro tag is {:.1} mDa heavier than the TMTpro mass \
         (304.2071 Da), so a 1000 Da peptide with three of them sits {:.1} ppm above its \
         candidate, outside --precursor-tol {tol}. Widen it to 10ppm or more, or the \
         deuterated channels' precursors are missed.",
        plex.name(),
        quant::isobaric::DEUTERATED_TAG_SHIFT_DA * 1e3,
        shift / 1000.0 * 1e6,
    ))
}

fn isobaric_label_of_plex(plex: quant::Plex) -> input::IsobaricLabel {
    match plex {
        quant::Plex::Tmt6 | quant::Plex::Tmt10 | quant::Plex::Tmt11 => input::IsobaricLabel::Tmt,
        quant::Plex::Tmt16 | quant::Plex::Tmt18 | quant::Plex::Tmt32 | quant::Plex::Tmt35 => {
            input::IsobaricLabel::TmtPro
        }
        quant::Plex::Itraq4 => input::IsobaricLabel::Itraq,
        quant::Plex::Itraq8 => input::IsobaricLabel::Itraq8,
    }
}

#[cfg(test)]
mod deuterated_tag_warning_tests {
    use super::deuterated_tag_warning;
    use model::Tolerance;
    use quant::Plex;

    #[test]
    fn warns_only_for_deuterated_kits_under_a_narrow_window() {
        let w = deuterated_tag_warning(Some(Plex::Tmt35), Tolerance::Ppm(5.0)).unwrap();
        assert!(w.contains("8.8 ppm") && w.contains("5ppm"), "{w}");
        assert!(deuterated_tag_warning(Some(Plex::Tmt32), Tolerance::Da(0.005)).is_some());
        assert!(deuterated_tag_warning(Some(Plex::Tmt35), Tolerance::Ppm(10.0)).is_none());
        assert!(deuterated_tag_warning(Some(Plex::Tmt35), Tolerance::Da(0.02)).is_none());
        assert!(deuterated_tag_warning(Some(Plex::Tmt18), Tolerance::Ppm(5.0)).is_none());
        assert!(deuterated_tag_warning(None, Tolerance::Ppm(5.0)).is_none());
    }
}

#[cfg(test)]
mod retrieval_choice_tests {
    use super::{retrieval_choice, FragmentIndexFlag};
    use model::tolerance::Tolerance;

    const HI: Tolerance = Tolerance::Ppm(20.0);
    const LO: Tolerance = Tolerance::Da(0.5);

    #[test]
    fn high_res_uses_the_index() {
        let c = retrieval_choice(HI, None, FragmentIndexFlag::Auto);
        assert!(c.use_index);
        assert!(c.report.unwrap().contains("high-resolution"));
    }

    #[test]
    fn low_res_uses_the_capped_index() {
        let c = retrieval_choice(LO, None, FragmentIndexFlag::Auto);
        assert!(c.use_index);
        assert!(c.report.unwrap().contains("most intense peaks vote"));
    }

    #[test]
    fn off_is_honoured_on_low_res() {
        let c = retrieval_choice(LO, None, FragmentIndexFlag::Off);
        assert!(!c.use_index);
    }

    #[test]
    fn glyco_never_uses_the_fragment_index() {
        const GLYCO: Option<&str> = Some("--glyco looks candidates up by backbone mass");
        for flag in [FragmentIndexFlag::Auto, FragmentIndexFlag::On] {
            let c = retrieval_choice(HI, GLYCO, flag);
            assert!(!c.use_index, "glyco looks candidates up by backbone mass");
        }
        let c = retrieval_choice(HI, GLYCO, FragmentIndexFlag::On);
        assert_eq!(
            c.refused_because,
            Some("--glyco looks candidates up by backbone mass")
        );
    }

    #[test]
    fn off_is_honoured_where_the_index_would_otherwise_run() {
        let c = retrieval_choice(HI, None, FragmentIndexFlag::Off);
        assert!(!c.use_index);
        assert!(c.report.unwrap().contains("switched off"));
    }

    #[test]
    fn activations_outside_the_indexed_ion_series_keep_enumeration() {
        use model::activation::ActivationMethod;
        for a in [
            ActivationMethod::CID,
            ActivationMethod::HCD,
            ActivationMethod::PQD,
            ActivationMethod::EThcD,
        ] {
            assert!(search::fragment_index::covers_activation(a), "{a:?}");
        }
        for a in [ActivationMethod::ETD, ActivationMethod::UVPD] {
            assert!(!search::fragment_index::covers_activation(a), "{a:?}");
        }
        const ETD: Option<&str> = Some("b/y only");
        for tol in [HI, LO] {
            assert!(!retrieval_choice(tol, ETD, FragmentIndexFlag::Auto).use_index);
            let forced = retrieval_choice(tol, ETD, FragmentIndexFlag::On);
            assert!(!forced.use_index);
            assert_eq!(forced.refused_because, ETD);
        }
    }

    #[test]
    fn a_tight_dalton_tolerance_counts_as_high_resolution() {
        let tight = retrieval_choice(Tolerance::Da(0.02), None, FragmentIndexFlag::Auto);
        assert!(tight.report.unwrap().contains("high-resolution"));
        let wide = retrieval_choice(Tolerance::Da(0.2), None, FragmentIndexFlag::Auto);
        assert!(wide.report.unwrap().contains("low-resolution"));
    }
}

#[cfg(test)]
mod rescore_pin_destination_tests {
    use super::rescore_pin_destination;
    use std::path::{Path, PathBuf};

    #[test]
    fn explicit_output_pin_is_pin_and_base() {
        let mut tmp = None;
        let (pin, base) = rescore_pin_destination(
            Some(Path::new("r/x.pin")),
            None,
            None,
            false,
            Path::new("."),
            &mut tmp,
        )
        .unwrap();
        assert_eq!(pin, PathBuf::from("r/x.pin"));
        assert_eq!(base, PathBuf::from("r/x.pin"));
        assert!(tmp.is_none());
    }

    #[test]
    fn keep_pin_writes_next_to_report_base() {
        let mut tmp = None;
        let (pin, base) = rescore_pin_destination(
            None,
            Some(Path::new("r/out.tsv")),
            Some(Path::new("r/out.idparquet")),
            true,
            Path::new("."),
            &mut tmp,
        )
        .unwrap();
        assert_eq!(pin, PathBuf::from("r/out.pin"));
        assert_eq!(base, PathBuf::from("r/out.tsv"));
        assert!(tmp.is_none());

        let (pin, _) = rescore_pin_destination(
            None,
            None,
            Some(Path::new("q.idparquet")),
            true,
            Path::new("."),
            &mut tmp,
        )
        .unwrap();
        assert_eq!(pin, PathBuf::from("q.pin"));
        let (pin, base) =
            rescore_pin_destination(None, None, None, true, Path::new("."), &mut tmp).unwrap();
        assert_eq!(pin, PathBuf::from("andes.pin"));
        assert_eq!(base, PathBuf::from("andes.pin"));
        assert!(tmp.is_none());
    }

    #[test]
    fn no_keep_pin_is_created_in_the_scratch_dir() {
        let scratch = tempfile::tempdir().unwrap();
        let mut tmp = None;
        let (pin, base) = rescore_pin_destination(
            None,
            Some(Path::new("out.tsv")),
            None,
            false,
            scratch.path(),
            &mut tmp,
        )
        .unwrap();
        let dir = tmp.as_ref().expect("temp dir held").path().to_path_buf();
        assert!(dir.starts_with(scratch.path()));
        assert!(pin.starts_with(&dir));
        assert_eq!(base, PathBuf::from("out.tsv"));
        drop(tmp);
        assert_eq!(std::fs::read_dir(scratch.path()).unwrap().count(), 0);
    }

    #[test]
    fn unwritable_scratch_dir_is_an_error_naming_the_flag() {
        let dir = tempfile::tempdir().unwrap();
        let blocker = dir.path().join("blocker");
        std::fs::write(&blocker, b"x").unwrap();
        let mut tmp = None;
        let err = rescore_pin_destination(
            None,
            Some(Path::new("out.tsv")),
            None,
            false,
            &blocker.join("sub"),
            &mut tmp,
        )
        .unwrap_err();
        assert!(err.to_string().contains("--temp-dir"), "{err}");
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

#[cfg(test)]
mod isobaric_peak_filter_tests {
    use super::isobaric_peak_filter;

    #[test]
    fn high_res_turns_the_isobaric_filter_off_unless_set() {
        assert_eq!(isobaric_peak_filter(None, false), None);
        assert_eq!(isobaric_peak_filter(None, true), Some((0.0, 0)));
        assert_eq!(
            isobaric_peak_filter(Some((100.0, 20)), true),
            Some((100.0, 20))
        );
        assert_eq!(
            isobaric_peak_filter(Some((50.0, 10)), false),
            Some((50.0, 10))
        );
    }
}

#[cfg(all(test, feature = "timstof"))]
mod d_archive_tests {
    use super::check_distinct_d_archive_runs;
    use std::path::PathBuf;

    #[test]
    fn archives_of_one_run_name_are_refused_before_unpacking() {
        let ok: Vec<PathBuf> = ["a.d.zip", "b.d.tar", "c.mzML"].map(PathBuf::from).to_vec();
        assert!(check_distinct_d_archive_runs(&ok).is_ok());
        let dup: Vec<PathBuf> = ["x/run.d.zip", "y/run.d.tar"].map(PathBuf::from).to_vec();
        let err = check_distinct_d_archive_runs(&dup).unwrap_err();
        assert!(
            err.contains("x/run.d.zip") && err.contains("y/run.d.tar"),
            "{err}"
        );
    }
}

fn with_label_mods(
    aa: &model::AminoAcidSet,
    spec: &quant::label::LabelSpec,
) -> Result<model::AminoAcidSet, String> {
    let (mods, channels) = crate::quant_run::label_run::search_label_mods(spec);
    let declared = aa.distinct_mods();
    let mut b = AminoAcidSetBuilder::new_standard();
    for m in &declared {
        if m.fixed {
            b = b.add_fixed_mod((*m).clone());
        } else {
            b = b.add_variable_mod((*m).clone());
        }
    }
    for m in mods {
        let dup = declared.iter().any(|d| {
            d.name == m.name && d.residue == m.residue && d.location == m.location && !d.fixed
        });
        if !dup {
            b = b.add_variable_mod(m);
        }
    }
    b.with_label_channels(channels)
        .build()
        .map_err(|e| format!("--label: building the amino-acid set: {e}"))
}
