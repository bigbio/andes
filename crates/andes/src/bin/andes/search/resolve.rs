//! The parameter-resolution stages of `andes search`: flag checks, scratch paths, the database,
//! input detection, the scoring model, the search parameters and the resolved-parameter log.
//! Each stage returns a typed result that `run` destructures into its working locals.

mod checks;
mod database;
mod inputs;
mod scorer;

use std::path::PathBuf;

use model::activation::ActivationMethod;
use model::{AminoAcidSetBuilder, InstrumentType, PrecursorTolerance};
use scoring_crate::{Param, RankScorer};
use search::precursor::monoisotopic::MonoParams;
use search::{PrecursorCalMode, SearchIndex, SearchParams};

use crate::acquisition::Acquisition;
use crate::cli::{Cli, Fragmentation, PrecursorMonoFlag, Protocol, ScoreFlag};
use crate::spectra::input_format_flags;
use crate::{arg_present, EXPLICIT_MISSED_CLEAVAGES};
pub(super) use checks::check_flag_combinations;
use checks::warn_restricted_artefacts;
pub(super) use database::load_search_database;
use database::{max_missed_cleavages, report_annotated_site_mods};
pub(super) use inputs::{detect_inputs, prepare_run};
use scorer::{isobaric_peak_filter, retrieval_choice};
pub(super) use scorer::{resolve_scoring_model, truncate_gbdt};

/// The temporary directories a run owns. They are removed when this value is dropped, in the
/// order of the fields.
pub(super) struct TempDirs {
    pub(super) _rescore: Option<tempfile::TempDir>,
    pub(super) _mbr: Option<tempfile::TempDir>,
    pub(super) _d_archives: Option<tempfile::TempDir>,
}

/// The directories and output paths a run works with.
pub(super) struct RunPaths {
    pub(super) scratch_dir: PathBuf,
    pub(super) original_spectrum_paths: Vec<PathBuf>,
    pub(super) output_pin_path: PathBuf,
    pub(super) report_base: PathBuf,
    pub(super) temp_dirs: TempDirs,
}

/// The loaded protein database, the search index built from it and the digestion settings.
pub(super) struct SearchDatabase {
    pub(super) target_db: model::ProteinDb,
    pub(super) idx: SearchIndex,
    pub(super) search_enzyme: model::enzyme::Enzyme,
    pub(super) extra_enzymes: Vec<model::enzyme::Enzyme>,
    pub(super) proteases: search::decoy::Proteases,
    pub(super) decoy_strategy: search::decoy::DecoyStrategy,
}

/// What the first spectrum file and the flags say about the input.
pub(super) struct RunInputs {
    pub(super) is_mzml: bool,
    pub(super) is_raw: bool,
    pub(super) is_d: bool,
    pub(super) is_mgf: bool,
    pub(super) detected_activation: Option<ActivationMethod>,
    pub(super) detected_instrument: Option<InstrumentType>,
    pub(super) instrument_was_detected: bool,
    pub(super) acquisition: Acquisition,
    pub(super) label: Option<input::IsobaricLabel>,
    pub(super) route_protocol: Protocol,
}

/// The scoring model and scorer chosen for the run.
pub(super) struct ScoringModel {
    pub(super) param: Param,
    pub(super) scorer: RankScorer,
    pub(super) routed_activation: ActivationMethod,
    pub(super) routed_instrument: Option<InstrumentType>,
    pub(super) resolved_model_id: String,
    pub(super) gbdt_k: usize,
    pub(super) frag_tol_override: Option<model::tolerance::Tolerance>,
}

/// The search parameters and the values derived while resolving them.
pub(super) struct ResolvedSearch {
    pub(super) params: SearchParams,
    pub(super) chimeric_active: bool,
    pub(super) mono_active: bool,
    pub(super) mono_params: MonoParams,
    pub(super) restricted_site_masses: Vec<f64>,
    pub(super) fixed_terminal_mods: Vec<(model::ModLocation, f64)>,
}

/// Builds the search parameters from the flags, the input and the chosen model.
pub(super) fn resolve_search_params(
    cli: &Cli,
    mut aa: model::AminoAcidSet,
    inputs: &RunInputs,
    model: &ScoringModel,
    db: &SearchDatabase,
    standalone: bool,
) -> Result<ResolvedSearch, Box<dyn std::error::Error>> {
    let spectrum_paths = &cli.spectrum;
    let spectrum_path: PathBuf = spectrum_paths[0].clone();
    let (is_mzml, is_raw, is_d) = (inputs.is_mzml, inputs.is_raw, inputs.is_d);
    let acquisition = &inputs.acquisition;
    let param = &model.param;
    let scorer = &model.scorer;
    let target_db = &db.target_db;
    let idx = &db.idx;
    let search_enzyme = db.search_enzyme;
    let extra_enzymes = &db.extra_enzymes;
    let proteases = &db.proteases;
    let decoy_strategy = db.decoy_strategy;
    report_annotated_site_mods(&aa, target_db);
    warn_restricted_artefacts(&aa, cli.db.uniprot_variants);

    let restricted_site_masses: Vec<f64> = aa
        .distinct_mods()
        .iter()
        .filter_map(|m| m.annotated_site_mass)
        .collect();
    let label_spec = match &cli.label {
        Some(text) => Some(crate::quant::run::label_spec(cli, text)?),
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
        log::warn!(target: "search",
            "--chimeric requires MS1 data (mzML, Thermo .raw or Bruker .d); the input \
             is MGF, so the co-isolation cascade is disabled and the search runs normally."
        );
    }
    if chimeric_active && cli.ms_level != 2 {
        log::warn!(target: "search",
            "--ms-level={} is ignored under --chimeric; the cascade always searches MS2.",
            cli.ms_level
        );
    }
    params.chimeric = chimeric_active;
    if cli.dia {
        if !(is_mzml || is_raw || is_d) {
            return Err(
                "--dia needs mzML, Thermo .raw or Bruker .d input (isolation windows and MS1)"
                    .into(),
            );
        }
        if spectrum_paths.iter().any(|p| input_format_flags(p).2) {
            log::warn!(
                "--dia on Bruker .d (diaPASEF) is not FDR-calibrated yet. \
                 On an entrapment test, q <= 0.01 measured a 5-7% true FDP. \
                 Do not use these q-values as a 1% FDR."
            );
        }
        params.dia = true;
        params.dia_ms1_tol_ppm = cli.dia_ms1_tol_ppm;
    }
    params.require_nxst_sequon = cli.glyco && cli.glyco_index_sequon_only;
    if params.require_nxst_sequon {
        log::info!(
            target: "candidate-index",
            "N-X-S/T sequon-bearing peptides only (--glyco-index-sequon-only)"
        );
    }
    let mono_active =
        cli.precursor_mono == PrecursorMonoFlag::Auto && (cli.glyco || cli.open_search) && has_ms1;
    if cli.precursor_mono == PrecursorMonoFlag::Auto
        && (cli.glyco || cli.open_search)
        && !mono_active
    {
        log::warn!(target: "search",
            "--precursor-mono auto needs an input with MS1 (mzML, Thermo .raw or \
             Bruker .d); \
             {} has none — precursors are left as recorded.",
            spectrum_path.display()
        );
    }
    let mono_params = MonoParams::default();
    if mono_active {
        log::info!(
            target: "precursor-mono",
            "auto (max shift {}, min fit {:.2}, min gain {:.2}, min SNR {:.1}, \
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
    if cli.dia && cli.predict.contains(&crate::cli::PredictFeature::Ms2) {
        params.top_n_psms_per_spectrum = cli.top_n.max(crate::dia::RERANK_POOL);
    }
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
        log::info!(target: "search",
            "isobaric windowed peak filter: off (high-resolution MS2; --peak-filter WINDOW:K turns it on)"
        );
    }
    let dia_d = cli.dia && cli.spectrum.iter().all(|p| input_format_flags(p).2);
    let vote_filter = peak_filter.filter(|&(w, k)| dia_d && w > 0.0 && k > 0);
    scoring_crate::psm::init_scoring_settings(scoring_crate::psm::ScoringSettings {
        peak_filter: if vote_filter.is_some() {
            Some((0.0, 0))
        } else {
            peak_filter
        },
    });
    params.dia_vote_filter = vote_filter;

    if cli.glyco
        && *EXPLICIT_MISSED_CLEAVAGES.get().unwrap_or(&false)
        && cli.max_missed_cleavages < 3
    {
        log::info!(target: "search",
            "--glyco normally raises missed cleavages to 3 (sequon-bearing \
             tryptic peptides often need a third to reach a scoreable length); \
             honouring your explicit --max-missed-cleavages {}. Expect fewer \
             glycopeptide IDs in exchange for a smaller candidate index.",
            cli.max_missed_cleavages
        );
    }
    params.max_missed_cleavages = max_missed_cleavages(cli);
    params.min_peaks = cli.min_peaks;
    params.min_length = cli.min_length;
    params.max_length = cli.max_length;
    params.max_variable_mods_per_peptide = cli.max_mods;
    if let Some(spec) = &label_spec {
        let enzymes: Vec<_> = std::iter::once(params.enzyme)
            .chain(params.extra_enzymes.iter().copied())
            .collect();
        let extra = crate::quant::run::labels::max_label_sites(
            spec,
            &enzymes,
            params.num_tolerable_termini,
            params.max_missed_cleavages,
            &idx.db,
            params.max_length,
        );
        params.max_label_sites = extra;
        log::info!(target: "search",
            "--label: channels {}; label mods added to the search (up to {extra} labelled sites \
             per peptide, apart from the {} other variable mods)",
            spec.describe(),
            params.max_variable_mods_per_peptide
        );
    }
    params.precursor_cal_mode = cli.precursor_cal;
    if cli.dia && params.precursor_cal_mode != PrecursorCalMode::Off {
        log::warn!(target: "search", "--precursor-cal is ignored under --dia (a DIA scan has no selected precursor)");
        params.precursor_cal_mode = PrecursorCalMode::Off;
    }
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
    } else if !search::fragments::index::covers_activation(param.data_type.activation) {
        Some(
            "it holds only b/y ions, and this activation method fragments into \
             c/z or a/x ions",
        )
    } else {
        None
    };
    if let (false, Some(reason)) = (standalone, retrieval_excluded) {
        log::info!(target: "search",
            "{} fragmentation -> per-spectrum candidate enumeration (the fragment-ion \
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
        search::params::FRAGMENT_INDEX_TOP_K
    } else {
        0
    };
    if cli.open_search {
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
        let open = search::params::OpenSearch {
            delta_min,
            delta_max,
            top_k,
        };
        params.precursor_tolerance = open.precursor_tolerance();
        params.open_search = Some(open);
        params.isotope_error_range = 0..=2;
        params.fragment_index_top_k = top_k;
        params.precursor_cal_mode = PrecursorCalMode::Off;
        log::info!(
            target: "open-search",
            "mass-offset window {delta_min:+}..{delta_max:+} Da, {top_k} candidates \
             per spectrum (fragment-ion shortlist), isotope range 0..=2 (reported offset is \
             isotope-corrected)"
        );
    }
    if let Some(refusal) = choice.refused_because {
        log::warn!(target: "search", "the fragment-ion index cannot be used here ({refusal}); using enumeration.");
    }
    if let Some(reason) = &choice.report {
        if choice.use_index {
            log::info!(
                target: "candidate-retrieval",
                "fragment-ion index ({reason}; {} candidates per spectrum, \
                 >= {} matched b/y ions)",
                params.fragment_index_top_k,
                search::params::FRAGMENT_INDEX_MIN_MATCHED
            );
        } else {
            log::info!(target: "candidate-retrieval", "per-spectrum enumeration ({reason})");
        }
    }

    if !cli.glyco && arg_present("--glyco-hcd-pair") {
        log::warn!(target: "search",
            "--glyco-hcd-pair has no effect without --glyco (it only drives \
             paired-scan candidate generation in glyco mode); ignoring it."
        );
    }
    if cli.glyco {
        if !matches!(
            decoy_strategy,
            search::decoy::DecoyStrategy::None | search::decoy::DecoyStrategy::Sequon
        ) {
            log::warn!(target: "search",
                "--glyco with --decoy-strategy {ds} GENERATES decoys ({ds} of the \
                 target proteome), which does NOT preserve N-X-S/T sequon density that \
                 glyco scoring gates on — the resulting .glyco.pin FDR is \
                 ANTI-CONSERVATIVE and its q-values should NOT be trusted. Use \
                 `--decoy-strategy sequon`, or supply an EXTERNAL target+decoy \
                 FASTA with `--decoy-strategy none --decoy-prefix <PREFIX>`.",
                ds = crate::cli::decoy_strategy_name(decoy_strategy)
            );
        }
        if cli.glyco_hcd_pair && spectrum_paths.len() > 1 {
            log::warn!(target: "search",
                "--glyco-hcd-pair is disabled for {} --spectrum inputs: paired-scan \
                 generation pairs ETD<->HCD over the concatenated spectra and would cross \
                 file boundaries. Run one file at a time to use paired generation.",
                spectrum_paths.len()
            );
        }
    }
    if params.chimeric && cli.refine {
        log::warn!(target: "search",
            "--refine + --chimeric do not fully stack in this release — chimeric \
             already explains many of the spectra refinement would otherwise rescue, so \
             refinement adds much less on top of chimeric than it does alone. See \
             docs/benchmarks/README.md for the measured combination."
        );
    }
    if params.score_mode == search::ScoreMode::Strong {
        log::info!(target: "search", "score mode: strong (ranking + PIN RawScore use StrongScore)");
    }

    Ok(ResolvedSearch {
        params,
        chimeric_active,
        mono_active,
        mono_params,
        restricted_site_masses,
        fixed_terminal_mods,
    })
}

/// Logs the resolved parameters of the run.
pub(super) fn log_resolved_parameters(
    cli: &Cli,
    inputs: &RunInputs,
    model: &ScoringModel,
    db: &SearchDatabase,
    params: &SearchParams,
) {
    let spectrum_path: PathBuf = cli.spectrum[0].clone();
    let detected_activation = inputs.detected_activation;
    let acquisition = &inputs.acquisition;
    let param = &model.param;
    let resolved_model_id = &model.resolved_model_id;
    let routed_activation = model.routed_activation;
    let routed_instrument = model.routed_instrument;
    let search_enzyme = db.search_enzyme;
    let extra_enzymes = &db.extra_enzymes;
    let decoy_strategy = db.decoy_strategy;
    log::info!(target: "run", "──────── andes resolved parameters ────────");
    log::info!(target: "run", "spectra        : {}", spectrum_path.display());
    let model_source = if cli.model_id_override.is_some() {
        "--model"
    } else {
        "selected from store"
    };
    log::info!(target: "run", "model          : {resolved_model_id} [{model_source}]");
    let activation_source = if cli.fragmentation != Fragmentation::Auto {
        "--fragmentation"
    } else if detected_activation.is_some() {
        "detected"
    } else {
        "default"
    };
    log::info!(
        target: "run",
        "activation     : {} [{activation_source}] (model trained on {})",
        routed_activation.name(),
        param.data_type.activation.name()
    );
    let instrument_source = acquisition.resolution_source.tag();
    log::info!(
        target: "run",
        "instrument     : {} [{instrument_source}] (model trained on {})",
        routed_instrument.map_or("LowRes", |i| i.name()),
        param.data_type.instrument.name()
    );
    log::info!(target: "run", "acquisition    : {}", acquisition.describe());
    log::info!(target: "run", "protocol       : {:?}", param.data_type.protocol);
    if extra_enzymes.is_empty() {
        log::info!(
            target: "run",
            "enzyme         : {} ({:?} termini, <={} missed cleavages)",
            search_enzyme.name(),
            cli.enzyme_specificity,
            params.max_missed_cleavages
        );
    } else {
        let extras: Vec<&str> = extra_enzymes.iter().map(|e| e.name()).collect();
        log::info!(target: "run", "enzyme         : {} (primary) + {} (multi-protease union) ({:?} termini, <={} missed cleavages)",
                  search_enzyme.name(), extras.join(","), cli.enzyme_specificity, params.max_missed_cleavages);
    }
    let mod_names = |fixed: bool| {
        let names: Vec<String> = params
            .aa_set
            .distinct_mods()
            .iter()
            .filter(|m| m.fixed == fixed && m.parts.is_empty())
            .map(|m| model::modification::spec::describe(m))
            .collect();
        if names.is_empty() {
            "none".to_string()
        } else {
            names.join(", ")
        }
    };
    log::info!(target: "run", "fixed mods     : {}", mod_names(true));
    log::info!(target: "run", "variable mods  : {}", mod_names(false));
    log::info!(
        target: "run",
        "max var-mods   : {} per peptide",
        params.max_variable_mods_per_peptide
    );
    log::info!(
        target: "run",
        "peptide length : {}-{}",
        params.min_length, params.max_length
    );
    log::info!(
        target: "run",
        "precursor tol  : {:?} (calibration: {:?})",
        params.precursor_tolerance, params.precursor_cal_mode
    );
    log::info!(
        target: "run",
        "charge range   : {}-{}",
        params.charge_range.start(),
        params.charge_range.end()
    );
    log::info!(
        target: "run",
        "isotope errors : {}..={}",
        params.isotope_error_range.start(),
        params.isotope_error_range.end()
    );
    log::info!(
        target: "run",
        "decoy          : {:?} (prefix {})   chimeric: {}",
        decoy_strategy, cli.db.decoy_prefix, params.chimeric
    );
}

/// The names of the flags in `flags` that are set.
fn set_flags<'a>(flags: &[(bool, &'a str)]) -> Vec<&'a str> {
    flags
        .iter()
        .filter_map(|&(on, flag)| on.then_some(flag))
        .collect()
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

fn with_label_mods(
    aa: &model::AminoAcidSet,
    spec: &quant::label::LabelSpec,
) -> Result<model::AminoAcidSet, String> {
    let (mods, channels) = crate::quant::run::labels::search_label_mods(spec);
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

#[cfg(test)]
#[path = "resolve/tests/deuterated.rs"]
mod deuterated_tag_warning_tests;

#[cfg(test)]
#[path = "resolve/tests/retrieval.rs"]
mod retrieval_choice_tests;

#[cfg(test)]
#[path = "resolve/tests/destination.rs"]
mod rescore_pin_destination_tests;

#[cfg(test)]
#[path = "resolve/tests/filter.rs"]
mod isobaric_peak_filter_tests;

#[cfg(all(test, feature = "timstof"))]
#[path = "resolve/tests/archives.rs"]
mod d_archive_tests;

#[cfg(test)]
#[path = "resolve/tests/flags.rs"]
mod flag_combination_tests;
