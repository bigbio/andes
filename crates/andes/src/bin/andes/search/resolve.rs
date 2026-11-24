//! The parameter-resolution stages of `andes search`: flag checks, scratch paths, the database,
//! input detection, the scoring model, the search parameters and the resolved-parameter log.
//! Each stage returns a typed result that `run` destructures into its working locals.

use std::path::{Path, PathBuf};
use std::sync::Arc;

use super::{load_claim_database, scratch_tempdir, ClaimDatabase};
use crate::acquisition::Acquisition;
use crate::cli::{
    Cli, EnzymeSpecificity, EthcdActivationFlag, FragmentIndexFlag, Fragmentation,
    PrecursorMonoFlag, Protocol, ScoreFlag,
};
use crate::models::{
    cli_fragment_tol_override, cli_fragmentation_to_activation_opt, label_from_mods,
    load_param_from_store, parse_enzymes, resolve_metadataless_selection, routing_protocol,
    search_aa_set, warn_if_universal_protease_combo,
};
use crate::quant::run::QuantSettings;
use crate::spectra::{
    compare_with_first_input, detect_acquisition_metadata, detect_isobaric_sampled,
    input_format_flags, sample_ms2_meta, sample_spectra,
};
use crate::{arg_present, EXPLICIT_MISSED_CLEAVAGES};
use model::{
    activation::ActivationMethod, AminoAcidSetBuilder, InstrumentType, PrecursorTolerance,
};
use scoring_crate::{Param, RankScorer};
use search::precursor::monoisotopic::MonoParams;
use search::{PrecursorCalMode, SearchIndex, SearchParams};

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

/// Rejects flag combinations that cannot work, using the flags alone.
pub(super) fn check_flag_combinations(cli: &Cli) -> Result<(), String> {
    if cli.top_n == 0 {
        return Err("--top-n must be at least 1".into());
    }
    if cli.dia_mobility_slices == 0 {
        return Err("--dia-mobility-slices must be at least 1".into());
    }
    if !cli.dia_ms1_tol_ppm.is_finite() || cli.dia_ms1_tol_ppm <= 0.0 {
        return Err("--dia-ms1-tol-ppm must be a positive number".into());
    }
    if cli.localize.is_some() && !cli.rescore {
        return Err("--localize ipf needs --rescore".into());
    }
    if cli.xl.is_none() && !cli.external_percolator() {
        if !cli.percolator_args.trim().is_empty() {
            return Err(
                "native rescoring takes no percolator_args; add --percolator-bin or --percolator-docker to run an external Percolator"
                    .into(),
            );
        }
        if cli.percolator_image != output::DEFAULT_PERCOLATOR_IMAGE {
            return Err(
                "--percolator-image has no effect without --percolator-docker; add --percolator-docker to run Percolator from that image"
                    .into(),
            );
        }
    }
    if cli.dia {
        if cli.chimeric {
            return Err("--dia and --chimeric cannot be combined: --dia searches every scan across its whole window".into());
        }
        if cli.open_search {
            return Err("--dia and --open-search cannot be combined: a DIA scan has no selected precursor to offset from".into());
        }
        if cli.refine {
            return Err(
                "--dia and --refine cannot be combined: the refinement pass has no MS1 gate".into(),
            );
        }
        if cli.glyco {
            return Err("--dia and --glyco cannot be combined: the glycopeptide search has no isolation-window or MS1 gate".into());
        }
        if cli.xl.is_some() {
            return Err("--dia and --xl cannot be combined: the cross-link search has no isolation-window or MS1 gate".into());
        }
    }
    if cli.glyco_index_sequon_only && !cli.glyco {
        return Err(
            "--glyco-index-sequon-only requires --glyco (the standard peptide \
                    search scores every candidate)"
                .into(),
        );
    }
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
    }
    if cli.xl.is_some() {
        let unsupported = set_flags(&[
            (cli.output_tsv.is_some(), "--output-tsv"),
            (cli.output_parquet.is_some(), "--output-parquet"),
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
        ]);
        if !unsupported.is_empty() {
            return Err(format!(
                "--xl does not support: {} (crosslink mode writes its own .xl.pin and \
                 tables, scored by the pair score)",
                unsupported.join(", ")
            ));
        }
    }
    if cli.glyco {
        let unsupported = set_flags(&[
            (cli.output_tsv.is_some(), "--output-tsv"),
            (cli.output_parquet.is_some(), "--output-parquet"),
            (cli.rescore, "--rescore"),
            (cli.refine, "--refine"),
            (cli.psm_fdr.is_some(), "--psm-fdr"),
            (cli.pep.is_some(), "--pep"),
        ]);
        if !unsupported.is_empty() {
            return Err(format!(
                "--glyco does not support: {} (glyco mode writes a standalone \
                 .glyco.pin and skips the standard PIN/rescore/TSV/Parquet/refine \
                 pipeline; run Percolator on the .glyco.pin separately)",
                unsupported.join(", ")
            ));
        }
        if cli.mbr && !cli.lfq {
            return Err(
                "--glyco --mbr needs --lfq: match-between-runs transfers the identified intact \
                 glycopeptides into the other runs"
                    .into(),
            );
        }
        if cli.optimize {
            return Err("--optimize is not supported with --glyco".into());
        }
    }
    Ok(())
}

/// Creates the scratch directory, unpacks archived inputs, validates the quantification
/// arguments and decides where the PIN goes.
pub(super) fn prepare_run(cli: &mut Cli) -> Result<RunPaths, Box<dyn std::error::Error>> {
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
    let _d_archives = unpack_d_archives(&mut cli.spectrum, &scratch_dir)?;
    let mut _mbr_tmp: Option<tempfile::TempDir> = None;
    if cli.mbr && !cli.glyco {
        check_mbr_args(cli)?;
        crate::quant::run::check_unique_stems(&cli.spectrum)?;
        if cli.quant_dir.is_none() && cli.label.is_none() {
            let dir = scratch_tempdir("andes-quant-", &scratch_dir)?;
            cli.quant_dir = Some(dir.path().to_path_buf());
            _mbr_tmp = Some(dir);
        }
    }
    if cli.quant_dir.is_some() {
        crate::quant::run::check_unique_stems(&cli.spectrum)?;
    }
    QuantSettings::check_args(cli)?;
    if let Some(warning) = deuterated_tag_warning(cli.tmt, cli.precursor_tol) {
        warn_search(&warning);
    }
    let mut _rescore_tmp: Option<tempfile::TempDir> = None;
    let (output_pin_path, report_base) = rescore_pin_destination(
        cli.output_pin.as_deref(),
        cli.output_tsv.as_deref(),
        cli.output_parquet.as_deref(),
        cli.keep_pin,
        &scratch_dir,
        &mut _rescore_tmp,
    )?;
    if cli.spectrum.len() > 1 {
        log::info!(target: "search",
            "Multi-spectrum search: {} inputs → one PIN",
            cli.spectrum.len()
        );
    }
    Ok(RunPaths {
        scratch_dir,
        original_spectrum_paths,
        output_pin_path,
        report_base,
        temp_dirs: TempDirs {
            _rescore: _rescore_tmp,
            _mbr: _mbr_tmp,
            _d_archives,
        },
    })
}

/// Missed cleavages per peptide: --max-missed-cleavages, raised to 3 under --glyco unless it was
/// given explicitly.
fn max_missed_cleavages(cli: &Cli) -> u32 {
    let explicit = *EXPLICIT_MISSED_CLEAVAGES.get().unwrap_or(&false);
    if cli.glyco && !explicit {
        cli.max_missed_cleavages.max(3)
    } else {
        cli.max_missed_cleavages
    }
}

/// Parses the enzymes, loads the protein database and builds the search index.
pub(super) fn load_search_database(
    cli: &Cli,
    database_path: &Path,
) -> Result<SearchDatabase, Box<dyn std::error::Error>> {
    let (search_enzyme, extra_enzymes) = parse_enzymes(&cli.db.enzyme)?;
    warn_if_universal_protease_combo(search_enzyme, &extra_enzymes);
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
    } = load_claim_database(database_path, &cli.db.annotations, &uniprot_opts)?;
    let target_db = loaded.db;
    if let Some(peff) = &loaded.peff {
        log_claims("PEFF", &peff.claims);
        if peff.sites_without_mass + peff.processing_not_chain > 0 {
            log::info!(target: "search",
                "PEFF: {} entries; left out {} sites with no known mass and {} processed spans \
                 that are not mature chains",
                peff.entries,
                peff.sites_without_mass,
                peff.processing_not_chain
            );
        }
    }
    if let Some(stats) = &annotations {
        log_claims("Annotations", stats);
    }
    if site_knockoffs.knockoffs > 0 {
        log::info!(
            target: "site-claims",
            "{} annotated sites, {} knockoff sites ({} sites with no position for \
             one)",
            site_knockoffs.sites, site_knockoffs.knockoffs, site_knockoffs.untested
        );
    }
    log::info!(
        target: "database",
        "Loaded {} target proteins from {} in {:.2} s",
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
                    log::warn!(target: "search",
                        "{flag} needs a UniProt XML database (*.xml[.gz]); {} is FASTA, \
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
        min_length: cli.min_length,
        max_length: cli.max_length,
        max_missed_cleavages: max_missed_cleavages(cli),
    };
    let prebuilt = crate::database::prebuilt_settings(
        database_path,
        &crate::database::SearchDecoys {
            strategy: decoy_strategy,
            seed: cli.decoy_seed,
            prefix: &cli.db.decoy_prefix,
            suffix: cli.db.decoy_suffix.as_deref(),
            digest: &proteases,
        },
    )?;
    if let Some(settings) = prebuilt {
        if !cli.db.annotations.is_empty() {
            return Err(format!(
                "--annotations cannot be used with {}, a database written by andes database:                  annotated sites change the decoys. Search the target FASTA instead.",
                database_path.display()
            )
            .into());
        }
        let t_phase = std::time::Instant::now();
        let idx = SearchIndex::from_prebuilt(
            target_db,
            &cli.db.decoy_prefix,
            cli.db.decoy_suffix.as_deref(),
        );
        let target_db = idx.target_db();
        let decoys = idx.db.proteins.len() - target_db.proteins.len();
        if (target_db.proteins.len() as u64, decoys as u64) != (settings.targets, settings.decoys)
        {
            return Err(format!(
                "{}: {} targets and {decoys} decoys, but andes database wrote {} and {}",
                database_path.display(),
                target_db.proteins.len(),
                settings.targets,
                settings.decoys
            )
            .into());
        }
        log::info!(target: "database",
            "{} targets and {decoys} decoys ({}) written by andes database, searched as written; \
             decoy key {}",
            target_db.proteins.len(),
            settings.decoy_strategy,
            settings.decoy_key
        );
        log::info!(
            target: "search_index_build",
            "done in {:.2} s",
            t_phase.elapsed().as_secs_f64()
        );
        return Ok(SearchDatabase {
            target_db,
            idx,
            search_enzyme,
            extra_enzymes,
            proteases,
            decoy_strategy,
        });
    }
    if decoy_strategy == search::decoy::DecoyStrategy::Mutate
        && (proteases.is_non_specific() || proteases.num_tolerable_termini < 2)
    {
        return Err(
            "--decoy-strategy mutate needs fully specific digestion: a semi- or non-specific \
             search reads sub-sequences of each decoy segment, and the ones without the changed \
             residue are also target peptides. Use --decoy-strategy reverse or shuffle."
                .into(),
        );
    }
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
            warn_search(&warning);
        }
        let suffix = cli.db.decoy_suffix.as_deref().filter(|s| !s.is_empty());
        log::info!(target: "database",
            "decoy key {} (the key of these decoys; andes database writes the same decoys \
             for these settings)",
            search::decoy::decoy_key(
                &search::decoy::targets_sha256(&target_db, &cli.db.decoy_prefix, suffix),
                decoy_strategy,
                cli.decoy_seed,
                &cli.db.decoy_prefix,
                suffix,
                std::slice::from_ref(&proteases),
            )
        );
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
    log::info!(
        target: "search_index_build",
        "done in {:.2} s",
        t_phase.elapsed().as_secs_f64()
    );
    Ok(SearchDatabase {
        target_db,
        idx,
        search_enzyme,
        extra_enzymes,
        proteases,
        decoy_strategy,
    })
}

/// Detects the input format, the acquisition metadata and any isobaric label.
pub(super) fn detect_inputs(
    cli: &Cli,
    aa: &model::AminoAcidSet,
) -> Result<RunInputs, Box<dyn std::error::Error>> {
    let spectrum_paths = &cli.spectrum;
    let spectrum_path: PathBuf = spectrum_paths[0].clone();
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
            warn_search(&warning);
        }
    }

    let peak_evidence = if !instrument_was_detected
        && cli.fragment_tol_ppm.is_none()
        && cli.fragment_tol_da.is_none()
    {
        let ev = sample_spectra(&spectrum_path, is_mzml, is_mgf, is_raw, cli.ms_level as u32)
            .map(|sample| input::infer_fragment_resolution(&sample));
        if let Some(e) = &ev {
            log::info!(
                target: "acquisition",
                "no analyzer metadata; fragment resolution from the peaks: {} -> {}",
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
            log::warn!(target: "search",
                "the SDRF dissociation method is {} but {} records {}; andes follows the file",
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
            label_from_mods(aa)
        } else {
            None
        }
    });
    let route_protocol = routing_protocol(cli.protocol, label);
    Ok(RunInputs {
        is_mzml,
        is_raw,
        is_d,
        is_mgf,
        detected_activation,
        detected_instrument,
        instrument_was_detected,
        acquisition,
        label,
        route_protocol,
    })
}

/// Chooses the scoring model for the detected activation and instrument, and builds the scorer.
pub(super) fn resolve_scoring_model(
    cli: &mut Cli,
    inputs: &RunInputs,
    aa: &mut model::AminoAcidSet,
    search_enzyme: model::enzyme::Enzyme,
) -> Result<ScoringModel, Box<dyn std::error::Error>> {
    let spectrum_path: PathBuf = cli.spectrum[0].clone();
    let detected_activation = inputs.detected_activation;
    let detected_instrument = inputs.detected_instrument;
    let instrument_was_detected = inputs.instrument_was_detected;
    let acquisition = &inputs.acquisition;
    let label = inputs.label;
    let route_protocol = inputs.route_protocol;
    let t_phase = std::time::Instant::now();
    let (mut param, routed_activation, routed_instrument, resolved_model_id) = {
        let detected = (detected_activation, detected_instrument);
        let (activation, instrument_opt): (ActivationMethod, Option<InstrumentType>) =
            match detected {
                (Some(method), Some(inst)) => {
                    let chosen =
                        cli_fragmentation_to_activation_opt(cli.fragmentation).unwrap_or(method);
                    if chosen != method {
                        log::info!(
                            target: "model",
                            "auto-detected activation = {} (instrument = {}) from {}, \
                             OVERRIDDEN by --fragmentation {}",
                            method.name(), inst.name(), spectrum_path.display(), chosen.name()
                        );
                    } else {
                        log::info!(
                            target: "model",
                            "auto-detected activation = {} (instrument = {}) from {}",
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
                    log::info!(
                        target: "model",
                        "instrument = {} from {}, no activation recorded; \
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
        log::info!(target: "model", "{model_id} (from store)");
        (p, activation, instrument_opt, model_id)
    };
    let gbdt_k = cli.gbdt_max_trees.unwrap_or(if cli.glyco {
        0
    } else {
        crate::cli::DEFAULT_GBDT_MAX_TREES
    });
    if gbdt_k > 0 {
        for (slot, name) in [
            (&mut param.frag_intensity_model, "frag-intensity"),
            (&mut param.rich_ion_model, "rich-ion"),
        ] {
            if let Some(before) = truncate_gbdt(slot, gbdt_k) {
                log::info!(target: "search",
                    "--gbdt-max-trees {gbdt_k}: {name} ensemble {before} -> {gbdt_k} trees \
                     (changes predictions; not byte-identical)"
                );
            }
        }
    }
    match cli.protocol {
        Protocol::Tmt => param.data_type.protocol = model::protocol::Protocol::TMT,
        Protocol::Itraq => param.data_type.protocol = model::protocol::Protocol::ITRAQ,
        Protocol::ItraqPhospho => {
            param.data_type.protocol = model::protocol::Protocol::ITRAQPhospho
        }
        Protocol::Auto => {
            if let Some(l) = label {
                log::info!(
                    target: "protocol",
                    "{} {l:?} → isobaric protocol (windowed peak filter on low-res MS2)",
                    if cli.tmt.is_some() {
                        "--tmt"
                    } else {
                        "detected"
                    }
                );
                param.data_type.protocol = match l {
                    input::IsobaricLabel::Tmt | input::IsobaricLabel::TmtPro => {
                        model::protocol::Protocol::TMT
                    }
                    input::IsobaricLabel::Itraq | input::IsobaricLabel::Itraq8 => {
                        model::protocol::Protocol::ITRAQ
                    }
                };
            }
        }
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
            *aa = search_aa_set(&cli.db, Some((name, mass)))?;
            log::info!(
                target: "protocol",
                "injected {name} fixed mod (+{mass:.4} on K + peptide N-term) \
                 into the candidate set (no --fixed-mods given)"
            );
        }
    }
    let mut scorer = RankScorer::new(&param);
    scorer.set_acquisition_high_res(Some(acquisition.high_res));
    let dia_tol_default = cli.dia
        && acquisition.analyzer == Some(InstrumentType::OrbitrapAstral)
        && cli.fragment_tol_ppm.is_none()
        && cli.fragment_tol_da.is_none();
    if dia_tol_default {
        cli.fragment_tol_ppm = Some(crate::dia::FRAGMENT_TOL_PPM);
        log::info!(target: "search",
            "--dia: fragment tolerance {} ppm (Orbitrap Astral MS2)",
            crate::dia::FRAGMENT_TOL_PPM
        );
    }
    let frag_tol_override = cli_fragment_tol_override(cli.fragment_tol_ppm, cli.fragment_tol_da);
    if let Some(t) = frag_tol_override {
        if instrument_was_detected && !dia_tol_default {
            log::info!(target: "search",
                "--fragment-tol-* {t:?} overrides the feature tolerance of the detected {} analyzer",
                acquisition.analyzer.map_or("", |i| i.name())
            );
        }
        scorer.set_fragment_tol_override(frag_tol_override);
    }
    log::info!(
        target: "param_and_scorer",
        "done in {:.2} s",
        t_phase.elapsed().as_secs_f64()
    );
    Ok(ScoringModel {
        param,
        scorer,
        routed_activation,
        routed_instrument,
        resolved_model_id,
        gbdt_k,
        frag_tol_override,
    })
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
            let report = match search::fragments::index::vote_peaks_for(fragment_tol) {
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

/// Logs a pre-formatted warning under the `search` target without its `WARN: ` prefix.
fn warn_search(warning: &str) {
    log::warn!(target: "search", "{}", warning.strip_prefix("WARN: ").unwrap_or(warning));
}

/// The names of the flags in `flags` that are set.
fn set_flags<'a>(flags: &[(bool, &'a str)]) -> Vec<&'a str> {
    flags
        .iter()
        .filter_map(|&(on, flag)| on.then_some(flag))
        .collect()
}

/// Shortens the GBDT ensemble in `slot` to `k` trees and returns its former tree count, or
/// `None` when it is absent or already no longer than `k`.
pub(super) fn truncate_gbdt(
    slot: &mut Option<Arc<scoring_crate::gbdt::GbdtPeakModel>>,
    k: usize,
) -> Option<usize> {
    let g = slot.as_ref()?;
    let before = g.trees.len();
    if k >= before {
        return None;
    }
    let mut t = (**g).clone();
    t.trees.truncate(k);
    *slot = Some(Arc::new(t));
    Some(before)
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
            log::info!(target: "search",
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
    log::info!(
        target: source,
        "{} rows{files}: {} variant, {} knockoff and {} chain entries, {} sites on \
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
        log::info!(
            target: source,
            "{} single-residue variants with no position for a knockoff",
            stats.variants_without_knockoff
        );
    }
    if stats.knockoffs_dropped > 0 {
        log::info!(
            target: source,
            "{} earlier knockoff entries dropped, within the knockoff distance of a \
             variant these files list",
            stats.knockoffs_dropped
        );
    }
}

fn log_uniprot_load(stats: &input::UniProtLoadStats, opts: &input::UniProtOptions) {
    log::info!(
        target: "uniprot",
        "{} entries, {} annotated modification sites, {} variant entries \
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
        log::info!(
            target: "uniprot",
            "{} knockoff variant entries ({} single-residue variants with no \
             position for a knockoff)",
            stats.knockoff_proteins, stats.variants_without_knockoff,
        );
    }
    if opts.isoforms || opts.chains {
        log::info!(
            target: "uniprot",
            "{} isoform entries ({} isoforms not rebuildable), {} mature-chain entries",
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
        log::info!(
            target: "uniprot",
            "{total} modification features have no single known mass and are not \
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
            log::warn!(target: "search",
                "{} ({:+.4} Da) is restricted to annotated sites while \
                 --uniprot-variants is on: variants of the same mass lose the \
                 competing explanation at unannotated residues.",
                m.name,
                m.mass_delta
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
            log::warn!(target: "search",
                "{name} ({mass:+.4} Da) is restricted to annotated sites (restrict_to: annotated), \
                 but no target protein annotates a site of that mass, so it is never placed. \
                 Annotated sites come from a UniProt XML database (--database *.xml[.gz]); \
                 --uniprot-evidence sets which are used."
            );
        } else {
            log::info!(target: "search",
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
    if !cli.rescore {
        return Err(
            "--mbr needs --rescore: the library ranks precursors by the rescored PSM scores".into(),
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
    crate::quant::evidence_params(cli.lfq_tol)?;
    Ok(())
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
            assert!(search::fragments::index::covers_activation(a), "{a:?}");
        }
        for a in [ActivationMethod::ETD, ActivationMethod::UVPD] {
            assert!(!search::fragments::index::covers_activation(a), "{a:?}");
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

#[cfg(test)]
mod flag_combination_tests {
    use super::*;
    use clap::Parser;

    #[derive(Parser, Debug)]
    struct Wrap {
        #[command(flatten)]
        args: Cli,
    }

    fn args(extra: &[&str]) -> Cli {
        let argv = std::iter::once("andes").chain(extra.iter().copied());
        Wrap::try_parse_from(argv).unwrap().args
    }

    fn accepted(extra: &[&str]) {
        if let Err(e) = check_flag_combinations(&args(extra)) {
            panic!("{extra:?} was rejected: {e}");
        }
    }

    fn rejected(extra: &[&str], fragment: &str) {
        match check_flag_combinations(&args(extra)) {
            Ok(()) => panic!("{extra:?} was accepted"),
            Err(e) => assert!(e.contains(fragment), "{extra:?}: {e}"),
        }
    }

    #[test]
    fn typical_runs_are_accepted() {
        accepted(&[]);
        accepted(&["--rescore"]);
        accepted(&["--rescore", "--localize", "ipf"]);
        accepted(&["--rescore", "--chimeric", "--refine"]);
        accepted(&["--dia"]);
        accepted(&["--open-search"]);
        accepted(&["--lfq", "--mbr", "--rescore"]);
        accepted(&["--glyco", "--glyco-index-sequon-only"]);
        accepted(&["--glyco", "--lfq", "--mbr"]);
        accepted(&["--xl", "DSSO"]);
    }

    #[test]
    fn top_n_must_be_positive() {
        rejected(&["--top-n", "0"], "--top-n must be at least 1");
    }

    #[test]
    fn localize_needs_rescoring() {
        rejected(&["--localize", "ipf"], "--localize ipf needs --rescore");
    }

    #[test]
    fn rescore_alone_selects_the_native_engine() {
        let cli = args(&["--rescore"]);
        assert!(cli.rescore && !cli.external_percolator());
    }

    #[test]
    fn a_percolator_binary_or_docker_selects_the_external_percolator() {
        let cli = args(&["--rescore", "--percolator-bin", "/opt/percolator"]);
        assert!(cli.rescore && cli.external_percolator());
        accepted(&["--rescore", "--percolator-bin", "/opt/percolator"]);
        let cli = args(&["--rescore", "--percolator-docker"]);
        assert!(cli.rescore && cli.external_percolator());
        accepted(&[
            "--rescore",
            "--percolator-docker",
            "--percolator-image",
            "other:1",
        ]);
        accepted(&[
            "--rescore",
            "--percolator-bin",
            "/opt/percolator",
            "--percolator-args=--maxiter 3",
        ]);
    }

    #[test]
    fn percolator_options_need_an_external_percolator() {
        rejected(
            &["--rescore", "--percolator-args=--maxiter 3"],
            "native rescoring takes no percolator_args",
        );
        rejected(
            &["--percolator-args=--maxiter 3"],
            "native rescoring takes no percolator_args",
        );
        rejected(
            &["--rescore", "--percolator-image", "other:1"],
            "--percolator-image has no effect without --percolator-docker",
        );
        accepted(&["--xl", "DSSO", "--percolator-args=--maxiter 3"]);
    }

    #[test]
    fn rescore_native_is_an_unknown_flag() {
        let argv = ["andes", "--rescore-native"];
        let err = Wrap::try_parse_from(argv).unwrap_err();
        assert_eq!(err.kind(), clap::error::ErrorKind::UnknownArgument);
    }

    #[test]
    fn dia_is_not_combined_with_other_modes() {
        rejected(&["--dia", "--chimeric"], "--dia and --chimeric");
        rejected(&["--dia", "--open-search"], "--dia and --open-search");
        rejected(&["--dia", "--refine"], "--dia and --refine");
        rejected(&["--dia", "--glyco"], "--dia and --glyco");
        rejected(&["--dia", "--xl", "DSSO"], "--dia and --xl");
    }

    #[test]
    fn sequon_only_needs_glyco() {
        rejected(
            &["--glyco-index-sequon-only"],
            "--glyco-index-sequon-only requires --glyco",
        );
    }

    #[test]
    fn open_search_is_not_combined_with_other_modes() {
        rejected(&["--open-search", "--glyco"], "not supported with --glyco");
        rejected(&["--open-search", "--chimeric"], "with --chimeric");
        rejected(&["--open-search", "--label", "x"], "with --label");
        rejected(&["--open-search", "--tmt", "tmt10plex"], "with --tmt");
        rejected(&["--open-search", "--lfq"], "with --lfq or --mbr");
        rejected(&["--open-search", "--mbr"], "with --lfq or --mbr");
        rejected(&["--open-search", "--refine"], "with --refine");
        rejected(&["--open-search", "--xl", "DSSO"], "with --xl");
        rejected(
            &["--open-search", "--isotope-error", "0..2"],
            "remove --isotope-error",
        );
    }

    #[test]
    fn xl_rejects_the_flags_it_does_not_support() {
        rejected(
            &[
                "--xl",
                "DSSO",
                "--output-tsv",
                "out.tsv",
                "--lfq",
                "--refine",
            ],
            "--xl does not support: --output-tsv, --refine, --lfq",
        );
        rejected(
            &["--xl", "DSSO", "--tmt", "tmt10plex"],
            "--xl does not support: --tmt",
        );
    }

    #[test]
    fn glyco_rejects_the_flags_it_does_not_support() {
        rejected(
            &[
                "--glyco",
                "--output-parquet",
                "out",
                "--rescore",
                "--pep",
                "0.01",
            ],
            "--glyco does not support: --output-parquet, --rescore, --pep",
        );
        rejected(
            &["--glyco", "--refine"],
            "--glyco does not support: --refine",
        );
    }

    #[test]
    fn glyco_match_between_runs_needs_lfq() {
        rejected(&["--glyco", "--mbr"], "--glyco --mbr needs --lfq");
    }

    #[test]
    fn optimize_is_not_supported_with_glyco() {
        rejected(
            &["--glyco", "--optimize"],
            "--optimize is not supported with --glyco",
        );
    }
}
