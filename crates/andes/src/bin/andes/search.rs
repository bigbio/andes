//! The standard database-search driver.

use std::fs::File;
use std::path::{Path, PathBuf};
use std::sync::mpsc::sync_channel;
use std::sync::Arc;
use std::thread;

use crate::cli::{
    CandidateIndexFlag, Cli, EnzymeSpecificity, EthcdActivationFlag, Fragmentation,
    GlycoIsotopeFlag, Protocol, ScoreFlag,
};
use crate::glyco_run::run_glyco;
use crate::model_select::{
    cli_fragment_tol_override, default_aa_set_with_tag, load_param_from_store, parse_enzymes,
    resolve_metadataless_selection, warn_if_universal_protease_combo,
};
use crate::rescore;
use crate::spectra::{
    detect_dominant_activation, detect_instrument_type_for_path, detect_isobaric_sampled,
    input_format_flags, merge_parse_stats, prefix_spectrum_titles, run_precursor_calibration,
    send_chunks, title_prefix_for, tolerance_ppm_display, warn_if_index_will_not_fit, ParseStats,
};
use crate::{
    arg_present, available_memory_bytes, log_rss, report_search_progress,
    EXPLICIT_MISSED_CLEAVAGES, RSS_PROBE,
};
use input::{FastaReader, MgfReader, Ms1Link, MzMLReader};
use model::{
    activation::ActivationMethod, AminoAcidSetBuilder, InstrumentType, PrecursorTolerance, Spectrum,
};
use scoring_crate::RankScorer;
use search::candidate_index::index_cache_path;
use search::{
    apply_shift_for_mode, apply_tightened_precursor_tolerance, PrecursorCalMode, PreparedSearch,
    SearchIndex, SearchParams, TopNQueue,
};

pub(crate) fn run(cli: Cli) -> Result<(), Box<dyn std::error::Error>> {
    // These three were validated as Some(..) by main() before calling run().
    if cli.spectrum.is_empty() {
        return Err("no --spectrum inputs".into());
    }
    // Parse the digestion enzyme(s) once. `--enzyme` accepts a comma-separated
    // list for a multi-protease digest; the FIRST entry is the primary (drives
    // model selection via build_selection_key + the cleavage-credit feature),
    // the rest widen candidate enumeration (params.extra_enzymes). The common
    // single-enzyme case yields an empty extras list ⇒ digestion bit-identical.
    let (search_enzyme, extra_enzymes) = parse_enzymes(&cli.enzyme)?;
    warn_if_universal_protease_combo(search_enzyme, &extra_enzymes);
    let spectrum_paths = &cli.spectrum;
    let spectrum_path: PathBuf = spectrum_paths[0].clone();
    let database_path: PathBuf = cli.database.clone().expect("database validated in main");
    // PIN destination. Normally `--output-pin`. Under `--rescore` without an
    // explicit path we write to a temp dir (kept alive in `_rescore_tmp` until
    // the rescore phase has parsed it). `--keep-pin false` also routes the PIN
    // to a temp dir so it is removed when the dir drops at function exit.
    let mut _rescore_tmp: Option<tempfile::TempDir> = None;
    let output_pin_path: PathBuf = match cli.output_pin.clone() {
        Some(p) if cli.keep_pin || !cli.rescore => p,
        // --rescore with --keep-pin false but an explicit path: still honor the
        // explicit path (user asked for it); keep-pin only governs the temp case.
        Some(p) => p,
        None => {
            // --rescore without --output-pin → temp PIN.
            let dir = tempfile::tempdir()
                .map_err(|e| format!("creating temp dir for rescore PIN: {e}"))?;
            let path = dir.path().join("andes.pin");
            _rescore_tmp = Some(dir);
            path
        }
    };
    // DURABLE base for run artifacts that must survive function exit — the
    // filtered q-value TSV and `statistics.log`. When the PIN is a real path
    // these sit next to it; but when the PIN is a temp file (`--rescore`/`--fdr`
    // without `--output-pin`) deriving them from `output_pin_path` would write
    // them INTO the temp dir, which is deleted at exit — so the user would be
    // told the TSV was written and then find it gone. Fall back to an explicit
    // `--output-tsv`/`--output-parquet` location, else the current directory.
    let report_base: PathBuf = if _rescore_tmp.is_some() {
        cli.output_tsv
            .clone()
            .or_else(|| cli.output_parquet.clone())
            .unwrap_or_else(|| PathBuf::from("andes.pin"))
    } else {
        output_pin_path.clone()
    };
    if spectrum_paths.len() > 1 {
        eprintln!(
            "Multi-spectrum search: {} inputs → one PIN",
            spectrum_paths.len()
        );
    }

    let _ = RSS_PROBE.set(cli.rss_probe);
    log_rss("startup");
    let t_total = std::time::Instant::now();
    let t_phase = std::time::Instant::now();
    // ── 1. Load FASTA target database ────────────────────────────────────────
    let target_db = FastaReader::load_all(input::open_buf_maybe_gz(&database_path)?)?;
    eprintln!(
        "Loaded {} target proteins from {} [PHASE fasta_load: {:.2}s]",
        target_db.proteins.len(),
        database_path.display(),
        t_phase.elapsed().as_secs_f64()
    );
    log_rss("after_fasta_load");

    // ── 2. Build SearchIndex (targets + strategy-generated decoys) ────────────
    let decoy_strategy =
        search::decoy::DecoyStrategy::from_name(&cli.decoy_strategy).ok_or_else(|| {
            format!(
                "unknown --decoy-strategy '{}' (expected reverse/shuffle/sequon-reverse/none)",
                cli.decoy_strategy
            )
        })?;
    let t_phase = std::time::Instant::now();
    let idx = SearchIndex::from_target_db_with_strategy(
        &target_db,
        &cli.decoy_prefix,
        cli.decoy_suffix.as_deref(),
        decoy_strategy,
        cli.decoy_seed,
    );
    eprintln!(
        "[PHASE search_index_build: {:.2}s]",
        t_phase.elapsed().as_secs_f64()
    );
    log_rss("after_search_index_build");

    // ── 3. Build AminoAcidSet ────────────────────────────────────────────────
    //
    // If --mod is given, parse the mods.txt file. Otherwise
    // fall back to andes's historical defaults (CAM fixed on C,
    // Oxidation variable on M + protein-N-term Acetyl variable).
    //
    // `num_mods_from_file` is populated only when --mod is given and the
    // file contains a `NumMods=N` line; it overrides the default
    // `max_variable_mods_per_peptide` (3) below.
    let (mut aa, num_mods_from_file) = match &cli.mods {
        Some(path) => {
            let n = AminoAcidSetBuilder::parse_num_mods_from_file(path)
                .map_err(|e| format!("parsing NumMods= from {}: {e}", path.display()))?;
            let set = AminoAcidSetBuilder::new_standard()
                .add_mods_from_file(path)
                .map_err(|e| format!("loading mods from {}: {e}", path.display()))?
                .build()
                .map_err(|e| format!("building amino-acid set from {}: {e}", path.display()))?;
            eprintln!(
                "Loaded modifications from {} (NumMods={})",
                path.display(),
                n.map(|v| v.to_string()).unwrap_or_else(|| "default".into()),
            );
            (set, n)
        }
        // No --mods: andes defaults (CAM-C fixed, Ox-M variable). The isobaric
        // tag (TMT/iTRAQ) is injected later, after protocol detection (C1).
        None => (default_aa_set_with_tag(None)?, None),
    };

    // ── 4. Load Param scoring model ───────────────────────────────────────────
    //
    // The model is selected from the canonical Parquet store: for mzML/.raw/.d
    // the activation+analyzer are auto-detected from metadata; for MGF
    // (metadata-less) the `--fragmentation` / `--fragment-tol-*` flags drive
    // `resolve_metadataless_selection`.
    // Detect the format from the underlying extension, stripping a trailing
    // `.gz` first so `spectra.mzML.gz` routes as mzML (not mis-routed to MGF by
    // the bare `.gz`); input::open_maybe_gz then reads it transparently. Native
    // `.raw` (Thermo, .NET 8) and `.d` (Bruker timsTOF, pure Rust) are
    // binary/directory and never gzipped. Anything else is treated as MGF.
    let (is_mzml, is_raw, is_d, is_mgf) = input_format_flags(&spectrum_path);

    // Detect (activation, instrument) from the input for auto-routing.
    // mzML peeks the file; Thermo `.raw` reads vendor metadata; Bruker `.d`
    // is always CID/TimsTOF (DDA-PASEF). Detection runs for every metadata-
    // bearing format and always wins over the MGF-only `--fragmentation` /
    // `--fragment-tol-*` flags (which carry no metadata of their own).
    let auto_route_eligible = is_mzml || is_raw || is_d;
    let detected_activation_instrument: Option<(ActivationMethod, Option<InstrumentType>)> =
        if !auto_route_eligible {
            None
        } else if is_mzml {
            detect_dominant_activation(&spectrum_path)
                .map(|m| (m, detect_instrument_type_for_path(&spectrum_path)))
        } else if is_raw {
            #[cfg(feature = "thermo")]
            {
                input::thermo::detect_activation_instrument(&spectrum_path, 64)
            }
            #[cfg(not(feature = "thermo"))]
            {
                None
            }
        } else {
            // is_d — timsTOF DDA-PASEF: CID fragmentation on a TOF analyzer.
            Some((ActivationMethod::CID, Some(InstrumentType::TimsTOF)))
        };
    // Pre-compute before the routing match consumes `detected_activation_instrument`.
    let instrument_was_detected = detected_activation_instrument
        .map(|(_, inst)| inst.is_some())
        .unwrap_or(false);

    let t_phase = std::time::Instant::now();
    let mut param = {
        // ── Resolve (activation, instrument) for the Parquet model store. ─────
        //
        // Metadata-first precedence: a fully detected (activation, instrument)
        // wins outright. When only the activation method is detected (analyzer
        // unknown), or nothing is detected (MGF / metadata-less mzML/.raw), the
        // metadata-less resolver folds in the MGF-only `--fragmentation` and
        // `--fragment-tol-*` flags. Default when no metadata: CID / low-res
        // (cid_lowres_tryp).
        let (activation, instrument_opt): (ActivationMethod, Option<InstrumentType>) =
            match detected_activation_instrument {
                // An EXPLICIT --fragmentation must win over detection. It previously
                // only applied when detection returned no instrument, so on any normal
                // mzML it was silently discarded — including the EThcD case, where the
                // reader relabels ETD to HCD and the warning tells the user to "pass
                // --fragmentation to override". It could not override anything.
                Some((method, Some(inst))) => {
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
                Some((method, None)) => resolve_metadataless_selection(
                    Some(method),
                    cli.fragmentation,
                    cli.fragment_tol_ppm,
                    cli.fragment_tol_da,
                ),
                None => resolve_metadataless_selection(
                    None,
                    cli.fragmentation,
                    cli.fragment_tol_ppm,
                    cli.fragment_tol_da,
                ),
            };

        let (model_id, p) = load_param_from_store(
            activation,
            instrument_opt,
            cli.protocol,
            search_enzyme,
            cli.model_store.as_deref(),
            cli.model_id_override.as_deref(),
        )?;
        eprintln!("Param model: {model_id} (from store)");
        // E5/E12: loud-fail the silent enzyme fallback. If the user explicitly
        // chose a non-Trypsin enzyme but no enzyme-matching model exists for the
        // detected activation/instrument, selection backs off to a Trypsin/generic
        // model whose cleavage prior + PIN enzymatic features are for the wrong
        // protease. Warn unmissably (skip when --model pins the choice on purpose).
        if cli.model_id_override.is_none() && search_enzyme != model::enzyme::Enzyme::Trypsin {
            if let Some(selected) = p.data_type.enzyme {
                if selected != search_enzyme {
                    eprintln!(
                        "WARN: --enzyme {} has no matching model for the detected \
                         activation/instrument; fell back to '{model_id}' (enzyme {:?}). \
                         Scores will use the wrong protease's cleavage prior + PIN features. \
                         Train a matching model or pass --model to choose explicitly.",
                        search_enzyme.name(),
                        selected
                    );
                }
            }
        }
        p
    };
    // Optional GBDT truncation (`--gbdt-max-trees`). Applied to the loaded model
    // before any scoring so every code path sees the same ensemble.
    //
    // BOTH ensembles are truncated. The shipped models carry a frag-intensity
    // GBDT (~799 KB) *and* a rich-ion GBDT (~842 KB) of comparable size, each
    // walked per fragment per candidate, so truncating only the first would
    // leave about half the ensemble cost untouched.
    // Glyco keeps the full ensembles: the glyco path is candidate-generation bound
    // (the ensembles profile at ~0.3% of its wall time), so truncating there would
    // give up prediction fidelity for no speed. An explicit flag still applies.
    let gbdt_k = if cli.glyco && !arg_present("--gbdt-max-trees") {
        0
    } else {
        cli.gbdt_max_trees
    };
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
    // Stamp the requested isobaric protocol onto the loaded model so the dense-
    // spectrum windowed peak filter (ScoredSpectrum) engages on TMT/iTRAQ
    // searches even when model selection fell back to a non-isobaric table
    // (there is no bundled CID-TMT model, so `--protocol TMT` resolves to
    // `cid_lowres_tryp`, whose stored protocol is Standard).
    // An explicit `--protocol` wins outright. When left at `auto` (the default),
    // auto-detect TMT/iTRAQ from MS2 reporter ions (mzML/MGF) so the dense-peak
    // windowed filter engages with zero config — the same path `--protocol TMT`
    // takes today. Detection returns None for label-free data, so non-isobaric
    // runs stay byte-identical.
    match cli.protocol {
        Protocol::Tmt => param.data_type.protocol = model::protocol::Protocol::TMT,
        Protocol::Itraq => param.data_type.protocol = model::protocol::Protocol::ITRAQ,
        Protocol::ItraqPhospho => {
            param.data_type.protocol = model::protocol::Protocol::ITRAQPhospho
        }
        Protocol::Auto => {
            let high_res = param.data_type.instrument.is_high_resolution();
            match detect_isobaric_sampled(
                &spectrum_path,
                is_mzml,
                is_mgf,
                cli.ms_level as u32,
                high_res,
            ) {
                Some(input::IsobaricLabel::Tmt) => {
                    eprintln!("Protocol resolver: auto-detected TMT reporter ions → engaging isobaric windowed peak filter");
                    param.data_type.protocol = model::protocol::Protocol::TMT;
                }
                Some(input::IsobaricLabel::Itraq) => {
                    eprintln!("Protocol resolver: auto-detected iTRAQ reporter ions → engaging isobaric windowed peak filter");
                    param.data_type.protocol = model::protocol::Protocol::ITRAQ;
                }
                None => {}
            }
        }
        _ => {}
    }
    // C1: parameter-free path only (no explicit --mods). When the protocol
    // resolves to TMT/iTRAQ, inject the tag as a fixed mod on K + peptide
    // N-term so labeled peptides match their precursor mass — otherwise the
    // reporter filter engages but every labeled candidate is +tag Da off and
    // misses. With explicit --mods the user owns the mod set (they may already
    // supply the tag), so those runs stay byte-identical.
    if cli.mods.is_none() {
        let tag = match param.data_type.protocol {
            model::protocol::Protocol::TMT => Some(("TMT6plex", 229.162932_f64)),
            model::protocol::Protocol::ITRAQ | model::protocol::Protocol::ITRAQPhospho => {
                Some(("iTRAQ4plex", 144.102063_f64))
            }
            _ => None,
        };
        if let Some((name, mass)) = tag {
            aa = default_aa_set_with_tag(Some((name, mass)))?;
            eprintln!(
                "Protocol resolver: injected {name} fixed mod (+{mass:.4} on K + peptide N-term) \
                 into the candidate set (no --mods given)"
            );
        }
    }
    let mut scorer = RankScorer::new(&param);
    // Fragment-tol override applies to metadata-less (MGF) input only. For
    // mzML/.raw/.d the analyzer is auto-detected, so the override is ignored.
    let frag_tol_override = cli_fragment_tol_override(cli.fragment_tol_ppm, cli.fragment_tol_da);
    if frag_tol_override.is_some() {
        if instrument_was_detected {
            eprintln!("WARN: --fragment-tol-* ignored — instrument auto-detected from metadata (use --fragment-tol-ppm/-da with MGF input only)");
        } else {
            scorer.set_fragment_tol_override(frag_tol_override);
        }
    }
    eprintln!(
        "[PHASE param_and_scorer: {:.2}s]",
        t_phase.elapsed().as_secs_f64()
    );

    // ── 5. Build SearchParams ─────────────────────────────────────────────────
    let mut params = SearchParams::default_tryptic(aa);
    params.precursor_tolerance = PrecursorTolerance::symmetric(cli.precursor_tol);
    // Ranges are validated (min <= max) by the clap value parsers.
    let (charge_min, charge_max) = cli.charge;
    params.charge_range = charge_min..=charge_max;
    // Round-8: resolve the default by MODE, not by sniffing the value. An explicit
    // `--isotope-error` (any range, including -1..2) is always honoured verbatim.
    // Unset under --glyco defaults to 0..=2: an MS1 envelope audit found the firmware
    // mis-picks the monoisotopic peak only ever too HIGH (0 of 23,907 scans needed a
    // negative shift), while the iso=-1 arm emitted 28.5% of all candidate rows for
    // 0.29% of the correct answers at a ~53:47 target:decoy ratio - pure FDR dilution.
    // Dropping it measured +81 backbone-correct @1%. ANDES_GLYCO_ISO_NEG=1 restores it.
    let iso_default = if cli.glyco && cli.glyco_isotope_error != GlycoIsotopeFlag::Negative {
        (0, 2)
    } else {
        (-1, 2)
    };
    let (iso_min, iso_max) = cli.isotope_error.unwrap_or(iso_default);
    params.isotope_error_range = iso_min..=iso_max;
    // Glyco high-mass precursors (backbone + multi-kDa glycan) frequently have the
    // monoisotopic peak mis-picked several 13C low, so the true neutral mass falls
    // outside the default -1..=2 sweep. Widen the upper bound for glyco so that
    // candidate mass is reachable. A/B-gated: ANDES_GLYCO_ISO_WIDE only.
    if cli.glyco && cli.glyco_isotope_error == GlycoIsotopeFlag::Wide {
        params.isotope_error_range = iso_min..=iso_max.max(5);
    }
    // Pass 2 co-isolation requires MS1 scans, captured by the mzML and Thermo
    // `.raw` readers. MGF (no MS1) and the Bruker `.d` reader (DDA MS2 only;
    // chimeric on `.d` is out of scope) make `--chimeric` inert, so keep
    // `params.chimeric` FALSE to turn the ENTIRE chimeric path off (Pass 2, PIN
    // column/SpecId gates, top-N forcing) — the run is then identical to a normal search.
    let chimeric_active = cli.chimeric && (is_mzml || is_raw);
    if cli.chimeric && !(is_mzml || is_raw) {
        eprintln!(
            "WARN: --chimeric requires MS1 data (mzML or Thermo .raw); the input is {} \
             so the co-isolation cascade is disabled and the search runs normally.",
            if is_d {
                "Bruker .d (DDA MS2 only),"
            } else {
                "MGF,"
            }
        );
    }
    // The cascade pairs MS2 with its preceding MS1 — it is MS2-only by
    // construction. Ignore a non-2 `--ms-level` under `--chimeric` so MS3+
    // (e.g. TMT SPS-MS3) can never enter the search on any input format.
    if chimeric_active && cli.ms_level != 2 {
        eprintln!(
            "WARN: --ms-level={} is ignored under --chimeric; the cascade always searches MS2.",
            cli.ms_level
        );
    }
    params.chimeric = chimeric_active;
    // Fallback isolation half-width (Da) used only when the file omits per-scan
    // isolation offsets; a fixed sensible default (was the --isolation-halfwidth flag).
    params.chimeric_isolation_halfwidth_da = 1.5;
    params.chimeric_max_coisolated = cli.chimeric_max_coisolated;
    params.chimeric_max_kl = cli.chimeric_max_kl;
    // FORCE top-1 under the cascade: Pass 1 emits only the best primary per scan;
    // secondaries come from Pass 2. The default top_n (10) would make Pass 1 emit
    // blind multi-emission per scan = inflated FDR.
    params.top_n_psms_per_spectrum = if chimeric_active { 1 } else { cli.top_n };
    params.enzyme = search_enzyme;
    params.extra_enzymes = extra_enzymes.clone();
    params.num_tolerable_termini = match cli.enzyme_specificity {
        EnzymeSpecificity::Fully => 2,
        EnzymeSpecificity::Semi => 1,
        EnzymeSpecificity::NonSpecific => 0,
    };
    // Glycopeptides frequently carry >=2 missed cleavages (the glycan sterically
    // protects the nearby cleavage site), so the general default of 1 silently
    // drops ~6% of true glycopeptides from the digest — concentrated at high
    // charge. Validated: raising to 2 for --glyco gave z5 22->54 (+2.5x) and
    // +116 backbone-correct @1% on the pooled AI-ETD benchmark. Floor glyco at 2.
    // Round-8: raised again 2 -> 3. A sequon-bearing tryptic peptide often needs a
    // third missed cleavage to reach a length the glyco path can score, and the
    // reference identification set for this benchmark was itself produced with 3.
    // Validated on the pooled AI-ETD benchmark: +44 backbone-correct @1%, +91
    // glycoPSMs, +19 unique glycopeptides, +19 glycosites — and an ENTRAPMENT
    // measurement (yeast/E.coli, where a glyco ID is false by construction) puts the
    // true error at 0.14%, i.e. it buys IDs without spending error budget.
    // Install engine-wide scoring settings from validated CLI values. Done once, before
    // any spectrum is scored, so no scoring code reads the environment.
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
    scoring_crate::scoring::init_scoring_settings(scoring_crate::scoring::ScoringSettings {
        peak_filter,
        precursor_offset_clamp: cli.precursor_offset_clamp,
        density_on_active_list: cli.density_on_active_list,
    });
    input::mzml::init_ethcd_as_etd(cli.ethcd_activation == EthcdActivationFlag::Etd);
    params.chimeric_allow_overlap = cli.chimeric_allow_overlap;

    params.max_missed_cleavages = if cli.glyco {
        // The floor is a default, not a mandate. Raising it to 3 grows the candidate
        // index by ~40% (13.2M -> 18.8M candidates on a 20k-protein human FASTA,
        // +4.4 GB resident), which is the difference between running and being
        // OOM-killed on a 16 GB machine. A user who asked for fewer gets fewer.
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
    if let Some(n) = num_mods_from_file {
        params.max_variable_mods_per_peptide = n; // NumMods= in --mods overrides --max-mods
    }
    params.precursor_cal_mode = cli.precursor_cal;
    // params.cal_min_spec_keys keeps its SearchParams default
    // (MIN_SPECKEYS_FOR_PREPASS); it is an internal threshold, no longer a flag.
    params.precursor_mass_shift_ppm = 0.0;
    params.refine_select_psm_fdr = cli.refine_select_psm_fdr;
    params.score_mode = match cli.score {
        ScoreFlag::Rank => search::ScoreMode::Rank,
        ScoreFlag::Strong => search::ScoreMode::Strong,
        // Auto: strong re-ranking pays on high-res (sharp fragment peaks → the
        // intensity / rich-ion GBDTs separate target/decoy); on low-res it
        // collapses the candidate pool → use rank. Keyed off the resolved
        // model's instrument resolution.
        ScoreFlag::Auto => {
            if param.data_type.instrument.is_high_resolution() {
                search::ScoreMode::Strong
            } else {
                search::ScoreMode::Rank
            }
        }
    };
    params.candidate_index = match cli.candidate_index {
        CandidateIndexFlag::Ram => search::CandidateIndexMode::Ram,
        CandidateIndexFlag::Mmap => search::CandidateIndexMode::Mmap,
        CandidateIndexFlag::Auto => {
            // mmap is not compatible with the chimeric / refine / glyco in-RAM
            // passes (handled below); those are not the OOM-prone giant-mod-space
            // case, so `auto` simply keeps them on RAM.
            if params.chimeric || cli.refine || cli.glyco {
                search::CandidateIndexMode::Ram
            } else {
                match available_memory_bytes() {
                    Some(avail) => {
                        // Budget 60% of available memory for the candidate index;
                        // the rest covers spectra, the model, scoring scratch, etc.
                        let budget = (avail as f64 * 0.60) as u64;
                        if search::candidate_index::ram_candidate_index_fits(
                            &idx,
                            &params,
                            &cli.decoy_prefix,
                            budget,
                        ) {
                            search::CandidateIndexMode::Ram
                        } else {
                            eprintln!(
                                "[auto] in-RAM candidate index would exceed the ~{} GiB budget \
                                 (60% of {} GiB available) → using out-of-core mmap \
                                 (force with --candidate-index ram)",
                                budget >> 30,
                                avail >> 30
                            );
                            search::CandidateIndexMode::Mmap
                        }
                    }
                    // Can't read available memory: keep the byte-identical RAM
                    // default; the user can force --candidate-index mmap.
                    None => search::CandidateIndexMode::Ram,
                }
            }
        }
    };
    // The out-of-core (`Mmap`) candidate path materializes candidates lazily and
    // only syncs them into `prepared.candidates` AFTER the scan. The chimeric
    // Pass 2 and the refinement cascade both read `prepared.candidates` /
    // `bucket_index` DURING scanning, so they are not supported together with
    // `--candidate-index mmap` in this phase (fail loud rather than silently
    // produce wrong results).
    if params.candidate_index == search::CandidateIndexMode::Mmap {
        if params.chimeric {
            return Err(
                "--candidate-index mmap is not yet compatible with --chimeric \
                        (the chimeric Pass 2 needs the in-RAM candidate index)"
                    .into(),
            );
        }
        if cli.refine {
            return Err(
                "--candidate-index mmap is not yet compatible with --refine \
                        (the refinement cascade needs the in-RAM candidate index)"
                    .into(),
            );
        }
        if cli.glyco {
            return Err("--candidate-index mmap is not yet compatible with --glyco \
                        (glyco_search_run needs the in-RAM candidate index/bucket_index)"
                .into());
        }
    }
    // --glyco-hcd-pair only takes effect inside the glyco driver; without --glyco
    // it is silently inert, which would mislead a user who set it on purpose (code
    // review). Warn rather than error so it stays a no-op knob outside glyco mode.
    if cli.glyco_hcd_pair && !cli.glyco {
        eprintln!(
            "WARN: --glyco-hcd-pair has no effect without --glyco (it only drives \
             paired-scan candidate generation in glyco mode); ignoring it."
        );
    }
    // --glyco is a standalone driver (see the `if cli.glyco` early-return block
    // below): it writes its own `.glyco.pin` and skips the standard PIN/rescore/
    // TSV/Parquet/refine machinery entirely. Silently ignoring those flags would
    // mislead a user who expects them to apply, so fail fast instead.
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
        if cli.fdr.is_some() {
            unsupported.push("--fdr");
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
        // FDR-TRUST GUARD (adversarial review): reverse/shuffle decoys are generated from
        // the TARGET proteome and do NOT preserve N-X-S/T sequon DENSITY, but glyco
        // scoring gates BOTH targets and decoys on that sequon — so a generated
        // decoy search space is systematically different and the glyco FDR is
        // ANTI-CONSERVATIVE (q-values users should not trust). Trustworthy glyco FDR
        // needs an EXTERNAL target+decoy FASTA consumed with `--decoy-strategy none`.
        // `sequon-reverse` was added precisely to preserve N-X-S/T sequon density, so
        // it does NOT trip this warning; only the sequon-depleting strategies do.
        if !cli.decoy_strategy.eq_ignore_ascii_case("none")
            && !cli.decoy_strategy.eq_ignore_ascii_case("sequon-reverse")
            && !cli.decoy_strategy.eq_ignore_ascii_case("sequon")
        {
            eprintln!(
                "WARN: --glyco with --decoy-strategy {ds} GENERATES decoys ({ds} of the \
                 target proteome), which does NOT preserve N-X-S/T sequon density that \
                 glyco scoring gates on — the resulting .glyco.pin FDR is \
                 ANTI-CONSERVATIVE and its q-values should NOT be trusted. Use \
                 `--decoy-strategy sequon-reverse`, or supply an EXTERNAL target+decoy \
                 FASTA with `--decoy-strategy none --decoy-prefix <PREFIX>`.",
                ds = cli.decoy_strategy
            );
        }
        // Multi-file paired-scan guard (code review): --glyco-hcd-pair pairs each
        // ETD spectrum to a nearby HCD spectrum by position + precursor m/z over the
        // concatenated `all_spectra`. Across multiple --spectrum inputs the pairing
        // window would straddle file boundaries and silently mis-pair, so disable it
        // (rather than warn-and-continue) for multi-file runs. Validated per-file.
        if cli.glyco_hcd_pair && spectrum_paths.len() > 1 {
            eprintln!(
                "WARN: --glyco-hcd-pair is disabled for {} --spectrum inputs: paired-scan \
                 generation pairs ETD<->HCD over the concatenated spectra and would cross \
                 file boundaries. Run one file at a time to use paired generation.",
                spectrum_paths.len()
            );
        }
    }
    // --refine + --chimeric run together correctly but do NOT fully stack.
    //
    // Mechanism: chimeric secondary PSMs do not collapse the refinement's
    // confident-anchor set (the anchor set GROWS when both are on). What shrinks
    // is the pool of spectra left for refinement to rescue: chimeric already
    // explains many of them, so refinement has less work available, not worse
    // anchors. See the docs for the measured numbers.
    //
    // Warn so the user isn't surprised that the combination is much closer to
    // chimeric alone than to the sum.
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

    // ── Resolved-parameter banner (reanalysis auditability) ───────────────────
    // One consolidated record of every resolved search parameter, so a
    // (zero-config) run is fully reproducible/auditable from its log. Values that
    // were auto-detected from the data/store are tagged [detected].
    eprintln!("──────── andes resolved parameters ────────");
    eprintln!("  spectra        : {}", spectrum_path.display());
    eprintln!("  model          : (see 'Param model:' line above) [detected]");
    eprintln!(
        "  activation     : {:?} [detected]",
        param.data_type.activation
    );
    eprintln!(
        "  instrument     : {:?} [detected]",
        param.data_type.instrument
    );
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
    eprintln!(
        "  mods           : {}",
        if cli.mods.is_some() {
            "from --mods file"
        } else {
            "defaults (Cam-C fixed, Ox-M variable) + isobaric tag if detected"
        }
    );
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
        decoy_strategy, cli.decoy_prefix, params.chimeric
    );
    eprintln!("───────────────────────────────────────────");

    // ── 6+7. Stream-load + chunked search ─────────────────────────────────
    //
    // Spectra are parsed and scored in chunks of CHUNK_SIZE. Each chunk's
    // peak data lives in RAM only for the time it takes to score the chunk,
    // then is dropped before the next chunk is read. The Vec<Spectrum> that
    // survives into the PIN/TSV writers retains scan/title/precursor_mz/scan
    // (the only fields the writers read) but has empty peaks.
    //
    // This bounds peak-data memory to ~CHUNK_SIZE × per-spectrum peak size
    // regardless of dataset size — fixes the Astral-scale OOM where loading
    // all 123k spectra at once pushed RSS to 28 GB on a 31 GB VM.
    const CHUNK_SIZE: usize = 5000;

    let t_phase = std::time::Instant::now();

    // Configure the global Rayon worker pool BEFORE we build PreparedSearch
    // or run any chunks. `build_global()` panics if called twice; guard with
    // `OnceLock` so repeated CLI invocations within a single test process
    // don't blow up.
    static POOL_INIT: std::sync::OnceLock<()> = std::sync::OnceLock::new();
    POOL_INIT.get_or_init(|| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(cli.threads)
            .build_global()
            .expect("build_global");
    });
    eprintln!("Using {} worker threads", cli.threads);

    // Fragment tolerance of 0.5 Da is the canonical low-res HCD default.
    let fragment_tol_da = 0.5_f64;

    let bench_cap = if cli.max_spectra > 0 {
        cli.max_spectra
    } else {
        usize::MAX
    };
    let ms_level_u32 = cli.ms_level as u32;
    // Native `.raw`/`.d` readers search MS2 (identification) scans only. A non-2
    // `--ms-level` would otherwise make the Thermo iterator emit MS3 scans (which
    // carry a precursor and would be searched — e.g. TMT SPS-MS3 reporter scans)
    // or MS1 (no precursor → an empty run). Force MS2 + warn. `--ms-level` still
    // applies to mzML.
    if (is_raw || is_d) && cli.ms_level != 2 {
        eprintln!(
            "WARN: --ms-level={} is ignored for native .raw/.d input; these formats \
             search MS2 (identification) scans only.",
            cli.ms_level
        );
    }

    // The precursor-calibration pre-pass currently reads only mzML/MGF. For a
    // Thermo `.raw` or Bruker `.d` it would be misread as MGF, so skip
    // calibration and warn (native-format calibration support is a follow-up).
    if (is_raw || is_d) && params.precursor_cal_mode != PrecursorCalMode::Off {
        let fmt = if is_raw { "Thermo .raw" } else { "Bruker .d" };
        eprintln!(
            "WARN: --precursor-cal is not yet supported for {fmt} input; \
             proceeding without calibration."
        );
        params.precursor_cal_mode = PrecursorCalMode::Off;
    }

    // Calibration pre-pass. Candidate enumeration is precursor-tolerance
    // independent, so keep the cal pass's `PreparedParts` and reuse them for the
    // main pass instead of re-enumerating all 16.8M candidates (~15s saved on
    // Astral). `into_parts()` runs BEFORE tightening so the owned parts outlive
    // the `params` borrow the cal `PreparedSearch` held.
    let reuse_parts = if params.precursor_cal_mode != PrecursorCalMode::Off {
        let cal_prepared =
            PreparedSearch::prepare(&idx, &params, &scorer, fragment_tol_da, &cli.decoy_prefix);
        let cal_stats = run_precursor_calibration(
            &spectrum_path,
            is_mzml,
            ms_level_u32,
            bench_cap,
            &params,
            &cal_prepared,
        )?;
        let parts = cal_prepared.into_parts();
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
        Some(parts)
    } else {
        None
    };

    let intensity_model: Option<Arc<scoring_crate::IntensityModel>> = cli
        .intensity_model
        .as_ref()
        .map(|path| {
            eprintln!("loading intensity model from {} ...", path.display());
            scoring_crate::IntensityModel::load(path)
                .map(Arc::new)
                .map_err(|e| format!("intensity model {}: {e}", path.display()))
        })
        .transpose()?;

    let mut prepared = match (reuse_parts, params.candidate_index) {
        // Calibration reuse always takes the in-RAM parts (calibration is RAM-only).
        // Warn if the user explicitly requested mmap so they know it was not applied.
        (Some(parts), _) => {
            if cli.candidate_index == CandidateIndexFlag::Mmap {
                eprintln!(
                    "WARN: --candidate-index mmap is ignored when precursor calibration \
                     reuse is active; running in-RAM for this search."
                );
            }
            PreparedSearch::from_parts(&idx, &params, &scorer, fragment_tol_da, parts)
        }
        (None, search::CandidateIndexMode::Mmap) => {
            // Use a content-addressed cache path so repeated searches over the
            // same FASTA + params reuse the index without rebuilding.
            let path = index_cache_path(&idx, &params);
            PreparedSearch::prepare_mmap(
                &idx,
                &params,
                &scorer,
                fragment_tol_da,
                &cli.decoy_prefix,
                &path,
            )
            .map_err(|e| format!("build out-of-core candidate index: {e}"))?
        }
        (None, search::CandidateIndexMode::Ram) => {
            PreparedSearch::prepare(&idx, &params, &scorer, fragment_tol_da, &cli.decoy_prefix)
        }
    }
    .with_intensity_model(intensity_model);
    log_rss("after_prepared_search");
    match params.candidate_index {
        search::CandidateIndexMode::Ram => {
            eprintln!(
                "PreparedSearch: {} candidates, {} mass buckets (candidate-index: ram)",
                prepared.candidates.len(),
                prepared.bucket_index.len(),
            );
            warn_if_index_will_not_fit(prepared.candidates.len(), cli.glyco);
        }
        search::CandidateIndexMode::Mmap => eprintln!(
            "PreparedSearch: out-of-core candidate-index: mmap \
             (base peptides resolved lazily per spectrum)"
        ),
    }

    let bench_mode = cli.max_spectra > 0;

    let mut all_spectra: Vec<Spectrum> = Vec::new();
    let mut all_queues: Vec<TopNQueue> = Vec::new();

    let t_search_start = std::time::Instant::now();

    // Pipeline mzML/MGF parsing with Rayon scoring via a bounded sync_channel.
    // The parser runs on a dedicated thread and pushes
    // CHUNK_SIZE-sized `Vec<Spectrum>` payloads through the channel; the main
    // thread (this one) drains the channel and calls `prepared.run_chunk` on
    // each chunk (which is itself Rayon-parallel internally). With capacity 2
    // the parser stays at most one chunk ahead of the scorer, overlapping
    // parse-of-chunk-(N+1) with score-of-chunk-N — so parse time (Astral
    // ~2-3s per chunk) overlaps scoring instead of running serially.
    // MGF carries no MS level (always treated as MS2). (Native `.raw`/`.d` are
    // warned separately above.)
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
    let _ = mzml_warn_ms_level_emitted; // silenced — unused for now.

    // Under `--chimeric` (mzML or Thermo `.raw`), stream MS2 in CHUNK_SIZE batches,
    // each paired with a bounded per-chunk MS1 link (`read_with_ms1_chunked`). Pass 1 + Pass 2
    // run per chunk on this thread and peaks are dropped immediately, so RSS stays
    // bounded to ~CHUNK_SIZE spectra (NOT the whole file). The parser runs on a
    // dedicated thread so chunk N+1 parses while chunk N scores. The streaming
    // pipeline in the `else` branch handles the (default) non-chimeric path.
    let chimeric_input = chimeric_active;
    let mut parse_stats = ParseStats::default();

    for (file_idx, input_path) in spectrum_paths.iter().enumerate() {
        if bench_mode && all_spectra.len() >= bench_cap {
            break;
        }
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
        let title_prefix = title_prefix_for(spectrum_paths.len(), file_stem);

        if spectrum_paths.len() > 1 {
            eprintln!(
                "=== spectrum [{}/{}] {} ===",
                file_idx + 1,
                spectrum_paths.len(),
                input_path.display()
            );
        }

        if chimeric_input && !(file_is_mzml || file_is_raw) {
            return Err(format!(
                "--chimeric only supports mzML/.raw inputs, got {}",
                input_path.display()
            )
            .into());
        }

        let file_stats = if chimeric_input {
            let (tx, rx) = sync_channel::<(Vec<Spectrum>, Ms1Link)>(2);
            let spectrum_path = input_path.clone();
            let cap = remaining_cap;
            // The cascade is MS2-only by construction (MS2 paired with its preceding
            // MS1); hardcode MS2 so `--ms-level 3` can never widen the mzML reader's
            // range to admit MS3 (the .raw chunked reader is already MS2-only).
            let mslevel = 2;
            let parser_handle = thread::spawn(move || -> Result<(usize, Vec<String>), String> {
                if file_is_mzml {
                    let reader = MzMLReader::new(
                        input::open_buf_maybe_gz(&spectrum_path)
                            .map_err(|e| format!("open mzML: {e}"))?,
                    )
                    .with_ms_level_range(mslevel, mslevel)
                    .with_ms1_capture(true);
                    let (errc, errs) =
                        reader.read_with_ms1_chunked(CHUNK_SIZE, cap, |chunk, link| {
                            let _ = tx.send((chunk, link));
                        });
                    Ok((errc, errs))
                } else {
                    #[cfg(feature = "thermo")]
                    {
                        let reader = input::ThermoRawReader::open(&spectrum_path)
                            .map_err(|e| format!("open Thermo .raw: {e}"))?;
                        let (errc, errs) =
                            reader.read_with_ms1_chunked(CHUNK_SIZE, cap, |chunk, link| {
                                let _ = tx.send((chunk, link));
                            });
                        Ok((errc, errs))
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
            for (mut chunk_spectra, chunk_link) in rx {
                if let Some(prefix) = &title_prefix {
                    prefix_spectrum_titles(&mut chunk_spectra, prefix);
                }
                let offset = all_spectra.len();
                let mut queues = prepared.run_chunk(&chunk_spectra, offset);
                search::match_engine::run_pass2_coisolation(
                    &prepared,
                    &chunk_spectra,
                    &mut queues,
                    &params,
                    &chunk_link,
                    offset,
                );
                file_offset += chunk_spectra.len();
                ms1_linked += chunk_link.ms1_peaks.len();
                all_queues.extend(queues);
                for mut spec in chunk_spectra.into_iter() {
                    // Peaks are normally dropped post-scoring to bound memory
                    // (only the metadata is needed downstream). Under `--refine`
                    // or `--glyco` we RETAIN them: refine's Pass-2 re-scores
                    // unidentified spectra; glyco_search_run needs the full
                    // peak lists for oxonium ion detection. Memory cost: full
                    // peak buffer stays resident (acceptable; both are opt-in).
                    if !cli.refine && !cli.glyco {
                        spec.peaks = Vec::new();
                    }
                    all_spectra.push(spec);
                }
                report_search_progress(all_spectra.len(), t_search_start);
            }
            let (err_count, first_errors) = parser_handle
                .join()
                .map_err(|_| "chimeric parser thread panicked".to_string())??;
            eprintln!(
                "chimeric mode: streamed {} MS2 spectra ({} MS1 scans linked) from {}",
                file_offset,
                ms1_linked,
                input_path.display()
            );
            log_rss("after_chimeric_stream_search");
            ParseStats {
                error_count: err_count,
                first_errors,
            }
        } else {
            let (tx, rx) = sync_channel::<Vec<Spectrum>>(2);
            let spectrum_path = input_path.clone();
            let parser_handle = thread::spawn(
                move || -> Result<ParseStats, Box<dyn std::error::Error + Send + Sync>> {
                    if file_is_mzml {
                        let reader = MzMLReader::new(
                            input::open_buf_maybe_gz(&spectrum_path)
                                .map_err(|e| format!("open mzML: {e}"))?,
                        )
                        .with_ms_level_range(ms_level_u32, ms_level_u32);
                        Ok(send_chunks(reader, CHUNK_SIZE, remaining_cap, tx))
                    } else if file_is_raw {
                        #[cfg(feature = "thermo")]
                        {
                            let reader = input::ThermoRawReader::open(&spectrum_path)
                                .map_err(|e| format!("open Thermo .raw: {e}"))?
                                .with_ms_level(Some(2));
                            Ok(send_chunks(reader, CHUNK_SIZE, remaining_cap, tx))
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
                            let reader = input::TimsTofReader::open(&spectrum_path)
                                .map_err(|e| format!("open Bruker .d: {e}"))?;
                            Ok(send_chunks(reader, CHUNK_SIZE, remaining_cap, tx))
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
                        Ok(send_chunks(reader, CHUNK_SIZE, remaining_cap, tx))
                    }
                },
            );

            log_rss("after_parser_thread_spawn");

            for mut chunk in rx {
                if chunk.is_empty() {
                    continue;
                }
                if let Some(prefix) = &title_prefix {
                    prefix_spectrum_titles(&mut chunk, prefix);
                }
                // SPEED (--glyco): the per-spectrum PEPTIDE search is pure waste in
                // glyco mode — glyco_search_run re-derives its own candidates from
                // `prepared` + the retained peaks and writes its own `.glyco.pin`
                // (the standard PIN is skipped), so these PSM queues are discarded.
                // Skip the scoring (parsing + peak retention below still run); ~16%
                // faster on glyco. Non-glyco path is unchanged.
                if !cli.glyco {
                    let offset = all_spectra.len();
                    let queues = prepared.run_chunk(&chunk, offset);
                    all_queues.extend(queues);
                }
                for mut spec in chunk.into_iter() {
                    // See the chimeric-loop note above: normally peaks are
                    // dropped post-scoring to bound memory, but `--refine` or
                    // `--glyco` needs the full peak lists retained.
                    if !cli.refine && !cli.glyco {
                        spec.peaks = Vec::new();
                    }
                    all_spectra.push(spec);
                }
                report_search_progress(all_spectra.len(), t_search_start);
                log_rss(&format!("after_chunk_{:06}_specs", all_spectra.len()));
            }

            match parser_handle.join() {
                Ok(Ok(stats)) => stats,
                Ok(Err(e)) => return Err(format!("parser thread error: {e}").into()),
                Err(_) => return Err("parser thread panicked".into()),
            }
        };
        merge_parse_stats(&mut parse_stats, file_stats);
    }

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

    log_rss("after_all_spectra");
    let search_elapsed = t_search_start.elapsed();
    eprintln!(
        "Loaded+scored {} spectra from {} in chunks of {} [PHASE stream_search: {:.2}s]",
        all_spectra.len(),
        spectrum_path.display(),
        CHUNK_SIZE,
        t_phase.elapsed().as_secs_f64()
    );
    if bench_mode {
        eprintln!("Bench mode: capped at {} spectra", cli.max_spectra);
    }

    // `Mmap` mode: drain the per-spectrum materialized candidates accumulated
    // during the scan into `prepared.candidates` so the PIN/TSV writers resolve
    // every PSM's `candidate_idxs` against a real candidate slice (no-op + cheap
    // in the default `Ram` mode, where candidates already live there).
    prepared.sync_materialized_candidates();

    // Downstream code uses these names.
    let spectra = all_spectra;
    let mut queues = all_queues;

    if cli.glyco {
        eprintln!(
            "Standard peptide search skipped in --glyco mode (PSMs unused; glyco driver \
             derives its own candidates); {} spectra loaded ({:.2}s)",
            spectra.len(),
            search_elapsed.as_secs_f64()
        );
    } else {
        let non_empty = queues.iter().filter(|q| !q.is_empty()).count();
        eprintln!(
            "Search complete: {non_empty} / {} spectra have PSMs (match_spectra wall: {:.2}s)",
            spectra.len(),
            search_elapsed.as_secs_f64()
        );
    }

    // ── 7a. Glyco mode: run glyco scoring and write .glyco.pin, then return ──
    // When --glyco is active, we run the glyco-PSM scoring driver over ALL
    // accumulated spectra (using the PreparedSearch from the standard search)
    // and write a separate `.glyco.pin` file.  The standard PIN is skipped.
    if cli.glyco {
        run_glyco(
            &cli,
            &spectra,
            &prepared,
            &params,
            &idx,
            &output_pin_path,
            spectrum_paths,
            &target_db,
            detected_activation_instrument,
            t_total,
        )?;
        return Ok(());
    }

    // ── 7b. PTM-refinement cascade (Pass-2) ───────────────────────────────────
    // Opt-in (`--refine`). Runs a scoped Pass-2 over the unidentified spectra
    // against the confidently-identified proteins, with the refinement
    // variable-mod tier applied. Computed BEFORE the Pass-1 PIN write (it reads
    // `&queues`/`&spectra` immutably). The Pass-2 winners ARE merged into
    // `queues` (force_pushed per scan by `merge_into_pass1`) and written to a
    // SINGLE unified PIN, so a scan's unmodified Pass-1 and modified Pass-2 PSMs
    // compete in one report. Collapse best-per-scan downstream; optionally split
    // by mokapot --group-column on IsRefinement/RefinementModClass.
    // `--refine-debug-split-pin` instead emits Pass-2 to a separate PIN.
    let refine_output = if cli.refine {
        // Refinement config: explicit YAML or the built-in 5-mod default tier.
        let base_cfg = match &cli.refine_config {
            Some(p) => search::RefineConfig::from_yaml_str(&std::fs::read_to_string(p)?)
                .map_err(|e| format!("parsing --refine-config {}: {e}", p.display()))?,
            None => search::RefineConfig::default_tier(),
        };
        // Max variable mods for refinement comes from the refine config/tier
        // (set it in the `--refine-config` YAML; the former CLI override was removed).
        let cfg = search::RefineConfig {
            max_mods: base_cfg.max_mods,
            ..base_cfg
        };

        // High-res signal: the resolved model's instrument class. High-res
        // instruments fragment-match in ppm (20 ppm vs 0.5 Da ion-trap), which is
        // exactly the regime where the near-isobaric refinement deltas (e.g.
        // deamidation +0.984 vs a C13 isotope error) are resolvable.
        let high_res = param.data_type.instrument.is_high_resolution();

        // `high_res_only` comes from the refine config/tier (set it in the
        // `--refine-config` YAML; the former CLI overrides were removed).
        let high_res_only = cfg.high_res_only;
        if high_res_only && !high_res {
            eprintln!(
                "WARN: refine is high-res-only and the data is low-res; skipping refinement."
            );
            None
        } else {
            // Target-only db recovered from the Pass-1 combined index; the
            // cascade regenerates its own decoys for the scoped Pass-2.
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
                0.01, // report_q = 1% FDR report threshold (spec)
                high_res,
                fragment_tol_da,
                &cli.decoy_prefix,
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

    // ── 8. Write PIN — single unified list (Pass-1 ⊕ merged refine) ──────────
    // Bench mode still writes PIN (so we can diff against the reference
    // fixture) but skips TSV.
    let t_phase = std::time::Instant::now();
    // `pin_candidates`/`pin_index` own the list written to PIN/TSV: the merged
    // Pass-1 ⊕ Pass-2 candidates when --refine is active, else the Pass-1 pool
    // moved straight through. Owned (not borrowed) so they outlive the writes.
    // Owned PIN candidate list + index. In the refine arm we MOVE
    // `prepared.candidates` into the merge (which extends it in place); in the
    // non-refine arm we move it straight into the tuple. Either way the Pass-1
    // pool is never duplicated. `prepared.candidates`/`idx` are not used past this
    // point (verified), so the move is safe.
    let refine_merged = refine_output.is_some();
    let pin_candidates;
    let pin_index;
    if let Some(out) = refine_output {
        let merged =
            search::refinement::merge_into_pass1(&mut queues, prepared.candidates, &idx, out);
        pin_candidates = merged.candidates;
        pin_index = merged.index;
    } else {
        pin_candidates = prepared.candidates;
        pin_index = idx;
    }

    // Engine-wide retention-time features (additive): after all PSMs are scored
    // and the final PIN candidate pool is fixed, run per-run RT self-calibration
    // and populate DeltaRT/AbsDeltaRT/DeltaRTNorm (+ predicted_rt for the QPX)
    // onto each PSM. Neutral 0.0 (baseline-identical) when observed RT is
    // unavailable or calibration cannot be fit. See
    // `docs/plans/glyco/50-roadmap/rt-prediction-design.md` (Commit 1).
    output::populate_rt_features(&spectra, &mut queues, &pin_candidates);

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
    log_rss("after_pin_write");
    if cli.refine {
        if refine_merged {
            eprintln!("Refinement PSMs merged into unified PIN (Pass-1 ⊕ Pass-2).");
        } else {
            eprintln!("Refinement produced no Pass-2 PSMs; unified PIN contains Pass-1 only.");
        }
    }

    // ── Rescore: run Percolator (or the native rescorer) on the PIN, join ─────
    // PEP/q-value downstream. Sits BETWEEN the PIN write and the QPX/TSV writes
    // so the parsed `SpecId → PercolatorPsm` map can be threaded into
    // `write_qpx` and the filtered TSV. `None` keeps every downstream output
    // identical to a non-rescore run.
    //
    // Backend choice: rescoring runs ONLY when explicitly requested — `--rescore`
    // forces Percolator, `--rescore-native` forces the native rescorer. `--fdr` /
    // `--pep` are thresholds applied by such a run; they never start one on their
    // own (that would silently launch a backend / the non-production native path).
    let fdr_threshold = cli.fdr.unwrap_or(0.01);
    let (use_percolator, use_native) = if cli.rescore {
        (true, false)
    } else if cli.rescore_native {
        (false, true)
    } else {
        (false, false)
    };
    let rescore_map: Option<std::collections::HashMap<String, output::PercolatorPsm>> =
        if use_percolator {
            let backend = output::resolve_backend(
                cli.percolator_bin.as_deref(),
                cli.percolator_docker,
                &cli.percolator_image,
            )
            .map_err(|e| format!("{e}"))?;
            eprintln!("Rescore: Percolator backend = {}", backend.describe());
            let extra: Vec<String> = cli
                .percolator_args
                .split_whitespace()
                .map(|s| s.to_string())
                .collect();
            let t_resc = std::time::Instant::now();
            let map = output::run_percolator(&backend, &output_pin_path, &extra)
                .map_err(|e| format!("{e}"))?;
            eprintln!(
                "Rescore: Percolator produced {} target PSMs [{:.2}s]",
                map.len(),
                t_resc.elapsed().as_secs_f64()
            );
            Some(map)
        } else if use_native {
            eprintln!(
                "Rescore: NATIVE GBDT rescorer (no Percolator; leakage-safe 3-fold target-decoy CV). \
                 Fallback — use --rescore (Percolator) for production-grade FDR."
            );
            let t_resc = std::time::Instant::now();
            let pin_text = std::fs::read_to_string(&output_pin_path)
                .map_err(|e| format!("reading PIN for native rescore: {e}"))?;
            let map = rescore::native_rescore_pin(&pin_text, 42)
                .map_err(|e| format!("native rescore: {e}"))?;
            eprintln!(
                "Rescore: native produced {} PSMs (q-value+PEP) [{:.2}s]",
                map.len(),
                t_resc.elapsed().as_secs_f64()
            );
            Some(map)
        } else {
            None
        };

    // Filtered TSV: target PSMs whose Percolator q-value ≤ --fdr. Written next to
    // the PIN as `<stem>.q<fdr>.tsv` (rescore only). One row per accepted PSM with
    // the join key + q/PEP + peptide/proteins from the Percolator result.
    if let Some(ref map) = rescore_map {
        let fdr = fdr_threshold;
        let stem = report_base
            .file_stem()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| "andes".to_string());
        let dir = report_base.parent().unwrap_or_else(|| Path::new("."));
        let q_tag = format!("{fdr}").replace('.', "p");
        let filtered_path = dir.join(format!("{stem}.q{q_tag}.tsv"));
        // q-value is the primary set-level FDR control; an optional --pep ANDs a
        // per-PSM local-FDR cap on top.
        let pep_cap = cli.pep;
        let mut accepted: Vec<&output::PercolatorPsm> = map
            .values()
            .filter(|p| p.q_value <= fdr && pep_cap.is_none_or(|t| p.pep <= t))
            .collect();
        accepted.sort_by(|a, b| {
            a.q_value
                .partial_cmp(&b.q_value)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let pep_note = pep_cap
            .map(|t| format!(" and PEP<={t}"))
            .unwrap_or_default();
        match write_filtered_tsv(&filtered_path, &accepted) {
            Ok(()) => eprintln!(
                "Wrote filtered TSV: {} ({} PSMs at q<={}{})",
                filtered_path.display(),
                accepted.len(),
                fdr,
                pep_note
            ),
            Err(e) => eprintln!("WARN: could not write {}: {e}", filtered_path.display()),
        }
    }

    // ── Run summary: final tolerances + per-modification PSM tally ────────────
    // andes auto-resolves the scoring model and the precursor/fragment
    // tolerances from the data, so the FINAL search parameters differ from
    // whatever the user passed on the CLI. Emit a summary to stderr and write a
    // `statistics.log` next to the PIN so a run's true parameters (and which
    // PTMs were identified, with PSM counts) are always recoverable.
    let run_stats = output::RunStatistics::compute(&queues, &pin_candidates, &params, &param.mme);
    eprint!("{}", run_stats.render());
    let stats_path = report_base
        .parent()
        .map(|d| d.join("statistics.log"))
        .unwrap_or_else(|| std::path::PathBuf::from("statistics.log"));
    match output::write_statistics_log(&run_stats, &stats_path) {
        Ok(()) => eprintln!("Wrote statistics: {}", stats_path.display()),
        Err(e) => eprintln!("WARN: could not write {}: {e}", stats_path.display()),
    }

    // ── QPX `.idparquet/` bundle (optional, OpenMS-compatible) ───────────────
    // Mirrors the PIN write from the SAME unified candidate/index pair. Sourced
    // here (not in the bench/TSV arms) so it is produced for every run that asks
    // for it, including bench mode. Run identifier = the spectrum file stem so
    // it is stable across re-runs of the same input.
    if let Some(ref parquet_dir) = cli.output_parquet {
        let run_id = spectrum_path
            .file_stem()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| "andes_run".to_string());
        let primary_paths: Vec<String> = spectrum_paths
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
            &param.mme,
            &run_id,
            &primary_paths,
            rescore_map.as_ref(),
        ) {
            Ok(()) => eprintln!("Wrote QPX bundle: {}", parquet_dir.display()),
            // Finding 3.9: failure to write an explicitly-requested output is a
            // hard error, not a WARN-and-exit-0. The user asked for this
            // artifact; silently exiting 0 without it hides data loss from any
            // calling pipeline.
            Err(e) => {
                return Err(format!(
                    "failed to write requested --output-parquet {}: {e}",
                    parquet_dir.display()
                )
                .into());
            }
        }
    }

    if bench_mode {
        eprintln!("Bench mode: skipping TSV write.");
        return Ok(());
    }

    // ── 9. Write TSV (optional) ───────────────────────────────────────────────
    if let Some(ref tsv_path) = cli.output_tsv {
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
        // Use the SAME merged candidate/index pair the PIN write used: after
        // `merge_into_pass1`, `queues` holds candidate indices into the merged
        // candidate list (Pass-1 ⊕ offset Pass-2), so resolving them against the
        // un-merged `prepared.candidates`/`idx` would be wrong/out-of-bounds.
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

    Ok(())
}

/// Write the FDR-filtered rescore TSV: one row per accepted target PSM (already
/// filtered to q ≤ fdr and sorted by q ascending). Columns: the Percolator join
/// key (`PSMId`, == PIN SpecId), `q-value`, `posterior_error_probability`,
/// `peptide`, `proteinIds`.
pub(crate) fn write_filtered_tsv(
    path: &Path,
    psms: &[&output::PercolatorPsm],
) -> std::io::Result<()> {
    use std::io::Write;
    let mut w = std::io::BufWriter::new(File::create(path)?);
    writeln!(
        w,
        "PSMId\tq-value\tposterior_error_probability\tpeptide\tproteinIds"
    )?;
    for p in psms {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}",
            p.psm_id, p.q_value, p.pep, p.peptide, p.proteins
        )?;
    }
    w.flush()
}

// ── Training pipeline ─────────────────────────────────────────────────────────
