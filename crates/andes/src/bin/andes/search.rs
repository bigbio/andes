//! The standard database-search driver.

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
    cli_fragment_tol_override, default_aa_set_with_tag, load_param_from_store, parse_enzymes,
    resolve_metadataless_selection, warn_if_universal_protease_combo,
};
use crate::mono::{correct_chunk, MonoDump, MonoStats};
use crate::quant_run::{QuantCollector, QuantSettings, RunScans};
use crate::rescore;
use crate::spectra::{
    add_mzml_skipped, detect_dominant_activation, detect_instrument_type_for_path,
    detect_isobaric_sampled, input_format_flags, merge_parse_stats, ms3_carries_reporters,
    prefix_spectrum_titles, run_precursor_calibration, sample_ms2_meta, sample_spectra,
    send_chunks, title_prefix_for, tolerance_ppm_display, ParseStats,
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

/// What [`retrieval_choice`] decided, and what to tell the user about it.
pub(crate) struct RetrievalChoice {
    /// Score each spectrum against a fragment-ion-index shortlist rather than
    /// every peptidoform in its precursor windows.
    pub use_index: bool,
    /// One human-readable clause naming the reason, or `None` when there is
    /// nothing to report (`--glyco`, which looks candidates up by backbone mass).
    pub report: Option<String>,
    /// Set when the user forced the index somewhere it cannot be used.
    pub refused_because: Option<&'static str>,
}

/// Pick the candidate-retrieval strategy: the fragment-ion index, because it
/// was MEASURED to win at both resolutions:
///
///  - **High-resolution fragments**: +23% PSMs on Astral at equal entrapment FDP
///    over per-spectrum enumeration.
///  - **Low-resolution fragments**: only the most intense peaks vote
///    (`fragment_index::vote_peaks_for`); with every peak voting, the wide bins
///    let noise swamp the count (TMT 12,281 -> 4,368 PSMs @1%). Capped, the
///    index matches enumeration on TMT and UPS1 at ~1/19 and ~1/11 of its CPU
///    time (TMT 289 s vs 5,404 s).
///
/// `excluded_because` names a search the index cannot serve, which keeps
/// enumeration: `--glyco` looks candidates up by backbone mass (precursor minus
/// glycan) in its own sequon-only candidate table, and ETD/UVPD spectra fragment
/// into ion series the b/y-only index does not hold.
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

/// One MS1 scan's (m/z, intensity) peaks, m/z-sorted.
type Ms1Peaks = Vec<(f64, f32)>;

/// Mass span (Da) of one fragment-index chunk.
///
/// A chunk's index holds every peptidoform in its precursor mass window, so
/// chunks are bounded by MASS SPAN, not only by spectrum count: 1,000 high-mass
/// phospho spectra spanned 450 Da and their window held 39.7M forms / 3.3G ion
/// entries. With the allocation-free index build a 150 Da slice at ~120k phospho
/// forms per Da is ~18M forms / ~5 GB, and each record is expanded once per slice
/// instead of once per 15 Da. The slice width follows the memory budget:
/// measured on the phospho space, a 150 Da slice peaks at ~25 GB RSS, about
/// 0.13 GB per Da over a ~5 GB base, so a 31 GB machine gets ~60 Da and a 96 GB
/// allowance the full 150 Da. `flag` is `--fragment-index-slice-da`.
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

/// Whether a spectrum path is gzip-compressed (`.gz`).
fn is_gzip_path(path: &Path) -> bool {
    path.extension()
        .and_then(|e| e.to_str())
        .is_some_and(|e| e.eq_ignore_ascii_case("gz"))
}

/// Sort one file's spectra by precursor mass and score them in fragment-index
/// chunks, appending to `all_spectra` / `all_queues` with peaks dropped.
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
            // `--refine` re-scores unidentified spectra in its Pass 2, so their
            // peaks stay resident.
            if !keep_peaks {
                spec.peaks = Vec::new();
            }
            all_spectra.push(spec);
        }
        report_search_progress(all_spectra.len(), t_search_start);
    }
}

/// Measures each input file's fragment mass error once, before its spectra
/// are scored, and installs it on the scorer (see
/// [`search::fragment_error`]). Active on high-resolution acquisitions without
/// a `--fragment-tol-*` override.
struct FragmentErrorMeasurer {
    enabled: bool,
    done: bool,
}

impl FragmentErrorMeasurer {
    /// Forget the previous file's measurement.
    fn start_file(&mut self, prepared: &PreparedSearch) {
        self.done = false;
        if self.enabled {
            prepared.scorer.set_fragment_error(None);
        }
    }

    /// Measure on `spectra` (a whole file, or its first chunk) unless this
    /// file was measured already. Scores up to
    /// [`search::fragment_error::SAMPLE_SPECTRA`] spectra around the median
    /// precursor mass, in fragment-index slices when the index is on.
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
                prepared.scorer.set_fragment_error(Some(m));
                eprintln!(
                    "Fragment error: offset {:+.2} ppm, spread {:.2} ppm from {} ions of {} PSMs \
                     at 1% q in {} sampled spectra -> PIN feature tolerance {:.1} ppm, \
                     ppm-Gaussian centre {:+.2} width {:.2} ppm ({:.1}s)",
                    m.offset_ppm,
                    m.sigma_ppm,
                    m.ions,
                    m.psms,
                    sample.len(),
                    m.feature_tolerance_ppm(),
                    m.offset_ppm,
                    m.gaussian_sigma_ppm(),
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

/// The windowed peak filter setting of the search. The isobaric filter (100 Da
/// windows, top 20) trims the dense tails of ion-trap CID-MS2; on high-res
/// HCD-MS2 TMT it measured neutral (TKO hrMS2: 13,390 PSMs off vs 13,387 on),
/// and on high-res label-free Astral data the same filter lost ~14% of PSMs,
/// so a high-res acquisition turns it off. An explicit `--peak-filter` wins.
fn isobaric_peak_filter(explicit: Option<(f64, usize)>, high_res: bool) -> Option<(f64, usize)> {
    explicit.or(high_res.then_some((0.0, 0)))
}

/// Resolve the PIN path and the durable report base (the filtered q-value
/// TSV and `<stem>.statistics.log` are named from it).
///
/// An explicit `--output-pin` is both. Otherwise the base is `--output-tsv`,
/// else `--output-parquet`, else `./andes.pin`; with `keep_pin` the PIN is
/// `<base>.pin`, without it the PIN goes into a temp dir stored in `tmp`.
/// The output base: `--output-pin`, else `--output-tsv`, else
/// `--output-parquet`, else `./andes.pin`. Report files and the default scratch
/// and cache directories sit next to it.
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

/// A temporary directory inside `scratch_dir` (created if missing), removed
/// when dropped; an error names `--temp-dir`.
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

pub(crate) fn run(cli: Cli) -> Result<(), Box<dyn std::error::Error>> {
    // These three were validated as Some(..) by main() before calling run().
    if cli.spectrum.is_empty() {
        return Err("no --spectrum inputs".into());
    }
    // Configure the global Rayon worker pool BEFORE anything runs on it: the
    // acquisition sample of a `.d` opens the timsTOF reader, which runs on the
    // global pool and would otherwise create it at the default size.
    // `build_global()` panics if called twice; guard with `OnceLock` so repeated
    // CLI invocations within a single test process don't blow up.
    static POOL_INIT: std::sync::OnceLock<()> = std::sync::OnceLock::new();
    POOL_INIT.get_or_init(|| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(cli.threads)
            .build_global()
            .expect("build_global");
    });
    let mut cli = cli;
    if cli.top_n == 0 {
        return Err("--top-n must be at least 1".into());
    }
    // Run-private temporary files go to `--temp-dir`, by default next to the
    // output; the quantification settings read the resolved directory.
    let scratch_dir = crate::cli::dir_or_output_dir(
        cli.temp_dir.as_deref(),
        &output_base(
            cli.output_pin.as_deref(),
            cli.output_tsv.as_deref(),
            cli.output_parquet.as_deref(),
        ),
    );
    cli.temp_dir = Some(scratch_dir.clone());
    // `--mbr`: the per-run artifacts go to `--quant-dir`, or to a temporary
    // directory that lives until the end of the run.
    let mut _mbr_tmp: Option<tempfile::TempDir> = None;
    if cli.mbr {
        check_mbr_args(&cli)?;
        if cli.quant_dir.is_none() {
            let dir = scratch_tempdir("andes-quant-", &scratch_dir)?;
            cli.quant_dir = Some(dir.path().to_path_buf());
            _mbr_tmp = Some(dir);
        }
    }
    // The per-run artifacts are named by the input's stem.
    if cli.quant_dir.is_some() {
        crate::quant_run::check_unique_stems(&cli.spectrum)?;
    }
    QuantSettings::check_args(&cli)?;
    if let Some(warning) = deuterated_tag_warning(cli.tmt, cli.precursor_tol) {
        eprintln!("{warning}");
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
    // PIN destination. Normally `--output-pin`. Under `--rescore` /
    // `--rescore-native` without an explicit path, the PIN (and Percolator's
    // result files, which sit next to it) goes next to the report base —
    // `--output-tsv`, else `--output-parquet`, else `./andes` — as
    // `<base>.pin` when `--keep-pin` is true (the default), or to a temp dir
    // (kept alive in `_rescore_tmp`, removed at function exit) when it is false.
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
    // ── 1. Load the target database (FASTA or UniProt XML) ───────────────────
    let uniprot_opts = input::UniProtOptions {
        variants: cli.uniprot_variants,
        isoforms: cli.uniprot_isoforms,
        chains: cli.uniprot_chains,
        evidence: cli.uniprot_evidence,
        knockoffs: Some(input::KnockoffOptions {
            per_variant: cli.knockoffs,
            enzyme: search_enzyme,
        }),
    };
    let loaded = input::load_protein_database(&database_path, &uniprot_opts)?;
    let mut target_db = loaded.db;
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
    if !cli.annotations.is_empty() {
        let stats =
            input::load_annotations(&cli.annotations, &mut target_db, uniprot_opts.knockoffs)?;
        log_claims("Annotations", &stats);
    }
    // Knockoff sites for every annotated site, before decoys are built so they
    // follow to the decoys like real sites.
    if let Some(ko) = uniprot_opts.knockoffs {
        let stats = input::add_site_knockoffs(&mut target_db, ko);
        if stats.knockoffs > 0 {
            eprintln!(
                "Site claims: {} annotated sites, {} knockoff sites ({} sites with no position \
                 for one)",
                stats.sites, stats.knockoffs, stats.untested
            );
        }
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
                (cli.uniprot_variants, "--uniprot-variants"),
                (cli.uniprot_isoforms, "--uniprot-isoforms"),
                (cli.uniprot_chains, "--uniprot-chains"),
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

    // ── 2. Build SearchIndex (targets + strategy-generated decoys) ────────────
    let decoy_strategy =
        search::decoy::DecoyStrategy::from_name(&cli.decoy_strategy).ok_or_else(|| {
            format!(
                "unknown --decoy-strategy '{}' (expected shuffle/reverse/sequon/none)",
                cli.decoy_strategy
            )
        })?;
    let proteases = search::decoy::Proteases {
        primary: search_enzyme,
        extra: extra_enzymes.clone(),
        num_tolerable_termini: match cli.enzyme_specificity {
            EnzymeSpecificity::Fully => 2,
            EnzymeSpecificity::Semi => 1,
            EnzymeSpecificity::NonSpecific => 0,
        },
    };
    // Isoform and chain entries share most peptides with their canonical
    // entry. A non-specific search shuffles each entry's whole sequence on its
    // own, so the decoys share none, and the provenance columns (and the
    // decoy search space) would differ between targets and decoys.
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
    let t_phase = std::time::Instant::now();
    let idx = SearchIndex::from_target_db_for(
        &target_db,
        &cli.decoy_prefix,
        cli.decoy_suffix.as_deref(),
        decoy_strategy,
        cli.decoy_seed,
        &proteases,
    );
    eprintln!(
        "[PHASE search_index_build: {:.2}s]",
        t_phase.elapsed().as_secs_f64()
    );

    // ── 3. Build AminoAcidSet ────────────────────────────────────────────────
    //
    // If --mods is given, parse the mods.txt file. Otherwise use the
    // defaults: CAM fixed on C, Oxidation variable on M and Acetyl
    // variable on the protein N-term.
    //
    // `num_mods_from_file` is populated only when --mods is given and the
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
        // No --mods: andes defaults (CAM-C fixed, Ox-M variable, protein-N-term
        // Acetyl variable). The isobaric
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

    // MUST precede every spectrum read. The EThcD labelling policy is consulted
    // by the mzML reader itself, and the activation pre-pass below reads 64
    // spectra: installed any later, that pre-pass sees the default (`false`),
    // relabels EThcD as HCD, and the dominant activation that drives model
    // selection is wrong no matter what `--ethcd-activation` said. The one-shot
    // warning would also fire there under the wrong policy and suppress the
    // correct message for the rest of the run.
    input::mzml::init_ethcd_as_etd(cli.ethcd_activation == EthcdActivationFlag::Etd);

    // Detect (activation, instrument) from the input for auto-routing.
    // mzML peeks the file; Thermo `.raw` reads vendor metadata; Bruker `.d`
    // is always CID/TimsTOF (DDA-PASEF). Detection runs for every metadata-
    // bearing format; an explicit `--fragmentation` still overrides the
    // detected activation (see the resolution below).
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

    // The acquisition context, carried separately from the model chosen below.
    // Without analyzer metadata and without a `--fragment-tol-*` unit, the
    // fragment resolution is read from the peaks of the first spectra.
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
        detected_activation_instrument,
        cli.fragment_tol_ppm,
        cli.fragment_tol_da,
        peak_evidence,
        sample_ms2_meta(&spectrum_path, is_mzml, is_raw, is_d),
    );

    let t_phase = std::time::Instant::now();
    // The block also returns the (activation, instrument) it routed on and the
    // chosen model id, so the resolved-parameter banner can report what was
    // detected separately from what the substituted model was trained on.
    let (mut param, routed_activation, routed_instrument, resolved_model_id) = {
        // ── Resolve (activation, instrument) for the Parquet model store. ─────
        //
        // Metadata-first precedence: a fully detected (activation, instrument)
        // wins outright. When only the activation method is detected (analyzer
        // unknown), or nothing is detected (MGF / metadata-less mzML/.raw), the
        // metadata-less resolver folds in `--fragmentation` (which wins over the
        // detected activation) and the `--fragment-tol-*` flags. Default when no metadata: CID / low-res
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
                    acquisition.inferred_high_res(),
                ),
                None => resolve_metadataless_selection(
                    None,
                    cli.fragmentation,
                    cli.fragment_tol_ppm,
                    cli.fragment_tol_da,
                    acquisition.inferred_high_res(),
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
        (p, activation, instrument_opt, model_id)
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
    // The label the run carries: `--tmt <plex>` names it; otherwise it is read
    // from the reporter ions (also under an explicit `--protocol`, so the
    // injected tag below gets the kit's mass: TMTpro and iTRAQ 8-plex are
    // 304.2 Da, not the 229.2 / 144.1 Da of TMT 6–11 and iTRAQ 4-plex).
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
    match cli.protocol {
        Protocol::Tmt => param.data_type.protocol = model::protocol::Protocol::TMT,
        Protocol::Itraq => param.data_type.protocol = model::protocol::Protocol::ITRAQ,
        Protocol::ItraqPhospho => {
            param.data_type.protocol = model::protocol::Protocol::ITRAQPhospho
        }
        // `--protocol auto` (the default) engages the isobaric windowed peak
        // filter when the run carries a label. Label-free runs are unchanged.
        Protocol::Auto => match label {
            Some(l @ (input::IsobaricLabel::Tmt | input::IsobaricLabel::TmtPro)) => {
                eprintln!(
                    "Protocol resolver: {} {l:?} → isobaric protocol (windowed peak filter on low-res MS2)",
                    if cli.tmt.is_some() {
                        "--tmt"
                    } else {
                        "auto-detected"
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
                        "auto-detected"
                    }
                );
                param.data_type.protocol = model::protocol::Protocol::ITRAQ;
            }
            None => {}
        },
        _ => {}
    }
    // C1: parameter-free path only (no explicit --mods). When the protocol
    // resolves to TMT/iTRAQ, or `--tmt` names the kit, inject the tag as a fixed mod on K + peptide
    // N-term so labeled peptides match their precursor mass — otherwise the
    // reporter filter engages but every labeled candidate is +tag Da off and
    // misses. The kit comes from `label` when known, else the protocol's
    // smaller kit. With explicit --mods the user owns the mod set (they may
    // already supply the tag), so those runs stay byte-identical.
    if cli.mods.is_none() {
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
            // `--tmt PLEX` names the kit whatever the protocol (e.g. phospho).
            (_, Some(l)) if cli.tmt.is_some() => Some(l.tag()),
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
    // Feature tolerances and the fragment-index vote mode follow the
    // acquisition, not the model's instrument.
    scorer.set_acquisition_high_res(Some(acquisition.high_res));
    // An explicit `--fragment-tol-*` sets the feature tolerance on every input
    // format, detected analyzer or not, and replaces the per-run measurement.
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
    warn_restricted_artefacts(&aa, cli.uniprot_variants);

    // ── 5. Build SearchParams ─────────────────────────────────────────────────
    // The `sites=annotated` mods, whose knockoff-site placements exist only for
    // the site claims and are not reported.
    let restricted_site_masses: Vec<f64> = aa
        .distinct_mods()
        .iter()
        .filter_map(|m| m.annotated_site_mass)
        .collect();
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
    // Dropping it measured +81 backbone-correct @1%.
    let iso_default = if cli.glyco { (0, 2) } else { (-1, 2) };
    let (iso_min, iso_max) = cli.isotope_error.unwrap_or(iso_default);
    params.isotope_error_range = iso_min..=iso_max;
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
    // `--glyco-index-sequon-only`: glyco-only index shrink (see the flag's doc).
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
    // `--precursor-mono auto`: MS1 isotope-envelope precursor correction. Glyco
    // only (the standard path's -1..2 sweep and PIN are untouched), and it needs
    // MS1, so mzML or Thermo `.raw`; anywhere else it is inert and says so. When
    // active the input is streamed with its MS1 link exactly like `--chimeric`
    // and each chunk's precursors are corrected BEFORE they reach the search.
    let mono_active =
        cli.precursor_mono == PrecursorMonoFlag::Auto && cli.glyco && (is_mzml || is_raw);
    // `auto` is the default and the correction is glyco-only, so a search
    // without --glyco is not told about it; a glyco search without MS1 is.
    if cli.precursor_mono == PrecursorMonoFlag::Auto && cli.glyco && !mono_active {
        eprintln!(
            "WARN: --precursor-mono auto needs an input with MS1 (mzML or Thermo .raw); \
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
    // Fallback isolation half-width (Da) used only when the file omits per-scan
    // isolation offsets; a fixed sensible default (was the --isolation-halfwidth flag).
    params.chimeric_isolation_halfwidth_da = 1.5;
    params.chimeric_max_coisolated = cli.chimeric_max_coisolated;
    params.chimeric_max_kl = cli.chimeric_max_kl;
    // `--chimeric` keeps the same top-N candidates per scan as a plain search.
    // Pass 2 takes the best one as the primary peptide; the alternatives stay in
    // the PIN, so a rescorer can promote a lower-ranked candidate. Measured on
    // the Astral entrapment database, keeping them does not raise the true FDP
    // at Percolator's nominal 1%.
    params.top_n_psms_per_spectrum = cli.top_n;
    params.enzyme = search_enzyme;
    params.extra_enzymes = extra_enzymes.clone();
    params.num_tolerable_termini = proteases.num_tolerable_termini;
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
    if let Some(n) = cli.cal_min_spec_keys {
        params.cal_min_spec_keys = n;
    }
    params.precursor_mass_shift_ppm = 0.0;
    params.refine_select_psm_fdr = cli.refine_select_psm_fdr;
    params.score_mode = match cli.score {
        ScoreFlag::Rank => search::ScoreMode::Rank,
        ScoreFlag::Strong => search::ScoreMode::Strong,
        // Auto: strong re-ranking pays on high-res (sharp fragment peaks → the
        // intensity / rich-ion GBDTs separate target/decoy); on low-res it
        // collapses the candidate pool → use rank. Keyed off the acquisition's
        // resolution class.
        ScoreFlag::Auto => {
            if acquisition.high_res {
                search::ScoreMode::Strong
            } else {
                search::ScoreMode::Rank
            }
        }
    };
    // Every search retrieves its candidates from the out-of-core base-peptide
    // index; the glyco driver materialises its own sequon-only table from it.
    // The retrieval strategy is chosen by `retrieval_choice` and reported once;
    // `--fragment-index` overrides it. Searches the fragment-ion index cannot
    // serve keep per-spectrum enumeration.
    let retrieval_excluded: Option<&'static str> = if cli.glyco {
        Some("--glyco looks candidates up by backbone mass")
    } else if !search::fragment_index::covers_activation(param.data_type.activation) {
        Some(
            "it holds only b/y ions, and this activation method fragments into \
             c/z or a/x ions",
        )
    } else {
        None
    };
    if let (false, Some(reason)) = (cli.glyco, retrieval_excluded) {
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

    // --glyco-hcd-pair only takes effect inside the glyco driver; without --glyco
    // it is silently inert, which would mislead a user who set it on purpose (code
    // review). Warn rather than error so it stays a no-op knob outside glyco mode.
    // The flag defaults on, so only a value typed on the command line warns.
    if !cli.glyco && arg_present("--glyco-hcd-pair") {
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
        // `sequon` restores the N-X-S/T sequon density, so it does NOT trip this
        // warning; only the sequon-depleting strategies do.
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
    // The activation/instrument lines report what the RUN was routed on and say
    // where each value came from; the model's own training regime is shown next
    // to it, so a substituted model (e.g. a high-res HCD model serving a TOF
    // run) is never presented as the detected instrument.
    let model_source = if cli.model_id_override.is_some() {
        "--model"
    } else {
        "selected from store"
    };
    eprintln!("  model          : {resolved_model_id} [{model_source}]");
    let activation_source = if cli.fragmentation != Fragmentation::Auto {
        "--fragmentation"
    } else if detected_activation_instrument.is_some() {
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

    eprintln!("Using {} worker threads", cli.threads);

    // `PreparedSearch::fragment_tolerance_da`: the Da window of the glyco
    // peptide-first fragment index on low-res data (high-res data uses the
    // ppm retrieval window). Node scoring takes the argument but does not
    // read it, and site evidence uses the scorer's PIN tolerance.
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

    // Calibration pre-pass on the out-of-core index. The index is
    // content-addressed and precursor-tolerance independent, so the main pass
    // reopens the same cached file the pre-pass built.
    let cache_dir = crate::cli::dir_or_output_dir(cli.index_cache_dir.as_deref(), &report_base);
    let index_path = index_cache_path(&idx, &params, &cache_dir);
    if params.precursor_cal_mode != PrecursorCalMode::Off {
        let t_cal = std::time::Instant::now();
        let mut cal_prepared = PreparedSearch::prepare_mmap(
            &idx,
            &params,
            &scorer,
            fragment_tol_da,
            &cli.decoy_prefix,
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

    // Content-addressed cache path, so repeated searches over the same FASTA +
    // params (and the calibration pre-pass just above) share one on-disk index
    // without rebuilding.
    let t_index = std::time::Instant::now();
    let mut prepared = PreparedSearch::prepare_mmap(
        &idx,
        &params,
        &scorer,
        fragment_tol_da,
        &cli.decoy_prefix,
        &index_path,
    )
    .map_err(|e| format!("build out-of-core candidate index: {e}"))?;
    eprintln!(
        "[PHASE candidate_index_open: {:.2}s]",
        t_index.elapsed().as_secs_f64()
    );
    eprintln!(
        "PreparedSearch: out-of-core candidate-index: mmap \
         (base peptides resolved lazily per spectrum)"
    );

    let bench_mode = cli.max_spectra > 0;

    // ── Quantification setup (`--tmt` / `--lfq`) ──────────────────────────────
    // Reporter ions are read from each MS2 while its peaks are still in memory;
    // MS1 scans (purity, label-free chromatograms) and MS3 scans (SPS-MS3
    // reporters) are captured by the readers next to the MS2 stream and handed
    // to the collector per file. Quantification itself runs after rescoring.
    let high_res = acquisition.high_res;
    // Without --tmt-level, SPS-MS3 reporters are read from the MS3 when the
    // low-resolution MS2 cannot separate the plex's channels.
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
    let mut quant: Option<QuantCollector> = QuantSettings::from_cli(&cli, high_res, tmt_level)?
        .map(|mut settings| {
            // Recorded in each run's MS2 table, where the library matches
            // its fragment spectra and extract scores transfers.
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
        if is_d && (s.needs_ms1_series() || s.needs_ms3()) {
            return Err(
                "--lfq, --quant-dir and --tmt-level 3 are not yet supported for Bruker .d input"
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
            return Err("--tmt/--lfq/--quant-dir are not available in --glyco mode".into());
        }
        if !cli.rescore && !cli.rescore_native && (s.plex.is_some() || s.lfq) {
            eprintln!(
                "WARN: --tmt/--lfq without --rescore or --rescore-native: every rank-1 target \
                 PSM is quantified; filter the tables by the q-values you compute downstream, \
                 or add --rescore so only PSMs at q<=--fdr are quantified."
            );
        }
        eprintln!("Quantification: {}", s.describe());
    }
    // `--predict`: load the predictors now so a missing model fails before the
    // search. Glyco mode writes its own `.glyco.pin` and takes none.
    let predictors = if cli.glyco {
        None
    } else {
        crate::predicted::Predictors::from_cli(&cli, &acquisition, param.data_type.instrument)?
    };
    if let Some(p) = &predictors {
        eprintln!("Predict: {}", p.describe());
    }
    // Peaks normally go once a spectrum is scored; `--refine` re-scores in its
    // Pass 2 and `--predict ms2` compares them with predicted spectra.
    let keep_peaks = cli.refine || predictors.as_ref().is_some_and(|p| p.needs_peaks());
    let mut all_spectra: Vec<Spectrum> = Vec::new();
    let mut all_queues: Vec<TopNQueue> = Vec::new();
    // Exclusive end index in `all_spectra` of each input, in input order.
    let mut run_ends: Vec<usize> = Vec::with_capacity(spectrum_paths.len());

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
    let ms1_input = chimeric_active || mono_active;
    let mut parse_stats = ParseStats::default();
    // `--precursor-mono`: one entry per accumulated spectrum (aligned with
    // `all_spectra`); stays empty when the flag is off.
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

    // The per-run fragment mass error (high-res acquisitions without an
    // explicit `--fragment-tol-*`; the glyco driver keeps its own tolerances).
    let mut measurer = FragmentErrorMeasurer {
        enabled: acquisition.high_res && frag_tol_override.is_none() && !cli.glyco,
        done: false,
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

        if ms1_input && !(file_is_mzml || file_is_raw) {
            return Err(format!(
                "{} only supports mzML/.raw inputs, got {}",
                if chimeric_active {
                    "--chimeric"
                } else {
                    "--precursor-mono"
                },
                input_path.display()
            )
            .into());
        }

        let mut run_scans = RunScans::default();
        let file_stats = if ms1_input {
            let (tx, rx) = sync_channel::<(Vec<Spectrum>, Ms1Link)>(2);
            let spectrum_path = input_path.clone();
            let cap = remaining_cap;
            // The cascade is MS2-only by construction (MS2 paired with its preceding
            // MS1); hardcode MS2 so `--ms-level 3` can never widen the mzML reader's
            // range to admit MS3 (the .raw chunked reader is already MS2-only).
            let mslevel = 2;
            // Label-free quantification reads the run's MS1 scans and each MS2's
            // context from this path too (the reader keeps them next to the
            // per-chunk co-isolation link).
            let capture_ms1 = quant.as_ref().is_some_and(|q| q.settings.needs_ms1());
            let capture_context = quant.is_some();
            let capture_noise = quant.as_ref().is_some_and(|q| q.settings.plex.is_some());
            type Ms1Parsed = (usize, Vec<String>, RunScans);
            let parser_handle = thread::spawn(move || -> Result<Ms1Parsed, String> {
                if file_is_mzml {
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
            // Fragment-index retrieval scores MASS-ORDERED chunks (see the
            // non-MS1 branch below). Each spectrum keeps the part of its linked
            // MS1 scan Pass 2 reads, so a mass-ordered chunk can carry its own
            // MS1 link.
            let index_ordered = chimeric_active && params.fragment_index_top_k > 0 && !cli.glyco;
            let mut pending: Vec<(Spectrum, Option<Ms1Peaks>)> = Vec::new();
            for (mut chunk_spectra, chunk_link) in rx {
                if let Some(prefix) = &title_prefix {
                    prefix_spectrum_titles(&mut chunk_spectra, prefix);
                }
                // `--precursor-mono`: correct precursors from the linked MS1
                // BEFORE anything scores the chunk, so every consumer (the glyco
                // driver, the PIN's ExpMass) sees the corrected m/z.
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
                } else if !cli.glyco {
                    // Unreachable today (`mono_active` implies `--glyco`), kept so
                    // an MS1-linked stream never silently skips the peptide search.
                    let queues = prepared.run_chunk(&chunk_spectra, offset);
                    all_queues.extend(queues);
                }
                file_offset += chunk_spectra.len();
                ms1_linked += chunk_link.ms1_peaks.len();
                for mut spec in chunk_spectra.into_iter() {
                    // Peaks are normally dropped post-scoring to bound memory
                    // (only the metadata is needed downstream). Under `--refine`
                    // or `--glyco` we RETAIN them: refine's Pass-2 re-scores
                    // unidentified spectra; glyco_search_run needs the full
                    // peak lists for oxonium ion detection. Memory cost: full
                    // peak buffer stays resident (acceptable; both are opt-in).
                    if let Some(q) = quant.as_mut() {
                        q.on_spectrum(&spec);
                    }
                    if !keep_peaks && !cli.glyco {
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
                // Pass 1 per mass-ordered chunk; Pass 2 once over the whole file
                // afterwards. A chunk's co-isolated precursors are spread over the
                // whole mass range, so a per-chunk Pass 2 would expand most of the
                // database once per chunk; one pass expands it once. Peaks stay
                // resident until then (they already are, for the mass sort).
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
            && !cli.glyco
            && file_is_mzml
            && !is_gzip_path(input_path)
            && !quant
                .as_ref()
                .is_some_and(|q| q.settings.needs_ms1() || q.settings.needs_ms3())
        {
            // Fragment-index mode scores mass-ordered chunks, so nothing is
            // scored until the whole file is read; read it with every thread.
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
                            let reader = input::TimsTofReader::open(&spectrum_path)
                                .map_err(|e| format!("open Bruker .d: {e}"))?;
                            Ok((
                                send_chunks(reader, CHUNK_SIZE, remaining_cap, tx),
                                RunScans::default(),
                            ))
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

            // Fragment-index mode (issue #76) scores MASS-ORDERED chunks: the
            // per-chunk index covers the peptidoforms of the chunk's precursor
            // windows, and a chunk in file order spans the whole mass range
            // (its index is the whole database; 100 GB on phospho). Collect the
            // stream, sort by precursor mass, then chunk. Peaks for every
            // spectrum are held until scored (~1 GB for 100k high-res MS2).
            let mass_ordered = params.fragment_index_top_k > 0 && !cli.glyco;
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
                // SPEED (--glyco): the per-spectrum PEPTIDE search is pure waste in
                // glyco mode — glyco_search_run re-derives its own candidates from
                // `prepared` + the retained peaks and writes its own `.glyco.pin`
                // (the standard PIN is skipped), so these PSM queues are discarded.
                // Skip the scoring (parsing + peak retention below still run); ~16%
                // faster on glyco. Non-glyco path is unchanged.
                if !cli.glyco {
                    measurer.measure(&prepared, chunk.iter(), None);
                    let offset = all_spectra.len();
                    let queues = prepared.run_chunk(&chunk, offset);
                    all_queues.extend(queues);
                }
                for mut spec in chunk.into_iter() {
                    // See the chimeric-loop note above: normally peaks are
                    // dropped post-scoring to bound memory, but `--refine` or
                    // `--glyco` needs the full peak lists retained.
                    if let Some(q) = quant.as_mut() {
                        q.on_spectrum(&spec);
                    }
                    if !keep_peaks && !cli.glyco {
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
    // Inputs skipped by the bench cap contribute no spectra.
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
    if bench_mode {
        eprintln!("Bench mode: capped at {} spectra", cli.max_spectra);
    }

    // Drain the per-spectrum materialized candidates accumulated during the scan
    // into `prepared.candidates` so the PIN/TSV writers resolve every PSM's
    // `candidate_idxs` against a real candidate slice.
    prepared.sync_materialized_candidates();

    // The glyco PIN carries the Mono columns whenever the correction ran, fitted or
    // not, so its header follows the command line and the input format, never the
    // data: the PINs of separately searched runs share one header.
    let mono_columns = mono_active;
    // No MS1 reached the correction (an MS2-only mzML, or no charges): nothing was
    // fitted, so the search runs as with the flag off and the Mono columns are 0.
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
    // Corrected precursors search a 0..=1 window: the `+2` step only carries the
    // Hex+Fuc/NeuGc composition degeneracy once the monoisotope is right
    // (bigbio/andes#64 arm F; rationale on `coupled_isotope_window`). Explicit
    // window flags win, and nothing changes unless an MS2 was actually fitted.
    // The glyco scorer reads its window from the prepared index's params, which
    // are borrowed for the rest of the run, so the narrowed window travels as a
    // `GlycoConfig` override instead of a mutation.
    let (window, narrowed) = crate::mono::coupled_isotope_window(
        mono_active,
        cli.isotope_error,
        params.isotope_error_range.clone(),
    );
    let glyco_isotope_override = if narrowed {
        // Only a spectrum whose envelope fitted, with its best fit inside the
        // narrowed window, has a verified monoisotope; the rest keep the wide
        // window they would have had without the flag.
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
            "Search complete: {non_empty} / {} spectra have PSMs (search wall: {:.2}s)",
            spectra.len(),
            search_elapsed.as_secs_f64()
        );
    }

    // ── 7a. Glyco mode: run glyco scoring and write .glyco.pin, then return ──
    // When --glyco is active, we run the glyco-PSM scoring driver over ALL
    // accumulated spectra (candidates from a sequon-only table built from the
    // out-of-core index) and write a separate `.glyco.pin` file.  The standard PIN is skipped.
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
    // compete in one report. Pass-2 rows carry `IsRefinement = 1`: one threshold
    // over the whole PIN lets that small group run at about twice the nominal FDR
    // (Astral entrapment: ~2% at a merged 1%), so for a 1% FDR on modified PSMs,
    // threshold the `IsRefinement = 1` rows on their own downstream.
    let refine_output = if cli.refine {
        // Refinement config: explicit YAML or the built-in 5-mod default tier.
        let base_cfg = match &cli.refine_config {
            Some(p) => search::RefineConfig::from_yaml_str(&std::fs::read_to_string(p)?)
                .map_err(|e| format!("parsing --refine-config {}: {e}", p.display()))?,
            None => search::RefineConfig::default_tier(),
        };
        // Max variable mods and the entrapment diagnostic for refinement come from
        // the refine config/tier (set them in the `--refine-config` YAML).
        let cfg = base_cfg;

        // High-res signal: the acquisition's resolution class. High-res
        // fragments match in ppm (20 ppm vs 0.5 Da ion-trap), which is exactly
        // the regime where the near-isobaric refinement deltas (e.g.
        // deamidation +0.984 vs a C13 isotope error) are resolvable.
        let high_res = acquisition.high_res;

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

    // Engine-wide retention-time features (additive): after all PSMs are scored
    // and the final PIN candidate pool is fixed, run per-run RT self-calibration
    // and populate DeltaRT/AbsDeltaRT/DeltaRTNorm (+ predicted_rt for the QPX)
    // onto each PSM. Neutral 0.0 (baseline-identical) when observed RT is
    // unavailable or calibration cannot be fit.
    output::populate_rt_features(&spectra, &mut queues, &pin_candidates);

    // Provenance (additive): entry class, gap to the best Reference candidate,
    // isotope-mimic flag. All 0 on a FASTA search except IsotopeMimic with
    // isotope-like variable mods.
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
    let proteoform_fdr = cli.proteoform_fdr.unwrap_or(fdr_threshold.max(0.05));
    let (use_percolator, use_native) = if cli.rescore {
        (true, false)
    } else if cli.rescore_native {
        (false, true)
    } else {
        // Say so. A threshold that silently does nothing is the same defect shape
        // as a flag that is parsed and never read.
        if cli.fdr.is_some() || cli.pep.is_some() {
            let which = match (cli.fdr.is_some(), cli.pep.is_some()) {
                (true, true) => "--fdr and --pep",
                (true, false) => "--fdr",
                _ => "--pep",
            };
            eprintln!(
                "WARN: {which} set without --rescore or --rescore-native, so no rescoring \
                 runs and no filtered .tsv is written. andes does not compute FDR on its \
                 own; the PIN is written for Percolator. Add --rescore to run Percolator \
                 in-process."
            );
        }
        (false, false)
    };
    // `--mbr`: the decoy PSMs with their scores, which the target map leaves
    // out, and where the rescoring result came from.
    let mut mbr_decoys: std::collections::HashMap<String, output::PercolatorPsm> =
        std::collections::HashMap::new();
    let mut mbr_source = String::new();
    // The rescorer's decoy-result ids: with target-decoy competition (Percolator
    // `-Y`, and the native rescorer) they let the QPX join tell a target that
    // lost its scan/ExpMass group from one whose join drifted.
    let mut rescore_decoys: Option<std::collections::HashSet<String>> = None;
    // The rescorer's decoy PSMs with their scores, for FDR within each class.
    let mut rescore_decoy_psms: std::collections::HashMap<String, output::PercolatorPsm> =
        std::collections::HashMap::new();
    let mut rescore_map: Option<std::collections::HashMap<String, output::PercolatorPsm>> =
        if use_percolator {
            let backend = output::resolve_backend(
                cli.percolator_bin.as_deref(),
                cli.percolator_docker,
                &cli.percolator_image,
            )
            .map_err(|e| format!("{e}"))?;
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
            if cli.mbr {
                mbr_source = format!("the Percolator results of {}", output_pin_path.display());
            }
            Some(map)
        } else if use_native {
            eprintln!(
                "Rescore: native rescorer (no Percolator; 3-fold cross-validated semi-supervised \
                 linear model, target-decoy competition per scan and precursor mass). \
                 Use --rescore (Percolator) for production FDR."
            );
            let t_resc = std::time::Instant::now();
            let pin_text = std::fs::read_to_string(&output_pin_path)
                .map_err(|e| format!("reading PIN for native rescore: {e}"))?;
            let all = rescore::native_rescore_pin_all(&pin_text, 42)
                .map_err(|e| format!("native rescore: {e}"))?;
            if all.untrained_folds > 0 {
                eprintln!(
                    "WARN: native rescore: {} of 3 cross-validation folds had too few confident \
                     targets or decoys to train on and ranked by the best single PIN feature. \
                     The q-values stay valid but identify fewer PSMs.",
                    all.untrained_folds
                );
            }
            // Like Percolator with -Y: the winning decoys tell the QPX join
            // which target rows lost their scan/ExpMass group.
            rescore_decoys = Some(all.decoys.keys().cloned().collect());
            rescore_decoy_psms = all.decoys;
            if cli.mbr {
                mbr_source = "the native rescoring result".to_string();
            }
            let map = all.targets;
            eprintln!(
                "Rescore: native produced {} target PSMs (q-value+PEP; {} rows lost \
                 target-decoy competition) [{:.2}s]",
                map.len(),
                all.competed_away,
                t_resc.elapsed().as_secs_f64()
            );
            Some(map)
        } else {
            None
        };

    // FDR within each evidence class (canonical, isoform, variant, chain): a
    // target's q-value becomes the larger of the global and its class q-value,
    // so the filtered TSV, the QPX and quantification accept a PSM only if it
    // passes both. A search with one class (any FASTA) is left unchanged.
    // Knockoff rows: in the PIN as targets; left out of the QPX bundle, the
    // filtered TSV, quantification and every other output.
    let mut knockoff_rows: std::collections::HashSet<String> = std::collections::HashSet::new();
    if let Some(ref mut map) = rescore_map {
        // A variant peptide a reference protein holds once I and L are taken
        // as equal identifies that reference too: share it before classes are
        // assigned. Only rows that could be reported are checked.
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
        // Knockoff variants: the knockoff PSMs accepted per substitution class
        // estimate the false real variant PSMs, and give them a q-value too.
        let outcome =
            output::apply_knockoff_fdr(map, &pin_index, proteoform_fdr, cli.knockoffs > 0);
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
        // Site claims: annotated sites against their knockoff sites; each
        // accepted site also gets its occupancy and identity from the
        // posterior of every row among its spectrum's rows.
        let has_site_claims = pin_index.iter_target_proteins().any(|p| {
            p.annotation
                .as_deref()
                .is_some_and(|a| !a.site_claims.is_empty())
        });
        let posteriors = has_site_claims.then(|| {
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
            output::site_identity::row_posteriors(map, &rescore_decoy_psms, &rows)
        });
        let sites = output::site_claims(
            map,
            &pin_index,
            fdr_threshold,
            proteoform_fdr,
            &restricted_site_masses,
            posteriors.as_ref(),
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
    }
    // Protein groups: parsimony per protein family and picked target-decoy
    // q-values, from every rescored row. Independent of quantification.
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
    // The cutoff each row is accepted at: --proteoform-fdr for proteoform
    // claims (variant-class rows), --fdr for the rest.
    let cutoffs = rescore_map
        .as_ref()
        .map(|m| output::Cutoffs::new(m, &pin_index, fdr_threshold, proteoform_fdr));

    // Filtered TSV: target PSMs whose q-value passes their cutoff. Written next
    // to the PIN as `<stem>.q<fdr>.tsv` (rescore only). One row per accepted PSM
    // with the join key + q/PEP + peptide/proteins from the Percolator result.
    if let (Some(map), Some(cutoffs)) = (rescore_map.as_ref(), cutoffs.as_ref()) {
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
        accepted.sort_by(|a, b| {
            a.q_value
                .partial_cmp(&b.q_value)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
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
    }

    // ── Run summary: final tolerances + per-modification PSM tally ────────────
    // andes auto-resolves the scoring model and the precursor/fragment
    // tolerances from the data, so the FINAL search parameters differ from
    // whatever the user passed on the CLI. Emit a summary to stderr and write a
    // `<stem>.statistics.log` next to the PIN so a run's true parameters (and which
    // PTMs were identified, with PSM counts) are always recoverable.
    let run_stats = output::RunStatistics::compute(&queues, &pin_candidates, &params, &param.mme);
    eprint!("{}", run_stats.render());
    // Named after the run like the other reports, so searches writing into
    // one directory do not overwrite each other's summary.
    let stats_path = {
        let stem = report_base
            .file_stem()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| "andes".to_string());
        let dir = report_base.parent().unwrap_or_else(|| Path::new("."));
        dir.join(format!("{stem}.statistics.log"))
    };
    match output::write_statistics_log(&run_stats, &stats_path) {
        Ok(()) => eprintln!("Wrote statistics: {}", stats_path.display()),
        Err(e) => eprintln!("WARN: could not write {}: {e}", stats_path.display()),
    }

    // ── QPX `.idparquet/` bundle (optional, quantms-compatible) ──────────────
    // Mirrors the PIN write from the SAME unified candidate/index pair. Sourced
    // here (not in the bench/TSV arms) so it is produced for every run that asks
    // for it, including bench mode. Run identifier = the first spectrum file
    // stem so it is stable across re-runs of the same input; each row's
    // `reference_file_name` is its own input's path.
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
            &run_ends,
            rescore_map.as_ref(),
            rescore_decoys.as_ref(),
            &knockoff_rows,
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

    if let Some(ref mut map) = rescore_map {
        map.retain(|id, _| !knockoff_rows.contains(id));
    }

    // ── Quantification (`--tmt` / `--lfq`) ────────────────────────────────────
    // After rescoring so the targets are the confident PSMs; after the QPX
    // write so the feature parquet joins a bundle that already exists.
    // The QPX views keep the protein groups at `--protein-fdr`, 0.01 unless
    // given.
    let qpx_groups = protein_groups
        .as_ref()
        .map(|inf| crate::quant_cmd::ProteinGroups {
            inf,
            max_q: cli.protein_fdr.unwrap_or(0.01),
        });
    let quant_summary = match &quant {
        Some(q) => Some(crate::quant_run::run_quant(
            q,
            &spectra,
            &queues,
            &pin_candidates,
            &pin_index,
            rescore_map.as_ref(),
            cutoffs.as_ref(),
            &report_base,
            cli.output_parquet.as_deref(),
            qpx_groups.as_ref(),
        )?),
        None => None,
    };

    // Proteoform table: per protein family, the canonical entry and each
    // isoform, variant and chain entry with accepted evidence, its confidence
    // and quantity. Written with rescoring when the database holds such entries.
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

    // `--mbr` inputs; the stages run after the last search output is
    // written and the search state is released.
    let mbr_job = if cli.mbr {
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
            .map(|q| (q.settings.quant_fdr, q.settings.pep_cap))
            .ok_or("--mbr needs --lfq")?;
        Some((dir, stems, scores, ep, donor_q, donor_pep))
    } else {
        None
    };
    drop(mbr_decoys);

    // ── 9. Write TSV (optional) ───────────────────────────────────────────────
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

    // ── Match-between-runs (`--mbr`): library, extract, finalize ────────────
    // Every search output is written; the search state is freed first, so
    // the stages hold one run's MS1 at a time next to the library.
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
        cli.mbr_library_fdr,
        donor_q,
        donor_pep,
        cli.mbr_fdr,
        cli.quant_feature_fdr.0,
        &report_base.with_extension(""),
        cli.output_parquet.as_deref(),
        qpx_groups.as_ref(),
    )?;
    Ok(())
}

/// Log what a UniProt XML load kept: annotated sites, variant entries, and
/// the features the evidence level or the canonical sequence ruled out.
/// Log what an annotation source (`--annotations` files or a PEFF header)
/// added to the database.
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

/// Say how many annotated sites each `sites=annotated` mod of `aa` can use,
/// and warn about one the database annotates no site of its mass for: it is
/// never placed.
/// Artefact masses that mimic common substitutions: deamidation (N→D, Q→E),
/// oxidation (A→S, P→T…), carbamylation, methylation (G→A, S→T, D→E…).
const ARTEFACT_MASSES: [f64; 4] = [0.984016, 15.994915, 43.005814, 14.015650];

/// With natural variants searched, restricting one of these mods to annotated
/// sites removes the canonical explanation that competes with a variant of the
/// same mass.
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
            // Annotated sites only, not the knockoff sites added for them.
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
                "WARN: {name} ({mass:+.4} Da) is restricted to annotated sites (sites=annotated), \
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

/// `--mbr` needs two or more inputs, `--lfq` and a rescoring run (the
/// library ranks precursors by the rescored PSM scores), and is label-free
/// only.
fn check_mbr_args(cli: &Cli) -> Result<(), String> {
    if cli.spectrum.len() < 2 {
        return Err(format!(
            "--mbr transfers identifications between runs and needs two or more --spectrum inputs (got {})",
            cli.spectrum.len()
        ));
    }
    if !cli.lfq {
        return Err("--mbr needs --lfq".into());
    }
    if !cli.rescore && !cli.rescore_native {
        return Err(
            "--mbr needs --rescore (or --rescore-native): the library ranks precursors by the rescored PSM scores"
                .into(),
        );
    }
    if cli.tmt.is_some() {
        return Err("--mbr is label-free and cannot be combined with --tmt".into());
    }
    crate::quant_cmd::evidence_params(cli.lfq_tol)?;
    Ok(())
}

/// The candidate pool is millions of small heap allocations and the process
/// exits right after the run; freeing them one by one took about a second of
/// wall time, so the pool is left for the operating system to reclaim.
fn release_at_exit(candidates: Vec<search::candidate_gen::Candidate>) {
    std::mem::forget(candidates);
}

/// Write the FDR-filtered rescore TSV: one row per accepted target PSM (already
/// filtered to q ≤ fdr and sorted by q ascending). Columns: the Percolator join
/// key (`PSMId`, == PIN SpecId), `q-value`, `posterior_error_probability`,
/// `peptide`, `proteinIds`.
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

/// TMTpro 32/35-plex peptides are searched with the TMTpro tag mass, but each
/// deuterated tag is 2.9 mDa heavier. Below three such tags on a 1000 Da
/// peptide (8.8 ppm), the precursor window misses the deuterated precursors.
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

/// The detected-label class of a `--tmt` kit.
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

// ── Training pipeline ─────────────────────────────────────────────────────────

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

    /// Low-res searches use the index with the vote cap: enumeration
    /// there took TMT 5,404 s of CPU against 289 s for the capped index, at the
    /// same PSMs.
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

    /// The index holds b/y ions only: ETD and UVPD keep enumeration at either
    /// resolution, and forcing the index there is refused with the reason.
    #[test]
    fn activations_outside_the_indexed_ion_series_keep_enumeration() {
        use model::activation::ActivationMethod;
        for a in [
            ActivationMethod::CID,
            ActivationMethod::HCD,
            ActivationMethod::PQD,
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

    /// A tolerance tight enough to bin usefully counts as high-resolution (no
    /// vote cap) even when it is expressed in Da.
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
        let (pin, base) = rescore_pin_destination(None, None, None, true, Path::new("."), &mut tmp).unwrap();
        assert_eq!(pin, PathBuf::from("andes.pin"));
        assert_eq!(base, PathBuf::from("andes.pin"));
        assert!(tmp.is_none());
    }

    /// Without `--keep-pin` the PIN lives in a directory under the scratch
    /// directory that was passed, never in the system temp directory.
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

    /// A scratch directory that cannot be created is an error naming the flag.
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

    /// The filtered TSV is written through `<path>.partial` and renamed, so a
    /// run that dies mid-write never leaves a truncated table under the name.
    #[test]
    fn filtered_tsv_is_written_through_a_partial_file() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("run.q0p01.tsv");
        let psm = output::PercolatorPsm {
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
        // Low-res: left unset, so the protocol engages the 100 Da / top-20 filter.
        assert_eq!(isobaric_peak_filter(None, false), None);
        // High-res: a zero window disables it.
        assert_eq!(isobaric_peak_filter(None, true), Some((0.0, 0)));
        // An explicit --peak-filter wins on either class.
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
