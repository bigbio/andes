//! msgf-rust: end-to-end MS-GF+ search.
//!
//! Loads an MGF or mzML spectrum file and a FASTA target database, runs a
//! tryptic database search with default MS-GF+ parameters, and writes output
//! in Percolator `.pin` format (and optionally `.tsv` format).
//!
//! Format dispatch: if `--spectrum` ends in `.mzML` or `.mzml`, `MzMLReader`
//! is used; otherwise `MgfReader` is used (default / backwards-compatible).

use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;
use std::process::ExitCode;

use clap::Parser;
use model::{AminoAcidSetBuilder, ModLocation, Modification, PrecursorTolerance, ResidueSpec, Spectrum, Tolerance};
use scoring_crate::{Param, RankScorer};
use search::{PreparedSearch, SearchIndex, SearchParams, TopNQueue};
use input::{FastaReader, MgfReader, MzMLReader};

#[derive(Parser, Debug)]
#[command(
    name = "msgf-rust",
    about = "MS-GF+ Rust port: database search of MGF/mzML spectra against FASTA"
)]
struct Cli {
    /// Input spectrum file (MGF or mzML). Format is auto-detected by extension:
    /// `.mzML`/`.mzml` → MzMLReader; anything else → MgfReader.
    #[arg(long)]
    spectrum: PathBuf,

    /// Input FASTA database (target sequences only; decoys are generated automatically).
    #[arg(long)]
    database: PathBuf,

    /// Output Percolator PIN file path.
    #[arg(long)]
    output_pin: PathBuf,

    /// Output TSV file path (optional).
    #[arg(long)]
    output_tsv: Option<PathBuf>,

    /// Decoy prefix used when generating reversed decoy sequences.
    #[arg(long, default_value = "XXX_")]
    decoy_prefix: String,

    /// Minimum isotope error offset to try (default -1).
    #[arg(long, default_value = "-1")]
    isotope_error_min: i8,

    /// Maximum isotope error offset to try (default 2).
    #[arg(long, default_value = "2")]
    isotope_error_max: i8,

    /// Precursor mass tolerance in ppm (default 20.0).
    #[arg(long, default_value = "20.0")]
    precursor_tol_ppm: f64,

    /// Minimum precursor charge to try when not specified in the spectrum.
    #[arg(long, default_value = "2")]
    charge_min: u8,

    /// Maximum precursor charge to try when not specified in the spectrum.
    #[arg(long, default_value = "3")]
    charge_max: u8,

    /// Maximum number of PSMs to retain per spectrum.
    #[arg(long, default_value = "10")]
    top_n: u32,

    /// Number of Tolerable Termini.
    ///
    /// Controls enzymatic-cleavage enforcement at span boundaries:
    ///   2 (default): both termini must be cleavage sites (strict / fully specific).
    ///   1: at least one terminus must be a cleavage site (semi-specific).
    ///   0: neither terminus needs to be a cleavage site (non-specific).
    #[arg(long, default_value = "2")]
    ntt: u8,

    /// Maximum number of missed cleavages per peptide (default 1).
    #[arg(long, default_value = "1")]
    max_missed_cleavages: u32,

    /// Minimum number of peaks required in an MS2 spectrum to attempt scoring.
    ///
    /// Spectra with fewer peaks are skipped (default 10).
    #[arg(long, default_value = "10")]
    min_peaks: u32,

    /// Minimum peptide length (in residues) to consider during the search.
    /// Default 6.
    #[arg(long, default_value = "6")]
    min_length: u32,

    /// Maximum peptide length (in residues) to consider during the search.
    /// Default 40.
    #[arg(long, default_value = "40")]
    max_length: u32,

    /// Path to the .param scoring model file.
    ///
    /// If not supplied, a bundled file under
    /// `src/main/resources/ionstat/` is selected from
    /// `(--fragmentation, --instrument, --protocol)` (default
    /// `HCD_QExactive_Tryp.param`). When running the binary outside the source
    /// tree this path may not exist; supply --param-file explicitly in that
    /// case.
    #[arg(long)]
    param_file: Option<PathBuf>,

    /// Path to a Java-format mods.txt file describing fixed and variable
    /// modifications. Format: each non-comment line is
    /// `<mass>,<aa>,<fix|opt>,<location>,<name>`, where:
    ///   - `<mass>` is a numeric monoisotopic mass delta (Da). Composition
    ///     strings (e.g. `C2H3N1O1`) are **not** yet supported.
    ///   - `<aa>` is a single uppercase letter or `*` (wildcard).
    ///   - `<location>` is one of `any|N-term|C-term|Prot-N-term|Prot-C-term`.
    /// A single `NumMods=N` line sets the max variable mods per peptide.
    /// Inline `#`-comments are stripped. Blank lines and full-line `#`-comments
    /// are ignored. When omitted, the binary uses its built-in defaults
    /// (Carbamidomethyl-C fixed, Oxidation-M variable).
    #[arg(long = "mod", value_name = "MODFILE")]
    mod_file: Option<PathBuf>,

    /// Fragmentation method index (Java's `-m`):
    ///   0=Auto/CID (default), 1=CID, 2=ETD, 3=HCD, 4=UVPD.
    /// Used to choose the bundled .param file when --param-file is not given.
    #[arg(long, value_name = "ID")]
    fragmentation: Option<u8>,

    /// Instrument type index (Java's `-inst`):
    ///   0=LowRes (default), 1=HighRes, 2=TOF, 3=QExactive.
    /// Used to choose the bundled .param file when --param-file is not given.
    #[arg(long, value_name = "ID")]
    instrument: Option<u8>,

    /// Protocol index (Java's `-protocol`):
    ///   0=Automatic (default), 1=Phosphorylation, 2=iTRAQ,
    ///   3=iTRAQPhospho, 4=TMT, 5=Standard.
    /// Used to choose the bundled .param file when --param-file is not given.
    #[arg(long, value_name = "ID")]
    protocol: Option<u8>,

    /// Number of worker threads for the search loop. Defaults to logical CPU count.
    #[arg(long, default_value_t = num_cpus::get())]
    threads: usize,

    /// Bench mode: process only the first N MS2 spectra and skip writing
    /// PIN/TSV. Use for fast Fix B iteration (1k-2k spectra ≈ <1 min vs
    /// 70 min on full PXD001819). When 0 (default) the full input is used.
    #[arg(long, default_value = "0")]
    max_spectra: usize,

    /// MS level to search. Default 2 (MS2). MS1 spectra (and any other levels)
    /// in the input file are filtered out at load time so they never enter
    /// the search loop or consume RAM. Only meaningful for mzML inputs — MGF
    /// files do not encode MS level and are treated as MS2 regardless.
    #[arg(long, default_value = "2")]
    ms_level: u8,
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    match run(cli) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("msgf-rust: {e}");
            ExitCode::from(1)
        }
    }
}

fn run(cli: Cli) -> Result<(), Box<dyn std::error::Error>> {
    let t_total = std::time::Instant::now();
    let t_phase = std::time::Instant::now();
    // ── 1. Load FASTA target database ────────────────────────────────────────
    let target_db =
        FastaReader::load_all(BufReader::new(File::open(&cli.database)?))?;
    eprintln!(
        "Loaded {} target proteins from {} [PHASE fasta_load: {:.2}s]",
        target_db.proteins.len(),
        cli.database.display(),
        t_phase.elapsed().as_secs_f64()
    );

    // ── 2. Build SearchIndex (target + reversed decoys) ───────────────────────
    let t_phase = std::time::Instant::now();
    let idx = SearchIndex::from_target_db(&target_db, &cli.decoy_prefix);
    eprintln!("[PHASE search_index_build: {:.2}s]", t_phase.elapsed().as_secs_f64());

    // ── 3. Build AminoAcidSet ────────────────────────────────────────────────
    //
    // If --mod is given, parse the Java-format mods.txt file. Otherwise
    // fall back to msgf-rust's historical defaults (CAM fixed on C,
    // Oxidation variable on M) so existing tests keep their behaviour.
    //
    // `num_mods_from_file` is populated only when --mod is given and the
    // file contains a `NumMods=N` line; it overrides the default
    // `max_variable_mods_per_peptide` (3) below.
    let (aa, num_mods_from_file) = match &cli.mod_file {
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
        None => {
            let cam = Modification {
                name: "Carbamidomethyl".into(),
                mass_delta: 57.02146,
                residue: ResidueSpec::Specific(b'C'),
                location: ModLocation::Anywhere,
                fixed: true,
                accession: None,
            };
            let ox = Modification {
                name: "Oxidation".into(),
                mass_delta: 15.99491,
                residue: ResidueSpec::Specific(b'M'),
                location: ModLocation::Anywhere,
                fixed: false,
                accession: None,
            };
            let set = AminoAcidSetBuilder::new_standard()
                .add_fixed_mod(cam)
                .add_variable_mod(ox)
                .build()?;
            (set, None)
        }
    };

    // ── 4. Load Param scoring model ───────────────────────────────────────────
    let param_path = match cli.param_file.clone() {
        Some(p) => p,
        None    => resolve_bundled_param(cli.fragmentation, cli.instrument, cli.protocol)?,
    };
    eprintln!("Param file: {}", param_path.display());

    let t_phase = std::time::Instant::now();
    let param = Param::load_from_file(&param_path)
        .map_err(|e| format!("loading param file {}: {e}", param_path.display()))?;
    let scorer = RankScorer::new(&param);
    eprintln!("[PHASE param_and_scorer: {:.2}s]", t_phase.elapsed().as_secs_f64());

    // ── 5. Build SearchParams ─────────────────────────────────────────────────
    let mut params = SearchParams::default_tryptic(aa);
    params.precursor_tolerance =
        PrecursorTolerance::symmetric(Tolerance::Ppm(cli.precursor_tol_ppm));
    params.charge_range = cli.charge_min..=cli.charge_max;
    params.isotope_error_range = cli.isotope_error_min..=cli.isotope_error_max;
    params.top_n_psms_per_spectrum = cli.top_n;
    params.num_tolerable_termini = cli.ntt;
    params.max_missed_cleavages = cli.max_missed_cleavages;
    params.min_peaks = cli.min_peaks;
    params.min_length = cli.min_length;
    params.max_length = cli.max_length;
    if let Some(n) = num_mods_from_file {
        params.max_variable_mods_per_peptide = n;
    }

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

    // Fragment tolerance of 0.5 Da matches the gf_bsa_parity integration test
    // (and the canonical HCD default).
    let fragment_tol_da = 0.5_f64;
    let prepared = PreparedSearch::prepare(
        &idx,
        &params,
        &scorer,
        fragment_tol_da,
        &cli.decoy_prefix,
    );

    let ext = cli.spectrum
        .extension()
        .and_then(|e| e.to_str())
        .map(|s| s.to_lowercase());
    let ms_level_u32 = cli.ms_level as u32;
    let bench_mode = cli.max_spectra > 0;
    let bench_cap = if bench_mode { cli.max_spectra } else { usize::MAX };

    let mut all_spectra: Vec<Spectrum> = Vec::new();
    let mut all_queues: Vec<TopNQueue> = Vec::new();
    let mut chunk: Vec<Spectrum> = Vec::with_capacity(CHUNK_SIZE);
    let mut error_count = 0usize;
    let mut first_errors: Vec<String> = Vec::with_capacity(3);

    let flush_chunk = |chunk: &mut Vec<Spectrum>,
                           all_spectra: &mut Vec<Spectrum>,
                           all_queues: &mut Vec<TopNQueue>| {
        if chunk.is_empty() {
            return;
        }
        let offset = all_spectra.len();
        let queues = prepared.run_chunk(chunk, offset);
        all_queues.extend(queues);
        for mut spec in chunk.drain(..) {
            spec.peaks = Vec::new();
            all_spectra.push(spec);
        }
    };

    let t_search_start = std::time::Instant::now();

    match ext.as_deref() {
        Some("mzml") => {
            let f = File::open(&cli.spectrum)?;
            let reader = MzMLReader::new(BufReader::new(f))
                .with_ms_level_range(ms_level_u32, ms_level_u32);
            for result in reader {
                if all_spectra.len() + chunk.len() >= bench_cap {
                    break;
                }
                match result {
                    Ok(s) => {
                        chunk.push(s);
                        if chunk.len() >= CHUNK_SIZE {
                            flush_chunk(&mut chunk, &mut all_spectra, &mut all_queues);
                        }
                    }
                    Err(e) => {
                        error_count += 1;
                        if error_count <= 3 {
                            eprintln!("WARN: mzML parse: {e}");
                        }
                    }
                }
            }
            if bench_mode && all_spectra.len() + chunk.len() > bench_cap {
                let keep = bench_cap.saturating_sub(all_spectra.len());
                chunk.truncate(keep);
            }
            flush_chunk(&mut chunk, &mut all_spectra, &mut all_queues);
            if error_count > 0 {
                eprintln!(
                    "WARN: {} mzML spectra failed to parse",
                    error_count
                );
            }
            eprintln!(
                "MS-level filter: {} (only MS{} spectra entered the search)",
                cli.ms_level, cli.ms_level
            );
        }
        _ => {
            // MGF (default / backwards-compatible). MGF files do not encode
            // MS level — they are treated as MS2 by convention. Warn if the
            // user requested a non-default level so the mismatch is visible.
            if cli.ms_level != 2 {
                eprintln!(
                    "WARN: --ms-level={} requested for an MGF input; MGF files \
                     do not record MS level (treated as MS2). The flag has \
                     no effect on this input.",
                    cli.ms_level
                );
            }
            let f = File::open(&cli.spectrum)?;
            for result in MgfReader::new(BufReader::new(f)) {
                if all_spectra.len() + chunk.len() >= bench_cap {
                    break;
                }
                match result {
                    Ok(s) => {
                        chunk.push(s);
                        if chunk.len() >= CHUNK_SIZE {
                            flush_chunk(&mut chunk, &mut all_spectra, &mut all_queues);
                        }
                    }
                    Err(e) => {
                        error_count += 1;
                        if first_errors.len() < 3 {
                            first_errors.push(format!("{e}"));
                        }
                    }
                }
            }
            if bench_mode && all_spectra.len() + chunk.len() > bench_cap {
                let keep = bench_cap.saturating_sub(all_spectra.len());
                chunk.truncate(keep);
            }
            flush_chunk(&mut chunk, &mut all_spectra, &mut all_queues);
            if error_count > 0 {
                eprintln!(
                    "WARN: {} MGF spectra failed to parse (first {} errors):",
                    error_count,
                    first_errors.len()
                );
                for e in &first_errors {
                    eprintln!("  - {e}");
                }
            }
        }
    }

    if all_spectra.is_empty() {
        return Err(format!(
            "no spectra parsed from {}",
            cli.spectrum.display()
        )
        .into());
    }

    let search_elapsed = t_search_start.elapsed();
    eprintln!(
        "Loaded+scored {} spectra from {} in chunks of {} [PHASE stream_search: {:.2}s]",
        all_spectra.len(),
        cli.spectrum.display(),
        CHUNK_SIZE,
        t_phase.elapsed().as_secs_f64()
    );
    if bench_mode {
        eprintln!("Bench mode: capped at {} spectra", cli.max_spectra);
    }

    // Downstream code uses these names.
    let spectra = all_spectra;
    let queues = all_queues;

    let non_empty = queues.iter().filter(|q| !q.is_empty()).count();
    eprintln!(
        "Search complete: {non_empty} / {} spectra have PSMs (match_spectra wall: {:.2}s)",
        spectra.len(),
        search_elapsed.as_secs_f64()
    );

    // ── 8. Write PIN ─────────────────────────────────────────────────────────
    // Bench mode still writes PIN (so we can diff against the reference
    // fixture) but skips TSV.
    let t_phase = std::time::Instant::now();
    output::write_pin(&cli.output_pin, &spectra, &queues, &params, &idx, &cli.decoy_prefix)?;
    eprintln!(
        "Wrote PIN: {} [PHASE pin_write: {:.2}s] [PHASE TOTAL: {:.2}s]",
        cli.output_pin.display(),
        t_phase.elapsed().as_secs_f64(),
        t_total.elapsed().as_secs_f64()
    );

    if bench_mode {
        eprintln!("Bench mode: skipping TSV write.");
        return Ok(());
    }

    // ── 9. Write TSV (optional) ───────────────────────────────────────────────
    if let Some(ref tsv_path) = cli.output_tsv {
        let spec_file_name = cli
            .spectrum
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_else(|| cli.spectrum.display().to_string());
        output::write_tsv(tsv_path, &spectra, &queues, &params, &idx, &spec_file_name, true)?;
        eprintln!("Wrote TSV: {}", tsv_path.display());
    }

    Ok(())
}

/// Translate `(--fragmentation, --instrument, --protocol)` into a bundled
/// `.param` filename and resolve it under
/// `src/main/resources/ionstat/` relative to the cargo manifest dir.
///
/// CLI indices match Java's:
/// - fragmentation: 0=Auto/CID, 1=CID, 2=ETD, 3=HCD, 4=UVPD
/// - instrument:    0=LowRes,   1=HighRes, 2=TOF, 3=QExactive
/// - protocol:      0=Automatic,1=Phosphorylation, 2=iTRAQ,
///                  3=iTRAQPhospho, 4=TMT, 5=Standard
///
/// When all three are `None`, the historical default
/// `HCD_QExactive_Tryp.param` is returned (preserving existing tests'
/// behaviour). Only Tryp is supported as the enzyme component for now;
/// other enzymes require the user to pass `--param-file` directly.
///
/// Returns an error if the resolved filename does not exist on disk.
fn resolve_bundled_param(
    fragmentation: Option<u8>,
    instrument:    Option<u8>,
    protocol:      Option<u8>,
) -> Result<PathBuf, String> {
    // Default file when no flags are given — preserves the previous
    // hard-coded behaviour.
    if fragmentation.is_none() && instrument.is_none() && protocol.is_none() {
        return canonicalize_bundled("HCD_QExactive_Tryp.param");
    }

    let frag = match fragmentation.unwrap_or(0) {
        // 0 is "Auto/CID" in Java's `-m` semantics. We don't yet implement
        // activation-method inference per spectrum, so we map 0 → CID for
        // the .param-file picker.
        0 | 1 => "CID",
        2     => "ETD",
        3     => "HCD",
        4     => "UVPD",
        n     => return Err(format!(
            "invalid --fragmentation {n}: valid range is 0..=4 \
             (0=Auto/CID, 1=CID, 2=ETD, 3=HCD, 4=UVPD)"
        )),
    };
    let inst = match instrument.unwrap_or(0) {
        0 => "LowRes",
        1 => "HighRes",
        2 => "TOF",
        3 => "QExactive",
        n => return Err(format!(
            "invalid --instrument {n}: valid range is 0..=3 \
             (0=LowRes, 1=HighRes, 2=TOF, 3=QExactive)"
        )),
    };
    let prot_suffix: &str = match protocol.unwrap_or(0) {
        // Automatic/Standard: no suffix.
        0 | 5 => "",
        1     => "_Phosphorylation",
        2     => "_iTRAQ",
        3     => "_iTRAQPhospho",
        4     => "_TMT",
        n     => return Err(format!(
            "invalid --protocol {n}: valid range is 0..=5 \
             (0=Automatic, 1=Phosphorylation, 2=iTRAQ, \
              3=iTRAQPhospho, 4=TMT, 5=Standard)"
        )),
    };

    let filename = format!("{frag}_{inst}_Tryp{prot_suffix}.param");
    canonicalize_bundled(&filename)
}

/// Resolve a bundled `.param` filename under
/// `src/main/resources/ionstat/` relative to the crate's cargo manifest
/// dir (set at compile time). Returns a helpful error if the file does
/// not exist.
fn canonicalize_bundled(filename: &str) -> Result<PathBuf, String> {
    let candidate = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .join("src/main/resources/ionstat")
        .join(filename);
    candidate.canonicalize().map_err(|e| format!(
        "bundled param file not found at `{}`: {e}\n\
         Hint: not every (fragmentation, instrument, protocol) combination \
         has a bundled .param file. Supply --param-file <PATH> to specify \
         the scoring model explicitly, or list available files under \
         `src/main/resources/ionstat/`.",
        candidate.display()
    ))
}

#[cfg(test)]
mod param_resolver_tests {
    use super::*;

    #[test]
    fn default_resolves_to_hcd_qexactive_tryp() {
        // No flags → existing default.
        let p = resolve_bundled_param(None, None, None).unwrap();
        let s = p.to_string_lossy();
        assert!(
            s.ends_with("HCD_QExactive_Tryp.param"),
            "expected HCD_QExactive_Tryp.param, got {s}"
        );
    }

    #[test]
    fn hcd_qexactive_tmt_combo_resolves() {
        // (HCD, QExactive, TMT) → bundled HCD_QExactive_Tryp_TMT.param.
        let p = resolve_bundled_param(Some(3), Some(3), Some(4)).unwrap();
        let s = p.to_string_lossy();
        assert!(
            s.ends_with("HCD_QExactive_Tryp_TMT.param"),
            "expected HCD_QExactive_Tryp_TMT.param, got {s}"
        );
    }

    #[test]
    fn cid_lowres_tryp_resolves() {
        // (CID, LowRes, Standard) → CID_LowRes_Tryp.param.
        let p = resolve_bundled_param(Some(1), Some(0), Some(5)).unwrap();
        let s = p.to_string_lossy();
        assert!(
            s.ends_with("CID_LowRes_Tryp.param"),
            "expected CID_LowRes_Tryp.param, got {s}"
        );
    }

    #[test]
    fn missing_combo_errors_with_helpful_hint() {
        // (CID, HighRes, TMT) — not in the bundle. Must surface a
        // helpful "supply --param-file" hint.
        let err = resolve_bundled_param(Some(1), Some(1), Some(4)).unwrap_err();
        assert!(
            err.contains("supply --param-file") || err.contains("--param-file"),
            "expected hint about --param-file, got: {err}"
        );
    }

    #[test]
    fn rejects_out_of_range_fragmentation() {
        let err = resolve_bundled_param(Some(99), None, None).unwrap_err();
        assert!(err.contains("--fragmentation"));
    }

    #[test]
    fn rejects_out_of_range_instrument() {
        let err = resolve_bundled_param(None, Some(99), None).unwrap_err();
        assert!(err.contains("--instrument"));
    }

    #[test]
    fn rejects_out_of_range_protocol() {
        let err = resolve_bundled_param(None, None, Some(99)).unwrap_err();
        assert!(err.contains("--protocol"));
    }
}
