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
    /// If not supplied, the bundled `HCD_QExactive_Tryp.param` file from the
    /// MS-GF+ source tree is used (resolved relative to the Cargo manifest
    /// directory at compile time). When running the binary outside the source
    /// tree this path may not exist; supply --param-file explicitly in that
    /// case.
    #[arg(long)]
    param_file: Option<PathBuf>,

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

    // ── 3. Build AminoAcidSet (default mods: CAM fixed, Oxidation M variable) ─
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
    let aa = AminoAcidSetBuilder::new_standard()
        .add_fixed_mod(cam)
        .add_variable_mod(ox)
        .build()?;

    // ── 4. Load Param scoring model ───────────────────────────────────────────
    let param_path = match cli.param_file {
        Some(p) => p,
        None => {
            // Resolve bundled param relative to the source tree.  This works
            // when the binary is run from the workspace (e.g., `cargo run`
            // or `cargo test --release`).  The path is embedded at compile
            // time, so it is stable across platforms as long as the directory
            // structure is unchanged.
            let candidate = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join("../../..")
                .join("src/main/resources/ionstat/HCD_QExactive_Tryp.param");
            match candidate.canonicalize() {
                Ok(p) => p,
                Err(e) => {
                    return Err(format!(
                        "bundled param file not found at `{}`: {e}\n\
                         Hint: supply --param-file <PATH> to specify the \
                         scoring model explicitly.",
                        candidate.display()
                    )
                    .into());
                }
            }
        }
    };

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
