//! msgf-rust: end-to-end MS-GF+ search.
//!
//! Loads an MGF spectrum file and a FASTA target database, runs a tryptic
//! database search with default MS-GF+ parameters, and writes output in
//! Percolator `.pin` format (and optionally `.tsv` format).

use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;
use std::process::ExitCode;

use clap::Parser;
use model::{AminoAcidSetBuilder, ModLocation, Modification, PrecursorTolerance, ResidueSpec, Tolerance};
use scoring_crate::{Param, RankScorer};
use search::{match_spectra, SearchIndex, SearchParams};
use input::{FastaReader, MgfReader};

#[derive(Parser, Debug)]
#[command(
    name = "msgf-rust",
    about = "MS-GF+ Rust port: database search of MGF spectra against FASTA"
)]
struct Cli {
    /// Input MGF spectrum file.
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

    /// Minimum isotope error offset to try (default -1, matching Java `-ti -1,2`).
    #[arg(long, default_value = "-1")]
    isotope_error_min: i8,

    /// Maximum isotope error offset to try (default 2, matching Java `-ti -1,2`).
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

    /// Path to the .param scoring model file.
    ///
    /// If not supplied, the bundled `HCD_QExactive_Tryp.param` file from the
    /// MS-GF+ source tree is used (resolved relative to the Cargo manifest
    /// directory at compile time). When running the binary outside the source
    /// tree this path may not exist; supply --param-file explicitly in that
    /// case.
    #[arg(long)]
    param_file: Option<PathBuf>,
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
    // ── 1. Load FASTA target database ────────────────────────────────────────
    let target_db =
        FastaReader::load_all(BufReader::new(File::open(&cli.database)?))?;
    eprintln!(
        "Loaded {} target proteins from {}",
        target_db.proteins.len(),
        cli.database.display()
    );

    // ── 2. Build SearchIndex (target + reversed decoys) ───────────────────────
    let idx = SearchIndex::from_target_db(&target_db, &cli.decoy_prefix);

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

    let param = Param::load_from_file(&param_path)
        .map_err(|e| format!("loading param file {}: {e}", param_path.display()))?;
    let scorer = RankScorer::new(&param);

    // ── 5. Build SearchParams ─────────────────────────────────────────────────
    let mut params = SearchParams::default_tryptic(aa);
    params.precursor_tolerance =
        PrecursorTolerance::symmetric(Tolerance::Ppm(cli.precursor_tol_ppm));
    params.charge_range = cli.charge_min..=cli.charge_max;
    params.isotope_error_range = cli.isotope_error_min..=cli.isotope_error_max;
    params.top_n_psms_per_spectrum = cli.top_n;

    // ── 6. Load MGF spectra ───────────────────────────────────────────────────
    let mgf_file = File::open(&cli.spectrum)?;
    let mut spectra: Vec<_> = Vec::new();
    let mut error_count = 0usize;
    let mut first_errors: Vec<String> = Vec::with_capacity(3);
    for result in MgfReader::new(BufReader::new(mgf_file)) {
        match result {
            Ok(s) => spectra.push(s),
            Err(e) => {
                error_count += 1;
                if first_errors.len() < 3 {
                    first_errors.push(format!("{e}"));
                }
            }
        }
    }
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

    if spectra.is_empty() {
        return Err(format!(
            "no spectra parsed from {}",
            cli.spectrum.display()
        )
        .into());
    }
    eprintln!(
        "Loaded {} spectra from {}",
        spectra.len(),
        cli.spectrum.display()
    );

    // ── 7. Run match_spectra ──────────────────────────────────────────────────
    // Fragment tolerance of 0.5 Da matches the gf_bsa_parity integration test
    // and the Java MS-GF+ default for HCD data.
    let fragment_tol_da = 0.5_f64;
    let queues = match_spectra(
        &spectra,
        &idx,
        &params,
        &scorer,
        fragment_tol_da,
        &cli.decoy_prefix,
    );

    let non_empty = queues.iter().filter(|q| !q.is_empty()).count();
    eprintln!("Search complete: {non_empty} / {} spectra have PSMs", spectra.len());

    // ── 8. Write PIN ─────────────────────────────────────────────────────────
    output::write_pin(&cli.output_pin, &spectra, &queues, &params, &idx, &cli.decoy_prefix)?;
    eprintln!("Wrote PIN: {}", cli.output_pin.display());

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
