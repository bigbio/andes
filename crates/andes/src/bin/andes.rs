//! andes: end-to-end peptide-spectrum database search.
//!
//! Loads an MGF or mzML spectrum file and a FASTA target database, runs a
//! tryptic database search and writes output
//! in Percolator `.pin` format (and optionally `.tsv` format).
//!
//! Format dispatch by `--spectrum` extension: `.mzML`/`.mzml` → `MzMLReader`;
//! `.d` → `TimsTofReader` (native Bruker timsTOF, only under `--features
//! timstof`); otherwise `MgfReader` (default reader).

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::mpsc::{sync_channel, SyncSender};
use std::sync::Arc;
use std::thread;

#[path = "../rescore.rs"]
mod rescore;
#[path = "../glyco_seeds.rs"]
mod glyco_seeds;
#[path = "../config.rs"]
mod config;

use clap::{Args, Parser, Subcommand, ValueEnum};
use model::{
    activation::ActivationMethod, AminoAcidSetBuilder, InstrumentType, ModLocation, Modification,
    PrecursorTolerance, ResidueSpec, Spectrum, Tolerance,
};
use model_train::{
    ModelStore,
    accumulate::{merge, StatsAccumulator},
    counts::CountStats,
    estimate::{Estimator, EstimatorConfig},
    gate::evaluate_candidate,
    geometry::{corpus_charge_masses, derive_geometry, GeometryConfig},
    labeled::bootstrap_labels,
    select::{select, select_nearest, SelectionKey},
    protocol_to_experiment_class as store_protocol_to_experiment_class,
    store::{
        SourceLedger,
        update_add, update_remove, update_reweight, update_decay, commit_update,
        write_all_models_with_sources_pub, write_all_models_with_sources_and_gbdt_pub,
    },
};
use scoring_crate::{Param, RankScorer};
use search::{
    apply_shift_for_mode, apply_tightened_precursor_tolerance, build_spec_keys,
    learn_calibration_stats, CalibrationStats,
    PreparedSearch, PrecursorCalMode, SearchIndex, SearchParams, SpecKey, TopNQueue,
};
use search::candidate_index::index_cache_path;
use search::precursor_cal::{constants as cal_constants, sample_every_nth};
use input::{detect_instrument_type, FastaReader, MgfReader, Ms1Link, MzMLReader};

// Type alias to reduce clippy type_complexity warnings in the train path.
type ModelEntryOwned = (String, Param, Vec<(SourceLedger, CountStats)>);

/// Fragmentation method. `Auto` detects from the mzML's activation block and
/// falls back to the bundled `hcd_qexactive_tryp` model when nothing is detected.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum Fragmentation {
    #[clap(name = "auto")] Auto,
    #[clap(name = "CID")]  Cid,
    #[clap(name = "ETD")]  Etd,
    #[clap(name = "HCD")]  Hcd,
    #[clap(name = "UVPD")] Uvpd,
}

/// Search protocol: sample labeling or enrichment strategy applied during the experiment.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum Protocol {
    #[clap(name = "auto")]          Auto,
    #[clap(name = "phospho")]       Phospho,
    #[clap(name = "iTRAQ")]         Itraq,
    #[clap(name = "iTRAQ-phospho")] ItraqPhospho,
    #[clap(name = "TMT")]           Tmt,
    #[clap(name = "standard")]      Standard,
}

/// Enzymatic-cleavage enforcement at peptide span boundaries:
/// 2=fully, 1=semi, 0=non-specific.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum EnzymeSpecificity {
    #[clap(name = "non-specific")] NonSpecific,
    #[clap(name = "semi")]         Semi,
    #[clap(name = "fully")]        Fully,
}

/// Primary ranking mode: inherited RawScore (`rank`) or fused strong score (`strong`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum)]
enum ScoreFlag {
    /// Pick by the resolved model's instrument: `strong` for high-res, `rank` for low-res.
    #[default]
    Auto,
    Rank,
    Strong,
}

/// Candidate-resolution backing: in-RAM (`ram`, default) or out-of-core mmap
/// base-peptide index with lazy mod enumeration (`mmap`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum)]
enum CandidateIndexFlag {
    /// Pick automatically: estimate the in-RAM candidate index size against
    /// available memory and use out-of-core mmap only if it would not fit
    /// (default). Errs toward RAM (byte-identical to prior releases) and only
    /// falls back to mmap when RAM would risk an OOM.
    #[default]
    Auto,
    /// Force the in-RAM candidate index (advanced; may OOM on very large mod
    /// searches — that is what `auto` protects against).
    Ram,
    /// Force the out-of-core mmap'd base-peptide index (advanced; lower peak RAM,
    /// result-equivalent but not byte-identical to `ram`).
    Mmap,
}

/// Available system memory in bytes, from Linux `/proc/meminfo` (`MemAvailable`).
/// Returns `None` when it cannot be determined (non-Linux, restricted sandbox);
/// callers should then keep the conservative RAM default.
fn available_memory_bytes() -> Option<u64> {
    let meminfo = std::fs::read_to_string("/proc/meminfo").ok()?;
    for line in meminfo.lines() {
        if let Some(rest) = line.strip_prefix("MemAvailable:") {
            // e.g. "MemAvailable:   12345678 kB"
            let kb: u64 = rest.split_whitespace().next()?.parse().ok()?;
            return Some(kb.saturating_mul(1024));
        }
    }
    None
}

/// Emit a one-line search-progress update after each scored chunk: cumulative
/// spectra scored, throughput, and elapsed time. The streaming search has no
/// upfront total for mzML/MGF, so this reports a running count + rate rather than
/// a percentage (use it to see the search is live and how fast it is going).
fn report_search_progress(scored: usize, start: std::time::Instant) {
    let secs = start.elapsed().as_secs_f64();
    let rate = if secs > 0.0 { scored as f64 / secs } else { 0.0 };
    eprintln!("[search] {scored} spectra scored (~{rate:.0}/s, {secs:.0}s elapsed)");
}

/// Search arguments (shared by the default search path and exposed as a
/// flat arg group so that `andes --spectrum X --database Y --output-pin Z`
/// keeps working unchanged).
///
/// Note: `spectrum`, `database`, and `output_pin` are declared `Option<PathBuf>`
/// at the clap level so that they are not required when a subcommand (e.g.
/// `train`) is given.  When no subcommand is present, `run()` validates them
/// manually and returns an early error if they are missing.
#[derive(Args, Debug)]
struct SearchArgs {
    /// YAML run-configuration file. Any parameter can be set here (grouped by
    /// experiment: io/search/scoring/decoys/chimeric/refine/rescoring/glyco; see
    /// DOCS §1b). An explicit CLI flag always overrides the config value.
    #[arg(long = "config", value_name = "FILE")]
    config: Option<PathBuf>,

    /// Input spectrum file(s). Repeat `--spectrum` for multiple inputs (one PIN).
    /// Format is auto-detected per file by extension.
    #[arg(long)]
    spectrum: Vec<PathBuf>,

    /// Input FASTA database (target sequences only; decoys are generated automatically).
    #[arg(long)]
    database: Option<PathBuf>,

    /// Output Percolator PIN file path.
    #[arg(long)]
    output_pin: Option<PathBuf>,

    /// Output TSV file path (optional).
    #[arg(long)]
    output_tsv: Option<PathBuf>,

    /// Output QPX `.idparquet/` bundle directory (optional; OpenMS-compatible).
    /// Writes `psms.parquet` + `proteins.parquet` + `search_params.parquet`.
    #[arg(long)]
    output_parquet: Option<PathBuf>,

    /// Decoy prefix used when generating reversed decoy sequences.
    #[arg(long, default_value = "XXX_")]
    decoy_prefix: String,

    /// Decoy-accession SUFFIX used to RECOGNIZE pre-built decoys in the input
    /// FASTA (e.g. `rev` for quantms/OpenMS `<orig>_rev` decoys). When set, a
    /// protein is a decoy iff its accession starts with `<decoy-prefix>_` OR ends
    /// with this suffix. Typically paired with `--decoy-strategy none` so andes
    /// consumes an externally-built target+decoy database instead of generating
    /// its own decoys (which would double-decoy and bias FDR).
    #[arg(long = "decoy-suffix")]
    decoy_suffix: Option<String>,

    /// How to generate decoys: `reverse` (default; reverse each sequence),
    /// `shuffle` (seeded reproducible shuffle), `sequon-reverse` (reverse but
    /// restore each N-X-S/T sequon at its mirrored position — RECOMMENDED with
    /// `--glyco`: plain reversal maps N-X-S/T to S/T-X-N, so reversed decoys reach
    /// the glyco sequon gate at a lower rate than targets and the resulting
    /// q-values are anti-conservative), or `none` (no decoys — for a FASTA that
    /// already contains decoys, or external FDR). `none` with a target-only FASTA
    /// leaves the search without decoys (FDR can't be estimated) and warns.
    #[arg(long = "decoy-strategy", default_value = "reverse")]
    decoy_strategy: String,

    /// Seed for `--decoy-strategy shuffle` (reproducible decoys). Ignored by
    /// reverse/none.
    #[arg(long = "decoy-seed", hide = true, default_value_t = search::decoy::DEFAULT_DECOY_SEED)]
    decoy_seed: u64,

    /// Isotope-error offset range to try, as `MIN..MAX` (also accepts `MIN-MAX`).
    /// Negative offsets allowed. Unset defaults to `-1..2`, or `0..2` under `--glyco`
    /// (see the resolution site). Left as an `Option` so an EXPLICIT `-1..2` is
    /// distinguishable from the default and is never silently overridden.
    #[arg(long = "isotope-error", hide = true, value_parser = parse_isotope_error_range)]
    isotope_error: Option<(i8, i8)>,

    /// Precursor-mass calibration: `off`, `auto`, or `on`. `auto`/`on` learn a
    /// systematic ppm shift from confident PSMs in a pre-pass and tighten the
    /// precursor tolerance for the main search; `auto` skips the correction when
    /// the sample is too small to be reliable.
    #[arg(long = "precursor-cal", default_value = "auto", value_parser = parse_precursor_cal)]
    precursor_cal: PrecursorCalMode,

    /// Precursor mass tolerance as `VALUE+unit`. Accepts ppm (e.g. `20ppm`,
    /// high-res) or Da (e.g. `0.02da`/`0.02Da`, low-res precursor selection).
    /// Default `20ppm`.
    #[arg(long = "precursor-tol", default_value = "20ppm", value_parser = parse_precursor_tol)]
    precursor_tol: Tolerance,

    /// Precursor charge range to try when not specified in the spectrum, as
    /// `MIN..MAX` (also accepts `MIN-MAX`). Default `2..5`.
    #[arg(long = "charge", hide = true, default_value = "2..5", value_parser = parse_charge_range)]
    charge: (u8, u8),

    /// Maximum number of PSMs to retain per spectrum.
    #[arg(long, hide = true, default_value = "10")]
    top_n: u32,

    /// Number of Tolerable Termini (enzymatic-cleavage enforcement at span
    /// boundaries). `fully`: both termini must be cleavage sites (strict).
    /// `semi`: at least one terminus must be a cleavage site. `non-specific`:
    /// neither terminus needs to be a cleavage site.
    #[arg(long = "enzyme-specificity", alias = "ntt",
          hide = true, default_value = "fully", value_parser = parse_enzyme_specificity)]
    enzyme_specificity: EnzymeSpecificity,

    /// Proteolytic enzyme for in-silico digestion. Named values: trypsin
    /// (default), chymotrypsin, lysc, aspn, gluc, lysn, argc, alphalp,
    /// nocleavage, nonspecific. A wrong enzyme yields ~no PSMs (fails loud,
    /// not silent).
    ///
    /// Multi-protease digest: pass a `,`- or `+`-separated list (e.g.
    /// `--enzyme gluc,trypsin`) to accept peptides cleaved by ANY listed
    /// protease (the union of their cleavage rules). The FIRST entry is the
    /// primary — it drives model selection and the cleavage-credit feature; if
    /// no model matches that enzyme for the data, andes WARNs and you should
    /// pass `--model`. Combining a universal protease (nonspecific/alphalp)
    /// with specific ones makes the whole digest non-specific (warned).
    #[arg(long, default_value = "trypsin")]
    enzyme: String,

    /// Maximum number of missed cleavages per peptide.
    #[arg(long, hide = true, default_value = "1")]
    max_missed_cleavages: u32,

    /// Minimum number of peaks an MS2 spectrum must have to be scored; spectra
    /// with fewer peaks are skipped.
    #[arg(long, hide = true, default_value = "10")]
    min_peaks: u32,

    /// Minimum peptide length, in residues.
    #[arg(long, hide = true, default_value = "6")]
    min_length: u32,

    /// Maximum peptide length, in residues. (50 matches the reference engine/a comparison engine defaults;
    /// 40 dropped long tryptic peptides.)
    #[arg(long, hide = true, default_value = "50")]
    max_length: u32,

    /// Maximum number of variable modifications per peptide. A `NumMods=N` line
    /// in a --mods file overrides this.
    #[arg(long = "max-mods", hide = true, default_value = "3")]
    max_mods: u32,

    /// Path to a mods.txt file describing fixed and variable modifications.
    /// Format: each non-comment line is
    /// `<mass>,<aa>,<fix|opt>,<location>,<name>`, where:
    ///   - `<mass>` is a numeric monoisotopic mass delta (Da). Composition
    ///     strings (e.g. `C2H3N1O1`) are **not** yet supported.
    ///   - `<aa>` is a single uppercase letter or `*` (wildcard).
    ///   - `<location>` is one of `any|N-term|C-term|Prot-N-term|Prot-C-term`.
    ///
    /// A single `NumMods=N` line sets the max variable mods per peptide.
    /// Inline `#`-comments are stripped. Blank lines and full-line `#`-comments
    /// are ignored. When omitted, the binary uses its built-in defaults
    /// (Carbamidomethyl-C fixed, Oxidation-M + protein-N-term-Acetyl variable).
    #[arg(long = "mods", value_name = "MODFILE")]
    mods: Option<PathBuf>,

    /// Fragmentation/activation method for MGF input only. mzML/.raw/.d
    /// auto-detect this. Named values: auto, CID, ETD, HCD, UVPD.
    #[arg(long, hide = true, default_value = "auto", value_parser = parse_fragmentation)]
    fragmentation: Fragmentation,

    /// Search protocol. Named values: auto, phospho, iTRAQ, iTRAQ-phospho, TMT, standard.
    #[arg(long, hide = true, default_value = "auto", value_parser = parse_protocol)]
    protocol: Protocol,

    /// Fragment-matching tolerance in ppm for **MGF input only** (high-resolution
    /// MS/MS). Has no effect on mzML/.raw/.d (analyzer auto-detected). Mutually
    /// exclusive with `--fragment-tol-da`.
    #[arg(long = "fragment-tol-ppm", hide = true, conflicts_with = "fragment_tol_da", value_parser = parse_positive_tol)]
    fragment_tol_ppm: Option<f64>,

    /// Fragment-matching tolerance in Da for **MGF input only** (low-resolution
    /// ion-trap MS/MS). Has no effect on mzML/.raw/.d. Mutually exclusive with
    /// `--fragment-tol-ppm`.
    #[arg(long = "fragment-tol-da", hide = true, conflicts_with = "fragment_tol_ppm", value_parser = parse_positive_tol)]
    fragment_tol_da: Option<f64>,

    /// Number of worker threads for the search loop. Defaults to logical CPU count.
    #[arg(long, default_value_t = num_cpus::get())]
    threads: usize,

    /// Debug/benchmark cap: process only the first N spectra (0 = no cap).
    #[arg(long, hide = true, default_value = "0")]
    max_spectra: usize,

    /// MS level to search. Defaults to MS2 (identification); MS1 and any higher
    /// levels (e.g. TMT SPS-MS3 reporter-quant scans) are filtered out at load
    /// time so they never enter the search loop. Override only if you explicitly
    /// want a different level. Applies to mzML and Thermo `.raw`; MGF files do
    /// not encode MS level and are always treated as MS2. The chimeric cascade
    /// always searches MS2 (it pairs MS2 with its preceding MS1).
    #[arg(long, hide = true, default_value = "2")]
    ms_level: u8,

    /// Enable the two-pass chimeric cascade for co-isolated (co-fragmented)
    /// peptides. Pass 1 is the normal top-1 search; Pass 2 detects co-isolated
    /// precursors in each scan's MS1 isolation window and runs a targeted search
    /// for the second peptide on the residual spectrum, emitting it as an extra
    /// PSM. Requires mzML (MS1 scans); has no effect on MGF input.
    #[arg(long, default_value = "false")]
    chimeric: bool,

    /// Chimeric mode: max co-isolated SECONDARY peptides to search per scan (the
    /// chimeric-N lever). Default 4 = the measured Astral sweet spot (+1.4% PSMs
    /// vs N=2 at flat FDP; saturates by N=4). Set 2 for the original behavior.
    #[arg(long = "chimeric-max-coisolated", hide = true, default_value = "4")]
    chimeric_max_coisolated: usize,

    /// Chimeric mode: averagine-envelope KL gate for accepting a co-isolated MS1
    /// envelope (lower = stricter/cleaner; fewer spurious secondaries).
    #[arg(long = "chimeric-max-kl", hide = true, default_value = "0.3")]
    chimeric_max_kl: f32,

    /// Path to a Parquet model store to use instead of the bundled
    /// `resources/models.parquet`. When set, model selection reads from
    /// this store; when unset, the bundled store is used.
    #[arg(long = "model-store", hide = true)]
    model_store: Option<PathBuf>,

    /// Exact model ID to load from the model store (bundled or `--model-store`).
    /// When set, skips automatic selection (metadata detection / `--fragmentation`
    /// / `--protocol`) and loads this ID directly. Useful after `andes train`
    /// to search with the freshly-trained model.
    #[arg(long = "model", hide = true)]
    model_id_override: Option<String>,

    /// Evaluate only the first N trees of each GBDT ensemble (0 = all).
    ///
    /// Applies to BOTH shipped ensembles — fragment-intensity and rich-ion — which are
    /// of comparable size and are each walked per fragment per candidate.
    ///
    /// `Tree::eval` on these ensembles is the single hottest operation in the search.
    /// Truncating the ensemble trades prediction fidelity for speed: the GBDT is
    /// additive, so the first trees carry the bulk of the signal and later ones
    /// refine it. UNLIKE the per-candidate de-duplication (which is byte-identical),
    /// this CHANGES the predicted intensities and therefore the emitted PIN feature
    /// values, so it can move identifications. Leave at 0 unless you have measured
    /// the identification cost on your own data.
    #[arg(long = "gbdt-max-trees", default_value_t = 0usize)]
    gbdt_max_trees: usize,

    /// Path to a trained intensity model parquet (`andes train-intensity` output).
    /// Populates the additive `IntensitySignal` PIN column; ranking stays on RawScore
    /// until `--score strong` is enabled in a later phase. When unset, the column is 0.0.
    #[arg(long = "intensity-model", hide = true)]
    intensity_model: Option<PathBuf>,

    /// Ranking / PIN RawScore source: `auto` (default — `strong` for high-res
    /// instruments, `rank` for low-res), `rank`, or `strong` (fused intensity +
    /// competition score from S1–S3).
    #[arg(long = "score", default_value = "auto")]
    score: ScoreFlag,

    /// Candidate-index backing: `auto` (default — automatically use out-of-core
    /// mmap only when the in-RAM candidate index would not fit available memory;
    /// otherwise RAM, byte-identical to prior releases), or force `ram` / `mmap`
    /// (advanced overrides). `mmap` lowers peak RAM with lazy per-spectrum mod
    /// enumeration (result-equivalent PSMs, not byte-identical).
    #[arg(long = "candidate-index", hide = true, default_value = "auto")]
    candidate_index: CandidateIndexFlag,

    /// Glycopeptide search mode: enumerate hybrid backbone candidates (DB + de-novo
    /// Y-ladder), filter by N-X-S/T sequon, score bare backbones, and write a
    /// `.glyco.pin` file instead of the standard PIN. Default off.
    #[arg(long = "glyco", default_value_t = false)]
    glyco: bool,

    /// Maximum backbone candidates per spectrum in glyco mode (DB + de-novo
    /// combined, after union-dedup). Hidden advanced knob; default 50.
    /// Raised from 20: core-Y evidence ranking means the cap now cuts fewer
    /// true positives, so more headroom is inexpensive and safe.
    #[arg(long = "glyco-backbone-top-k", hide = true, default_value_t = 150usize)]
    glyco_backbone_top_k: usize,


    /// Cap the peaks the glyco GENERATION stage considers, keeping the most
    /// intense N. The backbone solver is superlinear in peak count, so an
    /// uncentroided profile scan or a very dense wide-window scan can take tens of
    /// seconds while a normal scan takes milliseconds — the run looks hung.
    /// Scoring always reads the full spectrum, so a generated candidate is never
    /// scored on truncated evidence. Default 0 = no cap; 300-500 is a reasonable
    /// value if you hit this. Changing it changes results.
    #[arg(long = "glyco-max-peaks", default_value_t = 0usize)]
    glyco_max_peaks: usize,

    /// Maximum c/z fragment charge to probe in `--glyco` ETD scoring. Unset derives it
    /// from whether the spectrum was deconvoluted, which is correct in almost all cases:
    /// after deconvolution multiply-charged fragments have already been moved to 1+.
    /// Set this only for data known to carry unresolved high-charge c/z ions.
    #[arg(long = "glyco-cz-max-charge")]
    glyco_cz_max_charge: Option<u8>,

    /// Weight the explained-c/z terms by observed peak intensity instead of treating a
    /// match as presence-only. Off by default: it measured -48 backbone-correct @1% on
    /// the benchmark, though presence-only scoring is a known weakness for large glycans.
    #[arg(long = "glyco-cz-intensity", default_value_t = false)]
    glyco_cz_intensity: bool,

    /// Maximum glycan-Y fragment charge. Default 3; raising it reaches 4+/5+ Y ions on
    /// highly-charged precursors at the cost of more chance matches.
    #[arg(long = "glyco-y-max-charge", default_value_t = 3u8)]
    glyco_y_max_charge: u8,

    /// Choose the glycosite by c/z evidence when a peptide carries more than one
    /// N-X-S/T sequon (~8% of tryptic N-glycopeptides). Off by default: the default
    /// positional convention is decoy-symmetric, and enabling this is gated on a
    /// decoy-controlled A/B that would surface any sequon-count asymmetry.
    #[arg(long = "glyco-cz-multisite", default_value_t = false)]
    glyco_cz_multisite: bool,

    /// Windowed peak filtering as `WINDOW_DA:PEAKS` (e.g. `100:20`). Unset uses the
    /// protocol default — on for isobaric-labelled data, off otherwise. A window of 0
    /// forces it off.
    #[arg(long = "peak-filter")]
    peak_filter: Option<String>,

    /// Clamp the precursor-offset lookup to the nearest available charge when the exact
    /// charge is missing from the model, rather than dropping the correction.
    #[arg(long = "precursor-offset-clamp", default_value_t = true, action = clap::ArgAction::Set)]
    precursor_offset_clamp: bool,

    /// Measure local peak density on the active (deconvoluted) peak list rather than the
    /// raw list.
    #[arg(long = "density-on-active-list", default_value_t = true, action = clap::ArgAction::Set)]
    density_on_active_list: bool,

    /// Serve high-resolution models at the 20 ppm window their rank tables were TRAINED
    /// with instead of the model's stored tolerance (0.5 Da for every bundled model, so
    /// ~50x wider at m/z 500).
    ///
    /// MEASURED AND NOT RECOMMENDED. It is a real train/serve mismatch and it does help
    /// glyco (+4.9% on the AI-ETD benchmark), but on ordinary peptide search it is
    /// catastrophic: Astral fell 36,719 -> 28,894 identifications at 1% FDR, a 21% loss.
    /// The wide window is evidently load-bearing for the rank tables as trained, so
    /// closing the mismatch requires retraining the models, not re-serving them. Kept as
    /// a flag only so the experiment is repeatable.
    #[arg(long = "tight-highres-scoring", default_value_t = false)]
    tight_highres_scoring: bool,

    /// Allow a Pass-2 co-isolated candidate to overlap the primary's matched peaks.
    /// Off by default: the residual spectrum has the primary's peaks removed, and
    /// permitting overlap lets the same evidence support two PSMs.
    #[arg(long = "chimeric-allow-overlap", default_value_t = false)]
    chimeric_allow_overlap: bool,

    /// How to label EThcD/ETciD spectra (electron transfer with a supplemental
    /// collisional term). `hcd` is the default and is what model routing expects, since
    /// no EThcD model exists; `etd` labels them ETD so the c/z scoring path engages.
    #[arg(long = "ethcd-activation", value_enum, default_value_t = EthcdActivationFlag::Hcd)]
    ethcd_activation: EthcdActivationFlag,

    /// Diagnostic: restrict `--glyco` scoring to the scan numbers in this file, one per
    /// line. Makes a `--debug-glyco` dump of a chosen set of scans affordable.
    #[arg(long = "glyco-scans")]
    glyco_scans: Option<PathBuf>,

    /// Diagnostic: log resident set size at each phase boundary.
    #[arg(long = "rss-probe", default_value_t = false)]
    rss_probe: bool,

    /// Glycan composition list for `--glyco`.
    ///
    /// `common` (~600, default) is the MEASURED-BEST list and what the benchmarks were
    /// run with. `reference-human` (~2,300) reaches high-antennary glycans `common`
    /// cannot name -- 100% of a curated 160-composition human reference vs `common`'s
    /// 68% -- but measured WORSE overall on human plasma with an entrapment database
    /// (228 glycoPSMs at 0.00% entrapment FDP, vs 365 at 0.55%). `full` (~4,034) is
    /// wider still and was measured to raise entrapment error 5.4x on a benchmark where
    /// it looked like a gain on yield alone.
    ///
    /// Bigger is not better here: a larger candidate space gives decoys more places to
    /// fit, which tightens Percolator's threshold and leaves real identifications
    /// behind. Prefer `common` unless you have measured otherwise on your own data.
    #[arg(long = "glyco-glycan-list", value_enum, default_value_t = GlycanListFlag::Common)]
    glyco_glycan_list: GlycanListFlag,

    /// Exclude NeuGc (N-glycolylneuraminic acid) glycans from the search list.
    ///
    /// Humans do not synthesise NeuGc — CMAH is inactivated in the human lineage, so
    /// NeuGc in a human sample is trace dietary only. Most other mammals (mouse
    /// included) DO make it, which is why the mouse-developed glyco benchmarks never
    /// surfaced this.
    ///
    /// NeuGc is also the ENTIRE source of isobaric ambiguity in the default list:
    /// `NeuGc - NeuAc = 15.994914` and `Hex - Fuc = 15.994915`, so a NeuGc composition
    /// is mass-degenerate with a NeuAc one. Measured on the default list: 600
    /// compositions over only 460 distinct masses, 140 masses (30%) carrying more than
    /// one composition, and 100% of those collisions involve NeuGc. Excluding it gives
    /// 360 compositions over 360 masses — zero collisions, by construction — and a 40%
    /// smaller list to search.
    ///
    /// Use for human samples. Leave off for mouse and other CMAH-competent species.
    #[arg(long = "glyco-no-neugc", default_value_t = false)]
    glyco_no_neugc: bool,

    /// Glycan biology to assume for the search space.
    ///
    /// `auto` (default) surveys the NeuGc/NeuAc oxonium ratio across the run, and treats
    /// the FASTA's `OX=` taxon ids as a VETO rather than a second vote: it narrows the
    /// list when the spectra show no NeuGc, unless the database is a CMAH-competent
    /// organism. A database with no `OX=` headers, or a mixed one, does not block
    /// narrowing. It prints both signals and the decision. `human` forces NeuGc out,
    /// `mammal` forces it in.
    ///
    /// `--glyco-no-neugc` is the explicit override and wins over this.
    #[arg(long = "glyco-taxon", value_enum, default_value_t = GlycoTaxonFlag::Auto)]
    glyco_taxon: GlycoTaxonFlag,


    /// Isotope-error range for `--glyco`. `default` uses 0..=2 — the -1 offset costs
    /// 0.29% of correct answers at a ~53:47 target:decoy ratio (pure FDR dilution),
    /// and dropping it measured +81 backbone-correct @1%. `negative` restores
    /// -1..=2; `wide` extends the upper bound to 5 for heavily-labelled precursors.
    #[arg(long = "glyco-isotope-error", value_enum, default_value_t = GlycoIsotopeFlag::Default)]
    glyco_isotope_error: GlycoIsotopeFlag,

    /// Fragment tolerance (ppm) for the glyco-specific matching: oxonium ions,
    /// the core-Y ladder, backbone mass search, and c/z. Default 20 ppm, which
    /// suits Orbitrap MS2. **Raise this for low-resolution (ion-trap) MS2** —
    /// at 20 ppm a 0.3-0.5 Da ion-trap peak never matches, so the oxonium gate
    /// never fires and glyco IDs collapse to near zero. This is separate from
    /// `--fragment-tol-ppm`, which the scoring model owns.
    #[arg(long = "glyco-tol-ppm", default_value_t = 20.0f64)]
    glyco_tol_ppm: f64,

    /// `gp` fused-selector ladder weight K (`rank + K·ladder + J·core_y + H·hyper`).
    /// Hidden tuning knob; default 10 (lowered from 50 in round-2 — K·ladder is
    /// per-backbone and non-discriminating between isobaric peptides; see
    /// GLYCO_GP_K_DEFAULT).
    #[arg(long = "glyco-gp-k", hide = true, default_value_t = andes_glyco::glyco_psm::GLYCO_GP_K_DEFAULT)]
    glyco_gp_k: f32,

    /// `gp` fused-selector core-Y hit-count weight J. Hidden tuning knob; default 5.
    #[arg(long = "glyco-gp-j", hide = true, default_value_t = andes_glyco::glyco_psm::GLYCO_GP_J_DEFAULT)]
    glyco_gp_j: f32,

    /// `gp` fused-selector hyperscore weight H (0 disables). Hidden tuning knob; default 1.
    #[arg(long = "glyco-gp-h", hide = true, default_value_t = andes_glyco::glyco_psm::GLYCO_GP_H_DEFAULT)]
    glyco_gp_h: f32,

    /// `gp` selector ETD c/z-hyperscore weight (added ONLY on ETD/AI-ETD spectra;
    /// inert on HCD). Hidden knob; default 15 (raised from 5 in round-2 — c/z is
    /// the only per-candidate discriminator on ETD). 0 disables ETD c/z selection.
    #[arg(long = "glyco-gp-cz", hide = true, default_value_t = andes_glyco::glyco_psm::GLYCO_GP_CZ_DEFAULT)]
    glyco_gp_cz: f32,

    /// `gp` selector weight on the COUNT of matched b/y ions. The collapse runs before
    /// feature extraction, so it cannot see the strong score — the engine's best
    /// discriminator. Measured over a benchmark's reference identifications, the terms it
    /// can see rank the correct candidate at median 15 (rank) and median 44 (ladder, the
    /// heaviest weight), while this count ranks it at median 1-2. It is free: the count
    /// falls out of the hyperscore the selector already computes per candidate.
    /// Default 0 reproduces the previous selector exactly.
    #[arg(long = "glyco-gp-m", hide = true, default_value_t = andes_glyco::glyco_psm::GLYCO_GP_M_DEFAULT)]
    glyco_gp_m: f32,

    /// Resolve isobaric glycan-composition collisions on Y-ladder evidence rather
    /// than sort order. Two compositions can be isobaric to ~1 uDa (Hex-Fuc and
    /// NeuGc-NeuAc both = 15.9949), and above 2000 Da more than half of the default
    /// list has such a twin; without this the survivor is chosen by `to_bits()`
    /// ordering, so the SAME glycan mass can be annotated with DIFFERENT compositions
    /// on different spectra. Measured on PXD030622 plasma: andes emitted 131
    /// composition strings over 53 distinct masses (~2.5 per mass) where Byonic was
    /// 1.0. Off by default: the original A/B scored peptide YIELD (-8), which cannot
    /// see this. Judge it on compositions-per-mass, not on ID count.
    #[arg(long = "glyco-isobar-rep", default_value_t = false)]
    glyco_isobar_rep: bool,

    /// Keep backbones by GLYCAN-Y evidence as well as by peptide b/y ("two-axis
    /// retention"), and enable glycan-Y-first candidate generation.
    ///
    /// Without this, backbone truncation retains on peptide b/y rank (axis 1), c/z
    /// (axis 4, ETD-ONLY) and transfer (axis 3, off by default). On an HCD-only run --
    /// the human plasma regime -- peptide b/y is therefore the ONLY surviving axis, and
    /// it is the weakest one for large glycopeptides: a backbone anchored by a strong
    /// core-Y ladder but with few b/y ions is truncated before the fused selector ever
    /// sees it. The glycan-Y evidence is already computed and is otherwise used only as
    /// a tiebreak.
    ///
    /// Default off (the validated baseline). Costs a second top-k retention pass.
    #[arg(long = "glyco-y-index", default_value_t = false)]
    glyco_y_index: bool,

    /// Compute the PIN feature vector against the GLYCAN-DECORATED backbone instead of
    /// the bare deglycosylated peptide.
    ///
    /// By default the ~40 feature columns are computed on the bare backbone, so every
    /// glycosite-spanning fragment sits at the wrong theoretical mass -- roughly half the
    /// b/y ladder of a glycopeptide. IntensitySignal, MatchedIonRatio,
    /// ExplainedIonCurrentRatio, LongestComplementaryLadder and strong_score therefore
    /// describe a molecule that was never in the tube, and those are the columns
    /// Percolator weights most heavily.
    ///
    /// MEASURED AT -41% AND LEFT OFF. On PXD030622 plasma with an E. coli entrapment
    /// database this took 365 glycoPSMs @0.55% FDP down to 215 @0.00%. The premise was
    /// wrong: under HCD a glycopeptide fragments at the GLYCOSIDIC bonds first, so b/y
    /// ions come from the backbone AFTER the glycan is lost -- the BARE backbone is the
    /// correct theoretical ladder, and decorating moves half the predicted ladder to
    /// masses with no peaks. Consistent with the -16 measured on the scoring peptide
    #[arg(long = "glyco-decorated-features", default_value_t = false)]
    glyco_decorated_features: bool,

    /// Require a matching sialic OXONIUM ion before a glycan composition may claim
    /// NeuAc or NeuGc, as a fraction of base-peak intensity. 0 disables the gate.
    ///
    /// NeuAc and NeuGc are indistinguishable by precursor mass when traded against
    /// Hex/Fuc -- Hex1NeuAc1 and Fuc1NeuGc1 are the SAME elemental formula -- but they
    /// are distinguishable in oxonium ions: NeuAc gives m/z 274.092/292.103, NeuGc gives
    /// 290.087/308.098. Gating on those is how pGlyco3 breaks the degeneracy, and it is
    /// the evidence-based alternative to excluding NeuGc by species
    /// (`--glyco-taxon` / `--glyco-no-neugc`), so it also works where NeuGc is real.
    ///
    /// Deliberately a threshold, not a presence test: Chalkley & Baker (MCP 2025) found
    /// ~70% of spectra carrying a NeuGc oxonium contained no NeuGc, from co-isolation, so
    /// a binary test admits almost everything.
    ///
    /// MEASURED on PXD030622 plasma with an E. coli entrapment database: it fixes
    /// CALIBRATION, not yield. 2% gives 267 glycoPSMs @0.00% entrapment FDP and 5% gives
    /// 241 @0.00%, against an ungated 268 @1.87% -- so it flips the verdict from
    /// OPTIMISTIC to CONSERVATIVE at no yield cost, but buys no identifications, and an
    /// FDP pinned at 0.00% means the threshold has tightened past the useful point.
    /// Species exclusion (`--glyco-no-neugc`) still wins on yield there: 365 @0.55%.
    /// If you tune this, go LOOSER (0.005-0.01), not stricter.
    ///
    /// Gates SIALIC only, never fucose -- PTM-Shepherd's published hit/miss ratios weight
    /// absence of a fucose oxonium 10x weaker than absence of a sialic one
    #[arg(long = "glyco-sialic-oxonium-min-frac", default_value_t = 0.0f32,
          value_parser = parse_unit_fraction_f32)]
    glyco_sialic_oxonium_min_frac: f32,

    /// Minimum trimannosyl-core Y ions required before `--glyco` reports a PSM for a
    /// scan. andes historically reported a best guess for every scan clearing the
    /// oxonium gate, so most reported rows had no glycan evidence at all. Every other
    /// engine requires this: pGlyco3 and O-Pair require 2 core Y ions, Glyco-Decipher 3
    /// with Y1 mandatory. Reads only spectral evidence, so it applies equally to target
    /// and decoy scans and cannot skew the target/decoy ratio. 0 disables.
    ///
    /// MEASURED TRADE-OFF, which is why the default is 0. On a pooled human plasma set
    /// (stepped-collision HCD) `2` took verified-correct identifications from 0 to 87 at
    /// 1% FDR with a measured 0.75% false-discovery proportion — the ungated run
    /// accepted 102 PSMs of which NONE were correct. On a mouse AI-ETD benchmark the
    /// same value cost 161 of 707 identifications. Set it for collision-dominant data;
    /// leave it off for electron-transfer data. Exempting ETD scans automatically was
    /// tried and made both regimes worse, because real files are mixed.
    #[arg(long = "glyco-min-core-y", default_value_t = 0u32)]
    glyco_min_core_y: u32,

    /// Minimum winner RawScore for a `--glyco` scan to emit a PIN row at all.
    /// Unset = emit a best guess for every gated scan (historical behaviour).
    ///
    /// Measured on plasma (2026-08-28): 90.5% of emitted rows sit on scans with no
    /// glycopeptide in them (median RawScore −2.5 vs +9.4 on real glyco scans);
    /// that stratum is what Percolator trains on. At 3, it removes 83% of those
    /// rows while keeping every measured agreement with an external engine.
    /// Label-blind: reads only the winner's spectral match quality.
    #[arg(long = "glyco-min-raw-score")]
    glyco_min_raw_score: Option<f32>,

    /// Run-ADAPTIVE emission floor: drop scans whose winner scores below this
    /// quantile of the run's own decoy winners (e.g. 0.95). Self-calibrating --
    /// unlike an absolute --glyco-min-raw-score, it transfers across datasets,
    /// instruments and models, because the decoy winners ARE the run's null.
    /// The derived threshold is printed and applied identically to target and
    /// decoy scans. Mutually exclusive with --glyco-min-raw-score.
    #[arg(long = "glyco-min-raw-score-quantile")]
    glyco_min_raw_score_quantile: Option<f64>,

    /// Diagnostic TSV of per-candidate split evidence with sampled shifted-ladder
    /// nulls (the LLR-calibration probe). Requires --debug-glyco; never affects
    /// the PIN.
    #[arg(long = "glyco-diag-splits", requires = "debug_glyco")]
    glyco_diag_splits: Option<std::path::PathBuf>,

    /// Minimum matched b/y sequence ions required before `--glyco` reports a PSM.
    /// MSFragger's equivalents are 4 matched fragments with at least 2 non-Y. 0 disables.
    #[arg(long = "glyco-min-matched-ions", default_value_t = 0u32)]
    glyco_min_matched_ions: u32,

    /// c/z truncation gate: keep the top-k backbones by glycosite-spanning c/z
    /// evidence (AXIS 4) so high-charge ETD glycopeptides supported mainly by c/z
    /// survive Phase-1 truncation. Default ON; ETD-only (inert on HCD/CID). Pass
    /// `--glyco-cz-gate false` to disable. `action = Set` so the bool takes an
    /// explicit value (a bare bool arg would be an un-disableable set-true flag).
    #[arg(long = "glyco-cz-gate", hide = true, default_value_t = true, action = clap::ArgAction::Set)]
    glyco_cz_gate: bool,

    /// Charge states indexed by the peptide-first fragment index (b/y at 1..=N,
    /// clamped 1..=3); targets high-charge glycopeptides. Hidden knob; default 2.
    #[arg(long = "glyco-pf-charge", hide = true, default_value_t = 2u8)]
    glyco_pf_charge: u8,

    /// Max peptide-first candidates per spectrum. Hidden knob; default 1024.
    #[arg(long = "glyco-max-pf", hide = true, default_value_t = 1024usize)]
    glyco_max_pf: usize,

    /// MEASURED AND NOT RECOMMENDED: dispatching ETD scans to `etd_highres_tryp`
    /// instead of the file's dominant HCD model LOSES identifications — mouse frac2
    /// 707 -> 692 at 1% FDR. The bundled ETD model is evidently a poorer fit for these
    /// spectra than the HCD model, so the long-standing "ETD scans are scored by an HCD
    /// model" behaviour is benign, not the bug it looked like.
    /// Diagnostic glyco mode: emit ALL candidate rows per scan (including de-novo
    /// mass-residual hits) and print transfer diagnostics. The resulting PIN is for
    /// inspection ONLY and must never be fed to an FDR tool. Hidden dev flag.
    #[arg(long = "debug-glyco", hide = true, default_value_t = false)]
    debug_glyco: bool,

    /// Emit paired glycan-axis decoy rows for 2D (peptide × glycan) FDR
    /// discrimination (experimental). Off by default. Hidden.
    #[arg(long = "glyco-decoy", hide = true, default_value_t = false)]
    glyco_decoy: bool,

    /// On ETD/AI-ETD spectra, generate candidate backbones from the paired HCD scan
    /// (same precursor) while scoring c/z on the ETD scan — targets high-charge
    /// glycopeptides (validated +153 backbone-correct @1%). DEFAULT ON; scans with no
    /// HCD partner (and multi-file runs) fall back to unpaired automatically. Disable
    /// with `--glyco-hcd-pair false`. `action = Set` so the bool takes an explicit value.
    #[arg(long = "glyco-hcd-pair", default_value_t = true, action = clap::ArgAction::Set)]
    glyco_hcd_pair: bool,

    /// BUG2 fix, EXPERIMENTAL: on ETD/AI-ETD spectra, score the rank/edge/
    /// hyperscore path (RawScore, EdgeScore, hyperscore, RankScoreFloat) against a
    /// peptide clone carrying the intact glycan on its glycosite instead of the
    /// bare backbone, so glycosite-spanning c/z fragments are computed at the real
    /// (glycan-carrying) mass. DEFAULT ON (round-6: validated +33 backbone-correct @1%,
    /// decoy-safe); inert on HCD/CID. Disable with `--glyco-etd-rank-glycan false`.
    #[arg(long = "glyco-etd-rank-glycan", default_value_t = true, action = clap::ArgAction::Set)]
    glyco_etd_rank_glycan: bool,

    /// BUG5, EXPERIMENTAL: per-spectrum-activation model dispatch. andes normally
    /// selects ONE scoring model for the whole file by majority vote over the
    /// first 64 mzML spectra (`detect_dominant_activation`), so on a mixed
    /// HCD/ETD (EThcD / AI-ETD) file every scan is scored with the dominant
    /// model regardless of its own activation. When this is on, the glyco driver
    /// ALSO loads the ETD-family model (same instrument/protocol/enzyme lookup,
    /// activation forced to ETD) and dispatches each spectrum's own rank/edge/
    /// hyperscore/RankScoreFloat scoring to whichever model matches ITS OWN
    /// activation method — an ETD/AI-ETD scan scores against the ETD model, an
    /// HCD scan (including an HCD partner read for `--glyco-hcd-pair`
    /// candidate generation) still scores against the file's dominant model.
    /// Off by default (byte-identical to the pre-dispatch single-model path).
    /// Falls back to the single-model behavior with a WARN if no ETD-family
    /// model exists in the store. Hidden.
    #[arg(long = "glyco-per-spectrum-model", hide = true, default_value_t = false)]
    glyco_per_spectrum_model: bool,

    /// Cross-spectrum transfer: q-value threshold for confident donor seeds
    /// (--glyco-transfer only). Hidden; default 0.05 (native GBDT q is conservative).
    #[arg(long = "glyco-transfer-seed-fdr", hide = true, default_value_t = 0.05f64)]
    glyco_transfer_seed_fdr: f64,

    /// Cross-spectrum transfer: RT co-elution window in seconds. Hidden; default 1800.
    #[arg(long = "glyco-rt-window", hide = true, default_value_t = 1800.0f32)]
    glyco_rt_window: f32,

    /// Cross-spectrum transfer: skip the RT co-elution gate (unsafe research opt-in,
    /// transfers across the whole run when RT is missing). Hidden; default off.
    #[arg(long = "glyco-transfer-ungated", hide = true, default_value_t = false)]
    glyco_transfer_ungated: bool,

    /// Cross-spectrum transfer: minimum independent-donor graph support to inject a
    /// transferred backbone. Hidden; default 1 (no gate).
    #[arg(long = "glyco-transfer-min-support", hide = true, default_value_t = 1u32)]
    glyco_transfer_min_support: u32,

    /// Cross-spectrum transfer: acceptor-side core-Y quorum (incl. mandatory Y1) to
    /// accept a transfer. Hidden; default 3 (0 disables the gate).
    #[arg(long = "glyco-transfer-core-y", hide = true, default_value_t = 3u8)]
    glyco_transfer_core_y: u8,

    /// Enable cross-spectrum backbone transfer (single-invocation two-pass;
    /// glyco mode only). Pass-1 glyco PSMs are native-GBDT-rescored in-process,
    /// 1%-FDR confident backbones (target+decoy) are propagated to co-eluting
    /// sibling spectra via a glycan-delta graph, and accepted transfers are
    /// re-scored as `Source::Transferred` Pass-2 candidates before the final
    /// `.glyco.pin` is written. Off by default — baseline output is unchanged.
    #[arg(long = "glyco-transfer", default_value_t = false)]
    glyco_transfer: bool,

    /// Enable the PTM-refinement cascade (Pass-2 over confident proteins). Default off.
    #[arg(long = "refine", default_value_t = false)]
    refine: bool,

    /// YAML refinement config; omit to use the built-in 5-mod DEFAULT tier.
    #[arg(long = "refine-config", hide = true)]
    refine_config: Option<std::path::PathBuf>,

    /// Confident-anchor SCOPING FDR (not a reported FDR). Default 0.01 — the same
    /// internal TDC q used for calibration/training/report. A looser gate (e.g.
    /// 0.10) admits low-confidence anchors that leak into the entrapment-FDP
    /// (b1931: 0.10 → 4.86% vs 0.01 → 0.29% true FDP). Hidden power-user knob;
    /// leave at the default unless you have a measured reason to widen it.
    #[arg(long = "refine-select-psm-fdr", default_value_t = 0.01, hide = true, value_parser = parse_unit_fraction)]
    refine_select_psm_fdr: f64,

    /// Run Percolator on the PIN after the search and join its PEP/q-value back
    /// into the outputs (QPX `posterior_error_probability` + a `q-value` score,
    /// and a filtered `<stem>.q<fdr>.tsv`). Needs a Percolator backend (see
    /// `--percolator-bin` / `--percolator-docker`). Default off.
    #[arg(long = "rescore", default_value_t = false)]
    rescore: bool,

    /// Rescore with the built-in NATIVE GBDT rescorer instead of Percolator (no
    /// Percolator backend needed). Leakage-safe 3-fold target-decoy cross-
    /// validation over the PIN features → q-value + PEP. A self-contained
    /// FALLBACK for benchmarking / offline use — NOT production-grade FDR; prefer
    /// `--rescore` (Percolator) for production. Writes the same QPX q-value/PEP +
    /// filtered `<stem>.q<fdr>.tsv` outputs. Ignored if `--rescore` is also set.
    #[arg(long = "rescore-native", hide = true, default_value_t = false)]
    rescore_native: bool,

    /// FDR (q-value) threshold for the filtered `<stem>.q<fdr>.tsv` output
    /// (target PSMs at q ≤ this). Setting it EXPLICITLY without `--rescore` /
    /// `--rescore-native` TRIGGERS rescoring and auto-picks the backend:
    /// Percolator if one is available, otherwise the built-in native rescorer.
    /// When rescoring runs, the threshold defaults to 0.01 if unset.
    #[arg(long = "fdr", value_parser = parse_unit_fraction)]
    fdr: Option<f64>,

    /// Optional per-PSM PEP (posterior error probability / local FDR) cap,
    /// applied IN ADDITION to `--fdr` (a PSM must pass both q ≤ `--fdr` AND
    /// PEP ≤ `--pep`). The q-value stays the primary set-level FDR control;
    /// `--pep` is a supplementary per-PSM gate. Like `--fdr`, setting it
    /// explicitly triggers rescoring. Default: no PEP cap.
    #[arg(long = "pep", hide = true, value_parser = parse_unit_fraction)]
    pep: Option<f64>,

    /// Explicit path to a Percolator binary (highest-priority backend). When
    /// omitted, `percolator` on `$PATH` is used, else the docker fallback.
    #[arg(long = "percolator-bin", hide = true)]
    percolator_bin: Option<std::path::PathBuf>,

    /// Force the Percolator docker fallback (the pinned biocontainers image)
    /// instead of looking for a native binary. Requires the `docker` CLI.
    #[arg(long = "percolator-docker", hide = true, default_value_t = false)]
    percolator_docker: bool,

    /// Percolator docker image tag for the docker fallback (power-user override).
    #[arg(long = "percolator-image", hide = true, default_value = output::DEFAULT_PERCOLATOR_IMAGE)]
    percolator_image: String,

    /// Extra arguments passed verbatim to Percolator (after the fixed flags,
    /// before the PIN path). e.g. `--percolator-args "--testFDR 0.05"`.
    #[arg(long = "percolator-args", hide = true, default_value = "")]
    percolator_args: String,

    /// Keep the PIN file after rescoring. With `--rescore` and no `--output-pin`,
    /// a temporary PIN is used and deleted unless this is true. Default true.
    #[arg(long = "keep-pin", hide = true, default_value_t = true, action = clap::ArgAction::Set)]
    keep_pin: bool,
}

/// Training arguments for `andes train-from-search`.
#[derive(Args, Debug)]
struct TrainFromSearchArgs {

    /// Reuse the seed model's geometry instead of deriving one from the corpus.
    /// Deriving (the default) fits segments and mass tiers to the training data.
    #[arg(long = "seed-geometry", default_value_t = false)]
    seed_geometry: bool,
    /// Input spectrum file (training data). Same format dispatch as for search:
    /// `.mzML`/`.mzml` → mzML reader; anything else → MGF reader.
    ///
    /// Required for initial training.  In `--update` mode with `--remove-source`
    /// or `--reweight` / `--decay`, `--spectra` is only required when
    /// `--validate` is also given (to run the acceptance gate).
    #[arg(long)]
    spectra: Option<PathBuf>,

    /// Input FASTA target database (decoys are generated automatically).
    ///
    /// Required for initial training and for `--update --add`.
    /// In `--update` mode without `--add`, only required when `--validate` is
    /// given.
    #[arg(long)]
    database: Option<PathBuf>,

    /// EXTERNAL LABELS (SP-B glyco training): a TSV with columns `scan`,
    /// `peptide`, `charge`. When given, labels come from this file instead of a
    /// seed search (`bootstrap_labels`), and `--database` is not required. The
    /// `peptide` is the BARE backbone sequence (the glycan is stripped — the
    /// rank model scores backbone b/y ions, which in HCD come from the
    /// deglycosylated backbone); Cam-C and any `--mods` are applied via the
    /// AminoAcidSet. Rows whose scan is absent from `--spectra` or whose peptide
    /// fails to parse are skipped. Use to train a glyco-regime rank model from
    /// the reference engine/a glyco search engine backbone IDs (see docs/plans/glyco/50-roadmap/spb-design.md).
    #[arg(long = "labels")]
    labels: Option<PathBuf>,

    /// Seed model: slug from the bundled store (e.g. `hcd_qexactive_tryp`).
    /// When omitted, the bundled `hcd_qexactive_tryp` model is used as the seed.
    #[arg(long = "seed-model")]
    seed_model: Option<String>,

    /// Target-decoy q-value threshold for accepting PSMs as confident training
    /// labels. Use a lenient value (e.g. 0.1 or 0.5) for small fixtures.
    #[arg(long = "train-fdr", default_value = "0.01", value_parser = parse_unit_fraction)]
    train_fdr: f64,

    /// Instrument tag to embed in the trained model's metadata. Default: `QExactive`.
    #[arg(long, default_value = "QExactive")]
    instrument: String,

    /// Experiment-class / protocol tag (e.g. `Automatic`, `TMT`). Default: `Automatic`.
    #[arg(long, default_value = "Automatic")]
    protocol: String,

    /// Path to the Parquet model store to write (created if absent, appended
    /// otherwise). REQUIRED.
    #[arg(long = "out-store")]
    out_store: PathBuf,

    /// Model ID written into the store. Default: `trained_<instrument>_<protocol>`.
    #[arg(long = "model-id")]
    model_id: Option<String>,

    /// Path to a mods.txt file (same format as `--mods` for search). When
    /// omitted, uses built-in defaults (Carbamidomethyl-C fixed, Oxidation-M + protein-N-term-Acetyl
    /// variable).
    #[arg(long)]
    mods: Option<PathBuf>,

    /// Number of worker threads. Defaults to logical CPU count.
    #[arg(long, default_value_t = num_cpus::get())]
    threads: usize,

    /// ISO 8601 date string (e.g. `2026-01-01`) recorded in the source ledger.
    /// When omitted, the current date is used for initial training; empty string
    /// is stored when `--date ""` is explicitly passed.
    #[arg(long)]
    date: Option<String>,

    // ── Update mode ──────────────────────────────────────────────────────────

    /// Switch to incremental update mode for this model ID.
    /// When set, one of `--add`, `--remove-source`, `--reweight`, or `--decay`
    /// must be provided.
    #[arg(long = "update", value_name = "MODEL_ID")]
    update_model: Option<String>,

    /// (Update mode) Add a new source from `--spectra`.
    /// Requires `--source-id` and `--database`.
    #[arg(long, requires = "update_model")]
    add: bool,

    /// (Update mode) Source identifier for the new source being added
    /// (used with `--add`).
    #[arg(long = "source-id", requires = "add", value_name = "ID")]
    source_id: Option<String>,

    /// (Update mode) Remove the source with this ID from the model.
    #[arg(long = "remove-source", requires = "update_model", value_name = "ID")]
    remove_source: Option<String>,

    /// (Update mode) Set a source's weight.  Format: `<source-id>=<weight>`,
    /// e.g. `--reweight s0=0.5`.
    #[arg(long = "reweight", requires = "update_model", value_name = "ID=W")]
    reweight: Option<String>,

    /// (Update mode) Apply exponential age-decay to all sources with this
    /// half-life in days.
    #[arg(long = "decay", requires = "update_model", value_name = "DAYS")]
    decay: Option<f32>,

    /// (Update mode) Held-out validation spectra for the acceptance gate.
    /// When omitted the gate is skipped (a warning is printed).
    #[arg(long = "validate", requires = "update_model")]
    validate: Option<PathBuf>,

    /// (Update mode) Commit the update even if the acceptance gate fails.
    #[arg(long, requires = "update_model")]
    force: bool,
}

/// Training arguments for `andes train`.
///
/// Trains a scoring model directly from externally-labeled, high-confidence
/// PSMs supplied as a "flat training parquet" (one row per PSM, each carrying
/// the spectrum peaks + identified peptide + resolved mod mass-deltas). This
/// bypasses the bootstrap-search label step entirely: every input row is a
/// label. The seed model supplies only structural hyperparameters (`mme`,
/// deconvolution, segments, frag/precursor offset tables, max_rank); all
/// learned distributions come from the input data.
#[derive(Args, Debug)]
struct TrainArgs {

    /// Reuse the seed model's geometry instead of deriving one from the corpus.
    /// Deriving (the default) fits segments and mass tiers to the training data.
    #[arg(long = "seed-geometry", default_value_t = false)]
    seed_geometry: bool,
    /// Input flat training parquet(s). Repeatable; stats accumulate across all
    /// inputs into a single model.
    #[arg(long = "in", required = true)]
    inputs: Vec<PathBuf>,

    /// Path to the Parquet model store to write (created if absent; existing
    /// models are preserved and re-written alongside the new one). REQUIRED.
    #[arg(long = "out-store")]
    out_store: PathBuf,

    /// Model ID written into the store. Default: `default`.
    #[arg(long = "model-id", default_value = "default")]
    model_id: String,

    /// Seed model: slug from the bundled store (e.g. `hcd_qexactive_tryp`).
    /// Supplies structural hyperparameters only.
    #[arg(long = "seed-model", default_value = "hcd_qexactive_tryp")]
    seed_model: String,

    /// Override the trained model's activation in the store `data_type`
    /// (e.g. `CID`, `HCD`, `ETD`, `UVPD`). Defaults to the seed's value.
    /// Together with `--instrument/--enzyme/--protocol` this lets a new slug
    /// carry the correct selection columns even when seeded from a related model.
    #[arg(long = "activation")]
    activation: Option<String>,

    /// Override the trained model's instrument in the store `data_type`
    /// (e.g. `LowRes`, `HighRes`, `QExactive`, `TOF`). Defaults to the seed's value.
    #[arg(long = "instrument")]
    instrument: Option<String>,

    /// Override the trained model's enzyme in the store `data_type`
    /// (e.g. `Trypsin`, `LysC`, `LysN`). Defaults to the seed's value.
    #[arg(long = "enzyme")]
    enzyme: Option<String>,

    /// Override the trained model's protocol in the store `data_type`
    /// (e.g. `TMT`, `iTRAQ`, `Phosphorylation`, `Automatic`). Drives
    /// `experiment_class` model selection. Defaults to the seed's value.
    #[arg(long = "protocol")]
    protocol: Option<String>,

    /// Fragment match tolerance in ppm. Overwrites the seed model's `mme`
    /// before training. Mutually exclusive with `--fragment-tol-da`. When
    /// neither is given, the seed model's `mme` is kept.
    #[arg(long = "fragment-tol-ppm", conflicts_with = "fragment_tol_da")]
    fragment_tol_ppm: Option<f64>,

    /// Fragment match tolerance in Da. Overwrites the seed model's `mme`
    /// before training. Mutually exclusive with `--fragment-tol-ppm`.
    #[arg(long = "fragment-tol-da")]
    fragment_tol_da: Option<f64>,

    /// Number of worker threads. Defaults to logical CPU count.
    #[arg(long, default_value_t = num_cpus::get())]
    threads: usize,

    /// Laplace pseudo-count for rank/error tables (lower = sharper; default 1.0).
    #[arg(long = "train-pseudo", default_value_t = 1.0)]
    train_pseudo: f32,

    /// Laplace pseudo-count for the NOISE rank distribution (lower = sharper).
    /// Noise is abundant and concentrated, so it needs far less smoothing than
    /// signal ions; the signal `--train-pseudo` over-flattens it. Default 0.05.
    #[arg(long = "train-noise-pseudo", default_value_t = 0.05)]
    train_noise_pseudo: f32,

    /// Partition backoff prior weight (lower = less smoothing toward parent; default 20).
    #[arg(long = "train-backoff-weight", default_value_t = 20.0)]
    train_backoff_weight: f32,

    /// Minimum partition count before backoff blending (default 50).
    #[arg(long = "train-min-count", default_value_t = 50)]
    train_min_count: u64,

    /// Optional path to an independent prior model store. Sparse partitions in
    /// the trained model shrink toward the matching prior model instead of the
    /// corpus-internal pool. Must be own-data (NOT a bundled seed model) to stay
    /// relicense-safe.
    #[arg(long)]
    prior_model_store: Option<PathBuf>,

    /// Model id to load from `--prior-model-store` (defaults to the trained
    /// model id when omitted).
    #[arg(long)]
    prior_model: Option<String>,

    /// Apply widening rank-window smoothing to signal rank distributions
    /// (Kim et al., Nat Commun 5:5277, 2014).
    #[arg(long)]
    rank_smoothing: bool,

    /// Source identifier for the source ledger. Defaults to "msnet".
    #[arg(long, default_value = "msnet")]
    source: String,

    /// Whether to also train and embed a GBDT peak model. `on` (default) trains
    /// GBDT and writes the blob; `off` writes rank-core only (byte-identical to
    /// the pre-GBDT path).
    #[arg(long, default_value = "on")]
    gbdt: GbdtMode,

    /// Opt-in fallback (finding 3.6): downgrade a failed GBDT quality gate to a
    /// warning and embed the degenerate model anyway. Default off.
    #[arg(long = "allow-degenerate-model", hide = true, default_value_t = false)]
    allow_degenerate_model: bool,
}

/// Training arguments for `andes train-intensity`.
///
/// Merges one or more partial intensity aggregation parquets (from
/// `msnet_intensity_agg.py`) into a finalized `intensity_model.parquet` with
/// `mean_log_rel` / `var_log_rel` columns for runtime lookup.
#[derive(Args, Debug)]
struct TrainIntensityArgs {
    /// Input partial or finalized intensity parquets. Repeatable; stats merge
    /// across all inputs.
    #[arg(long = "in", required = true)]
    inputs: Vec<PathBuf>,

    /// Output path for the finalized intensity model parquet.
    #[arg(long = "out", required = true)]
    out: PathBuf,
}

/// GBDT training mode for `andes train`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum)]
enum GbdtMode {
    /// Train and embed a GBDT peak model (default).
    #[default]
    #[clap(name = "on")]
    On,
    /// Skip GBDT; write rank-core only (byte-identical to pre-GBDT path).
    #[clap(name = "off")]
    Off,
}

/// Training arguments for `andes train-intensity-gbdt`.
///
/// Reads flat training parquets (same schema as `train-from-msnet`) and fits a
/// GBDT fragment-intensity regressor (`v3 frag model`).  The trained model is
/// written into `--out-store` alongside any existing models under `--model-id`.
#[derive(Args, Debug)]
struct TrainIntensityGbdtArgs {
    /// Input flat training parquet(s). Repeatable; data accumulate across all
    /// inputs into a single frag-intensity model.
    #[arg(long = "in", required = true)]
    inputs: Vec<PathBuf>,

    /// Path to the Parquet model store to write (created if absent; existing
    /// models are preserved and re-written alongside the new one). REQUIRED.
    #[arg(long = "out-store", required = true)]
    out_store: PathBuf,

    /// Model ID written into the store. Default: `default`.
    #[arg(long = "model-id", default_value = "default")]
    model_id: String,

    /// Seed model: slug from the bundled store (e.g. `hcd_qexactive_tryp`).
    /// Supplies structural hyperparameters
    /// (fragment tolerance, charge range) used when building the frag dataset.
    #[arg(long = "seed-model", default_value = "hcd_qexactive_tryp")]
    seed_model: String,

    /// Number of worker threads (Rayon). Default: 8.
    #[arg(long, default_value_t = 8usize)]
    threads: usize,

    /// Opt-in fallback (finding 3.6): when set, a failed GBDT quality gate
    /// (too few rows / no held-out signal / empty ensemble) is downgraded from a
    /// hard error to a warning and the degenerate model is still written. Default
    /// off — gate failures abort with a non-zero exit. Intended for small
    /// synthetic fixtures / benchmarking only.
    #[arg(long = "allow-degenerate-model", hide = true, default_value_t = false)]
    allow_degenerate_model: bool,
}

/// Training arguments for `andes train-rich-ion-llr`.
///
/// Reads flat training parquets (same schema as `train-from-msnet`) and fits a
/// GBDT rich-ion LLR classifier (logistic; decoy-aware).  The trained model is
/// written into `--out-store` alongside any existing models under `--model-id`.
#[derive(Args, Debug)]
struct TrainRichIonLlrArgs {
    /// Input flat training parquet(s). Repeatable; data accumulate across all
    /// inputs into a single rich-ion model.
    #[arg(long = "in", required = true)]
    inputs: Vec<PathBuf>,

    /// Path to the Parquet model store to write (created if absent; existing
    /// models are preserved and re-written alongside the new one). REQUIRED.
    #[arg(long = "out-store", required = true)]
    out_store: PathBuf,

    /// Model ID written into the store. Default: `default`.
    #[arg(long = "model-id", default_value = "default")]
    model_id: String,

    /// Seed model: slug from the bundled store (e.g. `hcd_qexactive_tryp`).
    /// Supplies structural hyperparameters
    /// (fragment tolerance, charge range) used when building the ion dataset.
    #[arg(long = "seed-model", default_value = "hcd_qexactive_tryp")]
    seed_model: String,

    /// Number of worker threads (Rayon). Default: 8.
    #[arg(long, default_value_t = 8usize)]
    threads: usize,

    /// Opt-in fallback (finding 3.6): downgrade a failed GBDT quality gate to a
    /// warning and write the degenerate model anyway. Default off.
    #[arg(long = "allow-degenerate-model", hide = true, default_value_t = false)]
    allow_degenerate_model: bool,
}

/// Available subcommands.
#[derive(clap::Args, Debug)]
struct RescorePinArgs {
    /// Input PIN.
    #[arg(long = "in")]
    input: PathBuf,
    /// Output target PSMs (Percolator `.psms` shape: PSMId, score, q-value, PEP).
    #[arg(long = "out-psms")]
    out_psms: PathBuf,
    /// Output decoy PSMs.
    #[arg(long = "out-dpsms")]
    out_dpsms: PathBuf,
    /// Cross-validation seed.
    #[arg(long = "seed", default_value_t = 42u64)]
    seed: u64,
}

fn run_rescore_pin(args: RescorePinArgs) -> Result<(), Box<dyn std::error::Error>> {
    let pin_text = std::fs::read_to_string(&args.input)
        .map_err(|e| format!("reading {}: {e}", args.input.display()))?;
    let rows = rescore::native_rescore_qvalues(&pin_text, args.seed)?;
    let mut t = std::io::BufWriter::new(std::fs::File::create(&args.out_psms)?);
    let mut d = std::io::BufWriter::new(std::fs::File::create(&args.out_dpsms)?);
    use std::io::Write;
    for w in [&mut t, &mut d] {
        writeln!(w, "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds")?;
    }
    let (mut nt, mut nd) = (0usize, 0usize);
    for (id, is_decoy, q, score) in &rows {
        let w: &mut dyn Write = if *is_decoy { &mut d } else { &mut t };
        writeln!(w, "{id}\t{score}\t{q}\t{q}\t-\t-")?;
        if *is_decoy { nd += 1 } else { nt += 1 }
    }
    eprintln!("rescore-pin: {nt} target and {nd} decoy rows written");
    Ok(())
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Rescore an existing PIN with the built-in rescorer and write Percolator-shaped
    /// `.psms` / `.dpsms` files.
    ///
    /// Exists so two rescorers can be compared on IDENTICAL input: without it the only
    /// way to try the native rescorer was to re-run the whole search, which changes the
    /// PIN as well as the rescoring and confounds the comparison.
    #[command(name = "rescore-pin", hide = true)]
    RescorePin(RescorePinArgs),

    /// Train a scoring model directly from externally-labeled, high-confidence
    /// PSMs supplied as flat training parquet(s), bypassing the bootstrap
    /// search. This is the primary training path for the Phase-3 "own models".
    ///
    /// Boxed to keep the `Command` enum compact (clippy `large_enum_variant`).
    #[command(hide = true)]
    Train(Box<TrainArgs>),

    /// Train a scoring model from spectra and a FASTA database, writing the
    /// result to a Parquet model store.
    ///
    /// Boxed to keep the `Command` enum compact.
    #[command(name = "train-from-search", hide = true)]
    TrainFromSearch(Box<TrainFromSearchArgs>),

    /// Merge MSNet intensity aggregation parquets into a finalized intensity
    /// model for the strong-score numerator.
    #[command(name = "train-intensity", hide = true)]
    TrainIntensity(Box<TrainIntensityArgs>),

    /// Train a v3 GBDT fragment-intensity regressor from flat training
    /// parquets and embed it in a Parquet model store alongside existing
    /// rank-core models.
    ///
    /// Boxed to keep the `Command` enum compact.
    #[command(name = "train-intensity-gbdt", hide = true)]
    TrainIntensityGbdt(Box<TrainIntensityGbdtArgs>),

    /// Train a GBDT rich-ion LLR classifier (logistic; decoy-aware) from flat
    /// training parquets and embed it in a Parquet model store alongside
    /// existing rank-core models.
    ///
    /// Boxed to keep the `Command` enum compact.
    #[command(name = "train-rich-ion-llr", hide = true)]
    TrainRichIonLlr(Box<TrainRichIonLlrArgs>),
}

/// Top-level CLI.  When no subcommand is given, the flattened `SearchArgs`
/// drive the existing search path (byte-identical to the pre-subcommand
/// behaviour).
#[derive(Parser, Debug)]
#[command(
    name = "andes",
    // `--version` / `-V` from the crate version. Every CLI is expected to answer it,
    // and a search result is only reproducible if you can record which build produced
    // it — so this is provenance, not just convention.
    version = env!("CARGO_PKG_VERSION"),
    about = "andes: database search of MGF/mzML spectra against FASTA",
    allow_hyphen_values = true,
)]
struct TopCli {
    #[command(subcommand)]
    command: Option<Command>,

    #[command(flatten)]
    search: SearchArgs,
}

// Alias used internally for the search-args type.
type Cli = SearchArgs;

fn main() -> ExitCode {
    #[cfg(feature = "thermo")]
    configure_bundled_dotnet();
    // Parse via get_matches so we can query each flag's ValueSource (for the
    // --config merge: an explicit CLI flag must override the YAML value).
    let matches = <TopCli as clap::CommandFactory>::command().get_matches();
    let mut top = <TopCli as clap::FromArgMatches>::from_arg_matches(&matches)
        .unwrap_or_else(|e| e.exit());
    // Record whether the user typed --max-missed-cleavages. `--glyco` raises the
    // floor to 3, which costs real memory (measured +4.4 GB on a 20k-protein human
    // FASTA), so an explicit lower value must be allowed to win. The default (1) is
    // itself a plausible explicit value, so comparing against it cannot distinguish
    // the two cases — only the ValueSource can.
    let _ = EXPLICIT_MISSED_CLEAVAGES.set(matches!(
        matches.value_source("max_missed_cleavages"),
        Some(clap::parser::ValueSource::CommandLine)
    ));
    let result = match top.command.take() {
        Some(Command::RescorePin(args)) => run_rescore_pin(args),
        Some(Command::Train(args)) => run_train(*args),
        Some(Command::TrainFromSearch(args)) => run_train_from_search(*args),
        Some(Command::TrainIntensity(args)) => run_train_intensity(*args),
        Some(Command::TrainIntensityGbdt(args)) => run_train_intensity_gbdt(*args),
        Some(Command::TrainRichIonLlr(args)) => run_train_rich_ion_llr(*args),
        None => {
            // --config: fill any parameter the user did not type on the CLI.
            if let Some(cfg_path) = top.search.config.clone() {
                if let Err(e) = config::RunConfig::load(&cfg_path)
                    .and_then(|c| config::apply(c, &mut top.search, &matches))
                {
                    eprintln!("error: {e}");
                    return ExitCode::from(2);
                }
            }
            // Validate required search args that are Option<> at the clap level.
            let search = top.search;
            if search.spectrum.is_empty() {
                eprintln!("error: --spectrum is required for search (or use `andes train`)");
                return ExitCode::from(2);
            }
            let database = match search.database {
                Some(p) => p,
                None => {
                    eprintln!("error: --database is required for search");
                    return ExitCode::from(2);
                }
            };
            // --output-pin is required UNLESS --rescore is set: rescore can route
            // the search through a temporary PIN (deleted afterwards when
            // --keep-pin false), so a PIN path is not mandatory in that mode.
            if search.output_pin.is_none() && !search.rescore && !search.rescore_native {
                eprintln!("error: --output-pin is required for search (or use --rescore)");
                return ExitCode::from(2);
            }
            // --fdr / --pep are thresholds APPLIED BY a rescoring run; they do not
            // start one on their own. Warn (don't silently launch a backend) if set
            // without an explicit --rescore / --rescore-native.
            if (search.fdr.is_some() || search.pep.is_some())
                && !search.rescore
                && !search.rescore_native
            {
                eprintln!(
                    "warning: --fdr/--pep set without --rescore or --rescore-native; \
                     they are ignored (add --rescore to run Percolator, or feed the PIN \
                     to Percolator yourself)."
                );
            }
            run(Cli {
                database: Some(database),
                ..search
            })
        }
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("andes: {e}");
            ExitCode::from(1)
        }
    }
}

/// Make Thermo `.raw` reading work with zero setup when the runtime is bundled.
///
/// If a .NET runtime ships next to the executable (`<exe_dir>/dotnet`, as the
/// release archives do), point `DOTNET_ROOT` at it so opening a `.raw` "just
/// works". An existing `DOTNET_ROOT` or a system-wide .NET install is left
/// untouched (it takes precedence). No effect on mzML/MGF, which never load .NET.
#[cfg(feature = "thermo")]
fn configure_bundled_dotnet() {
    if std::env::var_os("DOTNET_ROOT").is_some() {
        return;
    }
    if let Ok(exe) = std::env::current_exe() {
        if let Some(dir) = exe.parent() {
            let bundled = dir.join("dotnet");
            if bundled.join("shared").join("Microsoft.NETCore.App").is_dir() {
                std::env::set_var("DOTNET_ROOT", &bundled);
            }
        }
    }
}

/// Print VmRSS for the current process when `ANDES_RSS_PROBE=1`. No-op
/// otherwise and a no-op on non-Linux platforms regardless of the env var.
///
/// We gate behind an env var so production runs stay quiet; flip the var on
/// when debugging memory regressions.
fn log_rss(tag: &str) {
    let probe_set = RSS_PROBE.get().copied().unwrap_or(false);
    if !probe_set {
        return;
    }
    #[cfg(target_os = "linux")]
    {
        if let Ok(s) = std::fs::read_to_string("/proc/self/status") {
            for line in s.lines() {
                if line.starts_with("VmRSS:") {
                    eprintln!(
                        "[RSS {tag}] {}",
                        line.trim_start_matches("VmRSS:").trim()
                    );
                    return;
                }
            }
        }
    }
    #[cfg(not(target_os = "linux"))]
    {
        let _ = tag;
    }
}

/// Statistics returned by the parser-thread helper.
#[derive(Debug, Default)]
struct ParseStats {
    error_count: usize,
    first_errors: Vec<String>,
}

/// Lowercased spectrum-file extension with a trailing `.gz` stripped, so
/// `run.mzML.gz` reports `mzml` rather than `gz`.
///
/// `Path::extension` returns `gz` for a double extension, which silently
/// defeated the `== "mzml"` guards on the metadata-detection helpers below:
/// a gzipped mzML skipped instrument and activation detection entirely and
/// fell back to the low-res default model. The readers those guards protect
/// all use `open_buf_maybe_gz`, so the guard -- not the reader -- was the
/// limitation. Mirrors the `.gz` handling in `input_format_flags`.
fn spectrum_ext_lower(path: &std::path::Path) -> Option<String> {
    let is_gz = path
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.eq_ignore_ascii_case("gz"))
        .unwrap_or(false);
    let effective: std::path::PathBuf = if is_gz {
        path.file_stem()
            .map(std::path::PathBuf::from)
            .unwrap_or_else(|| path.to_path_buf())
    } else {
        path.to_path_buf()
    };
    effective
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.to_ascii_lowercase())
}

fn input_format_flags(path: &Path) -> (bool, bool, bool, bool) {
    // Strip a trailing `.gz` so the format is detected from the underlying
    // extension (`spectra.mzML.gz` → mzML, `spectra.mgf.gz` → MGF). `.raw`/`.d`
    // are binary/directory and never gzipped, so this only affects mzML/MGF.
    let is_gz = path
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.eq_ignore_ascii_case("gz"))
        .unwrap_or(false);
    let effective: std::path::PathBuf = if is_gz {
        path.file_stem().map(std::path::PathBuf::from).unwrap_or_else(|| path.to_path_buf())
    } else {
        path.to_path_buf()
    };
    let ext = effective
        .extension()
        .and_then(|e| e.to_str())
        .map(|s| s.to_lowercase());
    let is_mzml = matches!(ext.as_deref(), Some("mzml"));
    let is_raw = matches!(ext.as_deref(), Some("raw"));
    let is_d = matches!(ext.as_deref(), Some("d"));
    let is_mgf = !is_mzml && !is_raw && !is_d;
    (is_mzml, is_raw, is_d, is_mgf)
}

/// Prefix spectrum titles so pooled multi-file PIN SpecIds stay unique.
/// Decide the per-file SpecId/title prefix. Returns `None` for a single-file
/// search so its PIN output stays byte-identical to the pre-multi-file path;
/// returns `Some("<stem>/")` only when disambiguating across multiple inputs.
fn title_prefix_for(num_files: usize, file_stem: &str) -> Option<String> {
    (num_files > 1).then(|| format!("{file_stem}/"))
}

fn prefix_spectrum_titles(chunk: &mut [Spectrum], prefix: &str) {
    for spec in chunk.iter_mut() {
        if spec.title.is_empty() {
            spec.title = format!("{prefix}scan={}", spec.scan.unwrap_or(0));
        } else {
            spec.title = format!("{prefix}{}", spec.title);
        }
    }
}

/// Build the geometry-derivation [`GeometryConfig`], honouring `ANDES_GEO_*`
/// env overrides so the structural knobs can be swept before settling on fixed
/// defaults. Unset vars fall back to the validated defaults.
/// Derived-geometry parameters for model training.
///
/// Every field has the default the training pipeline shipped with; the `train`
/// subcommands expose them as flags so a run can be reproduced from its command line
/// rather than from an environment the command line does not record.
#[derive(Debug, Clone, Copy, clap::Args)]
pub struct GeometryArgs {
    /// Number of score segments.
    #[arg(long = "geo-segments", default_value_t = 2i32)]
    pub segments: i32,
    /// Maximum peak rank considered.
    #[arg(long = "geo-max-rank", default_value_t = 150i32)]
    pub max_rank: i32,
    /// Target peptides per mass tier.
    #[arg(long = "geo-occupancy", default_value_t = 2500usize)]
    pub occupancy: usize,
    /// Maximum number of mass tiers.
    #[arg(long = "geo-max-tiers", default_value_t = 33usize)]
    pub max_tiers: usize,
    /// Maximum fragment charge modelled.
    #[arg(long = "geo-max-fragment-charge", default_value_t = 3i32)]
    pub max_fragment_charge: i32,
}

impl Default for GeometryArgs {
    fn default() -> Self {
        Self { segments: 2, max_rank: 150, occupancy: 2500, max_tiers: 33, max_fragment_charge: 3 }
    }
}

impl GeometryArgs {
    fn to_config(self) -> GeometryConfig {
        GeometryConfig {
            num_segments: self.segments.max(1),
            max_rank: self.max_rank.max(1),
            mass_tier_occupancy: self.occupancy.max(1),
            max_mass_tiers: self.max_tiers.max(1),
            max_fragment_charge: self.max_fragment_charge.max(1),
        }
    }
}

fn merge_parse_stats(acc: &mut ParseStats, part: ParseStats) {
    acc.error_count += part.error_count;
    for e in part.first_errors {
        if acc.first_errors.len() < 10 {
            acc.first_errors.push(e);
        }
    }
}

/// Producer helper: drains `reader` into fixed-size chunks of `Spectrum`
/// and sends them through `tx`. Stops at `bench_cap` total spectra (or
/// `usize::MAX` for unbounded). Parse errors are counted and the first few
/// captured for downstream reporting; the channel is closed when the
/// reader is exhausted or the consumer hangs up.
///
/// Generic over the reader's error type so the same helper serves both
/// MGF and mzML.
///
/// Runs on a dedicated thread so chunk N+1 is PARSED while chunk N is SCORED.
/// Channel capacity is 2 (one in-flight + one queued) so the producer stays at
/// most one chunk ahead.
fn send_chunks<R, E>(
    reader: R,
    chunk_size: usize,
    bench_cap: usize,
    tx: SyncSender<Vec<Spectrum>>,
) -> ParseStats
where
    R: Iterator<Item = Result<Spectrum, E>>,
    E: std::fmt::Display,
{
    let mut stats = ParseStats::default();
    let mut chunk: Vec<Spectrum> = Vec::with_capacity(chunk_size);
    let mut total = 0usize;
    for result in reader {
        if total >= bench_cap {
            break;
        }
        match result {
            Ok(s) => {
                chunk.push(s);
                total += 1;
                if chunk.len() >= chunk_size {
                    // If the consumer hung up, stop. Sender is moved into the
                    // function, so dropping returns `Err(SendError(chunk))`.
                    let payload = std::mem::replace(&mut chunk, Vec::with_capacity(chunk_size));
                    if tx.send(payload).is_err() {
                        return stats;
                    }
                }
            }
            Err(e) => {
                stats.error_count += 1;
                if stats.first_errors.len() < 3 {
                    stats.first_errors.push(format!("{e}"));
                }
            }
        }
    }
    if !chunk.is_empty() {
        let _ = tx.send(chunk);
    }
    stats
}

/// Lightweight metadata collected in one linear file scan for precursorCal.
#[derive(Debug, Clone)]
struct SpectrumMeta {
    precursor_mz: f64,
    precursor_charge: Option<i32>,
    num_peaks: usize,
}

fn scan_spectrum_metadata(
    path: &Path,
    is_mzml: bool,
    ms_level: u32,
    bench_cap: usize,
) -> Result<Vec<SpectrumMeta>, Box<dyn std::error::Error>> {
    let mut out = Vec::new();
    if is_mzml {
        let reader = MzMLReader::new(input::open_buf_maybe_gz(path)?)
            .with_ms_level_range(ms_level, ms_level);
        for result in reader {
            if out.len() >= bench_cap {
                break;
            }
            let spec = result.map_err(|e| format!("mzML parse: {e}"))?;
            out.push(SpectrumMeta {
                precursor_mz: spec.precursor_mz,
                precursor_charge: spec.precursor_charge,
                num_peaks: spec.peaks.len(),
            });
        }
    } else {
        let reader = MgfReader::new(input::open_buf_maybe_gz(path)?);
        for result in reader {
            if out.len() >= bench_cap {
                break;
            }
            let spec = result.map_err(|e| format!("MGF parse: {e}"))?;
            out.push(SpectrumMeta {
                precursor_mz: spec.precursor_mz,
                precursor_charge: spec.precursor_charge,
                num_peaks: spec.peaks.len(),
            });
        }
    }
    Ok(out)
}

fn build_spec_keys_from_metadata(
    meta: &[SpectrumMeta],
    charge_range: std::ops::RangeInclusive<u8>,
    min_peaks: u32,
) -> Vec<SpecKey> {
    let spectra: Vec<Spectrum> = meta
        .iter()
        .map(|m| Spectrum {
            title: String::new(),
            precursor_mz: m.precursor_mz,
            precursor_intensity: None,
            precursor_charge: m.precursor_charge,
            rt_seconds: None,
            scan: None,
            peaks: vec![(0.0, 0.0); m.num_peaks],
            activation_method: None,
            isolation_lower_offset: None,
            isolation_upper_offset: None,
        })
        .collect();
    build_spec_keys(&spectra, &charge_range, min_peaks)
}

fn load_spectra_by_index(
    path: &Path,
    is_mzml: bool,
    ms_level: u32,
    indices: &HashSet<usize>,
    bench_cap: usize,
) -> Result<HashMap<usize, Spectrum>, Box<dyn std::error::Error>> {
    let mut loaded = HashMap::new();
    if indices.is_empty() {
        return Ok(loaded);
    }
    if is_mzml {
        let reader = MzMLReader::new(input::open_buf_maybe_gz(path)?)
            .with_ms_level_range(ms_level, ms_level);
        for (idx, result) in reader.enumerate() {
            if idx >= bench_cap {
                break;
            }
            if !indices.contains(&idx) {
                continue;
            }
            let spec = result.map_err(|e| format!("mzML parse: {e}"))?;
            loaded.insert(idx, spec);
            if loaded.len() == indices.len() {
                break;
            }
        }
    } else {
        let reader = MgfReader::new(input::open_buf_maybe_gz(path)?);
        for (idx, result) in reader.enumerate() {
            if idx >= bench_cap {
                break;
            }
            if !indices.contains(&idx) {
                continue;
            }
            let spec = result.map_err(|e| format!("MGF parse: {e}"))?;
            loaded.insert(idx, spec);
            if loaded.len() == indices.len() {
                break;
            }
        }
    }
    Ok(loaded)
}

/// Auto-detect an isobaric label (TMT/iTRAQ) by sampling the first `SAMPLE_N`
/// MS2 spectra and inspecting their reporter-ion region. Used only when
/// `--protocol auto` is left at its default, to engage the isobaric windowed
/// peak filter with zero config.
///
/// Returns `None` for `.raw`/`.d` (the sampling reader here is mzML/MGF only —
/// the protocol then stays as-is, byte-identical) and for label-free data, so
/// non-isobaric runs are unchanged. The mzML benchmark datasets (Astral, UPS1,
/// TMT) all flow through the mzML branch.
fn detect_isobaric_sampled(
    path: &Path,
    is_mzml: bool,
    is_mgf: bool,
    ms_level: u32,
    high_res: bool,
) -> Option<input::IsobaricLabel> {
    const SAMPLE_N: usize = 1000;
    if !(is_mzml || is_mgf) {
        return None;
    }
    let indices: HashSet<usize> = (0..SAMPLE_N).collect();
    let loaded = load_spectra_by_index(path, is_mzml, ms_level, &indices, usize::MAX).ok()?;
    let sample: Vec<Spectrum> = loaded.into_values().collect();
    input::detect_isobaric(&sample, high_res)
}

fn tolerance_ppm_display(t: Tolerance) -> Option<f64> {
    match t {
        Tolerance::Ppm(v) => Some(v),
        Tolerance::Da(_) => None,
    }
}

fn run_precursor_calibration(
    spectrum_path: &Path,
    is_mzml: bool,
    ms_level: u32,
    bench_cap: usize,
    params: &SearchParams,
    prepared: &PreparedSearch<'_>,
) -> Result<CalibrationStats, Box<dyn std::error::Error>> {
    if params.precursor_cal_mode == PrecursorCalMode::Off {
        return Ok(CalibrationStats::default());
    }

    let t_cal = std::time::Instant::now();
    let meta = scan_spectrum_metadata(spectrum_path, is_mzml, ms_level, bench_cap)?;
    let spec_keys = build_spec_keys_from_metadata(&meta, params.charge_range.clone(), params.min_peaks);

    if spec_keys.len() < params.cal_min_spec_keys {
        eprintln!(
            "Precursor mass calibration skipped ({} SpecKeys < {} threshold; elapsed: {:.2}s). \
             The sample is too small for a reliable calibration pre-pass.",
            spec_keys.len(),
            params.cal_min_spec_keys,
            t_cal.elapsed().as_secs_f64()
        );
        return Ok(CalibrationStats::default());
    }

    let sampled = sample_every_nth(
        &spec_keys,
        cal_constants::SAMPLING_STRIDE,
        cal_constants::MAX_SAMPLED,
    );
    let needed: HashSet<usize> = sampled.iter().map(|k| k.spectrum_idx).collect();
    let originals = load_spectra_by_index(spectrum_path, is_mzml, ms_level, &needed, bench_cap)?;

    let stats = learn_calibration_stats(&spec_keys, &originals, prepared, params);

    if stats.has_reliable_stats() {
        eprintln!(
            "Precursor mass shift learned: {:.3} ppm from {} confident PSMs (robust sigma {:.3} ppm; elapsed: {:.2}s)",
            stats.shift_ppm,
            stats.confident_psm_count,
            stats.robust_sigma_ppm,
            t_cal.elapsed().as_secs_f64()
        );
    } else {
        eprintln!(
            "Precursor mass calibration skipped (insufficient confident PSMs: {} with PSMs, {} below RawScore floor, {} failed |residual|>50ppm; elapsed: {:.2}s)",
            stats.queues_with_psm,
            stats.rejected_low_score,
            stats.rejected_residual,
            t_cal.elapsed().as_secs_f64()
        );
    }
    Ok(stats)
}

#[cfg(test)]
mod format_routing_tests {
    use super::input_format_flags;
    use std::path::Path;

    // A gzipped spectrum is read transparently (input::open_maybe_gz) and must be
    // routed by its UNDERLYING extension, not the bare `.gz` (finding 2.5).
    #[test]
    fn gz_is_routed_by_the_underlying_extension() {
        // (is_mzml, is_raw, is_d, is_mgf)
        assert_eq!(input_format_flags(Path::new("x/foo.mzML.gz")), (true, false, false, false));
        assert_eq!(input_format_flags(Path::new("foo.MGF.GZ")), (false, false, false, true));
        assert_eq!(input_format_flags(Path::new("foo.mzML")), (true, false, false, false));
        assert_eq!(input_format_flags(Path::new("foo.raw")), (false, true, false, false));
        assert_eq!(input_format_flags(Path::new("run.d")), (false, false, true, false));
        assert_eq!(input_format_flags(Path::new("foo.mgf")), (false, false, false, true));
    }
}

/// Warn BEFORE scoring when the in-RAM candidate index is unlikely to fit.
///
/// Without this, a large database (a whole human proteome at three missed cleavages)
/// runs for half an hour and is then killed by the OOM killer with no message from
/// andes at all — the user sees only a dead process and no output. Measured on a
/// 20k-protein human FASTA: ~0.65 KB resident per candidate for a plain search and
/// ~0.92 KB under `--glyco`, which holds per-spectrum glyco state on top.
///
/// This warns rather than aborts: the estimate is a linear fit, machines differ, and
/// refusing to start a run that would have succeeded is worse than a noisy warning.
/// Only Linux exposes MemAvailable cheaply; elsewhere the check is skipped.
fn warn_if_index_will_not_fit(n_candidates: usize, glyco: bool) {
    const BYTES_PER_CANDIDATE_PLAIN: f64 = 665.0;
    const BYTES_PER_CANDIDATE_GLYCO: f64 = 940.0;

    let available = match std::fs::read_to_string("/proc/meminfo") {
        Ok(text) => text
            .lines()
            .find_map(|l| l.strip_prefix("MemAvailable:"))
            .and_then(|v| v.split_whitespace().next()?.parse::<u64>().ok())
            .map(|kb| kb * 1024),
        Err(_) => None,
    };
    let Some(available) = available else { return };

    let per = if glyco { BYTES_PER_CANDIDATE_GLYCO } else { BYTES_PER_CANDIDATE_PLAIN };
    let estimate = (n_candidates as f64 * per) as u64;
    if estimate <= available {
        return;
    }
    let gb = |b: u64| b as f64 / (1024.0 * 1024.0 * 1024.0);
    eprintln!(
        "WARNING: this search needs roughly {:.1} GB for the in-RAM candidate index \
         ({} candidates) but only {:.1} GB is available. The process is likely to be \
         killed by the operating system partway through, with no result written.",
        gb(estimate), n_candidates, gb(available)
    );
    if glyco {
        eprintln!(
            "  --glyco cannot use the out-of-core index yet (--candidate-index mmap is \
             rejected in glyco mode), so reduce the search space instead: pass \
             --max-missed-cleavages 1 or 2 (glyco defaults to 3), restrict the FASTA to \
             the proteins of interest, or split the database and merge the .glyco.pin \
             files afterwards."
        );
    } else {
        eprintln!("  Re-run with --candidate-index mmap to page the index from disk instead.");
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, clap::ValueEnum)]
enum EthcdActivationFlag {
    /// Label EThcD/ETciD as HCD (default; matches model routing).
    Hcd,
    /// Label them ETD so the c/z scoring path engages.
    Etd,
}

#[derive(Copy, Clone, PartialEq, Eq, Debug, clap::ValueEnum)]
enum GlycoTaxonFlag {
    /// Decide from the data: the NeuGc/NeuAc oxonium ratio across the run, cross-checked
    /// against `OX=` taxon ids in the FASTA. Conservative — only narrows the list when
    /// the evidence supports it, and always says what it decided.
    Auto,
    /// CMAH-inactivated (human and the great apes): exclude NeuGc.
    Human,
    /// CMAH-competent (mouse, rat, pig, bovine, CHO...): keep NeuGc.
    Mammal,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, clap::ValueEnum)]
enum GlycanListFlag {
    /// Reference-fitted human list (HexNAc up to 11, high-antennary reachable). Covers
    /// 100% of a curated 160-composition human reference vs `common`'s 68% -- but MEASURED
    /// WORSE overall on human plasma (228 glycoPSMs @0.00% entrapment FDP vs `common`'s
    /// 365 @0.55%), because the larger space tightens Percolator's threshold. Use only
    /// when the sample genuinely carries high-antennary glycans, and measure.
    ReferenceHuman,
    /// ~600 compositions. The measured-best default; what the benchmarks used.
    Common,
    /// The full ~4,034-composition list. Widest coverage, worst error control.
    Full,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, clap::ValueEnum)]
enum GlycoIsotopeFlag {
    /// 0..=2 — drops the -1 offset, which is pure FDR dilution for glyco.
    Default,
    /// -1..=2 — the pre-round-6 behaviour.
    Negative,
    /// 0..=5 — reaches candidates far above the monoisotopic peak.
    Wide,
}

/// Diagnostic RSS logging, installed from `--rss-probe`.
static RSS_PROBE: std::sync::OnceLock<bool> = std::sync::OnceLock::new();

/// Set in `main` from clap's `ValueSource`: did the user type
/// `--max-missed-cleavages` themselves? See the glyco floor below.
static EXPLICIT_MISSED_CLEAVAGES: std::sync::OnceLock<bool> = std::sync::OnceLock::new();

fn run(cli: Cli) -> Result<(), Box<dyn std::error::Error>> {
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
    let database_path: PathBuf = cli.database.expect("database validated in main");
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

    log_rss("startup");
    let t_total = std::time::Instant::now();
    let t_phase = std::time::Instant::now();
    // ── 1. Load FASTA target database ────────────────────────────────────────
    let target_db =
        FastaReader::load_all(input::open_buf_maybe_gz(&database_path)?)?;
    eprintln!(
        "Loaded {} target proteins from {} [PHASE fasta_load: {:.2}s]",
        target_db.proteins.len(),
        database_path.display(),
        t_phase.elapsed().as_secs_f64()
    );
    log_rss("after_fasta_load");

    // ── 2. Build SearchIndex (targets + strategy-generated decoys) ────────────
    let decoy_strategy = search::decoy::DecoyStrategy::from_name(&cli.decoy_strategy)
        .ok_or_else(|| format!(
            "unknown --decoy-strategy '{}' (expected reverse/shuffle/sequon-reverse/none)",
            cli.decoy_strategy
        ))?;
    let t_phase = std::time::Instant::now();
    let idx = SearchIndex::from_target_db_with_strategy(
        &target_db,
        &cli.decoy_prefix,
        cli.decoy_suffix.as_deref(),
        decoy_strategy,
        cli.decoy_seed,
    );
    eprintln!("[PHASE search_index_build: {:.2}s]", t_phase.elapsed().as_secs_f64());
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
                    Some(method), cli.fragmentation, cli.fragment_tol_ppm, cli.fragment_tol_da,
                ),
                None => resolve_metadataless_selection(
                    None, cli.fragmentation, cli.fragment_tol_ppm, cli.fragment_tol_da,
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
                        search_enzyme.name(), selected
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
    if cli.gbdt_max_trees > 0 {
        let k = cli.gbdt_max_trees;
        let truncate = |slot: &mut Option<Arc<scoring_crate::gbdt_eval::GbdtPeakModel>>, name: &str| {
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
        Protocol::ItraqPhospho => param.data_type.protocol = model::protocol::Protocol::ITRAQPhospho,
        Protocol::Auto => {
            let high_res = param.data_type.instrument.is_high_resolution();
            match detect_isobaric_sampled(&spectrum_path, is_mzml, is_mgf, cli.ms_level as u32, high_res) {
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
    eprintln!("[PHASE param_and_scorer: {:.2}s]", t_phase.elapsed().as_secs_f64());

    // BUG5 (per-spectrum model dispatch prototype, `--glyco-per-spectrum-model`):
    // load a SECOND model, forcing activation=ETD (same instrument/protocol/
    // enzyme lookup that just resolved the primary `param`/`scorer`), so the
    // glyco driver can dispatch each spectrum's rank/edge/hyperscore scoring to
    // whichever model matches ITS OWN activation, instead of the single
    // whole-file majority model every scan uses today (`detect_dominant_
    // activation` picks one model for the whole file by majority vote over the
    // first 64 spectra). `None` (the default — flag off, or no ETD-family model
    // in the store) reproduces the prior single-model behavior exactly; only the
    // glyco path reads this (see `GlycoScoreCtx::etd_scorer`).
    let etd_scorer_owned: Option<RankScorer> = if cli.glyco && cli.glyco_per_spectrum_model {
        let etd_instrument = detected_activation_instrument.and_then(|(_, inst)| inst);
        match load_param_from_store(
            ActivationMethod::ETD,
            etd_instrument,
            cli.protocol,
            search_enzyme,
            cli.model_store.as_deref(),
            None,
        ) {
            Ok((etd_model_id, mut etd_param)) => {
                // Reuse the resolved protocol (the auto-detected TMT/iTRAQ handling
                // applied to the primary `param` above) rather than raw `cli.protocol`,
                // so the ETD scorer's isobaric handling matches the primary scorer.
                etd_param.data_type.protocol = param.data_type.protocol;
                eprintln!(
                    "[glyco] --glyco-per-spectrum-model: loaded ETD model '{etd_model_id}' — \
                     ETD/AI-ETD spectra dispatch to this model for their own rank/edge/\
                     hyperscore/RankScoreFloat scoring; HCD spectra (and the HCD partner read \
                     for --glyco-hcd-pair candidate generation) keep the file's dominant model"
                );
                let mut etd_scorer = RankScorer::new(&etd_param);
                // Apply the same MGF fragment-tol override the primary scorer got
                // (ignored when the instrument was auto-detected from metadata).
                if frag_tol_override.is_some() && !instrument_was_detected {
                    etd_scorer.set_fragment_tol_override(frag_tol_override);
                }
                Some(etd_scorer)
            }
            Err(e) => {
                eprintln!(
                    "WARN: --glyco-per-spectrum-model: failed to load an ETD-family model \
                     ({e}) — falling back to the whole-file dominant model for every spectrum \
                     (no per-spectrum dispatch)"
                );
                None
            }
        }
    } else {
        None
    };

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
            if is_d { "Bruker .d (DDA MS2 only)," } else { "MGF," }
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
            let w: f64 = w.trim().parse().map_err(|_| format!("--peak-filter window `{w}` is not a number"))?;
            let k: usize = k.trim().parse().map_err(|_| format!("--peak-filter count `{k}` is not an integer"))?;
            Some((w, k))
        }
    };
    let _ = RSS_PROBE.set(cli.rss_probe);
    scoring_crate::scoring::init_scoring_settings(scoring_crate::scoring::ScoringSettings {
        peak_filter,
        precursor_offset_clamp: cli.precursor_offset_clamp,
        density_on_active_list: cli.density_on_active_list,
        tight_highres_scoring: cli.tight_highres_scoring,
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
                            &idx, &params, &cli.decoy_prefix, budget,
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
            return Err("--candidate-index mmap is not yet compatible with --chimeric \
                        (the chimeric Pass 2 needs the in-RAM candidate index)"
                .into());
        }
        if cli.refine {
            return Err("--candidate-index mmap is not yet compatible with --refine \
                        (the refinement cascade needs the in-RAM candidate index)"
                .into());
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
        // Multi-file transfer guard (code review): the cross-spectrum transfer
        // join is scan-keyed, and scan numbering conventionally restarts per file, so
        // multiple --spectrum inputs would collide and the fail-loud duplicate-scan
        // check would abort — but only AFTER the full (expensive) search. Warn up front.
        if cli.glyco_transfer && spectrum_paths.len() > 1 {
            eprintln!(
                "WARN: --glyco-transfer with {} --spectrum inputs: the transfer join is \
                 scan-keyed and per-file scan numbering (common in mzML/.raw) will \
                 collide and ABORT the run after the full search. Run --glyco-transfer \
                 per file, or pre-disambiguate scan ids across files.",
                spectrum_paths.len()
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
        // Paired-scan + transfer guard (code review): --glyco-hcd-pair pairing is a
        // Pass-1 property; the --glyco-transfer Pass-2 rescore rebuilds its context
        // WITHOUT the partner map, so an ETD acceptor spectrum would be re-scored
        // unpaired and lose its Pass-1 paired-generation result. Warn so the user
        // knows the combination does not stack for acceptor spectra.
        if cli.glyco_hcd_pair && cli.glyco_transfer {
            eprintln!(
                "WARN: --glyco-hcd-pair with --glyco-transfer: paired-scan generation is \
                 Pass-1 only; the transfer Pass-2 rescores acceptor spectra unpaired, so \
                 pairing does not carry through transfer for those spectra."
            );
        }
    }
    // --refine + --chimeric run together correctly but do NOT currently STACK:
    // the chimeric secondary (co-isolated) PSMs collapse the refinement's
    // confident-anchor set, so refinement adds little on top of chimeric. Warn so
    // the user isn't surprised that the combination ≈ chimeric alone.
    if params.chimeric && cli.refine {
        eprintln!(
            "WARN: --refine + --chimeric do not stack in this release — the chimeric \
             secondary PSMs shrink the refinement's confident-anchor set, so refinement \
             contributes little on top of chimeric. Consider running them separately."
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
    eprintln!("  activation     : {:?} [detected]", param.data_type.activation);
    eprintln!("  instrument     : {:?} [detected]", param.data_type.instrument);
    eprintln!("  protocol       : {:?}", param.data_type.protocol);
    if extra_enzymes.is_empty() {
        eprintln!("  enzyme         : {} ({:?} termini, <={} missed cleavages)",
                  search_enzyme.name(), cli.enzyme_specificity, params.max_missed_cleavages);
    } else {
        let extras: Vec<&str> = extra_enzymes.iter().map(|e| e.name()).collect();
        eprintln!("  enzyme         : {} (primary) + {} (multi-protease union) ({:?} termini, <={} missed cleavages)",
                  search_enzyme.name(), extras.join(","), cli.enzyme_specificity, params.max_missed_cleavages);
    }
    eprintln!("  mods           : {}",
              if cli.mods.is_some() { "from --mods file" }
              else { "defaults (Cam-C fixed, Ox-M variable) + isobaric tag if detected" });
    eprintln!("  max var-mods   : {} per peptide", params.max_variable_mods_per_peptide);
    eprintln!("  peptide length : {}-{}", params.min_length, params.max_length);
    eprintln!("  precursor tol  : {:?} (calibration: {:?})", params.precursor_tolerance, params.precursor_cal_mode);
    eprintln!("  charge range   : {}-{}", params.charge_range.start(), params.charge_range.end());
    eprintln!("  isotope errors : {}..={}", params.isotope_error_range.start(), params.isotope_error_range.end());
    eprintln!("  decoy          : {:?} (prefix {})   chimeric: {}",
              decoy_strategy, cli.decoy_prefix, params.chimeric);
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
        let cal_prepared = PreparedSearch::prepare(
            &idx,
            &params,
            &scorer,
            fragment_tol_da,
            &cli.decoy_prefix,
        );
        let cal_stats = run_precursor_calibration(
            &spectrum_path,
            is_mzml,
            ms_level_u32,
            bench_cap,
            &params,
            &cal_prepared,
        )?;
        let parts = cal_prepared.into_parts();
        params.precursor_mass_shift_ppm = apply_shift_for_mode(params.precursor_cal_mode, cal_stats);
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
                &idx, &params, &scorer, fragment_tol_da, &cli.decoy_prefix, &path,
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

        let (file_is_mzml, file_is_raw, file_is_d, _file_is_mgf) =
            input_format_flags(input_path);
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
                        input::open_buf_maybe_gz(&spectrum_path).map_err(|e| format!("open mzML: {e}"))?,
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
                file_offset, ms1_linked, input_path.display()
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
        let t_glyco = std::time::Instant::now();
        // Use the curated common list (~600 glycans) so that ALL backbone candidates
        // can be b/y-scored in phase-1 (avoids the Y-ladder pre-filter ceiling). The
        // full ~4034-entry list (n_glycan_list()) was A/B-refuted at 1% FDR when
        // the cz collapse scoring was buggy (long/decoy candidates won). The bug
        // hunt (2026-07-16) showed the default ~612 list MISSES the mouse-brain
        // glycome at high charge (z5 69%/z6 38% coverage) — a generation ceiling.
        // ANDES_GLYCO_FULL_GLYCANS re-tests the full list now that cz is fixable.
        // Install scoring settings that reach hot inner functions. Done once here, from
        // validated CLI values, so no scoring code has to read the environment.
        scoring_crate::scoring::init_cz_settings(scoring_crate::scoring::CzSettings {
            zmax_override: cli.glyco_cz_max_charge,
            use_intensity: cli.glyco_cz_intensity,
        });
        andes_glyco::backbone::init_y_max_charge(cli.glyco_y_max_charge);

        let mut glycan_list = match cli.glyco_glycan_list {
            GlycanListFlag::Full => andes_glyco::glycan_db::n_glycan_list(),
            GlycanListFlag::ReferenceHuman => {
                andes_glyco::glycan_db::n_glycan_list_reference_human()
            }
            GlycanListFlag::Common => andes_glyco::glycan_db::n_glycan_list_common(),
        };
        // Decide whether NeuGc belongs in the search space. NeuGc is the sole source of
        // isobaric mass degeneracy in this list (Fuc+NeuGc and Hex+NeuAc are the SAME
        // elemental formula), so getting this right is worth more than any downstream
        // scoring fix -- you cannot resolve from fragments what should not be enumerated.
        let drop_neugc = if cli.glyco_no_neugc {
            eprintln!("--glyco-no-neugc: NeuGc excluded (explicit).");
            true
        } else {
            match cli.glyco_taxon {
                GlycoTaxonFlag::Human => {
                    eprintln!("--glyco-taxon human: NeuGc excluded (CMAH-inactivated lineage).");
                    true
                }
                GlycoTaxonFlag::Mammal => {
                    eprintln!("--glyco-taxon mammal: NeuGc kept (CMAH-competent lineage).");
                    false
                }
                GlycoTaxonFlag::Auto => {
                    // Signal 1: the spectra themselves. Stronger than the FASTA, because it
                    // measures what is in the tube rather than what was searched.
                    let survey = andes_glyco::oxonium::survey_sialic_oxonium(
                        spectra.iter().map(|s| s.peaks.as_slice()),
                        cli.glyco_tol_ppm,
                        0.10,
                    );
                    // Signal 2: OX= taxon ids in the database.
                    let (taxon, n_hu, n_nh, n_ox) = andes_glyco::glycan_db::taxon_from_headers(
                        target_db.proteins.iter().map(|p| p.description.as_str()),
                    );
                    eprintln!(
                        "--glyco-taxon auto: sialylated spectra {} ({} with NeuGc oxonium >=10% of \
                         NeuAc = {:.2}%, {}), FASTA OX= {:?} ({} CMAH-null / {} competent of {})",
                        survey.neuac_spectra,
                        survey.neugc_spectra,
                        100.0 * survey.neugc_fraction,
                        if survey.conclusive { "conclusive" } else { "INCONCLUSIVE" },
                        taxon,
                        n_hu,
                        n_nh,
                        n_ox
                    );
                    // Narrow only when BOTH signals agree, and never on an inconclusive
                    // survey. A CMAH-competent FASTA vetoes: mouse genuinely has NeuGc, and
                    // recombinant human protein from a murine/CHO host does too.
                    let spectra_say_no = survey.conclusive && survey.neugc_fraction < 0.01;
                    let fasta_objects = taxon == andes_glyco::glycan_db::Taxon::CmahCompetent;
                    let decide = spectra_say_no && !fasta_objects;
                    eprintln!(
                        "--glyco-taxon auto: {} (override with --glyco-taxon human|mammal or \
                         --glyco-no-neugc)",
                        if decide {
                            "NeuGc EXCLUDED - no NeuGc oxonium evidence in this run"
                        } else if !survey.conclusive {
                            "NeuGc kept - too little sialic signal to judge"
                        } else if fasta_objects {
                            "NeuGc kept - FASTA is a CMAH-competent organism"
                        } else {
                            "NeuGc kept - NeuGc oxonium evidence present"
                        }
                    );
                    decide
                }
            }
        };
        if drop_neugc {
            let before = glycan_list.len();
            glycan_list.retain(|g| g.neugc == 0);
            eprintln!(
                "glycan list: {} -> {} compositions (NeuGc removed; it is the only source of \
                 isobaric mass degeneracy here)",
                before,
                glycan_list.len()
            );
        }
        let glycan_list = glycan_list;
        let glyco_tol_ppm = cli.glyco_tol_ppm;
        // `!(x > 0.0)` rather than `x <= 0.0` so NaN is rejected too: every
        // comparison against NaN is false, so a NaN tolerance would sail past a
        // `<= 0.0` check and then match nothing, silently, for the whole run.
        if !(glyco_tol_ppm.is_finite() && glyco_tol_ppm > 0.0) {
            eprintln!("error: --glyco-tol-ppm must be a finite value > 0 (got {glyco_tol_ppm})");
            std::process::exit(2);
        }
        if glyco_tol_ppm < 20.0 {
            eprintln!(
                "warning: --glyco-tol-ppm {glyco_tol_ppm} is tighter than the 20 ppm the \
                 glyco defaults were validated at; oxonium, core-Y and c/z matching may \
                 under-fire"
            );
        }
        // Dev cap on glyco scoring uses the global --max-spectra (was the
        // redundant --glyco-max-spectra).
        let spectra_for_glyco: &[_] = if cli.max_spectra > 0 {
            &spectra[..spectra.len().min(cli.max_spectra)]
        } else {
            &spectra
        };
        let glyco_cfg = search::glyco_search::GlycoConfig {
            gp_k: cli.glyco_gp_k,
            gp_j: cli.glyco_gp_j,
            gp_h: cli.glyco_gp_h,
            gp_cz: cli.glyco_gp_cz,
            gp_m: cli.glyco_gp_m,
            min_core_y: cli.glyco_min_core_y,
            min_raw_score: cli.glyco_min_raw_score,
            diag_splits: cli.glyco_diag_splits.clone(),
            min_matched_by: cli.glyco_min_matched_ions,
            max_gen_peaks: cli.glyco_max_peaks,
            cz_multisite: cli.glyco_cz_multisite,
            isobar_rep: cli.glyco_isobar_rep,
            y_index: cli.glyco_y_index,
            decorated_features: cli.glyco_decorated_features,
            sialic_oxonium_min_frac: cli.glyco_sialic_oxonium_min_frac,
            scan_filter_path: cli.glyco_scans.clone(),
            pf_charge: cli.glyco_pf_charge,
            max_pf: cli.glyco_max_pf,
            debug: cli.debug_glyco,
            glyco_decoy: cli.glyco_decoy,
            // Single-file only: cross-file pairing is unsound (see guard above).
            hcd_pair: cli.glyco_hcd_pair && spectrum_paths.len() == 1,
            etd_rank_glycan: cli.glyco_etd_rank_glycan,
            cz_gate: cli.glyco_cz_gate,
        };
        let pass1 = search::glyco_search::glyco_search_run(
            spectra_for_glyco,
            &prepared,
            &glycan_list,
            glyco_tol_ppm,
            cli.glyco_backbone_top_k,
            glyco_cfg.clone(),
            etd_scorer_owned.as_ref(),
        );
        let total_pass1_rows: usize = pass1.iter().map(|r| r.hits.len()).sum();
        eprintln!(
            "[glyco] scored {} spectra → {} glyco-PSM rows [{:.2}s]",
            pass1.len(),
            total_pass1_rows,
            t_glyco.elapsed().as_secs_f64()
        );

        // Derive the glyco PIN path: strip a trailing `.pin` (only), then append
        // `.glyco.pin`.
        //
        // The previous form was `with_extension("").with_extension("glyco.pin")`,
        // which strips at the LAST dot TWICE — so `PXD011533.Frac1.pin` and
        // `PXD011533.Frac2.pin` both collapsed to `PXD011533.glyco.pin` and the
        // second run silently overwrote the first. Pooling fractions is the
        // recommended practice for stable glyco FDR, so `dataset.FracN.pin` is
        // exactly the naming users reach for, and the data loss was silent.
        let glyco_pin_path = {
            let mut s = if output_pin_path
                .extension()
                .is_some_and(|e| e.eq_ignore_ascii_case("pin"))
            {
                output_pin_path.with_extension("").into_os_string()
            } else {
                output_pin_path.clone().into_os_string()
            };
            s.push(".glyco.pin");
            std::path::PathBuf::from(s)
        };
        eprintln!("Glyco PIN will be written to: {}", glyco_pin_path.display());
        // G3: opt-in glycan-axis decoy rows (2D-FDR discrimination on the glycan
        // axis), from --glyco-decoy. Default off — no change to the shipping PIN.
        let emit_glycan_decoy = cli.glyco_decoy;

        // ── Task 8d: single-invocation cross-spectrum backbone transfer ──
        // Off by default (`--glyco-transfer`); when off, behave EXACTLY as
        // before (write Pass-1 unchanged) — this is the gate.
        let glyco_results = if cli.glyco_transfer {
            let t_xfer = std::time::Instant::now();

            // FDR-SOUND as of 2026-07-11 (all 5 design-doc soundness bugs fixed +
            // Codex-reviewed): transferred hits are LOCKED to the seed peptide and
            // carry the seed's target/decoy label (decoy seeds emit decoy rows →
            // symmetric graph), the seed/node joins fail loud on duplicate scans,
            // dedup preserves transfer provenance, the RT gate rejects missing RT by
            // default, and the glycan tolerance is per-acceptor. Validated on
            // PXD025455 Fc3_r1: decoys@1% unchanged (honest), but NET-NEUTRAL — the
            // b/y-dominated winner selection does not promote weak-b/y transferred
            // backbones to top-1 (the conversion wall). See the design doc.
            eprintln!(
                "[glyco-transfer] FDR-sound (peptide-locked, decoy-symmetric). \
                 Net-neutral on Fc3_r1: transfer solves candidate existence but the \
                 selector does not promote weak-b/y transferred backbones. \
                 ANDES_GLYCO_TRANSFER_UNGATED=1 disables the RT gate (unsafe)."
            );

            // Step 2: write Pass-1 to an in-memory PIN, then (step 3)
            // native-GBDT-rescore it in-process to get target+decoy q-values.
            let mut buf: Vec<u8> = Vec::new();
            output::glyco_pin::write_glyco_pin_to(
                &mut buf, &spectra, &pass1, &prepared.candidates, &params, &idx, false, false,
            )?;
            let pin_text = String::from_utf8(buf)
                .map_err(|e| format!("Pass-1 glyco PIN is not valid UTF-8: {e}"))?;
            let q_rows = rescore::native_rescore_qvalues(&pin_text, 42)?;

            // Step 4: seed lookup — reproduce write_glyco_pin's top-1-per-scan
            // collapse winner EXACTLY (same comparator, same enumerated-only
            // gate) so `scan -> peptide_idx/backbone_mass` matches the row the
            // PIN (and thus the q-values above) actually describe.
            // Transfer seed-selection proxy (ladder-primary collapse_cmp).
            let y_primary = true;
            // scan (as emitted into ScanNr/SpecId, i.e. spec.scan.unwrap_or(0))
            // -> (peptide_idx, backbone_mass, rt_seconds, spec_idx).
            let mut seed_lookup: std::collections::BTreeMap<u32, (u32, f64, Option<f64>, usize)> =
                std::collections::BTreeMap::new();
            for r in &pass1 {
                if r.spectrum_idx >= spectra.len() || r.hits.is_empty() {
                    continue;
                }
                let winner = (0..r.hits.len())
                    .max_by(|&a, &b| {
                        andes_glyco::glyco_psm::collapse_cmp(
                            r.hits[a].psm.rank_score,
                            r.hits[a].glycan_key.y_ladder_intensity_score,
                            r.hits[b].psm.rank_score,
                            r.hits[b].glycan_key.y_ladder_intensity_score,
                            y_primary,
                        )
                        .then(b.cmp(&a))
                    })
                    .expect("non-empty hits");
                let hit = &r.hits[winner];
                // enumerated-only gate: a de-novo (unenumerated) winner is not
                // an emitted PIN row (see select_emitted_hits GI-1), so it
                // cannot have a q-value to seed from either.
                if hit.glycan_key.glycan.is_none() {
                    continue;
                }
                let spec = &spectra[r.spectrum_idx];
                // FDR-soundness (design bug #2): the seed→acceptor join is
                // scan-keyed, so a duplicate scan would silently overwrite a seed
                // and map it to the wrong spectrum. Fail loud on missing/duplicate.
                let scan = spec.scan.ok_or_else(|| {
                    format!(
                        "glyco-transfer: Pass-1 winner spectrum {} has no scan number; \
                         transfer requires unique scan ids (seed join is scan-keyed)",
                        r.spectrum_idx
                    )
                })? as u32;
                if seed_lookup
                    .insert(
                        scan,
                        (
                            hit.psm.primary_candidate_idx(),
                            hit.glycan_key.backbone_mass,
                            spec.rt_seconds,
                            r.spectrum_idx,
                        ),
                    )
                    .is_some()
                {
                    return Err(format!(
                        "glyco-transfer: duplicate scan {scan} among Pass-1 glyco winners; \
                         the scan-keyed seed join is not safe with duplicate scans."
                    )
                    .into());
                }
            }

            // Step 5: SeedRow extraction. native_rescore_qvalues already
            // derives is_decoy from the PIN's own fail-loud Label parse (see
            // rescore::parse_pin), so ambiguity is impossible here — a
            // malformed Label would already have failed loud inside
            // native_rescore_qvalues above.
            let rows: Vec<glyco_seeds::SeedRow> = q_rows
                .iter()
                .map(|(spec_id, is_decoy, q, score)| {
                    let scan = glyco_seeds::extract_scan(spec_id).ok_or_else(|| {
                        format!("Task 8d: could not extract scan from rescored SpecId {spec_id:?}")
                    })?;
                    Ok(glyco_seeds::SeedRow { scan, is_decoy: *is_decoy, q_value: *q, score: *score })
                })
                .collect::<Result<Vec<_>, String>>()?;

            // Seed FDR threshold on the NATIVE-GBDT q-value. The in-process
            // native rescorer ranks glyco PSMs as well as Percolator (≈258
            // targets @ native-q≤0.05 vs Percolator's 253 @ q≤0.01 on the same
            // PIN) but its plain target-decoy q-value is more CONSERVATIVE than
            // Percolator's π₀/mix-max-corrected q (its best achievable q floors
            // near ~0.028), so a 0.01 gate yields ZERO seeds even though the
            // confident set exists. 0.05 recovers the Percolator-equivalent
            // confident seed set; final FDR is still Percolator on the merged
            // PIN (the symmetric decoy graph keeps that honest). Tunable for A/B.
            let seed_q: f64 = cli.glyco_transfer_seed_fdr.clamp(f64::MIN_POSITIVE, 1.0);
            let seeds = glyco_seeds::seeds_at_fdr(&rows, seed_q, |scan| {
                seed_lookup.get(&scan).map(|&(pep_idx, bb, rt, _spec_idx)| (pep_idx, bb, rt))
            });

            // Step 7: oxonium-positive spectra as graph nodes, sorted by scan.
            let mut nodes: Vec<andes_glyco::crossspectrum::GlycoNode> = Vec::new();
            let mut scan_to_spec_idx: std::collections::HashMap<u32, usize> =
                std::collections::HashMap::new();
            for (spec_idx, spec) in spectra_for_glyco.iter().enumerate() {
                if spec.peaks.len() < params.min_peaks as usize {
                    continue;
                }
                if !andes_glyco::oxonium::oxonium_gate(&spec.peaks, 0.10, glyco_tol_ppm).fired {
                    continue;
                }
                let z = match spec.precursor_charge {
                    Some(z) if z > 0 => z as f64,
                    _ => continue,
                };
                let precursor_neutral = (spec.precursor_mz - model::mass::PROTON) * z - model::mass::H2O;
                if precursor_neutral <= 0.0 {
                    continue;
                }
                // FDR-soundness (design bug #2): the transfer graph joins seeds and
                // acceptors by scan number. If two spectra share a scan (multi-file
                // input, or MGF without SCANS where scan defaults to 0), the
                // scan→spec_idx map would silently last-wins and seed the WRONG
                // spectrum. Fail loud instead of transferring onto a mislabelled scan.
                let scan = spec.scan.ok_or_else(|| {
                    format!(
                        "glyco-transfer: spectrum {spec_idx} has no scan number; \
                         transfer requires unique scan ids (join is scan-keyed)"
                    )
                })? as u32;
                if let Some(prev) = scan_to_spec_idx.insert(scan, spec_idx) {
                    return Err(format!(
                        "glyco-transfer: duplicate scan {scan} (spectra {prev} and {spec_idx}); \
                         the scan-keyed transfer join is not safe with duplicate scans. \
                         Use single-file input with unique scan ids, or disable --glyco-transfer."
                    )
                    .into());
                }
                nodes.push(andes_glyco::crossspectrum::GlycoNode {
                    scan,
                    precursor_neutral,
                    rt_seconds: spec.rt_seconds,
                });
            }
            nodes.sort_by_key(|n| n.scan);

            // Step 8/9: propagate transfers over the glycan-delta graph, then
            // group provenance-bearing BackboneHits by acceptor spec_idx.
            let glycan_sorted: Vec<(f64, usize)> = {
                let mut v: Vec<(f64, usize)> =
                    glycan_list.iter().enumerate().map(|(i, g)| (g.mass, i)).collect();
                v.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
                v
            };
            let rt_window: f32 = cli.glyco_rt_window;
            const MIN_GLYCAN: f64 = 406.0;
            // FDR-soundness (design bug #4): require an RT co-elution check on both
            // ends by default. `ANDES_GLYCO_TRANSFER_UNGATED=1` is the explicit
            // unsafe opt-in that transfers across the whole run when RT is missing
            // (research only). `propagate_transfers` now scales the glycan-match
            // tolerance PER ACCEPTOR from `glyco_tol_ppm` (design bug #5), so no
            // fixed representative-mass tolerance is passed.
            let require_rt = !cli.glyco_transfer_ungated;
            let transferred = andes_glyco::crossspectrum::propagate_transfers(
                &seeds, &nodes, &glycan_sorted, &glycan_list, rt_window, MIN_GLYCAN,
                glyco_tol_ppm, require_rt,
            );

            // Min graph-support injection gate: only inject a transferred
            // backbone corroborated by >= this many co-eluting, glycan-delta-
            // linked sibling spectra (a real glycoform ladder), cutting the
            // mass-coincidence singletons that the wide "any glycan-delta" edge
            // otherwise floods in. Default 1 (no gate); tune via env for A/B.
            let min_support: u32 = cli.glyco_transfer_min_support;
            // Acceptor-side core-Y acceptance gate (a published spectrum-expansion
            // requirement): a transfer is accepted only if the acceptor spectrum
            // PHYSICALLY shows >= this many core-Y ions WITH Y1 (peptide+HexNAc)
            // present — otherwise transfer floods mass-coincidence candidates onto
            // spectra with no glycan evidence. Default 3; 0 disables the gate.
            let transfer_core_y: u8 = cli.glyco_transfer_core_y;
            let mut gated_out = 0usize;
            let mut injected: std::collections::BTreeMap<usize, Vec<andes_glyco::hybrid::BackboneHit>> =
                std::collections::BTreeMap::new();
            for tc in &transferred {
                if tc.graph_support < min_support {
                    continue;
                }
                let Some(&spec_idx) = scan_to_spec_idx.get(&tc.acceptor_scan) else {
                    continue;
                };
                let charge = spectra_for_glyco[spec_idx]
                    .precursor_charge
                    .filter(|&z| z > 0)
                    .map(|z| z as u8)
                    .unwrap_or(*params.charge_range.start());
                // Core-Y acceptance gate on the acceptor spectrum (bb NEUTRAL = residue + H2O).
                if transfer_core_y > 0 {
                    let spec = &spectra_for_glyco[spec_idx];
                    let stats = andes_glyco::backbone::SpectrumStats::new(&spec.peaks);
                    if !andes_glyco::backbone::acceptor_core_y_gate(
                        &spec.peaks,
                        &stats,
                        tc.backbone_mass + model::mass::H2O,
                        glyco_tol_ppm,
                        charge,
                        transfer_core_y,
                    ) {
                        gated_out += 1;
                        continue;
                    }
                }
                injected.entry(spec_idx).or_default().push(andes_glyco::hybrid::BackboneHit {
                    backbone_mass: tc.backbone_mass,
                    glycan: Some(tc.glycan.clone()),
                    source: andes_glyco::hybrid::Source::Transferred,
                    charge,
                    isotope_offset: 0,
                    glycan_mass_residual: tc.glycan.mass,
                    is_transferred: true,
                    transfer_graph_support: tc.graph_support,
                    transfer_seed_score: tc.seed_score as f32,
                    transfer_rt_delta: tc.rt_delta as f32,
                    transfer_ungated: tc.ungated,
                    // FDR-soundness (design bug #1): lock the transferred hit to the
                    // EXACT seed peptide + carry the seed's decoy label, so Pass-2
                    // scores only that peptide and a decoy seed emits a decoy row.
                    transfer_peptide_idx: Some(tc.peptide_idx),
                    transfer_seed_is_decoy: tc.is_decoy,
                });
            }

            let injected_cands: usize = injected.values().map(|v| v.len()).sum();
            eprintln!(
                "[glyco-transfer] {} Pass-1 rows rescored, {} seeds @{:.1}% native-q ({} decoy), {} nodes, {} transferred candidates -> {} injected (min_support>={}, core_y_gate>={} dropped {}) onto {} acceptor spectra [{:.2}s]",
                q_rows.len(),
                seeds.len(),
                seed_q * 100.0,
                seeds.iter().filter(|s| s.is_decoy).count(),
                nodes.len(),
                transferred.len(),
                injected_cands,
                min_support,
                transfer_core_y,
                gated_out,
                injected.len(),
                t_xfer.elapsed().as_secs_f64()
            );

            // Step 10: Pass-2 re-score only the acceptor spectra, superseding
            // their Pass-1 entry; everything else keeps its Pass-1 result.
            search::glyco_search::glyco_transfer_pass2(
                &spectra,
                &prepared,
                &glycan_list,
                glyco_tol_ppm,
                cli.glyco_backbone_top_k,
                glyco_cfg,
                pass1,
                &injected,
                etd_scorer_owned.as_ref(),
            )
        } else {
            pass1
        };
        let mut glyco_results = glyco_results;
        let total_glyco_rows: usize = glyco_results.iter().map(|r| r.hits.len()).sum();

        // Populate glyco RT PIN features (DeltaRT/AbsDeltaRT/DeltaRTNorm +
        // predicted_rt_min) in place on each hit, using the engine-wide backbone
        // RT index + per-monosaccharide offset + per-run self-calibration. The
        // glyco PIN writer then also appends the within-scan DeltaRTRank. Neutral
        // 0.0 without observed RT / <MIN_CALIBRATION_ANCHORS anchors (baseline-safe).
        if let Some(q) = cli.glyco_min_raw_score_quantile {
            if cli.glyco_min_raw_score.is_some() {
                return Err("--glyco-min-raw-score and --glyco-min-raw-score-quantile are \
                            mutually exclusive: one is an absolute floor, the other derives \
                            the floor from this run's decoy winners"
                    .into());
            }
            let cands = &prepared.candidates;
            let is_decoy = |h: &search::glyco_search::FullGlycoPsm| -> bool {
                h.psm
                    .candidate_idxs
                    .first()
                    .map(|&i| cands[i as usize].is_decoy)
                    .unwrap_or(false)
            };
            match search::glyco_search::apply_adaptive_emission_floor(
                &mut glyco_results,
                &is_decoy,
                q,
            ) {
                Some((floor, before, kept)) => eprintln!(
                    "--glyco-min-raw-score-quantile {q}: derived RawScore floor {floor:.3} \
                     from this run's decoy winners; scans {before} -> {kept}"
                ),
                None => eprintln!(
                    "WARN: --glyco-min-raw-score-quantile {q} did nothing: this run has no \
                     decoy winners to calibrate on (tiny input?); emitting ungated"
                ),
            }
        }

        output::populate_glyco_rt_features(
            &spectra,
            &mut glyco_results,
            &prepared.candidates,
            &glycan_list,
        );

        output::write_glyco_pin(
            &glyco_pin_path,
            &spectra,
            &glyco_results,
            &prepared.candidates,
            &params,
            &idx,
            emit_glycan_decoy,
            cli.debug_glyco,
        )?;
        eprintln!(
            "Wrote glyco PIN: {} ({} PSM rows) [PHASE TOTAL: {:.2}s]",
            glyco_pin_path.display(),
            total_glyco_rows,
            t_total.elapsed().as_secs_f64()
        );
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
            eprintln!("WARN: refine is high-res-only and the data is low-res; skipping refinement.");
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
            eprintln!("[PHASE refinement: {:.2}s]", t_refine.elapsed().as_secs_f64());
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

    output::write_pin(&output_pin_path, &spectra, &queues, &pin_candidates, &params, &pin_index)?;
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
            let extra: Vec<String> =
                cli.percolator_args.split_whitespace().map(|s| s.to_string()).collect();
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
            a.q_value.partial_cmp(&b.q_value).unwrap_or(std::cmp::Ordering::Equal)
        });
        let pep_note = pep_cap.map(|t| format!(" and PEP<={t}")).unwrap_or_default();
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
    let run_stats =
        output::RunStatistics::compute(&queues, &pin_candidates, &params, &param.mme);
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
        let primary_paths: Vec<String> =
            spectrum_paths.iter().map(|p| p.display().to_string()).collect();
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
        output::write_tsv(tsv_path, &spectra, &queues, &pin_candidates, &params, &pin_index, &spec_file_name, is_mgf)?;
        eprintln!("Wrote TSV: {}", tsv_path.display());
    }

    Ok(())
}

/// Write the FDR-filtered rescore TSV: one row per accepted target PSM (already
/// filtered to q ≤ fdr and sorted by q ascending). Columns: the Percolator join
/// key (`PSMId`, == PIN SpecId), `q-value`, `posterior_error_probability`,
/// `peptide`, `proteinIds`.
fn write_filtered_tsv(
    path: &Path,
    psms: &[&output::PercolatorPsm],
) -> std::io::Result<()> {
    use std::io::Write;
    let mut w = std::io::BufWriter::new(File::create(path)?);
    writeln!(w, "PSMId\tq-value\tposterior_error_probability\tpeptide\tproteinIds")?;
    for p in psms {
        writeln!(w, "{}\t{}\t{}\t{}\t{}", p.psm_id, p.q_value, p.pep, p.peptide, p.proteins)?;
    }
    w.flush()
}

// ── Training pipeline ─────────────────────────────────────────────────────────

/// Load all MS2 spectra from a path using the same format-dispatch logic as
/// the search path (mzML by extension, otherwise MGF).
fn load_spectra_for_train(
    path: &Path,
) -> Result<Vec<Spectrum>, Box<dyn std::error::Error>> {
    let ext_lower = spectrum_ext_lower(path);
    let mut spectra = Vec::new();
    match ext_lower.as_deref() {
        Some("mzml") => {
            let reader = MzMLReader::new(input::open_buf_maybe_gz(path)?)
                .with_ms_level_range(2, 2);
            for item in reader {
                match item {
                    Ok(s) => spectra.push(s),
                    Err(e) => eprintln!("WARN: mzML parse: {e}"),
                }
            }
        }
        // Native Thermo `.raw` — MS2 only, same Spectrum stream as the search
        // path. Requires building with `--features thermo`.
        Some("raw") => {
            #[cfg(feature = "thermo")]
            {
                let reader = input::ThermoRawReader::open(path)
                    .map_err(|e| format!("open Thermo .raw {}: {e}", path.display()))?
                    .with_ms_level(Some(2));
                for item in reader {
                    match item {
                        Ok(s) => spectra.push(s),
                        Err(e) => eprintln!("WARN: .raw parse: {e}"),
                    }
                }
            }
            #[cfg(not(feature = "thermo"))]
            {
                return Err(format!(
                    "native Thermo `.raw` training input requires building with \
                     `--features thermo` (and the .NET 8 runtime): {}",
                    path.display()
                )
                .into());
            }
        }
        // Native Bruker timsTOF `.d` (a directory). Requires `--features timstof`.
        Some("d") => {
            #[cfg(feature = "timstof")]
            {
                let reader = input::TimsTofReader::open(path)
                    .map_err(|e| format!("open Bruker .d {}: {e}", path.display()))?;
                for item in reader {
                    match item {
                        Ok(s) => spectra.push(s),
                        Err(e) => eprintln!("WARN: .d parse: {e}"),
                    }
                }
            }
            #[cfg(not(feature = "timstof"))]
            {
                return Err(format!(
                    "native Bruker `.d` training input requires building with \
                     `--features timstof`: {}",
                    path.display()
                )
                .into());
            }
        }
        // MGF (default reader).
        _ => {
            let reader = MgfReader::new(input::open_buf_maybe_gz(path)?);
            for item in reader {
                match item {
                    Ok(s) => spectra.push(s),
                    Err(e) => eprintln!("WARN: MGF parse: {e}"),
                }
            }
        }
    }
    Ok(spectra)
}

/// Build the `SearchParams` used by every training mode (initial bootstrap,
/// `--add`, and the acceptance gate) in one place, so they stay consistent with
/// each other and with the production search:
///   - charge span `2..=5` (the search binary's default; `default_tryptic` alone
///     is the narrow `2..=3`, which drops z=4/5 labels common on Astral/timsTOF),
///   - the `NumMods=` variable-mod limit from the `--mods` file when present
///     (the search path applies it too).
fn build_train_search_params(
    mods: &Option<PathBuf>,
) -> Result<SearchParams, Box<dyn std::error::Error>> {
    let aa = build_aa_set(mods)?;
    let mut params = SearchParams::default_tryptic(aa);
    params.charge_range = 2..=5;
    if let Some(path) = mods {
        if let Some(n) = AminoAcidSetBuilder::parse_num_mods_from_file(path)
            .map_err(|e| format!("parsing NumMods= from {}: {e}", path.display()))?
        {
            params.max_variable_mods_per_peptide = n;
        }
    }
    Ok(params)
}

/// Run the full training pipeline and write a model to a Parquet store.
///
/// When `args.update_model` is set, runs in incremental update mode (Part D).
/// Otherwise runs the standard initial-training pipeline (Part A).
/// Load external training labels from a `scan\tpeptide\tcharge` TSV (SP-B glyco
/// training). Column order is discovered from the header (case-insensitive), so
/// extra columns are ignored. The peptide is parsed with `aa_set` (Cam-C + any
/// `--mods` applied). Rows with an unknown scan or an unparseable peptide/charge
/// are skipped and counted. `confidence` is set to 0.0 (labels are pre-filtered
/// by the external engine; the value is unused by the accumulator).
/// First N-X-S/T sequon position (0-based index of the N; X != P), or None.
/// Used to place the glycan mass for training an ETD c/z model: the glycan rides
/// on glycosite-spanning c/z fragments, so it must be baked onto the sequon N.
fn first_nglyco_site(seq: &[u8]) -> Option<usize> {
    (0..seq.len()).find(|&i| {
        seq[i] == b'N'
            && i + 2 < seq.len()
            && seq[i + 1] != b'P'
            && (seq[i + 2] == b'S' || seq[i + 2] == b'T')
    })
}

fn load_labels_from_tsv(
    path: &std::path::Path,
    spectra: &[model::spectrum::Spectrum],
    aa_set: &model::aa_set::AminoAcidSet,
) -> Result<Vec<model_train::labeled::LabeledMatch>, Box<dyn std::error::Error>> {
    use model::peptide::Peptide;
    use std::collections::HashMap;

    let mut scan_to_idx: HashMap<i32, usize> = HashMap::new();
    for (i, s) in spectra.iter().enumerate() {
        if let Some(sc) = s.scan {
            scan_to_idx.entry(sc).or_insert(i);
        }
    }

    let text = std::fs::read_to_string(path)
        .map_err(|e| format!("reading labels {}: {e}", path.display()))?;
    let mut lines = text.lines();
    let header = lines.next().ok_or("labels file is empty")?;
    let cols: Vec<String> = header.split('\t').map(|c| c.trim().to_ascii_lowercase()).collect();
    let col = |names: &[&str]| -> Option<usize> {
        cols.iter().position(|c| names.iter().any(|n| c == n))
    };
    let scan_c = col(&["scan", "scannr", "scan_number"]).ok_or("labels: no scan column")?;
    let pep_c = col(&["peptide", "backbone", "sequence", "peptide_sequence"])
        .ok_or("labels: no peptide column")?;
    let chg_c = col(&["charge", "z", "precursor_charge"]).ok_or("labels: no charge column")?;
    // Optional glyco columns: when present, the glycan mass is baked onto the
    // glycosite N so glycosite-spanning c/z fragments carry the intact glycan
    // (required to train a correct ETD c/z model — without it the model learns
    // that glycosite c/z ions are "usually missing", which is false).
    let gly_c = col(&["glycan_mass", "glycan", "glycan_neutral"]);
    let site_c = col(&["glycosite", "site", "glyco_site"]);
    let max_c = scan_c.max(pep_c).max(chg_c);

    // Fixed-mod deltas (e.g. Cam-C 57.02146) from the aa_set: the label peptides
    // are BARE sequences, so we annotate each fixed-mod residue with `+delta` and
    // wrap in `-.SEQ.-` flanking so `Peptide::from_str` mass-matches the fixed
    // variant (a bare `C` would otherwise parse as UNMODIFIED — wrong b/y masses).
    // Variable mods (Ox-M) are left off: the label doesn't say which are modified,
    // and the unmodified form is the correct default for the rank corpus.
    let fixed_deltas: std::collections::HashMap<u8, f64> =
        aa_set.fixed_mod_deltas().into_iter().collect();
    let decorate = |seq: &str| -> String {
        let mut d = String::with_capacity(seq.len() + 4);
        d.push_str("-.");
        for &b in seq.as_bytes() {
            d.push(b as char);
            if let Some(delta) = fixed_deltas.get(&b) {
                d.push_str(&format!("+{:.5}", delta));
            }
        }
        d.push_str(".-");
        d
    };

    let mut labels = Vec::new();
    // ONE label per spectrum: the accumulator must not tally the same spectrum
    // multiple times (it would over-weight that spectrum's rank/edge/charge
    // histograms). External glyco exports can carry rank alternatives or
    // glycoform/site alternatives for one scan, so reject duplicate scans loudly
    // rather than silently double-counting (Codex adversarial-review finding).
    let mut seen_scans: std::collections::HashSet<i32> = std::collections::HashSet::new();
    let (mut miss_scan, mut miss_pep, mut miss_other) = (0usize, 0usize, 0usize);
    let (mut dup_scan, mut charge_mismatch) = (0usize, 0usize);
    for line in lines {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() <= max_c {
            miss_other += 1;
            continue;
        }
        let scan: i32 = match f[scan_c].trim().parse() {
            Ok(v) => v,
            Err(_) => { miss_other += 1; continue; }
        };
        let charge: u8 = match f[chg_c].trim().parse() {
            Ok(v) => v,
            Err(_) => { miss_other += 1; continue; }
        };
        let idx = match scan_to_idx.get(&scan) {
            Some(&v) => v,
            None => { miss_scan += 1; continue; }
        };
        // Parse the peptide BEFORE the dedup guard: an unparseable/rank-alternative
        // row must not "claim" the scan and cause a later VALID row for the same
        // scan to be dropped as a duplicate (Codex + code-review finding — external
        // exports can order a weaker alternative first).
        let seq = f[pep_c].trim();
        // Resolve the glycan placement (if this is a glyco corpus): glycan mass
        // from the column, glycosite from an explicit column or the first sequon.
        let glyco: Option<(usize, f64)> = gly_c.and_then(|gc| {
            let gmass: f64 = f.get(gc)?.trim().parse().ok()?;
            if gmass <= 0.0 {
                return None;
            }
            let site = site_c
                .and_then(|sc| f.get(sc)?.trim().parse::<usize>().ok())
                .or_else(|| first_nglyco_site(seq.as_bytes()))?;
            (site < seq.len()).then_some((site, gmass))
        });
        let mut peptide = match Peptide::from_str(&decorate(seq), aa_set) {
            Ok(p) => p,
            Err(_) => { miss_pep += 1; continue; }
        };
        // Bake the glycan directly onto the glycosite residue. We do NOT put it in
        // the peptide string: Peptide::from_str only accepts `+delta` mods that
        // match a REGISTERED variant, and per-PSM glycan masses aren't registered
        // (that dropped ~96% of glyco labels as "unparseable"). The glycan delta on
        // the N flows into every glycosite-spanning c/z fragment mass generically.
        if let Some((site, gmass)) = glyco {
            if let Some(res) = peptide.residues.get_mut(site) {
                let base = res.mod_.as_ref().map_or(0.0, |m| m.mass_delta);
                res.mod_ = Some(std::sync::Arc::new(model::modification::Modification {
                    name: "Glycan".to_string(),
                    mass_delta: base + gmass,
                    residue: model::modification::ResidueSpec::Specific(res.residue),
                    location: model::modification::ModLocation::Anywhere,
                    fixed: false,
                    accession: None,
                    neutral_losses: Vec::new(),
                    loss_class: 1,
                }));
            }
        }
        if !seen_scans.insert(scan) {
            // A second VALID row for a scan already labeled: skip (keep the first).
            dup_scan += 1;
            continue;
        }
        // Charge cross-check: a label charge that disagrees with the spectrum's
        // own precursor charge signals a stale annotation or a scan-mapping error
        // — count it for visibility (the label charge is used for scoring, so a
        // systematic mismatch means the corpus is built on the wrong spectra).
        if let Some(spec_z) = spectra[idx].precursor_charge {
            if spec_z != charge as i32 {
                charge_mismatch += 1;
            }
        }
        labels.push(model_train::labeled::LabeledMatch {
            spectrum_index: idx,
            peptide,
            charge,
            confidence: 0.0,
        });
    }
    eprintln!(
        "train: loaded {} labels from {} ({} skipped: {} scan-not-found, {} dup-scan, {} unparseable-peptide, {} malformed; {} charge-mismatch vs mzML)",
        labels.len(),
        path.display(),
        miss_scan + dup_scan + miss_pep + miss_other,
        miss_scan,
        dup_scan,
        miss_pep,
        miss_other,
        charge_mismatch,
    );
    if charge_mismatch * 5 > labels.len().max(1) {
        eprintln!(
            "train: WARNING — {charge_mismatch} labels ({}%) disagree with the mzML precursor \
             charge; the scan->spectrum mapping or the label charges may be wrong.",
            charge_mismatch * 100 / labels.len().max(1),
        );
    }
    Ok(labels)
}

fn run_train_from_search(args: TrainFromSearchArgs) -> Result<(), Box<dyn std::error::Error>> {
    let t0 = std::time::Instant::now();

    // ── 1. Configure Rayon thread pool ────────────────────────────────────────
    static POOL_INIT: std::sync::OnceLock<()> = std::sync::OnceLock::new();
    POOL_INIT.get_or_init(|| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .expect("build_global");
    });

    if let Some(ref update_model_id) = args.update_model.clone() {
        return run_train_update(args, update_model_id, t0);
    }

    // ── Standard training path ────────────────────────────────────────────────

    // ── 2. Load spectra ───────────────────────────────────────────────────────
    let spectra_path = args.spectra.clone().ok_or("--spectra is required for initial training")?;
    eprintln!("train: loading spectra from {} ...", spectra_path.display());
    let spectra = load_spectra_for_train(&spectra_path)?;
    eprintln!("train: loaded {} spectra", spectra.len());

    // ── 3. Load seed Param + RankScorer ──────────────────────────────────────
    let (seed_model_id, seed_param): (String, Param) = load_seed_param(&args.seed_model)?;
    eprintln!("train: seed model = {seed_model_id}");
    let seed_scorer = RankScorer::new(&seed_param);

    // ── 4-5. Labels: external TSV (SP-B) OR a seed search (bootstrap) ──────────
    let labels = if let Some(ref labels_path) = args.labels {
        // SP-B glyco path: labels come from an external engine's backbone IDs.
        // No seed search / database needed — the peptides are given directly.
        let aa_set = build_aa_set(&args.mods)?;
        load_labels_from_tsv(labels_path, &spectra, &aa_set)?
    } else {
        let database = args
            .database
            .clone()
            .ok_or("--database is required for initial training (or pass --labels)")?;
        let search_params = build_train_search_params(&args.mods)?;
        eprintln!(
            "train: running seed search (train-fdr = {}) ...",
            args.train_fdr
        );
        bootstrap_labels(
            &spectra,
            &database,
            &seed_scorer,
            &search_params,
            args.train_fdr,
        )
        .map_err(|e| format!("bootstrap_labels: {e}"))?
    };
    eprintln!("train: {} confident labels", labels.len());

    if labels.is_empty() {
        return Err(format!(
            "no confident labels found at train-fdr={} — try a higher --train-fdr",
            args.train_fdr
        )
        .into());
    }

    // ── 5b. Geometry template: own-derived (DEFAULT) or seed (opt-out) ─────────
    // andes derives the partition/segment geometry from THIS corpus by default
    // (own geometry — no MS-GF+ partition structure); the seed supplies only
    // non-geometry metadata. Opt out with ANDES_SEED_GEOMETRY=1 to reuse the
    // seed's geometry (e.g. to reproduce a legacy model). Own geometry is
    // entrapment-FDP-validated to beat seed geometry on honest PSMs AND speed
    // across Astral (+57%), UPS1 (+15%) and TMT (+50%).
    let use_seed_geometry = args.seed_geometry;
    let template: Param = if !use_seed_geometry {
        eprintln!(
            "train: deriving own partition geometry from {} PSMs (set ANDES_SEED_GEOMETRY=1 to reuse seed geometry)",
            labels.len()
        );
        let corpus = corpus_charge_masses(&labels);
        let geo_cfg = GeometryArgs::default().to_config();
        derive_geometry(&corpus, &seed_param, &geo_cfg)
    } else {
        eprintln!("train: ANDES_SEED_GEOMETRY set — reusing seed partition geometry");
        seed_param.clone()
    };
    let accum_scorer = RankScorer::new(&template);

    // ── 6. Accumulate stats ───────────────────────────────────────────────────
    eprintln!("train: accumulating ion-match statistics ...");
    let accumulator = StatsAccumulator::new(&accum_scorer);
    let mut stats = CountStats::new();
    for label in &labels {
        let spec = &spectra[label.spectrum_index];
        accumulator.accumulate(&mut stats, spec, &label.peptide, label.charge);
    }
    let stats = merge(vec![stats]);

    // ── 7. Estimate model ─────────────────────────────────────────────────────
    eprintln!("train: estimating model parameters ...");
    let cfg = EstimatorConfig::default();
    let estimator = Estimator::new(cfg);
    let trained_param = estimator.estimate(&stats, &template);
    let n_partitions = trained_param.partitions.len();
    eprintln!("train: trained model has {} partitions", n_partitions);

    // ── 8. Determine model ID ─────────────────────────────────────────────────
    let model_id = args
        .model_id
        .clone()
        .unwrap_or_else(|| format!("trained_{}_{}", args.instrument, args.protocol));
    eprintln!("train: model ID = {model_id}");

    // ── 9. Build source ledger (Part A) ───────────────────────────────────────
    // Determine date: use args.date if provided, else today's date in ISO 8601.
    let date_str = args.date.clone().unwrap_or_else(|| {
        // Format today as YYYY-MM-DD using std::time.
        format_today_iso8601()
    });
    let spectra_filename = spectra_path
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_else(|| spectra_path.display().to_string());
    let source_id = format!("bootstrap_{model_id}");
    let ledger = SourceLedger {
        source_id: source_id.clone(),
        dataset: spectra_filename,
        n_psms: labels.len() as i64,
        date: date_str,
        weight: 1.0,
        train_fdr: args.train_fdr as f32,
        instrument: args.instrument.clone(),
        experiment_class: args.protocol.clone(),
    };

    // ── 10. Write to store with source tracking ────────────────────────────────
    // Read existing OTHER models from the store (preserve them on append).
    let store_path = &args.out_store;
    let mut existing_other: Vec<ModelEntryOwned> = Vec::new();
    if store_path.exists() {
        let store = ModelStore::open(store_path)
            .map_err(|e| format!("opening existing store {}: {e}", store_path.display()))?;
        for id in store.model_ids() {
            if id == model_id {
                eprintln!("train: overwriting existing model '{id}' in store");
                continue;
            }
            let p = store.load_param(&id)
                .map_err(|e| format!("reading model '{id}': {e}"))?;
            let src_ledgers = store.load_sources(&id)
                .unwrap_or_default();
            let mut src = Vec::new();
            for l in src_ledgers {
                if let Ok(s) = store.load_source_stats(&id, &l.source_id) {
                    src.push((l, s));
                }
            }
            existing_other.push((id, p, src));
        }
    }

    // Combine all models: write the trained model + existing others together.
    let mut all_entries: Vec<ModelEntryOwned> = Vec::new();
    all_entries.push((model_id.clone(), trained_param.clone(), vec![(ledger, stats)]));
    for (id, p, src) in existing_other {
        all_entries.push((id, p, src));
    }

    write_all_models_with_sources_pub(
        store_path,
        &all_entries.iter()
            .map(|(id, p, s)| (id.as_str(), p, s.as_slice()))
            .collect::<Vec<_>>(),
    )
    .map_err(|e| format!("writing model store {}: {e}", store_path.display()))?;

    eprintln!(
        "train: wrote model '{model_id}' to {} (source '{source_id}') [{:.2}s]",
        store_path.display(),
        t0.elapsed().as_secs_f64(),
    );
    eprintln!(
        "train: summary — labels={}, partitions={}, store={}",
        labels.len(),
        n_partitions,
        store_path.display(),
    );

    Ok(())
}

// ════════════════════════════════════════════════════════════════════════════
// train-from-msnet: supervised training from externally-labeled PSM parquets
// ════════════════════════════════════════════════════════════════════════════

/// One confident, externally-labeled PSM read from a flat training parquet.
///
/// This is the in-memory form of one parquet row (see the column contract in
/// the `train-from-msnet` CLI docs). The peaks are stored already sorted
/// ascending by m/z (the parquet stores acquisition order; the reader sorts).
struct MsnetPsm {
    spectrum: Spectrum,
    peptide: model::peptide::Peptide,
    charge: u8,
}

/// Build a [`model::peptide::Peptide`] from a bare uppercase sequence plus
/// resolved modification mass-deltas.
///
/// Modifications are applied two ways, matching how the scoring code computes
/// peptide mass (`Peptide::new` sums `aa.mass + mod_.mass_delta` per residue):
///
/// - **Residue mods** (`res_mod_pos` / `res_mod_delta`, 1-based positions):
///   attached to the residue at that position via `AminoAcid::with_mod` with a
///   `ModLocation::Anywhere`, `ResidueSpec::Specific(residue)` modification.
/// - **Terminal mods** (`nterm_delta` / `cterm_delta`): folded onto the first /
///   last residue's mass-delta. Because `Peptide::new` only sums per-residue
///   `mod_.mass_delta`, a terminal delta must be carried by a residue's `mod_`
///   to be counted. If a residue already carries a residue mod, the terminal
///   delta is *added* to that residue's existing delta (a single combined
///   `Modification`); otherwise a fresh terminal `Modification` is attached.
///   This keeps `peptide.mass()` correct regardless of overlap.
///
/// Returns an error if `seq` contains a non-standard residue or the mod arrays
/// are misaligned.
fn build_msnet_peptide(
    seq: &str,
    res_mod_pos: &[i32],
    res_mod_delta: &[f64],
    nterm_delta: f64,
    cterm_delta: f64,
) -> Result<model::peptide::Peptide, Box<dyn std::error::Error>> {
    use std::sync::Arc;

    if res_mod_pos.len() != res_mod_delta.len() {
        return Err(format!(
            "res_mod_pos ({}) and res_mod_delta ({}) length mismatch",
            res_mod_pos.len(),
            res_mod_delta.len()
        )
        .into());
    }
    let bytes = seq.as_bytes();
    if bytes.is_empty() {
        return Err("empty peptide sequence".into());
    }
    let n = bytes.len();

    // Accumulate the total mod delta to apply to each residue (1-based -> 0-based).
    // Residue mods first, then terminal deltas folded onto the end residues.
    let mut residue_delta = vec![0.0f64; n];
    let mut residue_modded = vec![false; n];
    for (&pos1, &delta) in res_mod_pos.iter().zip(res_mod_delta.iter()) {
        if pos1 < 1 || (pos1 as usize) > n {
            return Err(format!(
                "res_mod_pos {pos1} out of range for sequence of length {n}"
            )
            .into());
        }
        let idx = (pos1 - 1) as usize;
        residue_delta[idx] += delta;
        residue_modded[idx] = true;
    }
    if nterm_delta != 0.0 {
        residue_delta[0] += nterm_delta;
        residue_modded[0] = true;
    }
    if cterm_delta != 0.0 {
        residue_delta[n - 1] += cterm_delta;
        residue_modded[n - 1] = true;
    }

    let mut residues = Vec::with_capacity(n);
    for (i, &r) in bytes.iter().enumerate() {
        let aa = model::AminoAcid::standard(r)
            .ok_or_else(|| format!("non-standard residue {:?} at position {}", r as char, i + 1))?;
        if residue_modded[i] {
            let m = Modification {
                name: "msnet".to_string(),
                mass_delta: residue_delta[i],
                residue: ResidueSpec::Specific(r),
                location: ModLocation::Anywhere,
                fixed: false,
                accession: None,
                neutral_losses: Vec::new(),
                loss_class: 0,
            };
            residues.push(aa.with_mod(Arc::new(m)));
        } else {
            residues.push(aa);
        }
    }

    // Flanking residues per the spec: pre=`_`, post=`-`.
    Ok(model::peptide::Peptide::new(residues, b'_', b'-'))
}

/// Read one flat training parquet into a vector of [`MsnetPsm`].
///
/// Reads via the workspace `parquet`/`arrow` crates in record-batch chunks.
/// List columns (`res_mod_pos`, `res_mod_delta`, `mz`, `intensity`) are
/// decoded per-row from their `ListArray` offsets. Peaks are sorted ascending
/// by m/z (the parquet stores acquisition order).
fn read_msnet_parquet(path: &Path) -> Result<Vec<MsnetPsm>, Box<dyn std::error::Error>> {
    use arrow::array::{Array, Float32Array, Float64Array, Int32Array, ListArray, StringArray};
    use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

    let file = File::open(path)
        .map_err(|e| format!("opening {}: {e}", path.display()))?;
    let reader = ParquetRecordBatchReaderBuilder::try_new(file)
        .map_err(|e| format!("parquet reader for {}: {e}", path.display()))?
        .build()
        .map_err(|e| format!("building parquet reader for {}: {e}", path.display()))?;

    let mut out = Vec::new();
    for batch in reader {
        let batch = batch.map_err(|e| format!("reading batch from {}: {e}", path.display()))?;

        let col = |name: &str| -> Result<&std::sync::Arc<dyn Array>, Box<dyn std::error::Error>> {
            batch
                .column_by_name(name)
                .ok_or_else(|| format!("missing column '{name}' in {}", path.display()).into())
        };

        let seq = col("seq")?.as_any().downcast_ref::<StringArray>()
            .ok_or("column 'seq' is not a STRING column")?;
        let charge = col("charge")?.as_any().downcast_ref::<Int32Array>()
            .ok_or("column 'charge' is not an INT32 column")?;
        let prec_mz = col("prec_mz")?.as_any().downcast_ref::<Float64Array>()
            .ok_or("column 'prec_mz' is not a DOUBLE column")?;
        let res_mod_pos = col("res_mod_pos")?.as_any().downcast_ref::<ListArray>()
            .ok_or("column 'res_mod_pos' is not a LIST column")?;
        let res_mod_delta = col("res_mod_delta")?.as_any().downcast_ref::<ListArray>()
            .ok_or("column 'res_mod_delta' is not a LIST column")?;
        let nterm = col("nterm_delta")?.as_any().downcast_ref::<Float64Array>()
            .ok_or("column 'nterm_delta' is not a DOUBLE column")?;
        let cterm = col("cterm_delta")?.as_any().downcast_ref::<Float64Array>()
            .ok_or("column 'cterm_delta' is not a DOUBLE column")?;
        let mz = col("mz")?.as_any().downcast_ref::<ListArray>()
            .ok_or("column 'mz' is not a LIST column")?;
        let intensity = col("intensity")?.as_any().downcast_ref::<ListArray>()
            .ok_or("column 'intensity' is not a LIST column")?;

        // Helper to pull a Vec<i32> out of one ListArray row.
        let list_i32 = |list: &ListArray, i: usize| -> Result<Vec<i32>, Box<dyn std::error::Error>> {
            if list.is_null(i) {
                return Ok(Vec::new());
            }
            let v = list.value(i);
            let a = v.as_any().downcast_ref::<Int32Array>()
                .ok_or("list element is not INT32")?;
            Ok((0..a.len()).map(|j| a.value(j)).collect())
        };
        let list_f64 = |list: &ListArray, i: usize| -> Result<Vec<f64>, Box<dyn std::error::Error>> {
            if list.is_null(i) {
                return Ok(Vec::new());
            }
            let v = list.value(i);
            let a = v.as_any().downcast_ref::<Float64Array>()
                .ok_or("list element is not DOUBLE")?;
            Ok((0..a.len()).map(|j| a.value(j)).collect())
        };
        let list_f32 = |list: &ListArray, i: usize| -> Result<Vec<f32>, Box<dyn std::error::Error>> {
            if list.is_null(i) {
                return Ok(Vec::new());
            }
            let v = list.value(i);
            let a = v.as_any().downcast_ref::<Float32Array>()
                .ok_or("list element is not FLOAT")?;
            Ok((0..a.len()).map(|j| a.value(j)).collect())
        };

        for i in 0..batch.num_rows() {
            let seq_s = seq.value(i);
            let ch = charge.value(i);
            if !(1..=255).contains(&ch) {
                return Err(format!("invalid charge {ch} at row {i} of {}", path.display()).into());
            }
            let charge_u8 = ch as u8;

            let positions = list_i32(res_mod_pos, i)?;
            let deltas = list_f64(res_mod_delta, i)?;
            let peptide = build_msnet_peptide(
                seq_s,
                &positions,
                &deltas,
                nterm.value(i),
                cterm.value(i),
            )?;

            let mzs = list_f32(mz, i)?;
            let ints = list_f32(intensity, i)?;
            if mzs.len() != ints.len() {
                return Err(format!(
                    "mz ({}) and intensity ({}) length mismatch at row {i} of {}",
                    mzs.len(),
                    ints.len(),
                    path.display()
                )
                .into());
            }
            let mut peaks: Vec<(f64, f32)> =
                mzs.iter().zip(ints.iter()).map(|(&m, &it)| (m as f64, it)).collect();
            // Input is acquisition order; the scoring path requires ascending m/z.
            peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

            let spectrum = Spectrum {
                title: format!("row{i}"),
                precursor_mz: prec_mz.value(i),
                precursor_intensity: None,
                precursor_charge: Some(charge_u8 as i32),
                rt_seconds: None,
                scan: None,
                peaks,
                activation_method: None,
                isolation_lower_offset: None,
                isolation_upper_offset: None,
            };

            out.push(MsnetPsm { spectrum, peptide, charge: charge_u8 });
        }
    }
    Ok(out)
}

// train-intensity: merge partial intensity stats into a finalized model parquet
// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct IntensityAggKey {
    ion_type: String,
    flank_n: String,
    flank_c: String,
    pos_bin: i32,
    charge: i32,
    nce_bin: String,
}

#[derive(Debug, Clone, Default)]
struct IntensityAggStats {
    count: i64,
    sum_log_rel: f64,
    sum_log_rel_sq: f64,
}

fn read_intensity_partial(path: &Path) -> Result<Vec<(IntensityAggKey, IntensityAggStats)>, Box<dyn std::error::Error>> {
    use arrow::array::{Array, Float64Array, Int32Array, Int64Array, StringArray};
    use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

    let file = std::fs::File::open(path)
        .map_err(|e| format!("open {}: {e}", path.display()))?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file)
        .map_err(|e| format!("parquet reader for {}: {e}", path.display()))?;
    let has_mean = builder.schema().field_with_name("mean_log_rel").is_ok();
    let mut rows = Vec::new();

    for batch_result in builder.build().map_err(|e| format!("build reader: {e}"))? {
        let batch = batch_result?;
        let ion_col = batch
            .column_by_name("ion_type")
            .and_then(|c| c.as_any().downcast_ref::<StringArray>())
            .ok_or("missing ion_type")?;
        let flank_n_col = batch
            .column_by_name("flank_n")
            .and_then(|c| c.as_any().downcast_ref::<StringArray>())
            .ok_or("missing flank_n")?;
        let flank_c_col = batch
            .column_by_name("flank_c")
            .and_then(|c| c.as_any().downcast_ref::<StringArray>())
            .ok_or("missing flank_c")?;
        let pos_col = batch
            .column_by_name("pos_bin")
            .and_then(|c| c.as_any().downcast_ref::<Int32Array>())
            .ok_or("missing pos_bin")?;
        let charge_col = batch
            .column_by_name("charge")
            .and_then(|c| c.as_any().downcast_ref::<Int32Array>())
            .ok_or("missing charge")?;
        let nce_col = batch
            .column_by_name("nce_bin")
            .and_then(|c| c.as_any().downcast_ref::<StringArray>())
            .ok_or("missing nce_bin")?;
        let count_col = batch
            .column_by_name("count")
            .and_then(|c| c.as_any().downcast_ref::<Int64Array>())
            .ok_or("missing count")?;

        let sum_col = if has_mean {
            None
        } else {
            Some(
                batch
                    .column_by_name("sum_log_rel")
                    .and_then(|c| c.as_any().downcast_ref::<Float64Array>())
                    .ok_or("missing sum_log_rel")?,
            )
        };
        let sum_sq_col = if has_mean {
            None
        } else {
            Some(
                batch
                    .column_by_name("sum_log_rel_sq")
                    .and_then(|c| c.as_any().downcast_ref::<Float64Array>())
                    .ok_or("missing sum_log_rel_sq")?,
            )
        };
        let mean_col = if has_mean {
            Some(
                batch
                    .column_by_name("mean_log_rel")
                    .and_then(|c| c.as_any().downcast_ref::<Float64Array>())
                    .ok_or("missing mean_log_rel")?,
            )
        } else {
            None
        };
        let var_col = if has_mean {
            Some(
                batch
                    .column_by_name("var_log_rel")
                    .and_then(|c| c.as_any().downcast_ref::<Float64Array>())
                    .ok_or("missing var_log_rel")?,
            )
        } else {
            None
        };

        for i in 0..batch.num_rows() {
            let count = count_col.value(i);
            let (sum, sum_sq) = if let (Some(sum_c), Some(sq_c)) = (sum_col, sum_sq_col) {
                (sum_c.value(i), sq_c.value(i))
            } else {
                let mean = mean_col
                    .ok_or("intensity partial: mean_log_rel column missing while in mean/var path")?
                    .value(i);
                let var = var_col
                    .ok_or("intensity partial: var_log_rel column missing while in mean/var path")?
                    .value(i);
                (mean * count as f64, (var + mean * mean) * count as f64)
            };
            let key = IntensityAggKey {
                ion_type: ion_col.value(i).to_string(),
                flank_n: flank_n_col.value(i).to_string(),
                flank_c: flank_c_col.value(i).to_string(),
                pos_bin: pos_col.value(i),
                charge: charge_col.value(i),
                nce_bin: nce_col.value(i).to_string(),
            };
            rows.push((
                key,
                IntensityAggStats {
                    count,
                    sum_log_rel: sum,
                    sum_log_rel_sq: sum_sq,
                },
            ));
        }
    }
    Ok(rows)
}

/// Finalize one aggregation cell into `(mean_log_rel, var_log_rel)`.
/// Returns `None` for `count <= 0` so empty cells (e.g. from a partial
/// aggregation parquet) are dropped instead of writing NaN. Variance is
/// clamped at 0 to absorb floating-point round-off in `E[x²] − E[x]²`.
fn finalize_intensity_stats(sum_log_rel: f64, sum_log_rel_sq: f64, count: i64) -> Option<(f64, f64)> {
    if count <= 0 {
        return None;
    }
    let n = count as f64;
    let mean = sum_log_rel / n;
    let var = (sum_log_rel_sq / n - mean * mean).max(0.0);
    Some((mean, var))
}

/// Sentinel flank residue for flank-marginalized backoff rows (matches the
/// `b'*'` key built by `IntensityModel::predict_log_rel`'s backoff).
const ANY_FLANK: &str = "*";
/// Sentinel nce bin for nce-marginalized backoff rows (matches `"__any__"`).
const ANY_NCE: &str = "__any__";

/// Emit backoff marginal rows so `IntensityModel::predict_log_rel`'s documented
/// sparse-key backoff (drop `nce_bin`, then flank residues) actually finds rows.
///
/// The aggregator only writes exact `(ion, flank_n, flank_c, pos_bin, charge,
/// nce_bin)` cells with real flanks and a real (numeric/`"unknown"`) nce bin.
/// The backoff probes `nce_bin="__any__"` and `flank=b'*'` keys that the
/// aggregator never produces, so without these marginals every lookup whose
/// exact key is absent falls through to the single global mean — making the
/// trained per-context intensities unused. In particular the inference path
/// (`compute_psm_features`) passes `nce_bin="unknown"`, which matches no trained
/// numeric bin, so *every* inference lookup hit the global mean and the
/// `IntensitySignal` cosine became model-value-independent.
///
/// We accumulate each real cell into three marginals: drop nce, drop flanks,
/// and drop both. Marginal keys are disjoint from real keys (real flanks are
/// residues, real nce bins are numeric or `"unknown"`), so snapshotting the
/// real cells first avoids double-counting.
fn add_backoff_marginals(merged: &mut rustc_hash::FxHashMap<IntensityAggKey, IntensityAggStats>) {
    let base: Vec<(IntensityAggKey, IntensityAggStats)> =
        merged.iter().map(|(k, s)| (k.clone(), s.clone())).collect();

    fn bump(
        merged: &mut rustc_hash::FxHashMap<IntensityAggKey, IntensityAggStats>,
        key: IntensityAggKey,
        s: &IntensityAggStats,
    ) {
        let slot = merged.entry(key).or_default();
        slot.count += s.count;
        slot.sum_log_rel += s.sum_log_rel;
        slot.sum_log_rel_sq += s.sum_log_rel_sq;
    }

    for (k, s) in &base {
        // Drop nce: keep flank/pos/charge structure, marginalize over collision energy.
        bump(
            merged,
            IntensityAggKey { nce_bin: ANY_NCE.to_string(), ..k.clone() },
            s,
        );
        // Drop flanks: keep pos/charge/nce, marginalize over the residue pair.
        bump(
            merged,
            IntensityAggKey {
                flank_n: ANY_FLANK.to_string(),
                flank_c: ANY_FLANK.to_string(),
                ..k.clone()
            },
            s,
        );
        // Drop both flanks and nce.
        bump(
            merged,
            IntensityAggKey {
                flank_n: ANY_FLANK.to_string(),
                flank_c: ANY_FLANK.to_string(),
                nce_bin: ANY_NCE.to_string(),
                ..k.clone()
            },
            s,
        );
    }
}

fn write_intensity_model(
    path: &Path,
    merged: &rustc_hash::FxHashMap<IntensityAggKey, IntensityAggStats>,
) -> Result<(), Box<dyn std::error::Error>> {
    use arrow::array::{Float64Array, Int32Array, Int64Array, StringArray};
    use arrow::datatypes::{DataType, Field, Schema};
    use arrow::record_batch::RecordBatch;
    use parquet::arrow::ArrowWriter;

    let mut keys: Vec<_> = merged.iter().collect();
    // Drop empty cells so they never write NaN mean/var (a partial aggregation
    // parquet can carry count==0 rows). finalize_intensity_stats gates on this.
    keys.retain(|(_, s)| s.count > 0);
    keys.sort_by(|a, b| {
        (&a.0.ion_type, &a.0.flank_n, &a.0.flank_c, a.0.pos_bin, a.0.charge, &a.0.nce_bin)
            .cmp(&(
                &b.0.ion_type,
                &b.0.flank_n,
                &b.0.flank_c,
                b.0.pos_bin,
                b.0.charge,
                &b.0.nce_bin,
            ))
    });

    let ion: Vec<_> = keys.iter().map(|(k, _)| k.ion_type.as_str()).collect();
    let flank_n: Vec<_> = keys.iter().map(|(k, _)| k.flank_n.as_str()).collect();
    let flank_c: Vec<_> = keys.iter().map(|(k, _)| k.flank_c.as_str()).collect();
    let pos_bin: Vec<_> = keys.iter().map(|(k, _)| k.pos_bin).collect();
    let charge: Vec<_> = keys.iter().map(|(k, _)| k.charge).collect();
    let nce: Vec<_> = keys.iter().map(|(k, _)| k.nce_bin.as_str()).collect();
    let count: Vec<_> = keys.iter().map(|(_, s)| s.count).collect();
    // Safe to unwrap: zero-count keys were retained out above.
    let mean: Vec<_> = keys
        .iter()
        .map(|(_, s)| finalize_intensity_stats(s.sum_log_rel, s.sum_log_rel_sq, s.count).unwrap().0)
        .collect();
    let var: Vec<_> = keys
        .iter()
        .map(|(_, s)| finalize_intensity_stats(s.sum_log_rel, s.sum_log_rel_sq, s.count).unwrap().1)
        .collect();

    let schema = Schema::new(vec![
        Field::new("ion_type", DataType::Utf8, false),
        Field::new("flank_n", DataType::Utf8, false),
        Field::new("flank_c", DataType::Utf8, false),
        Field::new("pos_bin", DataType::Int32, false),
        Field::new("charge", DataType::Int32, false),
        Field::new("nce_bin", DataType::Utf8, false),
        Field::new("count", DataType::Int64, false),
        Field::new("mean_log_rel", DataType::Float64, false),
        Field::new("var_log_rel", DataType::Float64, false),
    ]);
    let batch = RecordBatch::try_new(
        std::sync::Arc::new(schema.clone()),
        vec![
            std::sync::Arc::new(StringArray::from(ion)),
            std::sync::Arc::new(StringArray::from(flank_n)),
            std::sync::Arc::new(StringArray::from(flank_c)),
            std::sync::Arc::new(Int32Array::from(pos_bin)),
            std::sync::Arc::new(Int32Array::from(charge)),
            std::sync::Arc::new(StringArray::from(nce)),
            std::sync::Arc::new(Int64Array::from(count)),
            std::sync::Arc::new(Float64Array::from(mean)),
            std::sync::Arc::new(Float64Array::from(var)),
        ],
    )?;
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let file = std::fs::File::create(path)?;
    let mut writer = ArrowWriter::try_new(file, std::sync::Arc::new(schema), None)?;
    writer.write(&batch)?;
    writer.close()?;
    Ok(())
}

fn sanity_check_intensity_model(model: &scoring_crate::IntensityModel) -> Result<(), Box<dyn std::error::Error>> {
    use scoring_crate::IntensityIonType;

    // y after K/R should be brighter than b1 at N-terminus (when keys exist).
    let (y_kr, _) = model.predict_log_rel(IntensityIonType::Y, b'K', b'R', 5, 2, "25");
    let (b1, _) = model.predict_log_rel(IntensityIonType::B, b'A', b'L', 1, 2, "25");
    if y_kr <= b1 {
        eprintln!(
            "train-intensity: warning: y(K|R) mean {y_kr:.3} not above b1 mean {b1:.3} (sparse training data?)"
        );
    } else {
        eprintln!("train-intensity: sanity OK: y(K|R)={y_kr:.3} > b1={b1:.3}");
    }
    Ok(())
}

/// `andes train-intensity`: merge partial intensity aggregation parquets.
fn run_train_intensity(args: TrainIntensityArgs) -> Result<(), Box<dyn std::error::Error>> {
    use rustc_hash::FxHashMap;
    use scoring_crate::IntensityModel;

    let t0 = std::time::Instant::now();
    let mut merged: FxHashMap<IntensityAggKey, IntensityAggStats> = FxHashMap::default();
    let mut rows_read = 0usize;

    for input in &args.inputs {
        eprintln!("train-intensity: reading {} ...", input.display());
        let part = read_intensity_partial(input)?;
        rows_read += part.len();
        for (key, stats) in part {
            let slot = merged.entry(key).or_default();
            slot.count += stats.count;
            slot.sum_log_rel += stats.sum_log_rel;
            slot.sum_log_rel_sq += stats.sum_log_rel_sq;
        }
        eprintln!("train-intensity:   {} key rows", rows_read);
    }
    if merged.is_empty() {
        return Err("no intensity key rows read from any --in parquet".into());
    }

    let exact_keys = merged.len();
    add_backoff_marginals(&mut merged);
    eprintln!(
        "train-intensity: {} exact keys + {} backoff marginals = {} total",
        exact_keys,
        merged.len() - exact_keys,
        merged.len()
    );

    write_intensity_model(&args.out, &merged)?;
    eprintln!(
        "train-intensity: wrote {} keys -> {}",
        merged.len(),
        args.out.display()
    );

    let model = IntensityModel::load(&args.out)?;
    sanity_check_intensity_model(&model)?;
    eprintln!("train-intensity: done in {:.1}s", t0.elapsed().as_secs_f64());
    Ok(())
}

/// `andes train-intensity-gbdt`: fit a v3 GBDT fragment-intensity regressor
/// from externally-labeled PSM parquets and embed it in a Parquet model store.
///
/// The function reuses `read_msnet_parquet` / `load_seed_param` / `RankScorer`
/// from the `run_train` path and delegates the store write to
/// `write_all_models_with_sources_and_gbdt_pub`, preserving all other models.
fn run_train_intensity_gbdt(
    args: TrainIntensityGbdtArgs,
) -> Result<(), Box<dyn std::error::Error>> {
    use std::sync::Arc;
    use model_train::gbdt::dataset::PsmRow;
    use model_train::gbdt::frag_dataset::build_frag_dataset;
    use model_train::gbdt::train::{train_gbdt_regression, TrainParams};

    let n_files = args.inputs.len();
    let model_id = args.model_id.clone();
    let seed = args.seed_model.clone();
    let t = args.threads;
    eprintln!("train-intensity-gbdt: in={n_files} model_id={model_id} seed={seed} threads={t}");

    let t0 = std::time::Instant::now();

    // ── 1. Configure Rayon thread pool ────────────────────────────────────────
    static POOL_INIT_FRAG: std::sync::OnceLock<()> = std::sync::OnceLock::new();
    POOL_INIT_FRAG.get_or_init(|| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .expect("build_global");
    });

    // ── 2. Read all input parquets ────────────────────────────────────────────
    let mut psms: Vec<MsnetPsm> = Vec::new();
    for input in &args.inputs {
        eprintln!("train-intensity-gbdt: reading {} ...", input.display());
        let part = read_msnet_parquet(input)?;
        eprintln!("train-intensity-gbdt:   {} PSM rows", part.len());
        psms.extend(part);
    }
    if psms.is_empty() {
        return Err("no PSM rows read from any --in parquet".into());
    }
    eprintln!(
        "train-intensity-gbdt: {} total PSM rows across {} file(s)",
        psms.len(),
        args.inputs.len()
    );

    // ── 3. Load seed Param and build the scorer ───────────────────────────────
    let (seed_model_id, seed_param) = load_seed_param(&Some(args.seed_model.clone()))?;
    eprintln!("train-intensity-gbdt: seed model = {seed_model_id}");
    let seed_scorer = RankScorer::new(&seed_param);

    // ── 4. Build frag-intensity regression dataset ────────────────────────────
    eprintln!(
        "train-intensity-gbdt: building frag dataset from {} PSMs ...",
        psms.len()
    );
    let rows: Vec<PsmRow<'_>> = psms
        .iter()
        .map(|psm| PsmRow {
            spectrum: &psm.spectrum,
            peptide: &psm.peptide,
            charge: psm.charge,
        })
        .collect();
    let ds = build_frag_dataset(&rows, &seed_scorer);
    eprintln!(
        "train-intensity-gbdt: dataset: {} ion rows, {} features",
        ds.y.len(),
        ds.n_features,
    );

    // ── 5. Train the GBDT regressor ───────────────────────────────────────────
    // Hard-error if the trainer's quality gate fails (finding 3.6): a training
    // subcommand must not silently emit a non-deployable model, unless the
    // operator opts into the degenerate fallback.
    let train_params = TrainParams { allow_degenerate: args.allow_degenerate_model, ..TrainParams::default() };
    let trained_frag = train_gbdt_regression(&ds, &train_params, 42)?;
    eprintln!(
        "train-intensity-gbdt: trained frag model: {} trees",
        trained_frag.trees.len()
    );

    // ── 6. Embed the trained model in the seed Param ──────────────────────────
    // The seed Param supplies selection columns (activation/instrument/enzyme/
    // protocol) so model routing works; we stamp the frag model onto it.
    let mut out_param = seed_param;
    out_param.frag_intensity_model = Some(Arc::new(trained_frag));

    // ── 7. Write to store, preserving existing models ─────────────────────────
    let store_path = &args.out_store;

    {
        let mut existing_other: Vec<ModelEntryOwned> = Vec::new();
        let mut existing_blobs: Vec<Option<Vec<u8>>> = Vec::new();
        if store_path.exists() {
            let store = ModelStore::open(store_path)
                .map_err(|e| format!("opening existing store {}: {e}", store_path.display()))?;
            for id in store.model_ids() {
                if id == args.model_id {
                    eprintln!("train-intensity-gbdt: overwriting existing model '{id}' in store");
                    continue;
                }
                let p = store
                    .load_param(&id)
                    .map_err(|e| format!("reading model '{id}': {e}"))?;
                let blob = p.gbdt_peak_model.as_ref().map(|m| m.to_bytes());
                let src_ledgers = store.load_sources(&id).unwrap_or_default();
                let mut src = Vec::new();
                for l in src_ledgers {
                    if let Ok(s) = store.load_source_stats(&id, &l.source_id) {
                        src.push((l, s));
                    }
                }
                existing_other.push((id, p, src));
                existing_blobs.push(blob);
            }
        }

        // New model has no rank-core sources; sources slice is empty.
        let mut all_entries: Vec<ModelEntryOwned> = Vec::new();
        all_entries.push((args.model_id.clone(), out_param, vec![]));
        for (id, p, src) in existing_other {
            all_entries.push((id, p, src));
        }

        // New model carries no separate GBDT peak-model blob (the frag-intensity
        // model is embedded directly on Param.frag_intensity_model, not in the
        // gbdt_model_bytes column).  Existing models' blobs are preserved.
        let mut all_blobs: Vec<Option<Vec<u8>>> = vec![None];
        all_blobs.extend(existing_blobs);

        write_all_models_with_sources_and_gbdt_pub(
            store_path,
            &all_entries
                .iter()
                .map(|(id, p, s)| (id.as_str(), p, s.as_slice()))
                .collect::<Vec<_>>(),
            &all_blobs,
        )
        .map_err(|e| format!("writing model store {}: {e}", store_path.display()))?;
    }

    eprintln!(
        "train-intensity-gbdt: wrote model '{model_id}' to {} [{:.2}s]",
        store_path.display(),
        t0.elapsed().as_secs_f64(),
    );

    Ok(())
}

/// `andes train-rich-ion-llr`: fit a GBDT rich-ion LLR classifier (logistic;
/// decoy-aware) from externally-labeled PSM parquets and embed it in a Parquet
/// model store.
///
/// The function reuses `read_msnet_parquet` / `load_seed_param` / `RankScorer`
/// from the `run_train` path and delegates the store write to
/// `write_all_models_with_sources_and_gbdt_pub`, preserving all other models.
fn run_train_rich_ion_llr(
    args: TrainRichIonLlrArgs,
) -> Result<(), Box<dyn std::error::Error>> {
    use std::sync::Arc;
    use model_train::gbdt::dataset::PsmRow;
    use model_train::gbdt::ion_dataset::build_ion_dataset;
    use model_train::gbdt::train::{train_gbdt, TrainParams};

    let n_files = args.inputs.len();
    let model_id = args.model_id.clone();
    let seed = args.seed_model.clone();
    let t = args.threads;
    eprintln!("train-rich-ion-llr: in={n_files} model_id={model_id} seed={seed} threads={t}");

    let t0 = std::time::Instant::now();

    // ── 1. Configure Rayon thread pool ────────────────────────────────────────
    static POOL_INIT_RICH_ION: std::sync::OnceLock<()> = std::sync::OnceLock::new();
    POOL_INIT_RICH_ION.get_or_init(|| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .expect("build_global");
    });

    // ── 2. Read all input parquets ────────────────────────────────────────────
    let mut psms: Vec<MsnetPsm> = Vec::new();
    for input in &args.inputs {
        eprintln!("train-rich-ion-llr: reading {} ...", input.display());
        let part = read_msnet_parquet(input)?;
        eprintln!("train-rich-ion-llr:   {} PSM rows", part.len());
        psms.extend(part);
    }
    if psms.is_empty() {
        return Err("no PSM rows read from any --in parquet".into());
    }
    eprintln!(
        "train-rich-ion-llr: {} total PSM rows across {} file(s)",
        psms.len(),
        args.inputs.len()
    );

    // ── 3. Load seed Param and build the scorer ───────────────────────────────
    let (seed_model_id, seed_param) = load_seed_param(&Some(args.seed_model.clone()))?;
    eprintln!("train-rich-ion-llr: seed model = {seed_model_id}");
    let seed_scorer = RankScorer::new(&seed_param);

    // ── 4. Build rich-ion classification dataset ──────────────────────────────
    eprintln!(
        "train-rich-ion-llr: building ion dataset from {} PSMs ...",
        psms.len()
    );
    let rows: Vec<PsmRow<'_>> = psms
        .iter()
        .map(|psm| PsmRow {
            spectrum: &psm.spectrum,
            peptide: &psm.peptide,
            charge: psm.charge,
        })
        .collect();
    let ds = build_ion_dataset(&rows, &seed_scorer);
    eprintln!(
        "train-rich-ion-llr: dataset: {} ion rows, {} features",
        ds.y.len(),
        ds.n_features,
    );

    // ── 5. Train the GBDT classifier (logits held-out AUC) ────────────────────
    // Hard-error on quality-gate failure (finding 3.6) unless opted out.
    let train_params = TrainParams { allow_degenerate: args.allow_degenerate_model, ..TrainParams::default() };
    let trained_rich_ion = train_gbdt(&ds, &train_params, 42)?;
    eprintln!(
        "train-rich-ion-llr: trained rich-ion model: {} trees",
        trained_rich_ion.trees.len()
    );

    // ── 6. Embed the trained model in the seed Param ──────────────────────────
    // The seed Param supplies selection columns (activation/instrument/enzyme/
    // protocol) so model routing works; we stamp the rich-ion model onto it.
    let mut out_param = seed_param;
    out_param.rich_ion_model = Some(Arc::new(trained_rich_ion));

    // ── 7. Write to store, preserving existing models ─────────────────────────
    let store_path = &args.out_store;

    {
        let mut existing_other: Vec<ModelEntryOwned> = Vec::new();
        let mut existing_blobs: Vec<Option<Vec<u8>>> = Vec::new();
        if store_path.exists() {
            let store = ModelStore::open(store_path)
                .map_err(|e| format!("opening existing store {}: {e}", store_path.display()))?;
            for id in store.model_ids() {
                if id == args.model_id {
                    eprintln!("train-rich-ion-llr: overwriting existing model '{id}' in store");
                    continue;
                }
                let p = store
                    .load_param(&id)
                    .map_err(|e| format!("reading model '{id}': {e}"))?;
                let blob = p.gbdt_peak_model.as_ref().map(|m| m.to_bytes());
                let src_ledgers = store.load_sources(&id).unwrap_or_default();
                let mut src = Vec::new();
                for l in src_ledgers {
                    if let Ok(s) = store.load_source_stats(&id, &l.source_id) {
                        src.push((l, s));
                    }
                }
                existing_other.push((id, p, src));
                existing_blobs.push(blob);
            }
        }

        // New model has no rank-core sources; sources slice is empty.
        let mut all_entries: Vec<ModelEntryOwned> = Vec::new();
        all_entries.push((args.model_id.clone(), out_param, vec![]));
        for (id, p, src) in existing_other {
            all_entries.push((id, p, src));
        }

        // New model carries no separate GBDT peak-model blob (the rich-ion
        // model is embedded directly on Param.rich_ion_model, not in the
        // gbdt_model_bytes column).  Existing models' blobs are preserved.
        let mut all_blobs: Vec<Option<Vec<u8>>> = vec![None];
        all_blobs.extend(existing_blobs);

        write_all_models_with_sources_and_gbdt_pub(
            store_path,
            &all_entries
                .iter()
                .map(|(id, p, s)| (id.as_str(), p, s.as_slice()))
                .collect::<Vec<_>>(),
            &all_blobs,
        )
        .map_err(|e| format!("writing model store {}: {e}", store_path.display()))?;
    }

    eprintln!(
        "train-rich-ion-llr: wrote model '{model_id}' to {} [{:.2}s]",
        store_path.display(),
        t0.elapsed().as_secs_f64(),
    );

    Ok(())
}

/// `andes train`: train a scoring model directly from externally-labeled PSM
/// parquets, reusing the existing accumulate → estimate → store machinery but
/// bypassing the bootstrap search.
fn run_train(
    args: TrainArgs,
) -> Result<(), Box<dyn std::error::Error>> {
    use rayon::prelude::*;

    let t0 = std::time::Instant::now();

    // ── 1. Configure Rayon thread pool ────────────────────────────────────────
    static POOL_INIT: std::sync::OnceLock<()> = std::sync::OnceLock::new();
    POOL_INIT.get_or_init(|| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .expect("build_global");
    });

    // ── 2. Read all input parquets ────────────────────────────────────────────
    let mut psms: Vec<MsnetPsm> = Vec::new();
    let mut rows_read = 0usize;
    for input in &args.inputs {
        eprintln!("train: reading {} ...", input.display());
        let part = read_msnet_parquet(input)?;
        rows_read += part.len();
        eprintln!("train:   {} PSM rows", part.len());
        psms.extend(part);
    }
    if psms.is_empty() {
        return Err("no PSM rows read from any --in parquet".into());
    }
    eprintln!("train: {rows_read} total PSM rows across {} file(s)", args.inputs.len());

    // ── 3. Load seed Param and apply the fragment-tolerance override ──────────
    let (seed_model_id, mut seed_param): (String, Param) =
        load_seed_param(&Some(args.seed_model.clone()))?;
    eprintln!("train: seed model = {seed_model_id}");
    if let Some(ppm) = args.fragment_tol_ppm {
        seed_param.mme = Tolerance::Ppm(ppm);
        eprintln!("train: fragment tolerance overridden to {ppm} ppm");
    } else if let Some(da) = args.fragment_tol_da {
        seed_param.mme = Tolerance::Da(da);
        eprintln!("train: fragment tolerance overridden to {da} Da");
    } else {
        eprintln!("train: using seed fragment tolerance {:?}", seed_param.mme);
    }

    // ── 3b. Geometry template: own-derived (DEFAULT) or seed (opt-out) ─────────
    // andes derives the partition/segment geometry from THIS corpus by default
    // (own geometry — no MS-GF+ partition structure); the seed supplies only
    // non-geometry metadata. Opt out with ANDES_SEED_GEOMETRY=1 to reuse the
    // seed's geometry (e.g. to reproduce a legacy model). Own geometry is
    // entrapment-FDP-validated to beat seed geometry on honest PSMs AND speed
    // across Astral (+57%), UPS1 (+15%) and TMT (+50%).
    let use_seed_geometry = args.seed_geometry;
    let template: Param = if !use_seed_geometry {
        eprintln!(
            "train: deriving own partition geometry from {} PSMs (set ANDES_SEED_GEOMETRY=1 to reuse seed geometry)",
            psms.len()
        );
        let corpus: Vec<(i32, f32)> = psms
            .iter()
            .map(|p| (p.charge as i32, p.peptide.mass() as f32))
            .collect();
        let geo_cfg = GeometryArgs::default().to_config();
        derive_geometry(&corpus, &seed_param, &geo_cfg)
    } else {
        eprintln!("train: ANDES_SEED_GEOMETRY set — reusing seed partition geometry");
        seed_param.clone()
    };

    // Build the scorer AFTER the tolerance override so accumulation uses it.
    let seed_scorer = RankScorer::new(&template);

    // ── 4. Accumulate ion-match statistics (parallel; per-worker CountStats) ──
    eprintln!("train: accumulating ion-match statistics ...");
    let stats = psms
        .par_iter()
        .fold(
            CountStats::new,
            |mut acc, psm| {
                let accumulator = StatsAccumulator::new(&seed_scorer);
                accumulator.accumulate(&mut acc, &psm.spectrum, &psm.peptide, psm.charge);
                acc
            },
        )
        .collect::<Vec<_>>();
    let stats = merge(stats);
    eprintln!("train: accumulated {} PSMs", psms.len());

    // ── 5. Estimate the model (replaces all learned tables in the seed) ───────
    eprintln!("train: estimating model parameters ...");
    let cfg = EstimatorConfig {
        pseudo: args.train_pseudo,
        noise_pseudo: args.train_noise_pseudo,
        min_count: args.train_min_count,
        backoff_weight: args.train_backoff_weight,
        error_scaling_factor_override: None,
        rank_smoothing: args.rank_smoothing,
    };
    eprintln!(
        "train: estimator pseudo={} noise_pseudo={} backoff_weight={} min_count={}",
        cfg.pseudo, cfg.noise_pseudo, cfg.backoff_weight, cfg.min_count
    );
    let estimator = Estimator::new(cfg);

    // Optional independent prior: sparse partitions shrink toward this model
    // (Level 0 of the backoff hierarchy) instead of the corpus-internal pool.
    // When `--prior-model` is omitted, default to the trained model id. The
    // selection columns passed to `load_param_from_store` are inert here because
    // `model_id_override` is `Some` (it loads that exact id).
    let prior_param: Option<Param> = match &args.prior_model_store {
        Some(store_path) => {
            let prior_id = args.prior_model.clone().unwrap_or_else(|| args.model_id.clone());
            // `load_param_from_store`'s activation/instrument/protocol are only
            // consulted for automatic selection; passing an explicit
            // `model_id_override` makes them inert, so `Protocol::Auto` (the CLI
            // enum the signature expects) is a harmless placeholder.
            let (_pid, p) = load_param_from_store(
                seed_param.data_type.activation,
                Some(seed_param.data_type.instrument),
                Protocol::Auto,
                model::enzyme::Enzyme::Trypsin, // inert (model_id_override below makes selection columns unused)
                Some(store_path.as_path()),
                Some(&prior_id),
            )
            .map_err(|e| format!("loading --prior-model '{prior_id}': {e}"))?;
            eprintln!("train: prior model = {prior_id} (from {})", store_path.display());
            Some(p)
        }
        None => None,
    };

    let mut trained_param =
        estimator.estimate_with_prior(&stats, &template, prior_param.as_ref());
    let n_partitions = trained_param.partitions.len();
    eprintln!("train: trained model has {n_partitions} partitions");

    // ── 5b. Override the selection-relevant data_type from flags ──────────────
    // The trained model inherits the seed's data_type; minting a NEW slug whose
    // (activation, instrument, enzyme, protocol) differs from the seed requires
    // overriding those columns explicitly, otherwise model selection (which keys
    // on these columns, not the model_id string) would never route to it.
    if let Some(act) = &args.activation {
        trained_param.data_type.activation = ActivationMethod::from_name(act)
            .ok_or_else(|| format!("unknown --activation '{act}' (expected CID/HCD/ETD/UVPD/PQD)"))?;
    }
    if let Some(inst) = &args.instrument {
        trained_param.data_type.instrument = InstrumentType::from_name(inst)
            .ok_or_else(|| format!("unknown --instrument '{inst}' (expected LowRes/HighRes/QExactive/TOF)"))?;
    }
    if let Some(enz) = &args.enzyme {
        trained_param.data_type.enzyme = Some(
            model::enzyme::Enzyme::from_name(enz)
                .ok_or_else(|| format!("unknown --enzyme '{enz}' (e.g. Trypsin/LysC/LysN/AspN/GluC/ArgC)"))?,
        );
    }
    if let Some(prot) = &args.protocol {
        trained_param.data_type.protocol = model::protocol::Protocol::from_name(prot)
            .ok_or_else(|| format!("unknown --protocol '{prot}' (expected Automatic/TMT/iTRAQ/iTRAQPhospho/Phosphorylation/Standard)"))?;
    }
    eprintln!(
        "train: model data_type = {:?}/{:?}/{:?}/{:?}",
        trained_param.data_type.activation,
        trained_param.data_type.instrument,
        trained_param.data_type.enzyme,
        trained_param.data_type.protocol,
    );

    // ── 6. Build the source ledger (sentinel train_fdr; pre-labeled input) ────
    let dataset = args
        .inputs
        .first()
        .and_then(|p| p.file_name())
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_else(|| "msnet".to_string());
    let ledger = SourceLedger {
        source_id: args.source.clone(),
        dataset,
        n_psms: psms.len() as i64,
        date: format_today_iso8601(),
        weight: 1.0,
        train_fdr: 1.0, // sentinel: input is pre-labeled, no q-value filtering here
        instrument: trained_param.data_type.instrument.name().to_string(),
        experiment_class: trained_param.data_type.protocol.name().to_string(),
    };

    // ── 6b. Optionally train a GBDT peak model ────────────────────────────────
    let gbdt_blob: Option<Vec<u8>> = if args.gbdt != GbdtMode::Off {
        use model_train::gbdt::dataset::{build_dataset, PsmRow};
        use model_train::gbdt::train::{train_gbdt, TrainParams};

        eprintln!("train: building GBDT dataset from {} PSMs ...", psms.len());
        let gbdt_rows: Vec<PsmRow<'_>> = psms.iter().map(|psm| {
            // Pass the full mod-carrying peptide: labels are mod-aware so b/y
            // ions over Cam-C/TMT/Ox-M land at the correct (shifted) m/z.
            PsmRow {
                spectrum: &psm.spectrum,
                peptide: &psm.peptide,
                charge: psm.charge,
            }
        }).collect();

        let gbdt_dataset = build_dataset(&gbdt_rows, &seed_scorer);
        eprintln!(
            "train: GBDT dataset: {} peak rows, {} positives",
            gbdt_dataset.y.len(),
            gbdt_dataset.y.iter().filter(|&&l| l == 1).count(),
        );

        let gbdt_params = TrainParams { allow_degenerate: args.allow_degenerate_model, ..TrainParams::default() };
        // Hard-error on quality-gate failure (finding 3.6) instead of writing a
        // degenerate model into the store (unless opted out).
        let trained_gbdt = train_gbdt(&gbdt_dataset, &gbdt_params, 42)?;
        eprintln!(
            "train: GBDT trained: {} trees",
            trained_gbdt.trees.len(),
        );
        Some(trained_gbdt.to_bytes())
    } else {
        None
    };

    // ── 7. Write to store, preserving any other existing models ───────────────
    // Unified path: ALWAYS reads existing models from the store, builds the
    // combined entries list (new model + existing others), and writes once via
    // write_all_models_with_sources_and_gbdt_pub.  The GBDT blob (Some or None)
    // flows through for the new model; existing models' blobs are re-serialised
    // from their loaded gbdt_peak_model so no blob is ever lost on re-write.
    let store_path = &args.out_store;
    let model_id = args.model_id.clone();

    {
        let mut existing_other: Vec<ModelEntryOwned> = Vec::new();
        let mut existing_blobs: Vec<Option<Vec<u8>>> = Vec::new();
        if store_path.exists() {
            let store = ModelStore::open(store_path)
                .map_err(|e| format!("opening existing store {}: {e}", store_path.display()))?;
            for id in store.model_ids() {
                if id == model_id {
                    eprintln!("train: overwriting existing model '{id}' in store");
                    continue;
                }
                let p = store.load_param(&id)
                    .map_err(|e| format!("reading model '{id}': {e}"))?;
                // Preserve the GBDT blob (if any) for the existing model.
                let blob = p.gbdt_peak_model.as_ref().map(|m| m.to_bytes());
                let src_ledgers = store.load_sources(&id).unwrap_or_default();
                let mut src = Vec::new();
                for l in src_ledgers {
                    if let Ok(s) = store.load_source_stats(&id, &l.source_id) {
                        src.push((l, s));
                    }
                }
                existing_other.push((id, p, src));
                existing_blobs.push(blob);
            }
        }

        // Build the combined entries list: new model first, then existing others.
        let mut all_entries: Vec<ModelEntryOwned> = Vec::new();
        all_entries.push((model_id.clone(), trained_param, vec![(ledger, stats)]));
        for (id, p, src) in existing_other {
            all_entries.push((id, p, src));
        }

        // Parallel blobs: new model gets gbdt_blob (Some or None); existing models
        // get their preserved blobs.
        let mut all_blobs: Vec<Option<Vec<u8>>> = vec![gbdt_blob];
        all_blobs.extend(existing_blobs);

        write_all_models_with_sources_and_gbdt_pub(
            store_path,
            &all_entries.iter()
                .map(|(id, p, s)| (id.as_str(), p, s.as_slice()))
                .collect::<Vec<_>>(),
            &all_blobs,
        )
        .map_err(|e| format!("writing model store {}: {e}", store_path.display()))?;
    }

    eprintln!(
        "train: wrote model '{model_id}' to {} (source '{}', {} PSMs, {n_partitions} partitions) [{:.2}s]",
        store_path.display(),
        args.source,
        psms.len(),
        t0.elapsed().as_secs_f64(),
    );

    Ok(())
}

/// Incremental update mode (Part D): `--update <MODEL_ID>` plus one of
/// `--add`, `--remove-source`, `--reweight`, `--decay`.
fn run_train_update(
    args: TrainFromSearchArgs,
    model_id: &str,
    t0: std::time::Instant,
) -> Result<(), Box<dyn std::error::Error>> {
    let store_path = &args.out_store;
    let cfg = EstimatorConfig::default();

    // ── Dispatch to the right update operation ────────────────────────────────
    let (candidate, new_sources) = if args.add {
        // --add mode: search spectra, accumulate stats, call update_add.
        let spectra_path = args.spectra.clone()
            .ok_or("--spectra is required with --add")?;
        let database = args.database.clone()
            .ok_or("--database is required with --add")?;
        let source_id = args.source_id.clone()
            .ok_or("--source-id is required with --add")?;

        eprintln!("train update: loading spectra from {} ...", spectra_path.display());
        let spectra = load_spectra_for_train(&spectra_path)?;
        eprintln!("train update: loaded {} spectra", spectra.len());

        // Load the current stored model as the seed.
        let store = ModelStore::open(store_path)
            .map_err(|e| format!("opening store {}: {e}", store_path.display()))?;
        let current_param = store.load_param(model_id)
            .map_err(|e| format!("loading model '{model_id}': {e}"))?;
        let current_scorer = RankScorer::new(&current_param);

        let search_params = build_train_search_params(&args.mods)?;

        eprintln!("train update: running seed search (train-fdr={}) ...", args.train_fdr);
        let labels = bootstrap_labels(
            &spectra,
            &database,
            &current_scorer,
            &search_params,
            args.train_fdr,
        )
        .map_err(|e| format!("bootstrap_labels: {e}"))?;
        eprintln!("train update: {} confident labels", labels.len());

        if labels.is_empty() {
            return Err(format!(
                "no confident labels at train-fdr={} — try a higher --train-fdr",
                args.train_fdr
            ).into());
        }

        let accumulator = StatsAccumulator::new(&current_scorer);
        let mut stats = CountStats::new();
        for label in &labels {
            accumulator.accumulate(&mut stats, &spectra[label.spectrum_index], &label.peptide, label.charge);
        }
        let stats = merge(vec![stats]);

        let spectra_filename = spectra_path
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_else(|| spectra_path.display().to_string());
        let date_str = args.date.clone().unwrap_or_else(format_today_iso8601);
        let ledger = SourceLedger {
            source_id: source_id.clone(),
            dataset: spectra_filename,
            n_psms: labels.len() as i64,
            date: date_str,
            weight: 1.0,
            train_fdr: args.train_fdr as f32,
            instrument: args.instrument.clone(),
            experiment_class: args.protocol.clone(),
        };

        update_add(store_path, model_id, ledger, stats, cfg)
            .map_err(|e| format!("update_add: {e}"))?

    } else if let Some(ref sid) = args.remove_source.clone() {
        update_remove(store_path, model_id, sid, cfg)
            .map_err(|e| format!("update_remove: {e}"))?

    } else if let Some(ref spec) = args.reweight.clone() {
        // Parse "source-id=weight"
        let (sid, weight) = parse_reweight_spec(spec)?;
        update_reweight(store_path, model_id, &sid, weight, cfg)
            .map_err(|e| format!("update_reweight: {e}"))?

    } else if let Some(half_life) = args.decay {
        update_decay(store_path, model_id, half_life, cfg)
            .map_err(|e| format!("update_decay: {e}"))?

    } else {
        return Err(
            "update mode requires one of: --add, --remove-source, --reweight, --decay".into()
        );
    };

    // ── Acceptance gate (Part D) ──────────────────────────────────────────────
    let commit = if let Some(ref validate_path) = args.validate.clone() {
        let database = args.database.clone()
            .ok_or("--database is required with --validate")?;

        eprintln!("train update: running acceptance gate on {} ...", validate_path.display());
        let val_spectra = load_spectra_for_train(validate_path)?;

        let store = ModelStore::open(store_path)
            .map_err(|e| format!("opening store for gate: {e}"))?;
        let current_param = store.load_param(model_id)
            .map_err(|e| format!("loading current model for gate: {e}"))?;
        let current_scorer = RankScorer::new(&current_param);
        let candidate_scorer = RankScorer::new(&candidate);

        let search_params = build_train_search_params(&args.mods)?;

        let delta = evaluate_candidate(
            &val_spectra,
            &database,
            &current_scorer,
            &candidate_scorer,
            &search_params,
            args.train_fdr,
        )
        .map_err(|e| format!("evaluate_candidate: {e}"))?;

        eprintln!(
            "train update: gate — current={} PSMs, candidate={} PSMs at FDR={}",
            delta.current_count, delta.candidate_count, args.train_fdr
        );

        if delta.is_accepted() {
            eprintln!("train update: ACCEPTED (candidate >= current)");
            true
        } else {
            eprintln!("train update: REJECTED (candidate < current)");
            if args.force {
                eprintln!("train update: --force set, committing anyway");
                true
            } else {
                eprintln!("train update: skipping commit (use --force to override)");
                false
            }
        }
    } else {
        eprintln!("train update: no --validate dataset; skipping acceptance gate");
        if args.force {
            eprintln!("train update: --force set, committing unconditionally");
        }
        // Without --validate, commit unless user explicitly uses --force to control.
        // Default: commit (no gate run = no evidence of regression).
        true
    };

    if commit {
        commit_update(store_path, model_id, &candidate, &new_sources)
            .map_err(|e| format!("commit_update: {e}"))?;
        eprintln!(
            "train update: committed model '{model_id}' to {} [{:.2}s]",
            store_path.display(),
            t0.elapsed().as_secs_f64(),
        );
    }

    Ok(())
}

/// Parse `"source-id=weight"` from a `--reweight` argument.
fn parse_reweight_spec(spec: &str) -> Result<(String, f32), Box<dyn std::error::Error>> {
    let pos = spec.rfind('=').ok_or_else(|| {
        format!("--reweight value must be <source-id>=<weight>, got '{spec}'")
    })?;
    let sid = spec[..pos].to_string();
    let weight: f32 = spec[pos + 1..].parse()
        .map_err(|e| format!("invalid weight in --reweight '{spec}': {e}"))?;
    Ok((sid, weight))
}

/// Load the seed Param from the optional seed model specifier.
fn load_seed_param(seed_model: &Option<String>) -> Result<(String, Param), Box<dyn std::error::Error>> {
    match seed_model {
        None => {
            let store_path = bundled_store_path();
            let store = ModelStore::open(&store_path)
                .map_err(|e| format!("opening bundled store: {e}"))?;
            let p = store.load_param("hcd_qexactive_tryp")
                .map_err(|e| format!("loading seed model: {e}"))?;
            Ok(("hcd_qexactive_tryp".to_string(), p))
        }
        Some(seed) => {
            // Seed by slug from the canonical Parquet store.
            let store_path = bundled_store_path();
            let store = ModelStore::open(&store_path)
                .map_err(|e| format!("opening bundled store: {e}"))?;
            let p = store.load_param(seed)
                .map_err(|e| format!("loading seed model '{seed}': {e}"))?;
            Ok((seed.clone(), p))
        }
    }
}

/// Build an `AminoAcidSet` from an optional mods file, defaulting to
/// Carbamidomethyl-C fixed + Oxidation-M variable + protein-N-term Acetyl variable.
fn build_aa_set(
    mods: &Option<PathBuf>,
) -> Result<model::AminoAcidSet, Box<dyn std::error::Error>> {
    match mods {
        Some(path) => {
            let set = AminoAcidSetBuilder::new_standard()
                .add_mods_from_file(path)
                .map_err(|e| format!("loading mods from {}: {e}", path.display()))?
                .build()
                .map_err(|e| format!("building amino-acid set: {e}"))?;
            Ok(set)
        }
        None => {
            let cam = Modification {
                name: "Carbamidomethyl".into(),
                mass_delta: 57.02146,
                residue: ResidueSpec::Specific(b'C'),
                location: ModLocation::Anywhere,
                fixed: true,
                accession: None,
                neutral_losses: Vec::new(),
                loss_class: 0,
            };
            let ox = Modification {
                name: "Oxidation".into(),
                mass_delta: 15.99491,
                residue: ResidueSpec::Specific(b'M'),
                location: ModLocation::Anywhere,
                fixed: false,
                accession: None,
                neutral_losses: Vec::new(),
                loss_class: 0,
            };
            Ok(AminoAcidSetBuilder::new_standard()
                .add_fixed_mod(cam)
                .add_variable_mod(ox)
                .build()?)
        }
    }
}

/// Build the default `AminoAcidSet` (CAM-C fixed + Ox-M variable), optionally
/// with an isobaric tag (TMT/iTRAQ) as a FIXED mod on K + peptide N-term.
///
/// Used by the no-`--mods` (parameter-free) path: when the protocol resolves to
/// TMT/iTRAQ the tag MUST be in the candidate set, or every labeled peptide is
/// `+tag·(nK+1)` Da off at the precursor and silently misses (C1).
fn default_aa_set_with_tag(
    tag: Option<(&str, f64)>,
) -> Result<model::AminoAcidSet, Box<dyn std::error::Error>> {
    let cam = Modification {
        name: "Carbamidomethyl".into(),
        mass_delta: 57.02146,
        residue: ResidueSpec::Specific(b'C'),
        location: ModLocation::Anywhere,
        fixed: true,
        accession: None,
        neutral_losses: Vec::new(),
        loss_class: 0,
    };
    let ox = Modification {
        name: "Oxidation".into(),
        mass_delta: 15.99491,
        residue: ResidueSpec::Specific(b'M'),
        location: ModLocation::Anywhere,
        fixed: false,
        accession: None,
        neutral_losses: Vec::new(),
        loss_class: 0,
    };
    // Protein N-terminal acetylation (+42.010565) — a near-universal default in
    // the field (the reference engine / a comparison engine / a quantitation tool). Restricted to the PROTEIN N-term
    // (one site per protein, after Met cleavage), so it is combinatorially cheap
    // (not every peptide N-term) and biologically correct.
    let acetyl = Modification {
        name: "Acetyl".into(),
        mass_delta: 42.010565,
        residue: ResidueSpec::Wildcard,
        location: ModLocation::ProtNTerm,
        fixed: false,
        accession: None,
        neutral_losses: Vec::new(),
        loss_class: 0,
    };
    let mut b = AminoAcidSetBuilder::new_standard()
        .add_fixed_mod(cam)
        .add_variable_mod(ox)
        .add_variable_mod(acetyl);
    if let Some((name, mass)) = tag {
        let tag_k = Modification {
            name: name.into(),
            mass_delta: mass,
            residue: ResidueSpec::Specific(b'K'),
            location: ModLocation::Anywhere,
            fixed: true,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
        };
        let tag_nterm = Modification {
            name: name.into(),
            mass_delta: mass,
            residue: ResidueSpec::Wildcard,
            location: ModLocation::NTerm,
            fixed: true,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
        };
        b = b.add_fixed_mod(tag_k).add_fixed_mod(tag_nterm);
    }
    Ok(b.build()?)
}

/// Format today's date as `YYYY-MM-DD` using `std::time::SystemTime`.
fn format_today_iso8601() -> String {
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    // Simple Gregorian calendar conversion from Unix timestamp (days since epoch).
    let days = secs / 86400;
    unix_days_to_iso8601(days)
}

fn unix_days_to_iso8601(days: u64) -> String {
    // Algorithm: Gregorian calendar from Julian Day Number.
    // JDN for 1970-01-01 = 2440588.
    let jdn = days as i64 + 2_440_588;
    let a = jdn + 32044;
    let b = (4 * a + 3) / 146097;
    let c = a - (146097 * b) / 4;
    let d = (4 * c + 3) / 1461;
    let e = c - (1461 * d) / 4;
    let m = (5 * e + 2) / 153;
    let day = e - (153 * m + 2) / 5 + 1;
    let month = m + 3 - 12 * (m / 10);
    let year = 100 * b + d - 4800 + m / 10;
    format!("{:04}-{:02}-{:02}", year, month, day)
}

/// Convert the CLI `Fragmentation` enum to `Option<ActivationMethod>`.
///
/// `Fragmentation::Auto` returns `None` (no activation explicitly requested);
/// every concrete variant maps to its `ActivationMethod`. Used by
/// [`resolve_metadataless_selection`] so that an unset `--fragmentation`
/// defers to detection or the class-consistent default.
fn cli_fragmentation_to_activation_opt(f: Fragmentation) -> Option<ActivationMethod> {
    match f {
        Fragmentation::Auto => None,
        Fragmentation::Cid  => Some(ActivationMethod::CID),
        Fragmentation::Etd  => Some(ActivationMethod::ETD),
        Fragmentation::Hcd  => Some(ActivationMethod::HCD),
        Fragmentation::Uvpd => Some(ActivationMethod::UVPD),
    }
}

/// Resolve the CLI fragment-tolerance override (MGF only) into a `Tolerance`.
/// `--fragment-tol-ppm` ⇒ `Ppm`; `--fragment-tol-da` ⇒ `Da`; none ⇒ `None`.
fn cli_fragment_tol_override(
    fragment_tol_ppm: Option<f64>,
    fragment_tol_da: Option<f64>,
) -> Option<model::tolerance::Tolerance> {
    use model::tolerance::Tolerance;
    fragment_tol_ppm
        .map(Tolerance::Ppm)
        .or_else(|| fragment_tol_da.map(Tolerance::Da))
}

/// Parse the `--enzyme` value into `(primary, extras)`.
///
/// Accepts `,` or `+` separators (the latter matching a quantitation tool/the reference engine
/// docs), trims whitespace, drops empty tokens, and de-duplicates by enzyme
/// while preserving order. The first surviving entry is the primary (drives
/// model selection and the cleavage-credit PIN feature); the rest only widen
/// candidate enumeration. Errors on an unknown name or empty input.
fn parse_enzymes(
    spec: &str,
) -> Result<(model::enzyme::Enzyme, Vec<model::enzyme::Enzyme>), Box<dyn std::error::Error>> {
    use model::enzyme::Enzyme;
    let mut all: Vec<Enzyme> = Vec::new();
    for tok in spec.split([',', '+']).map(str::trim).filter(|s| !s.is_empty()) {
        let e = Enzyme::from_name(tok).ok_or_else(|| format!(
            "unknown --enzyme '{tok}' (expected trypsin/chymotrypsin/lysc/aspn/gluc/lysn/argc/\
             alphalp/nocleavage/nonspecific/elastase)"
        ))?;
        if !all.contains(&e) {
            all.push(e); // dedup, keep order (first = primary)
        }
    }
    let (primary, extras) = all.split_first().ok_or("empty --enzyme")?;
    Ok((*primary, extras.to_vec()))
}

/// Warn when a universal protease (NonSpecific, or AlphaLP which cleaves every
/// bond) is combined with specific protease(s): the cleavage-site union then
/// covers EVERY position, so the digest is effectively non-specific — usually
/// not what a user adding it to a specific list intends.
fn warn_if_universal_protease_combo(
    primary: model::enzyme::Enzyme,
    extras: &[model::enzyme::Enzyme],
) {
    use model::enzyme::Enzyme;
    let is_universal = |e: &Enzyme| matches!(e, Enzyme::NonSpecific | Enzyme::AlphaLP);
    let is_specific =
        |e: &Enzyme| !matches!(e, Enzyme::NonSpecific | Enzyme::AlphaLP | Enzyme::NoCleavage);
    let all: Vec<Enzyme> = std::iter::once(primary).chain(extras.iter().copied()).collect();
    if all.len() > 1 && all.iter().any(is_universal) && all.iter().any(is_specific) {
        eprintln!(
            "WARN: --enzyme combines a universal protease (NonSpecific/AlphaLP) with specific \
             one(s); their union cleaves at EVERY position, so this is effectively a \
             non-specific digest (much larger search space). Drop the universal protease to \
             keep the specific cleavage rules."
        );
    }
}

/// Resolve (activation, instrument) for model selection on metadata-less input
/// (MGF, or mzML/.raw with no analyzer metadata). Resolution class comes from
/// the `--fragment-tol-*` unit; activation from detected method, else
/// `--fragmentation`, else the class-consistent default. When nothing
/// disambiguates, defaults to CID / LowRes (→ `cid_lowres_tryp`) + a warning.
fn resolve_metadataless_selection(
    detected_activation: Option<ActivationMethod>,
    fragmentation: Fragmentation,
    fragment_tol_ppm: Option<f64>,
    fragment_tol_da: Option<f64>,
) -> (ActivationMethod, Option<InstrumentType>) {
    let instrument: Option<InstrumentType> = if fragment_tol_ppm.is_some() {
        Some(InstrumentType::QExactive)
    } else if fragment_tol_da.is_some() {
        Some(InstrumentType::LowRes)
    } else {
        None
    };
    let explicit = cli_fragmentation_to_activation_opt(fragmentation);
    // Class-consistent default when neither detection nor `--fragmentation`
    // names an activation: high-res classes imply HCD, otherwise CID.
    let class_default = match instrument {
        Some(InstrumentType::QExactive)
        | Some(InstrumentType::HighRes)
        | Some(InstrumentType::TOF) => ActivationMethod::HCD,
        _ => ActivationMethod::CID,
    };
    let activation = detected_activation.or(explicit).unwrap_or(class_default);
    if detected_activation.is_none() && explicit.is_none() && instrument.is_none() {
        eprintln!(
            "WARN: MGF input with no --fragmentation/--fragment-tol; assuming \
             CID / low-res / 0.5 Da. Pass --fragmentation and --fragment-tol-ppm/-da \
             to override."
        );
    }
    (activation, instrument)
}

/// Peek the spectrum file and return the dominant
/// `ActivationMethod` across the first several MS2 spectra.
///
/// Reads up to `MAX_PEEK` spectra (early-exit) and tallies a histogram of
/// activation methods. Returns the most-common method, or `None` when no
/// spectra carry an activation cvParam (older mzMLs, MGF, etc.).
///
/// Currently only mzML files (`.mzml` / `.mzML` extension) carry an
/// `<activation>` block. For anything else (MGF, unknown extension) we
/// return `None` and the caller falls back to the historical default.
///
/// When multiple activation methods are present, prints a single
/// `eprintln!` warning naming the runner-up and its count.
fn detect_dominant_activation(spectrum_path: &std::path::Path) -> Option<ActivationMethod> {
    // Only mzML carries `<activation>`. Other formats: caller falls back.
    let ext_lower = spectrum_ext_lower(spectrum_path);
    if ext_lower.as_deref() != Some("mzml") {
        return None;
    }

    const MAX_PEEK: usize = 64;

    let reader = MzMLReader::new(input::open_buf_maybe_gz(spectrum_path).ok()?);

    // Tally counts keyed by ActivationMethod variant.
    let mut counts: std::collections::HashMap<ActivationMethod, usize> =
        std::collections::HashMap::new();
    for (seen, item) in reader.enumerate() {
        if seen >= MAX_PEEK {
            break;
        }
        if let Ok(spec) = item {
            if let Some(m) = spec.activation_method {
                *counts.entry(m).or_insert(0) += 1;
            }
        }
    }

    if counts.is_empty() {
        return None;
    }

    // Find the dominant method. Ties are broken by ActivationMethod's
    // declaration order via match below, which is stable.
    let dominant = counts
        .iter()
        // Deterministic on ties: HashMap iteration order is randomised per
        // process, so a bare `max_by_key(count)` picks an ARBITRARY maximum. A
        // 1:1 interleaved HCD/ETD acquisition ties exactly at 32/32 over the
        // peeked window and flipped the selected model run-to-run on identical
        // input. Tie-break on the discriminant so the choice is reproducible.
        .max_by_key(|(&m, &n)| (n, std::cmp::Reverse(m as u8)))
        .map(|(&m, _)| m)?;

    // Warn on mixed activation. The dominant method still wins; this is
    // purely informational so the user can spot heterogeneous mzMLs.
    if counts.len() > 1 {
        let mut other_pairs: Vec<(ActivationMethod, usize)> = counts
            .iter()
            .filter(|(&m, _)| m != dominant)
            .map(|(&m, &n)| (m, n))
            .collect();
        other_pairs.sort_by(|a, b| b.1.cmp(&a.1));
        let total: usize = counts.values().sum();
        let dominant_count = counts[&dominant];
        eprintln!(
            "Param resolver: mixed activation methods in input ({} different methods \
             across {} peeked MS2 spectra). Using dominant = {} ({}/{}); other methods \
             present: {}",
            counts.len(),
            total,
            dominant.name(),
            dominant_count,
            total,
            other_pairs
                .iter()
                .map(|(m, n)| format!("{}={}", m.name(), n))
                .collect::<Vec<_>>()
                .join(", "),
        );
    }

    Some(dominant)
}

/// Helper to call `input::detect_instrument_type` on an mzML path.
///
/// Mirrors the structure of `detect_dominant_activation` so the two
/// detection passes look symmetric at the call site. Returns `None` for
/// non-mzML inputs or when the mzML has no recoverable instrument metadata.
fn detect_instrument_type_for_path(spectrum_path: &std::path::Path) -> Option<InstrumentType> {
    let ext_lower = spectrum_ext_lower(spectrum_path);
    if ext_lower.as_deref() != Some("mzml") {
        return None;
    }

    detect_instrument_type(input::open_buf_maybe_gz(spectrum_path).ok()?)
}

/// Resolve the path to the bundled model store.
///
/// The store may ship either as a single `resources/models.parquet` file or as
/// a per-protocol partitioned directory `resources/models/` (Hive-style
/// `protocol=<P>/models.parquet`). [`ModelStore::open`] accepts both, so the
/// resolver prefers the partitioned **directory** when present and falls back
/// to the single file.
///
/// A packaged release ships `resources/` next to the binary, so prefer
/// `<exe_dir>/resources/...` when it exists — that makes an installed binary
/// self-contained regardless of where it runs. Fall back to the compile-time
/// source tree (`CARGO_MANIFEST_DIR`) for `cargo run` / tests.
fn bundled_store_path() -> PathBuf {
    // Search roots in priority order: next to the binary, then the source tree.
    let mut roots: Vec<PathBuf> = Vec::new();
    if let Ok(exe) = std::env::current_exe() {
        if let Some(dir) = exe.parent() {
            roots.push(dir.join("resources"));
        }
    }
    roots.push(PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../resources"));

    for root in &roots {
        let partitioned = root.join("models");
        if partitioned.is_dir() {
            return partitioned;
        }
        let single = root.join("models.parquet");
        if single.exists() {
            return single;
        }
    }

    // Last-resort default (source tree single file) for error messages.
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../resources/models.parquet")
}

/// Build a [`SelectionKey`] from `(activation, instrument, protocol)` applying
/// all old-ladder normalizations. This is the new entry point used by the
/// search binary, replacing the former filename-based resolution ladder.
///
/// `activation`: the detected or explicitly set `ActivationMethod`.
/// `instrument`: the detected or explicitly set `InstrumentType` (None = undetected → LowRes).
/// `protocol`:   the CLI `Protocol` value.
fn build_selection_key(
    activation: ActivationMethod,
    instrument: Option<InstrumentType>,
    protocol: Protocol,
    enzyme: model::enzyme::Enzyme,
) -> SelectionKey {
    use std::collections::BTreeSet;

    // 1. PQD → CID for model routing.
    let act_str: &str = match activation {
        ActivationMethod::PQD => "CID",
        other                 => other.name(),
    };
    // 2. Apply family fallback (OrbitrapAstral → QExactive, TimsTOF → TOF).
    let inst_after_family: &str = match instrument {
        Some(i) => i.family_fallback().name(),
        None    => "LowRes",
    };

    // 3. Apply old-ladder (activation, instrument) normalization.
    //    Because `normalize_for_store` returns `&'static str` only for the
    //    normalizing arms (to avoid lifetime issues), we handle identity
    //    inline here.
    let (final_act, final_inst, drop_protocol): (&str, &str, bool) =
        match (act_str, inst_after_family) {
            // H5: low-res (ion-trap) HCD. Routing this to the high-res
            // QExactive model would match 0.5-Da peaks at 20 ppm and lose
            // ~18% of PSMs silently. No hcd_lowres model exists, so route to the
            // low-res b/y model (cid_lowres_tryp) instead — correct fragment
            // tolerance and ion series. Pinned by store_selection_equivalence.rs.
            ("HCD", "LowRes")    => {
                eprintln!(
                    "WARN: low-res (ion-trap) HCD detected — no hcd_lowres model exists; \
                     routing to cid_lowres_tryp (low-res b/y, 0.5-Da tolerance) rather than \
                     the high-res QExactive model. Pass --model to override."
                );
                ("CID", "LowRes", true)
            }
            ("HCD", "TOF")       => ("CID", "TOF",       true),
            ("CID", "QExactive") => ("CID", "LowRes",    true),
            ("ETD", i) if !matches!(i, "LowRes" | "HighRes") => ("ETD", "LowRes", true),
            ("UVPD", i) if i != "QExactive" => ("CID", "LowRes", true),
            _ => (act_str, inst_after_family, false),
        };

    // 4. Build experiment_class from protocol (unless the final fallback dropped it).
    //    Protocol → experiment_class mapping matches the parquet's `protocol` column.
    let protocol_for_store: &str = match protocol {
        Protocol::Auto | Protocol::Standard => "Automatic",
        Protocol::Tmt          => "TMT",
        Protocol::Phospho      => "Phosphorylation",
        Protocol::Itraq        => "iTRAQ",
        Protocol::ItraqPhospho => "iTRAQPhospho",
    };
    let experiment_class: BTreeSet<String> = if drop_protocol {
        BTreeSet::new()
    } else {
        store_protocol_to_experiment_class(protocol_for_store)
    };

    SelectionKey {
        activation: final_act.to_string(),
        instrument: final_inst.to_string(),
        // Parquet stores the enzyme as its `Enzyme::name()` ("Trypsin", "LysC", ...).
        enzyme: enzyme.name().to_string(),
        experiment_class,
    }
}

/// Load the scoring [`Param`] from the bundled Parquet store for the given
/// `(activation, instrument, protocol)` combination.
///
/// The `model_id` selected from the store is guaranteed to match the
/// reference resolution ladder by the equivalence gate test
/// `store_selection_matches_old_ladder_for_all_combos`.
///
/// `custom_store_path`: when `Some`, use that Parquet file instead of the
/// bundled `resources/models.parquet` (honours `--model-store`).
///
/// `model_id_override`: when `Some`, skip automatic selection and load this
/// exact model ID (honours `--model`).
fn load_param_from_store(
    activation: ActivationMethod,
    instrument: Option<InstrumentType>,
    protocol: Protocol,
    enzyme: model::enzyme::Enzyme,
    custom_store_path: Option<&Path>,
    model_id_override: Option<&str>,
) -> Result<(String, Param), Box<dyn std::error::Error>> {
    let store_path = custom_store_path
        .map(|p| p.to_owned())
        .unwrap_or_else(bundled_store_path);
    let store = ModelStore::open(&store_path)
        .map_err(|e| format!("opening model store {}: {e}", store_path.display()))?;

    let model_id: String = if let Some(id) = model_id_override {
        id.to_string()
    } else {
        let entries = store.selection_entries();
        let key = build_selection_key(activation, instrument, protocol, enzyme);

        // Forward-compat: `build_selection_key` collapses instruments with a real
        // family fallback (OrbitrapAstral → QExactive, TimsTOF → TOF) so the
        // bundled models resolve correctly. But that also hides a model trained
        // for the EXACT instrument (e.g. a user-trained OrbitrapAstral model).
        // Try the exact detected instrument FIRST; only when no such model exists
        // (the bundled case) do we fall through to the normalized family ladder —
        // so bundled selection (and the equivalence gate) is unchanged.
        let exact_id: Option<String> = match instrument {
            Some(i) if i.family_fallback().name() != i.name() => {
                let raw_key = SelectionKey {
                    instrument: i.name().to_string(),
                    ..key.clone()
                };
                select(&entries, &raw_key, |s| s.to_string(), None).map(|s| s.to_string())
            }
            _ => None,
        };

        exact_id.unwrap_or_else(|| {
            // `build_selection_key` already applies family fallback + all
            // normalizations, so the family_fn here is the identity. When the
            // exact ladder misses (e.g. a protocol the own-only store doesn't
            // carry), select_nearest routes to the CLOSEST own model — relaxing
            // the enzyme (keeping the activation and instrument) — and only
            // resolves to the standard base as a last resort, WARNing which model
            // it substituted so the user can pin one with --model.
            let (id, substituted) =
                select_nearest(&entries, &key, |i| i.to_string(), "hcd_qexactive_tryp");
            if substituted {
                eprintln!(
                    "WARN: no model matched (activation={}, instrument={}, enzyme={}, class={:?}) \
                     — using the nearest available model '{}'; scores may be mis-calibrated for \
                     this data. Pin a model with --model if this is wrong.",
                    key.activation, key.instrument, key.enzyme, key.experiment_class, id
                );
            }
            id.to_string()
        })
    };

    let param = store.load_param(&model_id)
        .map_err(|e| format!("loading model '{model_id}' from store: {e}"))?;

    Ok((model_id, param))
}

/// Parse `--fragmentation` value. Accepts named values (case-insensitive: auto,
/// CID, ETD, HCD, UVPD).
fn parse_fragmentation(s: &str) -> Result<Fragmentation, String> {
    <Fragmentation as ValueEnum>::from_str(s, true).map_err(|_| {
        format!("invalid fragmentation `{s}`: expected auto|CID|ETD|HCD|UVPD")
    })
}

/// Parse `--protocol` value. Accepts named values only.
fn parse_protocol(s: &str) -> Result<Protocol, String> {
    <Protocol as ValueEnum>::from_str(s, true).map_err(|_| {
        format!(
            "invalid --protocol `{s}`: expected \
             auto|phospho|iTRAQ|iTRAQ-phospho|TMT|standard"
        )
    })
}

/// In a `MIN-MAX` dash range, the separator is the first `-` whose preceding
/// char is an ASCII digit: any earlier `-` is MIN's sign and any `-` right after
/// it is MAX's sign. This makes signed endpoints unambiguous — `2-5`→(2,5),
/// `-1-2`→(-1,2), `-3--1`→(-3,-1) — which a naive `split`/`rsplit` on `-` cannot.
fn dash_sep_index(s: &str) -> Option<usize> {
    let b = s.as_bytes();
    (1..b.len()).find(|&i| b[i] == b'-' && b[i - 1].is_ascii_digit())
}

/// Parse a `MIN..MAX` (or `MIN-MAX`) range into a `(min, max)` pair, generic
/// over the integer type so it serves both `--charge` (u8) and
/// `--isotope-error` (i8, negatives allowed). The `-` separator is tried only
/// when the value does not parse as `..`; signed endpoints are supported in both
/// forms (`-1..2`, `-3--1`).
fn parse_int_range<T>(s: &str, label: &str) -> Result<(T, T), String>
where
    T: std::str::FromStr + PartialOrd + std::fmt::Display + Copy,
{
    let trimmed = s.trim();
    let (lo_s, hi_s) = if let Some((a, b)) = trimmed.split_once("..") {
        (a.trim(), b.trim())
    } else if let Some(idx) = dash_sep_index(trimmed) {
        (trimmed[..idx].trim(), trimmed[idx + 1..].trim())
    } else {
        return Err(format!("invalid {label} `{s}`: expected MIN..MAX (or MIN-MAX)"));
    };
    let lo: T = lo_s
        .parse()
        .map_err(|_| format!("invalid {label} `{s}`: bad MIN `{lo_s}`"))?;
    let hi: T = hi_s
        .parse()
        .map_err(|_| format!("invalid {label} `{s}`: bad MAX `{hi_s}`"))?;
    if lo > hi {
        return Err(format!("invalid {label} `{s}`: MIN {lo} > MAX {hi}"));
    }
    Ok((lo, hi))
}

/// Largest precursor charge state andes will search. Beyond this the candidate
/// space is dominated by noise and the fragment model is untrained.
const MAX_SUPPORTED_CHARGE: u8 = 50;

/// Parse `--charge MIN..MAX` (also `MIN-MAX`) into `(u8, u8)`.
///
/// Domain-validates (finding 3.8): the minimum charge must be >= 1 (charge 0 is
/// not a real precursor) and the maximum must be within the supported bound.
fn parse_charge_range(s: &str) -> Result<(u8, u8), String> {
    let (lo, hi) = parse_int_range::<u8>(s, "charge")?;
    if lo < 1 {
        return Err(format!(
            "invalid charge `{s}`: minimum charge must be >= 1 (got {lo})"
        ));
    }
    if hi > MAX_SUPPORTED_CHARGE {
        return Err(format!(
            "invalid charge `{s}`: maximum charge {hi} exceeds supported bound {MAX_SUPPORTED_CHARGE}"
        ));
    }
    Ok((lo, hi))
}

/// Parse `--isotope-error MIN..MAX` (also `MIN-MAX`) into `(i8, i8)`; negatives allowed.
fn parse_isotope_error_range(s: &str) -> Result<(i8, i8), String> {
    parse_int_range::<i8>(s, "isotope-error")
}

/// Parse `--precursor-tol VALUE+unit` (e.g. `20ppm`, `0.02da`/`0.02Da`).
fn parse_precursor_tol(s: &str) -> Result<Tolerance, String> {
    let t = s.trim();
    let lower = t.to_ascii_lowercase();
    let (num_str, is_ppm) = if let Some(n) = lower.strip_suffix("ppm") {
        (n, true)
    } else if let Some(n) = lower.strip_suffix("da") {
        (n, false)
    } else {
        return Err(format!(
            "invalid --precursor-tol `{s}`: expected VALUE+unit, e.g. 20ppm or 0.02da"
        ));
    };
    let v: f64 = num_str
        .trim()
        .parse()
        .map_err(|_| format!("invalid --precursor-tol `{s}`: bad number `{num_str}`"))?;
    // Domain check (finding 3.8): a tolerance must be a finite positive width.
    if !v.is_finite() || v <= 0.0 {
        return Err(format!(
            "invalid --precursor-tol `{s}`: value must be finite and > 0 (got {v})"
        ));
    }
    Ok(if is_ppm { Tolerance::Ppm(v) } else { Tolerance::Da(v) })
}

/// f32 companion to [`parse_unit_fraction`], for CLI fractions stored as `f32`.
/// Rejects NaN, negatives and values above 1 at PARSE time rather than letting a nonsense
/// threshold silently disable or invert a gate.
fn parse_unit_fraction_f32(s: &str) -> Result<f32, String> {
    parse_unit_fraction(s).map(|v| v as f32)
}

/// Parse a probability-domain CLI value (FDR / PEP / refine-FDR) — must be a
/// finite number in `[0, 1]` (finding 3.8). Used as a clap `value_parser` so a
/// bad value is rejected at parse time with a clear message.
fn parse_unit_fraction(s: &str) -> Result<f64, String> {
    let v: f64 = s
        .trim()
        .parse()
        .map_err(|_| format!("invalid value `{s}`: expected a number in [0, 1]"))?;
    if !v.is_finite() || !(0.0..=1.0).contains(&v) {
        return Err(format!(
            "invalid value `{s}`: must be finite and within [0, 1] (got {v})"
        ));
    }
    Ok(v)
}

/// Parse a fragment-tolerance scalar (`--fragment-tol-ppm` / `--fragment-tol-da`)
/// — must be a finite, strictly-positive width (finding 3.8).
fn parse_positive_tol(s: &str) -> Result<f64, String> {
    let v: f64 = s
        .trim()
        .parse()
        .map_err(|_| format!("invalid tolerance `{s}`: expected a positive number"))?;
    if !v.is_finite() || v <= 0.0 {
        return Err(format!(
            "invalid tolerance `{s}`: must be finite and > 0 (got {v})"
        ));
    }
    Ok(v)
}

/// Parse `--precursor-cal` value. Accepts auto|on|off.
fn parse_precursor_cal(s: &str) -> Result<PrecursorCalMode, String> {
    match s.to_ascii_lowercase().as_str() {
        "auto" => Ok(PrecursorCalMode::Auto),
        "on" => Ok(PrecursorCalMode::On),
        "off" => Ok(PrecursorCalMode::Off),
        _ => Err(format!(
            "invalid precursor-cal `{s}`: expected auto|on|off (Java -precursorCal)"
        )),
    }
}

fn parse_enzyme_specificity(s: &str) -> Result<EnzymeSpecificity, String> {
    <EnzymeSpecificity as ValueEnum>::from_str(s, true).map_err(|_| {
        format!("invalid enzyme specificity `{s}`: expected non-specific|semi|fully")
    })
}

#[cfg(test)]
mod cli_domain_validator_tests {
    use super::*;

    #[test]
    fn charge_range_rejects_zero_and_out_of_bounds() {
        assert!(parse_charge_range("2..5").is_ok());
        assert!(parse_charge_range("1..1").is_ok());
        assert!(parse_charge_range("0..5").is_err(), "charge 0 must be rejected");
        assert!(parse_charge_range("0..0").is_err(), "charge 0..0 must be rejected");
        assert!(parse_charge_range("2..60").is_err(), "charge above bound must be rejected");
    }

    #[test]
    fn precursor_tol_rejects_nonpositive_and_nonfinite() {
        assert!(parse_precursor_tol("20ppm").is_ok());
        assert!(parse_precursor_tol("0.02da").is_ok());
        assert!(parse_precursor_tol("0ppm").is_err(), "zero tol rejected");
        assert!(parse_precursor_tol("-5ppm").is_err(), "negative tol rejected");
        assert!(parse_precursor_tol("infda").is_err(), "inf tol rejected");
        assert!(parse_precursor_tol("nanppm").is_err(), "nan tol rejected");
    }

    #[test]
    fn unit_fraction_rejects_out_of_range() {
        assert_eq!(parse_unit_fraction("0.01").unwrap(), 0.01);
        assert_eq!(parse_unit_fraction("0").unwrap(), 0.0);
        assert_eq!(parse_unit_fraction("1").unwrap(), 1.0);
        assert!(parse_unit_fraction("-0.1").is_err());
        assert!(parse_unit_fraction("1.5").is_err());
        assert!(parse_unit_fraction("nan").is_err());
        assert!(parse_unit_fraction("inf").is_err());
    }

    #[test]
    fn positive_tol_rejects_nonpositive_and_nonfinite() {
        assert_eq!(parse_positive_tol("20").unwrap(), 20.0);
        assert!(parse_positive_tol("0").is_err());
        assert!(parse_positive_tol("-1").is_err());
        assert!(parse_positive_tol("nan").is_err());
        assert!(parse_positive_tol("inf").is_err());
    }
}

#[cfg(test)]
mod enzyme_cli_tests {
    use super::parse_enzymes;
    use model::enzyme::Enzyme;

    #[test]
    fn single_enzyme_has_no_extras() {
        let (primary, extras) = parse_enzymes("trypsin").unwrap();
        assert_eq!(primary, Enzyme::Trypsin);
        assert!(extras.is_empty(), "single enzyme ⇒ empty extras (bit-identical path)");
    }

    #[test]
    fn comma_and_plus_separators_both_parse() {
        let (p1, e1) = parse_enzymes("gluc,trypsin").unwrap();
        let (p2, e2) = parse_enzymes("gluc+trypsin").unwrap();
        assert_eq!((p1, e1.as_slice()), (Enzyme::GluC, [Enzyme::Trypsin].as_slice()));
        assert_eq!((p2, e2.as_slice()), (Enzyme::GluC, [Enzyme::Trypsin].as_slice()));
    }

    #[test]
    fn first_token_is_primary_order_matters_for_selection() {
        assert_eq!(parse_enzymes("trypsin,gluc").unwrap().0, Enzyme::Trypsin);
        assert_eq!(parse_enzymes("gluc,trypsin").unwrap().0, Enzyme::GluC);
    }

    #[test]
    fn duplicate_tokens_are_deduped_preserving_order() {
        // "trypsin,trypsin,gluc" and the alias "tryp" collapse to [Trypsin, GluC].
        let (primary, extras) = parse_enzymes("trypsin, tryp , gluc, gluc").unwrap();
        assert_eq!(primary, Enzyme::Trypsin);
        assert_eq!(extras, vec![Enzyme::GluC]);
    }

    #[test]
    fn whitespace_and_empties_tolerated() {
        let (primary, extras) = parse_enzymes("  lysc , , trypsin ").unwrap();
        assert_eq!(primary, Enzyme::LysC);
        assert_eq!(extras, vec![Enzyme::Trypsin]);
    }

    #[test]
    fn unknown_enzyme_errors_and_empty_errors() {
        assert!(parse_enzymes("bogus").is_err());
        assert!(parse_enzymes("").is_err());
        assert!(parse_enzymes("  , ").is_err());
    }

    #[test]
    fn elastase_aliases_to_nonspecific() {
        assert_eq!(parse_enzymes("elastase").unwrap().0, Enzyme::NonSpecific);
    }
}

#[cfg(test)]
mod param_resolver_tests {
    use super::*;

    // ── Model resolution is store-based: all bundled models live in
    //    resources/models.parquet and are selected by `model_id`.
    //    The store_selection_equivalence integration test covers the
    //    activation/instrument/protocol → model_id selection invariant.

    #[test]
    fn parse_fragmentation_accepts_named_rejects_numeric() {
        assert_eq!(parse_fragmentation("HCD").unwrap(), Fragmentation::Hcd);
        assert_eq!(parse_fragmentation("auto").unwrap(), Fragmentation::Auto);
        // Legacy numeric forms are no longer accepted.
        assert!(parse_fragmentation("3").is_err());
        assert!(parse_fragmentation("99").is_err());
    }

    #[test]
    fn parse_protocol_accepts_named_rejects_numeric() {
        assert_eq!(parse_protocol("TMT").unwrap(), Protocol::Tmt);
        assert_eq!(parse_protocol("auto").unwrap(), Protocol::Auto);
        assert!(parse_protocol("4").is_err());
        assert!(parse_protocol("99").is_err());
    }

    #[test]
    fn parse_enzyme_specificity_named_only() {
        assert_eq!(parse_enzyme_specificity("fully").unwrap(), EnzymeSpecificity::Fully);
        assert_eq!(parse_enzyme_specificity("semi").unwrap(), EnzymeSpecificity::Semi);
        assert!(parse_enzyme_specificity("2").is_err());
    }

    #[test]
    fn parse_charge_range_forms() {
        assert_eq!(parse_charge_range("2..5").unwrap(), (2, 5));
        assert_eq!(parse_charge_range("2-5").unwrap(), (2, 5));
        assert_eq!(parse_charge_range("3..3").unwrap(), (3, 3));
        assert!(parse_charge_range("5..2").is_err());
        assert!(parse_charge_range("x..5").is_err());
    }

    #[test]
    fn parse_isotope_error_range_allows_negatives() {
        assert_eq!(parse_isotope_error_range("-1..2").unwrap(), (-1, 2));
        assert_eq!(parse_isotope_error_range("-1-2").unwrap(), (-1, 2));
        assert_eq!(parse_isotope_error_range("0..0").unwrap(), (0, 0));
        assert!(parse_isotope_error_range("2..-1").is_err());
    }

    #[test]
    fn parse_precursor_tol_units() {
        assert_eq!(parse_precursor_tol("20ppm").unwrap(), Tolerance::Ppm(20.0));
        assert_eq!(parse_precursor_tol("0.02da").unwrap(), Tolerance::Da(0.02));
        assert_eq!(parse_precursor_tol("0.02Da").unwrap(), Tolerance::Da(0.02));
        assert!(parse_precursor_tol("20").is_err());
        assert!(parse_precursor_tol("xppm").is_err());
    }

    #[test]
    fn parse_precursor_cal_accepts_named_modes() {
        assert_eq!(parse_precursor_cal("auto").unwrap(), PrecursorCalMode::Auto);
        assert_eq!(parse_precursor_cal("OFF").unwrap(), PrecursorCalMode::Off);
        assert_eq!(parse_precursor_cal("on").unwrap(), PrecursorCalMode::On);
        assert!(parse_precursor_cal("bogus").is_err());
    }

    #[test]
    fn finalize_intensity_stats_drops_zero_count() {
        // A count==0 key (e.g. from a partial aggregation parquet) must not
        // produce NaN mean/var in the finalized model — it carries no signal.
        assert_eq!(finalize_intensity_stats(0.0, 0.0, 0), None);
    }

    #[test]
    fn finalize_intensity_stats_computes_mean_and_clamped_var() {
        // sum=6, sum_sq=14, count=2 -> mean=3, var=max(0, 14/2 - 9)= max(0,-2)=0.
        let (mean, var) = finalize_intensity_stats(6.0, 14.0, 2).unwrap();
        assert_eq!(mean, 3.0);
        assert_eq!(var, 0.0);
    }

    #[test]
    fn single_file_search_does_not_prefix_titles() {
        // Regression (default-path PIN parity): with one input file the
        // SpecId/title must stay unprefixed so PIN output is byte-identical to
        // the single-file behavior that predates multi-`--spectrum` support.
        assert_eq!(title_prefix_for(1, "myfile"), None);
    }

    #[test]
    fn multi_file_search_prefixes_titles_with_file_stem() {
        assert_eq!(title_prefix_for(2, "myfile").as_deref(), Some("myfile/"));
    }
}

#[cfg(test)]
mod train_from_msnet_tests {
    use super::*;
    use model::mass::H2O;

    /// Reference: unmodified peptide mass = sum of residue masses + H2O.
    fn unmod_mass(seq: &[u8]) -> f64 {
        let rsum: f64 = seq
            .iter()
            .map(|&r| model::AminoAcid::standard(r).unwrap().mass)
            .sum();
        rsum + H2O
    }

    #[test]
    fn unmodified_peptide_mass_matches_reference() {
        let p = build_msnet_peptide("PEPTIDEK", &[], &[], 0.0, 0.0).unwrap();
        let expected = unmod_mass(b"PEPTIDEK");
        assert_eq!(p.mass().to_bits(), expected.to_bits());
        assert!(p.residues.iter().all(|aa| !aa.is_modified()));
    }

    /// Oxidation on residue 1 of "MPEPTIDE" must add exactly +15.994915 Da to
    /// the unmodified mass — verifying mods are actually applied, not dropped.
    #[test]
    fn residue_mod_mass_is_correct() {
        const OX: f64 = 15.994915;
        let p = build_msnet_peptide("MPEPTIDE", &[1], &[OX], 0.0, 0.0).unwrap();
        let expected = unmod_mass(b"MPEPTIDE") + OX;
        // Exact: Peptide::new sums aa.mass + mod_.mass_delta, so the delta is
        // added with the same arithmetic as the reference.
        assert_eq!(p.mass().to_bits(), expected.to_bits());
        assert!(p.residues[0].is_modified(), "residue 1 should carry the mod");
        assert_eq!(p.residues[0].mod_.as_ref().unwrap().mass_delta, OX);
    }

    /// N-terminal Acetyl folds onto the first residue's mass-delta.
    #[test]
    fn nterm_mod_mass_is_correct() {
        const ACETYL: f64 = 42.010565;
        let p = build_msnet_peptide("PEPTIDEK", &[], &[], ACETYL, 0.0).unwrap();
        let expected = unmod_mass(b"PEPTIDEK") + ACETYL;
        assert_eq!(p.mass().to_bits(), expected.to_bits());
        assert!(p.residues[0].is_modified());
    }

    /// A residue mod and a terminal mod on the SAME residue must sum (not
    /// clobber) so the total mass stays correct.
    #[test]
    fn overlapping_residue_and_nterm_mods_sum() {
        const OX: f64 = 15.994915;
        const ACETYL: f64 = 42.010565;
        // Oxidation on residue 1 (M) AND N-term acetyl on the same first residue.
        let p = build_msnet_peptide("MPEPTIDE", &[1], &[OX], ACETYL, 0.0).unwrap();
        let expected = unmod_mass(b"MPEPTIDE") + OX + ACETYL;
        assert_eq!(p.mass().to_bits(), expected.to_bits());
        assert_eq!(p.residues[0].mod_.as_ref().unwrap().mass_delta, OX + ACETYL);
    }

    #[test]
    fn mismatched_mod_arrays_error() {
        let r = build_msnet_peptide("PEPTIDE", &[1, 2], &[1.0], 0.0, 0.0);
        assert!(r.is_err());
    }

    #[test]
    fn nonstandard_residue_errors() {
        let r = build_msnet_peptide("PEPTBDE", &[], &[], 0.0, 0.0);
        assert!(r.is_err());
    }
}

#[cfg(test)]
mod intensity_marginal_tests {
    use super::*;
    use rustc_hash::FxHashMap;

    fn cell(ion: &str, fn_: &str, fc: &str, nce: &str, count: i64, sum: f64) -> (IntensityAggKey, IntensityAggStats) {
        (
            IntensityAggKey {
                ion_type: ion.to_string(),
                flank_n: fn_.to_string(),
                flank_c: fc.to_string(),
                pos_bin: 5,
                charge: 2,
                nce_bin: nce.to_string(),
            },
            IntensityAggStats { count, sum_log_rel: sum, sum_log_rel_sq: sum * sum / count as f64 },
        )
    }

    /// The finalizer must emit the backoff marginals that `predict_log_rel`
    /// probes, or the trained per-context intensities go unused at inference
    /// (which passes `nce_bin="unknown"`).
    #[test]
    fn backoff_marginals_are_emitted_with_summed_stats() {
        let mut merged: FxHashMap<IntensityAggKey, IntensityAggStats> = FxHashMap::default();
        // Two real cells: same flank/pos/charge, different numeric nce bins.
        for (k, s) in [
            cell("y", "K", "R", "30", 100, -20.0),
            cell("y", "K", "R", "40", 50, -15.0),
        ] {
            merged.insert(k, s);
        }
        let exact = merged.len();
        add_backoff_marginals(&mut merged);

        // nce-marginal (`__any__`) keeps flanks, sums across both nce bins.
        let any_nce = merged
            .get(&cell("y", "K", "R", ANY_NCE, 0, 0.0).0)
            .expect("nce-marginal must exist");
        assert_eq!(any_nce.count, 150, "summed over the two nce bins");
        assert!((any_nce.sum_log_rel - (-35.0)).abs() < 1e-9);

        // flank-marginal (`*`/`*`) keeps each nce bin separately.
        assert!(merged.contains_key(&cell("y", ANY_FLANK, ANY_FLANK, "30", 0, 0.0).0));
        // both-marginal exists too.
        assert!(merged.contains_key(&cell("y", ANY_FLANK, ANY_FLANK, ANY_NCE, 0, 0.0).0));

        // Real cells are untouched (no double counting).
        assert_eq!(merged.get(&cell("y", "K", "R", "30", 0, 0.0).0).unwrap().count, 100);
        assert!(merged.len() > exact, "marginals added new keys");
    }
}
