//! Command-line argument definitions and value parsers for the search subcommand.

use std::path::PathBuf;

use clap::{Args, ValueEnum};
use model::Tolerance;
use model_train::geometry::GeometryConfig;
use search::PrecursorCalMode;

/// Fragmentation method. `Auto` detects from the mzML's activation block and
/// falls back to the bundled `hcd_qexactive_tryp` model when nothing is detected.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum Fragmentation {
    #[clap(name = "auto")]
    Auto,
    #[clap(name = "CID")]
    Cid,
    #[clap(name = "ETD")]
    Etd,
    #[clap(name = "HCD")]
    Hcd,
    #[clap(name = "UVPD")]
    Uvpd,
}

/// Search protocol: sample labeling or enrichment strategy applied during the experiment.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum Protocol {
    #[clap(name = "auto")]
    Auto,
    #[clap(name = "phospho")]
    Phospho,
    #[clap(name = "iTRAQ")]
    Itraq,
    #[clap(name = "iTRAQ-phospho")]
    ItraqPhospho,
    #[clap(name = "TMT")]
    Tmt,
    #[clap(name = "standard")]
    Standard,
}

/// Enzymatic-cleavage enforcement at peptide span boundaries:
/// 2=fully, 1=semi, 0=non-specific.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum EnzymeSpecificity {
    #[clap(name = "non-specific")]
    NonSpecific,
    #[clap(name = "semi")]
    Semi,
    #[clap(name = "fully")]
    Fully,
}

/// Primary ranking mode: inherited RawScore (`rank`) or fused strong score (`strong`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum)]
pub(crate) enum ScoreFlag {
    /// Pick by the resolved model's instrument: `strong` for high-res, `rank` for low-res.
    #[default]
    Auto,
    Rank,
    Strong,
}

/// Fragment-ion index for out-of-core searches (issue #76).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum)]
pub(crate) enum FragmentIndexFlag {
    /// On when the candidate index is out-of-core and fragment matching is
    /// high-resolution (default).
    #[default]
    Auto,
    /// On in out-of-core mode; a warning (and enumeration) when the candidate
    /// index is in RAM.
    On,
    /// Per-spectrum enumeration even out-of-core.
    Off,
}

/// Candidate-resolution backing: in-RAM (`ram`) or out-of-core mmap base-peptide
/// index with lazy mod enumeration (`mmap`); `auto` (default) picks one.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum)]
pub(crate) enum CandidateIndexFlag {
    /// Pick automatically (default): out-of-core mmap with fragment-ion retrieval
    /// for high-resolution fragments; otherwise RAM, unless the in-RAM index
    /// would not fit available memory. `--refine` and `--glyco` stay on RAM.
    #[default]
    Auto,
    /// Force the in-RAM candidate index (advanced; may OOM on very large mod
    /// searches — that is what `auto` protects against).
    Ram,
    /// Force the out-of-core mmap'd base-peptide index (advanced; lower peak RAM,
    /// result-equivalent but not byte-identical to `ram`).
    Mmap,
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
#[command(group(
    clap::ArgGroup::new("glycan_db")
        .args(["glyco_glycan_gdb", "glyco_species"])
        .multiple(false)
        .required(false),
))]
pub(crate) struct SearchArgs {
    /// YAML run-configuration file. Any parameter can be set here (grouped by
    /// experiment: io/search/scoring/decoys/chimeric/refine/rescoring/glyco; see
    /// DOCS §1b). An explicit CLI flag always overrides the config value.
    #[arg(long = "config", value_name = "FILE")]
    pub(crate) config: Option<PathBuf>,

    /// Input spectrum file(s). Repeat `--spectrum` for multiple inputs (one PIN).
    /// Format is auto-detected per file by extension.
    #[arg(long)]
    pub(crate) spectrum: Vec<PathBuf>,

    /// Input FASTA database (target sequences only; decoys are generated automatically).
    #[arg(long)]
    pub(crate) database: Option<PathBuf>,

    /// Output Percolator PIN file path.
    #[arg(long)]
    pub(crate) output_pin: Option<PathBuf>,

    /// Output TSV file path (optional).
    #[arg(long)]
    pub(crate) output_tsv: Option<PathBuf>,

    /// Output QPX `.idparquet/` bundle directory (optional; OpenMS-compatible).
    /// Writes `psms.parquet` + `proteins.parquet` + `search_params.parquet`.
    #[arg(long)]
    pub(crate) output_parquet: Option<PathBuf>,

    /// Decoy prefix used when generating reversed decoy sequences.
    #[arg(long, default_value = "XXX_")]
    pub(crate) decoy_prefix: String,

    /// Decoy-accession SUFFIX used to RECOGNIZE pre-built decoys in the input
    /// FASTA (e.g. `rev` for quantms/OpenMS `<orig>_rev` decoys). When set, a
    /// protein is a decoy iff its accession starts with `<decoy-prefix>_` OR ends
    /// with this suffix. Typically paired with `--decoy-strategy none` so andes
    /// consumes an externally-built target+decoy database instead of generating
    /// its own decoys (which would double-decoy and bias FDR).
    #[arg(long = "decoy-suffix")]
    pub(crate) decoy_suffix: Option<String>,

    /// How to generate decoys: `reverse` (default; reverse each sequence),
    /// `shuffle` (seeded reproducible shuffle), `sequon-reverse` (reverse but
    /// restore each N-X-S/T sequon at its mirrored position — RECOMMENDED with
    /// `--glyco`: plain reversal maps N-X-S/T to S/T-X-N, so reversed decoys reach
    /// the glyco sequon gate at a lower rate than targets and the resulting
    /// q-values are anti-conservative), or `none` (no decoys — for a FASTA that
    /// already contains decoys, or external FDR). `none` with a target-only FASTA
    /// leaves the search without decoys (FDR can't be estimated) and warns.
    #[arg(long = "decoy-strategy", default_value = "reverse")]
    pub(crate) decoy_strategy: String,

    /// Seed for `--decoy-strategy shuffle` (reproducible decoys). Ignored by
    /// reverse/none.
    #[arg(long = "decoy-seed", hide = true, default_value_t = search::decoy::DEFAULT_DECOY_SEED)]
    pub(crate) decoy_seed: u64,

    /// Isotope-error offset range to try, as `MIN..MAX` (also accepts `MIN-MAX`).
    /// Negative offsets allowed. Unset defaults to `-1..2`, or `0..2` under `--glyco`
    /// (see the resolution site). Left as an `Option` so an EXPLICIT `-1..2` is
    /// distinguishable from the default and is never silently overridden.
    #[arg(long = "isotope-error", hide = true, value_parser = parse_isotope_error_range)]
    pub(crate) isotope_error: Option<(i8, i8)>,

    /// Precursor-mass calibration: `off` (default), `auto`, or `on`. `auto`/`on`
    /// learn a systematic ppm shift from confident PSMs in a pre-pass and tighten
    /// the precursor tolerance for the main search; `auto` skips the correction
    /// when the sample is too small to be reliable. Off by default: on the
    /// benchmark sets the pre-pass cost 9-17% of the run and fired on none of
    /// Astral, TMT, UPS1 or phospho; where a larger sample made it fire (UPS1),
    /// the tightened window lost 7.6% of PSMs.
    #[arg(long = "precursor-cal", default_value = "off", value_parser = parse_precursor_cal)]
    pub(crate) precursor_cal: PrecursorCalMode,

    /// Minimum number of (spectrum, charge) keys before the calibration pre-pass
    /// runs; below it calibration is skipped as unreliable. Default 10000. Hidden:
    /// exists so the out-of-core equivalence gate can force the pre-pass on the
    /// small in-repo fixture; production runs should leave it alone.
    #[arg(long = "cal-min-spec-keys", hide = true)]
    pub(crate) cal_min_spec_keys: Option<usize>,

    /// Candidate retrieval for out-of-core searches. `auto` (the default, and
    /// what every normal run uses) selects the fragment-ion index whenever the
    /// candidate index is out-of-core; at a low-resolution fragment tolerance
    /// only the 150 most intense peaks vote. `off` forces per-spectrum
    /// enumeration and `on` is the same as `auto` — both are escape hatches for
    /// reproducing a measurement, not tuning knobs. Hidden on purpose: the
    /// engine reports which strategy it chose.
    #[arg(long = "fragment-index", hide = true, default_value = "auto")]
    pub(crate) fragment_index: FragmentIndexFlag,

    /// Fragment-index mode: precursor-mass width of one index slice in Da.
    /// Default: derived from the memory budget (10–150 Da; ~0.13 GB per Da on
    /// a phospho search).
    #[arg(long = "fragment-index-slice-da", hide = true)]
    pub(crate) fragment_index_slice_da: Option<f64>,

    /// Precursor mass tolerance as `VALUE+unit`. Accepts ppm (e.g. `20ppm`,
    /// high-res) or Da (e.g. `0.02da`/`0.02Da`, low-res precursor selection).
    /// Default `20ppm`.
    #[arg(long = "precursor-tol", default_value = "20ppm", value_parser = parse_precursor_tol)]
    pub(crate) precursor_tol: Tolerance,

    /// Precursor charge range to try when not specified in the spectrum, as
    /// `MIN..MAX` (also accepts `MIN-MAX`). Default `2..5`.
    #[arg(long = "charge", hide = true, default_value = "2..5", value_parser = parse_charge_range)]
    pub(crate) charge: (u8, u8),

    /// Maximum number of PSMs to retain per spectrum.
    #[arg(long, hide = true, default_value = "10")]
    pub(crate) top_n: u32,

    /// Number of Tolerable Termini (enzymatic-cleavage enforcement at span
    /// boundaries). `fully`: both termini must be cleavage sites (strict).
    /// `semi`: at least one terminus must be a cleavage site. `non-specific`:
    /// neither terminus needs to be a cleavage site.
    #[arg(long = "enzyme-specificity",
          hide = true, default_value = "fully", value_parser = parse_enzyme_specificity)]
    pub(crate) enzyme_specificity: EnzymeSpecificity,

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
    pub(crate) enzyme: String,

    /// Maximum number of missed cleavages per peptide.
    #[arg(long, hide = true, default_value = "1")]
    pub(crate) max_missed_cleavages: u32,

    /// Minimum number of peaks an MS2 spectrum must have to be scored; spectra
    /// with fewer peaks are skipped.
    #[arg(long, hide = true, default_value = "10")]
    pub(crate) min_peaks: u32,

    /// Minimum peptide length, in residues.
    #[arg(long, hide = true, default_value = "6")]
    pub(crate) min_length: u32,

    /// Maximum peptide length, in residues. (50 matches the reference engine/a comparison engine defaults;
    /// 40 dropped long tryptic peptides.)
    #[arg(long, hide = true, default_value = "50")]
    pub(crate) max_length: u32,

    /// Maximum number of variable modifications per peptide. A `NumMods=N` line
    /// in a --mods file overrides this.
    #[arg(long = "max-mods", hide = true, default_value = "3")]
    pub(crate) max_mods: u32,

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
    pub(crate) mods: Option<PathBuf>,

    /// Fragmentation/activation method for MGF input only. mzML/.raw/.d
    /// auto-detect this. Named values: auto, CID, ETD, HCD, UVPD.
    #[arg(long, hide = true, default_value = "auto", value_parser = parse_fragmentation)]
    pub(crate) fragmentation: Fragmentation,

    /// Search protocol. Named values: auto, phospho, iTRAQ, iTRAQ-phospho, TMT, standard.
    #[arg(long, hide = true, default_value = "auto", value_parser = parse_protocol)]
    pub(crate) protocol: Protocol,

    /// Fragment-matching tolerance in ppm for **MGF input only** (high-resolution
    /// MS/MS). Has no effect on mzML/.raw/.d (analyzer auto-detected). Mutually
    /// exclusive with `--fragment-tol-da`.
    #[arg(long = "fragment-tol-ppm", hide = true, conflicts_with = "fragment_tol_da", value_parser = parse_positive_tol)]
    pub(crate) fragment_tol_ppm: Option<f64>,

    /// Fragment-matching tolerance in Da for **MGF input only** (low-resolution
    /// ion-trap MS/MS). Has no effect on mzML/.raw/.d. Mutually exclusive with
    /// `--fragment-tol-ppm`.
    #[arg(long = "fragment-tol-da", hide = true, conflicts_with = "fragment_tol_ppm", value_parser = parse_positive_tol)]
    pub(crate) fragment_tol_da: Option<f64>,

    /// Number of worker threads for the search loop. Defaults to logical CPU count.
    #[arg(long, default_value_t = num_cpus::get())]
    pub(crate) threads: usize,

    /// Debug/benchmark cap: process only the first N spectra (0 = no cap).
    #[arg(long, hide = true, default_value = "0")]
    pub(crate) max_spectra: usize,

    /// MS level to search. Defaults to MS2 (identification); MS1 and any higher
    /// levels (e.g. TMT SPS-MS3 reporter-quant scans) are filtered out at load
    /// time so they never enter the search loop. Override only if you explicitly
    /// want a different level. Applies to mzML and Thermo `.raw`; MGF files do
    /// not encode MS level and are always treated as MS2. The chimeric cascade
    /// always searches MS2 (it pairs MS2 with its preceding MS1).
    #[arg(long, hide = true, default_value = "2")]
    pub(crate) ms_level: u8,

    /// Enable the two-pass chimeric cascade for co-isolated (co-fragmented)
    /// peptides. Pass 1 is the normal top-1 search; Pass 2 detects co-isolated
    /// precursors in each scan's MS1 isolation window and runs a targeted search
    /// for the second peptide on the residual spectrum, emitting it as an extra
    /// PSM. Requires mzML (MS1 scans); has no effect on MGF input.
    #[arg(long, default_value = "false")]
    pub(crate) chimeric: bool,

    /// Chimeric mode: max co-isolated SECONDARY peptides to search per scan (the
    /// chimeric-N lever). Default 4 = the measured Astral sweet spot (+1.4% PSMs
    /// vs N=2 at flat FDP; saturates by N=4). Set 2 for the original behavior.
    #[arg(long = "chimeric-max-coisolated", hide = true, default_value = "4")]
    pub(crate) chimeric_max_coisolated: usize,

    /// Chimeric mode: averagine-envelope KL gate for accepting a co-isolated MS1
    /// envelope (lower = stricter/cleaner; fewer spurious secondaries).
    #[arg(long = "chimeric-max-kl", hide = true, default_value = "0.3")]
    pub(crate) chimeric_max_kl: f32,

    /// Path to a Parquet model store (a single file or a partitioned directory)
    /// to use instead of the bundled `resources/models/`. When set, model selection reads from
    /// this store; when unset, the bundled store is used.
    #[arg(long = "model-store", hide = true)]
    pub(crate) model_store: Option<PathBuf>,

    /// Exact model ID to load from the model store (bundled or `--model-store`).
    /// When set, skips automatic selection (metadata detection / `--fragmentation`
    /// / `--protocol`) and loads this ID directly. Useful after `andes train`
    /// to search with the freshly-trained model.
    #[arg(long = "model", hide = true)]
    pub(crate) model_id_override: Option<String>,

    /// Evaluate only the first N trees of each GBDT ensemble (0 = all trees).
    ///
    /// Applies to BOTH shipped ensembles — fragment-intensity and rich-ion — which are
    /// of comparable size and are each walked per fragment per candidate. `Tree::eval`
    /// on them is the single hottest operation in a standard search (77% of self time
    /// in a native profile). The GBDT is additive, so the early trees carry the signal
    /// and later ones only refine it.
    ///
    /// WHAT TRUNCATION ACTUALLY AFFECTS, by score mode -- these differ, and conflating
    /// them is a documented past error:
    ///   `rank` (low-res default): the ensembles feed PIN feature columns only, never
    ///     the ranking score, so the emitted PSM row SET is byte-identical at every K
    ///     (verified: UPS1 emitted 382,703 rows at K=0, 100 and 25 alike).
    ///   `strong` (high-res default): `reorder_by_strong_score` ranks by StrongScore,
    ///     which consumes `intensity_signal` and therefore the frag-intensity ensemble.
    ///     Truncation there CHANGES WHICH PEPTIDE WINS, not just feature values
    ///     (verified: Astral emitted 1,214,417 / 1,214,757 / 1,214,892 / 1,215,148 rows
    ///     at K=50 / 25 / 10 / 1). The identification counts below are still flat at
    ///     K=100, but that is an empirical result about yield, NOT evidence that the
    ///     model is inert to selection.
    ///
    /// DEFAULT 100 for standard search, measured, not guessed. Five Percolator seeds
    /// per point, identical data and protocol:
    ///   Astral high-res: 300 trees 38,437 PSMs @1% FDR; 100 trees 38,444; 50 trees
    ///     38,357; 25 trees 38,094; 10 trees 37,769; 1 tree 37,533.
    ///   UPS1 low-res (yeast target + E. coli entrapment): 300 trees 15,813 PSMs at
    ///     3.39% true FDP; 100 trees 15,832 at 3.44%; 25 trees 15,779 at 3.61%.
    ///     (FDP = entrapment_hits * (1 + T/E) / n with T/E = 2.419 measured over the
    ///     searchable tryptic space, 734,280 yeast vs 303,537 E. coli peptides. An
    ///     earlier revision of this comment printed ~2% here from a 1:1 assumption
    ///     that does not hold for this FASTA; the ARM-TO-ARM comparison is unaffected
    ///     because the scale factor is common to all three.)
    /// So the last 200 trees buy nothing on either resolution class while costing
    /// 33–41% of wall time, and true error is unchanged. Below ~50 the loss becomes
    /// real. (For scale: even ONE tree still beats Comet by ~19% on Astral, so the
    /// whole ensemble is worth ~2% of identifications for 77% of the runtime.)
    ///
    /// GLYCO IS EXEMPT: under `--glyco` the ensembles are ~0.3% of wall time, so
    /// truncation would trade identification fidelity for nothing. Pass the flag
    /// explicitly to truncate there anyway.
    ///
    /// This CHANGES predicted intensities and therefore the emitted PIN feature
    /// values; `--gbdt-max-trees 0` restores the full ensembles.
    #[arg(long = "gbdt-max-trees", hide = true, default_value_t = 100usize)]
    pub(crate) gbdt_max_trees: usize,

    /// Ranking / PIN RawScore source: `auto` (default — `strong` for high-res
    /// instruments, `rank` for low-res), `rank`, or `strong` (fused intensity +
    /// competition score from S1–S3).
    #[arg(long = "score", default_value = "auto")]
    pub(crate) score: ScoreFlag,

    /// Candidate-index backing: `auto` (default — out-of-core mmap with
    /// fragment-ion retrieval for high-resolution fragments, cached in the temp
    /// directory; otherwise RAM unless it would not fit available memory), or
    /// force `ram` / `mmap` (advanced overrides). `mmap` lowers peak RAM with lazy per-spectrum mod
    /// enumeration (result-equivalent PSMs, not byte-identical).
    #[arg(long = "candidate-index", hide = true, default_value = "auto")]
    pub(crate) candidate_index: CandidateIndexFlag,

    /// Glycopeptide search mode: retrieve glycans from the loaded database, solve
    /// the backbone as `precursor − glycan`, filter by N-X-S/T sequon, and write a
    /// `.glyco.pin` file instead of the standard PIN. Default off.
    ///
    /// Requires a glycan database: `--glyco-species` for a bundled one, or
    /// `--glyco-glycan-gdb` for your own file. The requirement is enforced at
    /// argument-parsing time, so a missing database fails immediately rather than
    /// after the digest, index build and spectrum load.
    #[arg(long = "glyco", default_value_t = false, requires = "glycan_db")]
    pub(crate) glyco: bool,

    /// Cap the peaks the glyco GENERATION stage considers, keeping the most
    /// intense N. The backbone solver is superlinear in peak count, so an
    /// uncentroided profile scan or a very dense wide-window scan can take tens of
    /// seconds while a normal scan takes milliseconds — the run looks hung.
    /// Scoring always reads the full spectrum, so a generated candidate is never
    /// scored on truncated evidence. Default 0 = no cap; 300-500 is a reasonable
    /// value if you hit this. Changing it changes results.
    #[arg(long = "glyco-max-peaks", hide = true, default_value_t = 0usize)]
    pub(crate) glyco_max_peaks: usize,

    /// Maximum c/z fragment charge to probe in `--glyco` ETD scoring. Unset derives it
    /// from whether the spectrum was deconvoluted, which is correct in almost all cases:
    /// after deconvolution multiply-charged fragments have already been moved to 1+.
    /// Set this only for data known to carry unresolved high-charge c/z ions.
    #[arg(long = "glyco-cz-max-charge", hide = true)]
    pub(crate) glyco_cz_max_charge: Option<u8>,

    /// Maximum glycan-Y fragment charge. Default 3; raising it reaches 4+/5+ Y ions on
    /// highly-charged precursors at the cost of more chance matches.
    #[arg(long = "glyco-y-max-charge", hide = true, default_value_t = 3u8)]
    pub(crate) glyco_y_max_charge: u8,

    /// Windowed peak filtering as `WINDOW_DA:PEAKS` (e.g. `100:20`). Unset uses the
    /// protocol default — on for isobaric-labelled data, off otherwise. A window of 0
    /// forces it off.
    #[arg(long = "peak-filter", hide = true)]
    pub(crate) peak_filter: Option<String>,

    /// How to label EThcD/ETciD spectra (electron transfer with a supplemental
    /// collisional term). `hcd` is the default and is what model routing expects, since
    /// no EThcD model exists; `etd` labels them ETD so the c/z scoring path engages.
    #[arg(long = "ethcd-activation", hide = true, value_enum, default_value_t = EthcdActivationFlag::Hcd)]
    pub(crate) ethcd_activation: EthcdActivationFlag,

    /// Diagnostic: restrict `--glyco` scoring to the scan numbers in this file, one per
    /// line. Makes a `--debug-glyco` dump of a chosen set of scans affordable.
    #[arg(long = "glyco-scans", hide = true)]
    pub(crate) glyco_scans: Option<PathBuf>,

    /// Memory: enumerate only N-X-S/T-bearing peptides into the `--glyco` candidate
    /// index. Glyco scoring never reads any other candidate (the peptide-first
    /// fragment index and the backbone→candidate loop both gate on the same sequon
    /// test), so this shrinks the in-RAM index roughly 10x. Measured on a
    /// 1,500-protein mouse subset: identical rows, peptides and labels; 16 of 7,113
    /// rows (0.2%) differ in `RawScore`/`CandidateRankEntropy` only (per-spectrum
    /// candidate-distribution terms). Use on hosts where the full index does not fit
    /// (the mouse entrapment recipe needs ~27 GB without it, 3.4 GB with it).
    /// Hidden; requires `--glyco`.
    #[arg(long = "glyco-index-sequon-only", hide = true, default_value_t = false)]
    pub(crate) glyco_index_sequon_only: bool,

    /// Load an external pGlyco-style `.gdb` glycan database: the glycan-first search
    /// space. Each canonical string is parsed with its tree structure preserved
    /// (core- vs antenna-fucose), and the file's NeuGc content is used as-is (no
    /// species filtering). Use a real pGlyco `.gdb` (e.g. `pGlyco-N-Human.gdb`) for a
    /// like-for-like comparison against pGlyco. Takes precedence over `--glyco-species`.
    #[arg(long = "glyco-glycan-gdb")]
    pub(crate) glyco_glycan_gdb: Option<PathBuf>,

    /// Select a bundled species-specific N-glycan database instead of an external
    /// `--glyco-glycan-gdb` file. `--glyco` requires one of the two; an explicit
    /// `--glyco-glycan-gdb` wins. Databases are pGlyco's (see
    /// `crates/andes-glyco/glycan-db/`).
    #[arg(long = "glyco-species", value_enum)]
    pub(crate) glyco_species: Option<GlycoSpeciesFlag>,

    /// Correct each precursor to its monoisotopic peak from the preceding MS1
    /// isotope envelope BEFORE searching (`--glyco`; needs MS1, so mzML or Thermo
    /// `.raw`). `auto` fits the observed envelope at the reported charge against a
    /// glycopeptide isotope model under the hypotheses "the recorded precursor is
    /// the M+k peak", k = 0..6, and moves the precursor down when a k > 0 clearly
    /// wins; every spectrum whose envelope was fitted is then searched with a
    /// `0..1` isotope window instead of the glyco default `0..2` (the `+2` step
    /// only carries the Hex+Fuc/NeuGc composition degeneracy once the monoisotope
    /// is verified; bigbio/andes#64 arm F), while spectra with no linked MS1 or
    /// charge keep the default window. An explicit `--isotope-error` is honoured.
    /// `auto` is the default; `off` leaves every precursor as recorded. Motivation: on pGlyco2
    /// mouse liver the firmware records 3-6 Da above the monoisotope on a class of
    /// wide, high-mass glycopeptide envelopes, and widening the isotope window to
    /// reach them is mass-degenerate with a glycan composition change
    /// (bigbio/andes#64). Byte-identical output with `off`, or without MS1.
    #[arg(long = "precursor-mono", value_enum, default_value_t = PrecursorMonoFlag::Auto)]
    pub(crate) precursor_mono: PrecursorMonoFlag,

    /// Diagnostic: write one TSV row per MS2 with the recorded/corrected precursor,
    /// the applied shift, and the envelope fit under every hypothesis. Hidden.
    #[arg(long = "precursor-mono-dump", hide = true)]
    pub(crate) precursor_mono_dump: Option<PathBuf>,

    /// Fragment tolerance (ppm) for the glyco-specific matching: oxonium ions,
    /// the core-Y ladder, backbone mass search, and c/z. Default 20 ppm, which
    /// suits Orbitrap MS2. **Raise this for low-resolution (ion-trap) MS2** —
    /// at 20 ppm a 0.3-0.5 Da ion-trap peak never matches, so the oxonium gate
    /// never fires and glyco IDs collapse to near zero. This is separate from
    /// `--fragment-tol-ppm`, which the scoring model owns.
    #[arg(long = "glyco-tol-ppm", hide = true, default_value_t = 20.0f64, value_parser = parse_positive_tol)]
    pub(crate) glyco_tol_ppm: f64,

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
    pub(crate) glyco_min_core_y: u32,

    /// Peptide-first candidate RETRIEVAL tolerance in ppm. Default: --glyco-tol-ppm
    /// on high-resolution MS2, the rank model's 0.5 Da window on low-resolution.
    /// Retrieval only; the rank scorer and its tolerance are unchanged. Measured
    /// 7x faster than 0.5 Da on high-res data with identifications neutral.
    #[arg(long = "glyco-retrieval-tol-ppm", hide = true, value_parser = parse_positive_tol)]
    pub(crate) glyco_retrieval_tol_ppm: Option<f64>,

    /// Use the peptide-first b/y fragment-index fallback instead of the default
    /// mass-driven full-glycan-list branch. On a mouse-liver fraction it finds ~2%
    /// more glycoPSMs at ~13x the run time and ~1.5x the peak memory.
    #[arg(long = "glyco-peptide-first", hide = true, default_value_t = false)]
    pub(crate) glyco_peptide_first: bool,

    /// Diagnostic glyco mode: emit ALL candidate rows per scan (including de-novo
    /// mass-residual hits). The resulting PIN is for inspection ONLY and must never
    /// be fed to an FDR tool. Hidden dev flag.
    #[arg(long = "debug-glyco", hide = true, default_value_t = false)]
    pub(crate) debug_glyco: bool,

    /// On ETD/AI-ETD spectra, generate candidate backbones from the paired HCD scan
    /// (same precursor) while scoring c/z on the ETD scan — targets high-charge
    /// glycopeptides (validated +153 backbone-correct @1%). DEFAULT ON; scans with no
    /// HCD partner (and multi-file runs) fall back to unpaired automatically. Disable
    /// with `--glyco-hcd-pair false`. `action = Set` so the bool takes an explicit value.
    #[arg(long = "glyco-hcd-pair", default_value_t = true, action = clap::ArgAction::Set)]
    pub(crate) glyco_hcd_pair: bool,

    /// Enable the PTM-refinement cascade (Pass-2 over confident proteins). Default off.
    #[arg(long = "refine", default_value_t = false)]
    pub(crate) refine: bool,

    /// YAML refinement config; omit to use the built-in 5-mod DEFAULT tier.
    #[arg(long = "refine-config", hide = true)]
    pub(crate) refine_config: Option<std::path::PathBuf>,

    /// Confident-anchor SCOPING FDR (not a reported FDR). Default 0.01 — the same
    /// internal TDC q used for calibration/training/report. A looser gate (e.g.
    /// 0.10) admits low-confidence anchors that leak into the entrapment-FDP
    /// (b1931: 0.10 → 4.86% vs 0.01 → 0.29% true FDP). Hidden power-user knob;
    /// leave at the default unless you have a measured reason to widen it.
    #[arg(long = "refine-select-psm-fdr", default_value_t = 0.01, hide = true, value_parser = parse_unit_fraction)]
    pub(crate) refine_select_psm_fdr: f64,

    /// Run Percolator on the PIN after the search and join its PEP/q-value back
    /// into the outputs (QPX `posterior_error_probability` + a `q-value` score,
    /// and a filtered `<stem>.q<fdr>.tsv`). Needs a Percolator backend (see
    /// `--percolator-bin` / `--percolator-docker`). Default off.
    #[arg(long = "rescore", default_value_t = false)]
    pub(crate) rescore: bool,

    /// Rescore with the built-in NATIVE GBDT rescorer instead of Percolator (no
    /// Percolator backend needed). Leakage-safe 3-fold target-decoy cross-
    /// validation over the PIN features → q-value + PEP. A self-contained
    /// FALLBACK for benchmarking / offline use — NOT production-grade FDR; prefer
    /// `--rescore` (Percolator) for production. Writes the same QPX q-value/PEP +
    /// filtered `<stem>.q<fdr>.tsv` outputs. Ignored if `--rescore` is also set.
    #[arg(long = "rescore-native", hide = true, default_value_t = false)]
    pub(crate) rescore_native: bool,

    /// FDR (q-value) threshold for the filtered `<stem>.q<fdr>.tsv` output
    /// (target PSMs at q ≤ this).
    ///
    /// This is a THRESHOLD, not a switch: it does nothing on its own. andes does
    /// not compute FDR unless you ask it to, because doing so would silently
    /// launch Percolator or the non-production native rescorer. Pair it with
    /// `--rescore` (Percolator) or `--rescore-native`. Setting it alone warns and
    /// writes no filtered output. Defaults to 0.01 when rescoring does run.
    #[arg(long = "fdr", value_parser = parse_unit_fraction)]
    pub(crate) fdr: Option<f64>,

    /// Optional per-PSM PEP (posterior error probability / local FDR) cap,
    /// applied IN ADDITION to `--fdr` (a PSM must pass both q ≤ `--fdr` AND
    /// PEP ≤ `--pep`). The q-value stays the primary set-level FDR control;
    /// `--pep` is a supplementary per-PSM gate. Like `--fdr` it is a threshold
    /// only: it does nothing without `--rescore` / `--rescore-native`.
    /// Default: no PEP cap.
    #[arg(long = "pep", hide = true, value_parser = parse_unit_fraction)]
    pub(crate) pep: Option<f64>,

    /// Explicit path to a Percolator binary (highest-priority backend). When
    /// omitted, `percolator` on `$PATH` is used, else the docker fallback.
    #[arg(long = "percolator-bin", hide = true)]
    pub(crate) percolator_bin: Option<std::path::PathBuf>,

    /// Force the Percolator docker fallback (the pinned biocontainers image)
    /// instead of looking for a native binary. Requires the `docker` CLI.
    #[arg(long = "percolator-docker", hide = true, default_value_t = false)]
    pub(crate) percolator_docker: bool,

    /// Percolator docker image tag for the docker fallback (power-user override).
    #[arg(long = "percolator-image", hide = true, default_value = output::DEFAULT_PERCOLATOR_IMAGE)]
    pub(crate) percolator_image: String,

    /// Extra arguments passed verbatim to Percolator (after the fixed flags,
    /// before the PIN path). e.g. `--percolator-args "--testFDR 0.05"`.
    #[arg(long = "percolator-args", hide = true, default_value = "")]
    pub(crate) percolator_args: String,

    /// Keep the PIN file after rescoring. With `--rescore` and no `--output-pin`,
    /// a temporary PIN is used and deleted unless this is true. Default true.
    #[arg(long = "keep-pin", hide = true, default_value_t = true, action = clap::ArgAction::Set)]
    pub(crate) keep_pin: bool,
}

// Alias used internally for the search-args type.
pub(crate) type Cli = SearchArgs;

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
        Self {
            segments: 2,
            max_rank: 150,
            occupancy: 2500,
            max_tiers: 33,
            max_fragment_charge: 3,
        }
    }
}

impl GeometryArgs {
    pub(crate) fn to_config(self) -> GeometryConfig {
        GeometryConfig {
            num_segments: self.segments.max(1),
            max_rank: self.max_rank.max(1),
            mass_tier_occupancy: self.occupancy.max(1),
            max_mass_tiers: self.max_tiers.max(1),
            max_fragment_charge: self.max_fragment_charge.max(1),
        }
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, clap::ValueEnum)]
pub(crate) enum EthcdActivationFlag {
    /// Label EThcD/ETciD as HCD (default; matches model routing).
    Hcd,
    /// Label them ETD so the c/z scoring path engages.
    Etd,
}

/// Bundled species-specific N-glycan databases for `--glyco-species`.
#[derive(Copy, Clone, Debug, PartialEq, Eq, clap::ValueEnum)]
pub(crate) enum GlycoSpeciesFlag {
    /// Human N-glycans.
    Human,
    /// Human, multi-antennary.
    HumanMulti,
    /// Mouse N-glycans.
    Mouse,
    /// Mouse, extended.
    MouseLarge,
    /// High-mannose N-glycans (species-independent).
    HighMannose,
}

impl GlycoSpeciesFlag {
    /// The kebab-case key `load_species_glycan_db` expects.
    pub(crate) fn key(self) -> &'static str {
        match self {
            GlycoSpeciesFlag::Human => "human",
            GlycoSpeciesFlag::HumanMulti => "human-multi",
            GlycoSpeciesFlag::Mouse => "mouse",
            GlycoSpeciesFlag::MouseLarge => "mouse-large",
            GlycoSpeciesFlag::HighMannose => "high-mannose",
        }
    }
}

/// `--precursor-mono`: MS1 isotope-envelope precursor correction.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub(crate) enum PrecursorMonoFlag {
    /// Leave every precursor m/z as recorded (default).
    Off,
    /// Move a precursor down to the monoisotope its MS1 envelope supports.
    Auto,
}

/// Parse `--fragmentation` value. Accepts named values (case-insensitive: auto,
/// CID, ETD, HCD, UVPD).
pub(crate) fn parse_fragmentation(s: &str) -> Result<Fragmentation, String> {
    <Fragmentation as ValueEnum>::from_str(s, true)
        .map_err(|_| format!("invalid fragmentation `{s}`: expected auto|CID|ETD|HCD|UVPD"))
}

/// Parse `--protocol` value. Accepts named values only.
pub(crate) fn parse_protocol(s: &str) -> Result<Protocol, String> {
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
pub(crate) fn dash_sep_index(s: &str) -> Option<usize> {
    let b = s.as_bytes();
    (1..b.len()).find(|&i| b[i] == b'-' && b[i - 1].is_ascii_digit())
}

/// Parse a `MIN..MAX` (or `MIN-MAX`) range into a `(min, max)` pair, generic
/// over the integer type so it serves both `--charge` (u8) and
/// `--isotope-error` (i8, negatives allowed). The `-` separator is tried only
/// when the value does not parse as `..`; signed endpoints are supported in both
/// forms (`-1..2`, `-3--1`).
pub(crate) fn parse_int_range<T>(s: &str, label: &str) -> Result<(T, T), String>
where
    T: std::str::FromStr + PartialOrd + std::fmt::Display + Copy,
{
    let trimmed = s.trim();
    let (lo_s, hi_s) = if let Some((a, b)) = trimmed.split_once("..") {
        (a.trim(), b.trim())
    } else if let Some(idx) = dash_sep_index(trimmed) {
        (trimmed[..idx].trim(), trimmed[idx + 1..].trim())
    } else {
        return Err(format!(
            "invalid {label} `{s}`: expected MIN..MAX (or MIN-MAX)"
        ));
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
pub(crate) const MAX_SUPPORTED_CHARGE: u8 = 50;

/// Parse `--charge MIN..MAX` (also `MIN-MAX`) into `(u8, u8)`.
///
/// Domain-validates (finding 3.8): the minimum charge must be >= 1 (charge 0 is
/// not a real precursor) and the maximum must be within the supported bound.
pub(crate) fn parse_charge_range(s: &str) -> Result<(u8, u8), String> {
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
pub(crate) fn parse_isotope_error_range(s: &str) -> Result<(i8, i8), String> {
    parse_int_range::<i8>(s, "isotope-error")
}

/// Parse `--precursor-tol VALUE+unit` (e.g. `20ppm`, `0.02da`/`0.02Da`).
pub(crate) fn parse_precursor_tol(s: &str) -> Result<Tolerance, String> {
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
    Ok(if is_ppm {
        Tolerance::Ppm(v)
    } else {
        Tolerance::Da(v)
    })
}

/// Parse a probability-domain CLI value (FDR / PEP / refine-FDR) — must be a
/// finite number in `[0, 1]` (finding 3.8). Used as a clap `value_parser` so a
/// bad value is rejected at parse time with a clear message.
pub(crate) fn parse_unit_fraction(s: &str) -> Result<f64, String> {
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
pub(crate) fn parse_positive_tol(s: &str) -> Result<f64, String> {
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
pub(crate) fn parse_precursor_cal(s: &str) -> Result<PrecursorCalMode, String> {
    match s.to_ascii_lowercase().as_str() {
        "auto" => Ok(PrecursorCalMode::Auto),
        "on" => Ok(PrecursorCalMode::On),
        "off" => Ok(PrecursorCalMode::Off),
        _ => Err(format!(
            "invalid precursor-cal `{s}`: expected auto|on|off (Java -precursorCal)"
        )),
    }
}

pub(crate) fn parse_enzyme_specificity(s: &str) -> Result<EnzymeSpecificity, String> {
    <EnzymeSpecificity as ValueEnum>::from_str(s, true)
        .map_err(|_| format!("invalid enzyme specificity `{s}`: expected non-specific|semi|fully"))
}

#[cfg(test)]
mod cli_domain_validator_tests {
    use super::*;

    #[test]
    fn charge_range_rejects_zero_and_out_of_bounds() {
        assert!(parse_charge_range("2..5").is_ok());
        assert!(parse_charge_range("1..1").is_ok());
        assert!(
            parse_charge_range("0..5").is_err(),
            "charge 0 must be rejected"
        );
        assert!(
            parse_charge_range("0..0").is_err(),
            "charge 0..0 must be rejected"
        );
        assert!(
            parse_charge_range("2..60").is_err(),
            "charge above bound must be rejected"
        );
    }

    #[test]
    fn precursor_tol_rejects_nonpositive_and_nonfinite() {
        assert!(parse_precursor_tol("20ppm").is_ok());
        assert!(parse_precursor_tol("0.02da").is_ok());
        assert!(parse_precursor_tol("0ppm").is_err(), "zero tol rejected");
        assert!(
            parse_precursor_tol("-5ppm").is_err(),
            "negative tol rejected"
        );
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
    use crate::model_select::parse_enzymes;
    use model::enzyme::Enzyme;

    #[test]
    fn single_enzyme_has_no_extras() {
        let (primary, extras) = parse_enzymes("trypsin").unwrap();
        assert_eq!(primary, Enzyme::Trypsin);
        assert!(
            extras.is_empty(),
            "single enzyme ⇒ empty extras (bit-identical path)"
        );
    }

    #[test]
    fn comma_and_plus_separators_both_parse() {
        let (p1, e1) = parse_enzymes("gluc,trypsin").unwrap();
        let (p2, e2) = parse_enzymes("gluc+trypsin").unwrap();
        assert_eq!(
            (p1, e1.as_slice()),
            (Enzyme::GluC, [Enzyme::Trypsin].as_slice())
        );
        assert_eq!(
            (p2, e2.as_slice()),
            (Enzyme::GluC, [Enzyme::Trypsin].as_slice())
        );
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
