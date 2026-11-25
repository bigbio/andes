//! The search arguments (`andes` with no subcommand): database, search, FDR, rescoring, quant and
//! runtime flags, and the geometry options of the trained models.

use std::path::PathBuf;

use clap::Args;
use model::Tolerance;
use model_train::geometry::GeometryConfig;
use search::decoy::DecoyStrategy;
use search::PrecursorCalMode;

use super::flags::{
    parse_decoy_strategy, parse_enzyme_list, parse_fixed_mods, parse_plex, parse_precursor_cal,
    parse_uniprot_evidence, parse_variable_mods, EnzymeSpecificity, EthcdActivationFlag,
    FragmentIndexFlag, Fragmentation, GlycoCoreFlag, GlycoSpeciesFlag, Listed, ModList,
    PrecursorMonoFlag, PredictFeature, PredictorBackend, Protocol, ScoreFlag, TmtUnmix, XlSet,
};
use super::{
    ADVANCED_FDR, COARSE_QUANT_FDR, DECOY_STRATEGIES, DEFAULT_DECOY_STRATEGY, DEFAULT_FDR, ENZYMES,
    EVIDENCE_LEVELS, H_DATABASE, H_FDR, H_GLYCO, H_IO, H_ISOBARIC, H_LABEL, H_LFQ, H_MODES,
    H_PREDICTION, H_QUANT, H_RESCORING, H_RUNTIME, H_SEARCH, PLEXES, PRECURSOR_CAL_MODES,
    PROTEOFORM_FDR_FLOOR,
};
use crate::cli::parse::{
    parse_charge_range, parse_dia_interference, parse_isotope_error_range, parse_non_negative,
    parse_non_negative_f32, parse_open_search_delta, parse_positive_tol, parse_precursor_tol,
    parse_unit_fraction, parse_xl,
};

/// The protein database options shared by the search, `andes inference` and
/// `andes quant finalize`: how the database is digested and expanded, and how
/// a decoy is recognised. Inference and finalize take the search's values.
#[derive(Args, Debug, Clone)]
pub(crate) struct DatabaseArgs {
    /// Protease; a comma list (gluc,trypsin) cleaves at the sites of each one
    ///
    /// The order of a list does not matter. The primary protease, which selects the model
    /// and the cleavage feature, is trypsin whenever the list has it; otherwise the one whose
    /// rule sets the C-terminal residue of most peptides (a C-terminal protease with more
    /// cleavage residues). A non-specific protease in a list makes the whole digest
    /// non-specific. Also accepted: tryp, chymo, lys-c, asp-n, glu-c, lys-n, arg-c, alp,
    /// alpha-lp, none (nocleavage), all and elastase (nonspecific), in any case.
    #[arg(long, value_name = "ENZYME", default_value = "trypsin", help_heading = H_SEARCH,
          value_parser = Listed::new(ENZYMES, parse_enzyme_list))]
    pub(crate) enzyme: String,

    /// Fixed modifications, a list of Name (sites) [default: Carbamidomethyl (C)]
    ///
    /// Comma-separated, e.g. "Carbamidomethyl (C), TMT6plex (K), TMT6plex (N-term)". The name
    /// is a UNIMOD name or accession (Phospho, UNIMOD:21; GlyGly, Deamidation, TMT10plex and
    /// TMT16plex are also accepted), which gives the mass. The sites are residues (STY),
    /// N-term, C-term, Protein N-term, Protein C-term, or a terminus followed by residues
    /// (N-term Q). The list replaces the default; with the default, a run found to carry a
    /// TMT or iTRAQ tag also gets the tag on K and the peptide N-term. Neutral losses, loss
    /// class, cleavage blocking, diagnostic ions, annotated-only sites and custom masses are
    /// set per modification in the config file (search.fixed_mods). Format and examples:
    /// DOCS.md §2.
    #[arg(long, value_name = "LIST", help_heading = H_SEARCH, value_parser = parse_fixed_mods)]
    pub(crate) fixed_mods: Option<ModList>,

    /// Variable modifications, a list of Name (sites) [default: Oxidation (M), Acetyl (Protein N-term)]
    ///
    /// Written as --fixed-mods, e.g. "Oxidation (M), Phospho (STY), Gln->pyro-Glu (N-term Q)";
    /// "" searches none. The list replaces the default, and --max-mods bounds the variable
    /// modifications per peptide. The per-modification settings are config-file keys
    /// (search.variable_mods). For `andes inference` and `andes quant finalize`, give the
    /// search's lists or its --config: PSMs that place a `restrict_to: annotated`
    /// modification on a knockoff site count for no protein group, as in the search.
    #[arg(long, value_name = "LIST", help_heading = H_SEARCH, value_parser = parse_variable_mods)]
    pub(crate) variable_mods: Option<ModList>,

    /// Accession prefix that marks a decoy, given to the generated decoys
    #[arg(long, value_name = "PREFIX", default_value = "XXX_", hide_short_help = true,
          help_heading = H_SEARCH)]
    pub(crate) decoy_prefix: String,

    /// Accession suffix that also marks a decoy (e.g. rev for <accession>_rev decoys)
    ///
    /// A protein is a decoy when its accession starts with the --decoy-prefix or ends with
    /// this suffix. For a database that holds its decoys, with --decoy-strategy none, so
    /// that andes adds none of its own.
    #[arg(long, value_name = "SUFFIX", hide_short_help = true, help_heading = H_SEARCH)]
    pub(crate) decoy_suffix: Option<String>,

    /// UniProt XML: also search each natural sequence variant as its own entry
    ///
    /// Only the peptides carrying the variant are added, under an accession naming the
    /// change (sp|P04637_VAR_005851_D7H|P53_HUMAN), and decoys mirror them. Grows the
    /// search space.
    #[arg(long, hide_short_help = true, help_heading = H_DATABASE)]
    pub(crate) uniprot_variants: bool,

    /// UniProt XML: also search each described isoform as its own entry
    ///
    /// An isoform (sp|P04637-2|P53_HUMAN) is rebuilt from the canonical sequence and its
    /// splice variants; it carries its own annotated sites, the canonical ones on the
    /// residues it shares and, with --uniprot-variants, its own natural variants. Grows the
    /// search space.
    #[arg(long, hide_short_help = true, help_heading = H_DATABASE)]
    pub(crate) uniprot_isoforms: bool,

    /// UniProt XML: also search each mature chain or peptide whose ends differ from the precursor's
    ///
    /// A chain is what remains after a signal peptide, propeptide or transit peptide is
    /// removed. Only the peptides at the new termini are added; they are protein-terminal
    /// there, so protein-terminal modifications apply. Grows the search space.
    #[arg(long, hide_short_help = true, help_heading = H_DATABASE)]
    pub(crate) uniprot_chains: bool,

    /// UniProt XML: least evidence of a variant, annotated site or chain that is used
    #[arg(long, value_name = "LEVEL", default_value = "experimental", hide_short_help = true,
          help_heading = H_DATABASE,
          value_parser = Listed::new(EVIDENCE_LEVELS, parse_uniprot_evidence))]
    pub(crate) uniprot_evidence: input::EvidenceLevel,

    /// Annotation file (<name>.annotations.parquet) of variants, sites and chains; repeatable
    ///
    /// Claims about the database's reference proteins from any source, one row per claim:
    /// sequence edits become variant entries, chains become chain entries and sites are
    /// added to the reference, as from a UniProt XML database. Several files are merged by
    /// content. Columns: DOCS.md.
    #[arg(long, value_name = "FILE", hide_short_help = true, help_heading = H_DATABASE)]
    pub(crate) annotations: Vec<PathBuf>,

    /// Knockoff entries per single-residue variant, which estimate the variant FDR (0: none)
    ///
    /// For --uniprot-variants and --annotations variants. A knockoff is the same
    /// substitution at another position holding the same residue, in a matched sequence
    /// context and at least 10 residues from every known variant: almost certainly absent
    /// from the sample, it fits a mass coincidence or an isotope error as well as the real
    /// variant. The knockoffs accepted per substitution class estimate how many accepted
    /// variant PSMs are false, and with --rescore variant PSMs get a q-value from that
    /// estimate too. A variant PSM is accepted only when its variant has a knockoff;
    /// knockoff PSMs are never reported. With 0, variant PSMs get the class FDR only.
    #[arg(long, value_name = "N", default_value_t = 4, hide_short_help = true,
          help_heading = H_DATABASE)]
    pub(crate) knockoffs: usize,
}

#[derive(Args, Debug)]
#[command(group(
    clap::ArgGroup::new("glycan_db")
        .args(["glyco_glycan_gdb", "glyco_glycan_wurcs", "glyco_species"])
        .multiple(false)
        .required(false),
))]
pub(crate) struct SearchArgs {
    /// YAML file of options (config.example.yaml); a flag on the command line wins
    ///
    /// Sections and keys: DOCS.md §1b.
    #[arg(long, value_name = "FILE", help_heading = H_IO)]
    pub(crate) config: Option<PathBuf>,

    /// Write the resolved settings as a config file that reproduces the run, then exit
    ///
    /// FILE, or `-` for stdout. Nothing is searched. Replay with --config FILE.
    #[arg(long, value_name = "FILE", help_heading = H_IO, hide_short_help = true)]
    pub(crate) print_config: Option<PathBuf>,

    /// SDRF-Proteomics file: data files, protease, modifications, tolerance and labels
    ///
    /// Reads `comment[data file]` (looked up next to the SDRF), the cleavage agents, the
    /// modification parameters, the precursor mass tolerance and `comment[label]`
    /// (label free -> --lfq, TMT / iTRAQ -> --tmt with the smallest plex holding the
    /// labels, SILAC and dimethyl -> --label). A value on the command line or in --config
    /// wins. Rows must agree on these settings. The fragment tolerance is not applied
    /// (andes measures each run's), and the dissociation method is checked against the
    /// files. DOCS.md §1b, SDRF input.
    #[arg(long, value_name = "FILE", help_heading = H_IO)]
    pub(crate) sdrf: Option<PathBuf>,

    /// The SDRF dissociation method, checked against the files (set from --sdrf).
    #[arg(skip)]
    pub(crate) sdrf_dissociation: Option<model::activation::ActivationMethod>,

    /// Spectrum file: mzML, MGF, Thermo .raw or Bruker .d; repeat for several runs
    ///
    /// The reader follows the file extension. Several files are searched into one PIN.
    #[arg(long, value_name = "FILE", help_heading = H_IO)]
    pub(crate) spectrum: Vec<PathBuf>,

    /// Protein database: FASTA, or UniProt XML (.xml); either may be gzipped
    ///
    /// Target sequences only: andes generates the decoys (--decoy-strategy), unless it is the
    /// target + decoy FASTA of `andes database`, whose decoys are searched as written. A UniProt XML
    /// database searches like the UniProt FASTA and also carries annotated modification
    /// sites (used by `restrict_to: annotated` modifications) and the entries of the Database
    /// variants options.
    #[arg(long, value_name = "FILE", help_heading = H_IO)]
    pub(crate) database: Option<PathBuf>,

    /// Percolator PIN output (required unless --rescore is given)
    #[arg(long, value_name = "FILE", help_heading = H_IO)]
    pub(crate) output_pin: Option<PathBuf>,

    /// Tab-separated PSM report
    #[arg(long, value_name = "FILE", help_heading = H_IO)]
    pub(crate) output_tsv: Option<PathBuf>,

    /// QPX .idparquet bundle directory (quantms-compatible)
    ///
    /// Writes psms.parquet, proteins.parquet, protein_groups.parquet (empty) and
    /// search_params.parquet, and with --tmt, --lfq or --label quantms.feature.parquet.
    #[arg(long, value_name = "DIR", help_heading = H_IO)]
    pub(crate) output_parquet: Option<PathBuf>,

    /// Precursor mass tolerance in ppm or Da, e.g. 20ppm or 0.02da
    #[arg(long, value_name = "TOL", default_value = "20ppm", value_parser = parse_precursor_tol,
          help_heading = H_SEARCH)]
    pub(crate) precursor_tol: Tolerance,

    #[command(flatten)]
    pub(crate) db: DatabaseArgs,

    /// How decoys are generated
    ///
    /// Decoys are peptide-level: each segment between two cleavage sites is permuted with
    /// its cleavage residues, its cleavage-side terminal residue and an initiator Met in
    /// place, so every enzymatic peptide has a decoy of the same composition and C-terminal
    /// residue (a non-specific search permutes whole sequences). Every decoy peptide the
    /// search enumerates (missed-cleavage joins, semi-specific and non-specific sub-spans
    /// included) is checked against the target peptides with I and L folded; one that
    /// equals a target is redrawn, and the few that still do are scored as targets. With
    /// none and a target-only database the search has no decoys and warns.
    #[arg(long, value_name = "STRATEGY", default_value = DEFAULT_DECOY_STRATEGY, help_heading = H_SEARCH,
          value_parser = Listed::new(DECOY_STRATEGIES, parse_decoy_strategy))]
    pub(crate) decoy_strategy: DecoyStrategy,

    /// Worker threads [default: the number of logical CPUs]
    #[arg(long, value_name = "N", default_value_t = num_cpus::get(), hide_default_value = true,
          help_heading = H_SEARCH)]
    pub(crate) threads: usize,

    /// Precursor charges tried when the spectrum records none, as MIN..MAX (or MIN-MAX)
    #[arg(long, value_name = "MIN..MAX", default_value = "2..5", value_parser = parse_charge_range,
          hide_short_help = true, help_heading = H_SEARCH)]
    pub(crate) charge: (u8, u8),

    /// Isotope-error offsets tried, as MIN..MAX (or MIN-MAX) [default: -1..2, 0..2 with --glyco or --xl]
    #[arg(long, value_name = "MIN..MAX", value_parser = parse_isotope_error_range,
          hide_short_help = true, help_heading = H_SEARCH)]
    pub(crate) isotope_error: Option<(i8, i8)>,

    /// Open (mass-offset) search: any peptide within --open-search-delta of the precursor
    #[arg(long, hide_short_help = true, help_heading = H_SEARCH)]
    pub(crate) open_search: bool,

    /// Mass-offset window for --open-search, as MIN..MAX in Da [default: -100..300]
    #[arg(long, value_name = "MIN..MAX", default_value = "-100..300", allow_hyphen_values = true,
          value_parser = parse_open_search_delta, hide_short_help = true, help_heading = H_SEARCH)]
    pub(crate) open_search_delta: (f64, f64),

    /// Candidates scored per spectrum under --open-search [default: 100]
    #[arg(long, value_name = "N", default_value = "100", hide_short_help = true,
          help_heading = H_SEARCH)]
    pub(crate) open_search_top_k: u32,

    /// Parameter search: grid precursor tolerance, top-N and charge on a sample, then run
    #[arg(long, hide_short_help = true, help_heading = H_SEARCH)]
    pub(crate) optimize: bool,

    /// Spectra the --optimize grid scores [default: 5000]
    #[arg(long, value_name = "N", default_value = "5000", hide = true)]
    pub(crate) optimize_sample: usize,

    /// Missed cleavages per peptide (not applied to a non-specific digest)
    ///
    /// --glyco raises the value to 3 unless it is given on the command line or in --config.
    #[arg(long, value_name = "N", default_value = "2", hide_short_help = true,
          help_heading = H_SEARCH)]
    pub(crate) max_missed_cleavages: u32,

    /// Peptide ends that must be cleavage sites
    #[arg(long, value_name = "SPECIFICITY", value_enum, ignore_case = true,
          default_value = "fully", hide_short_help = true, help_heading = H_SEARCH)]
    pub(crate) enzyme_specificity: EnzymeSpecificity,

    /// Shortest peptide, in residues
    #[arg(long, value_name = "N", default_value = "6", hide_short_help = true,
          help_heading = H_SEARCH)]
    pub(crate) min_length: u32,

    /// Longest peptide, in residues
    #[arg(long, value_name = "N", default_value = "50", hide_short_help = true,
          help_heading = H_SEARCH)]
    pub(crate) max_length: u32,

    /// Variable modifications per peptide
    #[arg(long, value_name = "N", default_value = "3", hide_short_help = true,
          help_heading = H_SEARCH)]
    pub(crate) max_mods: u32,

    /// Fewest peaks of a scored MS2 spectrum; spectra with fewer are skipped
    #[arg(long, value_name = "N", default_value = "10", hide_short_help = true,
          help_heading = H_SEARCH)]
    pub(crate) min_peaks: u32,

    /// PSMs kept per spectrum
    #[arg(long, value_name = "N", default_value = "10", hide_short_help = true,
          help_heading = H_SEARCH)]
    pub(crate) top_n: u32,

    /// Seed of the decoy shuffles (shuffle, sequon, and the redraw of a segment
    /// whose decoy reads as a target peptide).
    #[arg(long, hide = true, default_value_t = search::decoy::DEFAULT_DECOY_SEED)]
    pub(crate) decoy_seed: u64,

    /// PSM-level FDR: the q-value threshold of the reported PSMs [default: 0.01]
    ///
    /// It sets the filtered <stem>.q<fdr>.tsv output (target PSMs at q <= this), the PSMs
    /// that quantification measures, the match-between-runs donors and the confident
    /// anchors of --refine. andes computes q-values only with --rescore; without it,
    /// setting it warns and writes no filtered output.
    #[arg(long, value_name = "FDR", value_parser = parse_unit_fraction, help_heading = H_FDR)]
    pub(crate) psm_fdr: Option<f64>,

    /// Protein-group q-value threshold of the filtered output [default: none]
    ///
    /// With --rescore, protein groups and their q-values (<stem>.protein_groups.tsv, picked
    /// target-decoy competition per protein family) are written either way; with this
    /// threshold the filtered TSV keeps only PSMs whose peptide counts for a group at or
    /// below it, and the proteoform table only the families of such groups. With --mbr and
    /// --output-parquet, the QPX pg and feature views keep the groups at or below it, 0.01
    /// when it is not given.
    #[arg(long, value_name = "FDR", value_parser = parse_unit_fraction, help_heading = H_FDR)]
    pub(crate) protein_fdr: Option<f64>,

    /// Quantification FDR: q-value threshold of values measured without an identification [default: 0.01]
    ///
    /// It covers the match-between-runs transfers of --lfq, --tmt, --label and --glyco
    /// --lfq, the --label channels matched without an identification and, with --lfq, the
    /// feature q-value of `accepted` in lfq.parquet. A quantification-quality q-value, not
    /// an identification FDR. When neither this nor the value's advanced flag is given, the
    /// --label channels and transfers, the lfq.parquet `accepted` flag and the --glyco --lfq
    /// --mbr transfers use 0.05: their q-values come from one target/decoy competition per
    /// run with few decoys.
    #[arg(long, value_name = "FDR", value_parser = parse_unit_fraction, help_heading = H_FDR)]
    pub(crate) quant_fdr: Option<f64>,

    /// Rescore the PIN with Percolator and write q-values, PEPs and the PSMs at --psm-fdr
    ///
    /// Percolator's PEP and q-value are joined into the outputs (QPX
    /// posterior_error_probability and a q-value score) and the filtered <stem>.q<fdr>.tsv
    /// is written. Every Percolator pass (the main run, the per-run q-values, the PEP
    /// repair and the scores of every PIN row) runs in process with the native engine, which
    /// implements the Percolator algorithm (Käll et al., Nat. Methods 2007), tracks upstream
    /// Percolator master (commit a9783ee8) and writes Percolator's output files; no Percolator binary or Docker is needed. --percolator-bin or
    /// --percolator-docker runs an external Percolator instead. --xl does not use the native
    /// engine: its Percolator passes always run an external Percolator, --percolator-bin,
    /// else `percolator` on PATH, else the Docker image.
    #[arg(long, help_heading = H_FDR)]
    pub(crate) rescore: bool,

    /// Also search co-isolated peptides: a second PSM per scan from its MS1 isolation window
    ///
    /// Pass 1 is the normal search and keeps the same --top-n candidates per scan; pass 2
    /// takes the best one as the primary peptide, detects co-isolated precursors in the
    /// scan's MS1 isolation window and searches the residual spectrum for a second peptide,
    /// written as an extra PSM. Needs MS1 scans (mzML, .raw); no effect on MGF input.
    #[arg(long, help_heading = H_MODES)]
    pub(crate) chimeric: bool,

    /// Data-independent acquisition: search every MS2 scan across its isolation window
    ///
    /// Each scan is searched against every peptide whose precursor m/z lies in its
    /// isolation window and has an MS1 peak (spectrum-centric DIA, as MSFragger-DIA); the PIN
    /// gets fragment and MS1 co-elution columns from the neighbouring cycles of the same
    /// window. With --rescore, Percolator runs with --trainFDR 0.05 and again after a per-run
    /// RT refinement. On Orbitrap Astral the fragment tolerance is 7 ppm unless
    /// --fragment-tol-* is set. Needs mzML or .raw with isolation windows and MS1 scans;
    /// replaces --chimeric and --precursor-cal; not with --open-search, --refine, --glyco or --xl.
    #[arg(long, help_heading = H_MODES)]
    pub(crate) dia: bool,

    /// diaPASEF (.d) only: split each isolation window's mobility range into N slices
    ///
    /// Each slice spans two N+1ths of the window's TIMS scan range, and neighbouring slices
    /// overlap by half. Each slice is searched as its own spectrum and traced on its own across
    /// cycles. A slice spectrum mixes fewer co-isolated peptides than the whole window. 1 keeps
    /// one spectrum per window. Ignored for mzML, .raw and ddaPASEF.
    #[arg(long, value_name = "N", default_value = "1", hide_short_help = true,
          help_heading = H_MODES)]
    pub(crate) dia_mobility_slices: usize,

    /// --dia MS1 gate: ppm tolerance for a candidate's monoisotopic peak in the survey scans
    ///
    /// A DIA candidate is scored only when its monoisotopic m/z has an MS1 peak within this
    /// tolerance in the survey scan nearest the MS2 scan or one either side. timsTOF MS1 is less
    /// accurate than Orbitrap MS1, so diaPASEF runs need a wider gate (8 on the ProteoBench
    /// low-input runs, after andes recentres their MS1 m/z).
    #[arg(long, value_name = "PPM", default_value = "5", hide_short_help = true,
          help_heading = H_MODES)]
    pub(crate) dia_ms1_tol_ppm: f64,

    /// diaPASEF MS1 gate in mobility (experiment): 1/K0 margin around the scan's mobility slice
    ///
    /// With a value, a candidate passes the MS1 gate only when its MS1 peak lies inside the
    /// MS2 scan's mobility slice widened by this margin (Vs/cm2) on each side; `auto` learns
    /// the margin from the run's RT calibration anchors.
    #[arg(long, value_name = "K0|auto", hide = true, value_parser = parse_mobility_margin)]
    pub(crate) dia_ms1_mobility_margin: Option<f64>,

    /// --dia --rescore: drop precursors whose fragments a stronger co-eluting precursor explains
    ///
    /// `F:DRT`, e.g. `0.4:4`, or `auto`. Before the precursor-level q, a precursor leaves the
    /// target-decoy competition (its rows get q = 1) when a higher-scoring precursor of any label
    /// has its best row in the same isolation window and an overlapping mobility slice within DRT
    /// seconds and holds at least F of its b/y fragment m/z. Only precursors at precursor q <= 0.2
    /// take part, on both sides. Labels are not used. `auto` estimates F and DRT from each run:
    /// pairs 60-200 s apart give the chance fragment share, DRT is where the excess share of closer
    /// pairs fades, and F is the smallest
    /// share whose chance pairs are at most 1% of those observed within DRT. When co-eluting pairs
    /// share fragments little more than chance, `auto` removes nothing. Off by default; measured
    /// on diaPASEF (docs/benchmarks/dia/README.md).
    #[arg(long, value_name = "F:DRT|auto", value_parser = parse_dia_interference,
          hide = true, help_heading = H_MODES)]
    pub(crate) dia_interference: Option<crate::dia::Interference>,

    /// Second pass for further modifications (oxidation, deamidation, ...) over the confident proteins
    ///
    /// The modifications are a built-in tier, or the tier of --refine-config.
    #[arg(long, help_heading = H_MODES)]
    pub(crate) refine: bool,

    /// Site localization of variable modifications (`ipf`); needs --rescore
    ///
    /// `ipf`: inference of peptidoforms (the OpenSWATH/PyProphet hierarchical Bayesian
    /// model) over every placement of each identified PSM's variable mods, after
    /// rescoring. Writes <base>.ipf.tsv, and the site probabilities into the QPX
    /// psms.parquet with --output-parquet. Experimental.
    #[arg(long = "localize", value_name = "METHOD", value_parser = ["ipf"],
          hide_short_help = true, help_heading = H_MODES)]
    pub(crate) localize: Option<String>,

    /// With --localize ipf: add a forms_debug column to <base>.ipf.tsv (diagnostics)
    ///
    /// Every placement as `sites|search score|IPF log-evidence`.
    #[arg(long = "ipf-debug", hide = true, requires = "localize", help_heading = H_MODES)]
    pub(crate) ipf_debug: bool,

    /// Intact glycopeptide search into a .glyco.pin (needs a glycan database, e.g. --glyco-species)
    ///
    /// Glycans come from the glycan database (--glyco-species, --glyco-glycan-gdb or
    /// --glyco-glycan-wurcs); the backbone is solved as precursor minus glycan and filtered
    /// by the N-X-S/T sequon (any S/T with --glyco-core o). With --lfq the glycopeptides
    /// are quantified from MS1 as intact glycopeptides (peptide, glycan composition and
    /// charge) into <base>.glyco.lfq.tsv (per feature), .glyco.quant.tsv (per
    /// glycopeptide), .glyco.site.tsv (glycoform shares per site), .glyco.protein.tsv and
    /// the cross-run .glyco.matrix.tsv and .glyco.protein_matrix.tsv; --lfq-tol,
    /// --lfq-rt-window, --lfq-min-cosine and --quant-transfer-fdr apply, and --mbr is the
    /// glycopeptide match between runs. The isobaric options do not apply.
    #[arg(long, help_heading = H_MODES)]
    pub(crate) glyco: bool,

    /// Crosslink search into <base>.xl.tsv, linker as Name (sites), e.g. "DSSO"; several as "DSSO; DSBU" (experimental)
    ///
    /// The linker is written like a modification: a reagent name (DSSO, DSBU or BuUrBu,
    /// DSBSO, CDI; non-cleavable DSS or BS3, BS2G or DSG, EGS, DST, DSP or DTSSP, PhoX or
    /// DSPP, DMP, DTBP, DFDNB, BMOE, SMCC, SDA, EDC or DMTMM), its XLMOD accession
    /// (XLMOD:02126), or the UNIMOD name or accession of its intact cross-link
    /// (Xlink:DSS[138], UNIMOD:1898), then optionally
    /// its sites in parentheses: "DSSO (K, S, T, Y, Protein N-term)". Without sites the
    /// reagent's usual ones apply (Lys and the protein N-terminus for NHS esters); N-term
    /// adds every peptide N-terminus (synthetic peptide libraries), and a
    /// heterobifunctional reagent separates its two ends with /. The reagents and their
    /// masses (UNIMOD records, or cited compositions for reagents UNIMOD lacks) are the
    /// entries of resources/crosslinkers.yaml; a reagent the table lacks is an entry of the
    /// same schema in the config file (search.xl). An --sdrf comment[cross-linker] column
    /// sets the linker when neither the command line nor the config does. Monolinks are
    /// variable modifications, by UNIMOD
    /// record or as reagent:quencher (--variable-mods "Xlink:DSSO[176] (K), DSBSO:water (K)").
    ///
    /// For an MS-cleavable linker each MS2 is searched for the reporter doublets; the doublet
    /// peptide masses and the precursor give the candidate peptide pairs. For a non-cleavable
    /// linker the fragment-ion index retrieves one peptide by its fragments that do not
    /// carry the link, and its partner is looked up by the remaining precursor mass. Each
    /// peptide is scored with the standard scorer carrying the cleaved stub or its partner on
    /// the link residue.
    /// Writes the best crosslink-spectrum match per spectrum into a pair-aware PIN
    /// (<base>.xl.pin; both link classes are scored by the pair score, and --rescore does
    /// not change them) and the q-values, from the
    /// target-target, target-decoy and decoy-decoy counts split by intra-protein, homotypic
    /// (one peptide sequence on both sides, or two peptides linked at the same or the next
    /// residue of one protein) and inter-protein links and pooled over the runs,
    /// of the matches (<base>.xl.tsv), residue pairs
    /// (<base>.xl.residue_pairs.tsv) and protein pairs (<base>.xl.ppi.tsv), with the
    /// doublet summary of every spectrum (<base>.xl.doublets.tsv).
    #[arg(long, value_name = "LINKER", value_parser = parse_xl, hide_short_help = true,
          help_heading = H_MODES, conflicts_with = "glyco")]
    pub(crate) xl: Option<XlSet>,

    /// With --xl: also run the exhaustive pair search (no doublet gate) on the first N spectra
    #[arg(
        long = "xl-oracle",
        value_name = "N",
        default_value = "0",
        hide = true,
        requires = "xl"
    )]
    pub(crate) xl_oracle: usize,

    /// With --xl: search only linear and monolinked peptides and write their protein groups
    ///
    /// The search runs with the variable modifications --xl implies (the default ones plus
    /// the linker's water-quenched monolink) and no cross-link search; with --rescore it
    /// writes <stem>.protein_groups.tsv.
    #[arg(
        long = "xl-linear-only",
        requires = "xl",
        hide_short_help = true,
        help_heading = H_MODES
    )]
    pub(crate) xl_linear_only: bool,

    /// With --xl: rescore the intra-link CSMs with Percolator (comparison only)
    ///
    /// Trains a Percolator model on the intra-link rows alone (<base>.xl.intra.pin) and
    /// scores them with it; without it the intra-links keep the pair score.
    #[arg(long = "xl-rescore-intra", hide = true, requires = "xl")]
    pub(crate) xl_rescore_intra: bool,

    /// With --xl: also rescore the inter-link CSMs with Percolator (comparison only)
    ///
    /// Trains a Percolator model on the inter-link rows alone (<base>.xl.inter.pin) and
    /// scores them with it; without it the inter-links keep the pair score.
    #[arg(long = "xl-rescore-inter", hide = true, requires = "xl")]
    pub(crate) xl_rescore_inter: bool,

    /// Bundled N-glycan database of --glyco
    ///
    /// The databases and their attribution: crates/glyco/glycan-db/.
    #[arg(long, value_name = "SPECIES", value_enum, help_heading = H_MODES)]
    pub(crate) glyco_species: Option<GlycoSpeciesFlag>,

    /// Most co-isolated peptides searched per scan with --chimeric
    #[arg(long, value_name = "N", default_value = "4", hide_short_help = true,
          help_heading = H_MODES)]
    pub(crate) chimeric_max_coisolated: usize,

    /// Largest averagine-envelope KL divergence of an accepted co-isolated precursor (--chimeric)
    ///
    /// Lower is stricter: fewer spurious second peptides.
    #[arg(long, value_name = "KL", default_value = "0.3", value_parser = parse_non_negative_f32,
          hide_short_help = true, help_heading = H_MODES)]
    pub(crate) chimeric_max_kl: f32,

    /// YAML file of the --refine modification tier [default: the built-in tier]
    #[arg(long, value_name = "FILE", hide_short_help = true, help_heading = H_MODES)]
    pub(crate) refine_config: Option<PathBuf>,

    /// Isobaric reporter-ion quantification with this labelling kit
    ///
    /// Writes <stem>.tmt.parquet and, with --output-parquet, the QPX
    /// quantms.feature.parquet. With --rescore only the PSMs at q <= --psm-fdr are
    /// quantified. Also accepted: the plex spellings (tmt10plex) and tmtpro, tmtpro16,
    /// tmtpro18, tmtpro32, tmtpro35, in any case.
    #[arg(long, value_name = "PLEX", help_heading = H_QUANT,
          value_parser = Listed::new(PLEXES, parse_plex))]
    pub(crate) tmt: Option<quant::Plex>,

    /// Label-free MS1 quantification of the identified precursors
    ///
    /// Isotope chromatograms are extracted around each identified peptidoform's retention
    /// time and integrated. Writes <stem>.lfq.parquet and, with --output-parquet, the QPX
    /// quantms.feature.parquet. mzML, Thermo .raw and Bruker .d input (a .d's MS1 frames
    /// are summed over ion mobility). With --glyco it quantifies intact glycopeptides
    /// instead (see --glyco).
    #[arg(long, help_heading = H_QUANT)]
    pub(crate) lfq: bool,

    /// MS1 labelled-pair quantification (SILAC, dimethyl): a preset or one bracket per channel
    ///
    /// Presets: silac2, silac3, dimethyl2, dimethyl3. Channels are brackets in increasing
    /// mass order, as OpenMS FeatureFinderMultiplex writes them ([][Lys8,Arg10],
    /// [Dimethyl0][Dimethyl4][Dimethyl8]); a label is one of Lys4 Lys6 Lys8 Arg6 Arg10 Leu3
    /// Dimethyl0 Dimethyl4 Dimethyl6 Dimethyl8, or <site>+<mass>[:[chem]<n>C<n>N<n>D] with
    /// the label's 13C, 15N and 2H atoms (chem: a chemical label that adds atoms; without
    /// atoms: a chemical shift of natural atoms). The labels are searched as variable
    /// modifications that every labelable site of a peptide carries alike, and each
    /// identified peptide is measured in every channel by one joint fit of their MS1
    /// isotope traces; channels without an identification are matched at their exact
    /// position (match between channels), and with --mbr peptides are transferred between
    /// runs. Writes <stem>.label.parquet. mzML or .raw with high-resolution MS1; not with
    /// --lfq or --tmt.
    #[arg(long, value_name = "SPEC", help_heading = H_QUANT)]
    pub(crate) label: Option<String>,

    /// Match between runs across two or more inputs (needs --rescore, and --lfq, --tmt or --label)
    ///
    /// Label-free: after the per-run quantification, builds the precursor library, aligns
    /// retention times across runs, measures every library precursor in every run and
    /// writes <base>.mbr.parquet and <base>.mbr_validation.parquet (with --output-parquet,
    /// quantms.feature.parquet holds these rows instead of the per-run features, without
    /// the transfers above the threshold). Uses --quant-dir for the per-run tables, or a
    /// temporary directory. With --rescore (Percolator) it gives the same result as `andes
    /// quant library`, `extract` and `finalize` on the same files and Percolator results.
    ///
    /// With --tmt: the reporter ions of unidentified MS2 scans at an identified peptide's
    /// m/z and aligned retention time are transferred into the runs that lack it, each with
    /// an evidence type and a quantification-quality q-value (`quant_q`) in tmt.parquet and
    /// the QPX feature view; high-resolution MS2 reporters only (not SPS-MS3 or ion-trap
    /// MS2); identified rows are unchanged.
    ///
    /// With --glyco --lfq: every accepted glycopeptide is measured in each run that did not
    /// identify it, around its predicted apex on a retention-time map fitted to the first
    /// input on the glycopeptides both runs identified, and the transfers are scored in
    /// their own target/decoy competition. Transferred rows say so in the Evidence column of
    /// .glyco.lfq.tsv. Details: DOCS.md §9.
    #[arg(long, help_heading = H_QUANT)]
    pub(crate) mbr: bool,

    /// Directory of the per-run tables of multi-run quantification
    ///
    /// For each input: <stem>.ms1.parquet (its MS1 scans, every FAIMS voltage),
    /// <stem>.ms2.parquet (its MS2 scans, reduced) and <stem>.quant.parquet (every rank-1
    /// PSM, targets and decoys, unfiltered; with --chimeric also the best PSM of each
    /// co-isolated precursor). Input as for --lfq (mzML, .raw or .d); not with MGF or
    /// --glyco.
    #[arg(long, value_name = "DIR", help_heading = H_QUANT)]
    pub(crate) quant_dir: Option<PathBuf>,

    /// Fixed storage budget of the native Percolator engine in GB, not a process-memory cap: PIN features, PSM texts, normalised training rows and score buffers spill to --temp-dir once they would take more than their share of it; per-PSM metadata (about 160 bytes a PSM) is charged to it but always stays in memory [default: what the PIN needs to stay in memory, at most 16 GB and at most the memory free to the job or machine less a margin]
    #[arg(long, value_name = "GB",
          value_parser = clap::value_parser!(u64).range(1..), hide = true,
          help_heading = H_RESCORING)]
    pub(crate) rescore_memory_gb: Option<u64>,

    /// PEP threshold applied with --psm-fdr: a PSM must pass both [default: none]
    #[arg(long, value_name = "PEP", value_parser = parse_unit_fraction, hide_short_help = true,
          help_heading = H_RESCORING)]
    pub(crate) pep: Option<f64>,

    /// Run --rescore with this external Percolator binary instead of the native engine
    ///
    /// With --xl the default is `percolator` on PATH, else the Docker image.
    #[arg(long, value_name = "FILE", hide_short_help = true, help_heading = H_RESCORING)]
    pub(crate) percolator_bin: Option<PathBuf>,

    /// Run --rescore with Percolator from its Docker image instead of the native engine
    ///
    /// Needs the docker command.
    #[arg(long, hide_short_help = true, help_heading = H_RESCORING)]
    pub(crate) percolator_docker: bool,

    /// Docker image of Percolator, with --percolator-docker or --xl
    #[arg(long, value_name = "IMAGE", default_value = output::DEFAULT_PERCOLATOR_IMAGE,
          hide_short_help = true, help_heading = H_RESCORING)]
    pub(crate) percolator_image: String,

    /// Extra arguments for an external Percolator, e.g. "--testFDR 0.05"
    ///
    /// Only with --percolator-bin, --percolator-docker or --xl; the native engine takes none.
    #[arg(long, value_name = "ARGS", default_value = "", hide_default_value = true,
          hide_short_help = true, help_heading = H_RESCORING)]
    pub(crate) percolator_args: String,

    /// Keep the PIN after rescoring
    ///
    /// Without --output-pin, true writes the PIN (and Percolator's result files) as
    /// <base>.pin, <base> being --output-tsv, else --output-parquet, else ./andes; false
    /// uses a temporary PIN deleted at exit. An explicit --output-pin is always kept.
    #[arg(long, value_name = "BOOL", default_value_t = true, action = clap::ArgAction::Set,
          hide_short_help = true, help_heading = H_RESCORING)]
    pub(crate) keep_pin: bool,

    /// Precursor mass calibration pre-pass
    ///
    /// auto and on learn a systematic ppm shift from confident PSMs in a pre-pass and
    /// tighten the precursor tolerance of the main search. Not available for Thermo .raw
    /// input (warns).
    #[arg(long, value_name = "MODE", default_value = "off", hide_short_help = true,
          help_heading = H_RUNTIME,
          value_parser = Listed::new(PRECURSOR_CAL_MODES, parse_precursor_cal))]
    pub(crate) precursor_cal: PrecursorCalMode,

    /// Fewest (spectrum, charge) keys for the calibration pre-pass to run; a
    /// test hook that forces the pre-pass on small fixtures. [default: 10000]
    #[arg(long, hide = true)]
    pub(crate) cal_min_spec_keys: Option<usize>,

    /// What ranks the candidates and fills the PIN RawScore
    #[arg(long, value_name = "SCORE", default_value = "auto", hide_short_help = true,
          help_heading = H_RUNTIME)]
    pub(crate) score: ScoreFlag,

    /// Fragmentation of MGF spectra (mzML, .raw and .d record it)
    #[arg(long, value_name = "METHOD", value_enum, ignore_case = true, default_value = "auto",
          hide_short_help = true, help_heading = H_RUNTIME)]
    pub(crate) fragmentation: Fragmentation,

    /// Labelling or enrichment protocol, which selects the protocol model
    ///
    /// auto keeps the selected model; when it finds TMT or iTRAQ reporter ions it turns on
    /// the isobaric peak filter and, without --fixed-mods, adds the tag as a fixed modification.
    #[arg(long, value_name = "PROTOCOL", value_enum, ignore_case = true,
          default_value = "auto", hide_short_help = true, help_heading = H_RUNTIME)]
    pub(crate) protocol: Protocol,

    /// Fragment tolerance in ppm of the PIN features, replacing the per-run measurement
    ///
    /// On MGF it also selects the high-resolution model class. Not with --fragment-tol-da.
    #[arg(long, value_name = "PPM", conflicts_with = "fragment_tol_da",
          value_parser = parse_positive_tol, hide_short_help = true, help_heading = H_RUNTIME)]
    pub(crate) fragment_tol_ppm: Option<f64>,

    /// Fragment tolerance in Da of the PIN features, replacing the per-run measurement
    ///
    /// On MGF it also selects the low-resolution model class. Not with --fragment-tol-ppm.
    #[arg(long, value_name = "DA", conflicts_with = "fragment_tol_ppm",
          value_parser = parse_positive_tol, hide_short_help = true, help_heading = H_RUNTIME)]
    pub(crate) fragment_tol_da: Option<f64>,

    /// Peak filter WINDOW_DA:PEAKS, e.g. 100:20 [default: 100:20 for isobaric data, else off]
    ///
    /// Keeps the PEAKS most intense peaks in every WINDOW_DA window; a window of 0 turns
    /// the filter off.
    #[arg(long, value_name = "WINDOW_DA:PEAKS", hide_short_help = true, help_heading = H_RUNTIME)]
    pub(crate) peak_filter: Option<String>,

    /// How EThcD and ETciD spectra are labelled
    #[arg(long, value_name = "LABEL", value_enum, default_value_t = EthcdActivationFlag::Ethcd,
          hide_short_help = true, help_heading = H_RUNTIME)]
    pub(crate) ethcd_activation: EthcdActivationFlag,

    /// MS level searched in mzML (MGF, .raw, .d and --chimeric search MS2)
    #[arg(long, value_name = "LEVEL", default_value = "2", hide_short_help = true,
          help_heading = H_RUNTIME)]
    pub(crate) ms_level: u8,

    /// Candidate retrieval
    ///
    /// At a low-resolution fragment tolerance only the 150 most intense peaks vote in the
    /// fragment-ion index. The search log reports the retrieval used.
    #[arg(long, value_name = "MODE", default_value = "auto", hide_short_help = true,
          help_heading = H_RUNTIME)]
    pub(crate) fragment_index: FragmentIndexFlag,

    /// Precursor-mass width of one fragment-index slice, in Da [default: from the memory budget]
    #[arg(long, value_name = "DA", hide_short_help = true, help_heading = H_RUNTIME)]
    pub(crate) fragment_index_slice_da: Option<f64>,

    /// Directory of the candidate-index cache [default: the directory of the outputs]
    ///
    /// The cache (andes-candidx-<key>.bin) is created if missing. A shared directory builds
    /// the index once for several tasks, a node-local one gives a cache per node. A
    /// directory that cannot be written is an error.
    #[arg(long, value_name = "DIR", hide_short_help = true, help_heading = H_RUNTIME)]
    pub(crate) index_cache_dir: Option<PathBuf>,

    /// Directory of the run's temporary files [default: the directory of the outputs]
    ///
    /// Created if missing and emptied at exit. It holds the MS1 scans kept between
    /// searching and quantification, the directory of --mbr without --quant-dir, and the
    /// PIN of --rescore with --keep-pin false. A directory that cannot be written is an
    /// error.
    #[arg(long, value_name = "DIR", hide_short_help = true, help_heading = H_RUNTIME)]
    pub(crate) temp_dir: Option<PathBuf>,

    /// Parquet model store (a file or a partitioned directory) instead of the bundled one
    #[arg(long, value_name = "PATH", hide_short_help = true, help_heading = H_RUNTIME)]
    pub(crate) model_store: Option<PathBuf>,

    /// Model id to load from the model store, skipping the automatic selection
    #[arg(long = "model", value_name = "ID", hide_short_help = true, help_heading = H_RUNTIME)]
    pub(crate) model_id_override: Option<String>,

    /// Trees evaluated per GBDT ensemble, 0 for all [default: 100; all with --glyco]
    ///
    /// Applies to the fragment-intensity and rich-ion ensembles. Fewer trees change the
    /// predicted intensities, so the PIN feature values and, with --score strong, which
    /// candidate ranks first.
    #[arg(long, value_name = "N", hide_short_help = true, help_heading = H_RUNTIME)]
    pub(crate) gbdt_max_trees: Option<usize>,

    /// Search scoring without the fragment-intensity and rich-ion GBDT ensembles
    ///
    /// Experimental speed switch: the search's intensity features (IntensitySignal,
    /// FragPred*, RichIonLLR) become 0 and the strong score ranks without them.
    #[arg(long, hide = true)]
    pub(crate) no_search_gbdt: bool,

    /// Search scoring without the rich-ion GBDT ensemble (experiment): RichIonLLR becomes 0
    #[arg(long, hide = true)]
    pub(crate) no_search_rich_ion: bool,

    /// Process only the first N spectra (0: all); a benchmarking cap.
    #[arg(long, hide = true, default_value = "0")]
    pub(crate) max_spectra: usize,

    /// Predicted features added to the PIN, comma-separated
    ///
    /// Predicted retention times are calibrated to the observed ones per input file.
    /// peptdeep-ms2, peptdeep-rt and deeplc-rt load exported ONNX models from
    /// --predictor-models and fail at start-up when a model file is missing.
    #[arg(long, value_name = "FEATURES", value_delimiter = ',', default_value = "rt",
          hide_short_help = true, help_heading = H_PREDICTION)]
    pub(crate) predict: Vec<PredictFeature>,

    /// Directory of the exported AlphaPeptDeep and DeepLC models
    ///
    /// peptdeep_mods.json, peptdeep_{ms2,rt}_dynamic.{json,onnx} and deeplc_rt.{json,onnx},
    /// as the onnx-infer export scripts write them.
    #[arg(long, value_name = "DIR", hide_short_help = true, help_heading = H_PREDICTION)]
    pub(crate) predictor_models: Option<PathBuf>,

    /// ONNX engine of the AlphaPeptDeep and DeepLC models
    #[arg(long, value_name = "ENGINE", default_value = "tract", hide_short_help = true,
          help_heading = H_PREDICTION)]
    pub(crate) predictor_backend: PredictorBackend,

    /// Collision energy given to the MS2 model [default: 28 on Q Exactive, else 30]
    #[arg(long, value_name = "NCE", hide_short_help = true, help_heading = H_PREDICTION)]
    pub(crate) predictor_nce: Option<f32>,

    /// Instrument given to the MS2 model, e.g. QE, Lumos, timsTOF [default: from the instrument]
    ///
    /// One of the model's instrument names. Unset: QE for a Q Exactive, timsTOF for a
    /// timsTOF, SciexTOF for another TOF, otherwise Lumos.
    #[arg(long, value_name = "NAME", hide_short_help = true, help_heading = H_PREDICTION)]
    pub(crate) predictor_instrument: Option<String>,

    /// Peptides per predictor batch
    #[arg(long, value_name = "N", default_value_t = 256, hide_short_help = true,
          help_heading = H_PREDICTION)]
    pub(crate) predictor_batch: usize,

    /// MS level of the reporter ions: 2 (the identified MS2) or 3 (its SPS-MS3 scan) [default: 3 for SPS-MS3 data, else 2]
    ///
    /// 3 is the default when the MS2 is low-resolution, cannot separate the plex's channels
    /// and the MS3 scans carry reporters. MS3 needs mzML or .raw (a .d has no MS3 scans).
    #[arg(long, value_name = "LEVEL", hide_short_help = true, help_heading = H_ISOBARIC)]
    pub(crate) tmt_level: Option<u8>,

    /// Reporter-ion m/z tolerance, e.g. 20ppm or 0.003da [default: 20ppm (10ppm for tmt32/tmt35), 0.3da on ion-trap MS2]
    ///
    /// 20 ppm (10 ppm for tmt32 and tmt35) applies to high-resolution MS2 and to
    /// --tmt-level 3, whose reporters are read in the Orbitrap. A window reaching half the
    /// spacing of the kit's closest channels is refused.
    #[arg(long, value_name = "TOL", value_parser = parse_precursor_tol, hide_short_help = true,
          help_heading = H_ISOBARIC)]
    pub(crate) tmt_tol: Option<Tolerance>,

    /// Isotope impurity matrix of the kit lot, in quantms text form [default: no correction]
    ///
    /// One line per channel with the -2/-1/+1/+2 percentages (eight columns for TMTpro
    /// 16/18). Corrected intensities are solved by non-negative least squares; the raw
    /// values are kept alongside. Not for tmt32 and tmt35.
    #[arg(long, value_name = "FILE", hide_short_help = true, help_heading = H_ISOBARIC)]
    pub(crate) tmt_correction: Option<PathBuf>,

    /// Quantification gate: least precursor isolation purity (0: off)
    ///
    /// The purity is the share of the MS1 isolation window that belongs to the precursor.
    /// PSMs that fail the gate are kept in tmt.parquet and the QPX parquet with quant_weight
    /// 0 (quant_gate_passed false in tmt.parquet).
    #[arg(long, value_name = "PURITY", default_value = "0.7", value_parser = parse_unit_fraction,
          hide_short_help = true, help_heading = H_ISOBARIC)]
    pub(crate) tmt_min_purity: f64,

    /// Quantification gate: least summed reporter signal-to-noise (0: off)
    ///
    /// Applies where the reporter scan records noise (Thermo .raw, or mzML converted with
    /// noise arrays); other PSMs are not held to it.
    #[arg(long, value_name = "SN", default_value = "0", value_parser = parse_non_negative,
          hide_short_help = true, help_heading = H_ISOBARIC)]
    pub(crate) tmt_min_sn: f64,

    /// Remove the co-isolation interference of high-resolution MS2 reporters
    ///
    /// Writes the reporter intensities with the interference removed as an extra unmix_
    /// channel block and an interference column in tmt.parquet (the other columns are
    /// unchanged). A per-run regression explains each scan's summed reporter signal by the
    /// precursor's, the co-isolated ions' and the sub-threshold ion current in the
    /// isolation window of the bordering MS1 scans; the interference is removed with the
    /// run's mean reporter profile. The auto strength is the one under which a protein's
    /// interference-heavy PSMs best agree with its clean PSMs (in the tmt.parquet
    /// metadata). auto and on apply to --tmt-level 2 only; on needs --tmt.
    #[arg(long, value_name = "MODE", value_enum, default_value = "auto", hide_short_help = true,
          help_heading = H_ISOBARIC)]
    pub(crate) tmt_unmix: TmtUnmix,

    /// With --tmt --mbr: RMS log2 deviation from its reference scan above which a transferred profile counts as wrong
    ///
    /// The definition of an error used to learn the transfer risk, not an acceptance
    /// threshold.
    #[arg(long = "mbr-wrong-log2", value_name = "LOG2", default_value = "0.5", value_parser = parse_positive_tol,
          hide_short_help = true, help_heading = H_ISOBARIC)]
    pub(crate) mbr_wrong_log2: f64,

    /// MS1 m/z tolerance of the isotope chromatograms, e.g. 10ppm [default: 10ppm]
    #[arg(long, value_name = "TOL", value_parser = parse_precursor_tol, hide_short_help = true,
          help_heading = H_LFQ)]
    pub(crate) lfq_tol: Option<Tolerance>,

    /// Retention-time half-window around the identifying PSM, in seconds [default: 60; 120 with --glyco]
    #[arg(long, value_name = "SECONDS", value_parser = parse_positive_tol, hide_short_help = true,
          help_heading = H_LFQ)]
    pub(crate) lfq_rt_window: Option<f64>,

    /// Experimental: Gaussian peak width for peak detection, in seconds (eight sigma) [default: Savitzky-Golay smoothing]
    ///
    /// Areas still use the raw signal. Without it, peaks are detected on five-scan
    /// Savitzky-Golay smoothed traces. Check the feature q-value coverage on the
    /// acquisition before relying on it.
    #[arg(long, value_name = "SECONDS", value_parser = parse_positive_tol, hide_short_help = true,
          help_heading = H_LFQ)]
    pub(crate) lfq_gaussian_width: Option<f64>,

    /// Least isotope-envelope cosine of a feature marked accepted in lfq.parquet
    ///
    /// lfq.parquet keeps every feature with its cosine.
    #[arg(long, value_name = "COSINE", default_value = "0.7", value_parser = parse_unit_fraction,
          hide_short_help = true, help_heading = H_LFQ)]
    pub(crate) lfq_min_cosine: f64,

    /// Reference channel of the --label ratios, 1-based
    #[arg(long, value_name = "N", default_value_t = 1, hide_short_help = true,
          help_heading = H_LABEL)]
    pub(crate) label_reference: usize,

    /// SILAC: fit the heavy proline made from heavy arginine and add it back to its channel
    ///
    /// Arg6 gives 13C5 Pro and Arg10 13C5 15N1 Pro, fitted as a satellite of each
    /// proline-containing peptide's arginine-labelled channels. label.parquet reports the
    /// converted fraction.
    #[arg(long, hide_short_help = true, help_heading = H_LABEL)]
    pub(crate) label_arg_pro: bool,

    /// 2H enrichment of metabolic deuterium labels (Lys4, Leu3)
    #[arg(long = "label-2h-enrichment-silac", value_name = "FRACTION", default_value = "0.96", value_parser = parse_unit_fraction,
          hide_short_help = true, help_heading = H_LABEL)]
    pub(crate) label_2h_enrichment_silac: f64,

    /// 2H enrichment of chemical deuterium labels (dimethyl)
    #[arg(long = "label-2h-enrichment-chemical", value_name = "FRACTION", default_value = "0.99", value_parser = parse_unit_fraction,
          hide_short_help = true, help_heading = H_LABEL)]
    pub(crate) label_2h_enrichment_chemical: f64,

    /// 13C and 15N enrichment of heavy labels, per atom
    #[arg(long, value_name = "FRACTION", default_value = "0.99", value_parser = parse_unit_fraction,
          hide_short_help = true, help_heading = H_LABEL)]
    pub(crate) label_heavy_enrichment: f64,

    /// External .gdb glycan database of --glyco (instead of --glyco-species)
    ///
    /// Each canonical string keeps its tree structure (core- or antenna-fucose), and the
    /// file's NeuGc content is used as it is. Not with --glyco-glycan-wurcs or
    /// --glyco-species.
    #[arg(long, value_name = "FILE", hide_short_help = true, help_heading = H_GLYCO)]
    pub(crate) glyco_glycan_gdb: Option<PathBuf>,

    /// WURCS 2.0 glycan database CSV of --glyco (instead of --glyco-species)
    ///
    /// Columns GlycanID, Mass, GlycanComposition and StructureInformation. Each
    /// StructureInformation is reduced to the five-symbol tree of a .gdb line (anomericity,
    /// ring and substituents dropped); mass and composition are computed from the tree. Not
    /// with --glyco-glycan-gdb or --glyco-species.
    #[arg(long, value_name = "FILE", hide_short_help = true, help_heading = H_GLYCO)]
    pub(crate) glyco_glycan_wurcs: Option<PathBuf>,

    /// Glycan core of --glyco
    ///
    /// An o run attaches glycans at any Ser or Thr (no sequon) and needs an O-glycan .gdb.
    #[arg(long, value_name = "CORE", value_enum, default_value_t = GlycoCoreFlag::N,
          hide_short_help = true, help_heading = H_GLYCO)]
    pub(crate) glyco_core: GlycoCoreFlag,

    /// Correct each precursor to the monoisotopic peak of its MS1 envelope (--glyco)
    ///
    /// auto fits the MS1 envelope at the reported charge against a glycopeptide isotope
    /// model under the hypotheses "the recorded precursor is the M+k peak", k = 0..6, and
    /// moves the precursor down when a k > 0 clearly wins. A fitted spectrum is searched
    /// with isotope errors 0..1 instead of the glyco default 0..2; spectra without a linked
    /// MS1 or charge keep the default, and an explicit --isotope-error is used as given.
    /// Needs MS1 (mzML, .raw or .d); without it the search is as with off.
    #[arg(long, value_name = "MODE", value_enum, default_value_t = PrecursorMonoFlag::Auto,
          hide_short_help = true, help_heading = H_GLYCO)]
    pub(crate) precursor_mono: PrecursorMonoFlag,

    /// Write one TSV row per MS2 with the recorded and corrected precursor and
    /// the envelope fit under every hypothesis; a diagnostic.
    #[arg(long, hide = true)]
    pub(crate) precursor_mono_dump: Option<PathBuf>,

    /// Fragment tolerance in ppm of the glyco matching: oxonium ions, core Y ions, backbone, c/z
    ///
    /// Unset on low-resolution (ion-trap) MS2, it is 500 ppm (about 0.5 Da at m/z 1000).
    /// It is separate from --fragment-tol-ppm; the glycan-by-subtraction precursor matching
    /// uses --precursor-tol.
    #[arg(long, value_name = "PPM", default_value_t = 20.0f64, value_parser = parse_positive_tol,
          hide_short_help = true, help_heading = H_GLYCO)]
    pub(crate) glyco_tol_ppm: f64,

    /// Fewest trimannosyl-core Y ions of a reported --glyco PSM (0: off)
    ///
    /// Reads spectral evidence only, so targets and decoys are held to it alike. 2 suits
    /// collision-dominated (HCD) data; electron-transfer data lose identifications with it
    /// (DOCS.md §9). With --glyco-core o the core ladder is Y0 and Y1 only, so a value above
    /// 2 rejects every O-glycopeptide.
    #[arg(long, value_name = "N", default_value_t = 0u32, hide_short_help = true,
          help_heading = H_GLYCO)]
    pub(crate) glyco_min_core_y: u32,

    /// On ETD scans, take the candidate backbones from the paired HCD scan of the same precursor
    ///
    /// c/z ions are still scored on the ETD scan. Scans without an HCD partner, and
    /// multi-file runs, are searched unpaired.
    #[arg(long, value_name = "BOOL", default_value_t = true, action = clap::ArgAction::Set,
          hide_short_help = true, help_heading = H_GLYCO)]
    pub(crate) glyco_hcd_pair: bool,

    /// Most intense peaks considered when generating --glyco candidates (0: all)
    ///
    /// Scoring always reads the whole spectrum. 300 to 500 keeps profile or very dense
    /// scans from taking seconds each. Changing it changes results.
    #[arg(long, value_name = "N", default_value_t = 0usize, hide_short_help = true,
          help_heading = H_GLYCO)]
    pub(crate) glyco_max_peaks: usize,

    /// Largest c/z fragment charge of --glyco ETD scoring [default: from whether the spectrum is deconvoluted]
    ///
    /// Set it only for data known to carry unresolved high-charge c/z ions.
    #[arg(long, value_name = "Z", hide_short_help = true, help_heading = H_GLYCO)]
    pub(crate) glyco_cz_max_charge: Option<u8>,

    /// Largest glycan Y fragment charge of --glyco
    #[arg(long, value_name = "Z", default_value_t = 3u8, hide_short_help = true,
          help_heading = H_GLYCO)]
    pub(crate) glyco_y_max_charge: u8,

    /// Candidate retrieval tolerance in ppm of --glyco [default: --glyco-tol-ppm on high-res MS2, 0.5 Da on low-res]
    ///
    /// Retrieval only; the scoring tolerance is unchanged.
    #[arg(long, value_name = "PPM", value_parser = parse_positive_tol, hide_short_help = true,
          help_heading = H_GLYCO)]
    pub(crate) glyco_retrieval_tol_ppm: Option<f64>,

    /// Peptide-first b/y fragment-index retrieval instead of the glycan-mass branch (slower)
    #[arg(long, hide_short_help = true, help_heading = H_GLYCO)]
    pub(crate) glyco_peptide_first: bool,

    /// Index only N-X-S/T peptides for --glyco, a much smaller index
    ///
    /// Glyco scoring reads only sequon-bearing candidates, so the results are those of the
    /// full index except for per-spectrum candidate-distribution features (RawScore,
    /// CandidateRankEntropy) on a few PSMs.
    #[arg(long, hide_short_help = true, help_heading = H_GLYCO)]
    pub(crate) glyco_index_sequon_only: bool,

    /// Restrict --glyco scoring to the scan numbers in this file, one per line;
    /// a diagnostic.
    #[arg(long, hide = true)]
    pub(crate) glyco_scans: Option<PathBuf>,

    /// Write every candidate row per scan, de-novo mass-residual hits included;
    /// a diagnostic PIN that must not be given to an FDR tool.
    #[arg(long, hide = true)]
    pub(crate) debug_glyco: bool,

    /// Write the training table of the --glyco --lfq --mbr match scorer to
    /// FILE (one row per scored candidate window: label, the ten match-model
    /// features and the peptide as group id), for retraining
    /// resources/quant/glyco_mbr_model.txt; the run is unchanged.
    #[arg(long, value_name = "FILE", hide = true)]
    pub(crate) glyco_mbr_match_train: Option<PathBuf>,

    /// q-value threshold of proteoform claims [default: the larger of --psm-fdr and 0.05]
    ///
    /// With rescoring, the filtered TSV, the proteoform table and quantification accept
    /// proteoform claims at this level instead of --psm-fdr. A proteoform claim is a PSM
    /// whose every protein is a UniProt variant entry; its q-value includes the class FDR
    /// and, with knockoffs, the knockoff estimate. Isoform and chain PSMs keep --psm-fdr.
    #[arg(long, value_name = "FDR", value_parser = parse_unit_fraction, hide_short_help = true,
          help_heading = ADVANCED_FDR)]
    pub(crate) proteoform_fdr: Option<f64>,

    /// Precursor q-value threshold of a match-between-runs library donor [default: inherits --psm-fdr]
    ///
    /// The best-PSM target/decoy q-value of the precursor.
    #[arg(long, value_name = "FDR", value_parser = parse_unit_fraction, hide_short_help = true,
          help_heading = ADVANCED_FDR)]
    pub(crate) mbr_donor_fdr: Option<f64>,

    /// q-value threshold of a reported transfer [default: inherits --quant-fdr]
    ///
    /// Transfers at or below it are reported as transfer_ms1 in <base>.mbr.parquet; with
    /// --tmt, transferred reporter rows at or below it pass the quantification gate; with
    /// --label --mbr it is the label transfer threshold (the identity-decoy false-donor
    /// rate added); with --glyco --lfq --mbr the glycopeptide transfer threshold. When
    /// --quant-fdr is not given either, --label and --glyco use 0.05.
    #[arg(long, value_name = "FDR", value_parser = parse_unit_fraction, hide_short_help = true,
          help_heading = ADVANCED_FDR)]
    pub(crate) quant_transfer_fdr: Option<f64>,

    /// --label: q-value threshold of a channel matched without an identification [default: inherits --quant-fdr]
    ///
    /// From the target/decoy competition, within its run, of the channel's fit against the
    /// same fit at non-label mass shifts. When --quant-fdr is not given either, 0.05.
    #[arg(long, value_name = "FDR", value_parser = parse_unit_fraction, hide_short_help = true,
          help_heading = ADVANCED_FDR)]
    pub(crate) quant_channel_fdr: Option<f64>,
}

impl SearchArgs {
    pub(crate) fn psm_fdr(&self) -> f64 {
        self.psm_fdr.unwrap_or(DEFAULT_FDR)
    }

    /// True when --rescore runs an external Percolator (--percolator-bin or
    /// --percolator-docker) instead of the native engine.
    pub(crate) fn external_percolator(&self) -> bool {
        self.percolator_bin.is_some() || self.percolator_docker
    }

    pub(crate) fn proteoform_fdr(&self) -> f64 {
        proteoform_fdr(self.proteoform_fdr, self.psm_fdr())
    }

    pub(crate) fn mbr_donor_fdr(&self) -> f64 {
        self.mbr_donor_fdr.unwrap_or_else(|| self.psm_fdr())
    }

    pub(crate) fn quant_transfer_fdr(&self) -> f64 {
        self.quant_transfer_fdr
            .or(self.quant_fdr)
            .unwrap_or(DEFAULT_FDR)
    }

    pub(crate) fn coarse_transfer_fdr(&self) -> f64 {
        self.quant_transfer_fdr
            .or(self.quant_fdr)
            .unwrap_or(COARSE_QUANT_FDR)
    }

    pub(crate) fn quant_channel_fdr(&self) -> f64 {
        self.quant_channel_fdr
            .or(self.quant_fdr)
            .unwrap_or(COARSE_QUANT_FDR)
    }

    pub(crate) fn lfq_feature_fdr(&self) -> f64 {
        self.quant_fdr.unwrap_or(COARSE_QUANT_FDR)
    }
}

pub(super) fn proteoform_fdr(explicit: Option<f64>, psm_fdr: f64) -> f64 {
    explicit.unwrap_or(psm_fdr.max(PROTEOFORM_FDR_FLOOR))
}

/// Derived-geometry parameters for model training.
///
/// Every field has the default the training pipeline shipped with; the `train`
/// subcommands expose them as flags so a run can be reproduced from its command line
/// rather than from an environment the command line does not record.
#[derive(Debug, Clone, Copy, clap::Args)]
pub(crate) struct GeometryArgs {
    /// Number of score segments.
    #[arg(long = "geo-segments", default_value_t = 2i32)]
    pub(crate) segments: i32,
    /// Maximum peak rank considered.
    #[arg(long = "geo-max-rank", default_value_t = 150i32)]
    pub(crate) max_rank: i32,
    /// Target peptides per mass tier.
    #[arg(long = "geo-occupancy", default_value_t = 2500usize)]
    pub(crate) occupancy: usize,
    /// Maximum number of mass tiers.
    #[arg(long = "geo-max-tiers", default_value_t = 33usize)]
    pub(crate) max_tiers: usize,
    /// Maximum fragment charge modelled.
    #[arg(long = "geo-max-fragment-charge", default_value_t = 3i32)]
    pub(crate) max_fragment_charge: i32,
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

/// `--dia-ms1-mobility-margin`: a non-negative 1/K0 margin, or `auto` (stored as -1).
fn parse_mobility_margin(v: &str) -> Result<f64, String> {
    if v == "auto" {
        return Ok(-1.0);
    }
    match v.parse::<f64>() {
        Ok(m) if m.is_finite() && m >= 0.0 => Ok(m),
        _ => Err(format!("expected a non-negative 1/K0 margin or `auto`, got `{v}`")),
    }
}
