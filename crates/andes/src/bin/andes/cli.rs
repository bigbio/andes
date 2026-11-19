use std::path::{Path, PathBuf};

use clap::builder::{PossibleValue, TypedValueParser};
use clap::{Args, ValueEnum};
use model::Tolerance;
use model_train::geometry::GeometryConfig;
use search::decoy::DecoyStrategy;
use search::PrecursorCalMode;

pub(crate) const H_IO: &str = "Input and output";
pub(crate) const H_SEARCH: &str = "Search";
pub(crate) const H_FDR: &str = "FDR";
pub(crate) const H_MODES: &str = "Modes";
pub(crate) const H_QUANT: &str = "Quantification";
pub(crate) const H_DATABASE: &str = "Database variants";
pub(crate) const H_RESCORING: &str = "Rescoring";
pub(crate) const H_RUNTIME: &str = "Calibration and runtime";
pub(crate) const H_PREDICTION: &str = "Prediction";
pub(crate) const H_ISOBARIC: &str = "Isobaric (TMT, iTRAQ)";
pub(crate) const H_LFQ: &str = "Label-free";
pub(crate) const H_LABEL: &str = "Labelled pairs";
pub(crate) const H_GLYCO: &str = "Glycopeptides";
pub(crate) const ADVANCED_FDR: &str = "Advanced FDR";
pub(crate) const H_ADVANCED: &str = "Advanced";

/// Fragmentation method of MGF spectra.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub(crate) enum Fragmentation {
    #[value(name = "auto")]
    Auto,
    #[value(name = "CID")]
    Cid,
    #[value(name = "ETD")]
    Etd,
    #[value(name = "HCD")]
    Hcd,
    #[value(name = "UVPD")]
    Uvpd,
}

/// Search protocol: the labelling or enrichment of the sample, which selects
/// the protocol model.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub(crate) enum Protocol {
    #[value(name = "auto")]
    Auto,
    #[value(name = "phospho")]
    Phospho,
    #[value(name = "iTRAQ")]
    Itraq,
    #[value(name = "iTRAQ-phospho")]
    ItraqPhospho,
    #[value(name = "TMT")]
    Tmt,
    #[value(name = "standard")]
    Standard,
}

/// `--tmt-unmix`: whether and how strongly the co-isolation interference of
/// high-resolution MS2 reporter ions is removed.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum)]
pub(crate) enum TmtUnmix {
    /// Remove it at a strength learned per run (a no-op without --tmt, and
    /// for low-resolution MS2 or MS3 reporters)
    #[default]
    #[value(name = "auto")]
    Auto,
    /// Write no interference estimate
    #[value(name = "off")]
    Off,
    /// Remove all of the estimated interference
    #[value(name = "on")]
    On,
}

/// `--enzyme-specificity`: the peptide ends that must be cleavage sites.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub(crate) enum EnzymeSpecificity {
    /// Neither end
    #[value(name = "non-specific")]
    NonSpecific,
    /// At least one end
    #[value(name = "semi")]
    Semi,
    /// Both ends
    #[value(name = "fully")]
    Fully,
}

/// `--score`: what ranks the candidates and fills the PIN `RawScore`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum)]
pub(crate) enum ScoreFlag {
    /// strong for a high-resolution model, rank for a low-resolution one
    #[default]
    Auto,
    /// The rank score of the model
    Rank,
    /// The fused intensity and competition score
    Strong,
}

pub(crate) const DEFAULT_GBDT_MAX_TREES: usize = 100;

pub(crate) const DEFAULT_FDR: f64 = 0.01;

pub(crate) const PROTEOFORM_FDR_FLOOR: f64 = 0.05;

pub(crate) const COARSE_QUANT_FDR: f64 = 0.05;

/// `--fragment-index`: candidate retrieval.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum)]
pub(crate) enum FragmentIndexFlag {
    /// The fragment-ion index (--glyco looks candidates up by backbone mass)
    #[default]
    Auto,
    /// As auto; with --glyco, a warning and per-spectrum enumeration
    On,
    /// Score every candidate in each precursor window
    Off,
}

/// One `--predict` entry: a predicted-feature block added to the PIN.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, ValueEnum)]
pub(crate) enum PredictFeature {
    /// Retention-time deviation from the bundled andes-rt model
    #[value(name = "rt")]
    Rt,
    /// Spectral similarity to the fragment intensities of the bundled andes-ms2 model
    #[value(name = "ms2")]
    Ms2,
    /// Spectral similarity to the fragment intensities of AlphaPeptDeep
    #[value(name = "peptdeep-ms2")]
    PeptDeepMs2,
    /// Retention-time deviation from AlphaPeptDeep
    #[value(name = "peptdeep-rt")]
    PeptDeepRt,
    /// Retention-time deviation from DeepLC
    #[value(name = "deeplc-rt")]
    DeepLcRt,
    /// No predicted features
    #[value(name = "none")]
    None,
}

/// `--predictor-backend`: ONNX engine of the third-party predictors.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub(crate) enum PredictorBackend {
    /// tract (pure Rust, built in)
    Tract,
    /// ONNX Runtime (a build with --features ort, plus the runtime library)
    Ort,
}

#[derive(Clone)]
pub(crate) struct Listed<T> {
    values: &'static [(&'static str, &'static str)],
    parse: fn(&str) -> Result<T, String>,
}

impl<T> Listed<T> {
    pub(crate) const fn new(
        values: &'static [(&'static str, &'static str)],
        parse: fn(&str) -> Result<T, String>,
    ) -> Self {
        Self { values, parse }
    }
}

impl<T: Clone + Send + Sync + 'static> TypedValueParser for Listed<T> {
    type Value = T;

    fn parse_ref(
        &self,
        cmd: &clap::Command,
        arg: Option<&clap::Arg>,
        value: &std::ffi::OsStr,
    ) -> Result<T, clap::Error> {
        self.parse.parse_ref(cmd, arg, value)
    }

    fn possible_values(&self) -> Option<Box<dyn Iterator<Item = PossibleValue> + '_>> {
        Some(Box::new(self.values.iter().map(|&(name, help)| {
            let value = PossibleValue::new(name);
            if help.is_empty() {
                value
            } else {
                value.help(help)
            }
        })))
    }
}

pub(crate) const ENZYMES: &[(&str, &str)] = &[
    ("trypsin", ""),
    ("chymotrypsin", ""),
    ("lysc", ""),
    ("aspn", ""),
    ("gluc", ""),
    ("lysn", ""),
    ("argc", ""),
    ("alphalp", ""),
    ("nocleavage", ""),
    ("nonspecific", ""),
];

pub(crate) const DECOY_STRATEGIES: &[(&str, &str)] = &[
    ("shuffle", "seeded shuffle of each segment"),
    (
        "reverse",
        "pseudo-reverse each segment; on ETD spectra a pseudo-reversed K-terminal decoy \
         shares c/z ions with its target at every I/L",
    ),
    (
        "sequon",
        "shuffle, then restore the segment's N-X-S/T sequon count; use it with --glyco",
    ),
    (
        "none",
        "no decoys: the database already holds them (--decoy-prefix, --decoy-suffix), or \
         the FDR is computed elsewhere",
    ),
];

pub(crate) const PLEXES: &[(&str, &str)] = &[
    ("tmt6", ""),
    ("tmt10", ""),
    ("tmt11", ""),
    ("tmt16", ""),
    ("tmt18", ""),
    ("tmt32", ""),
    ("tmt35", ""),
    ("itraq4", ""),
    ("itraq8", ""),
];

pub(crate) const EVIDENCE_LEVELS: &[(&str, &str)] = &[
    (
        "experimental",
        "experimental and large-scale proteomics evidence",
    ),
    (
        "curated",
        "also manual assertions without experimental support (e.g. by similarity); \
         Swiss-Prot chains without evidence codes count here",
    ),
    (
        "any",
        "also automatic annotation and features without evidence",
    ),
];

pub(crate) const PRECURSOR_CAL_MODES: &[(&str, &str)] = &[
    ("off", ""),
    ("auto", "skip the correction when the sample is too small"),
    ("on", ""),
];

#[derive(Debug, Clone, Default)]
pub(crate) struct ModList(pub(crate) Vec<model::Modification>);

pub(crate) fn parse_fixed_mods(s: &str) -> Result<ModList, String> {
    model::mod_spec::parse_mod_list(s, true)
        .map(ModList)
        .map_err(|e| e.to_string())
}

pub(crate) fn parse_variable_mods(s: &str) -> Result<ModList, String> {
    model::mod_spec::parse_mod_list(s, false)
        .map(ModList)
        .map_err(|e| e.to_string())
}

#[derive(Debug, Clone)]
pub(crate) struct XlSet(pub(crate) Vec<xl::crosslinker::Crosslinker>);

impl XlSet {
    /// One or more crosslinkers with distinct names; `key` names the option in an error.
    pub(crate) fn new(
        key: &str,
        linkers: Vec<xl::crosslinker::Crosslinker>,
    ) -> Result<Self, String> {
        if linkers.is_empty() {
            return Err(format!("{key}: no crosslinker"));
        }
        for (i, a) in linkers.iter().enumerate() {
            if linkers[..i].iter().any(|b| b.name == a.name) {
                return Err(format!("{key}: {} is given twice", a.name));
            }
        }
        Ok(Self(linkers))
    }
}

pub(crate) fn parse_xl_one(s: &str) -> Result<xl::crosslinker::Crosslinker, String> {
    xl::spec::XlDecl::parse(s)
        .and_then(|d| d.resolve(model::unimod::bundled()))
        .map_err(|e| e.to_string())
}

pub(crate) fn parse_xl(s: &str) -> Result<XlSet, String> {
    let linkers = s
        .split(';')
        .map(str::trim)
        .filter(|x| !x.is_empty())
        .map(parse_xl_one)
        .collect::<Result<Vec<_>, _>>()?;
    XlSet::new("--xl", linkers)
}

/// The protein database options shared by the search, `andes inference` and
/// `andes quant finalize`: how the database is digested and expanded, and how
/// a decoy is recognised. Inference and finalize take the search's values.
#[derive(Args, Debug, Clone)]
pub(crate) struct DatabaseArgs {
    /// Protease; a list (gluc,trypsin or gluc+trypsin) cleaves at the sites of each one
    ///
    /// The first protease of a list is the primary one: it selects the model and the
    /// cleavage feature. A non-specific protease in a list makes the whole digest
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
    /// Target sequences only: andes generates the decoys (--decoy-strategy). A UniProt XML
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
    /// residue (a non-specific search permutes whole sequences). With none and a
    /// target-only database the search has no decoys and warns.
    #[arg(long, value_name = "STRATEGY", default_value = "shuffle", help_heading = H_SEARCH,
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
    /// anchors of --refine. andes computes q-values only with --rescore or
    /// --rescore-native; without either, setting it warns and writes no filtered output.
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

    /// Run Percolator on the PIN and write q-values, PEPs and the PSMs at --psm-fdr
    ///
    /// Percolator's PEP and q-value are joined into the outputs (QPX
    /// posterior_error_probability and a q-value score) and the filtered <stem>.q<fdr>.tsv
    /// is written. Percolator is --percolator-bin, else `percolator` on PATH, else the
    /// Docker image.
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

    /// Second pass for further modifications (oxidation, deamidation, ...) over the confident proteins
    ///
    /// The modifications are a built-in tier, or the tier of --refine-config.
    #[arg(long, help_heading = H_MODES)]
    pub(crate) refine: bool,

    /// Site localization of variable modifications (`ipf`); needs --rescore or --rescore-native
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
    /// (one peptide sequence on both sides) and inter-protein links and pooled over the runs,
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

    /// Rescore as --rescore, with the native Percolator engine instead of the binary
    ///
    /// Every Percolator pass of --rescore (the main run, the per-run q-values, the PEP
    /// repair and the scores of every PIN row) runs in process with the algorithm of
    /// upstream Percolator master (commit a9783ee8) and Percolator's output files; no
    /// Percolator binary or Docker is needed. Takes no --percolator-args, --percolator-bin,
    /// --percolator-docker or --percolator-image.
    #[arg(long, hide_short_help = true, help_heading = H_RESCORING)]
    pub(crate) rescore_native: bool,

    /// Storage budget of the native Percolator engine in GB, not a process-memory cap: PIN features, PSM texts, normalised training rows and score buffers spill to --temp-dir once they would take more than their share of it; per-PSM metadata (about 160 bytes a PSM) is charged to it but always stays in memory
    #[arg(long, value_name = "GB", default_value_t = 16,
          value_parser = clap::value_parser!(u64).range(1..), hide = true,
          help_heading = H_RESCORING)]
    pub(crate) rescore_memory_gb: u64,

    /// PEP threshold applied with --psm-fdr: a PSM must pass both [default: none]
    #[arg(long, value_name = "PEP", value_parser = parse_unit_fraction, hide_short_help = true,
          help_heading = H_RESCORING)]
    pub(crate) pep: Option<f64>,

    /// Percolator binary [default: percolator on PATH, else the Docker image]
    #[arg(long, value_name = "FILE", hide_short_help = true, help_heading = H_RESCORING)]
    pub(crate) percolator_bin: Option<PathBuf>,

    /// Run Percolator from its Docker image (needs the docker command)
    #[arg(long, hide_short_help = true, help_heading = H_RESCORING)]
    pub(crate) percolator_docker: bool,

    /// Docker image of Percolator
    #[arg(long, value_name = "IMAGE", default_value = output::DEFAULT_PERCOLATOR_IMAGE,
          hide_short_help = true, help_heading = H_RESCORING)]
    pub(crate) percolator_image: String,

    /// Extra arguments for Percolator, e.g. "--testFDR 0.05"
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

    /// Quantification gate: largest Percolator PEP (1: off; needs --rescore)
    #[arg(long, value_name = "PEP", default_value = "0.05", value_parser = parse_unit_fraction,
          hide_short_help = true, help_heading = H_ISOBARIC)]
    pub(crate) tmt_max_pep: f64,

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

#[derive(Args, Debug)]
pub(crate) struct PercolatorArgs {
    /// Storage budget of the engine in GB, not a process-memory cap: PIN features, PSM texts, normalised training rows and score buffers spill to --temp-dir once they would take more than their share of it; per-PSM metadata (about 160 bytes a PSM) is charged to it but always stays in memory
    #[arg(long, value_name = "GB", default_value_t = 16,
          value_parser = clap::value_parser!(u64).range(1..))]
    pub(crate) rescore_memory_gb: u64,

    /// Directory of the spill files, created if missing
    #[arg(long, value_name = "DIR", required = true)]
    pub(crate) temp_dir: PathBuf,

    /// Percolator arguments as given to the percolator binary, the PIN last
    #[arg(value_name = "PERCOLATOR_ARGS", required = true, num_args = 1..,
          trailing_var_arg = true, allow_hyphen_values = true)]
    pub(crate) args: Vec<String>,
}

/// `andes inference`: protein groups over the rescored PSMs of one or more runs.
#[derive(Args, Debug)]
pub(crate) struct InferenceArgs {
    /// Percolator target PSM results (--results-psms) of every run
    ///
    /// Their scores must be on one scale: one Percolator run over all the runs (one file,
    /// or one file per run from that run).
    #[arg(long = "psms", value_name = "TSV", num_args = 1.., required = true)]
    pub(crate) psms: Vec<PathBuf>,

    /// Percolator decoy PSM results (--decoy-results-psms) of the same Percolator run
    #[arg(long = "decoy-psms", value_name = "TSV", num_args = 1.., required = true)]
    pub(crate) decoy_psms: Vec<PathBuf>,

    /// The protein database of the search, with its options below
    ///
    /// With the search's --enzyme, --fixed-mods, --variable-mods, decoy affixes,
    /// --annotations, --uniprot-* and --knockoffs, protein families, knockoff entries, knockoff sites and decoy pairs are
    /// the search's.
    #[arg(long, value_name = "FILE")]
    pub(crate) database: PathBuf,

    #[command(flatten)]
    pub(crate) proteins: ProteinArgs,

    /// Run names of the groups' grouped_runs [default: one per --psms file, named after it]
    #[arg(long = "runs", value_name = "NAME", num_args = 1..)]
    pub(crate) runs: Vec<String>,

    /// Output prefix: writes <PREFIX>.pg.parquet (the QPX protein-group view)
    #[arg(short = 'o', long = "output", value_name = "PREFIX")]
    pub(crate) output: PathBuf,
}

/// How protein inference reads the database and which PSMs count as evidence:
/// the search's database options and the PSM thresholds.
#[derive(Args, Debug, Clone)]
pub(crate) struct ProteinArgs {
    #[command(flatten)]
    pub(crate) db: DatabaseArgs,

    /// The search's --config file, for its database and modification keys
    ///
    /// Its enzyme, modification lists (with their config-only settings), decoy affixes,
    /// UniProt options, annotations and knockoffs apply here, below the flags given on the
    /// command line. Its other keys are not read.
    #[arg(long, value_name = "FILE", help_heading = H_IO)]
    pub(crate) config: Option<PathBuf>,

    /// PSM q-value threshold of the PSMs and peptides counted per group
    #[arg(long = "psm-fdr", value_name = "FDR", default_value = "0.01",
          value_parser = parse_unit_fraction, help_heading = H_FDR)]
    pub(crate) psm_fdr: f64,

    /// PSM q-value threshold of proteoform claims (variant-class PSMs) [default: the larger of --psm-fdr and 0.05]
    #[arg(long = "proteoform-fdr", value_name = "FDR", value_parser = parse_unit_fraction,
          hide_short_help = true, help_heading = ADVANCED_FDR)]
    pub(crate) proteoform_fdr: Option<f64>,

    /// Protein q-value threshold of the reported groups (1 keeps every group)
    ///
    /// The pg view keeps the target and decoy groups at or below it (the decoys show the
    /// FDR at the threshold), and `quant finalize`'s feature view the features of those
    /// target groups.
    #[arg(long = "protein-fdr", value_name = "FDR", default_value = "0.01",
          value_parser = parse_unit_fraction, help_heading = H_FDR)]
    pub(crate) protein_fdr: f64,
}

/// `andes quant`: the match-between-runs stages over a quant directory
/// written by a search with `--quant-dir`.
#[derive(clap::Subcommand, Debug)]
pub(crate) enum CohortCommand {
    /// Check a cohort and write its plan: every run's search settings and the FDR groups
    ///
    /// Searches nothing. Prints the groups and every problem found, and writes the plan
    /// only when there is no error.
    Plan(CohortPlanArgs),
    /// Crosslink FDR over the CSM tables of several runs, as one search would compute it
    ///
    /// Reads each run's <base>.xl.tsv, groups all CSMs into residue pairs and protein pairs,
    /// and gives every level its q-values once over all runs. Writes <output>.xl.tsv,
    /// <output>.xl.residue_pairs.tsv and <output>.xl.ppi.tsv.
    #[command(name = "xl-fdr", hide = true)]
    XlFdr(CohortXlFdrArgs),
    /// Plan the cohort, search every run, then run each group's stages
    ///
    /// Each run's outputs go to <output_dir>/runs/<run>/, each group's to
    /// <output_dir>/groups/<group>/. A run whose search has finished is kept, so a rerun
    /// searches only new or failed runs.
    Run(CohortRunArgs),
}

#[derive(clap::Args, Debug)]
pub(crate) struct CohortRunArgs {
    /// Config file with a `cohort` section
    #[arg(value_name = "CONFIG")]
    pub(crate) config: PathBuf,
    /// Worker threads of each search [default: the number of logical CPUs]
    #[arg(long, value_name = "N")]
    pub(crate) threads: Option<usize>,
}

#[derive(clap::Args, Debug)]
pub(crate) struct CohortXlFdrArgs {
    /// The runs' CSM tables (<base>.xl.tsv from `andes --xl`)
    #[arg(value_name = "CSM_TABLE", required = true)]
    pub(crate) tables: Vec<PathBuf>,
    /// Output base name: <output>.xl.tsv, <output>.xl.residue_pairs.tsv, <output>.xl.ppi.tsv
    #[arg(long, short, value_name = "BASE")]
    pub(crate) output: PathBuf,
}

#[derive(clap::Args, Debug)]
pub(crate) struct CohortPlanArgs {
    /// Config file with a `cohort` section; its other sections are the settings every
    /// run shares
    #[arg(value_name = "CONFIG")]
    pub(crate) config: PathBuf,
    /// Where to write the plan [default: <output_dir>/plan.json]
    #[arg(long, short, value_name = "FILE")]
    pub(crate) output: Option<PathBuf>,
}

#[derive(clap::Subcommand, Debug)]
pub(crate) enum QuantCommand {
    /// Build the precursor library and the joint retention-time alignment
    ///
    /// From every run's PSM table and the Percolator results of those PSMs. Writes
    /// library.parquet and alignment.parquet into the quant directory.
    Library(QuantLibraryArgs),
    /// Measure every library precursor in one run from its MS1 table
    ///
    /// Writes <stem>.extract.parquet.
    Extract(QuantExtractArgs),
    /// Transfer q-values and the match-between-runs tables from every run's extract table
    ///
    /// Writes <base>.mbr.parquet and <base>.mbr_validation.parquet.
    Finalize(Box<QuantFinalizeArgs>),
}

impl QuantCommand {
    pub(crate) fn threads(&self) -> Option<usize> {
        match self {
            QuantCommand::Library(a) => a.threads,
            QuantCommand::Extract(a) => a.threads,
            QuantCommand::Finalize(a) => a.threads,
        }
    }
}

#[derive(Args, Debug)]
pub(crate) struct QuantLibraryArgs {
    /// Directory with each run's <stem>.quant.parquet and <stem>.ms1.parquet
    #[arg(long = "quant-dir", value_name = "DIR")]
    pub(crate) quant_dir: PathBuf,

    /// Run stems to include [default: every <stem>.quant.parquet in the directory]
    ///
    /// A run whose PSM or MS1 table is missing is an error.
    #[arg(long = "runs", value_name = "STEM", num_args = 1..)]
    pub(crate) runs: Vec<String>,

    /// Percolator target PSM results (--results-psms) with their score column
    ///
    /// One file over every run (SpecIds distinct across runs), or <stem>=<file> once per
    /// run. PSMs are joined to the run tables on (run, SpecId).
    #[arg(long = "psms", value_name = "TSV|STEM=TSV", num_args = 1.., required = true)]
    pub(crate) psms: Vec<String>,

    /// Percolator decoy PSM results (--decoy-results-psms), in the form of --psms
    #[arg(long = "decoy-psms", value_name = "TSV|STEM=TSV", num_args = 1.., required = true)]
    pub(crate) decoy_psms: Vec<String>,

    /// PSM q-value threshold of a donor: the search's --psm-fdr gives the result of --mbr
    ///
    /// A target PSM at or below it identifies its precursor in its run: those runs are the
    /// precursor's donors and get identified rows; the other runs get transfers.
    #[arg(long = "psm-fdr", value_name = "FDR", default_value = "0.01",
          value_parser = parse_unit_fraction, help_heading = H_FDR)]
    pub(crate) psm_fdr: f64,

    /// Precursor q-value threshold of a library entry [default: inherits --psm-fdr]
    ///
    /// The best-PSM target/decoy q-value of the precursor. The search's --mbr-donor-fdr
    /// gives the result of --mbr.
    #[arg(long = "mbr-donor-fdr", value_name = "FDR", value_parser = parse_unit_fraction,
          hide_short_help = true, help_heading = ADVANCED_FDR)]
    pub(crate) mbr_donor_fdr: Option<f64>,

    /// PEP threshold of a donor PSM, applied with --psm-fdr [default: none]
    ///
    /// The Percolator posterior_error_prob column. The search's --pep gives the result of
    /// --mbr.
    #[arg(long = "pep", value_name = "PEP", value_parser = parse_unit_fraction,
          help_heading = H_FDR)]
    pub(crate) pep: Option<f64>,

    /// Worker threads [default: the number of logical CPUs]
    #[arg(long, value_name = "N")]
    pub(crate) threads: Option<usize>,
}

#[derive(Args, Debug)]
pub(crate) struct QuantExtractArgs {
    /// Directory with the library and the run's <stem>.ms1.parquet
    #[arg(long = "quant-dir", value_name = "DIR")]
    pub(crate) quant_dir: PathBuf,

    /// Stem of the run to measure
    #[arg(long = "run", value_name = "STEM")]
    pub(crate) run: String,

    /// MS1 m/z tolerance of the template fit in ppm, e.g. 10ppm [default: 10ppm]
    ///
    /// The search's --lfq-tol gives the result of --mbr.
    #[arg(long = "lfq-tol", value_name = "TOL", value_parser = parse_precursor_tol)]
    pub(crate) lfq_tol: Option<Tolerance>,

    /// Measure only the transfers that finalize can accept
    ///
    /// true measures a transfer only when the entry's best donor PSM passes the transfer
    /// donor threshold and it has at least --mbr-min-donor-fraction of the other aligned
    /// runs as donors; a skipped pair gets neither its target nor its decoy rows. false
    /// measures every transfer pair. The search's --mbr measures with true.
    #[arg(long = "mbr-preselect", value_name = "BOOL", default_value_t = true,
          action = clap::ArgAction::Set, hide_short_help = true, help_heading = H_ADVANCED)]
    pub(crate) mbr_preselect: bool,

    /// Fewest other aligned donor runs of a measured transfer, as a fraction of the aligned runs (0: no minimum)
    ///
    /// With --mbr-preselect true: the fraction of the library's aligned runs other than the
    /// acceptor, rounded up. A minimum leaves out some transfer fits (and their CPU time),
    /// including a few transfers that would have been accepted.
    #[arg(long = "mbr-min-donor-fraction", value_name = "FRACTION",
          default_value_t = crate::quant_mbr::DEFAULT_MIN_DONOR_FRACTION,
          value_parser = parse_unit_fraction, hide_short_help = true, help_heading = H_ADVANCED)]
    pub(crate) mbr_min_donor_fraction: f64,

    /// Worker threads [default: the number of logical CPUs]
    #[arg(long, value_name = "N")]
    pub(crate) threads: Option<usize>,
}

#[derive(Args, Debug)]
pub(crate) struct QuantFinalizeArgs {
    /// Directory with the library and every run's <stem>.extract.parquet
    #[arg(long = "quant-dir", value_name = "DIR")]
    pub(crate) quant_dir: PathBuf,

    /// Output base: writes <base>.mbr.parquet and <base>.mbr_validation.parquet
    #[arg(short = 'o', long = "output", value_name = "BASE")]
    pub(crate) output: PathBuf,

    /// Quantification FDR: q-value threshold of the transfers
    ///
    /// A quantification-quality q-value, not an identification FDR. The search's
    /// --quant-fdr gives the result of --mbr.
    #[arg(long = "quant-fdr", value_name = "FDR", default_value = "0.01",
          value_parser = parse_unit_fraction, help_heading = H_FDR)]
    pub(crate) quant_fdr: f64,

    /// q-value threshold of a transfer reported as transfer_ms1 [default: inherits --quant-fdr]
    #[arg(long = "quant-transfer-fdr", value_name = "FDR", value_parser = parse_unit_fraction,
          hide_short_help = true, help_heading = ADVANCED_FDR)]
    pub(crate) quant_transfer_fdr: Option<f64>,

    /// Also write the QPX views into this directory (needs --database, --psms and --decoy-psms)
    ///
    /// The protein groups of every run's PSMs as quantms.pg.parquet, and the features of
    /// the groups at --protein-fdr as quantms.feature.parquet, each with its group.
    #[arg(long = "output-parquet", value_name = "DIR",
          requires_all = ["database", "inference_psms", "inference_decoy_psms"])]
    pub(crate) output_parquet: Option<PathBuf>,

    /// The protein database of the search, with its options (as for andes inference)
    #[arg(long, value_name = "FILE")]
    pub(crate) database: Option<PathBuf>,

    /// Percolator target PSM results of every run (as for andes inference)
    ///
    /// One file from one Percolator run over all the runs, or one per run.
    #[arg(id = "inference_psms", long = "psms", value_name = "TSV", num_args = 1..)]
    pub(crate) psms: Vec<PathBuf>,

    /// Percolator decoy PSM results of the same Percolator run
    #[arg(id = "inference_decoy_psms", long = "decoy-psms", value_name = "TSV", num_args = 1..)]
    pub(crate) decoy_psms: Vec<PathBuf>,

    #[command(flatten)]
    pub(crate) proteins: ProteinArgs,

    /// Training rows per transfer or identified-quality classifier
    ///
    /// Above it the classifier is trained on the rows of a hashed sample of (precursor,
    /// run) pairs (every row is still scored); below it every row trains.
    #[arg(long = "max-train-rows", value_name = "N",
          default_value_t = crate::quant_finalize::DEFAULT_MAX_TRAIN_ROWS,
          hide_short_help = true, help_heading = H_ADVANCED)]
    pub(crate) max_train_rows: usize,

    /// Members of the bagged transfer classifier (1 to 64)
    ///
    /// Each fold's transfer score is the mean of this many models, each trained on its own
    /// sample of the training rows with its own seed, which keeps the accepted transfers
    /// stable across seeds at hundreds of runs; training and scoring time grow with it.
    #[arg(long = "transfer-models", value_name = "N",
          default_value_t = crate::quant_finalize::DEFAULT_TRANSFER_MODELS,
          value_parser = parse_transfer_models, hide_short_help = true, help_heading = H_ADVANCED)]
    pub(crate) transfer_models: usize,

    /// Reference cells of the run normalisation (about 150 bytes each)
    ///
    /// Above it the fit uses a hashed sample of whole precursors, which coarsens the
    /// RT-local offsets, so keep it above the cells of the study when memory allows.
    #[arg(long = "max-norm-cells", value_name = "N",
          default_value_t = crate::quant_finalize::DEFAULT_MAX_NORM_CELLS,
          hide_short_help = true, help_heading = H_ADVANCED)]
    pub(crate) max_norm_cells: usize,

    /// Donor fits summarised in one pass over the extract tables (about 90 bytes each)
    ///
    /// With more, the precursors are summarised in groups, one pass each. The results do
    /// not depend on it.
    #[arg(long = "max-donor-fits", value_name = "N",
          default_value_t = crate::quant_finalize::DEFAULT_MAX_DONOR_FITS,
          hide_short_help = true, help_heading = H_ADVANCED)]
    pub(crate) max_donor_fits: usize,

    /// Held-out classifier scores kept in memory (4 bytes each)
    ///
    /// Above it they are sorted and written to temporary tables. The results do not depend
    /// on it.
    #[arg(long = "max-held-scores", value_name = "N",
          default_value_t = crate::quant_finalize::DEFAULT_MAX_HELD_SCORES,
          hide_short_help = true, help_heading = H_ADVANCED)]
    pub(crate) max_held_scores: usize,

    /// Directory of the temporary tables, e.g. a local disk [default: the directory of the output base]
    ///
    /// The tables are removed at the end.
    #[arg(long = "temp-dir", value_name = "DIR")]
    pub(crate) temp_dir: Option<PathBuf>,

    /// Seed of the transfer FDR folds and models and of the samples
    #[arg(long = "seed", value_name = "N", default_value_t = crate::quant_cmd::MBR_SEED,
          hide_short_help = true, help_heading = H_ADVANCED)]
    pub(crate) seed: u64,

    /// Share of the precursors with an identified fit that report the trace area instead of the template area
    ///
    /// The trace area (trace_quantity of the extract) is the baseline-subtracted isotope
    /// traces integrated around the apex; the template area the Gaussian amplitude times
    /// the template width. The share, in [0, 1], goes to the brightest precursors by the
    /// median amplitude signal-to-noise of their identified fits, in every run: 0 reports
    /// the template area everywhere, 1 switches every such precursor. quant_q and the
    /// transfer q-values do not depend on it.
    #[arg(long = "quant-trace-fraction", value_name = "FRACTION",
          default_value_t = crate::quant_finalize::DEFAULT_TRACE_FRACTION,
          value_parser = parse_unit_fraction, hide_short_help = true, help_heading = H_ADVANCED)]
    pub(crate) quant_trace_fraction: f64,

    /// Worker threads [default: the number of logical CPUs]
    #[arg(long, value_name = "N")]
    pub(crate) threads: Option<usize>,
}

pub(crate) type Cli = SearchArgs;

impl SearchArgs {
    pub(crate) fn psm_fdr(&self) -> f64 {
        self.psm_fdr.unwrap_or(DEFAULT_FDR)
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

impl ProteinArgs {
    pub(crate) fn proteoform_fdr(&self) -> f64 {
        proteoform_fdr(self.proteoform_fdr, self.psm_fdr)
    }
}

impl QuantLibraryArgs {
    pub(crate) fn mbr_donor_fdr(&self) -> f64 {
        self.mbr_donor_fdr.unwrap_or(self.psm_fdr)
    }
}

impl QuantFinalizeArgs {
    pub(crate) fn quant_transfer_fdr(&self) -> f64 {
        self.quant_transfer_fdr.unwrap_or(self.quant_fdr)
    }
}

fn proteoform_fdr(explicit: Option<f64>, psm_fdr: f64) -> f64 {
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

#[derive(Copy, Clone, Debug, PartialEq, Eq, clap::ValueEnum)]
pub(crate) enum EthcdActivationFlag {
    /// Label EThcD/ETciD as EThcD (default; routes to the EThcD model).
    Ethcd,
    /// Label them ETD so the c/z scoring path engages.
    Etd,
}

/// Bundled species-specific N-glycan databases for `--glyco-species`.
#[derive(Copy, Clone, Debug, PartialEq, Eq, clap::ValueEnum)]
pub(crate) enum GlycoSpeciesFlag {
    /// Human N-glycans.
    Human,
    /// Mouse N-glycans.
    Mouse,
    /// Mouse, extended.
    MouseLarge,
    /// High-mannose N-glycans (species-independent).
    HighMannose,
}

impl GlycoSpeciesFlag {
    pub(crate) fn key(self) -> &'static str {
        match self {
            GlycoSpeciesFlag::Human => "human",
            GlycoSpeciesFlag::Mouse => "mouse",
            GlycoSpeciesFlag::MouseLarge => "mouse-large",
            GlycoSpeciesFlag::HighMannose => "high-mannose",
        }
    }
}

/// Glycan core type for `--glyco-core`: N- vs O-glycan search mode.
#[derive(Copy, Clone, Debug, PartialEq, Eq, clap::ValueEnum)]
pub(crate) enum GlycoCoreFlag {
    /// N-glycan: N-X-S/T sequon gate, trimannosyl core.
    N,
    /// O-glycan: Ser/Thr gate, single reducing-end GalNAc core.
    O,
}

impl From<GlycoCoreFlag> for glyco::glycan_first::GlycanCore {
    fn from(value: GlycoCoreFlag) -> Self {
        match value {
            GlycoCoreFlag::N => glyco::glycan_first::GlycanCore::NGlycan,
            GlycoCoreFlag::O => glyco::glycan_first::GlycanCore::OGlycan,
        }
    }
}

/// `--precursor-mono`: MS1 isotope-envelope precursor correction.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub(crate) enum PrecursorMonoFlag {
    /// Leave every precursor m/z as recorded.
    Off,
    /// Move a precursor down to the monoisotope its MS1 envelope supports (default).
    Auto,
}

pub(crate) fn parse_enum<E: ValueEnum>(what: &str, s: &str, ignore_case: bool) -> Result<E, String> {
    <E as ValueEnum>::from_str(s, ignore_case).map_err(|_| {
        let names: Vec<String> = E::value_variants()
            .iter()
            .filter_map(|v| v.to_possible_value())
            .map(|p| p.get_name().to_string())
            .collect();
        format!("invalid {what} `{s}`: expected {}", names.join("|"))
    })
}

pub(crate) fn parse_decoy_strategy(s: &str) -> Result<DecoyStrategy, String> {
    DecoyStrategy::from_name(s).ok_or_else(|| {
        format!("unknown --decoy-strategy '{s}' (expected shuffle/reverse/sequon/none)")
    })
}

pub(crate) fn decoy_strategy_name(strategy: DecoyStrategy) -> &'static str {
    match strategy {
        DecoyStrategy::Shuffle => "shuffle",
        DecoyStrategy::Reverse => "reverse",
        DecoyStrategy::Sequon => "sequon",
        DecoyStrategy::None => "none",
    }
}

pub(crate) fn parse_enzyme_list(s: &str) -> Result<String, String> {
    crate::model_select::parse_enzymes(s)
        .map(|_| s.to_string())
        .map_err(|e| e.to_string())
}

pub(crate) fn parse_plex(s: &str) -> Result<quant::Plex, String> {
    quant::Plex::parse(s)
}

pub(crate) fn dash_sep_index(s: &str) -> Option<usize> {
    let b = s.as_bytes();
    (1..b.len()).find(|&i| b[i] == b'-' && b[i - 1].is_ascii_digit())
}

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

pub(crate) const MAX_SUPPORTED_CHARGE: u8 = 50;

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

pub(crate) fn parse_isotope_error_range(s: &str) -> Result<(i8, i8), String> {
    parse_int_range::<i8>(s, "isotope-error")
}

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

pub(crate) fn parse_open_search_delta(s: &str) -> Result<(f64, f64), String> {
    let t = s.trim();
    let (min_s, max_s) = t.split_once("..").ok_or_else(|| {
        format!("invalid --open-search-delta `{s}`: expected MIN..MAX in Da, e.g. -100..300")
    })?;
    let parse = |v: &str| -> Result<f64, String> {
        v.trim()
            .parse::<f64>()
            .map_err(|_| format!("invalid --open-search-delta `{s}`: bad number `{}`", v.trim()))
    };
    let (lo, hi) = (parse(min_s)?, parse(max_s)?);
    if !lo.is_finite() || !hi.is_finite() || lo >= hi {
        return Err(format!(
            "invalid --open-search-delta `{s}`: need finite MIN < MAX (got {lo}..{hi})"
        ));
    }
    if lo > 0.0 || hi < 0.0 {
        return Err(format!(
            "invalid --open-search-delta `{s}`: the window must contain 0 (got {lo}..{hi})"
        ));
    }
    if hi - lo > 2000.0 {
        return Err(format!(
            "invalid --open-search-delta `{s}`: width {} Da exceeds the 2000 Da limit",
            hi - lo
        ));
    }
    Ok((lo, hi))
}

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

pub(crate) fn parse_non_negative(s: &str) -> Result<f64, String> {
    let v: f64 = s
        .trim()
        .parse()
        .map_err(|_| format!("invalid value `{s}`: expected a number >= 0"))?;
    if !v.is_finite() || v < 0.0 {
        return Err(format!(
            "invalid value `{s}`: must be finite and >= 0 (got {v})"
        ));
    }
    Ok(v)
}

pub(crate) fn parse_non_negative_f32(s: &str) -> Result<f32, String> {
    let v = parse_non_negative(s)?;
    if v > f32::MAX as f64 {
        return Err(format!("invalid value `{s}`: out of range"));
    }
    Ok(v as f32)
}

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

pub(crate) fn dir_or_output_dir(flag: Option<&Path>, output: &Path) -> PathBuf {
    match flag {
        Some(dir) => dir.to_path_buf(),
        None => match output.parent() {
            Some(p) if !p.as_os_str().is_empty() => p.to_path_buf(),
            _ => PathBuf::from("."),
        },
    }
}

pub(crate) fn ensure_scratch_dir(dir: &Path) -> Result<(), String> {
    std::fs::create_dir_all(dir).map_err(|e| scratch_dir_error(dir, &e))
}

pub(crate) fn scratch_dir_error(dir: &Path, e: &dyn std::fmt::Display) -> String {
    format!(
        "cannot create temporary files in {} ({e}); pass --temp-dir with a writable directory",
        dir.display()
    )
}

pub(crate) fn parse_precursor_cal(s: &str) -> Result<PrecursorCalMode, String> {
    match s.to_ascii_lowercase().as_str() {
        "auto" => Ok(PrecursorCalMode::Auto),
        "on" => Ok(PrecursorCalMode::On),
        "off" => Ok(PrecursorCalMode::Off),
        _ => Err(format!(
            "invalid precursor-cal `{s}`: expected off|auto|on"
        )),
    }
}

pub(crate) fn parse_uniprot_evidence(s: &str) -> Result<input::EvidenceLevel, String> {
    input::EvidenceLevel::from_name(s).ok_or_else(|| {
        format!("invalid UniProt evidence level `{s}`: expected experimental|curated|any")
    })
}

fn parse_transfer_models(s: &str) -> Result<usize, String> {
    match s.parse::<usize>() {
        Ok(n) if (1..=64).contains(&n) => Ok(n),
        _ => Err(format!("{s}: expected 1 to 64 transfer classifier members")),
    }
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
    fn open_search_delta_parses_a_range() {
        assert_eq!(parse_open_search_delta("-100..300"), Ok((-100.0, 300.0)));
        assert_eq!(parse_open_search_delta(" -20.5..79.966 "), Ok((-20.5, 79.966)));
        assert!(parse_open_search_delta("-100-300").is_err());
    }

    #[test]
    fn open_search_delta_rejects_bad_windows() {
        assert!(parse_open_search_delta("300..-100").is_err(), "min >= max");
        assert!(parse_open_search_delta("10..300").is_err(), "0 not in window");
        assert!(parse_open_search_delta("-300..-1").is_err(), "0 not in window");
        assert!(parse_open_search_delta("0..0").is_err(), "empty window");
        assert!(parse_open_search_delta("nan..300").is_err(), "non-finite");
        assert!(parse_open_search_delta("-1200..1200").is_err(), "width cap");
        assert!(parse_open_search_delta("abc").is_err(), "no separator");
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
    fn non_negative_f32_rejects_negative_and_nonfinite() {
        assert_eq!(parse_non_negative_f32("0.3").unwrap(), 0.3);
        assert_eq!(parse_non_negative_f32("0").unwrap(), 0.0);
        assert!(parse_non_negative_f32("-0.1").is_err());
        assert!(parse_non_negative_f32("nan").is_err());
        assert!(parse_non_negative_f32("inf").is_err());
        assert!(parse_non_negative_f32("1e300").is_err());
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
mod fdr_parameter_tests {
    use super::*;
    use clap::{CommandFactory, Parser};

    #[derive(Parser, Debug)]
    struct SearchWrap {
        #[command(flatten)]
        args: SearchArgs,
    }

    #[derive(Parser, Debug)]
    struct LibraryWrap {
        #[command(flatten)]
        args: QuantLibraryArgs,
    }

    #[derive(Parser, Debug)]
    struct FinalizeWrap {
        #[command(flatten)]
        args: QuantFinalizeArgs,
    }

    fn search(extra: &[&str]) -> SearchArgs {
        let argv = std::iter::once("andes").chain(extra.iter().copied());
        SearchWrap::try_parse_from(argv).unwrap().args
    }

    #[test]
    fn defaults_keep_the_previous_thresholds() {
        let s = search(&[]);
        assert_eq!(s.psm_fdr(), 0.01);
        assert_eq!(s.proteoform_fdr(), 0.05);
        assert_eq!(s.mbr_donor_fdr(), 0.01);
        assert_eq!(s.quant_transfer_fdr(), 0.01);
        assert_eq!(s.coarse_transfer_fdr(), 0.05);
        assert_eq!(s.quant_channel_fdr(), 0.05);
        assert_eq!(s.lfq_feature_fdr(), 0.05);
        assert_eq!(s.protein_fdr, None);
    }

    #[test]
    fn advanced_flags_inherit_their_main_flag() {
        let s = search(&["--psm-fdr", "0.02", "--quant-fdr", "0.03"]);
        assert_eq!(s.psm_fdr(), 0.02);
        assert_eq!(s.proteoform_fdr(), 0.05, "at least 0.05");
        assert_eq!(s.mbr_donor_fdr(), 0.02);
        assert_eq!(s.quant_transfer_fdr(), 0.03);
        assert_eq!(s.coarse_transfer_fdr(), 0.03);
        assert_eq!(s.quant_channel_fdr(), 0.03);
        assert_eq!(s.lfq_feature_fdr(), 0.03);
        assert_eq!(search(&["--psm-fdr", "0.1"]).proteoform_fdr(), 0.1);
    }

    #[test]
    fn advanced_flags_set_win_over_their_main_flag() {
        let s = search(&[
            "--psm-fdr",
            "0.02",
            "--quant-fdr",
            "0.03",
            "--proteoform-fdr",
            "0.04",
            "--mbr-donor-fdr",
            "0.005",
            "--quant-transfer-fdr",
            "0.06",
            "--quant-channel-fdr",
            "0.07",
        ]);
        assert_eq!(s.proteoform_fdr(), 0.04);
        assert_eq!(s.mbr_donor_fdr(), 0.005);
        assert_eq!(s.quant_transfer_fdr(), 0.06);
        assert_eq!(s.coarse_transfer_fdr(), 0.06);
        assert_eq!(s.quant_channel_fdr(), 0.07);
        assert_eq!(s.lfq_feature_fdr(), 0.03);
    }

    #[test]
    fn quant_stages_inherit_the_same_way() {
        let l = LibraryWrap::try_parse_from([
            "andes",
            "--quant-dir",
            "q",
            "--psms",
            "t",
            "--decoy-psms",
            "d",
        ])
        .unwrap()
        .args;
        assert_eq!((l.psm_fdr, l.mbr_donor_fdr()), (0.01, 0.01));
        let l = LibraryWrap::try_parse_from([
            "andes",
            "--quant-dir",
            "q",
            "--psms",
            "t",
            "--decoy-psms",
            "d",
            "--psm-fdr",
            "0.02",
        ])
        .unwrap()
        .args;
        assert_eq!(l.mbr_donor_fdr(), 0.02);
        let f = FinalizeWrap::try_parse_from(["andes", "--quant-dir", "q", "-o", "x"])
            .unwrap()
            .args;
        assert_eq!(f.quant_transfer_fdr(), 0.01);
        assert_eq!(f.proteins.psm_fdr, 0.01);
        assert_eq!(f.proteins.proteoform_fdr(), 0.05);
        let f = FinalizeWrap::try_parse_from([
            "andes",
            "--quant-dir",
            "q",
            "-o",
            "x",
            "--quant-fdr",
            "0.02",
        ])
        .unwrap()
        .args;
        assert_eq!(f.quant_transfer_fdr(), 0.02);
        let f = FinalizeWrap::try_parse_from([
            "andes",
            "--quant-dir",
            "q",
            "-o",
            "x",
            "--quant-fdr",
            "0.02",
            "--quant-transfer-fdr",
            "0.05",
        ])
        .unwrap()
        .args;
        assert_eq!(f.quant_transfer_fdr(), 0.05);
    }

    #[test]
    fn removed_flags_are_rejected() {
        for flag in [
            "--fdr",
            "--mbr-fdr",
            "--mbr-library-fdr",
            "--label-fdr",
            "--label-mbr-fdr",
            "--lfq-feature-fdr",
            "--quant-feature-fdr",
            "--refine-select-psm-fdr",
        ] {
            assert!(
                SearchWrap::try_parse_from(["andes", flag, "0.01"]).is_err(),
                "{flag}"
            );
        }
    }

    #[test]
    fn short_help_omits_the_advanced_flags() {
        let advanced = [
            "--proteoform-fdr",
            "--mbr-donor-fdr",
            "--quant-transfer-fdr",
            "--quant-channel-fdr",
        ];
        let mut cmd = SearchWrap::command();
        let short = cmd.render_help().to_string();
        let long = cmd.render_long_help().to_string();
        let line = |flag: &str| format!("{flag} <");
        for main in ["--psm-fdr", "--protein-fdr", "--quant-fdr"] {
            assert!(short.contains(&line(main)), "-h misses {main}");
            assert!(long.contains(&line(main)), "--help misses {main}");
        }
        for flag in advanced {
            assert!(!short.contains(&line(flag)), "-h shows {flag}");
            assert!(long.contains(&line(flag)), "--help misses {flag}");
        }
        assert!(!short.contains(ADVANCED_FDR));
        assert!(long.contains(ADVANCED_FDR));
        let mut cmd = FinalizeWrap::command();
        let short = cmd.render_help().to_string();
        assert!(short.contains(&line("--quant-fdr")));
        assert!(!short.contains(&line("--quant-transfer-fdr")));
        assert!(!short.contains(&line("--proteoform-fdr")));
        let long = cmd.render_long_help().to_string();
        assert!(long.contains(&line("--quant-transfer-fdr")));
        assert!(long.contains(&line("--proteoform-fdr")));
    }
}

#[cfg(test)]
mod dir_or_output_dir_tests {
    use super::*;

    #[test]
    fn the_flag_wins() {
        assert_eq!(
            dir_or_output_dir(Some(Path::new("/shared")), Path::new("out/run.pin")),
            PathBuf::from("/shared")
        );
    }

    #[test]
    fn the_default_is_the_output_directory() {
        assert_eq!(
            dir_or_output_dir(None, Path::new("results/run.pin")),
            PathBuf::from("results")
        );
        assert_eq!(
            dir_or_output_dir(None, Path::new("/data/out/run")),
            PathBuf::from("/data/out")
        );
    }

    #[test]
    fn a_bare_output_name_means_the_working_directory() {
        assert_eq!(
            dir_or_output_dir(None, Path::new("run.pin")),
            PathBuf::from(".")
        );
    }
}

#[cfg(test)]
mod spelling_tests {
    use super::*;
    use clap::Parser;

    #[derive(Parser, Debug)]
    struct Wrap {
        #[command(flatten)]
        args: SearchArgs,
    }

    fn parse(flag: &str, value: &str) -> Result<SearchArgs, clap::Error> {
        Wrap::try_parse_from(["andes", flag, value]).map(|w| w.args)
    }

    fn accepts(flag: &str, values: &[&str]) {
        for v in values {
            if let Err(e) = parse(flag, v) {
                panic!("{flag} {v:?} rejected: {e}");
            }
        }
    }

    fn rejects(flag: &str, values: &[&str]) {
        for v in values {
            assert!(parse(flag, v).is_err(), "{flag} {v:?} accepted");
        }
    }

    #[test]
    fn listed_flags_accept_the_spellings_of_their_parsers() {
        accepts(
            "--decoy-strategy",
            &[
                "shuffle", "Shuffle", "shuf", "SHUF", "reverse", "rev", "REV", " reverse ",
                "sequon", "SEQUON", "none", "off", "NONE",
            ],
        );
        rejects("--decoy-strategy", &["bogus", "", "random"]);
        assert_eq!(
            parse("--decoy-strategy", "rev").unwrap().decoy_strategy,
            DecoyStrategy::Reverse
        );
        assert_eq!(
            parse("--decoy-strategy", "off").unwrap().decoy_strategy,
            DecoyStrategy::None
        );

        accepts(
            "--enzyme",
            &[
                "trypsin", "Trypsin", "tryp", "chymotrypsin", "chymo", "lysc", "lys-c", "aspn",
                "Asp-N", "gluc", "glu-c", "lysn", "lys-n", "argc", "arg-c", "alp", "alpha-lp",
                "alphalp", "nocleavage", "none", "nonspecific", "all", "elastase",
                "gluc,trypsin", "gluc+trypsin", "  lysc , , trypsin ",
            ],
        );
        rejects("--enzyme", &["bogus", "", " , ", "trypsin,bogus"]);
        assert_eq!(
            parse("--enzyme", "gluc+trypsin").unwrap().db.enzyme,
            "gluc+trypsin",
            "the value is kept as given"
        );

        accepts(
            "--tmt",
            &[
                "tmt6", "TMT10", "tmt10plex", "TMT11plex", "tmt16", "tmtpro", "tmtpro16",
                "tmt18", "tmtpro18", "tmt32", "tmtpro32", "tmt35", "tmtpro35", "itraq4",
                "iTRAQ8plex",
            ],
        );
        rejects("--tmt", &["tmt9", "", "tmt"]);

        accepts("--uniprot-evidence", &["experimental", "Curated", "ANY"]);
        rejects("--uniprot-evidence", &["strong", ""]);

        accepts("--precursor-cal", &["off", "OFF", "Auto", "on"]);
        rejects("--precursor-cal", &["bogus", "1", ""]);
        assert_eq!(
            parse("--precursor-cal", "On").unwrap().precursor_cal,
            PrecursorCalMode::On
        );
    }

    #[test]
    fn case_insensitive_enum_flags() {
        accepts(
            "--fragmentation",
            &["auto", "AUTO", "CID", "cid", "ETD", "etd", "HCD", "Hcd", "UVPD", "uvpd"],
        );
        rejects("--fragmentation", &["3", "99", "", "EThcD"]);
        assert_eq!(
            parse("--fragmentation", "hcd").unwrap().fragmentation,
            Fragmentation::Hcd
        );
        accepts(
            "--protocol",
            &[
                "auto", "phospho", "PHOSPHO", "iTRAQ", "itraq", "iTRAQ-phospho", "ITRAQ-PHOSPHO",
                "TMT", "tmt", "standard",
            ],
        );
        rejects("--protocol", &["4", "42", "", "silac"]);
        accepts(
            "--enzyme-specificity",
            &["non-specific", "NON-SPECIFIC", "semi", "Semi", "fully", "Fully"],
        );
        rejects("--enzyme-specificity", &["2", "nonspecific", ""]);
    }

    #[test]
    fn exact_enum_flags_are_unchanged() {
        accepts("--score", &["auto", "rank", "strong"]);
        rejects("--score", &["Rank", "STRONG"]);
        accepts("--fragment-index", &["auto", "on", "off"]);
        rejects("--fragment-index", &["OFF"]);
        accepts("--tmt-unmix", &["auto", "off", "on"]);
        rejects("--tmt-unmix", &["Auto"]);
        accepts("--predict", &["rt", "ms2", "peptdeep-ms2", "peptdeep-rt", "deeplc-rt", "none", "rt,ms2"]);
        rejects("--predict", &["RT", "ms3"]);
        accepts("--predictor-backend", &["tract", "ort"]);
        rejects("--predictor-backend", &["ORT"]);
        accepts("--glyco-core", &["n", "o"]);
        rejects("--glyco-core", &["N"]);
        accepts("--glyco-species", &["human", "mouse", "mouse-large", "high-mannose"]);
        rejects("--glyco-species", &["Human", "rat"]);
        accepts("--precursor-mono", &["auto", "off"]);
        rejects("--precursor-mono", &["on"]);
        accepts("--ethcd-activation", &["ethcd", "etd"]);
        rejects("--ethcd-activation", &["hcd"]);
    }

    #[test]
    fn config_enum_parser_names_the_values() {
        let e = parse_enum::<PredictFeature>("predict entry", "ms3", true).unwrap_err();
        assert_eq!(
            e,
            "invalid predict entry `ms3`: expected rt|ms2|peptdeep-ms2|peptdeep-rt|deeplc-rt|none"
        );
        assert_eq!(
            parse_enum::<PredictorBackend>("predictor_backend", "ORT", true).unwrap(),
            PredictorBackend::Ort
        );
        assert!(parse_enum::<ScoreFlag>("score", "Rank", false).is_err());
    }
}

#[cfg(test)]
mod help_layout_tests {
    use super::*;
    use clap::{CommandFactory, Parser};
    use std::collections::BTreeSet;

    #[derive(Parser, Debug)]
    struct SearchWrap {
        #[command(flatten)]
        args: SearchArgs,
    }

    #[derive(Parser, Debug)]
    struct InferenceWrap {
        #[command(flatten)]
        args: InferenceArgs,
    }

    #[derive(Parser, Debug)]
    struct FinalizeWrap {
        #[command(flatten)]
        args: QuantFinalizeArgs,
    }

    fn listed(help: &str) -> Vec<String> {
        help.lines()
            .filter(|l| l.starts_with("      --"))
            .map(|l| {
                l.trim_start()
                    .split([' ', '='])
                    .next()
                    .unwrap()
                    .to_string()
            })
            .collect()
    }

    const MAIN: &[&str] = &[
        "--config",
        "--sdrf",
        "--spectrum",
        "--database",
        "--output-pin",
        "--output-tsv",
        "--output-parquet",
        "--precursor-tol",
        "--enzyme",
        "--fixed-mods",
        "--variable-mods",
        "--decoy-strategy",
        "--threads",
        "--psm-fdr",
        "--protein-fdr",
        "--quant-fdr",
        "--rescore",
        "--chimeric",
        "--dia",
        "--refine",
        "--glyco",
        "--glyco-species",
        "--tmt",
        "--lfq",
        "--label",
        "--mbr",
        "--quant-dir",
    ];

    const HIDDEN: &[&str] = &[
        "decoy_seed",
        "cal_min_spec_keys",
        "max_spectra",
        "precursor_mono_dump",
        "glyco_scans",
        "debug_glyco",
        "glyco_mbr_match_train",
        "ipf_debug",
        "xl_oracle",
        "xl_rescore_intra",
        "xl_rescore_inter",
        "optimize_sample",
        "rescore_memory_gb",
    ];

    #[test]
    fn short_help_lists_the_main_options_under_five_headings() {
        let mut cmd = SearchWrap::command();
        let short = cmd.render_help().to_string();
        assert_eq!(listed(&short), MAIN, "{short}");
        for heading in [H_IO, H_SEARCH, H_FDR, H_MODES, H_QUANT] {
            assert!(short.contains(&format!("{heading}:")), "-h misses {heading}");
        }
        for heading in [
            H_DATABASE,
            H_RESCORING,
            H_RUNTIME,
            H_PREDICTION,
            H_ISOBARIC,
            H_LFQ,
            H_LABEL,
            H_GLYCO,
            ADVANCED_FDR,
        ] {
            assert!(!short.contains(&format!("{heading}:")), "-h shows {heading}");
        }
    }

    #[test]
    fn long_help_lists_every_option_but_the_hidden_ones() {
        let mut cmd = SearchWrap::command();
        let long = cmd.render_long_help().to_string();
        let listed: BTreeSet<String> = listed(&long).into_iter().collect();
        let mut hidden = BTreeSet::new();
        for arg in cmd.get_arguments() {
            let id = arg.get_id().as_str();
            if matches!(id, "help" | "version") {
                continue;
            }
            assert!(
                arg.is_hide_set() || arg.get_help_heading().is_some(),
                "{id} has no heading"
            );
            let help = arg.get_help().map(|h| h.to_string()).unwrap_or_default();
            assert!(!help.trim().is_empty(), "{id} has no help");
            let flag = format!("--{}", arg.get_long().unwrap());
            if arg.is_hide_set() {
                hidden.insert(id.to_string());
                assert!(!listed.contains(&flag), "--help lists hidden {flag}");
            } else {
                assert!(listed.contains(&flag), "--help misses {flag}");
            }
        }
        let expected: BTreeSet<String> = HIDDEN.iter().map(|s| s.to_string()).collect();
        assert_eq!(hidden, expected);
        for text in [
            "[default: 0.01]",
            "[default: inherits --psm-fdr]",
            "[default: inherits --quant-fdr]",
        ] {
            assert!(long.contains(text), "--help misses {text}");
        }
    }

    #[test]
    fn help_text_is_plain() {
        for mut cmd in [
            SearchWrap::command(),
            InferenceWrap::command(),
            FinalizeWrap::command(),
        ] {
            let long = cmd.render_long_help().to_string();
            for bad in ["**", "#64", "bigbio/", "@1%", "action = Set", "Every other engine"] {
                assert!(!long.contains(bad), "help contains {bad:?}");
            }
        }
    }

    #[test]
    fn modification_lists_parse_and_suggest() {
        let a = SearchWrap::try_parse_from([
            "andes",
            "--fixed-mods",
            "Carbamidomethyl (C), TMT10plex (K), TMT10plex (N-term)",
            "--variable-mods",
            "Oxidation (M), UNIMOD:21 (STY), Gln->pyro-Glu (N-term Q)",
        ])
        .unwrap()
        .args;
        let fixed = a.db.fixed_mods.unwrap().0;
        assert_eq!(fixed.len(), 3);
        assert_eq!(fixed[1].name, "TMT6plex");
        assert_eq!(a.db.variable_mods.unwrap().0.len(), 5);
        let none = SearchWrap::try_parse_from(["andes", "--variable-mods", ""]).unwrap();
        assert!(none.args.db.variable_mods.unwrap().0.is_empty());
        let e = SearchWrap::try_parse_from(["andes", "--variable-mods", "Oxidaton (M)"])
            .unwrap_err()
            .to_string();
        assert!(e.contains("did you mean Oxidation?"), "{e}");
    }

    #[test]
    fn database_options_have_the_same_help_everywhere() {
        let help_of = |cmd: &clap::Command, id: &str| {
            cmd.get_arguments()
                .find(|a| a.get_id().as_str() == id)
                .and_then(|a| a.get_help().map(|h| h.to_string()))
                .unwrap_or_default()
        };
        let search = SearchWrap::command();
        for cmd in [InferenceWrap::command(), FinalizeWrap::command()] {
            for arg in cmd.get_arguments() {
                let id = arg.get_id().as_str();
                if matches!(id, "help" | "version") {
                    continue;
                }
                let help = help_of(&cmd, id);
                assert!(!help.trim().is_empty(), "{id} has no help");
            }
            for id in [
                "enzyme",
                "fixed_mods",
                "variable_mods",
                "decoy_prefix",
                "decoy_suffix",
                "uniprot_variants",
                "uniprot_isoforms",
                "uniprot_chains",
                "uniprot_evidence",
                "annotations",
                "knockoffs",
            ] {
                assert_eq!(help_of(&cmd, id), help_of(&search, id), "{id}");
            }
        }
    }
}
