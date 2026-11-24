//! The search command's parameter table: what clap cannot hold about each argument.
//!
//! clap declares the flag, its value type, parser, default and help text. This table adds the
//! config-file section and key, the unit, the accepted range as short text, what an unset value
//! means when clap has no default, the legacy form of the flag, and an optional example. One row
//! exists for every argument of the search command, and the tests in
//! `tests/params_table.rs` keep the table and the CLI in step. `params_render` reads the table
//! through [`spec`] and [`PARAMS`] to write `config.example.yaml`, the DOCS.md flag tables and
//! `--print-config`.

pub(crate) mod render;
pub(crate) mod sources;

/// The `RunConfig` section a parameter belongs to; the YAML name is the field name of
/// `RunConfig`. `CliOnly` arguments have no config key.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Section {
    Io,
    Search,
    Scoring,
    Decoys,
    Chimeric,
    Refine,
    Rescoring,
    Glyco,
    Quant,
    CliOnly,
}

impl Section {
    /// The YAML section name, or `None` for a command-line-only argument.
    pub fn yaml_name(self) -> Option<&'static str> {
        match self {
            Section::Io => Some("io"),
            Section::Search => Some("search"),
            Section::Scoring => Some("scoring"),
            Section::Decoys => Some("decoys"),
            Section::Chimeric => Some("chimeric"),
            Section::Refine => Some("refine"),
            Section::Rescoring => Some("rescoring"),
            Section::Glyco => Some("glyco"),
            Section::Quant => Some("quant"),
            Section::CliOnly => None,
        }
    }
}

/// The unit of a numeric value, or the kind of a path; `None` for a name, a list or a switch.
/// `Count` is an integer, `Number` a plain real number without a unit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Unit {
    None,
    Ppm,
    Da,
    PpmOrDa,
    Seconds,
    Gb,
    Fraction,
    Count,
    Number,
    Charge,
    Path,
}

/// What an argument means when no value is given and clap holds no default.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Unset {
    /// clap holds the default; `--help` shows it.
    ClapDefault,
    /// The value is worked out at run time; the text says how, for documentation.
    Auto(&'static str),
    /// The feature is off.
    Off,
    /// The run cannot start without a value.
    Required,
}

/// One search argument: the clap id and what clap cannot declare about it.
#[derive(Debug, Clone, Copy)]
pub struct ParamSpec {
    pub id: &'static str,
    pub section: Section,
    /// The YAML key inside the section, usually the id.
    pub key: &'static str,
    pub unit: Unit,
    /// The accepted values as short text, such as `">= 1"` or `"0..1"`; `None` for an
    /// argument with clap possible values, which the renderers list from clap.
    pub range: Option<&'static str>,
    pub unset: Unset,
    /// The Java-form text of the DOCS.md "Legacy form" column.
    pub legacy: Option<&'static str>,
    /// A worked example value for documentation.
    pub example: Option<&'static str>,
    /// Comment lines written under the row in `config.example.yaml`, without the `# ` prefix:
    /// the worked multi-line forms of a config-only setting.
    pub example_block: Option<&'static str>,
}

const fn row(
    id: &'static str,
    section: Section,
    key: &'static str,
    unit: Unit,
    range: Option<&'static str>,
    unset: Unset,
    legacy: Option<&'static str>,
) -> ParamSpec {
    ParamSpec {
        id,
        section,
        key,
        unit,
        range,
        unset,
        legacy,
        example: None,
        example_block: None,
    }
}

const fn with_example(mut p: ParamSpec, example: &'static str) -> ParamSpec {
    p.example = Some(example);
    p
}

const fn with_block(mut p: ParamSpec, block: &'static str) -> ParamSpec {
    p.example_block = Some(block);
    p
}

/// The row of the argument with clap id `id`.
pub fn spec(id: &str) -> Option<&'static ParamSpec> {
    PARAMS.iter().find(|p| p.id == id)
}

const MODS_BLOCK: &str = "\
A list item is the command-line string or a map that adds the config-only settings (example):
variable_mods:
  - Oxidation (M)
  - {name: Phospho, sites: STY, neutral_loss: [97.976896], loss_class: phospho}
  - {name: GlyGly, sites: K, cleavage: blocked, diagnostic_ions: [215.1502, 198.1237]}
  - {name: Phospho, sites: Y, restrict_to: annotated}   # only at sites a UniProt XML database annotates
  - {name: MyTag, mass: 123.4567, sites: N-term}         # a custom mass, no UNIMOD lookup";

const XL_BLOCK: &str = "\
A reagent the table (resources/crosslinkers.yaml) lacks is an entry of that file's schema (example):
xl:
  name: MyLinker
  sites: K, Protein N-term       # one end, or \"end / end\" for a heterobifunctional linker
  intact: {mass: 158.003765}     # or {composition: \"H(6) C(6) O(3) S\"} or {unimod: \"Xlink:DSSO[158]\"}
  stubs: [{name: alkene, mass: 54.010565}, {name: thiol, mass: 85.982635}]   # MS-cleavage stubs
  doublets: [[alkene, thiol]]    # reporter doublets, [light, heavy] stub names
  monolinks: [{name: water, mass: 176.01433}]";

/// Every argument of the search command, grouped by config section.
#[rustfmt::skip]
pub const PARAMS: &[ParamSpec] = &[
    row("uniprot_variants", Section::Io, "uniprot_variants", Unit::None, None, Unset::Off, None),
    row("uniprot_isoforms", Section::Io, "uniprot_isoforms", Unit::None, None, Unset::Off, None),
    row("uniprot_chains", Section::Io, "uniprot_chains", Unit::None, None, Unset::Off, None),
    row("uniprot_evidence", Section::Io, "uniprot_evidence", Unit::None, None, Unset::ClapDefault, None),
    with_example(row("annotations", Section::Io, "annotations", Unit::Path, None, Unset::Off, None), "db.annotations.parquet"),
    row("knockoffs", Section::Io, "knockoffs", Unit::Count, None, Unset::ClapDefault, None),
    with_example(row("sdrf", Section::Io, "sdrf", Unit::Path, None, Unset::Off, None), "design.sdrf.tsv"),
    with_example(row("spectrum", Section::Io, "spectrum", Unit::Path, None, Unset::Required, Some("Java -s <FILE>")), "sample.mzML"),
    with_example(row("database", Section::Io, "database", Unit::Path, None, Unset::Required, Some("Java -d <FILE>")), "human.fasta"),
    with_example(row("output_pin", Section::Io, "output_pin", Unit::Path, None, Unset::Auto("required unless --rescore"), Some("Java -o <FILE> (when -outputFormat pin)")), "out.pin"),
    with_example(row("output_tsv", Section::Io, "output_tsv", Unit::Path, None, Unset::Off, Some("Java -outputFormat 1 with output path")), "out.tsv"),
    with_example(row("output_parquet", Section::Io, "output_parquet", Unit::Path, None, Unset::Off, None), "out.idparquet"),
    with_example(row("temp_dir", Section::Io, "temp_dir", Unit::Path, None, Unset::Auto("the directory of the outputs"), None), "tmp"),
    row("enzyme", Section::Search, "enzyme", Unit::None, None, Unset::ClapDefault, Some("Java -e")),
    with_example(row("fixed_mods", Section::Search, "fixed_mods", Unit::None, Some("list of Name (sites)"), Unset::Auto("Carbamidomethyl (C)"), Some("Java -mod <FILE>")), "Carbamidomethyl (C)"),
    with_block(with_example(row("variable_mods", Section::Search, "variable_mods", Unit::None, Some("list of Name (sites)"), Unset::Auto("Oxidation (M), Acetyl (Protein N-term)"), Some("Java -mod <FILE>")), "Oxidation (M)"), MODS_BLOCK),
    row("precursor_tol", Section::Search, "precursor_tol", Unit::PpmOrDa, Some("> 0, VALUE+unit, e.g. 20ppm or 0.02da"), Unset::ClapDefault, Some("Java -t 20ppm")),
    row("charge", Section::Search, "charge", Unit::Charge, Some("MIN..MAX, 1..50, MIN <= MAX"), Unset::ClapDefault, Some("no direct Java flag; set via param file in Java")),
    with_example(row("isotope_error", Section::Search, "isotope_error", Unit::None, Some("MIN..MAX, MIN <= MAX"), Unset::Auto("-1..2, 0..2 with --glyco or --xl"), Some("Java -ti -1,2")), "0..1"),
    row("open_search", Section::Search, "open_search", Unit::None, None, Unset::Off, None),
    row("open_search_delta", Section::Search, "open_search_delta", Unit::Da, Some("MIN..MAX in Da, contains 0, width <= 2000"), Unset::ClapDefault, None),
    row("open_search_top_k", Section::Search, "open_search_top_k", Unit::Count, None, Unset::ClapDefault, None),
    row("optimize", Section::Search, "optimize", Unit::None, None, Unset::Off, None),
    row("optimize_sample", Section::Search, "optimize_sample", Unit::Count, None, Unset::ClapDefault, None),
    row("max_missed_cleavages", Section::Search, "max_missed_cleavages", Unit::Count, None, Unset::ClapDefault, Some("Java -maxMissedCleavages 2")),
    row("enzyme_specificity", Section::Search, "enzyme_specificity", Unit::None, None, Unset::ClapDefault, Some("Java -ntt")),
    row("min_length", Section::Search, "min_length", Unit::Count, None, Unset::ClapDefault, Some("Java -minLength 6")),
    row("max_length", Section::Search, "max_length", Unit::Count, None, Unset::ClapDefault, Some("Java -maxLength 40")),
    row("max_mods", Section::Search, "max_mods", Unit::Count, None, Unset::ClapDefault, None),
    row("min_peaks", Section::Search, "min_peaks", Unit::Count, None, Unset::ClapDefault, Some("Java -minNumPeaks 10")),
    row("top_n", Section::Search, "top_n", Unit::Count, None, Unset::ClapDefault, Some("Java -n 10")),
    row("dia", Section::Search, "dia", Unit::None, None, Unset::Off, None),
    row("dia_mobility_slices", Section::Search, "dia_mobility_slices", Unit::Count, None, Unset::ClapDefault, None),
    row("dia_ms1_tol_ppm", Section::Search, "dia_ms1_tol_ppm", Unit::Ppm, None, Unset::ClapDefault, None),
    with_example(row("localize", Section::Search, "localize", Unit::None, None, Unset::Off, None), "ipf"),
    with_block(with_example(row("xl", Section::Search, "xl", Unit::None, None, Unset::Off, None), "DSSO (K, Protein N-term)"), XL_BLOCK),
    row("xl_linear_only", Section::Search, "xl_linear_only", Unit::None, None, Unset::Off, None),
    row("fragment_index", Section::Search, "fragment_index", Unit::None, None, Unset::ClapDefault, None),
    with_example(row("fragment_index_slice_da", Section::Search, "fragment_index_slice_da", Unit::Da, None, Unset::Auto("from the memory budget"), None), "80"),
    with_example(row("index_cache_dir", Section::Search, "index_cache_dir", Unit::Path, None, Unset::Auto("the directory of the outputs"), None), "cache"),
    with_example(row("threads", Section::Scoring, "threads", Unit::Count, None, Unset::Auto("the number of logical CPUs"), Some("Java -thread N")), "8"),
    row("precursor_cal", Section::Scoring, "precursor_cal", Unit::None, None, Unset::ClapDefault, Some("Java -precursorCal auto|on|off")),
    row("score", Section::Scoring, "score", Unit::None, None, Unset::ClapDefault, None),
    row("fragmentation", Section::Scoring, "fragmentation", Unit::None, None, Unset::ClapDefault, Some("Java -m")),
    row("protocol", Section::Scoring, "protocol", Unit::None, None, Unset::ClapDefault, Some("Java -protocol")),
    with_example(row("fragment_tol_ppm", Section::Scoring, "fragment_tol_ppm", Unit::Ppm, Some("> 0"), Unset::Auto("measured per run"), None), "20"),
    with_example(row("fragment_tol_da", Section::Scoring, "fragment_tol_da", Unit::Da, Some("> 0"), Unset::Auto("measured per run"), None), "0.5"),
    with_example(row("peak_filter", Section::Scoring, "peak_filter", Unit::None, Some("WINDOW_DA:PEAKS"), Unset::Auto("100:20 for isobaric data, else off"), None), "100:20"),
    row("ethcd_activation", Section::Scoring, "ethcd_activation", Unit::None, None, Unset::ClapDefault, None),
    row("ms_level", Section::Scoring, "ms_level", Unit::Count, None, Unset::ClapDefault, None),
    with_example(row("model_store", Section::Scoring, "model_store", Unit::Path, None, Unset::Auto("the bundled store"), None), "models.parquet"),
    with_example(row("model_id_override", Section::Scoring, "model", Unit::None, None, Unset::Auto("selected from the data"), None), "hcd_qexactive_tryp"),
    with_example(row("gbdt_max_trees", Section::Scoring, "gbdt_max_trees", Unit::Count, None, Unset::Auto("100; all with --glyco"), None), "100"),
    row("max_spectra", Section::Scoring, "max_spectra", Unit::Count, None, Unset::ClapDefault, None),
    row("precursor_mono", Section::Scoring, "precursor_mono", Unit::None, None, Unset::ClapDefault, None),
    row("decoy_prefix", Section::Decoys, "prefix", Unit::None, None, Unset::ClapDefault, None),
    with_example(row("decoy_suffix", Section::Decoys, "suffix", Unit::None, None, Unset::Off, None), "rev"),
    row("decoy_strategy", Section::Decoys, "strategy", Unit::None, None, Unset::ClapDefault, None),
    row("decoy_seed", Section::Decoys, "seed", Unit::Count, None, Unset::ClapDefault, None),
    row("chimeric", Section::Chimeric, "enabled", Unit::None, None, Unset::Off, None),
    row("chimeric_max_coisolated", Section::Chimeric, "max_coisolated", Unit::Count, None, Unset::ClapDefault, None),
    row("chimeric_max_kl", Section::Chimeric, "max_kl", Unit::Number, Some(">= 0"), Unset::ClapDefault, None),
    row("refine", Section::Refine, "enabled", Unit::None, None, Unset::Off, None),
    with_example(row("refine_config", Section::Refine, "config", Unit::Path, None, Unset::Auto("the built-in tier"), None), "refine_tiers.yaml"),
    with_example(row("psm_fdr", Section::Rescoring, "psm_fdr", Unit::Fraction, Some("0..1"), Unset::Auto("0.01"), None), "0.01"),
    with_example(row("protein_fdr", Section::Rescoring, "protein_fdr", Unit::Fraction, Some("0..1"), Unset::Off, None), "0.01"),
    row("rescore", Section::Rescoring, "enabled", Unit::None, None, Unset::Off, None),
    with_example(row("rescore_memory_gb", Section::Rescoring, "memory_gb", Unit::Gb, Some(">= 1"), Unset::Auto("what the PIN needs to stay in memory, at most 16 GB and at most the memory free to the job or machine less a margin"), None), "16"),
    with_example(row("pep", Section::Rescoring, "pep", Unit::Fraction, Some("0..1"), Unset::Off, None), "0.05"),
    with_example(row("percolator_bin", Section::Rescoring, "percolator_bin", Unit::Path, None, Unset::Off, None), "percolator"),
    row("percolator_docker", Section::Rescoring, "percolator_docker", Unit::None, None, Unset::Off, None),
    row("percolator_image", Section::Rescoring, "percolator_image", Unit::None, None, Unset::ClapDefault, None),
    with_example(row("percolator_args", Section::Rescoring, "percolator_args", Unit::None, None, Unset::Auto("no extra arguments"), None), "--testFDR 0.05"),
    row("keep_pin", Section::Rescoring, "keep_pin", Unit::None, None, Unset::ClapDefault, None),
    row("predict", Section::Rescoring, "predict", Unit::None, None, Unset::ClapDefault, None),
    with_example(row("predictor_models", Section::Rescoring, "predictor_models", Unit::Path, None, Unset::Off, None), "models/"),
    row("predictor_backend", Section::Rescoring, "predictor_backend", Unit::None, None, Unset::ClapDefault, None),
    with_example(row("predictor_nce", Section::Rescoring, "predictor_nce", Unit::Number, None, Unset::Auto("28 on Q Exactive, else 30"), None), "28"),
    with_example(row("predictor_instrument", Section::Rescoring, "predictor_instrument", Unit::None, None, Unset::Auto("from the instrument"), None), "QE"),
    row("predictor_batch", Section::Rescoring, "predictor_batch", Unit::Count, None, Unset::ClapDefault, None),
    with_example(row("proteoform_fdr", Section::Rescoring, "proteoform_fdr", Unit::Fraction, Some("0..1"), Unset::Auto("the larger of --psm-fdr and 0.05"), None), "0.05"),
    row("glyco", Section::Glyco, "enabled", Unit::None, None, Unset::Off, None),
    with_example(row("glyco_species", Section::Glyco, "species", Unit::None, None, Unset::Off, None), "human"),
    with_example(row("glyco_glycan_gdb", Section::Glyco, "glycan_gdb", Unit::Path, None, Unset::Off, None), "glycans.gdb"),
    with_example(row("glyco_glycan_wurcs", Section::Glyco, "glycan_wurcs", Unit::Path, None, Unset::Off, None), "glycans.csv"),
    row("glyco_core", Section::Glyco, "core", Unit::None, None, Unset::ClapDefault, None),
    row("glyco_tol_ppm", Section::Glyco, "tol_ppm", Unit::Ppm, Some("> 0"), Unset::ClapDefault, None),
    row("glyco_min_core_y", Section::Glyco, "min_core_y", Unit::Count, None, Unset::ClapDefault, None),
    row("glyco_hcd_pair", Section::Glyco, "hcd_pair", Unit::None, None, Unset::ClapDefault, None),
    row("glyco_max_peaks", Section::Glyco, "max_peaks", Unit::Count, None, Unset::ClapDefault, None),
    with_example(row("glyco_cz_max_charge", Section::Glyco, "cz_max_charge", Unit::Charge, None, Unset::Auto("from whether the spectrum is deconvoluted"), None), "1"),
    row("glyco_y_max_charge", Section::Glyco, "y_max_charge", Unit::Charge, None, Unset::ClapDefault, None),
    with_example(row("glyco_retrieval_tol_ppm", Section::Glyco, "retrieval_tol_ppm", Unit::Ppm, Some("> 0"), Unset::Auto("--glyco-tol-ppm on high-res MS2, 0.5 Da on low-res"), None), "20"),
    row("glyco_peptide_first", Section::Glyco, "peptide_first", Unit::None, None, Unset::Off, None),
    row("glyco_index_sequon_only", Section::Glyco, "index_sequon_only", Unit::None, None, Unset::Off, None),
    row("debug_glyco", Section::Glyco, "debug", Unit::None, None, Unset::Off, None),
    with_example(row("quant_fdr", Section::Quant, "quant_fdr", Unit::Fraction, Some("0..1"), Unset::Auto("0.01"), None), "0.01"),
    with_example(row("tmt", Section::Quant, "tmt", Unit::None, None, Unset::Off, None), "tmt10"),
    row("lfq", Section::Quant, "lfq", Unit::None, None, Unset::Off, None),
    with_example(row("label", Section::Quant, "label", Unit::None, Some("preset or one bracket per channel"), Unset::Off, None), "silac2"),
    row("mbr", Section::Quant, "mbr", Unit::None, None, Unset::Off, None),
    with_example(row("quant_dir", Section::Quant, "dir", Unit::Path, None, Unset::Off, None), "quant_out"),
    with_example(row("tmt_level", Section::Quant, "tmt_level", Unit::Count, Some("2 or 3"), Unset::Auto("3 for SPS-MS3 data, else 2"), None), "3"),
    with_example(row("tmt_tol", Section::Quant, "tmt_tol", Unit::PpmOrDa, Some("> 0, VALUE+unit, e.g. 20ppm or 0.003da"), Unset::Auto("20ppm (10ppm for tmt32/tmt35), 0.3da on ion-trap MS2"), None), "20ppm"),
    with_example(row("tmt_correction", Section::Quant, "tmt_correction", Unit::Path, None, Unset::Off, None), "lot.txt"),
    row("tmt_min_purity", Section::Quant, "tmt_min_purity", Unit::Fraction, Some("0..1"), Unset::ClapDefault, None),
    row("tmt_min_sn", Section::Quant, "tmt_min_sn", Unit::Number, Some(">= 0"), Unset::ClapDefault, None),
    row("tmt_unmix", Section::Quant, "tmt_unmix", Unit::None, None, Unset::ClapDefault, None),
    row("mbr_wrong_log2", Section::Quant, "mbr_wrong_log2", Unit::Number, Some("> 0"), Unset::ClapDefault, None),
    with_example(row("lfq_tol", Section::Quant, "lfq_tol", Unit::PpmOrDa, Some("> 0, VALUE+unit, e.g. 20ppm or 0.02da"), Unset::Auto("10ppm"), None), "10ppm"),
    with_example(row("lfq_rt_window", Section::Quant, "lfq_rt_window", Unit::Seconds, Some("> 0"), Unset::Auto("60; 120 with --glyco"), None), "60"),
    with_example(row("lfq_gaussian_width", Section::Quant, "lfq_gaussian_width", Unit::Seconds, Some("> 0"), Unset::Auto("Savitzky-Golay smoothing"), None), "12"),
    row("lfq_min_cosine", Section::Quant, "lfq_min_cosine", Unit::Fraction, Some("0..1"), Unset::ClapDefault, None),
    row("label_reference", Section::Quant, "label_reference", Unit::Count, Some(">= 1"), Unset::ClapDefault, None),
    row("label_arg_pro", Section::Quant, "label_arg_pro", Unit::None, None, Unset::Off, None),
    row("label_2h_enrichment_silac", Section::Quant, "label_2h_enrichment_silac", Unit::Fraction, Some("0..1"), Unset::ClapDefault, None),
    row("label_2h_enrichment_chemical", Section::Quant, "label_2h_enrichment_chemical", Unit::Fraction, Some("0..1"), Unset::ClapDefault, None),
    row("label_heavy_enrichment", Section::Quant, "label_heavy_enrichment", Unit::Fraction, Some("0..1"), Unset::ClapDefault, None),
    with_example(row("mbr_donor_fdr", Section::Quant, "mbr_donor_fdr", Unit::Fraction, Some("0..1"), Unset::Auto("inherits --psm-fdr"), None), "0.01"),
    with_example(row("quant_transfer_fdr", Section::Quant, "quant_transfer_fdr", Unit::Fraction, Some("0..1"), Unset::Auto("inherits --quant-fdr"), None), "0.01"),
    with_example(row("quant_channel_fdr", Section::Quant, "quant_channel_fdr", Unit::Fraction, Some("0..1"), Unset::Auto("inherits --quant-fdr"), None), "0.05"),
    row("config", Section::CliOnly, "config", Unit::Path, None, Unset::Off, None),
    row("print_config", Section::CliOnly, "print_config", Unit::Path, None, Unset::Off, None),
    row("ipf_debug", Section::CliOnly, "ipf_debug", Unit::None, None, Unset::Off, None),
    row("xl_oracle", Section::CliOnly, "xl_oracle", Unit::Count, None, Unset::ClapDefault, None),
    row("xl_rescore_intra", Section::CliOnly, "xl_rescore_intra", Unit::None, None, Unset::Off, None),
    row("xl_rescore_inter", Section::CliOnly, "xl_rescore_inter", Unit::None, None, Unset::Off, None),
    row("cal_min_spec_keys", Section::CliOnly, "cal_min_spec_keys", Unit::Count, None, Unset::Auto("10000"), None),
    row("precursor_mono_dump", Section::CliOnly, "precursor_mono_dump", Unit::Path, None, Unset::Off, None),
    row("glyco_scans", Section::CliOnly, "glyco_scans", Unit::Path, None, Unset::Off, None),
    row("glyco_mbr_match_train", Section::CliOnly, "glyco_mbr_match_train", Unit::Path, None, Unset::Off, None),
    row("quiet", Section::CliOnly, "quiet", Unit::None, None, Unset::Off, None),
    row("verbose", Section::CliOnly, "verbose", Unit::None, None, Unset::Off, None),
];
