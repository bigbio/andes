//! The search command's parameter table: what clap cannot hold about each argument.
//!
//! clap declares the flag, its value type, parser, default and help text. This table adds the
//! config-file section and key, the unit, the accepted range as short text, what an unset value
//! means when clap has no default, the legacy form of the flag, and an optional example. One row
//! exists for every argument of the search command, and the tests in
//! `tests/params_table.rs` keep the table and the CLI in step. The table is read by tests today;
//! later consumers (documentation, `--print-config`) read it through [`spec`] and [`PARAMS`].

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
    Minutes,
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
    /// The accepted values as short text, such as `">= 1"` or `"0..1"`.
    pub range: Option<&'static str>,
    pub unset: Unset,
    /// The Java-form text of the DOCS.md "Legacy form" column.
    pub legacy: Option<&'static str>,
    /// A worked example value for documentation.
    pub example: Option<&'static str>,
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
    }
}

/// The row of the argument with clap id `id`.
pub fn spec(id: &str) -> Option<&'static ParamSpec> {
    PARAMS.iter().find(|p| p.id == id)
}

/// Every argument of the search command, grouped by config section.
#[rustfmt::skip]
pub const PARAMS: &[ParamSpec] = &[
    row("uniprot_variants", Section::Io, "uniprot_variants", Unit::None, None, Unset::Off, None),
    row("uniprot_isoforms", Section::Io, "uniprot_isoforms", Unit::None, None, Unset::Off, None),
    row("uniprot_chains", Section::Io, "uniprot_chains", Unit::None, None, Unset::Off, None),
    row("uniprot_evidence", Section::Io, "uniprot_evidence", Unit::None, Some("experimental|curated|any"), Unset::ClapDefault, None),
    row("annotations", Section::Io, "annotations", Unit::Path, None, Unset::Off, None),
    row("knockoffs", Section::Io, "knockoffs", Unit::Count, None, Unset::ClapDefault, None),
    row("sdrf", Section::Io, "sdrf", Unit::Path, None, Unset::Off, None),
    row("spectrum", Section::Io, "spectrum", Unit::Path, None, Unset::Required, Some("Java -s <FILE>")),
    row("database", Section::Io, "database", Unit::Path, None, Unset::Required, Some("Java -d <FILE>")),
    row("output_pin", Section::Io, "output_pin", Unit::Path, None, Unset::Auto("required unless --rescore"), Some("Java -o <FILE> (when -outputFormat pin)")),
    row("output_tsv", Section::Io, "output_tsv", Unit::Path, None, Unset::Off, Some("Java -outputFormat 1 with output path")),
    row("output_parquet", Section::Io, "output_parquet", Unit::Path, None, Unset::Off, None),
    row("temp_dir", Section::Io, "temp_dir", Unit::Path, None, Unset::Auto("the directory of the outputs"), None),
    row("enzyme", Section::Search, "enzyme", Unit::None, Some("comma list of protease names"), Unset::ClapDefault, Some("Java -e")),
    row("fixed_mods", Section::Search, "fixed_mods", Unit::None, Some("list of Name (sites)"), Unset::Auto("Carbamidomethyl (C)"), Some("Java -mod <FILE>")),
    row("variable_mods", Section::Search, "variable_mods", Unit::None, Some("list of Name (sites)"), Unset::Auto("Oxidation (M), Acetyl (Protein N-term)"), Some("Java -mod <FILE>")),
    row("precursor_tol", Section::Search, "precursor_tol", Unit::PpmOrDa, Some("> 0, VALUE+unit, e.g. 20ppm or 0.02da"), Unset::ClapDefault, Some("Java -t 20ppm")),
    row("charge", Section::Search, "charge", Unit::Charge, Some("MIN..MAX, 1..50, MIN <= MAX"), Unset::ClapDefault, Some("no direct Java flag; set via param file in Java")),
    row("isotope_error", Section::Search, "isotope_error", Unit::None, Some("MIN..MAX, MIN <= MAX"), Unset::Auto("-1..2, 0..2 with --glyco or --xl"), Some("Java -ti -1,2")),
    row("open_search", Section::Search, "open_search", Unit::None, None, Unset::Off, None),
    row("open_search_delta", Section::Search, "open_search_delta", Unit::Da, Some("MIN..MAX in Da, contains 0, width <= 2000"), Unset::ClapDefault, None),
    row("open_search_top_k", Section::Search, "open_search_top_k", Unit::Count, None, Unset::ClapDefault, None),
    row("optimize", Section::Search, "optimize", Unit::None, None, Unset::Off, None),
    row("optimize_sample", Section::Search, "optimize_sample", Unit::Count, None, Unset::ClapDefault, None),
    row("max_missed_cleavages", Section::Search, "max_missed_cleavages", Unit::Count, None, Unset::ClapDefault, Some("Java -maxMissedCleavages 2")),
    row("enzyme_specificity", Section::Search, "enzyme_specificity", Unit::None, Some("enum"), Unset::ClapDefault, Some("Java -ntt")),
    row("min_length", Section::Search, "min_length", Unit::Count, None, Unset::ClapDefault, Some("Java -minLength 6")),
    row("max_length", Section::Search, "max_length", Unit::Count, None, Unset::ClapDefault, Some("Java -maxLength 40")),
    row("max_mods", Section::Search, "max_mods", Unit::Count, None, Unset::ClapDefault, None),
    row("min_peaks", Section::Search, "min_peaks", Unit::Count, None, Unset::ClapDefault, Some("Java -minNumPeaks 10")),
    row("top_n", Section::Search, "top_n", Unit::Count, None, Unset::ClapDefault, Some("Java -n 10")),
    row("dia", Section::Search, "dia", Unit::None, None, Unset::Off, None),
    row("dia_mobility_slices", Section::Search, "dia_mobility_slices", Unit::Count, None, Unset::ClapDefault, None),
    row("dia_ms1_tol_ppm", Section::Search, "dia_ms1_tol_ppm", Unit::Ppm, None, Unset::ClapDefault, None),
    row("localize", Section::Search, "localize", Unit::None, Some("ipf"), Unset::Off, None),
    row("xl", Section::Search, "xl", Unit::None, None, Unset::Off, None),
    row("fragment_index", Section::Search, "fragment_index", Unit::None, Some("enum"), Unset::ClapDefault, None),
    row("fragment_index_slice_da", Section::Search, "fragment_index_slice_da", Unit::Da, None, Unset::Auto("from the memory budget"), None),
    row("index_cache_dir", Section::Search, "index_cache_dir", Unit::Path, None, Unset::Auto("the directory of the outputs"), None),
    row("threads", Section::Scoring, "threads", Unit::Count, None, Unset::ClapDefault, Some("Java -thread N")),
    row("precursor_cal", Section::Scoring, "precursor_cal", Unit::None, Some("off|auto|on"), Unset::ClapDefault, Some("Java -precursorCal auto|on|off")),
    row("score", Section::Scoring, "score", Unit::None, Some("enum"), Unset::ClapDefault, None),
    row("fragmentation", Section::Scoring, "fragmentation", Unit::None, Some("enum"), Unset::ClapDefault, Some("Java -m")),
    row("protocol", Section::Scoring, "protocol", Unit::None, Some("enum"), Unset::ClapDefault, Some("Java -protocol")),
    row("fragment_tol_ppm", Section::Scoring, "fragment_tol_ppm", Unit::Ppm, Some("> 0"), Unset::Auto("measured per run"), None),
    row("fragment_tol_da", Section::Scoring, "fragment_tol_da", Unit::Da, Some("> 0"), Unset::Auto("measured per run"), None),
    row("peak_filter", Section::Scoring, "peak_filter", Unit::None, Some("WINDOW_DA:PEAKS"), Unset::Auto("100:20 for isobaric data, else off"), None),
    row("ethcd_activation", Section::Scoring, "ethcd_activation", Unit::None, Some("enum"), Unset::ClapDefault, None),
    row("ms_level", Section::Scoring, "ms_level", Unit::Count, None, Unset::ClapDefault, None),
    row("model_store", Section::Scoring, "model_store", Unit::Path, None, Unset::Auto("the bundled store"), None),
    row("model_id_override", Section::Scoring, "model", Unit::None, None, Unset::Auto("selected from the data"), None),
    row("gbdt_max_trees", Section::Scoring, "gbdt_max_trees", Unit::Count, None, Unset::Auto("100; all with --glyco"), None),
    row("max_spectra", Section::Scoring, "max_spectra", Unit::Count, None, Unset::ClapDefault, None),
    row("precursor_mono", Section::Scoring, "precursor_mono", Unit::None, Some("enum"), Unset::ClapDefault, None),
    row("decoy_prefix", Section::Decoys, "prefix", Unit::None, None, Unset::ClapDefault, None),
    row("decoy_suffix", Section::Decoys, "suffix", Unit::None, None, Unset::Off, None),
    row("decoy_strategy", Section::Decoys, "strategy", Unit::None, Some("shuffle|reverse|sequon|mutate|none"), Unset::ClapDefault, None),
    row("decoy_seed", Section::Decoys, "seed", Unit::Count, None, Unset::ClapDefault, None),
    row("chimeric", Section::Chimeric, "enabled", Unit::None, None, Unset::Off, None),
    row("chimeric_max_coisolated", Section::Chimeric, "max_coisolated", Unit::Count, None, Unset::ClapDefault, None),
    row("chimeric_max_kl", Section::Chimeric, "max_kl", Unit::Number, Some(">= 0"), Unset::ClapDefault, None),
    row("refine", Section::Refine, "enabled", Unit::None, None, Unset::Off, None),
    row("refine_config", Section::Refine, "config", Unit::Path, None, Unset::Auto("the built-in tier"), None),
    row("psm_fdr", Section::Rescoring, "psm_fdr", Unit::Fraction, Some("0..1"), Unset::Auto("0.01; no filtered output without --rescore"), None),
    row("protein_fdr", Section::Rescoring, "protein_fdr", Unit::Fraction, Some("0..1"), Unset::Off, None),
    row("rescore", Section::Rescoring, "enabled", Unit::None, None, Unset::Off, None),
    row("rescore_memory_gb", Section::Rescoring, "memory_gb", Unit::Gb, Some(">= 1"), Unset::Auto("what the PIN needs, at most 16 GB and the free memory less a margin"), None),
    row("pep", Section::Rescoring, "pep", Unit::Fraction, Some("0..1"), Unset::Off, None),
    row("percolator_bin", Section::Rescoring, "percolator_bin", Unit::Path, None, Unset::Off, None),
    row("percolator_docker", Section::Rescoring, "percolator_docker", Unit::None, None, Unset::Off, None),
    row("percolator_image", Section::Rescoring, "percolator_image", Unit::None, None, Unset::ClapDefault, None),
    row("percolator_args", Section::Rescoring, "percolator_args", Unit::None, None, Unset::ClapDefault, None),
    row("keep_pin", Section::Rescoring, "keep_pin", Unit::None, None, Unset::ClapDefault, None),
    row("predict", Section::Rescoring, "predict", Unit::None, Some("list of predicted features"), Unset::ClapDefault, None),
    row("predictor_models", Section::Rescoring, "predictor_models", Unit::Path, None, Unset::Off, None),
    row("predictor_backend", Section::Rescoring, "predictor_backend", Unit::None, Some("enum"), Unset::ClapDefault, None),
    row("predictor_nce", Section::Rescoring, "predictor_nce", Unit::Number, None, Unset::Auto("28 on Q Exactive, else 30"), None),
    row("predictor_instrument", Section::Rescoring, "predictor_instrument", Unit::None, None, Unset::Auto("from the instrument"), None),
    row("predictor_batch", Section::Rescoring, "predictor_batch", Unit::Count, None, Unset::ClapDefault, None),
    row("proteoform_fdr", Section::Rescoring, "proteoform_fdr", Unit::Fraction, Some("0..1"), Unset::Auto("the larger of --psm-fdr and 0.05"), None),
    row("glyco", Section::Glyco, "enabled", Unit::None, None, Unset::Off, None),
    row("glyco_species", Section::Glyco, "species", Unit::None, Some("enum"), Unset::Off, None),
    row("glyco_glycan_gdb", Section::Glyco, "glycan_gdb", Unit::Path, None, Unset::Off, None),
    row("glyco_glycan_wurcs", Section::Glyco, "glycan_wurcs", Unit::Path, None, Unset::Off, None),
    row("glyco_core", Section::Glyco, "core", Unit::None, Some("enum"), Unset::ClapDefault, None),
    row("glyco_tol_ppm", Section::Glyco, "tol_ppm", Unit::Ppm, Some("> 0"), Unset::ClapDefault, None),
    row("glyco_min_core_y", Section::Glyco, "min_core_y", Unit::Count, None, Unset::ClapDefault, None),
    row("glyco_hcd_pair", Section::Glyco, "hcd_pair", Unit::None, None, Unset::ClapDefault, None),
    row("glyco_max_peaks", Section::Glyco, "max_peaks", Unit::Count, None, Unset::ClapDefault, None),
    row("glyco_cz_max_charge", Section::Glyco, "cz_max_charge", Unit::Charge, None, Unset::Auto("from whether the spectrum is deconvoluted"), None),
    row("glyco_y_max_charge", Section::Glyco, "y_max_charge", Unit::Charge, None, Unset::ClapDefault, None),
    row("glyco_retrieval_tol_ppm", Section::Glyco, "retrieval_tol_ppm", Unit::Ppm, Some("> 0"), Unset::Auto("--glyco-tol-ppm on high-res MS2, 0.5 Da on low-res"), None),
    row("glyco_peptide_first", Section::Glyco, "peptide_first", Unit::None, None, Unset::Off, None),
    row("glyco_index_sequon_only", Section::Glyco, "index_sequon_only", Unit::None, None, Unset::Off, None),
    row("debug_glyco", Section::Glyco, "debug", Unit::None, None, Unset::Off, None),
    row("quant_fdr", Section::Quant, "quant_fdr", Unit::Fraction, Some("0..1"), Unset::Auto("0.01; 0.05 for --label and transfers"), None),
    row("tmt", Section::Quant, "tmt", Unit::None, Some("plex name, e.g. tmt10"), Unset::Off, None),
    row("lfq", Section::Quant, "lfq", Unit::None, None, Unset::Off, None),
    row("label", Section::Quant, "label", Unit::None, Some("preset or one bracket per channel"), Unset::Off, None),
    row("mbr", Section::Quant, "mbr", Unit::None, None, Unset::Off, None),
    row("quant_dir", Section::Quant, "dir", Unit::Path, None, Unset::Off, None),
    row("tmt_level", Section::Quant, "tmt_level", Unit::Count, Some("2 or 3"), Unset::Auto("3 for SPS-MS3 data, else 2"), None),
    row("tmt_tol", Section::Quant, "tmt_tol", Unit::PpmOrDa, Some("> 0, VALUE+unit, e.g. 20ppm or 0.02da"), Unset::Auto("20ppm (10ppm for tmt32/tmt35), 0.3da on ion-trap MS2"), None),
    row("tmt_correction", Section::Quant, "tmt_correction", Unit::Path, None, Unset::Off, None),
    row("tmt_min_purity", Section::Quant, "tmt_min_purity", Unit::Fraction, Some("0..1"), Unset::ClapDefault, None),
    row("tmt_min_sn", Section::Quant, "tmt_min_sn", Unit::Number, Some(">= 0"), Unset::ClapDefault, None),
    row("tmt_unmix", Section::Quant, "tmt_unmix", Unit::None, Some("auto|off|on"), Unset::ClapDefault, None),
    row("mbr_wrong_log2", Section::Quant, "mbr_wrong_log2", Unit::Number, Some("> 0"), Unset::ClapDefault, None),
    row("lfq_tol", Section::Quant, "lfq_tol", Unit::PpmOrDa, Some("> 0, VALUE+unit, e.g. 20ppm or 0.02da"), Unset::Auto("10ppm"), None),
    row("lfq_rt_window", Section::Quant, "lfq_rt_window", Unit::Seconds, Some("> 0"), Unset::Auto("60 s, 120 s with --glyco"), None),
    row("lfq_gaussian_width", Section::Quant, "lfq_gaussian_width", Unit::Seconds, Some("> 0"), Unset::Auto("Savitzky-Golay smoothing"), None),
    row("lfq_min_cosine", Section::Quant, "lfq_min_cosine", Unit::Fraction, Some("0..1"), Unset::ClapDefault, None),
    row("label_reference", Section::Quant, "label_reference", Unit::Count, Some(">= 1"), Unset::ClapDefault, None),
    row("label_arg_pro", Section::Quant, "label_arg_pro", Unit::None, None, Unset::Off, None),
    row("label_2h_enrichment_silac", Section::Quant, "label_2h_enrichment_silac", Unit::Fraction, Some("0..1"), Unset::ClapDefault, None),
    row("label_2h_enrichment_chemical", Section::Quant, "label_2h_enrichment_chemical", Unit::Fraction, Some("0..1"), Unset::ClapDefault, None),
    row("label_heavy_enrichment", Section::Quant, "label_heavy_enrichment", Unit::Fraction, Some("0..1"), Unset::ClapDefault, None),
    row("mbr_donor_fdr", Section::Quant, "mbr_donor_fdr", Unit::Fraction, Some("0..1"), Unset::Auto("inherits --psm-fdr"), None),
    row("quant_transfer_fdr", Section::Quant, "quant_transfer_fdr", Unit::Fraction, Some("0..1"), Unset::Auto("inherits --quant-fdr"), None),
    row("quant_channel_fdr", Section::Quant, "quant_channel_fdr", Unit::Fraction, Some("0..1"), Unset::Auto("inherits --quant-fdr"), None),
    row("config", Section::CliOnly, "config", Unit::Path, None, Unset::Off, None),
    row("ipf_debug", Section::CliOnly, "ipf_debug", Unit::None, None, Unset::Off, None),
    row("xl_oracle", Section::CliOnly, "xl_oracle", Unit::Count, None, Unset::ClapDefault, None),
    row("xl_linear_only", Section::CliOnly, "xl_linear_only", Unit::None, None, Unset::Off, None),
    row("xl_rescore_intra", Section::CliOnly, "xl_rescore_intra", Unit::None, None, Unset::Off, None),
    row("xl_rescore_inter", Section::CliOnly, "xl_rescore_inter", Unit::None, None, Unset::Off, None),
    row("cal_min_spec_keys", Section::CliOnly, "cal_min_spec_keys", Unit::Count, None, Unset::Auto("10000"), None),
    row("precursor_mono_dump", Section::CliOnly, "precursor_mono_dump", Unit::Path, None, Unset::Off, None),
    row("glyco_scans", Section::CliOnly, "glyco_scans", Unit::Path, None, Unset::Off, None),
    row("glyco_mbr_match_train", Section::CliOnly, "glyco_mbr_match_train", Unit::Path, None, Unset::Off, None),
    row("quiet", Section::CliOnly, "quiet", Unit::None, None, Unset::Off, None),
    row("verbose", Section::CliOnly, "verbose", Unit::None, None, Unset::Off, None),
];
