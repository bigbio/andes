//! The command line: help headings, shared defaults and value lists, and the small path
//! helpers; the flag types, the search arguments and the subcommand arguments are in child modules.

mod commands;
mod flags;
mod options;
pub(crate) mod parse;

use std::path::{Path, PathBuf};

pub(crate) use crate::cli::parse::{
    parse_charge_range, parse_isotope_error_range, parse_open_search_delta, parse_precursor_tol,
    parse_unit_fraction, parse_xl,
};
pub(crate) use commands::{
    default_value, CheckDigest, CohortCommand, CohortPlanArgs, CohortRunArgs, CohortXlFdrArgs,
    DatabaseCmdArgs, InferenceArgs, ProteinArgs, QuantCommand, QuantExtractArgs, QuantFinalizeArgs,
    QuantLibraryArgs, RescoreArgs,
};
pub(crate) use flags::{
    decoy_strategy_name, parse_decoy_strategy, parse_enum, parse_enzyme_list, parse_fixed_mods,
    parse_plex, parse_precursor_cal, parse_uniprot_evidence, parse_variable_mods,
    EnzymeSpecificity, EthcdActivationFlag, FragmentIndexFlag, Fragmentation, GlycoCoreFlag,
    GlycoSpeciesFlag, ModList, PrecursorMonoFlag, PredictFeature, PredictorBackend, Protocol,
    ScoreFlag, TmtUnmix, XlSet,
};
pub(crate) use options::{DatabaseArgs, GeometryArgs, SearchArgs};

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

pub(crate) const DEFAULT_GBDT_MAX_TREES: usize = 100;

pub(crate) const DEFAULT_FDR: f64 = 0.01;

pub(crate) const PROTEOFORM_FDR_FLOOR: f64 = 0.05;

pub(crate) const COARSE_QUANT_FDR: f64 = 0.05;

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

/// The `--decoy-strategy` of the search and of `andes database`.
pub(crate) const DEFAULT_DECOY_STRATEGY: &str = "shuffle";

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
        "mutate",
        "substitute one interior residue per tryptic segment of 7 or more (mass change > 1.1 Da), \
         and one in each peptide that still equals a target, such as a join of shorter segments; \
         near-variant decoys that share most fragments with their target, for --dia on \
         wide-window data such as diaPASEF",
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

/// Highest storage budget of the native Percolator engine when --rescore-memory-gb is not
/// given, in GB: the engine takes what the PIN needs, up to this and the memory free.
pub(crate) const DEFAULT_RESCORE_MEMORY_GB: u64 = 16;

pub(crate) type Cli = SearchArgs;

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

#[cfg(test)]
#[path = "cli/tests/validators.rs"]
mod cli_domain_validator_tests;

#[cfg(test)]
#[path = "cli/tests/enzymes.rs"]
mod enzyme_cli_tests;

#[cfg(test)]
#[path = "cli/tests/fdr.rs"]
mod fdr_parameter_tests;

#[cfg(test)]
#[path = "cli/tests/dirs.rs"]
mod dir_or_output_dir_tests;

#[cfg(test)]
#[path = "cli/tests/spelling.rs"]
mod spelling_tests;

#[cfg(test)]
#[path = "cli/tests/help.rs"]
mod help_layout_tests;
