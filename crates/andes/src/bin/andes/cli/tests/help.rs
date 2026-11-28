//! Tests of the help layout: headings, order and advanced flags.

use super::commands::{InferenceArgs, QuantFinalizeArgs};
use super::options::SearchArgs;

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
        .map(|l| l.trim_start().split([' ', '=']).next().unwrap().to_string())
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
    "no_search_gbdt",
    "no_search_rich_ion",
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
    "dia_interference",
];

#[test]
fn experimental_entrapment_flags_are_advanced_diagnostics() {
    let mut cmd = SearchWrap::command();
    let short = cmd.render_help().to_string();
    let long = cmd.render_long_help().to_string();
    for flag in ["--add-entrapment", "--entrapment-seed"] {
        assert!(!listed(&short).iter().any(|listed| listed == flag));
        assert!(listed(&long).iter().any(|listed| listed == flag));
    }
    assert!(long.contains("experimental"));
    assert!(long.contains("have not been validated"));
    for id in ["add_entrapment", "entrapment_seed"] {
        let arg = cmd.get_arguments().find(|arg| arg.get_id() == id).unwrap();
        assert_eq!(arg.get_help_heading(), Some(ADVANCED_FDR));
    }
}

#[test]
fn short_help_lists_the_main_options_under_five_headings() {
    let mut cmd = SearchWrap::command();
    let short = cmd.render_help().to_string();
    assert_eq!(listed(&short), MAIN, "{short}");
    for heading in [H_IO, H_SEARCH, H_FDR, H_MODES, H_QUANT] {
        assert!(
            short.contains(&format!("{heading}:")),
            "-h misses {heading}"
        );
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
        assert!(
            !short.contains(&format!("{heading}:")),
            "-h shows {heading}"
        );
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
        for bad in [
            "**",
            "#64",
            "bigbio/",
            "@1%",
            "action = Set",
            "Every other engine",
        ] {
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
