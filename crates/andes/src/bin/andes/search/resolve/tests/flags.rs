//! Tests of incompatible flag combinations.

use super::checks::check_flag_combinations;

use super::*;
use clap::Parser;

#[derive(Parser, Debug)]
struct Wrap {
    #[command(flatten)]
    args: Cli,
}

fn args(extra: &[&str]) -> Cli {
    let argv = std::iter::once("andes").chain(extra.iter().copied());
    Wrap::try_parse_from(argv).unwrap().args
}

fn accepted(extra: &[&str]) {
    if let Err(e) = check_flag_combinations(&args(extra)) {
        panic!("{extra:?} was rejected: {e}");
    }
}

fn rejected(extra: &[&str], fragment: &str) {
    match check_flag_combinations(&args(extra)) {
        Ok(()) => panic!("{extra:?} was accepted"),
        Err(e) => assert!(e.contains(fragment), "{extra:?}: {e}"),
    }
}

#[test]
fn typical_runs_are_accepted() {
    accepted(&[]);
    accepted(&["--rescore"]);
    accepted(&["--rescore", "--localize", "ipf"]);
    accepted(&["--rescore", "--chimeric", "--refine"]);
    accepted(&["--dia"]);
    accepted(&["--open-search"]);
    accepted(&["--lfq", "--mbr", "--rescore"]);
    accepted(&["--glyco", "--glyco-index-sequon-only"]);
    accepted(&["--glyco", "--lfq", "--mbr"]);
    accepted(&["--xl", "DSSO"]);
}

#[test]
fn top_n_must_be_positive() {
    rejected(&["--top-n", "0"], "--top-n must be at least 1");
}

#[test]
fn localize_needs_rescoring() {
    rejected(&["--localize", "ipf"], "--localize ipf needs --rescore");
}

#[test]
fn rescore_alone_selects_the_native_engine() {
    let cli = args(&["--rescore"]);
    assert!(cli.rescore && !cli.external_percolator());
}

#[test]
fn a_percolator_binary_or_docker_selects_the_external_percolator() {
    let cli = args(&["--rescore", "--percolator-bin", "/opt/percolator"]);
    assert!(cli.rescore && cli.external_percolator());
    accepted(&["--rescore", "--percolator-bin", "/opt/percolator"]);
    let cli = args(&["--rescore", "--percolator-docker"]);
    assert!(cli.rescore && cli.external_percolator());
    accepted(&[
        "--rescore",
        "--percolator-docker",
        "--percolator-image",
        "other:1",
    ]);
    accepted(&[
        "--rescore",
        "--percolator-bin",
        "/opt/percolator",
        "--percolator-args=--maxiter 3",
    ]);
}

#[test]
fn percolator_options_need_an_external_percolator() {
    rejected(
        &["--rescore", "--percolator-args=--maxiter 3"],
        "native rescoring takes no percolator_args",
    );
    rejected(
        &["--percolator-args=--maxiter 3"],
        "native rescoring takes no percolator_args",
    );
    rejected(
        &["--rescore", "--percolator-image", "other:1"],
        "--percolator-image has no effect without --percolator-docker",
    );
    accepted(&["--xl", "DSSO", "--percolator-args=--maxiter 3"]);
}

#[test]
fn rescore_native_is_an_unknown_flag() {
    let argv = ["andes", "--rescore-native"];
    let err = Wrap::try_parse_from(argv).unwrap_err();
    assert_eq!(err.kind(), clap::error::ErrorKind::UnknownArgument);
}

#[test]
fn dia_is_not_combined_with_other_modes() {
    rejected(&["--dia", "--chimeric"], "--dia and --chimeric");
    rejected(&["--dia", "--open-search"], "--dia and --open-search");
    rejected(&["--dia", "--refine"], "--dia and --refine");
    rejected(&["--dia", "--glyco"], "--dia and --glyco");
    rejected(&["--dia", "--xl", "DSSO"], "--dia and --xl");
}

#[test]
fn sequon_only_needs_glyco() {
    rejected(
        &["--glyco-index-sequon-only"],
        "--glyco-index-sequon-only requires --glyco",
    );
}

#[test]
fn open_search_is_not_combined_with_other_modes() {
    rejected(&["--open-search", "--glyco"], "not supported with --glyco");
    rejected(&["--open-search", "--chimeric"], "with --chimeric");
    rejected(&["--open-search", "--label", "x"], "with --label");
    rejected(&["--open-search", "--tmt", "tmt10plex"], "with --tmt");
    rejected(&["--open-search", "--lfq"], "with --lfq or --mbr");
    rejected(&["--open-search", "--mbr"], "with --lfq or --mbr");
    rejected(&["--open-search", "--refine"], "with --refine");
    rejected(&["--open-search", "--xl", "DSSO"], "with --xl");
    rejected(
        &["--open-search", "--isotope-error", "0..2"],
        "remove --isotope-error",
    );
}

#[test]
fn xl_rejects_the_flags_it_does_not_support() {
    rejected(
        &[
            "--xl",
            "DSSO",
            "--output-tsv",
            "out.tsv",
            "--lfq",
            "--refine",
        ],
        "--xl does not support: --output-tsv, --refine, --lfq",
    );
    rejected(
        &["--xl", "DSSO", "--tmt", "tmt10plex"],
        "--xl does not support: --tmt",
    );
}

#[test]
fn glyco_rejects_the_flags_it_does_not_support() {
    rejected(
        &[
            "--glyco",
            "--output-parquet",
            "out",
            "--rescore",
            "--pep",
            "0.01",
        ],
        "--glyco does not support: --output-parquet, --rescore, --pep",
    );
    rejected(
        &["--glyco", "--refine"],
        "--glyco does not support: --refine",
    );
}

#[test]
fn glyco_match_between_runs_needs_lfq() {
    rejected(&["--glyco", "--mbr"], "--glyco --mbr needs --lfq");
}

#[test]
fn optimize_is_not_supported_with_glyco() {
    rejected(
        &["--glyco", "--optimize"],
        "--optimize is not supported with --glyco",
    );
}
