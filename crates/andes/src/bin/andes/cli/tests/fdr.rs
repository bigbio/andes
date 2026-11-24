//! Tests of the FDR flags and their defaults.

use super::commands::{QuantFinalizeArgs, QuantLibraryArgs};
use super::options::SearchArgs;

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
