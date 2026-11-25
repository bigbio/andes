//! Tests of the run report built from the parameter sources and a search's results.

use clap::{CommandFactory, FromArgMatches};

use super::build::{run_report, Results, Trace};
use super::{InputHasher, StageClock};
use crate::params::sources::from_command_line;

fn matches(extra: &[&str]) -> clap::ArgMatches {
    let base = ["andes", "--spectrum", "x.mzML", "--database", "y.fasta", "--output-pin", "o.pin"];
    crate::TopCli::command()
        .try_get_matches_from(base.iter().chain(extra))
        .expect("parse args")
}

#[test]
fn a_flag_set_to_its_default_is_reported_as_command_line() {
    let sources = from_command_line(&matches(&["--precursor-tol", "20ppm"]));
    let rows = crate::params::render::report_rows(&sources);
    let row = rows.iter().find(|r| r.id == "precursor_tol").unwrap();
    assert_eq!((row.value.as_str(), row.source.as_str()), ("20ppm", "command line"));
    let missed = rows.iter().find(|r| r.id == "max_missed_cleavages").unwrap();
    assert_eq!(missed.source, "default");
    assert!(missed.number.is_some(), "{missed:?}");
    assert!(rows.iter().all(|r| r.id != "print_config"));
}

#[test]
fn a_run_without_rescoring_has_no_fdr_counts_and_no_percolator_sentence() {
    let m = matches(&[]);
    let cli = crate::SearchArgs::from_arg_matches(&m).unwrap();
    let aa = model::amino_acid::set::AminoAcidSetBuilder::new_standard().build().unwrap();
    let params = search::SearchParams::default_tryptic(aa);
    let mut clock = StageClock::new();
    clock.stage("search");
    let trace = Trace {
        clock,
        hasher: InputHasher::start(vec![("spectrum".to_string(), "x.mzML".into())]),
        sources: from_command_line(&m),
        command_line: "andes --spectrum x.mzML".to_string(),
        models: vec![("scoring".to_string(), "hcd_highres_tryp".to_string())],
        decoy_strategy: search::decoy::DecoyStrategy::Shuffle,
        decoy_key: None,
    };
    let report = run_report(
        trace,
        &Results {
            cli: &cli,
            params: &params,
            fragment_tolerance: model::tolerance::Tolerance::Ppm(20.0),
            spectra: 10,
            rank1_target_psms: 7,
            rescored: None,
            protein_groups: None,
        },
        "o.statistics.parquet",
    );
    assert!(report.counts.iter().all(|c| c.fdr.is_none()), "{:?}", report.counts);
    assert!(!report.methods.contains("Percolator"), "{}", report.methods);
    assert!(report.methods.contains("precursor tolerance of 20 ppm"), "{}", report.methods);
    let stages: Vec<&str> = report.stages.iter().map(|s| s.name.as_str()).collect();
    assert_eq!(stages, ["search", "hash inputs"]);
    assert!(report.inputs[0].sha256.is_err());
    assert_eq!(report.references.len(), 2);
    let text = report.render();
    assert!(text.contains("──────── Methods ────────"), "{text}");
}
