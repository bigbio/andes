//! Tests of the cross-link quantification tables.

use quant::lfq::{FeatureQuant, LfqParams};

use super::super::library::{Precursor, Side};
use super::super::Evidence;
use super::*;

fn settings() -> Settings {
    Settings {
        params: LfqParams::default(),
        rt_window_s: None,
        xl_fdr: 0.01,
        transfer_fdr: 0.05,
        mbr: true,
    }
}

fn side(peptide: &str, site: usize) -> Side {
    Side {
        peptide: peptide.into(),
        bare: peptide.into(),
        site,
        site_kind: "K".into(),
    }
}

fn precursor(charge: u8, residue_pair: &str) -> Precursor {
    Precursor {
        key: format!("DSS|AKA@2|GKG@2/{charge}"),
        linker: "DSS".into(),
        sides: [side("AKA", 2), side("GKG", 2)],
        charge,
        mono_mz: 500.0,
        envelope: vec![0.4, 0.3, 0.2, 0.1],
        isotope_share: 0.95,
        residue_pair: residue_pair.into(),
        protein_pair: "P -- P".into(),
        link: "intra",
        q: 0.004,
    }
}

fn row(precursor: usize, run: usize, evidence: Evidence, area: f64, q: Option<f64>) -> Row {
    Row {
        precursor,
        run,
        evidence,
        anchor_rt: 100.0,
        anchor_sd_s: 5.0,
        window_s: 40.0,
        scan: 17,
        csm_score: 12.5,
        n_csms: 2,
        feature: Some(FeatureQuant {
            apex_rt: 101.0,
            rt_start: 95.0,
            rt_stop: 108.0,
            area,
            apex_intensity: 1e5,
            isotope_areas: vec![4.0, 3.0, 2.0, 1.0],
            cosine: 0.99,
            n_isotopes_found: 4,
            rt_delta_s: 1.0,
            n_scans: 30,
        }),
        decoy: None,
        score: Some(0.9),
        q,
    }
}

#[test]
fn a_residue_pair_is_the_sum_of_its_accepted_precursors_per_run() {
    let lib = Library {
        precursors: vec![
            precursor(3, "P:2 -- P:9"),
            precursor(4, "P:2 -- P:9"),
            precursor(3, "P:5 -- P:9"),
        ],
        identifications: Vec::new(),
    };
    let rows = vec![
        row(0, 0, Evidence::Identified, 100.0, Some(0.0)),
        row(1, 0, Evidence::Identified, 50.0, Some(0.0)),
        row(0, 1, Evidence::Transferred, 30.0, Some(0.01)),
        row(1, 1, Evidence::Transferred, 999.0, Some(0.5)),
        row(2, 1, Evidence::Transferred, 999.0, None),
    ];
    let dir = tempfile::tempdir().expect("tempdir");
    let path = dir.path().join("out.xl.residue_pairs.quant.tsv");
    let n = write_residue_pairs(&path, &lib, &rows, &["a", "b"], &settings()).expect("written");
    assert_eq!(n, 2);
    let text = std::fs::read_to_string(&path).expect("read");
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(
        lines[0],
        "residue_pair\tprotein_pair\tlink\tq_value\tn_precursors\ta\tb"
    );
    assert_eq!(
        lines[1],
        "P:2 -- P:9\tP -- P\tintra\t0.004000\t2\t150.0\t30.0"
    );
    assert_eq!(lines[2], "P:5 -- P:9\tP -- P\tintra\t0.004000\t0\t\t");
}

#[test]
fn a_feature_row_carries_the_precursor_and_the_measurement() {
    let lib = Library {
        precursors: vec![precursor(3, "P:2 -- P:9")],
        identifications: Vec::new(),
    };
    let s = settings();
    let accepted = feature_row(
        &row(0, 0, Evidence::Transferred, 30.0, Some(0.01)),
        &lib,
        "a",
        &s,
    );
    assert_eq!(accepted.run, "a");
    assert_eq!(accepted.precursor, "DSS|AKA@2|GKG@2/3");
    assert_eq!(accepted.evidence, "transferred");
    assert!(accepted.accepted);
    assert_eq!(accepted.quantity, Some(30.0));
    assert_eq!(
        (accepted.alpha_peptide.as_str(), accepted.alpha_site),
        ("AKA", 2)
    );
    assert_eq!((accepted.charge, accepted.n_isotopes), (3, 4));
    assert_eq!(accepted.anchor_sd_s, 5.0);
    let rejected = feature_row(
        &row(0, 0, Evidence::Transferred, 30.0, Some(0.5)),
        &lib,
        "a",
        &s,
    );
    assert!(!rejected.accepted);
    assert_eq!(
        rejected.quantity,
        Some(30.0),
        "the value stays in the table"
    );
}

#[test]
fn the_tables_are_named_after_the_report_without_a_pin_extension() {
    assert_eq!(
        sibling(Path::new("out/run.pin"), ".xl.features.parquet"),
        PathBuf::from("out/run.xl.features.parquet")
    );
    assert_eq!(
        sibling(Path::new("out/run.v2"), ".xl.features.parquet"),
        PathBuf::from("out/run.v2.xl.features.parquet")
    );
}
