//! Shared test fixtures for the search crate's integration tests.
//!
//! Used via `mod common; use common::*;` in each integration test file.
//! Cargo treats `tests/common/mod.rs` as a non-test module per
//! https://doc.rust-lang.org/cargo/guide/tests.html#integration-tests.

#![allow(dead_code)] // some helpers are used by only a subset of tests

use std::path::PathBuf;

use model::{AminoAcidSetBuilder, ModLocation, Modification, ResidueSpec};
use scoring_crate::{Param, RankScorer};

/// Resolve a path relative to the workspace root (CARGO_MANIFEST_DIR/../../..).
///
/// Pass the full path from the repo root, e.g.
/// `fixture("src/test/resources/BSA.fasta")`.
pub fn fixture(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .join(rel)
        .canonicalize()
        .unwrap_or_else(|e| panic!("canonicalize {rel}: {e}"))
}

/// Standard BSA-search aa_set: Carbamidomethyl-C fixed + Oxidation-M variable.
pub fn aa_set() -> model::AminoAcidSet {
    let cam = Modification {
        name: "Carbamidomethyl".into(),
        mass_delta: 57.02146,
        residue: ResidueSpec::Specific(b'C'),
        location: ModLocation::Anywhere,
        fixed: true,
        accession: None,
    };
    let ox = Modification {
        name: "Oxidation".into(),
        mass_delta: 15.99491,
        residue: ResidueSpec::Specific(b'M'),
        location: ModLocation::Anywhere,
        fixed: false,
        accession: None,
    };
    AminoAcidSetBuilder::new_standard()
        .add_fixed_mod(cam)
        .add_variable_mod(ox)
        .build()
        .unwrap()
}

/// Load the bundled `HCD_QExactive_Tryp.param` and construct a RankScorer.
pub fn rank_scorer() -> RankScorer {
    let param_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .join("src/main/resources/ionstat/HCD_QExactive_Tryp.param")
        .canonicalize()
        .unwrap_or_else(|e| panic!("canonicalize HCD_QExactive_Tryp.param: {e}"));
    let param = Param::load_from_file(&param_path)
        .unwrap_or_else(|e| panic!("load HCD_QExactive_Tryp.param: {e}"));
    RankScorer::new(&param)
}
