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

/// Strip Percolator flanking (`X.PEPTIDE.Y`) and mod-mass tokens like
/// `+57.021` / `-18.0` from a `.pin`-format peptide string. Returns the
/// residue-only sequence in uppercase.
///
/// Implementation note: a naive `split('.').nth(1)` is WRONG for any peptide
/// containing a mod-mass (e.g. `K.GAC+57.021LLPK.E` → buggy parser yields
/// `"GAC+57"` → `"GAC"`). The flanking dots are at fixed byte positions
/// (1 and len-2) when the flanking residue is a single character (always
/// the case in `.pin` output). Mod-mass dots lie strictly inside that
/// middle range. We extract the middle and strip mod-mass tokens
/// (`[+-]\d+(\.\d+)?`) explicitly.
pub fn strip_flanking_and_mods(pin_pep: &str) -> String {
    let bytes = pin_pep.as_bytes();
    if bytes.len() < 5 {
        return String::new();
    }
    if bytes[1] != b'.' || bytes[bytes.len() - 2] != b'.' {
        return String::new();
    }
    let middle = &pin_pep[2..pin_pep.len() - 2];
    let mut out = String::with_capacity(middle.len());
    let mut chars = middle.chars().peekable();
    while let Some(c) = chars.next() {
        if c == '+' || c == '-' {
            // Consume mod-mass tail: digits, optional dot, optional digits.
            while let Some(&nc) = chars.peek() {
                if nc.is_ascii_digit() || nc == '.' {
                    chars.next();
                } else {
                    break;
                }
            }
        } else if c.is_ascii_uppercase() {
            out.push(c);
        }
    }
    out
}

#[cfg(test)]
mod parser_tests {
    use super::strip_flanking_and_mods;

    #[test]
    fn strips_flanking_only() {
        assert_eq!(strip_flanking_and_mods("R.PEPTIDE.K"), "PEPTIDE");
    }

    #[test]
    fn strips_one_mod_mass() {
        assert_eq!(strip_flanking_and_mods("K.PEPTM+15.995DE.R"), "PEPTMDE");
    }

    #[test]
    fn strips_multiple_mod_masses() {
        // Regression: the case that broke the prior naive parser.
        assert_eq!(
            strip_flanking_and_mods("K.GAC+57.021LLPKIETM+15.995R.E"),
            "GACLLPKIETMR"
        );
    }

    #[test]
    fn strips_negative_mod_mass() {
        assert_eq!(strip_flanking_and_mods("K.PEPM-18.0R.E"), "PEPMR");
    }

    #[test]
    fn handles_protein_terminal_dash_flanking() {
        assert_eq!(strip_flanking_and_mods("-.PEPTIDE.R"), "PEPTIDE");
        assert_eq!(strip_flanking_and_mods("R.PEPTIDE.-"), "PEPTIDE");
    }
}
