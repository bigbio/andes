//! Integration tests for [`model_train::labeled::bootstrap_labels`].
//!
//! Uses the BSA fixture (test-fixtures/test.mgf + test-fixtures/BSA.fasta)
//! and the HCD_QExactive_Tryp.param seed model.

use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use input::MgfReader;
use model::{AminoAcidSetBuilder, ModLocation, Modification, ResidueSpec};
use scoring_crate::{Param, RankScorer};
use search::SearchParams;

use model_train::labeled::bootstrap_labels;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Resolve path relative to the workspace root (CARGO_MANIFEST_DIR/../../..).
fn fixture(rel: &str) -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(rel)
        .canonicalize()
        .unwrap_or_else(|e| panic!("canonicalize {rel}: {e}"))
}

/// Standard BSA-search AminoAcidSet: Carbamidomethyl-C fixed + Oxidation-M variable.
fn bsa_aa_set() -> model::AminoAcidSet {
    let cam = Modification {
        name: "Carbamidomethyl".into(),
        mass_delta: 57.02146,
        residue: ResidueSpec::Specific(b'C'),
        location: ModLocation::Anywhere,
        fixed: true,
        accession: None,
        neutral_losses: Vec::new(),
        loss_class: 0,
    };
    let ox = Modification {
        name: "Oxidation".into(),
        mass_delta: 15.99491,
        residue: ResidueSpec::Specific(b'M'),
        location: ModLocation::Anywhere,
        fixed: false,
        accession: None,
        neutral_losses: Vec::new(),
        loss_class: 0,
    };
    AminoAcidSetBuilder::new_standard()
        .add_fixed_mod(cam)
        .add_variable_mod(ox)
        .build()
        .unwrap()
}

/// Load HCD_QExactive_Tryp.param from the model-train fixtures directory.
fn load_hcd_scorer() -> RankScorer {
    let param_path = Path::new(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/HCD_QExactive_Tryp.param"
    ));
    let param = Param::load_from_file(param_path).expect("load HCD_QExactive_Tryp.param");
    RankScorer::new(&param)
}

/// Load spectra from BSA test.mgf fixture.
fn load_bsa_spectra() -> Vec<model::Spectrum> {
    let path = fixture("test-fixtures/test.mgf");
    let file = File::open(&path).expect("open test.mgf");
    MgfReader::new(BufReader::new(file))
        .filter_map(|r| r.ok())
        .collect()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

/// Basic smoke test: bootstrap_labels runs without error, returns a Vec,
/// is deterministic, and every returned match has confidence <= train_fdr.
///
/// We use a lenient FDR of 0.5 to maximise the chance of a non-empty result
/// on the small BSA fixture.
#[test]
fn bootstrap_labels_bsa_smoke() {
    let spectra = load_bsa_spectra();
    assert!(!spectra.is_empty(), "test.mgf must contain at least one spectrum");

    let scorer = load_hcd_scorer();
    let db_path = fixture("test-fixtures/BSA.fasta");
    let params = SearchParams::default_tryptic(bsa_aa_set());
    let train_fdr = 0.5_f64;

    let result = bootstrap_labels(&spectra, &db_path, &scorer, &params, train_fdr);
    assert!(result.is_ok(), "bootstrap_labels must succeed: {:?}", result.err());

    let matches = result.unwrap();

    // Every returned match must satisfy its claimed confidence threshold.
    for m in &matches {
        assert!(
            m.confidence <= train_fdr,
            "match confidence {} > train_fdr {}",
            m.confidence,
            train_fdr
        );
        // spectrum_index must be in range.
        assert!(
            m.spectrum_index < spectra.len(),
            "spectrum_index {} out of range (spectra.len()={})",
            m.spectrum_index,
            spectra.len()
        );
    }

    // Determinism: running twice must produce the same count.
    let result2 = bootstrap_labels(&spectra, &db_path, &scorer, &params, train_fdr)
        .expect("second bootstrap_labels run must succeed");
    assert_eq!(
        matches.len(),
        result2.len(),
        "bootstrap_labels must be deterministic: first run returned {}, second returned {}",
        matches.len(),
        result2.len()
    );

    // Log the count for informational purposes (not a hard failure if 0).
    eprintln!(
        "bootstrap_labels BSA @ train_fdr={}: {} confident target PSMs",
        train_fdr,
        matches.len()
    );
}

/// Strict test with train_fdr = 0.05 — fewer matches expected.
/// Every confidence value must still be <= 0.05.
#[test]
fn bootstrap_labels_bsa_strict_fdr() {
    let spectra = load_bsa_spectra();
    let scorer = load_hcd_scorer();
    let db_path = fixture("test-fixtures/BSA.fasta");
    let params = SearchParams::default_tryptic(bsa_aa_set());
    let train_fdr = 0.05_f64;

    let matches = bootstrap_labels(&spectra, &db_path, &scorer, &params, train_fdr)
        .expect("bootstrap_labels must succeed at strict FDR");

    for m in &matches {
        assert!(
            m.confidence <= train_fdr,
            "match confidence {} > train_fdr {}",
            m.confidence,
            train_fdr
        );
    }

    eprintln!(
        "bootstrap_labels BSA @ train_fdr={}: {} confident target PSMs",
        train_fdr,
        matches.len()
    );
}

/// FDR monotonicity: lenient FDR must produce at least as many matches as strict.
#[test]
fn bootstrap_labels_bsa_lenient_ge_strict() {
    let spectra = load_bsa_spectra();
    let scorer = load_hcd_scorer();
    let db_path = fixture("test-fixtures/BSA.fasta");
    let params = SearchParams::default_tryptic(bsa_aa_set());

    let lenient = bootstrap_labels(&spectra, &db_path, &scorer, &params, 0.5)
        .expect("lenient bootstrap");
    let strict = bootstrap_labels(&spectra, &db_path, &scorer, &params, 0.05)
        .expect("strict bootstrap");

    assert!(
        lenient.len() >= strict.len(),
        "lenient FDR=0.5 should give >= matches as strict FDR=0.05: \
         lenient={}, strict={}",
        lenient.len(),
        strict.len()
    );
}
