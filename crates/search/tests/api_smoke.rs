//! Smoke test exercising the re-exported public API end-to-end. If this
//! compiles and passes, downstream crates can import the same types
//! without touching submodule paths.

use model::{
    AminoAcid, AminoAcidSetBuilder, Enzyme, ModLocation, Modification,
    Peptide, PrecursorTolerance, ResidueSpec, Tolerance, H2O, PROTON,
};

#[test]
fn build_set_and_peptide_via_public_api() {
    let cam = Modification {
        name: "Carbamidomethyl".to_string(),
        mass_delta: 57.02146,
        residue: ResidueSpec::Specific(b'C'),
        location: ModLocation::Anywhere,
        fixed: true,
        accession: None,
        neutral_losses: Vec::new(),
    };
    let set = AminoAcidSetBuilder::new_standard()
        .add_fixed_mod(cam)
        .build()
        .unwrap();

    let residues: Vec<AminoAcid> = b"PEPTIDE".iter()
        .map(|&r| AminoAcid::standard(r).unwrap())
        .collect();
    let p = Peptide::new(residues, b'_', b'-').with_charge(2);

    assert_eq!(p.length(), 7);
    assert_eq!(p.charge, Some(2));
    assert_eq!(p.to_string(), "_.PEPTIDE.-");

    let p2 = Peptide::from_str("_.PEPTIDE.-", &set).unwrap();
    assert_eq!(p2.to_string(), p.to_string());
}

#[test]
fn enzyme_and_tolerance_via_public_api() {
    assert!(Enzyme::Trypsin.is_cleavable_after(b'K'));
    let t = Tolerance::Ppm(10.0);
    assert_eq!(t.as_da(1000.0), 0.01);
    let pt = PrecursorTolerance::symmetric(t);
    assert_eq!(pt.left.as_da(1000.0), pt.right.as_da(1000.0));
}

#[test]
fn chemistry_constants_via_public_api() {
    // Just confirm they're reachable through the re-export and have
    // sensible (non-zero) values; bit-exact pinning lives in the
    // chemistry parity test.
    assert_eq!(PROTON, 1.00727649);
    // Compare via a runtime binding to avoid clippy::assertions_on_constants.
    let h2o: f64 = H2O;
    assert!(h2o > 18.0 && h2o < 18.1);
}
