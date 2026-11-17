use model::{
    AminoAcid, AminoAcidSetBuilder, Enzyme, ModLocation, Modification, Peptide, PrecursorTolerance,
    ResidueSpec, Tolerance, H2O, PROTON,
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
        loss_class: 0,
        annotated_site_mass: None,
        blocks_cleavage: false,
        parts: Vec::new(),
        diagnostic_ions: Vec::new(),
    };
    let set = AminoAcidSetBuilder::new_standard()
        .add_fixed_mod(cam)
        .build()
        .unwrap();

    let residues: Vec<AminoAcid> = b"PEPTIDE"
        .iter()
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
    assert_eq!(PROTON, 1.00727649);
    let h2o: f64 = H2O;
    assert!(h2o > 18.0 && h2o < 18.1);
}
