//! Tests of the cross-link precursor library: envelopes, formulas and keys.

use model::amino_acid::AminoAcid;
use xl::crosslinker::Crosslinker;

use super::*;

fn peptide(sequence: &[u8]) -> Peptide {
    let residues = sequence
        .iter()
        .map(|&r| AminoAcid::standard(r).expect("standard residue"))
        .collect();
    Peptide::new(residues, b'K', b'A')
}

fn linker(spec: &str) -> Crosslinker {
    xl::spec::XlDecl::parse(spec)
        .and_then(|d| d.resolve(model::unimod::bundled()))
        .expect("bundled reagent")
}

#[test]
fn the_traced_envelope_holds_the_coverage_and_sums_to_one() {
    let full = [0.1, 0.2, 0.3, 0.2, 0.1, 0.05, 0.03, 0.02];
    let (env, share) = traced_envelope(&full);
    assert_eq!(env.len(), 6);
    assert!((env.iter().sum::<f64>() - 1.0).abs() < 1e-12);
    assert!((share - 0.95).abs() < 1e-12);
    assert!((env[2] - 0.3 / 0.95).abs() < 1e-12);
}

#[test]
fn the_traced_envelope_keeps_between_the_fewest_and_the_most_isotopes() {
    let (short, share) = traced_envelope(&[0.9, 0.08, 0.02]);
    assert_eq!(short.len(), 3);
    assert!((share - 1.0).abs() < 1e-12);
    let narrow = [0.97, 0.01, 0.01, 0.005, 0.005];
    assert_eq!(traced_envelope(&narrow).0.len(), MIN_ISOTOPES);
    let flat = [1.0 / 30.0; 30];
    assert_eq!(traced_envelope(&flat).0.len(), MAX_ISOTOPES);
    assert_eq!(traced_envelope(&[]), (vec![1.0], 1.0));
    assert_eq!(traced_envelope(&[0.0, 0.0]), (vec![1.0], 1.0));
}

#[test]
fn the_pair_formula_is_both_peptides_and_the_link() {
    let (a, b) = (peptide(b"AAK"), peptide(b"GK"));
    let f = pair_formula(&a, &b, &linker("DSS"));
    let mut expected = Formula::peptide(b"AAK".iter().map(|&r| (r, 0.0)));
    expected.add(Formula::peptide(b"GK".iter().map(|&r| (r, 0.0))));
    assert_eq!(f.c, expected.c + 8);
    assert_eq!(f.h, expected.h + 10);
    assert_eq!(f.o, expected.o + 2);
    assert_eq!((f.n, f.s), (expected.n, expected.s));
    let neutral = a.mass() + b.mass() + linker("DSS").intact;
    assert!(
        (f.monoisotopic_mass() - neutral).abs() < 1e-6,
        "{} vs {neutral}",
        f.monoisotopic_mass()
    );
}

#[test]
fn a_zero_length_link_removes_water() {
    let (a, b) = (peptide(b"AAK"), peptide(b"GDK"));
    let f = pair_formula(&a, &b, &linker("EDC"));
    let neutral = a.mass() + b.mass() + linker("EDC").intact;
    assert!((f.monoisotopic_mass() - neutral).abs() < 1e-6);
}

#[test]
fn a_cross_linked_pair_has_a_wider_envelope_than_one_of_its_peptides() {
    let a = peptide(b"LVNEVTEFAKTCVADESAENCDK");
    let b = peptide(b"ALVLIAFAQYLQQCPFEDHVK");
    let (pair, _) =
        traced_envelope(&pair_formula(&a, &b, &linker("DSS")).envelope(ENVELOPE_ISOTOPES));
    let (single, _) = traced_envelope(&peptide_formula(&a).envelope(ENVELOPE_ISOTOPES));
    assert!(
        pair.len() > single.len(),
        "{} vs {}",
        pair.len(),
        single.len()
    );
    assert!(pair[0] < single[0]);
    let top = pair.iter().copied().fold(0.0, f64::max);
    assert!(
        pair[0] < 0.5 * top,
        "the monoisotopic peak is weak: {pair:?}"
    );
}

#[test]
fn the_peptidoform_has_no_flanking_residues() {
    assert_eq!(peptidoform(&peptide(b"AAK")), "AAK");
}

#[test]
fn a_side_key_marks_a_terminal_link() {
    let side = |kind: &str| Side {
        peptide: "AAK".into(),
        bare: "AAK".into(),
        site: 1,
        site_kind: kind.into(),
    };
    assert_eq!(side_key(&side("A")), "AAK@1");
    assert_eq!(side_key(&side("nterm")), "AAK@1n");
    assert_eq!(side_key(&side("cterm")), "AAK@1c");
}
