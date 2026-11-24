//! Tests of training from MSnet PSM tables.

use super::inputs::run_stem;
use super::msnet::build_msnet_peptide;
use std::path::Path;

use model::mass::H2O;

#[test]
fn run_stem_strips_spectrum_extensions() {
    assert_eq!(run_stem(Path::new("/x/a.mzML.gz")), "a");
    assert_eq!(run_stem(Path::new("/x/PXD1_run.d")), "PXD1_run");
    assert_eq!(run_stem(Path::new("/x/PXD1_run.d/")), "PXD1_run");
    assert_eq!(run_stem(Path::new("b.MGF")), "b");
    assert_eq!(run_stem(Path::new("c.raw")), "c");
    assert_eq!(run_stem(Path::new(".d")), ".d");
}

fn unmod_mass(seq: &[u8]) -> f64 {
    let rsum: f64 = seq
        .iter()
        .map(|&r| model::AminoAcid::standard(r).unwrap().mass)
        .sum();
    rsum + H2O
}

#[test]
fn unmodified_peptide_mass_matches_reference() {
    let p = build_msnet_peptide("PEPTIDEK", &[], &[], 0.0, 0.0).unwrap();
    let expected = unmod_mass(b"PEPTIDEK");
    assert_eq!(p.mass().to_bits(), expected.to_bits());
    assert!(p.residues.iter().all(|aa| !aa.is_modified()));
}

#[test]
fn residue_mod_mass_is_correct() {
    const OX: f64 = 15.994915;
    let p = build_msnet_peptide("MPEPTIDE", &[1], &[OX], 0.0, 0.0).unwrap();
    let expected = unmod_mass(b"MPEPTIDE") + OX;
    assert_eq!(p.mass().to_bits(), expected.to_bits());
    assert!(
        p.residues[0].is_modified(),
        "residue 1 should carry the mod"
    );
    assert_eq!(p.residues[0].mod_.as_ref().unwrap().mass_delta, OX);
}

#[test]
fn nterm_mod_mass_is_correct() {
    const ACETYL: f64 = 42.010565;
    let p = build_msnet_peptide("PEPTIDEK", &[], &[], ACETYL, 0.0).unwrap();
    let expected = unmod_mass(b"PEPTIDEK") + ACETYL;
    assert_eq!(p.mass().to_bits(), expected.to_bits());
    assert!(p.residues[0].is_modified());
}

#[test]
fn overlapping_residue_and_nterm_mods_sum() {
    const OX: f64 = 15.994915;
    const ACETYL: f64 = 42.010565;
    let p = build_msnet_peptide("MPEPTIDE", &[1], &[OX], ACETYL, 0.0).unwrap();
    let expected = unmod_mass(b"MPEPTIDE") + OX + ACETYL;
    assert_eq!(p.mass().to_bits(), expected.to_bits());
    assert_eq!(p.residues[0].mod_.as_ref().unwrap().mass_delta, OX + ACETYL);
}

#[test]
fn mismatched_mod_arrays_error() {
    let r = build_msnet_peptide("PEPTIDE", &[1, 2], &[1.0], 0.0, 0.0);
    assert!(r.is_err());
}

#[test]
fn nonstandard_residue_errors() {
    let r = build_msnet_peptide("PEPTBDE", &[], &[], 0.0, 0.0);
    assert!(r.is_err());
}
