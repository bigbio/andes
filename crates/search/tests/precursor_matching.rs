//! Precursor-mass tolerance tests.

use model::{AminoAcid, Peptide, PrecursorTolerance, Spectrum, Tolerance, PROTON};
use search::{matches_precursor};

fn make_peptide(seq: &[u8]) -> Peptide {
    let residues: Vec<AminoAcid> = seq.iter().map(|&r| AminoAcid::standard(r).unwrap()).collect();
    Peptide::new(residues, b'_', b'-')
}

fn make_spectrum(precursor_mz: f64, charge: Option<i32>) -> Spectrum {
    Spectrum {
        title: "test".into(),
        precursor_mz,
        precursor_intensity: None,
        precursor_charge: charge,
        rt_seconds: None,
        scan: None,
        peaks: vec![],
        activation_method: None,
        isolation_lower_offset: None,
        isolation_upper_offset: None,
    }
}

#[test]
fn exact_mass_match() {
    let peptide = make_peptide(b"AR");
    let mass = peptide.mass();
    let charge = 2u8;
    let mz = (mass + charge as f64 * PROTON) / charge as f64;
    let spec = make_spectrum(mz, Some(charge as i32));
    let tol = PrecursorTolerance::symmetric(Tolerance::Ppm(20.0));
    let err = matches_precursor(&spec, &peptide, charge, 0, &tol, 0.0).expect("should match");
    assert!(err.mass_error_ppm.abs() < 0.001, "error too large: {}", err.mass_error_ppm);
}

#[test]
fn within_tolerance() {
    let peptide = make_peptide(b"AR");
    let mass = peptide.mass();
    let charge = 2u8;
    let drift = mass * 5e-6;
    let mz_drifted = (mass + drift + charge as f64 * PROTON) / charge as f64;
    let spec = make_spectrum(mz_drifted, Some(charge as i32));
    let tol = PrecursorTolerance::symmetric(Tolerance::Ppm(20.0));
    assert!(matches_precursor(&spec, &peptide, charge, 0, &tol, 0.0).is_some());
}

#[test]
fn outside_tolerance() {
    let peptide = make_peptide(b"AR");
    let mass = peptide.mass();
    let charge = 2u8;
    let drift = mass * 50e-6;
    let mz_drifted = (mass + drift + charge as f64 * PROTON) / charge as f64;
    let spec = make_spectrum(mz_drifted, Some(charge as i32));
    let tol = PrecursorTolerance::symmetric(Tolerance::Ppm(20.0));
    assert!(matches_precursor(&spec, &peptide, charge, 0, &tol, 0.0).is_none());
}

#[test]
fn da_tolerance() {
    let peptide = make_peptide(b"AR");
    let mass = peptide.mass();
    let charge = 2u8;
    let mz_drifted = (mass + 0.005 + charge as f64 * PROTON) / charge as f64;
    let spec = make_spectrum(mz_drifted, Some(charge as i32));
    let tol = PrecursorTolerance::symmetric(Tolerance::Da(0.01));
    assert!(matches_precursor(&spec, &peptide, charge, 0, &tol, 0.0).is_some());
    let tol_tight = PrecursorTolerance::symmetric(Tolerance::Da(0.001));
    assert!(matches_precursor(&spec, &peptide, charge, 0, &tol_tight, 0.0).is_none());
}

#[test]
fn asymmetric_tolerance_rejects_excessive_negative_drift() {
    let peptide = make_peptide(b"AR");
    let mass = peptide.mass();
    let charge = 2u8;
    // Construct a spectrum where peptide is 15 ppm LIGHTER (negative error).
    let drift = mass * 15e-6;
    // spectrum implies a NEUTRAL mass of `mass + drift`. peptide_mass < spectrum mass.
    let spec_neutral = mass + drift;
    let mz_drifted = (spec_neutral + charge as f64 * PROTON) / charge as f64;
    let spec = make_spectrum(mz_drifted, Some(charge as i32));
    // Asymmetric: 5 ppm left (negative), 20 ppm right (positive). 15 ppm > 5 → reject.
    let tol = PrecursorTolerance::asymmetric(Tolerance::Ppm(5.0), Tolerance::Ppm(20.0));
    let result = matches_precursor(&spec, &peptide, charge, 0, &tol, 0.0);
    assert!(result.is_none(), "expected no match (15 ppm > 5 ppm left tolerance)");
}

#[test]
fn positive_shift_compensates_observed_bias() {
    let peptide = make_peptide(b"AR");
    let mass = peptide.mass();
    let charge = 2u8;
    // Spectrum reports +5 ppm heavy observed mass.
    let drift = mass * 5e-6;
    let mz_heavy = (mass + drift + charge as f64 * PROTON) / charge as f64;
    let spec = make_spectrum(mz_heavy, Some(charge as i32));
    let tol = PrecursorTolerance::symmetric(Tolerance::Ppm(1.0));
    assert!(
        matches_precursor(&spec, &peptide, charge, 0, &tol, 0.0).is_none(),
        "without shift, +5 ppm drift should miss 1 ppm tolerance"
    );
    let err = matches_precursor(&spec, &peptide, charge, 0, &tol, 5.0)
        .expect("shift should cancel +5 ppm bias");
    assert!(
        err.mass_error_ppm.abs() < 0.01,
        "residual ppm after shift: {}",
        err.mass_error_ppm
    );
}

#[test]
fn charge_zero_returns_none() {
    let peptide = make_peptide(b"AR");
    let spec = make_spectrum(100.0, Some(2));
    let tol = PrecursorTolerance::symmetric(Tolerance::Ppm(20.0));
    assert!(matches_precursor(&spec, &peptide, 0, 0, &tol, 0.0).is_none());
}
