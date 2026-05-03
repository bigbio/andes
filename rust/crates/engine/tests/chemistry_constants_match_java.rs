//! Value-table test pinning `engine::mass` constants to Java's
//! `edu.ucsd.msjava.msutil.Composition`. If a Java constant ever changes
//! upstream, or a Rust re-typed value drifts, this test fails loudly.

use engine::mass::{H2O, PROTON};

/// Bit-equality on f64 — masses must match Java to the full mantissa.
fn bit_eq(a: f64, b: f64) -> bool {
    a.to_bits() == b.to_bits()
}

#[test]
fn h2o_matches_java() {
    // Java: Composition.H2O = 18.010565
    assert!(
        bit_eq(H2O, 18.010565),
        "H2O drifted: rust={H2O}, expected 18.010565"
    );
}

#[test]
fn proton_matches_java() {
    // Java: IsotopeMass.getProtonMass() = 1.007276466
    assert!(
        bit_eq(PROTON, 1.007276466),
        "PROTON drifted: rust={PROTON}, expected 1.007276466"
    );
}

#[test]
fn nominal_from_well_known_values() {
    // Reference values computed with Java's `(int)(mass * 0.9995f + 0.5)`.
    use engine::mass::nominal_from;
    assert_eq!(nominal_from(0.0), 0);
    assert_eq!(nominal_from(57.02146), 57);   // Gly
    assert_eq!(nominal_from(71.03711), 71);   // Ala
    assert_eq!(nominal_from(186.07931), 186); // Trp
    assert_eq!(nominal_from(1000.0), 1000);
}
