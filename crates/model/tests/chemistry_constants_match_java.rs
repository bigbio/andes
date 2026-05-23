//! Value-table test pinning `model::mass` constants to Java's
//! `edu.ucsd.msjava.msutil.Composition` and `Constants`. References are
//! the actual IEEE 754 bit patterns Java produces — verified against
//! the Java source, not against the same Rust literals.

use model::mass::{nominal_from, C, H, H2O, INTEGER_MASS_SCALER, N, O, PROTON, S};

/// Bit-equality on f64 — masses must match Java to the full mantissa.
fn bit_eq(a: f64, b: f64) -> bool {
    a.to_bits() == b.to_bits()
}

#[test]
fn h_o_match_java_literals() {
    // Source: src/main/java/edu/ucsd/msjava/msutil/Composition.java
    //   public static final double H = 1.007825035;
    //   public static final double O = 15.99491463;
    assert_eq!(H.to_bits(),  1.007825035_f64.to_bits());
    assert_eq!(O.to_bits(), 15.99491463_f64.to_bits());
}

#[test]
fn h2o_matches_java_computed() {
    // Java: public static final double H2O = H * 2 + O;
    // The IEEE 754 result is 18.0105647... (bit pattern 0x403202b45e40fdf7).
    // The naive literal 18.010565 is NOT bit-equal — that drift is what
    // this test exists to catch.
    assert_eq!(H2O.to_bits(), 0x403202b45e40fdf7);
    assert!(
        bit_eq(H2O, 1.007825035_f64 * 2.0 + 15.99491463_f64),
        "H2O drifted from H*2+O: rust=0x{:016x}", H2O.to_bits()
    );
}

#[test]
fn proton_matches_java() {
    // Source: Composition.java line 30: public static final double PROTON = 1.00727649;
    assert_eq!(PROTON.to_bits(), 1.00727649_f64.to_bits());
}

#[test]
fn integer_mass_scaler_matches_java() {
    // Source: Constants.java line 13:
    //   public static final float INTEGER_MASS_SCALER = 0.999497f;
    assert_eq!(INTEGER_MASS_SCALER.to_bits(), 0.999497_f32.to_bits());
}

#[test]
fn nominal_from_matches_java_aminoacid_constructor() {
    // Reference values: each computed by `Math.round(INTEGER_MASS_SCALER * (float) mass)`
    // exactly as Java AminoAcid.java:33 does it.
    assert_eq!(nominal_from(0.0), 0);
    assert_eq!(nominal_from(57.02146), 57);   // Gly
    assert_eq!(nominal_from(71.03711), 71);   // Ala
    assert_eq!(nominal_from(113.08406), 113); // Leu/Ile
    assert_eq!(nominal_from(186.07931), 186); // Trp
    assert_eq!(nominal_from(1000.0), 999);    // boundary anchoring f32 scaler
}

#[test]
fn c_n_s_match_java_literals() {
    // Source: Composition.java
    //   public static final double C = 12.0;
    //   public static final double N = 14.003074;
    //   public static final double S = 31.9720707;
    assert_eq!(C.to_bits(), 12.0_f64.to_bits());
    assert_eq!(N.to_bits(), 14.003074_f64.to_bits());
    assert_eq!(S.to_bits(), 31.9720707_f64.to_bits());
}
