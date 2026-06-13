//! Value-table test pinning `model::mass` chemistry constants to IUPAC
//! monoisotopic atomic masses and the andes nominal-mass scaler.

use model::mass::{nominal_from, C, H, H2O, INTEGER_MASS_SCALER, N, O, PROTON, S};

/// Bit-equality on f64 — masses must match to the full mantissa.
fn bit_eq(a: f64, b: f64) -> bool {
    a.to_bits() == b.to_bits()
}

#[test]
fn h_o_match_iupac_literals() {
    // IUPAC monoisotopic masses (Audi et al., Pure Appl. Chem.).
    assert_eq!(H.to_bits(),  1.007825035_f64.to_bits());
    assert_eq!(O.to_bits(), 15.99491463_f64.to_bits());
}

#[test]
fn h2o_matches_composition() {
    // H2O = 2×H + O; the naive literal 18.010565 is NOT bit-equal.
    assert_eq!(H2O.to_bits(), 0x403202b45e40fdf7);
    assert!(
        bit_eq(H2O, 1.007825035_f64 * 2.0 + 15.99491463_f64),
        "H2O drifted from H*2+O: rust=0x{:016x}", H2O.to_bits()
    );
}

#[test]
fn proton_matches_iupac() {
    assert_eq!(PROTON.to_bits(), 1.00727649_f64.to_bits());
}

#[test]
fn integer_mass_scaler_matches_residue_table_mean() {
    // Re-derive INTEGER_MASS_SCALER from the 20 standard residues:
    //   ratio = integer_formula_mass / monoisotopic_mass
    //   mean(ratio) ≈ 0.999497 (stored as f32).
    let compositions: &[(u32, u32, u32, u32, u32, f64)] = &[
        (2,  3, 1, 1, 0, 57.02146),   // G
        (3,  5, 1, 1, 0, 71.03711),   // A
        (3,  5, 1, 2, 0, 87.03203),   // S
        (5,  7, 1, 1, 0, 97.05276),   // P
        (5,  9, 1, 1, 0, 99.06841),   // V
        (4,  7, 1, 2, 0, 101.04768),  // T
        (3,  5, 1, 1, 1, 103.00919),  // C
        (6, 11, 1, 1, 0, 113.08406),  // L
        (6, 11, 1, 1, 0, 113.08406),  // I
        (4,  6, 2, 2, 0, 114.04293),  // N
        (4,  5, 1, 3, 0, 115.02694),  // D
        (5,  8, 2, 2, 0, 128.05858),  // Q
        (6, 12, 2, 1, 0, 128.09496),  // K
        (5,  7, 1, 3, 0, 129.04259),  // E
        (5,  9, 1, 1, 1, 131.04049),  // M
        (6,  7, 3, 1, 0, 137.05891),  // H
        (9,  9, 1, 1, 0, 147.06841),  // F
        (6, 12, 4, 1, 0, 156.10111),  // R
        (9,  9, 1, 2, 0, 163.06333),  // Y
        (11, 10, 2, 1, 0, 186.07931), // W
    ];
    let mut sum = 0.0f64;
    for &(c, h, n, o, s, mono) in compositions {
        let nominal = (c * 12 + h * 1 + n * 14 + o * 16 + s * 32) as f64;
        sum += nominal / mono;
    }
    let mean = sum / compositions.len() as f64;
    assert!(
        (mean - INTEGER_MASS_SCALER as f64).abs() < 1e-6,
        "INTEGER_MASS_SCALER drift: table mean={mean}, constant={INTEGER_MASS_SCALER}"
    );
    assert_eq!(INTEGER_MASS_SCALER.to_bits(), 0.999497_f32.to_bits());
}

#[test]
fn nominal_from_matches_f32_round_boundary() {
    assert_eq!(nominal_from(0.0), 0);
    assert_eq!(nominal_from(57.02146), 57);   // Gly
    assert_eq!(nominal_from(71.03711), 71);   // Ala
    assert_eq!(nominal_from(113.08406), 113); // Leu/Ile
    assert_eq!(nominal_from(186.07931), 186); // Trp
    assert_eq!(nominal_from(1000.0), 999);    // boundary anchoring f32 scaler
}

#[test]
fn c_n_s_match_iupac_literals() {
    assert_eq!(C.to_bits(), 12.0_f64.to_bits());
    assert_eq!(N.to_bits(), 14.003074_f64.to_bits());
    assert_eq!(S.to_bits(), 31.9720707_f64.to_bits());
}
