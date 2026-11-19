use model::mass::{
    nominal_from, C, C13, CO, H, H2, H2O, INTEGER_MASS_SCALER, ISOTOPE, N, N15, NH3, O, P, PROTON,
    S, Z_DOT_OFFSET,
};

fn bit_eq(a: f64, b: f64) -> bool {
    a.to_bits() == b.to_bits()
}

#[test]
fn h_o_match_iupac_literals() {
    assert_eq!(H.to_bits(), 1.00782503223_f64.to_bits());
    assert_eq!(O.to_bits(), 15.99491461957_f64.to_bits());
}

#[test]
fn h2o_matches_composition() {
    assert_eq!(
        H2O.to_bits(),
        (1.00782503223_f64 * 2.0 + 15.99491461957_f64).to_bits()
    );
    assert!(
        bit_eq(H2O, 1.00782503223_f64 * 2.0 + 15.99491461957_f64),
        "H2O drifted from H*2+O: rust=0x{:016x}",
        H2O.to_bits()
    );
}

#[test]
fn proton_matches_iupac() {
    assert_eq!(PROTON.to_bits(), 1.0072764665789_f64.to_bits());
}

#[test]
fn integer_mass_scaler_matches_residue_table_mean() {
    let compositions: &[(u32, u32, u32, u32, u32, f64)] = &[
        (2, 3, 1, 1, 0, 57.02146),
        (3, 5, 1, 1, 0, 71.03711),
        (3, 5, 1, 2, 0, 87.03203),
        (5, 7, 1, 1, 0, 97.05276),
        (5, 9, 1, 1, 0, 99.06841),
        (4, 7, 1, 2, 0, 101.04768),
        (3, 5, 1, 1, 1, 103.00919),
        (6, 11, 1, 1, 0, 113.08406),
        (6, 11, 1, 1, 0, 113.08406),
        (4, 6, 2, 2, 0, 114.04293),
        (4, 5, 1, 3, 0, 115.02694),
        (5, 8, 2, 2, 0, 128.05858),
        (6, 12, 2, 1, 0, 128.09496),
        (5, 7, 1, 3, 0, 129.04259),
        (5, 9, 1, 1, 1, 131.04049),
        (6, 7, 3, 1, 0, 137.05891),
        (9, 9, 1, 1, 0, 147.06841),
        (6, 12, 4, 1, 0, 156.10111),
        (9, 9, 1, 2, 0, 163.06333),
        (11, 10, 2, 1, 0, 186.07931),
    ];
    let mut sum = 0.0f64;
    for &(c, h, n, o, s, mono) in compositions {
        let nominal = (c * 12 + h + n * 14 + o * 16 + s * 32) as f64;
        sum += nominal / mono;
    }
    let mean = sum / compositions.len() as f64;
    assert!(
        (mean - INTEGER_MASS_SCALER as f64).abs() < 1e-3,
        "INTEGER_MASS_SCALER sanity: table mean={mean}, constant={INTEGER_MASS_SCALER} \
         (loose bound; exact value pinned by assert_eq below)"
    );
    assert_eq!(INTEGER_MASS_SCALER.to_bits(), 0.999497_f32.to_bits());
}

#[test]
fn nominal_from_matches_f32_round_boundary() {
    assert_eq!(nominal_from(0.0), 0);
    assert_eq!(nominal_from(57.02146), 57);
    assert_eq!(nominal_from(71.03711), 71);
    assert_eq!(nominal_from(113.08406), 113);
    assert_eq!(nominal_from(186.07931), 186);
    assert_eq!(nominal_from(1000.0), 999);
}

#[test]
fn c_n_s_match_iupac_literals() {
    assert_eq!(C.to_bits(), 12.0_f64.to_bits());
    assert_eq!(N.to_bits(), 14.00307400443_f64.to_bits());
    assert_eq!(S.to_bits(), 31.9720711744_f64.to_bits());
    assert_eq!(P.to_bits(), 30.97376199842_f64.to_bits());
    assert_eq!(C13.to_bits(), 13.00335483507_f64.to_bits());
    assert_eq!(N15.to_bits(), 15.00010889888_f64.to_bits());
    assert_eq!(H2.to_bits(), 2.01410177812_f64.to_bits());
}

#[test]
fn derived_masses_follow_their_elements() {
    assert_eq!(NH3.to_bits(), (N + 3.0 * H).to_bits());
    assert_eq!(CO.to_bits(), (C + O).to_bits());
    assert_eq!(ISOTOPE.to_bits(), (C13 - C).to_bits());
    assert_eq!(Z_DOT_OFFSET.to_bits(), (-(N + 2.0 * H)).to_bits());
    assert!((NH3 - 17.0265491015).abs() < 1e-9);
    assert!((H2O - 18.0105646837).abs() < 1e-9);
    assert!((Z_DOT_OFFSET + 16.0187240).abs() < 1e-6);
}
