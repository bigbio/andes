//! Physical and chemical constants of andes: the one place where element, isotope and
//! small-molecule masses are defined. Every crate imports them from here; none types its own
//! copy (`crates/model/tests/single_mass_constants.rs` checks this). Monoisotopic masses in
//! daltons; compounds are derived from the elements, so they cannot disagree.

/// Monoisotopic mass of ¹H, in Da: 1.00782503223(9), NIST Atomic Weights and Isotopic
/// Compositions (AME2020).
pub const H: f64 = 1.00782503223;

/// Monoisotopic mass of ²H (deuterium), in Da: 2.01410177812(12), NIST (AME2020).
pub const H2: f64 = 2.01410177812;

/// Monoisotopic mass of ¹⁶O, in Da: 15.99491461957(17), NIST (AME2020).
pub const O: f64 = 15.99491461957;

/// Monoisotopic mass of ¹²C, in Da (exact by definition of the dalton).
pub const C: f64 = 12.0;

/// Monoisotopic mass of ¹⁴N, in Da: 14.00307400443(20), NIST (AME2020).
pub const N: f64 = 14.00307400443;

/// Monoisotopic mass of ¹⁵N, in Da: 15.00010889888(64), NIST (AME2020).
pub const N15: f64 = 15.00010889888;

/// Monoisotopic mass of ³²S, in Da: 31.9720711744(14), NIST (AME2020).
pub const S: f64 = 31.9720711744;

/// Monoisotopic mass of ³¹P, in Da: 30.97376199842(70), NIST (AME2020).
pub const P: f64 = 30.97376199842;

/// Monoisotopic mass of water, in Da.
pub const H2O: f64 = H * 2.0 + O;

/// Monoisotopic mass of ammonia, in Da.
pub const NH3: f64 = N + H * 3.0;

/// Monoisotopic mass of carbon monoxide (b → a ion), in Da.
pub const CO: f64 = C + O;

/// Mass of a proton, in Da: 1.0072764665789(83), CODATA 2022 (proton mass in u).
pub const PROTON: f64 = 1.0072764665789;

/// Monoisotopic mass of ¹³C, in Da: 13.00335483507(23), NIST (AME2020).
pub const C13: f64 = 13.00335483507;

/// Spacing of the ¹³C isotope peaks, in Da (divide by the charge for m/z).
pub const ISOTOPE: f64 = C13 - C;

/// Offset of a z• ion from the y ion of the same residues, in Da: the y ion loses NH₃ and
/// gains a hydrogen atom, which is `-(N + 2H)`.
pub const Z_DOT_OFFSET: f64 = -(N + H * 2.0);

pub const INTEGER_MASS_SCALER: f32 = 0.999497;

pub fn nominal_from(mass: f64) -> i32 {
    (INTEGER_MASS_SCALER * mass as f32).round() as i32
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nominal_from_zero() {
        assert_eq!(nominal_from(0.0), 0);
    }

    #[test]
    fn nominal_from_glycine() {
        assert_eq!(nominal_from(57.02146), 57);
    }

    #[test]
    fn nominal_from_alanine() {
        assert_eq!(nominal_from(71.03711), 71);
    }

    #[test]
    fn nominal_from_tryptophan() {
        assert_eq!(nominal_from(186.07931), 186);
    }

    #[test]
    fn nominal_from_h2o() {
        assert_eq!(nominal_from(18.010565), 18);
    }

    #[test]
    fn nominal_from_one_kilodalton() {
        assert_eq!(nominal_from(1000.0), 999);
    }
}
