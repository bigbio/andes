//! Physical and chemical constants of andes: the one place where element, isotope and
//! small-molecule masses are defined. Every crate imports them from here; none types its own
//! copy (`crates/model/tests/single_mass_constants.rs` checks this). Monoisotopic masses in
//! daltons; compounds are derived from the elements, so they cannot disagree.

/// Monoisotopic mass of ¹H, in Da.
pub const H: f64 = 1.007825035;

/// Monoisotopic mass of ¹⁶O, in Da.
pub const O: f64 = 15.99491463;

/// Monoisotopic mass of ¹²C, in Da (exact by definition of the dalton).
pub const C: f64 = 12.0;

/// Monoisotopic mass of ¹⁴N, in Da.
pub const N: f64 = 14.003074;

/// Monoisotopic mass of ³²S, in Da.
pub const S: f64 = 31.9720707;

/// Monoisotopic mass of ³¹P, in Da.
pub const P: f64 = 30.973761998;

/// Monoisotopic mass of water, in Da.
pub const H2O: f64 = H * 2.0 + O;

/// Monoisotopic mass of ammonia, in Da. Kept as the typed value rather than derived as
/// `N + 3.0 * H` (17.026549105, 1.05e-7 Da higher) so that defining the constants in one
/// place leaves search results unchanged; moving to the derived value is its own measured
/// change.
pub const NH3: f64 = 17.026549;

/// Monoisotopic mass of carbon monoxide (b → a ion), in Da.
pub const CO: f64 = C + O;

/// Mass of a proton, in Da.
pub const PROTON: f64 = 1.00727649;

/// Monoisotopic mass of ¹³C, in Da.
pub const C13: f64 = 13.00335483;

/// Spacing of the ¹³C isotope peaks, in Da (divide by the charge for m/z).
pub const ISOTOPE: f64 = C13 - C;

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
