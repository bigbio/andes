//! Chemistry constants and mass utilities. See
//! `tests/chemistry_constants_match_java.rs` for the parity gate.

/// Monoisotopic mass of hydrogen.
pub const H: f64 = 1.007825035;

/// Monoisotopic mass of oxygen.
pub const O: f64 = 15.99491463;

/// Monoisotopic mass of carbon-12.
pub const C: f64 = 12.0;

/// Monoisotopic mass of nitrogen-14.
pub const N: f64 = 14.003074;

/// Monoisotopic mass of sulfur-32.
pub const S: f64 = 31.9720707;

/// Monoisotopic mass of H2O, computed as `H * 2 + O` so the IEEE 754
/// rounding matches the canonical bit pattern. The literal `18.010565`
/// is *not* bit-equal (mantissa drifts by 0x05).
pub const H2O: f64 = H * 2.0 + O;

/// Proton mass used as the default charge carrier.
pub const PROTON: f64 = 1.00727649;

/// Monoisotopic mass of carbon-13.
pub const C13: f64 = 13.00335483;

/// Mass difference between carbon-13 and carbon-12, used as the unit
/// step for isotope-error tolerance.
pub const ISOTOPE: f64 = C13 - C;

/// Single-precision integer-mass scaler. Used in `nominal_from` via
/// float-domain arithmetic; the multiply must happen in f32 (single
/// precision) before rounding to preserve the rounding boundary.
pub const INTEGER_MASS_SCALER: f32 = 0.999497;

/// Convert a monoisotopic mass to the integer "nominal" mass that
/// indexes MS-GF+'s scoring DP table.
///
/// The multiply happens in f32 (single precision) before rounding —
/// this is the rounding boundary the DP table is built against.
/// For non-negative inputs this matches `f32::round()` (round half-up).
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
        // 0.999497f * 57.02146f = 57.0228... → round → 57
        assert_eq!(nominal_from(57.02146), 57);
    }

    #[test]
    fn nominal_from_alanine() {
        // 0.999497f * 71.03711f = 71.001... → round → 71
        assert_eq!(nominal_from(71.03711), 71);
    }

    #[test]
    fn nominal_from_tryptophan() {
        // 0.999497f * 186.07931f = 185.9857... → round → 186
        assert_eq!(nominal_from(186.07931), 186);
    }

    #[test]
    fn nominal_from_h2o() {
        // 0.999497f * 18.010565f = 18.0014... → round → 18
        assert_eq!(nominal_from(18.010565), 18);
    }

    #[test]
    fn nominal_from_one_kilodalton() {
        // 0.999497f * 1000.0f = 999.497 → round → 999 (NOT 1000)
        // Anchors that the f32 scaler is in use; the f64 literal 0.9995
        // would give 1000 here.
        assert_eq!(nominal_from(1000.0), 999);
    }
}
