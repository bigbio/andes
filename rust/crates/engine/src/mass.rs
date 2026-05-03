//! Chemistry constants and mass utilities. Values pinned to Java
//! `edu.ucsd.msjava.msutil.Composition` and `Constants` — see
//! `tests/chemistry_constants_match_java.rs` for the parity gate.

/// Monoisotopic mass of hydrogen. Java: `Composition.H = 1.007825035`.
pub const H: f64 = 1.007825035;

/// Monoisotopic mass of oxygen. Java: `Composition.O = 15.99491463`.
pub const O: f64 = 15.99491463;

/// Monoisotopic mass of H2O, computed as `H * 2 + O` so the IEEE 754
/// rounding matches Java's `Composition.H2O` to the bit. The literal
/// `18.010565` is *not* bit-equal (mantissa drifts by 0x05).
pub const H2O: f64 = H * 2.0 + O;

/// Proton mass used as the default charge carrier. Java:
/// `Composition.PROTON = 1.00727649`.
pub const PROTON: f64 = 1.00727649;

/// Single-precision integer-mass scaler. Java:
/// `Constants.INTEGER_MASS_SCALER = 0.999497f`. Used in `nominal_from`
/// via float-domain arithmetic to mirror Java's
/// `AminoAcid.java:33: Math.round(INTEGER_MASS_SCALER * (float) mass)`.
pub const INTEGER_MASS_SCALER: f32 = 0.999497;

/// Convert a monoisotopic mass to the integer "nominal" mass that
/// indexes MS-GF+'s scoring DP table. Mirrors Java
/// `AminoAcid.java:33`: `Math.round(INTEGER_MASS_SCALER * (float) mass)`
/// — the multiply happens in f32 (single precision) before rounding.
/// Java's `Math.round(float)` is `floor(x + 0.5)`; for non-negative
/// inputs this matches Rust's `f32::round()`.
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
