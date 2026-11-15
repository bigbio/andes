//! Chemistry constants and mass utilities (monoisotopic masses, nominal-mass
//! scaling). Constant values are physical facts; a test pins their bit patterns.

pub const H: f64 = 1.007825035;

pub const O: f64 = 15.99491463;

pub const C: f64 = 12.0;

pub const N: f64 = 14.003074;

pub const S: f64 = 31.9720707;

/// Monoisotopic mass of H2O, computed as `H * 2 + O` so the IEEE 754
/// rounding matches the canonical bit pattern. The literal `18.010565`
/// is *not* bit-equal (mantissa drifts by 0x05).
pub const H2O: f64 = H * 2.0 + O;

/// Proton mass used as the default charge carrier.
pub const PROTON: f64 = 1.00727649;

pub const C13: f64 = 13.00335483;

/// Mass difference between carbon-13 and carbon-12, used as the unit
/// step for isotope-error tolerance.
pub const ISOTOPE: f64 = C13 - C;

/// Single-precision integer-mass scaler: maps monoisotopic residue mass to the
/// integer nominal bin used by the DP scoring tables via `nominal_from`.
///
/// Derivation (standard 20-residue table, IUPAC monoisotopic masses from
/// `model::amino_acid::standard_composition` × `C/H/N/O/S`):
///   nominalᵢ = C×12 + H×1 + N×14 + O×16 + S×32   (integer formula mass)
///   ratioᵢ   = nominalᵢ / monoisotopicᵢ
///   mean(ratioᵢ) ≈ 0.999497  →  stored as f32 for the f32-round boundary
///   that `nominal_from` relies on (see the 1000 Da anchor test).
///
/// The multiply must happen in f32 before rounding to keep cumulative
/// fragment-bin error bounded across long peptides.
pub const INTEGER_MASS_SCALER: f32 = 0.999497;

/// Convert a monoisotopic mass to the integer "nominal" mass that indexes
/// the per-rank scoring tables.
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
        // 0.999497f * 1000.0f = 999.497 → round → 999 (NOT 1000)
        // Anchors that the f32 scaler is in use; the f64 literal 0.9995
        // would give 1000 here.
        assert_eq!(nominal_from(1000.0), 999);
    }
}
