//! Chemistry constants and mass utilities. Values pinned to Java
//! `edu.ucsd.msjava.msutil.Composition` — see the value-table test in
//! `tests/chemistry_constants_match_java.rs`.

/// Monoisotopic mass of H2O. Java: `Composition.H2O`. Used as the
/// neutral-mass adjustment in `Peptide::mass`.
pub const H2O: f64 = 18.010565;

/// Monoisotopic mass of a proton (charge carrier). Java:
/// `IsotopeMass.getProtonMass()` (1.007276466).
pub const PROTON: f64 = 1.007276466;

/// Convert a monoisotopic mass to the integer "nominal" mass used by
/// MS-GF+'s scoring DP table. Magic constant 0.9995 mirrors Java
/// `Constants.NOMINAL_MASS_FACTOR`.
pub fn nominal_from(mass: f64) -> i32 {
    (mass * 0.9995 + 0.5) as i32
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nominal_from_zero() {
        assert_eq!(nominal_from(0.0), 0);
    }

    #[test]
    fn nominal_from_small_positive() {
        // 1.0 Da → 1 (1.0 * 0.9995 + 0.5 = 1.4995 → trunc to 1)
        assert_eq!(nominal_from(1.0), 1);
    }

    #[test]
    fn nominal_from_typical_aa() {
        // Glycine residue mass 57.02146 → nominal 57
        // 57.02146 * 0.9995 + 0.5 = 57.49294... → trunc to 57
        assert_eq!(nominal_from(57.02146), 57);
    }

    #[test]
    fn nominal_from_typical_peptide() {
        // ~1000 Da peptide → 1000
        // 1000.0 * 0.9995 + 0.5 = 1000.0 → trunc to 1000
        assert_eq!(nominal_from(1000.0), 1000);
    }

    #[test]
    fn nominal_from_h2o() {
        // 18.010565 * 0.9995 + 0.5 = 18.50156... → trunc 18
        assert_eq!(nominal_from(18.010565), 18);
    }
}
