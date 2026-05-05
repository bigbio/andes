//! Mass tolerances. Mirrors Java `edu.ucsd.msjava.msgf.Tolerance`.

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Tolerance {
    Ppm(f64),
    Da(f64),
}

impl Tolerance {
    /// Convert this tolerance to absolute Daltons relative to a target mass.
    /// For `Da`, returns the constant; for `Ppm`, returns `mass * ppm * 1e-6`.
    pub fn as_da(&self, mass: f64) -> f64 {
        match self {
            Tolerance::Ppm(ppm) => mass * ppm * 1e-6,
            Tolerance::Da(da)   => *da,
        }
    }

    /// Return the raw numeric value stored in the tolerance — NOT converted to Da.
    ///
    /// For `Ppm(20.0)` this returns `20.0`; for `Da(0.5)` it returns `0.5`.
    /// Mirrors Java's `Tolerance.getValue()` which returns the stored primitive.
    pub fn raw_value(&self) -> f64 {
        match self {
            Tolerance::Ppm(v) => *v,
            Tolerance::Da(v)  => *v,
        }
    }
}

/// Asymmetric precursor mass tolerance. Phase B's calibrator produces
/// asymmetric `(left, right)` pairs; symmetric tolerances are a special
/// case constructed via `symmetric`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PrecursorTolerance {
    pub left:  Tolerance,
    pub right: Tolerance,
}

impl PrecursorTolerance {
    pub fn symmetric(t: Tolerance) -> Self {
        Self { left: t, right: t }
    }

    pub fn asymmetric(left: Tolerance, right: Tolerance) -> Self {
        Self { left, right }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ppm_at_1000_da() {
        // 10 ppm of 1000 Da = 0.01 Da
        let t = Tolerance::Ppm(10.0);
        assert_eq!(t.as_da(1000.0), 0.01);
    }

    #[test]
    fn ppm_at_500_da() {
        // 20 ppm of 500 Da = 0.01 Da
        let t = Tolerance::Ppm(20.0);
        assert_eq!(t.as_da(500.0), 0.01);
    }

    #[test]
    fn da_is_constant_under_mass() {
        let t = Tolerance::Da(0.5);
        assert_eq!(t.as_da(100.0), 0.5);
        assert_eq!(t.as_da(1000.0), 0.5);
        assert_eq!(t.as_da(0.0), 0.5);
    }

    #[test]
    fn precursor_symmetric_left_eq_right() {
        let t = PrecursorTolerance::symmetric(Tolerance::Ppm(10.0));
        assert_eq!(t.left.as_da(1000.0), t.right.as_da(1000.0));
    }

    #[test]
    fn precursor_asymmetric() {
        let t = PrecursorTolerance::asymmetric(Tolerance::Ppm(5.0), Tolerance::Ppm(20.0));
        assert_eq!(t.left.as_da(1000.0), 0.005);
        assert_eq!(t.right.as_da(1000.0), 0.02);
    }
}
