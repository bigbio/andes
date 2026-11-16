#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Tolerance {
    Ppm(f64),
    Da(f64),
}

impl Tolerance {
    pub fn as_da(&self, mass: f64) -> f64 {
        match self {
            Tolerance::Ppm(ppm) => mass * ppm * 1e-6,
            Tolerance::Da(da) => *da,
        }
    }

    pub fn raw_value(&self) -> f64 {
        match self {
            Tolerance::Ppm(v) => *v,
            Tolerance::Da(v) => *v,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PrecursorTolerance {
    pub left: Tolerance,
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
        let t = Tolerance::Ppm(10.0);
        assert_eq!(t.as_da(1000.0), 0.01);
    }

    #[test]
    fn ppm_at_500_da() {
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
