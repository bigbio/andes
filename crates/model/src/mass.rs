pub const H: f64 = 1.007825035;

pub const O: f64 = 15.99491463;

pub const C: f64 = 12.0;

pub const N: f64 = 14.003074;

pub const S: f64 = 31.9720707;

pub const H2O: f64 = H * 2.0 + O;

pub const PROTON: f64 = 1.00727649;

pub const C13: f64 = 13.00335483;

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
