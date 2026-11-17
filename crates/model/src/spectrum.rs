use crate::activation::ActivationMethod;

pub const MAX_PRECURSOR_CHARGE: i32 = u8::MAX as i32;

pub fn known_precursor_charge(z: i64) -> Option<i32> {
    (1..=i64::from(MAX_PRECURSOR_CHARGE))
        .contains(&z)
        .then_some(z as i32)
}

#[derive(Debug, Clone)]
pub struct Spectrum {
    pub title: String,
    pub precursor_mz: f64,
    pub precursor_intensity: Option<f32>,
    pub precursor_charge: Option<i32>,
    pub rt_seconds: Option<f64>,
    pub scan: Option<i32>,
    pub peaks: Vec<(f64, f32)>,
    pub activation_method: Option<ActivationMethod>,
    pub isolation_lower_offset: Option<f64>,
    pub isolation_upper_offset: Option<f64>,
    pub isolation_target_mz: Option<f64>,
}

impl Spectrum {
    pub fn len(&self) -> usize {
        self.peaks.len()
    }
    pub fn is_empty(&self) -> bool {
        self.peaks.is_empty()
    }
}

impl Default for Spectrum {
    fn default() -> Self {
        Spectrum {
            title: String::new(),
            precursor_mz: 0.0,
            precursor_intensity: None,
            precursor_charge: None,
            rt_seconds: None,
            scan: None,
            peaks: Vec::new(),
            activation_method: None,
            isolation_lower_offset: None,
            isolation_upper_offset: None,
            isolation_target_mz: None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn known_precursor_charge_bounds() {
        assert_eq!(known_precursor_charge(2), Some(2));
        assert_eq!(known_precursor_charge(255), Some(255));
        assert_eq!(known_precursor_charge(256), None);
        assert_eq!(known_precursor_charge(300), None);
        assert_eq!(known_precursor_charge(0), None);
        assert_eq!(known_precursor_charge(-2), None);
        assert_eq!(known_precursor_charge(i64::MAX), None);
    }

    fn make_spectrum() -> Spectrum {
        Spectrum {
            title: "Scan 100".to_string(),
            precursor_mz: 500.123,
            precursor_intensity: Some(1234.5),
            precursor_charge: Some(2),
            rt_seconds: Some(123.45),
            scan: Some(100),
            peaks: vec![(100.0, 1.0), (200.0, 2.0), (300.0, 3.0)],
            activation_method: None,
            isolation_lower_offset: None,
            isolation_upper_offset: None,
            isolation_target_mz: None,
        }
    }

    #[test]
    fn spectrum_isolation_offsets_default_none() {
        let s = Spectrum::default();
        assert_eq!(s.isolation_lower_offset, None);
        assert_eq!(s.isolation_upper_offset, None);
        assert_eq!(s.isolation_target_mz, None);
    }

    #[test]
    fn len_returns_peak_count() {
        let s = make_spectrum();
        assert_eq!(s.len(), 3);
    }

    #[test]
    fn is_empty_false_with_peaks() {
        let s = make_spectrum();
        assert!(!s.is_empty());
    }

    #[test]
    fn is_empty_true_no_peaks() {
        let s = Spectrum {
            title: "x".into(),
            precursor_mz: 0.0,
            precursor_intensity: None,
            precursor_charge: None,
            rt_seconds: None,
            scan: None,
            peaks: vec![],
            activation_method: None,
            isolation_lower_offset: None,
            isolation_upper_offset: None,
            isolation_target_mz: None,
        };
        assert!(s.is_empty());
        assert_eq!(s.len(), 0);
    }
}
