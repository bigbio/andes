//! Spectrum — a single tandem MS scan. Lives in `input` for Phase 3a;
//! relocates to `engine` when Phase 4 search lands and needs to consume it.

#[derive(Debug, Clone)]
pub struct Spectrum {
    /// MGF `TITLE=` value (or `<spectrumRef>` for mzML in Phase 3b).
    /// Used as the PSM `SpecID` column in `.pin` output.
    pub title: String,
    /// `PEPMASS=` first value: precursor m/z.
    pub precursor_mz: f64,
    /// `PEPMASS=` second value (optional): precursor intensity.
    pub precursor_intensity: Option<f32>,
    /// `CHARGE=` value, e.g. `2+`. None when absent.
    pub precursor_charge: Option<i32>,
    /// `RTINSECONDS=` value. None when absent.
    pub rt_seconds: Option<f64>,
    /// `SCANS=` value (scan number). None when absent.
    pub scan: Option<i32>,
    /// Peak list: (m/z f64, intensity f32). Sorted ascending by m/z by
    /// the parser.
    pub peaks: Vec<(f64, f32)>,
}

impl Spectrum {
    pub fn len(&self) -> usize { self.peaks.len() }
    pub fn is_empty(&self) -> bool { self.peaks.is_empty() }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_spectrum() -> Spectrum {
        Spectrum {
            title: "Scan 100".to_string(),
            precursor_mz: 500.123,
            precursor_intensity: Some(1234.5),
            precursor_charge: Some(2),
            rt_seconds: Some(123.45),
            scan: Some(100),
            peaks: vec![(100.0, 1.0), (200.0, 2.0), (300.0, 3.0)],
        }
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
            title: "x".into(), precursor_mz: 0.0, precursor_intensity: None,
            precursor_charge: None, rt_seconds: None, scan: None,
            peaks: vec![],
        };
        assert!(s.is_empty());
        assert_eq!(s.len(), 0);
    }
}
