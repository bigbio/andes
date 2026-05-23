//! Spectrum — a single tandem MS scan.

use crate::activation::ActivationMethod;

#[derive(Debug, Clone)]
pub struct Spectrum {
    /// MGF `TITLE=` value (or `<spectrumRef>` for mzML).
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
    /// Activation method recorded in the source file (mzML `<activation>`
    /// cvParam, or `ACTIVATION=` in MGF). `None` when the source doesn't
    /// record one. This is *informational* — used by the CLI binary to
    /// auto-route to the matching bundled `.param` file when the user
    /// hasn't overridden `--param-file`/`--fragmentation`/`--instrument`.
    /// It is NOT used by the scoring loop directly.
    pub activation_method: Option<ActivationMethod>,
}

impl Spectrum {
    pub fn len(&self) -> usize { self.peaks.len() }
    pub fn is_empty(&self) -> bool { self.peaks.is_empty() }
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
        }
    }
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
            activation_method: None,
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
            activation_method: None,
        };
        assert!(s.is_empty());
        assert_eq!(s.len(), 0);
    }
}
