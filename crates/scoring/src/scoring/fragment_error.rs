//! A run's measured fragment mass error and the tolerances derived from it.

use std::sync::atomic::{AtomicU64, Ordering};

/// Spreads added to the offset for the PIN feature tolerance.
pub const TOL_SIGMAS: f64 = 4.0;
/// Narrowest derived feature tolerance (ppm).
pub const MIN_FEATURE_TOL_PPM: f64 = 5.0;
/// Widest derived feature tolerance (ppm).
pub const MAX_FEATURE_TOL_PPM: f64 = 40.0;
/// Width of the ppm-Gaussian term relative to the measured spread. A width of
/// 1.5-1.8 times the fragment-error spread scored best on Astral.
pub const GAUSSIAN_SIGMA_SCALE: f64 = 1.5;
/// ppm-Gaussian width when no measurement exists.
pub const DEFAULT_GAUSSIAN_SIGMA_PPM: f64 = 7.0;
/// Narrowest ppm-Gaussian width (ppm).
pub const MIN_GAUSSIAN_SIGMA_PPM: f64 = 1.0;

/// Fragment mass error of a run, measured on its confident PSMs.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FragmentErrorModel {
    /// Median signed error (observed - theoretical), ppm.
    pub offset_ppm: f64,
    /// Robust spread (normal-consistent MAD) of the errors, ppm.
    pub sigma_ppm: f64,
    /// Ions the estimate rests on.
    pub ions: usize,
    /// Confident PSMs the ions came from.
    pub psms: usize,
}

impl FragmentErrorModel {
    pub fn new(offset_ppm: f64, sigma_ppm: f64, ions: usize, psms: usize) -> Self {
        Self {
            offset_ppm,
            sigma_ppm,
            ions,
            psms,
        }
    }

    /// Symmetric feature tolerance covering the error distribution:
    /// `|offset| + 4 sigma`, clamped to 5-40 ppm.
    pub fn feature_tolerance_ppm(&self) -> f64 {
        (self.offset_ppm.abs() + TOL_SIGMAS * self.sigma_ppm)
            .clamp(MIN_FEATURE_TOL_PPM, MAX_FEATURE_TOL_PPM)
    }

    /// Width of the ppm-Gaussian term.
    pub fn gaussian_sigma_ppm(&self) -> f64 {
        (GAUSSIAN_SIGMA_SCALE * self.sigma_ppm).max(MIN_GAUSSIAN_SIGMA_PPM)
    }
}

/// A [`FragmentErrorModel`] shared by the scoring threads: set between files,
/// read per PSM without locking. Only the offset and spread are kept.
#[derive(Debug)]
pub(crate) struct SharedFragmentError {
    offset_bits: AtomicU64,
    sigma_bits: AtomicU64,
}

impl SharedFragmentError {
    pub(crate) fn new() -> Self {
        Self {
            offset_bits: AtomicU64::new(f64::NAN.to_bits()),
            sigma_bits: AtomicU64::new(f64::NAN.to_bits()),
        }
    }

    pub(crate) fn set(&self, m: Option<FragmentErrorModel>) {
        let (o, s) = m.map_or((f64::NAN, f64::NAN), |m| (m.offset_ppm, m.sigma_ppm));
        self.offset_bits.store(o.to_bits(), Ordering::Relaxed);
        self.sigma_bits.store(s.to_bits(), Ordering::Relaxed);
    }

    pub(crate) fn get(&self) -> Option<FragmentErrorModel> {
        let o = f64::from_bits(self.offset_bits.load(Ordering::Relaxed));
        let s = f64::from_bits(self.sigma_bits.load(Ordering::Relaxed));
        (o.is_finite() && s.is_finite()).then(|| FragmentErrorModel::new(o, s, 0, 0))
    }
}

impl Clone for SharedFragmentError {
    fn clone(&self) -> Self {
        let c = Self::new();
        c.set(self.get());
        c
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tolerance_covers_offset_and_four_spreads_within_bounds() {
        // Astral-like: -2.5 ppm, spread 1.6 ppm.
        let astral = FragmentErrorModel::new(-2.5, 1.6, 5000, 800);
        assert!((astral.feature_tolerance_ppm() - 8.9).abs() < 1e-9);
        assert!((astral.gaussian_sigma_ppm() - 2.4).abs() < 1e-9);
        // Very tight data is floored, very loose data capped.
        assert_eq!(
            FragmentErrorModel::new(0.1, 0.3, 1, 1).feature_tolerance_ppm(),
            MIN_FEATURE_TOL_PPM
        );
        assert_eq!(
            FragmentErrorModel::new(12.0, 9.0, 1, 1).feature_tolerance_ppm(),
            MAX_FEATURE_TOL_PPM
        );
        assert_eq!(
            FragmentErrorModel::new(0.0, 0.2, 1, 1).gaussian_sigma_ppm(),
            MIN_GAUSSIAN_SIGMA_PPM
        );
    }

    #[test]
    fn shared_value_round_trips_and_clears() {
        let s = SharedFragmentError::new();
        assert_eq!(s.get(), None);
        s.set(Some(FragmentErrorModel::new(4.2, 3.0, 10, 2)));
        let got = s.get().unwrap();
        assert_eq!((got.offset_ppm, got.sigma_ppm), (4.2, 3.0));
        assert_eq!(s.clone().get().map(|m| m.offset_ppm), Some(4.2));
        s.set(None);
        assert_eq!(s.get(), None);
    }
}
