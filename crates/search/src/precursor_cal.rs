//! Two-pass precursor mass calibration (Java P2-cal / `MassCalibrator`).
//!
//! Phase 0–1: helpers + mode wiring. The pre-pass calibrator (Phase 3) learns
//! a file-wide ppm shift; the main pass applies it to observed neutral masses
//! without mutating [`model::Spectrum`] objects.

/// Java `-precursorCal` modes.
///
/// `Default` is `Off` until the G1 ship gate closes — matches the CLI default
/// and `SearchParams::default_tryptic`, so library consumers that derive
/// `Default` on a struct containing this field cannot silently enable the
/// pre-pass.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PrecursorCalMode {
    /// Run pre-pass; apply shift only when ≥200 confident PSMs yield non-zero median.
    Auto,
    /// Run pre-pass; always store learned shift (still no-op when shift is 0.0).
    On,
    /// Skip pre-pass; `precursor_mass_shift_ppm` stays 0.0 (bit-identical baseline).
    #[default]
    Off,
}

/// Sample every `stride`-th element, capped at `cap`. Mirrors Java
/// `MassCalibrator.sampleEveryNth`.
pub fn sample_every_nth<T: Clone>(source: &[T], stride: usize, cap: usize) -> Vec<T> {
    if source.is_empty() || stride == 0 || cap == 0 {
        return Vec::new();
    }
    let mut out = Vec::with_capacity(cap.min(source.len() / stride + 1));
    for (i, item) in source.iter().enumerate() {
        if i % stride != 0 {
            continue;
        }
        out.push(item.clone());
        if out.len() >= cap {
            break;
        }
    }
    out
}

/// Residual in ppm: `(observed - theoretical) / theoretical * 1e6`.
pub fn residual_ppm(observed_mass: f64, theoretical_mass: f64) -> f64 {
    (observed_mass - theoretical_mass) / theoretical_mass * 1e6
}

/// Median of `values`. Empty input returns 0.0 (Java "no shift" contract).
pub fn median(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let mut copy: Vec<f64> = values.to_vec();
    copy.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = copy.len();
    if n % 2 == 1 {
        copy[n / 2]
    } else {
        (copy[n / 2 - 1] + copy[n / 2]) / 2.0
    }
}

/// Apply learned calibration to an observed neutral peptide mass.
///
/// When `shift_ppm == 0.0`, returns `raw_neutral` unchanged — the mandatory
/// fast path for `-precursorCal off` and the pre-calibration baseline.
pub fn adjusted_observed_neutral_mass(raw_neutral: f64, shift_ppm: f64) -> f64 {
    if shift_ppm == 0.0 {
        raw_neutral
    } else {
        raw_neutral * (1.0 - shift_ppm * 1e-6)
    }
}

/// Pre-pass tuning constants (Java `MassCalibrator`).
pub mod constants {
    pub const SAMPLING_STRIDE: usize = 10;
    pub const MAX_SAMPLED: usize = 500;
    pub const MIN_CONFIDENT_PSMS: usize = 200;
    pub const MAX_SPEC_EVALUE: f64 = 1e-6;
    pub const MIN_SPECKEYS_FOR_PREPASS: usize = 10_000;

    /// Java `DEFAULT_TIGHTENED_WINDOW_*` — post-cal main-pass tolerance tightening.
    pub const TIGHTENED_WINDOW_FLOOR_PPM: f64 = 2.0;
    pub const TIGHTENED_WINDOW_MARGIN_PPM: f64 = 0.5;
    pub const TIGHTENED_WINDOW_SIGMA_MULTIPLIER: f64 = 3.0;
    pub const MAD_TO_SIGMA_SCALE: f64 = 1.4826;
}

/// Median absolute deviation from `center`.
pub fn median_absolute_deviation(values: &[f64], center: f64) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let deviations: Vec<f64> = values.iter().map(|v| (v - center).abs()).collect();
    median(&deviations)
}

/// Robust Gaussian-equivalent sigma (MAD × 1.4826).
pub fn robust_sigma_ppm(residuals: &[f64], center: f64) -> f64 {
    constants::MAD_TO_SIGMA_SCALE * median_absolute_deviation(residuals, center)
}

/// Conservative tightened ppm half-window for a calibrated main pass.
pub fn tightened_tolerance_ppm(
    user_ppm: f64,
    robust_sigma_ppm: f64,
    sigma_multiplier: f64,
    floor_ppm: f64,
    margin_ppm: f64,
) -> f64 {
    if user_ppm <= 0.0 {
        return user_ppm;
    }
    let tightened = floor_ppm.max(sigma_multiplier * robust_sigma_ppm + margin_ppm);
    user_ppm.min(tightened)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn median_odd_even_empty_and_unsorted() {
        assert_eq!(median(&[]), 0.0);
        assert_eq!(median(&[1.0, 3.0, 5.0]), 3.0);
        assert_eq!(median(&[1.0, 2.0, 3.0, 4.0]), 2.5);
        assert_eq!(median(&[5.0, 1.0, 3.0]), 3.0);
        assert_eq!(median(&[1.0, 2.0, 3.0, 4.0, 1000.0]), 3.0);
        assert_eq!(median(&[7.5]), 7.5);
    }

    #[test]
    fn residual_ppm_sign_convention() {
        assert!((residual_ppm(1001.0, 1000.0) - 1000.0).abs() < 0.5);
        assert!((residual_ppm(999.0, 1000.0) + 1000.0).abs() < 0.5);
        assert_eq!(residual_ppm(1000.0, 1000.0), 0.0);
        let observed = 1000.0 + 1000.0 * 5e-6;
        assert!((residual_ppm(observed, 1000.0) - 5.0).abs() < 1e-6);
    }

    #[test]
    fn sample_every_nth_stride_and_cap() {
        let source: Vec<i32> = (0..100).collect();
        let sampled = sample_every_nth(&source, 10, 500);
        assert_eq!(sampled.len(), 10);
        assert_eq!(sampled[0], 0);
        assert_eq!(sampled[9], 90);

        let big: Vec<i32> = (0..10_000).collect();
        assert_eq!(sample_every_nth(&big, 10, 500).len(), 500);
        assert!(sample_every_nth::<i32>(&[], 10, 500).is_empty());
        assert_eq!(sample_every_nth(&[0, 1, 2], 10, 500), vec![0]);
    }

    #[test]
    fn adjusted_mass_zero_shift_is_identity() {
        assert_eq!(adjusted_observed_neutral_mass(1234.567, 0.0), 1234.567);
        let adj = adjusted_observed_neutral_mass(1000.0, 5.0);
        assert!((adj - (1000.0 * (1.0 - 5e-6))).abs() < 1e-9);
    }

    #[test]
    fn robust_sigma_matches_java_mad_scale() {
        let residuals = vec![9.0, 10.0, 11.0];
        assert!((robust_sigma_ppm(&residuals, 10.0) - 1.4826).abs() < 1e-6);
    }

    #[test]
    fn tightened_tolerance_respects_floor_and_cap() {
        assert!((tightened_tolerance_ppm(10.0, 0.2, 3.0, 2.0, 0.5) - 2.0).abs() < 1e-6);
        assert!((tightened_tolerance_ppm(1.5, 0.2, 3.0, 2.0, 0.5) - 1.5).abs() < 1e-6);
        assert!((tightened_tolerance_ppm(12.0, 1.0, 3.0, 2.0, 0.5) - 3.5).abs() < 1e-6);
    }
}
