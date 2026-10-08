//! Precursor isolation purity: the share of the MS1 ion current inside the
//! isolation window that belongs to the selected precursor's isotope envelope.
//!
//! Follows the OpenMS `PrecursorPurity` / `IsobaricChannelExtractor` rule:
//! from the MS1 peak nearest the precursor m/z, walk the isotope ladder
//! (`ISOTOPE / z`) outward in both directions, matching each expected position
//! to the most intense peak within a ppm tolerance, and divide the summed
//! envelope intensity by everything in the window. Low purity means reporter
//! ions of co-isolated peptides compress the ratios of this PSM.

use model::mass::ISOTOPE;

/// Purity in `0..=1`, or `None` when the window is degenerate or empty.
///
/// `ms1` is m/z-sorted. `lower_off` / `upper_off` are the isolation offsets
/// below and above `precursor_mz` (mzML `isolation window lower/upper offset`).
/// `charge` 0 is treated as 1. `iso_tol_ppm` matches isotope peaks (OpenMS
/// default 10 ppm).
pub fn precursor_purity(
    ms1: &[(f64, f32)],
    precursor_mz: f64,
    charge: u8,
    lower_off: f64,
    upper_off: f64,
    iso_tol_ppm: f64,
) -> Option<f32> {
    let lo = precursor_mz - lower_off;
    let hi = precursor_mz + upper_off;
    if hi.partial_cmp(&lo) != Some(std::cmp::Ordering::Greater) || ms1.is_empty() {
        return None;
    }
    let start = ms1.partition_point(|&(mz, _)| mz < lo);
    let window = &ms1[start..];
    let end = window.partition_point(|&(mz, _)| mz <= hi);
    let window = &window[..end];
    let total: f64 = window.iter().map(|&(_, i)| i as f64).sum();
    if total <= 0.0 {
        return None;
    }
    let z = charge.max(1) as f64;
    let spacing = ISOTOPE / z;
    let tol_da = |mz: f64| (mz * iso_tol_ppm * 1e-6).max(0.001);

    // Anchor: most intense peak within tolerance of the precursor m/z. Without
    // one the precursor contributes nothing (purity 0), not "unknown".
    let anchor = best_within(window, precursor_mz, tol_da(precursor_mz));
    let Some((anchor_mz, anchor_int)) = anchor else {
        return Some(0.0);
    };
    let mut envelope = anchor_int as f64;
    for dir in [-1.0f64, 1.0] {
        let mut expected = anchor_mz + dir * spacing;
        let mut k = 0usize;
        while expected >= lo && expected <= hi && k < 10 {
            if let Some((mz, inten)) = best_within(window, expected, tol_da(expected)) {
                envelope += inten as f64;
                expected = mz + dir * spacing;
            } else {
                expected += dir * spacing;
            }
            k += 1;
        }
    }
    Some((envelope / total).min(1.0) as f32)
}

fn best_within(peaks: &[(f64, f32)], target: f64, tol: f64) -> Option<(f64, f32)> {
    let lo = target - tol;
    let hi = target + tol;
    let start = peaks.partition_point(|&(mz, _)| mz < lo);
    let mut best: Option<(f64, f32)> = None;
    for &(mz, inten) in &peaks[start..] {
        if mz > hi {
            break;
        }
        if best.is_none_or(|(_, b)| inten > b) {
            best = Some((mz, inten));
        }
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pure_precursor_is_one() {
        let mz = 500.25;
        let ms1 = vec![
            (mz, 1000.0f32),
            (mz + ISOTOPE / 2.0, 500.0),
            (mz + ISOTOPE, 150.0),
        ];
        let p = precursor_purity(&ms1, mz, 2, 0.7, 0.7, 10.0).unwrap();
        assert!((p - 1.0).abs() < 1e-6, "{p}");
    }

    #[test]
    fn interference_lowers_purity() {
        let mz = 500.25;
        let ms1 = vec![
            (mz - 0.3, 1000.0f32), // co-isolated, not on the ladder
            (mz, 1000.0),
            (mz + ISOTOPE / 2.0, 500.0),
        ];
        let p = precursor_purity(&ms1, mz, 2, 0.7, 0.7, 10.0).unwrap();
        assert!((p - 0.6).abs() < 1e-6, "{p}");
    }

    #[test]
    fn missing_precursor_peak_is_zero_and_empty_window_is_none() {
        let ms1 = vec![(400.0, 10.0f32), (600.0, 10.0)];
        assert_eq!(precursor_purity(&ms1, 500.0, 2, 0.7, 0.7, 10.0), None);
        let ms1 = vec![(500.3, 10.0f32)];
        assert_eq!(precursor_purity(&ms1, 500.0, 2, 0.7, 0.7, 10.0), Some(0.0));
    }

    #[test]
    fn walks_past_a_missing_isotope() {
        let mz = 800.0;
        // M+1 missing, M+2 present: still counted as envelope.
        let ms1 = vec![
            (mz, 100.0f32),
            (mz + 2.0 * ISOTOPE / 2.0, 50.0),
            (mz + 0.9, 50.0),
        ];
        let p = precursor_purity(&ms1, mz, 2, 1.5, 1.5, 10.0).unwrap();
        assert!((p - 0.75).abs() < 1e-6, "{p}");
    }
}
