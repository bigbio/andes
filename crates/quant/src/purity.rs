//! Precursor isolation purity: the share of the MS1 ion current inside the
//! isolation window that belongs to the selected precursor's isotope envelope.
//!
//! Follows the reference purity-filter rule:
//! from the MS1 peak nearest the precursor m/z, walk the isotope ladder
//! (`ISOTOPE / z`) outward in both directions, matching each expected position
//! to the most intense peak within a ppm tolerance, and divide the summed
//! envelope intensity inside the window by everything in the window. Low
//! purity means reporter ions of co-isolated peptides compress the ratios of
//! this PSM.
//!
//! The window is centred on the isolation target m/z the instrument used, not
//! on the precursor m/z: a monoisotopic precursor assigned after acquisition
//! can sit below the window when the instrument isolated a heavier isotope.
//! The ladder still starts at the precursor m/z, so isotopes of the precursor
//! that fall inside the window count wherever the anchor is.

use model::mass::ISOTOPE;

/// Isotope steps counted inside the window, per direction.
const MAX_WINDOW_STEPS: usize = 10;
/// Ladder steps per direction overall, including those outside the window.
const MAX_LADDER_STEPS: usize = 64;

/// Purity in `0..=1`, or `None` when the window is degenerate or empty.
///
/// `ms1` is m/z-sorted. The window is `[centre - lower_off, centre + upper_off]`
/// where `centre` is `target_mz` (mzML `isolation window target m/z`) or, when
/// the source records no target, `precursor_mz`. `lower_off` / `upper_off` are
/// the mzML `isolation window lower/upper offset`. The isotope ladder is
/// anchored on `precursor_mz`. `charge` 0 is treated as 1. `iso_tol_ppm`
/// matches isotope peaks (reference default 10 ppm).
pub fn precursor_purity(
    ms1: &[(f64, f32)],
    precursor_mz: f64,
    target_mz: Option<f64>,
    charge: u8,
    lower_off: f64,
    upper_off: f64,
    iso_tol_ppm: f64,
) -> Option<f32> {
    let centre = target_mz.unwrap_or(precursor_mz);
    let lo = centre - lower_off;
    let hi = centre + upper_off;
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
    let in_window = |mz: f64| (lo..=hi).contains(&mz);

    // Anchor: most intense MS1 peak within tolerance of the precursor m/z.
    // Without one the precursor contributes nothing (purity 0), not "unknown".
    let Some((anchor_mz, anchor_int)) = best_within(ms1, precursor_mz, tol_da(precursor_mz)) else {
        return Some(0.0);
    };
    let mut envelope = if in_window(anchor_mz) {
        anchor_int as f64
    } else {
        0.0
    };
    for (dir, edge) in [(-1.0f64, lo), (1.0, hi)] {
        let mut expected = anchor_mz + dir * spacing;
        let (mut inside, mut steps) = (0usize, 0usize);
        // Walk while the next rung has not passed the window edge on this side.
        while dir * (edge - expected) >= 0.0
            && inside < MAX_WINDOW_STEPS
            && steps < MAX_LADDER_STEPS
        {
            let was_inside = in_window(expected);
            if let Some((mz, inten)) = best_within(ms1, expected, tol_da(expected)) {
                if in_window(mz) {
                    envelope += inten as f64;
                }
                expected = mz + dir * spacing;
            } else {
                expected += dir * spacing;
            }
            inside += usize::from(was_inside);
            steps += 1;
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
        let p = precursor_purity(&ms1, mz, None, 2, 0.7, 0.7, 10.0).unwrap();
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
        let p = precursor_purity(&ms1, mz, None, 2, 0.7, 0.7, 10.0).unwrap();
        assert!((p - 0.6).abs() < 1e-6, "{p}");
    }

    #[test]
    fn missing_precursor_peak_is_zero_and_empty_window_is_none() {
        let ms1 = vec![(400.0, 10.0f32), (600.0, 10.0)];
        assert_eq!(precursor_purity(&ms1, 500.0, None, 2, 0.7, 0.7, 10.0), None);
        let ms1 = vec![(500.3, 10.0f32)];
        assert_eq!(
            precursor_purity(&ms1, 500.0, None, 2, 0.7, 0.7, 10.0),
            Some(0.0)
        );
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
        let p = precursor_purity(&ms1, mz, None, 2, 1.5, 1.5, 10.0).unwrap();
        assert!((p - 0.75).abs() < 1e-6, "{p}");
    }

    /// The window is centred on the isolation target, not on the monoisotopic
    /// precursor: the instrument isolated M+1 (650.83) of a z=2 precursor whose
    /// monoisotopic m/z is 650.33. Inside `[650.48, 651.18]` sit M+1 and an
    /// interferent at 650.95, so purity is 1.1 / 2.1.
    #[test]
    fn window_is_centred_on_the_isolation_target() {
        let ms1 = vec![
            (650.33, 1.0e6f32),
            (650.83, 1.1e6),
            (650.95, 1.0e6),
            (651.33, 0.6e6),
        ];
        let p = precursor_purity(&ms1, 650.33, Some(650.83), 2, 0.35, 0.35, 10.0).unwrap();
        assert!((p - 1.1 / 2.1).abs() < 1e-3, "{p}");
        // Without a recorded target the window falls back to the precursor
        // m/z, which here excludes the interferent entirely.
        let p = precursor_purity(&ms1, 650.33, None, 2, 0.35, 0.35, 10.0).unwrap();
        assert!((p - 1.0).abs() < 1e-6, "{p}");
    }

    /// A target equal to the precursor m/z gives the precursor-centred result.
    #[test]
    fn target_at_the_precursor_matches_the_fallback() {
        let mz = 500.25;
        let ms1 = vec![
            (mz - 0.3, 1000.0f32),
            (mz, 1000.0),
            (mz + ISOTOPE / 2.0, 500.0),
        ];
        let with = precursor_purity(&ms1, mz, Some(mz), 2, 0.7, 0.7, 10.0);
        let without = precursor_purity(&ms1, mz, None, 2, 0.7, 0.7, 10.0);
        assert_eq!(with, without);
    }
}
