use model::mass::ISOTOPE;

const MAX_WINDOW_STEPS: usize = 10;
const MAX_LADDER_STEPS: usize = 64;

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
            (mz - 0.3, 1000.0f32),
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
        let ms1 = vec![
            (mz, 100.0f32),
            (mz + 2.0 * ISOTOPE / 2.0, 50.0),
            (mz + 0.9, 50.0),
        ];
        let p = precursor_purity(&ms1, mz, None, 2, 1.5, 1.5, 10.0).unwrap();
        assert!((p - 0.75).abs() < 1e-6, "{p}");
    }

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
        let p = precursor_purity(&ms1, 650.33, None, 2, 0.35, 0.35, 10.0).unwrap();
        assert!((p - 1.0).abs() < 1e-6, "{p}");
    }

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
