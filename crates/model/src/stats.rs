//! Median, median absolute deviation and quantiles shared by the quantification, input,
//! search and output crates.
//!
//! Every function orders values with `total_cmp`, so a NaN sorts last and the result does not
//! depend on the algorithm used to find the middle element. Callers choose the even-length rule
//! explicitly with [`Even`] and keep their own convention for an empty slice through the
//! returned `Option`.

/// Factor that turns a median absolute deviation into a standard deviation for normal data.
pub const MAD_TO_SD: f64 = 1.4826;

/// `x`, or negative infinity for a NaN, so that a NaN score ranks last when sorting by score.
pub fn nan_low(x: f64) -> f64 {
    if x.is_nan() {
        f64::NEG_INFINITY
    } else {
        x
    }
}

/// [`nan_low`] for `f32`.
pub fn nan_low_f32(x: f32) -> f32 {
    if x.is_nan() {
        f32::NEG_INFINITY
    } else {
        x
    }
}

/// How the median of an even number of values is taken.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Even {
    /// `0.5 * (lower middle + upper middle)`.
    Mean,
    /// The upper of the two middle values, `sorted[n / 2]`.
    Upper,
}

/// Median of `v` in `f64`; `v` is reordered. `None` when `v` is empty.
pub fn median_mut(v: &mut [f64], even: Even) -> Option<f64> {
    let n = v.len();
    if n == 0 {
        return None;
    }
    let (left, &mut upper, _) = v.select_nth_unstable_by(n / 2, f64::total_cmp);
    if n % 2 == 1 || even == Even::Upper {
        return Some(upper);
    }
    let lower = left.iter().copied().max_by(f64::total_cmp).unwrap_or(upper);
    Some(0.5 * (lower + upper))
}

/// Median of `v` in `f32` arithmetic; `v` is reordered. `None` when `v` is empty.
pub fn median_f32_mut(v: &mut [f32], even: Even) -> Option<f32> {
    let n = v.len();
    if n == 0 {
        return None;
    }
    let (left, &mut upper, _) = v.select_nth_unstable_by(n / 2, f32::total_cmp);
    if n % 2 == 1 || even == Even::Upper {
        return Some(upper);
    }
    let lower = left.iter().copied().max_by(f32::total_cmp).unwrap_or(upper);
    Some(0.5 * (lower + upper))
}

/// Median of `f32` values with the two middle values widened to `f64` before they are averaged
/// (`Even::Mean`); `v` is reordered. `None` when `v` is empty.
pub fn median_f32_as_f64_mut(v: &mut [f32]) -> Option<f64> {
    let n = v.len();
    if n == 0 {
        return None;
    }
    let (left, &mut upper, _) = v.select_nth_unstable_by(n / 2, f32::total_cmp);
    if n % 2 == 1 {
        return Some(upper as f64);
    }
    let lower = left.iter().copied().max_by(f32::total_cmp).unwrap_or(upper);
    Some(0.5 * (lower as f64 + upper as f64))
}

/// The finite values of `v`, in input order.
pub fn finite_vec(v: &[f64]) -> Vec<f64> {
    v.iter().copied().filter(|x| x.is_finite()).collect()
}

/// Median of `v` and the raw median absolute deviation around that median.
///
/// `v` is overwritten with the absolute deviations. Multiply the second value by [`MAD_TO_SD`]
/// for a standard deviation. `None` when `v` is empty.
pub fn median_mad_mut(v: &mut [f64], even: Even) -> Option<(f64, f64)> {
    let m = median_mut(v, even)?;
    for x in v.iter_mut() {
        *x = (*x - m).abs();
    }
    let mad = median_mut(v, even)?;
    Some((m, mad))
}

/// [`median_mad_mut`] on a copy of `v`.
pub fn median_mad(v: &[f64], even: Even) -> Option<(f64, f64)> {
    median_mad_mut(&mut v.to_vec(), even)
}

/// Quantile `q` of the ascending `sorted` values by linear interpolation between the two
/// neighbouring ranks; NaN when `sorted` is empty.
///
/// With `clamp`, `q` is first limited to `[0, 1]`. Without it, `q` must already lie in `[0, 1]`.
/// At `q = 0.5` and an even length this is `a + (b - a) * 0.5`, which is not always bit-equal to
/// the [`Even::Mean`] median `0.5 * (a + b)`.
pub fn quantile_sorted(sorted: &[f64], q: f64, clamp: bool) -> f64 {
    if sorted.is_empty() {
        return f64::NAN;
    }
    let q = if clamp { q.clamp(0.0, 1.0) } else { q };
    let pos = q * (sorted.len() - 1) as f64;
    let lo = pos.floor() as usize;
    let hi = pos.ceil() as usize;
    sorted[lo] + (sorted[hi] - sorted[lo]) * (pos - lo as f64)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nan_scores_rank_below_every_number() {
        assert_eq!(nan_low(f64::NAN), f64::NEG_INFINITY);
        assert_eq!(nan_low(-3.5), -3.5);
        assert_eq!(nan_low_f32(f32::NAN), f32::NEG_INFINITY);
        assert_eq!(nan_low_f32(2.0), 2.0);
    }

    #[test]
    fn median_odd_even_and_empty() {
        assert_eq!(median_mut(&mut [], Even::Mean), None);
        assert_eq!(median_mut(&mut [], Even::Upper), None);
        assert_eq!(median_mut(&mut [3.0, 1.0, 2.0], Even::Mean), Some(2.0));
        assert_eq!(median_mut(&mut [3.0, 1.0, 2.0], Even::Upper), Some(2.0));
        assert_eq!(median_mut(&mut [4.0, 1.0, 3.0, 2.0], Even::Mean), Some(2.5));
        assert_eq!(
            median_mut(&mut [4.0, 1.0, 3.0, 2.0], Even::Upper),
            Some(3.0)
        );
        assert_eq!(median_mut(&mut [7.5], Even::Mean), Some(7.5));
        assert_eq!(
            median_mut(&mut [1.0, 2.0, 3.0, 4.0, 1000.0], Even::Mean),
            Some(3.0)
        );
    }

    #[test]
    fn median_agrees_with_a_full_sort() {
        let mut state = 0x9E37_79B9_7F4A_7C15u64;
        for n in 1..40usize {
            let v: Vec<f64> = (0..n)
                .map(|_| {
                    state = state
                        .wrapping_mul(6364136223846793005)
                        .wrapping_add(1442695040888963407);
                    ((state >> 33) % 17) as f64 / 3.0 - 2.0
                })
                .collect();
            let mut sorted = v.clone();
            sorted.sort_by(f64::total_cmp);
            let upper = sorted[n / 2];
            let mean = if n % 2 == 1 {
                sorted[n / 2]
            } else {
                0.5 * (sorted[n / 2 - 1] + sorted[n / 2])
            };
            assert_eq!(median_mut(&mut v.clone(), Even::Mean), Some(mean));
            assert_eq!(median_mut(&mut v.clone(), Even::Upper), Some(upper));
        }
    }

    #[test]
    fn nan_sorts_last_and_negative_zero_before_positive_zero() {
        assert_eq!(median_mut(&mut [f64::NAN, 1.0, 2.0], Even::Mean), Some(2.0));
        assert!(median_mut(&mut [f64::NAN, f64::NAN, 1.0], Even::Mean)
            .unwrap()
            .is_nan());
        let m = median_mut(&mut [0.0, -0.0, 5.0], Even::Mean).unwrap();
        assert_eq!(m.to_bits(), 0.0f64.to_bits());
        let m = median_mut(&mut [0.0, -0.0], Even::Upper).unwrap();
        assert_eq!(m.to_bits(), 0.0f64.to_bits());
    }

    #[test]
    fn f32_median_uses_f32_arithmetic() {
        assert_eq!(median_f32_mut(&mut [], Even::Mean), None);
        assert_eq!(median_f32_mut(&mut [3.0, 1.0, 2.0], Even::Mean), Some(2.0));
        assert_eq!(
            median_f32_mut(&mut [4.0, 1.0, 3.0, 2.0], Even::Mean),
            Some(2.5)
        );
        assert_eq!(
            median_f32_mut(&mut [4.0, 1.0, 3.0, 2.0], Even::Upper),
            Some(3.0)
        );
        let (a, b) = (16_777_216.0f32, 16_777_218.0f32);
        let in_f32 = 0.5 * (a + b);
        assert_eq!(median_f32_mut(&mut [a, b], Even::Mean), Some(in_f32));
        assert_ne!(in_f32 as f64, 0.5 * (a as f64 + b as f64));
    }

    #[test]
    fn f32_median_widened_to_f64() {
        assert_eq!(median_f32_as_f64_mut(&mut []), None);
        assert_eq!(median_f32_as_f64_mut(&mut [3.0, 1.0, 2.0]), Some(2.0));
        assert_eq!(median_f32_as_f64_mut(&mut [4.0, 1.0, 3.0, 2.0]), Some(2.5));
        let (a, b) = (16_777_216.0f32, 16_777_218.0f32);
        assert_eq!(median_f32_as_f64_mut(&mut [a, b]), Some(16_777_217.0));
        let (a, b) = (1.0f32, 1.0f32 + f32::EPSILON);
        let wide = 0.5 * (a as f64 + b as f64);
        assert_eq!(median_f32_as_f64_mut(&mut [a, b]), Some(wide));
    }

    #[test]
    fn finite_vec_drops_non_finite_values() {
        let v = [1.0, f64::NAN, f64::INFINITY, -2.0, f64::NEG_INFINITY];
        assert_eq!(finite_vec(&v), vec![1.0, -2.0]);
    }

    #[test]
    fn median_and_mad_around_the_median() {
        assert_eq!(median_mad(&[], Even::Mean), None);
        let (m, mad) = median_mad(&[9.0, 10.0, 11.0], Even::Mean).unwrap();
        assert_eq!((m, mad), (10.0, 1.0));
        assert!((MAD_TO_SD * mad - 1.4826).abs() < 1e-12);
        let (m, mad) = median_mad(&[1.0, 2.0, 3.0, 4.0], Even::Mean).unwrap();
        assert_eq!((m, mad), (2.5, 1.0));
        let (m, mad) = median_mad(&[1.0, 2.0, 3.0, 4.0], Even::Upper).unwrap();
        assert_eq!((m, mad), (3.0, 1.0));
    }

    #[test]
    fn median_mad_mut_leaves_the_deviations() {
        let mut v = vec![1.0, 2.0, 4.0];
        let (m, mad) = median_mad_mut(&mut v, Even::Mean).unwrap();
        assert_eq!((m, mad), (2.0, 1.0));
        v.sort_by(f64::total_cmp);
        assert_eq!(v, vec![0.0, 1.0, 2.0]);
    }

    #[test]
    fn quantile_interpolates_between_ranks() {
        let s = [1.0, 2.0, 3.0, 4.0, 5.0];
        assert!(quantile_sorted(&[], 0.5, true).is_nan());
        assert_eq!(quantile_sorted(&[7.0], 0.9, true), 7.0);
        assert_eq!(quantile_sorted(&s, 0.0, true), 1.0);
        assert_eq!(quantile_sorted(&s, 0.5, true), 3.0);
        assert_eq!(quantile_sorted(&s, 1.0, true), 5.0);
        assert_eq!(quantile_sorted(&s, 0.25, true), 2.0);
        assert_eq!(quantile_sorted(&s, 0.125, true), 1.5);
        assert_eq!(quantile_sorted(&[1.0, 2.0, 3.0, 4.0], 0.5, false), 2.5);
    }

    #[test]
    fn quantile_clamps_only_when_asked() {
        let s = [1.0, 2.0, 3.0];
        assert_eq!(quantile_sorted(&s, -1.0, true), 1.0);
        assert_eq!(quantile_sorted(&s, 2.0, true), 3.0);
        assert_eq!(quantile_sorted(&s, 1.0, false), 3.0);
    }
}
