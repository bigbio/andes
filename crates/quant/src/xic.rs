//! Chromatographic peak handling on extracted-ion chromatograms: smoothing,
//! apex and boundary detection, integration, envelope similarity.

/// Five-point quadratic Savitzky–Golay smoothing (`[-3, 12, 17, 12, -3] / 35`).
/// The two values at either end are copied through; negative results are
/// clamped to zero (intensities).
pub fn savitzky_golay5(x: &[f32]) -> Vec<f32> {
    let n = x.len();
    if n < 5 {
        return x.to_vec();
    }
    let mut out = x.to_vec();
    for i in 2..n - 2 {
        let v = (-3.0 * x[i - 2] + 12.0 * x[i - 1] + 17.0 * x[i] + 12.0 * x[i + 1]
            - 3.0 * x[i + 2])
            / 35.0;
        out[i] = v.max(0.0);
    }
    out
}

/// A chromatographic peak on a trace: apex scan and inclusive boundaries.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PeakBounds {
    pub apex: usize,
    pub left: usize,
    pub right: usize,
}

/// Fraction of the apex below which a boundary scan is reached.
const BOUNDARY_FRACTION: f32 = 0.02;
/// A rise above the running minimum by this factor ends the peak (valley).
const VALLEY_RISE: f32 = 1.10;

/// Locate the peak that contains `anchor`: climb from the anchor to the
/// nearest local maximum (at most `max_climb` scans away), then extend the
/// boundaries outward while the trace keeps falling, stopping at a valley
/// (the trace rises again above the running minimum by [`VALLEY_RISE`]) or
/// at [`BOUNDARY_FRACTION`] of the apex. One missing scan (zero) inside a
/// peak is tolerated. `None` when the trace has no signal at the anchor's
/// peak.
pub fn find_peak(trace: &[f32], anchor: usize, max_climb: usize) -> Option<PeakBounds> {
    let n = trace.len();
    if n == 0 || anchor >= n {
        return None;
    }
    // Climb both ways; take the higher local maximum.
    let climb = |dir: isize| -> usize {
        let mut i = anchor;
        let mut steps = 0usize;
        loop {
            let next = i as isize + dir;
            if next < 0 || next >= n as isize || steps >= max_climb {
                break;
            }
            let next = next as usize;
            // Allow stepping through a single zero scan to reach a higher value.
            let peek = (next as isize + dir).clamp(0, n as isize - 1) as usize;
            if trace[next] > trace[i]
                || (trace[next] == 0.0 && trace[peek] > trace[i] && peek != next)
            {
                i = next;
                steps += 1;
            } else {
                break;
            }
        }
        // If we stopped on a zero bridge, step back to the real maximum.
        if trace[i] == 0.0 {
            i = anchor;
        }
        i
    };
    let l = climb(-1);
    let r = climb(1);
    let apex = if trace[r] >= trace[l] { r } else { l };
    let apex_int = trace[apex];
    if apex_int <= 0.0 {
        return None;
    }
    let floor = apex_int * BOUNDARY_FRACTION;

    let extend = |dir: isize| -> usize {
        let mut i = apex;
        let mut running_min = apex_int;
        let mut gap = false;
        loop {
            let next = i as isize + dir;
            if next < 0 || next >= n as isize {
                break;
            }
            let next = next as usize;
            let v = trace[next];
            if v <= floor {
                // Include the first below-floor scan as the boundary unless it
                // is a one-scan gap followed by real signal.
                let peek = next as isize + dir;
                if !gap
                    && peek >= 0
                    && peek < n as isize
                    && trace[peek as usize] > floor
                    && trace[peek as usize] <= running_min
                {
                    gap = true;
                    i = next;
                    continue;
                }
                i = next;
                break;
            }
            if v > running_min * VALLEY_RISE {
                break; // rising again: the valley was the previous scan
            }
            running_min = running_min.min(v);
            i = next;
        }
        i
    };
    let left = extend(-1);
    let right = extend(1);
    Some(PeakBounds { apex, left, right })
}

/// Trapezoidal area of `trace` between `left..=right` over `rt` (seconds).
/// A single-scan peak returns its height (as if one second wide).
pub fn integrate(rt: &[f64], trace: &[f32], left: usize, right: usize) -> f64 {
    if trace.is_empty() || left > right || right >= trace.len() {
        return 0.0;
    }
    if left == right {
        return trace[left] as f64;
    }
    let mut area = 0.0;
    for i in left..right {
        let dt = (rt[i + 1] - rt[i]).max(0.0);
        area += 0.5 * (trace[i] as f64 + trace[i + 1] as f64) * dt;
    }
    area
}

/// Cosine similarity of two non-negative vectors (0 when either is empty or
/// all zero).
pub fn cosine(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len().min(b.len());
    if n == 0 {
        return 0.0;
    }
    let (mut dot, mut na, mut nb) = (0.0, 0.0, 0.0);
    for i in 0..n {
        dot += a[i] * b[i];
        na += a[i] * a[i];
        nb += b[i] * b[i];
    }
    if na <= 0.0 || nb <= 0.0 {
        0.0
    } else {
        (dot / (na.sqrt() * nb.sqrt())).clamp(0.0, 1.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn gaussian(n: usize, center: f64, sigma: f64, height: f32) -> Vec<f32> {
        (0..n)
            .map(|i| {
                let d = (i as f64 - center) / sigma;
                (height as f64 * (-0.5 * d * d).exp()) as f32
            })
            .collect()
    }

    #[test]
    fn smoothing_preserves_a_flat_trace_and_clamps() {
        let flat = vec![5.0f32; 9];
        assert_eq!(savitzky_golay5(&flat), flat);
        let spike = vec![0.0, 0.0, 100.0, 0.0, 0.0, 0.0, 0.0];
        let s = savitzky_golay5(&spike);
        assert!(s.iter().all(|&v| v >= 0.0));
        assert!(s[2] > s[3]);
    }

    #[test]
    fn finds_apex_and_bounds_of_a_gaussian() {
        let t = gaussian(41, 20.0, 3.0, 1000.0);
        let p = find_peak(&t, 17, 10).unwrap();
        assert_eq!(p.apex, 20);
        // 2 % of the apex is ~2.8 sigma away → bounds near 11 / 29.
        assert!(p.left >= 10 && p.left <= 12, "{p:?}");
        assert!(p.right >= 28 && p.right <= 30, "{p:?}");
    }

    #[test]
    fn stops_at_a_valley_between_two_peaks() {
        let mut t = gaussian(61, 20.0, 3.0, 1000.0);
        let second = gaussian(61, 34.0, 3.0, 800.0);
        for (a, b) in t.iter_mut().zip(second) {
            *a += b;
        }
        let p = find_peak(&t, 21, 10).unwrap();
        assert_eq!(p.apex, 20);
        assert!(p.right >= 25 && p.right <= 28, "{p:?}");
    }

    #[test]
    fn empty_or_zero_trace_is_none() {
        assert_eq!(find_peak(&[], 0, 5), None);
        assert_eq!(find_peak(&[0.0, 0.0, 0.0], 1, 5), None);
        assert_eq!(find_peak(&[1.0], 3, 5), None);
    }

    #[test]
    fn integrates_trapezoids() {
        let rt = [0.0, 1.0, 2.0, 3.0];
        let t = [0.0, 10.0, 10.0, 0.0];
        assert!((integrate(&rt, &t, 0, 3) - 20.0).abs() < 1e-9);
        assert_eq!(integrate(&rt, &t, 1, 1), 10.0);
        assert_eq!(integrate(&rt, &t, 2, 1), 0.0);
    }

    #[test]
    fn cosine_of_proportional_vectors_is_one() {
        assert!((cosine(&[1.0, 2.0, 3.0], &[2.0, 4.0, 6.0]) - 1.0).abs() < 1e-12);
        assert_eq!(cosine(&[1.0, 0.0], &[0.0, 1.0]), 0.0);
        assert_eq!(cosine(&[], &[1.0]), 0.0);
    }
}
