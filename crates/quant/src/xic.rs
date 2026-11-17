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

pub fn gaussian_smooth(rt: &[f64], x: &[f32], width_s: f64) -> Vec<f32> {
    assert_eq!(rt.len(), x.len());
    if x.len() < 2 || !width_s.is_finite() || width_s <= 0.0 {
        return x.to_vec();
    }
    let sigma = width_s / 8.0;
    let radius = width_s / 2.0;
    let weights: Vec<f64> = (0..rt.len())
        .map(|i| {
            let left = if i == 0 { rt[i] } else { rt[i - 1] };
            let right = if i + 1 == rt.len() { rt[i] } else { rt[i + 1] };
            ((right - left) / 2.0).max(0.0)
        })
        .collect();
    let mut left = 0;
    let mut right = 0;
    rt.iter()
        .enumerate()
        .map(|(i, &center)| {
            while left < rt.len() && rt[left] < center - radius {
                left += 1;
            }
            while right < rt.len() && rt[right] <= center + radius {
                right += 1;
            }
            let mut total = 0.0;
            let mut signal = 0.0;
            for j in left..right {
                let d = (rt[j] - center) / sigma;
                let weight = (-0.5 * d * d).exp() * weights[j];
                total += weight;
                signal += weight * x[j] as f64;
            }
            if total > 0.0 {
                (signal / total) as f32
            } else {
                x[i]
            }
        })
        .collect()
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PeakBounds {
    pub apex: usize,
    pub left: usize,
    pub right: usize,
}

const HALF_WIDTHS: f32 = 2.0;
const DEEP_VALLEY: f32 = 0.5;
const VALLEY_RISE: f32 = 1.5;

pub fn find_peak(trace: &[f32], anchor: usize, max_climb: usize) -> Option<PeakBounds> {
    let n = trace.len();
    if n == 0 || anchor >= n {
        return None;
    }
    let climb = |dir: isize| -> usize {
        let mut i = anchor;
        let mut steps = 0usize;
        loop {
            let next = i as isize + dir;
            if next < 0 || next >= n as isize || steps >= max_climb {
                break;
            }
            let next = next as usize;
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
    let extend = |dir: isize| -> usize {
        let limit = (HALF_WIDTHS * half_width(trace, apex, dir)).round() as usize;
        let mut i = apex;
        let mut running_min = apex_int;
        for _ in 0..limit {
            let next = i as isize + dir;
            if next < 0 || next >= n as isize {
                break;
            }
            let v = trace[next as usize];
            if running_min < DEEP_VALLEY * apex_int && v > VALLEY_RISE * running_min {
                break;
            }
            running_min = running_min.min(v);
            i = next as usize;
        }
        i
    };
    let left = extend(-1);
    let right = extend(1);
    Some(PeakBounds { apex, left, right })
}

fn half_width(trace: &[f32], apex: usize, dir: isize) -> f32 {
    let half = trace[apex] / 2.0;
    let mut i = apex;
    loop {
        let next = i as isize + dir;
        if next < 0 || next >= trace.len() as isize {
            return i.abs_diff(apex) as f32;
        }
        let next = next as usize;
        if trace[next] <= half {
            let drop = (trace[i] - trace[next]).max(f32::EPSILON);
            return i.abs_diff(apex) as f32 + (trace[i] - half) / drop;
        }
        i = next;
    }
}

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
    fn gaussian_smoothing_weights_irregular_samples_in_time() {
        let rt = [0.0, 0.2, 0.4, 2.0, 4.0, 7.0, 10.0];
        let flat = gaussian_smooth(&rt, &[5.0; 7], 8.0);
        assert!(flat.iter().all(|v| (*v - 5.0).abs() < 1e-6));
        let spike = gaussian_smooth(&rt, &[0.0, 0.0, 0.0, 100.0, 0.0, 0.0, 0.0], 8.0);
        assert!(spike[3] > spike[2] && spike[3] > spike[4]);
        assert_eq!(spike[6], 0.0);
        assert_eq!(gaussian_smooth(&[1.0], &[9.0], 8.0), vec![9.0]);
        let rt: Vec<f64> = (0..100)
            .map(|i| i as f64 * 0.05)
            .chain((0..=10).map(|i| 5.0 + i as f64 * 0.5))
            .collect();
        let ramp: Vec<f32> = rt.iter().map(|&t| t as f32).collect();
        assert!((gaussian_smooth(&rt, &ramp, 8.0)[100] - 5.0).abs() < 0.02);
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
        assert_eq!((p.left, p.right), (13, 27), "{p:?}");
    }

    #[test]
    fn boundaries_ignore_noise_on_a_broad_peak() {
        let t: Vec<f32> = gaussian(121, 60.0, 12.0, 1000.0)
            .into_iter()
            .enumerate()
            .map(|(i, v)| v * if i % 2 == 0 { 1.15 } else { 0.85 })
            .collect();
        let s = savitzky_golay5(&t);
        let p = find_peak(&s, 58, 12).unwrap();
        assert!((p.apex as i64 - 60).abs() <= 2, "{p:?}");
        assert!(p.right - p.left >= 50, "{p:?}");
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
