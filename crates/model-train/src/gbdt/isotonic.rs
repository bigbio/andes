pub fn pava(x: &[f32], y: &[f32]) -> (Vec<f32>, Vec<f32>) {
    assert_eq!(x.len(), y.len(), "pava: x and y must be equal length");
    if x.is_empty() {
        return (Vec::new(), Vec::new());
    }
    let mut means: Vec<f32> = Vec::new();
    let mut weights: Vec<f32> = Vec::new();
    let mut right_x: Vec<f32> = Vec::new();
    for i in 0..x.len() {
        means.push(y[i]);
        weights.push(1.0);
        right_x.push(x[i]);
        while means.len() > 1 && means[means.len() - 2] > means[means.len() - 1] {
            let n = means.len();
            let w = weights[n - 2] + weights[n - 1];
            let m = (means[n - 2] * weights[n - 2] + means[n - 1] * weights[n - 1]) / w;
            means[n - 2] = m;
            weights[n - 2] = w;
            right_x[n - 2] = right_x[n - 1];
            means.pop();
            weights.pop();
            right_x.pop();
        }
    }
    let mut bx = vec![x[0]];
    let mut by = vec![means[0]];
    for b in 0..means.len() {
        bx.push(right_x[b]);
        by.push(means[b]);
    }
    (bx, by)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pava_is_monotone_and_fits_steps() {
        let x = vec![0.1_f32, 0.2, 0.3, 0.4, 0.5];
        let y = vec![0.0_f32, 1.0, 0.0, 1.0, 1.0];
        let (bx, by) = pava(&x, &y);
        assert_eq!(bx.len(), by.len());
        for w in by.windows(2) {
            assert!(w[1] >= w[0] - 1e-6, "not monotone: {by:?}");
        }
        assert!((bx[0] - 0.1).abs() < 1e-6 && (*bx.last().unwrap() - 0.5).abs() < 1e-6);
    }

    #[test]
    fn pava_already_monotone_is_unchanged_mean() {
        let x = vec![0.0_f32, 0.5, 1.0];
        let y = vec![0.1_f32, 0.5, 0.9];
        let (_bx, by) = pava(&x, &y);
        assert!((by[by.len() - 1] - 0.9).abs() < 1e-6);
        for w in by.windows(2) {
            assert!(w[1] >= w[0] - 1e-6);
        }
    }

    #[test]
    fn pava_empty_is_empty() {
        let (bx, by) = pava(&[], &[]);
        assert!(bx.is_empty() && by.is_empty());
    }
}
