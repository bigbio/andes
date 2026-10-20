//! Per-run retention-time calibration: predicted index -> observed minutes.
//!
//! An isotonic (non-decreasing) fit by pool-adjacent-violators on anchor pairs
//! (predicted, observed), applied to every row by linear interpolation between
//! the fitted knots and clamped at the ends. With fewer than
//! [`MIN_ISOTONIC_ANCHORS`] anchors a least-squares line is used instead.

/// Anchors below which the calibration is a least-squares line.
pub const MIN_ISOTONIC_ANCHORS: usize = 20;

/// A fitted calibration curve.
#[derive(Clone, Debug, PartialEq)]
pub enum RtCalibration {
    /// `y = slope * x + intercept`.
    Linear {
        /// Slope.
        slope: f64,
        /// Intercept.
        intercept: f64,
    },
    /// Knots `(x, y)`, `x` strictly increasing and `y` non-decreasing.
    Isotonic {
        /// Knot positions.
        x: Vec<f64>,
        /// Knot values.
        y: Vec<f64>,
    },
}

impl RtCalibration {
    /// Fit on `(predicted, observed)` anchors; non-finite pairs are ignored.
    /// `None` with fewer than two usable anchors or a degenerate line fit.
    pub fn fit(anchors: &[(f64, f64)]) -> Option<Self> {
        let mut a: Vec<(f64, f64)> = anchors
            .iter()
            .copied()
            .filter(|(x, y)| x.is_finite() && y.is_finite())
            .collect();
        if a.len() < 2 {
            return None;
        }
        if a.len() < MIN_ISOTONIC_ANCHORS {
            let n = a.len() as f64;
            let mx = a.iter().map(|p| p.0).sum::<f64>() / n;
            let my = a.iter().map(|p| p.1).sum::<f64>() / n;
            let sxx: f64 = a.iter().map(|p| (p.0 - mx) * (p.0 - mx)).sum();
            let sxy: f64 = a.iter().map(|p| (p.0 - mx) * (p.1 - my)).sum();
            if sxx <= 0.0 {
                return None;
            }
            let slope = sxy / sxx;
            return Some(RtCalibration::Linear {
                slope,
                intercept: my - slope * mx,
            });
        }
        a.sort_by(|p, q| p.0.total_cmp(&q.0));
        // Tied x collapse to their mean y, weighted by count.
        let mut ux: Vec<f64> = Vec::new();
        let mut uy: Vec<f64> = Vec::new();
        let mut w: Vec<f64> = Vec::new();
        for (x, y) in a {
            if ux.last() == Some(&x) {
                let k = ux.len() - 1;
                uy[k] = (uy[k] * w[k] + y) / (w[k] + 1.0);
                w[k] += 1.0;
            } else {
                ux.push(x);
                uy.push(y);
                w.push(1.0);
            }
        }
        let y = pav(&uy, &w);
        Some(RtCalibration::Isotonic { x: ux, y })
    }

    /// Calibrated value at `x`; non-finite `x` gives NaN.
    pub fn apply(&self, x: f64) -> f64 {
        if !x.is_finite() {
            return f64::NAN;
        }
        match self {
            RtCalibration::Linear { slope, intercept } => slope * x + intercept,
            RtCalibration::Isotonic { x: kx, y: ky } => {
                if x <= kx[0] {
                    return ky[0];
                }
                let last = kx.len() - 1;
                if x >= kx[last] {
                    return ky[last];
                }
                let j = kx.partition_point(|v| *v <= x);
                let (x0, x1, y0, y1) = (kx[j - 1], kx[j], ky[j - 1], ky[j]);
                y0 + (y1 - y0) * (x - x0) / (x1 - x0)
            }
        }
    }
}

/// Weighted non-decreasing least-squares fit of `y` (ordered by x).
fn pav(y: &[f64], w: &[f64]) -> Vec<f64> {
    let mut val: Vec<f64> = Vec::with_capacity(y.len());
    let mut wt: Vec<f64> = Vec::with_capacity(y.len());
    let mut cnt: Vec<usize> = Vec::with_capacity(y.len());
    for (&yi, &wi) in y.iter().zip(w) {
        val.push(yi);
        wt.push(wi);
        cnt.push(1);
        while val.len() > 1 && val[val.len() - 2] > val[val.len() - 1] {
            let (v1, w1, c1) = (val.pop().unwrap(), wt.pop().unwrap(), cnt.pop().unwrap());
            let k = val.len() - 1;
            let tw = wt[k] + w1;
            val[k] = (val[k] * wt[k] + v1 * w1) / tw;
            wt[k] = tw;
            cnt[k] += c1;
        }
    }
    val.iter()
        .zip(&cnt)
        .flat_map(|(v, c)| std::iter::repeat_n(*v, *c))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn isotonic_is_monotone_and_robust_to_an_outlier_block() {
        // y = 10 + 100 x with a deterministic wobble and one off-trend block.
        let n = 2000;
        let mut a: Vec<(f64, f64)> = (0..n)
            .map(|i| {
                let x = i as f64 / n as f64;
                let wobble = ((i * 7919) % 13) as f64 / 13.0 - 0.5;
                (x, 10.0 + 100.0 * x + 2.0 * wobble)
            })
            .collect();
        for p in a.iter_mut().take(640).skip(600) {
            p.1 = 90.0;
        }
        let c = RtCalibration::fit(&a).unwrap();
        assert!(matches!(c, RtCalibration::Isotonic { .. }));
        let mut prev = f64::NEG_INFINITY;
        for i in 0..=200 {
            let v = c.apply(i as f64 / 200.0);
            assert!(v >= prev - 1e-9);
            prev = v;
        }
        let mut err: Vec<f64> = a
            .iter()
            .map(|(x, _)| (c.apply(*x) - (10.0 + 100.0 * x)).abs())
            .collect();
        err.sort_by(f64::total_cmp);
        assert!(
            err[err.len() / 2] < 2.0,
            "median error {}",
            err[err.len() / 2]
        );
    }

    #[test]
    fn few_anchors_fit_a_line() {
        let a = [(0.0, 1.0), (1.0, 3.0), (2.0, 5.0)];
        let c = RtCalibration::fit(&a).unwrap();
        assert!(matches!(c, RtCalibration::Linear { .. }));
        assert!((c.apply(3.0) - 7.0).abs() < 1e-12);
    }

    #[test]
    fn degenerate_inputs_give_none() {
        assert!(RtCalibration::fit(&[]).is_none());
        assert!(RtCalibration::fit(&[(1.0, 2.0)]).is_none());
        assert!(RtCalibration::fit(&[(1.0, 2.0), (1.0, 3.0)]).is_none());
        assert!(RtCalibration::fit(&[(f64::NAN, 2.0), (1.0, f64::INFINITY)]).is_none());
    }

    #[test]
    fn isotonic_clamps_outside_the_anchors_and_ties_average() {
        let mut a: Vec<(f64, f64)> = (0..30).map(|i| (i as f64, i as f64)).collect();
        a.push((5.0, 7.0));
        let c = RtCalibration::fit(&a).unwrap();
        assert_eq!(c.apply(-10.0), 0.0);
        assert_eq!(c.apply(100.0), 29.0);
        assert!((c.apply(5.0) - 6.0).abs() < 1e-12);
        assert!(c.apply(f64::NAN).is_nan());
    }
}
