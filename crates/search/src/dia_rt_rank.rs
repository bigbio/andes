//! Retention-time term of the DIA candidate rank (`--dia` only).
//!
//! A DIA scan's precursor window holds about a thousand forms with three or more fragment
//! votes per charge, and the fragment-index cut keeps the best hundred by votes. Forms whose
//! predicted retention time is far from the scan's own retention time are rarely the eluting
//! precursor, so the cut ranks by
//! `votes - RT_PENALTY_WEIGHT * min(|scan RT - f(raw RT)| / sigma, RT_Z_CAP_SIGMAS)^2`,
//! with the vote intensity as tie-break. `raw RT` is the andes-rt prediction of the form's own
//! sequence (targets and decoys alike), `f` a monotone piecewise-linear map from raw
//! predictions to the run's minutes and `sigma` the robust spread of confident matches around
//! it; both are fitted per run by [`DiaRtPrior::fit`]. Beyond the end knots `f` continues
//! with the map's overall slope: a form predicted to elute before the first or after the last
//! anchor is calibrated outside the anchors' RT range, not pinned at its edge (pinned, every
//! such form would agree with the scans at the edge, and on a short run or an RT slice that is
//! about half of the forms).

use model::peptide::Peptide;
use model::stats::MAD_TO_SD;

/// Weight of the squared RT z-score against one fragment vote.
pub const RT_PENALTY_WEIGHT: f64 = 0.5;
/// RT z-scores are capped here, so a form far off its predicted RT loses at most
/// `RT_PENALTY_WEIGHT * RT_Z_CAP_SIGMAS^2` votes; a form without a prediction gets the cap.
pub const RT_Z_CAP_SIGMAS: f64 = 6.0;
/// Fewest calibration anchors for a prior; with fewer the cut keeps the plain vote rank.
pub const MIN_ANCHORS: usize = 50;
/// Most knots of the monotone map.
const MAX_KNOTS: usize = 40;
/// Fewest knots of the monotone map.
const MIN_KNOTS: usize = 6;
/// Anchors per knot when the anchors are too few for `MAX_KNOTS`.
const ANCHORS_PER_KNOT: usize = 25;
/// Anchors at which the monotone map reaches its full resolution (`MAX_KNOTS` knots); the
/// calibration sample grows until it holds this many or the run is used up.
pub const FULL_MAP_ANCHORS: usize = MAX_KNOTS * ANCHORS_PER_KNOT;
/// Rounds of residual trimming in the fit.
const TRIM_ROUNDS: usize = 5;
/// Residuals beyond this many sigmas leave the next fitting round.
const TRIM_SIGMAS: f64 = 3.0;
/// Predicts the raw andes-rt value of one peptidoform; NaN when the model cannot encode it.
pub type FormRtModel = dyn Fn(&Peptide) -> f32 + Send + Sync;

/// The run's map from raw RT predictions to minutes and the spread around it.
#[derive(Clone, Debug, PartialEq)]
pub struct DiaRtPrior {
    knots_raw: Vec<f64>,
    knots_min: Vec<f64>,
    sigma_min: f64,
    anchors: usize,
    kept: usize,
}

impl DiaRtPrior {
    /// Fits the prior on `(raw prediction, observed minutes)` anchors: per round, the kept
    /// anchors are split into equal-count bins of the raw prediction, each bin contributes the
    /// knot (median raw, median minutes), knot minutes are made non-decreasing, sigma is the
    /// scaled median absolute deviation of the kept residuals, and anchors beyond
    /// `TRIM_SIGMAS` sigmas leave the next round. `None` with fewer than [`MIN_ANCHORS`]
    /// finite anchors or a degenerate spread.
    pub fn fit(anchors: &[(f64, f64)]) -> Option<Self> {
        let pts: Vec<(f64, f64)> = anchors
            .iter()
            .copied()
            .filter(|(x, y)| x.is_finite() && y.is_finite())
            .collect();
        let n = pts.len();
        if n < MIN_ANCHORS {
            return None;
        }
        let n_bins = (n / ANCHORS_PER_KNOT).clamp(MIN_KNOTS, MAX_KNOTS);
        let mut keep = vec![true; n];
        let mut knots_raw: Vec<f64> = Vec::new();
        let mut knots_min: Vec<f64> = Vec::new();
        let mut sigma = f64::NAN;
        let mut kept = n;
        for _ in 0..TRIM_ROUNDS {
            let mut xs: Vec<f64> = pts
                .iter()
                .zip(&keep)
                .filter(|(_, k)| **k)
                .map(|(p, _)| p.0)
                .collect();
            if xs.len() < MIN_ANCHORS {
                return None;
            }
            xs.sort_by(f64::total_cmp);
            let edges: Vec<f64> = (0..=n_bins)
                .map(|j| quantile_sorted(&xs, j as f64 / n_bins as f64))
                .collect();
            let mut bins: Vec<(Vec<f64>, Vec<f64>)> = vec![(Vec::new(), Vec::new()); n_bins];
            for (&(x, y), _) in pts.iter().zip(&keep).filter(|(_, k)| **k) {
                let b = edges
                    .partition_point(|&e| e <= x)
                    .saturating_sub(1)
                    .min(n_bins - 1);
                bins[b].0.push(x);
                bins[b].1.push(y);
            }
            knots_raw.clear();
            knots_min.clear();
            for (bx, by) in bins.iter_mut().filter(|(bx, _)| !bx.is_empty()) {
                let (mx, my) = (median_in_place(bx), median_in_place(by));
                let my = knots_min.last().map_or(my, |&last: &f64| last.max(my));
                if knots_raw.last() == Some(&mx) {
                    *knots_min.last_mut().expect("a knot exists") = my;
                } else {
                    knots_raw.push(mx);
                    knots_min.push(my);
                }
            }
            let resid: Vec<f64> = pts
                .iter()
                .map(|&(x, y)| y - interp(&knots_raw, &knots_min, x))
                .collect();
            let mut r_kept: Vec<f64> = resid
                .iter()
                .zip(&keep)
                .filter(|(_, k)| **k)
                .map(|(r, _)| *r)
                .collect();
            let centre = median_in_place(&mut r_kept);
            let mut dev: Vec<f64> = r_kept.iter().map(|r| (r - centre).abs()).collect();
            sigma = MAD_TO_SD * median_in_place(&mut dev);
            if !(sigma.is_finite() && sigma > 0.0) {
                return None;
            }
            for (k, r) in keep.iter_mut().zip(&resid) {
                *k = r.abs() < TRIM_SIGMAS * sigma;
            }
            kept = keep.iter().filter(|k| **k).count();
        }
        (knots_raw.len() >= 2).then_some(Self {
            knots_raw,
            knots_min,
            sigma_min: sigma,
            anchors: n,
            kept,
        })
    }

    /// Calibrated RT in minutes of a raw prediction (NaN stays NaN; beyond the end knots the
    /// map continues with its overall slope).
    pub fn apply(&self, raw: f64) -> f64 {
        if raw.is_finite() {
            interp(&self.knots_raw, &self.knots_min, raw)
        } else {
            f64::NAN
        }
    }

    pub fn sigma_min(&self) -> f64 {
        self.sigma_min
    }

    pub fn knots(&self) -> usize {
        self.knots_raw.len()
    }

    /// Finite anchors the fit started from.
    pub fn anchors(&self) -> usize {
        self.anchors
    }

    /// Anchors inside `TRIM_SIGMAS` sigmas after the last round.
    pub fn kept(&self) -> usize {
        self.kept
    }
}

/// The RT term of the rank key (`<= 0`) for an RT difference of `diff_min` minutes at spread
/// `sigma_min`; a NaN difference gets the capped penalty.
#[inline]
pub fn rt_penalty(diff_min: f64, sigma_min: f64) -> f64 {
    let z = if diff_min.is_finite() {
        (diff_min.abs() / sigma_min).min(RT_Z_CAP_SIGMAS)
    } else {
        RT_Z_CAP_SIGMAS
    };
    -RT_PENALTY_WEIGHT * z * z
}

fn quantile_sorted(xs: &[f64], q: f64) -> f64 {
    let pos = q * (xs.len() - 1) as f64;
    let lo = pos.floor() as usize;
    let hi = (lo + 1).min(xs.len() - 1);
    xs[lo] + (pos - lo as f64) * (xs[hi] - xs[lo])
}

fn median_in_place(v: &mut [f64]) -> f64 {
    v.sort_by(f64::total_cmp);
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2]
    } else {
        0.5 * (v[n / 2 - 1] + v[n / 2])
    }
}

/// Piecewise-linear interpolation through strictly increasing `xs`, extended beyond the ends
/// with the overall slope `(y_last - y_first) / (x_last - x_first)` (constant with one knot).
fn interp(xs: &[f64], ys: &[f64], x: f64) -> f64 {
    let last = xs.len() - 1;
    let slope = if last > 0 {
        (ys[last] - ys[0]) / (xs[last] - xs[0])
    } else {
        0.0
    };
    let i = xs.partition_point(|&k| k <= x);
    if i == 0 {
        return ys[0] + slope * (x - xs[0]);
    }
    if i == xs.len() {
        return ys[last] + slope * (x - xs[last]);
    }
    let (x0, x1, y0, y1) = (xs[i - 1], xs[i], ys[i - 1], ys[i]);
    y0 + (x - x0) * (y1 - y0) / (x1 - x0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lcg(seed: &mut u64) -> f64 {
        *seed = seed
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (*seed >> 11) as f64 / (1u64 << 53) as f64
    }

    fn normal(seed: &mut u64) -> f64 {
        let (u, v) = (lcg(seed).max(1e-12), lcg(seed));
        (-2.0 * u.ln()).sqrt() * (2.0 * std::f64::consts::PI * v).cos()
    }

    #[test]
    fn the_fit_recovers_a_monotone_map_and_its_spread_despite_wrong_anchors() {
        let mut s = 7u64;
        let truth = |x: f64| 2.0 + 10.0 * x + 3.0 * x * x;
        let mut pts = Vec::new();
        for _ in 0..1500 {
            let x = lcg(&mut s);
            pts.push((x, truth(x) + 0.3 * normal(&mut s)));
        }
        for _ in 0..300 {
            pts.push((lcg(&mut s), 15.0 * lcg(&mut s)));
        }
        let p = DiaRtPrior::fit(&pts).expect("enough anchors");
        assert_eq!(p.knots(), MAX_KNOTS);
        assert!(
            (p.sigma_min() - 0.3).abs() < 0.05,
            "sigma {}",
            p.sigma_min()
        );
        for x in [0.1, 0.5, 0.9] {
            assert!(
                (p.apply(x) - truth(x)).abs() < 0.15,
                "x {x}: {}",
                p.apply(x)
            );
        }
        assert_eq!(p.anchors(), 1800);
        assert!(p.kept() < 1800 && p.kept() > 1400, "kept {}", p.kept());
    }

    #[test]
    fn the_map_reaches_its_full_resolution_at_full_map_anchors() {
        let mut s = 11u64;
        let pts = |n: usize, s: &mut u64| -> Vec<(f64, f64)> {
            (0..n)
                .map(|_| {
                    let x = lcg(s);
                    (x, 8.0 * x + 0.2 * normal(s))
                })
                .collect()
        };
        let full = DiaRtPrior::fit(&pts(FULL_MAP_ANCHORS, &mut s)).unwrap();
        assert_eq!(full.knots(), MAX_KNOTS);
        let short = DiaRtPrior::fit(&pts(FULL_MAP_ANCHORS / 4, &mut s)).unwrap();
        assert!(short.knots() < MAX_KNOTS / 2, "knots {}", short.knots());
    }

    #[test]
    fn too_few_or_non_finite_anchors_give_no_prior() {
        let few: Vec<(f64, f64)> = (0..MIN_ANCHORS - 1).map(|i| (i as f64, i as f64)).collect();
        assert!(DiaRtPrior::fit(&few).is_none());
        let mut nan = few.clone();
        nan.push((f64::NAN, 1.0));
        assert!(
            DiaRtPrior::fit(&nan).is_none(),
            "a NaN anchor does not count"
        );
    }

    #[test]
    fn calibrated_rt_is_non_decreasing_and_extends_beyond_the_anchors() {
        let mut s = 3u64;
        let pts: Vec<(f64, f64)> = (0..400)
            .map(|_| {
                let x = lcg(&mut s);
                (x, 5.0 * x + 0.2 * normal(&mut s))
            })
            .collect();
        let p = DiaRtPrior::fit(&pts).unwrap();
        let mut last = f64::MIN;
        for i in 0..=100 {
            let y = p.apply(i as f64 / 100.0);
            assert!(y >= last);
            last = y;
        }
        assert!(p.apply(-1.0) < p.apply(0.0) - 4.0, "{}", p.apply(-1.0));
        assert!(p.apply(2.0) > p.apply(1.0) + 4.0, "{}", p.apply(2.0));
        let slope = (p.apply(3.0) - p.apply(2.0)) / 1.0;
        assert!((slope - 5.0).abs() < 0.3, "overall slope {slope}");
        assert!(p.apply(f64::NAN).is_nan());
    }

    #[test]
    fn the_penalty_is_capped_and_a_missing_prediction_gets_the_cap() {
        assert_eq!(rt_penalty(0.0, 0.4), 0.0);
        assert!((rt_penalty(0.8, 0.4) - -2.0).abs() < 1e-12);
        let cap = -RT_PENALTY_WEIGHT * RT_Z_CAP_SIGMAS * RT_Z_CAP_SIGMAS;
        assert_eq!(rt_penalty(100.0, 0.4), cap);
        assert_eq!(rt_penalty(f64::NAN, 0.4), cap);
        assert_eq!(rt_penalty(-0.8, 0.4), rt_penalty(0.8, 0.4));
    }
}
