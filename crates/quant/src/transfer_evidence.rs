//! Extra evidence for match-between-runs transfers, computed the same way for
//! a target probe and for each of its decoy probes.
//!
//! * Charge-state partner: the precursor's neutral mass at charge `z - 1` and
//!   `z + 1` (the probe's own m/z, so a mass-shift decoy's partner is at its
//!   shifted mass). Each partner is fitted twice: with the probe's own RT
//!   prior, and anchored on the probe's fitted apex (prior centred on `t0`
//!   with `s = w`, profile width `w`). [`PartnerEvidence`] keeps, per prior,
//!   the best partner LLR and its apex offset from the probe's apex. A peptide
//!   that is present often shows another charge state at the same apex; a
//!   decoy rarely does.
//! * Local RT: the alignment predicts `t_hat` from a smooth map per run. The
//!   [`LocalRt`] correction adds the median residual of the `k` identified
//!   precursors of the receiving run nearest in predicted RT (each residual
//!   against its own leave-this-run-out prediction), and replaces the run's
//!   residual width in `s` with the spread of those residuals. A decoy keeps
//!   its target's corrected prior.
//!
//! Neither uses a quantity measured in another run.

use model::mass::PROTON;

use crate::evidence::Evidence;

/// Highest partner charge fitted.
pub const MAX_PARTNER_CHARGE: u8 = 8;
/// Floor of the anchored partner's prior width and profile width, seconds.
pub const MIN_ANCHOR_S: f64 = 1.0;
/// Floor of a local RT uncertainty, seconds.
pub const MIN_LOCAL_S: f64 = 0.5;
/// Consistency factor of the median absolute deviation for a Gaussian.
const MAD_SCALE: f64 = 1.4826;

/// `(charge, m/z)` of the charge-state partners of a probe at monoisotopic
/// m/z `mz` and charge `charge`: `z - 1` (when at least 1) and `z + 1` (when
/// at most [`MAX_PARTNER_CHARGE`]).
pub fn partner_mzs(mz: f64, charge: u8) -> Vec<(u8, f64)> {
    let z = charge.max(1);
    let neutral = (mz - PROTON) * z as f64;
    [z.checked_sub(1), z.checked_add(1)]
        .into_iter()
        .flatten()
        .filter(|&c| (1..=MAX_PARTNER_CHARGE).contains(&c))
        .map(|c| (c, neutral / c as f64 + PROTON))
        .collect()
}

/// `(t_hat, s, w_prior)` of the anchored partner prior of a probe fit.
pub fn anchored_prior(main: &Evidence) -> (f64, f64, f64) {
    let w = main.w.max(MIN_ANCHOR_S);
    (main.t0, w, w)
}

/// The partner evidence of one probe.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PartnerEvidence {
    /// Best partner LLR with the probe's own prior (0 without a fit).
    pub llr: f64,
    /// That partner's apex minus the probe's apex, seconds (0 without a fit).
    pub dt: f64,
    /// Best partner LLR with the prior anchored on the probe's apex.
    pub coel_llr: f64,
    /// That partner's apex minus the probe's apex, seconds.
    pub coel_dt: f64,
}

/// `(llr, t0 - main.t0)` of the fit with the highest LLR among `fits`
/// (`(0, 0)` when none).
fn best(main: &Evidence, fits: &[Option<Evidence>]) -> (f64, f64) {
    fits.iter()
        .flatten()
        .filter(|e| e.llr.is_finite() && e.t0.is_finite())
        .max_by(|a, b| a.llr.total_cmp(&b.llr))
        .map_or((0.0, 0.0), |e| (e.llr, e.t0 - main.t0))
}

impl PartnerEvidence {
    /// From the partner fits with the probe's prior (`same`) and anchored on
    /// its apex (`anchored`).
    pub fn from_fits(
        main: &Evidence,
        same: &[Option<Evidence>],
        anchored: &[Option<Evidence>],
    ) -> Self {
        let (llr, dt) = best(main, same);
        let (coel_llr, coel_dt) = best(main, anchored);
        Self {
            llr,
            dt,
            coel_llr,
            coel_dt,
        }
    }
}

/// A locally corrected RT prior.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LocalPrior {
    pub t_hat: f64,
    pub s: f64,
}

/// The local RT correction of one receiving run.
#[derive(Debug, Clone)]
pub struct LocalRt {
    /// `(predicted RT, residual, precursor)` of the run's identified
    /// precursors, ascending by predicted RT.
    points: Vec<(f64, f64, u32)>,
    k: usize,
}

impl LocalRt {
    /// From `(precursor, observed RT, predicted RT)` of the run's identified
    /// precursors (the prediction without this run's observation); points
    /// that are not finite are dropped. `k` is the neighbour count.
    pub fn new(points: impl IntoIterator<Item = (u32, f64, f64)>, k: usize) -> Self {
        let mut points: Vec<(f64, f64, u32)> = points
            .into_iter()
            .filter(|(_, rt, t)| rt.is_finite() && t.is_finite())
            .map(|(p, rt, t)| (t, rt - t, p))
            .collect();
        points.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.2.cmp(&b.2)));
        Self { points, k: k.max(3) }
    }

    /// Number of reference points.
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// True without reference points.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }

    /// The `k` points nearest to `t_hat` in predicted RT, precursor
    /// `exclude` left out (ties broken towards earlier RT).
    fn neighbours(&self, t_hat: f64, exclude: u32) -> Vec<f64> {
        let pts = &self.points;
        let mut hi = pts.partition_point(|p| p.0 < t_hat);
        let mut lo = hi;
        let mut out = Vec::with_capacity(self.k);
        while out.len() < self.k && (lo > 0 || hi < pts.len()) {
            let left = lo.checked_sub(1).map(|i| t_hat - pts[i].0);
            let right = (hi < pts.len()).then(|| pts[hi].0 - t_hat);
            let take_left = match (left, right) {
                (Some(l), Some(r)) => l <= r,
                (Some(_), None) => true,
                _ => false,
            };
            let i = if take_left {
                lo -= 1;
                lo
            } else {
                hi += 1;
                hi - 1
            };
            if pts[i].2 != exclude {
                out.push(pts[i].1);
            }
        }
        out
    }

    /// The corrected prior of precursor `p` predicted at `(t_hat, s)` in this
    /// run, whose residual width there is `sigma_run`: `t_hat` plus the
    /// median neighbour residual, and `s^2 - sigma_run^2 + sigma_local^2`
    /// (the first difference floored at 0), with `sigma_local` the scaled
    /// MAD of the neighbour residuals, floored at [`MIN_LOCAL_S`]. `None`
    /// with fewer than `k` neighbours or a non-finite input.
    pub fn correct(&self, p: u32, t_hat: f64, s: f64, sigma_run: f64) -> Option<LocalPrior> {
        if !(t_hat.is_finite() && s.is_finite() && sigma_run.is_finite()) {
            return None;
        }
        let mut r = self.neighbours(t_hat, p);
        if r.len() < self.k {
            return None;
        }
        let delta = median(&mut r);
        let mut dev: Vec<f64> = r.iter().map(|v| (v - delta).abs()).collect();
        let sigma_local = MAD_SCALE * median(&mut dev);
        let rest = (s * s - sigma_run * sigma_run).max(0.0);
        Some(LocalPrior {
            t_hat: t_hat + delta,
            s: (rest + sigma_local * sigma_local).sqrt().max(MIN_LOCAL_S),
        })
    }
}

/// Median of a non-empty slice (sorted in place).
fn median(v: &mut [f64]) -> f64 {
    v.sort_by(f64::total_cmp);
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2]
    } else {
        0.5 * (v[n / 2 - 1] + v[n / 2])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ev(llr: f64, t0: f64, w: f64) -> Evidence {
        Evidence {
            amplitude: 1.0,
            amplitude_se: 0.1,
            quantity: 1.0,
            quantity_se: 0.1,
            llr,
            rt_z: 0.0,
            width_ratio: 1.0,
            envelope_cos: 1.0,
            coherence: 1.0,
            m1_ratio: 0.0,
            mass_err_ppm: 0.0,
            n_points: 10,
            isotopes_used: 0b111,
            t0,
            w,
            iso_heights: [0.0; 4],
            tailing: 0.0,
        }
    }

    #[test]
    fn partners_keep_the_neutral_mass() {
        let mz = 600.0;
        let m = (mz - PROTON) * 2.0;
        let ps = partner_mzs(mz, 2);
        assert_eq!(ps.len(), 2);
        for (z, pmz) in ps {
            assert!(((pmz - PROTON) * z as f64 - m).abs() < 1e-9);
        }
        assert_eq!(partner_mzs(900.0, 1).iter().map(|p| p.0).collect::<Vec<_>>(), [2]);
        assert_eq!(
            partner_mzs(500.0, MAX_PARTNER_CHARGE)
                .iter()
                .map(|p| p.0)
                .collect::<Vec<_>>(),
            [MAX_PARTNER_CHARGE - 1]
        );
    }

    #[test]
    fn partner_evidence_takes_the_best_fit_and_its_offset() {
        let main = ev(100.0, 1000.0, 5.0);
        let same = [Some(ev(20.0, 1030.0, 5.0)), Some(ev(50.0, 990.0, 5.0)), None];
        let anchored = [None, Some(ev(40.0, 1001.0, 5.0))];
        let p = PartnerEvidence::from_fits(&main, &same, &anchored);
        assert_eq!((p.llr, p.dt), (50.0, -10.0));
        assert_eq!((p.coel_llr, p.coel_dt), (40.0, 1.0));
        let none = PartnerEvidence::from_fits(&main, &[None], &[]);
        assert_eq!((none.llr, none.dt, none.coel_llr, none.coel_dt), (0.0, 0.0, 0.0, 0.0));
        assert_eq!(anchored_prior(&main), (1000.0, 5.0, 5.0));
        assert_eq!(anchored_prior(&ev(1.0, 10.0, 0.2)).1, MIN_ANCHOR_S);
    }

    #[test]
    fn local_correction_follows_the_neighbours_and_excludes_the_precursor() {
        // Residual +3 s around 1000 s, -3 s around 2000 s, small jitter.
        let pts = (0..200u32).map(|i| {
            let t = 500.0 + 10.0 * i as f64;
            let jitter = ((i * 7919) % 11) as f64 * 0.05 - 0.25;
            let r = if t < 1500.0 { 3.0 } else { -3.0 } + jitter;
            (i, t + r, t)
        });
        let l = LocalRt::new(pts, 20);
        assert_eq!(l.len(), 200);
        let a = l.correct(10_000, 1000.0, 5.0, 4.0).unwrap();
        assert!((a.t_hat - 1003.0).abs() < 0.3, "{a:?}");
        assert!(a.s < 5.0 && a.s >= MIN_LOCAL_S, "{a:?}");
        let b = l.correct(10_000, 2200.0, 5.0, 4.0).unwrap();
        assert!((b.t_hat - 2197.0).abs() < 0.3, "{b:?}");
        // The precursor itself is never its own neighbour: an outlier at
        // 1000 s does not move its own correction.
        let mut with_outlier: Vec<(u32, f64, f64)> =
            (0..200u32).map(|i| (i, 500.0 + 10.0 * i as f64, 500.0 + 10.0 * i as f64)).collect();
        with_outlier[50] = (50, 1000.0 + 60.0, 1000.0);
        let l = LocalRt::new(with_outlier, 3);
        assert_eq!(l.correct(50, 1000.0, 5.0, 4.0).unwrap().t_hat, 1000.0);
        // Too few points: no correction.
        assert!(LocalRt::new([(1, 10.0, 10.0)], 5).correct(2, 10.0, 1.0, 1.0).is_none());
    }
}
