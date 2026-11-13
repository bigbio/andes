//! Decoy hypotheses for match-between-runs.
//!
//! A decoy is a probe that is certainly false but goes through the same fit as
//! the target it shadows. Two families, one of each per target and run:
//!
//! * mass shift: the m/z moves by a whole number of nominal Da scaled by the
//!   average peptide mass-per-nominal-Da, so the decoy stays in the peptide
//!   mass-defect band and sees the same background density as a target;
//! * RT shift: the m/z is kept and the predicted apex moves far outside the
//!   target's search window.
//!
//! Both keep the target's `s` and `w_prior`, so the search breadth is equal.

/// What a probe asserts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Hypothesis {
    Target,
    /// Mass-shift decoy carrying the nominal shift `n` in Da.
    MassShiftDecoy(i8),
    RtShiftDecoy,
    /// Random-RT decoy `i` of [`random_rt_probes`]: the target's m/z at the
    /// predicted retention time of another library precursor.
    RandomRtDecoy(u8),
}

/// One (m/z, predicted apex) hypothesis to be fitted.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Probe {
    pub kind: Hypothesis,
    pub mz: f64,
    pub t_hat: f64,
    pub s: f64,
    pub w_prior: f64,
}

/// Nominal mass shifts (Da) of the mass-shift decoy family.
pub const MASS_SHIFT_UNITS: [i8; 4] = [-11, -7, 7, 11];
/// Average peptide mass per nominal Da.
pub const PEPTIDE_MASS_SLOPE: f64 = 1.000495;

fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    x ^ (x >> 31)
}

/// The target probe followed by its mass-shift decoy and, when a sign fits
/// inside `run_rt`, its RT-shift decoy.
///
/// `charge` must be positive (it divides the mass shift); a zero charge is a
/// caller bug and trips a debug assertion. `run_rt` is the (start, end) of the
/// run in the same RT unit as `t_hat`.
#[allow(clippy::too_many_arguments)]
pub fn probes(
    precursor: u32,
    run: u16,
    mz: f64,
    charge: u8,
    t_hat: f64,
    s: f64,
    w_prior: f64,
    run_rt: (f64, f64),
) -> Vec<Probe> {
    debug_assert!(charge > 0, "charge must be positive");
    let h = splitmix64(((precursor as u64) << 16) | run as u64);
    let n = MASS_SHIFT_UNITS[(h % MASS_SHIFT_UNITS.len() as u64) as usize];
    let positive = (h >> 2) & 1 == 1;

    let mut out = Vec::with_capacity(3);
    out.push(Probe {
        kind: Hypothesis::Target,
        mz,
        t_hat,
        s,
        w_prior,
    });
    out.push(Probe {
        kind: Hypothesis::MassShiftDecoy(n),
        mz: mz + n as f64 * PEPTIDE_MASS_SLOPE / charge as f64,
        t_hat,
        s,
        w_prior,
    });

    let d = (8.0 * s + 4.0 * w_prior).max(0.08 * (run_rt.1 - run_rt.0));
    let sign = if positive { 1.0 } else { -1.0 };
    let shifted = [t_hat + sign * d, t_hat - sign * d]
        .into_iter()
        .find(|t| *t >= run_rt.0 && *t <= run_rt.1);
    if let Some(t) = shifted {
        out.push(Probe {
            kind: Hypothesis::RtShiftDecoy,
            mz,
            t_hat: t,
            s,
            w_prior,
        });
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mass_shift_stays_in_mass_defect_band() {
        let ps = probes(5, 2, 700.0, 2, 1500.0, 10.0, 6.0, (0.0, 3600.0));
        let d = ps
            .iter()
            .find(|p| matches!(p.kind, Hypothesis::MassShiftDecoy(_)))
            .unwrap();
        let Hypothesis::MassShiftDecoy(n) = d.kind else {
            unreachable!()
        };
        assert!(((d.mz - 700.0) * 2.0 - n as f64 * PEPTIDE_MASS_SLOPE).abs() < 1e-9);
        assert_eq!(d.t_hat, 1500.0);
        assert_eq!((d.s, d.w_prior), (10.0, 6.0));
    }

    #[test]
    fn rt_shift_keeps_search_breadth_and_stays_in_run() {
        let run_rt = (0.0, 3600.0);
        for pre in 0..200u32 {
            let ps = probes(pre, 1, 600.0, 3, 1800.0, 10.0, 6.0, run_rt);
            assert_eq!(ps[0].kind, Hypothesis::Target);
            assert_eq!(
                (ps[0].mz, ps[0].t_hat, ps[0].s, ps[0].w_prior),
                (600.0, 1800.0, 10.0, 6.0)
            );
            let r = ps
                .iter()
                .find(|p| p.kind == Hypothesis::RtShiftDecoy)
                .unwrap();
            assert_eq!((r.mz, r.s, r.w_prior), (600.0, 10.0, 6.0));
            assert!(r.t_hat >= run_rt.0 && r.t_hat <= run_rt.1);
            assert!((r.t_hat - 1800.0).abs() >= 8.0 * 10.0 + 4.0 * 6.0);
        }
    }

    #[test]
    fn rt_shift_flips_sign_at_run_edge() {
        // Near the start only the positive direction fits, whatever the hash says.
        for pre in 0..50u32 {
            let ps = probes(pre, 0, 500.0, 2, 100.0, 10.0, 6.0, (0.0, 3600.0));
            let r = ps
                .iter()
                .find(|p| p.kind == Hypothesis::RtShiftDecoy)
                .unwrap();
            assert!(r.t_hat > 100.0 && r.t_hat <= 3600.0);
        }
    }

    #[test]
    fn deterministic_per_precursor_and_run() {
        let a = probes(9, 3, 650.5, 2, 900.0, 8.0, 5.0, (0.0, 3000.0));
        let b = probes(9, 3, 650.5, 2, 900.0, 8.0, 5.0, (0.0, 3000.0));
        assert_eq!(a, b);
        let ns: std::collections::HashSet<i8> = (0..64u16)
            .filter_map(
                |r| match probes(9, r, 650.5, 2, 900.0, 8.0, 5.0, (0.0, 3000.0))[1].kind {
                    Hypothesis::MassShiftDecoy(n) => Some(n),
                    _ => None,
                },
            )
            .collect();
        assert!(ns.len() > 1, "different runs should pick different shifts");
    }

    #[test]
    fn no_rt_decoy_when_run_too_short() {
        let ps = probes(5, 2, 700.0, 2, 1500.0, 10.0, 6.0, (1490.0, 1510.0));
        assert_eq!(ps.len(), 2);
        assert_eq!(ps[0].kind, Hypothesis::Target);
        assert!(matches!(ps[1].kind, Hypothesis::MassShiftDecoy(_)));
    }

    #[test]
    fn shifts_avoid_common_deltas() {
        let deltas = [
            0.984, 14.016, 15.995, 17.027, 18.011, 21.982, 27.995, 31.990, 37.956,
        ];
        for n in MASS_SHIFT_UNITS {
            for d in deltas {
                assert!(((n as f64).abs() * PEPTIDE_MASS_SLOPE - d).abs() > 0.9);
            }
        }
    }
}

/// Random-RT decoys per transfer target ([`random_rt_probes`]).
pub const RANDOM_RT_DECOYS: usize = 4;
/// Neutral-mass distance (Da) of the precursors whose predicted RTs the
/// random-RT decoys borrow: close enough to share the target's background,
/// far enough not to be an isotope or a common modification of it.
pub const RANDOM_RT_MASS: (f64, f64) = (5.0, 11.0);

/// The predicted retention times of a run's library precursors, by neutral
/// mass: the RTs the random-RT decoys of the run are drawn from.
#[derive(Debug, Clone, Default)]
pub struct RtPool {
    /// (neutral mass, predicted RT), ascending in mass.
    at: Vec<(f64, f64)>,
}

impl RtPool {
    /// `items` are (neutral mass, predicted RT) pairs in any order.
    pub fn new(mut items: Vec<(f64, f64)>) -> Self {
        items.retain(|(m, t)| m.is_finite() && t.is_finite());
        items.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
        Self { at: items }
    }

    pub fn len(&self) -> usize {
        self.at.len()
    }

    pub fn is_empty(&self) -> bool {
        self.at.is_empty()
    }
}

/// Up to `k` random-RT decoys of a target (`mz`, neutral mass `mass`,
/// predicted at `t_hat` with RT spread `s` and width prior `w_prior`), after
/// PIP-ECHO (Solivais et al., J Proteome Res 2025): each is the target's
/// own m/z at the predicted RT of another library precursor of the run whose
/// neutral mass is [`RANDOM_RT_MASS`] away and whose RT lies outside the
/// target's search window (more than `4 (s + w_prior)` from `t_hat`), inside
/// `run_rt`. The RTs follow where real peptides elute, so the decoys see the
/// background a wrong transfer sees, unlike an RT far outside the window.
/// The draw is deterministic in (precursor, run); fewer than `k` come back
/// when the pool has fewer candidates.
#[allow(clippy::too_many_arguments)]
pub fn random_rt_probes(
    precursor: u32,
    run: u16,
    mz: f64,
    mass: f64,
    t_hat: f64,
    s: f64,
    w_prior: f64,
    run_rt: (f64, f64),
    pool: &RtPool,
    k: usize,
) -> Vec<Probe> {
    let (lo, hi) = RANDOM_RT_MASS;
    let sep = 4.0 * (s + w_prior);
    let range = |a: f64, b: f64| {
        let i = pool.at.partition_point(|x| x.0 < a);
        let j = pool.at.partition_point(|x| x.0 <= b);
        i..j
    };
    let mut cand: Vec<f64> = range(mass - hi, mass - lo)
        .chain(range(mass + lo, mass + hi))
        .map(|i| pool.at[i].1)
        .filter(|&t| (t - t_hat).abs() > sep && t >= run_rt.0 && t <= run_rt.1)
        .collect();
    let mut out = Vec::with_capacity(k);
    let mut h = splitmix64(((precursor as u64) << 16) | run as u64 | (1u64 << 63));
    while out.len() < k && !cand.is_empty() {
        h = splitmix64(h);
        let t = cand.swap_remove((h % cand.len() as u64) as usize);
        out.push(Probe {
            kind: Hypothesis::RandomRtDecoy(out.len() as u8),
            mz,
            t_hat: t,
            s,
            w_prior,
        });
    }
    out
}

#[cfg(test)]
mod random_rt_tests {
    use super::*;

    #[test]
    fn random_rt_decoys_borrow_rts_of_nearby_masses_outside_the_window() {
        // Masses 1000 (the target's), 994, 1003 (too close), 1008, 1020 (too far).
        let pool = RtPool::new(vec![(994.0, 300.0), (1003.0, 900.0), (1008.0, 1500.0), (1020.0, 2000.0), (1008.0, 1010.0)]);
        let ps = random_rt_probes(7, 1, 501.0, 1000.0, 1000.0, 5.0, 5.0, (0.0, 3600.0), &pool, 4);
        // 994 -> 300 and 1008 -> 1500 qualify; 1008 -> 1010 is inside the window (sep 40 s).
        let mut ts: Vec<f64> = ps.iter().map(|p| p.t_hat).collect();
        ts.sort_by(f64::total_cmp);
        assert_eq!(ts, vec![300.0, 1500.0]);
        assert!(ps.iter().all(|p| p.mz == 501.0 && matches!(p.kind, Hypothesis::RandomRtDecoy(_))));
        // Deterministic in (precursor, run).
        let again = random_rt_probes(7, 1, 501.0, 1000.0, 1000.0, 5.0, 5.0, (0.0, 3600.0), &pool, 4);
        assert_eq!(again, ps);
        // An empty pool gives none.
        assert!(random_rt_probes(7, 1, 501.0, 1000.0, 1000.0, 5.0, 5.0, (0.0, 3600.0), &RtPool::default(), 4).is_empty());
    }
}
