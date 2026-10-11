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
