//! Same-run claims on an MS1 feature: which other probed target of the run
//! sits on the same feature.
//!
//! A target row (identified or transferred, fitted) is *claimed* by another
//! target row of the same run when the two have the same charge, neutral
//! masses within [`CLAIM_PPM`] of each other and apexes within
//! [`CLAIM_WINDOW_S`] seconds. The claim is *co-apex* when the apexes are
//! closer than the sum of the two fitted widths: the two probes then measure
//! one MS1 feature, whose quantity is inherently shared between the two
//! compounds (isobaric residue combinations, permutations) or belongs to the
//! claimant alone (a transfer into a run that does not contain the analyte).
//!
//! The claim is a fact about the run's probes, not a verdict: nothing here
//! removes or re-weights a measurement. The best claimant is the identified
//! one when there is any, then the one with the nearest apex.

/// Mass tolerance of a claim, parts per million of the neutral mass.
pub const CLAIM_PPM: f64 = 2.0;
/// Largest apex difference of a claim, seconds.
pub const CLAIM_WINDOW_S: f64 = 60.0;

/// One fitted target row of a run, as the claim search sees it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Probe {
    /// Precursor id (distinguishes the row from its own claimants).
    pub precursor: u32,
    pub charge: u8,
    /// Neutral mass, Da.
    pub mass: f64,
    /// Fitted apex and width, seconds.
    pub t0: f64,
    pub w: f64,
    /// The precursor is identified in this run.
    pub identified: bool,
}

/// The claims on one row.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Claim {
    /// Best claimant: identified first, then nearest apex; `None` without a
    /// claimant.
    pub claimant: Option<u32>,
    /// Apex of the best claimant minus the row's apex, seconds.
    pub dt_s: Option<f32>,
    /// The best claimant is identified in the run.
    pub identified: bool,
    /// The best claimant is co-apex.
    pub co_apex: bool,
    /// Claimants within the window.
    pub n: u16,
    /// Co-apex claimants.
    pub n_co_apex: u16,
    /// Co-apex claimants identified in the run.
    pub n_identified_co_apex: u16,
}

/// The claims on the rows `rows` (indices into `probes`) from all `probes`;
/// one [`Claim`] per index of `rows`, in order. Rows or probes without a
/// finite apex, width or mass take no part.
pub fn assign(probes: &[Probe], rows: &[usize]) -> Vec<Claim> {
    let mut order: Vec<usize> = (0..probes.len())
        .filter(|&i| {
            let p = &probes[i];
            p.mass.is_finite() && p.t0.is_finite() && p.w.is_finite()
        })
        .collect();
    order.sort_unstable_by(|&a, &b| {
        (probes[a].charge, probes[a].mass)
            .partial_cmp(&(probes[b].charge, probes[b].mass))
            .expect("finite masses")
    });
    rows.iter()
        .map(|&i| {
            let me = &probes[i];
            if !(me.mass.is_finite() && me.t0.is_finite() && me.w.is_finite()) {
                return Claim::default();
            }
            let tol = me.mass * CLAIM_PPM * 1e-6;
            let lo = order.partition_point(|&j| {
                (probes[j].charge, probes[j].mass) < (me.charge, me.mass - tol)
            });
            let mut claim = Claim::default();
            let mut best: Option<(bool, f64, &Probe)> = None;
            for &j in &order[lo..] {
                let p = &probes[j];
                if p.charge != me.charge || p.mass > me.mass + tol {
                    break;
                }
                let dt = p.t0 - me.t0;
                if p.precursor == me.precursor || dt.abs() > CLAIM_WINDOW_S {
                    continue;
                }
                claim.n += 1;
                if dt.abs() <= me.w + p.w {
                    claim.n_co_apex += 1;
                    claim.n_identified_co_apex += u16::from(p.identified);
                }
                let better = match best {
                    None => true,
                    Some((id, d, _)) => {
                        (p.identified, std::cmp::Reverse(ord(dt.abs())))
                            > (id, std::cmp::Reverse(ord(d)))
                    }
                };
                if better {
                    best = Some((p.identified, dt.abs(), p));
                }
            }
            if let Some((_, _, p)) = best {
                let dt = p.t0 - me.t0;
                claim.claimant = Some(p.precursor);
                claim.dt_s = Some(dt as f32);
                claim.identified = p.identified;
                claim.co_apex = dt.abs() <= me.w + p.w;
            }
            claim
        })
        .collect()
}

/// Total order key of a finite non-negative float.
fn ord(x: f64) -> u64 {
    x.to_bits()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn probe(precursor: u32, mass: f64, t0: f64, identified: bool) -> Probe {
        Probe {
            precursor,
            charge: 2,
            mass,
            t0,
            w: 4.0,
            identified,
        }
    }

    #[test]
    fn co_apex_identified_claimant_is_best() {
        let probes = [
            probe(0, 1000.0, 100.0, false),
            probe(1, 1000.0 + 1.0e-3, 130.0, true),
            probe(2, 1000.0 + 1.0e-3, 102.0, false),
            probe(3, 1000.0 + 0.1, 100.0, true),
        ];
        let c = assign(&probes, &[0]);
        assert_eq!(c[0].claimant, Some(1));
        assert_eq!(c[0].n, 2);
        assert_eq!(c[0].n_co_apex, 1);
        assert_eq!(c[0].n_identified_co_apex, 0);
        assert!(c[0].identified && !c[0].co_apex);
        assert_eq!(c[0].dt_s, Some(30.0));
    }

    #[test]
    fn nearest_apex_among_equals_and_window() {
        let probes = [
            probe(0, 1000.0, 100.0, false),
            probe(1, 1000.0, 105.0, false),
            probe(2, 1000.0, 99.0, false),
            probe(3, 1000.0, 400.0, true),
        ];
        let c = assign(&probes, &[0]);
        assert_eq!(c[0].claimant, Some(2));
        assert!(c[0].co_apex && !c[0].identified);
        assert_eq!(c[0].n, 2);
        assert_eq!(c[0].n_co_apex, 2);
    }

    #[test]
    fn charge_mass_and_self_are_excluded() {
        let mut other_charge = probe(1, 1000.0, 100.0, true);
        other_charge.charge = 3;
        let probes = [
            probe(0, 1000.0, 100.0, false),
            other_charge,
            probe(2, 1000.0 * (1.0 + 3.0e-6), 100.0, true),
        ];
        let c = assign(&probes, &[0]);
        assert_eq!(c[0], Claim::default());
    }

    #[test]
    fn unfitted_probes_take_no_part() {
        let mut nan = probe(1, 1000.0, f64::NAN, true);
        nan.w = f64::NAN;
        let probes = [probe(0, 1000.0, 100.0, false), nan];
        assert_eq!(assign(&probes, &[0, 1]), vec![Claim::default(); 2]);
    }
}
