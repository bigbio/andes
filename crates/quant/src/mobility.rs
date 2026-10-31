//! The ion-mobility window of a run's MS1 extraction.
//!
//! On a run whose survey scans are resolved by ion mobility
//! ([`crate::ms1_index::Ms1Scan::mobility`]), a precursor's traces are read
//! within a window of 1/K0 around its library mobility (the median precursor
//! 1/K0 of its donor PSMs), so a co-eluting ion at the same m/z but another
//! mobility does not enter its traces. The window is learned per run from
//! the run's identified precursors: for each, the 1/K0 of the most intense
//! MS1 centroid of its monoisotopic m/z near its PSM's retention time,
//! against its library mobility. The median of those differences is the
//! run's offset (the run's mobility scale against the library's) and their
//! robust standard deviation (1.4826 times the median absolute deviation) its
//! scale; the window is the library mobility plus the offset, within
//! [`GATE_SIGMAS`] scales. The differences are collected within
//! [`SEARCH_K0`] of the library mobility, then again within [`REFINE_SIGMAS`]
//! scales of the first offset, and the offset and scale refitted.

use crate::ms1_index::{MobilityGate, Ms1RunIndex};

/// Half-width (V·s/cm²) of the first search for each identified precursor's
/// centroid around its library mobility; it only bounds the pairs the
/// offset and scale are estimated from.
pub const SEARCH_K0: f64 = 0.1;
/// Pairs beyond this many scales of the first offset are left out of the
/// refit.
pub const REFINE_SIGMAS: f64 = 4.0;
/// Half-width of the extraction window, in learned scales.
pub const GATE_SIGMAS: f64 = 3.0;
/// Fewest pairs a window is learned from; with fewer, the run has none.
pub const MIN_PAIRS: usize = 50;
/// Floor of the learned scale (V·s/cm²), about one mobility scan.
pub const MIN_SIGMA: f64 = 1e-3;

/// One identified precursor of the run: monoisotopic m/z, PSM retention time
/// (seconds) and library mobility (1/K0).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MobilityObservation {
    pub mz: f64,
    pub rt: f64,
    pub k0: f32,
}

/// A run's learned mobility scale (see the module docs).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RunMobility {
    /// Median of (centroid 1/K0 − library 1/K0) over the pairs.
    pub offset: f64,
    /// Robust standard deviation of those differences, floored at
    /// [`MIN_SIGMA`].
    pub sigma: f64,
    /// Pairs in the final estimate.
    pub pairs: usize,
}

impl RunMobility {
    /// The extraction window of a precursor with library mobility `k0`.
    pub fn gate(&self, k0: f32) -> MobilityGate {
        let centre = k0 as f64 + self.offset;
        let half = GATE_SIGMAS * self.sigma;
        MobilityGate {
            lo: (centre - half) as f32,
            hi: (centre + half) as f32,
        }
    }
}

fn median(v: &mut [f64]) -> f64 {
    v.sort_by(f64::total_cmp);
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2]
    } else {
        0.5 * (v[n / 2 - 1] + v[n / 2])
    }
}

/// Median and robust standard deviation of `d` (not empty).
fn location_scale(d: &[f64]) -> (f64, f64) {
    let mut v = d.to_vec();
    let m = median(&mut v);
    let mut dev: Vec<f64> = d.iter().map(|x| (x - m).abs()).collect();
    (m, 1.4826 * median(&mut dev))
}

/// The 1/K0 differences of the observations whose most intense centroid in
/// the isotope-0 slice (`± tol_ppm`), over the scans within `rt_half_s` of
/// the PSM, lies within `half` of `k0 + centre`.
fn differences(
    ms1: &Ms1RunIndex,
    obs: &[MobilityObservation],
    tol_ppm: f64,
    rt_half_s: f64,
    centre: f64,
    half: f64,
) -> Vec<f64> {
    let mut out = Vec::new();
    for o in obs {
        if !(o.mz.is_finite() && o.rt.is_finite() && o.k0.is_finite()) {
            continue;
        }
        let tol = o.mz * tol_ppm * 1e-6;
        let (lo, hi) = (o.mz - tol, o.mz + tol);
        let target = o.k0 as f64 + centre;
        let mut best: Option<(f32, f32)> = None;
        for j in ms1.scan_window(o.rt - rt_half_s, o.rt + rt_half_s) {
            let peaks = ms1.peaks(j);
            let k0s = ms1.mobility(j);
            if k0s.is_empty() {
                continue;
            }
            let start = peaks.partition_point(|p| p.0 < lo);
            for (k, &(mz, inten)) in peaks.iter().enumerate().skip(start) {
                if mz > hi {
                    break;
                }
                if (k0s[k] as f64 - target).abs() > half {
                    continue;
                }
                if inten > best.map_or(0.0, |b| b.0) {
                    best = Some((inten, k0s[k]));
                }
            }
        }
        if let Some((_, k0)) = best {
            out.push(k0 as f64 - o.k0 as f64);
        }
    }
    out
}

/// Learn the run's mobility scale from its identified precursors `obs`
/// (see the module docs); `tol_ppm` is the extraction's m/z half-width and
/// `rt_half_s` the retention-time half-width searched around each PSM.
/// `None` when the run is not mobility-resolved or gives fewer than
/// [`MIN_PAIRS`] pairs.
pub fn learn(
    ms1: &Ms1RunIndex,
    obs: &[MobilityObservation],
    tol_ppm: f64,
    rt_half_s: f64,
) -> Option<RunMobility> {
    if !ms1.has_mobility() {
        return None;
    }
    let first = differences(ms1, obs, tol_ppm, rt_half_s, 0.0, SEARCH_K0);
    if first.len() < MIN_PAIRS {
        return None;
    }
    let (offset, sigma) = location_scale(&first);
    let half = REFINE_SIGMAS * sigma.max(MIN_SIGMA);
    let second = differences(ms1, obs, tol_ppm, rt_half_s, offset, half);
    if second.len() < MIN_PAIRS {
        return None;
    }
    let (offset, sigma) = location_scale(&second);
    Some(RunMobility {
        offset,
        sigma: sigma.max(MIN_SIGMA),
        pairs: second.len(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ms1_index::Ms1Scan;

    /// Scans every second: precursor `p` (m/z 500 + p) elutes at 100 + 10 p
    /// s with 1/K0 0.8 + 0.001 p + `shift`; an interferer at the same m/z,
    /// four times as intense, sits 0.2 V·s/cm² away (outside the search).
    fn run(n: usize, shift: f32) -> Ms1RunIndex {
        let scans = (0..(20 * n + 200))
            .map(|k| {
                let rt = k as f64 * 0.5;
                let mut rows: Vec<(f64, f32, f32)> = Vec::new();
                for p in 0..n {
                    let t0 = 100.0 + 10.0 * p as f64;
                    let a = (-0.5 * ((rt - t0) / 3.0).powi(2)).exp() as f32;
                    if a < 1e-3 {
                        continue;
                    }
                    let mz = 500.0 + p as f64;
                    let k0 = 0.8 + 0.001 * p as f32 + shift;
                    // Centroid mobility jitter of +-0.002.
                    let jitter = if (k + p) % 2 == 0 { 0.002 } else { -0.002 };
                    rows.push((mz, 1000.0 * a, k0 + jitter));
                    rows.push((mz + 1e-4, 4000.0 * a, k0 + 0.2));
                }
                rows.sort_by(|a, b| a.0.total_cmp(&b.0));
                Ms1Scan {
                    rt,
                    peaks: rows.iter().map(|r| (r.0, r.1)).collect(),
                    mobility: rows.iter().map(|r| r.2).collect(),
                    ..Ms1Scan::default()
                }
            })
            .collect();
        Ms1RunIndex::new(scans)
    }

    fn obs(n: usize) -> Vec<MobilityObservation> {
        (0..n)
            .map(|p| MobilityObservation {
                mz: 500.0 + p as f64,
                rt: 100.0 + 10.0 * p as f64,
                k0: 0.8 + 0.001 * p as f32,
            })
            .collect()
    }

    #[test]
    fn learns_the_run_offset_and_scale() {
        let m = learn(&run(200, 0.015), &obs(200), 10.0, 6.0).expect("enough pairs");
        assert!((m.offset - 0.015).abs() < 0.0025, "{m:?}");
        assert!(m.sigma <= 0.004, "{m:?}");
        assert_eq!(m.pairs, 200);
        let g = m.gate(0.9);
        assert!(g.contains(0.915) && !g.contains(0.9 + 0.2), "{g:?}");
    }

    #[test]
    fn no_window_without_mobility_or_pairs() {
        let plain = Ms1RunIndex::new(vec![Ms1Scan {
            rt: 1.0,
            peaks: vec![(500.0, 1.0)],
            ..Ms1Scan::default()
        }]);
        assert!(learn(&plain, &obs(10), 10.0, 6.0).is_none());
        assert!(learn(&run(20, 0.0), &obs(20), 10.0, 6.0).is_none());
    }
}
