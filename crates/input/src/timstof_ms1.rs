//! MS1 survey frames of a Bruker timsTOF `.d`: summed over ion mobility,
//! centroided along the time-of-flight axis, and converted to m/z with a
//! per-run calibration learned from the instrument's own precursor records.
//!
//! A TIMS MS1 frame holds one TOF spectrum per mobility scan. Quantification
//! and the MS1-linked search paths read survey scans as one centroid list per
//! retention time, so each frame is summed over its scans (per TOF index),
//! split into peaks at the valleys of the summed profile, and each peak is
//! reported at its intensity-weighted TOF position.
//!
//! The TOF → m/z conversion that `timsrust` derives from the acquisition
//! range (`sqrt(m/z)` linear in the TOF index) is an approximation of the
//! vendor calibration. The acquisition software records, for every PASEF
//! precursor, the m/z of its largest isotope peak (`Precursors.LargestPeakMz`)
//! computed with the full instrument calibration. [`fit_calibration`] matches
//! those m/z values to the centroids of their parent MS1 frames and fits
//! `sqrt(m/z)` as a cubic in the TOF index plus a linear term in each of the
//! two frame temperatures the TDF records (`Frames.T1`, `Frames.T2`), which is
//! what the vendor calibration compensates for. Pairs are matched in two
//! passes (a wide window to find the run's offset, a narrow one around it)
//! and the fit is refined by discarding residuals beyond four robust
//! standard deviations.

/// A centroid in the TOF domain: intensity-weighted TOF index, summed
/// intensity.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct TofCentroid {
    pub tof: f64,
    pub intensity: f64,
}

/// Sum the intensities of one frame's points per TOF index. `tof` and
/// `intensity` are parallel arrays over every mobility scan of the frame;
/// the result is sorted by TOF index.
pub(crate) fn sum_over_mobility(tof: &[u32], intensity: &[u32]) -> Vec<(u32, u64)> {
    let n = tof.len().min(intensity.len());
    let mut pairs: Vec<(u32, u32)> = (0..n).map(|i| (tof[i], intensity[i])).collect();
    pairs.sort_unstable_by_key(|p| p.0);
    let mut out: Vec<(u32, u64)> = Vec::with_capacity(pairs.len() / 2 + 1);
    for (t, i) in pairs {
        match out.last_mut() {
            Some(last) if last.0 == t => last.1 += u64::from(i),
            _ => out.push((t, u64::from(i))),
        }
    }
    out
}

/// Largest gap (in TOF indices) between two points of one peak. A missing
/// index between two recorded ones is a zero-count bin inside the peak.
const MAX_TOF_GAP: u32 = 2;

/// Two maxima separated by a valley deeper than this fraction of the lower
/// maximum are two peaks; a shallower dip is noise on one peak.
const VALLEY_SPLIT_FRACTION: f64 = 0.5;

/// Centroid a summed TOF profile (sorted by TOF index). Runs of points no
/// more than [`MAX_TOF_GAP`] apart are split at their local minima, and
/// neighbouring segments are merged again when the valley between them is
/// shallower than [`VALLEY_SPLIT_FRACTION`] of the lower segment maximum.
/// Each peak is reported at its intensity-weighted TOF index with its summed
/// intensity.
pub(crate) fn centroid_profile(profile: &[(u32, u64)]) -> Vec<TofCentroid> {
    let mut out = Vec::new();
    let mut start = 0;
    while start < profile.len() {
        let mut end = start + 1;
        while end < profile.len() && profile[end].0 - profile[end - 1].0 <= MAX_TOF_GAP {
            end += 1;
        }
        centroid_run(&profile[start..end], &mut out);
        start = end;
    }
    out
}

/// One peak segment: its maximum and its intensity and TOF-weighted sums.
#[derive(Debug, Clone, Copy)]
struct Segment {
    max: f64,
    sum: f64,
    weighted: f64,
}

impl Segment {
    fn of(points: &[(u32, u64)]) -> Self {
        let mut s = Self {
            max: 0.0,
            sum: 0.0,
            weighted: 0.0,
        };
        for &(t, i) in points {
            let i = i as f64;
            s.max = s.max.max(i);
            s.sum += i;
            s.weighted += i * f64::from(t);
        }
        s
    }

    fn merge(self, other: Self) -> Self {
        Self {
            max: self.max.max(other.max),
            sum: self.sum + other.sum,
            weighted: self.weighted + other.weighted,
        }
    }
}

/// Centroid one run of contiguous points (see [`centroid_profile`]).
fn centroid_run(run: &[(u32, u64)], out: &mut Vec<TofCentroid>) {
    // Split at every strict local minimum; the minimum opens the right
    // segment, and its intensity is the valley in front of that segment.
    let mut segments: Vec<(Segment, f64)> = Vec::new();
    let mut seg_start = 0;
    let mut valley = 0.0;
    for i in 1..run.len().saturating_sub(1) {
        let (prev, cur, next) = (run[i - 1].1, run[i].1, run[i + 1].1);
        if cur < prev && cur <= next {
            segments.push((Segment::of(&run[seg_start..i]), valley));
            seg_start = i;
            valley = cur as f64;
        }
    }
    segments.push((Segment::of(&run[seg_start..]), valley));
    let mut merged: Vec<Segment> = Vec::with_capacity(segments.len());
    for (seg, valley) in segments {
        match merged.last_mut() {
            Some(last) if valley >= VALLEY_SPLIT_FRACTION * last.max.min(seg.max) => {
                *last = last.merge(seg);
            }
            _ => merged.push(seg),
        }
    }
    for seg in merged {
        if seg.sum > 0.0 {
            out.push(TofCentroid {
                tof: seg.weighted / seg.sum,
                intensity: seg.sum,
            });
        }
    }
}

/// TOF → m/z conversion: `sqrt(m/z)` as a polynomial in the scaled TOF index
/// plus temperature terms. With all temperature coefficients zero it is the
/// acquisition-range model (`sqrt(m/z) = a + b * tof`).
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct MzCalibration {
    /// TOF indices are divided by this before the polynomial (numerical
    /// conditioning).
    scale: f64,
    /// Coefficients of `1, x, x^2, x^3` (x = tof / scale).
    poly: [f64; 4],
    /// Coefficients of `T1 - t1_ref` and `T2 - t2_ref`.
    temp: [f64; 2],
    t_ref: [f64; 2],
}

impl MzCalibration {
    /// The acquisition-range model `sqrt(m/z) = intercept + slope * tof`.
    pub fn linear(intercept: f64, slope: f64, scale: f64) -> Self {
        Self {
            scale,
            poly: [intercept, slope * scale, 0.0, 0.0],
            temp: [0.0; 2],
            t_ref: [0.0; 2],
        }
    }

    /// m/z of a (fractional) TOF index in a frame recorded at temperatures `t`.
    pub fn mz(&self, tof: f64, t: [f64; 2]) -> f64 {
        let x = tof / self.scale;
        let p = &self.poly;
        let root = p[0]
            + x * (p[1] + x * (p[2] + x * p[3]))
            + self.temp[0] * (t[0] - self.t_ref[0])
            + self.temp[1] * (t[1] - self.t_ref[1]);
        root * root
    }

    /// The TOF index of `mz` under the linear part of the model (the
    /// acquisition-range inverse, used to locate a vendor m/z among the
    /// centroids before the fit exists).
    pub fn linear_tof(&self, mz: f64) -> f64 {
        (mz.sqrt() - self.poly[0]) / self.poly[1] * self.scale
    }
}

/// A vendor precursor record matched against its parent frame.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct CalibrationPoint {
    /// Centroid TOF index in the parent MS1 frame.
    pub tof: f64,
    /// Vendor-calibrated m/z of the precursor's largest peak.
    pub mz: f64,
    /// Parent frame temperatures.
    pub t: [f64; 2],
}

/// Matching window (ppm) of the first pass, which finds the run's offset
/// between the acquisition-range model and the vendor m/z.
pub(crate) const WIDE_WINDOW_PPM: f64 = 200.0;
/// Half-width (ppm) of the second matching pass, around that offset.
pub(crate) const NARROW_WINDOW_PPM: f64 = 15.0;
/// Fewest matched pairs a fit is made from.
pub(crate) const MIN_CALIBRATION_POINTS: usize = 200;
/// Residuals beyond this many robust standard deviations are left out of the
/// refit.
const OUTLIER_SIGMAS: f64 = 4.0;

/// Pair each vendor m/z of a frame with the nearest centroid within
/// `offset_ppm ± window_ppm` (positions under the linear model `base`).
/// `centroids` are sorted by TOF.
pub(crate) fn match_precursors(
    centroids: &[TofCentroid],
    vendor_mz: &[f64],
    t: [f64; 2],
    base: &MzCalibration,
    offset_ppm: f64,
    window_ppm: f64,
    out: &mut Vec<CalibrationPoint>,
) {
    for &mz in vendor_mz {
        if !(mz.is_finite() && mz > 0.0) {
            continue;
        }
        let tof = base.linear_tof(mz);
        let i = centroids.partition_point(|c| c.tof < tof);
        let mut best: Option<(f64, f64)> = None; // (|ppm - offset|, tof)
        for c in centroids[i.saturating_sub(3)..(i + 3).min(centroids.len())].iter() {
            let ppm = (base.mz(c.tof, t) - mz) / mz * 1e6;
            let d = (ppm - offset_ppm).abs();
            if d <= window_ppm && best.is_none_or(|b| d < b.0) {
                best = Some((d, c.tof));
            }
        }
        if let Some((_, tof)) = best {
            out.push(CalibrationPoint { tof, mz, t });
        }
    }
}

/// Median of `v` (reorders it); `None` when empty.
fn median(v: &mut [f64]) -> Option<f64> {
    if v.is_empty() {
        return None;
    }
    let mid = v.len() / 2;
    let (_, m, _) = v.select_nth_unstable_by(mid, f64::total_cmp);
    Some(*m)
}

/// The run's median offset (ppm) of matched pairs from the linear model.
pub(crate) fn median_offset_ppm(points: &[CalibrationPoint], base: &MzCalibration) -> Option<f64> {
    let mut ppm: Vec<f64> = points
        .iter()
        .map(|p| (base.mz(p.tof, p.t) - p.mz) / p.mz * 1e6)
        .collect();
    median(&mut ppm)
}

/// Least-squares solution of `a x = b` for a small dense symmetric system
/// (Gaussian elimination with partial pivoting); `None` when singular.
fn solve(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let n = b.len();
    for col in 0..n {
        let pivot = (col..n).max_by(|&i, &j| a[i][col].abs().total_cmp(&a[j][col].abs()))?;
        if a[pivot][col].abs() < 1e-300 {
            return None;
        }
        a.swap(col, pivot);
        b.swap(col, pivot);
        for row in col + 1..n {
            let f = a[row][col] / a[col][col];
            for k in col..n {
                a[row][k] -= f * a[col][k];
            }
            b[row] -= f * b[col];
        }
    }
    let mut x = vec![0.0; n];
    for row in (0..n).rev() {
        let s: f64 = (row + 1..n).map(|k| a[row][k] * x[k]).sum();
        x[row] = (b[row] - s) / a[row][row];
    }
    x.iter().all(|v| v.is_finite()).then_some(x)
}

/// Spread (standard deviation) below which a temperature is treated as
/// constant over the run and gets no term.
const MIN_TEMPERATURE_SPREAD: f64 = 1e-3;

/// Summary of a calibration fit.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CalibrationFit {
    pub calibration: MzCalibration,
    /// Pairs used in the final fit.
    pub n_used: usize,
    /// Median and robust standard deviation (1.4826 MAD) of the final
    /// residuals against the vendor m/z, ppm.
    pub median_ppm: f64,
    pub sigma_ppm: f64,
    /// The same for the acquisition-range model on the same pairs.
    pub base_median_ppm: f64,
    pub base_sigma_ppm: f64,
}

fn residual_stats(ppm: &[f64]) -> (f64, f64) {
    let mut v = ppm.to_vec();
    let med = median(&mut v).unwrap_or(0.0);
    let mut dev: Vec<f64> = ppm.iter().map(|r| (r - med).abs()).collect();
    let mad = median(&mut dev).unwrap_or(0.0);
    (med, 1.4826 * mad)
}

/// Fit `sqrt(m/z)` against the matched pairs (see the module docs). `None`
/// with fewer than [`MIN_CALIBRATION_POINTS`] pairs or a singular system.
pub(crate) fn fit_calibration(
    points: &[CalibrationPoint],
    base: &MzCalibration,
) -> Option<CalibrationFit> {
    if points.len() < MIN_CALIBRATION_POINTS {
        return None;
    }
    let scale = base.scale;
    let mean_t = |k: usize| points.iter().map(|p| p.t[k]).sum::<f64>() / points.len() as f64;
    let t_ref = [mean_t(0), mean_t(1)];
    let spread = |k: usize| {
        (points
            .iter()
            .map(|p| (p.t[k] - t_ref[k]).powi(2))
            .sum::<f64>()
            / points.len() as f64)
            .sqrt()
    };
    let use_temp = [
        spread(0) > MIN_TEMPERATURE_SPREAD,
        spread(1) > MIN_TEMPERATURE_SPREAD,
    ];
    let features = |p: &CalibrationPoint| -> Vec<f64> {
        let x = p.tof / scale;
        let mut f = vec![1.0, x, x * x, x * x * x];
        for k in 0..2 {
            if use_temp[k] {
                f.push(p.t[k] - t_ref[k]);
            }
        }
        f
    };
    let base_ppm: Vec<f64> = points
        .iter()
        .map(|p| (base.mz(p.tof, p.t) - p.mz) / p.mz * 1e6)
        .collect();
    let (base_median_ppm, base_sigma_ppm) = residual_stats(&base_ppm);
    let mut keep: Vec<bool> = vec![true; points.len()];
    let mut result: Option<CalibrationFit> = None;
    for _ in 0..3 {
        let used: Vec<&CalibrationPoint> = points
            .iter()
            .zip(&keep)
            .filter(|(_, &k)| k)
            .map(|(p, _)| p)
            .collect();
        if used.len() < MIN_CALIBRATION_POINTS {
            break;
        }
        let n = features(used[0]).len();
        let mut ata = vec![vec![0.0; n]; n];
        let mut atb = vec![0.0; n];
        for p in &used {
            let f = features(p);
            let y = p.mz.sqrt();
            for i in 0..n {
                atb[i] += f[i] * y;
                for j in 0..n {
                    ata[i][j] += f[i] * f[j];
                }
            }
        }
        let coef = solve(ata, atb)?;
        let mut temp = [0.0; 2];
        let mut k = 4;
        for (t, &on) in temp.iter_mut().zip(&use_temp) {
            if on {
                *t = coef[k];
                k += 1;
            }
        }
        let calibration = MzCalibration {
            scale,
            poly: [coef[0], coef[1], coef[2], coef[3]],
            temp,
            t_ref,
        };
        let ppm: Vec<f64> = points
            .iter()
            .map(|p| (calibration.mz(p.tof, p.t) - p.mz) / p.mz * 1e6)
            .collect();
        let used_ppm: Vec<f64> = ppm
            .iter()
            .zip(&keep)
            .filter(|(_, &k)| k)
            .map(|(r, _)| *r)
            .collect();
        let (median_ppm, sigma_ppm) = residual_stats(&used_ppm);
        result = Some(CalibrationFit {
            calibration,
            n_used: used_ppm.len(),
            median_ppm,
            sigma_ppm,
            base_median_ppm,
            base_sigma_ppm,
        });
        let limit = OUTLIER_SIGMAS * sigma_ppm.max(1e-3);
        keep = ppm
            .iter()
            .map(|r| (r - median_ppm).abs() <= limit)
            .collect();
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn summing_merges_equal_tof_indices_across_scans() {
        let tof = [10, 12, 10, 11, 12];
        let inten = [1, 2, 3, 4, 5];
        assert_eq!(
            sum_over_mobility(&tof, &inten),
            vec![(10, 4), (11, 4), (12, 7)]
        );
    }

    #[test]
    fn one_peak_is_one_weighted_centroid() {
        let profile = vec![(100, 1), (101, 3), (102, 5), (103, 3), (104, 1)];
        let c = centroid_profile(&profile);
        assert_eq!(c.len(), 1);
        assert!((c[0].tof - 102.0).abs() < 1e-12);
        assert_eq!(c[0].intensity, 13.0);
    }

    #[test]
    fn a_deep_valley_splits_and_a_shallow_dip_does_not() {
        // Two maxima (10, 8) with a valley of 1: two peaks.
        let deep = vec![(1, 2), (2, 10), (3, 1), (4, 8), (5, 2)];
        assert_eq!(centroid_profile(&deep).len(), 2);
        // A dip to 9 between 10 and 10: one peak.
        let shallow = vec![(1, 2), (2, 10), (3, 9), (4, 10), (5, 2)];
        assert_eq!(centroid_profile(&shallow).len(), 1);
    }

    #[test]
    fn a_gap_wider_than_the_limit_separates_peaks() {
        let profile = vec![(1, 5), (2, 5), (3 + MAX_TOF_GAP + 1, 5)];
        assert_eq!(centroid_profile(&profile).len(), 2);
        let close = vec![(1, 5), (1 + MAX_TOF_GAP, 5)];
        assert_eq!(centroid_profile(&close).len(), 1);
    }

    /// A synthetic instrument: the vendor calibration is a cubic with a
    /// temperature term; the acquisition-range model is linear and off.
    fn truth(tof: f64, t: [f64; 2]) -> f64 {
        let x = tof / 400_000.0;
        let root = 10.0 + 31.0 * x + 0.02 * x * x - 0.01 * x * x * x + 2e-4 * (t[1] - 27.0);
        root * root
    }

    #[test]
    fn the_fit_recovers_a_cubic_temperature_dependent_calibration() {
        let base = MzCalibration::linear(10.0, 31.0 / 400_000.0, 400_000.0);
        let mut points = Vec::new();
        for i in 0..2000 {
            let tof = 20_000.0 + (i as f64) * 180.0;
            let t = [25.7, 27.0 + 0.3 * ((i % 7) as f64 / 6.0)];
            points.push(CalibrationPoint {
                tof,
                mz: truth(tof, t),
                t,
            });
        }
        let fit = fit_calibration(&points, &base).expect("enough pairs");
        assert!(fit.sigma_ppm < 0.01, "{fit:?}");
        assert!(fit.median_ppm.abs() < 0.01, "{fit:?}");
        assert!(fit.base_sigma_ppm > fit.sigma_ppm);
        for &(tof, t) in &[(50_000.0, [25.7, 27.1]), (300_000.0, [25.7, 27.25])] {
            let ppm = (fit.calibration.mz(tof, t) - truth(tof, t)) / truth(tof, t) * 1e6;
            assert!(ppm.abs() < 0.01, "{ppm}");
        }
    }

    #[test]
    fn outliers_are_left_out_of_the_refit() {
        let base = MzCalibration::linear(10.0, 31.0 / 400_000.0, 400_000.0);
        let mut points: Vec<CalibrationPoint> = (0..1000)
            .map(|i| {
                let tof = 20_000.0 + (i as f64) * 350.0;
                let t = [25.7, 27.0];
                CalibrationPoint {
                    tof,
                    mz: truth(tof, t),
                    t,
                }
            })
            .collect();
        // 5 % wrong matches 30 ppm off.
        for p in points.iter_mut().step_by(20) {
            p.mz *= 1.0 + 30e-6;
        }
        let fit = fit_calibration(&points, &base).unwrap();
        assert!(fit.median_ppm.abs() < 0.05, "{fit:?}");
        assert!(fit.n_used <= 950 + 1, "{fit:?}");
    }

    #[test]
    fn too_few_pairs_give_no_fit() {
        let base = MzCalibration::linear(10.0, 1e-4, 1.0);
        let points = vec![
            CalibrationPoint {
                tof: 1.0,
                mz: 100.0,
                t: [0.0; 2]
            };
            MIN_CALIBRATION_POINTS - 1
        ];
        assert!(fit_calibration(&points, &base).is_none());
    }

    #[test]
    fn matching_picks_the_nearest_centroid_inside_the_window() {
        let base = MzCalibration::linear(10.0, 31.0 / 400_000.0, 400_000.0);
        let t = [0.0; 2];
        let target = base.mz(100_000.0, t);
        let centroids = vec![
            TofCentroid {
                tof: 99_990.0,
                intensity: 1.0,
            },
            TofCentroid {
                tof: 100_000.4,
                intensity: 1.0,
            },
            TofCentroid {
                tof: 100_020.0,
                intensity: 1.0,
            },
        ];
        let mut out = Vec::new();
        match_precursors(&centroids, &[target], t, &base, 0.0, 20.0, &mut out);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].tof, 100_000.4);
        out.clear();
        // Nothing within 0.001 ppm of an offset of +500 ppm.
        match_precursors(&centroids, &[target], t, &base, 500.0, 0.001, &mut out);
        assert!(out.is_empty());
    }
}
