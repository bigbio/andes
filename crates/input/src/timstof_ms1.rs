//! MS1 survey frames of a Bruker timsTOF `.d`: summed over ion mobility,
//! centroided along the time-of-flight axis, and converted to m/z with a
//! per-run calibration learned from the instrument's own precursor records.
//!
//! A TIMS MS1 frame holds one TOF spectrum per mobility scan. Quantification
//! and the MS1-linked search paths read survey scans as one centroid list per
//! retention time, so each frame is summed over its scans (per TOF index),
//! split into peaks at the valleys of the lightly smoothed summed profile,
//! and each peak is reported at the intensity-weighted TOF position of its
//! half-height core.
//!
//! The TOF → m/z conversion that `timsrust` derives from the acquisition
//! range (`sqrt(m/z)` linear in the TOF index) is an approximation of the
//! vendor calibration. The acquisition software records, for every PASEF
//! precursor, its monoisotopic m/z (`Precursors.MonoisotopicMz`) computed with
//! the full instrument calibration; it agrees with the theoretical m/z of
//! identified precursors to about 1 ppm, while the m/z of the largest isotope
//! peak recorded beside it (`LargestPeakMz`, used only when the monoisotopic
//! peak was not determined) sits several ppm off. [`fit_calibration`] matches
//! those m/z values to the centroids of their parent MS1 frames and fits
//! `sqrt(m/z)` as a cubic in the TOF index plus a linear term in each of the
//! two frame temperatures the TDF records (`Frames.T1`, `Frames.T2`), which is
//! what the vendor calibration compensates for. Each vendor m/z is paired
//! with the most intense centroid in a window, in two passes: a wide window
//! whose pairs give the run's offset as the mode of their offsets, then a
//! narrow window around that offset. The fit is refined by discarding residuals beyond four robust
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

/// Two maxima of the smoothed profile separated by a valley deeper than this
/// fraction of the lower maximum are two peaks; a shallower dip is noise on
/// one peak.
const VALLEY_SPLIT_FRACTION: f64 = 0.5;

/// Binomial smoothing kernel (standard deviation one TOF bin). An MS1 peak
/// spans several TOF bins (its full width at half maximum is 3 to 6 bins over
/// the m/z range), while the counts of a summed frame fluctuate from bin to
/// bin; peaks are found on the smoothed profile so that this noise does not
/// split one peak into several.
const SMOOTHING: [f64; 5] = [1.0 / 16.0, 4.0 / 16.0, 6.0 / 16.0, 4.0 / 16.0, 1.0 / 16.0];

/// The part of a peak its position is taken from: the points whose smoothed
/// intensity is at least this fraction of the peak's smoothed maximum (the
/// half-height core, where the shape is symmetric and the tails of
/// neighbouring peaks weigh least).
const CORE_FRACTION: f64 = 0.5;

/// Centroid a summed TOF profile (sorted by TOF index). Runs of points no
/// more than [`MAX_TOF_GAP`] apart are laid out on a dense TOF grid, smoothed
/// with [`SMOOTHING`], and split at the minima of the smoothed profile;
/// neighbouring segments are merged again when the smoothed valley between
/// them is shallower than [`VALLEY_SPLIT_FRACTION`] of the lower segment
/// maximum. Each peak is reported at the intensity-weighted TOF index of its
/// core (see [`CORE_FRACTION`]) with the summed intensity of the whole peak.
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

/// Centroid one run of points no more than [`MAX_TOF_GAP`] apart (see
/// [`centroid_profile`]).
fn centroid_run(run: &[(u32, u64)], out: &mut Vec<TofCentroid>) {
    let first = run[0].0;
    let n = (run[run.len() - 1].0 - first) as usize + 1;
    let mut raw = vec![0.0f64; n];
    for &(t, i) in run {
        raw[(t - first) as usize] = i as f64;
    }
    let half = SMOOTHING.len() / 2;
    let smooth: Vec<f64> = (0..n)
        .map(|i| {
            SMOOTHING
                .iter()
                .enumerate()
                .filter_map(|(k, w)| {
                    let j = (i + k).checked_sub(half)?;
                    raw.get(j).map(|v| v * w)
                })
                .sum()
        })
        .collect();
    // Segments [lo, hi) split at the strict minima of the smoothed profile,
    // each with its smoothed maximum and the smoothed valley in front of it.
    let mut segments: Vec<(usize, usize, f64, f64)> = Vec::new();
    let mut lo = 0;
    let mut valley = 0.0;
    for i in 1..n.saturating_sub(1) {
        if smooth[i] < smooth[i - 1] && smooth[i] <= smooth[i + 1] {
            let max = smooth[lo..i].iter().copied().fold(0.0, f64::max);
            segments.push((lo, i, max, valley));
            lo = i;
            valley = smooth[i];
        }
    }
    let max = smooth[lo..].iter().copied().fold(0.0, f64::max);
    segments.push((lo, n, max, valley));
    let mut merged: Vec<(usize, usize, f64)> = Vec::with_capacity(segments.len());
    for (lo, hi, max, valley) in segments {
        match merged.last_mut() {
            Some(last) if valley >= VALLEY_SPLIT_FRACTION * last.2.min(max) => {
                last.1 = hi;
                last.2 = last.2.max(max);
            }
            _ => merged.push((lo, hi, max)),
        }
    }
    for (lo, hi, max) in merged {
        let total: f64 = raw[lo..hi].iter().sum();
        if total <= 0.0 {
            continue;
        }
        let (mut w, mut wt) = (0.0, 0.0);
        for i in lo..hi {
            if smooth[i] >= CORE_FRACTION * max {
                w += raw[i];
                wt += raw[i] * (first as f64 + i as f64);
            }
        }
        if w <= 0.0 {
            w = total;
            wt = (lo..hi).map(|i| raw[i] * (first as f64 + i as f64)).sum();
        }
        out.push(TofCentroid {
            tof: wt / w,
            intensity: total,
        });
    }
}

/// Conditions a centroid's m/z depends on besides its TOF index: the two
/// frame temperatures (`T1`, `T2`) and the natural log of the centroid's
/// summed intensity. Intense peaks sit at a slightly later TOF than weak ones
/// of the same m/z (several ppm across the intensity range of a run), which
/// the vendor conversion corrects and the TOF polynomial alone cannot.
pub(crate) type Conditions = [f64; 3];

/// The conditions of a centroid of `intensity` in a frame recorded at
/// temperatures `t`.
pub(crate) fn conditions(t: [f64; 2], intensity: f64) -> Conditions {
    [t[0], t[1], intensity.max(1.0).ln()]
}

/// TOF → m/z conversion: `sqrt(m/z)` as a polynomial in the scaled TOF index
/// plus a linear term in each of the [`Conditions`]. With all condition
/// coefficients zero it is the acquisition-range model
/// (`sqrt(m/z) = a + b * tof`).
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct MzCalibration {
    /// TOF indices are divided by this before the polynomial (numerical
    /// conditioning).
    scale: f64,
    /// Coefficients of `1, x, x^2, x^3` (x = tof / scale).
    poly: [f64; 4],
    /// Coefficients of the conditions, each centred on `c_ref`.
    cond: [f64; 3],
    c_ref: Conditions,
}

impl MzCalibration {
    /// The acquisition-range model `sqrt(m/z) = intercept + slope * tof`.
    pub fn linear(intercept: f64, slope: f64, scale: f64) -> Self {
        Self {
            scale,
            poly: [intercept, slope * scale, 0.0, 0.0],
            cond: [0.0; 3],
            c_ref: [0.0; 3],
        }
    }

    /// m/z of a (fractional) TOF index under conditions `c`.
    pub fn mz(&self, tof: f64, c: Conditions) -> f64 {
        let x = tof / self.scale;
        let p = &self.poly;
        let terms: f64 = (0..3).map(|k| self.cond[k] * (c[k] - self.c_ref[k])).sum();
        let root = p[0] + x * (p[1] + x * (p[2] + x * p[3])) + terms;
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
    /// Vendor-calibrated monoisotopic m/z of the precursor.
    pub mz: f64,
    /// Parent frame temperatures and the centroid's log intensity.
    pub t: Conditions,
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

/// Pair each vendor m/z of a frame with the most intense centroid within
/// `offset_ppm ± window_ppm` (positions under the linear model `base`). A
/// selected precursor is among the most intense ions of its frame, so among
/// the centroids near its m/z the most intense is its counterpart; the
/// nearest one is as often a neighbouring peak. `centroids` are sorted by TOF.
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
        // The window in the TOF domain; the base model is monotonic in TOF,
        // and the condition terms of `base` are zero.
        let at = |ppm: f64| base.linear_tof(mz * (1.0 + ppm * 1e-6));
        let (lo, hi) = (at(offset_ppm - window_ppm), at(offset_ppm + window_ppm));
        let i = centroids.partition_point(|c| c.tof < lo);
        let j = centroids.partition_point(|c| c.tof <= hi);
        let best = centroids[i..j.max(i)]
            .iter()
            .max_by(|a, b| a.intensity.total_cmp(&b.intensity));
        if let Some(c) = best {
            out.push(CalibrationPoint {
                tof: c.tof,
                mz,
                t: conditions(t, c.intensity),
            });
        }
    }
}

/// Width (ppm) of the offset histogram bins of [`mode_offset_ppm`].
const OFFSET_BIN_PPM: f64 = 1.0;

/// The run's offset (ppm) of the linear model from the vendor m/z: the mode
/// of the matched pairs' offsets (histogram of [`OFFSET_BIN_PPM`] bins,
/// smoothed over three bins). Correct pairs pile up at the offset while pairs
/// with a neighbouring peak spread over the whole window, so the mode finds
/// the offset where a median would be pulled toward the window centre.
pub(crate) fn mode_offset_ppm(points: &[CalibrationPoint], base: &MzCalibration) -> Option<f64> {
    if points.is_empty() {
        return None;
    }
    let half = (WIDE_WINDOW_PPM / OFFSET_BIN_PPM).ceil() as i64;
    let mut hist = vec![0.0f64; (2 * half + 1) as usize];
    for p in points {
        let ppm = (base.mz(p.tof, p.t) - p.mz) / p.mz * 1e6;
        let k = (ppm / OFFSET_BIN_PPM).round() as i64 + half;
        if let Some(h) = usize::try_from(k).ok().and_then(|k| hist.get_mut(k)) {
            *h += 1.0;
        }
    }
    let smooth = |k: usize| {
        hist[k.saturating_sub(1)..(k + 2).min(hist.len())]
            .iter()
            .sum::<f64>()
    };
    let k = (0..hist.len()).max_by(|&a, &b| smooth(a).total_cmp(&smooth(b)))?;
    // Refine within the peak: the mean offset of the pairs within one bin.
    let centre = (k as i64 - half) as f64 * OFFSET_BIN_PPM;
    let near: Vec<f64> = points
        .iter()
        .map(|p| (base.mz(p.tof, p.t) - p.mz) / p.mz * 1e6)
        .filter(|ppm| (ppm - centre).abs() <= 1.5 * OFFSET_BIN_PPM)
        .collect();
    Some(near.iter().sum::<f64>() / near.len().max(1) as f64)
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

/// Spread (standard deviation) below which a condition is treated as
/// constant over the run and gets no term.
const MIN_CONDITION_SPREAD: f64 = 1e-3;

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
    let t_ref = [mean_t(0), mean_t(1), mean_t(2)];
    let spread = |k: usize| {
        (points
            .iter()
            .map(|p| (p.t[k] - t_ref[k]).powi(2))
            .sum::<f64>()
            / points.len() as f64)
            .sqrt()
    };
    let use_temp = [0, 1, 2].map(|k| spread(k) > MIN_CONDITION_SPREAD);
    let features = |p: &CalibrationPoint| -> Vec<f64> {
        let x = p.tof / scale;
        let mut f = vec![1.0, x, x * x, x * x * x];
        for k in 0..3 {
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
        let mut temp = [0.0; 3];
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
            cond: temp,
            c_ref: t_ref,
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

    /// A Gaussian peak of `height` at TOF `centre` (standard deviation
    /// `sigma` bins) over `lo..hi`, with a deterministic count jitter of
    /// `jitter` (fraction of the value, alternating sign).
    fn gaussian(
        lo: u32,
        hi: u32,
        centre: f64,
        sigma: f64,
        height: f64,
        jitter: f64,
    ) -> Vec<(u32, u64)> {
        (lo..hi)
            .map(|t| {
                let x = (t as f64 - centre) / sigma;
                let sign = if t % 2 == 0 { 1.0 } else { -1.0 };
                let v = height * (-0.5 * x * x).exp() * (1.0 + sign * jitter);
                (t, v.round() as u64)
            })
            .filter(|p| p.1 > 0)
            .collect()
    }

    fn add(a: &[(u32, u64)], b: &[(u32, u64)]) -> Vec<(u32, u64)> {
        let mut m = std::collections::BTreeMap::new();
        for &(t, i) in a.iter().chain(b) {
            *m.entry(t).or_insert(0) += i;
        }
        m.into_iter().collect()
    }

    #[test]
    fn bin_to_bin_noise_does_not_split_one_peak() {
        // A 2-bin sigma peak with 30 % alternating count noise: one centroid
        // within a quarter bin (about 1 ppm) of the true centre.
        let p = gaussian(990, 1011, 1000.3, 2.0, 40.0, 0.3);
        let c = centroid_profile(&p);
        assert_eq!(c.len(), 1, "{c:?}");
        assert!((c[0].tof - 1000.3).abs() < 0.25, "{c:?}");
    }

    #[test]
    fn resolved_peaks_split_and_an_unresolved_shoulder_does_not() {
        // Two peaks 10 bins apart (5 sigma): two centroids at their centres.
        let p = add(
            &gaussian(990, 1011, 1000.0, 2.0, 100.0, 0.0),
            &gaussian(1000, 1021, 1010.0, 2.0, 60.0, 0.0),
        );
        let c = centroid_profile(&p);
        assert_eq!(c.len(), 2, "{c:?}");
        assert!(
            (c[0].tof - 1000.0).abs() < 0.2 && (c[1].tof - 1010.0).abs() < 0.2,
            "{c:?}"
        );
        // 3 bins apart: one peak.
        let q = add(
            &gaussian(990, 1011, 1000.0, 2.0, 100.0, 0.0),
            &gaussian(993, 1014, 1003.0, 2.0, 60.0, 0.0),
        );
        assert_eq!(centroid_profile(&q).len(), 1);
    }

    #[test]
    fn a_gap_wider_than_the_limit_separates_peaks() {
        let profile = vec![(1, 5), (2, 5), (3 + MAX_TOF_GAP + 1, 5)];
        assert_eq!(centroid_profile(&profile).len(), 2);
        let close = vec![(1, 5), (1 + MAX_TOF_GAP, 5)];
        assert_eq!(centroid_profile(&close).len(), 1);
    }

    /// A synthetic instrument: the vendor calibration is a cubic with a
    /// temperature and an intensity term; the acquisition-range model is
    /// linear and off.
    fn truth(tof: f64, t: Conditions) -> f64 {
        let x = tof / 400_000.0;
        let root = 10.0 + 31.0 * x + 0.02 * x * x - 0.01 * x * x * x
            + 2e-4 * (t[1] - 27.0)
            + 4e-5 * (t[2] - 5.0);
        root * root
    }

    #[test]
    fn the_fit_recovers_a_cubic_temperature_dependent_calibration() {
        let base = MzCalibration::linear(10.0, 31.0 / 400_000.0, 400_000.0);
        let mut points = Vec::new();
        for i in 0..2000 {
            let tof = 20_000.0 + (i as f64) * 180.0;
            let t = [
                25.7,
                27.0 + 0.3 * ((i % 7) as f64 / 6.0),
                3.0 + 0.5 * (i % 11) as f64,
            ];
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
        for &(tof, t) in &[
            (50_000.0, [25.7, 27.1, 4.0]),
            (300_000.0, [25.7, 27.25, 7.5]),
        ] {
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
                let t = [25.7, 27.0, 5.0];
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
                t: [0.0; 3]
            };
            MIN_CALIBRATION_POINTS - 1
        ];
        assert!(fit_calibration(&points, &base).is_none());
    }

    #[test]
    fn matching_picks_the_most_intense_centroid_inside_the_window() {
        let base = MzCalibration::linear(10.0, 31.0 / 400_000.0, 400_000.0);
        let t = [0.0; 2];
        let target = base.mz(100_000.0, [0.0; 3]);
        let c = |tof: f64, intensity: f64| TofCentroid { tof, intensity };
        // Nearest (100_000.4) is weak; the intense one 1 bin away is the peak.
        let centroids = vec![
            c(99_990.0, 50.0),
            c(100_000.4, 1.0),
            c(100_001.0, 20.0),
            c(100_200.0, 99.0),
        ];
        let mut out = Vec::new();
        match_precursors(&centroids, &[target], t, &base, 0.0, 20.0, &mut out);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].tof, 100_001.0);
        out.clear();
        // Nothing within 0.001 ppm of an offset of +500 ppm.
        match_precursors(&centroids, &[target], t, &base, 500.0, 0.001, &mut out);
        assert!(out.is_empty());
    }

    #[test]
    fn the_offset_is_the_mode_not_the_median() {
        // 300 pairs at +8 ppm, 700 spread uniformly over -150..150 ppm.
        let base = MzCalibration::linear(10.0, 31.0 / 400_000.0, 400_000.0);
        let t = [0.0; 3];
        let mut points = Vec::new();
        for i in 0..1000 {
            let tof = 50_000.0 + i as f64 * 200.0;
            let m = base.mz(tof, t);
            let ppm = if i < 300 {
                8.0
            } else {
                -150.0 + 300.0 * ((i * 37 % 700) as f64 / 700.0)
            };
            points.push(CalibrationPoint {
                tof,
                mz: m / (1.0 + ppm * 1e-6),
                t,
            });
        }
        let off = mode_offset_ppm(&points, &base).unwrap();
        assert!((off - 8.0).abs() < 0.2, "{off}");
    }
}
