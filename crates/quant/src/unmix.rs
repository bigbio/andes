//! Reporter-ion interference correction for isobaric MS2 scans.
//!
//! On high-resolution MS2 the reporter ions of a scan come from every ion the
//! quadrupole let through, not only from the identified precursor. Each scan
//! gets an isolation-window record from the MS1 scans that border it in time
//! (the ion current inside the window and the share of it on the precursor's
//! isotope envelope), and the run's summed reporter intensities are explained
//! by one regression:
//!
//! ```text
//! S_s = a_z · PIC_s · f_s  +  b · PIC_s · (1 − f_s)  +  c · H_s
//! ```
//!
//! * `S_s` is the scan's summed reporter intensity,
//! * `PIC_s` its peptide ion current (MS2 ion current above the reporter
//!   region),
//! * `f_s` the precursor's share of the window's MS1 ion current,
//! * `H_s` the ion current hidden below the MS1 detection limit inside the
//!   window: detection limit × centroid density around the precursor ×
//!   window width,
//! * `a_z` (per precursor charge class), `b` and `c` are non-negative and fitted
//!   per run by robust (bisquare-reweighted) non-negative least squares on
//!   relative residuals.
//!
//! The coefficients are the run's scale factors between MS1 ion current and
//! reporter yield, so no ground truth is needed. A scan's interference
//! fraction is the share of its predicted reporter signal that does not come
//! from the precursor, `1 − a_z·PIC·f / predicted`, capped at
//! [`MAX_INTERFERENCE`]. The interference is removed with the run's mean
//! reporter profile (each channel's share of all reporter signal): the
//! co-isolated ions are, on average, the run's average peptide.
//!
//! Joint non-negative least squares over scans linked by shared precursors
//! (every identified co-isolated precursor subtracted with its own reporter
//! profile) was measured against this and recovered less of the known
//! ratios, so only the regression is implemented.

use model::mass::ISOTOPE;

use crate::nnls::nnls;

/// Largest interference fraction removed from one scan.
pub const MAX_INTERFERENCE: f64 = 0.8;
/// A run needs at least this many usable scans to fit its model.
pub const MIN_FIT_SCANS: usize = 100;
/// Isotope peaks of the precursor are matched within this tolerance.
const ENVELOPE_PPM: f64 = 10.0;
/// Isotope positions of the precursor counted inside the window, relative to
/// the monoisotopic peak.
const ENVELOPE_ISOTOPES: std::ops::RangeInclusive<i32> = -1..=8;
/// Half width (Th) of the m/z range whose MS1 centroid density scales the
/// hidden ion current.
const DENSITY_HALF_WIDTH: f64 = 5.0;
/// Bisquare tuning constant and reweighting rounds of the robust fit.
const BISQUARE_K: f64 = 4.685;
const ROBUST_ROUNDS: usize = 8;
/// Regression columns: three precursor charge classes, co-isolated, hidden.
const N_COEF: usize = 5;
/// Floor of the robust scale of the relative residuals.
const MIN_MAD: f64 = 1e-6;

/// One quantified MS2 scan as the interference model sees it.
#[derive(Debug, Clone, PartialEq)]
pub struct ScanInput {
    /// Charge of the identified precursor.
    pub charge: u8,
    /// Reporter intensities, impurity-corrected when a sheet was applied.
    pub reporters: Vec<f32>,
    /// MS2 ion current above the reporter region.
    pub pic: f64,
    /// The isolation-window record, `None` without MS1 around the scan.
    pub window: Option<WindowRecord>,
}

/// The MS1 view of one isolation window at the time of the MS2 scan.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WindowRecord {
    /// Ion current inside the window, interpolated between the bordering MS1
    /// scans.
    pub total: f64,
    /// The part of `total` on the precursor's isotope envelope.
    pub precursor: f64,
    /// Ion current estimated below the MS1 detection limit inside the window.
    pub hidden: f64,
}

impl WindowRecord {
    /// The precursor's share of the window's visible ion current (1 when the
    /// window holds no centroid).
    pub fn precursor_fraction(&self) -> f64 {
        if self.total > 0.0 {
            (self.precursor / self.total).clamp(0.0, 1.0)
        } else {
            1.0
        }
    }
}

/// One bordering MS1 scan: its retention time, m/z-sorted centroids and
/// detection limit.
#[derive(Debug, Clone, Copy)]
pub struct Survey<'a> {
    pub rt: f64,
    pub peaks: &'a [(f64, f32)],
    pub detection_limit: f32,
}

/// Window intensity and precursor-envelope intensity of one MS1 scan.
fn window_sums(peaks: &[(f64, f32)], lo: f64, hi: f64, mono_mz: f64, charge: u8) -> (f64, f64) {
    let start = peaks.partition_point(|p| p.0 < lo);
    let end = start + peaks[start..].partition_point(|p| p.0 <= hi);
    let window = &peaks[start..end];
    let total: f64 = window.iter().map(|p| p.1 as f64).sum();
    let z = charge.max(1) as f64;
    let mut envelope = 0.0;
    for i in ENVELOPE_ISOTOPES {
        let expected = mono_mz + i as f64 * ISOTOPE / z;
        if expected < lo || expected > hi {
            continue;
        }
        let tol = expected * ENVELOPE_PPM * 1e-6;
        let a = window.partition_point(|p| p.0 < expected - tol);
        let best = window[a..]
            .iter()
            .take_while(|p| p.0 <= expected + tol)
            .map(|p| p.1 as f64)
            .fold(0.0, f64::max);
        envelope += best;
    }
    (total, envelope.min(total))
}

/// Centroids per Th within [`DENSITY_HALF_WIDTH`] of `centre`.
fn centroid_density(peaks: &[(f64, f32)], centre: f64) -> f64 {
    let lo = peaks.partition_point(|p| p.0 < centre - DENSITY_HALF_WIDTH);
    let hi = peaks.partition_point(|p| p.0 <= centre + DENSITY_HALF_WIDTH);
    (hi - lo) as f64 / (2.0 * DENSITY_HALF_WIDTH)
}

/// The isolation-window record of an MS2 scan at `rt` isolating
/// `[lo, hi]`, from the MS1 scans before (`prev`) and after (`next`) it.
/// Intensities are interpolated linearly in retention time; with one
/// bordering scan, that scan is used. The precursor envelope is the isotope
/// ladder of `mono_mz` at `charge`. The hidden ion current is the bordering
/// scans' detection limit times the centroid density around the window times
/// its width. `None` without a bordering scan or for an empty window range.
pub fn window_record(
    prev: Option<Survey<'_>>,
    next: Option<Survey<'_>>,
    rt: f64,
    lo: f64,
    hi: f64,
    mono_mz: f64,
    charge: u8,
) -> Option<WindowRecord> {
    if hi.partial_cmp(&lo) != Some(std::cmp::Ordering::Greater) {
        return None;
    }
    let weight_prev = match (prev, next) {
        (Some(p), Some(n)) if n.rt > p.rt => ((n.rt - rt) / (n.rt - p.rt)).clamp(0.0, 1.0),
        (Some(_), Some(_)) => 0.5,
        (Some(_), None) => 1.0,
        (None, Some(_)) => 0.0,
        (None, None) => return None,
    };
    let centre = 0.5 * (lo + hi);
    let mut rec = WindowRecord {
        total: 0.0,
        precursor: 0.0,
        hidden: 0.0,
    };
    for (scan, w) in [(prev, weight_prev), (next, 1.0 - weight_prev)] {
        let Some(s) = scan else { continue };
        if w <= 0.0 {
            continue;
        }
        let (total, precursor) = window_sums(s.peaks, lo, hi, mono_mz, charge);
        rec.total += w * total;
        rec.precursor += w * precursor;
        rec.hidden += w * s.detection_limit as f64 * centroid_density(s.peaks, centre) * (hi - lo);
    }
    Some(rec)
}

/// A run's fitted interference model.
#[derive(Debug, Clone, PartialEq)]
pub struct InterferenceModel {
    /// Reporter yield per unit of precursor ion current, by charge class
    /// (1–2, 3, 4 and more).
    pub precursor: [f64; 3],
    /// Reporter yield per unit of co-isolated (non-precursor) ion current.
    pub co_isolated: f64,
    /// Reporter yield per unit of hidden ion current.
    pub hidden: f64,
    /// Scans the model was fitted on.
    pub n: usize,
    /// Share of the variance of log summed reporter intensity the model
    /// explains on those scans.
    pub r2_log: f64,
}

fn charge_class(z: u8) -> usize {
    match z {
        0..=2 => 0,
        3 => 1,
        _ => 2,
    }
}

/// Regression row of a scan, `None` when the scan cannot enter the model.
fn design_row(s: &ScanInput) -> Option<[f64; N_COEF]> {
    let w = s.window?;
    if !(s.pic.is_finite() && s.pic > 0.0) {
        return None;
    }
    let f = w.precursor_fraction();
    let mut row = [0.0; N_COEF];
    row[charge_class(s.charge)] = s.pic * f;
    row[3] = s.pic * (1.0 - f);
    row[4] = if w.hidden.is_finite() { w.hidden } else { 0.0 };
    Some(row)
}

fn summed(r: &[f32]) -> f64 {
    r.iter().map(|&v| (v as f64).max(0.0)).sum()
}

fn median(mut v: Vec<f64>) -> f64 {
    if v.is_empty() {
        return 0.0;
    }
    let mid = v.len() / 2;
    *v.select_nth_unstable_by(mid, f64::total_cmp).1
}

impl InterferenceModel {
    /// Fit the model on a run's scans. `None` with fewer than
    /// [`MIN_FIT_SCANS`] usable scans or when every coefficient fits to 0.
    pub fn fit(scans: &[ScanInput]) -> Option<Self> {
        let (rows, y): (Vec<[f64; N_COEF]>, Vec<f64>) = scans
            .iter()
            .filter_map(|s| {
                let total = summed(&s.reporters);
                (total > 0.0).then_some(())?;
                Some((design_row(s)?, total))
            })
            .unzip();
        let m = rows.len();
        if m < MIN_FIT_SCANS {
            return None;
        }
        // Relative residuals: each scan counts by its relative error, so the
        // most intense scans do not dominate the fit.
        let mut weight = vec![1.0f64; m];
        let mut beta = vec![0.0f64; N_COEF];
        let mut a = vec![0.0f64; m * N_COEF];
        let mut b = vec![0.0f64; m];
        for _ in 0..ROBUST_ROUNDS {
            for i in 0..m {
                let s = weight[i] / y[i];
                for j in 0..N_COEF {
                    a[i * N_COEF + j] = rows[i][j] * s;
                }
                b[i] = weight[i];
            }
            beta = nnls(&a, m, N_COEF, &b, 50 * N_COEF);
            let resid: Vec<f64> = (0..m)
                .map(|i| (y[i] - dot(&rows[i], &beta)) / y[i])
                .collect();
            let centre = median(resid.clone());
            let mad =
                (1.4826 * median(resid.iter().map(|r| (r - centre).abs()).collect())).max(MIN_MAD);
            for (w, r) in weight.iter_mut().zip(&resid) {
                let u = r / (BISQUARE_K * mad);
                *w = if u.abs() < 1.0 { 1.0 - u * u } else { 0.0 };
            }
        }
        if beta.iter().all(|&v| v <= 0.0) {
            return None;
        }
        let logs: Vec<(f64, f64)> = rows
            .iter()
            .zip(&y)
            .map(|(r, &t)| (t.ln(), dot(r, &beta).max(f64::MIN_POSITIVE).ln()))
            .collect();
        let mean = logs.iter().map(|l| l.0).sum::<f64>() / m as f64;
        let ss_tot: f64 = logs.iter().map(|l| (l.0 - mean).powi(2)).sum();
        let ss_res: f64 = logs.iter().map(|l| (l.0 - l.1).powi(2)).sum();
        Some(Self {
            precursor: [beta[0], beta[1], beta[2]],
            co_isolated: beta[3],
            hidden: beta[4],
            n: m,
            r2_log: if ss_tot > 0.0 {
                1.0 - ss_res / ss_tot
            } else {
                0.0
            },
        })
    }

    /// The scan's interference fraction in `0..=`[`MAX_INTERFERENCE`], `None`
    /// when the scan cannot be modelled (no window record, no ion current).
    pub fn interference(&self, scan: &ScanInput) -> Option<f64> {
        let row = design_row(scan)?;
        let beta = self.coefficients();
        let predicted = dot(&row, &beta);
        if predicted <= 0.0 {
            return None;
        }
        let own = row[charge_class(scan.charge)] * beta[charge_class(scan.charge)];
        Some((1.0 - own / predicted).clamp(0.0, MAX_INTERFERENCE))
    }

    fn coefficients(&self) -> [f64; N_COEF] {
        [
            self.precursor[0],
            self.precursor[1],
            self.precursor[2],
            self.co_isolated,
            self.hidden,
        ]
    }
}

fn dot(row: &[f64; N_COEF], beta: &[f64]) -> f64 {
    row.iter().zip(beta).map(|(a, b)| a * b).sum()
}

/// Each channel's share of all reporter signal of the run's scans.
pub fn background_profile(scans: &[ScanInput]) -> Option<Vec<f64>> {
    let n = scans.first()?.reporters.len();
    let mut sum = vec![0.0f64; n];
    for s in scans.iter().filter(|s| s.reporters.len() == n) {
        for (acc, &v) in sum.iter_mut().zip(&s.reporters) {
            *acc += (v as f64).max(0.0);
        }
    }
    let total: f64 = sum.iter().sum();
    (total > 0.0).then(|| sum.iter().map(|v| v / total).collect())
}

/// One scan after interference removal.
#[derive(Debug, Clone, PartialEq)]
pub struct Unmixed {
    /// Estimated interference fraction of the scan's reporter signal.
    pub interference: f32,
    /// Reporter intensities with the interference removed (floored at 0).
    pub reporters: Vec<f32>,
}

/// Fit the run's model and remove each scan's interference. Returns the
/// model (`None` when the run cannot be fitted, and then every scan is
/// `None`) and one entry per input scan in order.
pub fn unmix_run(scans: &[ScanInput]) -> (Option<InterferenceModel>, Vec<Option<Unmixed>>) {
    let none = || vec![None; scans.len()];
    let Some(model) = InterferenceModel::fit(scans) else {
        return (None, none());
    };
    let Some(profile) = background_profile(scans) else {
        return (None, none());
    };
    let out = scans
        .iter()
        .map(|s| {
            let eil = model.interference(s)?;
            let total = summed(&s.reporters);
            let reporters = s
                .reporters
                .iter()
                .zip(&profile)
                .map(|(&v, &p)| ((v as f64) - eil * total * p).max(0.0) as f32)
                .collect();
            Some(Unmixed {
                interference: eil as f32,
                reporters,
            })
        })
        .collect();
    (Some(model), out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn survey(rt: f64, peaks: &[(f64, f32)]) -> Survey<'_> {
        Survey {
            rt,
            peaks,
            detection_limit: 10.0,
        }
    }

    #[test]
    fn window_record_interpolates_between_the_bordering_scans() {
        // z = 2 precursor at 500.0; one co-isolated peak at 500.2.
        let before = [
            (499.0, 5.0f32),
            (500.0, 100.0),
            (500.2, 50.0),
            (500.5017, 50.0),
        ];
        let after = [(500.0, 300.0f32), (500.2, 50.0), (500.5017, 150.0)];
        let rec = window_record(
            Some(survey(10.0, &before)),
            Some(survey(20.0, &after)),
            15.0,
            499.8,
            500.6,
            500.0,
            2,
        )
        .unwrap();
        // Halfway: precursor (100+50)/2 + (300+150)/2 = 300, interferent 50.
        assert!((rec.precursor - 300.0).abs() < 1e-9, "{rec:?}");
        assert!((rec.total - 350.0).abs() < 1e-9, "{rec:?}");
        assert!((rec.precursor_fraction() - 300.0 / 350.0).abs() < 1e-12);
        // Hidden: detection limit 10 × density × width 0.8.
        let dens_before = 4.0 / 10.0;
        let dens_after = 3.0 / 10.0;
        let expect = 0.5 * 10.0 * dens_before * 0.8 + 0.5 * 10.0 * dens_after * 0.8;
        assert!((rec.hidden - expect).abs() < 1e-9, "{rec:?}");
    }

    #[test]
    fn window_record_needs_a_bordering_scan_and_a_window() {
        assert_eq!(window_record(None, None, 1.0, 1.0, 2.0, 1.5, 2), None);
        let p = [(500.0, 1.0f32)];
        assert_eq!(
            window_record(Some(survey(0.0, &p)), None, 1.0, 2.0, 2.0, 1.5, 2),
            None
        );
        let one = window_record(Some(survey(0.0, &p)), None, 1.0, 499.0, 501.0, 500.0, 1).unwrap();
        assert_eq!(one.precursor_fraction(), 1.0);
    }

    /// Scans generated from the model itself: the fit recovers the
    /// coefficients and the interference of each scan.
    fn synthetic(n: usize) -> Vec<ScanInput> {
        let (a, b, c) = (1.0, 0.5, 2.0);
        (0..n)
            .map(|i| {
                let pic = 1e5 * (1.0 + (i % 7) as f64);
                let f = 0.3 + 0.7 * ((i * 37) % 100) as f64 / 100.0;
                let hidden = 1e4 * ((i * 13) % 5) as f64;
                let total = a * pic * f + b * pic * (1.0 - f) + c * hidden;
                let reporters = vec![(total / 4.0) as f32; 4];
                // The window record stores absolute currents; any scale works
                // since only the fraction enters the model.
                ScanInput {
                    charge: 2,
                    reporters,
                    pic,
                    window: Some(WindowRecord {
                        total: 1000.0,
                        precursor: 1000.0 * f,
                        hidden,
                    }),
                }
            })
            .collect()
    }

    #[test]
    fn fit_recovers_the_generating_model() {
        let scans = synthetic(400);
        let m = InterferenceModel::fit(&scans).unwrap();
        assert!((m.precursor[0] - 1.0).abs() < 1e-6, "{m:?}");
        assert!((m.co_isolated - 0.5).abs() < 1e-6, "{m:?}");
        assert!((m.hidden - 2.0).abs() < 1e-6, "{m:?}");
        assert!(m.r2_log > 0.999, "{m:?}");
        let s = &scans[5];
        let w = s.window.unwrap();
        let f = w.precursor_fraction();
        let own = s.pic * f;
        let total = own + 0.5 * s.pic * (1.0 - f) + 2.0 * w.hidden;
        let eil = m.interference(s).unwrap();
        assert!((eil - (1.0 - own / total).min(MAX_INTERFERENCE)).abs() < 1e-6);
    }

    #[test]
    fn too_few_scans_are_not_fitted() {
        let (model, out) = unmix_run(&synthetic(MIN_FIT_SCANS - 1));
        assert!(model.is_none());
        assert!(out.iter().all(Option::is_none));
    }

    #[test]
    fn interference_is_removed_with_the_run_profile() {
        // Background profile 1:1:1:1 (every scan is flat), except one scan
        // whose first channel holds the precursor's own signal.
        let mut scans = synthetic(400);
        scans[0].reporters = vec![400.0, 100.0, 100.0, 100.0];
        let (model, out) = unmix_run(&scans);
        let model = model.unwrap();
        let eil = model.interference(&scans[0]).unwrap();
        let profile = background_profile(&scans).unwrap();
        let u = out[0].as_ref().unwrap();
        assert!((u.interference as f64 - eil).abs() < 1e-6);
        for (k, (&v, &p)) in scans[0].reporters.iter().zip(&profile).enumerate() {
            let expect = ((v as f64) - eil * 700.0 * p).max(0.0);
            assert!((u.reporters[k] as f64 - expect).abs() < 1e-2, "{k}");
        }
        // Interference lowers the flat channels relative to the first.
        assert!(u.reporters[1] / u.reporters[0] < 0.25);
    }

    #[test]
    fn scans_without_a_window_record_are_left_out() {
        let mut scans = synthetic(400);
        scans[3].window = None;
        let (_, out) = unmix_run(&scans);
        assert!(out[3].is_none());
        assert!(out[4].is_some());
    }
}
