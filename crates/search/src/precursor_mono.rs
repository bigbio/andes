//! MS1 isotope-envelope precursor correction (`--precursor-mono`).
//!
//! WHY. On pGlyco2 mouse liver (`MouseLiver-Z-T-1`), 84 of the 515 reference
//! spectra andes misses have the true backbone in the candidate set, but the
//! precursor m/z the instrument recorded sits 3–6 Da (70 of them exactly
//! +4.01 Da) ABOVE the glycopeptide's monoisotopic mass: the firmware's
//! monoisotopic-peak pick failed on a wide, high-mass envelope and the scan
//! was recorded on an M+3..M+6 isotopologue. Widening `--isotope-error` to
//! reach them is not a fix — every extra offset `k` is mass-degenerate with a
//! glycan-composition change (`(k, X)` ≡ `(k−1, X + Hex + Fuc − NeuGc)` within
//! 20 ppm), so the search fills the largest allowed offset with hundreds of
//! PSMs whose composition is arbitrary (bigbio/andes#64, arms C and D).
//!
//! WHAT. For each MS2, read the preceding MS1, and score the observed isotope
//! envelope at the reported charge against the theoretical glycopeptide
//! envelope under every hypothesis "the recorded precursor is the M+k peak",
//! `k = 0..=max_shift`. If a `k > 0` hypothesis fits clearly better than the
//! recorded one, move the precursor m/z down by `k` isotopes and search it
//! with the DEFAULT narrow window. Nothing is widened.
//!
//! HOW THE FIT WORKS. Under hypothesis `k` the monoisotope sits at
//! `mz − k·ISO/z`; the observed vector is the max MS1 intensity within
//! tolerance of `mono + i·ISO/z` for `i = −1..n`, and the theoretical vector
//! is `[0, t_0, …, t_{n−1}]` from [`model::isotope::glycopeptide_isotope_envelope`].
//! The fit is their cosine. Two properties make it discriminative between
//! neighbouring `k` (which a KL on the observed distribution is not):
//!
//! * the `i = −1` slot, whose theoretical intensity is 0, penalises a
//!   hypothesis that leaves a real peak one isotope BELOW its monoisotope —
//!   exactly what the recorded (`k = 0`) hypothesis does on a mis-picked scan;
//! * a missing peak where theory expects one lowers the cosine, so a
//!   hypothesis whose monoisotope is not observed at all cannot win.
//!
//! WHEN IT APPLIES. A shift is applied only when the best `k > 0` clears an
//! absolute fit floor, beats the recorded hypothesis by a margin, and its
//! monoisotope is above the MS1 noise floor ([`MonoParams`]). Otherwise the
//! precursor is left as recorded, and with no MS1 (MGF, MS2-only mzML) or the
//! flag off nothing in the run changes.
//!
//! BACK-OFF BY ONE. Neighbouring hypotheses differ by one smooth step of a
//! wide envelope, so `k` versus `k−1` is the one confusion the fit cannot
//! always settle (synthetic z=3, 3.5 kDa: 1.00 vs 0.87). The search still
//! sweeps `--isotope-error 0..2` from the corrected mass, so an UNDERSHOOT by
//! one costs nothing while an OVERSHOOT by one puts the true mass out of
//! reach. The applied shift is therefore `best_k − backoff` (default 1): a
//! scan recorded on M+4 is searched from M+3 and the sweep's `+1` lands it,
//! and a scan the fit reads as M+1 is left alone — the sweep already covers
//! it. `MonoShift` in the PIN is the applied shift; the effective offset of a
//! PSM is `MonoShift + isotope_error`.

use model::isotope::glycopeptide_isotope_envelope;
use model::mass::{ISOTOPE, PROTON};

/// Tunables for [`correct_precursor`]. Defaults are what `--precursor-mono
/// auto` ships; every field has a hidden CLI override.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MonoParams {
    /// Largest number of isotopes the recorded precursor may be moved DOWN by.
    /// The pGlyco2 liver misses reach +6 (bigbio/andes#64 test set).
    pub max_shift: u8,
    /// MS1 peak-matching tolerance, ppm of the target m/z.
    pub tol_ppm: f64,
    /// Floor on the matching tolerance in Da (low m/z).
    pub tol_min_da: f64,
    /// Minimum cosine fit of the winning shifted hypothesis.
    pub min_fit: f32,
    /// Minimum cosine improvement of the winning shifted hypothesis over the
    /// recorded (`k = 0`) hypothesis.
    pub min_gain: f32,
    /// Minimum monoisotope intensity of the winning hypothesis, in units of
    /// the MS1 median non-zero intensity.
    pub min_snr: f32,
    /// Longest theoretical envelope considered.
    pub max_isotopes: usize,
    /// Isotopes held back from the best-fitting shift when applying it (see
    /// the module doc). 0 applies the best fit verbatim.
    pub backoff: u8,
}

impl Default for MonoParams {
    fn default() -> Self {
        Self {
            max_shift: 6,
            tol_ppm: 10.0,
            tol_min_da: 0.005,
            min_fit: 0.90,
            min_gain: 0.15,
            min_snr: 3.0,
            max_isotopes: 10,
            backoff: 1,
        }
    }
}

/// Outcome of [`correct_precursor`] for one MS2. `shift == 0` means the
/// precursor was left as recorded (the diagnostic fields still describe the
/// envelope at the recorded m/z).
#[derive(Debug, Clone, PartialEq)]
pub struct MonoCorrection {
    /// Precursor m/z as recorded by the instrument.
    pub recorded_mz: f64,
    /// Precursor m/z to search: `recorded_mz − shift · ISO / charge`.
    pub corrected_mz: f64,
    pub charge: u8,
    /// Isotopes subtracted from the recorded precursor (0 = unchanged). This
    /// is `best_shift − backoff` when the correction was accepted.
    pub shift: i8,
    /// Best-fitting `k` in `0..=max_shift` whether or not it was applied.
    pub best_shift: i8,
    /// Cosine fit of the recorded hypothesis (`k = 0`).
    pub fit_recorded: f32,
    /// Cosine fit of `best_shift`.
    pub fit_best: f32,
    /// Cosine fit of the "monoisotope one isotope ABOVE the recorded peak"
    /// hypothesis (`k = −1`). Diagnostic only, never applied: the firmware
    /// mis-picks only ever too high.
    pub fit_up: f32,
    /// Monoisotope intensity of the SEARCHED hypothesis over the MS1 median
    /// non-zero intensity.
    pub snr: f32,
    /// Number of theoretical isotopes fitted.
    pub n_iso: u8,
    /// Cosine fit per `k` in `0..=max_shift` (diagnostic dump).
    pub fits: Vec<f32>,
}

impl MonoCorrection {
    pub fn applied(&self) -> bool {
        self.shift != 0
    }
    /// Fit of the hypothesis the search actually uses.
    pub fn fit_searched(&self) -> f32 {
        if self.applied() {
            self.fit_best
        } else {
            self.fit_recorded
        }
    }
}

/// Number of isotopes to fit for a theoretical envelope: enough to cover
/// 98.5% of the distribution, at least 4, at most `max`.
fn envelope_length(theo: &[f64], max: usize) -> usize {
    let mut cum = 0.0;
    for (i, &t) in theo.iter().enumerate() {
        cum += t;
        if cum >= 0.985 {
            return (i + 1).clamp(4, max);
        }
    }
    theo.len().clamp(1, max)
}

/// Max MS1 intensity within `tol` of `target` (0.0 if none). `ms1` is
/// m/z-ascending.
fn peak_at(ms1: &[(f64, f32)], target: f64, tol: f64) -> f32 {
    let lo = target - tol;
    let hi = target + tol;
    let start = ms1.partition_point(|&(mz, _)| mz < lo);
    let mut best = 0.0f32;
    for &(mz, inten) in &ms1[start..] {
        if mz > hi {
            break;
        }
        if inten > best {
            best = inten;
        }
    }
    best
}

/// Cosine fit of the observed MS1 envelope with its monoisotope at `mono_mz`
/// against the theoretical glycopeptide envelope. Returns
/// `(fit, mono_intensity, n_iso)`; `fit` is 0.0 when nothing is observed at
/// any position.
pub fn envelope_fit(
    ms1: &[(f64, f32)],
    mono_mz: f64,
    charge: u8,
    params: &MonoParams,
) -> (f32, f32, usize) {
    let z = charge as f64;
    let neutral = (mono_mz - PROTON) * z;
    if neutral <= 0.0 {
        return (0.0, 0.0, 0);
    }
    let theo_full = glycopeptide_isotope_envelope(neutral, params.max_isotopes.max(1));
    let n = envelope_length(&theo_full, params.max_isotopes.max(1));
    let spacing = ISOTOPE / z;

    // Slot 0 is the i = −1 position (theory 0), slots 1..=n are M+0..M+(n−1).
    let mut dot = 0.0f64;
    let mut o2 = 0.0f64;
    let mut t2 = 0.0f64;
    let mut mono_intensity = 0.0f32;
    for i in -1..(n as i32) {
        let target = mono_mz + (i as f64) * spacing;
        let tol = (target * params.tol_ppm * 1e-6).max(params.tol_min_da);
        let obs = peak_at(ms1, target, tol) as f64;
        let theo = if i < 0 { 0.0 } else { theo_full[i as usize] };
        if i == 0 {
            mono_intensity = obs as f32;
        }
        dot += obs * theo;
        o2 += obs * obs;
        t2 += theo * theo;
    }
    if o2 <= 0.0 || t2 <= 0.0 {
        return (0.0, mono_intensity, n);
    }
    let fit = (dot / (o2.sqrt() * t2.sqrt())).clamp(0.0, 1.0) as f32;
    (fit, mono_intensity, n)
}

/// Median of the non-zero MS1 intensities (robust noise proxy); `None` when
/// there are none.
pub fn median_nonzero_intensity(ms1: &[(f64, f32)]) -> Option<f32> {
    let mut nz: Vec<f32> = ms1.iter().map(|&(_, i)| i).filter(|&i| i > 0.0).collect();
    if nz.is_empty() {
        return None;
    }
    let mid = nz.len() / 2;
    let (_, m, _) = nz.select_nth_unstable_by(mid, |a, b| a.total_cmp(b));
    let m = *m;
    if nz.len() % 2 == 1 {
        Some(m)
    } else {
        let (_, lo, _) = nz.select_nth_unstable_by(mid - 1, |a, b| a.total_cmp(b));
        Some(0.5 * (*lo + m))
    }
}

/// Decide whether the recorded precursor `mz` (charge `charge`) of an MS2 is
/// a higher isotopologue of its glycopeptide, from the preceding MS1 `ms1`
/// (m/z-ascending peaks). `median_noise` is
/// [`median_nonzero_intensity`] of that MS1 (computed once per MS1 by the
/// caller; `None` disables the SNR gate).
///
/// Returns `None` when there is nothing to fit (empty MS1, unknown charge,
/// non-positive m/z); otherwise a [`MonoCorrection`] whose `shift` is 0 when
/// the precursor is kept as recorded.
pub fn correct_precursor(
    ms1: &[(f64, f32)],
    mz: f64,
    charge: u8,
    median_noise: Option<f32>,
    params: &MonoParams,
) -> Option<MonoCorrection> {
    if ms1.is_empty() || charge == 0 || mz.is_nan() || mz <= 0.0 {
        return None;
    }
    let spacing = ISOTOPE / charge as f64;
    let max_shift = params.max_shift as i32;

    let mut fits: Vec<f32> = Vec::with_capacity(max_shift as usize + 1);
    let mut monos: Vec<f32> = Vec::with_capacity(max_shift as usize + 1);
    let mut n_iso = 0usize;
    for k in 0..=max_shift {
        let mono_mz = mz - (k as f64) * spacing;
        let (fit, mono, n) = envelope_fit(ms1, mono_mz, charge, params);
        fits.push(fit);
        monos.push(mono);
        if k == 0 {
            n_iso = n;
        }
    }
    let (fit_up, _, _) = envelope_fit(ms1, mz + spacing, charge, params);

    // argmax over k, the SMALLEST k taking a tie: the correction should move
    // a precursor no further than the envelope demands.
    let mut best_shift = 0usize;
    for (k, &f) in fits.iter().enumerate() {
        if f > fits[best_shift] {
            best_shift = k;
        }
    }
    let snr_of = |k: usize| -> f32 {
        match median_noise {
            Some(med) if med > 0.0 => monos[k] / med,
            _ => f32::INFINITY,
        }
    };
    let apply = best_shift > 0
        && fits[best_shift] >= params.min_fit
        && fits[best_shift] - fits[0] >= params.min_gain
        && snr_of(best_shift) >= params.min_snr;
    let shift = if apply {
        best_shift.saturating_sub(params.backoff as usize)
    } else {
        0
    };
    let snr = {
        let v = snr_of(shift);
        if v.is_finite() {
            v
        } else {
            0.0
        }
    };
    Some(MonoCorrection {
        recorded_mz: mz,
        corrected_mz: mz - (shift as f64) * spacing,
        charge,
        shift: shift as i8,
        best_shift: best_shift as i8,
        fit_recorded: fits[0],
        fit_best: fits[best_shift],
        fit_up,
        snr,
        n_iso: n_iso as u8,
        fits,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// An MS1 carrying one glycopeptide envelope with its monoisotope at
    /// `mono_mz`, on a floor of `noise` peaks away from it.
    fn synth_ms1(mono_mz: f64, charge: u8, n: usize, scale: f32, noise: f32) -> Vec<(f64, f32)> {
        let neutral = (mono_mz - PROTON) * charge as f64;
        let env = glycopeptide_isotope_envelope(neutral, n);
        let spacing = ISOTOPE / charge as f64;
        let mut peaks: Vec<(f64, f32)> = (0..40)
            .map(|i| (mono_mz - 60.0 + i as f64 * 1.37, noise))
            .collect();
        for (k, &e) in env.iter().enumerate() {
            peaks.push((mono_mz + k as f64 * spacing, scale * e as f32));
        }
        peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
        peaks
    }

    #[test]
    fn recorded_monoisotope_is_left_alone() {
        let mono = 1500.25;
        let ms1 = synth_ms1(mono, 2, 8, 1e6, 100.0);
        let med = median_nonzero_intensity(&ms1);
        let c = correct_precursor(&ms1, mono, 2, med, &MonoParams::default()).unwrap();
        assert_eq!(c.shift, 0, "{c:?}");
        assert_eq!(c.best_shift, 0, "{c:?}");
        assert!(c.fit_recorded > 0.99, "{c:?}");
        assert_eq!(c.corrected_mz, mono);
    }

    #[test]
    fn precursor_recorded_on_m_plus_4_is_moved_down_four() {
        // A 3 kDa glycopeptide at z=2 whose scan was recorded on M+4.
        let mono = 1500.25;
        let ms1 = synth_ms1(mono, 2, 8, 1e6, 100.0);
        let med = median_nonzero_intensity(&ms1);
        let recorded = mono + 4.0 * ISOTOPE / 2.0;
        let c = correct_precursor(&ms1, recorded, 2, med, &MonoParams::default()).unwrap();
        assert_eq!(c.best_shift, 4, "{c:?}");
        // Default back-off: searched from M+1, the 0..2 sweep takes the last step.
        assert_eq!(c.shift, 3, "{c:?}");
        assert!((c.corrected_mz - (mono + ISOTOPE / 2.0)).abs() < 1e-9);
        assert!(c.fit_best > 0.99 && c.fit_recorded < 0.8, "{c:?}");
        assert!(c.applied());
        assert!(c.snr > 100.0, "{c:?}");
        let exact = MonoParams {
            backoff: 0,
            ..MonoParams::default()
        };
        let c0 = correct_precursor(&ms1, recorded, 2, med, &exact).unwrap();
        assert_eq!(c0.shift, 4, "{c0:?}");
        assert!((c0.corrected_mz - mono).abs() < 1e-9);
    }

    #[test]
    fn every_shift_from_two_to_six_is_recovered_at_z3() {
        let mono = 1180.7;
        let ms1 = synth_ms1(mono, 3, 10, 5e5, 50.0);
        let med = median_nonzero_intensity(&ms1);
        for k in 2..=6i32 {
            let recorded = mono + k as f64 * ISOTOPE / 3.0;
            let c = correct_precursor(&ms1, recorded, 3, med, &MonoParams::default()).unwrap();
            assert_eq!(c.best_shift as i32, k, "k={k}: {c:?}");
            assert_eq!(c.shift as i32, k - 1, "k={k}: {c:?}");
        }
        // k = 1 is the by-one ambiguity the sweep owns: the fit still names it,
        // but nothing is applied (best 1 − backoff 1 = 0).
        let c =
            correct_precursor(&ms1, mono + ISOTOPE / 3.0, 3, med, &MonoParams::default()).unwrap();
        assert_eq!(c.best_shift, 1, "{c:?}");
        assert_eq!(c.shift, 0, "{c:?}");
        assert!(!c.applied());
    }

    #[test]
    fn no_envelope_near_precursor_keeps_recorded() {
        let ms1: Vec<(f64, f32)> = (0..50).map(|i| (400.0 + i as f64 * 7.3, 1000.0)).collect();
        let med = median_nonzero_intensity(&ms1);
        let c = correct_precursor(&ms1, 1500.0, 2, med, &MonoParams::default()).unwrap();
        assert_eq!(c.shift, 0);
        assert_eq!(c.fit_recorded, 0.0);
    }

    #[test]
    fn weak_monoisotope_below_snr_floor_is_not_applied() {
        let mono = 1500.25;
        // Envelope barely above a heavy noise floor: fits may look fine, but
        // the mono peak is ~1x the median, under the default SNR gate.
        let ms1 = synth_ms1(mono, 2, 8, 800.0, 200.0);
        let med = median_nonzero_intensity(&ms1);
        let recorded = mono + 4.0 * ISOTOPE / 2.0;
        let c = correct_precursor(&ms1, recorded, 2, med, &MonoParams::default()).unwrap();
        assert_eq!(c.best_shift, 4, "{c:?}");
        assert_eq!(c.shift, 0, "SNR gate must hold: {c:?}");
        let loose = MonoParams {
            min_snr: 0.0,
            ..MonoParams::default()
        };
        let c2 = correct_precursor(&ms1, recorded, 2, med, &loose).unwrap();
        assert_eq!(c2.shift, 3, "{c2:?}");
    }

    #[test]
    fn degenerate_inputs_return_none() {
        let ms1 = synth_ms1(1500.25, 2, 8, 1e6, 100.0);
        assert!(correct_precursor(&[], 1500.0, 2, None, &MonoParams::default()).is_none());
        assert!(correct_precursor(&ms1, 1500.0, 0, None, &MonoParams::default()).is_none());
        assert!(correct_precursor(&ms1, 0.0, 2, None, &MonoParams::default()).is_none());
    }

    #[test]
    fn median_nonzero_ignores_zero_intensities() {
        let ms1 = vec![(1.0, 0.0), (2.0, 5.0), (3.0, 1.0), (4.0, 3.0)];
        assert_eq!(median_nonzero_intensity(&ms1), Some(3.0));
        let even = vec![(1.0, 4.0), (2.0, 1.0), (3.0, 3.0), (4.0, 2.0)];
        assert_eq!(median_nonzero_intensity(&even), Some(2.5));
        assert_eq!(median_nonzero_intensity(&[(1.0, 0.0)]), None);
    }
}
