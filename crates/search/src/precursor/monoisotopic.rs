use model::isotope::glycopeptide_isotope_envelope;
use model::mass::{ISOTOPE, PROTON};
use model::stats::{median_f32_mut, Even};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MonoParams {
    pub max_shift: u8,
    pub tol_ppm: f64,
    pub tol_min_da: f64,
    pub min_fit: f32,
    pub min_gain: f32,
    pub min_snr: f32,
    pub max_isotopes: usize,
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

#[derive(Debug, Clone, PartialEq)]
pub struct MonoCorrection {
    pub recorded_mz: f64,
    pub corrected_mz: f64,
    pub charge: u8,
    pub shift: i8,
    pub best_shift: i8,
    pub fit_recorded: f32,
    pub fit_best: f32,
    pub fit_up: f32,
    pub snr: f32,
    pub n_iso: u8,
    pub fits: Vec<f32>,
}

impl MonoCorrection {
    pub fn applied(&self) -> bool {
        self.shift != 0
    }
    pub fn fit_searched(&self) -> f32 {
        if self.applied() {
            self.fit_best
        } else {
            self.fit_recorded
        }
    }

    pub fn verified_within(&self, min_fit: f32, max_offset: i8) -> bool {
        self.fit_best >= min_fit && self.best_shift - self.shift <= max_offset
    }
}

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

/// The m/z of the most intense MS1 peak within `ppm` of the reported precursor `mz`, when it
/// is at least `ratio` times as intense as the peak at `mz` itself (searched within
/// `ppm / 8`). The instrument at times reports a side peak of a Fourier-transform signal, a
/// fraction as intense and 15-25 ppm off the ion it belongs to. `None` keeps the reported m/z.
pub fn snap_precursor_mz(ms1: &[(f64, f32)], mz: f64, ppm: f64, ratio: f32) -> Option<f64> {
    if ms1.is_empty() || mz.is_nan() || mz <= 0.0 {
        return None;
    }
    let tol = mz * ppm * 1e-6;
    let own = peak_at(ms1, mz, tol / 8.0);
    let start = ms1.partition_point(|&(m, _)| m < mz - tol);
    let mut best: Option<(f64, f32)> = None;
    for &(m, inten) in &ms1[start..] {
        if m > mz + tol {
            break;
        }
        if best.is_none_or(|(_, b)| inten > b) {
            best = Some((m, inten));
        }
    }
    let (m, inten) = best?;
    ((m - mz).abs() > tol / 8.0 && inten >= ratio * own.max(f32::MIN_POSITIVE)).then_some(m)
}

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

pub fn median_nonzero_intensity(ms1: &[(f64, f32)]) -> Option<f32> {
    let mut nz: Vec<f32> = ms1.iter().map(|&(_, i)| i).filter(|&i| i > 0.0).collect();
    median_f32_mut(&mut nz, Even::Mean)
}

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

    #[test]
    fn a_weak_side_peak_moves_to_the_intense_ion_beside_it() {
        // The reported 619.8634 is a side peak of the ion at 619.8765 (+21 ppm, 9x as intense).
        let ms1 = vec![
            (619.4734, 4e4),
            (619.8634, 3.9e5),
            (619.8765, 3.6e6),
            (620.0771, 5.0e6),
        ];
        let got = snap_precursor_mz(&ms1, 619.8634, 25.0, 3.0).unwrap();
        assert!((got - 619.8765).abs() < 1e-9);
        // A neighbour under three times as intense, or beyond the window, leaves it.
        assert_eq!(snap_precursor_mz(&ms1, 619.8634, 25.0, 10.0), None);
        assert_eq!(snap_precursor_mz(&ms1, 619.8634, 15.0, 3.0), None);
        // The reported peak is itself the most intense one.
        assert_eq!(snap_precursor_mz(&ms1, 619.8765, 25.0, 3.0), None);
        // Nothing at the reported m/z: the most intense peak in the window.
        let got = snap_precursor_mz(&ms1, 619.8700, 25.0, 3.0).unwrap();
        assert!((got - 619.8765).abs() < 1e-9);
        assert_eq!(snap_precursor_mz(&[], 619.8634, 25.0, 3.0), None);
    }

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
    fn verified_within_requires_a_fitted_envelope_the_window_reaches() {
        let p = MonoParams::default();
        let mono = 1500.25;
        let ms1 = synth_ms1(mono, 2, 8, 1e6, 100.0);
        let c = correct_precursor(&ms1, mono, 2, median_nonzero_intensity(&ms1), &p).unwrap();
        assert!(c.verified_within(p.min_fit, 1));
        let m4 = mono + 4.0 * ISOTOPE / 2.0;
        let c = correct_precursor(&ms1, m4, 2, median_nonzero_intensity(&ms1), &p).unwrap();
        assert!(c.applied());
        assert!(c.verified_within(p.min_fit, 1));
        assert!(!c.verified_within(p.min_fit, 0));
        let empty = synth_ms1(mono + 300.0, 2, 8, 1e6, 100.0);
        let c = correct_precursor(&empty, mono, 2, median_nonzero_intensity(&empty), &p).unwrap();
        assert!(!c.verified_within(p.min_fit, 1));
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
        let mono = 1500.25;
        let ms1 = synth_ms1(mono, 2, 8, 1e6, 100.0);
        let med = median_nonzero_intensity(&ms1);
        let recorded = mono + 4.0 * ISOTOPE / 2.0;
        let c = correct_precursor(&ms1, recorded, 2, med, &MonoParams::default()).unwrap();
        assert_eq!(c.best_shift, 4, "{c:?}");
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
