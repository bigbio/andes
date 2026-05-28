//! Precursor isotope-envelope matching against the observed MS1.
//!
//! For a chimeric PSM, this scores how well the peptide's THEORETICAL
//! precursor isotope envelope (averagine + Poisson model) matches the
//! OBSERVED MS1 isotope envelope around the precursor m/z. A poor match
//! (high KL divergence) flags a likely spurious co-isolation; a low SNR
//! flags a precursor barely above MS1 noise. Both are emitted as additive
//! `PsmFeatures` columns so Percolator can learn weights without disturbing
//! existing feature distributions.
//!
//! NOTE (Task 3, commit 1): the helper is wired into the feature fill in
//! commit 2. Until then the module is `allow(dead_code)` so the
//! self-contained helper + tests can land first.
#![allow(dead_code)]

use model::isotope::averagine_isotope_envelope;
use model::mass::ISOTOPE;

/// Match a peptide's theoretical precursor isotope envelope against the
/// observed MS1 peaks. Returns `(kl_divergence, snr)`.
///
/// `ms1_peaks` MUST be sorted ascending by m/z (the [`crate`]'s
/// `Ms1Link::ms1_peaks` are emitted in mzML order, which is m/z-sorted).
/// For `k` in `0..n_isotopes`, the observed intensity is the
/// **maximum-intensity** MS1 peak whose m/z lies within `tol_da` of
/// `theo_mono_mz + k * ISOTOPE / charge` (0.0 if no peak falls in the
/// window). The observed vector is normalized to sum 1.0, then
/// `KL(observed || averagine(neutral_mass, n_isotopes))` is computed with an
/// epsilon floor on the theoretical probabilities.
///
/// SNR = monoisotope (`k = 0`) observed intensity / median of the nonzero
/// MS1 intensities (a robust noise proxy).
///
/// Degenerate cases:
/// * No peak found at ANY isotope position (observed sum is 0) → returns
///   `(10.0, 0.0)` — a large KL so Percolator can treat "no envelope at the
///   precursor m/z" as a strong negative signal.
/// * `ms1_peaks` empty / `n_isotopes == 0` → same `(10.0, 0.0)`.
pub(crate) fn precursor_isotope_match(
    ms1_peaks: &[(f64, f32)],
    theo_mono_mz: f64,
    charge: u8,
    neutral_mass: f64,
    tol_da: f64,
    n_isotopes: usize,
) -> (f32, f32) {
    // KL value returned when no observed envelope can be formed. Large so
    // Percolator reads "no envelope at precursor m/z" as strongly negative.
    const NO_ENVELOPE_KL: f32 = 10.0;
    // Floor on the theoretical (and observed, post-normalization) probability
    // inside the log, to keep KL finite.
    const EPS: f64 = 1e-9;

    if ms1_peaks.is_empty() || n_isotopes == 0 || charge == 0 {
        return (NO_ENVELOPE_KL, 0.0);
    }

    let charge_f = charge as f64;
    let spacing = ISOTOPE / charge_f;

    // Observed intensity at each isotope position: the max-intensity peak
    // within tol_da of the target m/z (binary search since ms1_peaks is
    // m/z-sorted), 0.0 if none.
    let mut observed = vec![0.0f64; n_isotopes];
    for (k, obs_k) in observed.iter_mut().enumerate() {
        let target = theo_mono_mz + (k as f64) * spacing;
        let lo = target - tol_da;
        let hi = target + tol_da;
        // First peak with m/z >= lo.
        let start = ms1_peaks.partition_point(|&(mz, _)| mz < lo);
        let mut best = 0.0f64;
        for &(mz, inten) in &ms1_peaks[start..] {
            if mz > hi {
                break;
            }
            if (inten as f64) > best {
                best = inten as f64;
            }
        }
        *obs_k = best;
    }

    let obs_sum: f64 = observed.iter().sum();
    if obs_sum <= 0.0 {
        // No peak at any isotope position — strong negative signal.
        return (NO_ENVELOPE_KL, 0.0);
    }

    // Normalize observed to a probability distribution summing to 1.0.
    for o in observed.iter_mut() {
        *o /= obs_sum;
    }

    // Theoretical averagine envelope (already normalized to sum 1.0).
    let theo = averagine_isotope_envelope(neutral_mass, n_isotopes);

    // KL(observed || theo) = Σ_k o_k * ln(o_k / max(t_k, eps)), over o_k > 0.
    let mut kl = 0.0f64;
    for (k, &o) in observed.iter().enumerate() {
        if o > 0.0 {
            let t = theo.get(k).copied().unwrap_or(0.0).max(EPS);
            kl += o * (o / t).ln();
        }
    }
    // KL is >= 0 in theory; clamp tiny negative rounding noise.
    let kl = kl.max(0.0) as f32;

    // SNR = monoisotope observed intensity / median nonzero MS1 intensity.
    // We recover the un-normalized monoisotope intensity as observed[0] *
    // obs_sum (observed[0] is the normalized share of the monoisotope).
    let mono_intensity = observed[0] * obs_sum;
    let snr = match median_nonzero_intensity(ms1_peaks) {
        Some(med) if med > 0.0 => (mono_intensity / med) as f32,
        _ => 0.0,
    };

    (kl, snr)
}

/// Median of the nonzero MS1 peak intensities (a robust noise proxy).
/// Returns `None` when there are no nonzero peaks.
fn median_nonzero_intensity(ms1_peaks: &[(f64, f32)]) -> Option<f64> {
    let mut nz: Vec<f64> = ms1_peaks
        .iter()
        .map(|&(_, i)| i as f64)
        .filter(|&i| i > 0.0)
        .collect();
    if nz.is_empty() {
        return None;
    }
    nz.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = nz.len();
    Some(if n % 2 == 1 {
        nz[n / 2]
    } else {
        0.5 * (nz[n / 2 - 1] + nz[n / 2])
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::isotope::averagine_isotope_envelope;
    use model::mass::ISOTOPE;

    /// Build an MS1 peak list whose envelope at `mono_mz` is shaped exactly
    /// like the averagine envelope for `neutral_mass`, scaled by `scale`,
    /// sitting on a noise floor of `noise` peaks placed away from the
    /// precursor.
    fn build_clean_envelope(
        mono_mz: f64,
        charge: u8,
        neutral_mass: f64,
        n: usize,
        scale: f32,
        noise: f32,
    ) -> Vec<(f64, f32)> {
        let env = averagine_isotope_envelope(neutral_mass, n);
        let spacing = ISOTOPE / charge as f64;
        let mut peaks: Vec<(f64, f32)> = Vec::new();
        // Noise peaks well below the precursor m/z (still m/z-sorted first).
        peaks.push((mono_mz - 50.0, noise));
        peaks.push((mono_mz - 40.0, noise));
        for (k, &e) in env.iter().enumerate() {
            let mz = mono_mz + k as f64 * spacing;
            peaks.push((mz, scale * e as f32));
        }
        // Noise peaks above the envelope.
        peaks.push((mono_mz + 50.0, noise));
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        peaks
    }

    #[test]
    fn precursor_isotope_match_clean_envelope_low_kl_high_snr() {
        let charge = 2u8;
        let neutral_mass = 1500.0;
        let mono_mz = (neutral_mass + charge as f64 * model::mass::PROTON) / charge as f64;
        let n = 4;
        // Strong precursor (scale 1000) over a low noise floor (1.0).
        let ms1 = build_clean_envelope(mono_mz, charge, neutral_mass, n, 1000.0, 1.0);

        let (kl, snr) =
            precursor_isotope_match(&ms1, mono_mz, charge, neutral_mass, 0.01, n);

        assert!(kl < 0.1, "clean averagine-shaped envelope should give small KL, got {kl}");
        assert!(snr > 1.0, "strong precursor over low noise should give SNR > 1, got {snr}");
    }

    #[test]
    fn precursor_isotope_match_no_envelope_large_kl_zero_snr() {
        let charge = 2u8;
        let neutral_mass = 1500.0;
        let mono_mz = (neutral_mass + charge as f64 * model::mass::PROTON) / charge as f64;
        let n = 4;
        // Peaks only at unrelated m/z values — nothing within tol of the
        // precursor isotope ladder.
        let ms1: Vec<(f64, f32)> = vec![
            (mono_mz - 30.0, 500.0),
            (mono_mz - 20.0, 800.0),
            (mono_mz + 25.0, 600.0),
            (mono_mz + 40.0, 700.0),
        ];

        let (kl, snr) =
            precursor_isotope_match(&ms1, mono_mz, charge, neutral_mass, 0.01, n);

        assert!(kl >= 1.0, "no envelope at precursor m/z should give large KL, got {kl}");
        assert!(snr.abs() < 1e-6, "no monoisotope peak should give SNR ~0, got {snr}");
    }

    #[test]
    fn precursor_isotope_match_empty_peaks_returns_large_kl_zero_snr() {
        let charge = 2u8;
        let neutral_mass = 1500.0;
        let mono_mz = (neutral_mass + charge as f64 * model::mass::PROTON) / charge as f64;
        let (kl, snr) = precursor_isotope_match(&[], mono_mz, charge, neutral_mass, 0.01, 4);
        assert_eq!(kl, 10.0, "empty MS1 should return the no-envelope KL sentinel");
        assert_eq!(snr, 0.0, "empty MS1 should return SNR 0.0");
    }

    #[test]
    fn precursor_isotope_match_picks_max_within_tolerance() {
        // Two peaks within tol of the monoisotope position; the max-intensity
        // one should drive the SNR.
        let charge = 1u8;
        let neutral_mass = 1000.0;
        let mono_mz = neutral_mass + model::mass::PROTON;
        let n = 3;
        let env = averagine_isotope_envelope(neutral_mass, n);
        let spacing = ISOTOPE / charge as f64;
        let mut ms1: Vec<(f64, f32)> = vec![
            (mono_mz - 0.005, 200.0), // within tol, smaller
            (mono_mz + 0.004, 900.0), // within tol, larger -> should win
            (mono_mz + spacing, 1000.0 * env[1] as f32),
            (mono_mz + 2.0 * spacing, 1000.0 * env[2] as f32),
            (mono_mz + 30.0, 1.0), // noise floor
        ];
        ms1.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let (_kl, snr) =
            precursor_isotope_match(&ms1, mono_mz, charge, neutral_mass, 0.01, n);
        // mono intensity 900 / median nonzero. Median of {200,900,~270,~36,1}
        // is the middle value -> well below 900, so SNR comfortably > 1.
        assert!(snr > 1.0, "max-within-tol monoisotope peak should give SNR > 1, got {snr}");
    }
}
