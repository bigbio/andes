//! Chimeric two-pass cascade: detect co-isolated precursors in an MS2 scan's MS1
//! isolation window (excluding the selected precursor), then run a targeted
//! second-peptide search at each. This is the speed-correct chimeric path: it
//! scores few candidates at MS1-confirmed masses instead of thousands across the
//! blind window (see docs/parity-analysis/notes/2026-05-30-chimeric-cost-profile.md).
//!
//! Task 1 (this commit) ships only the detector. The targeted second-peptide
//! search (`search_secondary`) and the binary-level driver land in Tasks 2/3,
//! which consume `CoIsolated` / `detect_coisolated`. Until then they are
//! unreferenced outside tests, so allow dead_code at the module level.
#![allow(dead_code)]

use crate::chimeric_features::precursor_isotope_match;
use model::mass::{ISOTOPE, PROTON};

/// A co-isolated precursor detected in the MS1 isolation window.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct CoIsolated {
    pub mono_mz: f64,
    pub charge: u8,
    pub neutral_mass: f64,
}

/// Detect co-isolated precursors in `ms1_peaks` (m/z-sorted) within the isolation
/// window `[win_lo, win_hi]`, EXCLUDING the envelope at `selected_mz` (the peptide
/// Pass 1 already searched). Tries charges in `charge_range`; accepts an envelope
/// whose averagine KL is below `max_kl`. Returns at most `max_n` (highest-intensity
/// monoisotopic peaks first).
#[allow(clippy::too_many_arguments)]
pub(crate) fn detect_coisolated(
    ms1_peaks: &[(f64, f32)],
    win_lo: f64,
    win_hi: f64,
    selected_mz: f64,
    charge_range: std::ops::RangeInclusive<u8>,
    tol_da: f64,
    max_kl: f32,
    max_n: usize,
) -> Vec<CoIsolated> {
    // Candidate monoisotopic peaks = peaks inside the window, sorted by intensity desc.
    let lo_idx = ms1_peaks.partition_point(|&(mz, _)| mz < win_lo);
    let mut cands: Vec<(f64, f32)> = ms1_peaks[lo_idx..]
        .iter()
        .take_while(|&&(mz, _)| mz <= win_hi)
        .copied()
        .collect();
    cands.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

    let mut out: Vec<CoIsolated> = Vec::new();
    for &(mz, _inten) in &cands {
        if (mz - selected_mz).abs() <= tol_da {
            continue; // skip the selected precursor (monoisotope)
        }
        // Skip the selected precursor's HIGHER isotope peaks too: a peak at
        // selected_mz + k*ISOTOPE/z (k >= 1) is part of the Pass-1 envelope, not
        // a distinct co-isolated species. The selected charge is unknown here, so
        // reject if the peak lines up with any isotope spacing in `charge_range`.
        if mz > selected_mz
            && charge_range.clone().filter(|&z| z != 0).any(|z| {
                let d = mz - selected_mz;
                (1..6).any(|k| (d - k as f64 * ISOTOPE / z as f64).abs() <= tol_da)
            })
        {
            continue;
        }
        // Don't re-report a peak that's an isotope of an already-accepted envelope.
        if out.iter().any(|c| {
            let d = (mz - c.mono_mz).abs();
            (0..6).any(|k| (d - k as f64 * ISOTOPE / c.charge as f64).abs() <= tol_da)
        }) {
            continue;
        }
        // Try charges; accept the lowest-KL charge under max_kl.
        let mut best: Option<(f32, CoIsolated)> = None;
        for z in charge_range.clone() {
            if z == 0 {
                continue;
            }
            let neutral = (mz - PROTON) * z as f64;
            let (kl, _snr) = precursor_isotope_match(ms1_peaks, mz, z, neutral, tol_da, 4);
            if kl <= max_kl && best.as_ref().is_none_or(|(bk, _)| kl < *bk) {
                best = Some((
                    kl,
                    CoIsolated {
                        mono_mz: mz,
                        charge: z,
                        neutral_mass: neutral,
                    },
                ));
            }
        }
        if let Some((_, c)) = best {
            out.push(c);
        }
        if out.len() >= max_n {
            break;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::isotope::averagine_isotope_envelope;

    /// Build a synthetic MS1 peak list (m/z-sorted) containing a 4-peak averagine
    /// envelope for `(mono_mz, charge, neutral_mass)` scaled by `scale`.
    fn envelope(mono_mz: f64, charge: u8, neutral: f64, scale: f32) -> Vec<(f64, f32)> {
        let env = averagine_isotope_envelope(neutral, 4);
        (0..4)
            .map(|k| {
                (
                    mono_mz + k as f64 * ISOTOPE / charge as f64,
                    (env[k] as f32) * scale,
                )
            })
            .collect()
    }

    #[test]
    fn detects_coisolated_excludes_selected() {
        let z = 2u8;
        let selected_mz = 600.0;
        let sel_neutral = (selected_mz - PROTON) * z as f64;
        let co_mz = 600.7; // a second precursor within a ~2 Da window
        let co_neutral = (co_mz - PROTON) * z as f64;
        let mut peaks = envelope(selected_mz, z, sel_neutral, 1000.0);
        peaks.extend(envelope(co_mz, z, co_neutral, 500.0));
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        let got = detect_coisolated(&peaks, 599.0, 601.5, selected_mz, 2..=3, 0.02, 0.5, 2);
        assert_eq!(got.len(), 1, "exactly one co-isolated (selected excluded)");
        assert!((got[0].mono_mz - co_mz).abs() < 0.02);
        assert_eq!(got[0].charge, z);
    }

    #[test]
    fn no_coisolation_when_only_selected_present() {
        let z = 2u8;
        let selected_mz = 600.0;
        let peaks = envelope(selected_mz, z, (selected_mz - PROTON) * z as f64, 1000.0);
        let got = detect_coisolated(&peaks, 599.0, 601.5, selected_mz, 2..=3, 0.02, 0.5, 2);
        assert!(got.is_empty(), "only the selected precursor -> no co-isolation");
    }
}
