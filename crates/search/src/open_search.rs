//! Localisation of an open-search mass offset to a residue (`--open-search`).
//!
//! An open-search match says *how much* mass is unaccounted for, not *where* it
//! sits. This module assigns the offset to one residue: for every candidate
//! position it shifts that residue's mass by the offset, generates the b/y ions
//! of the resulting form, and sums the intensity of the spectrum's matched
//! peaks. The best position wins; the margin over the runner-up is the
//! confidence.
//!
//! The machinery mirrors the site-determining-ion idea of
//! [`crate::site_evidence`] and [`scoring_crate::mod_site_features`], but it
//! works from the bare peptide + an external delta instead of a peptide that
//! already carries a declared modification.
//!
//! The offset inferred from the precursor carries that precursor's mass error
//! (tens of mDa on a long peptide), which is large against a ppm fragment
//! window. [`refine_offset`] re-centres it on the shifted ions of one residue,
//! and [`localize_and_refine`] uses that per residue: a coarse window to find
//! the shifted ions, then the tight fragment window to separate residues. The
//! open-search driver re-scores the shifted form at the winning residue
//! (localization-aware scoring, Yu et al., *Nat. Commun.* 2020).
//!
//! Cost is O(sites x ions x log peaks) per candidate (binary search into the
//! sorted peak list), so the driver bounds how many re-scores it runs per
//! spectrum.

use model::mass::{H2O, PROTON};
use model::peptide::Peptide;

/// One localisation result.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Localization {
    /// 1-based residue carrying the offset; 0 when no position separates from
    /// the runner-up (or no ions match at all).
    pub site: u16,
    /// Share of the spectrum's total peak intensity explained by the winning
    /// site's matched b/y ions (0.0..=1.0).
    pub score: f32,
    /// `(best - second) / best` over the candidate positions (0.0..=1.0).
    /// Larger means the spectrum pins the offset to one residue.
    pub margin: f32,
}

/// Below this margin the best position is not called and `site` stays 0.
/// 2% of the winning site's intensity is about the noise between two positions
/// that differ by one weak fragment.
pub const MIN_LOCALIZATION_MARGIN: f32 = 0.02;

/// Largest offset correction (Da) `refine_offset` will apply. The correction
/// exists to remove the precursor's own mass error, which is tens of mDa; a
/// larger one is a wrong-site artefact, not calibration.
pub const MAX_REFINE_DA: f64 = 0.1;

/// The peptide with `delta` added to residue `site0` (0-based), rebuilt so the
/// scorer predicts the **shifted** fragment ions (Yu et al., *Nat. Commun.*
/// 2020, localization-aware open search).
///
/// The offset is carried in the residue mass, not in a `Modification` handle:
/// it is a runtime value with no declared record, and this form is used for
/// scoring only — it is never printed.
pub fn shifted_peptide(peptide: &Peptide, site0: usize, delta: f64) -> Peptide {
    let mut residues = peptide.residues.clone();
    if let Some(aa) = residues.get_mut(site0) {
        aa.mass += delta;
    }
    Peptide::new(residues, peptide.pre, peptide.post)
}

/// Nearest peak to `mz` within `tol_da`, from a list sorted by m/z.
fn nearest_peak(peaks: &[(f64, f32)], mz: f64, tol_da: f64) -> Option<usize> {
    let i = peaks.partition_point(|p| p.0 < mz);
    [i.saturating_sub(1), i]
        .into_iter()
        .find(|&j| j < peaks.len() && (peaks[j].0 - mz).abs() <= tol_da)
}

/// Re-centre `delta` on the fragments that carry it: the median signed
/// `theoretical - observed` m/z over the shifted b/y ions of `site0`, converted
/// to neutral mass.
///
/// The offset inferred from the precursor carries that precursor's measurement
/// error (tens of mDa on a long peptide), which is large against a ppm fragment
/// window; the shifted ions of the true form would then miss by more than the
/// scorer's tolerance. This is the mass-calibration step of the method. Returns
/// `delta` unchanged when fewer than two shifted ions match.
pub fn refine_offset(
    peptide: &Peptide,
    site0: usize,
    delta: f64,
    peaks: &[(f64, f32)],
    tol_da: f64,
) -> f64 {
    let n = peptide.length();
    if !delta.is_finite() || site0 >= n || n < 2 {
        return delta;
    }
    let base: Vec<f64> = peptide
        .residues
        .iter()
        .map(|a| a.mass + a.mod_.as_ref().map_or(0.0, |m| m.mass_delta))
        .collect();
    let mut diffs: Vec<f64> = Vec::new();
    // b_{i+1} covers residues 0..=i and carries the offset iff i >= site0.
    let mut prefix = 0.0;
    for (i, &m) in base.iter().enumerate().take(n - 1) {
        prefix += m + if i == site0 { delta } else { 0.0 };
        if i < site0 {
            continue;
        }
        for z in 1..=2u8 {
            let zf = f64::from(z);
            let mz = (prefix + zf * PROTON) / zf;
            if let Some(j) = nearest_peak(peaks, mz, tol_da) {
                diffs.push((mz - peaks[j].0) * zf);
            }
        }
    }
    // y over residues i..=n-1 carries the offset iff i <= site0.
    let mut suffix = 0.0;
    for i in (1..n).rev() {
        suffix += base[i] + if i == site0 { delta } else { 0.0 };
        if i > site0 {
            continue;
        }
        for z in 1..=2u8 {
            let zf = f64::from(z);
            let mz = (suffix + H2O + zf * PROTON) / zf;
            if let Some(j) = nearest_peak(peaks, mz, tol_da) {
                diffs.push((mz - peaks[j].0) * zf);
            }
        }
    }
    if diffs.len() < 2 {
        return delta;
    }
    diffs.sort_by(|a, b| a.total_cmp(b));
    let mid = diffs.len() / 2;
    let median = if diffs.len() % 2 == 1 {
        diffs[mid]
    } else {
        0.5 * (diffs[mid - 1] + diffs[mid])
    };
    // theoretical - observed = delta_inferred - delta_true. A correction beyond
    // the precursor's own accuracy is a wrong-site artefact: keep the inferred
    // offset instead.
    if median.abs() > MAX_REFINE_DA {
        return delta;
    }
    delta - median
}

/// A localisation whose offset was re-centred on the winning site's ions.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RefinedLocalization {
    /// The called site (1-based), 0 when the margin is below the reporting
    /// threshold. Used for output only.
    pub site: u16,
    /// The highest-scoring site (1-based), always set. The re-score uses this
    /// one: whether a site is *reportable* is a separate, stricter decision
    /// than which form scores best.
    pub best_site: u16,
    /// The offset re-centred on the winning site (Da, observed - peptide).
    pub delta: f64,
    pub score: f32,
    pub margin: f32,
}

/// Two-stage localisation for an offset that carries the precursor's mass
/// error: for every residue, re-centre the offset on that residue's shifted
/// ions with the coarse window, then score the re-centred form with the tight
/// (fragment) window. The coarse pass alone cannot separate residues — a wide
/// window matches at every site and collapses the margin — so the winner is
/// chosen on the tight score.
pub fn localize_and_refine(
    peptide: &Peptide,
    peaks: &[(f64, f32)],
    charge: u8,
    delta: f64,
    coarse_tol: f64,
    fine_tol: f64,
) -> Option<RefinedLocalization> {
    let n = peptide.length();
    if n < 2 || charge == 0 || !delta.is_finite() || !fine_tol.is_finite() || fine_tol <= 0.0 {
        return None;
    }
    let total: f64 = peaks.iter().map(|&(_, i)| f64::from(i.max(0.0))).sum();
    if total <= 0.0 {
        return None;
    }
    let base: Vec<f64> = peptide
        .residues
        .iter()
        .map(|a| a.mass + a.mod_.as_ref().map_or(0.0, |m| m.mass_delta))
        .collect();
    let max_z = charge.clamp(1, 2);
    let mut used = vec![0u32; peaks.len()];
    let mut generation: u32 = 0;
    let mut best = (0usize, 0.0f64, delta);
    let mut second = 0.0f64;
    for site in 0..n {
        let d = refine_offset(peptide, site, delta, peaks, coarse_tol);
        generation = generation.wrapping_add(1);
        let s = site_score(&base, site, d, peaks, fine_tol, max_z, &mut used, generation);
        if s > best.1 {
            second = best.1;
            best = (site, s, d);
        } else if s > second {
            second = s;
        }
    }
    if best.1 <= 0.0 {
        return None;
    }
    let margin = ((best.1 - second) / best.1) as f32;
    Some(RefinedLocalization {
        site: if margin >= MIN_LOCALIZATION_MARGIN {
            (best.0 + 1) as u16
        } else {
            0
        },
        best_site: (best.0 + 1) as u16,
        delta: best.2,
        score: (best.1 / total).clamp(0.0, 1.0) as f32,
        margin: margin.clamp(0.0, 1.0),
    })
}

/// Localise `delta` (Da, added to one residue) on `peptide` from `peaks`
/// (m/z, intensity). `tol_da` is the fragment match window in Da; `charge` is
/// the PSM's charge, which caps the fragment charge at 2.
pub fn localize_offset(
    peptide: &Peptide,
    peaks: &[(f64, f32)],
    charge: u8,
    delta: f64,
    tol_da: f64,
) -> Option<Localization> {
    let sorted = sorted_peaks(peaks);
    if sorted.len() < 2 {
        return None;
    }
    localize_offset_sorted(peptide, &sorted, charge, delta, tol_da)
}

/// The peaks of one spectrum in m/z order, with non-positive / non-finite
/// entries dropped. Build once per spectrum and reuse for every candidate: the
/// per-candidate localisation re-sorts otherwise.
pub fn sorted_peaks(peaks: &[(f64, f32)]) -> Vec<(f64, f32)> {
    let mut sorted: Vec<(f64, f32)> = peaks
        .iter()
        .copied()
        .filter(|&(mz, i)| mz.is_finite() && mz > 0.0 && i.is_finite())
        .collect();
    sorted.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));
    sorted
}

/// [`localize_offset`] over an already-sorted peak list (see [`sorted_peaks`]).
pub fn localize_offset_sorted(
    peptide: &Peptide,
    sorted: &[(f64, f32)],
    charge: u8,
    delta: f64,
    tol_da: f64,
) -> Option<Localization> {
    let n = peptide.length();
    if n < 2 || charge == 0 || !delta.is_finite() || !tol_da.is_finite() || tol_da <= 0.0 {
        return None;
    }
    if sorted.len() < 2 {
        return None;
    }
    let total: f64 = sorted.iter().map(|&(_, i)| f64::from(i.max(0.0))).sum();
    if total <= 0.0 {
        return None;
    }

    // Residue masses of the unmodified-offset form; `delta` is added to the
    // candidate site only.
    let base: Vec<f64> = peptide
        .residues
        .iter()
        .map(|a| a.mass + a.mod_.as_ref().map_or(0.0, |m| m.mass_delta))
        .collect();
    let max_z = charge.clamp(1, 2);

    // Dedup peaks per site: a peak explains the site once, even when two ions
    // of the same site fall in the same window.
    let mut used = vec![0u32; sorted.len()];
    let mut generation: u32 = 0;

    let mut best = (0usize, 0.0f64);
    let mut second = 0.0f64;
    for site in 0..n {
        generation = generation.wrapping_add(1);
        let s = site_score(
            &base, site, delta, sorted, tol_da, max_z, &mut used, generation,
        );
        if s > best.1 {
            second = best.1;
            best = (site, s);
        } else if s > second {
            second = s;
        }
    }
    if best.1 <= 0.0 {
        return None;
    }
    let margin = ((best.1 - second) / best.1) as f32;
    Some(Localization {
        site: if margin >= MIN_LOCALIZATION_MARGIN {
            (best.0 + 1) as u16
        } else {
            0
        },
        score: (best.1 / total).clamp(0.0, 1.0) as f32,
        margin: margin.clamp(0.0, 1.0),
    })
}

/// Summed intensity of the peaks matched by the b/y ions of `base` with
/// `delta` added to residue index `site`.
#[allow(clippy::too_many_arguments)]
fn site_score(
    base: &[f64],
    site: usize,
    delta: f64,
    peaks: &[(f64, f32)],
    tol_da: f64,
    max_z: u8,
    used: &mut [u32],
    generation: u32,
) -> f64 {
    let n = base.len();
    let mut sum = 0.0f64;
    let mut add = |mz: f64, sum: &mut f64| {
        if !mz.is_finite() || mz <= 0.0 {
            return;
        }
        // Nearest peak to `mz`; sorted, so a binary search lands on the first
        // peak >= mz and the match is at `i` or `i - 1`.
        let i = peaks.partition_point(|p| p.0 < mz);
        for j in [i.saturating_sub(1), i] {
            if j >= peaks.len() {
                continue;
            }
            if (peaks[j].0 - mz).abs() <= tol_da {
                if used[j] != generation {
                    used[j] = generation;
                    *sum += f64::from(peaks[j].1.max(0.0));
                }
                break;
            }
        }
    };

    // b ions: residues 0..=i, i < n-1 (b_n is the precursor).
    let mut prefix = 0.0;
    for (i, &m) in base.iter().enumerate().take(n - 1) {
        prefix += m + if i == site { delta } else { 0.0 };
        for z in 1..=max_z {
            let zf = f64::from(z);
            add((prefix + zf * PROTON) / zf, &mut sum);
        }
    }
    // y ions: residues i..n-1, i > 0.
    let mut suffix = 0.0;
    for i in (1..n).rev() {
        suffix += base[i] + if i == site { delta } else { 0.0 };
        for z in 1..=max_z {
            let zf = f64::from(z);
            add((suffix + H2O + zf * PROTON) / zf, &mut sum);
        }
    }
    sum
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::aa_set::AminoAcidSetBuilder;

    fn peptide(seq: &str) -> Peptide {
        let aa = AminoAcidSetBuilder::new_standard().build().unwrap();
        Peptide::from_str(&format!("K.{seq}.R"), &aa).expect("valid peptide")
    }

    /// b/y ions (charges 1 and 2) of `seq` with `delta` on residue `site`,
    /// as a synthetic spectrum.
    fn spectrum(seq: &str, site: usize, delta: f64) -> Vec<(f64, f32)> {
        let aa = AminoAcidSetBuilder::new_standard().build().unwrap();
        let pep = Peptide::from_str(&format!("K.{seq}.R"), &aa).unwrap();
        let base: Vec<f64> = pep.residues.iter().map(|a| a.mass).collect();
        let n = base.len();
        let mut out = Vec::new();
        let mut prefix = 0.0;
        for (i, &m) in base.iter().enumerate().take(n - 1) {
            prefix += m + if i == site { delta } else { 0.0 };
            out.push((prefix + PROTON, 1000.0));
        }
        let mut suffix = 0.0;
        for (k, &m) in base.iter().enumerate().rev().take(n - 1) {
            suffix += m + if k == site { delta } else { 0.0 };
            out.push((suffix + H2O + PROTON, 1000.0));
        }
        out
    }

    #[test]
    fn localizes_an_offset_to_the_shifted_residue() {
        let seq = "PEPTIDEK";
        let pep = peptide(seq);
        let delta = 79.96633; // phosphorylation mass, as an unknown offset
        let peaks = spectrum(seq, 4, delta);
        let loc = localize_offset(&pep, &peaks, 2, delta, 0.02).expect("localised");
        assert_eq!(loc.site, 5, "1-based position of residue 4");
        // An adjacent site matches every ion except the two that bracket the
        // offset, so the honest margin on an 8-residue peptide is ~1/7.
        assert!(loc.margin > 0.1, "margin {}", loc.margin);
        assert!(loc.score > 0.5, "score {}", loc.score);
    }

    #[test]
    fn an_unshifted_form_has_no_winning_site() {
        let seq = "PEPTIDEK";
        let pep = peptide(seq);
        // delta 0 shifts no ion, so every site scores the same.
        let peaks = spectrum(seq, 4, 0.0);
        let loc = localize_offset(&pep, &peaks, 2, 0.0, 0.02).expect("scores exist");
        assert_eq!(loc.site, 0, "no position separates from the runner-up");
        assert!(
            loc.margin < MIN_LOCALIZATION_MARGIN,
            "margin {}",
            loc.margin
        );
    }

    #[test]
    fn no_matching_peaks_means_no_call() {
        let pep = peptide("PEPTIDEK");
        let peaks = vec![(2000.0, 10.0), (2100.0, 10.0)];
        assert!(localize_offset(&pep, &peaks, 2, 79.96633, 0.02).is_none());
    }

    #[test]
    fn empty_spectrum_is_refused() {
        let pep = peptide("PEPTIDEK");
        assert!(localize_offset(&pep, &[], 2, 79.96633, 0.02).is_none());
    }

    #[test]
    fn two_stage_localisation_survives_a_precursor_error() {
        let seq = "PEPTIDEK";
        let pep = peptide(seq);
        let true_delta = 79.96633;
        let inferred = true_delta + 0.025; // 25 mDa precursor measurement error
        let peaks = sorted_peaks(&spectrum(seq, 4, true_delta));
        let loc = localize_and_refine(&pep, &peaks, 2, inferred, 0.06, 0.02)
            .expect("localised");
        assert_eq!(loc.site, 5);
        assert!(
            (loc.delta - true_delta).abs() < 0.005,
            "refined delta {}",
            loc.delta
        );
        assert!(loc.margin > 0.1, "margin {}", loc.margin);
    }

    #[test]
    fn refine_offset_recovers_the_true_mass_from_the_shifted_ions() {
        let seq = "PEPTIDEK";
        let pep = peptide(seq);
        let true_delta = 15.99491;
        // The precursor-derived offset carries a 25 mDa measurement error.
        let inferred = true_delta + 0.025;
        let peaks = sorted_peaks(&spectrum(seq, 4, true_delta));
        let refined = refine_offset(&pep, 4, inferred, &peaks, 0.06);
        assert!(
            (refined - true_delta).abs() < 0.005,
            "refined {refined} vs {true_delta}"
        );
    }

    #[test]
    fn shifted_peptide_moves_only_the_named_residue() {
        let pep = peptide("PEPTIDEK");
        let delta = 79.96633;
        let shifted = shifted_peptide(&pep, 4, delta);
        assert!((shifted.mass() - pep.mass() - delta).abs() < 1e-9);
        assert!((shifted.residues[3].mass - pep.residues[3].mass).abs() < 1e-9);
        assert!((shifted.residues[4].mass - pep.residues[4].mass - delta).abs() < 1e-9);
    }
}
