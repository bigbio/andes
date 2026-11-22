use model::mass::{H2O, PROTON};
use model::peptide::Peptide;
use model::stats::{median_mut, Even};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Localization {
    pub site: u16,
    pub score: f32,
    pub margin: f32,
}

pub const MIN_LOCALIZATION_MARGIN: f32 = 0.02;

pub const MAX_REFINE_DA: f64 = 0.1;

pub fn shifted_peptide(peptide: &Peptide, site0: usize, delta: f64) -> Peptide {
    let mut residues = peptide.residues.clone();
    if let Some(aa) = residues.get_mut(site0) {
        aa.mass += delta;
    }
    Peptide::new(residues, peptide.pre, peptide.post)
}

fn nearest_peak(peaks: &[(f64, f32)], mz: f64, tol_da: f64) -> Option<usize> {
    let i = peaks.partition_point(|p| p.0 < mz);
    [i.saturating_sub(1), i]
        .into_iter()
        .find(|&j| j < peaks.len() && (peaks[j].0 - mz).abs() <= tol_da)
}

/// Refines a precursor offset from the localised fragment ions: the median mass error of the
/// ions that carry the offset, matched within `tol_da`, is subtracted, up to `MAX_REFINE_DA`.
/// A matching window wider than `MAX_REFINE_DA` (low-resolution ion-trap MS2) admits fragment
/// errors larger than any correction the refinement may apply, so the median is noise at that
/// scale and the offset measured on the precursor is kept.
pub fn refine_offset(
    peptide: &Peptide,
    site0: usize,
    delta: f64,
    peaks: &[(f64, f32)],
    tol_da: f64,
) -> f64 {
    let n = peptide.length();
    if !delta.is_finite() || site0 >= n || n < 2 || !(tol_da <= MAX_REFINE_DA) {
        return delta;
    }
    let base: Vec<f64> = peptide
        .residues
        .iter()
        .map(|a| a.mass + a.mod_.as_ref().map_or(0.0, |m| m.mass_delta))
        .collect();
    let mut diffs: Vec<f64> = Vec::new();
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
    let median = median_mut(&mut diffs, Even::Mean).unwrap_or(0.0);
    if median.abs() > MAX_REFINE_DA {
        return delta;
    }
    delta - median
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RefinedLocalization {
    pub site: u16,
    pub best_site: u16,
    pub delta: f64,
    pub score: f32,
    pub margin: f32,
}

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
        let s = site_score(
            &base, site, d, peaks, fine_tol, max_z, &mut used, generation,
        );
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

pub fn sorted_peaks(peaks: &[(f64, f32)]) -> Vec<(f64, f32)> {
    let mut sorted: Vec<(f64, f32)> = peaks
        .iter()
        .copied()
        .filter(|&(mz, i)| mz.is_finite() && mz > 0.0 && i.is_finite())
        .collect();
    sorted.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));
    sorted
}

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

    let base: Vec<f64> = peptide
        .residues
        .iter()
        .map(|a| a.mass + a.mod_.as_ref().map_or(0.0, |m| m.mass_delta))
        .collect();
    let max_z = charge.clamp(1, 2);

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

    let mut prefix = 0.0;
    for (i, &m) in base.iter().enumerate().take(n - 1) {
        prefix += m + if i == site { delta } else { 0.0 };
        for z in 1..=max_z {
            let zf = f64::from(z);
            add((prefix + zf * PROTON) / zf, &mut sum);
        }
    }
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
        let delta = 79.96633;
        let peaks = spectrum(seq, 4, delta);
        let loc = localize_offset(&pep, &peaks, 2, delta, 0.02).expect("localised");
        assert_eq!(loc.site, 5, "1-based position of residue 4");
        assert!(loc.margin > 0.1, "margin {}", loc.margin);
        assert!(loc.score > 0.5, "score {}", loc.score);
    }

    #[test]
    fn an_unshifted_form_has_no_winning_site() {
        let seq = "PEPTIDEK";
        let pep = peptide(seq);
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
        let inferred = true_delta + 0.025;
        let peaks = sorted_peaks(&spectrum(seq, 4, true_delta));
        let loc = localize_and_refine(&pep, &peaks, 2, inferred, 0.06, 0.02).expect("localised");
        assert_eq!(loc.site, 5);
        assert!(
            (loc.delta - true_delta).abs() < 0.005,
            "refined delta {}",
            loc.delta
        );
        assert!(loc.margin > 0.1, "margin {}", loc.margin);
    }

    #[test]
    fn a_low_resolution_fragment_window_keeps_the_precursor_offset() {
        let seq = "PEPTIDEK";
        let pep = peptide(seq);
        let true_delta = 15.99491;
        let inferred = true_delta + 0.025;
        let peaks = sorted_peaks(&spectrum(seq, 4, true_delta));
        assert_eq!(refine_offset(&pep, 4, inferred, &peaks, 0.5), inferred);
        assert_eq!(refine_offset(&pep, 4, inferred, &peaks, f64::NAN), inferred);
    }

    #[test]
    fn refine_offset_recovers_the_true_mass_from_the_shifted_ions() {
        let seq = "PEPTIDEK";
        let pep = peptide(seq);
        let true_delta = 15.99491;
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
