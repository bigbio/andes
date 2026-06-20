//! Theoretical-ion oracle for GBDT signal/noise labels. Uses the engine's
//! MOD-AWARE residue masses (taken from the parsed `Peptide`, i.e. base residue
//! mass + each residue's modification delta) so labels reflect the same ion
//! chemistry the engine scores — including Carbamidomethyl-C, TMT, Oxidation,
//! and terminal mods (which the peptide builder folds onto the end residues).
//! Enumerates b/y (+a, −H2O, −NH3) at charges 1..=max(1, z-1).

use model::mass::{H2O, PROTON};
use model::peptide::Peptide;
use model::Tolerance;

const NH3: f64 = 17.026549;
const CO: f64 = 27.994915;

/// Mod-aware monoisotopic residue masses for a parsed peptide. Each entry is
/// `base_residue_mass + modification_delta`, matching the engine's
/// `residue_mass_with_mod`. No string parsing (the residues are already parsed
/// `AminoAcid`s), so unknown-byte dropping can never misalign the masses.
fn residue_masses(peptide: &Peptide) -> Vec<f64> {
    peptide
        .residues
        .iter()
        .map(|aa| aa.mass + aa.mod_.as_ref().map_or(0.0, |m| m.mass_delta))
        .collect()
}

/// Sorted theoretical ion m/z for `peptide` at precursor charge `z`
/// (b/y + a + water/ammonia losses, charges 1..=max(1, z-1)). Mod-aware: the
/// prefix/total sums use modified residue masses, so b/y ions spanning a
/// modified residue land at the correct m/z (the v1 bug enumerated them at
/// unmodified masses and mislabeled them as noise).
pub fn theoretical_ion_mzs(peptide: &Peptide, z: u8) -> Vec<f64> {
    let res = residue_masses(peptide);
    let n = res.len();
    if n < 2 {
        return Vec::new();
    }
    let total: f64 = res.iter().sum::<f64>() + H2O;
    let max_z = z.max(1).saturating_sub(1).max(1); // 1..=max(1, z-1)
    let mut mzs: Vec<f64> = Vec::new();
    let mut prefix = 0.0;
    for i in 1..n {
        prefix += res[i - 1];
        let b = prefix + PROTON; // singly-charged b (mod-aware prefix)
        let y = (total - prefix) + PROTON; // singly-charged y (mod-aware suffix)
        for base in [b, y] {
            for zz in 1..=max_z {
                let zz = zz as f64;
                mzs.push((base + (zz - 1.0) * PROTON) / zz);
                mzs.push((base - H2O + (zz - 1.0) * PROTON) / zz);
                mzs.push((base - NH3 + (zz - 1.0) * PROTON) / zz);
            }
        }
        // a-ion = b - CO
        for zz in 1..=max_z {
            let zz = zz as f64;
            mzs.push((b - CO + (zz - 1.0) * PROTON) / zz);
        }
    }
    mzs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    mzs
}

/// 1 if a peak matches any theoretical ion of any confident peptide within the
/// per-peak `mme` tolerance (the SAME fragment tolerance the matcher uses — so
/// high-res runs no longer admit chance matches at the old hardcoded 0.5 Da),
/// else 0. `peptides` is the confident set (union for chimeric IDs).
pub fn label_peaks(peaks: &[(f64, f32)], peptides: &[&Peptide], z: u8, mme: &Tolerance) -> Vec<u8> {
    let mut theo: Vec<f64> = peptides.iter().flat_map(|p| theoretical_ion_mzs(p, z)).collect();
    theo.sort_by(|a, b| a.partial_cmp(b).unwrap());
    peaks
        .iter()
        .map(|&(mz, _)| {
            let tol = mme.as_da(mz);
            let lo = theo.partition_point(|&t| t < mz - tol);
            u8::from(lo < theo.len() && theo[lo] <= mz + tol)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::amino_acid::AminoAcid;

    /// Build an unmodified test peptide from an uppercase sequence.
    fn pep(seq: &str) -> Peptide {
        let residues: Vec<AminoAcid> =
            seq.bytes().map(|b| AminoAcid::standard(b).unwrap()).collect();
        Peptide::new(residues, b'K', b'R')
    }

    /// Build a test peptide applying a mass delta to the residue at 1-based `pos`.
    fn pep_with_mod(seq: &str, pos: usize, delta: f64) -> Peptide {
        use model::modification::{ModLocation, Modification, ResidueSpec};
        let mut residues: Vec<AminoAcid> =
            seq.bytes().map(|b| AminoAcid::standard(b).unwrap()).collect();
        let r = residues[pos - 1].residue;
        let m = Modification {
            name: "test".to_string(),
            mass_delta: delta,
            residue: ResidueSpec::Specific(r),
            location: ModLocation::Anywhere,
            fixed: false,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
        };
        residues[pos - 1] = residues[pos - 1].clone().with_mod(m);
        Peptide::new(residues, b'K', b'R')
    }

    #[test]
    fn b_ion_peak_labeled_signal() {
        // "PEPTIDE" b2 (P+E) singly charged ≈ 227.1026. Place a peak there and a
        // junk peak; expect [signal, noise].
        let p = pep("PEPTIDE");
        let theo = theoretical_ion_mzs(&p, 2);
        assert!(theo.iter().any(|m| (m - 227.1026).abs() < 0.02), "b2 not enumerated: {theo:?}");
        let peaks = [(227.1026_f64, 500.0_f32), (999.999, 10.0)];
        let labels = label_peaks(&peaks, &[&p], 2, &Tolerance::Da(0.02));
        assert_eq!(labels, vec![1u8, 0]);
    }

    #[test]
    fn chimeric_union_labels_either_peptide() {
        // A peak matching peptide B (but placed where A has no ion) is still
        // signal when both are confident.
        let a = pep("PEPTIDE");
        let b = pep("SAMPLER");
        let b_theo = theoretical_ion_mzs(&b, 2);
        let only_b = b_theo[0];
        let peaks = [(only_b, 100.0_f32)];
        let labels = label_peaks(&peaks, &[&a, &b], 2, &Tolerance::Da(0.02));
        assert_eq!(labels, vec![1u8]);
    }

    #[test]
    fn mod_shifts_theoretical_ion_mass() {
        // The mod-aware oracle must place b/y ions over a modified residue at the
        // SHIFTED mass (the v1 bug placed them at the unmodified mass). Add a
        // +57.02146 (Carbamidomethyl) to residue 2 of "PEPTIDE": every b-ion that
        // spans residue 2 (b2..b6) shifts up by 57.02146 vs the unmodified peptide.
        let plain = theoretical_ion_mzs(&pep("PEPTIDE"), 2);
        let modded = theoretical_ion_mzs(&pep_with_mod("PEPTIDE", 2, 57.02146), 2);
        // Unmodified b2 ≈ 227.1026; modified b2 ≈ 284.1240 (shifted by the delta).
        assert!(
            modded.iter().any(|m| (m - (227.1026 + 57.02146)).abs() < 0.01),
            "modified b2 not at shifted mass: {modded:?}"
        );
        // A peak at the UNMODIFIED b2 must NOT be labeled signal for the modified
        // peptide (the v1 mod-blind bug would have mislabeled it).
        let labels = label_peaks(
            &[(227.1026_f64, 100.0_f32)],
            &[&pep_with_mod("PEPTIDE", 2, 57.02146)],
            2,
            &Tolerance::Da(0.02),
        );
        assert_eq!(labels, vec![0u8], "unmodified b2 mass must be noise for the modified peptide");
        // Sanity: the plain peptide DOES label that peak as signal.
        let plain_labels =
            label_peaks(&[(227.1026_f64, 100.0_f32)], &[&pep("PEPTIDE")], 2, &Tolerance::Da(0.02));
        assert_eq!(plain_labels, vec![1u8]);
        assert_ne!(plain, modded, "mod must change the theoretical ion set");
    }

    #[test]
    fn ppm_tolerance_is_tight_on_high_res() {
        // With a 20 ppm tolerance (high-res), a peak 0.4 Da off b2 is noise —
        // the old hardcoded 0.5 Da would have mislabeled it as signal.
        let p = pep("PEPTIDE");
        let labels = label_peaks(&[(227.1026_f64 + 0.4, 100.0_f32)], &[&p], 2, &Tolerance::Ppm(20.0));
        assert_eq!(labels, vec![0u8], "0.4 Da off must be noise at 20 ppm");
    }
}
