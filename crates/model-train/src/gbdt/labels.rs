//! Theoretical-ion oracle for GBDT signal/noise labels. Uses the engine's
//! residue-mass table (the `model` crate) so labels are defined by the same ion
//! chemistry the engine scores. Enumerates b/y (+a, −H2O, −NH3) at charges
//! 1..=max(1, z-1).

use model::amino_acid::AminoAcid;
use model::mass::{H2O, PROTON};

const NH3: f64 = 17.026549;
const CO: f64 = 27.994915;

/// Monoisotopic residue masses for an uppercase AA string; unknown bytes skipped.
fn residue_masses(peptide: &str) -> Vec<f64> {
    peptide
        .bytes()
        .filter_map(|c| AminoAcid::standard(c).map(|aa| aa.mass))
        .collect()
}

/// Sorted theoretical ion m/z for `peptide` at precursor charge `z`
/// (b/y + a + water/ammonia losses, charges 1..=max(1, z-1)).
pub fn theoretical_ion_mzs(peptide: &str, z: u8) -> Vec<f64> {
    let res = residue_masses(peptide);
    let n = res.len();
    let total: f64 = res.iter().sum::<f64>() + H2O;
    let max_z = z.max(1).saturating_sub(1).max(1); // 1..=max(1, z-1)
    let mut mzs: Vec<f64> = Vec::new();
    let mut prefix = 0.0;
    for i in 1..n {
        prefix += res[i - 1];
        let b = prefix + PROTON;           // singly-charged b
        let y = (total - prefix) + PROTON; // singly-charged y
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

/// 1 if a peak matches any theoretical ion of any confident peptide within
/// `tol_da`, else 0. `peptides` is the confident set (union for chimeric IDs).
pub fn label_peaks(peaks: &[(f64, f32)], peptides: &[&str], z: u8, tol_da: f64) -> Vec<u8> {
    let mut theo: Vec<f64> = peptides.iter().flat_map(|p| theoretical_ion_mzs(p, z)).collect();
    theo.sort_by(|a, b| a.partial_cmp(b).unwrap());
    peaks
        .iter()
        .map(|&(mz, _)| {
            let lo = theo.partition_point(|&t| t < mz - tol_da);
            u8::from(lo < theo.len() && theo[lo] <= mz + tol_da)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn b_ion_peak_labeled_signal() {
        // "PEPTIDE" b2 (P+E) singly charged ≈ 227.1026. Place a peak there and a
        // junk peak; expect [signal, noise].
        let theo = theoretical_ion_mzs("PEPTIDE", 2);
        assert!(theo.iter().any(|m| (m - 227.1026).abs() < 0.02), "b2 not enumerated: {theo:?}");
        let peaks = [(227.1026_f64, 500.0_f32), (999.999, 10.0)];
        let labels = label_peaks(&peaks, &["PEPTIDE"], 2, 0.02);
        assert_eq!(labels, vec![1u8, 0]);
    }

    #[test]
    fn chimeric_union_labels_either_peptide() {
        // A peak matching peptide B (but placed where A has no ion) is still
        // signal when both are confident.
        let b = "SAMPLER";
        let b_theo = theoretical_ion_mzs(b, 2);
        let only_b = b_theo[0];
        let peaks = [(only_b, 100.0_f32)];
        let labels = label_peaks(&peaks, &["PEPTIDE", b], 2, 0.02);
        assert_eq!(labels, vec![1u8]);
    }
}
