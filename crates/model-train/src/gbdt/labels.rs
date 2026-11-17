use model::mass::{H2O, PROTON};
use model::peptide::Peptide;
use model::Tolerance;

const NH3: f64 = 17.026549;
const CO: f64 = 27.994915;

fn residue_masses(peptide: &Peptide) -> Vec<f64> {
    peptide
        .residues
        .iter()
        .map(|aa| aa.mass + aa.mod_.as_ref().map_or(0.0, |m| m.mass_delta))
        .collect()
}

pub fn theoretical_ion_mzs(peptide: &Peptide, z: u8) -> Vec<f64> {
    let res = residue_masses(peptide);
    let n = res.len();
    if n < 2 {
        return Vec::new();
    }
    let total: f64 = res.iter().sum::<f64>() + H2O;
    let max_z = z.max(1).saturating_sub(1).max(1);
    let mut mzs: Vec<f64> = Vec::new();
    let mut prefix = 0.0;
    for i in 1..n {
        prefix += res[i - 1];
        let b = prefix + PROTON;
        let y = (total - prefix) + PROTON;
        for base in [b, y] {
            for zz in 1..=max_z {
                let zz = zz as f64;
                mzs.push((base + (zz - 1.0) * PROTON) / zz);
                mzs.push((base - H2O + (zz - 1.0) * PROTON) / zz);
                mzs.push((base - NH3 + (zz - 1.0) * PROTON) / zz);
            }
        }
        for zz in 1..=max_z {
            let zz = zz as f64;
            mzs.push((b - CO + (zz - 1.0) * PROTON) / zz);
        }
    }
    mzs.sort_by(f64::total_cmp);
    mzs
}

pub fn label_peaks(peaks: &[(f64, f32)], peptides: &[&Peptide], z: u8, mme: &Tolerance) -> Vec<u8> {
    let mut theo: Vec<f64> = peptides
        .iter()
        .flat_map(|p| theoretical_ion_mzs(p, z))
        .collect();
    theo.sort_by(f64::total_cmp);
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

    fn pep(seq: &str) -> Peptide {
        let residues: Vec<AminoAcid> = seq
            .bytes()
            .map(|b| AminoAcid::standard(b).unwrap())
            .collect();
        Peptide::new(residues, b'K', b'R')
    }

    fn pep_with_mod(seq: &str, pos: usize, delta: f64) -> Peptide {
        use model::modification::{ModLocation, Modification, ResidueSpec};
        let mut residues: Vec<AminoAcid> = seq
            .bytes()
            .map(|b| AminoAcid::standard(b).unwrap())
            .collect();
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
            annotated_site_mass: None,
            blocks_cleavage: false,
            parts: Vec::new(),
            diagnostic_ions: Vec::new(),
        };
        residues[pos - 1] = residues[pos - 1].clone().with_mod(m);
        Peptide::new(residues, b'K', b'R')
    }

    #[test]
    fn b_ion_peak_labeled_signal() {
        let p = pep("PEPTIDE");
        let theo = theoretical_ion_mzs(&p, 2);
        assert!(
            theo.iter().any(|m| (m - 227.1026).abs() < 0.02),
            "b2 not enumerated: {theo:?}"
        );
        let peaks = [(227.1026_f64, 500.0_f32), (999.999, 10.0)];
        let labels = label_peaks(&peaks, &[&p], 2, &Tolerance::Da(0.02));
        assert_eq!(labels, vec![1u8, 0]);
    }

    #[test]
    fn chimeric_union_labels_either_peptide() {
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
        let plain = theoretical_ion_mzs(&pep("PEPTIDE"), 2);
        let modded = theoretical_ion_mzs(&pep_with_mod("PEPTIDE", 2, 57.02146), 2);
        assert!(
            modded
                .iter()
                .any(|m| (m - (227.1026 + 57.02146)).abs() < 0.01),
            "modified b2 not at shifted mass: {modded:?}"
        );
        let labels = label_peaks(
            &[(227.1026_f64, 100.0_f32)],
            &[&pep_with_mod("PEPTIDE", 2, 57.02146)],
            2,
            &Tolerance::Da(0.02),
        );
        assert_eq!(
            labels,
            vec![0u8],
            "unmodified b2 mass must be noise for the modified peptide"
        );
        let plain_labels = label_peaks(
            &[(227.1026_f64, 100.0_f32)],
            &[&pep("PEPTIDE")],
            2,
            &Tolerance::Da(0.02),
        );
        assert_eq!(plain_labels, vec![1u8]);
        assert_ne!(plain, modded, "mod must change the theoretical ion set");
    }

    #[test]
    fn ppm_tolerance_is_tight_on_high_res() {
        let p = pep("PEPTIDE");
        let labels = label_peaks(
            &[(227.1026_f64 + 0.4, 100.0_f32)],
            &[&p],
            2,
            &Tolerance::Ppm(20.0),
        );
        assert_eq!(labels, vec![0u8], "0.4 Da off must be noise at 20 ppm");
    }
}
