use model::peptide::Peptide;
use scoring_crate::ion_features::{extract_ion_features, ION_FEATURE_NAMES, N_ION_FEATURES};
use scoring_crate::scoring::fragment_ions::predict_by_ions;
use scoring_crate::scoring::scored_spectrum::ScoredSpectrum;
use scoring_crate::RankScorer;

use crate::gbdt::dataset::{group_id, PsmRow};
use crate::gbdt::train::Dataset;

fn reverse_peptide(p: &Peptide) -> Peptide {
    let mut residues = p.residues.clone();
    residues.reverse();
    Peptide::new(residues, p.pre, p.post)
}

pub fn build_ion_dataset(rows: &[PsmRow<'_>], scorer: &RankScorer) -> Dataset {
    let n_features = N_ION_FEATURES;
    let mut x: Vec<f32> = Vec::new();
    let mut y: Vec<u8> = Vec::new();
    let mut groups: Vec<u32> = Vec::new();

    let feat_tol = scorer.feature_match_tolerance();
    let is_ppm = matches!(feat_tol, model::tolerance::Tolerance::Ppm(_));
    let tol_val = feat_tol.raw_value();

    for row in rows {
        if row.peptide.length() < 2 {
            continue;
        }
        let ss = ScoredSpectrum::new(row.spectrum, scorer, row.charge);
        let seq: String = row
            .peptide
            .residues
            .iter()
            .map(|aa| aa.residue as char)
            .collect();
        let gid = group_id(&seq);
        let decoy = reverse_peptide(row.peptide);

        for (pep, label) in [(row.peptide, 1u8), (&decoy, 0u8)] {
            for ion in predict_by_ions(pep, 1..=2) {
                let td = if is_ppm {
                    ion.mz * tol_val / 1e6
                } else {
                    tol_val
                };
                if ss.nearest_peak_full(ion.mz, td).is_some() {
                    let f = extract_ion_features(
                        pep,
                        &ss,
                        ion.kind,
                        ion.position,
                        row.charge,
                        ion.charge,
                        tol_val,
                        is_ppm,
                    );
                    x.extend_from_slice(&f);
                    y.push(label);
                    groups.push(gid);
                }
            }
        }
    }

    Dataset {
        x,
        y,
        groups,
        n_features,
        feature_names: ION_FEATURE_NAMES.iter().map(|n| n.to_string()).collect(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ModelStore;
    use model::amino_acid::AminoAcid;
    use model::spectrum::Spectrum;
    use std::path::PathBuf;

    fn make_scorer() -> RankScorer {
        let store = ModelStore::open(
            &PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../resources/models"),
        )
        .unwrap();
        RankScorer::new(&store.load_param("hcd_qexactive_tryp").unwrap())
    }

    fn pep(seq: &str) -> Peptide {
        Peptide::new(
            seq.bytes()
                .map(|b| AminoAcid::standard(b).unwrap())
                .collect(),
            b'K',
            b'R',
        )
    }

    fn spectrum(mut peaks: Vec<(f64, f32)>, mz: f64, z: i32) -> Spectrum {
        peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
        Spectrum {
            title: "test".to_string(),
            precursor_mz: mz,
            precursor_intensity: None,
            precursor_charge: Some(z),
            rt_seconds: None,
            scan: None,
            peaks,
            activation_method: None,
            isolation_lower_offset: None,
            isolation_upper_offset: None,
            isolation_target_mz: None,
            mobility_window: None,
        }
    }

    #[test]
    fn builds_positive_and_negative_rows_with_n_ion_features() {
        let scorer = make_scorer();
        let p = pep("PEPTIDEK");
        let decoy = reverse_peptide(&p);
        let mut peaks: Vec<(f64, f32)> = Vec::new();
        for ion in predict_by_ions(&p, 1..=1) {
            peaks.push((ion.mz, 1000.0));
        }
        for ion in predict_by_ions(&decoy, 1..=1) {
            peaks.push((ion.mz, 500.0));
        }
        peaks.push((50.0, 5.0));
        let spec = spectrum(peaks, 500.0, 2);
        let row = PsmRow {
            spectrum: &spec,
            peptide: &p,
            charge: 2,
        };

        let ds = build_ion_dataset(&[row], &scorer);
        assert_eq!(ds.n_features, N_ION_FEATURES);
        assert!(ds.y.contains(&1), "target-matched positives present");
        assert!(ds.y.contains(&0), "decoy-matched negatives present");
        assert_eq!(ds.x.len(), ds.y.len() * ds.n_features);
        assert_eq!(ds.groups.len(), ds.y.len());
    }

    #[test]
    fn reverse_peptide_preserves_mass_and_length() {
        let p = pep("PEPTIDEK");
        let r = reverse_peptide(&p);
        assert_eq!(r.length(), p.length());
        assert!(
            (r.mass() - p.mass()).abs() < 1e-6,
            "reverse keeps precursor mass"
        );
    }
}
