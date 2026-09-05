//! Regression training dataset builder for per-fragment intensity prediction.
//!
//! Pairs per-fragment features (from [`scoring_crate::frag_features`]) with the
//! observed log-relative intensity, for every matched b/y ion in each gold PSM.
//! Only MATCHED ions contribute rows — unmatched ions are skipped so the model
//! trains only on actual observations.
//!
//! The target `y` is `ln(obs_intensity / base_peak)` — the SAME log-relative-
//! intensity space that `IntensityModel::predict_log_rel` uses, making the
//! trained regression GBDT a drop-in for that prediction path.

use scoring_crate::frag_features::{extract_frag_features, N_FRAG_FEATURES};
use scoring_crate::scoring::fragment_ions::predict_by_ions;
use scoring_crate::scoring::scored_spectrum::ScoredSpectrum;
use scoring_crate::RankScorer;

use crate::gbdt::dataset::{group_id, PsmRow};
use crate::gbdt::train::RegressionDataset;

/// Build a [`RegressionDataset`] from `rows` using `scorer` to rank peaks.
///
/// For each PSM:
/// 1. Build a [`ScoredSpectrum`] and compute the base peak (max active-peak
///    intensity). Skip if base_peak ≤ 0 or peptide length < 2.
/// 2. Enumerate charge-1 b/y ions via [`predict_by_ions`].
/// 3. For each ion with an observed peak within `scorer.feature_match_tolerance()`
///    (the SAME tolerance the strong-score cosine uses at serve), push one
///    regression row: features from [`extract_frag_features`] and target
///    `y = ln(obs_intensity / base_peak)`.
/// 4. Unmatched ions contribute no row.
pub fn build_frag_dataset(rows: &[PsmRow<'_>], scorer: &RankScorer) -> RegressionDataset {
    let n_features = N_FRAG_FEATURES;

    let mut x: Vec<f32> = Vec::new();
    let mut y: Vec<f32> = Vec::new();
    let mut groups: Vec<u32> = Vec::new();

    for row in rows {
        let ss = ScoredSpectrum::new(row.spectrum, scorer, row.charge);

        // Base peak: max intensity across all active peaks.
        let base_peak = ss
            .dump_active_peaks()
            .iter()
            .map(|(_, _, intensity)| *intensity as f64)
            .fold(0.0_f64, f64::max);

        if base_peak <= 0.0 || row.peptide.length() < 2 {
            continue;
        }

        // Match observed peaks at the SAME tolerance the strong-score cosine
        // uses at serve time (feature_match_tolerance: 20 ppm high-res / 0.5 Da
        // low-res), NOT the coarse param.mme rank-binning tolerance. On high-res
        // data a 0.5 Da window pairs spurious junk peaks to non-fragmenting
        // ions, corrupting the intensity targets and skewing train vs serve.
        let feat_tol = scorer.feature_match_tolerance();
        let seq: String = row
            .peptide
            .residues
            .iter()
            .map(|aa| aa.residue as char)
            .collect();
        let gid = group_id(&seq, row.charge);

        // 1+ AND 2+ fragments — MUST match the serve range in strong_score.rs.
        let predicted = predict_by_ions(row.peptide, 1..=2);

        for ion in &predicted {
            let tol_da = feat_tol.as_da(ion.mz);
            if let Some((_, obs_intensity, _)) = ss.nearest_peak_full(ion.mz, tol_da) {
                let log_rel = ((obs_intensity as f64 / base_peak).max(1e-12)).ln() as f32;

                let feats = extract_frag_features(
                    row.peptide,
                    ion.kind,
                    ion.position,
                    row.charge,
                    ion.charge,
                    0.0,
                );
                x.extend_from_slice(&feats);
                y.push(log_rel);
                groups.push(gid);
            }
        }
    }

    RegressionDataset {
        x,
        y,
        groups,
        n_features,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ModelStore;
    use model::amino_acid::AminoAcid;
    use model::peptide::Peptide;
    use model::spectrum::Spectrum;
    use scoring_crate::scoring::fragment_ions::predict_by_ions;
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

    fn spectrum(peaks: Vec<(f64, f32)>, mz: f64, z: i32) -> Spectrum {
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
        }
    }

    #[test]
    fn frag_dataset_emits_one_row_per_matched_ion() {
        let scorer = make_scorer();
        let p = pep("PEPTIDEK"); // n=8
        let ions = predict_by_ions(&p, 1..=1);
        // place observed peaks at the m/z of the first two b/y ions, plus a junk peak
        let mut peaks: Vec<(f64, f32)> = ions.iter().take(2).map(|i| (i.mz, 1000.0_f32)).collect();
        peaks.push((50.0, 10.0)); // junk + ensures >=... peaks
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let spec = spectrum(peaks, 500.0, 2);
        let row = PsmRow {
            spectrum: &spec,
            peptide: &p,
            charge: 2,
        };
        let ds = build_frag_dataset(&[row], &scorer);
        assert_eq!(ds.n_features, scoring_crate::frag_features::N_FRAG_FEATURES);
        assert!(
            ds.y.len() >= 2,
            "at least the 2 placed ions matched, got {}",
            ds.y.len()
        );
        assert_eq!(ds.x.len(), ds.y.len() * ds.n_features);
        assert_eq!(ds.groups.len(), ds.y.len());
        assert!(ds.y.iter().all(|v| v.is_finite()), "targets finite");
    }
}
