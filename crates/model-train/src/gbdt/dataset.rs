//! One-pass GBDT training dataset builder.
//!
//! Given a list of high-confidence PSMs (each carrying a reference spectrum,
//! peptide sequence string, and charge), build a row-major feature matrix, a
//! binary label vector, and a group-id vector suitable for [`super::train::Dataset`].
//!
//! No peptide identity is ever consulted during feature extraction — features
//! are spectrum-only.  The peptide sequence is only used for theoretical-ion
//! labeling (signal vs. noise ground truth).

use crate::gbdt::labels::label_peaks;
use crate::gbdt::train::Dataset;
use model::peptide::Peptide;
use model::spectrum::Spectrum;
use scoring_crate::peak_features::{extract_peak_features, PeakFeatureCtx, N_FEATURES};
use scoring_crate::scoring::scored_spectrum::ScoredSpectrum;
use scoring_crate::RankScorer;

/// One input PSM row for dataset building.
pub struct PsmRow<'a> {
    /// The spectrum associated with this PSM.
    pub spectrum: &'a Spectrum,
    /// The identified peptide, INCLUDING modifications. The mods are used only
    /// for mod-aware theoretical-ion labeling (signal vs. noise ground truth);
    /// feature extraction stays peptide-agnostic.
    pub peptide: &'a Peptide,
    /// Precursor charge state.
    pub charge: u8,
}

/// FNV-1a 32-bit hash for deterministic group-id generation.
fn fnv1a_32(data: &[u8]) -> u32 {
    let mut h: u32 = 2_166_136_261;
    for &b in data {
        h ^= b as u32;
        h = h.wrapping_mul(16_777_619);
    }
    h
}

/// Compute a stable group id from peptide sequence and charge.
///
/// Rows from the same (peptide, charge) pair share one group id, which
/// prevents leakage across the group-disjoint train/validation split in
/// [`super::train::train_gbdt`].
pub(crate) fn group_id(peptide_seq: &str, charge: u8) -> u32 {
    let mut buf = Vec::with_capacity(peptide_seq.len() + 1);
    buf.extend_from_slice(peptide_seq.as_bytes());
    buf.push(charge);
    fnv1a_32(&buf)
}

/// Build a GBDT [`Dataset`] from `rows` using `scorer` to rank peaks.
///
/// For each PSM:
/// 1. Build a [`ScoredSpectrum`] to get the ranked, filtered peak list.
/// 2. Compute per-spectrum feature context via [`PeakFeatureCtx::for_spectrum`].
/// 3. Extract per-peak features via [`extract_peak_features`].
/// 4. Label peaks as signal (1) or noise (0) via theoretical-ion matching.
/// 5. Accumulate into the flat `x`, `y`, `groups` vectors.
pub fn build_dataset(rows: &[PsmRow<'_>], scorer: &RankScorer) -> Dataset {
    let n_features = N_FEATURES;

    let mut x: Vec<f32> = Vec::new();
    let mut y: Vec<u8> = Vec::new();
    let mut groups: Vec<u32> = Vec::new();

    for row in rows {
        let ss = ScoredSpectrum::new(row.spectrum, scorer, row.charge);
        let (active_peaks, active_ranks) = ss.active_peaks_and_ranks();

        if active_peaks.is_empty() {
            continue;
        }

        let parent_mass = ss.parent_mass();
        let precursor_mz = row.spectrum.precursor_mz;
        let param = scorer.param();

        let ctx = PeakFeatureCtx::for_spectrum(
            precursor_mz,
            row.charge,
            parent_mass,
            active_peaks,
            &param.mme,
        );

        let features = extract_peak_features(active_peaks, active_ranks, &ctx);
        // Mod-aware labels at the SAME fragment tolerance the matcher uses
        // (param.mme), not the old hardcoded 0.5 Da.
        let labels = label_peaks(active_peaks, &[row.peptide], row.charge, &param.mme);

        // Group by unmodified sequence + charge (prevents train/val leakage of
        // the same backbone; conservative across peptidoforms).
        let seq: String = row
            .peptide
            .residues
            .iter()
            .map(|aa| aa.residue as char)
            .collect();
        let gid = group_id(&seq, row.charge);

        for (feat, &label) in features.iter().zip(labels.iter()) {
            x.extend_from_slice(feat);
            y.push(label);
            groups.push(gid);
        }
    }

    Dataset {
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
    use model::spectrum::Spectrum;
    use scoring_crate::RankScorer;
    use std::path::PathBuf;

    fn pep(seq: &str) -> Peptide {
        let residues: Vec<AminoAcid> = seq
            .bytes()
            .map(|b| AminoAcid::standard(b).unwrap())
            .collect();
        Peptide::new(residues, b'K', b'R')
    }

    fn bundled_store_path() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../resources/models")
    }

    fn make_scorer() -> RankScorer {
        let store = ModelStore::open(&bundled_store_path()).expect("bundled model store");
        let param = store
            .load_param("hcd_qexactive_tryp")
            .expect("hcd_qexactive_tryp");
        RankScorer::new(&param)
    }

    fn tiny_spectrum(peaks: Vec<(f64, f32)>, precursor_mz: f64, charge: i32) -> Spectrum {
        Spectrum {
            title: "test".to_string(),
            precursor_mz,
            precursor_intensity: None,
            precursor_charge: Some(charge),
            rt_seconds: None,
            scan: None,
            peaks,
            activation_method: None,
            isolation_lower_offset: None,
            isolation_upper_offset: None,
        }
    }

    #[test]
    fn build_dataset_empty_rows() {
        let scorer = make_scorer();
        let ds = build_dataset(&[], &scorer);
        assert_eq!(ds.x.len(), 0);
        assert_eq!(ds.y.len(), 0);
        assert_eq!(ds.groups.len(), 0);
        assert_eq!(ds.n_features, N_FEATURES);
    }

    #[test]
    fn build_dataset_single_row_produces_n_features_per_peak() {
        let scorer = make_scorer();

        // Tiny synthetic spectrum: 5 peaks in sorted order
        let peaks = vec![
            (100.0_f64, 10.0_f32),
            (200.0, 100.0),
            (300.0, 50.0),
            (400.0, 30.0),
            (500.0, 20.0),
        ];
        let spec = tiny_spectrum(peaks, 600.0, 2);

        let peptide = pep("PEPTIDE");
        let row = PsmRow {
            spectrum: &spec,
            peptide: &peptide,
            charge: 2,
        };

        let ds = build_dataset(&[row], &scorer);
        // Each active peak contributes N_FEATURES features.
        assert!(!ds.y.is_empty(), "must have at least one peak row");
        assert_eq!(ds.x.len(), ds.y.len() * N_FEATURES);
        assert_eq!(ds.groups.len(), ds.y.len());
        // All rows share the same group id.
        let first_gid = ds.groups[0];
        assert!(ds.groups.iter().all(|&g| g == first_gid));
    }

    #[test]
    fn build_dataset_labels_are_binary() {
        let scorer = make_scorer();

        let peaks = vec![(200.0_f64, 100.0_f32), (999.0, 5.0)];
        let spec = tiny_spectrum(peaks, 500.0, 2);

        let peptide = pep("PEPTIDE");
        let row = PsmRow {
            spectrum: &spec,
            peptide: &peptide,
            charge: 2,
        };

        let ds = build_dataset(&[row], &scorer);
        for &label in &ds.y {
            assert!(
                label == 0 || label == 1,
                "label must be 0 or 1, got {label}"
            );
        }
    }

    #[test]
    fn group_ids_differ_for_different_peptides() {
        let g1 = group_id("PEPTIDE", 2);
        let g2 = group_id("SAMPLER", 2);
        let g3 = group_id("PEPTIDE", 3);
        assert_ne!(
            g1, g2,
            "different sequences should give different group ids"
        );
        assert_ne!(g1, g3, "different charges should give different group ids");
    }
}
