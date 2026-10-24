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
use scoring_crate::peak_features::{extract_peak_features, PeakFeatureCtx, FEATURE_NAMES, N_FEATURES};
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

/// Compute a stable group id from a stripped peptide sequence.
///
/// Every row of the same backbone shares one group id regardless of charge
/// state or modifications, and I/L are treated as the same residue (they are
/// isobaric, so `PEPTIDE` and `PEPTLDE` are the same spectral evidence). This
/// keeps the group-disjoint train/validation split in
/// [`super::train::train_gbdt`] peptide-disjoint.
pub(crate) fn group_id(peptide_seq: &str) -> u32 {
    let buf: Vec<u8> = peptide_seq
        .bytes()
        .map(|b| if b == b'I' { b'L' } else { b })
        .collect();
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

        // Group by unmodified sequence (prevents train/val leakage of the same
        // backbone across charge states and peptidoforms).
        let seq: String = row
            .peptide
            .residues
            .iter()
            .map(|aa| aa.residue as char)
            .collect();
        let gid = group_id(&seq);

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
        feature_names: FEATURE_NAMES.iter().map(|n| n.to_string()).collect(),
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
            isolation_target_mz: None,
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
        let g1 = group_id("PEPTIDE");
        let g2 = group_id("SAMPLER");
        assert_ne!(
            g1, g2,
            "different sequences should give different group ids"
        );
    }

    #[test]
    fn group_ids_fold_isoleucine_to_leucine() {
        assert_eq!(group_id("PEPTIDEK"), group_id("PEPTLDEK"));
        assert_eq!(group_id("IIAK"), group_id("LLAK"));
    }

    #[test]
    fn build_dataset_groups_one_peptide_across_charge_states() {
        // The same backbone observed at 2+ and 3+ must share one group id so
        // the group-disjoint train/validation split cannot place it on both
        // sides.
        let scorer = make_scorer();
        let peaks = vec![
            (200.0_f64, 100.0_f32),
            (300.0, 50.0),
            (400.0, 30.0),
            (500.0, 20.0),
        ];
        let spec_z2 = tiny_spectrum(peaks.clone(), 600.0, 2);
        let spec_z3 = tiny_spectrum(peaks, 400.0, 3);
        let peptide = pep("PEPTIDE");
        let rows = [
            PsmRow {
                spectrum: &spec_z2,
                peptide: &peptide,
                charge: 2,
            },
            PsmRow {
                spectrum: &spec_z3,
                peptide: &peptide,
                charge: 3,
            },
        ];
        let ds = build_dataset(&rows, &scorer);
        let z2_only = build_dataset(&rows[..1], &scorer);
        assert!(
            ds.groups.len() > z2_only.groups.len(),
            "the 3+ row must contribute peaks"
        );
        let first = ds.groups[0];
        assert!(ds.groups.iter().all(|&g| g == first));
    }
}
