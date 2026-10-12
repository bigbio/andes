//! Always-on API / schema contracts for future consumers (andes, quantms-rescoring).
//!
//! These tests do **not** require `--features tract` or ONNX Runtime. They freeze
//! feature names, batching policy, PeptDeep I/O packing, and fixture manifests so
//! later wiring (#92 PIN columns, quantms idparquet) cannot silently drift.

use std::fs;
use std::path::PathBuf;

use onnx_infer::peptdeep::{
    encode_aa_sequence, pack_ms2_inputs, pack_rt_inputs, INSTRUMENT_DEFAULT, MOD_ELEMENT_COUNT,
    NCE_DEFAULT,
};
use onnx_infer::rescoring::{
    CandidateKey, FeatureRow, RescoringBatch, TopNPolicy, ALPHAPEPTDEEP_CORE_FEATURES,
    DEEPLC_CORE_FEATURES, SPECTRUM_CORE_FEATURES,
};
use onnx_infer::{ModelManifest, ModelTask};

fn fixtures() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("fixtures")
}

/// #92 / quantms-rescoring short names — do not rename without a migration note.
#[test]
fn alphapeptdeep_core_feature_names_frozen() {
    assert_eq!(
        ALPHAPEPTDEEP_CORE_FEATURES,
        &[
            "spec_pearson_norm",
            "ionb_pearson_norm",
            "iony_pearson_norm",
            "spec_mse_norm",
            "dotprod_norm",
            "cos_norm",
            "rt_diff",
            "predicted_retention_time",
        ]
    );
}

#[test]
fn deeplc_and_spectrum_feature_names_frozen() {
    assert_eq!(
        DEEPLC_CORE_FEATURES,
        &[
            "observed_retention_time",
            "predicted_retention_time",
            "rt_diff",
        ]
    );
    assert_eq!(
        SPECTRUM_CORE_FEATURES,
        &[
            "snr",
            "spectral_entropy",
            "tic_top_fraction",
            "mz_intensity_weighted_std",
        ]
    );
}

#[test]
fn top_n_policy_rejects_zero_and_filters_per_spectrum() {
    assert!(TopNPolicy::new(0).is_err());
    let policy = TopNPolicy::new(2).unwrap();
    let items = [("a", 1u32), ("a", 2), ("a", 3), ("b", 1)];
    let kept = policy.filter_sorted_by_rank(&items, |x| x.0);
    assert_eq!(kept.len(), 3);
    assert_eq!(kept[2].0, "b");
}

#[test]
fn rescoring_batch_column_order_matches_feature_names() {
    use std::collections::BTreeMap;
    let rows = vec![FeatureRow {
        key: CandidateKey {
            spectrum_ref: "scan=9".into(),
            rank: 1,
            charge: 3,
            peptide: "PEPTIDER".into(),
        },
        values: BTreeMap::from([("spec_pearson_norm".into(), 0.91), ("rt_diff".into(), -1.5)]),
    }];
    let names = ["spec_pearson_norm", "missing_col", "rt_diff"];
    let batch = RescoringBatch::from_rows(rows, &names, "features").unwrap();
    assert_eq!(batch.feature_names, names);
    assert_eq!(batch.features.shape, vec![1, 3]);
    assert_eq!(batch.features.data, vec![0.91, 0.0, -1.5]);
    assert_eq!(batch.keys[0].spectrum_ref, "scan=9");
}

#[test]
fn peptdeep_packing_matches_openms_export_tensor_names() {
    let rt = pack_rt_inputs(&["PEPTIDE"], 16).unwrap();
    assert_eq!(rt.input_sequences.name, "input_sequences");
    assert_eq!(rt.mod_x.name, "mod_x");
    assert_eq!(rt.mod_x.shape, vec![1, 18, MOD_ELEMENT_COUNT]);

    let ms2 = pack_ms2_inputs(
        &["PEPTIDE", "SAMPLER"],
        &[2, 3],
        20,
        NCE_DEFAULT,
        INSTRUMENT_DEFAULT,
    )
    .unwrap();
    assert_eq!(ms2.aa_indices.name, "aa_indices");
    assert_eq!(ms2.charges.name, "charges");
    assert_eq!(ms2.nce.name, "nce");
    assert_eq!(ms2.instrument_indices.name, "instrument_indices");
    assert_eq!(ms2.aa_indices.shape, vec![2, 22]);
    // peptdeep scales charge by 0.1 and NCE by 0.01 before the network.
    assert!((ms2.charges.data[0] - 0.2).abs() < 1e-6 && (ms2.charges.data[1] - 0.3).abs() < 1e-6);
    assert!((ms2.nce.data[0] - NCE_DEFAULT * 0.01).abs() < 1e-6);

    // ASCII indices (A=1 … Z=26) with a zero pad at each terminus.
    let idx = encode_aa_sequence("ACDEFGHIKLMNPQRSTVWY", 20).unwrap();
    assert_eq!(idx.len(), 22);
    assert_eq!((idx[0], idx[1], idx[2], idx[20], idx[21]), (0, 1, 3, 25, 0));
}

#[test]
fn all_fixture_manifests_validate_and_resolve_onnx() {
    let dir = fixtures();
    let mut found = 0usize;
    for entry in fs::read_dir(&dir).unwrap() {
        let path = entry.unwrap().path();
        if path.extension().and_then(|e| e.to_str()) != Some("json") {
            continue;
        }
        let man = ModelManifest::from_path(&path)
            .unwrap_or_else(|e| panic!("manifest {} invalid: {e}", path.display()));
        assert!(
            man.onnx_path().is_file(),
            "missing onnx for {}: {}",
            man.id,
            man.onnx_path().display()
        );
        assert!(!man.id.is_empty());
        found += 1;
    }
    assert!(
        found >= 4,
        "expected ≥4 fixture manifests (double, intensity_mlp, tree_leq, peak_gbdt), got {found}"
    );
}

#[test]
fn fixture_tasks_cover_rescoring_and_tree_ensemble() {
    let mlp = ModelManifest::from_path(fixtures().join("intensity_mlp.json")).unwrap();
    assert_eq!(mlp.task, ModelTask::Rescoring);

    let tree = ModelManifest::from_path(fixtures().join("peak_gbdt.json")).unwrap();
    assert_eq!(tree.task, ModelTask::TreeEnsemble);
    assert_eq!(
        tree.extras.get("n_features").and_then(|v| v.as_u64()),
        Some(3)
    );

    let agbd = fixtures().join("peak_gbdt.agbd");
    let bytes = fs::read(&agbd).unwrap();
    assert!(bytes.starts_with(b"AGBD"), "peak_gbdt.agbd must be AGBD v1");
    assert_eq!(&bytes[4..8], &1u32.to_le_bytes());
}

/// Placeholder gate for a future PeptDeep ONNX drop under `models/peptdeep/`.
/// Ignored in CI until weights are exported locally (see MODELS.md).
#[test]
#[ignore = "requires models/peptdeep/*.onnx — see crates/onnx-infer/MODELS.md"]
fn peptdeep_onnx_present_when_exported() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("models/peptdeep");
    let rt = root.join("peptdeep_rt_dynamic.onnx");
    assert!(
        rt.is_file(),
        "missing {rt:?}; run scripts/export_peptdeep_to_onnx.py"
    );
}
