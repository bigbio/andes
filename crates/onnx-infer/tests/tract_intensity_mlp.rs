//! Tract smoke for the toy intensity MLP fixture (MatMul+Add+ReLU).

#![cfg(feature = "tract")]

use onnx_infer::{BackendKind, NamedTensor, OnnxSession};

fn fixture(name: &str) -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("fixtures")
        .join(name)
}

#[test]
fn intensity_mlp_relu_scores() {
    let mut session = OnnxSession::from_manifest(fixture("intensity_mlp.json"), BackendKind::Tract)
        .expect("load intensity_mlp");
    // rows: [1,0,0,0] → 0.55; [0,0,1,0] → 0.80; [-1,-1,-1,-1] → 0 (ReLU)
    let input = NamedTensor::new(
        "features",
        vec![3, 4],
        vec![
            1., 0., 0., 0., //
            0., 0., 1., 0., //
            -1., -1., -1., -1.,
        ],
    )
    .unwrap();
    let out = session.run_one(input).expect("infer");
    assert_eq!(out.shape, vec![3, 1]);
    assert!((out.data[0] - 0.55).abs() < 1e-5, "got {}", out.data[0]);
    assert!((out.data[1] - 0.80).abs() < 1e-5, "got {}", out.data[1]);
    assert!(out.data[2].abs() < 1e-6, "got {}", out.data[2]);
}

#[test]
fn rescoring_batch_feeds_mlp() {
    use onnx_infer::rescoring::{CandidateKey, FeatureRow, RescoringBatch};
    use std::collections::BTreeMap;

    let names = ["f0", "f1", "f2", "f3"];
    let rows = vec![
        FeatureRow {
            key: CandidateKey {
                spectrum_ref: "s1".into(),
                rank: 1,
                charge: 2,
                peptide: "PEPTIDE".into(),
            },
            values: BTreeMap::from([
                ("f0".into(), 1.0),
                ("f1".into(), 0.0),
                ("f2".into(), 0.0),
                ("f3".into(), 0.0),
            ]),
        },
        FeatureRow {
            key: CandidateKey {
                spectrum_ref: "s1".into(),
                rank: 2,
                charge: 2,
                peptide: "PEPTIDEK".into(),
            },
            values: BTreeMap::from([
                ("f0".into(), 0.0),
                ("f1".into(), 0.0),
                ("f2".into(), 1.0),
                ("f3".into(), 0.0),
            ]),
        },
    ];
    let batch = RescoringBatch::from_rows(rows, &names, "features").unwrap();
    let mut session =
        OnnxSession::from_onnx(fixture("intensity_mlp.onnx"), BackendKind::Tract).unwrap();
    let out = session.run_one(batch.features).unwrap();
    assert_eq!(out.shape[0], 2);
    assert!((out.data[0] - 0.55).abs() < 1e-5);
    assert!((out.data[1] - 0.80).abs() < 1e-5);
}
