//! Integration test: tract backend × fixture `double.onnx` (y = 2x).

#![cfg(feature = "tract")]

use std::path::PathBuf;

use onnx_infer::{BackendKind, NamedTensor, OnnxSession};

fn fixture_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("fixtures")
}

#[test]
fn doubles_batched_input_from_onnx_path() {
    let path = fixture_dir().join("double.onnx");
    let mut session = OnnxSession::from_onnx(&path, BackendKind::Tract).expect("load");
    assert_eq!(session.backend_kind(), BackendKind::Tract);

    let input = NamedTensor::new(
        "x",
        vec![2, 4],
        vec![1., 2., 3., 4., 5., 6., 7., 8.],
    )
    .unwrap();
    let outs = session.run(&[input]).expect("run");
    assert_eq!(outs.len(), 1);
    let y = &outs[0];
    assert_eq!(y.shape, vec![2, 4]);
    let expected: Vec<f32> = [1., 2., 3., 4., 5., 6., 7., 8.]
        .into_iter()
        .map(|v| v * 2.0)
        .collect();
    for (a, b) in y.data.iter().zip(expected.iter()) {
        assert!((a - b).abs() < 1e-5, "{a} vs {b}");
    }
}

#[test]
fn loads_via_manifest() {
    let path = fixture_dir().join("double.json");
    let mut session = OnnxSession::from_manifest(&path, BackendKind::Tract).expect("manifest");
    let man = session.manifest().expect("sidecar");
    assert_eq!(man.id, "fixture-double-v1");

    let input = NamedTensor::new("x", vec![1, 4], vec![0.5, 1.0, 1.5, 2.0]).unwrap();
    let y = session.run_one(input).unwrap();
    assert_eq!(y.data, vec![1.0, 2.0, 3.0, 4.0]);
}
