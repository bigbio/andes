//! Integration test: tract reports ONNX output tensor names, not node names.

#![cfg(feature = "tract")]

use std::path::PathBuf;

use onnx_infer::{BackendKind, NamedTensor, OnnxSession};

fn fixture_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("fixtures")
}

/// `identity_named.onnx` is `Identity(x) -> y` in a node named `copy_node`.
#[test]
fn output_is_named_after_the_graph_output_not_the_node() {
    let path = fixture_dir().join("identity_named.onnx");
    let mut session = OnnxSession::from_onnx(&path, BackendKind::Tract).expect("load");
    assert_eq!(session.output_names(), ["y".to_string()]);

    let input = NamedTensor::new("x", vec![1, 1], vec![2.0]).unwrap();
    let outs = session.run(&[input]).expect("run");
    assert_eq!(outs.len(), 1);
    assert_eq!(outs[0].name, "y");
    assert_eq!(outs[0].data, vec![2.0]);
}
