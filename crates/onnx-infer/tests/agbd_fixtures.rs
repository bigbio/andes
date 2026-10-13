//! The committed `.agbd` fixtures decode with the native andes evaluator
//! (`scoring::gbdt_eval::GbdtPeakModel::from_bytes`) and give the values the
//! ONNX exports of the same trees are tested against.

use std::fs;
use std::path::PathBuf;

use scoring::gbdt_eval::GbdtPeakModel;

fn load(name: &str) -> GbdtPeakModel {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("fixtures")
        .join(name);
    let bytes = fs::read(&path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    let model = GbdtPeakModel::from_bytes(&bytes)
        .unwrap_or_else(|e| panic!("{name} does not decode natively: {e}"));
    assert_eq!(
        model.to_bytes(),
        bytes,
        "{name} is not byte-identical to its native re-serialisation"
    );
    model
}

#[test]
fn tree_leq_decodes_natively() {
    let m = load("tree_leq.agbd");
    assert_eq!(m.n_features, 1);
    assert_eq!(m.trees.len(), 1);
    assert!(!m.apply_sigmoid);
    assert!(m.iso_x.is_empty() && m.iso_y.is_empty());
    // Same rows and values as `tests/ort_backend.rs::tree_ensemble_regressor`.
    for (x, want) in [(0.2f32, 1.0f32), (0.5, 1.0), (0.9, 3.0)] {
        assert_eq!(m.predict_value(&[x]), want, "x={x}");
    }
}

#[test]
fn peak_gbdt_decodes_natively_with_calibration_footer() {
    let m = load("peak_gbdt.agbd");
    assert_eq!(m.n_features, 3);
    assert_eq!(m.trees.len(), 3);
    assert!(m.apply_sigmoid);
    assert_eq!(m.iso_x, vec![0.0, 0.5, 1.0]);
    assert_eq!(m.iso_y, vec![0.05, 0.4, 0.95]);

    // Raw tree sum (what the ONNX export computes), walked by hand:
    // tree 0: x0=0.2 <= 0.5 -> node 1; x1=0.0 > -2.0 -> leaf 4 = 2.5e-3
    // tree 1: x2=0.0 <= 1e-3 -> node 1; x0=0.2 <= 0.5 -> leaf 3 = 0.42
    // tree 2: leaf 0.125
    let x = [0.2f32, 0.0, 0.0];
    let raw = m.predict_value(&x);
    assert!((raw - (2.5e-3 + 0.42 + 0.125)).abs() < 1e-6, "raw={raw}");

    // The calibrated probability applies sigmoid then the isotonic map.
    let p = 1.0 / (1.0 + (-raw).exp());
    let want = 0.4 + (0.95 - 0.4) * (p - 0.5) / 0.5;
    assert!((m.predict_proba(&x) - want).abs() < 1e-6);
}
