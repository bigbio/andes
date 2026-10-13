//! The committed AGBD fixtures load with andes's native decoder
//! (`scoring::gbdt_eval::GbdtPeakModel`) and score the rows in
//! `fixtures/parity/agbd.json` to the expected values. `scripts/test_fixtures.py`
//! checks ONNX Runtime against the same table, so the Python parser and the
//! exported ONNX are both tied to the native format rather than to each other.

use std::path::PathBuf;

use scoring::gbdt_eval::GbdtPeakModel;
use serde_json::Value;

fn fixture_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("fixtures")
}

fn parity_table() -> Vec<Value> {
    let text = std::fs::read_to_string(fixture_dir().join("parity/agbd.json")).expect("read");
    serde_json::from_str::<Value>(&text)
        .expect("parse")
        .as_array()
        .expect("array")
        .clone()
}

fn load_native(name: &str) -> (Vec<u8>, GbdtPeakModel) {
    let bytes = std::fs::read(fixture_dir().join(name)).expect("read agbd");
    let model = GbdtPeakModel::from_bytes(&bytes)
        .unwrap_or_else(|e| panic!("native decoder rejects {name}: {e}"));
    (bytes, model)
}

#[test]
fn agbd_fixtures_decode_and_reserialize_byte_identically() {
    let table = parity_table();
    assert!(!table.is_empty());
    for entry in &table {
        let name = entry["agbd"].as_str().unwrap();
        let (bytes, model) = load_native(name);
        assert_eq!(
            model.to_bytes(),
            bytes,
            "{name}: native re-serialization differs"
        );
    }
}

#[test]
fn agbd_fixtures_score_parity_rows_natively() {
    for entry in &parity_table() {
        let name = entry["agbd"].as_str().unwrap();
        let (_, model) = load_native(name);
        let rows = entry["rows"].as_array().unwrap();
        let expected = entry["expected"].as_array().unwrap();
        assert_eq!(rows.len(), expected.len());
        for (row, want) in rows.iter().zip(expected) {
            let x: Vec<f32> = row
                .as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_f64().map_or(f32::NAN, |f| f as f32))
                .collect();
            assert_eq!(x.len(), model.n_features as usize, "{name}");
            let got = model.predict_value(&x);
            let want = want.as_f64().unwrap() as f32;
            assert!(
                (got - want).abs() <= 1e-6,
                "{name} row {x:?}: native {got} vs {want}"
            );
        }
    }
}
