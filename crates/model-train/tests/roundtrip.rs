use model_train::store::{write_models, ModelStore};
use scoring_crate::param_model::Param;
use std::path::Path;

fn fixture() -> Param {
    // Load from the local test fixtures directory (not the bundled resources).
    let param_path = Path::new(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/CID_TOF_aLP.param"
    ));
    Param::load_from_file(param_path).expect("load fixture CID_TOF_aLP.param")
}

fn fixture2() -> Param {
    let param_path = Path::new(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/HCD_TOF_aLP.param"
    ));
    Param::load_from_file(param_path).expect("load fixture HCD_TOF_aLP.param")
}

#[test]
fn write_creates_a_nonempty_parquet_file() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("models.parquet");
    write_models(&path, &[("cid_tof_alp".to_string(), &fixture())]).unwrap();
    assert!(std::fs::metadata(&path).unwrap().len() > 0);
}

#[test]
fn roundtrip_param_is_equal() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("models.parquet");
    let original = fixture();
    write_models(&path, &[("m".to_string(), &original)]).unwrap();
    let store = ModelStore::open(&path).unwrap();
    let loaded = store.load_param("m").unwrap();
    assert_eq!(loaded, original, "round-tripped Param must equal the original");
}

#[test]
fn roundtrip_two_models_isolated() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("models.parquet");
    let orig1 = fixture();
    let orig2 = fixture2();
    write_models(
        &path,
        &[
            ("model_a".to_string(), &orig1),
            ("model_b".to_string(), &orig2),
        ],
    )
    .unwrap();
    let store = ModelStore::open(&path).unwrap();
    assert_eq!(store.model_ids().len(), 2);
    let loaded1 = store.load_param("model_a").unwrap();
    let loaded2 = store.load_param("model_b").unwrap();
    assert_eq!(loaded1, orig1, "model_a round-trip failed");
    assert_eq!(loaded2, orig2, "model_b round-trip failed");
}
