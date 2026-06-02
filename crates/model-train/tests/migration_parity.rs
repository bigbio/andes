use model_train::store::{migrate_dir, ModelStore};
use scoring_crate::param_model::Param;
use std::path::Path;

#[test]
fn every_bundled_model_round_trips_via_store() {
    let dir = tempfile::tempdir().unwrap();
    let store_path = dir.path().join("models.parquet");

    // Integration tests run with the crate root as cwd, so use CARGO_MANIFEST_DIR
    // to build an absolute path to the resources/ionstat directory.
    let manifest = std::env!("CARGO_MANIFEST_DIR");
    let ionstat = Path::new(manifest).join("../../resources/ionstat");

    let ids = migrate_dir(&ionstat, &store_path).expect("migrate");
    assert!(ids.len() >= 39, "expected >=39 migrated models, got {}", ids.len());

    let store = ModelStore::open(&store_path).unwrap();
    for (id, file) in &ids {
        let from_binary = Param::load_from_file(file).unwrap();
        let from_store = store.load_param(id).unwrap();
        assert_eq!(from_store, from_binary, "model {id} differs after migration");
    }
}
