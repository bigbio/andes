use model_train::store::{migrate_dir, ModelStore};
use scoring_crate::param_model::Param;
use std::path::Path;

#[test]
fn every_bundled_model_round_trips_via_store() {
    let dir = tempfile::tempdir().unwrap();
    let store_path = dir.path().join("models.parquet");

    // Migrate the local test fixtures (3 representative .param files) instead of
    // the full resources/ionstat directory (which no longer ships .param files —
    // models are now bundled in a single resources/ionstat/models.parquet store).
    let fixtures = Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures"));

    let ids = migrate_dir(fixtures, &store_path).expect("migrate");
    assert!(ids.len() >= 3, "expected >=3 migrated fixture models, got {}", ids.len());

    let store = ModelStore::open(&store_path).unwrap();
    for (id, file) in &ids {
        let from_binary = Param::load_from_file(file).unwrap();
        let from_store = store.load_param(id).unwrap();
        assert_eq!(from_store, from_binary, "model {id} differs after migration");
    }
}
