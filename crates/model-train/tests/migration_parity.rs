use model_train::store::{migrate_dir, write_models, ModelStore};
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

/// Standing all-models fidelity gate.
///
/// The 39 original `.param` source files were migrated into the bundled
/// `resources/ionstat/models.parquet` and removed from the tree, so the
/// migration can no longer be re-validated against the `.param` ground truth in
/// CI. This test instead asserts that EVERY model shipped in the bundled store
/// survives a `Param -> write -> read -> Param` round-trip byte-identically —
/// the store-I/O self-consistency that the "byte-identical for all 39 models"
/// guarantee depends on going forward.
#[test]
fn all_bundled_store_models_round_trip_write_read() {
    let bundled = Path::new(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../resources/ionstat/models.parquet"
    ));
    let store = ModelStore::open(bundled).expect("open bundled models.parquet");
    let ids = store.model_ids();
    assert!(
        ids.len() >= 39,
        "expected the bundled store to ship >=39 models, got {}",
        ids.len()
    );

    let originals: Vec<(String, Param)> = ids
        .iter()
        .map(|id| (id.clone(), store.load_param(id).expect("load bundled model")))
        .collect();

    let dir = tempfile::tempdir().unwrap();
    let tmp = dir.path().join("roundtrip.parquet");
    let refs: Vec<(String, &Param)> = originals.iter().map(|(id, p)| (id.clone(), p)).collect();
    write_models(&tmp, &refs).expect("write round-trip store");

    let reopened = ModelStore::open(&tmp).unwrap();
    for (id, original) in &originals {
        let back = reopened.load_param(id).expect("reload model");
        assert_eq!(
            &back, original,
            "bundled model {id} differs after write -> read round-trip"
        );
    }
}
