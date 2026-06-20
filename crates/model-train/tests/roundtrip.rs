use model_train::store::{write_models, write_models_with_gbdt, write_all_models_with_sources_and_gbdt_pub, ModelStore};
use model_train::store::SourceLedger;
use model_train::counts::CountStats;
use rustc_hash::FxHashMap;
use scoring_crate::param_model::{FragmentOffsetFrequency, IonType, Param, Partition, SpecDataType};
use model::activation::ActivationMethod;
use model::instrument::InstrumentType;
use model::protocol::Protocol;
use model::tolerance::Tolerance;
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

/// Construct a minimal `Param` that contains loss-ion entries (loss_class != 0)
/// in both `rank_dist_table` and `frag_off_table`.
fn param_with_loss_ions() -> Param {
    let part = Partition { charge: 2, parent_mass: 1500.0, seg_num: 0 };

    // intact prefix ion (loss_class = 0)
    let intact = IonType::Prefix { charge: 1, offset_bits: 1.0_f32.to_bits(), loss_class: 0 };
    // phospho-loss prefix ion (loss_class = 2)
    let phospho_loss = IonType::Prefix { charge: 1, offset_bits: 2.0_f32.to_bits(), loss_class: 2 };
    // generic-loss suffix ion (loss_class = 255)
    let generic_loss = IonType::Suffix { charge: 1, offset_bits: 3.0_f32.to_bits(), loss_class: 255 };
    let noise = IonType::Noise;

    let mut ion_table: FxHashMap<IonType, Vec<f32>> = FxHashMap::default();
    ion_table.insert(intact, vec![0.6, 0.3, 0.05, 0.001]);
    ion_table.insert(phospho_loss, vec![0.4, 0.2, 0.03, 0.001]);
    ion_table.insert(generic_loss, vec![0.3, 0.15, 0.02, 0.001]);
    ion_table.insert(noise, vec![0.1, 0.2, 0.3, 0.4]);

    let mut rank_dist_table: FxHashMap<Partition, FxHashMap<IonType, Vec<f32>>> = FxHashMap::default();
    rank_dist_table.insert(part, ion_table);

    let mut frag_off_table: FxHashMap<Partition, Vec<FragmentOffsetFrequency>> = FxHashMap::default();
    frag_off_table.insert(part, vec![
        FragmentOffsetFrequency { ion_type: intact, frequency: 0.7 },
        FragmentOffsetFrequency { ion_type: phospho_loss, frequency: 0.5 },
        FragmentOffsetFrequency { ion_type: generic_loss, frequency: 0.3 },
    ]);

    let mut p = Param {
        version: 10001,
        data_type: SpecDataType {
            activation: ActivationMethod::HCD,
            instrument: InstrumentType::QExactive,
            enzyme: None,
            protocol: Protocol::Automatic,
        },
        mme: Tolerance::Ppm(20.0),
        apply_deconvolution: false,
        deconvolution_error_tolerance: 0.0,
        charge_hist: vec![(2, 100)],
        min_charge: 2,
        max_charge: 2,
        num_segments: 1,
        partitions: vec![part],
        num_precursor_off: 0,
        precursor_off_map: FxHashMap::default(),
        frag_off_table,
        max_rank: 3,
        rank_dist_table,
        error_scaling_factor: 0,
        ion_err_dist_table: FxHashMap::default(),
        noise_err_dist_table: FxHashMap::default(),
        ion_existence_table: FxHashMap::default(),
        partition_ion_types_cache: FxHashMap::default(),
        gbdt_peak_model: None,
        frag_intensity_model: None,
        rich_ion_model: None,
    };
    p.rebuild_cache();
    p
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

/// A model with loss ions (loss_class != 0) must survive a store round-trip
/// with loss_class preserved.
#[test]
fn loss_class_survives_store_round_trip() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("loss_models.parquet");
    let original = param_with_loss_ions();

    write_models(&path, &[("loss_model".to_string(), &original)]).unwrap();
    let store = ModelStore::open(&path).unwrap();
    let loaded = store.load_param("loss_model").unwrap();

    assert_eq!(loaded, original, "round-tripped loss-ion Param must equal the original");

    // Also directly verify that the specific loss_class values are preserved.
    let part = Partition { charge: 2, parent_mass: 1500.0, seg_num: 0 };

    // Check rank_dist_table ion keys
    let ion_map = loaded.rank_dist_table.get(&part)
        .expect("partition must exist in rank_dist_table");
    let phospho = IonType::Prefix { charge: 1, offset_bits: 2.0_f32.to_bits(), loss_class: 2 };
    let generic  = IonType::Suffix { charge: 1, offset_bits: 3.0_f32.to_bits(), loss_class: 255 };
    assert!(ion_map.contains_key(&phospho),
        "rank_dist_table must contain phospho-loss ion (loss_class=2); got keys: {:?}",
        ion_map.keys().collect::<Vec<_>>());
    assert!(ion_map.contains_key(&generic),
        "rank_dist_table must contain generic-loss ion (loss_class=255); got keys: {:?}",
        ion_map.keys().collect::<Vec<_>>());

    // Check frag_off_table ion types
    let frags = loaded.frag_off_table.get(&part)
        .expect("partition must exist in frag_off_table");
    let all_loss_classes: Vec<u8> = frags.iter().map(|f| f.ion_type.loss_class()).collect();
    assert!(all_loss_classes.contains(&2),
        "frag_off_table must contain an entry with loss_class=2; got {:?}", all_loss_classes);
    assert!(all_loss_classes.contains(&255),
        "frag_off_table must contain an entry with loss_class=255; got {:?}", all_loss_classes);
}

/// A GBDT blob stored via `write_models_with_gbdt` must survive a round-trip:
/// `load_param` should return a `Param` whose `gbdt_peak_model` is `Some(…)`
/// and whose `predict_logit` produces the expected value.
#[test]
fn gbdt_blob_roundtrips_through_store() {
    use scoring_crate::gbdt_eval::{GbdtPeakModel, Tree};

    // Minimal one-split tree: feature[0] <= 0.5 → leaf -1.0; else → leaf +2.0
    // apply_sigmoid=true, isotonic=identity ([0,1]→[0,1]).
    // predict_logit([1.0]) → sigmoid(2.0) ≈ 0.8808 → logit ≈ 2.0
    let model = GbdtPeakModel {
        n_features: 1,
        apply_sigmoid: true,
        trees: vec![Tree {
            feature: vec![0, -1, -1],
            threshold: vec![0.5, 0.0, 0.0],
            left: vec![1, -1, -1],
            right: vec![2, -1, -1],
            value: vec![0.0, -1.0, 2.0],
            default_left: vec![1, 1, 1],
        }],
        iso_x: vec![0.0, 1.0],
        iso_y: vec![0.0, 1.0],
    };
    let blob = model.to_bytes();

    let tmp = tempfile::tempdir().unwrap();
    let store_path = tmp.path().join("with_gbdt.parquet");

    // Re-use param_with_loss_ions() as a structurally valid Param; the GBDT
    // blob is passed separately via write_models_with_gbdt.
    let param = param_with_loss_ions();

    write_models_with_gbdt(
        &store_path,
        &[("toy", &param, Some(blob.clone()))],
    )
    .expect("write_models_with_gbdt failed");

    let loaded = ModelStore::open(&store_path)
        .unwrap()
        .load_param("toy")
        .expect("load 'toy' model");

    let gm = loaded
        .gbdt_peak_model
        .expect("gbdt_peak_model must be Some after round-trip");

    let logit = gm.predict_logit(&[1.0]);
    assert!(
        (logit - 2.0).abs() < 1e-4,
        "predict_logit([1.0]) expected ≈ 2.0, got {logit}"
    );
}

/// Reading the existing bundled store (written without the loss_class column)
/// must yield all ion types with loss_class == 0.
#[test]
fn old_store_without_loss_class_reads_as_zero() {
    let bundled = Path::new(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../resources/models.parquet"
    ));
    let store = ModelStore::open(bundled).expect("open bundled models.parquet");
    let ids = store.model_ids();
    assert!(ids.len() >= 39, "expected >=39 bundled models, got {}", ids.len());

    for id in &ids {
        let param = store.load_param(id).expect("load bundled model");

        // All rank_dist ion types must have loss_class == 0
        for (part, ion_map) in &param.rank_dist_table {
            for ion in ion_map.keys() {
                assert_eq!(
                    ion.loss_class(), 0,
                    "bundled model {id} partition {:?}: ion {:?} has non-zero loss_class",
                    part, ion
                );
            }
        }

        // All frag_off ion types must have loss_class == 0
        for (part, frags) in &param.frag_off_table {
            for fof in frags {
                assert_eq!(
                    fof.ion_type.loss_class(), 0,
                    "bundled model {id} partition {:?}: frag_off ion {:?} has non-zero loss_class",
                    part, fof.ion_type
                );
            }
        }
    }
}

/// Multi-model write with per-model GBDT blobs must round-trip correctly:
/// - model A has a GBDT blob → loaded param has Some(gbdt_peak_model)
/// - model B has no blob → loaded param has None
/// - both models' source ledgers are preserved
/// - both models can be loaded after write
#[test]
fn multi_model_write_preserves_blobs_and_sources() {
    use scoring_crate::gbdt_eval::{GbdtPeakModel, Tree};

    // Build a minimal GbdtPeakModel and serialise it to bytes (same as A2.4 test).
    let model = GbdtPeakModel {
        n_features: 1,
        apply_sigmoid: true,
        trees: vec![Tree {
            feature: vec![0, -1, -1],
            threshold: vec![0.5, 0.0, 0.0],
            left: vec![1, -1, -1],
            right: vec![2, -1, -1],
            value: vec![0.0, -1.0, 2.0],
            default_left: vec![1, 1, 1],
        }],
        iso_x: vec![0.0, 1.0],
        iso_y: vec![0.0, 1.0],
    };
    let blob_a = model.to_bytes();

    let param_a = param_with_loss_ions(); // model A: will get GBDT blob
    let param_b = fixture();              // model B: no GBDT blob

    // Build a minimal SourceLedger + empty CountStats for each model.
    let make_ledger = |id: &str| SourceLedger {
        source_id: id.to_string(),
        dataset: format!("dataset_{id}"),
        n_psms: 42,
        date: "2026-01-01".to_string(),
        weight: 1.0,
        train_fdr: 0.01,
        instrument: "QExactive".to_string(),
        experiment_class: "standard".to_string(),
    };
    let sources_a: Vec<(SourceLedger, CountStats)> = vec![(make_ledger("src_a"), CountStats::new())];
    let sources_b: Vec<(SourceLedger, CountStats)> = vec![(make_ledger("src_b"), CountStats::new())];

    let tmp = tempfile::tempdir().unwrap();
    let store_path = tmp.path().join("multi_blob.parquet");

    // Write both models in one call via the new function.
    let models = &[
        ("model_a", &param_a, sources_a.as_slice()),
        ("model_b", &param_b, sources_b.as_slice()),
    ];
    let blobs: &[Option<Vec<u8>>] = &[
        Some(blob_a.clone()), // model A has a GBDT blob
        None,                 // model B has no blob
    ];
    write_all_models_with_sources_and_gbdt_pub(&store_path, models, blobs)
        .expect("write_all_models_with_sources_and_gbdt_pub failed");

    // --- Round-trip assertions ---
    let store = ModelStore::open(&store_path).unwrap();
    let ids = store.model_ids();
    assert_eq!(ids.len(), 2, "must have 2 models; got {:?}", ids);
    assert!(ids.contains(&"model_a".to_string()), "model_a missing");
    assert!(ids.contains(&"model_b".to_string()), "model_b missing");

    // model A: blob present + correct prediction
    let loaded_a = store.load_param("model_a").expect("load model_a");
    let gm_a = loaded_a.gbdt_peak_model
        .expect("model_a must have gbdt_peak_model after round-trip");
    let logit = gm_a.predict_logit(&[1.0]);
    assert!(
        (logit - 2.0).abs() < 1e-4,
        "model_a: predict_logit([1.0]) expected ≈ 2.0, got {logit}"
    );

    // model B: no blob
    let loaded_b = store.load_param("model_b").expect("load model_b");
    assert!(
        loaded_b.gbdt_peak_model.is_none(),
        "model_b must have no gbdt_peak_model"
    );

    // Both models' source ledgers survived.
    let ledgers_a = store.load_sources("model_a").expect("load_sources model_a");
    assert_eq!(ledgers_a.len(), 1, "model_a must have 1 source ledger");
    assert_eq!(ledgers_a[0].source_id, "src_a");

    let ledgers_b = store.load_sources("model_b").expect("load_sources model_b");
    assert_eq!(ledgers_b.len(), 1, "model_b must have 1 source ledger");
    assert_eq!(ledgers_b[0].source_id, "src_b");
}

/// A `frag_intensity_model_bytes` blob stored alongside a manifest row must
/// survive a write→read round-trip: `load_param` returns a `Param` whose
/// `frag_intensity_model` is `Some(…)` and whose `predict_value` gives the
/// expected raw tree sum.
#[test]
fn frag_intensity_model_roundtrips_through_store() {
    use std::sync::Arc;
    use scoring_crate::gbdt_eval::{GbdtPeakModel, Tree};

    // Minimal model: single leaf → constant 3.5.
    let model = GbdtPeakModel {
        n_features: 1,
        apply_sigmoid: false,
        trees: vec![Tree {
            feature: vec![-1],
            threshold: vec![0.0],
            left: vec![-1],
            right: vec![-1],
            value: vec![3.5],
            default_left: vec![1],
        }],
        iso_x: vec![],
        iso_y: vec![],
    };
    let blob = model.to_bytes();

    let tmp = tempfile::tempdir().unwrap();
    let store_path = tmp.path().join("frag_intensity.parquet");

    let mut param = param_with_loss_ions();
    // Attach the model to the Param before writing, so the writer can read it from the field.
    param.frag_intensity_model = Some(Arc::new(model));

    // The write path reads frag_intensity_model from the Param and persists the blob.
    model_train::store::write_models(&store_path, &[("toy".to_string(), &param)])
        .expect("write_models failed");

    let loaded = ModelStore::open(&store_path)
        .unwrap()
        .load_param("toy")
        .expect("load 'toy' model");

    let fim = loaded
        .frag_intensity_model
        .expect("frag_intensity_model must be Some after round-trip");

    let v = fim.predict_value(&[0.0]);
    assert!(
        (v - 3.5).abs() < 1e-5,
        "predict_value expected 3.5, got {v}"
    );

    // Verify the blob byte-matches.
    assert_eq!(fim.to_bytes(), blob, "round-tripped blob must be byte-identical");
}

/// A `rich_ion_model_bytes` blob stored alongside a manifest row must
/// survive a write→read round-trip: `load_param` returns a `Param` whose
/// `rich_ion_model` is `Some(…)` and whose `predict_value` gives the
/// expected raw tree sum.
#[test]
fn rich_ion_model_round_trips_through_store() {
    use std::sync::Arc;
    use scoring_crate::gbdt_eval::{GbdtPeakModel, Tree};

    // Minimal model: single leaf → constant 2.71.
    let model = GbdtPeakModel {
        n_features: 1,
        apply_sigmoid: false,
        trees: vec![Tree {
            feature: vec![-1],
            threshold: vec![0.0],
            left: vec![-1],
            right: vec![-1],
            value: vec![2.71],
            default_left: vec![1],
        }],
        iso_x: vec![],
        iso_y: vec![],
    };
    let blob = model.to_bytes();

    let tmp = tempfile::tempdir().unwrap();
    let store_path = tmp.path().join("rich_ion.parquet");

    let mut param = param_with_loss_ions();
    // Attach the model to the Param before writing, so the writer can read it from the field.
    param.rich_ion_model = Some(Arc::new(model));

    // The write path reads rich_ion_model from the Param and persists the blob.
    model_train::store::write_models(&store_path, &[("toy".to_string(), &param)])
        .expect("write_models failed");

    let loaded = ModelStore::open(&store_path)
        .unwrap()
        .load_param("toy")
        .expect("load 'toy' model");

    let rim = loaded
        .rich_ion_model
        .expect("rich_ion_model must be Some after round-trip");

    let v = rim.predict_value(&[0.0]);
    assert!(
        (v - 2.71).abs() < 1e-5,
        "predict_value expected 2.71, got {v}"
    );

    // Verify the blob byte-matches.
    assert_eq!(rim.to_bytes(), blob, "round-tripped blob must be byte-identical");
}

/// A model written WITHOUT `rich_ion_model` (None) must load with
/// `rich_ion_model == None` — back-compat with old stores.
#[test]
fn missing_rich_ion_model_column_loads_as_none() {
    let tmp = tempfile::tempdir().unwrap();
    let store_path = tmp.path().join("no_rich_ion.parquet");

    let param = param_with_loss_ions(); // rich_ion_model is None
    model_train::store::write_models(&store_path, &[("m".to_string(), &param)])
        .expect("write_models failed");

    let loaded = ModelStore::open(&store_path)
        .unwrap()
        .load_param("m")
        .expect("load 'm' model");

    assert!(
        loaded.rich_ion_model.is_none(),
        "rich_ion_model must be None when not stored"
    );
}

/// A model written WITHOUT `frag_intensity_model` (None) must load with
/// `frag_intensity_model == None` — back-compat with old stores.
#[test]
fn missing_frag_intensity_model_column_loads_as_none() {
    let tmp = tempfile::tempdir().unwrap();
    let store_path = tmp.path().join("no_frag_intensity.parquet");

    let param = param_with_loss_ions(); // frag_intensity_model is None
    model_train::store::write_models(&store_path, &[("m".to_string(), &param)])
        .expect("write_models failed");

    let loaded = ModelStore::open(&store_path)
        .unwrap()
        .load_param("m")
        .expect("load 'm' model");

    assert!(
        loaded.frag_intensity_model.is_none(),
        "frag_intensity_model must be None when not stored"
    );
}
