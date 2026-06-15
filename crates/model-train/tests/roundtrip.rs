use model_train::store::{write_models, ModelStore};
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

/// Reading the existing bundled store (written without the loss_class column)
/// must yield all ion types with loss_class == 0.
#[test]
fn old_store_without_loss_class_reads_as_zero() {
    let bundled = Path::new(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../resources/ionstat/models.parquet"
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
