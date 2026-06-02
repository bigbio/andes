use model_train::counts::CountStats;
use model_train::store::{ModelStore, SourceLedger, write_model_with_sources};
use scoring_crate::param_model::{IonType, Param, Partition};
use std::path::Path;

fn fixture_param() -> Param {
    let param_path = Path::new(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/HCD_QExactive_Tryp.param"
    ));
    Param::load_from_file(param_path).expect("load fixture HCD_QExactive_Tryp.param")
}

fn test_partition_a() -> Partition {
    Partition { charge: 2, parent_mass: 800.0_f32, seg_num: 0 }
}

fn test_partition_b() -> Partition {
    Partition { charge: 3, parent_mass: 1200.0_f32, seg_num: 1 }
}

fn test_ion_prefix() -> IonType {
    IonType::Prefix { charge: 1, offset_bits: 1.007_f32.to_bits() }
}

fn test_ion_suffix() -> IonType {
    IonType::Suffix { charge: 1, offset_bits: 0.0_f32.to_bits() }
}

fn make_count_stats_0() -> CountStats {
    let pa = test_partition_a();
    let pb = test_partition_b();
    let ion_p = test_ion_prefix();
    let ion_s = test_ion_suffix();

    let mut s = CountStats::new();
    s.bump_rank(pa, ion_p, 0);
    s.bump_rank(pa, ion_p, 2);
    s.bump_rank(pa, ion_s, 1);
    s.bump_rank(pb, ion_p, 0);
    s.bump_error(pa, 5);
    s.bump_error(pa, 7);
    s.bump_error(pb, 3);
    s.bump_noise_error(pa, 2);
    s.bump_noise_error(pb, 4);
    s.bump_existence(pa, 0);
    s.bump_existence(pa, 0);
    s.bump_existence(pa, 1);
    s.bump_existence(pb, 2);
    s.bump_charge(2);
    s.bump_charge(2);
    s.bump_charge(3);
    s
}

fn make_count_stats_1() -> CountStats {
    let pa = test_partition_a();
    let pb = test_partition_b();
    let ion_p = test_ion_prefix();

    let mut s = CountStats::new();
    s.bump_rank(pa, ion_p, 1);
    s.bump_rank(pb, ion_p, 0);
    s.bump_rank(pb, ion_p, 0);
    s.bump_error(pa, 5);
    s.bump_noise_error(pb, 4);
    s.bump_existence(pa, 0);
    s.bump_existence(pb, 2);
    s.bump_charge(2);
    s.bump_charge(4);
    s
}

#[test]
fn per_source_stats_round_trip_and_sum() {
    let param = fixture_param();
    let s0 = make_count_stats_0();
    let s1 = make_count_stats_1();

    let ledger0 = SourceLedger {
        source_id: "s0".to_string(),
        dataset: "PXD001819".to_string(),
        n_psms: 14839,
        date: "2026-01-01".to_string(),
        weight: 1.0,
        train_fdr: 0.01,
        instrument: "QExactive".to_string(),
        experiment_class: "standard".to_string(),
    };
    let ledger1 = SourceLedger {
        source_id: "s1".to_string(),
        dataset: "PXD009630".to_string(),
        n_psms: 9788,
        date: "2026-02-01".to_string(),
        weight: 0.5,
        train_fdr: 0.01,
        instrument: "Fusion".to_string(),
        experiment_class: "tmt".to_string(),
    };

    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("models.parquet");

    write_model_with_sources(
        &path,
        "m",
        &param,
        &[(ledger0.clone(), s0.clone()), (ledger1.clone(), s1.clone())],
    )
    .unwrap();

    let store = ModelStore::open(&path).unwrap();

    // ledger round-trips
    let ledgers = store.load_sources("m").unwrap();
    assert_eq!(ledgers.len(), 2, "expected 2 source ledger entries");
    let found_s0 = ledgers.iter().any(|l| l.source_id == "s0");
    let found_s1 = ledgers.iter().any(|l| l.source_id == "s1");
    assert!(found_s0, "ledger s0 not found");
    assert!(found_s1, "ledger s1 not found");

    // Check ledger fields round-trip
    let l0 = ledgers.iter().find(|l| l.source_id == "s0").unwrap();
    assert_eq!(l0.dataset, "PXD001819");
    assert_eq!(l0.n_psms, 14839);
    assert_eq!(l0.date, "2026-01-01");
    assert!((l0.weight - 1.0_f32).abs() < 1e-6);
    assert!((l0.train_fdr - 0.01_f32).abs() < 1e-5);
    assert_eq!(l0.instrument, "QExactive");
    assert_eq!(l0.experiment_class, "standard");

    // each source's stats round-trip exactly
    let loaded_s0 = store.load_source_stats("m", "s0").unwrap();
    let loaded_s1 = store.load_source_stats("m", "s1").unwrap();
    assert_eq!(loaded_s0, s0, "CountStats for s0 did not round-trip exactly");
    assert_eq!(loaded_s1, s1, "CountStats for s1 did not round-trip exactly");

    // and load_param still works (search path unaffected)
    assert!(store.load_param("m").is_ok(), "load_param failed on store with source/stat rows");

    // sum of sources == an aggregate you can recompute
    let mut agg = s0.clone();
    agg.add(&s1);
    let mut from_store = store.load_source_stats("m", "s0").unwrap();
    from_store.add(&store.load_source_stats("m", "s1").unwrap());
    assert_eq!(from_store, agg, "sum of loaded sources must equal direct aggregate");
}

#[test]
fn load_param_ignores_source_stat_rows() {
    // Store with no sources must still work.
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("models.parquet");
    let param = fixture_param();

    write_model_with_sources(&path, "only_model", &param, &[]).unwrap();

    let store = ModelStore::open(&path).unwrap();
    let loaded = store.load_param("only_model").unwrap();
    assert_eq!(loaded, param, "load_param must still work with zero source rows");
    let sources = store.load_sources("only_model").unwrap();
    assert_eq!(sources.len(), 0);
}

#[test]
fn write_models_still_works_without_sources() {
    // The legacy write_models function must continue to work.
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("models.parquet");
    let param = fixture_param();

    model_train::store::write_models(&path, &[("m".to_string(), &param)]).unwrap();

    let store = ModelStore::open(&path).unwrap();
    assert!(store.load_param("m").is_ok());
    // No sources in a legacy store.
    let sources = store.load_sources("m").unwrap();
    assert_eq!(sources.len(), 0);
}
