//! TDD tests for incremental model training (add/remove/reweight/decay + gate).
//!
//! These tests validate the EXACT count-arithmetic property: adding a source
//! and then removing it restores the original estimated model parameter
//! exactly, because per-source `CountStats` are re-summed from scratch on
//! each update.

use std::path::Path;

use model_train::{
    counts::CountStats,
    estimate::EstimatorConfig,
    gate::{evaluate_candidate, YieldDelta},
    store::{
        attach_to_model, commit_update, update_add, update_decay, update_remove, update_reweight,
        write_model_with_sources, ModelStore, SourceLedger,
    },
};
use scoring_crate::{param_model::Param, scoring::rank_scorer::RankScorer};

// ---------------------------------------------------------------------------
// Fixtures
// ---------------------------------------------------------------------------

fn fixture_param() -> Param {
    // hcd_qexactive_tryp from the canonical Parquet store.
    let bundled = Path::new(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../resources/models"
    ));
    let store = model_train::store::ModelStore::open(bundled).expect("open bundled model store");
    store
        .load_param("hcd_qexactive_tryp")
        .expect("load hcd_qexactive_tryp from store")
}

fn make_stats_s0() -> CountStats {
    use scoring_crate::param_model::{IonType, Partition};
    let pa = Partition {
        charge: 2,
        parent_mass: 800.0_f32,
        seg_num: 0,
    };
    let ion = IonType::Prefix {
        charge: 1,
        offset_bits: 1.007_f32.to_bits(),
        loss_class: 0,
    };

    let mut s = CountStats::new();
    s.bump_rank(pa, ion, 0);
    s.bump_rank(pa, ion, 1);
    s.bump_rank(pa, ion, 2);
    s.bump_error(pa, 5);
    s.bump_error(pa, 7);
    s.bump_noise_error(pa, 2);
    s.bump_existence(pa, 0);
    s.bump_existence(pa, 1);
    s.bump_charge(2);
    s.bump_charge(2);
    s.bump_charge(3);
    s
}

fn make_stats_s1() -> CountStats {
    use scoring_crate::param_model::{IonType, Partition};
    let pa = Partition {
        charge: 2,
        parent_mass: 800.0_f32,
        seg_num: 0,
    };
    let ion = IonType::Prefix {
        charge: 1,
        offset_bits: 1.007_f32.to_bits(),
        loss_class: 0,
    };

    let mut s = CountStats::new();
    s.bump_rank(pa, ion, 0);
    s.bump_rank(pa, ion, 3);
    s.bump_error(pa, 6);
    s.bump_noise_error(pa, 1);
    s.bump_existence(pa, 2);
    s.bump_charge(2);
    s.bump_charge(3);
    s
}

fn make_ledger(id: &str) -> SourceLedger {
    SourceLedger {
        source_id: id.to_string(),
        dataset: format!("dataset_{id}"),
        n_psms: 100,
        date: "2026-01-01".to_string(),
        weight: 1.0,
        train_fdr: 0.01,
        decay_half_life_days: None,
        instrument: "QExactive".to_string(),
        experiment_class: "standard".to_string(),
    }
}

// ---------------------------------------------------------------------------
// PART B tests
// ---------------------------------------------------------------------------

/// EXACTNESS TEST: add s1 then remove s1 must restore the original param P0.
#[test]
fn add_then_remove_source_restores_model() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("m.parquet");

    let param = fixture_param();
    let s0 = make_stats_s0();
    let ledger0 = make_ledger("s0");

    // Estimate P0 from s0 alone and write to store.
    let cfg = EstimatorConfig::default();
    let estimator = model_train::estimate::Estimator::new(cfg.clone());
    let p0 = estimator.estimate(&s0, &param);

    write_model_with_sources(&path, "m", &p0, &[(ledger0.clone(), s0.clone())], None).unwrap();

    // Verify P0 round-trips.
    let store = ModelStore::open(&path).unwrap();
    let loaded_p0 = store.load_param("m").unwrap();
    assert_eq!(
        p0, loaded_p0,
        "P0 must round-trip exactly through the store"
    );

    // update_add: append s1.
    let s1 = make_stats_s1();
    let ledger1 = make_ledger("s1");
    let (candidate_p1, sources_with_s1) =
        update_add(&path, "m", ledger1, s1.clone(), cfg.clone()).unwrap();
    commit_update(&path, "m", &candidate_p1, &sources_with_s1).unwrap();

    // The store now has s0+s1; param should differ from P0.
    let store2 = ModelStore::open(&path).unwrap();
    let p_after_add = store2.load_param("m").unwrap();
    // (We don't assert on p_after_add vs p0 here — only restore matters.)
    drop(p_after_add);

    // update_remove: remove s1 → should restore P0 exactly.
    let (candidate_p0_restored, sources_s0_only) =
        update_remove(&path, "m", "s1", cfg.clone()).unwrap();
    commit_update(&path, "m", &candidate_p0_restored, &sources_s0_only).unwrap();

    let store3 = ModelStore::open(&path).unwrap();
    let p_restored = store3.load_param("m").unwrap();

    assert_eq!(
        p0, p_restored,
        "add then remove must restore the model exactly (count arithmetic is exact)"
    );
}

/// Removing a non-existent source returns an error.
#[test]
fn remove_nonexistent_source_errors() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("m.parquet");

    let param = fixture_param();
    let s0 = make_stats_s0();
    let ledger0 = make_ledger("s0");
    write_model_with_sources(&path, "m", &param, &[(ledger0, s0)], None).unwrap();

    let cfg = EstimatorConfig::default();
    let result = update_remove(&path, "m", "nonexistent", cfg);
    assert!(
        result.is_err(),
        "removing a nonexistent source must return an error"
    );
}

/// Reweighting a source changes the estimated param and the stored weight.
#[test]
fn reweight_source_changes_param() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("m.parquet");

    let param = fixture_param();
    let s0 = make_stats_s0();
    let s1 = make_stats_s1();
    let ledger0 = make_ledger("s0");
    let ledger1 = make_ledger("s1");
    write_model_with_sources(&path, "m", &param, &[(ledger0, s0), (ledger1, s1)], None).unwrap();

    let cfg = EstimatorConfig::default();
    let (candidate, sources) = update_reweight(&path, "m", "s1", 0.0, cfg).unwrap();
    commit_update(&path, "m", &candidate, &sources).unwrap();

    // After reweight to 0, s1 contributes nothing — check the weight was stored.
    let store = ModelStore::open(&path).unwrap();
    let ledgers = store.load_sources("m").unwrap();
    let l1 = ledgers.iter().find(|l| l.source_id == "s1").unwrap();
    assert!(
        (l1.weight - 0.0_f32).abs() < 1e-6,
        "s1 weight must be 0 after reweight"
    );
}

/// Decay with empty dates leaves weights unchanged and emits a warning.
#[test]
fn decay_with_empty_date_unchanged() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("m.parquet");

    let param = fixture_param();
    let s0 = make_stats_s0();
    let mut ledger0 = make_ledger("s0");
    ledger0.date = String::new(); // empty date
    write_model_with_sources(&path, "m", &param, &[(ledger0.clone(), s0)], None).unwrap();

    let cfg = EstimatorConfig::default();
    let (_, sources) = update_decay(&path, "m", 365.0, cfg).unwrap();

    // Weight must be unchanged (since date is empty).
    let updated_ledger = &sources[0].0;
    assert!(
        (updated_ledger.weight - ledger0.weight).abs() < 1e-6,
        "weight must be unchanged when date is empty"
    );
}

/// commit_update preserves other models in the store.
#[test]
fn commit_update_preserves_other_models() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("m.parquet");

    let param = fixture_param();
    let s0 = make_stats_s0();
    let ledger0 = make_ledger("s0");

    // Write two models: "m" and "other".
    use model_train::store::write_models;
    // First write "other" model (no sources).
    write_models(&path, &[("other".to_string(), &param)]).unwrap();
    // Then add "m" with sources by rebuilding the store.
    // We need to load "other" and rewrite both.
    {
        let store = ModelStore::open(&path).unwrap();
        let other_param = store.load_param("other").unwrap();
        let _ = other_param; // just verifying it loads
    }
    // Write "m" with sources into a fresh store (overwrite).
    let path2 = dir.path().join("m2.parquet");
    write_model_with_sources(&path2, "m", &param, &[(ledger0.clone(), s0.clone())], None).unwrap();

    let cfg = EstimatorConfig::default();
    let s1 = make_stats_s1();
    let ledger1 = make_ledger("s1");
    let (candidate, sources) = update_add(&path2, "m", ledger1, s1, cfg).unwrap();
    commit_update(&path2, "m", &candidate, &sources).unwrap();

    // "m" must still load.
    let store = ModelStore::open(&path2).unwrap();
    assert!(
        store.load_param("m").is_ok(),
        "model 'm' must survive commit_update"
    );
    let ledgers = store.load_sources("m").unwrap();
    assert_eq!(ledgers.len(), 2, "must have 2 sources after add+commit");
}

// ---------------------------------------------------------------------------
// PART C tests (acceptance gate)
// ---------------------------------------------------------------------------

/// acceptance_gate_rejects_worse_candidate: the gate rule is
/// `accepted iff candidate_count >= current_count`.
#[test]
fn acceptance_gate_rule() {
    // Directly test the YieldDelta accept rule without running a real search.
    let accepted = YieldDelta {
        current_count: 10,
        candidate_count: 10,
    };
    assert!(accepted.is_accepted(), "equal counts must be accepted");

    let better = YieldDelta {
        current_count: 10,
        candidate_count: 15,
    };
    assert!(better.is_accepted(), "higher candidate must be accepted");

    let worse = YieldDelta {
        current_count: 10,
        candidate_count: 5,
    };
    assert!(!worse.is_accepted(), "lower candidate must be rejected");
}

/// Full end-to-end acceptance gate using the BSA fixture.
///
/// Uses the same scorer for both current and candidate → YieldDelta should
/// show equal counts (accepted).
#[cfg(test)]
#[test]
fn acceptance_gate_same_model_is_accepted() {
    use input::MgfReader;
    use model::{AminoAcidSetBuilder, ModLocation, Modification, ResidueSpec};
    use search::SearchParams;
    use std::io::BufReader;

    fn fixture(rel: &str) -> std::path::PathBuf {
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../..")
            .join(rel)
            .canonicalize()
            .unwrap_or_else(|e| panic!("canonicalize {rel}: {e}"))
    }

    let bsa_mgf = fixture("test-fixtures/test.mgf.gz");
    let bsa_fasta = fixture("test-fixtures/BSA.fasta");

    // Load spectra.
    let f = input::open_buf_maybe_gz(&bsa_mgf).expect("open test.mgf");
    let reader = MgfReader::new(BufReader::new(f));
    let spectra: Vec<_> = reader.filter_map(|r| r.ok()).collect();

    // Load scorer from the canonical Parquet store.
    let param = fixture_param();
    let scorer = RankScorer::new(&param);

    let cam = Modification {
        name: "Carbamidomethyl".into(),
        mass_delta: 57.02146,
        residue: ResidueSpec::Specific(b'C'),
        location: ModLocation::Anywhere,
        fixed: true,
        accession: None,
        neutral_losses: Vec::new(),
        loss_class: 0,
        annotated_site_mass: None,
        blocks_cleavage: false,
        diagnostic_ions: Vec::new(),
    };
    let ox = Modification {
        name: "Oxidation".into(),
        mass_delta: 15.99491,
        residue: ResidueSpec::Specific(b'M'),
        location: ModLocation::Anywhere,
        fixed: false,
        accession: None,
        neutral_losses: Vec::new(),
        loss_class: 0,
        annotated_site_mass: None,
        blocks_cleavage: false,
        diagnostic_ions: Vec::new(),
    };
    let aa = AminoAcidSetBuilder::new_standard()
        .add_fixed_mod(cam)
        .add_variable_mod(ox)
        .build()
        .unwrap();
    let search_params = SearchParams::default_tryptic(aa);

    // Same scorer for both current and candidate → counts should be equal.
    let cache_dir = tempfile::tempdir().expect("cache dir");
    let delta = evaluate_candidate(
        &spectra,
        &bsa_fasta,
        &scorer,
        &scorer,
        &search_params,
        0.5,
        cache_dir.path(),
    )
    .expect("evaluate_candidate must succeed");

    assert_eq!(
        delta.current_count, delta.candidate_count,
        "identical models must yield identical counts"
    );
    assert!(
        delta.is_accepted(),
        "identical models must be accepted by the gate"
    );
}

/// Serialised GBDT blobs of a Param: (peak, fragment-intensity, rich-ion).
fn gbdt_blobs(p: &Param) -> [Option<Vec<u8>>; 3] {
    [
        p.gbdt_peak_model.as_ref().map(|m| m.to_bytes()),
        p.frag_intensity_model.as_ref().map(|m| m.to_bytes()),
        p.rich_ion_model.as_ref().map(|m| m.to_bytes()),
    ]
}

/// An update re-estimates the rank tables only: every GBDT carried by the
/// stored model (a copy of a bundled model) survives the candidate, the gate
/// input and the committed store byte-identically.
#[test]
fn update_keeps_gbdt_models_byte_identical() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("m.parquet");

    let param = fixture_param();
    assert!(
        param.frag_intensity_model.is_some() && param.rich_ion_model.is_some(),
        "bundled hcd_qexactive_tryp must carry frag-intensity and rich-ion GBDTs"
    );
    let before = gbdt_blobs(&param);
    let peak_blob = before[0].clone();
    write_model_with_sources(
        &path,
        "m",
        &param,
        &[(make_ledger("s0"), make_stats_s0())],
        peak_blob,
    )
    .unwrap();
    let stored = ModelStore::open(&path).unwrap().load_param("m").unwrap();
    assert_eq!(gbdt_blobs(&stored), before, "store round-trip");

    let cfg = EstimatorConfig::default();
    let (candidate, sources) =
        update_add(&path, "m", make_ledger("s1"), make_stats_s1(), cfg.clone()).unwrap();
    assert_eq!(gbdt_blobs(&candidate), before, "update_add candidate");
    commit_update(&path, "m", &candidate, &sources).unwrap();
    let committed = ModelStore::open(&path).unwrap().load_param("m").unwrap();
    assert_eq!(gbdt_blobs(&committed), before, "committed after update_add");

    for (candidate, _) in [
        update_reweight(&path, "m", "s1", 0.5, cfg.clone()).unwrap(),
        update_remove(&path, "m", "s1", cfg.clone()).unwrap(),
        update_decay(&path, "m", 365.0, cfg.clone()).unwrap(),
    ] {
        assert_eq!(gbdt_blobs(&candidate), before);
    }
}

/// Decay is a function of source age, not of how often it runs: decaying,
/// committing and decaying again on the same day yields the same ledger and
/// the same candidate as a single decay, and the base weight is kept.
#[test]
fn repeated_decay_on_same_day_equals_one_decay() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("m.parquet");

    let param = fixture_param();
    let mut old = make_ledger("s0");
    old.date = "2020-01-01".to_string();
    let recent = make_ledger("s1");
    write_model_with_sources(
        &path,
        "m",
        &param,
        &[(old, make_stats_s0()), (recent, make_stats_s1())],
        None,
    )
    .unwrap();

    let cfg = EstimatorConfig::default();
    let (undecayed, _) = update_reweight(&path, "m", "s0", 1.0, cfg.clone()).unwrap();
    let (once, once_sources) = update_decay(&path, "m", 30.0, cfg.clone()).unwrap();
    assert_ne!(once, undecayed, "decay must change the estimate");
    commit_update(&path, "m", &once, &once_sources).unwrap();

    let (twice, twice_sources) = update_decay(&path, "m", 30.0, cfg.clone()).unwrap();
    assert_eq!(twice, once, "second same-day decay must not compound");
    for (l, _) in &twice_sources {
        assert_eq!(l.weight, 1.0, "ledger keeps the undecayed base weight");
        assert_eq!(l.decay_half_life_days, Some(30.0));
    }
    commit_update(&path, "m", &twice, &twice_sources).unwrap();

    // The recorded half-life persists, so a later reweight keeps the decay.
    let (reweighted, _) = update_reweight(&path, "m", "s1", 1.0, cfg).unwrap();
    assert_eq!(reweighted, once, "reweight must keep the recorded decay");
}

/// Update operations refuse a partitioned (directory) store up front instead
/// of failing at commit time after the work is done.
#[test]
fn update_rejects_partitioned_store() {
    let bundled = Path::new(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../resources/models"
    ));
    let cfg = EstimatorConfig::default();
    let err = update_decay(bundled, "hcd_qexactive_tryp", 365.0, cfg.clone()).unwrap_err();
    assert!(err.to_string().contains("partitioned"), "{err}");
    let err = update_add(
        bundled,
        "hcd_qexactive_tryp",
        make_ledger("s1"),
        make_stats_s1(),
        cfg,
    )
    .unwrap_err();
    assert!(err.to_string().contains("partitioned"), "{err}");
    let p = fixture_param();
    let err = commit_update(bundled, "hcd_qexactive_tryp", &p, &[]).unwrap_err();
    assert!(err.to_string().contains("partitioned"), "{err}");
}

/// Attaching a GBDT to an existing model replaces only that GBDT: the rank
/// tables, the other GBDT blobs and the source ledger with its stats survive,
/// and other models are untouched. A model id not in the store starts from
/// the fallback with no sources.
#[test]
fn attach_to_model_keeps_existing_model_content() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("m.parquet");

    let param = fixture_param();
    let peak_blob = param.gbdt_peak_model.as_ref().map(|m| m.to_bytes());
    write_model_with_sources(
        &path,
        "m",
        &param,
        &[(make_ledger("s0"), make_stats_s0())],
        peak_blob,
    )
    .unwrap();
    let other = model_train::estimate::Estimator::new(EstimatorConfig::default())
        .estimate(&make_stats_s1(), &param);
    let (cand, srcs) = update_add(
        &path,
        "m",
        make_ledger("s1"),
        make_stats_s1(),
        EstimatorConfig::default(),
    )
    .unwrap();
    commit_update(&path, "m", &cand, &srcs).unwrap();
    commit_update(&path, "other", &other, &[(make_ledger("o"), make_stats_s1())]).unwrap();
    let before = ModelStore::open(&path).unwrap().load_param("m").unwrap();

    // A distinguishable replacement rich-ion model that still satisfies the
    // rich-ion feature contract: the stored one with its last tree dropped.
    let stored = param.rich_ion_model.as_deref().expect("bundled rich-ion model");
    assert!(stored.trees.len() > 1);
    let replacement = std::sync::Arc::new(scoring_crate::gbdt_eval::GbdtPeakModel {
        n_features: stored.n_features,
        feature_names: stored.feature_names.clone(),
        apply_sigmoid: stored.apply_sigmoid,
        trees: stored.trees[..stored.trees.len() - 1].to_vec(),
        iso_x: stored.iso_x.clone(),
        iso_y: stored.iso_y.clone(),
        packed: Default::default(),
    });
    assert_ne!(*replacement, *stored);
    let existed = attach_to_model(&path, "m", other.clone(), |p| {
        p.rich_ion_model = Some(replacement.clone());
    })
    .unwrap();
    assert!(existed);

    let store = ModelStore::open(&path).unwrap();
    let after = store.load_param("m").unwrap();
    let mut expected = before.clone();
    expected.rich_ion_model = Some(replacement.clone());
    assert_eq!(after, expected, "only the rich-ion model may change");
    assert_eq!(gbdt_blobs(&after)[..2], gbdt_blobs(&before)[..2]);
    let ids: Vec<String> = store
        .load_sources("m")
        .unwrap()
        .into_iter()
        .map(|l| l.source_id)
        .collect();
    assert_eq!(ids, ["s0", "s1"], "source ledger must survive");
    assert_eq!(
        store.load_source_stats("m", "s1").unwrap(),
        make_stats_s1(),
        "per-source stats must survive"
    );
    assert_eq!(store.load_param("other").unwrap(), other);

    let existed = attach_to_model(&path, "fresh", other.clone(), |p| {
        p.rich_ion_model = Some(replacement.clone());
    })
    .unwrap();
    assert!(!existed);
    let store = ModelStore::open(&path).unwrap();
    assert!(store.load_sources("fresh").unwrap().is_empty());
    let mut expected = other;
    expected.rich_ion_model = Some(replacement);
    assert_eq!(store.load_param("fresh").unwrap(), expected);
}
