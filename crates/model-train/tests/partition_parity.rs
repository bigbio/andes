//! Parity gate for protocol-partitioned model stores.
//!
//! A partitioned store (a directory of `protocol=<P>/models.parquet` files)
//! MUST load identically to the single-file store it was split from:
//! - same set of `model_id`s,
//! - same `selection_entries()`,
//! - per-model `load_param(id)` semantically identical (every field, every
//!   table, every blob).
//!
//! The bundled-store test partitions `resources/models.parquet` into a temp
//! dir and proves equivalence. A second test does the same for an external
//! v1 store if `ANDES_V1_STORE` points at one (used to build the shipped
//! `resources/models/` partition set).

use std::path::PathBuf;

use model_train::store::{split_store_by_protocol, ModelStore};
use scoring_crate::param_model::Param;

fn bundled_single_file() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../resources/models.parquet")
}

/// Assert two stores are equivalent: same model_ids, same selection_entries,
/// and per-model Param equality.
fn assert_stores_equivalent(single: &ModelStore, partitioned: &ModelStore) {
    // 1. Same model_id set.
    let mut a = single.model_ids();
    let mut b = partitioned.model_ids();
    a.sort();
    b.sort();
    assert_eq!(a, b, "model_id sets differ between single file and partitioned dir");

    // 2. Same selection_entries (as a sorted set).
    let mut sa: Vec<_> = single
        .selection_entries()
        .into_iter()
        .map(|e| {
            (
                e.model_id,
                e.activation,
                e.instrument,
                e.enzyme,
                e.experiment_class.into_iter().collect::<Vec<_>>().join("+"),
            )
        })
        .collect();
    let mut sb: Vec<_> = partitioned
        .selection_entries()
        .into_iter()
        .map(|e| {
            (
                e.model_id,
                e.activation,
                e.instrument,
                e.enzyme,
                e.experiment_class.into_iter().collect::<Vec<_>>().join("+"),
            )
        })
        .collect();
    sa.sort();
    sb.sort();
    assert_eq!(sa, sb, "selection_entries differ between single file and partitioned dir");

    // 3. Per-model Param equality.
    for id in &a {
        let pa: Param = single.load_param(id).expect("load_param single");
        let pb: Param = partitioned.load_param(id).expect("load_param partitioned");
        assert_params_equal(id, &pa, &pb);
    }

    // 4. Per-model source ledgers + per-source stats equivalence.
    for id in &a {
        let mut la = single.load_sources(id).expect("load_sources single");
        let mut lb = partitioned.load_sources(id).expect("load_sources partitioned");
        la.sort_by(|x, y| x.source_id.cmp(&y.source_id));
        lb.sort_by(|x, y| x.source_id.cmp(&y.source_id));
        assert_eq!(la, lb, "[{id}] load_sources differ between single and partitioned");

        // Each (model_id, source_id) that has stat rows must reconstruct
        // identical sufficient statistics from either layout (CountStats'
        // PartialEq ignores trailing histogram zeros). Sources without stat
        // rows must be absent (NoModel) in BOTH layouts.
        for src in &la {
            match (
                single.load_source_stats(id, &src.source_id),
                partitioned.load_source_stats(id, &src.source_id),
            ) {
                (Ok(sa), Ok(sb)) => assert_eq!(
                    sa, sb,
                    "[{id}/{}] load_source_stats differ between single and partitioned",
                    src.source_id
                ),
                (Err(_), Err(_)) => {}
                (a, b) => panic!(
                    "[{id}/{}] load_source_stats presence differs: single={:?} partitioned={:?}",
                    src.source_id,
                    a.is_ok(),
                    b.is_ok()
                ),
            }
        }
    }
}

/// Compare two reconstructed Params field-by-field (the cache field and the
/// Arc-wrapped models are compared on content, not pointer identity).
fn assert_params_equal(id: &str, a: &Param, b: &Param) {
    assert_eq!(a.version, b.version, "[{id}] version");
    assert_eq!(a.data_type, b.data_type, "[{id}] data_type");
    assert_eq!(a.mme, b.mme, "[{id}] mme");
    assert_eq!(a.apply_deconvolution, b.apply_deconvolution, "[{id}] apply_deconvolution");
    assert_eq!(
        a.deconvolution_error_tolerance.to_bits(),
        b.deconvolution_error_tolerance.to_bits(),
        "[{id}] deconvolution_error_tolerance"
    );
    assert_eq!(a.charge_hist, b.charge_hist, "[{id}] charge_hist");
    assert_eq!(a.min_charge, b.min_charge, "[{id}] min_charge");
    assert_eq!(a.max_charge, b.max_charge, "[{id}] max_charge");
    assert_eq!(a.num_segments, b.num_segments, "[{id}] num_segments");
    assert_eq!(a.partitions, b.partitions, "[{id}] partitions");
    assert_eq!(a.num_precursor_off, b.num_precursor_off, "[{id}] num_precursor_off");
    assert_eq!(a.precursor_off_map, b.precursor_off_map, "[{id}] precursor_off_map");
    assert_eq!(a.frag_off_table, b.frag_off_table, "[{id}] frag_off_table");
    assert_eq!(a.max_rank, b.max_rank, "[{id}] max_rank");
    assert_eq!(a.rank_dist_table, b.rank_dist_table, "[{id}] rank_dist_table");
    assert_eq!(a.error_scaling_factor, b.error_scaling_factor, "[{id}] error_scaling_factor");
    assert_eq!(a.ion_err_dist_table, b.ion_err_dist_table, "[{id}] ion_err_dist_table");
    assert_eq!(a.noise_err_dist_table, b.noise_err_dist_table, "[{id}] noise_err_dist_table");
    assert_eq!(a.ion_existence_table, b.ion_existence_table, "[{id}] ion_existence_table");

    // Model blobs: compare presence + serialized bytes.
    assert_eq!(
        a.frag_intensity_model.as_ref().map(|m| m.to_bytes()),
        b.frag_intensity_model.as_ref().map(|m| m.to_bytes()),
        "[{id}] frag_intensity_model"
    );
    assert_eq!(
        a.rich_ion_model.as_ref().map(|m| m.to_bytes()),
        b.rich_ion_model.as_ref().map(|m| m.to_bytes()),
        "[{id}] rich_ion_model"
    );
    assert_eq!(
        a.gbdt_peak_model.as_ref().map(|m| m.to_bytes()),
        b.gbdt_peak_model.as_ref().map(|m| m.to_bytes()),
        "[{id}] gbdt_peak_model"
    );
}

#[test]
fn partitioned_store_loads_identical_to_bundled_single_file() {
    let single_path = bundled_single_file();
    let single = ModelStore::open(&single_path).expect("open bundled single-file store");

    let tmp = tempfile::tempdir().expect("tempdir");
    let dir = tmp.path().join("models");
    let written = split_store_by_protocol(&single_path, &dir).expect("split store");
    assert!(!written.is_empty(), "split produced no partitions");

    let partitioned = ModelStore::open(&dir).expect("open partitioned dir");
    assert_stores_equivalent(&single, &partitioned);
}

/// Build + verify the partitioned layout for an external v1 store.
///
/// Set `ANDES_V1_STORE=/path/to/models_v1_own.parquet` and optionally
/// `ANDES_V1_OUT=/path/to/resources/models` to also write the partitions
/// to a fixed location (used to produce the shipped `resources/models/`).
/// Verify the *shipped* pyarrow-written partitions (`resources/models/`) load
/// identically to the v1 single-file store. Unlike the test below, this does
/// NOT run the Rust splitter — it reads the partitions exactly as they ship,
/// which are produced by `scripts/split_store_by_protocol.py` (pyarrow/zstd).
///
/// Set `ANDES_V1_STORE=/path/to/models_v1_own.parquet` to enable. Skipped
/// otherwise (the v1 single file is not committed to the repo).
#[test]
fn shipped_partitions_load_identical_to_v1_single_file() {
    let Ok(src) = std::env::var("ANDES_V1_STORE") else {
        eprintln!("ANDES_V1_STORE not set — skipping shipped-partition parity test");
        return;
    };
    let shipped =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../resources/models");
    if !shipped.is_dir() {
        eprintln!("resources/models/ not present — skipping shipped-partition parity test");
        return;
    }
    let single = ModelStore::open(&PathBuf::from(src)).expect("open v1 single-file store");
    let partitioned = ModelStore::open(&shipped).expect("open shipped resources/models");
    assert_stores_equivalent(&single, &partitioned);
    eprintln!(
        "SHIPPED PARTITION PARITY: PASS ({} models, {})",
        single.model_ids().len(),
        shipped.display()
    );
}

#[test]
fn partitioned_v1_store_loads_identical_when_present() {
    let Ok(src) = std::env::var("ANDES_V1_STORE") else {
        eprintln!("ANDES_V1_STORE not set — skipping v1 partition parity test");
        return;
    };
    let src = PathBuf::from(src);
    let single = ModelStore::open(&src).expect("open v1 single-file store");

    let out: PathBuf = match std::env::var("ANDES_V1_OUT") {
        Ok(p) => PathBuf::from(p),
        Err(_) => tempfile::tempdir().unwrap().keep().join("models"),
    };
    if out.exists() {
        std::fs::remove_dir_all(&out)
            .unwrap_or_else(|e| panic!("failed to clear {} before re-split: {e}", out.display()));
    }
    let written = split_store_by_protocol(&src, &out).expect("split v1 store");
    eprintln!("v1 partitions written to {}:", out.display());
    for (proto, path) in &written {
        let sz = std::fs::metadata(path).map(|m| m.len()).unwrap_or(0);
        eprintln!("  protocol={proto:<16} {sz:>12} bytes  {}", path.display());
    }

    let partitioned = ModelStore::open(&out).expect("open v1 partitioned dir");
    assert_stores_equivalent(&single, &partitioned);
    eprintln!("v1 PARITY: PASS ({} models)", single.model_ids().len());
}
