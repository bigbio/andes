use std::path::PathBuf;

use model_train::store::{split_store_by_protocol, ModelStore};
use scoring_crate::param_model::Param;

fn bundled_single_file() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../resources/models/protocol=Automatic/models.parquet")
}

fn assert_stores_equivalent(single: &ModelStore, partitioned: &ModelStore) {
    let mut a = single.model_ids();
    let mut b = partitioned.model_ids();
    a.sort();
    b.sort();
    assert_eq!(
        a, b,
        "model_id sets differ between single file and partitioned dir"
    );

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
    assert_eq!(
        sa, sb,
        "selection_entries differ between single file and partitioned dir"
    );

    for id in &a {
        let pa: Param = single.load_param(id).expect("load_param single");
        let pb: Param = partitioned.load_param(id).expect("load_param partitioned");
        assert_params_equal(id, &pa, &pb);
    }

    for id in &a {
        let mut la = single.load_sources(id).expect("load_sources single");
        let mut lb = partitioned
            .load_sources(id)
            .expect("load_sources partitioned");
        la.sort_by(|x, y| x.source_id.cmp(&y.source_id));
        lb.sort_by(|x, y| x.source_id.cmp(&y.source_id));
        assert_eq!(
            la, lb,
            "[{id}] load_sources differ between single and partitioned"
        );

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

fn assert_params_equal(id: &str, a: &Param, b: &Param) {
    assert_eq!(a.version, b.version, "[{id}] version");
    assert_eq!(a.data_type, b.data_type, "[{id}] data_type");
    assert_eq!(a.mme, b.mme, "[{id}] mme");
    assert_eq!(
        a.apply_deconvolution, b.apply_deconvolution,
        "[{id}] apply_deconvolution"
    );
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
    assert_eq!(
        a.num_precursor_off, b.num_precursor_off,
        "[{id}] num_precursor_off"
    );
    assert_eq!(
        a.precursor_off_map, b.precursor_off_map,
        "[{id}] precursor_off_map"
    );
    assert_eq!(a.frag_off_table, b.frag_off_table, "[{id}] frag_off_table");
    assert_eq!(a.max_rank, b.max_rank, "[{id}] max_rank");
    assert_eq!(
        a.rank_dist_table, b.rank_dist_table,
        "[{id}] rank_dist_table"
    );
    assert_eq!(
        a.error_scaling_factor, b.error_scaling_factor,
        "[{id}] error_scaling_factor"
    );
    assert_eq!(
        a.ion_err_dist_table, b.ion_err_dist_table,
        "[{id}] ion_err_dist_table"
    );
    assert_eq!(
        a.noise_err_dist_table, b.noise_err_dist_table,
        "[{id}] noise_err_dist_table"
    );
    assert_eq!(
        a.ion_existence_table, b.ion_existence_table,
        "[{id}] ion_existence_table"
    );

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
    let single =
        ModelStore::open(&single_path).expect("open bundled partition as single-file store");
    assert_eq!(
        single.model_ids().len(),
        11,
        "Automatic partition ships 11 models"
    );

    let tmp = tempfile::tempdir().expect("tempdir");
    let dir = tmp.path().join("models");
    let written = split_store_by_protocol(&single_path, &dir).expect("split store");
    assert!(!written.is_empty(), "split produced no partitions");

    let partitioned = ModelStore::open(&dir).expect("open partitioned dir");
    assert_stores_equivalent(&single, &partitioned);
}

#[test]
fn shipped_partitions_load_identical_to_v1_single_file() {
    let Ok(src) = std::env::var("ANDES_V1_STORE") else {
        eprintln!("ANDES_V1_STORE not set — skipping shipped-partition parity test");
        return;
    };
    let shipped = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../resources/models");
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

    let tmp = tempfile::tempdir().expect("tempdir");
    let out: PathBuf = tmp.path().join("models");
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
