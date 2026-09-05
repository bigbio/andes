//! Integration test for `andes train-intensity-gbdt`.
//!
//! Verifies the round-trip: build a tiny synthetic flat parquet → run
//! `train-intensity-gbdt` → reload the produced store → assert
//! `param.frag_intensity_model.is_some()` and `predict_value` is finite.

use std::path::Path;
use std::process::Command;
use std::sync::Arc;

use arrow::array::{
    ArrayRef, Float32Builder, Float64Array, Int32Array, Int32Builder, ListBuilder, StringArray,
};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;

use model_train::ModelStore;

/// One synthetic PSM row (same schema as the train-from-msnet parquet).
struct Row {
    seq: &'static str,
    charge: i32,
    prec_mz: f64,
    res_mod_pos: Vec<i32>,
    res_mod_delta: Vec<f64>,
    nterm_delta: f64,
    cterm_delta: f64,
    peaks: Vec<(f32, f32)>,
}

/// Write the flat training parquet at `path` from `rows`.
/// Schema mirrors `read_msnet_parquet` in andes.rs.
fn write_flat_parquet(path: &Path, rows: &[Row]) {
    let schema = Arc::new(Schema::new(vec![
        Field::new("seq", DataType::Utf8, false),
        Field::new("charge", DataType::Int32, false),
        Field::new("prec_mz", DataType::Float64, false),
        Field::new(
            "res_mod_pos",
            DataType::List(Arc::new(Field::new("item", DataType::Int32, true))),
            true,
        ),
        Field::new(
            "res_mod_delta",
            DataType::List(Arc::new(Field::new("item", DataType::Float64, true))),
            true,
        ),
        Field::new("nterm_delta", DataType::Float64, false),
        Field::new("cterm_delta", DataType::Float64, false),
        Field::new(
            "mz",
            DataType::List(Arc::new(Field::new("item", DataType::Float32, true))),
            true,
        ),
        Field::new(
            "intensity",
            DataType::List(Arc::new(Field::new("item", DataType::Float32, true))),
            true,
        ),
    ]));

    let seq: ArrayRef = Arc::new(StringArray::from(
        rows.iter().map(|r| r.seq).collect::<Vec<_>>(),
    ));
    let charge: ArrayRef = Arc::new(Int32Array::from(
        rows.iter().map(|r| r.charge).collect::<Vec<_>>(),
    ));
    let prec_mz: ArrayRef = Arc::new(Float64Array::from(
        rows.iter().map(|r| r.prec_mz).collect::<Vec<_>>(),
    ));
    let nterm: ArrayRef = Arc::new(Float64Array::from(
        rows.iter().map(|r| r.nterm_delta).collect::<Vec<_>>(),
    ));
    let cterm: ArrayRef = Arc::new(Float64Array::from(
        rows.iter().map(|r| r.cterm_delta).collect::<Vec<_>>(),
    ));

    let mut pos_b = ListBuilder::new(Int32Builder::new());
    for r in rows {
        for &p in &r.res_mod_pos {
            pos_b.values().append_value(p);
        }
        pos_b.append(true);
    }
    let res_mod_pos: ArrayRef = Arc::new(pos_b.finish());

    let mut delta_b = ListBuilder::new(arrow::array::Float64Builder::new());
    for r in rows {
        for &d in &r.res_mod_delta {
            delta_b.values().append_value(d);
        }
        delta_b.append(true);
    }
    let res_mod_delta: ArrayRef = Arc::new(delta_b.finish());

    let mut mz_b = ListBuilder::new(Float32Builder::new());
    for r in rows {
        for &(m, _) in &r.peaks {
            mz_b.values().append_value(m);
        }
        mz_b.append(true);
    }
    let mz: ArrayRef = Arc::new(mz_b.finish());

    let mut int_b = ListBuilder::new(Float32Builder::new());
    for r in rows {
        for &(_, it) in &r.peaks {
            int_b.values().append_value(it);
        }
        int_b.append(true);
    }
    let intensity: ArrayRef = Arc::new(int_b.finish());

    let batch = RecordBatch::try_new(
        schema.clone(),
        vec![
            seq,
            charge,
            prec_mz,
            res_mod_pos,
            res_mod_delta,
            nterm,
            cterm,
            mz,
            intensity,
        ],
    )
    .expect("build record batch");

    let file = std::fs::File::create(path).expect("create parquet");
    let mut writer = ArrowWriter::try_new(file, schema, None).expect("arrow writer");
    writer.write(&batch).expect("write batch");
    writer.close().expect("close writer");
}

/// A handful of fake peaks at plausible b/y-ion m/z values for PEPTIDEK.
fn fake_peaks() -> Vec<(f32, f32)> {
    vec![
        (98.06, 1000.0),  // b1 (P)
        (227.10, 2000.0), // b2 (PE)
        (324.15, 3000.0), // b3 (PEP)
        (147.11, 1500.0), // y1 (K)
        (276.15, 2500.0), // y2 (EK)
        (500.30, 500.0),
        (600.35, 800.0),
        (700.40, 700.0),
    ]
}

fn synthetic_rows() -> Vec<Row> {
    vec![
        Row {
            seq: "PEPTIDEK",
            charge: 2,
            prec_mz: 472.75,
            res_mod_pos: vec![],
            res_mod_delta: vec![],
            nterm_delta: 0.0,
            cterm_delta: 0.0,
            peaks: fake_peaks(),
        },
        Row {
            seq: "MPEPTIDER",
            charge: 2,
            prec_mz: 545.27,
            res_mod_pos: vec![1],
            res_mod_delta: vec![15.994915],
            nterm_delta: 0.0,
            cterm_delta: 0.0,
            peaks: fake_peaks(),
        },
        Row {
            seq: "SAMPLERK",
            charge: 3,
            prec_mz: 320.5,
            res_mod_pos: vec![],
            res_mod_delta: vec![],
            nterm_delta: 42.010565,
            cterm_delta: 0.0,
            peaks: fake_peaks(),
        },
    ]
}

/// Round-trip test: build flat parquet → `train-intensity-gbdt` CLI → reload
/// store → `frag_intensity_model.is_some()` + `predict_value` is finite.
#[test]
fn train_intensity_gbdt_roundtrip() {
    let dir = tempfile::tempdir().expect("tempdir");
    let in_parquet = dir.path().join("psms.parquet");
    let store_path = dir.path().join("models.parquet");

    write_flat_parquet(&in_parquet, &synthetic_rows());

    let status = Command::new(env!("CARGO_BIN_EXE_andes"))
        .arg("train-intensity-gbdt")
        .arg("--in")
        .arg(&in_parquet)
        .arg("--out-store")
        .arg(&store_path)
        .arg("--model-id")
        .arg("test_frag")
        .arg("--threads")
        .arg("1")
        // The synthetic fixture is tiny (a few PSMs) and intentionally below the
        // GBDT quality gate (finding 3.6); opt into the degenerate fallback so
        // this round-trip test still produces a model to inspect.
        .arg("--allow-degenerate-model")
        .status()
        .expect("run andes train-intensity-gbdt");

    assert!(
        status.success(),
        "train-intensity-gbdt should exit 0, got: {status}"
    );

    assert!(store_path.exists(), "store should be written");

    let store = ModelStore::open(&store_path).expect("open store");
    let ids = store.model_ids();
    assert!(
        ids.contains(&"test_frag".to_string()),
        "store should contain 'test_frag'; found: {ids:?}"
    );

    let param = store.load_param("test_frag").expect("load param");
    assert!(
        param.frag_intensity_model.is_some(),
        "frag_intensity_model must be present after train-intensity-gbdt"
    );

    // The model must produce a finite prediction.
    let fim = param.frag_intensity_model.as_ref().unwrap();
    let n_features = scoring_crate::frag_features::N_FRAG_FEATURES;
    let dummy = vec![0.0f32; n_features];
    let val = fim.predict_value(&dummy);
    assert!(
        val.is_finite(),
        "predict_value on zero features must be finite, got {val}"
    );
}

/// A second model pre-existing in the store must be preserved after
/// `train-intensity-gbdt` writes the new frag model.
#[test]
fn train_intensity_gbdt_preserves_existing_models() {
    let dir = tempfile::tempdir().expect("tempdir");
    let in_parquet = dir.path().join("psms.parquet");
    let store_path = dir.path().join("models.parquet");

    write_flat_parquet(&in_parquet, &synthetic_rows());

    // First: train a rank-core model into the store via train-from-msnet.
    let s1 = Command::new(env!("CARGO_BIN_EXE_andes"))
        .arg("train")
        .arg("--in")
        .arg(&in_parquet)
        .arg("--out-store")
        .arg(&store_path)
        .arg("--model-id")
        .arg("existing_model")
        .arg("--threads")
        .arg("1")
        .arg("--allow-degenerate-model")
        .status()
        .expect("run train-from-msnet");
    assert!(s1.success(), "train-from-msnet should exit 0, got {s1}");

    // Second: run train-intensity-gbdt; must preserve "existing_model".
    let s2 = Command::new(env!("CARGO_BIN_EXE_andes"))
        .arg("train-intensity-gbdt")
        .arg("--in")
        .arg(&in_parquet)
        .arg("--out-store")
        .arg(&store_path)
        .arg("--model-id")
        .arg("frag_model")
        .arg("--threads")
        .arg("1")
        .arg("--allow-degenerate-model")
        .status()
        .expect("run train-intensity-gbdt");
    assert!(s2.success(), "train-intensity-gbdt should exit 0, got {s2}");

    let store = ModelStore::open(&store_path).expect("open store");
    let ids = store.model_ids();
    assert!(
        ids.contains(&"existing_model".to_string()),
        "existing model must be preserved; found: {ids:?}"
    );
    assert!(
        ids.contains(&"frag_model".to_string()),
        "new frag model must be present; found: {ids:?}"
    );

    let param = store.load_param("frag_model").expect("load frag_model");
    assert!(
        param.frag_intensity_model.is_some(),
        "frag_intensity_model must be Some for the newly trained model"
    );
}
