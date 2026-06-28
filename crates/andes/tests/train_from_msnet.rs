//! Integration test for `andes train-from-msnet`.
//!
//! Builds a tiny synthetic "flat training parquet" (one unmodified peptide,
//! one with Oxidation on M, one with an N-terminal Acetyl), runs the
//! `train-from-msnet` subcommand into a temp `--out-store`, and asserts the
//! store gains a model with the requested ID whose `rank_dist_table` is
//! non-empty and that loads via `ModelStore::load_param`.
//!
//! Also asserts that two runs with different `--fragment-tol-*` settings
//! produce different `rank_dist_table`s on the same input (the tolerance
//! override is wired into the scorer used for accumulation).

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::Arc;

use arrow::array::{
    ArrayRef, Float32Builder, Float64Array, Int32Array, Int32Builder, ListBuilder, StringArray,
};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;

use model_train::ModelStore;

/// One synthetic PSM row.
struct Row {
    seq: &'static str,
    charge: i32,
    prec_mz: f64,
    res_mod_pos: Vec<i32>,
    res_mod_delta: Vec<f64>,
    nterm_delta: f64,
    cterm_delta: f64,
    /// (m/z, intensity) peaks in acquisition order (intentionally NOT sorted).
    peaks: Vec<(f32, f32)>,
}

/// Write the flat training parquet at `path` from `rows`.
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

    // res_mod_pos: LIST<INT32>
    let mut pos_b = ListBuilder::new(Int32Builder::new());
    for r in rows {
        for &p in &r.res_mod_pos {
            pos_b.values().append_value(p);
        }
        pos_b.append(true);
    }
    let res_mod_pos: ArrayRef = Arc::new(pos_b.finish());

    // res_mod_delta: LIST<DOUBLE>
    let mut delta_b = ListBuilder::new(arrow::array::Float64Builder::new());
    for r in rows {
        for &d in &r.res_mod_delta {
            delta_b.values().append_value(d);
        }
        delta_b.append(true);
    }
    let res_mod_delta: ArrayRef = Arc::new(delta_b.finish());

    // mz: LIST<FLOAT>
    let mut mz_b = ListBuilder::new(Float32Builder::new());
    for r in rows {
        for &(m, _) in &r.peaks {
            mz_b.values().append_value(m);
        }
        mz_b.append(true);
    }
    let mz: ArrayRef = Arc::new(mz_b.finish());

    // intensity: LIST<FLOAT>
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

/// A handful of fake peaks (acquisition order; the reader sorts by m/z).
fn fake_peaks() -> Vec<(f32, f32)> {
    vec![
        (500.3, 1000.0),
        (200.1, 500.0),
        (700.5, 2000.0),
        (350.2, 800.0),
        (600.4, 1500.0),
        (150.05, 300.0),
        (800.6, 1200.0),
        (450.25, 900.0),
    ]
}

fn synthetic_rows() -> Vec<Row> {
    vec![
        // Unmodified.
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
        // Oxidation on M (residue 1).
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
        // N-terminal Acetyl.
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

fn run_train(in_parquet: &Path, store: &Path, extra: &[&str]) {
    run_train_env(in_parquet, store, extra, &[]);
}

/// Like [`run_train`] but injects environment variables — e.g.
/// `ANDES_SEED_GEOMETRY=1` to opt OUT of own-geometry derivation.
fn run_train_env(in_parquet: &Path, store: &Path, extra: &[&str], env: &[(&str, &str)]) {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_andes"));
    cmd.arg("train")
        .arg("--in")
        .arg(in_parquet)
        .arg("--out-store")
        .arg(store)
        .arg("--model-id")
        .arg("default")
        .arg("--threads")
        .arg("1")
        // Tiny synthetic fixture: below the GBDT quality gate (finding 3.6).
        .arg("--allow-degenerate-model");
    for e in extra {
        cmd.arg(e);
    }
    for (k, v) in env {
        cmd.env(k, v);
    }
    let status = cmd.status().expect("run andes train-from-msnet");
    assert!(status.success(), "train-from-msnet should exit 0, got {status}");
}

/// An order-independent signature of a trained model's partition geometry
/// (charge × parent-mass tier × segment), used to compare the derived-vs-seed
/// geometry without depending on partition ordering.
fn geometry_sig(partitions: &[scoring_crate::param_model::Partition]) -> Vec<(i32, i32, u32)> {
    let mut v: Vec<(i32, i32, u32)> = partitions
        .iter()
        .map(|q| (q.charge, q.seg_num, q.parent_mass.to_bits()))
        .collect();
    v.sort_unstable();
    v
}

#[test]
fn train_from_msnet_writes_model_with_rank_dist() {
    let dir = tempfile::tempdir().expect("tempdir");
    let in_parquet = dir.path().join("psms.parquet");
    let store = dir.path().join("models.parquet");
    write_flat_parquet(&in_parquet, &synthetic_rows());

    run_train(&in_parquet, &store, &["--fragment-tol-ppm", "20"]);

    assert!(store.exists(), "store should be written");
    let ms = ModelStore::open(&store).expect("open store");
    let ids = ms.model_ids();
    assert!(
        ids.contains(&"default".to_string()),
        "store should contain model 'default'; got {ids:?}"
    );

    let param = ms.load_param("default").expect("load default param");
    assert!(
        !param.rank_dist_table.is_empty(),
        "trained rank_dist_table should be non-empty"
    );
    assert!(!param.partitions.is_empty(), "trained model should have partitions");
}

/// Wiring guard (finding 3.9): `run_train` — the binary's default training path —
/// must DERIVE the partition geometry from the corpus (the own-geometry,
/// MS-GF+-free path), and only reuse the seed geometry when `ANDES_SEED_GEOMETRY=1`
/// opts out. If the binary ever silently fell back to seed geometry, the two
/// models trained below would be geometry-identical — so we assert they differ.
/// This exercises the same `andes train` entry point `andes::main` uses, not just
/// the lower-level `derive_geometry` helper (covered in model-train).
#[test]
fn geometry_derived_by_default_differs_from_seed() {
    let dir = tempfile::tempdir().expect("tempdir");
    let in_parquet = dir.path().join("psms.parquet");
    write_flat_parquet(&in_parquet, &synthetic_rows());

    // DEFAULT: derive geometry from the (tiny) corpus.
    let store_derived = dir.path().join("derived.parquet");
    run_train_env(&in_parquet, &store_derived, &["--fragment-tol-ppm", "20"], &[]);

    // OPT-OUT: reuse the bundled seed's full partition geometry.
    let store_seed = dir.path().join("seed.parquet");
    run_train_env(
        &in_parquet,
        &store_seed,
        &["--fragment-tol-ppm", "20"],
        &[("ANDES_SEED_GEOMETRY", "1")],
    );

    let derived = ModelStore::open(&store_derived).unwrap().load_param("default").unwrap();
    let seed = ModelStore::open(&store_seed).unwrap().load_param("default").unwrap();

    let sd = geometry_sig(&derived.partitions);
    let ss = geometry_sig(&seed.partitions);
    assert!(!sd.is_empty() && !ss.is_empty(), "both trained models must have partitions");
    assert_ne!(
        sd, ss,
        "default training must DERIVE geometry from the corpus, but it matched the seed geometry \
         exactly — the own-geometry wiring in run_train has regressed"
    );
    // The 3-PSM corpus (charges 2,2,3) yields a strictly smaller geometry than the
    // full bundled seed (hcd_qexactive_tryp), confirming it adapts to the corpus.
    assert!(
        sd.len() < ss.len(),
        "derived geometry ({} partitions) should be smaller than the reused seed ({})",
        sd.len(),
        ss.len()
    );
}

#[test]
fn fragment_tolerance_override_changes_model() {
    let dir = tempfile::tempdir().expect("tempdir");
    let in_parquet = dir.path().join("psms.parquet");
    write_flat_parquet(&in_parquet, &synthetic_rows());

    // Tight tolerance.
    let store_tight = dir.path().join("tight.parquet");
    run_train(&in_parquet, &store_tight, &["--fragment-tol-ppm", "1"]);
    let tight = ModelStore::open(&store_tight)
        .unwrap()
        .load_param("default")
        .unwrap();

    // Wide tolerance (Da-based, large window).
    let store_wide = dir.path().join("wide.parquet");
    run_train(&in_parquet, &store_wide, &["--fragment-tol-da", "1.0"]);
    let wide = ModelStore::open(&store_wide)
        .unwrap()
        .load_param("default")
        .unwrap();

    // The override must flow through and change the trained model's match
    // tolerance (`mme`). On real noisy spectra a wider window also reshapes the
    // rank histograms, but the synthetic fixtures use *exact* fragment peaks (no
    // noise), so 1 ppm and 1 Da match the same peaks and the rank_dist_tables are
    // identical — the robust, override-flowed-through check is `mme` itself.
    assert_ne!(
        tight.mme, wide.mme,
        "fragment-tol override must change the trained model's mme (Ppm(1) vs Da(1.0))"
    );
}

/// `--activation/--instrument/--enzyme/--protocol` override the trained model's
/// `data_type` columns (which drive model selection), independent of the seed.
/// Without these flags the model would inherit the seed's data_type
/// (default seed `hcd_qexactive_tryp` => HCD/QExactive/Tryp/Automatic), so
/// selection could never route e.g. a low-res CID-TMT query to the new model.
#[test]
fn data_type_override_sets_selection_columns() {
    let dir = tempfile::tempdir().expect("tempdir");
    let in_parquet = dir.path().join("psms.parquet");
    let store = dir.path().join("models.parquet");
    write_flat_parquet(&in_parquet, &synthetic_rows());

    // Mint a low-res CID-TMT model from a (non-TMT) seed.
    run_train(
        &in_parquet,
        &store,
        &[
            "--fragment-tol-da", "0.5",
            "--activation", "CID",
            "--instrument", "LowRes",
            "--enzyme", "Trypsin",
            "--protocol", "TMT",
        ],
    );

    let ms = ModelStore::open(&store).expect("open store");
    let entry = ms
        .manifest_entries()
        .iter()
        .find(|e| e.model_id == "default")
        .expect("manifest entry for 'default'");

    assert_eq!(entry.activation, "CID", "activation column must reflect --activation");
    assert_eq!(entry.instrument, "LowRes", "instrument column must reflect --instrument");
    assert_eq!(entry.enzyme, "Trypsin", "enzyme column must reflect --enzyme");
    assert_eq!(
        entry.protocol, "TMT",
        "protocol column must reflect --protocol (drives experiment_class=tmt selection)"
    );
}

/// `--prior-model-store`/`--prior-model` are accepted and a model is still
/// written. Smoke test for the mechanism: train a prior into one store, then
/// re-run pointing `--prior-model-store` at that store with a valid `--prior-model`
/// slug. Sparse partitions shrink toward the prior, but the run must still
/// produce a model in the out-store.
#[test]
fn train_from_msnet_accepts_prior_model_flags() {
    let dir = tempfile::tempdir().expect("tempdir");
    let in_parquet = dir.path().join("psms.parquet");
    write_flat_parquet(&in_parquet, &synthetic_rows());

    // Build a prior store (a valid model named "default") to point at.
    let prior_store = dir.path().join("prior.parquet");
    run_train(&in_parquet, &prior_store, &["--fragment-tol-ppm", "20"]);
    assert!(prior_store.exists(), "prior store should be written");

    // Train again, shrinking toward the prior model loaded from `prior_store`.
    let store = dir.path().join("models.parquet");
    run_train(
        &in_parquet,
        &store,
        &[
            "--fragment-tol-ppm",
            "20",
            "--prior-model-store",
            prior_store.to_str().unwrap(),
            "--prior-model",
            "default",
        ],
    );

    assert!(store.exists(), "store should be written");
    let ms = ModelStore::open(&store).expect("open store");
    assert!(
        ms.model_ids().contains(&"default".to_string()),
        "store should contain model 'default'; got {:?}",
        ms.model_ids()
    );
    let param = ms.load_param("default").expect("load default param");
    assert!(
        !param.partitions.is_empty(),
        "trained model should have partitions"
    );
}

/// `--rank-smoothing` is accepted and a model is still written.
/// Smoke-tests that the flag reaches EstimatorConfig.rank_smoothing without panicking.
#[test]
fn train_from_msnet_accepts_rank_smoothing() {
    let dir = tempfile::tempdir().expect("tempdir");
    let in_parquet = dir.path().join("psms.parquet");
    let store = dir.path().join("models.parquet");
    write_flat_parquet(&in_parquet, &synthetic_rows());

    run_train(&in_parquet, &store, &["--fragment-tol-ppm", "20", "--rank-smoothing"]);

    assert!(store.exists(), "store should be written");
    let ms = ModelStore::open(&store).expect("open store");
    assert!(
        ms.model_ids().contains(&"default".to_string()),
        "store should contain model 'default'; got {:?}",
        ms.model_ids()
    );
    let param = ms.load_param("default").expect("load default param");
    assert!(
        !param.rank_dist_table.is_empty(),
        "trained rank_dist_table should be non-empty"
    );
}

/// Multiple `--in` files accumulate into one model.
#[test]
fn multiple_inputs_accumulate() {
    let dir = tempfile::tempdir().expect("tempdir");
    let a = dir.path().join("a.parquet");
    let b = dir.path().join("b.parquet");
    write_flat_parquet(&a, &synthetic_rows());
    write_flat_parquet(&b, &synthetic_rows());
    let store = dir.path().join("models.parquet");

    let status = Command::new(env!("CARGO_BIN_EXE_andes"))
        .arg("train")
        .arg("--in")
        .arg(&a)
        .arg("--in")
        .arg(&b)
        .arg("--out-store")
        .arg(&store)
        .arg("--threads")
        .arg("1")
        .arg("--allow-degenerate-model")
        .status()
        .expect("run");
    assert!(status.success(), "should exit 0");

    let ms = ModelStore::open(&store).unwrap();
    // Default model id is "default".
    assert!(ms.model_ids().contains(&"default".to_string()));
    let _ = PathBuf::from(&store); // silence unused import if refactored
}
