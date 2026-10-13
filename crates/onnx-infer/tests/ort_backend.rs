//! ONNX Runtime backend on the fixture graphs. Needs the runtime library at
//! test time (`ORT_DYLIB_PATH`, or `libonnxruntime` on the loader path).

#![cfg(feature = "ort")]

use std::path::PathBuf;

use onnx_infer::{BackendKind, BackendOptions, NamedTensor, OnnxSession};

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("fixtures")
        .join(name)
}

#[test]
fn doubles_from_path_and_manifest() {
    let mut s = OnnxSession::from_onnx(fixture("double.onnx"), BackendKind::Ort).expect("load");
    assert_eq!(s.backend_kind(), BackendKind::Ort);
    assert_eq!(s.input_names(), ["x"]);
    assert_eq!(s.output_names(), ["y"]);
    let input = NamedTensor::new("x", vec![2, 4], (1..=8).map(|v| v as f32).collect()).unwrap();
    let y = s.run_one(input).expect("run");
    assert_eq!(y.shape, vec![2, 4]);
    assert_eq!(y.data, (1..=8).map(|v| 2.0 * v as f32).collect::<Vec<_>>());

    let mut m =
        OnnxSession::from_manifest(fixture("double.json"), BackendKind::Ort).expect("manifest");
    assert_eq!(m.manifest().unwrap().id, "fixture-double-v1");
    let y = m
        .run_one(NamedTensor::new("x", vec![1, 4], vec![0.5, 1.0, 1.5, 2.0]).unwrap())
        .unwrap();
    assert_eq!(y.data, vec![1.0, 2.0, 3.0, 4.0]);
}

#[test]
fn reports_onnx_tensor_names() {
    let mut s =
        OnnxSession::from_manifest(fixture("intensity_mlp.json"), BackendKind::Ort).unwrap();
    assert_eq!(s.input_names(), ["features"]);
    assert_eq!(s.output_names(), ["score"]);
    let out = s
        .run_one(NamedTensor::new("features", vec![1, 4], vec![0.0; 4]).unwrap())
        .unwrap();
    assert_eq!(out.name, "score");
}

#[test]
fn rejects_unknown_input_name() {
    let mut s =
        OnnxSession::from_manifest(fixture("intensity_mlp.json"), BackendKind::Ort).unwrap();
    let a = NamedTensor::new("features", vec![1, 4], vec![0.0; 4]).unwrap();
    let b = NamedTensor::new("nope", vec![1, 4], vec![0.0; 4]).unwrap();
    assert!(s.run(&[a, b]).is_err());
}

#[test]
fn intensity_mlp_relu_scores() {
    let opts = BackendOptions { intra_threads: 2 };
    let mut s =
        OnnxSession::from_manifest_with(fixture("intensity_mlp.json"), BackendKind::Ort, &opts)
            .expect("load intensity_mlp");
    let input = NamedTensor::new(
        "features",
        vec![3, 4],
        vec![1., 0., 0., 0., 0., 0., 1., 0., -1., -1., -1., -1.],
    )
    .unwrap();
    let out = s.run_one(input).unwrap();
    assert_eq!(out.shape, vec![3, 1]);
    assert!((out.data[0] - 0.55).abs() < 1e-5, "got {}", out.data[0]);
    assert!((out.data[1] - 0.80).abs() < 1e-5, "got {}", out.data[1]);
    assert!(out.data[2].abs() < 1e-6, "got {}", out.data[2]);
}

/// `ai.onnx.ml` TreeEnsembleRegressor: feature0 <= 0.5 → 1.0 else 3.0.
#[test]
fn tree_ensemble_regressor() {
    let mut s = OnnxSession::from_manifest(fixture("tree_leq.json"), BackendKind::Ort)
        .expect("load tree_leq");
    let out = s
        .run_one(NamedTensor::new("features", vec![3, 1], vec![0.2, 0.5, 0.9]).unwrap())
        .unwrap();
    assert_eq!(out.shape, vec![3, 1]);
    assert_eq!(out.data, vec![1.0, 1.0, 3.0]);

    // Multi-tree ensemble: one finite score per row.
    let mut g = OnnxSession::from_manifest(fixture("peak_gbdt.json"), BackendKind::Ort)
        .expect("load peak_gbdt");
    let out = g
        .run_one(
            NamedTensor::new(
                "features",
                vec![4, 3],
                (0..12).map(|i| i as f32 / 4.0).collect(),
            )
            .unwrap(),
        )
        .unwrap();
    assert_eq!(out.shape, vec![4, 1]);
    assert!(out.data.iter().all(|v| v.is_finite()));
}

/// Shared handles run the same session concurrently and agree with a
/// single-threaded run.
#[test]
fn shared_handles_run_concurrently() {
    let master =
        OnnxSession::from_manifest(fixture("intensity_mlp.json"), BackendKind::Ort).unwrap();
    let rows = 64;
    let data: Vec<f32> = (0..rows * 4)
        .map(|i| ((i * 37 % 23) as f32 - 11.0) / 7.0)
        .collect();
    let want = master
        .share()
        .unwrap()
        .run_one(NamedTensor::new("features", vec![rows, 4], data.clone()).unwrap())
        .unwrap();
    std::thread::scope(|sc| {
        for _ in 0..8 {
            let mut s = master.share().unwrap();
            let (data, want) = (&data, &want);
            sc.spawn(move || {
                for _ in 0..200 {
                    let got = s
                        .run_one(NamedTensor::new("features", vec![rows, 4], data.clone()).unwrap())
                        .unwrap();
                    assert_eq!(got.data, want.data);
                }
            });
        }
    });
}

/// The same graph gives the same numbers on both engines (tract has no
/// `ai.onnx.ml` tree ensembles, so only the dense fixtures are compared).
#[cfg(feature = "tract")]
#[test]
fn matches_tract() {
    for (manifest, cols) in [("intensity_mlp.json", 4usize), ("double.json", 4)] {
        let rows = 50;
        let data: Vec<f32> = (0..rows * cols)
            .map(|i| ((i * 53 % 29) as f32 - 14.0) / 9.0)
            .collect();
        let x = || NamedTensor::new("in", vec![rows, cols], data.clone()).unwrap();
        let mut o = OnnxSession::from_manifest(fixture(manifest), BackendKind::Ort).unwrap();
        let mut t = OnnxSession::from_manifest(fixture(manifest), BackendKind::Tract).unwrap();
        let (a, b) = (o.run_one(x()).unwrap(), t.run_one(x()).unwrap());
        assert_eq!(a.shape, b.shape, "{manifest}");
        for (p, q) in a.data.iter().zip(&b.data) {
            assert!(
                (p - q).abs() <= 1e-5 * q.abs().max(1.0),
                "{manifest}: ort {p} tract {q}"
            );
        }
    }
}

/// A graph without nodes whose input `x` is also its output: the output
/// ONNX Runtime returns aliases the input buffer.
#[test]
fn output_aliasing_an_input() {
    let mut s = OnnxSession::from_onnx(fixture("identity.onnx"), BackendKind::Ort).expect("load");
    let y = s
        .run_one(NamedTensor::new("x", vec![1], vec![17.0]).unwrap())
        .unwrap();
    assert_eq!((y.shape, y.data), (vec![1], vec![17.0]));
    let x: Vec<f32> = (0..4096).map(|v| v as f32).collect();
    let y = s
        .run_one(NamedTensor::new("x", vec![4096], x.clone()).unwrap())
        .unwrap();
    assert_eq!(y.data, x);
}
