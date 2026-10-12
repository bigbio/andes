//! Parity of onnx-infer (tract) predictions with a model's reference Python
//! implementation, from the TSV its `scripts/<model>_reference.py` writes.
//!
//! ```text
//! cargo run --release -p onnx-infer --features tract --example model_parity -- \
//!   --model ms2pip --reference ms2pip_reference.tsv --models DIR
//! ```
//!
//! Prints the largest absolute and relative deviation per model variant.

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

use onnx_infer::{BackendKind, OnnxSession, ms2pip};

struct Args {
    model: String,
    reference: PathBuf,
    models: PathBuf,
}

fn parse_args() -> Args {
    let v: Vec<String> = std::env::args().skip(1).collect();
    let get = |k: &str| {
        v.iter()
            .position(|a| a == k)
            .and_then(|i| v.get(i + 1))
            .unwrap_or_else(|| panic!("missing {k}"))
            .clone()
    };
    Args {
        model: get("--model"),
        reference: get("--reference").into(),
        models: get("--models").into(),
    }
}

fn floats(s: &str) -> Vec<f32> {
    s.split(',').filter(|x| !x.is_empty()).map(|x| x.parse().expect("float")).collect()
}

#[derive(Default)]
struct Dev {
    abs: f64,
    rel: f64,
    n: usize,
}

impl Dev {
    fn add(&mut self, got: f32, want: f32) {
        let d = (got as f64 - want as f64).abs();
        self.abs = self.abs.max(d);
        self.rel = self.rel.max(d / (want as f64).abs().max(1e-6));
        self.n += 1;
    }
}

fn load(dir: &Path, stem: &str) -> OnnxSession {
    let t = Instant::now();
    let s = OnnxSession::from_manifest(dir.join(format!("{stem}.json")), BackendKind::Tract)
        .unwrap_or_else(|e| panic!("load {stem}: {e}"));
    println!("load {stem}: {:.2}s", t.elapsed().as_secs_f64());
    s
}

fn ms2pip_parity(a: &Args) {
    let text = fs::read_to_string(&a.reference).expect("reference TSV");
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let col = |n: &str| header.iter().position(|h| *h == n).unwrap();
    let mut by_model: BTreeMap<String, Vec<Vec<&str>>> = BTreeMap::new();
    for l in lines {
        let x: Vec<&str> = l.split('\t').collect();
        by_model.entry(x[col("model")].to_string()).or_default().push(x);
    }
    for (model, rows) in by_model {
        let mut sb = load(&a.models, &format!("ms2pip_{model}_b"));
        let mut sy = load(&a.models, &format!("ms2pip_{model}_y"));
        let (mut feat, mut pb, mut py) = (Dev::default(), Dev::default(), Dev::default());
        for x in &rows {
            let seq = x[col("sequence")];
            let z: u8 = x[col("charge")].parse().unwrap();
            let batch = ms2pip::pack_batch(&[(seq, z)]).unwrap();
            for (g, w) in batch.features.data.iter().zip(floats(x[col("features")])) {
                feat.add(*g, w);
            }
            let input = std::slice::from_ref(&batch.features);
            let b = ms2pip::log2_predictions(&sb.run(input).unwrap()[1]).unwrap();
            let y = ms2pip::log2_predictions(&sy.run(input).unwrap()[1]).unwrap();
            let (rb, ry) = (floats(x[col("b_log2")]), floats(x[col("y_log2")]));
            assert_eq!((b.len(), y.len()), (rb.len(), ry.len()), "{seq} row count");
            b.iter().zip(&rb).for_each(|(g, w)| pb.add(*g, *w));
            y.iter().zip(&ry).for_each(|(g, w)| py.add(*g, *w));
        }
        println!(
            "PARITY ms2pip {model}: {} peptides | features max|Δ| {:.3e} ({} values) | \
             b log2 max|Δ| {:.3e} rel {:.3e} | y log2 max|Δ| {:.3e} rel {:.3e}",
            rows.len(),
            feat.abs,
            feat.n,
            pb.abs,
            pb.rel,
            py.abs,
            py.rel
        );
    }
}

fn main() {
    let a = parse_args();
    match a.model.as_str() {
        "ms2pip" => ms2pip_parity(&a),
        other => panic!("unknown --model {other}"),
    }
}
