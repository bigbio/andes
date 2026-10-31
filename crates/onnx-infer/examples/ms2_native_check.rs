//! Parity and speed check for `ms2_native` against the PyTorch reference.
//!
//! ```text
//! cargo run --release -p onnx-infer --example ms2_native_check -- \
//!   --weights andes_ms2.ams2 --reference ms2_reference.tsv [--threads 8] [--repeat 5]
//! ```
//!
//! The reference TSV (from `scripts/student/ms2/ms2_export.py`) has columns
//! `format` (`proforma` | `andes`), `peptide`, `charge`, `nce` (`NA` = unknown),
//! `activation`, `instrument`, `values` and `raw` (comma-separated, site-major).
//! Prints the largest absolute and relative differences on both and the
//! single-thread / multi-thread throughput.

use std::time::Instant;

use onnx_infer::ms2_native::{Ms2Conditions, Ms2Model, Ms2Scratch};
use onnx_infer::rt_native::RtPeptide;

fn main() {
    let v: Vec<String> = std::env::args().skip(1).collect();
    let get = |k: &str| v.iter().position(|x| x == k).map(|i| v[i + 1].clone());
    let model = Ms2Model::load(get("--weights").expect("--weights")).expect("load weights");
    let threads: usize = get("--threads").map_or(8, |t| t.parse().unwrap());
    let repeat: usize = get("--repeat").map_or(5, |t| t.parse().unwrap());
    println!("model params {}", model.n_params());

    let text = std::fs::read_to_string(get("--reference").expect("--reference")).unwrap();
    let mut lines = text.lines();
    let h: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let col = |n: &str| {
        h.iter()
            .position(|x| *x == n)
            .unwrap_or_else(|| panic!("column {n}"))
    };
    let (cf, cp, cz, cn, ca, ci, cv, cr) = (
        col("format"),
        col("peptide"),
        col("charge"),
        col("nce"),
        col("activation"),
        col("instrument"),
        col("values"),
        col("raw"),
    );
    let parse = |s: &str| {
        s.split(',')
            .map(|x| x.parse::<f32>().unwrap())
            .collect::<Vec<_>>()
    };
    let mut items = Vec::new();
    let (mut abs_err, mut rel_err, mut n, mut n_unsup) = (0f32, 0f32, 0usize, 0usize);
    let mut worst = String::new();
    let mut s = Ms2Scratch::default();
    for line in lines {
        let x: Vec<&str> = line.split('\t').collect();
        let p = match x[cf] {
            "proforma" => RtPeptide::from_proforma(x[cp]),
            "andes" => RtPeptide::from_andes(x[cp]),
            other => panic!("format {other}"),
        };
        let c = Ms2Conditions::new(x[cz].parse().unwrap(), x[cn].parse().ok(), x[ca], x[ci]);
        let got = p.as_ref().and_then(|p| model.predict_with(p, &c, &mut s));
        match got {
            Some(got) => {
                assert!(x[cv] != "NA", "rust encodes {} but torch did not", x[cp]);
                let raw = model
                    .predict_raw_with(p.as_ref().unwrap(), &c, &mut s)
                    .unwrap();
                let (want, want_raw) = (parse(x[cv]), parse(x[cr]));
                assert_eq!(got.len(), want.len(), "{}: output length", x[cp]);
                for (g, w) in got.iter().zip(&want).chain(raw.iter().zip(&want_raw)) {
                    abs_err = abs_err.max((g - w).abs());
                    let rel = (g - w).abs() / (1.0 + w.abs());
                    if rel > rel_err {
                        worst = format!("{} z{} {} {} {}", x[cp], x[cz], x[cn], x[ca], x[ci]);
                    }
                    rel_err = rel_err.max(rel);
                }
                n += 1;
                items.push((p.unwrap(), c));
            }
            None => {
                assert!(x[cv] == "NA", "torch encodes {} but rust did not", x[cp]);
                n_unsup += 1;
            }
        }
    }
    println!("PARITY n={n} unsupported(both)={n_unsup} max|rust - torch| = {abs_err:.3e} (relative {rel_err:.3e})");
    println!("worst row: {worst}");

    let mean_len =
        items.iter().map(|(p, _)| p.sequence.len()).sum::<usize>() as f64 / items.len() as f64;
    let t = Instant::now();
    let mut acc = 0.0f32;
    for _ in 0..repeat {
        for (p, c) in &items {
            acc += model.predict_with(p, c, &mut s).unwrap()[0];
        }
    }
    let per = t.elapsed().as_secs_f64() / (repeat * items.len()) as f64;
    println!(
        "SPEED 1 thread: {:.1} us/peptide (mean length {mean_len:.1}) [{acc:.1}]",
        per * 1e6
    );

    let t = Instant::now();
    std::thread::scope(|sc| {
        let chunk = items.len().div_ceil(threads);
        for part in items.chunks(chunk) {
            let model = &model;
            sc.spawn(move || {
                let mut s = Ms2Scratch::default();
                let mut a = 0.0f32;
                for _ in 0..repeat {
                    for (p, c) in part {
                        a += model.predict_with(p, c, &mut s).unwrap()[0];
                    }
                }
                a
            });
        }
    });
    let wall = t.elapsed().as_secs_f64();
    println!(
        "SPEED {threads} threads: {:.0} peptides/s, {:.1} us/peptide/thread",
        (repeat * items.len()) as f64 / wall,
        wall * threads as f64 / (repeat * items.len()) as f64 * 1e6
    );
}
