//! Parity and speed check for `rt_native` against the PyTorch reference.
//!
//! ```text
//! cargo run --release -p onnx-infer --example rt_native_check -- \
//!   --weights andes_rt.artn --reference rt_reference.tsv [--threads 8] [--repeat 20]
//! ```
//!
//! The reference TSV (from `scripts/student/rt_export.py`) has columns
//! `format` (`proforma` | `andes`), `peptide`, `index`. Prints the largest
//! absolute difference and single-thread / multi-thread throughput.

use std::io::{BufRead, BufReader};
use std::time::Instant;

use onnx_infer::rt_native::{RtModel, RtPeptide, Scratch};

fn main() {
    let v: Vec<String> = std::env::args().skip(1).collect();
    let get = |k: &str| v.iter().position(|x| x == k).map(|i| v[i + 1].clone());
    let model = RtModel::load(get("--weights").expect("--weights")).expect("load weights");
    let threads: usize = get("--threads").map_or(8, |t| t.parse().unwrap());
    let repeat: usize = get("--repeat").map_or(20, |t| t.parse().unwrap());
    println!("model params {}", model.n_params());

    let f = BufReader::new(std::fs::File::open(get("--reference").expect("--reference")).unwrap());
    let mut lines = f.lines();
    let h: Vec<String> = lines
        .next()
        .unwrap()
        .unwrap()
        .split('\t')
        .map(String::from)
        .collect();
    let col = |n: &str| {
        h.iter()
            .position(|x| x == n)
            .unwrap_or_else(|| panic!("column {n}"))
    };
    let (cf, cp, ci) = (col("format"), col("peptide"), col("index"));
    let mut peps = Vec::new();
    let (mut max_err, mut n, mut n_unsup) = (0f64, 0usize, 0usize);
    for line in lines {
        let line = line.unwrap();
        let x: Vec<&str> = line.split('\t').collect();
        let p = match x[cf] {
            "proforma" => RtPeptide::from_proforma(x[cp]),
            "andes" => RtPeptide::from_andes(x[cp]),
            other => panic!("format {other}"),
        };
        let want: f64 = x[ci].parse().unwrap_or(f64::NAN);
        match p.as_ref().and_then(|p| model.predict(p)) {
            Some(got) => {
                assert!(want.is_finite(), "rust encodes {} but torch did not", x[cp]);
                max_err = max_err.max((got as f64 - want).abs());
                n += 1;
                peps.push(p.unwrap());
            }
            None => {
                assert!(
                    !want.is_finite(),
                    "torch encodes {} but rust did not",
                    x[cp]
                );
                n_unsup += 1;
            }
        }
    }
    println!("PARITY n={n} unsupported(both)={n_unsup} max|rust - torch| = {max_err:.3e}");

    let mean_len = peps.iter().map(|p| p.sequence.len()).sum::<usize>() as f64 / peps.len() as f64;
    let mut s = Scratch::default();
    let t = Instant::now();
    let mut acc = 0.0f32;
    for _ in 0..repeat {
        for p in &peps {
            acc += model.predict_with(p, &mut s).unwrap();
        }
    }
    let per = t.elapsed().as_secs_f64() / (repeat * peps.len()) as f64;
    println!(
        "SPEED 1 thread: {:.1} us/peptide (mean length {mean_len:.1}) [{acc:.1}]",
        per * 1e6
    );

    let t = Instant::now();
    std::thread::scope(|sc| {
        let chunk = peps.len().div_ceil(threads);
        for part in peps.chunks(chunk) {
            let model = &model;
            sc.spawn(move || {
                let mut s = Scratch::default();
                let mut a = 0.0f32;
                for _ in 0..repeat {
                    for p in part {
                        a += model.predict_with(p, &mut s).unwrap();
                    }
                }
                a
            });
        }
    });
    let wall = t.elapsed().as_secs_f64();
    println!(
        "SPEED {threads} threads: {:.0} peptides/s, {:.1} us/peptide/thread",
        (repeat * peps.len()) as f64 / wall,
        wall * threads as f64 / (repeat * peps.len()) as f64 * 1e6
    );
}
