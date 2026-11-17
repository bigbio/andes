use percolator_native::{run, MemLimit};
use std::fmt::Write as _;
use std::path::{Path, PathBuf};

const ROWS_PER_FILE: usize = 1000;
const FILES: usize = 50;
const NFEAT: usize = 100;

struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64) / ((1u64 << 53) as f64)
    }

    fn gauss(&mut self) -> f64 {
        let mut s = 0.0;
        for _ in 0..6 {
            s += self.next();
        }
        (s - 3.0) * 1.4142135623730951
    }
}

fn synthetic_pin() -> String {
    let mut rng = Lcg(20261116);
    let mut out = String::new();
    out.push_str("SpecId\tLabel\tScanNr\tExpMass\tFileName");
    for f in 0..NFEAT {
        write!(out, "\tf{f}").unwrap();
    }
    out.push_str("\tPeptide\tProteins\n");
    let mut id = 0usize;
    for file in 0..FILES {
        for k in 0..ROWS_PER_FILE / 2 {
            let scan = (k / 2 + 1) as u64;
            let mass = 800.0 + ((k * 37 + file * 11) % 4000) as f64 * 0.25;
            let true_hit = rng.next() < 0.55;
            for decoy in [false, true] {
                let shift = if !decoy && true_hit { 1.6 } else { 0.0 };
                let label = if decoy { -1 } else { 1 };
                write!(out, "r{id}\t{label}\t{scan}\t{mass:.4}\tfile{file}.mzML").unwrap();
                for f in 0..NFEAT {
                    let v = match f {
                        0..=5 => rng.gauss() + shift * (1.0 - f as f64 * 0.12),
                        6..=9 => ((rng.gauss() + shift).round()).clamp(-3.0, 3.0),
                        10..=19 => (rng.gauss() * 4.0).round() / 4.0,
                        20..=24 => 3.5,
                        25..=29 => {
                            if rng.next() < 0.001 {
                                1.0
                            } else {
                                0.0
                            }
                        }
                        30..=39 => (rng.next() * 10.0).floor(),
                        _ => rng.gauss() * (1.0 + (f % 7) as f64) + 100.0 * (f % 3) as f64,
                    };
                    write!(out, "\t{v:.5}").unwrap();
                }
                let pep = format!("K.PEP{}TIDE{}K.R", id % 9000, (id / 3) % 17);
                let prot = if decoy {
                    format!("decoy_P{}", id % 5000)
                } else {
                    format!("P{}", id % 5000)
                };
                writeln!(out, "\t{pep}\t{prot}").unwrap();
                id += 1;
            }
        }
    }
    out
}

struct Outputs {
    target: Vec<u8>,
    decoy: Vec<u8>,
    weights: Option<Vec<u8>>,
}

fn p(path: &Path) -> String {
    path.to_string_lossy().into_owned()
}

fn strs(v: &[String]) -> Vec<String> {
    v.to_vec()
}

fn training(pin: &Path, dir: &Path, threads: usize, limit: u64) -> Outputs {
    let t = dir.join("t.psms");
    let d = dir.join("d.psms");
    let w = dir.join("w.tsv");
    let args = vec![
        "--seed".to_string(),
        "42".to_string(),
        "--results-psms".to_string(),
        p(&t),
        "--decoy-results-psms".to_string(),
        p(&d),
        "--weights".to_string(),
        p(&w),
        "--only-psms".to_string(),
        "-Y".to_string(),
        "--num-threads".to_string(),
        threads.to_string(),
        p(pin),
    ];
    run(
        &strs(&args),
        MemLimit {
            bytes: limit,
            spill_dir: dir.to_path_buf(),
        },
    )
    .unwrap();
    Outputs {
        target: std::fs::read(&t).unwrap(),
        decoy: std::fs::read(&d).unwrap(),
        weights: Some(std::fs::read(&w).unwrap()),
    }
}

fn static_pass(
    pin: &Path,
    init: &Path,
    dir: &Path,
    threads: usize,
    limit: u64,
    with_weights: bool,
) -> Outputs {
    let t = dir.join("st.psms");
    let d = dir.join("sd.psms");
    let w = dir.join("sw.tsv");
    let mut args = vec![
        "--seed".to_string(),
        "42".to_string(),
        "--results-psms".to_string(),
        p(&t),
        "--decoy-results-psms".to_string(),
        p(&d),
        "--only-psms".to_string(),
        "--static".to_string(),
        "--init-weights".to_string(),
        p(init),
        "--override".to_string(),
        "-Y".to_string(),
        "--num-threads".to_string(),
        threads.to_string(),
        p(pin),
    ];
    if with_weights {
        args.insert(6, "--weights".to_string());
        args.insert(7, p(&w));
    }
    run(
        &strs(&args),
        MemLimit {
            bytes: limit,
            spill_dir: dir.to_path_buf(),
        },
    )
    .unwrap();
    Outputs {
        target: std::fs::read(&t).unwrap(),
        decoy: std::fs::read(&d).unwrap(),
        weights: with_weights.then(|| std::fs::read(&w).unwrap()),
    }
}

fn same(label: &str, a: &Outputs, b: &Outputs) {
    assert!(a.target.len() > 100, "{label}: empty target output");
    assert!(a.target == b.target, "{label}: target file differs");
    assert!(a.decoy == b.decoy, "{label}: decoy file differs");
    assert!(a.weights == b.weights, "{label}: weights file differs");
    assert!(a.weights.is_some() || b.weights.is_none());
}

fn write_pin(dir: &Path) -> PathBuf {
    let pin = dir.join("synthetic.pin");
    std::fs::write(&pin, synthetic_pin()).unwrap();
    pin
}

#[test]
fn spill_path_crosses_block_boundaries() {
    let rows = FILES * ROWS_PER_FILE;
    let block_rows = (8usize << 20) / (4 * NFEAT);
    let blocks = rows.div_ceil(block_rows);
    assert!(rows >= 50_000);
    assert!(blocks >= 3, "only {blocks} blocks");
    assert!(rows % block_rows != 0, "last block is full");
}

#[test]
fn output_identical_at_any_memory_limit() {
    let dir = tempfile::tempdir().unwrap();
    let pin = write_pin(dir.path());
    let mut runs = Vec::new();
    for (i, limit) in [64u64 << 30, 1 << 20, 1].into_iter().enumerate() {
        let sub = dir.path().join(format!("train{i}"));
        std::fs::create_dir(&sub).unwrap();
        runs.push((sub.clone(), training(&pin, &sub, 3, limit), limit));
    }
    for (_, o, limit) in &runs[1..] {
        same(&format!("training at {limit} bytes"), &runs[0].1, o);
    }
    let init = runs[0].0.join("w.tsv");
    for with_weights in [true, false] {
        let mut statics = Vec::new();
        for (i, limit) in [64u64 << 30, 1 << 20, 1].into_iter().enumerate() {
            let sub = dir.path().join(format!("static{with_weights}{i}"));
            std::fs::create_dir(&sub).unwrap();
            statics.push((
                static_pass(&pin, &init, &sub, 3, limit, with_weights),
                limit,
            ));
        }
        for (o, limit) in &statics[1..] {
            same(
                &format!("static (weights {with_weights}) at {limit} bytes"),
                &statics[0].0,
                o,
            );
        }
    }
}

#[test]
fn output_identical_at_1_and_16_threads() {
    let dir = tempfile::tempdir().unwrap();
    let pin = write_pin(dir.path());
    let mut runs = Vec::new();
    for (i, threads) in [1usize, 16].into_iter().enumerate() {
        let sub = dir.path().join(format!("train{i}"));
        std::fs::create_dir(&sub).unwrap();
        runs.push((sub.clone(), training(&pin, &sub, threads, 1 << 20)));
    }
    same("training at 16 threads", &runs[0].1, &runs[1].1);
    let init = runs[0].0.join("w.tsv");
    for with_weights in [true, false] {
        let mut statics = Vec::new();
        for (i, threads) in [1usize, 16].into_iter().enumerate() {
            let sub = dir.path().join(format!("static{with_weights}{i}"));
            std::fs::create_dir(&sub).unwrap();
            statics.push(static_pass(
                &pin,
                &init,
                &sub,
                threads,
                1 << 20,
                with_weights,
            ));
        }
        same(
            &format!("static (weights {with_weights}) at 16 threads"),
            &statics[0],
            &statics[1],
        );
    }
}
