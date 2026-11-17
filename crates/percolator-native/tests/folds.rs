use percolator_native::{
    count_at, folds_by_spectrum, init_direction, libstdcxx_sort_by, read_pin, score, stdv, tdc_q,
    MemLimit,
};
use std::path::PathBuf;

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

fn sorted(keys: Vec<(u32, i32)>) -> Vec<i32> {
    let mut v = keys;
    libstdcxx_sort_by(&mut v, |a, b| a.0 < b.0);
    v.iter().map(|p| p.1).collect()
}

#[test]
fn sort_matches_libstdcxx_std_sort() {
    let reference = std::fs::read_to_string(data("libstdcxx_sort_ref.txt")).unwrap();
    let mut want = reference.lines();
    let mut s = 7u64;
    let mut checked = 0;
    for n in [1, 2, 5, 16, 17, 33, 100, 1000, 5000] {
        for range in [3u64, 50, 1_000_000] {
            let keys: Vec<(u32, i32)> = (0..n)
                .map(|i| {
                    s = (s * 279_470_273) % 4_294_967_291;
                    ((s % range) as u32, i)
                })
                .collect();
            let line = want.next().unwrap();
            let (head, order) = line.split_once(':').unwrap();
            assert_eq!(head, format!("n {n} range {range}"));
            let order: Vec<i32> = order
                .split_whitespace()
                .map(|x| x.parse().unwrap())
                .collect();
            assert_eq!(sorted(keys), order, "{head}");
            checked += 1;
        }
    }
    let pattern: Vec<(u32, i32)> = (0..3000)
        .map(|i| {
            (
                if i % 7 == 0 {
                    5
                } else {
                    ((i * 7919) % 13) as u32
                },
                i,
            )
        })
        .collect();
    let mut descending: Vec<(u32, i32)> = (0..4000).map(|i| ((4000 - i) as u32, i)).collect();
    for p in descending.iter_mut().step_by(2) {
        p.0 = 7;
    }
    for (name, keys) in [
        ("n 3000 pattern", pattern),
        ("n 4000 descending", descending),
    ] {
        let line = want.next().unwrap();
        let (head, order) = line.split_once(':').unwrap();
        assert_eq!(head, name);
        let order: Vec<i32> = order
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        assert_eq!(sorted(keys), order, "{head}");
        checked += 1;
    }
    assert_eq!(checked, 29);
    assert!(want.next().is_none());
}

fn reference_counts(seed: u64) -> Vec<(usize, usize, usize)> {
    let text = std::fs::read_to_string(data("folds_ref_percolator_v3.txt")).unwrap();
    let prefix = format!("seed {seed}\t");
    let mut init = Vec::new();
    let mut training = Vec::new();
    for line in text.lines().filter_map(|l| l.strip_prefix(&prefix)) {
        let words: Vec<&str> = line.split_whitespace().collect();
        if let Some(k) = words.iter().position(|&w| w == "separate") {
            let feature: usize = words[words.iter().position(|&w| w == "feature").unwrap() + 1]
                .parse()
                .unwrap();
            init.push((feature, words[k + 1].parse::<usize>().unwrap()));
        }
        if let Some(k) = words.iter().position(|&w| w == "Training") {
            training.push((
                words[k + 2].parse::<usize>().unwrap(),
                words[k + 5].parse::<usize>().unwrap(),
            ));
        }
    }
    assert_eq!(init.len(), 3);
    assert_eq!(training.len(), 3);
    init.iter()
        .zip(&training)
        .map(|(&(feature, p), &(tp, neg))| {
            assert_eq!(p, tp);
            (feature, p, neg)
        })
        .collect()
}

#[test]
fn folds_reproduce_percolator_training_sets() {
    let d = tempfile::tempdir().unwrap();
    let pin = read_pin(
        &data("folds_ref.pin"),
        MemLimit {
            bytes: u64::MAX,
            spill_dir: d.path().to_path_buf(),
        },
    )
    .unwrap();
    let norm = stdv(&pin);
    for seed in [42u64, 7, 1] {
        let fold = folds_by_spectrum(&pin.meta, 3, seed);
        let got: Vec<(usize, usize, usize)> = (0..3u8)
            .map(|f| {
                let train: Vec<usize> = (0..pin.meta.len()).filter(|&i| fold[i] != f).collect();
                let decoy: Vec<bool> = train.iter().map(|&i| pin.meta[i].label < 0).collect();
                let w = init_direction(&pin, &norm, &train, 0.01, false).unwrap();
                let feature = w.iter().position(|&v| v != 0.0).unwrap() + 1;
                let s = score(&pin.rows, &train, &w, Some(&norm));
                let p = count_at(&tdc_q(&s, &decoy, true), &decoy, 0.01);
                (feature, p, decoy.iter().filter(|&&x| x).count())
            })
            .collect();
        assert_eq!(got, reference_counts(seed), "seed {seed}");
    }
}
