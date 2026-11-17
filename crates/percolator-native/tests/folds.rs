use percolator_native::{
    count_at, folds_by_spectrum, init_direction, read_pin, score, stdv, tdc_q, MemLimit,
};
use std::path::PathBuf;

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
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

type Split = (usize, usize, usize);

fn collision_reference() -> Vec<(String, Vec<Split>)> {
    std::fs::read_to_string(data("fold_collision_percolator_v3.txt"))
        .unwrap()
        .lines()
        .skip(1)
        .map(|l| {
            let mut cells = l.split('\t');
            let threads = cells.next().unwrap().to_string();
            let splits = cells
                .map(|c| {
                    let v: Vec<usize> = c.split(' ').map(|x| x.parse().unwrap()).collect();
                    (v[0], v[1], v[2])
                })
                .collect();
            (threads, splits)
        })
        .collect()
}

#[test]
fn colliding_scan_hashes_split_as_percolator_at_one_thread() {
    let d = tempfile::tempdir().unwrap();
    let pin = read_pin(
        &data("fold_collision.pin"),
        MemLimit {
            bytes: u64::MAX,
            spill_dir: d.path().to_path_buf(),
        },
    )
    .unwrap();
    let norm = stdv(&pin);
    let fold = folds_by_spectrum(&pin.meta, 3, 42);
    let got: Vec<Split> = (0..3u8)
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
    let reference = collision_reference();
    assert_eq!(reference[0].0, "threads 1");
    assert_eq!(got, reference[0].1);
    assert!(reference[1..].iter().any(|(_, r)| *r != got));
}
