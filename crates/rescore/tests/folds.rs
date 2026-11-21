use rescore::{
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
            sizing: rescore::Sizing::Fixed,
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
fn colliding_scan_hashes_split_as_percolator_at_any_thread_count() {
    let d = tempfile::tempdir().unwrap();
    let pin = read_pin(
        &data("fold_collision.pin"),
        MemLimit {
            bytes: u64::MAX,
            spill_dir: d.path().to_path_buf(),
            sizing: rescore::Sizing::Fixed,
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
    assert_eq!(reference.len(), 4);
    for (threads, splits) in &reference {
        assert_eq!(&got, splits, "{threads}");
    }
}

#[test]
fn seeds_reduce_modulo_the_generator_and_zero_acts_as_one() {
    let d = tempfile::tempdir().unwrap();
    let pin = read_pin(
        &data("folds_ref.pin"),
        MemLimit {
            bytes: u64::MAX,
            spill_dir: d.path().to_path_buf(),
            sizing: rescore::Sizing::Fixed,
        },
    )
    .unwrap();
    let one = folds_by_spectrum(&pin.meta, 3, 1);
    let m = 4_294_967_291u64;
    for seed in [0, m, 2 * m, m + 1, u64::MAX - u64::MAX % m] {
        assert_eq!(folds_by_spectrum(&pin.meta, 3, seed), one, "seed {seed}");
    }
    assert_eq!(
        folds_by_spectrum(&pin.meta, 3, m + 42),
        folds_by_spectrum(&pin.meta, 3, 42)
    );
    assert_eq!(
        folds_by_spectrum(&pin.meta, 3, u64::MAX),
        folds_by_spectrum(&pin.meta, 3, u64::MAX % m)
    );
}

fn letters(mut i: u64) -> String {
    const ALPHABET: &[u8] = b"ACDEFGHKMNPQRSTVWY";
    let mut s = String::new();
    loop {
        s.push(ALPHABET[(i % 18) as usize] as char);
        i /= 18;
        if i == 0 {
            return s;
        }
    }
}

fn row(file: u32, scan: u64, label: i8, peptide: &str) -> rescore::RowMeta {
    rescore::RowMeta {
        label,
        file,
        scan,
        exp_mass: 1000.0,
        sequence: rescore::sequence_key(peptide),
    }
}

#[test]
fn sequence_key_ignores_flanks_modifications_and_isoleucine() {
    use rescore::sequence_key;
    let base = sequence_key("K.PEPTLDEM.R");
    assert_eq!(sequence_key("PEPTLDEM"), base);
    assert_eq!(sequence_key("-.PEPTIDEM[15.9949].-"), base);
    assert_eq!(sequence_key("n[42.0106]PEPTIDEM(UNIMOD:35)"), base);
    assert_eq!(sequence_key("PEPTIDEM+15.995"), base);
    assert_ne!(sequence_key("PEPTIDEK"), base);
}

#[test]
fn spectra_with_the_same_best_peptide_share_a_fold() {
    let mut meta = Vec::new();
    for s in 0..600u64 {
        let best = format!("K.PEPTWDE{}K.R", letters(s / 6));
        meta.push(row((s % 2) as u32, s, 1, &best));
        meta.push(row((s % 2) as u32, s, -1, &format!("K.DCQY{}K.R", letters(s))));
        meta.push(row((s % 2) as u32, s, 1, &format!("K.QTHER{}K.R", letters(s % 7))));
    }
    for seed in [1u64, 42, 7] {
        let fold = rescore::folds_by_precursor(&meta, 3, seed);
        for s in 0..600usize {
            let f = fold[3 * s];
            assert_eq!(fold[3 * s + 1], f, "rows of one spectrum split");
            assert_eq!(fold[3 * s + 2], f, "rows of one spectrum split");
            let first = 3 * (s / 6 * 6);
            assert_eq!(fold[first], f, "spectra of one best peptide split");
        }
        let sizes: Vec<usize> = (0..3u8).map(|f| fold.iter().filter(|&&x| x == f).count()).collect();
        assert!(sizes.iter().all(|&c| c > 400), "unbalanced folds {sizes:?}");
    }
}

#[test]
fn precursor_folds_fill_every_fold_with_few_units() {
    let meta: Vec<_> = (0..3u64).map(|s| row(0, s, 1, &format!("PEPTWDE{}", letters(s)))).collect();
    for seed in [1u64, 2, 3, 42] {
        let fold = rescore::folds_by_precursor(&meta, 3, seed);
        let mut seen = [false; 3];
        for &f in &fold {
            seen[f as usize] = true;
        }
        assert!(seen.iter().all(|&s| s), "seed {seed}: {fold:?}");
    }
    let one = vec![row(0, 0, 1, "PEPTIDE"), row(0, 1, -1, "PEPTIDE")];
    assert_eq!(rescore::folds_by_precursor(&one, 3, 1).len(), 2);
    assert!(rescore::folds_by_precursor(&[], 3, 1).is_empty());
}
