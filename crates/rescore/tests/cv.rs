use rescore::{
    count_at, init_direction, l2_svm_mfn, read_pin, score, stdv, tdc_q, train, CvConfig, CvResult,
    MemLimit, MfnOpts, Pin, Rows,
};

fn lcg(state: &mut u64) -> f64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*state >> 11) as f64 + 0.5) / (1u64 << 53) as f64
}

fn noise(state: &mut u64) -> f64 {
    lcg(state) + lcg(state) + lcg(state) - 1.5
}

fn letters(mut i: usize) -> String {
    const ALPHABET: &[u8] = b"ACDEFGHKMNPQRSTVWY";
    let mut s = String::new();
    loop {
        s.push(ALPHABET[i % 18] as char);
        i /= 18;
        if i == 0 {
            return s;
        }
    }
}

fn pin_text(rows: &[(i8, u64, Vec<f64>)], default: Option<&[f64]>) -> String {
    let nf = rows[0].2.len();
    let names: Vec<String> = (1..=nf).map(|j| format!("f{j}")).collect();
    let mut text = format!(
        "SpecId\tLabel\tScanNr\tExpMass\tFileName\t{}\tPeptide\tProteins\n",
        names.join("\t")
    );
    if let Some(d) = default {
        let cells: Vec<String> = d.iter().map(|v| v.to_string()).collect();
        text.push_str(&format!(
            "DefaultDirection\t-\t-\t-\t-\t{}\t-\t-\n",
            cells.join("\t")
        ));
    }
    for (i, (label, scan, f)) in rows.iter().enumerate() {
        let cells: Vec<String> = f.iter().map(|v| format!("{v:.6}")).collect();
        text.push_str(&format!(
            "r{i}\t{label}\t{scan}\t1000.0\tfile0\t{}\tK.PEPTWDE{}.R\tP{i}\n",
            cells.join("\t"),
            letters(i)
        ));
    }
    text
}

fn load(rows: &[(i8, u64, Vec<f64>)], default: Option<&[f64]>) -> (tempfile::TempDir, Pin) {
    load_within(rows, default, u64::MAX)
}

fn load_within(
    rows: &[(i8, u64, Vec<f64>)],
    default: Option<&[f64]>,
    bytes: u64,
) -> (tempfile::TempDir, Pin) {
    let d = tempfile::tempdir().unwrap();
    let p = d.path().join("x.pin");
    std::fs::write(&p, pin_text(rows, default)).unwrap();
    let pin = read_pin(
        &p,
        MemLimit {
            bytes,
            spill_dir: d.path().to_path_buf(),
            sizing: rescore::Sizing::Fixed,
        },
    )
    .unwrap();
    (d, pin)
}

fn synthetic(n: usize, shifts: &[f64], seed: u64) -> Vec<(i8, u64, Vec<f64>)> {
    let mut s = seed;
    (0..n)
        .map(|i| {
            let target = lcg(&mut s) < 0.5;
            let f = shifts
                .iter()
                .map(|&sh| noise(&mut s) + if target { sh } else { 0.0 })
                .collect();
            (if target { 1 } else { -1 }, i as u64 + 1, f)
        })
        .collect()
}

fn pool(threads: usize) -> rayon::ThreadPool {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .unwrap()
}

fn rows_of(res: &CvResult, f: usize, train_rows: bool) -> Vec<rescore::RowId> {
    (0..res.fold.len() as rescore::RowId)
        .filter(|&i| (res.fold[i as usize] as usize == f) != train_rows)
        .collect()
}

fn decoys(pin: &Pin, rows: &[rescore::RowId]) -> Vec<bool> {
    rows.iter().map(|&i| pin.meta[i as usize].label < 0).collect()
}

fn bits(v: &[f64]) -> Vec<u64> {
    v.iter().map(|x| x.to_bits()).collect()
}

#[test]
fn cv_deterministic_across_threads() {
    let (_d, pin) = load(&synthetic(900, &[0.9, 0.6, 0.0, 0.3], 11), None);
    let norm = stdv(&pin);
    let cfg = CvConfig::default();
    let a = train(&pin, &norm, &cfg, &pool(1)).unwrap();
    let b = train(&pin, &norm, &cfg, &pool(8)).unwrap();
    assert_eq!(a.fold, b.fold);
    for (x, y) in a
        .w_norm
        .iter()
        .chain(&a.w_raw)
        .zip(b.w_norm.iter().chain(&b.w_raw))
    {
        assert_eq!(bits(x), bits(y));
    }
    assert_eq!(bits(&a.test_scores), bits(&b.test_scores));
    assert!(a.test_scores.iter().all(|s| s.is_finite()));
}

#[test]
fn cv_grid_tie_picks_later_pair() {
    let (_d, pin) = load(&synthetic(600, &[12.0, 0.0, 0.0], 5), None);
    let norm = stdv(&pin);
    let p = pool(4);
    let cfg = |cpos: Vec<f64>| CvConfig {
        iterations: 1,
        cpos,
        cfrac: vec![1.0],
        ..CvConfig::default()
    };
    let both = train(&pin, &norm, &cfg(vec![10.0, 0.1]), &p).unwrap();
    let early = train(&pin, &norm, &cfg(vec![10.0]), &p).unwrap();
    let late = train(&pin, &norm, &cfg(vec![0.1]), &p).unwrap();
    for f in 0..3 {
        let tr = rows_of(&both, f, true);
        let dec = decoys(&pin, &tr);
        let targets = dec.iter().filter(|&&d| !d).count();
        for r in [&early, &late] {
            let s = score(&pin.rows, &tr, &r.w_norm[f], Some(&norm));
            assert_eq!(count_at(&tdc_q(&s, &dec, true), &dec, 0.01), targets);
        }
        assert_ne!(bits(&early.w_norm[f]), bits(&late.w_norm[f]));
        assert_eq!(bits(&both.w_norm[f]), bits(&late.w_norm[f]));
    }
}

#[test]
fn cv_grid_with_no_budget_equals_unlimited() {
    let rows = synthetic(900, &[0.9, 0.6, 0.0, 0.3], 17);
    let (_a, big) = load(&rows, None);
    let (_b, small) = load_within(&rows, None, 1);
    assert!(small.rows.is_spilled());
    let cfg = CvConfig {
        iterations: 4,
        ..CvConfig::default()
    };
    let p = pool(8);
    let a = train(&big, &stdv(&big), &cfg, &p).unwrap();
    let b = train(&small, &stdv(&small), &cfg, &p).unwrap();
    assert_eq!(a.fold, b.fold);
    for (x, y) in a
        .w_norm
        .iter()
        .chain(&a.w_raw)
        .zip(b.w_norm.iter().chain(&b.w_raw))
    {
        assert_eq!(bits(x), bits(y));
    }
    assert_eq!(bits(&a.test_scores), bits(&b.test_scores));
}

fn twins(n: usize) -> Vec<(i8, u64, Vec<f64>)> {
    let mut s = 3u64;
    let mut rows = Vec::new();
    for i in 0..n {
        let f: Vec<f64> = (0..3).map(|_| noise(&mut s)).collect();
        rows.push((1, i as u64 + 1, f.clone()));
        rows.push((-1, i as u64 + 1, f));
    }
    rows
}

fn all_rows(pin: &Pin) -> Vec<rescore::RowId> {
    (0..pin.meta.len() as rescore::RowId).collect()
}

#[test]
fn pin_direction_has_priority_and_overrule_keeps_each_fold_trained() {
    let p = pool(2);
    let (_d, pin) = load(&twins(200), Some(&[0.0, 1.5, -2.0]));
    let norm = stdv(&pin);
    assert_eq!(
        init_direction(&pin, &norm, &all_rows(&pin), 0.01, false).unwrap(),
        vec![0.0, 1.5, -2.0, 0.0]
    );
    let res = train(&pin, &norm, &CvConfig::default(), &p).unwrap();
    for (wn, wr) in res.w_norm.iter().zip(&res.w_raw) {
        assert_eq!(*wn, vec![0.0, 1.5, -2.0, 0.0]);
        assert!(wr.iter().all(|v| v.is_finite()));
    }
    assert!(res.test_scores.iter().all(|s| s.is_finite()));

    let (_d2, pin) = load(&twins(200), None);
    let norm = stdv(&pin);
    assert_eq!(
        init_direction(&pin, &norm, &all_rows(&pin), 0.01, true).unwrap(),
        vec![-1.0, 0.0, 0.0, 0.0]
    );
    let cfg = CvConfig {
        overrule: true,
        ..CvConfig::default()
    };
    let res = train(&pin, &norm, &cfg, &p).unwrap();
    for wn in &res.w_norm {
        assert!(wn.iter().all(|v| v.is_finite()));
    }
    for i in 0..pin.meta.len() {
        let s = score(
            &pin.rows,
            &[i as rescore::RowId],
            &res.w_norm[res.fold[i] as usize],
            Some(&norm),
        )[0];
        assert_eq!(res.test_scores[i].to_bits(), s.to_bits());
    }
}

#[test]
fn no_separating_feature_is_an_error() {
    let (_d, pin) = load(&twins(200), None);
    let norm = stdv(&pin);
    let e = init_direction(&pin, &norm, &all_rows(&pin), 0.01, false).unwrap_err();
    assert!(e.contains("cannot find an initial direction"), "{e}");
    let Err(e) = train(&pin, &norm, &CvConfig::default(), &pool(2)) else {
        panic!("expected an error");
    };
    assert!(e.contains("cannot find an initial direction"), "{e}");
}

#[test]
fn zero_decoys_or_targets_is_an_error() {
    let rows = synthetic(300, &[2.0, 0.0, 0.0], 4);
    for (label, word) in [(1i8, "decoy"), (-1i8, "target")] {
        let one: Vec<(i8, u64, Vec<f64>)> = rows
            .iter()
            .map(|(_, scan, f)| (label, *scan, f.clone()))
            .collect();
        let (_d, pin) = load(&one, None);
        let norm = stdv(&pin);
        let Err(e) = train(&pin, &norm, &CvConfig::default(), &pool(2)) else {
            panic!("expected an error");
        };
        assert_eq!(e, format!("Error: no {word} PSMs were provided."));
    }
}

#[test]
fn init_direction_picks_best_feature_with_sign() {
    let (_d, pin) = load(&synthetic(600, &[0.0, -3.0, 1.0], 9), None);
    let norm = stdv(&pin);
    let all: Vec<rescore::RowId> = (0..pin.meta.len() as rescore::RowId).collect();
    assert_eq!(
        init_direction(&pin, &norm, &all, 0.01, false).unwrap(),
        vec![0.0, -1.0, 0.0, 0.0]
    );
}

#[test]
fn cv_improves_over_init() {
    let (_d, pin) = load(&synthetic(1200, &[2.5, 0.0, 0.0], 21), None);
    let norm = stdv(&pin);
    let cfg = CvConfig::default();
    let res = train(&pin, &norm, &cfg, &pool(4)).unwrap();
    let (init_pos, trained_pos, _) = init_and_trained(&pin, &res, &cfg);
    assert!(init_pos > 0);
    assert!(trained_pos >= init_pos, "{trained_pos} < {init_pos}");
}

#[test]
fn cv_training_combines_features() {
    let (_d, pin) = load(&synthetic(3000, &[1.2, 1.2, 0.0], 23), None);
    let norm = stdv(&pin);
    let cfg = CvConfig::default();
    let res = train(&pin, &norm, &cfg, &pool(4)).unwrap();
    let (init_pos, trained_pos, moved) = init_and_trained(&pin, &res, &cfg);
    assert!(trained_pos > init_pos, "{trained_pos} <= {init_pos}");
    assert!(moved);
}

fn init_and_trained(pin: &Pin, res: &CvResult, cfg: &CvConfig) -> (usize, usize, bool) {
    let norm = stdv(pin);
    let mut init_pos = 0;
    let mut trained_pos = 0;
    let mut moved = false;
    for f in 0..3 {
        let tr = rows_of(res, f, true);
        let te = rows_of(res, f, false);
        let dec = decoys(pin, &te);
        let w0 = init_direction(pin, &norm, &tr, cfg.initial_selection_fdr, false).unwrap();
        let s0 = score(&pin.rows, &te, &w0, Some(&norm));
        init_pos += count_at(&tdc_q(&s0, &dec, false), &dec, cfg.test_fdr);
        let s1: Vec<f64> = te.iter().map(|&i| res.test_scores[i as usize]).collect();
        trained_pos += count_at(&tdc_q(&s1, &dec, false), &dec, cfg.test_fdr);
        moved |= res.w_norm[f][..3] != w0[..3];
    }
    (init_pos, trained_pos, moved)
}

#[test]
fn mfn_finite_with_one_class() {
    let mut s = 17u64;
    let mut data = Vec::new();
    for _ in 0..40 {
        data.push(noise(&mut s));
        data.push(noise(&mut s));
        data.push(1.0);
    }
    let x = Rows::from_dense(data, 40, 3);
    let all: Vec<usize> = (0..40).collect();
    let neg_only = l2_svm_mfn(&x, &[], &all, 1.0, 1.0, &MfnOpts::default());
    let pos_only = l2_svm_mfn(&x, &all, &[], 1.0, 1.0, &MfnOpts::default());
    assert!(neg_only.iter().all(|v| v.is_finite()), "{neg_only:?}");
    assert!(pos_only.iter().all(|v| v.is_finite()), "{pos_only:?}");
    assert!(neg_only[2] < 0.0 && pos_only[2] > 0.0);
}
