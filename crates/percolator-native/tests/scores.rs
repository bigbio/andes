use percolator_native::{
    count_at, folds_by_spectrum, raw_weights, read_pin, score, stdv, tdc_q, tdc_winners,
    training_sets, MemLimit, Pin, RowMeta,
};
use std::collections::HashMap;

fn pin_from(rows: &[(i8, u32, u64, f64, [f32; 3])]) -> (tempfile::TempDir, Pin) {
    let d = tempfile::tempdir().unwrap();
    let mut text =
        String::from("SpecId\tLabel\tScanNr\tExpMass\tFileName\tf1\tf2\tf3\tPeptide\tProteins\n");
    for (i, (label, file, scan, mass, f)) in rows.iter().enumerate() {
        text.push_str(&format!(
            "r{i}\t{label}\t{scan}\t{mass}\tfile{file}\t{}\t{}\t{}\tK.A.R\tP\n",
            f[0], f[1], f[2]
        ));
    }
    let p = d.path().join("x.pin");
    std::fs::write(&p, text).unwrap();
    let pin = read_pin(
        &p,
        MemLimit {
            bytes: u64::MAX,
            spill_dir: d.path().to_path_buf(),
        },
    )
    .unwrap();
    (d, pin)
}

fn matrix() -> (tempfile::TempDir, Pin) {
    pin_from(&[
        (1, 0, 1, 1.0, [1.0, 5.0, 0.0]),
        (1, 0, 2, 1.0, [2.0, 5.0, 0.0]),
        (-1, 0, 3, 1.0, [3.0, 5.0, 0.0]),
        (-1, 0, 4, 1.0, [4.0, 5.0, 4.0]),
    ])
}

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-9
}

#[test]
fn stdv_matches_reference() {
    let (_d, pin) = matrix();
    let n = stdv(&pin);
    assert!(close(n.mean[0], 2.5));
    assert!(close(n.mean[1], 5.0));
    assert!(close(n.mean[2], 1.0));
    assert!(close(1.0 / n.inv_sd[0], 1.25f64.sqrt()));
    assert!(close(1.0 / n.inv_sd[2], 3.0f64.sqrt()));
}

#[test]
fn stdv_zero_variance_feature_gives_finite_weights() {
    let (_d, pin) = matrix();
    let n = stdv(&pin);
    assert_eq!(n.inv_sd[1], 1.0);
    let raw = raw_weights(&n, &[1.0, 2.0, 3.0, 0.5]);
    assert!(raw.iter().all(|x| x.is_finite()));
    let s = score(&pin.rows, &[0, 1, 2, 3], &[1.0, 2.0, 3.0, 0.5], Some(&n));
    assert!(s.iter().all(|x| x.is_finite()));
}

#[test]
fn raw_weights_roundtrip() {
    let (_d, pin) = matrix();
    let n = stdv(&pin);
    let w = [0.7, -1.3, 2.1, 0.4];
    let raw = raw_weights(&n, &w);
    let idx = [3, 0, 2, 1, 1];
    let a = score(&pin.rows, &idx, &w, Some(&n));
    let b = score(&pin.rows, &idx, &raw, None);
    for (x, y) in a.iter().zip(&b) {
        assert!(close(*x, *y), "{x} {y}");
    }
}

const SCORES: [f64; 10] = [10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0];
const DECOY: [bool; 10] = [
    false, false, true, false, false, true, false, true, false, false,
];

#[test]
fn tdc_q_skip_plus_one() {
    let skip = tdc_q(&SCORES, &DECOY, true);
    let t = 3.0 / 7.0;
    let want = [0.0, 0.0, 0.25, 0.25, 0.25, 0.4, 0.4, t, t, t];
    for (g, w) in skip.iter().zip(want) {
        assert!(close(*g, w), "{skip:?}");
    }
    let plus = tdc_q(&SCORES, &DECOY, false);
    let u = 4.0 / 7.0;
    let want = [0.5, 0.5, 0.5, 0.5, 0.5, u, u, u, u, u];
    for (g, w) in plus.iter().zip(want) {
        assert!(close(*g, w), "{plus:?}");
    }
    let shuffled_scores = [1.0, 10.0, 5.0, 9.0];
    let shuffled_decoy = [false, false, true, false];
    let q = tdc_q(&shuffled_scores, &shuffled_decoy, true);
    assert!(close(q[1], 0.0) && close(q[3], 0.0));
    assert!(close(q[2], 1.0 / 3.0) && close(q[0], 1.0 / 3.0));
    assert_eq!(count_at(&q, &shuffled_decoy, 0.5), 3);
    assert_eq!(count_at(&q, &shuffled_decoy, 1.0 / 3.0), 2);
}

#[test]
fn tdc_q_ties_are_one_group() {
    let q = tdc_q(&[1.0, 1.0, 0.0], &[false, true, false], true);
    assert!(q.iter().all(|&x| close(x, 0.5)), "{q:?}");
}

#[test]
fn training_sets_select_targets_at_q_and_all_decoys() {
    let q = tdc_q(&SCORES, &DECOY, true);
    let (pos, neg) = training_sets(&q, &DECOY, 0.25);
    assert_eq!(pos, vec![0, 1, 3, 4]);
    assert_eq!(neg, vec![2, 5, 7]);
}

fn meta(file: u32, scan: u64, mass: f64) -> RowMeta {
    RowMeta {
        label: 1,
        file,
        scan,
        exp_mass: mass,
    }
}

#[test]
fn folds_keep_spectrum_together() {
    let mut m = Vec::new();
    for file in 0..2u32 {
        for scan in 0..500u64 {
            for _ in 0..(1 + scan % 3) {
                m.push(meta(file, scan, 100.0));
            }
        }
    }
    let f = folds_by_spectrum(&m, 3, 17);
    assert_eq!(f, folds_by_spectrum(&m, 3, 17));
    assert_ne!(f, folds_by_spectrum(&m, 3, 18));
    let mut by: HashMap<(u32, u64), u8> = HashMap::new();
    for (r, &fold) in m.iter().zip(&f) {
        assert!(fold < 3);
        assert_eq!(*by.entry((r.file, r.scan)).or_insert(fold), fold);
    }
    let mut sizes = [0usize; 3];
    for &fold in &f {
        sizes[fold as usize] += 1;
    }
    let (lo, hi) = (sizes.iter().min().unwrap(), sizes.iter().max().unwrap());
    assert!(hi - lo <= 3, "{sizes:?}");
}

fn id(i: usize) -> Result<String, String> {
    Ok(format!("r{i}"))
}

#[test]
fn tdc_winner_per_file_scan_mass() {
    let m = vec![
        meta(0, 1, 100.0),
        meta(0, 1, 100.0),
        meta(0, 1, 200.0),
        meta(1, 1, 100.0),
        meta(0, 2, 100.0),
        meta(0, 2, 100.0),
    ];
    let s = [1.0, 3.0, 0.5, 0.1, 2.0, 2.0];
    assert_eq!(tdc_winners(&m, &s, id).unwrap(), vec![1, 2, 3, 4]);
}

#[test]
fn tdc_winner_tie_prefers_decoy_then_smaller_id() {
    let mut m = vec![
        meta(0, 1, 100.0),
        meta(0, 1, 100.0),
        meta(0, 2, 100.0),
        meta(0, 2, 100.0),
    ];
    m[0].label = -1;
    m[3].label = -1;
    let s = [2.0, 2.0, 1.0, 1.0];
    assert_eq!(tdc_winners(&m, &s, id).unwrap(), vec![0, 3]);
    let m = vec![meta(0, 1, 100.0), meta(0, 1, 100.0), meta(0, 1, 100.0)];
    let names = ["r2", "r10", "r10"];
    let by_name = |i: usize| Ok(names[i].to_string());
    assert_eq!(tdc_winners(&m, &[1.0, 1.0, 1.0], by_name).unwrap(), vec![1]);
}
