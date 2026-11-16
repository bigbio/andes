use percolator_native::{
    estimate_pi0, pep_is_degenerate, pep_mix_max, pep_tdc, pep_tdc_fit, tdc_q, too_good_separation,
};
use std::path::PathBuf;

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

fn numbers(name: &str) -> Vec<f64> {
    std::fs::read_to_string(data(name))
        .unwrap()
        .split_whitespace()
        .map(|v| v.parse().unwrap())
        .collect()
}

fn fixture(far: Option<f64>) -> (Vec<f64>, Vec<bool>) {
    let t = numbers("pep_ref/t.txt");
    let mut d = numbers("pep_ref/d.txt");
    if let Some(v) = far {
        let last = d.len() - 1;
        d[last] = v;
    }
    let decoy: Vec<bool> = t
        .iter()
        .map(|_| false)
        .chain(d.iter().map(|_| true))
        .collect();
    (t.into_iter().chain(d).collect(), decoy)
}

fn qvality(name: &str) -> Vec<(f64, f64)> {
    std::fs::read_to_string(data(name))
        .unwrap()
        .lines()
        .skip(1)
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            (c[0].parse().unwrap(), c[1].parse().unwrap())
        })
        .collect()
}

fn by_score_desc(scores: &[f64], pep: &[f64]) -> Vec<(f64, f64)> {
    let mut v: Vec<(f64, f64)> = scores.iter().copied().zip(pep.iter().copied()).collect();
    v.sort_by(|a, b| b.0.total_cmp(&a.0));
    v
}

fn max_rel_diff(ours: &[(f64, f64)], reference: &[(f64, f64)]) -> (f64, usize) {
    assert_eq!(ours.len(), reference.len());
    let mut worst = (0.0f64, 0usize);
    for (i, (a, b)) in ours.iter().zip(reference).enumerate() {
        assert!(
            (a.0 - b.0).abs() < 1e-9,
            "row {i}: score {} vs {}",
            a.0,
            b.0
        );
        let rel = (a.1 - b.1).abs() / b.1.abs().max(1e-300);
        if rel > worst.0 {
            worst = (rel, i);
        }
    }
    worst
}

fn assert_monotone_bounded(scores: &[f64], pep: &[f64]) {
    assert_eq!(scores.len(), pep.len());
    let v = by_score_desc(scores, pep);
    for (i, w) in v.windows(2).enumerate() {
        assert!(w[0].1 <= w[1].1, "row {i}: PEP {} above {}", w[0].1, w[1].1);
        if w[0].0 == w[1].0 {
            assert_eq!(w[0].1, w[1].1);
        }
    }
    for &p in pep {
        assert!(
            p.is_finite() && (0.0..=1.0).contains(&p),
            "PEP {p} outside [0, 1]"
        );
    }
}

#[test]
fn pep_matches_qvality_tdc() {
    let (s, d) = fixture(None);
    let pep = pep_tdc_fit(&s, &d).expect("fit");
    let reference = qvality("pep_ref/q_d.txt");
    let (worst, at) = max_rel_diff(&by_score_desc(&s, &pep), &reference);
    eprintln!("qvality -Y: max relative PEP difference {worst:e} at row {at}");
    assert!(
        worst < 1e-4,
        "max relative PEP difference {worst:e} at row {at}"
    );
    assert_eq!(pep, pep_tdc(&s, &d));
}

#[test]
fn estimate_pi0_matches_qvality() {
    let (s, d) = fixture(None);
    let pi0 = estimate_pi0(&s, &d);
    eprintln!("pi0 {pi0}");
    assert!((pi0 - 0.409259).abs() < 5e-7, "pi0 {pi0}");
}

#[test]
fn pep_matches_qvality_mix_max() {
    let (s, d) = fixture(None);
    let pep = pep_mix_max(&s, &d, estimate_pi0(&s, &d));
    let reference = qvality("pep_ref/qm_d.txt");
    let (worst, at) = max_rel_diff(&by_score_desc(&s, &pep), &reference);
    eprintln!("qvality mix-max: max relative PEP difference {worst:e} at row {at}");
    assert!(
        worst < 1e-4,
        "max relative PEP difference {worst:e} at row {at}"
    );
}

#[test]
fn pep_matches_percolator_svm_ref() {
    let mut rows: Vec<(f64, f64)> = Vec::new();
    let mut scores = Vec::new();
    let mut decoy = Vec::new();
    for (name, is_decoy) in [
        ("svm_ref_perc.target.psms", false),
        ("svm_ref_perc.decoy.psms", true),
    ] {
        for l in std::fs::read_to_string(data(name)).unwrap().lines().skip(1) {
            let c: Vec<&str> = l.split('\t').collect();
            let (s, p): (f64, f64) = (c[1].parse().unwrap(), c[3].parse().unwrap());
            rows.push((s, p));
            scores.push(s);
            decoy.push(is_decoy);
        }
    }
    let pep = pep_tdc(&scores, &decoy);
    let worst = rows
        .iter()
        .zip(&pep)
        .map(|((_, want), got)| (got.ln() - want.ln()).abs())
        .fold(0.0f64, f64::max);
    eprintln!("svm_ref: max |ln PEP - ln Percolator PEP| {worst:e}");
    for (i, ((s, want), got)) in rows.iter().zip(&pep).enumerate() {
        let close = (got - want).abs() <= 1e-6 || ((got.ln() - want.ln()).abs() < 0.05);
        assert!(close, "row {i} score {s}: PEP {got:e}, Percolator {want:e}");
    }
}

#[test]
fn pep_monotone_and_bounded() {
    for far in [None, Some(8.0), Some(12.0)] {
        let (s, d) = fixture(far);
        assert_monotone_bounded(&s, &pep_tdc(&s, &d));
        assert_monotone_bounded(&s, &pep_mix_max(&s, &d, estimate_pi0(&s, &d)));
    }
}

#[test]
fn pep_never_degenerate_on_far_decoy() {
    for far in [8.0, 10.0, 12.0] {
        let (s, d) = fixture(Some(far));
        let q = tdc_q(&s, &d, false);
        let raw = pep_tdc_fit(&s, &d).expect("fit");
        assert!(
            pep_is_degenerate(&raw, &q, &d),
            "far decoy {far}: the 3.7.1 fit should give PEP 1 to every accepted PSM"
        );
        let pep = pep_tdc(&s, &d);
        assert!(!pep_is_degenerate(&pep, &q, &d), "far decoy {far}");
        let mut accepted: Vec<f64> = pep
            .iter()
            .zip(&q)
            .zip(&d)
            .filter(|((_, &qi), &di)| !di && qi <= 0.01)
            .map(|((&p, _), _)| p)
            .collect();
        accepted.sort_by(f64::total_cmp);
        assert!(accepted.len() > 1000);
        assert!(
            accepted[accepted.len() / 2] < 0.05,
            "median accepted PEP {}",
            accepted[accepted.len() / 2]
        );
        let mm = pep_mix_max(&s, &d, estimate_pi0(&s, &d));
        assert!(!pep_is_degenerate(&mm, &q, &d), "mix-max, far decoy {far}");
    }
}

#[test]
fn healthy_fit_is_returned_unchanged() {
    for far in [None, Some(5.0)] {
        let (s, d) = fixture(far);
        assert_eq!(pep_tdc(&s, &d), pep_tdc_fit(&s, &d).expect("fit"));
    }
}

#[test]
fn too_few_bins_falls_back_to_bounded_pep() {
    let s = [3.0, 2.0, 1.0];
    let d = [false, false, true];
    assert!(pep_tdc_fit(&s, &d).is_none());
    assert_monotone_bounded(&s, &pep_tdc(&s, &d));
}

#[test]
fn separation_check_as_check_separation() {
    let (s, d) = fixture(None);
    assert!(!too_good_separation(&s, &d));
    let s: Vec<f64> = (0..600).map(|i| i as f64).collect();
    let d: Vec<bool> = (0..600).map(|i| i < 300).collect();
    assert!(too_good_separation(&s, &d));
    let d: Vec<bool> = (0..600).map(|i| i < 100).collect();
    assert!(!too_good_separation(&s, &d));
}
