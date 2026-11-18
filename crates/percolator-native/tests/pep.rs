use percolator_native::{estimate_pi0, mix_max_q, pep_from_q, tdc_q, too_good_separation};
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

fn fixture() -> (Vec<f64>, Vec<bool>) {
    let t = numbers("pep_ref/t.txt");
    let d = numbers("pep_ref/d.txt");
    let decoy: Vec<bool> = t
        .iter()
        .map(|_| false)
        .chain(d.iter().map(|_| true))
        .collect();
    (t.into_iter().chain(d).collect(), decoy)
}

fn order(scores: &[f64], decoy: &[bool]) -> Vec<usize> {
    let mut o: Vec<usize> = (0..scores.len()).collect();
    o.sort_by(|&a, &b| {
        scores[b]
            .total_cmp(&scores[a])
            .then(decoy[a].cmp(&decoy[b]))
            .then(b.cmp(&a))
    });
    o
}

fn qvality(name: &str) -> Vec<(f64, f64, f64)> {
    std::fs::read_to_string(data(name))
        .unwrap()
        .lines()
        .skip(1)
        .map(|l| {
            let c: Vec<f64> = l.split('\t').map(|v| v.parse().unwrap()).collect();
            (c[0], c[1], c[2])
        })
        .collect()
}

fn printed(a: f64, b: f64) -> bool {
    (a - b).abs() <= 5e-6 * b.abs().max(1e-300)
}

fn compare(scores: &[f64], decoy: &[bool], q: &[f64], reference: &str) {
    let o = order(scores, decoy);
    let pep = pep_from_q(&o, scores, q, decoy).unwrap();
    let want = qvality(reference);
    assert_eq!(want.len(), scores.len());
    for (row, (&k, &(s, p, qv))) in o.iter().zip(&want).enumerate() {
        assert!(
            printed(scores[k], s),
            "row {row}: score {} vs {s}",
            scores[k]
        );
        assert!(printed(q[k], qv), "row {row}: q {} vs {qv}", q[k]);
        assert!(
            printed(pep[k], p),
            "row {row}: PEP {} vs qvality {p}",
            pep[k]
        );
    }
}

#[test]
fn pep_matches_qvality_tdc() {
    let (s, d) = fixture();
    let q = tdc_q(&s, &d, false);
    compare(&s, &d, &q, "pep_ref/q_d.txt");
}

#[test]
fn estimate_pi0_matches_qvality() {
    let (s, d) = fixture();
    let pi0 = estimate_pi0(&s, &d);
    assert!((pi0 - 0.409259).abs() < 5e-7, "pi0 {pi0}");
}

#[test]
fn pep_matches_qvality_mix_max() {
    let (s, d) = fixture();
    let q = mix_max_q(&s, &d, estimate_pi0(&s, &d));
    compare(&s, &d, &q, "pep_ref/qm_d.txt");
}

#[test]
fn pep_matches_percolator_svm_ref() {
    let mut scores = Vec::new();
    let mut q = Vec::new();
    let mut want = Vec::new();
    let mut decoy = Vec::new();
    for (name, is_decoy) in [
        ("svm_ref_perc.target.psms", false),
        ("svm_ref_perc.decoy.psms", true),
    ] {
        for l in std::fs::read_to_string(data(name)).unwrap().lines().skip(1) {
            let c: Vec<&str> = l.split('\t').collect();
            scores.push(c[1].parse::<f64>().unwrap());
            q.push(c[2].parse::<f64>().unwrap());
            want.push(c[3].parse::<f64>().unwrap());
            decoy.push(is_decoy);
        }
    }
    let o = order(&scores, &decoy);
    let pep = pep_from_q(&o, &scores, &q, &decoy).unwrap();
    for (i, (got, want)) in pep.iter().zip(&want).enumerate() {
        assert!(
            (got - want).abs() <= 1e-5 * want.max(1e-6),
            "row {i} score {}: PEP {got:e}, Percolator {want:e}",
            scores[i]
        );
    }
}

#[test]
fn pep_monotone_and_bounded() {
    let (s, d) = fixture();
    for q in [
        tdc_q(&s, &d, false),
        mix_max_q(&s, &d, estimate_pi0(&s, &d)),
    ] {
        let o = order(&s, &d);
        let pep = pep_from_q(&o, &s, &q, &d).unwrap();
        for w in o.windows(2) {
            assert!(
                pep[w[0]] <= pep[w[1]] + 1e-12,
                "{} above {}",
                pep[w[0]],
                pep[w[1]]
            );
        }
        assert!(pep.iter().all(|p| (0.0..=1.0).contains(p)));
    }
}

#[test]
fn few_psms_give_bounded_pep() {
    let s = [3.0, 2.0, 1.0];
    let d = [false, false, true];
    let q = tdc_q(&s, &d, false);
    let pep = pep_from_q(&order(&s, &d), &s, &q, &d).unwrap();
    assert!(pep.iter().all(|p| p.is_finite() && (0.0..=1.0).contains(p)));
}

#[test]
fn infinite_score_span_is_an_error() {
    let s = [1e308, -1e308, 0.0];
    let d = [false, false, true];
    let q = tdc_q(&s, &d, false);
    let err = pep_from_q(&order(&s, &d), &s, &q, &d).unwrap_err();
    assert!(err.contains("span is not a finite number"), "{err}");
}

#[test]
fn separation_check_as_check_separation() {
    let (s, d) = fixture();
    assert!(!too_good_separation(&s, &d));
    let s: Vec<f64> = (0..600).map(|i| i as f64).collect();
    let d: Vec<bool> = (0..600).map(|i| i < 300).collect();
    assert!(too_good_separation(&s, &d));
    let d: Vec<bool> = (0..600).map(|i| i < 100).collect();
    assert!(!too_good_separation(&s, &d));
}
