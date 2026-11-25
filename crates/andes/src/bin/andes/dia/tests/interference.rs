use super::*;
use onnx_infer::predictors::form::parse_peptide;

#[test]
fn interference_removes_the_lower_coeluting_sharer_only() {
    let (_, m) = parse_peptide("K.PEPTIDEK.L").unwrap();
    let (_, other) = parse_peptide("K.WWHHWWR.L").unwrap();
    let f = fragment_bins(&m, 2);
    assert_eq!(f.len(), 14);
    assert_eq!(fragment_bins(&m, 3).len(), 28);
    let p = |key, score, rt, mob, frags: &Vec<i64>| Coeluting {
        key,
        score,
        rt,
        mob,
        seq: format!("S{key}"),
        frags: frags.clone(),
    };
    let g = fragment_bins(&other, 2);
    let mut ps = vec![
        p(0, 5.0, 100.0, (10, 30), &f),
        p(1, 4.0, 101.0, (20, 40), &f),
        p(2, 3.0, 110.0, (20, 40), &f),
        p(3, 2.0, 100.0, (30, 50), &f),
        p(4, 1.0, 100.5, (10, 30), &g),
    ];
    let hit = interfered(&mut ps, (0.5, 4.0));
    let by_key: Vec<bool> = (0..5)
        .map(|k| hit[ps.iter().position(|x| x.key == k).unwrap()])
        .collect();
    // 1 is beaten by 0 (1 s apart, slices overlap); 2 is 9 s from 1 and 10 s from 0; 3 only
    // touches 0's slice edge but overlaps 1; 4 shares no fragment with any.
    assert_eq!(by_key, [false, true, false, true, false]);
    let mut touching = vec![p(0, 5.0, 100.0, (10, 30), &f), p(1, 4.0, 100.0, (30, 50), &f)];
    assert_eq!(interfered(&mut touching, (0.5, 4.0)), [false, false], "touching slices");
}

#[test]
fn interference_estimate_finds_the_coelution_scale_and_a_safe_share() {
    // Chance pairs over 0-200 s: shares 0..0.3 spread evenly (the null). Co-eluting pairs within
    // 2 s: share 0.6. The excess fades at 2 s; 0.3 is the null's top, so F is the first grid
    // value above it.
    let null = |n: usize| -> Vec<(f64, f64)> {
        (0..n)
            .map(|i| (i as f64 * 200.0 / n as f64, (i % 300) as f64 / 1000.0))
            .collect()
    };
    let mut pairs = null(40_000);
    pairs.extend((0..400).map(|i| (i as f64 * 2.0 / 400.0, 0.6)));
    let e = estimate_from_pairs(&pairs);
    let (f, drt) = e.rule.unwrap_or_else(|| panic!("{e}"));
    assert_eq!(drt, 2.0, "{e}");
    assert!((0.3..=0.35).contains(&f), "{e}");
    // Too few pairs: off.
    assert!(estimate_from_pairs(&pairs[..500]).rule.is_none());
    // Co-eluting pairs no likelier to share than chance (no excess at 0 s): off.
    assert!(estimate_from_pairs(&null(40_000)).rule.is_none());
    // Two forms of one sequence 100 s apart are no null pair; two sequences are.
    let (_, m) = parse_peptide("K.PEPTIDEK.L").unwrap();
    let c = |key, rt, seq: &str| Coeluting {
        key,
        score: 1.0,
        rt,
        mob: (0, 10),
        seq: seq.to_string(),
        frags: fragment_bins(&m, 2),
    };
    let null_pairs = |b: &str| {
        let w = HashMap::from([((0, 0), vec![c(0, 0.0, "PEPTIDEK"), c(1, 100.0, b)])]);
        estimate_interference(&w).null_pairs
    };
    assert_eq!(null_pairs("PEPTIDEK"), 0);
    assert_eq!(null_pairs("PEPTLDEK"), 1);
}

#[test]
fn interference_never_removes_a_form_of_the_same_sequence() {
    let (_, m) = parse_peptide("K.PEPTIDEK.L").unwrap();
    let f = fragment_bins(&m, 2);
    let p = |key, score, seq: &str| Coeluting {
        key,
        score,
        rt: 100.0,
        mob: (10, 30),
        seq: seq.to_string(),
        frags: f.clone(),
    };
    let mut same = vec![p(0, 5.0, "PEPTIDEK"), p(1, 4.0, "PEPTIDEK")];
    assert_eq!(interfered(&mut same, (0.5, 4.0)), [false, false]);
    let mut other = vec![p(0, 5.0, "PEPTIDEK"), p(1, 4.0, "PEPTLDEK")];
    let hit = interfered(&mut other, (0.5, 4.0));
    assert!(hit[other.iter().position(|x| x.key == 1).unwrap()]);
}
