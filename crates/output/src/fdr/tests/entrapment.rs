use super::*;

fn d(key: &str, pair: Option<&str>, entrapment: bool, score: f64) -> Discovery {
    Discovery {
        key: key.to_string(),
        pair: pair.map(str::to_string),
        entrapment,
        score,
        q_value: 0.0,
    }
}

#[test]
fn known_unobserved_counterparts_are_counted_below_the_cut() {
    let mut rows = vec![
        d("T1", Some("A"), false, 10.0),
        d("E1", Some("A"), true, 9.0),
        d("E2", Some("B"), true, 8.0),
        d("T2", Some("B"), false, 7.0),
        d("E3", Some("C"), true, 6.0),
        d("T3", Some("D"), false, 5.0),
    ];
    let c = curve(&mut rows, 1.0);
    let last = c[c.len() - 1];
    assert_eq!((last.targets, last.entrapment), (3, 3));
    assert_eq!((last.pair_below_cut, last.pair_lower), (1, 1));
    assert!((last.lower_bound - 0.5).abs() < 1e-12);
    assert!((last.combined - 1.0).abs() < 1e-12);
    assert_eq!(last.paired, 1.0);
    assert_eq!(at_q(&rows, &c, 0.01).paired, 1.0);
    let at_e2 = c[2];
    assert_eq!((at_e2.pair_below_cut, at_e2.pair_lower), (1, 0));
}

#[test]
fn a_target_found_after_its_entrapment_moves_it_from_below_the_cut_to_lower() {
    let mut rows = vec![d("E", Some("A"), true, 5.0), d("T", Some("A"), false, 4.0)];
    let c = curve(&mut rows, 1.0);
    assert_eq!((c[0].pair_below_cut, c[0].pair_lower), (1, 0));
    assert_eq!((c[1].pair_below_cut, c[1].pair_lower), (0, 1));
}

#[test]
fn the_combined_estimate_scales_with_the_ratio() {
    let mut rows = vec![d("T", None, false, 2.0), d("E", None, true, 1.0)];
    let c = curve(&mut rows, 2.0);
    assert!((c[1].combined - 0.75).abs() < 1e-12);
}

#[test]
fn best_per_key_keeps_the_highest_score() {
    let best = best_per_key([
        d("P", None, false, 1.0),
        d("P", None, false, 3.0),
        d("Q", None, true, 2.0),
    ]);
    assert_eq!(best.len(), 2);
    assert_eq!(best[0].score, 3.0);
}

#[test]
fn at_q_counts_the_rows_within_the_threshold() {
    let mut rows = vec![d("A", None, false, 3.0), d("B", None, true, 2.0)];
    rows[1].q_value = 0.02;
    let c = curve(&mut rows, 1.0);
    assert_eq!(at_q(&rows, &c, 0.01).entrapment, 0);
    assert_eq!(at_q(&rows, &c, 0.05).entrapment, 1);
}

#[test]
fn tied_scores_never_count_as_strict_entrapment_wins() {
    let mut rows = vec![d("E", Some("A"), true, 5.0), d("T", Some("A"), false, 5.0)];
    let c = curve(&mut rows, 1.0);
    assert_eq!(c[0].pair_below_cut, 0);
    assert_eq!(c[1].pair_lower, 0);
    assert_eq!(c[1].paired, 0.5);
}

#[test]
fn at_q_counts_only_accepted_rows_and_accepted_pairs() {
    let mut rows = vec![
        d("T", Some("A"), false, 6.0),
        d("E", Some("A"), true, 5.0),
        d("U", Some("B"), false, 4.0),
    ];
    rows[0].q_value = 0.2;
    let c = curve(&mut rows, 2.0);
    let accepted = at_q(&rows, &c, 0.01);
    assert_eq!((accepted.targets, accepted.entrapment), (1, 1));
    assert_eq!(accepted.combined, 0.75);
    assert_eq!(accepted.pair_below_cut, 0);
    assert_eq!(accepted.paired, 0.5);
}

#[test]
fn ambiguous_missing_and_many_to_one_pairs_make_only_paired_unavailable() {
    for mut rows in [
        vec![d("E", None, true, 5.0), d("T", Some("A"), false, 4.0)],
        vec![
            d("E", Some("A"), true, 5.0),
            d("F", Some("A"), true, 4.0),
            d("T", Some("A"), false, 3.0),
        ],
        vec![
            d("E", Some("A"), true, 5.0),
            d("T", Some("A"), false, 4.0),
            d("U", Some("A"), false, 3.0),
        ],
    ] {
        let c = curve(&mut rows, 1.0);
        assert!(c.iter().all(|e| e.paired.is_nan()));
        assert!(c
            .iter()
            .all(|e| e.lower_bound.is_finite() && e.combined.is_finite()));
        assert!(at_q(&rows, &c, 0.01).paired.is_nan());
    }
}

#[test]
fn exact_best_score_ties_are_independent_of_input_order() {
    let rows = vec![d("P", Some("B"), true, 5.0), d("P", Some("A"), false, 5.0)];
    let forward = best_per_key(rows.clone());
    let reverse = best_per_key(rows.into_iter().rev());
    assert_eq!(forward, reverse);
}

#[test]
fn at_q_recomputes_pair_membership_instead_of_reusing_a_prefix() {
    let mut rows = vec![
        d("E", Some("A"), true, 6.0),
        d("T", Some("A"), false, 5.0),
        d("U", Some("B"), false, 4.0),
    ];
    rows[1].q_value = 0.2;
    let c = curve(&mut rows, 1.0);
    let accepted = at_q(&rows, &c, 0.01);
    assert_eq!((accepted.targets, accepted.entrapment), (1, 1));
    assert_eq!((accepted.pair_below_cut, accepted.pair_lower), (1, 0));
    assert_eq!(accepted.paired, 1.0);
}

#[test]
fn best_per_key_prefers_lower_q_then_lexicographic_pair_for_exact_score_ties() {
    let mut better_q = d("P", Some("Z"), true, 5.0);
    better_q.q_value = 0.01;
    let mut worse_q = d("P", Some("A"), false, 5.0);
    worse_q.q_value = 0.02;
    assert_eq!(best_per_key([worse_q, better_q.clone()]), vec![better_q]);
    let a = d("P", Some("A"), true, 5.0);
    let b = d("P", Some("B"), true, 5.0);
    assert_eq!(best_per_key([b.clone(), a.clone()]), vec![a.clone()]);
    assert_eq!(best_per_key([a.clone(), b]), vec![a]);
}

#[test]
fn empty_acceptance_has_zero_counts_without_losing_pair_applicability() {
    let mut rows = vec![d("E", Some("A"), true, 5.0), d("T", Some("A"), false, 4.0)];
    let c = curve(&mut rows, 0.0);
    let accepted = at_q(&rows, &c, -1.0);
    assert_eq!((accepted.targets, accepted.entrapment), (0, 0));
    assert_eq!(accepted.paired, 0.0);
    assert!(accepted.combined.is_nan());
    assert_eq!(at_q(&[], &[], 0.01).targets, 0);
}

#[test]
fn valid_one_to_one_pairing_keeps_both_strict_win_terms() {
    let mut rows = vec![
        d("T1", Some("A"), false, 10.0),
        d("E1", Some("A"), true, 9.0),
        d("E2", Some("B"), true, 8.0),
        d("T2", Some("B"), false, 7.0),
        d("E3", Some("C"), true, 6.0),
        d("T3", Some("C"), false, 5.0),
    ];
    let c = curve(&mut rows, 1.0);
    assert_eq!((c[4].pair_below_cut, c[4].pair_lower), (1, 1));
    assert_eq!(c[4].paired, 6.0 / 5.0);
    assert_eq!((c[5].pair_below_cut, c[5].pair_lower), (0, 2));
    assert_eq!(c[5].paired, 7.0 / 6.0);
}

#[test]
fn empty_pair_identifiers_and_nonfinite_scores_cannot_validate_pairing() {
    for mut rows in [
        vec![d("E", Some(""), true, 5.0), d("T", Some(""), false, 4.0)],
        vec![
            d("E", Some("A"), true, f64::NAN),
            d("T", Some("A"), false, 4.0),
        ],
    ] {
        let c = curve(&mut rows, 1.0);
        assert!(c.iter().all(|estimate| estimate.paired.is_nan()));
        assert!(c.iter().all(|estimate| estimate.combined.is_finite()));
    }
}

#[test]
fn a_known_unobserved_pair_is_distinct_from_an_unknown_pair_identity() {
    let mut known = vec![d("E", Some("source"), true, 5.0)];
    let known_curve = curve(&mut known, 1.0);
    assert_eq!(known_curve[0].pair_below_cut, 1);
    assert_eq!(known_curve[0].paired, 2.0);
    assert_eq!(at_q(&known, &known_curve, 0.01).paired, 2.0);
    let mut unknown = vec![d("E", None, true, 5.0)];
    let unknown_curve = curve(&mut unknown, 1.0);
    assert!(unknown_curve[0].paired.is_nan());
    assert!(at_q(&unknown, &unknown_curve, 0.01).paired.is_nan());
    assert_eq!(known_curve[0].lower_bound, unknown_curve[0].lower_bound);
    assert_eq!(known_curve[0].combined, unknown_curve[0].combined);
}

#[test]
fn repeated_entrapment_pairs_remain_ambiguous_even_without_observed_targets() {
    let mut rows = vec![
        d("E", Some("source"), true, 5.0),
        d("F", Some("source"), true, 4.0),
    ];
    let c = curve(&mut rows, 1.0);
    assert!(c.iter().all(|estimate| estimate.paired.is_nan()));
    assert!(at_q(&rows, &c, 0.01).paired.is_nan());
}
