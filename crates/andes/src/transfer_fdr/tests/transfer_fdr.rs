use super::*;

struct Normal(u64);

impl Normal {
    fn unit(&mut self) -> f64 {
        let z = model::hash::splitmix64_next(&mut self.0);
        ((z >> 11) as f64 + 0.5) / (1u64 << 53) as f64
    }
    fn draw(&mut self) -> f32 {
        let (u1, u2) = (self.unit(), self.unit());
        ((-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()) as f32
    }
    fn vector(&mut self, mean: f32) -> Vec<f32> {
        (0..TRANSFER_FEATURES.len())
            .map(|_| mean + self.draw())
            .collect()
    }
}

fn rows(n: u32, target_mean: f32, families: &[u8], seed: u64) -> Vec<FdrRow> {
    let mut g = Normal(seed);
    let mut out = Vec::new();
    for p in 0..n {
        out.push(FdrRow {
            precursor: p,
            is_decoy: false,
            family: 0,
            x: g.vector(target_mean),
        });
        for &f in families {
            out.push(FdrRow {
                precursor: p,
                is_decoy: true,
                family: f,
                x: g.vector(0.0),
            });
        }
    }
    out
}

fn accepted_targets(rows: &[FdrRow], q: &[Option<f64>]) -> usize {
    rows.iter()
        .zip(q)
        .filter(|(r, q)| !r.is_decoy && q.is_some_and(|q| q <= 0.01))
        .count()
}

#[test]
fn separable_targets_get_low_q() {
    let rows = rows(500, 3.0, &[1, 2], 11);
    let r = transfer_qvalues(&rows, 42);
    let n = accepted_targets(&rows, &r.q);
    eprintln!(
        "separable: {n} of 500 targets at q <= 0.01 (mass-only {}, rt-only {})",
        accepted_targets(&rows, &r.q_mass_only),
        accepted_targets(&rows, &r.q_rt_only)
    );
    assert!(n >= 400, "{n} targets at q <= 0.01");
    assert!(r.score.iter().all(|s| s.is_finite()));
}

#[test]
fn null_targets_are_not_accepted() {
    for seed in [3u64, 5, 9] {
        let rows = rows(500, 0.0, &[1, 2], seed);
        let r = transfer_qvalues(&rows, seed);
        let n = accepted_targets(&rows, &r.q);
        eprintln!("null seed {seed}: {n} targets at q <= 0.01");
        assert!(n <= 5, "seed {seed}: {n} null targets at q <= 0.01");
    }
}

#[test]
fn null_experiments_rarely_accept_a_target() {
    let mut g = Normal(77);
    let n_exp = 10_000;
    let mut with_false = 0usize;
    let mut q = vec![None; 60];
    for _ in 0..n_exp {
        let items: Vec<(usize, f32, u8)> = (0..60)
            .map(|i| (i, g.draw(), if i >= 20 { MASS_SHIFT } else { TARGET }))
            .collect();
        q.fill(None);
        kdecoy_q(&items, &mut q);
        if q[..20].iter().any(|q| q.is_some_and(|q| q <= 0.01)) {
            with_false += 1;
        }
    }
    let frac = with_false as f64 / n_exp as f64;
    eprintln!(
        "null experiments: {with_false} of {n_exp} ({frac:.4}) accept a target at q <= 0.01"
    );
    assert!(
        frac <= 0.015,
        "{frac} of the null experiments accept a target"
    );
}

#[test]
fn mixture_fdp_stays_near_the_estimate() {
    let mut g = Normal(21);
    let mut rows = Vec::new();
    let mut present = Vec::new();
    for p in 0..4000u32 {
        let real = p % 10 < 3;
        let mut x = g.vector(0.0);
        if real {
            for v in &mut x[..4] {
                *v += 1.5;
            }
        }
        rows.push(FdrRow {
            precursor: p,
            is_decoy: false,
            family: 0,
            x,
        });
        present.push(real);
        for f in [1, 2] {
            rows.push(FdrRow {
                precursor: p,
                is_decoy: true,
                family: f,
                x: g.vector(0.0),
            });
            present.push(false);
        }
    }
    let r = transfer_qvalues(&rows, 17);
    for thr in [0.01, 0.05, 0.1] {
        let (mut acc, mut false_acc) = (0usize, 0usize);
        for (i, row) in rows.iter().enumerate() {
            if !row.is_decoy && r.q[i].is_some_and(|q| q <= thr) {
                acc += 1;
                false_acc += usize::from(!present[i]);
            }
        }
        let fdp = false_acc as f64 / acc.max(1) as f64;
        eprintln!("mixture: q <= {thr}: {acc} accepted, {false_acc} absent, FDP {fdp:.4}");
        assert!(acc > 100, "q <= {thr}: {acc} accepted");
        assert!(fdp <= 2.0 * thr + 0.01, "q <= {thr}: FDP {fdp}");
    }
}

#[test]
fn folds_group_by_precursor() {
    let rows = rows(600, 3.0, &[1, 2], 1);
    let mut fold_of_precursor = std::collections::HashMap::new();
    let mut sizes = [0usize; N_FOLDS];
    for r in &rows {
        let f = fold_of(r.precursor, 7);
        assert!(f < N_FOLDS);
        let first = *fold_of_precursor.entry(r.precursor).or_insert(f);
        assert_eq!(first, f, "precursor {} spans folds", r.precursor);
    }
    for f in fold_of_precursor.values() {
        sizes[*f] += 1;
    }
    assert!(sizes.iter().all(|&s| s > 150), "fold sizes {sizes:?}");
    let moved = (0..600u32)
        .filter(|&p| fold_of(p, 7) != fold_of(p, 8))
        .count();
    assert!(moved > 200, "the seed changes {moved} of 600 folds");
}

#[test]
fn fdr_degenerate_inputs_yield_none() {
    let empty = transfer_qvalues(&[], 1);
    assert!(empty.q.is_empty() && empty.score.is_empty());

    let no_decoys = rows(300, 3.0, &[], 2);
    let r = transfer_fdr(&no_decoys, 1);
    assert!(r.result.q.iter().all(Option::is_none));
    assert_eq!(r.counts.failed_folds.len(), N_FOLDS);
    assert_eq!(r.counts.unscored, no_decoys.len());

    let only_decoys: Vec<FdrRow> = rows(300, 3.0, &[1, 2], 3)
        .into_iter()
        .filter(|r| r.is_decoy)
        .collect();
    let r = transfer_qvalues(&only_decoys, 1);
    assert!(r.q.iter().all(Option::is_none));

    let seed = 5;
    let mixed: Vec<FdrRow> = rows(600, 3.0, &[1, 2], 4)
        .into_iter()
        .filter(|r| !(r.is_decoy && fold_of(r.precursor, seed) == 0))
        .collect();
    let r = transfer_fdr(&mixed, seed);
    for (row, q) in mixed.iter().zip(&r.result.q) {
        if fold_of(row.precursor, seed) == 0 {
            assert!(q.is_none(), "one-class fold row got q {q:?}");
        }
    }
    assert!(accepted_targets(&mixed, &r.result.q) > 200);
    assert!(r.counts.no_q > 0);
}

#[test]
fn per_family_estimates_reported() {
    let rows = rows(500, 3.0, &[1], 6);
    let r = transfer_qvalues(&rows, 42);
    assert!(r.q_rt_only.iter().all(Option::is_none));
    assert!(accepted_targets(&rows, &r.q_mass_only) >= 400);
    assert_eq!(r.q, r.q_mass_only);
}

#[test]
fn q_values_are_monotone_in_score_within_a_fold() {
    let rows = rows(400, 1.0, &[1, 2], 8);
    let r = transfer_qvalues(&rows, 3);
    let mut by_fold: Vec<Vec<(f32, f64)>> = vec![Vec::new(); N_FOLDS];
    for (i, row) in rows.iter().enumerate() {
        if let Some(q) = r.q[i] {
            by_fold[fold_of(row.precursor, 3)].push((r.score[i], q));
        }
    }
    for f in by_fold.iter_mut() {
        assert!(!f.is_empty());
        f.sort_by(|a, b| a.0.total_cmp(&b.0));
        for w in f.windows(2) {
            assert!(w[1].1 <= w[0].1, "q rises with score: {w:?}");
        }
        assert!(f.iter().all(|&(_, q)| (0.0..=1.0).contains(&q)));
    }
}

#[test]
fn k_follows_the_decoy_to_target_ratio() {
    let mut rows = Vec::new();
    for p in 0..400u32 {
        let v = vec![p as f32; TRANSFER_FEATURES.len()];
        rows.push(FdrRow {
            precursor: p,
            is_decoy: false,
            family: 0,
            x: v.clone(),
        });
        for _ in 0..4 {
            rows.push(FdrRow {
                precursor: p,
                is_decoy: true,
                family: 1,
                x: v.clone(),
            });
        }
    }
    let mut q = vec![None; rows.len()];
    let scored: Vec<(usize, f32, u8)> = rows
        .iter()
        .enumerate()
        .map(|(i, r)| (i, r.x[0], r.family))
        .collect();
    kdecoy_q(&scored, &mut q);
    assert!(
        q.iter().all(|q| (q.unwrap() - 1.0).abs() < 1e-12),
        "{:?}",
        &q[..5]
    );
}

#[test]
fn pseudo_rows_take_the_q_of_the_nearest_lower_score() {
    let rows = rows(500, 3.0, &[1, 2], 12);
    let fdr = transfer_fdr(&rows, 4);
    for (i, row) in rows
        .iter()
        .enumerate()
        .filter(|(_, r)| !r.is_decoy)
        .take(50)
    {
        assert_eq!(fdr.pseudo_q(row.precursor, &row.x), fdr.result.q[i]);
    }
    let low = vec![-50.0f32; TRANSFER_FEATURES.len()];
    assert!(fdr.pseudo_q(0, &low).is_some_and(|q| q > 0.1));
}

fn with_identity(n: u32, n_id: u32, seed: u64) -> Vec<FdrRow> {
    let mut out = rows(n, 3.0, &[MASS_SHIFT, RT_SHIFT], seed);
    let mut g = Normal(seed ^ 0xabc);
    out.extend((0..n_id).map(|i| FdrRow {
        precursor: n + i,
        is_decoy: true,
        family: IDENTITY,
        x: g.vector(3.0),
    }));
    out
}

#[test]
fn identity_decoys_add_their_share_to_the_estimate() {
    let rows = with_identity(1000, 30, 31);
    {
        let r = transfer_fdr(&rows, 5).result;
        let wp = accepted_targets(&rows, &r.q_wrong_peak);
        let all = accepted_targets(&rows, &r.q);
        eprintln!("identity: wrong-peak {wp}, combined {all} at q <= 0.01");
        assert!(wp >= 800, "{wp}");
        assert!(5 * all < wp, "{all} of {wp}");
        let at10 = rows
            .iter()
            .zip(&r.q)
            .filter(|(row, q)| row.family == TARGET && q.is_some_and(|q| q <= 0.1))
            .count();
        assert!(at10 >= 800, "{at10}");
        assert!(rows
            .iter()
            .zip(&r.q)
            .filter(|(row, _)| row.family == IDENTITY)
            .all(|(_, q)| q.is_some()));
    }
}

#[test]
fn kdecoy_counts_identity_decoys_once_each() {
    let mut items: Vec<(usize, f32, u8)> =
        (0..100).map(|i| (i, 1000.0 - i as f32, TARGET)).collect();
    items.extend((100..150).map(|i| (i, 10.0 - i as f32 * 0.01, MASS_SHIFT)));
    items.extend([(150, 0.5, IDENTITY), (151, 0.4, IDENTITY)]);
    let mut q = vec![None; items.len()];
    kdecoy_q(&items, &mut q);
    assert!((q[0].unwrap() - 1.0 / 50.0).abs() < 1e-12, "{:?}", q[0]);
    assert_eq!(q[151], Some(1.0));
    items[150].1 = 2000.0;
    kdecoy_q(&items, &mut q);
    assert!((q[0].unwrap() - (1.0 / 50.0 + 1.0 / 100.0)).abs() < 1e-12);
}

#[test]
fn results_are_deterministic_for_a_seed() {
    let rows = rows(300, 1.5, &[1, 2], 13);
    assert_eq!(transfer_qvalues(&rows, 9), transfer_qvalues(&rows, 9));
}

#[test]
fn step_functions_match_the_per_row_estimate() {
    let mut g = Normal(99);
    for trial in 0..20u64 {
        let n = 50 + 37 * trial as usize;
        let items: Vec<(usize, f32, u8)> = (0..n)
            .map(|i| {
                let fam = match i % 7 {
                    0..=2 => TARGET,
                    3 | 4 => MASS_SHIFT,
                    5 => RT_SHIFT,
                    _ => IDENTITY,
                };
                let shift = if fam == TARGET { 1.0 } else { 0.0 };
                let mut s = ((g.draw() + shift) * 4.0).round() / 4.0;
                if s == 0.0 && i % 2 == 0 {
                    s = -0.0;
                }
                (i, s, fam)
            })
            .collect();
        let mut by_family: [Vec<f32>; 4] = Default::default();
        for &(_, s, f) in &items {
            by_family[f as usize].push(canonical(s));
        }
        for v in &mut by_family {
            v.sort_by(f32::total_cmp);
        }
        for e in ESTIMATES {
            let kept: Vec<(usize, f32, u8)> =
                items.iter().copied().filter(|it| e.keeps(it.2)).collect();
            let mut want = vec![None; n];
            kdecoy_q(&kept, &mut want);
            let steps = kdecoy_steps(&by_family, e);
            for &(i, s, f) in &items {
                let got = steps
                    .as_ref()
                    .filter(|_| e.keeps(f))
                    .map(|st| st.of_member(s));
                assert_eq!(got, want[i], "trial {trial} {e:?} row {i} score {s}");
            }
            if let (Some(st), Estimate::Pooled) = (&steps, e) {
                let mut table: Vec<(f32, f64)> = kept
                    .iter()
                    .filter_map(|it| want[it.0].map(|q| (it.1, q)))
                    .collect();
                table.sort_by(|a, b| a.0.total_cmp(&b.0).then(b.1.total_cmp(&a.1)));
                for probe in [-9.0f32, -0.3, -0.0, 0.0, 0.1, 0.25, 1.3, 9.0, f32::NAN] {
                    let at = table.partition_point(|&(s, _)| s <= probe);
                    let old = if at == 0 { table[0].1 } else { table[at - 1].1 };
                    assert_eq!(st.at_or_below(probe), old, "probe {probe}");
                }
            }
        }
    }
}

#[test]
fn threshold_steps_accept_as_the_whole_ones() {
    let mut g = Normal(7);
    let (mut whole, mut cut) = (FoldScores::new(3), FoldScores::new(3));
    let mut items = Vec::new();
    for i in 0..6000u32 {
        let fam = match i % 7 {
            0..=2 => TARGET,
            3 | 4 => MASS_SHIFT,
            5 => RT_SHIFT,
            _ => IDENTITY,
        };
        let shift = if fam == TARGET { 1.5 } else { 0.0 };
        let s = Some(((g.draw() + shift) * 8.0).round() / 8.0);
        let fam = u32::from(fam);
        items.push((i % 800, fam, s));
        whole.add(i % 800, fam, s);
        cut.add(i % 800, fam, s);
    }
    let ts = [0.01, 0.05, 0.3];
    let whole = whole.finish(&[]);
    let cut = cut.finish_with(&[], Some(&ts), |_, _| Vec::new());
    assert!(cut.steps() < whole.steps());
    assert_eq!(cut.pooled_steps(), whole.pooled_steps());
    for &(p, f, s) in &items {
        for e in ESTIMATES {
            for t in ts {
                assert_eq!(whole.accepts(e, p, f, s, t), cut.accepts(e, p, f, s, t));
            }
        }
        assert_eq!(
            whole.q(Estimate::Pooled, p, f, s),
            cut.q(Estimate::Pooled, p, f, s)
        );
        assert_eq!(whole.pseudo_q(p, 0, s), cut.pseudo_q(p, 0, s));
    }
}

fn bagged_q(rows: &[FdrRow], members: usize, seed: u64) -> Vec<Option<f64>> {
    let sets: Vec<&[FdrRow]> = (0..members).map(|_| rows).collect();
    let models = train_bagged_folds(&sets, seed, &TrainParams::default());
    assert_eq!(models.members(), members);
    let precursors: Vec<u32> = rows.iter().map(|r| r.precursor).collect();
    let xs: Vec<&[f32]> = rows.iter().map(|r| r.x.as_slice()).collect();
    let mut scores = Vec::new();
    models.score(&precursors, &xs, &mut scores);
    let mut fs = FoldScores::new(seed);
    for (r, s) in rows.iter().zip(&scores) {
        fs.add(r.precursor, r.family.into(), *s);
    }
    let q = fs.finish(&models.failed_folds);
    rows.iter()
        .zip(&scores)
        .map(|(r, s)| q.q(Estimate::Pooled, r.precursor, r.family.into(), *s))
        .collect()
}

#[test]
fn one_member_is_the_single_model() {
    let rows = rows(400, 1.0, &[1, 2], 23);
    assert_eq!(bagged_q(&rows, 1, 6), transfer_qvalues(&rows, 6).q);
}

#[test]
fn bagged_scores_are_the_mean_of_the_members() {
    let rows = rows(300, 1.0, &[1, 2], 29);
    let sets: Vec<&[FdrRow]> = vec![&rows; 3];
    let bagged = train_bagged_folds(&sets, 4, &TrainParams::default());
    let xs: Vec<&[f32]> = rows.iter().map(|r| r.x.as_slice()).collect();
    let precursors: Vec<u32> = rows.iter().map(|r| r.precursor).collect();
    let mut got = Vec::new();
    bagged.score(&precursors, &xs, &mut got);
    for (fold, members) in bagged.models.iter().enumerate() {
        assert_eq!(members.len(), 3);
        for (i, r) in rows.iter().enumerate() {
            if fold_of(r.precursor, 4) != fold {
                continue;
            }
            let sum: f32 = members.iter().map(|m| m.predict_value(&r.x)).sum();
            let want = sum / 3.0;
            let g = got[i].unwrap();
            assert!(
                (g - want).abs() <= 1e-5 * want.abs().max(1.0),
                "{g} vs {want}"
            );
        }
    }
    let m = &bagged.models[0];
    assert_ne!(
        m[0].predict_value(&rows[0].x),
        m[1].predict_value(&rows[0].x)
    );
}

#[test]
fn bagging_keeps_the_mixture_calibrated_and_steadier_across_seeds() {
    let mut g = Normal(41);
    let mut rows = Vec::new();
    let mut present = Vec::new();
    for p in 0..3000u32 {
        let real = p % 10 < 3;
        let mut x = g.vector(0.0);
        if real {
            for v in &mut x[..4] {
                *v += 1.2;
            }
        }
        rows.push(FdrRow {
            precursor: p,
            is_decoy: false,
            family: 0,
            x,
        });
        present.push(real);
        for f in [1, 2] {
            rows.push(FdrRow {
                precursor: p,
                is_decoy: true,
                family: f,
                x: g.vector(0.0),
            });
            present.push(false);
        }
    }
    let accepted = |q: &[Option<f64>]| -> std::collections::BTreeSet<usize> {
        (0..rows.len())
            .filter(|&i| !rows[i].is_decoy && q[i].is_some_and(|q| q <= 0.05))
            .collect()
    };
    let jaccard = |a: &std::collections::BTreeSet<usize>,
                   b: &std::collections::BTreeSet<usize>| {
        a.intersection(b).count() as f64 / a.union(b).count().max(1) as f64
    };
    let mut overlap = Vec::new();
    for members in [1, 8] {
        let (a, b) = (
            accepted(&bagged_q(&rows, members, 1)),
            accepted(&bagged_q(&rows, members, 2)),
        );
        for set in [&a, &b] {
            let false_acc = set.iter().filter(|&&i| !present[i]).count();
            let fdp = false_acc as f64 / set.len().max(1) as f64;
            assert!(set.len() > 100, "{members} members: {} accepted", set.len());
            assert!(fdp <= 0.11, "{members} members: FDP {fdp}");
        }
        let j = jaccard(&a, &b);
        eprintln!(
            "{members} members: accepted {} and {}, Jaccard {j:.4}",
            a.len(),
            b.len()
        );
        overlap.push(j);
    }
    assert!(overlap[1] > overlap[0], "{overlap:?}");
}

#[test]
fn training_sample_wants_what_push_keeps() {
    let mut a: TrainingSample<usize, u32> = TrainingSample::new(500, 9);
    let mut b: TrainingSample<usize, u32> = TrainingSample::new(500, 9);
    for i in 0..5000usize {
        let unit = (i / 3) as u64;
        a.push(i, unit, i as u32);
        if b.wants(unit) {
            b.push(i, unit, i as u32);
        } else {
            b.skip();
        }
    }
    assert_eq!(a.stats(), b.stats());
    assert_eq!(a.into_items(), b.into_items());
}

#[test]
fn training_sample_keeps_whole_precursors_whatever_the_order() {
    let rows = rows(3000, 1.0, &[1, 2], 5);
    let mut a: TrainingSample<usize> = TrainingSample::new(2500, 42);
    for (i, r) in rows.iter().enumerate() {
        a.push(i, r.precursor as u64, r.clone());
    }
    let mut b: TrainingSample<usize> = TrainingSample::new(2500, 42);
    for (i, r) in rows.iter().enumerate().rev() {
        b.push(i, r.precursor as u64, r.clone());
    }
    let (seen, kept, share) = a.stats();
    assert_eq!(seen, 9000);
    assert!(kept <= 2500 && kept > 1000, "{kept}");
    assert!(share < 1.0);
    let (ra, rb) = (a.into_items(), b.into_items());
    assert_eq!(
        ra, rb,
        "the kept set and its order do not depend on arrival"
    );
    let mut per: std::collections::HashMap<u32, usize> = std::collections::HashMap::new();
    for r in &ra {
        *per.entry(r.precursor).or_default() += 1;
    }
    assert!(per.values().all(|&c| c == 3));
    let mut all: TrainingSample<usize> = TrainingSample::new(10_000, 42);
    for (i, r) in rows.iter().enumerate() {
        all.push(i, r.precursor as u64, r.clone());
    }
    assert_eq!(all.into_items(), rows);
}

#[test]
fn the_transfer_classifier_differs_from_the_default_trainer_in_its_learning_rate() {
    let (t, d) = (transfer_train_params(), TrainParams::default());
    assert_eq!(t.learning_rate, TRANSFER_LEARNING_RATE);
    assert!(t.learning_rate > d.learning_rate);
    assert_eq!(
        (t.n_rounds, t.max_depth, t.n_bins, t.early_stop_rounds),
        (d.n_rounds, d.max_depth, d.n_bins, d.early_stop_rounds)
    );
    let rows = rows(400, 3.0, &[1, 2], 31);
    let sets: Vec<&[FdrRow]> = vec![&rows];
    let models = train_bagged_folds(&sets, 3, &t);
    assert!(models.failed_folds.is_empty());
    let precursors: Vec<u32> = rows.iter().map(|r| r.precursor).collect();
    let xs: Vec<&[f32]> = rows.iter().map(|r| r.x.as_slice()).collect();
    let mut scores = Vec::new();
    models.score(&precursors, &xs, &mut scores);
    let mean = |decoy: bool| {
        let v: Vec<f32> = rows
            .iter()
            .zip(&scores)
            .filter(|(r, _)| r.is_decoy == decoy)
            .filter_map(|(_, s)| *s)
            .collect();
        v.iter().sum::<f32>() / v.len() as f32
    };
    assert!(
        mean(false) > mean(true) + 1.0,
        "{} {}",
        mean(false),
        mean(true)
    );
}

fn two_strata_scores(strata: usize) -> FoldScores {
    let mut g = Normal(77);
    let mut fs = FoldScores::with_strata(9, strata);
    for p in 0..30_000u32 {
        let stratum = (p % 2) as u8;
        let packed = |family: u8| pack_family(family, stratum);
        fs.add(p, packed(TARGET), Some(g.draw() + 2.5));
        for f in [MASS_SHIFT, RT_SHIFT] {
            let hit = stratum == 1 && g.unit() < 0.05;
            fs.add(p, packed(f), Some(g.draw() + if hit { 2.5 } else { 0.0 }));
        }
    }
    fs
}

#[test]
fn strata_calibrate_each_stratum_on_its_own_decoys() {
    let pooled = two_strata_scores(1).finish(&[]);
    let strat = two_strata_scores(5).finish(&[]);
    let at = |q: &QValues, stratum: u8, s: f32| {
        q.q(
            Estimate::Pooled,
            (stratum as u32) + 2 * 3,
            pack_family(TARGET, stratum),
            Some(s),
        )
    };
    let (a, b) = (at(&strat, 0, 2.0).unwrap(), at(&strat, 1, 2.0).unwrap());
    assert!(b > 2.0 * a, "stratum 0 q {a}, stratum 1 q {b}");
    let p0 = at(&pooled, 0, 2.0).unwrap();
    let p1 = at(&pooled, 1, 2.0).unwrap();
    assert!((p0 - p1).abs() < 0.2 * p0.max(p1) + 1e-9, "{p0} {p1}");
    assert!(a < p0 && b > p0, "{a} < {p0} < {b}");
}

#[test]
fn prior_strata_follow_the_edges() {
    let e = [0.0, 0.05, 0.2, 0.5];
    assert_eq!(stratum_of(0.0, &e), 0);
    assert_eq!(stratum_of(0.01, &e), 1);
    assert_eq!(stratum_of(0.05, &e), 1);
    assert_eq!(stratum_of(0.3, &e), 3);
    assert_eq!(stratum_of(1.0, &e), 4);
    assert_eq!(stratum_of(f64::NAN, &e), 0);
    assert_eq!(unpack_family(pack_family(IDENTITY, 4u32)), (IDENTITY, 4));
}

#[test]
fn shrinkage_follows_the_between_stratum_spread() {
    let apart = [(5000.0, 5.0), (5000.0, 100.0), (5000.0, 400.0)];
    assert!(stratum_kappa(&apart, 1.0) < 50.0);
    let alike = [(5000.0, 52.0), (5000.0, 48.0), (5000.0, 50.0)];
    assert!(stratum_kappa(&alike, 1.0).is_infinite());
    assert!(stratum_kappa(&[(100.0, 1.0)], 1.0).is_infinite());
}
