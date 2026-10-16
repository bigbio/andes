//! Transfer FDR of match-between-runs: a cross-validated GBDT separates
//! transfer targets from their decoys, and the K-decoy estimator turns the
//! held-out scores into q-values.
//!
//! The transfer q-value is an internal quant-quality score that decides which
//! transferred MS1 quantities are reported; it is not an identification FDR.
//! The identified rows' quality q-value (`quant_q`) is the same estimate over
//! the identified targets and their decoys, fitted separately with
//! [`IDENTIFIED_FEATURES`].
//!
//! Folds: a hash of the precursor id (splitmix64, seeded) gives each precursor
//! one of three folds, so its targets and decoys in every run share a fold.
//! For each fold a classifier (target = 1, decoy = 0, validation groups = the
//! precursor) is trained on the other two folds and scores the held-out fold.
//! Scores are the raw boosted margins and are compared within a fold only.
//!
//! q-values, per fold: rank the fold's scored rows by score; at threshold `c`,
//! `FDR(c) = (D(c) + 1) / (K * max(T(c), 1))`, capped at 1, with `D`/`T` the
//! decoys/targets at or above `c` and `K` the fold's decoy-to-target row
//! ratio; q is the running minimum of `FDR` from the lowest score up. The `+1`
//! is the finite-sample correction: a target above every decoy is not given
//! q = 0, and no target reaches q <= 0.01 unless at least `1 / (0.01 K)`
//! targets score above every decoy. The accepted set at
//! any threshold is the union of the folds' accepted sets, so its estimated FDR
//! is at most the threshold. The same estimate is made with one decoy family
//! at a time (`K` = that family's ratio).
//!
//! `K` is counted over the rows passed in, which are the rows that have
//! features: a target or decoy without evidence takes no part. When decoys are
//! missing (no RT-shift decoy near a run edge, a window too short to fit), `K`
//! falls below the nominal number of decoys per target and each remaining
//! decoy counts for more targets, rather than the missing decoys being counted
//! as absent false hits; this errs on the conservative side.
//!
//! A row gets no score when its fold's training set is single-class or fails
//! the trainer's quality gate, and no q-value when it has no score or its fold
//! (restricted to the estimate's families) lacks targets or decoys.

use model_train::gbdt::train::{train_gbdt, Dataset, TrainParams};
use rayon::prelude::*;
use scoring_crate::gbdt_eval::GbdtPeakModel;

/// Cross-validation folds.
pub const N_FOLDS: usize = 3;

/// Feature names of the transfer classifier, in [`FdrRow::x`] order.
pub const TRANSFER_FEATURES: &[&str] = &[
    "llr",
    "rt_z",
    "abs_rt_z",
    "width_ratio",
    "envelope_cos",
    "coherence",
    "m1_ratio",
    "abs_mass_err_ppm",
    "n_points",
    "n_isotopes_used",
    "log10_amplitude_snr",
    "n_donor_runs",
    "library_q",
    "shape_log_width",
    "shape_tailing_diff",
    "shape_iso_cos",
];

/// Feature names of the identified-row quality classifier: the template-fit
/// evidence of [`TRANSFER_FEATURES`] without `n_donor_runs` and `library_q`,
/// which are the same for an identified target and its decoys, and without
/// the donor shape comparison (an identified row is itself a donor).
pub const IDENTIFIED_FEATURES: &[&str] = &[
    "llr",
    "rt_z",
    "abs_rt_z",
    "width_ratio",
    "envelope_cos",
    "coherence",
    "m1_ratio",
    "abs_mass_err_ppm",
    "n_points",
    "n_isotopes_used",
    "log10_amplitude_snr",
];

/// One target or decoy to be scored.
#[derive(Debug, Clone, PartialEq)]
pub struct FdrRow {
    pub precursor: u32,
    pub is_decoy: bool,
    /// 0 target, 1 mass-shift decoy, 2 RT-shift decoy.
    pub family: u8,
    pub x: Vec<f32>,
}

/// Per-row output of [`transfer_qvalues`], in input order.
#[derive(Debug, Clone, PartialEq)]
pub struct FdrResult {
    /// Both decoy families.
    pub q: Vec<Option<f64>>,
    /// Mass-shift decoys only; `None` for RT-shift decoy rows.
    pub q_mass_only: Vec<Option<f64>>,
    /// RT-shift decoys only; `None` for mass-shift decoy rows.
    pub q_rt_only: Vec<Option<f64>>,
    /// Held-out raw margin, comparable within a fold; NaN when unscored.
    pub score: Vec<f32>,
}

/// What the cross-validation could not score, for the caller to log.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct FdrCounts {
    /// One message per fold without a model.
    pub failed_folds: Vec<String>,
    /// Rows without a score.
    pub unscored: usize,
    /// Scored rows without a pooled q-value.
    pub no_q: usize,
}

/// The fitted fold models with their q-values, kept to score rows that took
/// no part in the fit.
pub struct TransferFdr {
    pub result: FdrResult,
    pub counts: FdrCounts,
    seed: u64,
    models: Vec<Option<GbdtPeakModel>>,
    /// Per fold, `(score, pooled q)` of its scored rows, ascending by score.
    tables: Vec<Vec<(f32, f64)>>,
}

/// The fold of a precursor (splitmix64 of its id).
pub fn fold_of(precursor: u32, seed: u64) -> usize {
    let mut z = (seed ^ (precursor as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15))
        .wrapping_add(0x9e37_79b9_7f4a_7c15);
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^= z >> 31;
    (z % N_FOLDS as u64) as usize
}

/// Cross-validated q-values of transfer targets and decoys.
#[allow(dead_code)]
pub fn transfer_qvalues(rows: &[FdrRow], seed: u64) -> FdrResult {
    transfer_fdr(rows, seed).result
}

/// Trains the classifier of `fold` on the rows of the other folds.
fn train_fold(
    rows: &[FdrRow],
    folds: &[usize],
    nf: usize,
    fold: usize,
    seed: u64,
) -> Result<GbdtPeakModel, String> {
    let mut x = Vec::new();
    let mut y = Vec::new();
    let mut groups = Vec::new();
    for (r, &f) in rows.iter().zip(folds) {
        if f == fold {
            continue;
        }
        x.extend_from_slice(&r.x);
        y.push(u8::from(!r.is_decoy));
        groups.push(r.precursor);
    }
    let positives = y.iter().filter(|&&v| v == 1).count();
    if positives == 0 || positives == y.len() {
        return Err(format!(
            "fold {fold}: training rows are single-class ({positives} targets of {})",
            y.len()
        ));
    }
    let ds = Dataset {
        x,
        y,
        groups,
        n_features: nf,
    };
    train_gbdt(
        &ds,
        &TrainParams::default(),
        seed.wrapping_add(fold as u64 + 1),
    )
    .map_err(|e| format!("fold {fold}: {e:?}"))
}

/// [`transfer_qvalues`] keeping the fold models.
///
/// # Panics
/// When the rows do not all have the same number of features.
pub fn transfer_fdr(rows: &[FdrRow], seed: u64) -> TransferFdr {
    let n = rows.len();
    let nf = rows.first().map_or(0, |r| r.x.len());
    assert!(
        rows.iter().all(|r| r.x.len() == nf),
        "transfer FDR rows have different feature counts"
    );
    let folds: Vec<usize> = rows.iter().map(|r| fold_of(r.precursor, seed)).collect();
    let mut result = FdrResult {
        q: vec![None; n],
        q_mass_only: vec![None; n],
        q_rt_only: vec![None; n],
        score: vec![f32::NAN; n],
    };
    let mut counts = FdrCounts::default();
    let fitted: Vec<Result<GbdtPeakModel, String>> = (0..N_FOLDS)
        .into_par_iter()
        .map(|f| train_fold(rows, &folds, nf, f, seed))
        .collect();
    let mut models = Vec::with_capacity(N_FOLDS);
    let mut tables = Vec::with_capacity(N_FOLDS);
    for (fold, fit) in fitted.into_iter().enumerate() {
        let members: Vec<usize> = (0..n).filter(|&i| folds[i] == fold).collect();
        let model = match fit {
            Ok(m) => m,
            Err(e) => {
                counts.failed_folds.push(e);
                counts.unscored += members.len();
                models.push(None);
                tables.push(Vec::new());
                continue;
            }
        };
        let xs: Vec<&[f32]> = members.iter().map(|&i| rows[i].x.as_slice()).collect();
        let mut scores = vec![0.0f32; members.len()];
        model.predict_value_batch(&xs, &mut scores);
        for (&i, &s) in members.iter().zip(&scores) {
            result.score[i] = s;
        }
        let estimate = |keep: &dyn Fn(&FdrRow) -> bool, q: &mut [Option<f64>]| {
            let items: Vec<(usize, f32, bool)> = members
                .iter()
                .filter(|&&i| result.score[i].is_finite() && keep(&rows[i]))
                .map(|&i| (i, result.score[i], rows[i].is_decoy))
                .collect();
            kdecoy_q(&items, q);
        };
        estimate(&|_| true, &mut result.q);
        estimate(&|r| !r.is_decoy || r.family == 1, &mut result.q_mass_only);
        estimate(&|r| !r.is_decoy || r.family == 2, &mut result.q_rt_only);
        let mut table: Vec<(f32, f64)> = members
            .iter()
            .filter_map(|&i| result.q[i].map(|q| (result.score[i], q)))
            .collect();
        table.sort_by(|a, b| a.0.total_cmp(&b.0).then(b.1.total_cmp(&a.1)));
        counts.no_q += members.len() - table.len();
        tables.push(table);
        models.push(Some(model));
    }
    TransferFdr {
        result,
        counts,
        seed,
        models,
        tables,
    }
}

impl TransferFdr {
    /// The q-value of a row outside the fit (a pseudo-transfer), scored by
    /// the model of its precursor's fold: the pooled q of the fold's scored
    /// row with the highest score at or below it, or of the fold's lowest row
    /// when it scores below all of them. `None` when the fold has no model or
    /// no q-values.
    pub fn pseudo_q(&self, precursor: u32, x: &[f32]) -> Option<f64> {
        let fold = fold_of(precursor, self.seed);
        let model = self.models.get(fold)?.as_ref()?;
        let table = &self.tables[fold];
        let first = table.first()?;
        let mut score = [0.0f32];
        model.predict_value_batch(&[x], &mut score);
        let at = table.partition_point(|&(s, _)| s <= score[0]);
        Some(if at == 0 { first.1 } else { table[at - 1].1 })
    }
}

/// K-decoy q-values of `(row, score, is_decoy)` into `q[row]`, with `K` the
/// decoy-to-target ratio of `items` and `FDR = (D + 1) / (K * max(T, 1))`. Leaves `q` untouched when `items` has no
/// targets or no decoys.
fn kdecoy_q(items: &[(usize, f32, bool)], q: &mut [Option<f64>]) {
    let targets = items.iter().filter(|it| !it.2).count();
    let decoys = items.len() - targets;
    if targets == 0 || decoys == 0 {
        return;
    }
    let k = decoys as f64 / targets as f64;
    let mut order: Vec<usize> = (0..items.len()).collect();
    order.sort_by(|&a, &b| items[b].1.total_cmp(&items[a].1));
    let mut fdr = vec![0.0f64; order.len()];
    let (mut t, mut d) = (0usize, 0usize);
    let mut start = 0;
    while start < order.len() {
        let score = items[order[start]].1;
        let mut end = start;
        while end < order.len() && (end == start || items[order[end]].1 == score) {
            if items[order[end]].2 {
                d += 1;
            } else {
                t += 1;
            }
            end += 1;
        }
        let v = ((d + 1) as f64 / (k * t.max(1) as f64)).min(1.0);
        fdr[start..end].fill(v);
        start = end;
    }
    let mut running = 1.0f64;
    for (pos, &j) in order.iter().enumerate().rev() {
        running = running.min(fdr[pos]);
        q[items[j].0] = Some(running);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Deterministic standard normal draws (splitmix64 + Box-Muller).
    struct Normal(u64);

    impl Normal {
        fn unit(&mut self) -> f64 {
            self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
            let mut z = self.0;
            z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
            z ^= z >> 31;
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

    /// `n` precursors, each with one target (mean `target_mean`) and one
    /// decoy per listed family (mean 0).
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
        // 10,000 experiments of 20 targets and 40 decoys, every score drawn
        // from one distribution: any accepted target is a false discovery,
        // so the share of experiments that accept one is their mean FDP.
        let mut g = Normal(77);
        let n_exp = 10_000;
        let mut with_false = 0usize;
        let mut q = vec![None; 60];
        for _ in 0..n_exp {
            let items: Vec<(usize, f32, bool)> = (0..60).map(|i| (i, g.draw(), i >= 20)).collect();
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
        // 0.01 plus about five binomial standard errors at n = 10,000.
        assert!(
            frac <= 0.015,
            "{frac} of the null experiments accept a target"
        );
    }

    #[test]
    fn mixture_fdp_stays_near_the_estimate() {
        // 30% of the targets are present (shifted on four features), the rest
        // are drawn like the decoys.
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

        // Fold 0 holds targets only: its rows have no decoys to estimate from.
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
        // Four decoys per target at the same scores as the targets: FDR = 1
        // only when K is the observed ratio 4, not an assumed 2.
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
        let scored: Vec<(usize, f32, bool)> = rows
            .iter()
            .enumerate()
            .map(|(i, r)| (i, r.x[0], r.is_decoy))
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
        // A pseudo row identical to a scored target gets that target's q.
        for (i, row) in rows
            .iter()
            .enumerate()
            .filter(|(_, r)| !r.is_decoy)
            .take(50)
        {
            assert_eq!(fdr.pseudo_q(row.precursor, &row.x), fdr.result.q[i]);
        }
        // Far below every row: the fold's largest q.
        let low = vec![-50.0f32; TRANSFER_FEATURES.len()];
        assert!(fdr.pseudo_q(0, &low).is_some_and(|q| q > 0.1));
    }

    #[test]
    fn results_are_deterministic_for_a_seed() {
        let rows = rows(300, 1.5, &[1, 2], 13);
        assert_eq!(transfer_qvalues(&rows, 9), transfer_qvalues(&rows, 9));
    }
}
