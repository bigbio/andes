//! Logistic gradient boosting trainer → [`GbdtPeakModel`].
//!
//! Implements Newton-step GBDT with logistic loss, quantile binning,
//! group-disjoint validation splitting, negative undersampling, early stopping,
//! and PAVA isotonic calibration on the hold-out validation set.

use scoring_crate::gbdt_eval::GbdtPeakModel;

use super::isotonic::pava;
use super::tree::{fit_tree, TreeParams};

// ---------------------------------------------------------------------------
// Public types
// ---------------------------------------------------------------------------

/// Row-major continuous feature matrix `x` (n × n_features), binary labels `y`
/// (0/1), and a `groups` id per row (run+peptide) for leakage-free validation
/// splitting.
pub struct Dataset {
    pub x: Vec<f32>,
    pub y: Vec<u8>,
    pub groups: Vec<u32>,
    pub n_features: usize,
}

/// Hyper-parameters for gradient boosting.
pub struct TrainParams {
    pub n_rounds: u32,
    pub learning_rate: f32,
    pub max_depth: u32,
    pub n_bins: usize,
    pub lambda: f32,
    pub min_hessian: f32,
    /// Number of negatives kept per positive during undersampling.
    pub neg_pos_ratio: f32,
    /// Fraction of groups assigned to validation.
    pub val_fraction: f32,
    /// Stop if validation logloss does not improve for this many rounds.
    pub early_stop_rounds: u32,
}

impl Default for TrainParams {
    fn default() -> Self {
        Self {
            n_rounds: 300,
            learning_rate: 0.05,
            max_depth: 6,
            n_bins: 64,
            lambda: 1.0,
            min_hessian: 5.0,
            neg_pos_ratio: 4.0,
            val_fraction: 0.2,
            early_stop_rounds: 30,
        }
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Splitmix64 — one step, returns a pseudo-random u64.
#[inline]
fn splitmix64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

/// Map a u64 into [0.0, 1.0).
#[inline]
fn to_unit(v: u64) -> f64 {
    (v >> 11) as f64 / (1u64 << 53) as f64
}

/// Sigmoid of x, numerically stable.
#[inline]
fn sigmoid(x: f32) -> f32 {
    1.0 / (1.0 + (-x).exp())
}

/// Log-loss over a validation set with raw scores (before sigmoid).
fn logloss(raw_val: &[f32], val_y: &[u8]) -> f32 {
    let mut sum = 0.0f64;
    for (&r, &y) in raw_val.iter().zip(val_y) {
        let p = sigmoid(r).clamp(1e-7, 1.0 - 1e-7) as f64;
        let label = y as f64;
        sum -= label * p.ln() + (1.0 - label) * (1.0 - p).ln();
    }
    if val_y.is_empty() { f32::INFINITY } else { (sum / val_y.len() as f64) as f32 }
}

/// ROC AUC of `scores` vs binary `labels` (Mann-Whitney U; ties get average
/// rank). Rank-based, so AUC of the raw GBDT scores equals AUC of the
/// calibrated P(signal). Returns 0.5 when either class is empty. This is the
/// offline gate metric: AUC ≈ 0.5 ⇒ the peptide-agnostic peak model has no
/// discriminative signal; higher ⇒ worth wiring as an additive feature.
fn auc(scores: &[f32], labels: &[u8]) -> f64 {
    let n = scores.len();
    let n_pos = labels.iter().filter(|&&y| y == 1).count();
    let n_neg = n - n_pos;
    if n_pos == 0 || n_neg == 0 {
        return 0.5;
    }
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&a, &b| scores[a].partial_cmp(&scores[b]).unwrap());
    let mut ranks = vec![0.0f64; n];
    let mut i = 0;
    while i < n {
        let mut j = i;
        while j + 1 < n && scores[idx[j + 1]] == scores[idx[i]] {
            j += 1;
        }
        let avg_rank = ((i + 1) + (j + 1)) as f64 / 2.0; // 1-based average rank for ties
        for &row in &idx[i..=j] {
            ranks[row] = avg_rank;
        }
        i = j + 1;
    }
    let sum_pos_ranks: f64 = (0..n).filter(|&r| labels[r] == 1).map(|r| ranks[r]).sum();
    (sum_pos_ranks - (n_pos * (n_pos + 1)) as f64 / 2.0) / (n_pos as f64 * n_neg as f64)
}

/// Build quantile bin upper edges for one feature column.
///
/// Returns up to `n_bins` distinct upper edges at evenly-spaced quantiles.
/// The result is ascending and deduplicated. If there are fewer distinct values
/// than `n_bins`, all distinct values are returned (each as its own upper edge).
fn compute_bin_uppers(col: &[f32], n_bins: usize) -> Vec<f32> {
    debug_assert!(n_bins <= 256, "n_bins must fit in u8");
    let mut sorted: Vec<f32> = col.iter().copied().filter(|v| v.is_finite()).collect();
    sorted.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
    if sorted.is_empty() {
        return vec![0.0];
    }
    let n = sorted.len();
    // Pick up to n_bins quantile positions (using 0-indexed boundary of each bin).
    let mut edges: Vec<f32> = Vec::with_capacity(n_bins);
    for b in 1..=n_bins {
        // Upper quantile of bin b: index = ceil(b * n / n_bins) - 1 (clamped).
        let idx = (b * n).div_ceil(n_bins).min(n) - 1;
        let v = sorted[idx];
        if edges.last().copied() != Some(v) {
            edges.push(v);
        }
    }
    if edges.is_empty() {
        edges.push(sorted[n - 1]);
    }
    edges
}

/// Map a continuous value to its bin index given the precomputed upper edges.
///
/// Bin = first index whose upper edge is >= v (i.e., `partition_point(|&u| u < v)`).
/// This is the EXACT inverse of how `fit_tree` uses thresholds: it stores
/// `threshold = bin_uppers[feat][bin]` and routes `x <= threshold` left, so a
/// value at or below the upper edge bins to `<= bin` ↔ goes left. Consistent.
///
/// Clamped to `uppers.len() - 1` so result always fits in u8 (since n_bins ≤ 256).
#[inline]
fn bin_of(v: f32, uppers: &[f32]) -> u8 {
    let b = uppers.partition_point(|&u| u < v);
    b.min(uppers.len() - 1) as u8
}

/// Fisher-Yates shuffle of indices (deterministic, seeded via splitmix64).
fn shuffle_indices(indices: &mut [usize], rng: &mut u64) {
    let n = indices.len();
    for i in (1..n).rev() {
        let j = (splitmix64(rng) as usize) % (i + 1);
        indices.swap(i, j);
    }
}

// ---------------------------------------------------------------------------
// Core training function
// ---------------------------------------------------------------------------

/// Fit a GBDT classifier with logistic loss and return a calibrated
/// [`GbdtPeakModel`].
///
/// All randomness is seeded from `seed`; two calls with identical arguments
/// and seed produce byte-for-byte identical models.
pub fn train_gbdt(ds: &Dataset, p: &TrainParams, seed: u64) -> GbdtPeakModel {
    let n_rows = ds.y.len();
    let nf = ds.n_features;

    // --- Guard: degenerate input → empty model ---------------------------------
    let n_pos = ds.y.iter().filter(|&&y| y == 1).count();
    if n_pos == 0 || n_rows == 0 || nf == 0 {
        return GbdtPeakModel {
            n_features: nf as u32,
            apply_sigmoid: true,
            trees: vec![],
            iso_x: vec![],
            iso_y: vec![],
        };
    }

    let n_bins = p.n_bins.min(256);

    // --- Step 1: Quantile-bin each feature ------------------------------------
    // bin_uppers[f] = ascending upper edges for feature f
    let bin_uppers: Vec<Vec<f32>> = (0..nf)
        .map(|f| {
            let col: Vec<f32> = (0..n_rows).map(|r| ds.x[r * nf + f]).collect();
            compute_bin_uppers(&col, n_bins)
        })
        .collect();

    // Build binned matrix (row-major, n_rows × nf, values in 0..n_bins as u8)
    let mut binned_all = vec![0u8; n_rows * nf];
    for r in 0..n_rows {
        for f in 0..nf {
            binned_all[r * nf + f] = bin_of(ds.x[r * nf + f], &bin_uppers[f]);
        }
    }

    // --- Step 2: Group-disjoint validation split (seeded) ---------------------
    // Collect distinct group ids in sorted order for determinism.
    let mut distinct_groups: Vec<u32> = ds.groups.clone();
    distinct_groups.sort_unstable();
    distinct_groups.dedup();

    // Hash each group id with the seed using splitmix64 to decide train/val.
    // A group is validation if hash → unit < val_fraction.
    let val_set: std::collections::BTreeSet<u32> = distinct_groups
        .iter()
        .filter(|&&g| {
            let mut state = seed ^ (g as u64).wrapping_mul(0x9e3779b97f4a7c15);
            let h = splitmix64(&mut state);
            to_unit(h) < p.val_fraction as f64
        })
        .copied()
        .collect();

    let mut train_rows: Vec<usize> = Vec::new();
    let mut val_rows: Vec<usize> = Vec::new();
    for r in 0..n_rows {
        if val_set.contains(&ds.groups[r]) {
            val_rows.push(r);
        } else {
            train_rows.push(r);
        }
    }

    // Guard: if no validation rows (e.g. empty val fraction), return empty model.
    if val_rows.is_empty() {
        return GbdtPeakModel {
            n_features: nf as u32,
            apply_sigmoid: true,
            trees: vec![],
            iso_x: vec![],
            iso_y: vec![],
        };
    }

    // --- Step 3: Negative undersample TRAIN rows (seeded) ---------------------
    let train_pos: Vec<usize> =
        train_rows.iter().copied().filter(|&r| ds.y[r] == 1).collect();
    let mut train_neg: Vec<usize> =
        train_rows.iter().copied().filter(|&r| ds.y[r] == 0).collect();

    let n_pos_train = train_pos.len();
    // Guard: if no positives in training set, return empty model.
    if n_pos_train == 0 {
        return GbdtPeakModel {
            n_features: nf as u32,
            apply_sigmoid: true,
            trees: vec![],
            iso_x: vec![],
            iso_y: vec![],
        };
    }

    let n_neg_keep =
        (p.neg_pos_ratio * n_pos_train as f32).ceil() as usize;
    let mut rng = seed;
    if train_neg.len() > n_neg_keep {
        shuffle_indices(&mut train_neg, &mut rng);
        train_neg.truncate(n_neg_keep);
        train_neg.sort_unstable(); // keep deterministic order
    }

    let mut active_train: Vec<usize> = train_pos;
    active_train.extend_from_slice(&train_neg);
    active_train.sort_unstable(); // deterministic order

    let n_active = active_train.len();

    // Build per-row compact maps for fast access.
    // active_train[i] → original row index.
    // We need: binned subset, x subset for tree eval, y subset.

    // Compact binned array for the active training rows.
    let mut train_binned = vec![0u8; n_active * nf];
    let mut train_y = vec![0u8; n_active];
    for (i, &r) in active_train.iter().enumerate() {
        train_y[i] = ds.y[r];
        for f in 0..nf {
            train_binned[i * nf + f] = binned_all[r * nf + f];
        }
    }

    let val_n = val_rows.len();
    let val_y: Vec<u8> = val_rows.iter().map(|&r| ds.y[r]).collect();

    // --- Step 4: Boosting loop ------------------------------------------------
    // raw scores in logit space; init 0 (sigmoid(0) = 0.5).
    let mut raw_train = vec![0.0f32; n_active];
    let mut raw_val = vec![0.0f32; val_n];

    let tree_params = TreeParams {
        max_depth: p.max_depth,
        n_bins,
        lambda: p.lambda,
        min_hessian: p.min_hessian,
        n_features: nf,
    };

    let mut trees: Vec<scoring_crate::gbdt_eval::Tree> = Vec::new();
    let mut best_loss = f32::INFINITY;
    let mut best_round: usize = 0;
    let mut no_improve = 0u32;

    for _round in 0..p.n_rounds {
        // Compute gradients and hessians for active training rows.
        let mut grad = vec![0.0f32; n_active];
        let mut hess = vec![0.0f32; n_active];
        for i in 0..n_active {
            let pi = sigmoid(raw_train[i]);
            grad[i] = pi - train_y[i] as f32;
            hess[i] = (pi * (1.0 - pi)).max(1e-6);
        }

        // Fit one tree on the active binned training data.
        let mut tree = fit_tree(&train_binned, &grad, &hess, &tree_params, &bin_uppers);

        // Apply learning rate shrinkage to all leaf values.
        // Scaling ALL values is safe: internal nodes have value 0 and are not
        // reached by eval; leaves carry the actual weight.
        for v in tree.value.iter_mut() {
            *v *= p.learning_rate;
        }

        // Update raw scores: training rows (use compact index i, continuous x).
        for i in 0..n_active {
            let r = active_train[i];
            let row = &ds.x[r * nf..(r + 1) * nf];
            raw_train[i] += tree.eval(row);
        }

        // Update raw scores: validation rows (continuous x).
        for (j, &r) in val_rows.iter().enumerate() {
            let row = &ds.x[r * nf..(r + 1) * nf];
            raw_val[j] += tree.eval(row);
        }

        trees.push(tree);

        // Evaluate validation logloss.
        let vloss = logloss(&raw_val, &val_y);
        if vloss < best_loss - 1e-8 {
            best_loss = vloss;
            best_round = trees.len(); // number of trees at best
            no_improve = 0;
        } else {
            no_improve += 1;
            if no_improve >= p.early_stop_rounds {
                break;
            }
        }
    }

    // Truncate to best round.
    trees.truncate(best_round);

    // --- Step 5: PAVA calibration on validation set ---------------------------
    // Recompute raw_val from the truncated tree set for correctness.
    let mut raw_val_final = vec![0.0f32; val_n];
    for (j, &r) in val_rows.iter().enumerate() {
        let row = &ds.x[r * nf..(r + 1) * nf];
        for t in &trees {
            raw_val_final[j] += t.eval(row);
        }
    }

    // OFFLINE GATE METRIC: held-out AUC of P(signal) on the group-disjoint
    // validation set. This is the review's free gate — AUC ≈ 0.5 means the
    // peptide-agnostic peak model does not discriminate signal from noise and
    // is not worth wiring (jump to the per-ion v3 design instead).
    let n_pos_val = val_y.iter().filter(|&&y| y == 1).count();
    let val_auc = auc(&raw_val_final, &val_y);
    eprintln!(
        "train-gbdt: validation AUC = {val_auc:.4}  (n_val={val_n}, pos={n_pos_val}, neg={})",
        val_n - n_pos_val
    );

    // Build (p_sigmoid, y) pairs sorted by p_sigmoid.
    let mut cal_pairs: Vec<(f32, f32)> = raw_val_final
        .iter()
        .zip(val_y.iter())
        .map(|(&r, &y)| (sigmoid(r), y as f32))
        .collect();
    cal_pairs.sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    let sorted_p: Vec<f32> = cal_pairs.iter().map(|&(p, _)| p).collect();
    let sorted_y: Vec<f32> = cal_pairs.iter().map(|&(_, y)| y).collect();
    let (iso_x, iso_y) = pava(&sorted_p, &sorted_y);

    // --- Return ---------------------------------------------------------------
    GbdtPeakModel {
        n_features: nf as u32,
        apply_sigmoid: true,
        trees,
        iso_x,
        iso_y,
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn auc_perfect_reversed_tied_and_empty() {
        // Perfect separation: all positives score above all negatives → 1.0.
        let scores = [0.1_f32, 0.2, 0.8, 0.9];
        let labels = [0u8, 0, 1, 1];
        assert!((auc(&scores, &labels) - 1.0).abs() < 1e-9);
        // Reversed labels (positives score lowest) → 0.0.
        let labels_rev = [1u8, 1, 0, 0];
        assert!((auc(&scores, &labels_rev) - 0.0).abs() < 1e-9);
        // All-tied scores → 0.5 (no separation).
        let tied = [0.5_f32, 0.5, 0.5, 0.5];
        assert!((auc(&tied, &labels) - 0.5).abs() < 1e-9);
        // One class empty → 0.5 by convention.
        assert!((auc(&scores, &[0u8, 0, 0, 0]) - 0.5).abs() < 1e-9);
    }

    fn lcg(state: &mut u64) -> f32 {
        *state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((*state >> 33) as f32) / ((1u64 << 31) as f32)
    }

    fn separable(nf: usize, n: usize, seed: u64) -> Dataset {
        let mut s = seed;
        let mut x = vec![0.0f32; n * nf];
        let mut y = vec![0u8; n];
        let mut groups = vec![0u32; n];
        for r in 0..n {
            for f in 0..nf {
                x[r * nf + f] = lcg(&mut s);
            }
            let label = if x[r * nf] + 0.2 * lcg(&mut s) > 0.6 { 1u8 } else { 0 };
            y[r] = label;
            groups[r] = (r as u32) % 50;
        }
        Dataset { x, y, groups, n_features: nf }
    }

    #[test]
    fn separable_data_trains_and_separates() {
        let ds = separable(6, 4000, 0x1234);
        let model = train_gbdt(&ds, &TrainParams::default(), 0xC0FFEE);
        assert!(model.apply_sigmoid);
        assert_eq!(model.n_features as usize, 6);
        assert!(!model.trees.is_empty(), "should fit at least one tree");
        let mut hi = vec![0.0f32; 6];
        hi[0] = 0.95;
        let mut lo = vec![0.0f32; 6];
        lo[0] = 0.05;
        assert!(
            model.predict_logit(&hi) > model.predict_logit(&lo) + 1.0,
            "hi={} lo={}",
            model.predict_logit(&hi),
            model.predict_logit(&lo)
        );
    }

    #[test]
    fn deterministic_for_same_seed() {
        let ds = separable(5, 1500, 7);
        let a = train_gbdt(&ds, &TrainParams::default(), 99);
        let b = train_gbdt(&ds, &TrainParams::default(), 99);
        assert_eq!(
            a.to_bytes(),
            b.to_bytes(),
            "same seed must yield identical model"
        );
    }

    #[test]
    fn no_positives_returns_empty_model() {
        let ds = Dataset {
            x: vec![0.1f32; 10 * 3],
            y: vec![0u8; 10],
            groups: (0..10).collect(),
            n_features: 3,
        };
        let m = train_gbdt(&ds, &TrainParams::default(), 1);
        assert!(m.trees.is_empty());
    }
}
