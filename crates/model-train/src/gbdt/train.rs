use model::hash::{splitmix64_next, SPLITMIX_GOLDEN};
use scoring_crate::gbdt_eval::GbdtPeakModel;

use super::isotonic::pava;
use super::tree::{fit_tree, fit_tree_with_leaves, TreeParams};
use crate::TrainError;

pub const MIN_TRAIN_ROWS: usize = 100;

pub const MIN_VAL_AUC: f64 = 0.55;

pub const MIN_VAL_PEARSON: f64 = 0.10;
pub const MIN_VAL_R2: f64 = 0.0;

pub struct Dataset {
    pub x: Vec<f32>,
    pub y: Vec<u8>,
    pub groups: Vec<u32>,
    pub n_features: usize,
    pub feature_names: Vec<String>,
}

pub struct TrainParams {
    pub n_rounds: u32,
    pub learning_rate: f32,
    pub max_depth: u32,
    pub n_bins: usize,
    pub lambda: f32,
    pub min_hessian: f32,
    pub neg_pos_ratio: f32,
    pub val_fraction: f32,
    pub early_stop_rounds: u32,
    pub pairwise: bool,

    pub allow_degenerate: bool,
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
            pairwise: false,
            allow_degenerate: false,
        }
    }
}

fn gate_or_fallback(
    allow_degenerate: bool,
    reason: String,
    fallback: impl FnOnce() -> GbdtPeakModel,
) -> Result<GbdtPeakModel, TrainError> {
    if allow_degenerate {
        log::warn!("gbdt quality gate failed ({reason}); returning degenerate model (--allow-degenerate-model)");
        Ok(fallback())
    } else {
        Err(TrainError::QualityGate(reason))
    }
}

#[inline]
fn to_unit(v: u64) -> f64 {
    (v >> 11) as f64 / (1u64 << 53) as f64
}

#[inline]
fn sigmoid(x: f32) -> f32 {
    1.0 / (1.0 + (-x).exp())
}

fn logloss(raw_val: &[f32], val_y: &[u8]) -> f32 {
    let mut sum = 0.0f64;
    for (&r, &y) in raw_val.iter().zip(val_y) {
        let p = sigmoid(r).clamp(1e-7, 1.0 - 1e-7) as f64;
        let label = y as f64;
        sum -= label * p.ln() + (1.0 - label) * (1.0 - p).ln();
    }
    if val_y.is_empty() {
        f32::INFINITY
    } else {
        (sum / val_y.len() as f64) as f32
    }
}

fn auc(scores: &[f32], labels: &[u8]) -> f64 {
    let n = scores.len();
    let n_pos = labels.iter().filter(|&&y| y == 1).count();
    let n_neg = n - n_pos;
    if n_pos == 0 || n_neg == 0 {
        return 0.5;
    }
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&a, &b| scores[a].total_cmp(&scores[b]));
    let mut ranks = vec![0.0f64; n];
    let mut i = 0;
    while i < n {
        let mut j = i;
        while j + 1 < n && scores[idx[j + 1]] == scores[idx[i]] {
            j += 1;
        }
        let avg_rank = ((i + 1) + (j + 1)) as f64 / 2.0;
        for &row in &idx[i..=j] {
            ranks[row] = avg_rank;
        }
        i = j + 1;
    }
    let sum_pos_ranks: f64 = (0..n).filter(|&r| labels[r] == 1).map(|r| ranks[r]).sum();
    (sum_pos_ranks - (n_pos * (n_pos + 1)) as f64 / 2.0) / (n_pos as f64 * n_neg as f64)
}

fn grouped_neg_auc(scores: &[f32], labels: &[u8], groups: &[u32]) -> Option<f32> {
    use std::collections::BTreeMap;
    let mut by_group: BTreeMap<u32, (Vec<f32>, Vec<u8>)> = BTreeMap::new();
    for i in 0..scores.len() {
        let e = by_group.entry(groups[i]).or_default();
        e.0.push(scores[i]);
        e.1.push(labels[i]);
    }
    let mut sum = 0.0f64;
    let mut n_groups = 0usize;
    for (_g, (s, y)) in by_group {
        let has_pos = y.contains(&1);
        let has_neg = y.contains(&0);
        if has_pos && has_neg {
            sum += auc(&s, &y);
            n_groups += 1;
        }
    }
    if n_groups == 0 {
        None
    } else {
        Some(-(sum / n_groups as f64) as f32)
    }
}

fn compute_bin_uppers(col: &[f32], n_bins: usize) -> Vec<f32> {
    debug_assert!(n_bins <= 256, "n_bins must fit in u8");
    let mut sorted: Vec<f32> = col.iter().copied().filter(|v| v.is_finite()).collect();
    sorted.sort_unstable_by(f32::total_cmp);
    if sorted.is_empty() {
        return vec![0.0];
    }
    let n = sorted.len();
    let mut edges: Vec<f32> = Vec::with_capacity(n_bins);
    for b in 1..=n_bins {
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

#[inline]
fn bin_of(v: f32, uppers: &[f32]) -> u8 {
    let b = uppers.partition_point(|&u| u < v);
    b.min(uppers.len() - 1) as u8
}

fn shuffle_indices(indices: &mut [usize], rng: &mut u64) {
    let n = indices.len();
    for i in (1..n).rev() {
        let j = (splitmix64_next(rng) as usize) % (i + 1);
        indices.swap(i, j);
    }
}

fn bin_uppers_per_feature(x: &[f32], n_rows: usize, nf: usize, n_bins: usize) -> Vec<Vec<f32>> {
    (0..nf)
        .map(|f| {
            let col: Vec<f32> = (0..n_rows).map(|r| x[r * nf + f]).collect();
            compute_bin_uppers(&col, n_bins)
        })
        .collect()
}

fn bin_rows(x: &[f32], nf: usize, rows: &[usize], bin_uppers: &[Vec<f32>]) -> Vec<u8> {
    let mut binned = vec![0u8; rows.len() * nf];
    for (i, &r) in rows.iter().enumerate() {
        for f in 0..nf {
            binned[i * nf + f] = bin_of(x[r * nf + f], &bin_uppers[f]);
        }
    }
    binned
}

fn split_rows_by_group(groups: &[u32], seed: u64, val_fraction: f32) -> (Vec<usize>, Vec<usize>) {
    let is_val = |g: u32| {
        let mut state = seed ^ (g as u64).wrapping_mul(SPLITMIX_GOLDEN);
        to_unit(splitmix64_next(&mut state)) < val_fraction as f64
    };
    (0..groups.len()).partition(|&r| !is_val(groups[r]))
}

fn gbdt_model(
    feature_names: &[String],
    apply_sigmoid: bool,
    trees: Vec<scoring_crate::gbdt_eval::Tree>,
) -> GbdtPeakModel {
    GbdtPeakModel {
        packed: Default::default(),
        n_features: feature_names.len() as u32,
        feature_names: feature_names.to_vec(),
        apply_sigmoid,
        trees,
        iso_x: vec![],
        iso_y: vec![],
    }
}

pub fn train_gbdt(ds: &Dataset, p: &TrainParams, seed: u64) -> Result<GbdtPeakModel, TrainError> {
    let n_rows = ds.y.len();
    let nf = ds.n_features;
    if ds.feature_names.len() != nf {
        return Err(TrainError::Other(format!(
            "{} feature names for {nf} features",
            ds.feature_names.len()
        )));
    }

    let empty_classifier = || gbdt_model(&ds.feature_names, true, vec![]);
    let n_pos = ds.y.iter().filter(|&&y| y == 1).count();
    if nf == 0 {
        return gate_or_fallback(
            p.allow_degenerate,
            "classifier: zero features".into(),
            empty_classifier,
        );
    }
    if n_rows < MIN_TRAIN_ROWS {
        return gate_or_fallback(
            p.allow_degenerate,
            format!("classifier: {n_rows} rows < minimum {MIN_TRAIN_ROWS}"),
            empty_classifier,
        );
    }
    if n_pos == 0 || n_pos == n_rows {
        return gate_or_fallback(
            p.allow_degenerate,
            "classifier: training labels are single-class (need both 0 and 1)".into(),
            empty_classifier,
        );
    }

    let n_bins = p.n_bins.min(256);

    let bin_uppers = bin_uppers_per_feature(&ds.x, n_rows, nf, n_bins);

    let (train_rows, val_rows) = split_rows_by_group(&ds.groups, seed, p.val_fraction);

    if val_rows.is_empty() {
        return gate_or_fallback(
            p.allow_degenerate,
            "classifier: validation split is empty (cannot evaluate held-out AUC)".into(),
            empty_classifier,
        );
    }

    let train_pos: Vec<usize> = train_rows
        .iter()
        .copied()
        .filter(|&r| ds.y[r] == 1)
        .collect();
    let mut train_neg: Vec<usize> = train_rows
        .iter()
        .copied()
        .filter(|&r| ds.y[r] == 0)
        .collect();

    let n_pos_train = train_pos.len();
    if n_pos_train == 0 {
        return gate_or_fallback(
            p.allow_degenerate,
            "classifier: training split has no positive rows".into(),
            empty_classifier,
        );
    }

    let n_neg_keep = (p.neg_pos_ratio * n_pos_train as f32).ceil() as usize;
    let mut rng = seed;
    if train_neg.len() > n_neg_keep {
        shuffle_indices(&mut train_neg, &mut rng);
        train_neg.truncate(n_neg_keep);
        train_neg.sort_unstable();
    }

    let mut active_train: Vec<usize> = train_pos;
    active_train.extend_from_slice(&train_neg);
    active_train.sort_unstable();

    let n_active = active_train.len();
    if n_active < MIN_TRAIN_ROWS {
        return gate_or_fallback(
            p.allow_degenerate,
            format!("classifier: {n_active} active training rows < minimum {MIN_TRAIN_ROWS}"),
            empty_classifier,
        );
    }

    let train_binned = bin_rows(&ds.x, nf, &active_train, &bin_uppers);
    let train_y: Vec<u8> = active_train.iter().map(|&r| ds.y[r]).collect();

    let val_n = val_rows.len();
    let val_y: Vec<u8> = val_rows.iter().map(|&r| ds.y[r]).collect();
    let val_groups: Vec<u32> = val_rows.iter().map(|&r| ds.groups[r]).collect();

    let pair_groups: Vec<(Vec<usize>, Vec<usize>)> = if p.pairwise {
        let mut pos_map: std::collections::HashMap<u32, Vec<usize>> =
            std::collections::HashMap::new();
        let mut neg_map: std::collections::HashMap<u32, Vec<usize>> =
            std::collections::HashMap::new();
        for i in 0..n_active {
            let g = ds.groups[active_train[i]];
            if train_y[i] == 1 {
                pos_map.entry(g).or_default().push(i);
            } else {
                neg_map.entry(g).or_default().push(i);
            }
        }
        let mut gs: Vec<u32> = pos_map
            .keys()
            .filter(|k| neg_map.contains_key(k))
            .copied()
            .collect();
        gs.sort_unstable();
        gs.into_iter()
            .map(|g| (pos_map[&g].clone(), neg_map[&g].clone()))
            .collect()
    } else {
        Vec::new()
    };
    if p.pairwise && pair_groups.is_empty() {
        return Err(TrainError::QualityGate(
            "pairwise (RankNet) training requested but no training group contains \
             both a positive and a negative row — no rank pairs can be formed"
                .into(),
        ));
    }

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
    let mut best_raw_val = vec![0.0f32; val_n];
    let mut grad = vec![0.0f32; n_active];
    let mut hess = vec![0.0f32; n_active];

    for _round in 0..p.n_rounds {
        if p.pairwise {
            grad.fill(0.0);
            hess.fill(0.0);
            for (pos, neg) in &pair_groups {
                for &i in pos {
                    for &j in neg {
                        let rho = sigmoid(-(raw_train[i] - raw_train[j]));
                        grad[i] -= rho;
                        grad[j] += rho;
                        let h = rho * (1.0 - rho);
                        hess[i] += h;
                        hess[j] += h;
                    }
                }
            }
            for h in hess.iter_mut() {
                *h = h.max(1e-6);
            }
        } else {
            for i in 0..n_active {
                let pi = sigmoid(raw_train[i]);
                grad[i] = pi - train_y[i] as f32;
                hess[i] = (pi * (1.0 - pi)).max(1e-6);
            }
        }

        let (mut tree, leaf_of) =
            fit_tree_with_leaves(&train_binned, &grad, &hess, &tree_params, &bin_uppers);

        for v in tree.value.iter_mut() {
            *v *= p.learning_rate;
        }

        for (raw, &leaf) in raw_train.iter_mut().zip(&leaf_of) {
            *raw += tree.value[leaf as usize];
        }

        for (j, &r) in val_rows.iter().enumerate() {
            let row = &ds.x[r * nf..(r + 1) * nf];
            raw_val[j] += tree.eval(row);
        }

        trees.push(tree);

        let vloss_opt = if p.pairwise {
            grouped_neg_auc(&raw_val, &val_y, &val_groups)
        } else {
            Some(logloss(&raw_val, &val_y))
        };
        match vloss_opt {
            Some(vloss) => {
                if vloss < best_loss - 1e-8 {
                    best_loss = vloss;
                    best_round = trees.len();
                    best_raw_val.copy_from_slice(&raw_val);
                    no_improve = 0;
                } else {
                    no_improve += 1;
                    if no_improve >= p.early_stop_rounds {
                        break;
                    }
                }
            }
            None => {
                best_round = trees.len();
                best_raw_val.copy_from_slice(&raw_val);
            }
        }
    }

    trees.truncate(best_round);

    let raw_val_final = best_raw_val;

    let n_pos_val = val_y.iter().filter(|&&y| y == 1).count();
    let val_auc = auc(&raw_val_final, &val_y);
    log::debug!(
        target: "train-gbdt",
        "validation AUC = {val_auc:.4}  (n_val={val_n}, pos={n_pos_val}, neg={})",
        val_n - n_pos_val
    );

    if trees.is_empty() {
        return gate_or_fallback(
            p.allow_degenerate,
            "classifier: empty tree ensemble (no useful split found)".into(),
            empty_classifier,
        );
    }
    if n_pos_val == 0 || n_pos_val == val_n {
        return gate_or_fallback(
            p.allow_degenerate,
            "classifier: validation split is single-class (AUC undefined)".into(),
            empty_classifier,
        );
    }
    if !val_auc.is_finite() || val_auc < MIN_VAL_AUC {
        let trees_for_fallback = trees.clone();
        return gate_or_fallback(
            p.allow_degenerate,
            format!("classifier: held-out AUC {val_auc:.4} < minimum {MIN_VAL_AUC}"),
            move || gbdt_model(&ds.feature_names, true, trees_for_fallback),
        );
    }

    let mut cal_pairs: Vec<(f32, f32)> = raw_val_final
        .iter()
        .zip(val_y.iter())
        .map(|(&r, &y)| (sigmoid(r), y as f32))
        .collect();
    cal_pairs.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));

    let sorted_p: Vec<f32> = cal_pairs.iter().map(|&(p, _)| p).collect();
    let sorted_y: Vec<f32> = cal_pairs.iter().map(|&(_, y)| y).collect();
    let (iso_x, iso_y) = pava(&sorted_p, &sorted_y);

    Ok(GbdtPeakModel {
        iso_x,
        iso_y,
        ..gbdt_model(&ds.feature_names, true, trees)
    })
}

pub struct RegressionDataset {
    pub x: Vec<f32>,
    pub y: Vec<f32>,
    pub groups: Vec<u32>,
    pub n_features: usize,
    pub feature_names: Vec<String>,
}

pub fn pearson_r2(pred: &[f32], y: &[f32]) -> (f64, f64) {
    let n = pred.len().min(y.len());
    if n == 0 {
        return (0.0, 0.0);
    }
    let mean_y: f64 = y[..n].iter().map(|&v| v as f64).sum::<f64>() / n as f64;
    let mean_p: f64 = pred[..n].iter().map(|&v| v as f64).sum::<f64>() / n as f64;

    let mut cov = 0.0f64;
    let mut var_y = 0.0f64;
    let mut var_p = 0.0f64;
    let mut ss_res = 0.0f64;

    for i in 0..n {
        let dy = y[i] as f64 - mean_y;
        let dp = pred[i] as f64 - mean_p;
        cov += dy * dp;
        var_y += dy * dy;
        var_p += dp * dp;
        let res = pred[i] as f64 - y[i] as f64;
        ss_res += res * res;
    }

    if var_y < f64::EPSILON || var_p < f64::EPSILON {
        return (0.0, 0.0);
    }

    let r = cov / (var_y.sqrt() * var_p.sqrt());
    let r2 = 1.0 - ss_res / var_y;
    (r, r2)
}

fn mse(raw_val: &[f32], val_y: &[f32]) -> f32 {
    if val_y.is_empty() {
        return f32::INFINITY;
    }
    let sum: f64 = raw_val
        .iter()
        .zip(val_y)
        .map(|(&p, &y)| {
            let d = (p - y) as f64;
            d * d
        })
        .sum();
    (sum / val_y.len() as f64) as f32
}

pub fn train_gbdt_regression(
    ds: &RegressionDataset,
    p: &TrainParams,
    seed: u64,
) -> Result<GbdtPeakModel, TrainError> {
    let n_rows = ds.y.len();
    let nf = ds.n_features;
    if ds.feature_names.len() != nf {
        return Err(TrainError::Other(format!(
            "{} feature names for {nf} features",
            ds.feature_names.len()
        )));
    }

    let empty_regressor = || gbdt_model(&ds.feature_names, false, vec![]);
    if nf == 0 {
        return gate_or_fallback(
            p.allow_degenerate,
            "regressor: zero features".into(),
            empty_regressor,
        );
    }
    if n_rows < MIN_TRAIN_ROWS {
        return gate_or_fallback(
            p.allow_degenerate,
            format!("regressor: {n_rows} rows < minimum {MIN_TRAIN_ROWS}"),
            empty_regressor,
        );
    }

    let n_bins = p.n_bins.min(256);

    let bin_uppers = bin_uppers_per_feature(&ds.x, n_rows, nf, n_bins);

    let (train_rows, val_rows) = split_rows_by_group(&ds.groups, seed, p.val_fraction);

    if val_rows.is_empty() || train_rows.is_empty() {
        return gate_or_fallback(
            p.allow_degenerate,
            "regressor: empty train or validation split".into(),
            empty_regressor,
        );
    }

    let n_active = train_rows.len();
    if n_active < MIN_TRAIN_ROWS {
        return gate_or_fallback(
            p.allow_degenerate,
            format!("regressor: {n_active} training rows after split < minimum {MIN_TRAIN_ROWS}"),
            empty_regressor,
        );
    }
    let train_binned = bin_rows(&ds.x, nf, &train_rows, &bin_uppers);
    let train_y: Vec<f32> = train_rows.iter().map(|&r| ds.y[r]).collect();

    let val_n = val_rows.len();
    let val_y: Vec<f32> = val_rows.iter().map(|&r| ds.y[r]).collect();

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
        let grad: Vec<f32> = (0..n_active).map(|i| raw_train[i] - train_y[i]).collect();
        let hess: Vec<f32> = vec![1.0f32; n_active];

        let mut tree = fit_tree(&train_binned, &grad, &hess, &tree_params, &bin_uppers);

        for v in tree.value.iter_mut() {
            *v *= p.learning_rate;
        }

        for i in 0..n_active {
            let r = train_rows[i];
            let row = &ds.x[r * nf..(r + 1) * nf];
            raw_train[i] += tree.eval(row);
        }

        for (j, &r) in val_rows.iter().enumerate() {
            let row = &ds.x[r * nf..(r + 1) * nf];
            raw_val[j] += tree.eval(row);
        }

        trees.push(tree);

        let vloss = mse(&raw_val, &val_y);
        if vloss < best_loss - 1e-8 {
            best_loss = vloss;
            best_round = trees.len();
            no_improve = 0;
        } else {
            no_improve += 1;
            if no_improve >= p.early_stop_rounds {
                break;
            }
        }
    }

    trees.truncate(best_round);

    let mut raw_val_final = vec![0.0f32; val_n];
    for (j, &r) in val_rows.iter().enumerate() {
        let row = &ds.x[r * nf..(r + 1) * nf];
        for t in &trees {
            raw_val_final[j] += t.eval(row);
        }
    }

    let (r, r2) = pearson_r2(&raw_val_final, &val_y);
    let n = val_n;
    log::debug!(target: "train-gbdt", "regression val Pearson r = {r:.4}  R2 = {r2:.4}  (n_val={n})");

    if trees.is_empty() {
        return gate_or_fallback(
            p.allow_degenerate,
            "regressor: empty tree ensemble (no useful split found)".into(),
            empty_regressor,
        );
    }
    if !r.is_finite() || !r2.is_finite() || r < MIN_VAL_PEARSON || r2 <= MIN_VAL_R2 {
        let trees_for_fallback = trees.clone();
        return gate_or_fallback(
            p.allow_degenerate,
            format!(
                "regressor: held-out Pearson r = {r:.4} (min {MIN_VAL_PEARSON}), \
                 R2 = {r2:.4} (min {MIN_VAL_R2})"
            ),
            move || gbdt_model(&ds.feature_names, false, trees_for_fallback),
        );
    }

    Ok(gbdt_model(&ds.feature_names, false, trees))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nan_scores_do_not_panic_the_sorts() {
        let auc_nan = auc(&[f32::NAN, 0.5, 0.2, 0.9], &[1, 0, 1, 0]);
        assert!((0.0..=1.0).contains(&auc_nan));
        let edges = compute_bin_uppers(&[0.3, f32::NAN, 0.1, f32::INFINITY, 0.2], 4);
        assert_eq!(edges, vec![0.1, 0.2, 0.3]);
    }

    #[test]
    fn auc_perfect_reversed_tied_and_empty() {
        let scores = [0.1_f32, 0.2, 0.8, 0.9];
        let labels = [0u8, 0, 1, 1];
        assert!((auc(&scores, &labels) - 1.0).abs() < 1e-9);
        let labels_rev = [1u8, 1, 0, 0];
        assert!((auc(&scores, &labels_rev) - 0.0).abs() < 1e-9);
        let tied = [0.5_f32, 0.5, 0.5, 0.5];
        assert!((auc(&tied, &labels) - 0.5).abs() < 1e-9);
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
            let label = if x[r * nf] + 0.2 * lcg(&mut s) > 0.6 {
                1u8
            } else {
                0
            };
            y[r] = label;
            groups[r] = (r as u32) % 50;
        }
        Dataset {
            x,
            y,
            groups,
            n_features: nf,
            feature_names: scoring_crate::gbdt_eval::indexed_feature_names(nf as u32),
        }
    }

    #[test]
    fn separable_data_trains_and_separates() {
        let ds = separable(6, 4000, 0x1234);
        let model =
            train_gbdt(&ds, &TrainParams::default(), 0xC0FFEE).expect("separable passes gate");
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
        let a = train_gbdt(&ds, &TrainParams::default(), 99).expect("separable passes gate");
        let b = train_gbdt(&ds, &TrainParams::default(), 99).expect("separable passes gate");
        assert_eq!(
            a.to_bytes(),
            b.to_bytes(),
            "same seed must yield identical model"
        );
    }

    #[test]
    fn no_positives_fails_quality_gate() {
        let ds = Dataset {
            x: vec![0.1f32; 10 * 3],
            y: vec![0u8; 10],
            groups: (0..10).collect(),
            n_features: 3,
            feature_names: scoring_crate::gbdt_eval::indexed_feature_names(3),
        };
        let m = train_gbdt(&ds, &TrainParams::default(), 1);
        assert!(matches!(m, Err(TrainError::QualityGate(_))), "got {m:?}");
    }

    #[test]
    fn too_few_rows_fails_quality_gate() {
        let n = MIN_TRAIN_ROWS - 1;
        let mut x = vec![0.0f32; n * 3];
        let mut y = vec![0u8; n];
        for r in 0..n {
            x[r * 3] = (r % 2) as f32;
            y[r] = (r % 2) as u8;
        }
        let ds = Dataset {
            x,
            y,
            groups: (0..n as u32).collect(),
            n_features: 3,
            feature_names: scoring_crate::gbdt_eval::indexed_feature_names(3),
        };
        let m = train_gbdt(&ds, &TrainParams::default(), 7);
        assert!(
            matches!(m, Err(TrainError::QualityGate(_))),
            "too few rows must fail gate, got {m:?}"
        );
    }

    #[test]
    fn allow_degenerate_returns_fallback_instead_of_error() {
        let ds = Dataset {
            x: vec![0.1f32; 10 * 3],
            y: vec![0u8; 10],
            groups: (0..10).collect(),
            n_features: 3,
            feature_names: scoring_crate::gbdt_eval::indexed_feature_names(3),
        };
        let p = TrainParams {
            allow_degenerate: true,
            ..TrainParams::default()
        };
        let m = train_gbdt(&ds, &p, 1).expect("fallback returns a model");
        assert!(
            m.trees.is_empty(),
            "fallback model is the degenerate empty model"
        );
    }

    #[test]
    fn pearson_r2_perfect_and_constant() {
        let y = [1.0_f32, 2.0, 3.0, 4.0];
        let (r, r2) = pearson_r2(&y, &y);
        assert!((r - 1.0).abs() < 1e-6 && (r2 - 1.0).abs() < 1e-6);
        let (_, r2c) = pearson_r2(&[2.5_f32; 4], &y);
        assert!(r2c <= 1e-6);
    }

    #[test]
    fn regression_fits_linear_target_monotone() {
        let mut st = 0x1234u64;
        let (n, nf) = (4000usize, 2usize);
        let (mut x, mut y, mut g) = (Vec::new(), Vec::new(), Vec::new());
        for i in 0..n {
            let x0 = lcg(&mut st);
            let x1 = lcg(&mut st);
            x.push(x0);
            x.push(x1);
            y.push(2.0 * x0 - 1.0);
            g.push(i as u32);
        }
        let ds = RegressionDataset {
            x,
            y,
            groups: g,
            n_features: nf,
            feature_names: scoring_crate::gbdt_eval::indexed_feature_names(nf as u32),
        };
        let p = TrainParams::default();
        let model = train_gbdt_regression(&ds, &p, 42).expect("linear target passes gate");
        let lo: f32 = model.trees.iter().map(|t| t.eval(&[0.0, 0.5])).sum();
        let hi: f32 = model.trees.iter().map(|t| t.eval(&[1.0, 0.5])).sum();
        assert!(hi > lo + 0.5, "must be monotone in x0: lo={lo} hi={hi}");
        assert!(
            !model.apply_sigmoid,
            "regression model must not apply sigmoid"
        );
        assert!(
            model.iso_x.is_empty(),
            "regression model must have no isotonic calibration"
        );
    }
}
