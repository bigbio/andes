use rayon::prelude::*;
use scoring_crate::gbdt::Tree;

/// Nodes with fewer rows than this build their histograms on the calling thread.
const PAR_MIN_ROWS: usize = 4096;

/// Most features per histogram task: each task reads a node's gradients once for its block.
/// Fewer features share a task when the blocks would not give every thread one.
const FEATURE_BLOCK: usize = 4;

pub struct TreeParams {
    pub max_depth: u32,
    pub n_bins: usize,
    pub lambda: f32,
    pub min_hessian: f32,
    pub n_features: usize,
}

#[derive(Clone)]
struct BinHist {
    g: Vec<f32>,
    h: Vec<f32>,
}

impl BinHist {
    fn new(n_bins: usize) -> Self {
        Self {
            g: vec![0.0; n_bins],
            h: vec![0.0; n_bins],
        }
    }

    fn accumulate(&mut self, bin: u8, grad: f32, hess: f32) {
        let b = bin as usize;
        self.g[b] += grad;
        self.h[b] += hess;
    }
}

struct Builder<'a> {
    params: &'a TreeParams,
    cols: &'a [u8],
    n_rows: usize,
    grad: &'a [f32],
    hess: &'a [f32],
    bin_uppers: &'a [Vec<f32>],
    feature: Vec<i32>,
    threshold: Vec<f32>,
    left: Vec<i32>,
    right: Vec<i32>,
    value: Vec<f32>,
    default_left: Vec<u8>,
    leaf_of: Vec<u32>,
}

impl<'a> Builder<'a> {
    fn new(
        params: &'a TreeParams,
        cols: &'a [u8],
        grad: &'a [f32],
        hess: &'a [f32],
        bin_uppers: &'a [Vec<f32>],
    ) -> Self {
        let capacity = 1 << (params.max_depth.min(20) + 1);
        Self {
            params,
            cols,
            n_rows: grad.len(),
            grad,
            hess,
            bin_uppers,
            feature: Vec::with_capacity(capacity),
            threshold: Vec::with_capacity(capacity),
            left: Vec::with_capacity(capacity),
            right: Vec::with_capacity(capacity),
            value: Vec::with_capacity(capacity),
            default_left: Vec::with_capacity(capacity),
            leaf_of: vec![0; grad.len()],
        }
    }

    fn reserve_node(&mut self) -> usize {
        let idx = self.feature.len();
        self.feature.push(0);
        self.threshold.push(0.0);
        self.left.push(-1);
        self.right.push(-1);
        self.value.push(0.0);
        self.default_left.push(1);
        idx
    }

    fn write_leaf(&mut self, idx: usize, weight: f32) {
        self.feature[idx] = -1;
        self.threshold[idx] = 0.0;
        self.left[idx] = -1;
        self.right[idx] = -1;
        self.value[idx] = weight;
        self.default_left[idx] = 1;
    }

    fn mark_leaf(&mut self, rows: &[usize], idx: usize) {
        for &r in rows {
            self.leaf_of[r] = idx as u32;
        }
    }

    fn write_internal(
        &mut self,
        idx: usize,
        feat: i32,
        thresh: f32,
        left_child: usize,
        right_child: usize,
    ) {
        self.feature[idx] = feat;
        self.threshold[idx] = thresh;
        self.left[idx] = left_child as i32;
        self.right[idx] = right_child as i32;
        self.value[idx] = 0.0;
        self.default_left[idx] = 1;
    }

    fn col(&self, feat: usize) -> &'a [u8] {
        &self.cols[feat * self.n_rows..(feat + 1) * self.n_rows]
    }

    /// Gradient and hessian histograms of every feature over `rows`. Each bin adds its rows in
    /// the order of `rows`, so the sums do not depend on how the features are split over threads.
    fn histograms(&self, rows: &[usize]) -> Vec<BinHist> {
        let nf = self.params.n_features;
        let nb = self.params.n_bins;
        let per_task = (nf / rayon::current_num_threads().max(1)).clamp(1, FEATURE_BLOCK);
        let block = |first: usize| -> Vec<BinHist> {
            let feats = first..(first + per_task).min(nf);
            let cols: Vec<&[u8]> = feats.clone().map(|f| self.col(f)).collect();
            let mut hists: Vec<BinHist> = feats.map(|_| BinHist::new(nb)).collect();
            for &r in rows {
                let (g, h) = (self.grad[r], self.hess[r]);
                for (hist, col) in hists.iter_mut().zip(&cols) {
                    hist.accumulate(col[r], g, h);
                }
            }
            hists
        };
        let firsts: Vec<usize> = (0..nf).step_by(per_task).collect();
        if rows.len() < PAR_MIN_ROWS {
            firsts.into_iter().flat_map(block).collect()
        } else {
            firsts.into_par_iter().map(block).collect::<Vec<_>>().concat()
        }
    }

    /// Sums of the gradients and hessians of `rows`, added in row order.
    fn totals(&self, rows: &[usize]) -> (f32, f32) {
        let g_total: f32 = rows.iter().map(|&r| self.grad[r]).sum();
        let h_total: f32 = rows.iter().map(|&r| self.hess[r]).sum();
        (g_total, h_total)
    }

    /// Rows with bin <= `bin` of feature `feat` and the others, each in the order of `rows`.
    fn split_rows(&self, rows: &[usize], feat: usize, bin: usize) -> (Vec<usize>, Vec<usize>) {
        let col = self.col(feat);
        let side = |c: &[usize]| -> (Vec<usize>, Vec<usize>) {
            c.iter().partition(|&&r| (col[r] as usize) <= bin)
        };
        if rows.len() < PAR_MIN_ROWS {
            return side(rows);
        }
        let parts: Vec<(Vec<usize>, Vec<usize>)> =
            rows.par_chunks(PAR_MIN_ROWS).map(side).collect();
        let n_left = parts.iter().map(|p| p.0.len()).sum();
        let mut left = Vec::with_capacity(n_left);
        let mut right = Vec::with_capacity(rows.len() - n_left);
        for (l, r) in parts {
            left.extend(l);
            right.extend(r);
        }
        (left, right)
    }

    fn build(&mut self, rows: &[usize], depth: u32) -> usize {
        let nb = self.params.n_bins;
        let lambda = self.params.lambda;
        let min_h = self.params.min_hessian;

        let leaf = depth >= self.params.max_depth || rows.len() < 2;
        let ((g_total, h_total), hists) = {
            let this = &*self;
            if leaf {
                (this.totals(rows), Vec::new())
            } else {
                rayon::join(|| this.totals(rows), || this.histograms(rows))
            }
        };
        let leaf_weight = -g_total / (h_total + lambda);

        let node_idx = self.reserve_node();

        if leaf {
            self.write_leaf(node_idx, leaf_weight);
            self.mark_leaf(rows, node_idx);
            return node_idx;
        }

        let base_score = g_total * g_total / (h_total + lambda);
        let mut best_gain = 0.0_f32;
        let mut best_feat: Option<usize> = None;
        let mut best_bin: usize = 0;

        for (feat, hist) in hists.iter().enumerate() {
            let mut gl = 0.0_f32;
            let mut hl = 0.0_f32;
            for bin in 0..(nb - 1) {
                gl += hist.g[bin];
                hl += hist.h[bin];
                let gr = g_total - gl;
                let hr = h_total - hl;
                if hl < min_h || hr < min_h {
                    continue;
                }
                let gain = gl * gl / (hl + lambda) + gr * gr / (hr + lambda) - base_score;
                if gain > best_gain {
                    best_gain = gain;
                    best_feat = Some(feat);
                    best_bin = bin;
                }
            }
        }

        match best_feat {
            None => {
                self.write_leaf(node_idx, leaf_weight);
                self.mark_leaf(rows, node_idx);
            }
            Some(feat) => {
                let (left_rows, right_rows) = self.split_rows(rows, feat, best_bin);

                let left_idx = self.build(&left_rows, depth + 1);
                let right_idx = self.build(&right_rows, depth + 1);

                let threshold = self.bin_uppers[feat][best_bin];
                self.write_internal(node_idx, feat as i32, threshold, left_idx, right_idx);
            }
        }

        node_idx
    }
}

/// Row-major bins (`binned[row * n_features + feature]`) to feature-major columns
/// (`cols[feature * n_rows + row]`), the layout of [`fit_tree_cols`].
pub fn to_columns(binned: &[u8], n_features: usize) -> Vec<u8> {
    let n_rows = binned.len() / n_features.max(1);
    let mut cols = vec![0u8; binned.len()];
    for (r, row) in binned.chunks_exact(n_features.max(1)).enumerate() {
        for (f, &b) in row.iter().enumerate() {
            cols[f * n_rows + r] = b;
        }
    }
    cols
}

pub fn fit_tree(
    binned: &[u8],
    grad: &[f32],
    hess: &[f32],
    params: &TreeParams,
    bin_uppers: &[Vec<f32>],
) -> Tree {
    fit_tree_with_leaves(binned, grad, hess, params, bin_uppers).0
}

pub fn fit_tree_with_leaves(
    binned: &[u8],
    grad: &[f32],
    hess: &[f32],
    params: &TreeParams,
    bin_uppers: &[Vec<f32>],
) -> (Tree, Vec<u32>) {
    let cols = to_columns(binned, params.n_features);
    fit_tree_cols(&cols, grad, hess, params, bin_uppers)
}

/// [`fit_tree_with_leaves`] on feature-major bins (`cols[feature * n_rows + row]`).
pub fn fit_tree_cols(
    cols: &[u8],
    grad: &[f32],
    hess: &[f32],
    params: &TreeParams,
    bin_uppers: &[Vec<f32>],
) -> (Tree, Vec<u32>) {
    debug_assert_eq!(grad.len(), hess.len(), "grad/hess length mismatch");
    debug_assert_eq!(
        cols.len(),
        grad.len() * params.n_features,
        "binned length mismatch"
    );

    let n_rows = grad.len();
    let all_rows: Vec<usize> = (0..n_rows).collect();

    let mut builder = Builder::new(params, cols, grad, hess, bin_uppers);
    builder.build(&all_rows, 0);

    let tree = Tree {
        feature: builder.feature,
        threshold: builder.threshold,
        left: builder.left,
        right: builder.right,
        value: builder.value,
        default_left: builder.default_left,
    };
    (tree, builder.leaf_of)
}

#[cfg(test)]
mod tests {
    use super::*;
    use scoring_crate::gbdt::Tree;

    #[test]
    fn leaves_match_eval_on_continuous_features() {
        let (n, nf, nb) = (500usize, 3usize, 8usize);
        let uppers: Vec<Vec<f32>> = (0..nf)
            .map(|f| (0..nb).map(|b| b as f32 + 0.25 * f as f32).collect())
            .collect();
        let mut state = 0x9e37_79b9_7f4a_7c15u64;
        let mut next = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (state >> 33) as u32
        };
        let mut x = vec![0.0f32; n * nf];
        let mut binned = vec![0u8; n * nf];
        let mut grad = vec![0.0f32; n];
        for r in 0..n {
            for f in 0..nf {
                let v = match next() % 10 {
                    0 => f32::NAN,
                    1 => uppers[f][(next() as usize) % nb],
                    _ => (next() % 1000) as f32 / 100.0 - 1.0,
                };
                x[r * nf + f] = v;
                let b = uppers[f].partition_point(|&u| u < v).min(nb - 1);
                binned[r * nf + f] = b as u8;
            }
            grad[r] = (next() % 200) as f32 / 100.0 - 1.0;
        }
        let hess = vec![1.0f32; n];
        let params = TreeParams {
            max_depth: 4,
            n_bins: nb,
            lambda: 1.0,
            min_hessian: 1.0,
            n_features: nf,
        };
        let (tree, leaf_of) = fit_tree_with_leaves(&binned, &grad, &hess, &params, &uppers);
        assert!(tree.feature.len() > 3, "the tree splits");
        for r in 0..n {
            let row = &x[r * nf..(r + 1) * nf];
            let leaf = leaf_of[r] as usize;
            assert!(tree.feature[leaf] < 0, "row {r} is assigned an internal node");
            assert_eq!(tree.value[leaf].to_bits(), tree.eval(row).to_bits(), "row {r}");
        }
    }

    #[test]
    fn tree_splits_on_discriminative_feature() {
        let n = 100usize;
        let n_features = 2usize;
        let n_bins = 4usize;
        let mut binned = vec![0u8; n * n_features];
        let mut grad = vec![0.0_f32; n];
        let hess = vec![1.0_f32; n];
        for r in 0..n {
            let hi = r % 2 == 0;
            binned[r * n_features] = if hi { 3 } else { 0 };
            binned[r * n_features + 1] = (r % n_bins) as u8;
            grad[r] = if hi { -1.0 } else { 1.0 };
        }
        let params = TreeParams {
            max_depth: 2,
            n_bins,
            lambda: 1.0,
            min_hessian: 1.0,
            n_features,
        };
        let bin_uppers = vec![vec![0.5_f32, 1.5, 2.5, 3.5]; n_features];
        let tree: Tree = fit_tree(&binned, &grad, &hess, &params, &bin_uppers);
        assert_eq!(
            tree.feature[0], 0,
            "root should split on feature 0, got {:?}",
            tree.feature
        );
        assert!(tree.feature.len() >= 3, "a split produces >= 3 nodes");
    }

    #[test]
    fn emitted_tree_passes_walker_decode() {
        use scoring_crate::gbdt::GbdtPeakModel;
        let n = 60usize;
        let nf = 1usize;
        let nb = 3usize;
        let mut binned = vec![0u8; n * nf];
        let mut grad = vec![0.0_f32; n];
        let hess = vec![1.0_f32; n];
        for r in 0..n {
            let hi = r < 30;
            binned[r] = if hi { 2 } else { 0 };
            grad[r] = if hi { -1.0 } else { 1.0 };
        }
        let params = TreeParams {
            max_depth: 3,
            n_bins: nb,
            lambda: 1.0,
            min_hessian: 1.0,
            n_features: nf,
        };
        let bin_uppers = vec![vec![0.5_f32, 1.5, 2.5]];
        let tree = fit_tree(&binned, &grad, &hess, &params, &bin_uppers);
        let model = GbdtPeakModel {
            packed: Default::default(),
            n_features: 1,
            feature_names: scoring_crate::gbdt::indexed_feature_names(1),
            apply_sigmoid: false,
            trees: vec![tree],
            iso_x: vec![],
            iso_y: vec![],
        };
        let back = GbdtPeakModel::from_bytes(&model.to_bytes());
        assert!(
            back.is_ok(),
            "emitted tree must pass walker validate-at-decode: {:?}",
            back.err()
        );
    }

    #[test]
    fn single_leaf_when_no_gain() {
        let n = 20usize;
        let nf = 1usize;
        let nb = 2usize;
        let binned = vec![0u8; n * nf];
        let grad = vec![0.5_f32; n];
        let hess = vec![1.0_f32; n];
        let params = TreeParams {
            max_depth: 3,
            n_bins: nb,
            lambda: 1.0,
            min_hessian: 1.0,
            n_features: nf,
        };
        let bin_uppers = vec![vec![0.5_f32, 1.5]];
        let tree = fit_tree(&binned, &grad, &hess, &params, &bin_uppers);
        assert_eq!(tree.feature, vec![-1], "expected a single leaf");
        let g: f32 = grad.iter().sum();
        let h: f32 = hess.iter().sum();
        assert!(
            (tree.value[0] - (-g / (h + 1.0))).abs() < 1e-4,
            "leaf weight = -G/(H+λ)"
        );
    }
}
