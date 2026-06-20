//! Histogram regression tree — the weak learner for one GBDT boosting round.
//!
//! Fits one depth-limited regression tree to per-row (gradient, hessian) pairs
//! using pre-binned integer features (0..n_bins). Emits a
//! `scoring_crate::gbdt_eval::Tree` in pre-order layout so that every child
//! index is strictly greater than its parent index — the invariant required by
//! `gbdt_eval::Tree::validate`.

use scoring_crate::gbdt_eval::Tree;

/// Hyper-parameters for one tree fit.
pub struct TreeParams {
    pub max_depth: u32,
    /// Number of bins per feature (bins are 0-indexed: 0..n_bins).
    pub n_bins: usize,
    /// L2 regularisation on leaf weights (λ in the XGBoost gain formula).
    pub lambda: f32,
    /// Minimum summed hessian required in each child after a split.
    pub min_hessian: f32,
    pub n_features: usize,
}

// ---------------------------------------------------------------------------
// Internal accumulator for one node's histogram over a single feature.
// ---------------------------------------------------------------------------

/// Per-bin sums of (gradient, hessian) for a single feature.
#[derive(Clone)]
struct BinHist {
    g: Vec<f32>,
    h: Vec<f32>,
}

impl BinHist {
    fn new(n_bins: usize) -> Self {
        Self { g: vec![0.0; n_bins], h: vec![0.0; n_bins] }
    }

    fn accumulate(&mut self, bin: u8, grad: f32, hess: f32) {
        let b = bin as usize;
        self.g[b] += grad;
        self.h[b] += hess;
    }
}

// ---------------------------------------------------------------------------
// Builder state — passed through the recursion so we can append nodes and
// back-patch child indices.
// ---------------------------------------------------------------------------

struct Builder<'a> {
    params: &'a TreeParams,
    binned: &'a [u8],  // row-major [row * n_features + feat]
    grad: &'a [f32],
    hess: &'a [f32],
    bin_uppers: &'a [Vec<f32>],
    // Growing SoA vecs — will be moved into the final Tree.
    feature: Vec<i32>,
    threshold: Vec<f32>,
    left: Vec<i32>,
    right: Vec<i32>,
    value: Vec<f32>,
    default_left: Vec<u8>,
}

impl<'a> Builder<'a> {
    fn new(
        params: &'a TreeParams,
        binned: &'a [u8],
        grad: &'a [f32],
        hess: &'a [f32],
        bin_uppers: &'a [Vec<f32>],
    ) -> Self {
        let capacity = 1 << (params.max_depth.min(20) + 1); // rough upper bound
        Self {
            params,
            binned,
            grad,
            hess,
            bin_uppers,
            feature: Vec::with_capacity(capacity),
            threshold: Vec::with_capacity(capacity),
            left: Vec::with_capacity(capacity),
            right: Vec::with_capacity(capacity),
            value: Vec::with_capacity(capacity),
            default_left: Vec::with_capacity(capacity),
        }
    }

    /// Append a placeholder node and return its index.
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

    /// Write a leaf into the slot at `idx`.
    fn write_leaf(&mut self, idx: usize, weight: f32) {
        self.feature[idx] = -1;
        self.threshold[idx] = 0.0;
        self.left[idx] = -1;
        self.right[idx] = -1;
        self.value[idx] = weight;
        self.default_left[idx] = 1;
    }

    /// Write an internal split node into the slot at `idx`.
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

    // -----------------------------------------------------------------------
    // Build a node recursively (pre-order: current node first, then children).
    // Returns the index of the node appended.
    // -----------------------------------------------------------------------
    fn build(&mut self, rows: &[usize], depth: u32) -> usize {
        let nf = self.params.n_features;
        let nb = self.params.n_bins;
        let lambda = self.params.lambda;
        let min_h = self.params.min_hessian;

        // Compute totals for this node.
        let g_total: f32 = rows.iter().map(|&r| self.grad[r]).sum();
        let h_total: f32 = rows.iter().map(|&r| self.hess[r]).sum();
        let leaf_weight = -g_total / (h_total + lambda);

        // Reserve the slot now so this node's index is lower than its children.
        let node_idx = self.reserve_node();

        // Check stopping conditions.
        if depth >= self.params.max_depth || rows.len() < 2 {
            self.write_leaf(node_idx, leaf_weight);
            return node_idx;
        }

        // Build per-feature histograms.
        let mut hists: Vec<BinHist> = (0..nf).map(|_| BinHist::new(nb)).collect();
        for &r in rows {
            let g = self.grad[r];
            let h = self.hess[r];
            for (feat, hist) in hists.iter_mut().enumerate() {
                let bin = self.binned[r * nf + feat];
                hist.accumulate(bin, g, h);
            }
        }

        // Find best split using prefix scan (gain = GL²/(HL+λ) + GR²/(HR+λ) − G²/(H+λ)).
        let base_score = g_total * g_total / (h_total + lambda);
        let mut best_gain = 0.0_f32;
        let mut best_feat: Option<usize> = None;
        let mut best_bin: usize = 0;

        for (feat, hist) in hists.iter().enumerate() {
            let mut gl = 0.0_f32;
            let mut hl = 0.0_f32;
            // Scan bins 0..nb-1 (split after bin b means left = bins 0..=b, right = b+1..nb-1).
            for bin in 0..(nb - 1) {
                gl += hist.g[bin];
                hl += hist.h[bin];
                let gr = g_total - gl;
                let hr = h_total - hl;
                if hl < min_h || hr < min_h {
                    continue;
                }
                let gain =
                    gl * gl / (hl + lambda) + gr * gr / (hr + lambda) - base_score;
                // Deterministic tie-break: prefer lower (feat, bin) — since we
                // iterate feat then bin in ascending order, strictly-greater
                // comparison keeps the first best found.
                if gain > best_gain {
                    best_gain = gain;
                    best_feat = Some(feat);
                    best_bin = bin;
                }
            }
        }

        match best_feat {
            None => {
                // No positive-gain split found — emit a leaf.
                self.write_leaf(node_idx, leaf_weight);
            }
            Some(feat) => {
                // Partition rows: left ≤ best_bin, right > best_bin.
                let (left_rows, right_rows): (Vec<usize>, Vec<usize>) = rows
                    .iter()
                    .partition(|&&r| (self.binned[r * nf + feat] as usize) <= best_bin);

                // Recurse — children are appended AFTER node_idx (pre-order).
                let left_idx = self.build(&left_rows, depth + 1);
                let right_idx = self.build(&right_rows, depth + 1);

                let threshold = self.bin_uppers[feat][best_bin];
                self.write_internal(node_idx, feat as i32, threshold, left_idx, right_idx);
            }
        }

        node_idx
    }
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Fit one regression tree to the given per-row gradients/hessians using
/// pre-binned features.
///
/// `binned` is row-major: `binned[row * n_features + feat]` holds the bin
/// index (0..n_bins) for that (row, feature) pair.
///
/// `bin_uppers[feat][bin]` is the continuous upper edge of `bin` for `feat`;
/// the emitted tree stores it as the split threshold so `x <= threshold`
/// reproduces the histogram split when the walker sees continuous features.
///
/// The returned `Tree` satisfies `gbdt_eval::Tree::validate` (pre-order,
/// strictly-increasing child indices, all SoA vecs share the same length).
pub fn fit_tree(
    binned: &[u8],
    grad: &[f32],
    hess: &[f32],
    params: &TreeParams,
    bin_uppers: &[Vec<f32>],
) -> Tree {
    debug_assert_eq!(grad.len(), hess.len(), "grad/hess length mismatch");
    debug_assert_eq!(
        binned.len(),
        grad.len() * params.n_features,
        "binned length mismatch"
    );

    let n_rows = grad.len();
    let all_rows: Vec<usize> = (0..n_rows).collect();

    let mut builder = Builder::new(params, binned, grad, hess, bin_uppers);
    builder.build(&all_rows, 0);

    Tree {
        feature: builder.feature,
        threshold: builder.threshold,
        left: builder.left,
        right: builder.right,
        value: builder.value,
        default_left: builder.default_left,
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use scoring_crate::gbdt_eval::Tree;

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
            binned[r * n_features] = if hi { 3 } else { 0 }; // feature 0 separates
            binned[r * n_features + 1] = (r % n_bins) as u8; // feature 1 noise
            grad[r] = if hi { -1.0 } else { 1.0 }; // opposite gradients
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
            tree.feature[0],
            0,
            "root should split on feature 0, got {:?}",
            tree.feature
        );
        assert!(tree.feature.len() >= 3, "a split produces >= 3 nodes");
    }

    #[test]
    fn emitted_tree_passes_walker_decode() {
        // Round-trip the tree through the AGBD blob to assert it satisfies
        // gbdt_eval's structural validation (pre-order, in-range children).
        use scoring_crate::gbdt_eval::GbdtPeakModel;
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
            n_features: 1,
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
        // All gradients equal → no split improves → a single leaf with weight −G/(H+λ).
        let n = 20usize;
        let nf = 1usize;
        let nb = 2usize;
        let binned = vec![0u8; n * nf]; // all same bin → no split possible anyway
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
