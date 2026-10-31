//! Hand-rolled struct-of-arrays GBDT evaluator (zero native deps).
//!
//! Decodes the `AGBD` v2 blob (produced offline by `training/gbdt/transcode.py`
//! from a LightGBM binary classifier + scikit-learn IsotonicRegression) and
//! evaluates it on a peptide-AGNOSTIC per-peak feature vector
//! (`crate::peak_features`). Output is `log(s/(1-s))`, the additive LLR term the
//! rank scorer folds in (`scored_spectrum`).

use std::io::{Cursor, Read};

use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};

use thiserror::Error;

const MAGIC: &[u8; 4] = b"AGBD";
const FORMAT_VERSION: u32 = 2;
const FLAG_SIGMOID: u32 = 1;
/// Upper bounds on the feature-name trailer, checked before allocating.
const MAX_FEATURES: u32 = 100_000;
const MAX_NAME_LEN: u32 = 1024;
/// Clamp for the calibrated probability before taking the logit, so the LLR
/// term stays finite (a degenerate 0 or 1 would give ±inf).
const PROB_EPS: f32 = 1e-6;

#[derive(Debug, Error)]
pub enum GbdtError {
    #[error("bad magic (not an AGBD blob)")]
    BadMagic,
    #[error("unsupported AGBD version {0}")]
    BadVersion(u32),
    #[error("truncated AGBD blob: {0}")]
    Io(#[from] std::io::Error),
    #[error("malformed AGBD blob: {0}")]
    Malformed(String),
    #[error("GBDT feature-contract violation for {kind} model: blob declares {got} features but the evaluator extracts {expected}")]
    FeatureCountMismatch {
        kind: &'static str,
        expected: u32,
        got: u32,
    },
    #[error("GBDT feature-contract violation for {kind} model: feature {position} is named '{got}' in the blob but the evaluator extracts '{expected}'")]
    FeatureNameMismatch {
        kind: &'static str,
        position: usize,
        expected: String,
        got: String,
    },
}

/// Reject implausibly large counts from a corrupt/version-skewed blob before
/// allocating, so decoding fails cleanly instead of OOM-aborting. Generous
/// bounds — a real per-peak GBDT is far smaller.
fn checked_count(n: u32, what: &str, max: u32) -> Result<usize, GbdtError> {
    if n > max {
        return Err(GbdtError::Malformed(format!(
            "{what} count {n} exceeds max {max}"
        )));
    }
    Ok(n as usize)
}

/// One regression tree in struct-of-arrays layout. All vecs have length
/// `n_nodes`; node 0 is the root.
#[derive(Debug, Clone, PartialEq)]
pub struct Tree {
    pub feature: Vec<i32>, // -1 => leaf
    pub threshold: Vec<f32>,
    pub left: Vec<i32>,
    pub right: Vec<i32>,
    pub value: Vec<f32>,       // leaf output
    pub default_left: Vec<u8>, // 1 => NaN feature descends left
}

impl Tree {
    /// Sum of the leaf value reached for `x`.
    pub fn eval(&self, x: &[f32]) -> f32 {
        // Assumes a structurally-validated tree (see `Tree::validate`, enforced
        // by `from_bytes`): feature/child arrays are consistent and internal
        // children are valid in-range indices, so direct indexing cannot panic
        // and traversal terminates. Models built in-process by the trainer
        // satisfy this by construction.
        let mut node = 0usize;
        loop {
            let feat = self.feature[node];
            if feat < 0 {
                return self.value[node];
            }
            let v = x.get(feat as usize).copied().unwrap_or(f32::NAN);
            let go_left = if v.is_nan() {
                self.default_left[node] == 1
            } else {
                v <= self.threshold[node]
            };
            node = if go_left {
                self.left[node]
            } else {
                self.right[node]
            } as usize;
        }
    }

    /// Validate structural invariants so `eval` can index children without
    /// bounds checks: all SoA vecs share one length n>=1; every node is either
    /// a leaf (feature < 0) or an internal node whose left/right children are
    /// valid indices in [0, n). (Leaf child links are not followed by `eval`,
    /// so they are unconstrained.)
    fn validate(&self) -> Result<(), GbdtError> {
        let n = self.feature.len();
        if n == 0 {
            return Err(GbdtError::Malformed("tree has zero nodes".into()));
        }
        if self.threshold.len() != n
            || self.left.len() != n
            || self.right.len() != n
            || self.value.len() != n
            || self.default_left.len() != n
        {
            return Err(GbdtError::Malformed(
                "tree SoA arrays have mismatched lengths".into(),
            ));
        }
        #[allow(clippy::needless_range_loop)] // indexing multiple parallel arrays by node index
        for node in 0..n {
            if self.feature[node] >= 0 {
                // internal node: children must be valid indices and strictly
                // greater (guarantees acyclicity for pre-order/DFS layout)
                for child in [self.left[node], self.right[node]] {
                    if child < 0 || (child as usize) >= n {
                        return Err(GbdtError::Malformed(format!(
                            "internal node {node} has out-of-range child {child} (n={n})"
                        )));
                    }
                    if (child as usize) <= node {
                        return Err(GbdtError::Malformed(format!(
                            "internal node {node} child {child} is not strictly greater (non-preorder/cyclic)"
                        )));
                    }
                }
            }
        }
        Ok(())
    }
}

#[derive(Debug, Clone)]
pub struct GbdtPeakModel {
    pub n_features: u32,
    /// Ordered feature names, one per feature (`feature_names.len() == n_features`);
    /// index `i` is the feature the tree splits call `i`. Written in the blob's
    /// trailer and checked against the extractor by [`Self::validate_features`].
    pub feature_names: Vec<String>,
    pub apply_sigmoid: bool,
    pub trees: Vec<Tree>,
    pub iso_x: Vec<f32>,
    pub iso_y: Vec<f32>,
    /// `trees` repacked for evaluation, built on first use. Derived from `trees`;
    /// construct with `Default::default()`.
    pub packed: std::sync::OnceLock<PackedTrees>,
}

/// Placeholder names `f0..f{n-1}` for a model that is not bound to an extractor.
pub fn indexed_feature_names(n: u32) -> Vec<String> {
    (0..n).map(|i| format!("f{i}")).collect()
}

/// Equality is over the model itself, not whether its evaluation cache is built.
impl PartialEq for GbdtPeakModel {
    fn eq(&self, other: &Self) -> bool {
        self.n_features == other.n_features
            && self.feature_names == other.feature_names
            && self.apply_sigmoid == other.apply_sigmoid
            && self.trees == other.trees
            && self.iso_x == other.iso_x
            && self.iso_y == other.iso_y
    }
}

/// High bit of [`PackedNode::right`]: a NaN feature value descends left.
const DEFAULT_LEFT_BIT: u32 = 1 << 31;

/// One tree node in 16 bytes, so a node visit touches one cache line instead of
/// the five arrays of [`Tree`]. `value` is the threshold of an internal node and
/// the output of a leaf (`feature < 0`).
#[derive(Debug, Clone, Copy)]
struct PackedNode {
    feature: i32,
    value: f32,
    left: u32,
    right: u32,
}

/// Every tree of a model as packed nodes; see [`GbdtPeakModel::packed`].
#[derive(Debug, Clone, Default)]
pub struct PackedTrees {
    trees: Vec<Box<[PackedNode]>>,
    /// The same trees in [`PerfectTree`] layout for the batched walk; `None` for a
    /// tree deeper than [`PERFECT_MAX_DEPTH`], which keeps the node walk.
    perfect: Vec<Option<PerfectTree>>,
    /// One more than the largest feature index any tree splits on (1 when none).
    n_features: usize,
}

impl PackedTrees {
    fn from_trees(trees: &[Tree]) -> Self {
        let packed: Vec<Box<[PackedNode]>> = trees
            .iter()
            .map(|t| {
                (0..t.feature.len())
                    .map(|i| {
                        let leaf = t.feature[i] < 0;
                        PackedNode {
                            feature: t.feature[i],
                            value: if leaf { t.value[i] } else { t.threshold[i] },
                            left: if leaf { 0 } else { t.left[i] as u32 },
                            right: if leaf {
                                0
                            } else {
                                t.right[i] as u32
                                    | if t.default_left[i] == 1 {
                                        DEFAULT_LEFT_BIT
                                    } else {
                                        0
                                    }
                            },
                        }
                    })
                    .collect()
            })
            .collect();
        let perfect = packed.iter().map(|t| PerfectTree::from_packed(t)).collect();
        let n_features = packed
            .iter()
            .flat_map(|t| t.iter())
            .filter(|n| n.feature >= 0)
            .map(|n| n.feature as usize + 1)
            .max()
            .unwrap_or(1);
        PackedTrees {
            trees: packed,
            perfect,
            n_features,
        }
    }
}

/// Deepest tree given the [`PerfectTree`] layout (4096 leaves); deeper trees keep
/// the node walk.
const PERFECT_MAX_DEPTH: u32 = 12;

/// Rows walked through a [`PerfectTree`] together.
const PERFECT_LANES: usize = 4;

/// One internal node of a [`PerfectTree`]: the split threshold and the feature
/// index, whose high bit ([`DEFAULT_LEFT_BIT`]) sends a NaN value left.
#[derive(Debug, Clone, Copy)]
#[repr(C)]
struct PerfectNode {
    threshold: f32,
    feature: u32,
}

impl PerfectNode {
    /// The comparison [`eval_packed`] makes at this node.
    #[inline(always)]
    fn goes_left(self, x: &[f32]) -> bool {
        let v = x
            .get((self.feature & !DEFAULT_LEFT_BIT) as usize)
            .copied()
            .unwrap_or(f32::NAN);
        (v <= self.threshold) | (v.is_nan() & (self.feature & DEFAULT_LEFT_BIT != 0))
    }
}

/// A tree re-laid out as a complete binary tree of its maximum depth, internal
/// nodes in breadth-first order (children of `k` at `2k+1` and `2k+2`), so a walk
/// is `depth` index updates with no data-dependent branch. A leaf above the
/// bottom level is pushed down: every node below it returns its value whichever
/// way it compares, so each row still lands on the value [`eval_packed`] returns.
#[derive(Debug, Clone)]
struct PerfectTree {
    depth: u32,
    nodes: Box<[PerfectNode]>,
    leaves: Box<[f32]>,
}

impl PerfectTree {
    fn from_packed(nodes: &[PackedNode]) -> Option<Self> {
        fn depth(nodes: &[PackedNode], i: usize) -> u32 {
            let n = &nodes[i];
            if n.feature < 0 {
                0
            } else {
                1 + depth(nodes, n.left as usize)
                    .max(depth(nodes, (n.right & !DEFAULT_LEFT_BIT) as usize))
            }
        }
        // Fill the subtree at `slot` (level `level`) from packed node `node`; a
        // leaf reached above the bottom level fills its whole subtree.
        fn fill(t: &mut PerfectTree, nodes: &[PackedNode], slot: usize, node: usize, level: u32) {
            let n = nodes[node];
            if level == t.depth {
                t.leaves[slot - t.nodes.len()] = n.value;
                return;
            }
            if n.feature < 0 {
                t.nodes[slot] = PerfectNode {
                    threshold: 0.0,
                    feature: 0,
                };
                fill(t, nodes, 2 * slot + 1, node, level + 1);
                fill(t, nodes, 2 * slot + 2, node, level + 1);
            } else {
                t.nodes[slot] = PerfectNode {
                    threshold: n.value,
                    feature: n.feature as u32 | (n.right & DEFAULT_LEFT_BIT),
                };
                fill(t, nodes, 2 * slot + 1, n.left as usize, level + 1);
                let right = (n.right & !DEFAULT_LEFT_BIT) as usize;
                fill(t, nodes, 2 * slot + 2, right, level + 1);
            }
        }
        let d = depth(nodes, 0);
        if d > PERFECT_MAX_DEPTH {
            return None;
        }
        let mut t = PerfectTree {
            depth: d,
            nodes: vec![
                PerfectNode {
                    threshold: 0.0,
                    feature: 0,
                };
                (1usize << d) - 1
            ]
            .into_boxed_slice(),
            leaves: vec![0.0; 1usize << d].into_boxed_slice(),
        };
        fill(&mut t, nodes, 0, 0, 0);
        Some(t)
    }

    /// Add this tree's leaf value to `out[r]` for every row `r`.
    ///
    /// [`PERFECT_LANES`] rows walk the tree together, one level per step, so their
    /// independent loads and comparisons overlap. The last group repeats the last
    /// row to fill its lanes and keeps only the real rows' results.
    fn add_to(&self, rows: &[&[f32]], out: &mut [f32]) {
        let n = rows.len();
        let first_leaf = self.nodes.len() as u32;
        for start in (0..n).step_by(PERFECT_LANES) {
            let lane_rows: [&[f32]; PERFECT_LANES] =
                std::array::from_fn(|l| rows[(start + l).min(n - 1)]);
            let mut idx = [0u32; PERFECT_LANES];
            for _ in 0..self.depth {
                for l in 0..PERFECT_LANES {
                    let go_left = self.nodes[idx[l] as usize].goes_left(lane_rows[l]);
                    idx[l] = 2 * idx[l] + 2 - u32::from(go_left);
                }
            }
            for (o, k) in out[start..].iter_mut().zip(idx) {
                *o += self.leaves[(k - first_leaf) as usize];
            }
        }
    }
}

/// The batched walk of [`PerfectTree`]s on AVX2: eight rows per vector, the
/// node and feature loads as gathers. The rows of a block are copied into a
/// feature-major table first, so a feature value of a row is one gather of
/// `feature * rows + row`. Every row still adds its leaf value to its sum in
/// tree order, so the sums are bit-identical to the scalar walk.
#[cfg(target_arch = "x86_64")]
mod simd {
    use super::{eval_packed, PackedTrees, PerfectTree, DEFAULT_LEFT_BIT};
    use std::arch::x86_64::*;

    /// Rows per block: the block's feature table (about 17 KB for 17
    /// features) and one tree stay in the L1 cache while the block's rows
    /// walk the tree.
    const BLOCK: usize = 256;
    /// Row vectors walked together, so their gathers overlap.
    const GROUPS: usize = 4;
    /// Rows per step.
    const STEP: usize = 8 * GROUPS;

    /// Adds every tree's value for every row to `out` (zeros on entry), in
    /// tree order.
    ///
    /// # Safety
    /// The CPU must support AVX2.
    #[target_feature(enable = "avx2")]
    pub(super) unsafe fn add_all(packed: &PackedTrees, rows: &[&[f32]], out: &mut [f32]) {
        let nf = packed.n_features;
        let mut table: Vec<f32> = Vec::new();
        let mut acc = [0.0f32; BLOCK];
        for (rows, out) in rows.chunks(BLOCK).zip(out.chunks_mut(BLOCK)) {
            let n = rows.len();
            let n_pad = n.div_ceil(STEP) * STEP;
            table.clear();
            table.resize(nf * n_pad, f32::NAN);
            for f in 0..nf {
                let col = &mut table[f * n_pad..f * n_pad + n_pad];
                for (r, c) in col.iter_mut().enumerate() {
                    *c = rows[r.min(n - 1)].get(f).copied().unwrap_or(f32::NAN);
                }
            }
            acc[..n_pad].fill(0.0);
            for (t, perfect) in packed.trees.iter().zip(&packed.perfect) {
                match perfect {
                    Some(p) => {
                        for row0 in (0..n_pad).step_by(STEP) {
                            walk(p, table.as_ptr(), n_pad as i32, row0, acc.as_mut_ptr());
                        }
                    }
                    None => {
                        for (a, x) in acc.iter_mut().zip(rows) {
                            *a += eval_packed(t, x);
                        }
                    }
                }
            }
            out.copy_from_slice(&acc[..n]);
        }
    }

    /// Walks rows `row0..row0 + STEP` of the feature table `xs` (`n_pad` rows
    /// per feature) down `p` and adds the leaf values to `acc[row0..]`.
    #[inline]
    #[target_feature(enable = "avx2")]
    unsafe fn walk(p: &PerfectTree, xs: *const f32, n_pad: i32, row0: usize, acc: *mut f32) {
        let nodes = p.nodes.as_ptr() as *const i32;
        let lane = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
        let two = _mm256_set1_epi32(2);
        let all_but_sign = _mm256_set1_epi32(!(DEFAULT_LEFT_BIT as i32));
        let stride = _mm256_set1_epi32(n_pad);
        let mut idx = [_mm256_setzero_si256(); GROUPS];
        let mut row = [_mm256_setzero_si256(); GROUPS];
        for (g, r) in row.iter_mut().enumerate() {
            *r = _mm256_add_epi32(_mm256_set1_epi32((row0 + 8 * g) as i32), lane);
        }
        for _ in 0..p.depth {
            for g in 0..GROUPS {
                // A node is `{ threshold: f32, feature: u32 }`: 8 bytes.
                let threshold = _mm256_i32gather_ps::<8>(nodes as *const f32, idx[g]);
                let feature = _mm256_i32gather_epi32::<8>(nodes.add(1), idx[g]);
                let at = _mm256_add_epi32(
                    _mm256_mullo_epi32(_mm256_and_si256(feature, all_but_sign), stride),
                    row[g],
                );
                let v = _mm256_i32gather_ps::<4>(xs, at);
                let le = _mm256_cmp_ps::<_CMP_LE_OQ>(v, threshold);
                let nan = _mm256_cmp_ps::<_CMP_UNORD_Q>(v, v);
                let nan_left =
                    _mm256_and_ps(nan, _mm256_castsi256_ps(_mm256_srai_epi32::<31>(feature)));
                let left = _mm256_castps_si256(_mm256_or_ps(le, nan_left));
                // 2 * idx + 2 - 1 for a left child (the mask is -1), 2 * idx + 2 otherwise.
                idx[g] = _mm256_add_epi32(
                    _mm256_add_epi32(idx[g], idx[g]),
                    _mm256_add_epi32(two, left),
                );
            }
        }
        let first_leaf = _mm256_set1_epi32(p.nodes.len() as i32);
        for (g, &i) in idx.iter().enumerate() {
            let leaf = _mm256_i32gather_ps::<4>(p.leaves.as_ptr(), _mm256_sub_epi32(i, first_leaf));
            let a = acc.add(row0 + 8 * g);
            _mm256_storeu_ps(a, _mm256_add_ps(_mm256_loadu_ps(a), leaf));
        }
    }
}

/// [`Tree::eval`] on packed nodes: the same comparisons in the same order, so the
/// result is bit-identical.
#[inline]
fn eval_packed(nodes: &[PackedNode], x: &[f32]) -> f32 {
    let mut i = 0usize;
    loop {
        let n = &nodes[i];
        if n.feature < 0 {
            return n.value;
        }
        let v = x.get(n.feature as usize).copied().unwrap_or(f32::NAN);
        let go_left = if v.is_nan() {
            n.right & DEFAULT_LEFT_BIT != 0
        } else {
            v <= n.value
        };
        i = if go_left {
            n.left
        } else {
            n.right & !DEFAULT_LEFT_BIT
        } as usize;
    }
}

impl GbdtPeakModel {
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, GbdtError> {
        let mut c = Cursor::new(bytes);
        let mut magic = [0u8; 4];
        c.read_exact(&mut magic)?;
        if &magic != MAGIC {
            return Err(GbdtError::BadMagic);
        }
        let version = c.read_u32::<LittleEndian>()?;
        if version != FORMAT_VERSION {
            return Err(GbdtError::BadVersion(version));
        }
        let n_features = c.read_u32::<LittleEndian>()?;
        let flags = c.read_u32::<LittleEndian>()?;
        let apply_sigmoid = flags & FLAG_SIGMOID != 0;
        let n_trees = checked_count(c.read_u32::<LittleEndian>()?, "n_trees", 1_000_000)?;

        let mut trees = Vec::with_capacity(n_trees);
        for _ in 0..n_trees {
            let n = checked_count(c.read_u32::<LittleEndian>()?, "n_nodes", 10_000_000)?;
            let feature = read_i32_vec(&mut c, n)?;
            let threshold = read_f32_vec(&mut c, n)?;
            let left = read_i32_vec(&mut c, n)?;
            let right = read_i32_vec(&mut c, n)?;
            let value = read_f32_vec(&mut c, n)?;
            let mut default_left = vec![0u8; n];
            c.read_exact(&mut default_left)?;
            let tree = Tree {
                feature,
                threshold,
                left,
                right,
                value,
                default_left,
            };
            tree.validate()?;
            trees.push(tree);
        }
        let n_iso = checked_count(c.read_u32::<LittleEndian>()?, "n_iso", 10_000_000)?;
        let iso_x = read_f32_vec(&mut c, n_iso)?;
        let iso_y = read_f32_vec(&mut c, n_iso)?;
        let n_names = c.read_u32::<LittleEndian>()?;
        if n_names != n_features {
            return Err(GbdtError::Malformed(format!(
                "{n_names} feature names for {n_features} features"
            )));
        }
        let n_names = checked_count(n_names, "feature name", MAX_FEATURES)?;
        let mut feature_names = Vec::with_capacity(n_names);
        for _ in 0..n_names {
            let len = checked_count(
                c.read_u32::<LittleEndian>()?,
                "feature name length",
                MAX_NAME_LEN,
            )?;
            let mut raw = vec![0u8; len];
            c.read_exact(&mut raw)?;
            feature_names.push(String::from_utf8(raw).map_err(|_| {
                GbdtError::Malformed("feature name is not valid UTF-8".to_string())
            })?);
        }
        if c.position() != bytes.len() as u64 {
            return Err(GbdtError::Malformed(format!(
                "{} trailing bytes after the feature names",
                bytes.len() as u64 - c.position()
            )));
        }
        Ok(Self {
            packed: Default::default(),
            n_features,
            feature_names,
            apply_sigmoid,
            trees,
            iso_x,
            iso_y,
        })
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        let mut b = Vec::new();
        b.extend_from_slice(MAGIC);
        b.write_u32::<LittleEndian>(FORMAT_VERSION).unwrap();
        b.write_u32::<LittleEndian>(self.n_features).unwrap();
        b.write_u32::<LittleEndian>(if self.apply_sigmoid { FLAG_SIGMOID } else { 0 })
            .unwrap();
        b.write_u32::<LittleEndian>(self.trees.len() as u32)
            .unwrap();
        for t in &self.trees {
            let n = t.feature.len();
            b.write_u32::<LittleEndian>(n as u32).unwrap();
            write_i32_vec(&mut b, &t.feature);
            write_f32_vec(&mut b, &t.threshold);
            write_i32_vec(&mut b, &t.left);
            write_i32_vec(&mut b, &t.right);
            write_f32_vec(&mut b, &t.value);
            b.extend_from_slice(&t.default_left);
        }
        b.write_u32::<LittleEndian>(self.iso_x.len() as u32)
            .unwrap();
        write_f32_vec(&mut b, &self.iso_x);
        write_f32_vec(&mut b, &self.iso_y);
        b.write_u32::<LittleEndian>(self.feature_names.len() as u32)
            .unwrap();
        for name in &self.feature_names {
            b.write_u32::<LittleEndian>(name.len() as u32).unwrap();
            b.extend_from_slice(name.as_bytes());
        }
        b
    }

    /// Enforce the feature contract: the blob's ordered feature names must equal
    /// the names of the in-process extractor for this model role, position by
    /// position. A skewed or reordered blob otherwise scores plausibly: `eval`
    /// reads a missing feature slot as `NaN` and a permuted feature as a
    /// different quantity, yielding a finite but meaningless score with no
    /// error. Call this once at load (see `model-train`'s store reader).
    ///
    /// A count difference is reported as [`GbdtError::FeatureCountMismatch`],
    /// otherwise the first differing position as
    /// [`GbdtError::FeatureNameMismatch`]. `kind` is a short role label used
    /// only in the error message (`"peak"`, `"frag-intensity"`, `"rich-ion"`).
    pub fn validate_features(
        &self,
        expected: &[&str],
        kind: &'static str,
    ) -> Result<(), GbdtError> {
        if self.n_features as usize != expected.len() || self.feature_names.len() != expected.len()
        {
            return Err(GbdtError::FeatureCountMismatch {
                kind,
                expected: expected.len() as u32,
                got: self.n_features,
            });
        }
        for (position, (want, got)) in expected.iter().zip(&self.feature_names).enumerate() {
            if want != got {
                return Err(GbdtError::FeatureNameMismatch {
                    kind,
                    position,
                    expected: (*want).to_string(),
                    got: got.clone(),
                });
            }
        }
        Ok(())
    }

    /// Raw regression output: sum of all tree leaf values, with NO sigmoid
    /// activation and NO isotonic calibration applied.
    ///
    /// Use this when the model was trained as a pure regressor (e.g. for
    /// fragment intensity prediction) and the caller wants the uncalibrated
    /// score directly. For the classifier path (signal/noise scoring) use
    /// [`predict_logit`] or [`predict_proba`] instead.
    fn packed(&self) -> &PackedTrees {
        self.packed
            .get_or_init(|| PackedTrees::from_trees(&self.trees))
    }

    pub fn predict_value(&self, x: &[f32]) -> f32 {
        self.packed().trees.iter().map(|t| eval_packed(t, x)).sum()
    }

    /// Batched `predict_value`: trees OUTER, rows INNER.
    ///
    /// The per-row loop re-streams the entire ensemble for every row, so a 300-tree
    /// ~790 KB model is pulled through L2/L3 once per fragment ion (~52 per PSM).
    /// Walking trees on the outside instead pulls each ~2.6 KB tree into L1 once and
    /// reuses it across every row, cutting cache traffic by roughly the row count.
    /// The K-sweep evidence says this walk is memory-residency-bound, not ALU-bound:
    /// dropping 300 -> 100 trees is a 3x work reduction but bought -41% wall.
    ///
    /// BIT-IDENTICAL to calling `predict_value` per row: `sum()` over f32 folds left
    /// from 0.0 in tree order, and accumulating `out[r] += tree.eval(row_r)` visits
    /// the trees in that same order, so every partial sum matches exactly.
    pub fn predict_value_batch(&self, rows: &[&[f32]], out: &mut [f32]) {
        assert_eq!(
            rows.len(),
            out.len(),
            "predict_value_batch: rows/out length mismatch"
        );
        out.fill(0.0);
        let packed = self.packed();
        #[cfg(target_arch = "x86_64")]
        if !rows.is_empty() && std::arch::is_x86_feature_detected!("avx2") {
            // SAFETY: the CPU supports AVX2, as just detected.
            unsafe { simd::add_all(packed, rows, out) };
            return;
        }
        for (t, perfect) in packed.trees.iter().zip(&packed.perfect) {
            match perfect {
                Some(p) => p.add_to(rows, out),
                None => {
                    for (o, x) in out.iter_mut().zip(rows.iter()) {
                        *o += eval_packed(t, x);
                    }
                }
            }
        }
    }

    /// Batched [`predict_logit`]; see [`Self::predict_value_batch`] for why this is
    /// faster and why it is bit-identical.
    pub fn predict_logit_batch(&self, rows: &[&[f32]], out: &mut [f32]) {
        self.predict_value_batch(rows, out);
        for o in out.iter_mut() {
            let raw = *o;
            let p = if self.apply_sigmoid {
                1.0 / (1.0 + (-raw).exp())
            } else {
                raw
            };
            let s = self.isotonic(p).clamp(PROB_EPS, 1.0 - PROB_EPS);
            *o = (s / (1.0 - s)).ln();
        }
    }

    /// Calibrated P(signal) in [0,1].
    pub fn predict_proba(&self, x: &[f32]) -> f32 {
        let raw: f32 = self.trees.iter().map(|t| t.eval(x)).sum();
        let p = if self.apply_sigmoid {
            1.0 / (1.0 + (-raw).exp())
        } else {
            raw
        };
        self.isotonic(p)
    }

    /// The additive LLR term `log(s/(1-s))` with s clamped to (eps, 1-eps).
    pub fn predict_logit(&self, x: &[f32]) -> f32 {
        let s = self.predict_proba(x).clamp(PROB_EPS, 1.0 - PROB_EPS);
        (s / (1.0 - s)).ln()
    }

    /// Piecewise-linear interpolation of the isotonic map. Empty map => identity.
    fn isotonic(&self, p: f32) -> f32 {
        let n = self.iso_x.len();
        if n == 0 {
            return p;
        }
        if p <= self.iso_x[0] {
            return self.iso_y[0];
        }
        if p >= self.iso_x[n - 1] {
            return self.iso_y[n - 1];
        }
        // binary search for the segment containing p
        let mut lo = 0usize;
        let mut hi = n - 1;
        while hi - lo > 1 {
            let mid = (lo + hi) / 2;
            if self.iso_x[mid] <= p {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let (x0, x1) = (self.iso_x[lo], self.iso_x[hi]);
        let (y0, y1) = (self.iso_y[lo], self.iso_y[hi]);
        if (x1 - x0).abs() < f32::EPSILON {
            return y0;
        }
        y0 + (y1 - y0) * (p - x0) / (x1 - x0)
    }
}

fn read_i32_vec(c: &mut Cursor<&[u8]>, n: usize) -> Result<Vec<i32>, GbdtError> {
    let mut v = Vec::with_capacity(n);
    for _ in 0..n {
        v.push(c.read_i32::<LittleEndian>()?);
    }
    Ok(v)
}
fn read_f32_vec(c: &mut Cursor<&[u8]>, n: usize) -> Result<Vec<f32>, GbdtError> {
    let mut v = Vec::with_capacity(n);
    for _ in 0..n {
        v.push(c.read_f32::<LittleEndian>()?);
    }
    Ok(v)
}
fn write_i32_vec(b: &mut Vec<u8>, v: &[i32]) {
    for &x in v {
        b.write_i32::<LittleEndian>(x).unwrap();
    }
}
fn write_f32_vec(b: &mut Vec<u8>, v: &[f32]) {
    for &x in v {
        b.write_f32::<LittleEndian>(x).unwrap();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// One tree, one split on feature 0 at threshold 0.5:
    ///   x0 <= 0.5 -> leaf value -1.0 ; else leaf value +2.0
    /// flags: sigmoid ON. isotonic: identity over [0,1] (two breakpoints).
    fn toy_model() -> GbdtPeakModel {
        GbdtPeakModel {
            packed: Default::default(),
            n_features: 1,
            feature_names: indexed_feature_names(1),
            apply_sigmoid: true,
            trees: vec![Tree {
                feature: vec![0, -1, -1],
                threshold: vec![0.5, 0.0, 0.0],
                left: vec![1, -1, -1],
                right: vec![2, -1, -1],
                value: vec![0.0, -1.0, 2.0],
                default_left: vec![1, 1, 1],
            }],
            iso_x: vec![0.0, 1.0],
            iso_y: vec![0.0, 1.0],
        }
    }

    #[test]
    fn packed_evaluation_is_bit_identical_to_tree_eval() {
        // Two-level trees over 3 features with both NaN directions, plus a leaf-only
        // tree; inputs cover each branch, NaN and a short row (missing feature).
        let mk = |f0: i32, t0: f32, f1: i32, t1: f32, dl: u8| Tree {
            feature: vec![f0, f1, -1, -1, -1],
            threshold: vec![t0, t1, 0.0, 0.0, 0.0],
            left: vec![1, 3, -1, -1, -1],
            right: vec![2, 4, -1, -1, -1],
            value: vec![0.0, 0.0, 0.37, -1.25, 2.5e-3],
            default_left: vec![dl, 1 - dl, 0, 0, 0],
        };
        let leaf_only = Tree {
            feature: vec![-1],
            threshold: vec![0.0],
            left: vec![-1],
            right: vec![-1],
            value: vec![0.125],
            default_left: vec![0],
        };
        let m = GbdtPeakModel {
            packed: Default::default(),
            n_features: 3,
            feature_names: indexed_feature_names(3),
            apply_sigmoid: false,
            trees: vec![mk(0, 0.5, 1, -2.0, 1), mk(2, 1e-3, 0, 0.5, 0), leaf_only],
            iso_x: vec![],
            iso_y: vec![],
        };
        let rows: Vec<Vec<f32>> = vec![
            vec![0.1, -3.0, 0.0],
            vec![0.9, 4.0, 2.0],
            vec![f32::NAN, f32::NAN, f32::NAN],
            vec![0.5, -2.0, 1e-3],
            vec![0.7],
        ];
        for x in &rows {
            let reference: f32 = m.trees.iter().map(|t| t.eval(x)).sum();
            assert_eq!(
                m.predict_value(x).to_bits(),
                reference.to_bits(),
                "row {x:?}"
            );
        }
        let refs: Vec<&[f32]> = rows.iter().map(|r| r.as_slice()).collect();
        let mut out = vec![0.0f32; refs.len()];
        m.predict_value_batch(&refs, &mut out);
        for (x, o) in rows.iter().zip(&out) {
            let reference: f32 = m.trees.iter().map(|t| t.eval(x)).sum();
            assert_eq!(o.to_bits(), reference.to_bits(), "batch row {x:?}");
        }
        // Enough rows to fill several lane groups plus a partial one, covering both
        // NaN directions and a short row; the batched walk must match the
        // single-row walk exactly.
        let vals = [0.1f32, 0.9, f32::NAN, 0.5, -2.0, 1e-3, 4.0, -3.0];
        let wide: Vec<Vec<f32>> = (0..19)
            .map(|i| {
                if i == 7 {
                    return vec![vals[i % vals.len()]];
                }
                vec![
                    vals[i % vals.len()],
                    vals[(i * 3 + 1) % vals.len()],
                    vals[(i * 5 + 2) % vals.len()],
                ]
            })
            .collect();
        let refs: Vec<&[f32]> = wide.iter().map(|r| r.as_slice()).collect();
        let mut out = vec![0.0f32; refs.len()];
        m.predict_value_batch(&refs, &mut out);
        for (x, o) in wide.iter().zip(&out) {
            let reference: f32 = m.trees.iter().map(|t| t.eval(x)).sum();
            assert_eq!(o.to_bits(), reference.to_bits(), "batch row {x:?}");
        }
        // An evaluated model still equals an unevaluated copy.
        let mut fresh = m.clone();
        fresh.packed = Default::default();
        assert_eq!(m, fresh);
    }

    /// A random tree of depth up to `depth` over `nf` features, as the
    /// trainer would write it: leaves may sit above the bottom level.
    fn random_tree(rng: &mut u64, depth: u32, nf: i32) -> Tree {
        let mut next = |m: u64| {
            *rng = rng
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (*rng >> 33) % m
        };
        let mut t = Tree {
            feature: vec![],
            threshold: vec![],
            left: vec![],
            right: vec![],
            value: vec![],
            default_left: vec![],
        };
        fn grow(
            t: &mut Tree,
            next: &mut dyn FnMut(u64) -> u64,
            level: u32,
            depth: u32,
            nf: i32,
        ) -> i32 {
            let at = t.feature.len() as i32;
            let leaf = level == depth || (level > 0 && next(4) == 0);
            t.feature
                .push(if leaf { -1 } else { next(nf as u64) as i32 });
            t.threshold.push((next(2000) as f32 - 1000.0) / 500.0);
            t.left.push(-1);
            t.right.push(-1);
            t.value.push((next(2001) as f32 - 1000.0) / 7.0 + 0.001);
            t.default_left.push(next(2) as u8);
            if !leaf {
                let l = grow(t, next, level + 1, depth, nf);
                let r = grow(t, next, level + 1, depth, nf);
                t.left[at as usize] = l;
                t.right[at as usize] = r;
            }
            at
        }
        grow(&mut t, &mut next, 0, depth, nf);
        t
    }

    #[test]
    fn simd_batches_equal_the_tree_walk_bit_for_bit() {
        let mut rng = 12345u64;
        let nf = 17;
        let mut trees: Vec<Tree> = (0..120)
            .map(|i| random_tree(&mut rng, 1 + (i % 8) as u32, nf))
            .collect();
        // A tree deeper than the perfect-tree limit keeps the node walk.
        trees.push(random_tree(&mut rng, 14, nf));
        trees.push(random_tree(&mut rng, 6, nf));
        let m = GbdtPeakModel {
            packed: Default::default(),
            n_features: nf as u32,
            feature_names: indexed_feature_names(nf as u32),
            apply_sigmoid: false,
            trees,
            iso_x: vec![],
            iso_y: vec![],
        };
        let mut next = move || {
            rng = rng
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (rng >> 33) as f32 / (1u64 << 31) as f32
        };
        // Row counts around the block (256) and step (16) sizes; some rows are
        // short or hold NaNs.
        for n in [1usize, 7, 8, 15, 16, 17, 255, 256, 257, 700] {
            let rows: Vec<Vec<f32>> = (0..n)
                .map(|i| {
                    let len = if i % 11 == 3 { 5 } else { nf as usize };
                    (0..len)
                        .map(|_| {
                            let v = next();
                            if v < 0.05 {
                                f32::NAN
                            } else {
                                v * 6.0 - 3.0
                            }
                        })
                        .collect()
                })
                .collect();
            let refs: Vec<&[f32]> = rows.iter().map(|r| r.as_slice()).collect();
            let mut out = vec![0.0f32; n];
            m.predict_value_batch(&refs, &mut out);
            for (x, o) in rows.iter().zip(&out) {
                let want = m.trees.iter().fold(0.0f32, |a, t| a + t.eval(x));
                assert_eq!(o.to_bits(), want.to_bits(), "n = {n}, row {x:?}");
            }
        }
    }

    #[test]
    #[ignore = "timing"]
    fn timing_of_the_batched_walk() {
        let mut rng = 99u64;
        let nf = 17;
        let trees: Vec<Tree> = (0..300).map(|_| random_tree(&mut rng, 6, nf)).collect();
        let m = GbdtPeakModel {
            packed: Default::default(),
            n_features: nf as u32,
            feature_names: indexed_feature_names(nf as u32),
            apply_sigmoid: false,
            trees,
            iso_x: vec![],
            iso_y: vec![],
        };
        let rows: Vec<Vec<f32>> = (0..4096)
            .map(|i| {
                (0..nf)
                    .map(|f| ((i * 31 + f * 17) % 997) as f32 / 150.0 - 3.0)
                    .collect()
            })
            .collect();
        let refs: Vec<&[f32]> = rows.iter().map(|r| r.as_slice()).collect();
        let mut out = vec![0.0f32; refs.len()];
        let reps = 20;
        let t = std::time::Instant::now();
        for _ in 0..reps {
            m.predict_value_batch(&refs, &mut out);
            std::hint::black_box(&out);
        }
        let batched = t.elapsed().as_secs_f64() / (reps * refs.len()) as f64 * 1e9;
        let t = std::time::Instant::now();
        for _ in 0..reps {
            let packed = m.packed();
            out.fill(0.0);
            for p in packed.perfect.iter().flatten() {
                p.add_to(&refs, &mut out);
            }
            std::hint::black_box(&out);
        }
        let scalar = t.elapsed().as_secs_f64() / (reps * refs.len()) as f64 * 1e9;
        eprintln!(
            "300 trees: batched walk {batched:.0} ns/row, scalar lane walk {scalar:.0} ns/row"
        );
    }

    #[test]
    fn roundtrip_bytes() {
        let m = toy_model();
        let bytes = m.to_bytes();
        let back = GbdtPeakModel::from_bytes(&bytes).expect("decode");
        assert_eq!(back, m, "round-trip must preserve every field");
    }

    #[test]
    fn predict_matches_manual() {
        let m = toy_model();
        // x0 = 0.0 -> leaf -1.0 -> raw=-1.0 -> sigmoid(-1)=0.26894 -> iso identity
        let s_lo = m.predict_proba(&[0.0]);
        assert!((s_lo - 0.2689414).abs() < 1e-5, "got {s_lo}");
        // x0 = 1.0 -> leaf +2.0 -> sigmoid(2)=0.880797
        let s_hi = m.predict_proba(&[1.0]);
        assert!((s_hi - 0.8807971).abs() < 1e-5, "got {s_hi}");
        // logit(s) recovers the raw sum for the identity isotonic map.
        let lg = m.predict_logit(&[1.0]);
        assert!((lg - 2.0).abs() < 1e-4, "logit got {lg}");
    }

    #[test]
    fn empty_iso_is_identity() {
        // No isotonic breakpoints -> calibrated prob == sigmoid(raw).
        let mut m = toy_model();
        m.iso_x.clear();
        m.iso_y.clear();
        let s = m.predict_proba(&[1.0]);
        assert!((s - 0.8807971).abs() < 1e-5);
    }

    #[test]
    fn predict_value_returns_raw_tree_sum() {
        // toy_model: x0 <= 0.5 → leaf -1.0; else → leaf +2.0
        // predict_value must return the raw tree sum without sigmoid or isotonic.
        let m = toy_model();
        // x0 = 0.0 → leaf -1.0; raw sum = -1.0 (no sigmoid/iso applied)
        let v_lo = m.predict_value(&[0.0]);
        assert!(
            (v_lo - (-1.0)).abs() < 1e-6,
            "predict_value([0.0]) expected -1.0, got {v_lo}"
        );
        // x0 = 1.0 → leaf +2.0; raw sum = 2.0
        let v_hi = m.predict_value(&[1.0]);
        assert!(
            (v_hi - 2.0).abs() < 1e-6,
            "predict_value([1.0]) expected 2.0, got {v_hi}"
        );
    }

    #[test]
    fn predict_value_ignores_sigmoid_flag() {
        // Even when apply_sigmoid=false, predict_value must still return raw sum (not apply isotonic).
        let mut m = toy_model();
        m.apply_sigmoid = false;
        // Set a non-identity isotonic map so we can verify it's NOT applied.
        m.iso_x = vec![0.0, 1.0];
        m.iso_y = vec![0.5, 0.5]; // constant 0.5 — if applied, predict_value would be 0.5
        let v = m.predict_value(&[1.0]);
        assert!(
            (v - 2.0).abs() < 1e-6,
            "predict_value must not apply isotonic map; expected 2.0, got {v}"
        );
    }

    fn named(names: &[&str]) -> GbdtPeakModel {
        let mut m = toy_model();
        m.n_features = names.len() as u32;
        m.feature_names = names.iter().map(|n| n.to_string()).collect();
        m
    }

    #[test]
    fn validate_features_accepts_the_same_names_in_the_same_order() {
        let m = named(&["a", "b", "c"]);
        assert!(m.validate_features(&["a", "b", "c"], "peak").is_ok());
    }

    #[test]
    fn validate_features_rejects_a_count_difference() {
        let m = named(&["a", "b", "c"]);
        assert!(matches!(
            m.validate_features(&["a", "b"], "peak"),
            Err(GbdtError::FeatureCountMismatch {
                kind: "peak",
                expected: 2,
                got: 3
            })
        ));
        assert!(matches!(
            m.validate_features(&[], "rich-ion"),
            Err(GbdtError::FeatureCountMismatch {
                expected: 0,
                got: 3,
                ..
            })
        ));
    }

    #[test]
    fn validate_features_names_the_first_differing_position() {
        // Positions 1 and 2 are swapped; the first one is reported.
        let m = named(&["a", "c", "b"]);
        match m.validate_features(&["a", "b", "c"], "frag-intensity") {
            Err(GbdtError::FeatureNameMismatch {
                kind,
                position,
                expected,
                got,
            }) => {
                assert_eq!(kind, "frag-intensity");
                assert_eq!(position, 1);
                assert_eq!(expected, "b");
                assert_eq!(got, "c");
            }
            other => panic!("expected a name mismatch, got {other:?}"),
        }
        let msg = m
            .validate_features(&["a", "b", "c"], "frag-intensity")
            .unwrap_err()
            .to_string();
        assert!(msg.contains("feature 1") && msg.contains("'c'") && msg.contains("'b'"));
    }

    #[test]
    fn feature_names_round_trip_through_the_blob() {
        let m = named(&["log_intensity", "mz", "mass_defect"]);
        let back = GbdtPeakModel::from_bytes(&m.to_bytes()).expect("round trip");
        assert_eq!(back.feature_names, m.feature_names);
        assert_eq!(back, m);
        assert_eq!(back.to_bytes(), m.to_bytes());
    }

    #[test]
    fn blob_without_the_name_trailer_is_rejected() {
        let m = toy_model();
        let mut bytes = m.to_bytes();
        // Drop the trailer (u32 count + one length-prefixed name): the blob now
        // ends after the calibration footer, as a version-1 blob would.
        let trailer = 4 + 4 + m.feature_names[0].len();
        bytes.truncate(bytes.len() - trailer);
        assert!(matches!(
            GbdtPeakModel::from_bytes(&bytes),
            Err(GbdtError::Io(_))
        ));
    }

    #[test]
    fn name_count_must_equal_feature_count() {
        let mut m = toy_model();
        m.feature_names.push("extra".to_string());
        assert!(matches!(
            GbdtPeakModel::from_bytes(&m.to_bytes()),
            Err(GbdtError::Malformed(_))
        ));
    }

    #[test]
    fn trailing_bytes_after_the_names_are_rejected() {
        let mut bytes = toy_model().to_bytes();
        bytes.push(0);
        assert!(matches!(
            GbdtPeakModel::from_bytes(&bytes),
            Err(GbdtError::Malformed(_))
        ));
    }

    #[test]
    fn previous_format_version_is_rejected() {
        let mut bytes = toy_model().to_bytes();
        bytes[4..8].copy_from_slice(&1u32.to_le_bytes());
        assert!(matches!(
            GbdtPeakModel::from_bytes(&bytes),
            Err(GbdtError::BadVersion(1))
        ));
    }

    #[test]
    fn malformed_child_index_is_rejected() {
        // An internal node (feature 0) whose left child points out of range.
        let bad = GbdtPeakModel {
            packed: Default::default(),
            n_features: 1,
            feature_names: indexed_feature_names(1),
            apply_sigmoid: true,
            trees: vec![Tree {
                feature: vec![0, -1], // node 0 internal, node 1 leaf
                threshold: vec![0.5, 0.0],
                left: vec![99, -1], // 99 is out of range (n=2)
                right: vec![1, -1],
                value: vec![0.0, 1.0],
                default_left: vec![1, 1],
            }],
            iso_x: vec![],
            iso_y: vec![],
        };
        let bytes = bad.to_bytes();
        assert!(
            matches!(
                GbdtPeakModel::from_bytes(&bytes),
                Err(GbdtError::Malformed(_))
            ),
            "out-of-range child must be rejected at decode"
        );
    }
}
