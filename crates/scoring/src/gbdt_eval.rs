//! Hand-rolled struct-of-arrays GBDT evaluator (zero native deps).
//!
//! Decodes the `AGBD` v1 blob (produced offline by `training/gbdt/transcode.py`
//! from a LightGBM binary classifier + scikit-learn IsotonicRegression) and
//! evaluates it on a peptide-AGNOSTIC per-peak feature vector
//! (`crate::peak_features`). Output is `log(s/(1-s))`, the additive LLR term the
//! rank scorer folds in (`scored_spectrum`).

use std::io::{Cursor, Read};

use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};

use thiserror::Error;

const MAGIC: &[u8; 4] = b"AGBD";
const FORMAT_VERSION: u32 = 1;
const FLAG_SIGMOID: u32 = 1;
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
}

/// Reject implausibly large counts from a corrupt/version-skewed blob before
/// allocating, so decoding fails cleanly instead of OOM-aborting. Generous
/// bounds — a real per-peak GBDT is far smaller.
fn checked_count(n: u32, what: &str, max: u32) -> Result<usize, GbdtError> {
    if n > max {
        return Err(GbdtError::Malformed(format!("{what} count {n} exceeds max {max}")));
    }
    Ok(n as usize)
}

/// One regression tree in struct-of-arrays layout. All vecs have length
/// `n_nodes`; node 0 is the root.
#[derive(Debug, Clone, PartialEq)]
pub struct Tree {
    pub feature: Vec<i32>,      // -1 => leaf
    pub threshold: Vec<f32>,
    pub left: Vec<i32>,
    pub right: Vec<i32>,
    pub value: Vec<f32>,        // leaf output
    pub default_left: Vec<u8>,  // 1 => NaN feature descends left
}

impl Tree {
    /// Sum of the leaf value reached for `x`.
    fn eval(&self, x: &[f32]) -> f32 {
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
            node = if go_left { self.left[node] } else { self.right[node] } as usize;
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
        if self.threshold.len() != n || self.left.len() != n || self.right.len() != n
            || self.value.len() != n || self.default_left.len() != n {
            return Err(GbdtError::Malformed("tree SoA arrays have mismatched lengths".into()));
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

#[derive(Debug, Clone, PartialEq)]
pub struct GbdtPeakModel {
    pub n_features: u32,
    pub apply_sigmoid: bool,
    pub trees: Vec<Tree>,
    pub iso_x: Vec<f32>,
    pub iso_y: Vec<f32>,
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
            let tree = Tree { feature, threshold, left, right, value, default_left };
            tree.validate()?;
            trees.push(tree);
        }
        let n_iso = checked_count(c.read_u32::<LittleEndian>()?, "n_iso", 10_000_000)?;
        let iso_x = read_f32_vec(&mut c, n_iso)?;
        let iso_y = read_f32_vec(&mut c, n_iso)?;
        Ok(Self { n_features, apply_sigmoid, trees, iso_x, iso_y })
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        let mut b = Vec::new();
        b.extend_from_slice(MAGIC);
        b.write_u32::<LittleEndian>(FORMAT_VERSION).unwrap();
        b.write_u32::<LittleEndian>(self.n_features).unwrap();
        b.write_u32::<LittleEndian>(if self.apply_sigmoid { FLAG_SIGMOID } else { 0 }).unwrap();
        b.write_u32::<LittleEndian>(self.trees.len() as u32).unwrap();
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
        b.write_u32::<LittleEndian>(self.iso_x.len() as u32).unwrap();
        write_f32_vec(&mut b, &self.iso_x);
        write_f32_vec(&mut b, &self.iso_y);
        b
    }

    /// Calibrated P(signal) in [0,1].
    pub fn predict_proba(&self, x: &[f32]) -> f32 {
        let raw: f32 = self.trees.iter().map(|t| t.eval(x)).sum();
        let p = if self.apply_sigmoid { 1.0 / (1.0 + (-raw).exp()) } else { raw };
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
            n_features: 1,
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
    fn malformed_child_index_is_rejected() {
        // An internal node (feature 0) whose left child points out of range.
        let bad = GbdtPeakModel {
            n_features: 1,
            apply_sigmoid: true,
            trees: vec![Tree {
                feature: vec![0, -1],     // node 0 internal, node 1 leaf
                threshold: vec![0.5, 0.0],
                left: vec![99, -1],       // 99 is out of range (n=2)
                right: vec![1, -1],
                value: vec![0.0, 1.0],
                default_left: vec![1, 1],
            }],
            iso_x: vec![],
            iso_y: vec![],
        };
        let bytes = bad.to_bytes();
        assert!(matches!(GbdtPeakModel::from_bytes(&bytes), Err(GbdtError::Malformed(_))),
            "out-of-range child must be rejected at decode");
    }
}
