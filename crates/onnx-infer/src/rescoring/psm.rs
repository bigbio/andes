//! Search-engine–agnostic PSM handles for rescoring batches.
//!
//! Deliberately thin: no dependency on andes `PsmMatch` or other engines' types.
//! Callers map their candidate structs into [`CandidateKey`] + feature rows.
//! This keeps `onnx-infer` extractable for quantms while staying compatible
//! with andes PIN / idparquet writers that already share `psm_feature_values`.

use std::collections::BTreeMap;

use crate::error::{Error, Result};
use crate::tensor::NamedTensor;

/// Stable identity for one PSM / candidate inside a rescoring batch.
///
/// Fields mirror the join keys used by quantms-rescoring idparquet + mzML
/// matching (`spectrum_ref` ≈ scan/native id; `rank` is 1-based engine rank).
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct CandidateKey {
    /// Spectrum native id / SpecId base (no rank suffix).
    pub spectrum_ref: String,
    /// 1-based rank among retained candidates for that spectrum.
    pub rank: u32,
    /// Precursor charge.
    pub charge: u8,
    /// Peptide sequence (ProForma or plain AA; caller-defined).
    pub peptide: String,
}

/// One row of named rescoring features for a candidate.
#[derive(Debug, Clone, PartialEq)]
pub struct FeatureRow {
    /// Candidate identity.
    pub key: CandidateKey,
    /// Feature name → value (PIN / idparquet column names).
    pub values: BTreeMap<String, f64>,
}

/// How many engine candidates per spectrum enter ONNX inference.
///
/// **Speed contract:** never score the full search lattice — only retained
/// top-N (andes default search `n` / Percolator PIN depth).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TopNPolicy {
    /// Maximum candidates per spectrum (must be ≥ 1).
    pub max_per_spectrum: usize,
}

impl TopNPolicy {
    /// Construct; rejects zero.
    pub fn new(max_per_spectrum: usize) -> Result<Self> {
        if max_per_spectrum == 0 {
            return Err(Error::InvalidArg(
                "TopNPolicy.max_per_spectrum must be ≥ 1".into(),
            ));
        }
        Ok(Self { max_per_spectrum })
    }

    /// Keep the first `max_per_spectrum` candidates per `spectrum_ref`.
    ///
    /// Input must already be sorted by spectrum then ascending rank (1,2,…).
    pub fn filter_sorted_by_rank<'a, T, F>(&self, items: &'a [T], spectrum_of: F) -> Vec<&'a T>
    where
        F: Fn(&T) -> &str,
    {
        let mut out = Vec::new();
        let mut cur_spec: Option<&str> = None;
        let mut kept = 0usize;
        for item in items {
            let spec = spectrum_of(item);
            if cur_spec != Some(spec) {
                cur_spec = Some(spec);
                kept = 0;
            }
            if kept < self.max_per_spectrum {
                out.push(item);
                kept += 1;
            }
        }
        out
    }
}

/// A dense batch ready for ONNX: keys aligned with matrix rows.
#[derive(Debug, Clone, PartialEq)]
pub struct RescoringBatch {
    /// Row i corresponds to `features` row i.
    pub keys: Vec<CandidateKey>,
    /// Contiguous `[n_rows, n_features]` feature matrix.
    pub features: NamedTensor,
    /// Column names matching `features` axis-1 (PIN feature names).
    pub feature_names: Vec<String>,
}

impl RescoringBatch {
    /// Pack rows into a single `[n, f]` tensor named `features` by default.
    pub fn from_rows(
        rows: Vec<FeatureRow>,
        feature_names: &[&str],
        tensor_name: impl Into<String>,
    ) -> Result<Self> {
        let n = rows.len();
        let f = feature_names.len();
        let mut data = vec![0f32; n * f];
        for (i, row) in rows.iter().enumerate() {
            for (j, name) in feature_names.iter().enumerate() {
                let v = row.values.get(*name).copied().unwrap_or(0.0) as f32;
                if !v.is_finite() {
                    return Err(Error::InvalidArg(format!(
                        "feature `{name}` of row {i} is not a finite f32 ({v}); \
                         impute or drop non-finite values before packing"
                    )));
                }
                data[i * f + j] = v;
            }
        }
        let features = NamedTensor::new(tensor_name, vec![n, f], data)?;
        Ok(Self {
            keys: rows.into_iter().map(|r| r.key).collect(),
            features,
            feature_names: feature_names.iter().map(|s| (*s).to_string()).collect(),
        })
    }

    /// Number of candidates in the batch.
    pub fn len(&self) -> usize {
        self.keys.len()
    }

    /// True when empty.
    pub fn is_empty(&self) -> bool {
        self.keys.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn top_n_keeps_first_k_per_spectrum() {
        let policy = TopNPolicy::new(2).unwrap();
        let items = vec![("s1", 1), ("s1", 2), ("s1", 3), ("s2", 1), ("s2", 2)];
        let kept = policy.filter_sorted_by_rank(&items, |x| x.0);
        assert_eq!(kept, vec![&("s1", 1), &("s1", 2), &("s2", 1), &("s2", 2)]);
    }

    #[test]
    fn batch_packs_missing_as_zero() {
        let rows = vec![FeatureRow {
            key: CandidateKey {
                spectrum_ref: "scan=1".into(),
                rank: 1,
                charge: 2,
                peptide: "PEPTIDE".into(),
            },
            values: BTreeMap::from([("spec_pearson_norm".into(), 0.9)]),
        }];
        let batch =
            RescoringBatch::from_rows(rows, &["spec_pearson_norm", "rt_diff"], "features").unwrap();
        assert_eq!(batch.features.shape, vec![1, 2]);
        assert_eq!(batch.features.data, vec![0.9, 0.0]);
    }

    #[test]
    fn non_finite_features_are_rejected() {
        for bad in [f64::NAN, f64::INFINITY, 1e300] {
            let rows = vec![FeatureRow {
                key: CandidateKey {
                    spectrum_ref: "scan=1".into(),
                    rank: 1,
                    charge: 2,
                    peptide: "PEPTIDE".into(),
                },
                values: BTreeMap::from([("rt_diff".into(), bad)]),
            }];
            assert!(
                RescoringBatch::from_rows(rows, &["rt_diff"], "features").is_err(),
                "{bad} must be rejected"
            );
        }
    }
}
