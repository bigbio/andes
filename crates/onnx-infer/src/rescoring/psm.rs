use std::collections::BTreeMap;

use crate::error::{Error, Result};
use crate::tensor::NamedTensor;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct CandidateKey {
    pub spectrum_ref: String,
    pub rank: u32,
    pub charge: u8,
    pub peptide: String,
}

#[derive(Debug, Clone, PartialEq)]
pub struct FeatureRow {
    pub key: CandidateKey,
    pub values: BTreeMap<String, f64>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TopNPolicy {
    pub max_per_spectrum: usize,
}

impl TopNPolicy {
    pub fn new(max_per_spectrum: usize) -> Result<Self> {
        if max_per_spectrum == 0 {
            return Err(Error::InvalidArg(
                "TopNPolicy.max_per_spectrum must be ≥ 1".into(),
            ));
        }
        Ok(Self { max_per_spectrum })
    }

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

#[derive(Debug, Clone, PartialEq)]
pub struct RescoringBatch {
    pub keys: Vec<CandidateKey>,
    pub features: NamedTensor,
    pub feature_names: Vec<String>,
}

impl RescoringBatch {
    pub fn from_rows(
        rows: Vec<FeatureRow>,
        feature_names: &[&str],
        tensor_name: impl Into<String>,
    ) -> Result<Self> {
        let n = rows.len();
        let f = feature_names.len();
        let mut data =
            vec![0f32; crate::tensor::checked_element_count(&[n, f], "rescoring batch")?];
        for (i, row) in rows.iter().enumerate() {
            for (j, name) in feature_names.iter().enumerate() {
                let v = row.values.get(*name).copied().ok_or_else(|| {
                    Error::InvalidArg(format!(
                        "row {i} ({}) has no feature `{name}`",
                        row.key.spectrum_ref
                    ))
                })? as f32;
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

    pub fn len(&self) -> usize {
        self.keys.len()
    }

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
    fn batch_rejects_a_missing_feature() {
        let rows = vec![FeatureRow {
            key: CandidateKey {
                spectrum_ref: "scan=1".into(),
                rank: 1,
                charge: 2,
                peptide: "PEPTIDE".into(),
            },
            values: BTreeMap::from([("spec_pearson_norm".into(), 0.9)]),
        }];
        let err =
            RescoringBatch::from_rows(rows.clone(), &["spec_pearson_norm", "rt_diff"], "features")
                .unwrap_err();
        assert!(err.to_string().contains("rt_diff"), "{err}");
        let batch = RescoringBatch::from_rows(rows, &["spec_pearson_norm"], "features").unwrap();
        assert_eq!(batch.features.shape, vec![1, 1]);
        assert_eq!(batch.features.data, vec![0.9]);
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
