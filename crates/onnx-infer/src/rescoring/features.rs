//! Feature name tables aligned with quantms-rescoring / ms2rescore.
//!
//! Names are the *short* forms written into PIN / idparquet `psm_metavalues`
//! (matching `quantmsrescore.constants` value side). Tool-prefixed aliases
//! (`MS2PIP:…`, `DeepLC:…`) live in the Python package; Rust consumers should
//! use these short names so andes PIN and quantms-rescoring never diverge.

/// Core AlphaPeptDeep / MS2PIP-style spectral-match features used by #92 A/B.
///
/// Subset of `quantmsrescore.alphapeptdeep.AlphaPeptDeepFeatureGenerator.feature_names`
/// — the five columns issue #92 plans to append first, plus a few high-signal
/// companions already emitted by the Python path.
pub const ALPHAPEPTDEEP_CORE_FEATURES: &[&str] = &[
    "spec_pearson_norm",
    "ionb_pearson_norm",
    "iony_pearson_norm",
    "spec_mse_norm",
    "dotprod_norm",
    "cos_norm",
    "rt_diff",
    "predicted_retention_time",
];

/// DeepLC retention-time features (`quantmsrescore.constants.DEEPLC_FEATURES` values).
pub const DEEPLC_CORE_FEATURES: &[&str] = &[
    "observed_retention_time",
    "predicted_retention_time",
    "rt_diff",
];

/// Spectrum-quality features from `rescoring spectrum2feature` / SNR path.
pub const SPECTRUM_CORE_FEATURES: &[&str] = &[
    "snr",
    "spectral_entropy",
    "tic_top_fraction",
    "mz_intensity_weighted_std",
];

/// AlphaPeptDeep core feature names (owned `String`s for PIN builders).
pub fn alphapeptdeep_feature_names() -> Vec<&'static str> {
    ALPHAPEPTDEEP_CORE_FEATURES.to_vec()
}

/// DeepLC core feature names.
pub fn deeplc_feature_names() -> Vec<&'static str> {
    DEEPLC_CORE_FEATURES.to_vec()
}

/// Spectrum-quality feature names.
pub fn spectrum_feature_names() -> Vec<&'static str> {
    SPECTRUM_CORE_FEATURES.to_vec()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn core_sets_are_non_empty_and_unique() {
        for set in [
            ALPHAPEPTDEEP_CORE_FEATURES,
            DEEPLC_CORE_FEATURES,
            SPECTRUM_CORE_FEATURES,
        ] {
            assert!(!set.is_empty());
            let mut v = set.to_vec();
            v.sort_unstable();
            v.dedup();
            assert_eq!(v.len(), set.len(), "duplicate feature name in {set:?}");
        }
    }
}
