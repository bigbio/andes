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

pub const DEEPLC_CORE_FEATURES: &[&str] = &[
    "observed_retention_time",
    "predicted_retention_time",
    "rt_diff",
];

pub const SPECTRUM_CORE_FEATURES: &[&str] = &[
    "snr",
    "spectral_entropy",
    "tic_top_fraction",
    "mz_intensity_weighted_std",
];

pub fn alphapeptdeep_feature_names() -> Vec<&'static str> {
    ALPHAPEPTDEEP_CORE_FEATURES.to_vec()
}

pub fn deeplc_feature_names() -> Vec<&'static str> {
    DEEPLC_CORE_FEATURES.to_vec()
}

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
