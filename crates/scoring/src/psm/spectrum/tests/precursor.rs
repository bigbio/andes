//! Tests of the precursor peak filter.

use super::*;
use crate::params::{Param, PrecursorOffsetFrequency, SpecDataType};
use model::activation::ActivationMethod;
use model::instrument::InstrumentType;
use model::protocol::Protocol;
use model::tolerance::Tolerance;
use rustc_hash::FxHashMap;

fn param_with_precursor_filter() -> Param {
    let mut precursor_off_map: FxHashMap<i32, Vec<PrecursorOffsetFrequency>> = FxHashMap::default();
    precursor_off_map.insert(
        2,
        vec![PrecursorOffsetFrequency {
            reduced_charge: 2,
            offset: 0.0,
            tolerance: Tolerance::Da(0.5),
            frequency: 1.0,
        }],
    );

    Param {
        version: 10001,
        data_type: SpecDataType {
            activation: ActivationMethod::HCD,
            instrument: InstrumentType::QExactive,
            enzyme: None,
            protocol: Protocol::Automatic,
        },
        mme: Tolerance::Ppm(20.0),
        apply_deconvolution: false,
        deconvolution_error_tolerance: 0.0,
        charge_hist: vec![(2, 100)],
        min_charge: 2,
        max_charge: 2,
        num_segments: 1,
        partitions: vec![],
        num_precursor_off: 1,
        precursor_off_map,
        frag_off_table: FxHashMap::default(),
        max_rank: 3,
        rank_dist_table: FxHashMap::default(),
        error_scaling_factor: 0,
        ion_err_dist_table: FxHashMap::default(),
        noise_err_dist_table: FxHashMap::default(),
        ion_existence_table: FxHashMap::default(),
        partition_ion_types_cache: FxHashMap::default(),
        gbdt_peak_model: None,
        frag_intensity_model: None,
        rich_ion_model: None,
    }
}

fn make_spec(precursor_mz: f64, peaks: &[(f64, f32)]) -> Spectrum {
    Spectrum {
        title: "test".into(),
        precursor_mz,
        precursor_intensity: None,
        precursor_charge: Some(2),
        rt_seconds: None,
        scan: None,
        peaks: peaks.to_vec(),
        activation_method: None,
        isolation_lower_offset: None,
        isolation_upper_offset: None,
        isolation_target_mz: None,
        mobility_window: None,
    }
}

fn param_with_precursor_filter_rc0() -> Param {
    let mut precursor_off_map: FxHashMap<i32, Vec<PrecursorOffsetFrequency>> = FxHashMap::default();
    precursor_off_map.insert(
        2,
        vec![PrecursorOffsetFrequency {
            reduced_charge: 0,
            offset: 0.0,
            tolerance: Tolerance::Da(0.5),
            frequency: 1.0,
        }],
    );

    Param {
        version: 10001,
        data_type: SpecDataType {
            activation: ActivationMethod::HCD,
            instrument: InstrumentType::QExactive,
            enzyme: None,
            protocol: Protocol::Automatic,
        },
        mme: Tolerance::Ppm(20.0),
        apply_deconvolution: false,
        deconvolution_error_tolerance: 0.0,
        charge_hist: vec![(2, 100)],
        min_charge: 2,
        max_charge: 2,
        num_segments: 1,
        partitions: vec![],
        num_precursor_off: 1,
        precursor_off_map,
        frag_off_table: FxHashMap::default(),
        max_rank: 3,
        rank_dist_table: FxHashMap::default(),
        error_scaling_factor: 0,
        ion_err_dist_table: FxHashMap::default(),
        noise_err_dist_table: FxHashMap::default(),
        ion_existence_table: FxHashMap::default(),
        partition_ion_types_cache: FxHashMap::default(),
        gbdt_peak_model: None,
        frag_intensity_model: None,
        rich_ion_model: None,
    }
}

#[test]
fn precursor_peak_is_filtered_out() {
    let s = make_spec(500.0, &[(100.0, 1.0), (300.0, 5.0), (500.0, 100.0)]);
    let param = param_with_precursor_filter_rc0();
    let scorer = RankScorer::new(&param);
    let ss = ScoredSpectrum::new(&s, &scorer, 2);

    assert!(
        ss.nearest_peak_rank(500.0, 0.1).is_none(),
        "precursor peak at 500.0 should be filtered, but a peak at that m/z was found"
    );

    assert_eq!(ss.nearest_peak_rank(300.0, 0.1), Some(1));
    assert_eq!(ss.nearest_peak_rank(100.0, 0.1), Some(2));
}

#[test]
fn non_precursor_peaks_kept() {
    let s = make_spec(500.0, &[(100.0, 1.0), (200.0, 50.0), (300.0, 5.0)]);
    let param = param_with_precursor_filter_rc0();
    let scorer = RankScorer::new(&param);
    let ss = ScoredSpectrum::new(&s, &scorer, 2);

    assert_eq!(ss.peak_count_after_filtering(), 3);
    assert_eq!(ss.nearest_peak_rank(200.0, 0.1), Some(1));
}

#[test]
fn missing_precursor_off_map_falls_back_to_unfiltered() {
    let mut param = param_with_precursor_filter_rc0();
    param.precursor_off_map.clear();
    let s = make_spec(500.0, &[(100.0, 1.0), (500.0, 100.0)]);
    let scorer = RankScorer::new(&param);
    let ss = ScoredSpectrum::new(&s, &scorer, 2);
    assert_eq!(ss.peak_count_after_filtering(), 2);
}

#[test]
fn invalid_reduced_charge_skipped() {
    let param = param_with_precursor_filter();
    let s = make_spec(500.0, &[(100.0, 1.0), (500.0, 100.0)]);
    let scorer = RankScorer::new(&param);
    let ss = ScoredSpectrum::new(&s, &scorer, 2);
    assert_eq!(ss.peak_count_after_filtering(), 2);
}
