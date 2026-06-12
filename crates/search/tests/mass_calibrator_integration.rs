//! Integration tests for the MassCalibrator pre-pass.
//!
//! Closes the "no isolated cal integration test" gap documented in
//! `docs/superpowers/specs/2026-05-25-precursor-cal-ship-design.md`.
//!
//! Asserts the contracts of each calibrator helper layer (threshold skip,
//! mode-aware shift application, mass-correction round-trip, tightening
//! bounds, SpecKey expansion). End-to-end benchmark validation remains the
//! harness's responsibility.

use std::collections::HashMap;
use rustc_hash::FxHashMap;

use model::{AminoAcidSetBuilder, Protein, ProteinDb, Spectrum};
use scoring_crate::param_model::{IonType, Partition, SpecDataType};
use scoring_crate::{Param, RankScorer};
use model::activation::ActivationMethod;
use model::instrument::InstrumentType;
use model::protocol::Protocol;
use model::tolerance::Tolerance;
use search::precursor_cal::{
    adjusted_observed_neutral_mass, constants as cal_constants, tightened_tolerance_ppm,
    PrecursorCalMode,
};
use search::{
    apply_shift_for_mode, build_spec_keys, learn_calibration_stats, CalibrationStats,
    PreparedSearch, SearchIndex, SearchParams, SpecKey,
};

fn tiny_scorer() -> RankScorer {
    let part = Partition { charge: 2, parent_mass: 1000.0, seg_num: 0 };
    let prefix1 = IonType::Prefix { charge: 1, offset_bits: 0.0_f32.to_bits() };
    let suffix1 = IonType::Suffix { charge: 1, offset_bits: 0.0_f32.to_bits() };
    let noise = IonType::Noise;

    let mut ion_table = FxHashMap::default();
    ion_table.insert(prefix1, vec![0.5_f32, 0.1, 0.05, 0.01]);
    ion_table.insert(suffix1, vec![0.5_f32, 0.1, 0.05, 0.01]);
    ion_table.insert(noise, vec![0.1_f32, 0.05, 0.02, 0.01]);

    let mut rank_dist_table = FxHashMap::default();
    rank_dist_table.insert(part, ion_table);

    let mut frag_off_table = FxHashMap::default();
    frag_off_table.insert(part, vec![]);

    let mut param = Param {
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
        partitions: vec![part],
        num_precursor_off: 0,
        precursor_off_map: FxHashMap::default(),
        frag_off_table,
        max_rank: 3,
        rank_dist_table,
        error_scaling_factor: 0,
        ion_err_dist_table: FxHashMap::default(),
        noise_err_dist_table: FxHashMap::default(),
        ion_existence_table: FxHashMap::default(),
        partition_ion_types_cache: FxHashMap::default(),
    };
    param.rebuild_cache();
    RankScorer::new(&param)
}

#[test]
fn threshold_skip_returns_empty_stats() {
    // Fewer than MIN_SPECKEYS_FOR_PREPASS spec keys → stats unreliable.
    let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
    let params = SearchParams::default_tryptic(aa_set);
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: String::new(),
            sequence: b"MKWVTFISLLR".to_vec(),
        }],
    };
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let scorer = tiny_scorer();
    let prepared = PreparedSearch::prepare(&idx, &params, &scorer, 0.5, "XXX");

    let too_few_keys: Vec<SpecKey> = (0..(cal_constants::MIN_SPECKEYS_FOR_PREPASS - 1))
        .map(|i| SpecKey { spectrum_idx: i, charge: 2 })
        .collect();

    let stats = learn_calibration_stats(
        &too_few_keys,
        &HashMap::new(),
        &prepared,
        &params,
    );
    assert!(!stats.has_reliable_stats());
    assert_eq!(stats.shift_ppm, 0.0);
    assert_eq!(stats.confident_psm_count, 0);
}

#[test]
fn apply_shift_for_mode_respects_mode_and_reliability() {
    let unreliable = CalibrationStats { shift_ppm: 5.0, ..Default::default() };
    assert_eq!(apply_shift_for_mode(PrecursorCalMode::Auto, unreliable), 0.0);
    assert_eq!(apply_shift_for_mode(PrecursorCalMode::On, unreliable), 5.0);
    assert_eq!(apply_shift_for_mode(PrecursorCalMode::Off, unreliable), 0.0);

    let reliable = CalibrationStats {
        shift_ppm: 5.0,
        robust_sigma_ppm: 0.5,
        confident_psm_count: 200,
        ..Default::default()
    };
    assert_eq!(apply_shift_for_mode(PrecursorCalMode::Auto, reliable), 5.0);
    assert_eq!(apply_shift_for_mode(PrecursorCalMode::On, reliable), 5.0);
    assert_eq!(apply_shift_for_mode(PrecursorCalMode::Off, reliable), 0.0);
}

#[test]
fn adjusted_observed_neutral_mass_round_trips_a_known_bias() {
    let theoretical = 1000.0_f64;
    let observed_with_bias = theoretical * (1.0 + 5e-6);
    let corrected = adjusted_observed_neutral_mass(observed_with_bias, 5.0);
    assert!(
        (corrected - theoretical).abs() < 1e-6,
        "after +5 ppm correction, observed should match theoretical; got delta {}",
        corrected - theoretical
    );
}

#[test]
fn tightened_tolerance_respects_user_floor_and_sigma() {
    assert!((tightened_tolerance_ppm(10.0, 0.2, 3.0, 2.0, 0.5) - 2.0).abs() < 1e-9);
    assert!((tightened_tolerance_ppm(10.0, 1.0, 3.0, 2.0, 0.5) - 3.5).abs() < 1e-9);
    assert!((tightened_tolerance_ppm(1.5, 1.0, 3.0, 2.0, 0.5) - 1.5).abs() < 1e-9);
}

#[test]
fn build_spec_keys_skips_below_min_peaks_and_expands_missing_charge() {
    fn spec_with_peaks(precursor_mz: f64, charge: Option<i32>, n_peaks: usize) -> Spectrum {
        Spectrum {
            title: "t".into(),
            precursor_mz,
            precursor_intensity: None,
            precursor_charge: charge,
            rt_seconds: None,
            scan: None,
            peaks: vec![(0.0, 0.0); n_peaks],
            activation_method: None,
            isolation_lower_offset: None,
            isolation_upper_offset: None,
        }
    }
    let spectra = vec![
        spec_with_peaks(500.0, Some(2), 50),
        spec_with_peaks(700.0, None, 50),
        spec_with_peaks(800.0, Some(2), 5),
    ];
    let keys = build_spec_keys(&spectra, &(2..=4), 10);
    assert_eq!(keys.len(), 4, "expected 4 keys, got {keys:?}");
    assert!(keys.iter().any(|k| k.spectrum_idx == 0 && k.charge == 2));
    assert!(keys.iter().any(|k| k.spectrum_idx == 1 && k.charge == 3));
    assert!(!keys.iter().any(|k| k.spectrum_idx == 2));
}
