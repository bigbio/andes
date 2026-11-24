//! Tests of spectrum scoring: ranks, node and edge scores, ion-match facts and deconvolution.

use super::deconvolve::main_ion_from_param;
use super::matches::{matched_peak_mz, nearest_peak_rank_and_mz_in, nearest_peak_rank_in};
use crate::params::Param;

use super::*;

#[test]
fn main_ion_deterministic_with_placeholder_zero_frequencies() {
    use crate::params::{FragmentOffsetFrequency, IonType, Partition};
    let part = Partition {
        charge: 2,
        parent_mass: 1500.0,
        seg_num: 0,
    };
    let zero = |ion| FragmentOffsetFrequency {
        ion_type: ion,
        frequency: 0.0,
    };
    let b1 = IonType::Prefix {
        charge: 1,
        offset_bits: 1.0_f32.to_bits(),
        loss_class: 0,
    };
    let y1 = IonType::Suffix {
        charge: 1,
        offset_bits: 19.0_f32.to_bits(),
        loss_class: 0,
    };
    let b2 = IonType::Prefix {
        charge: 2,
        offset_bits: 1.0_f32.to_bits(),
        loss_class: 0,
    };
    let y2 = IonType::Suffix {
        charge: 2,
        offset_bits: 19.0_f32.to_bits(),
        loss_class: 0,
    };

    let mut p = crate::testutil::tiny_param();
    p.num_segments = 1;
    p.frag_off_table.clear();
    p.frag_off_table.insert(
        part,
        vec![zero(b1), zero(y1), zero(b2), zero(y2), zero(IonType::Noise)],
    );
    p.rebuild_cache();

    let a = main_ion_from_param(&p, part);
    let b = main_ion_from_param(&p, part);
    assert_eq!(a, b, "main ion must be deterministic across calls");
    assert_eq!(
        a, y1,
        "tied-zero frequencies must resolve to y1 (Suffix, charge 1), not a hash-order ion"
    );
}

/// The scorer blends with σ equal to the full fragment tolerance (`soft_blend(.., tol_da)`),
/// so a peak at the window edge keeps exp(-1/2), about 61 %, of the matched credit.
const TOL_DA: f32 = 0.5;

#[test]
fn soft_match_centre_gives_full_credit() {
    let matched = 4.0_f32;
    let missing = -2.0_f32;
    let s = soft_blend(matched, missing, 0.0, TOL_DA);
    assert!(
        (s - matched).abs() < 1e-6,
        "centre peak must keep full matched credit, got {s}"
    );
}

#[test]
fn soft_match_edge_keeps_most_of_the_matched_credit() {
    let matched = 4.0_f32;
    let missing = -2.0_f32;
    let s = soft_blend(matched, missing, TOL_DA, TOL_DA);
    assert!(
        s > 0.5 * (matched + missing),
        "with σ = tolerance an edge peak leans toward matched, got {s}"
    );
    let w = (-0.5f32).exp();
    let expected = w * matched + (1.0 - w) * missing;
    assert!(
        (s - expected).abs() < 1e-5,
        "edge blend math mismatch: {s} vs {expected}"
    );
}

#[test]
fn soft_match_monotonic_in_abs_dm() {
    let matched = 5.0_f32;
    let missing = -1.0_f32;
    let mut prev = f32::INFINITY;
    for k in 0..=25 {
        let dm = k as f32 * 0.02;
        let s = soft_blend(matched, missing, dm, TOL_DA);
        assert!(
            s <= prev + 1e-6,
            "blend must be non-increasing in |Δm|: dm={dm} s={s} prev={prev}"
        );
        prev = s;
    }
    let plus = soft_blend(matched, missing, 0.1, TOL_DA);
    let minus = soft_blend(matched, missing, -0.1, TOL_DA);
    assert!(
        (plus - minus).abs() < 1e-6,
        "blend must be symmetric in Δm sign"
    );
}

#[test]
fn soft_match_large_factor_approaches_hard() {
    let matched = 3.0_f32;
    let missing = -2.0_f32;
    let huge_sigma = 1000.0_f32;
    let s = soft_blend(matched, missing, 0.4, huge_sigma);
    assert!(
        (s - matched).abs() < 1e-3,
        "σ ≫ tol must degenerate to hard matching, got {s}"
    );
}

#[test]
fn nearest_peak_rank_and_mz_single_scan_matches_split_helpers() {
    let peaks = vec![(99.95, 10.0_f32), (100.05, 50.0_f32), (100.5, 5.0_f32)];
    let ranks = vec![2_u32, 1, 3];
    let tol = 0.1;
    let combined = nearest_peak_rank_and_mz_in(&peaks, &ranks, 100.0, tol);
    let rank = nearest_peak_rank_in(&peaks, &ranks, 100.0, tol);
    let mz = matched_peak_mz(&peaks, &ranks, 100.0, tol);
    assert_eq!(
        combined.map(|(r, _)| r),
        rank,
        "single-scan rank must match"
    );
    assert_eq!(combined.map(|(_, m)| m), mz, "single-scan m/z must match");
    assert_eq!(combined, Some((1, 100.05)));
    assert_eq!(
        nearest_peak_rank_and_mz_in(&peaks, &ranks, 200.0, tol),
        None
    );
}

use crate::params::{IonType, Partition, SpecDataType};
use crate::psm::rank::RankScorer;
use crate::testutil::tiny_param_with_ions;

fn spec(peaks: &[(f64, f32)]) -> Spectrum {
    Spectrum {
        title: "test".into(),
        precursor_mz: 500.0,
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

#[test]
fn isobaric_windowed_peak_filter_caps_dense_window() {
    use model::protocol::Protocol;
    let peaks: Vec<(f64, f32)> = (0..30)
        .map(|i| (200.0 + i as f64 * 0.9, (100 - i) as f32))
        .collect();
    let s = spec(&peaks);

    let mut p = tiny_param_with_ions();
    p.data_type.protocol = Protocol::Automatic;
    let scorer = RankScorer::new(&p);
    assert_eq!(
        ScoredSpectrum::new(&s, &scorer, 2)
            .dump_active_peaks()
            .len(),
        30,
        "non-isobaric must keep all peaks"
    );

    let mut pt = tiny_param_with_ions();
    pt.data_type.protocol = Protocol::TMT;
    let scorer_t = RankScorer::new(&pt);
    let kept = ScoredSpectrum::new(&s, &scorer_t, 2).dump_active_peaks();
    assert_eq!(kept.len(), 20, "TMT must cap to top-20 per 100-Da window");
    let min_kept = kept
        .iter()
        .map(|&(_, _, it)| it)
        .fold(f32::INFINITY, f32::min);
    assert!(
        min_kept >= 81.0,
        "kept peaks must be the most intense (got min {min_kept})"
    );
}

#[test]
fn loss_node_score_scores_loss_shifted_peak_via_trained_table() {
    let part = Partition {
        charge: 2,
        parent_mass: 1000.0,
        seg_num: 0,
    };
    let loss_ion = IonType::Prefix {
        charge: 1,
        offset_bits: 0.0_f32.to_bits(),
        loss_class: 1,
    };

    let mut p = tiny_param_with_ions();
    p.rank_dist_table
        .get_mut(&part)
        .unwrap()
        .insert(loss_ion, vec![0.6_f32, 0.3, 0.05, 0.001]);

    let nominal = 500.0_f64;
    let loss_mass = 162.0528_f64;
    let theo_mz = loss_ion.mz(nominal) - loss_mass;

    let s = spec(&[(theo_mz, 100.0)]);
    let scorer = RankScorer::new(&p);
    let ss = ScoredSpectrum::new(&s, &scorer, 2);

    let active = [(loss_mass, 1u8)];
    let with_loss = ss.loss_node_score(nominal, 0.0, &scorer, 2, 1000.0, &active, &[]);
    assert!(
        with_loss > 0,
        "loss ion at rank 1 should score positive, got {with_loss}"
    );

    assert_eq!(
        ss.loss_node_score(nominal, 0.0, &scorer, 2, 1000.0, &[], &[]),
        0
    );

    assert_eq!(
        ss.loss_node_score(nominal, 0.0, &scorer, 2, 1000.0, &[(loss_mass, 2u8)], &[]),
        0
    );

    let plain = RankScorer::new(&tiny_param_with_ions());
    let ss2 = ScoredSpectrum::new(&s, &plain, 2);
    assert_eq!(
        ss2.loss_node_score(nominal, 0.0, &plain, 2, 1000.0, &active, &[]),
        0
    );
}

#[test]
fn ion_match_facts_emits_loss_facts_for_loss_peptide() {
    use model::amino_acid::AminoAcid;
    use model::modification::{ModLocation, Modification, ResidueSpec};
    use model::peptide::Peptide;

    let p = tiny_param_with_ions();
    let scorer = RankScorer::new(&p);

    let m = Modification {
        name: "L".into(),
        mass_delta: 0.0,
        residue: ResidueSpec::Specific(b'T'),
        location: ModLocation::Anywhere,
        fixed: false,
        accession: None,
        neutral_losses: vec![162.0528],
        loss_class: 1,
        annotated_site_mass: None,
        blocks_cleavage: false,
        parts: Vec::new(),
        diagnostic_ions: Vec::new(),
    };
    let t = AminoAcid::standard(b'T').unwrap().with_mod(m);
    let w1 = AminoAcid::standard(b'W').unwrap();
    let w2 = AminoAcid::standard(b'W').unwrap();
    let pep = Peptide::new(vec![t, w1, w2], b'_', b'-');

    let intact = IonType::Prefix {
        charge: 1,
        offset_bits: 0.0_f32.to_bits(),
        loss_class: 0,
    };
    let prefix_nominal = nominal_from(pep.residues[0].mass + pep.residues[1].mass);
    let theo = intact.mz(prefix_nominal as f64) - 162.0528;

    let s = Spectrum {
        title: "loss".into(),
        precursor_mz: 600.0,
        precursor_intensity: None,
        precursor_charge: Some(2),
        rt_seconds: None,
        scan: None,
        peaks: vec![(theo, 100.0)],
        activation_method: None,
        isolation_lower_offset: None,
        isolation_upper_offset: None,
        isolation_target_mz: None,
        mobility_window: None,
    };
    let ss = ScoredSpectrum::new(&s, &scorer, 2);

    let facts = ss.ion_match_facts(&pep, &scorer);
    assert!(
        facts
            .iter()
            .any(|f| f.ion_type.loss_class() == 1 && f.rank.is_some()),
        "expected a matched glyco loss fact (rank Some, class 1)"
    );

    let plain = Peptide::new(
        vec![
            AminoAcid::standard(b'W').unwrap(),
            AminoAcid::standard(b'W').unwrap(),
            AminoAcid::standard(b'A').unwrap(),
        ],
        b'_',
        b'-',
    );
    assert!(
        ss.ion_match_facts(&plain, &scorer)
            .iter()
            .all(|f| !f.ion_type.is_loss()),
        "plain peptide must emit no loss facts"
    );
}

#[test]
fn prob_peak_uses_raw_mme_value_not_da_converted() {
    use crate::params::SpecDataType;
    use model::activation::ActivationMethod;
    use model::instrument::InstrumentType;
    use model::protocol::Protocol;
    use model::tolerance::Tolerance;
    use rustc_hash::FxHashMap;

    let precursor_mz = 501.007_276_49_f64;
    let s = Spectrum {
        title: "prob_peak_test".into(),
        precursor_mz,
        precursor_intensity: None,
        precursor_charge: Some(2),
        rt_seconds: None,
        scan: None,
        peaks: vec![(100.0, 1.0), (200.0, 2.0), (300.0, 3.0)],
        activation_method: None,
        isolation_lower_offset: None,
        isolation_upper_offset: None,
        isolation_target_mz: None,
        mobility_window: None,
    };

    let param = Param {
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
        num_precursor_off: 0,
        precursor_off_map: FxHashMap::default(),
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
    };

    let scorer = RankScorer::new(&param);
    let ss = ScoredSpectrum::new(&s, &scorer, 2);

    let parent_mass = (precursor_mz - PROTON) * 2.0;
    let raw_mme = 20.0_f64;
    let approx_num_bins = parent_mass / (raw_mme * 2.0);
    let expected_prob_peak = (3.0_f64 / approx_num_bins.max(1.0)) as f32;

    let wrong_approx_num_bins = parent_mass / (parent_mass * 20e-6 * 2.0);
    let wrong_prob_peak = (3.0_f64 / wrong_approx_num_bins.max(1.0)) as f32;

    assert!(
        (expected_prob_peak - wrong_prob_peak).abs() > 0.001,
        "test precondition failed: Ppm raw vs Da-converted did not produce different prob_peak values"
    );

    assert!(
        (ss.prob_peak - expected_prob_peak).abs() < 1e-5,
        "prob_peak={} but expected={} (raw-mme formula). Wrong Da-converted value would be {}",
        ss.prob_peak,
        expected_prob_peak,
        wrong_prob_peak
    );
}

fn deconv_param(apply: bool) -> Param {
    use model::activation::ActivationMethod;
    use model::instrument::InstrumentType;
    use model::protocol::Protocol;
    use model::tolerance::Tolerance;
    use rustc_hash::FxHashMap;
    Param {
        version: 10001,
        data_type: SpecDataType {
            activation: ActivationMethod::HCD,
            instrument: InstrumentType::QExactive,
            enzyme: None,
            protocol: Protocol::Automatic,
        },
        mme: Tolerance::Ppm(20.0),
        apply_deconvolution: apply,
        deconvolution_error_tolerance: 0.05,
        charge_hist: vec![(2, 100)],
        min_charge: 2,
        max_charge: 4,
        num_segments: 1,
        partitions: vec![],
        num_precursor_off: 0,
        precursor_off_map: FxHashMap::default(),
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
fn deconv_active_for_charge_2_produces_input_equivalent_peaks() {
    let s = Spectrum {
        title: "deconv_test".into(),
        precursor_mz: 501.007_276_49_f64,
        precursor_intensity: None,
        precursor_charge: Some(2),
        rt_seconds: None,
        scan: None,
        peaks: vec![(100.0, 1.0), (200.0, 2.0), (300.0, 3.0)],
        activation_method: None,
        isolation_lower_offset: None,
        isolation_upper_offset: None,
        isolation_target_mz: None,
        mobility_window: None,
    };
    let param = deconv_param(true);
    let scorer = RankScorer::new(&param);
    let ss = ScoredSpectrum::new(&s, &scorer, 2);

    let parent_mass = (s.precursor_mz - PROTON) * 2.0;
    let approx = parent_mass / (20.0_f64 * 2.0);
    let expected = (3.0_f64 / approx.max(1.0)) as f32;
    assert!(
        (ss.prob_peak - expected).abs() < 1e-5,
        "charge=2 deconv-active spectrum: prob_peak={} expected={} (active_count=3)",
        ss.prob_peak,
        expected
    );
}

#[test]
fn deconv_active_for_charge_3_uses_post_deconv_peak_count_for_prob_peak() {
    use model::mass::ISOTOPE;
    let p1 = 100.0;
    let p2 = p1 + ISOTOPE / 2.0;
    let p3 = 500.0;
    let s = Spectrum {
        title: "deconv_charge3".into(),
        precursor_mz: 401.0,
        precursor_intensity: None,
        precursor_charge: Some(3),
        rt_seconds: None,
        scan: None,
        peaks: vec![(p1, 10.0), (p2, 5.0), (p3, 1.0)],
        activation_method: None,
        isolation_lower_offset: None,
        isolation_upper_offset: None,
        isolation_target_mz: None,
        mobility_window: None,
    };
    let param = deconv_param(true);
    let scorer = RankScorer::new(&param);
    let ss = ScoredSpectrum::new(&s, &scorer, 3);

    let active_count = ss.deconv_peaks.as_ref().map(|p| p.len()).unwrap_or(0);
    assert!(
        active_count >= 1,
        "deconv_peaks should be populated for charge=3 + apply_deconvolution=true"
    );
    let parent_mass = (401.0 - PROTON) * 3.0;
    let approx = parent_mass / (20.0_f64 * 2.0);
    let expected = (active_count as f64 / approx.max(1.0)) as f32;
    assert!(
        (ss.prob_peak - expected).abs() < 1e-5,
        "charge=3 deconv-active spectrum: prob_peak={} expected={} (post-deconv count={})",
        ss.prob_peak,
        expected,
        active_count
    );
}

#[test]
fn deconv_off_uses_kept_count_for_prob_peak() {
    let s = Spectrum {
        title: "no_deconv".into(),
        precursor_mz: 501.007_276_49_f64,
        precursor_intensity: None,
        precursor_charge: Some(2),
        rt_seconds: None,
        scan: None,
        peaks: vec![(100.0, 1.0), (200.0, 2.0), (300.0, 3.0), (400.0, 4.0)],
        activation_method: None,
        isolation_lower_offset: None,
        isolation_upper_offset: None,
        isolation_target_mz: None,
        mobility_window: None,
    };
    let param = deconv_param(false);
    let scorer = RankScorer::new(&param);
    let ss = ScoredSpectrum::new(&s, &scorer, 2);

    let parent_mass = (s.precursor_mz - PROTON) * 2.0;
    let approx = parent_mass / (20.0_f64 * 2.0);
    let expected = (4.0_f64 / approx.max(1.0)) as f32;
    assert!(
        (ss.prob_peak - expected).abs() < 1e-5,
        "deconv-off: prob_peak={} expected={} (kept_count=4)",
        ss.prob_peak,
        expected
    );
    assert!(
        ss.deconv_peaks.is_none(),
        "deconv_peaks must be None when apply_deconvolution=false"
    );
}

#[test]
fn observed_node_mass_picks_highest_intensity_peak_in_window() {
    use model::mass::INTEGER_MASS_SCALER;
    let node_nominal = 100_i32;
    let theo_mz = node_nominal as f64 / INTEGER_MASS_SCALER as f64;
    let closer_mz = theo_mz + 0.09;
    let farther_mz = theo_mz + 0.39;
    let s = spec(&[(closer_mz, 1.0), (farther_mz, 100.0)]);
    let param = tiny_param_with_ions();
    let scorer = RankScorer::new(&param);
    let ss = ScoredSpectrum::new_without_filtering(&s);
    let result = ss.observed_node_mass(node_nominal, &scorer, 2, 1000.0);
    let result_mass = result.expect("should find a peak in the window");
    let expected_mass = farther_mz;
    let wrong_mass = closer_mz;
    assert!(
        (result_mass - expected_mass).abs() < 1e-6,
        "expected highest-intensity (farther) peak mass {expected_mass:.6}, \
             got {result_mass:.6} (closest/wrong would be {wrong_mass:.6})"
    );
}

#[test]
fn node_score_does_not_panic_on_empty_spectrum() {
    let s = spec(&[]);
    let param = tiny_param_with_ions();
    let scorer = RankScorer::new(&param);
    let ss = ScoredSpectrum::new_without_filtering(&s);
    let n = ss.node_score(100.0, 900.0, &scorer, 2, 1000.0, 0.5);
    assert!(
        n <= 0,
        "missing-ion score on empty spectrum should be non-positive, got {n}"
    );
}

#[test]
fn node_score_nonzero_when_peak_matches_prefix_ion() {
    use model::mass::INTEGER_MASS_SCALER;
    let nominal = 100.0_f64;
    let b1_mz = nominal / INTEGER_MASS_SCALER as f64;
    let s = spec(&[(50.0, 1.0), (b1_mz, 100.0), (200.0, 2.0)]);
    let param = tiny_param_with_ions();
    let scorer = RankScorer::new(&param);
    let ss = ScoredSpectrum::new_without_filtering(&s);
    let n = ss.node_score(nominal, 900.0, &scorer, 2, 1000.0, 0.5);
    assert!(
        n > 0,
        "expected positive node_score when b-ion peak present, got {n}"
    );
}

#[test]
fn node_score_prefix_only_match() {
    use model::mass::INTEGER_MASS_SCALER;
    let nominal = 57.0_f64;
    let mz = nominal / INTEGER_MASS_SCALER as f64;
    let s = spec(&[(mz, 50.0), (300.0, 1.0)]);
    let param = tiny_param_with_ions();
    let scorer = RankScorer::new(&param);
    let ss = ScoredSpectrum::new_without_filtering(&s);
    let n = ss.node_score(nominal, 900.0, &scorer, 2, 1000.0, 0.5);
    assert!(n > 0, "prefix-only match: expected positive score, got {n}");
}

#[test]
fn node_score_no_matching_ions_returns_negative_or_zero() {
    let s = spec(&[(5000.0, 100.0)]);
    let param = tiny_param_with_ions();
    let scorer = RankScorer::new(&param);
    let ss = ScoredSpectrum::new_without_filtering(&s);
    let n = ss.node_score(100.0, 900.0, &scorer, 2, 1000.0, 0.5);
    assert!(
        n <= 0,
        "missing ion should produce non-positive score, got {n}"
    );
}

#[test]
fn node_score_nominal_mass_zero_prefix_returns_zero() {
    let s = spec(&[]);
    let param = tiny_param_with_ions();
    let scorer = RankScorer::new(&param);
    let ss = ScoredSpectrum::new_without_filtering(&s);
    let n = ss.node_score(0.0, 1000.0, &scorer, 2, 1000.0, 0.5);
    assert!(
        n <= 0,
        "source-node score with empty spectrum should be non-positive, got {n}"
    );
}

#[test]
fn edge_score_returns_zero_when_table_empty() {
    let s = spec(&[(100.0, 1.0)]);
    let mut param = tiny_param_with_ions();
    param.ion_existence_table.clear();
    let scorer = RankScorer::new(&param);
    let ss = ScoredSpectrum::new_without_filtering(&s);
    let e = ss.edge_score(150, 100, 50.0, &scorer, 2, 1000.0);
    assert_eq!(e, 0);
}

#[test]
fn edge_score_returns_zero_when_error_scaling_factor_zero() {
    let s = spec(&[(100.0, 1.0)]);
    let param = tiny_param_with_ions();
    assert_eq!(param.error_scaling_factor, 0);
    let scorer = RankScorer::new(&param);
    let ss = ScoredSpectrum::new_without_filtering(&s);
    let e = ss.edge_score(150, 100, 50.0, &scorer, 2, 1000.0);
    assert_eq!(e, 0);
}

#[test]
fn edge_score_nonzero_with_existence_table() {
    use crate::params::{FragmentOffsetFrequency, SpecDataType};
    use model::activation::ActivationMethod;
    use model::instrument::InstrumentType;
    use model::protocol::Protocol;
    use model::tolerance::Tolerance;
    use rustc_hash::FxHashMap;

    let part = Partition {
        charge: 2,
        parent_mass: 1000.0,
        seg_num: 0,
    };
    let prefix1 = IonType::Prefix {
        charge: 1,
        offset_bits: 0.0_f32.to_bits(),
        loss_class: 0,
    };
    let noise = IonType::Noise;

    let ion_freqs = vec![0.6_f32, 0.3, 0.05, 0.001];
    let noise_freqs = vec![0.1_f32, 0.2, 0.3, 0.4];

    let mut ion_table: FxHashMap<IonType, Vec<f32>> = FxHashMap::default();
    ion_table.insert(prefix1, ion_freqs);
    ion_table.insert(noise, noise_freqs);

    let mut rank_dist_table: FxHashMap<Partition, FxHashMap<IonType, Vec<f32>>> =
        FxHashMap::default();
    rank_dist_table.insert(part, ion_table);

    let mut frag_off_table = FxHashMap::default();
    frag_off_table.insert(
        part,
        vec![FragmentOffsetFrequency {
            ion_type: prefix1,
            frequency: 0.7,
        }],
    );

    let error_scaling_factor = 2_i32;
    let dist_len = (error_scaling_factor as usize) * 2 + 1;

    let mut ion_err_dist_table: FxHashMap<Partition, Vec<f32>> = FxHashMap::default();
    ion_err_dist_table.insert(part, vec![0.1_f32, 0.2, 0.4, 0.2, 0.1]);

    let mut noise_err_dist_table: FxHashMap<Partition, Vec<f32>> = FxHashMap::default();
    noise_err_dist_table.insert(part, vec![0.05_f32, 0.1, 0.7, 0.1, 0.05]);

    let mut ion_existence_table: FxHashMap<Partition, Vec<f32>> = FxHashMap::default();
    ion_existence_table.insert(part, vec![0.1_f32, 0.3, 0.3, 0.5]);

    let _ = dist_len;

    let mut param = Param {
        version: 10001,
        data_type: SpecDataType {
            activation: ActivationMethod::HCD,
            instrument: InstrumentType::QExactive,
            enzyme: None,
            protocol: Protocol::Automatic,
        },
        mme: Tolerance::Da(0.5),
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
        error_scaling_factor,
        ion_err_dist_table,
        noise_err_dist_table,
        ion_existence_table,
        partition_ion_types_cache: FxHashMap::default(),
        gbdt_peak_model: None,
        frag_intensity_model: None,
        rich_ion_model: None,
    };
    param.rebuild_cache();

    let s = spec(&[]);
    let scorer = RankScorer::new(&param);
    let ss = ScoredSpectrum::new_without_filtering(&s);
    let e = ss.edge_score(150, 100, 50.0, &scorer, 2, 1000.0);
    assert_ne!(
        e, 0,
        "edge_score should be nonzero with populated existence table"
    );
}

#[test]
fn empty_spectrum_yields_no_ranks() {
    let s = spec(&[]);
    let ss = ScoredSpectrum::new_without_filtering(&s);
    assert_eq!(ss.peak_count(), 0);
    assert!(ss.nearest_peak_rank(500.0, 0.1).is_none());
}

#[test]
fn highest_intensity_gets_rank_1() {
    let s = spec(&[(100.0, 1.0), (200.0, 5.0), (300.0, 3.0)]);
    let ss = ScoredSpectrum::new_without_filtering(&s);
    assert_eq!(ss.peak_count(), 3);
    assert_eq!(ss.nearest_peak_rank(200.0, 0.1), Some(1));
    assert_eq!(ss.nearest_peak_rank(300.0, 0.1), Some(2));
    assert_eq!(ss.nearest_peak_rank(100.0, 0.1), Some(3));
}

#[test]
fn nearest_peak_within_tolerance() {
    let s = spec(&[(100.0, 1.0), (200.5, 5.0), (300.0, 3.0)]);
    let ss = ScoredSpectrum::new_without_filtering(&s);
    assert_eq!(ss.nearest_peak_rank(200.4, 0.2), Some(1));
    assert_eq!(ss.nearest_peak_rank(200.5, 0.001), Some(1));
    assert_eq!(ss.nearest_peak_rank(200.4, 0.05), None);
}

#[test]
fn ties_broken_deterministically() {
    let s = spec(&[(100.0, 5.0), (200.0, 5.0)]);
    let ss = ScoredSpectrum::new_without_filtering(&s);
    let r1 = ss.nearest_peak_rank(100.0, 0.1).unwrap();
    let r2 = ss.nearest_peak_rank(200.0, 0.1).unwrap();
    assert_ne!(r1, r2);
    assert!(r1 == 1 || r2 == 1);
    assert!(r1 == 2 || r2 == 2);
}

#[test]
fn nearest_peak_full_picks_max_intensity_within_tolerance() {
    let s = spec(&[(99.6, 5.0), (100.0, 9.0), (100.5, 1.0)]);
    let ss = ScoredSpectrum::new_without_filtering(&s);

    let (rank, intensity, mz) = ss.nearest_peak_full(100.45, 0.6).unwrap();
    assert_eq!(
        rank, 1,
        "max-intensity peak (100.0) must win, not nearest-mz (100.5)"
    );
    assert_eq!(mz, 100.0);
    assert_eq!(intensity, 9.0);

    assert_eq!(ss.nearest_peak_rank(100.45, 0.6), Some(1));
}

#[test]
fn nearest_peak_rank_matches_linear_scan_on_many_peaks() {
    let mut peaks: Vec<(f64, f32)> = (0..100)
        .map(|i| (i as f64 * 10.0 + 0.5, (100 - i) as f32))
        .collect();
    peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let s = Spectrum {
        title: "many".into(),
        precursor_mz: 500.0,
        precursor_intensity: None,
        precursor_charge: Some(2),
        rt_seconds: None,
        scan: None,
        peaks,
        activation_method: None,
        isolation_lower_offset: None,
        isolation_upper_offset: None,
        isolation_target_mz: None,
        mobility_window: None,
    };
    let ss = ScoredSpectrum::new_without_filtering(&s);

    for target in [50.5, 100.0, 250.0, 333.7, 500.5, 750.5, 999.5] {
        let tol = 5.0_f64;
        let bs_result = ss.nearest_peak_rank(target, tol);
        let bf_result = {
            let mut best: Option<(usize, f64)> = None;
            for (i, &(mz, _)) in s.peaks.iter().enumerate() {
                if (mz - target).abs() <= tol
                    && best.as_ref().is_none_or(|(_, d)| (mz - target).abs() < *d)
                {
                    best = Some((i, (mz - target).abs()));
                }
            }
            best.map(|(i, _)| ss.peak_rank_at(i).unwrap_or(u32::MAX))
        };
        assert_eq!(
            bs_result, bf_result,
            "binary search and linear scan differ at target {target}"
        );
    }
}

#[test]
fn gbdt_term_raises_node_score_for_signal_peak() {
    use crate::gbdt::{GbdtPeakModel, Tree};
    use model::mass::INTEGER_MASS_SCALER;

    let gbdt = GbdtPeakModel {
        packed: Default::default(),
        n_features: 1,
        feature_names: crate::gbdt::indexed_feature_names(1),
        apply_sigmoid: true,
        trees: vec![Tree {
            feature: vec![-1],
            threshold: vec![0.0],
            left: vec![-1],
            right: vec![-1],
            value: vec![4.0],
            default_left: vec![1],
        }],
        iso_x: vec![],
        iso_y: vec![],
    };

    const PROTON_F64: f64 = model::mass::PROTON;
    let nominal = 100.0_f64;
    let b1_mz = nominal / INTEGER_MASS_SCALER as f64;
    let precursor_mz = (1000.0 + 2.0 * PROTON_F64) / 2.0;

    let peaks_vec = vec![(b1_mz, 1000.0_f32)];
    let spec_plain = Spectrum {
        title: "gbdt_test".into(),
        precursor_mz,
        precursor_intensity: None,
        precursor_charge: Some(2),
        rt_seconds: None,
        scan: None,
        peaks: peaks_vec.clone(),
        activation_method: None,
        isolation_lower_offset: None,
        isolation_upper_offset: None,
        isolation_target_mz: None,
        mobility_window: None,
    };

    let param_plain = tiny_param_with_ions();
    let scorer_plain = RankScorer::new(&param_plain);
    let s_plain = ScoredSpectrum::new(&spec_plain, &scorer_plain, 2);
    let parent_mass_plain = s_plain.parent_mass;
    let base =
        s_plain.directional_node_score(nominal, true, &scorer_plain, 2, parent_mass_plain, 0.5);

    let mut param_gbdt = tiny_param_with_ions();
    param_gbdt.gbdt_peak_model = Some(gbdt);
    let scorer_gbdt = RankScorer::new(&param_gbdt);
    let spec_gbdt = Spectrum {
        title: "gbdt_test".into(),
        precursor_mz,
        precursor_intensity: None,
        precursor_charge: Some(2),
        rt_seconds: None,
        scan: None,
        peaks: peaks_vec,
        activation_method: None,
        isolation_lower_offset: None,
        isolation_upper_offset: None,
        isolation_target_mz: None,
        mobility_window: None,
    };
    let s_gbdt = ScoredSpectrum::new(&spec_gbdt, &scorer_gbdt, 2);
    let parent_mass_gbdt = s_gbdt.parent_mass;
    let with = s_gbdt.directional_node_score(nominal, true, &scorer_gbdt, 2, parent_mass_gbdt, 0.5);

    assert!(
        with > base + 3.0,
        "expected GBDT bump ≥ 3.0: base={base}, with={with}"
    );
}
