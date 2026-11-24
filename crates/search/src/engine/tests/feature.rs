use super::*;
use model::activation::ActivationMethod;
use model::amino_acid::AminoAcid;
use model::instrument::InstrumentType;
use model::mass::PROTON;
use model::peptide::Peptide;
use model::protocol::Protocol;
use model::spectrum::Spectrum;
use model::tolerance::Tolerance;
use rustc_hash::FxHashMap;
use scoring_crate::params::{FragmentOffsetFrequency, IonType, Partition, SpecDataType};
use scoring_crate::psm::ions::predict_by_ions;
use scoring_crate::psm::ScoredSpectrum;

fn make_scorer(tol_da: f64) -> RankScorer {
    use model::mass::{H2O, PROTON};
    let part = Partition {
        charge: 2,
        parent_mass: 0.0,
        seg_num: 0,
    };
    let prefix1 = IonType::Prefix {
        charge: 1,
        offset_bits: (PROTON as f32).to_bits(),
        loss_class: 0,
    };
    let suffix1 = IonType::Suffix {
        charge: 1,
        offset_bits: ((H2O + PROTON) as f32).to_bits(),
        loss_class: 0,
    };
    let noise = IonType::Noise;
    let mut ion_table = FxHashMap::default();
    ion_table.insert(prefix1, vec![0.6_f32, 0.3, 0.05, 0.001]);
    ion_table.insert(suffix1, vec![0.6_f32, 0.3, 0.05, 0.001]);
    ion_table.insert(noise, vec![0.1_f32, 0.2, 0.3, 0.4]);
    let mut rank_dist_table = FxHashMap::default();
    rank_dist_table.insert(part, ion_table);
    let mut frag_off_table = FxHashMap::default();
    frag_off_table.insert(
        part,
        vec![
            FragmentOffsetFrequency {
                ion_type: prefix1,
                frequency: 0.7,
            },
            FragmentOffsetFrequency {
                ion_type: suffix1,
                frequency: 0.7,
            },
        ],
    );
    let mut param = scoring_crate::Param {
        version: 10001,
        data_type: SpecDataType {
            activation: ActivationMethod::HCD,
            instrument: InstrumentType::QExactive,
            enzyme: None,
            protocol: Protocol::Automatic,
        },
        mme: Tolerance::Da(tol_da),
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
        gbdt_peak_model: None,
        frag_intensity_model: None,
        rich_ion_model: None,
    };
    param.rebuild_cache();
    RankScorer::new(&param)
}

fn ala_peptide(len: usize) -> Peptide {
    let aa = AminoAcid::standard(b'A').unwrap();
    Peptide::new(vec![aa; len], b'_', b'-')
}

fn make_spectrum(peaks: Vec<(f64, f32)>) -> Spectrum {
    Spectrum {
        title: "test".into(),
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
    }
}

#[test]
fn recentred_spectra_removes_the_measured_offset() {
    let spec = make_spectrum(vec![(300.0, 1.0), (1000.0, 2.0)]);
    assert!(recentred_spectra(std::slice::from_ref(&spec), 0.0).is_none());
    let theo = [300.0, 1000.0];
    let low = make_spectrum(theo.iter().map(|m| (m * (1.0 - 8e-6), 1.0)).collect());
    let c = recentred_spectra(std::slice::from_ref(&low), -8.0).unwrap();
    for (p, t) in c[0].peaks.iter().zip(theo) {
        assert!((p.0 - t).abs() / t < 1e-9, "{} vs {t}", p.0);
    }
    assert_eq!(low.peaks[0].0, 300.0 * (1.0 - 8e-6));
    assert!(c[0].peaks.windows(2).all(|w| w[0].0 < w[1].0));
}

#[test]
fn compute_psm_features_top7_error_stats_zero_when_no_matches() {
    let pep = ala_peptide(4);
    let spec = make_spectrum(vec![]);
    let ss = ScoredSpectrum::new_without_filtering(&spec);
    let f = compute_psm_features(&ss, &pep, &make_scorer(0.5), 2);
    assert_eq!(
        f.mean_error_top7, 0.0,
        "mean_error_top7 should be 0 with no matches"
    );
    assert_eq!(
        f.stdev_error_top7, 0.0,
        "stdev_error_top7 should be 0 with no matches"
    );
    assert_eq!(
        f.mean_rel_error_top7, 0.0,
        "mean_rel_error_top7 should be 0 with no matches"
    );
    assert_eq!(
        f.stdev_rel_error_top7, 0.0,
        "stdev_rel_error_top7 should be 0 with no matches"
    );
    assert_eq!(
        f.explained_ion_current_ratio, 0.0,
        "ratio should be 0 with no peaks"
    );
    assert_eq!(
        f.ms2_ion_current, 0.0,
        "ms2_ion_current should be 0 with no peaks"
    );
}

#[test]
fn compute_psm_features_populates_ion_current_ratios() {
    let pep = ala_peptide(3);
    let predicted = predict_by_ions(&pep, 1..=1);

    let mut peaks: Vec<(f64, f32)> = predicted
        .iter()
        .enumerate()
        .map(|(i, p)| (p.mz, (i + 1) as f32 * 10.0))
        .collect();
    peaks.push((1500.0, 5.0));
    peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    let spec = make_spectrum(peaks);
    let ss = ScoredSpectrum::new_without_filtering(&spec);
    let f = compute_psm_features(&ss, &pep, &make_scorer(0.01), 2);

    assert!(
        f.explained_ion_current_ratio > 0.0,
        "explained_ion_current_ratio should be > 0 when ions match, got {}",
        f.explained_ion_current_ratio
    );
    assert!(
        f.n_term_ion_current_ratio > 0.0,
        "n_term_ion_current_ratio should be > 0 when b-ions match"
    );
    assert!(
        f.c_term_ion_current_ratio > 0.0,
        "c_term_ion_current_ratio should be > 0 when y-ions match"
    );

    let sum = f.n_term_ion_current_ratio + f.c_term_ion_current_ratio;
    assert!(
        (sum - f.explained_ion_current_ratio).abs() < 1e-5,
        "n_term + c_term should == explained ({} + {} != {})",
        f.n_term_ion_current_ratio,
        f.c_term_ion_current_ratio,
        f.explained_ion_current_ratio
    );

    let total: f32 = ss.total_intensity() as f32;
    assert!(
        (f.ms2_ion_current - total).abs() < 1.0,
        "ms2_ion_current {} should match total spectrum intensity {}",
        f.ms2_ion_current,
        total
    );

    assert_eq!(f.isolation_window_efficiency, 0.0);
}

#[test]
fn compute_psm_features_unique_match_fraction() {
    let scorer = make_scorer(0.01);

    let pep = ala_peptide(5);
    let predicted = predict_by_ions(&pep, 1..=1);
    let mut peaks: Vec<(f64, f32)> = predicted
        .iter()
        .enumerate()
        .map(|(i, p)| (p.mz, (i + 1) as f32 * 10.0))
        .collect();
    peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let f_unique = compute_psm_features(
        &ScoredSpectrum::new_without_filtering(&make_spectrum(peaks)),
        &pep,
        &scorer,
        2,
    );
    assert!(
        f_unique.unique_match_fraction > 0.7,
        "separated ions should yield a high fraction, got {}",
        f_unique.unique_match_fraction
    );

    const H2O: f64 = model::mass::H2O;
    const NH3: f64 = model::mass::NH3;
    let mut contested_case: Option<(Peptide, f64, f64)> = None;
    'search: for seq in ["PFAH", "PFSR", "PFLM", "PFML"] {
        let residues: Vec<AminoAcid> = seq
            .bytes()
            .map(|r| AminoAcid::standard(r).unwrap())
            .collect();
        let pep2 = Peptide::new(residues, b'_', b'-');
        let pred2 = predict_by_ions(&pep2, 1..=1);
        let mut theo: Vec<f64> = Vec::new();
        for p in &pred2 {
            theo.push(p.mz);
            theo.push(p.mz - H2O);
            theo.push(p.mz - NH3);
        }
        for p in &pred2 {
            let tol = p.mz * 20e-6;
            let near_real = pred2.iter().filter(|q| (q.mz - p.mz).abs() < 0.05).count();
            let others: Vec<f64> = theo
                .iter()
                .copied()
                .filter(|t| (t - p.mz).abs() > 1e-6 && (t - p.mz).abs() <= tol)
                .collect();
            if near_real == 1 && others.len() == 1 {
                contested_case = Some((pep2, p.mz, others[0]));
                break 'search;
            }
        }
    }
    let (pep2, peak_mz, other_mz) =
        contested_case.expect("should find a contested theo pair among the candidate peptides");
    let tol = peak_mz * 20e-6;
    assert!(
        (other_mz - peak_mz).abs() <= tol,
        "contested pair should be within feature tolerance"
    );
    let contested = vec![(peak_mz, 100.0_f32)];
    let f_contested = compute_psm_features(
        &ScoredSpectrum::new_without_filtering(&make_spectrum(contested)),
        &pep2,
        &scorer,
        2,
    );
    assert!(
        f_contested.unique_match_fraction < f_unique.unique_match_fraction,
        "contested peak should lower fraction (unique={}, contested={})",
        f_unique.unique_match_fraction,
        f_contested.unique_match_fraction
    );
    assert!(
        (f_contested.unique_match_fraction - 0.5).abs() < 0.15,
        "two ions on one peak → fraction ~0.5, got {}",
        f_contested.unique_match_fraction
    );
}

#[test]
fn compute_psm_features_doubly_charged_matched_ion_count() {
    let pep = ala_peptide(4);
    let predicted_z2 = predict_by_ions(&pep, 2..=2);
    assert!(
        !predicted_z2.is_empty(),
        "charge-2 ions should exist for 4-mer"
    );

    let f_empty = compute_psm_features(
        &ScoredSpectrum::new_without_filtering(&make_spectrum(vec![])),
        &pep,
        &make_scorer(0.01),
        3,
    );
    assert_eq!(f_empty.doubly_charged_matched_ion_count, 0);

    let mut peaks: Vec<(f64, f32)> = predicted_z2
        .iter()
        .enumerate()
        .map(|(i, p)| (p.mz, (i + 1) as f32 * 10.0))
        .collect();
    peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let f = compute_psm_features(
        &ScoredSpectrum::new_without_filtering(&make_spectrum(peaks)),
        &pep,
        &make_scorer(0.01),
        3,
    );
    assert!(
        f.doubly_charged_matched_ion_count > 0,
        "charge-2 peaks should match, got {}",
        f.doubly_charged_matched_ion_count
    );
    assert_eq!(
        f.doubly_charged_matched_ion_count,
        predicted_z2.len() as u32,
        "all charge-2 ions placed → count should equal {}",
        predicted_z2.len()
    );
}

#[test]
fn compute_psm_features_mean_matched_intensity_rank() {
    let pep = ala_peptide(3);
    let predicted = predict_by_ions(&pep, 1..=1);

    let mut top: Vec<(f64, f32)> = predicted
        .iter()
        .enumerate()
        .map(|(i, p)| (p.mz, 1000.0 - i as f32))
        .collect();
    top.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let f_top = compute_psm_features(
        &ScoredSpectrum::new_without_filtering(&make_spectrum(top)),
        &pep,
        &make_scorer(0.01),
        2,
    );
    assert!(
        (f_top.mean_matched_intensity_rank - 2.5).abs() < 0.1,
        "ions as the only peaks → mean rank should be ~2.5, got {}",
        f_top.mean_matched_intensity_rank
    );

    let mut buried: Vec<(f64, f32)> = predicted
        .iter()
        .enumerate()
        .map(|(i, p)| (p.mz, 10.0 + i as f32))
        .collect();
    for j in 0..10 {
        buried.push((1500.0 + j as f64 * 10.0, 500.0 + j as f32 * 50.0));
    }
    buried.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let f_buried = compute_psm_features(
        &ScoredSpectrum::new_without_filtering(&make_spectrum(buried)),
        &pep,
        &make_scorer(0.01),
        2,
    );
    assert!(
        f_buried.mean_matched_intensity_rank > f_top.mean_matched_intensity_rank,
        "buried ions should have higher mean rank (top={}, buried={})",
        f_top.mean_matched_intensity_rank,
        f_buried.mean_matched_intensity_rank
    );
}

#[test]
fn compute_psm_features_longest_complementary_ladder() {
    let pep = ala_peptide(5);
    let predicted = predict_by_ions(&pep, 1..=1);
    let n = pep.length();

    let mut all: Vec<(f64, f32)> = predicted
        .iter()
        .enumerate()
        .map(|(i, p)| (p.mz, (i + 1) as f32 * 10.0))
        .collect();
    all.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let f_all = compute_psm_features(
        &ScoredSpectrum::new_without_filtering(&make_spectrum(all)),
        &pep,
        &make_scorer(0.01),
        2,
    );
    let expected = (n - 1) as u32;
    assert_eq!(
        f_all.longest_complementary_ladder, expected,
        "all ions matched → ladder should equal n-1={expected}, got {}",
        f_all.longest_complementary_ladder
    );

    let mut split: Vec<(f64, f32)> = predicted
        .iter()
        .enumerate()
        .filter(|(_, p)| !(matches!(p.kind, IonKind::Y) && p.position == 3))
        .map(|(i, p)| (p.mz, (i + 1) as f32 * 10.0))
        .collect();
    split.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let f_split = compute_psm_features(
        &ScoredSpectrum::new_without_filtering(&make_spectrum(split)),
        &pep,
        &make_scorer(0.01),
        2,
    );
    assert_eq!(
        f_split.longest_complementary_ladder, 2,
        "missing middle bond should leave a length-2 ladder, got {}",
        f_split.longest_complementary_ladder
    );
    assert!(
        f_split.longest_complementary_ladder < f_all.longest_complementary_ladder,
        "split ladder ({}) should be shorter than full ({})",
        f_split.longest_complementary_ladder,
        f_all.longest_complementary_ladder
    );
}

#[test]
fn compute_psm_features_complementary_ion_balance() {
    let pep = ala_peptide(5);
    let predicted = predict_by_ions(&pep, 1..=1);

    let mut all: Vec<(f64, f32)> = predicted
        .iter()
        .enumerate()
        .map(|(i, p)| (p.mz, (i + 1) as f32 * 10.0))
        .collect();
    all.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let f_all = compute_psm_features(
        &ScoredSpectrum::new_without_filtering(&make_spectrum(all)),
        &pep,
        &make_scorer(0.01),
        2,
    );
    assert!(
        f_all.complementary_ion_balance > 0.0,
        "all ions matched → balance should be > 0, got {}",
        f_all.complementary_ion_balance
    );

    let mut b_only: Vec<(f64, f32)> = predicted
        .iter()
        .enumerate()
        .filter(|(_, p)| matches!(p.kind, IonKind::B))
        .map(|(i, p)| (p.mz, (i + 1) as f32 * 10.0))
        .collect();
    b_only.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let f_b = compute_psm_features(
        &ScoredSpectrum::new_without_filtering(&make_spectrum(b_only)),
        &pep,
        &make_scorer(0.01),
        2,
    );
    assert_eq!(
        f_b.complementary_ion_balance, 0.0,
        "b-only spectrum → no complementary pair → balance 0, got {}",
        f_b.complementary_ion_balance
    );
}

#[test]
fn compute_psm_features_neutral_loss_ion_count() {
    let pep = ala_peptide(3);
    let predicted = predict_by_ions(&pep, 1..=1);
    let b2 = predicted
        .iter()
        .find(|p| matches!(p.kind, IonKind::B) && p.position == 2)
        .expect("b2 should exist");

    let mut peaks = vec![(b2.mz, 100.0_f32), (b2.mz - 18.0106, 50.0_f32)];
    peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let f = compute_psm_features(
        &ScoredSpectrum::new_without_filtering(&make_spectrum(peaks)),
        &pep,
        &make_scorer(0.01),
        2,
    );
    assert_eq!(
        f.neutral_loss_ion_count, 1,
        "b2 + b2−H2O partner should yield neutral_loss_ion_count=1, got {}",
        f.neutral_loss_ion_count
    );
}

#[test]
fn compute_psm_features_strong_score_stage1_bolt_ons() {
    let pep = ala_peptide(3);
    let predicted = predict_by_ions(&pep, 1..=1);

    let mut exact: Vec<(f64, f32)> = predicted
        .iter()
        .enumerate()
        .map(|(i, p)| (p.mz, (i + 1) as f32 * 10.0))
        .collect();
    exact.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let f_exact = compute_psm_features(
        &ScoredSpectrum::new_without_filtering(&make_spectrum(exact)),
        &pep,
        &make_scorer(0.01),
        2,
    );

    assert!(
        (f_exact.ppm_gaussian_score - 4.0).abs() < 0.05,
        "ppm_gaussian_score should be ~4 for 4 exact matches, got {}",
        f_exact.ppm_gaussian_score
    );
    assert_eq!(
        f_exact.longest_complementary_ladder, 2,
        "both cleavage sites should form a length-2 ladder, got {}",
        f_exact.longest_complementary_ladder
    );

    let mut shifted: Vec<(f64, f32)> = predicted
        .iter()
        .enumerate()
        .map(|(i, p)| (p.mz * (1.0 + 15e-6), (i + 1) as f32 * 10.0))
        .collect();
    shifted.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let f_shift = compute_psm_features(
        &ScoredSpectrum::new_without_filtering(&make_spectrum(shifted)),
        &pep,
        &make_scorer(0.01),
        2,
    );
    assert!(
        f_shift.ppm_gaussian_score < f_exact.ppm_gaussian_score * 0.5,
        "ppm_gaussian_score should fall sharply when matches scatter (exact={}, shifted={})",
        f_exact.ppm_gaussian_score,
        f_shift.ppm_gaussian_score
    );
}

#[test]
fn compute_psm_features_error_stats_nonzero_when_ions_match_with_offset() {
    let pep = ala_peptide(5);
    let predicted = predict_by_ions(&pep, 1..=1);

    let offset_da = 0.0005_f64;
    let mut peaks: Vec<(f64, f32)> = predicted
        .iter()
        .enumerate()
        .map(|(i, p)| (p.mz + offset_da, (i + 1) as f32 * 10.0))
        .collect();
    peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    let spec = make_spectrum(peaks);
    let ss = ScoredSpectrum::new_without_filtering(&spec);
    let f = compute_psm_features(&ss, &pep, &make_scorer(0.05), 2);

    assert!(
        f.mean_error_top7 > 0.0,
        "mean_error_top7 should be > 0 when peaks are systematically offset, got {}",
        f.mean_error_top7
    );
    assert!(
        f.stdev_error_top7 < 20.0,
        "stdev_error_top7 should be small (single-digit ppm) for identical-Da offset, got {}",
        f.stdev_error_top7
    );
    assert!(
        f.mean_rel_error_top7 != 0.0,
        "mean_rel_error_top7 should be nonzero when peaks are offset"
    );
}

#[test]
fn ms2_ion_current_equals_total_intensity() {
    let pep = ala_peptide(3);
    let peaks = vec![(100.0, 50.0_f32), (200.0, 30.0), (300.0, 20.0)];
    let spec = make_spectrum(peaks.clone());
    let ss = ScoredSpectrum::new_without_filtering(&spec);
    let f = compute_psm_features(&ss, &pep, &make_scorer(0.5), 2);

    let expected: f32 = peaks.iter().map(|&(_, i)| i).sum();
    assert_eq!(
        f.ms2_ion_current, expected,
        "ms2_ion_current {} should equal sum of peak intensities {}",
        f.ms2_ion_current, expected
    );
}

#[test]
fn b1_mz_for_alanine_is_proton_plus_residue_mass() {
    use model::amino_acid::AminoAcid;
    let aa = AminoAcid::standard(b'A').unwrap();
    let residue_mass = aa.mass;
    let expected_b1_mz = residue_mass + PROTON;
    let pep = ala_peptide(2);
    let predicted = predict_by_ions(&pep, 1..=1);
    let b1 = predicted
        .iter()
        .find(|p| matches!(p.kind, IonKind::B) && p.position == 1)
        .expect("b1 ion should exist");
    assert!(
        (b1.mz - expected_b1_mz).abs() < 1e-6,
        "b1 mz {} expected {}",
        b1.mz,
        expected_b1_mz
    );
}
