use super::*;
use crate::index::SearchIndex;
use model::isotope::averagine_isotope_envelope;

fn envelope(mono_mz: f64, charge: u8, neutral: f64, scale: f32) -> Vec<(f64, f32)> {
    let env = averagine_isotope_envelope(neutral, 4);
    (0..4)
        .map(|k| {
            (
                mono_mz + k as f64 * ISOTOPE / charge as f64,
                (env[k] as f32) * scale,
            )
        })
        .collect()
}

#[test]
fn detects_coisolated_excludes_selected() {
    let z = 2u8;
    let selected_mz = 600.0;
    let sel_neutral = (selected_mz - PROTON) * z as f64;
    let co_mz = 600.7;
    let co_neutral = (co_mz - PROTON) * z as f64;
    let mut peaks = envelope(selected_mz, z, sel_neutral, 1000.0);
    peaks.extend(envelope(co_mz, z, co_neutral, 500.0));
    peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    let got = detect_coisolated(&peaks, 599.0, 601.5, selected_mz, 2..=3, 0.02, 0.5, 2);
    assert_eq!(got.len(), 1, "exactly one co-isolated (selected excluded)");
    assert!((got[0].mono_mz - co_mz).abs() < 0.02);
    assert_eq!(got[0].charge, z);
}

#[test]
fn heavy_coisolated_envelope_is_reported_at_its_monoisotope() {
    let z = 2u8;
    let selected_mz = 1200.0;
    let sel_neutral = (selected_mz - PROTON) * z as f64;
    let co_mz = 1501.2;
    let co_neutral = (co_mz - PROTON) * z as f64;
    let env = averagine_isotope_envelope(co_neutral, 4);
    assert!(env[1] > env[0], "fixture needs M+1 > M");
    let mut peaks = envelope(selected_mz, z, sel_neutral, 1000.0);
    peaks.extend(envelope(co_mz, z, co_neutral, 500.0));
    peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    let got = detect_coisolated(&peaks, 1499.0, 1504.0, selected_mz, 2..=3, 0.01, 0.5, 4);
    assert_eq!(got.len(), 1, "one co-isolated envelope: {got:?}");
    assert!(
        (got[0].mono_mz - co_mz).abs() < 0.01,
        "reported at {} instead of the monoisotope {co_mz}",
        got[0].mono_mz
    );
    assert_eq!(got[0].charge, z);
    assert!((got[0].neutral_mass - co_neutral).abs() < 0.02);
}

#[test]
fn max_n_zero_is_a_hard_disable() {
    let z = 2u8;
    let selected_mz = 600.0;
    let sel_neutral = (selected_mz - PROTON) * z as f64;
    let co_mz = 600.7;
    let co_neutral = (co_mz - PROTON) * z as f64;
    let mut peaks = envelope(selected_mz, z, sel_neutral, 1000.0);
    peaks.extend(envelope(co_mz, z, co_neutral, 500.0));
    peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    let got = detect_coisolated(&peaks, 599.0, 601.5, selected_mz, 2..=3, 0.02, 0.5, 0);
    assert!(
        got.is_empty(),
        "max_n=0 must emit zero secondaries, got {}",
        got.len()
    );
    let one = detect_coisolated(&peaks, 599.0, 601.5, selected_mz, 2..=3, 0.02, 0.5, 1);
    assert_eq!(one.len(), 1, "max_n=1 caps at one");
}

#[test]
fn no_coisolation_when_only_selected_present() {
    let z = 2u8;
    let selected_mz = 600.0;
    let peaks = envelope(selected_mz, z, (selected_mz - PROTON) * z as f64, 1000.0);
    let got = detect_coisolated(&peaks, 599.0, 601.5, selected_mz, 2..=3, 0.02, 0.5, 2);
    assert!(
        got.is_empty(),
        "only the selected precursor -> no co-isolation"
    );
}

use crate::PreparedSearch;
use model::activation::ActivationMethod;
use model::instrument::InstrumentType;
use model::protocol::Protocol;
use model::Tolerance;
use model::{AminoAcid, AminoAcidSetBuilder, Protein, ProteinDb};
use rustc_hash::FxHashMap;
use scoring_crate::params::{IonType, Partition, SpecDataType};
use scoring_crate::psm::ions::predict_by_ions;
use scoring_crate::Param;

pub(crate) fn tiny_scorer() -> RankScorer {
    let part = Partition {
        charge: 2,
        parent_mass: 500.0,
        seg_num: 0,
    };
    let prefix1 = IonType::Prefix {
        charge: 1,
        offset_bits: 0.0_f32.to_bits(),
        loss_class: 0,
    };
    let suffix1 = IonType::Suffix {
        charge: 1,
        offset_bits: 0.0_f32.to_bits(),
        loss_class: 0,
    };
    let noise = IonType::Noise;

    let mut ion_table = FxHashMap::default();
    ion_table.insert(prefix1, vec![0.5_f32, 0.1, 0.05, 0.01]);
    ion_table.insert(suffix1, vec![0.5_f32, 0.1, 0.05, 0.01]);
    ion_table.insert(noise, vec![0.05_f32, 0.05, 0.05, 0.05]);

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
        gbdt_peak_model: None,
        frag_intensity_model: None,
        rich_ion_model: None,
    };
    param.rebuild_cache();
    RankScorer::new(&param)
}

fn ion_scorer() -> RankScorer {
    use scoring_crate::params::FragmentOffsetFrequency;
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
    let mut ion_table = FxHashMap::default();
    ion_table.insert(prefix1, vec![0.6_f32, 0.3, 0.05, 0.001]);
    ion_table.insert(suffix1, vec![0.6_f32, 0.3, 0.05, 0.001]);
    ion_table.insert(IonType::Noise, vec![0.1_f32, 0.2, 0.3, 0.4]);
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
    let mut param = Param {
        version: 10001,
        data_type: SpecDataType {
            activation: ActivationMethod::HCD,
            instrument: InstrumentType::QExactive,
            enzyme: None,
            protocol: Protocol::Automatic,
        },
        mme: Tolerance::Da(0.02),
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

fn whole_db_pool(prepared: &PreparedSearch<'_>) -> (Vec<Candidate>, BTreeMap<i32, Vec<usize>>) {
    let pool = prepared.candidates_where(|_, _| true);
    let mut bucket_index: BTreeMap<i32, Vec<usize>> = BTreeMap::new();
    for (i, c) in pool.iter().enumerate() {
        bucket_index
            .entry(c.peptide.nominal_residue_mass())
            .or_default()
            .push(i);
    }
    (pool, bucket_index)
}

fn peptide(residues: &[u8], pre: u8, post: u8) -> Peptide {
    let aas: Vec<AminoAcid> = residues
        .iter()
        .map(|&r| AminoAcid::standard(r).unwrap())
        .collect();
    Peptide::new(aas, pre, post)
}

#[test]
fn secondary_search_finds_planted_peptide() {
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "".into(),
            sequence: b"MKWVTFISLLRKDAFLGSFLYEYSRK".to_vec(),
            annotation: None,
        }],
    };
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
    let params = SearchParams::default_tryptic(aa_set);
    let scorer = tiny_scorer();
    let frag_tol = 0.02_f64;

    let prepared =
        PreparedSearch::prepare_in_memory(&idx, &params, &scorer, frag_tol, "XXX").unwrap();
    let (candidates, bucket_index) = whole_db_pool(&prepared);

    let planted = peptide(b"DAFLGSFLYEYSR", b'K', b'K');
    let z = 2u8;
    let co_neutral = planted.mass();
    let co_mz = (co_neutral + z as f64 * PROTON) / z as f64;
    let co = CoIsolated {
        mono_mz: co_mz,
        charge: z,
        neutral_mass: co_neutral,
    };

    let peaks: Vec<(f64, f32)> = predict_by_ions(&planted, 1..=1)
        .into_iter()
        .map(|p| (p.mz, 100.0_f32))
        .collect();
    let mut peaks = peaks;
    peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    let primary = peptide(b"WVTFISLLR", b'K', b'K');
    let prim_mz = (primary.mass() + z as f64 * PROTON) / z as f64;
    let spec = Spectrum {
        title: "chimeric-secondary".into(),
        precursor_mz: prim_mz,
        precursor_intensity: None,
        precursor_charge: Some(z as i32),
        rt_seconds: None,
        scan: None,
        peaks,
        activation_method: None,
        isolation_lower_offset: None,
        isolation_upper_offset: None,
        isolation_target_mz: None,
        mobility_window: None,
    };

    let got = search_secondary(
        &spec,
        &primary,
        &std::collections::HashSet::new(),
        co,
        &candidates,
        &bucket_index,
        &scorer,
        &prepared.aa_set_for_scoring,
        &params,
        frag_tol,
    );

    let (psm, winner_claimed) =
        got.expect("secondary search should return a PSM at the co-isolated mass");
    let found = &candidates[psm.primary_candidate_idx() as usize].peptide;
    assert_eq!(
        found.residues, planted.residues,
        "secondary PSM should resolve to the planted peptide DAFLGSFLYEYSR"
    );

    assert!(
        !winner_claimed.is_empty(),
        "winner should claim its matched peaks"
    );
    let again = search_secondary(
        &spec,
        &primary,
        &winner_claimed,
        co,
        &candidates,
        &bucket_index,
        &scorer,
        &prepared.aa_set_for_scoring,
        &params,
        frag_tol,
    );
    let matched_after = again
        .as_ref()
        .map(|(p, _)| p.features.num_matched_main_ions)
        .unwrap_or(0);
    assert!(
        matched_after < psm.features.num_matched_main_ions,
        "with its peaks already claimed, the secondary must match fewer ions \
         (was {}, now {})",
        psm.features.num_matched_main_ions,
        matched_after
    );
    assert!(
        psm.rank_score > 0.0,
        "secondary PSM rank_score (RawScore) must be positive, got {}",
        psm.rank_score
    );
}

fn planted_secondary_spectrum(planted: &Peptide, primary: &Peptide) -> (Spectrum, CoIsolated) {
    let z = 2u8;
    let mut peaks: Vec<(f64, f32)> = predict_by_ions(planted, 1..=1)
        .into_iter()
        .map(|p| (p.mz, 100.0_f32))
        .collect();
    peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let spec = Spectrum {
        title: "chimeric-planted".into(),
        precursor_mz: (primary.mass() + z as f64 * PROTON) / z as f64,
        precursor_intensity: None,
        precursor_charge: Some(z as i32),
        rt_seconds: None,
        scan: None,
        peaks,
        activation_method: None,
        isolation_lower_offset: None,
        isolation_upper_offset: None,
        isolation_target_mz: None,
        mobility_window: None,
    };
    let co = CoIsolated {
        mono_mz: (planted.mass() + z as f64 * PROTON) / z as f64,
        charge: z,
        neutral_mass: planted.mass(),
    };
    (spec, co)
}

#[test]
fn secondary_listwise_gap_ignores_duplicate_occurrences_of_the_winner() {
    let target = ProteinDb {
        proteins: vec![
            Protein {
                accession: "P1".into(),
                description: "".into(),
                sequence: b"MKWVTFISLLRKDAFLGSFLYEYSRK".to_vec(),
                annotation: None,
            },
            Protein {
                accession: "P2".into(),
                description: "".into(),
                sequence: b"GGGKDAFLGSFLYEYSRKGGG".to_vec(),
                annotation: None,
            },
            Protein {
                accession: "P3".into(),
                description: "".into(),
                sequence: b"GGGKADFGLSFLYEYSRKGGG".to_vec(),
                annotation: None,
            },
        ],
    };
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
    let params = SearchParams::default_tryptic(aa_set);
    let scorer = ion_scorer();
    let frag_tol = 0.02_f64;
    let prepared =
        PreparedSearch::prepare_in_memory(&idx, &params, &scorer, frag_tol, "XXX").unwrap();
    let (candidates, bucket_index) = whole_db_pool(&prepared);

    let planted = peptide(b"DAFLGSFLYEYSR", b'K', b'K');
    let occurrences = candidates
        .iter()
        .filter(|c| !c.is_decoy && c.peptide.residues == planted.residues)
        .count();
    assert!(
        occurrences >= 2,
        "fixture must hold the winner twice, got {occurrences}"
    );

    let primary = peptide(b"WVTFISLLR", b'K', b'K');
    let (spec, co) = planted_secondary_spectrum(&planted, &primary);
    let (psm, _) = search_secondary(
        &spec,
        &primary,
        &std::collections::HashSet::new(),
        co,
        &candidates,
        &bucket_index,
        &scorer,
        &prepared.aa_set_for_scoring,
        &params,
        frag_tol,
    )
    .expect("planted secondary must be found");
    let winner = &candidates[psm.primary_candidate_idx() as usize];
    assert!(!winner.is_decoy);
    assert_eq!(winner.peptide.residues, planted.residues);
    assert!(
        psm.features.listwise_score_gap > 0.0,
        "gap must be measured against a different peptidoform, got {}",
        psm.features.listwise_score_gap
    );
}

#[test]
fn secondary_features_share_the_primary_scales() {
    let target = ProteinDb {
        proteins: vec![
            Protein {
                accession: "P1".into(),
                description: "".into(),
                sequence: b"MKWVTFISLLRKDAFLGSFLYEYSRK".to_vec(),
                annotation: None,
            },
            Protein {
                accession: "P2".into(),
                description: "".into(),
                sequence: b"GGGKADFGLSFLYEYSRKGGG".to_vec(),
                annotation: None,
            },
        ],
    };
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
    let scorer = ion_scorer();
    let frag_tol = 0.02_f64;
    let planted = peptide(b"DAFLGSFLYEYSR", b'K', b'K');
    let primary = peptide(b"WVTFISLLR", b'K', b'K');
    let (spec, co) = planted_secondary_spectrum(&planted, &primary);

    let run = |mode: crate::params::ScoreMode| {
        let mut params = SearchParams::default_tryptic(aa_set.clone());
        params.score_mode = mode;
        let prepared =
            PreparedSearch::prepare_in_memory(&idx, &params, &scorer, frag_tol, "XXX").unwrap();
        let (candidates, bucket_index) = whole_db_pool(&prepared);
        search_secondary(
            &spec,
            &primary,
            &std::collections::HashSet::new(),
            co,
            &candidates,
            &bucket_index,
            &scorer,
            &prepared.aa_set_for_scoring,
            &params,
            frag_tol,
        )
        .expect("planted secondary must be found")
        .0
    };
    let rank = run(crate::params::ScoreMode::Rank);
    let f = &rank.features;
    assert_ne!(f.rank_score_float, 0.0);
    assert!(
        (f.rank_score_float - rank.score).abs() < 0.5 * f.num_matched_main_ions as f32 + 1.0,
        "RankScoreFloat {} vs RawScore {}",
        f.rank_score_float,
        rank.score
    );
    assert!(f.listwise_score_gap > 0.0);
    assert!(f.strong_score_cal.is_finite());
    assert_ne!(f.strong_score_cal, f.strong_score);

    let strong = run(crate::params::ScoreMode::Strong);
    assert_eq!(strong.score, rank.score);
    assert_eq!(strong.rank_score, strong.features.strong_score);
}

#[test]
fn secondary_search_applies_precursor_calibration_shift() {
    let shift_ppm = 30.0_f64;
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "".into(),
            sequence: b"MKWVTFISLLRKDAFLGSFLYEYSRK".to_vec(),
            annotation: None,
        }],
    };
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
    let mut params = SearchParams::default_tryptic(aa_set);
    params.precursor_mass_shift_ppm = shift_ppm;
    let scorer = tiny_scorer();
    let frag_tol = 0.02_f64;
    let prepared =
        PreparedSearch::prepare_in_memory(&idx, &params, &scorer, frag_tol, "XXX").unwrap();
    let (candidates, bucket_index) = whole_db_pool(&prepared);

    let planted = peptide(b"DAFLGSFLYEYSR", b'K', b'K');
    let z = 2u8;
    let true_mass = planted.mass();
    let raw_neutral = true_mass / (1.0 - shift_ppm * 1e-6);
    let co_mz = (raw_neutral + z as f64 * PROTON) / z as f64;
    let co = CoIsolated {
        mono_mz: co_mz,
        charge: z,
        neutral_mass: raw_neutral,
    };

    let mut peaks: Vec<(f64, f32)> = predict_by_ions(&planted, 1..=1)
        .into_iter()
        .map(|p| (p.mz, 100.0_f32))
        .collect();
    peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let primary = peptide(b"WVTFISLLR", b'K', b'K');
    let prim_mz = (primary.mass() + z as f64 * PROTON) / z as f64;
    let spec = Spectrum {
        title: "chimeric-cal".into(),
        precursor_mz: prim_mz,
        precursor_intensity: None,
        precursor_charge: Some(z as i32),
        rt_seconds: None,
        scan: None,
        peaks,
        activation_method: None,
        isolation_lower_offset: None,
        isolation_upper_offset: None,
        isolation_target_mz: None,
        mobility_window: None,
    };

    let got = search_secondary(
        &spec,
        &primary,
        &std::collections::HashSet::new(),
        co,
        &candidates,
        &bucket_index,
        &scorer,
        &prepared.aa_set_for_scoring,
        &params,
        frag_tol,
    );
    let (psm, _claimed) =
        got.expect("secondary must be found at the calibration-adjusted mass");
    let found = &candidates[psm.primary_candidate_idx() as usize].peptide;
    assert_eq!(
        found.residues, planted.residues,
        "should resolve to the planted peptide"
    );
    assert!(
        psm.mass_error_ppm.abs() < 5.0,
        "mass error must be reported on the calibrated scale (~0), got {} ppm",
        psm.mass_error_ppm
    );
}

#[test]
fn secondary_window_uses_left_for_lighter_and_right_for_heavier() {
    use model::tolerance::{PrecursorTolerance, Tolerance};
    let aa_set = model::AminoAcidSetBuilder::new_standard().build().unwrap();
    let mut params = SearchParams::default_tryptic(aa_set);
    params.precursor_tolerance =
        PrecursorTolerance::asymmetric(Tolerance::Da(0.05), Tolerance::Da(0.5));
    let co = CoIsolated {
        mono_mz: 600.0,
        charge: 2,
        neutral_mass: 1198.0,
    };
    let (center, left, right) = secondary_mass_window(co, &params);
    assert_eq!((left, right), (0.05, 0.5));
    assert!(in_secondary_mass_window(center + 0.3, center, left, right));
    assert!(!in_secondary_mass_window(center - 0.3, center, left, right));
    assert!(in_secondary_mass_window(center - 0.04, center, left, right));
    assert!(!in_secondary_mass_window(center + 0.6, center, left, right));
    let (lo, hi) = secondary_nominal_window(co, &params);
    assert!(lo <= nominal_from(center - left - H2O));
    assert!(hi >= nominal_from(center + right - H2O));
}
