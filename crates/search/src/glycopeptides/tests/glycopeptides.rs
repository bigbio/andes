#[test]
fn glycan_subtraction_uses_the_precursor_tolerance() {
    use model::tolerance::{PrecursorTolerance, Tolerance};
    let aa = model::AminoAcidSetBuilder::new_standard().build().unwrap();
    let mut params = crate::params::SearchParams::default_tryptic(aa);
    params.precursor_tolerance = PrecursorTolerance::symmetric(Tolerance::Ppm(10.0));
    assert!((super::precursor_window_da(&params, 3000.0) - 0.03).abs() < 1e-12);
    params.precursor_tolerance =
        PrecursorTolerance::asymmetric(Tolerance::Ppm(5.0), Tolerance::Da(0.05));
    assert!((super::precursor_window_da(&params, 3000.0) - 0.05).abs() < 1e-12);
}

#[test]
fn o_glycan_floor_admits_a_light_tn_precursor() {
    use super::{min_glycan_floor, GlycanCore};
    use glyco::glycan::mass::HEXNAC;
    assert!(HEXNAC - 0.01 >= min_glycan_floor(GlycanCore::OGlycan));
    assert!(HEXNAC < min_glycan_floor(GlycanCore::NGlycan));
}

#[test]
fn o_glycan_hex_gate_drops_hex_without_evidence() {
    use super::{hex_gate_passes, GlycanComp, GlycanCore};
    let hex_glycan = GlycanComp::from_counts(1, 1, 0, 0, 0, 0);
    let no_hex = GlycanComp::from_counts(1, 0, 0, 0, 0, 0);
    assert!(hex_gate_passes(
        GlycanCore::NGlycan,
        Some(&hex_glycan),
        false
    ));
    assert!(hex_gate_passes(GlycanCore::OGlycan, Some(&no_hex), false));
    assert!(hex_gate_passes(GlycanCore::OGlycan, None, false));
    assert!(!hex_gate_passes(
        GlycanCore::OGlycan,
        Some(&hex_glycan),
        false
    ));
    assert!(hex_gate_passes(
        GlycanCore::OGlycan,
        Some(&hex_glycan),
        true
    ));
}

#[test]
fn peptide_first_backbones_are_exempt_from_the_core_y_prefilter() {
    use super::{dedup_backbone_hits, has_backbone_evidence, BackboneHit, Source};
    let mass_only = BackboneHit {
        backbone_mass: 1500.0,
        glycan: None,
        source: Source::Db,
        charge: 3,
        isotope_offset: 0,
        glycan_mass_residual: 1800.0,
        peptide_first: false,
    };
    let peptide_first = BackboneHit {
        peptide_first: true,
        ..mass_only.clone()
    };
    assert!(!has_backbone_evidence(&mass_only, 0));
    assert!(has_backbone_evidence(&mass_only, 1));
    assert!(has_backbone_evidence(&peptide_first, 0));
    let merged = dedup_backbone_hits(vec![mass_only.clone(), peptide_first.clone()], 20.0);
    assert_eq!(merged.len(), 1);
    assert!(merged[0].peptide_first);
    let merged = dedup_backbone_hits(vec![peptide_first, mass_only], 20.0);
    assert_eq!(merged.len(), 1);
    assert!(merged[0].peptide_first);
}

#[test]
fn isotope_window_for_applies_the_override_per_spectrum() {
    use super::isotope_window_for;
    let default = 0..=2;
    let narrow = 0..=1;
    assert_eq!(isotope_window_for(0, &default, None, None), 0..=2);
    assert_eq!(isotope_window_for(0, &default, None, Some(&[true])), 0..=2);
    assert_eq!(isotope_window_for(7, &default, Some(&narrow), None), 0..=1);
    let mask = [true, false, true];
    assert_eq!(
        isotope_window_for(0, &default, Some(&narrow), Some(&mask)),
        0..=1
    );
    assert_eq!(
        isotope_window_for(1, &default, Some(&narrow), Some(&mask)),
        0..=2
    );
    assert_eq!(
        isotope_window_for(2, &default, Some(&narrow), Some(&mask)),
        0..=1
    );
    assert_eq!(
        isotope_window_for(3, &default, Some(&narrow), Some(&mask)),
        0..=2
    );
}

use super::*;

#[test]
fn build_hcd_partners_pairs_etd_to_precursor_matched_hcd() {
    use model::activation::ActivationMethod;
    let sp = |mz: f64, am: ActivationMethod| Spectrum {
        title: String::new(),
        precursor_mz: mz,
        precursor_intensity: None,
        precursor_charge: Some(3),
        rt_seconds: None,
        scan: None,
        peaks: vec![],
        activation_method: Some(am),
        isolation_lower_offset: None,
        isolation_upper_offset: None,
        isolation_target_mz: None,
        mobility_window: None,
    };
    use ActivationMethod::{ETD, HCD};
    let spectra = vec![
        sp(500.0, HCD),
        sp(500.0, ETD),
        sp(600.0, HCD),
        sp(600.005, ETD),
        sp(999.0, ETD),
    ];
    let p = build_hcd_partners(&spectra);
    assert_eq!(p[0], None, "HCD spectra never get a partner");
    assert_eq!(p[1], Some(0), "ETD@500 pairs to HCD@500");
    assert_eq!(p[2], None);
    assert_eq!(p[3], Some(2), "ETD@600.005 pairs to HCD@600 within m/z tol");
    assert_eq!(p[4], None, "ETD@999 has no precursor-matched HCD partner");
}

#[test]
fn order_peptide_first_is_deterministic_across_input_orders() {
    let permutations = vec![
        vec![(8u32, 9u32), (5, 7), (2, 5), (3, 9), (1, 7)],
        vec![(1, 7), (3, 9), (2, 5), (8, 9), (5, 7)],
        vec![(2, 5), (1, 7), (5, 7), (8, 9), (3, 9)],
        vec![(5, 7), (8, 9), (3, 9), (1, 7), (2, 5)],
    ];

    let mut reference = permutations[0].clone();
    order_peptide_first(&mut reference);

    for perm in &permutations {
        let mut got = perm.clone();
        order_peptide_first(&mut got);
        assert_eq!(
            &got, &reference,
            "input order {perm:?} must collapse to the same deterministic order"
        );
    }

    for w in reference.windows(2) {
        assert!(
            w[0].1 >= w[1].1,
            "counts must be non-increasing: {reference:?}"
        );
    }
    assert_eq!(reference[0].1, 9, "strongest count first");
    assert_eq!(reference.last().unwrap().1, 5, "weakest count last");
}

#[test]
fn isotope_sweep_annotates_glycan_only_at_correct_offset() {
    use glyco::glycan::db::GlycanComp;
    use glyco::glycan::mass::{CORE_Y_STEPS, PROTON as GLY_PROTON};
    use glyco::hybrid::hybrid_candidates_with_isotope;

    let true_backbone_residue = 1500.0_f64;
    let glycan = GlycanComp::from_counts(2, 3, 0, 0, 0, 0);
    let true_precursor_neutral = true_backbone_residue + glycan.mass;

    let observed_neutral = true_precursor_neutral + ISOTOPE;

    let glycans = vec![glycan];

    let y0_neutral = true_backbone_residue + H2O;
    let mut peaks: Vec<(f64, f32)> = vec![(204.08665, 200.0), (138.05496, 120.0)];
    peaks.push((y0_neutral + GLY_PROTON, 150.0));
    for &s in CORE_Y_STEPS.iter() {
        peaks.push((y0_neutral + s + GLY_PROTON, 100.0));
    }

    let hits0 =
        hybrid_candidates_with_isotope(&peaks, observed_neutral, 2, 0, &glycans, 20.0, 5);
    let matching0: Vec<_> = hits0
        .iter()
        .filter(|h| (h.backbone_mass - true_backbone_residue).abs() < 0.05)
        .collect();
    assert!(
        !matching0.is_empty(),
        "backbone is read from the ladder → recovered even at the wrong isotope offset"
    );
    assert!(
        matching0.iter().any(|h| h.source == Source::DeNovo),
        "at the wrong offset the glycan is ~1 ISOTOPE off and must NOT annotate to a known composition"
    );
    assert!(
        matching0.iter().all(|h| h.source != Source::Db),
        "wrong isotope offset must not produce a DB annotation for the true backbone"
    );

    let precursor_neutral_iso1 = observed_neutral - ISOTOPE;
    let hits1 =
        hybrid_candidates_with_isotope(&peaks, precursor_neutral_iso1, 2, 1, &glycans, 20.0, 5);
    let hit = hits1
        .iter()
        .find(|h| {
            (h.backbone_mass - true_backbone_residue).abs() < 0.05 && h.source == Source::Db
        })
        .expect("offset +1 must recover AND annotate the backbone via the corrected precursor");
    assert_eq!(
        hit.isotope_offset, 1,
        "recovered hit must record isotope_offset=1"
    );
    assert_eq!(
        hit.charge, 2,
        "recovered hit must record the charge it was matched at"
    );
}

#[test]
fn charges_to_try_is_the_reported_charge_or_the_range() {
    let range = 2u8..=3u8;
    assert_eq!(glyco_charges_to_try(Some(5), &range), vec![5]);
    assert_eq!(glyco_charges_to_try(Some(2), &range), vec![2]);
    assert_eq!(glyco_charges_to_try(None, &range), vec![2, 3]);
    assert_eq!(glyco_charges_to_try(Some(0), &range), vec![2, 3]);
    assert_eq!(glyco_charges_to_try(Some(-1), &range), vec![2, 3]);
}

#[test]
fn dedup_preserves_distinct_isotope_residuals_for_novel_glycans() {
    use glyco::hybrid::BackboneHit;
    let mk = |iso: i8, residual: f64| BackboneHit {
        backbone_mass: 1500.0,
        glycan: None,
        source: Source::DeNovo,
        charge: 3,
        isotope_offset: iso,
        glycan_mass_residual: residual,
        peptide_first: false,
    };
    let out = dedup_backbone_hits(vec![mk(0, 892.317), mk(1, 891.313)], 20.0);
    assert_eq!(
        out.len(),
        2,
        "distinct isotope residuals must not be merged: {out:?}"
    );
    let residuals: Vec<f64> = out.iter().map(|h| h.glycan_mass_residual).collect();
    assert!(residuals.iter().any(|r| (r - 892.317).abs() < 1e-6));
    assert!(residuals.iter().any(|r| (r - 891.313).abs() < 1e-6));

    let dup = dedup_backbone_hits(vec![mk(0, 892.317), mk(0, 892.317)], 20.0);
    assert_eq!(dup.len(), 1, "same-hypothesis duplicates must merge");
}

#[test]
fn full_glyco_psm_is_clone() {
    let psm = PsmMatch {
        spectrum_idx: 0,
        candidate_idxs: smallvec::smallvec![0],
        charge_used: 2,
        mass_error_ppm: 0.0,
        score: 0.0,
        rank_score: 0.0,
        edge_score: 0,
        activation_method: None,
        features: PsmFeatures::default(),
        isotope_offset: 0,
        precursor_mz_override: None,
        site_evidence: None,
        tie_key: 0,
    };
    let key = GlycoPsmKey {
        spectrum_idx: 0,
        glycan: None,
        glycan_source: Source::Db,
        oxonium_summed_frac: 0.0,
        n_core_oxonium_ions: 0,
        y_ladder_intensity_score: 0.0,
        y_hit_frac: 0.0,
        partial_glycan_by: 0.0,
        y0y1_anchor_score: 0.0,
        sialic_consistency: 0.0,
        core_y_hits: 0,
        glycan_mass: 0.0,
        backbone_mass: 0.0,
        cz_hyperscore: 0.0,
        cz_intensity: 0.0,
        cz_explained: 0.0,
        cz_chance_llr: 0.0,
    };
    let hit = FullGlycoPsm {
        glycan_key: key,
        psm,
    };
    let cloned = hit.clone();
    assert_eq!(cloned.psm.spectrum_idx, 0);
}

#[test]
fn glyco_spectrum_result_is_clone() {
    let result = GlycoSpectrumResult {
        spectrum_idx: 7,
        hits: vec![],
    };
    let c = result.clone();
    assert_eq!(c.spectrum_idx, 7);
}

#[test]
fn core_y_ranking_promotes_supported_backbone_over_unsupported() {
    use glyco::backbone::{count_core_y_hits, SpectrumStats};
    use glyco::glycan::mass::{CORE_Y_STEPS, PROTON};

    let true_bb = 1100.0_f64;

    let noise_bb = 2800.0_f64;

    let mut peaks: Vec<(f64, f32)> = vec![
        (true_bb + PROTON, 500.0),
        (true_bb + PROTON + CORE_Y_STEPS[0], 400.0),
        (true_bb + PROTON + CORE_Y_STEPS[1], 350.0),
        (true_bb + PROTON + CORE_Y_STEPS[2], 300.0),
        (true_bb + PROTON + CORE_Y_STEPS[3], 250.0),
        (true_bb + PROTON + CORE_Y_STEPS[4], 200.0),
        (900.0, 10.0),
        (1050.0, 10.0),
    ];
    peaks.push((noise_bb + PROTON + 5.0, 50.0));

    let tol_ppm = 20.0;

    let stats = SpectrumStats::new(&peaks);
    let true_hits = count_core_y_hits(&peaks, &stats, true_bb, tol_ppm, 3, GlycanCore::NGlycan);
    let noise_hits =
        count_core_y_hits(&peaks, &stats, noise_bb, tol_ppm, 3, GlycanCore::NGlycan);

    assert_eq!(
        true_hits, 6,
        "expected all 6 core-Y hits for true_bb, got {true_hits}",
    );
    assert_eq!(
        noise_hits, 0,
        "expected 0 core-Y hits for noise_bb, got {noise_hits}",
    );

    let mut candidates = [
        (noise_bb, noise_hits),
        (true_bb, true_hits),
    ];
    candidates.sort_by(|&(am, ay), &(bm, by)| {
        by.cmp(&ay)
            .then_with(|| bm.partial_cmp(&am).unwrap_or(std::cmp::Ordering::Equal))
    });

    assert!(
        (candidates[0].0 - true_bb).abs() < 0.01,
        "expected true_bb ranked first after core-Y sort, got backbone_mass={}",
        candidates[0].0
    );
    assert!(
        (candidates[1].0 - noise_bb).abs() < 0.01,
        "expected noise_bb ranked second, got backbone_mass={}",
        candidates[1].0
    );
}

#[test]
fn by_rank_promotes_byone_matching_backbone_over_y_ladder_backbone() {
    let true_bb_mass = 1500.0_f64;
    let true_bb_best_rank: f32 = 10.0;
    let true_bb_core_y: u8 = 0;

    let noise_bb_mass = 2000.0_f64;
    let noise_bb_best_rank: f32 = 2.0;
    let noise_bb_core_y: u8 = 6;

    let backbones = [
        (noise_bb_mass, noise_bb_best_rank, noise_bb_core_y),
        (true_bb_mass, true_bb_best_rank, true_bb_core_y),
    ];
    let mut order: Vec<usize> = (0..backbones.len()).collect();
    order.sort_by(|&ai, &bi| {
        backbones[bi]
            .1
            .partial_cmp(&backbones[ai].1)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| backbones[bi].2.cmp(&backbones[ai].2))
            .then_with(|| {
                backbones[bi]
                    .0
                    .partial_cmp(&backbones[ai].0)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    });

    assert_eq!(
        order[0], 1,
        "expected true_bb (idx=1) ranked first by b/y rank, got idx={}",
        order[0]
    );
    assert_eq!(
        order[1], 0,
        "expected noise_bb (idx=0) ranked second, got idx={}",
        order[1]
    );
}
