//! PSM features written to the PIN: enzyme-cleavage credit, fragment-prediction features and
//! the per-PSM feature vector of a kept match. Part of the match engine.

use super::*;

#[inline(always)]
#[allow(
    clippy::too_many_arguments,
    reason = "hot-loop driver; all args are orthogonal cleavage parameters"
)]
pub(crate) fn enzyme_cleavage_credit(
    cand: &Candidate,
    enz: Enzyme,
    enz_is_c_term: bool,
    enz_is_n_term: bool,
    credit_neighboring: i32,
    penalty_neighboring: i32,
    credit_peptide: i32,
    penalty_peptide: i32,
) -> i32 {
    let mut score: i32 = 0;
    let pre = cand.peptide.pre;
    let post = cand.peptide.post;
    if enz_is_c_term {
        score += if cand.is_protein_n_term || enz.is_cleavable(pre) {
            credit_neighboring
        } else {
            penalty_neighboring
        };
        let last = match cand.peptide.residues.last() {
            Some(aa) => aa.residue,
            None => 0,
        };
        score += if enz.is_cleavable(last) {
            credit_peptide
        } else {
            penalty_peptide
        };
    } else if enz_is_n_term {
        let first = match cand.peptide.residues.first() {
            Some(aa) => aa.residue,
            None => 0,
        };
        score += if enz.is_cleavable(first) {
            credit_peptide
        } else {
            penalty_peptide
        };
        score += if cand.is_protein_c_term || enz.is_cleavable(post) {
            credit_neighboring
        } else {
            penalty_neighboring
        };
    }
    score
}

pub(crate) fn cleavage_credit_for(cand: &Candidate, enz: Enzyme, aa_set: &AminoAcidSet) -> i32 {
    enzyme_cleavage_credit(
        cand,
        enz,
        enz.is_c_term(),
        enz.is_n_term(),
        aa_set.neighboring_aa_cleavage_credit(),
        aa_set.neighboring_aa_cleavage_penalty(),
        aa_set.peptide_cleavage_credit(),
        aa_set.peptide_cleavage_penalty(),
    )
}

pub(crate) fn compute_psm_features(
    scored_spec: &ScoredSpectrum<'_>,
    peptide: &Peptide,
    scorer: &RankScorer,
    charge: u8,
) -> PsmFeatures {
    compute_psm_features_with_edge(scored_spec, peptide, scorer, charge, None, true, None)
}

/// [`predict_frag_intensities`], through `cache` when one is given.
pub(super) fn predict_frag_with(
    model: &scoring_crate::gbdt_eval::GbdtPeakModel,
    peptide: &Peptide,
    charge: u8,
    cache: Option<&crate::fragments::predictions::FragPredCache>,
) -> Vec<f64> {
    match cache {
        Some(c) => c.intensities(model, peptide, charge),
        None => predict_frag_intensities(model, peptide, charge),
    }
}

/// The PSM features that only a full feature pass computes: the fragment-intensity LLR
/// battery, the rich-ion LLR and the modification-site features.
#[derive(Default)]
pub(super) struct KeptPsmFeatureValues {
    frag_pred_explained: f32,
    frag_pred_chance_llr: f32,
    frag_topk_observed: f32,
    rich_ion_llr: f32,
    mod_site: [f32; scoring_crate::mod_site_features::N_MOD_SITE_FEATURES],
}

/// Computes the full-pass-only features of `peptide` at `charge`. `frag_pred` is the
/// fragment-intensity prediction of the scorer's model for this peptide and charge (`None`
/// without a model).
pub(super) fn kept_psm_feature_values(
    scored_spec: &ScoredSpectrum<'_>,
    peptide: &Peptide,
    scorer: &RankScorer,
    charge: u8,
    frag_pred: Option<&[f64]>,
) -> KeptPsmFeatureValues {
    let model_tol = scorer.feature_match_tolerance();
    let model_tol_is_ppm = matches!(model_tol, model::tolerance::Tolerance::Ppm(_));
    let model_tol = model_tol.raw_value();
    let (frag_pred_explained, frag_pred_chance_llr, frag_topk_observed) = frag_llr_battery(
        scorer.param().frag_intensity_model.as_deref(),
        scored_spec,
        peptide,
        charge,
        model_tol,
        model_tol_is_ppm,
        frag_pred,
    );
    let rich_ion_llr = rich_ion_llr(
        scorer.param().rich_ion_model.as_deref(),
        scored_spec,
        peptide,
        charge,
        model_tol,
        model_tol_is_ppm,
    );
    let mod_site = mod_site_features(peptide, scored_spec, model_tol, model_tol_is_ppm);
    KeptPsmFeatureValues {
        frag_pred_explained,
        frag_pred_chance_llr,
        frag_topk_observed,
        rich_ion_llr,
        mod_site,
    }
}

/// Completes `features`, computed by a pass with `full = false` for `peptide` at `charge`,
/// with the full-pass-only features, giving the same values as a full pass. `frag_cache`
/// shares fragment-intensity predictions within a chunk.
pub(crate) fn complete_psm_features(
    features: &mut PsmFeatures,
    scored_spec: &ScoredSpectrum<'_>,
    peptide: &Peptide,
    scorer: &RankScorer,
    charge: u8,
    frag_cache: Option<&crate::fragments::predictions::FragPredCache>,
) {
    if peptide.length() < 2 {
        return;
    }
    let frag_pred: Option<Vec<f64>> = scorer
        .param()
        .frag_intensity_model
        .as_deref()
        .map(|g| predict_frag_with(g, peptide, charge, frag_cache));
    let kept = kept_psm_feature_values(scored_spec, peptide, scorer, charge, frag_pred.as_deref());
    features.frag_pred_explained = kept.frag_pred_explained;
    features.frag_pred_chance_llr = kept.frag_pred_chance_llr;
    features.frag_topk_observed = kept.frag_topk_observed;
    features.rich_ion_llr = kept.rich_ion_llr;
    features.mod_site_shifted_matched = kept.mod_site[FEAT_SITE_SHIFTED_MATCHED];
    features.mod_site_shifted_frac = kept.mod_site[FEAT_SITE_SHIFTED_FRAC];
    features.mod_site_intens_frac = kept.mod_site[FEAT_SITE_INTENS_FRAC];
    features.mod_site_localized = kept.mod_site[FEAT_SITE_LOCALIZED];
    features.mod_site_det_count = kept.mod_site[FEAT_SITE_DET_COUNT];
}

pub(crate) fn compute_psm_features_with_edge(
    scored_spec: &ScoredSpectrum<'_>,
    peptide: &Peptide,
    scorer: &RankScorer,
    charge: u8,
    known_edge_score: Option<i32>,
    full: bool,
    frag_cache: Option<&crate::fragments::predictions::FragPredCache>,
) -> PsmFeatures {
    let n = peptide.length();
    if n < 2 {
        return PsmFeatures::default();
    }

    let edge_score =
        known_edge_score.unwrap_or_else(|| psm_edge_score(scored_spec, peptide, scorer, charge));

    let predict_losses = scorer
        .param()
        .data_type
        .activation
        .predicts_neutral_losses();
    let predicted = predict_by_ions_with_losses(peptide, 1..=1, predict_losses);
    let mut b_matched: SmallVec<[bool; 64]> = smallvec![false; n - 1];
    let mut y_matched: SmallVec<[bool; 64]> = smallvec![false; n - 1];
    let mut b_rank: SmallVec<[u32; 64]> = smallvec![u32::MAX; n - 1];
    let mut y_rank: SmallVec<[u32; 64]> = smallvec![u32::MAX; n - 1];

    let mut matched_ions: SmallVec<[MatchedIon; 96]> = SmallVec::new();
    let mut matched_ion_ranks: SmallVec<[u32; 96]> = SmallVec::new();

    let feat_tol = scorer.pin_match_tolerance();
    let feature_tol_is_ppm = matches!(feat_tol, model::tolerance::Tolerance::Ppm(_));
    let feature_tol = feat_tol.raw_value();

    let model_tol = scorer.feature_match_tolerance();
    let model_tol_is_ppm = matches!(model_tol, model::tolerance::Tolerance::Ppm(_));
    let model_tol = model_tol.raw_value();

    const H2O_LOSS_DA: f64 = model::mass::H2O;
    const NH3_LOSS_DA: f64 = model::mass::NH3;
    let mut neutral_loss_ion_count: u32 = 0;

    for p in &predicted {
        let tol_da = if feature_tol_is_ppm {
            p.mz * feature_tol / 1e6
        } else {
            feature_tol
        };
        if let Some((rank, intensity, peak_mz)) = scored_spec.nearest_peak_full(p.mz, tol_da) {
            let is_b = matches!(p.kind, IonKind::B | IonKind::C);
            matched_ions.push((intensity, peak_mz, p.mz, is_b, FragmentIon::of_predicted(p)));
            matched_ion_ranks.push(rank);

            let pos = (p.position - 1) as usize;
            match p.kind {
                IonKind::B | IonKind::C => {
                    if pos < b_matched.len() {
                        b_matched[pos] = true;
                        if !p.is_loss() && rank < b_rank[pos] {
                            b_rank[pos] = rank;
                        }
                    }
                }
                IonKind::Y | IonKind::Z => {
                    if pos < y_matched.len() {
                        y_matched[pos] = true;
                        if !p.is_loss() && rank < y_rank[pos] {
                            y_rank[pos] = rank;
                        }
                    }
                }
            }

            let ion_charge = 1.0_f64;
            let h2o_target = p.mz - H2O_LOSS_DA / ion_charge;
            let nh3_target = p.mz - NH3_LOSS_DA / ion_charge;
            if scored_spec.nearest_peak_full(h2o_target, tol_da).is_some()
                || scored_spec.nearest_peak_full(nh3_target, tol_da).is_some()
            {
                neutral_loss_ion_count += 1;
            }
        }
    }

    let num_matched: u32 = (b_matched.iter().filter(|&&m| m).count()
        + y_matched.iter().filter(|&&m| m).count()) as u32;

    let predicted_z2 = predict_by_ions_with_losses(peptide, 2..=2, predict_losses);
    let mut doubly_charged_matched_ion_count: u32 = 0;
    for p in &predicted_z2 {
        let tol_da = if feature_tol_is_ppm {
            p.mz * feature_tol / 1e6
        } else {
            feature_tol
        };
        if scored_spec.nearest_peak_full(p.mz, tol_da).is_some() {
            doubly_charged_matched_ion_count += 1;
        }
    }

    fn longest_run(matched: &[bool]) -> u32 {
        let mut best = 0u32;
        let mut cur = 0u32;
        for &m in matched {
            if m {
                cur += 1;
                if cur > best {
                    best = cur;
                }
            } else {
                cur = 0;
            }
        }
        best
    }

    let parent_mass = scored_spec.parent_mass();
    let num_segments = scorer.param().num_segments.max(1) as usize;

    let mut b_any_matched: SmallVec<[bool; 64]> = smallvec![false; n - 1];
    let mut y_any_matched: SmallVec<[bool; 64]> = smallvec![false; n - 1];
    let mut sum_prefix_intensity: f64 = 0.0;
    let mut sum_suffix_intensity: f64 = 0.0;

    let mut prm_accurate: f64 = 0.0;
    let mut srm_accurate: f64 = 0.0;

    let segment_ions: SmallVec<[&[scoring_crate::param_model::IonType]; 8]> = (0..num_segments)
        .map(|seg| {
            scorer
                .param()
                .ion_types_for_partition_slice(charge, parent_mass, seg)
        })
        .collect();

    for i in 0..(n - 1) {
        let aa_n = &peptide.residues[i];
        let aa_c = &peptide.residues[n - 1 - i];
        prm_accurate += aa_n.mass + aa_n.mod_.as_ref().map_or(0.0, |m| m.mass_delta);
        srm_accurate += aa_c.mass + aa_c.mod_.as_ref().map_or(0.0, |m| m.mass_delta);

        let mut b_any_this = false;
        let mut y_any_this = false;

        for seg in 0..num_segments {
            let ions = segment_ions[seg];
            for &ion in ions {
                let (is_prefix, residue_mass) = match ion {
                    scoring_crate::param_model::IonType::Prefix {
                        charge: ic,
                        offset_bits,
                        ..
                    } => {
                        let offset = f32::from_bits(offset_bits) as f64;
                        let z = ic as f64;
                        (true, (prm_accurate / z + offset, ion))
                    }
                    scoring_crate::param_model::IonType::Suffix {
                        charge: ic,
                        offset_bits,
                        ..
                    } => {
                        let offset = f32::from_bits(offset_bits) as f64;
                        let z = ic as f64;
                        (false, (srm_accurate / z + offset, ion))
                    }
                    scoring_crate::param_model::IonType::Noise => continue,
                };
                let theo_mz = residue_mass.0;
                if scorer.param().segment_num(theo_mz, parent_mass) != seg {
                    continue;
                }
                let tol_da = if feature_tol_is_ppm {
                    theo_mz * feature_tol / 1e6
                } else {
                    feature_tol
                };
                if let Some((_rank, intensity, _peak_mz)) =
                    scored_spec.nearest_peak_full(theo_mz, tol_da)
                {
                    if is_prefix {
                        sum_prefix_intensity += intensity as f64;
                        b_any_this = true;
                    } else {
                        sum_suffix_intensity += intensity as f64;
                        y_any_this = true;
                    }
                }
            }
        }

        b_any_matched[i] = b_any_this;
        y_any_matched[i] = y_any_this;
    }

    let longest_b = longest_run(&b_any_matched);
    let longest_y = longest_run(&y_any_matched);

    let total_intensity = scored_spec.total_intensity();
    let matched_b_intensity: f64 = sum_prefix_intensity;
    let matched_y_intensity: f64 = sum_suffix_intensity;
    let matched_total = matched_b_intensity + matched_y_intensity;

    let safe_div = |num: f64, denom: f64| -> f32 {
        if denom > 0.0 {
            (num / denom) as f32
        } else {
            0.0
        }
    };

    let explained_ion_current_ratio = safe_div(matched_total, total_intensity);
    let n_term_ion_current_ratio = safe_div(matched_b_intensity, total_intensity);
    let c_term_ion_current_ratio = safe_div(matched_y_intensity, total_intensity);
    let ms2_ion_current = if total_intensity > 0.0 {
        total_intensity as f32
    } else {
        0.0
    };
    let isolation_window_efficiency = 0.0_f32;

    matched_ions.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
    let top7 = &matched_ions[..matched_ions.len().min(7)];

    let abs_ppm_errors: Vec<f64> = top7
        .iter()
        .filter(|&&(_, _, pred, _, _)| pred > 0.0)
        .map(|&(_, obs, pred, _, _)| ((obs - pred) / pred * 1e6).abs())
        .collect();
    let rel_ppm_errors: Vec<f64> = top7
        .iter()
        .filter(|&&(_, _, pred, _, _)| pred > 0.0)
        .map(|&(_, obs, pred, _, _)| (obs - pred) / pred * 1e6)
        .collect();

    fn mean_and_pop_stdev(values: &[f64]) -> (f32, f32) {
        if values.is_empty() {
            return (0.0, 0.0);
        }
        let n = values.len() as f64;
        let mean = values.iter().sum::<f64>() / n;
        let sum_sq: f64 = values.iter().map(|v| v * v).sum();
        let var = (sum_sq / n - mean * mean).max(0.0);
        (mean as f32, var.sqrt() as f32)
    }

    let (mean_error_top7, stdev_error_top7) = mean_and_pop_stdev(&abs_ppm_errors);
    let (mean_rel_error_top7, stdev_rel_error_top7) = mean_and_pop_stdev(&rel_ppm_errors);

    let (ppm_centre, ppm_sigma) = scorer.ppm_gaussian();
    let ppm_gaussian_score: f32 = matched_ions
        .iter()
        .filter(|&&(_, _, pred, _, _)| pred > 0.0)
        .map(|&(_, obs, pred, _, _)| {
            let z = ((obs - pred) / pred * 1e6 - ppm_centre) / ppm_sigma;
            (-0.5 * z * z).exp()
        })
        .sum::<f64>() as f32;

    let chance_match_surprise: f32 = matched_ions
        .iter()
        .filter(|&&(_, _, pred, _, _)| pred > 0.0)
        .map(|&(_, obs, pred, _, _)| {
            let tol_da = if feature_tol_is_ppm {
                pred * feature_tol / 1e6
            } else {
                feature_tol
            };
            let rho = scored_spec.local_peak_density(obs, DENSITY_HW);
            let p_chance = rho * (2.0 * tol_da);
            if p_chance > 0.0 {
                (-p_chance.ln()).max(0.0)
            } else {
                0.0
            }
        })
        .sum::<f64>() as f32;

    let mut complementary_sites: SmallVec<[bool; 64]> = smallvec![false; n - 1];
    for i in 1..n {
        complementary_sites[i - 1] = b_any_matched[i - 1] && y_any_matched[n - i - 1];
    }
    let longest_complementary_ladder = longest_run(&complementary_sites);

    let mut complementary_ion_balance = 0.0_f32;
    for i in 1..n {
        let rb = b_rank[i - 1];
        let ry = y_rank[n - i - 1];
        if rb != u32::MAX && ry != u32::MAX {
            let diff = (rb as i64 - ry as i64).unsigned_abs() as f32;
            complementary_ion_balance += 1.0 / (1.0 + diff);
        }
    }

    let mean_matched_intensity_rank = if matched_ion_ranks.is_empty() {
        0.0
    } else {
        matched_ion_ranks.iter().map(|&r| r as f32).sum::<f32>() / matched_ion_ranks.len() as f32
    };

    let charge1_predicted = predict_by_ions(peptide, 1..=1);
    let theoretical_ions = theoretical_fragment_ions(&charge1_predicted);
    let unique_match_fraction = if num_matched == 0 {
        0.0
    } else {
        let evidence_sum: f64 = matched_ions
            .iter()
            .map(|&(_, obs, _, _, ion)| {
                let ambiguity = fragment_ambiguity(
                    &theoretical_ions,
                    obs,
                    feature_tol,
                    feature_tol_is_ppm,
                    ion,
                );
                1.0 / (1.0 + ambiguity as f64)
            })
            .sum();
        (evidence_sum / num_matched as f64) as f32
    };
    let mass_competition_evidence_val = mass_competition_evidence(
        scored_spec,
        &matched_ions,
        &theoretical_ions,
        feature_tol,
        feature_tol_is_ppm,
    );

    let frag_intensity_model = scorer.param().frag_intensity_model.as_deref();
    let frag_pred: Option<Vec<f64>> = frag_intensity_model
        .filter(|_| peptide.length() >= 2)
        .map(|g| predict_frag_with(g, peptide, charge, frag_cache));
    let frag_pred_slice = frag_pred.as_deref();
    let intensity_signal_val = intensity_signal(
        frag_intensity_model,
        scored_spec,
        peptide,
        charge,
        model_tol,
        model_tol_is_ppm,
        frag_pred_slice,
    );

    let kept = if full {
        kept_psm_feature_values(scored_spec, peptide, scorer, charge, frag_pred_slice)
    } else {
        KeptPsmFeatureValues::default()
    };
    let KeptPsmFeatureValues {
        frag_pred_explained,
        frag_pred_chance_llr,
        frag_topk_observed,
        rich_ion_llr: rich_ion_llr_val,
        mod_site,
    } = kept;

    PsmFeatures {
        num_matched_main_ions: num_matched,
        longest_b,
        longest_y,
        longest_y_pct: longest_y as f32 / n as f32,
        matched_ion_ratio: num_matched as f32 / n as f32,
        explained_ion_current_ratio,
        n_term_ion_current_ratio,
        c_term_ion_current_ratio,
        ms2_ion_current,
        isolation_window_efficiency,
        mean_error_top7,
        stdev_error_top7,
        mean_rel_error_top7,
        stdev_rel_error_top7,
        edge_score,
        precursor_isotope_kl: 0.0,
        precursor_snr: 0.0,
        delta_raw_score: 0.0,
        tailor_score: 0.0,
        rank_score_float: 0.0,
        ppm_gaussian_score,
        longest_complementary_ladder,
        complementary_ion_balance,
        neutral_loss_ion_count,
        mean_matched_intensity_rank,
        doubly_charged_matched_ion_count,
        chance_match_surprise,
        unique_match_fraction,
        intensity_signal: intensity_signal_val,
        frag_pred_explained,
        frag_pred_chance_llr,
        frag_topk_observed,
        rich_ion_llr: rich_ion_llr_val,
        is_refinement: 0,
        num_mods: 0,
        refine_mod_class: 0,
        mod_site_shifted_matched: mod_site[FEAT_SITE_SHIFTED_MATCHED],
        mod_site_shifted_frac: mod_site[FEAT_SITE_SHIFTED_FRAC],
        mod_site_intens_frac: mod_site[FEAT_SITE_INTENS_FRAC],
        mod_site_localized: mod_site[FEAT_SITE_LOCALIZED],
        mod_site_det_count: mod_site[FEAT_SITE_DET_COUNT],
        mass_competition_evidence: mass_competition_evidence_val,
        candidate_rank_entropy: 0.0,
        listwise_score_gap: 0.0,
        strong_score: 0.0,
        strong_score_cal: 0.0,
        delta_rt: 0.0,
        abs_delta_rt: 0.0,
        delta_rt_norm: 0.0,
        ref_competitor: 0,
        ref_best_raw: 0.0,
        is_variant: 0,
        is_isoform: 0,
        is_chain: 0,
        delta_ref_score: 0.0,
        isotope_mimic: 0,
        predicted_rt_min: 0.0,
        isobaric_rt_margin: 0.0,
        open_site: 0,
        open_loc_score: 0.0,
        open_loc_margin: 0.0,
        open_delta: 0.0,
        open_used_shift: false,
    }
}
