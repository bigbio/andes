use model::spectrum::Spectrum;
use scoring_crate::rt_model::{calibrate, RtIndexModel};
use search::candidates::Candidate;
use search::glycopeptides::GlycoSpectrumResult;

use glyco::glycan_db::GlycanComp;
use glyco::rt_offset::{glyco_index, GlycanRtCoeffs};

const SECONDS_PER_MINUTE: f64 = 60.0;

fn predict_glyco_index(
    model: &RtIndexModel,
    peptide: &model::peptide::Peptide,
    glycan: Option<&GlycanComp>,
    coeffs: &GlycanRtCoeffs,
) -> f64 {
    let backbone = model.predict_peptide(peptide);
    match glycan {
        Some(comp) => glyco_index(backbone, comp, coeffs),
        None => backbone,
    }
}

pub fn populate_glyco_rt_features(
    spectra: &[Spectrum],
    results: &mut [GlycoSpectrumResult],
    candidates: &[Candidate],
    glycan_list: &[GlycanComp],
) {
    let model = RtIndexModel::seed();
    let coeffs = GlycanRtCoeffs::default();
    const ISOBARIC_TOL_DA: f64 = 0.05;

    let mut anchors: Vec<(f64, f64)> = Vec::new();
    for result in results.iter() {
        let spec_idx = result.spectrum_idx;
        if spec_idx >= spectra.len() {
            continue;
        }
        let obs_rt_min = match spectra[spec_idx].rt_seconds {
            Some(s) if s.is_finite() => s / SECONDS_PER_MINUTE,
            _ => continue,
        };
        for hit in result.hits.iter() {
            let cand_idx = hit.psm.primary_candidate_idx() as usize;
            if cand_idx >= candidates.len() {
                continue;
            }
            if candidates[cand_idx].is_decoy {
                continue;
            }
            let index = predict_glyco_index(
                &model,
                &candidates[cand_idx].peptide,
                hit.glycan_key.glycan.as_ref(),
                &coeffs,
            );
            anchors.push((index, obs_rt_min));
            break;
        }
    }

    let (slope, intercept) = match calibrate(&anchors) {
        Some(si) => si,
        None => return,
    };

    let mut min_rt = f64::INFINITY;
    let mut max_rt = f64::NEG_INFINITY;
    for result in results.iter() {
        if result.hits.is_empty() {
            continue;
        }
        if let Some(s) = spectra
            .get(result.spectrum_idx)
            .and_then(|sp| sp.rt_seconds)
        {
            if s.is_finite() {
                let m = s / SECONDS_PER_MINUTE;
                min_rt = min_rt.min(m);
                max_rt = max_rt.max(m);
            }
        }
    }
    let span = max_rt - min_rt;

    for result in results.iter_mut() {
        let spec_idx = result.spectrum_idx;
        let obs_rt_min = spectra
            .get(spec_idx)
            .and_then(|sp| sp.rt_seconds)
            .filter(|s| s.is_finite())
            .map(|s| s / SECONDS_PER_MINUTE);
        let obs = match obs_rt_min {
            Some(o) => o,
            None => continue,
        };
        for hit in result.hits.iter_mut() {
            let cand_idx = hit.psm.primary_candidate_idx() as usize;
            if cand_idx >= candidates.len() {
                continue;
            }
            let backbone = model.predict_peptide(&candidates[cand_idx].peptide);
            let index = match hit.glycan_key.glycan.as_ref() {
                Some(comp) => glyco_index(backbone, comp, &coeffs),
                None => backbone,
            };
            let pred = slope * index + intercept;
            let delta = obs - pred;
            if !pred.is_finite() || !delta.is_finite() {
                continue;
            }
            let norm = if span > 0.0 { delta / span } else { 0.0 };
            hit.psm.features.predicted_rt_min = pred as f32;
            hit.psm.features.delta_rt = delta as f32;
            hit.psm.features.abs_delta_rt = delta.abs() as f32;
            hit.psm.features.delta_rt_norm = norm as f32;

            if let Some(g) = hit.glycan_key.glycan.as_ref() {
                let assigned_res = delta.abs();
                let mut best_alt_res = f64::INFINITY;
                for alt in glycan_list {
                    if (alt.mass - g.mass).abs() > ISOBARIC_TOL_DA {
                        continue;
                    }
                    if alt.hexnac == g.hexnac
                        && alt.hex == g.hex
                        && alt.fuc == g.fuc
                        && alt.neuac == g.neuac
                        && alt.neugc == g.neugc
                    {
                        continue;
                    }
                    let alt_pred = slope * glyco_index(backbone, alt, &coeffs) + intercept;
                    let alt_res = (obs - alt_pred).abs();
                    if alt_res.is_finite() && alt_res < best_alt_res {
                        best_alt_res = alt_res;
                    }
                }
                if best_alt_res.is_finite() {
                    hit.psm.features.isobaric_rt_margin = (best_alt_res - assigned_res) as f32;
                }
            }
        }
    }
}

pub fn delta_rt_rank(hits: &[search::glycopeptides::FullGlycoPsm], emitted_hit_idx: usize) -> f64 {
    let n = hits.len();
    if n < 2 {
        return 0.0;
    }
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| {
        let ka = hits[a].psm.features.abs_delta_rt;
        let kb = hits[b].psm.features.abs_delta_rt;
        ka.total_cmp(&kb)
            .then_with(|| glycan_tiebreak_key(&hits[a]).cmp(&glycan_tiebreak_key(&hits[b])))
            .then(a.cmp(&b))
    });
    match order.iter().position(|&i| i == emitted_hit_idx) {
        Some(pos) => (pos + 1) as f64 / n as f64,
        None => 0.0,
    }
}

fn glycan_tiebreak_key(hit: &search::glycopeptides::FullGlycoPsm) -> (u8, u8, u8, u8, u8, u8) {
    match &hit.glycan_key.glycan {
        Some(g) => (0, g.hexnac, g.hex, g.fuc, g.neuac, g.neugc),
        None => (1, 0, 0, 0, 0, 0),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use glyco::glyco_psm::GlycoPsmKey;
    use glyco::hybrid::Source;
    use glyco::rt_offset::glycan_index_offset;
    use model::amino_acid::AminoAcid;
    use model::peptide::Peptide;
    use scoring_crate::rt_model::MIN_CALIBRATION_ANCHORS;
    use search::glycopeptides::FullGlycoPsm;
    use search::psm::{PsmFeatures, PsmMatch};

    fn pep(seq: &[u8]) -> Peptide {
        let residues: Vec<_> = seq
            .iter()
            .map(|&r| AminoAcid::standard(r).unwrap())
            .collect();
        Peptide::new(residues, b'K', b'-')
    }

    fn cand(peptide: Peptide, is_decoy: bool) -> Candidate {
        Candidate {
            peptide,
            protein_index: 0,
            start_offset_in_protein: 0,
            is_decoy,
            is_protein_n_term: false,
            is_protein_c_term: false,
        }
    }

    fn spec(scan: i32, rt_seconds: Option<f64>) -> Spectrum {
        Spectrum {
            title: format!("scan={scan}"),
            precursor_mz: 500.0,
            precursor_intensity: None,
            precursor_charge: Some(2),
            rt_seconds,
            scan: Some(scan),
            peaks: vec![],
            activation_method: None,
            isolation_lower_offset: None,
            isolation_upper_offset: None,
            isolation_target_mz: None,
            mobility_window: None,
        }
    }

    fn glycan(hexnac: u8, hex: u8, fuc: u8, neuac: u8, neugc: u8) -> GlycanComp {
        GlycanComp::from_counts(hexnac, hex, fuc, neuac, neugc, 0)
    }

    fn key(glycan: Option<GlycanComp>) -> GlycoPsmKey {
        GlycoPsmKey {
            spectrum_idx: 0,
            glycan_mass: glycan.as_ref().map(|g| g.mass).unwrap_or(0.0),
            glycan,
            glycan_source: Source::Db,
            oxonium_summed_frac: 0.0,
            n_core_oxonium_ions: 0,
            y_ladder_intensity_score: 0.0,
            y_hit_frac: 0.0,
            partial_glycan_by: 0.0,
            y0y1_anchor_score: 0.0,
            sialic_consistency: 0.0,
            core_y_hits: 0,
            backbone_mass: 0.0,
            cz_hyperscore: 0.0,
            cz_intensity: 0.0,
            cz_explained: 0.0,
            cz_chance_llr: 0.0,
        }
    }

    fn psm(candidate_idx: u32, rank_score: f32) -> PsmMatch {
        PsmMatch {
            spectrum_idx: 0,
            candidate_idxs: vec![candidate_idx].into(),
            charge_used: 2,
            mass_error_ppm: 0.0,
            score: rank_score,
            rank_score,
            edge_score: 0,
            activation_method: None,
            features: PsmFeatures::default(),
            isotope_offset: 0,
            precursor_mz_override: None,
            site_evidence: None,
            tie_key: 0,
        }
    }

    fn hit(candidate_idx: u32, glycan: Option<GlycanComp>) -> FullGlycoPsm {
        FullGlycoPsm {
            glycan_key: key(glycan),
            psm: psm(candidate_idx, 10.0),
        }
    }

    #[test]
    fn glyco_anchor_skips_decoy_first_hit_to_the_target() {
        let n = MIN_CALIBRATION_ANCHORS + 2;
        let mut candidates = Vec::new();
        let mut spectra = Vec::new();
        let mut results = Vec::new();
        for i in 0..n {
            let d = candidates.len() as u32;
            candidates.push(cand(pep(b"PEPTIDEDK"), true));
            let mut seq = b"TARGETPEP".to_vec();
            seq.push(b"ACDEFGHIKLMNPQRSTVWY"[i % 20]);
            let t = candidates.len() as u32;
            candidates.push(cand(pep(&seq), false));
            spectra.push(spec(i as i32, Some((i as f64 + 10.0) * 60.0)));
            results.push(GlycoSpectrumResult {
                spectrum_idx: i,
                hits: vec![
                    hit(d, Some(glycan(2, 3, 0, 0, 0))),
                    hit(t, Some(glycan(2, 3, 0, 0, 0))),
                ],
            });
        }
        populate_glyco_rt_features(&spectra, &mut results, &candidates, &[]);
        let populated = results.iter().any(|r| {
            r.hits
                .iter()
                .any(|h| h.psm.features.predicted_rt_min != 0.0)
        });
        assert!(
            populated,
            "decoy-first scans must still calibrate via the next target hit"
        );
    }

    #[test]
    fn delta_rt_rank_assigns_normalized_ranks_by_abs_delta() {
        let mut hits = vec![
            hit(0, Some(glycan(2, 3, 0, 0, 0))),
            hit(0, Some(glycan(2, 4, 0, 0, 0))),
            hit(0, Some(glycan(2, 5, 0, 0, 0))),
        ];
        hits[0].psm.features.abs_delta_rt = 5.0;
        hits[1].psm.features.abs_delta_rt = 1.0;
        hits[2].psm.features.abs_delta_rt = 3.0;

        assert!((delta_rt_rank(&hits, 1) - 1.0 / 3.0).abs() < 1e-9);
        assert!((delta_rt_rank(&hits, 2) - 2.0 / 3.0).abs() < 1e-9);
        assert!((delta_rt_rank(&hits, 0) - 3.0 / 3.0).abs() < 1e-9);
    }

    #[test]
    fn delta_rt_rank_deterministic_tie_break_by_glycan_key() {
        let mut hits = vec![
            hit(0, Some(glycan(2, 5, 0, 0, 0))),
            hit(0, Some(glycan(2, 3, 0, 0, 0))),
        ];
        hits[0].psm.features.abs_delta_rt = 2.0;
        hits[1].psm.features.abs_delta_rt = 2.0;
        assert!((delta_rt_rank(&hits, 1) - 0.5).abs() < 1e-9);
        assert!((delta_rt_rank(&hits, 0) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn delta_rt_rank_neutral_when_fewer_than_two_candidates() {
        let hits = vec![hit(0, Some(glycan(2, 3, 0, 0, 0)))];
        assert_eq!(delta_rt_rank(&hits, 0), 0.0);
        let empty: Vec<FullGlycoPsm> = vec![];
        assert_eq!(delta_rt_rank(&empty, 0), 0.0);
    }

    #[test]
    fn delta_rt_populates_with_glycan_offset_and_neuac_opposite_sign() {
        let model = RtIndexModel::seed();
        let coeffs = GlycanRtCoeffs::default();
        let (slope, intercept) = (2.0_f64, 10.0_f64);

        let seqs: Vec<&[u8]> = vec![
            b"LLIIVV", b"PEPTIDE", b"AAAAAA", b"WFWFWF", b"DEKRDE", b"GGGGGG", b"STSTST", b"MMMMMM",
        ];
        let mut candidates = Vec::new();
        let mut spectra = Vec::new();
        let mut results = Vec::new();
        for (i, seq) in seqs.iter().enumerate() {
            let p = pep(seq);
            let gl = glycan(2, 3, 0, (i % 3) as u8, 0);
            let idx = glyco_index(model.predict_peptide(&p), &gl, &coeffs);
            let rt_min = slope * idx + intercept;
            candidates.push(cand(p, false));
            spectra.push(spec(i as i32, Some(rt_min * SECONDS_PER_MINUTE)));
            results.push(GlycoSpectrumResult {
                spectrum_idx: i,
                hits: vec![hit(i as u32, Some(gl))],
            });
        }

        populate_glyco_rt_features(&spectra, &mut results, &candidates, &[]);

        for r in &results {
            let f = &r.hits[0].psm.features;
            assert!(
                f.delta_rt.abs() < 0.5,
                "exact-fit glyco DeltaRT should be ≈ 0, got {}",
                f.delta_rt
            );
            assert!(f.predicted_rt_min != 0.0, "predicted_rt_min must be set");
            assert!((f.abs_delta_rt - f.delta_rt.abs()).abs() < 1e-6);
        }

        let neutral = glycan_index_offset(&glycan(2, 3, 0, 0, 0), &coeffs);
        let sialic = glycan_index_offset(&glycan(2, 3, 0, 1, 0), &coeffs);
        assert!(
            sialic > neutral,
            "NeuAc must raise the index (opposite sign)"
        );
    }

    #[test]
    fn delta_rt_equals_obs_minus_pred_on_perturbed_scan() {
        let model = RtIndexModel::seed();
        let coeffs = GlycanRtCoeffs::default();
        let (slope, intercept) = (2.0_f64, 10.0_f64);
        let seqs: Vec<&[u8]> = vec![
            b"LLIIVV", b"PEPTIDE", b"AAAAAA", b"WFWFWF", b"DEKRDE", b"GGGGGG", b"STSTST", b"MMMMMM",
        ];
        let mut candidates = Vec::new();
        let mut spectra = Vec::new();
        let mut results = Vec::new();
        for (i, seq) in seqs.iter().enumerate() {
            let p = pep(seq);
            let gl = glycan(2, 3, 1, 0, 0);
            let idx = glyco_index(model.predict_peptide(&p), &gl, &coeffs);
            let rt_min = slope * idx + intercept;
            candidates.push(cand(p, false));
            spectra.push(spec(i as i32, Some(rt_min * SECONDS_PER_MINUTE)));
            results.push(GlycoSpectrumResult {
                spectrum_idx: i,
                hits: vec![hit(i as u32, Some(gl))],
            });
        }
        let gl0 = glycan(2, 3, 1, 0, 0);
        let idx0 = glyco_index(model.predict_peptide(&candidates[0].peptide), &gl0, &coeffs);
        spectra[0].rt_seconds = Some((slope * idx0 + intercept + 5.0) * SECONDS_PER_MINUTE);

        populate_glyco_rt_features(&spectra, &mut results, &candidates, &[]);

        let f = &results[0].hits[0].psm.features;
        assert!(
            (f.delta_rt - 5.0).abs() < 0.5,
            "DeltaRT should be ≈ +5 on the perturbed scan, got {}",
            f.delta_rt
        );
    }

    #[test]
    fn isobaric_rt_margin_positive_when_assigned_fits_better() {
        let model = RtIndexModel::seed();
        let coeffs = GlycanRtCoeffs::default();
        let (slope, intercept) = (2.0_f64, 10.0_f64);
        let assigned = glycan(2, 4, 0, 1, 0);
        let alt = glycan(2, 3, 1, 0, 1);
        let glycan_list = vec![assigned.clone(), alt.clone()];
        assert!(
            (glyco_index(1.0, &assigned, &coeffs) - glyco_index(1.0, &alt, &coeffs)).abs() > 1e-9,
            "isobaric pair must have distinct RT offsets for the test to be meaningful"
        );
        let seqs: Vec<&[u8]> = vec![
            b"LLIIVV", b"PEPTIDE", b"AAAAAA", b"WFWFWF", b"DEKRDE", b"GGGGGG", b"STSTST", b"MMMMMM",
        ];
        let (mut candidates, mut spectra, mut results) = (Vec::new(), Vec::new(), Vec::new());
        for (i, seq) in seqs.iter().enumerate() {
            let p = pep(seq);
            let idx = glyco_index(model.predict_peptide(&p), &assigned, &coeffs);
            let rt_min = slope * idx + intercept;
            candidates.push(cand(p, false));
            spectra.push(spec(i as i32, Some(rt_min * SECONDS_PER_MINUTE)));
            results.push(GlycoSpectrumResult {
                spectrum_idx: i,
                hits: vec![hit(i as u32, Some(assigned.clone()))],
            });
        }
        populate_glyco_rt_features(&spectra, &mut results, &candidates, &glycan_list);
        let m = results[0].hits[0].psm.features.isobaric_rt_margin;
        assert!(
            m > 0.0,
            "assigned fits RT better than the isobaric alternative → margin>0, got {m}"
        );
    }

    #[test]
    fn neutral_when_rt_unavailable() {
        let mut candidates = Vec::new();
        let mut spectra = Vec::new();
        let mut results = Vec::new();
        for i in 0..8 {
            candidates.push(cand(pep(b"PEPTIDE"), false));
            spectra.push(spec(i, None));
            results.push(GlycoSpectrumResult {
                spectrum_idx: i as usize,
                hits: vec![hit(i as u32, Some(glycan(2, 3, 0, 0, 0)))],
            });
        }
        populate_glyco_rt_features(&spectra, &mut results, &candidates, &[]);
        for r in &results {
            let f = &r.hits[0].psm.features;
            assert_eq!(f.delta_rt, 0.0);
            assert_eq!(f.abs_delta_rt, 0.0);
            assert_eq!(f.delta_rt_norm, 0.0);
            assert_eq!(f.predicted_rt_min, 0.0);
        }
    }

    #[test]
    fn neutral_when_too_few_anchors() {
        let n = MIN_CALIBRATION_ANCHORS - 1;
        let mut candidates = Vec::new();
        let mut spectra = Vec::new();
        let mut results = Vec::new();
        for i in 0..n {
            let mut seq = b"PEPTIDE".to_vec();
            seq.push(b"ACDEFGHIKLMNPQRSTVWY"[i % 20]);
            candidates.push(cand(pep(&seq), false));
            spectra.push(spec(i as i32, Some((i as f64 + 1.0) * 60.0)));
            results.push(GlycoSpectrumResult {
                spectrum_idx: i,
                hits: vec![hit(i as u32, Some(glycan(2, 3, 0, 0, 0)))],
            });
        }
        populate_glyco_rt_features(&spectra, &mut results, &candidates, &[]);
        for r in &results {
            let f = &r.hits[0].psm.features;
            assert_eq!(f.delta_rt, 0.0);
            assert_eq!(f.abs_delta_rt, 0.0);
            assert_eq!(f.delta_rt_norm, 0.0);
        }
    }

    #[test]
    fn decoy_rank1_hits_are_not_anchors() {
        let seqs: Vec<&[u8]> = vec![
            b"LLIIVV", b"PEPTIDE", b"AAAAAA", b"WFWFWF", b"DEKRDE", b"GGGGGG", b"STSTST", b"MMMMMM",
        ];
        let mut candidates = Vec::new();
        let mut spectra = Vec::new();
        let mut results = Vec::new();
        for (i, seq) in seqs.iter().enumerate() {
            candidates.push(cand(pep(seq), true));
            spectra.push(spec(i as i32, Some((i as f64 + 1.0) * 60.0)));
            results.push(GlycoSpectrumResult {
                spectrum_idx: i,
                hits: vec![hit(i as u32, Some(glycan(2, 3, 0, 0, 0)))],
            });
        }
        populate_glyco_rt_features(&spectra, &mut results, &candidates, &[]);
        for r in &results {
            assert_eq!(
                r.hits[0].psm.features.delta_rt, 0.0,
                "decoy-only ⇒ no calibration"
            );
        }
    }
}
