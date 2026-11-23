use model::spectrum::Spectrum;
use scoring_crate::rt_model::{calibrate, RtIndexModel};
use search::candidates::Candidate;
use search::psm::TopNQueue;

use crate::rows::iter_ranked_by_rank_score;

const SECONDS_PER_MINUTE: f64 = 60.0;

pub fn populate_rt_features_per_run(
    spectra: &[Spectrum],
    queues: &mut [TopNQueue],
    candidates: &[Candidate],
    run_ends: &[usize],
) {
    let mut start = 0;
    for &end in run_ends {
        let end = end.min(spectra.len()).max(start);
        populate_rt_features(&spectra[start..end], &mut queues[start..end], candidates);
        start = end;
    }
    if start < spectra.len() {
        populate_rt_features(&spectra[start..], &mut queues[start..], candidates);
    }
}

pub fn populate_rt_features(
    spectra: &[Spectrum],
    queues: &mut [TopNQueue],
    candidates: &[Candidate],
) {
    let model = RtIndexModel::seed();

    let mut anchors: Vec<(f64, f64)> = Vec::new();
    for (spec_idx, queue) in queues.iter().enumerate() {
        if queue.is_empty() {
            continue;
        }
        let obs_rt_min = match spectra[spec_idx].rt_seconds {
            Some(s) => s / SECONDS_PER_MINUTE,
            None => continue,
        };
        let sorted = queue.clone().into_rank_sorted_vec();
        for (rank, psm) in iter_ranked_by_rank_score(&sorted) {
            if rank != 1 {
                break;
            }
            let cand = &candidates[psm.primary_candidate_idx() as usize];
            if cand.is_decoy {
                continue;
            }
            let index = model.predict_peptide(&cand.peptide);
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
    for (spec_idx, queue) in queues.iter().enumerate() {
        if queue.is_empty() {
            continue;
        }
        if let Some(s) = spectra[spec_idx].rt_seconds {
            if s.is_finite() {
                let m = s / SECONDS_PER_MINUTE;
                min_rt = min_rt.min(m);
                max_rt = max_rt.max(m);
            }
        }
    }
    let span = max_rt - min_rt;

    for (spec_idx, queue) in queues.iter_mut().enumerate() {
        if queue.is_empty() {
            continue;
        }
        let obs_rt_min = spectra[spec_idx]
            .rt_seconds
            .filter(|s| s.is_finite())
            .map(|s| s / SECONDS_PER_MINUTE);
        queue.fill_post_topn(|psm| {
            let obs = match obs_rt_min {
                Some(o) => o,
                None => return,
            };
            let cand = &candidates[psm.primary_candidate_idx() as usize];
            let index = model.predict_peptide(&cand.peptide);
            let pred = slope * index + intercept;
            let delta = obs - pred;
            if !pred.is_finite() || !delta.is_finite() {
                return;
            }
            let norm = if span > 0.0 { delta / span } else { 0.0 };
            psm.features.predicted_rt_min = pred as f32;
            psm.features.delta_rt = delta as f32;
            psm.features.abs_delta_rt = delta.abs() as f32;
            psm.features.delta_rt_norm = norm as f32;
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::amino_acid::AminoAcid;
    use model::peptide::Peptide;
    use scoring_crate::rt_model::MIN_CALIBRATION_ANCHORS;
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

    fn queue_with(psms: Vec<PsmMatch>) -> TopNQueue {
        let mut q = TopNQueue::new(10);
        for p in psms {
            q.push(p);
        }
        q
    }

    fn linear_fixture(
        n: usize,
        slope: f64,
        intercept: f64,
    ) -> (Vec<Spectrum>, Vec<TopNQueue>, Vec<Candidate>, RtIndexModel) {
        let model = RtIndexModel::seed();
        let seqs: Vec<&[u8]> = vec![
            b"LLIIVV", b"PEPTIDE", b"AAAAAA", b"WFWFWF", b"DEKRDE", b"GGGGGG", b"STSTST",
            b"MMMMMM", b"YYYYYY", b"NQNQNQ",
        ];
        let mut candidates = Vec::new();
        let mut spectra = Vec::new();
        let mut queues = Vec::new();
        for i in 0..n {
            let p = pep(seqs[i % seqs.len()]);
            let index = model.predict_peptide(&p);
            let rt_min = slope * index + intercept;
            candidates.push(cand(p, false));
            spectra.push(spec(i as i32, Some(rt_min * SECONDS_PER_MINUTE)));
            queues.push(queue_with(vec![psm(i as u32, 100.0 - i as f32)]));
        }
        (spectra, queues, candidates, model)
    }

    #[test]
    fn emits_delta_abs_norm_with_known_calibration() {
        let (mut spectra, mut queues, candidates, model) = linear_fixture(8, 2.0, 10.0);

        let idx0 = model.predict_peptide(&candidates[0].peptide);
        let perturbed_rt_min = 2.0 * idx0 + 10.0 + 5.0;
        spectra[0].rt_seconds = Some(perturbed_rt_min * SECONDS_PER_MINUTE);

        populate_rt_features(&spectra, &mut queues, &candidates);

        let f = queues[0]
            .iter_psms()
            .next()
            .expect("psm present")
            .features
            .clone();
        assert!(
            (f.delta_rt - 5.0).abs() < 0.5,
            "delta_rt {} should be ≈ +5",
            f.delta_rt
        );
        assert!(
            (f.abs_delta_rt - f.delta_rt.abs()).abs() < 1e-6,
            "abs_delta_rt {} must equal |delta_rt|",
            f.abs_delta_rt
        );
        let mins: Vec<f64> = spectra
            .iter()
            .map(|s| s.rt_seconds.unwrap() / 60.0)
            .collect();
        let span = mins.iter().cloned().fold(f64::MIN, f64::max)
            - mins.iter().cloned().fold(f64::MAX, f64::min);
        let expected_norm = f.delta_rt as f64 / span;
        assert!(
            (f.delta_rt_norm as f64 - expected_norm).abs() < 1e-4,
            "delta_rt_norm {} should be delta/span {}",
            f.delta_rt_norm,
            expected_norm
        );
    }

    #[test]
    fn each_run_gets_its_own_calibration() {
        let (mut spectra, mut queues, mut candidates, _) = linear_fixture(8, 1.5, 20.0);
        let (s2, q2, c2, _) = linear_fixture(8, 4.0, 2.0);
        let shift = candidates.len() as u32;
        spectra.extend(s2);
        candidates.extend(c2);
        for q in q2 {
            let p = q.iter_psms().next().unwrap().clone();
            queues.push(queue_with(vec![psm(p.candidate_idxs[0] + shift, p.rank_score)]));
        }
        let worst = |queues: &[TopNQueue]| {
            queues
                .iter()
                .map(|q| q.iter_psms().next().unwrap().features.delta_rt.abs())
                .fold(0.0f32, f32::max)
        };
        let mut per_run = queues.clone();
        populate_rt_features_per_run(&spectra, &mut per_run, &candidates, &[8, 16]);
        assert!(worst(&per_run) < 0.5, "per-run DeltaRT {}", worst(&per_run));
        let mut pooled = queues.clone();
        populate_rt_features(&spectra, &mut pooled, &candidates);
        assert!(worst(&pooled) > 1.0, "pooled DeltaRT {}", worst(&pooled));
    }

    #[test]
    fn delta_rt_equals_obs_minus_pred_exactly() {
        let (spectra, mut queues, candidates, _model) = linear_fixture(8, 1.5, 20.0);
        populate_rt_features(&spectra, &mut queues, &candidates);
        for q in &queues {
            let f = q.iter_psms().next().unwrap().features.clone();
            assert!(
                f.delta_rt.abs() < 0.5,
                "exact-fit DeltaRT should be ≈ 0, got {}",
                f.delta_rt
            );
            assert!(f.predicted_rt_min != 0.0, "predicted_rt_min should be set");
        }
    }

    #[test]
    fn neutral_when_rt_unavailable() {
        let mut candidates = Vec::new();
        let mut spectra = Vec::new();
        let mut queues = Vec::new();
        for i in 0..8 {
            candidates.push(cand(pep(b"PEPTIDE"), false));
            spectra.push(spec(i, None));
            queues.push(queue_with(vec![psm(i as u32, 100.0)]));
        }
        populate_rt_features(&spectra, &mut queues, &candidates);
        for q in &queues {
            let f = q.iter_psms().next().unwrap().features.clone();
            assert_eq!(f.delta_rt, 0.0);
            assert_eq!(f.abs_delta_rt, 0.0);
            assert_eq!(f.delta_rt_norm, 0.0);
            assert_eq!(f.predicted_rt_min, 0.0);
        }
    }

    #[test]
    fn non_finite_rt_never_writes_inf_or_nan() {
        let mut candidates = Vec::new();
        let mut spectra = Vec::new();
        let mut queues = Vec::new();
        for i in 0..8 {
            let mut seq = b"PEPTIDE".to_vec();
            seq.push(b"ACDEFGHIKLMNPQRSTVWY"[i % 20]);
            candidates.push(cand(pep(&seq), false));
            spectra.push(spec(i as i32, Some((i as f64 + 10.0) * 60.0)));
            queues.push(queue_with(vec![psm(i as u32, 100.0)]));
        }
        let bad = candidates.len() as u32;
        candidates.push(cand(pep(b"PEPTIDEK"), false));
        spectra.push(spec(100, Some(f64::INFINITY)));
        queues.push(queue_with(vec![psm(bad, 100.0)]));
        candidates.push(cand(pep(b"PEPTIDER"), false));
        spectra.push(spec(101, Some(f64::NAN)));
        queues.push(queue_with(vec![psm(bad + 1, 100.0)]));

        populate_rt_features(&spectra, &mut queues, &candidates);

        for q in &queues {
            let f = q.iter_psms().next().unwrap().features.clone();
            assert!(
                f.delta_rt.is_finite()
                    && f.abs_delta_rt.is_finite()
                    && f.delta_rt_norm.is_finite()
                    && f.predicted_rt_min.is_finite(),
                "non-finite RT feature leaked: {:?}",
                (
                    f.delta_rt,
                    f.abs_delta_rt,
                    f.delta_rt_norm,
                    f.predicted_rt_min
                )
            );
        }
        for qi in [8usize, 9] {
            let f = queues[qi].iter_psms().next().unwrap().features.clone();
            assert_eq!(f.delta_rt, 0.0);
            assert_eq!(f.abs_delta_rt, 0.0);
            assert_eq!(f.delta_rt_norm, 0.0);
            assert_eq!(f.predicted_rt_min, 0.0);
        }
    }

    #[test]
    fn neutral_when_too_few_anchors() {
        let n = MIN_CALIBRATION_ANCHORS - 1;
        let (spectra, mut queues, candidates, _m) = linear_fixture(n, 2.0, 10.0);
        populate_rt_features(&spectra, &mut queues, &candidates);
        for q in &queues {
            let f = q.iter_psms().next().unwrap().features.clone();
            assert_eq!(f.delta_rt, 0.0);
            assert_eq!(f.abs_delta_rt, 0.0);
            assert_eq!(f.delta_rt_norm, 0.0);
        }
    }

    #[test]
    fn decoy_rank1_psms_are_not_anchors() {
        let (spectra, mut q_dec, mut candidates, _m) = linear_fixture(8, 2.0, 10.0);
        for c in &mut candidates {
            c.is_decoy = true;
        }
        populate_rt_features(&spectra, &mut q_dec, &candidates);
        for q in &q_dec {
            let f = q.iter_psms().next().unwrap().features.clone();
            assert_eq!(f.delta_rt, 0.0, "decoy-only run must not calibrate");
        }
    }

    #[test]
    fn anchor_gathering_is_deterministic() {
        let (spectra, mut a, candidates, _m) = linear_fixture(8, 2.0, 10.0);
        let (_s2, mut b, _c2, _m2) = linear_fixture(8, 2.0, 10.0);
        populate_rt_features(&spectra, &mut a, &candidates);
        populate_rt_features(&spectra, &mut b, &candidates);
        for (qa, qb) in a.iter().zip(b.iter()) {
            let fa = qa.iter_psms().next().unwrap().features.clone();
            let fb = qb.iter_psms().next().unwrap().features.clone();
            assert_eq!(fa.delta_rt, fb.delta_rt);
            assert_eq!(fa.abs_delta_rt, fb.abs_delta_rt);
            assert_eq!(fa.delta_rt_norm, fb.delta_rt_norm);
            assert_eq!(fa.predicted_rt_min, fb.predicted_rt_min);
        }
    }
}
