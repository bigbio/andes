//! Glyco retention-time (RT) PIN-feature wiring (design "Commit 2").
//!
//! This is the glyco counterpart of [`crate::rt_wiring`]. It reuses the SAME
//! engine-wide backbone RT-index predictor (`scoring::rt_model`) and per-run
//! linear calibration, then adds the glyco-only pieces:
//!
//! 1. The predicted index is `glyco_index(backbone_index, glycan_comp, coeffs)`
//!    — the backbone index PLUS the per-monosaccharide RP offset
//!    (`andes_glyco::rt_offset`), where NeuAc/NeuGc carry the OPPOSITE sign to
//!    the neutral hexoses (they retain the glycopeptide LATER on RP).
//! 2. The three shared PIN features (`DeltaRT`/`AbsDeltaRT`/`DeltaRTNorm`) are
//!    populated in place on each glyco hit's `PsmMatch::features` so the glyco
//!    PIN — which reuses `pin::psm_feature_values` — emits them (they are 0.0
//!    today because the glyco path never calls `rt_wiring`).
//! 3. A glyco-only `DeltaRTRank` feature (computed by the glyco PIN writer, not
//!    here) ranks a scan's competing glycoform candidates by `AbsDeltaRT` — the
//!    GlycReSoft isobaric-glycoform-disambiguation trick.
//!
//! # Baseline safety (non-negotiable)
//! With observed RT unavailable (`Spectrum.rt_seconds` is `None`) or the per-run
//! calibration unfittable (fewer than [`MIN_CALIBRATION_ANCHORS`] anchors /
//! degenerate indices), all RT features stay at their `Default` `0.0` — the
//! NEUTRAL value — so the glyco PIN is byte-identical to the pre-RT baseline on
//! those columns.
//!
//! # Determinism (sacred in this codebase)
//! Anchors are gathered by iterating `results` in slice order and, per result,
//! the scan's hit order. No `HashMap`/`HashSet` on any path that feeds the
//! calibration or the emitted feature values.

use model::spectrum::Spectrum;
use scoring_crate::rt_model::{calibrate, RtIndexModel};
use search::candidate_gen::Candidate;
use search::glyco_search::GlycoSpectrumResult;

use andes_glyco::glycan_db::GlycanComp;
use andes_glyco::rt_offset::{glyco_index, GlycanRtCoeffs};

/// Seconds-per-minute conversion for `rt_seconds` → minutes.
const SECONDS_PER_MINUTE: f64 = 60.0;

/// Predicted glyco RT INDEX for one hit: backbone index (Kyte-Doolittle seed)
/// plus the per-monosaccharide offset. De-novo hits with no resolved glycan add
/// a zero offset (so they reduce to the bare backbone index).
fn predict_glyco_index(
    model: &RtIndexModel,
    peptide: &model::peptide::Peptide,
    glycan: Option<&GlycanComp>,
    coeffs: &GlycanRtCoeffs,
) -> f64 {
    let backbone = model.predict_peptide(peptide);
    match glycan {
        Some(comp) => glyco_index(backbone, comp, coeffs),
        // No composition ⇒ empty offset; same code path serves both.
        None => backbone,
    }
}

/// Populate the glyco RT PIN features (`DeltaRT`/`AbsDeltaRT`/`DeltaRTNorm` and
/// `predicted_rt_min`) on every glyco hit in `results`, in place.
///
/// Steps (mirror `rt_wiring::populate_rt_features`, glyco index substituted):
/// 1. Gather calibration anchors `(predicted_index, observed_rt_min)` from
///    RANK-1 TARGET glyco hits (the scan's best hit, `is_decoy == false`), in a
///    stable (result-order, then hit-order) order. No hashing.
/// 2. Fit the per-run line `rt_obs = slope·index + intercept` via [`calibrate`]
///    (anchors in MINUTES ⇒ the line yields minutes directly). On failure leave
///    every RT feature at the neutral `0.0` and return.
/// 3. Gradient span = max−min observed RT (minutes) over emitting spectra.
/// 4. For each hit: predict its glyco index, map to `predicted_rt_min`, set
///    `delta_rt = obs − pred`, `abs_delta_rt`, `delta_rt_norm = delta / span`.
///    Hits on spectra without an observed RT keep the neutral `0.0`.
pub fn populate_glyco_rt_features(
    spectra: &[Spectrum],
    results: &mut [GlycoSpectrumResult],
    candidates: &[Candidate],
) {
    let model = RtIndexModel::seed();
    let coeffs = GlycanRtCoeffs::default();

    // ── Step 1: gather calibration anchors from rank-1 TARGET glyco hits ─────
    // The scan's rank-1 hit is `hits[0]` (glyco_search emits per-scan hits in
    // descending rank order). Stable order: results in slice order, then the
    // scan's first hit. No hashing.
    let mut anchors: Vec<(f64, f64)> = Vec::new();
    for result in results.iter() {
        let spec_idx = result.spectrum_idx;
        if spec_idx >= spectra.len() {
            continue;
        }
        let obs_rt_min = match spectra[spec_idx].rt_seconds {
            Some(s) => s / SECONDS_PER_MINUTE,
            None => continue,
        };
        let Some(hit) = result.hits.first() else {
            continue;
        };
        let cand_idx = hit.psm.primary_candidate_idx() as usize;
        if cand_idx >= candidates.len() {
            continue;
        }
        let cand = &candidates[cand_idx];
        if cand.is_decoy {
            continue; // targets only
        }
        let index =
            predict_glyco_index(&model, &cand.peptide, hit.glycan_key.glycan.as_ref(), &coeffs);
        anchors.push((index, obs_rt_min));
    }

    // ── Step 2: per-run calibration (anchors in minutes ⇒ line in minutes) ──
    let (slope, intercept) = match calibrate(&anchors) {
        Some(si) => si,
        None => return, // too few anchors / degenerate ⇒ leave neutral 0.0
    };

    // ── Step 3: gradient span = max−min observed RT (minutes) over the run ──
    let mut min_rt = f64::INFINITY;
    let mut max_rt = f64::NEG_INFINITY;
    for result in results.iter() {
        if result.hits.is_empty() {
            continue;
        }
        if let Some(s) = spectra.get(result.spectrum_idx).and_then(|sp| sp.rt_seconds) {
            let m = s / SECONDS_PER_MINUTE;
            min_rt = min_rt.min(m);
            max_rt = max_rt.max(m);
        }
    }
    let span = max_rt - min_rt; // ≤ 0 (single RT / none) ⇒ norm = 0.0

    // ── Step 4: write features onto each glyco hit (neutral when no RT) ─────
    for result in results.iter_mut() {
        let spec_idx = result.spectrum_idx;
        let obs_rt_min = spectra
            .get(spec_idx)
            .and_then(|sp| sp.rt_seconds)
            .map(|s| s / SECONDS_PER_MINUTE);
        let obs = match obs_rt_min {
            Some(o) => o,
            None => continue, // no observed RT ⇒ leave every hit neutral 0.0
        };
        for hit in result.hits.iter_mut() {
            let cand_idx = hit.psm.primary_candidate_idx() as usize;
            if cand_idx >= candidates.len() {
                continue;
            }
            let index = predict_glyco_index(
                &model,
                &candidates[cand_idx].peptide,
                hit.glycan_key.glycan.as_ref(),
                &coeffs,
            );
            let pred = slope * index + intercept;
            let delta = obs - pred;
            let norm = if span > 0.0 { delta / span } else { 0.0 };
            hit.psm.features.predicted_rt_min = pred as f32;
            hit.psm.features.delta_rt = delta as f32;
            hit.psm.features.abs_delta_rt = delta.abs() as f32;
            hit.psm.features.delta_rt_norm = norm as f32;
        }
    }
}

/// Compute the glyco-only `DeltaRTRank` for one hit within its scan.
///
/// Within a scan, rank the scan's candidate hits by `AbsDeltaRT` ascending
/// (the best RT match ranks first). The returned value is the 1-based rank
/// NORMALIZED by the candidate count (`rank / n`), so it lies in `(0, 1]` and is
/// scale-free — the GlycReSoft isobaric-glycoform disambiguation feature.
///
/// Neutral `0.0` when the scan has fewer than 2 candidate hits (no competition)
/// or when the emitted hit is not found in the candidate list.
///
/// Deterministic tie-break: hits with equal `AbsDeltaRT` are ordered by their
/// glycan key (HexNAc,Hex,Fuc,NeuAc,NeuGc then de-novo-last), then by original
/// hit index, so equal RT deltas produce a stable rank.
pub fn delta_rt_rank(hits: &[search::glyco_search::FullGlycoPsm], emitted_hit_idx: usize) -> f64 {
    let n = hits.len();
    if n < 2 {
        return 0.0; // no competition ⇒ neutral
    }
    // Sort a (abs_delta_rt, glycan_tiebreak_key, orig_idx) index list ascending.
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| {
        let ka = hits[a].psm.features.abs_delta_rt;
        let kb = hits[b].psm.features.abs_delta_rt;
        ka.total_cmp(&kb)
            .then_with(|| glycan_tiebreak_key(&hits[a]).cmp(&glycan_tiebreak_key(&hits[b])))
            .then(a.cmp(&b))
    });
    // 1-based rank of the emitted hit within the sorted order.
    match order.iter().position(|&i| i == emitted_hit_idx) {
        Some(pos) => (pos + 1) as f64 / n as f64,
        None => 0.0,
    }
}

/// Deterministic secondary tie-break key for `DeltaRTRank`: the glycan
/// composition tuple (de-novo hits with no glycan sort last via a max sentinel).
fn glycan_tiebreak_key(hit: &search::glyco_search::FullGlycoPsm) -> (u8, u8, u8, u8, u8, u8) {
    match &hit.glycan_key.glycan {
        // Leading 0 flags "has glycan" so composed glycans precede de-novo.
        Some(g) => (0, g.hexnac, g.hex, g.fuc, g.neuac, g.neugc),
        None => (1, 0, 0, 0, 0, 0),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use andes_glyco::glyco_psm::GlycoPsmKey;
    use andes_glyco::hybrid::Source;
    use andes_glyco::rt_offset::glycan_index_offset;
    use model::amino_acid::AminoAcid;
    use model::peptide::Peptide;
    use scoring_crate::rt_model::MIN_CALIBRATION_ANCHORS;
    use search::glyco_search::FullGlycoPsm;
    use search::psm::{PsmFeatures, PsmMatch};

    fn pep(seq: &[u8]) -> Peptide {
        let residues: Vec<_> = seq.iter().map(|&r| AminoAcid::standard(r).unwrap()).collect();
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
        }
    }

    fn glycan(hexnac: u8, hex: u8, fuc: u8, neuac: u8, neugc: u8) -> GlycanComp {
        GlycanComp { hexnac, hex, fuc, neuac, neugc, mass: 0.0 }
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
            y_ladder_decoy_score: 0.0,
            y0y1_anchor_score: 0.0,
            sialic_consistency: 0.0,
            core_y_hits: 0,
            backbone_mass: 0.0,
            is_transferred: false,
            transfer_graph_support: 0,
            transfer_seed_score: 0.0,
            transfer_rt_delta: 0.0,
            transfer_ungated: false,
        }
    }

    fn psm(candidate_idx: u32, rank_score: f32) -> PsmMatch {
        PsmMatch {
            spectrum_idx: 0,
            candidate_idxs: vec![candidate_idx],
            charge_used: 2,
            mass_error_ppm: 0.0,
            score: rank_score,
            rank_score,
            edge_score: 0,
            activation_method: None,
            features: PsmFeatures::default(),
            isotope_offset: 0,
            precursor_mz_override: None,
        }
    }

    fn hit(candidate_idx: u32, glycan: Option<GlycanComp>) -> FullGlycoPsm {
        FullGlycoPsm { glycan_key: key(glycan), psm: psm(candidate_idx, 10.0) }
    }

    // ── DeltaRTRank: known AbsDeltaRT → correct normalized ranks ──────────────

    #[test]
    fn delta_rt_rank_assigns_normalized_ranks_by_abs_delta() {
        // 3 candidate hits with distinct AbsDeltaRT: 5.0, 1.0, 3.0.
        // Ascending order: idx1 (1.0) < idx2 (3.0) < idx0 (5.0).
        let mut hits = vec![
            hit(0, Some(glycan(2, 3, 0, 0, 0))),
            hit(0, Some(glycan(2, 4, 0, 0, 0))),
            hit(0, Some(glycan(2, 5, 0, 0, 0))),
        ];
        hits[0].psm.features.abs_delta_rt = 5.0;
        hits[1].psm.features.abs_delta_rt = 1.0;
        hits[2].psm.features.abs_delta_rt = 3.0;

        // rank/n: idx1 → 1/3, idx2 → 2/3, idx0 → 3/3.
        assert!((delta_rt_rank(&hits, 1) - 1.0 / 3.0).abs() < 1e-9);
        assert!((delta_rt_rank(&hits, 2) - 2.0 / 3.0).abs() < 1e-9);
        assert!((delta_rt_rank(&hits, 0) - 3.0 / 3.0).abs() < 1e-9);
    }

    #[test]
    fn delta_rt_rank_deterministic_tie_break_by_glycan_key() {
        // Two hits with IDENTICAL AbsDeltaRT ⇒ tie broken by glycan key.
        // glycan(2,3,...) < glycan(2,5,...), so the (2,3) hit ranks first.
        let mut hits = vec![
            hit(0, Some(glycan(2, 5, 0, 0, 0))), // idx0
            hit(0, Some(glycan(2, 3, 0, 0, 0))), // idx1 (smaller glycan key)
        ];
        hits[0].psm.features.abs_delta_rt = 2.0;
        hits[1].psm.features.abs_delta_rt = 2.0;
        // idx1 (smaller key) → rank 1/2; idx0 → 2/2.
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

    // ── DeltaRT: known calibration + known glycan offset ─────────────────────

    /// With observed RT = slope·glyco_index + intercept for every scan, the
    /// calibration recovers (slope, intercept) and DeltaRT ≈ 0. The glyco index
    /// must include the per-monosaccharide offset (NeuAc opposite sign): a
    /// sialylated glycoform's predicted index is HIGHER than the same neutral
    /// composition's, so if the offset were dropped DeltaRT would be non-zero.
    #[test]
    fn delta_rt_populates_with_glycan_offset_and_neuac_opposite_sign() {
        let model = RtIndexModel::seed();
        let coeffs = GlycanRtCoeffs::default();
        let (slope, intercept) = (2.0_f64, 10.0_f64);

        // 8 distinct peptides × distinct glycans so indices spread apart.
        let seqs: Vec<&[u8]> = vec![
            b"LLIIVV", b"PEPTIDE", b"AAAAAA", b"WFWFWF", b"DEKRDE", b"GGGGGG", b"STSTST", b"MMMMMM",
        ];
        let mut candidates = Vec::new();
        let mut spectra = Vec::new();
        let mut results = Vec::new();
        for (i, seq) in seqs.iter().enumerate() {
            let p = pep(seq);
            let gl = glycan(2, 3, 0, (i % 3) as u8, 0); // vary NeuAc 0..2
            let idx = glyco_index(model.predict_peptide(&p), &gl, &coeffs);
            let rt_min = slope * idx + intercept;
            candidates.push(cand(p, false));
            spectra.push(spec(i as i32, Some(rt_min * SECONDS_PER_MINUTE)));
            results.push(GlycoSpectrumResult {
                spectrum_idx: i,
                hits: vec![hit(i as u32, Some(gl))],
            });
        }

        populate_glyco_rt_features(&spectra, &mut results, &candidates);

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

        // Sanity: the NeuAc offset is positive (opposite the neutral sugars), so
        // adding a NeuAc raises the predicted index — the mechanism DeltaRT
        // depends on. (Guards against silently dropping the glycan offset.)
        let neutral = glycan_index_offset(&glycan(2, 3, 0, 0, 0), &coeffs);
        let sialic = glycan_index_offset(&glycan(2, 3, 0, 1, 0), &coeffs);
        assert!(sialic > neutral, "NeuAc must raise the index (opposite sign)");
    }

    /// A hand-checkable single perturbed scan: DeltaRT = observed − predicted.
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
        // Perturb scan 0 by +5 min.
        let gl0 = glycan(2, 3, 1, 0, 0);
        let idx0 = glyco_index(model.predict_peptide(&candidates[0].peptide), &gl0, &coeffs);
        spectra[0].rt_seconds = Some((slope * idx0 + intercept + 5.0) * SECONDS_PER_MINUTE);

        populate_glyco_rt_features(&spectra, &mut results, &candidates);

        let f = &results[0].hits[0].psm.features;
        assert!(
            (f.delta_rt - 5.0).abs() < 0.5,
            "DeltaRT should be ≈ +5 on the perturbed scan, got {}",
            f.delta_rt
        );
    }

    // ── Neutral paths ────────────────────────────────────────────────────────

    #[test]
    fn neutral_when_rt_unavailable() {
        let mut candidates = Vec::new();
        let mut spectra = Vec::new();
        let mut results = Vec::new();
        for i in 0..8 {
            candidates.push(cand(pep(b"PEPTIDE"), false));
            spectra.push(spec(i, None)); // rt_seconds = None
            results.push(GlycoSpectrumResult {
                spectrum_idx: i as usize,
                hits: vec![hit(i as u32, Some(glycan(2, 3, 0, 0, 0)))],
            });
        }
        populate_glyco_rt_features(&spectra, &mut results, &candidates);
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
            // DISTINCT peptides (varying trailing residue) so the predicted
            // indices differ — the calibration is non-degenerate on every axis
            // EXCEPT anchor count, isolating the `< MIN_CALIBRATION_ANCHORS`
            // threshold as the sole reason the RT features stay neutral (CodeRabbit).
            let mut seq = b"PEPTIDE".to_vec();
            seq.push(b"ACDEFGHIKLMNPQRSTVWY"[i % 20]);
            candidates.push(cand(pep(&seq), false));
            spectra.push(spec(i as i32, Some((i as f64 + 1.0) * 60.0)));
            results.push(GlycoSpectrumResult {
                spectrum_idx: i,
                hits: vec![hit(i as u32, Some(glycan(2, 3, 0, 0, 0)))],
            });
        }
        populate_glyco_rt_features(&spectra, &mut results, &candidates);
        for r in &results {
            let f = &r.hits[0].psm.features;
            assert_eq!(f.delta_rt, 0.0);
            assert_eq!(f.abs_delta_rt, 0.0);
            assert_eq!(f.delta_rt_norm, 0.0);
        }
    }

    #[test]
    fn decoy_rank1_hits_are_not_anchors() {
        // All rank-1 hits are decoys ⇒ zero anchors ⇒ calibration fails ⇒ neutral.
        let seqs: Vec<&[u8]> = vec![
            b"LLIIVV", b"PEPTIDE", b"AAAAAA", b"WFWFWF", b"DEKRDE", b"GGGGGG", b"STSTST", b"MMMMMM",
        ];
        let mut candidates = Vec::new();
        let mut spectra = Vec::new();
        let mut results = Vec::new();
        for (i, seq) in seqs.iter().enumerate() {
            candidates.push(cand(pep(seq), true)); // decoy
            spectra.push(spec(i as i32, Some((i as f64 + 1.0) * 60.0)));
            results.push(GlycoSpectrumResult {
                spectrum_idx: i,
                hits: vec![hit(i as u32, Some(glycan(2, 3, 0, 0, 0)))],
            });
        }
        populate_glyco_rt_features(&spectra, &mut results, &candidates);
        for r in &results {
            assert_eq!(r.hits[0].psm.features.delta_rt, 0.0, "decoy-only ⇒ no calibration");
        }
    }
}
