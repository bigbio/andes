//! Engine-wide retention-time (RT) PIN-feature wiring.
//!
//! This is the *wiring* layer for the pure, tested RT-prediction core in
//! `scoring::rt_model`. It runs ONCE per search, after all PSMs are scored and
//! before the PIN/QPX writers, and populates the additive RT feature fields on
//! each PSM (`PsmFeatures::{delta_rt, abs_delta_rt, delta_rt_norm,
//! predicted_rt_min}`) so the SHARED `pin::psm_feature_values` vector — the
//! single source of truth for both the PIN and the QPX `.idparquet` — can emit
//! them without any extra plumbing.
//!
//! Design reference: `docs/plans/glyco/50-roadmap/rt-prediction-design.md`
//! ("Commit 1 (engine-wide)"). This implements the engine-wide backbone RT
//! index + per-run linear calibration ONLY; the glyco per-monosaccharide offset
//! and `DeltaRTRank` (Commit 2) are a separate later task.
//!
//! # Baseline safety (non-negotiable)
//! When observed RT is unavailable (`Spectrum.rt_seconds` is `None`) or the
//! per-run calibration cannot be fit (fewer than
//! [`MIN_CALIBRATION_ANCHORS`](scoring_crate::rt_model::MIN_CALIBRATION_ANCHORS)
//! confident anchors, or degenerate indices), the three PIN features stay at
//! their `Default` `0.0` — the NEUTRAL value — so the PIN is byte-identical to
//! the pre-RT baseline on these columns.
//!
//! # Determinism (sacred in this codebase)
//! Anchors are gathered by iterating `spectra`/`queues` in index order and, per
//! spectrum, the queue's rank-sorted PSM order — no `HashMap`/`HashSet` on any
//! path that feeds the calibration or the emitted feature values.

use model::spectrum::Spectrum;
use scoring_crate::rt_model::{calibrate, RtIndexModel};
use search::candidate_gen::Candidate;
use search::psm::TopNQueue;

use crate::row_context::iter_ranked_by_rank_score;

/// Seconds-per-minute conversion for `rt_seconds` → minutes (features are in
/// minutes, keeping f64 precision until the final `f32` store).
const SECONDS_PER_MINUTE: f64 = 60.0;

/// Populate the engine-wide RT PIN features on every PSM in `queues`.
///
/// Steps (design "Commit 1"):
/// 1. Collect calibration anchors `(predicted_index, observed_rt_min)` from all
///    RANK-1 TARGET PSMs, in a stable (spectrum-index, then rank-sorted) order.
/// 2. Fit the per-run line `rt_obs = slope·index + intercept` via
///    [`calibrate`]. On failure (too few anchors / degenerate), leave every RT
///    feature at the neutral `0.0` and return.
/// 3. Compute the run's gradient span = max−min observed RT (minutes) over all
///    emitting spectra with an RT.
/// 4. For each PSM: predict its index, map to `predicted_rt_min` via the
///    calibration, and set `delta_rt = obs − pred`, `abs_delta_rt`,
///    `delta_rt_norm = delta_rt / span` (0.0 when span ≤ 0). PSMs on spectra
///    without an observed RT keep the neutral `0.0`.
///
/// `spectra` and `queues` are parallel slices (`queues[i]` holds the top-N PSMs
/// for `spectra[i]`); `candidates` is the pool the PSMs' `primary_candidate_idx`
/// resolves against (used for the peptide and the target/decoy label).
pub fn populate_rt_features(
    spectra: &[Spectrum],
    queues: &mut [TopNQueue],
    candidates: &[Candidate],
) {
    // The seed RT-index model (Kyte-Doolittle prior). The calibration line
    // absorbs the arbitrary index scale, so the seed is sufficient for a
    // per-run relative RT ordering; a fitted/GBDT model can drop in later.
    let model = RtIndexModel::seed();

    // ── Step 1: gather calibration anchors from rank-1 target PSMs ──────────
    // Stable order: spectra in index order, PSMs in rank-sorted order. No hashing.
    let mut anchors: Vec<(f64, f64)> = Vec::new();
    for (spec_idx, queue) in queues.iter().enumerate() {
        if queue.is_empty() {
            continue;
        }
        let obs_rt_min = match spectra[spec_idx].rt_seconds {
            Some(s) => s / SECONDS_PER_MINUTE,
            None => continue, // no observed RT ⇒ can't anchor on this spectrum
        };
        let sorted = queue.clone().into_rank_sorted_vec();
        for (rank, psm) in iter_ranked_by_rank_score(&sorted) {
            if rank != 1 {
                break; // only the rank-1 row anchors; queue is rank-sorted
            }
            let cand = &candidates[psm.primary_candidate_idx() as usize];
            if cand.is_decoy {
                continue; // skip a decoy rank-1, try the next rank-1 target
            }
            let index = model.predict_peptide(&cand.peptide);
            anchors.push((index, obs_rt_min));
            // ONE anchor per spectrum: tied rank-1 targets (equal rank_score,
            // distinct peptides) would otherwise inject several points sharing a
            // single observed RT — a vertical cluster that biases the OLS slope
            // toward 0 for that RT (review finding). The first rank-1 target is
            // the representative anchor.
            break;
        }
    }

    // ── Step 2: per-run calibration. `calibrate` is unit-agnostic — it fits
    // `y = slope·x + intercept` on whatever `y` it is given. We pass `y` in
    // MINUTES (not the seconds `scoring::rt_model::predicted_rt` assumes), so
    // the fitted line yields minutes directly and we map with `slope·index +
    // intercept` below (no /60). ─────────────────────────────────────────────
    let (slope, intercept) = match calibrate(&anchors) {
        Some(si) => si,
        None => return, // too few anchors / degenerate ⇒ leave neutral 0.0
    };

    // ── Step 3: gradient span = max−min observed RT (minutes) over the run ──
    let mut min_rt = f64::INFINITY;
    let mut max_rt = f64::NEG_INFINITY;
    for (spec_idx, queue) in queues.iter().enumerate() {
        if queue.is_empty() {
            continue;
        }
        // Guard against a non-finite observed RT poisoning the span (Codex review):
        // an inf/NaN `rt_seconds` was excluded from the calibration fit but would
        // still make max/min — and hence every DeltaRTNorm — non-finite here.
        if let Some(s) = spectra[spec_idx].rt_seconds {
            if s.is_finite() {
                let m = s / SECONDS_PER_MINUTE;
                min_rt = min_rt.min(m);
                max_rt = max_rt.max(m);
            }
        }
    }
    let span = max_rt - min_rt; // may be ≤ 0 (single RT / none) ⇒ norm = 0.0

    // ── Step 4: write features onto each PSM (neutral when no observed RT) ──
    for (spec_idx, queue) in queues.iter_mut().enumerate() {
        if queue.is_empty() {
            continue;
        }
        // Only a FINITE observed RT anchors a feature; a non-finite `rt_seconds`
        // leaves the row neutral (0.0) rather than writing inf/NaN into the PIN.
        let obs_rt_min = spectra[spec_idx]
            .rt_seconds
            .filter(|s| s.is_finite())
            .map(|s| s / SECONDS_PER_MINUTE);
        queue.fill_post_topn(|psm| {
            let obs = match obs_rt_min {
                Some(o) => o,
                None => return, // no / non-finite observed RT ⇒ leave neutral 0.0
            };
            let cand = &candidates[psm.primary_candidate_idx() as usize];
            let index = model.predict_peptide(&cand.peptide);
            // Anchors were fit in minutes, so the line yields minutes directly.
            let pred = slope * index + intercept;
            let delta = obs - pred;
            // Belt-and-braces: never emit a non-finite RT feature (leave neutral).
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

    /// Build an unmodified peptide from a byte sequence.
    fn pep(seq: &[u8]) -> Peptide {
        let residues: Vec<_> = seq
            .iter()
            .map(|&r| AminoAcid::standard(r).unwrap())
            .collect();
        Peptide::new(residues, b'K', b'-')
    }

    /// A target/decoy candidate carrying `peptide`.
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

    /// A spectrum with an optional RT in seconds.
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

    /// A rank-1 PSM pointing at `candidate_idx` with the given ranking score.
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

    fn queue_with(psms: Vec<PsmMatch>) -> TopNQueue {
        let mut q = TopNQueue::new(10);
        for p in psms {
            q.push(p);
        }
        q
    }

    /// Build `n` distinct target candidates whose predicted indices span a range,
    /// plus spectra whose observed RT is an exact linear function of that index —
    /// so the calibration recovers a known (slope, intercept) and DeltaRT is
    /// computable in closed form.
    fn linear_fixture(
        n: usize,
        slope: f64,
        intercept: f64,
    ) -> (Vec<Spectrum>, Vec<TopNQueue>, Vec<Candidate>, RtIndexModel) {
        let model = RtIndexModel::seed();
        // A set of chemically-distinct peptides to spread indices apart.
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
        // Observed RT = 2.0·index + 10.0 (minutes) exactly for every spectrum,
        // so a perfect fit ⇒ DeltaRT ≈ 0 on the anchors themselves. To get a
        // NON-zero, checkable DeltaRT, offset ONE spectrum's observed RT.
        let (mut spectra, mut queues, candidates, model) = linear_fixture(8, 2.0, 10.0);

        // Perturb spectrum 0's observed RT by +5 min so its DeltaRT is known.
        // (One perturbed point barely moves the trimmed fit.)
        let idx0 = model.predict_peptide(&candidates[0].peptide);
        let perturbed_rt_min = 2.0 * idx0 + 10.0 + 5.0;
        spectra[0].rt_seconds = Some(perturbed_rt_min * SECONDS_PER_MINUTE);

        populate_rt_features(&spectra, &mut queues, &candidates);

        // The calibration should recover ~ (slope=2, intercept=10). PSM 0's
        // predicted RT ≈ 2·idx0 + 10, observed = that + 5 ⇒ DeltaRT ≈ +5.
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
        // DeltaRTNorm = DeltaRT / gradient span. Span = max−min observed RT.
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
    fn delta_rt_equals_obs_minus_pred_exactly() {
        // Exact-fit fixture: DeltaRT on each anchor should be ≈ 0 because
        // observed RT is an exact linear function of the predicted index.
        let (spectra, mut queues, candidates, _model) = linear_fixture(8, 1.5, 20.0);
        populate_rt_features(&spectra, &mut queues, &candidates);
        for q in &queues {
            let f = q.iter_psms().next().unwrap().features.clone();
            assert!(
                f.delta_rt.abs() < 0.5,
                "exact-fit DeltaRT should be ≈ 0, got {}",
                f.delta_rt
            );
            // predicted_rt_min must be populated (non-zero for a real fit).
            assert!(f.predicted_rt_min != 0.0, "predicted_rt_min should be set");
        }
    }

    #[test]
    fn neutral_when_rt_unavailable() {
        // No spectrum has an observed RT ⇒ calibration cannot run ⇒ all three
        // features stay at the neutral 0.0 (baseline-identical).
        let mut candidates = Vec::new();
        let mut spectra = Vec::new();
        let mut queues = Vec::new();
        for i in 0..8 {
            candidates.push(cand(pep(b"PEPTIDE"), false));
            spectra.push(spec(i, None)); // rt_seconds = None
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
        // A run with enough FINITE anchors to calibrate, PLUS spectra whose
        // rt_seconds is inf / NaN. The fit succeeds on the finite anchors, but the
        // bad-RT rows must stay neutral 0.0 — a non-finite rt_seconds must NEVER
        // leak into DeltaRT/AbsDeltaRT/DeltaRTNorm/predicted_rt (Codex review).
        let mut candidates = Vec::new();
        let mut spectra = Vec::new();
        let mut queues = Vec::new();
        for i in 0..8 {
            // Distinct peptides → non-degenerate calibration.
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

        // No feature anywhere may be inf/NaN.
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
        // The two bad-RT rows are exactly neutral.
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
        // Fewer than MIN_CALIBRATION_ANCHORS target anchors ⇒ calibrate returns
        // None ⇒ neutral 0.0 everywhere, no crash.
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
        // Two runs, identical except one uses decoy candidates for the SAME
        // peptides. The decoy run has zero valid anchors ⇒ calibration fails ⇒
        // neutral, proving decoys are excluded from the anchor set.
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
        // Same PSMs in ⇒ identical emitted features across repeated runs.
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
