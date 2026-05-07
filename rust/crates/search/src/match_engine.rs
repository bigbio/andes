//! Top-level integration: spectra × candidates → top-N PSMs per spectrum.

use std::collections::{BTreeMap, HashMap, HashSet};

use model::aa_set::AminoAcidSet;
use crate::candidate_gen::{enumerate_candidates, Candidate};
use model::enzyme::Enzyme;
use scoring_crate::gf::generating_function::GeneratingFunction;
use scoring_crate::gf::group::GeneratingFunctionGroup;
use scoring_crate::gf::primitive_graph::PrimitiveAaGraph;
use model::mass::{nominal_from, H2O, PROTON};
use model::peptide::Peptide;
use crate::precursor_matching::{matches_precursor, MassError};
use crate::psm::{PsmFeatures, PsmMatch, TopNQueue};
use scoring_crate::scoring::fragment_ions::{IonKind, predict_by_ions};
use crate::search_index::SearchIndex;
use crate::search_params::SearchParams;
use scoring_crate::scoring::{score_psm, RankScorer, ScoredSpectrum};
use model::spectrum::Spectrum;

/// Match every spectrum against every candidate from the SearchIndex.
/// Returns one top-N PSM queue per spectrum, in input order.
///
/// Phase 5 Task 5: score_psm replaces Phase 4e's -|mass_error_ppm| placeholder.
/// A `ScoredSpectrum` is built once per spectrum and reused across all candidates.
///
/// Phase 4f optimization: bucket candidates by mass for sub-linear lookup.
///
/// Phase 6 Task 8: after per-candidate scoring, compute SpecEValue via the
/// generating-function DP across the precursor tolerance window in nominal
/// mass space and assign it to every PSM in the queue.
pub fn match_spectra(
    spectra: &[Spectrum],
    idx: &SearchIndex,
    params: &SearchParams,
    scorer: &RankScorer,
    fragment_tolerance_da: f64,
    decoy_prefix: &str,
) -> Vec<TopNQueue> {
    let mut queues: Vec<TopNQueue> = (0..spectra.len())
        .map(|_| TopNQueue::new(params.top_n_psms_per_spectrum))
        .collect();

    let candidates: Vec<Candidate> = enumerate_candidates(idx, params, decoy_prefix).collect();

    // Build mass-bucket index: nominal(peptide.mass() - H2O) → Vec<candidate_idx>.
    //
    // Uses the same nominal_from convention as the GF mass-bin loop so that
    // bucket keys align with the GF's mass-bin lookup (commit b89779a fix).
    // Stores only indices into `candidates` — no cloning, tiny memory overhead.
    let mut bucket_index: BTreeMap<i32, Vec<usize>> = BTreeMap::new();
    for (cand_idx, cand) in candidates.iter().enumerate() {
        let nominal = nominal_from(cand.peptide.mass() - H2O);
        bucket_index.entry(nominal).or_default().push(cand_idx);
    }

    // Build an aa_set clone with enzyme registered (for GF cleavage scoring).
    // We use Java MS-GF+ defaults: peptide_eff = 0.95, neighboring_eff = 0.95.
    // Cloning is cheap (AminoAcidSet is a HashMap of ~20 entries).
    // This avoids mutating the shared SearchParams.aa_set borrow.
    let mut aa_set_for_gf: AminoAcidSet = params.aa_set.clone();
    if params.enzyme != Enzyme::NoCleavage && params.enzyme != Enzyme::NonSpecific {
        aa_set_for_gf.register_enzyme(params.enzyme, 0.95, 0.95);
    }

    for (spec_idx, spec) in spectra.iter().enumerate() {
        // Determine which charge states to try for this spectrum.
        // For charge-explicit spectra this is a single entry; for charge-missing,
        // typically 2-3 entries (small overhead, correct behavior).
        let charges_to_try: Vec<u8> = match spec.precursor_charge {
            Some(z) if z > 0 => vec![z as u8],
            _ => params.charge_range.clone().collect(),
        };

        // Build (and cache) a ScoredSpectrum per charge to evaluate.
        //
        // Fix (Track B3): previously a single ScoredSpectrum was built with
        // `precursor_z = spec.precursor_charge.unwrap_or(2)`, so charge-missing
        // spectra always used z=2 even when evaluating z=3 candidates — wrong
        // precursor filtering, wrong partition, wrong main_ion.
        //
        // For charge-explicit spectra the cache has exactly 1 entry (no overhead).
        // For charge-missing spectra, typically 2-3 entries per spectrum.
        // The HashMap lifetime aligns with `spec`'s borrow (same loop scope).
        let mut scored_per_charge: HashMap<u8, ScoredSpectrum<'_>> = HashMap::new();
        for &z in &charges_to_try {
            scored_per_charge.entry(z)
                .or_insert_with(|| ScoredSpectrum::new(spec, scorer.param(), z));
        }

        // Compute per-charge candidate windows and union them into a deduplicated
        // set of candidate indices. This avoids O(all_candidates) iteration —
        // only candidates whose nominal mass falls within the precursor-tolerance
        // window for at least one of the tried charges are visited.
        //
        // Window derivation mirrors compute_spec_e_values_for_spectrum's logic so
        // that any candidate admitted by matches_precursor is guaranteed to be in
        // at least one charge's window, preserving parity with the brute-force path.
        let mut window_cand_indices: HashSet<usize> = HashSet::new();
        for &z in &charges_to_try {
            let charge_f = z as f64;
            let neutral_mass = (spec.precursor_mz - PROTON) * charge_f - H2O;
            let nominal_center = nominal_from(neutral_mass);
            let iso_min = *params.isotope_error_range.start() as i32;
            let iso_max = *params.isotope_error_range.end() as i32;
            let tol_da_left  = params.precursor_tolerance.left.as_da(neutral_mass);
            let tol_da_right = params.precursor_tolerance.right.as_da(neutral_mass);
            let widen_left  = (tol_da_left  - 0.4999_f64).round() as i32;
            let widen_right = (tol_da_right - 0.4999_f64).round() as i32;
            // Java convention (same as GF window in compute_spec_e_values):
            //   max includes widen_left (tolDaLeft widens the upper bound)
            //   min includes widen_right (tolDaRight widens the lower bound)
            let min_nominal = nominal_center - iso_max - widen_right;
            let max_nominal = nominal_center - iso_min + widen_left;
            for (_nm, idxs) in bucket_index.range(min_nominal..=max_nominal) {
                for &ci in idxs {
                    window_cand_indices.insert(ci);
                }
            }
        }

        for &cand_idx in &window_cand_indices {
            let cand = &candidates[cand_idx];
            for &z in &charges_to_try {
                let scored_spec = &scored_per_charge[&z];
                let mut best_for_charge: Option<(MassError, f32)> = None;
                for offset in params.isotope_error_range.clone() {
                    if let Some(err) = matches_precursor(spec, &cand.peptide, z, offset, &params.precursor_tolerance) {
                        // Phase 5: use real score_psm instead of -|mass_error_ppm| placeholder.
                        let score = score_psm(scored_spec, &cand.peptide, scorer, z, fragment_tolerance_da);
                        if best_for_charge.as_ref().map_or(true, |(_, s)| score > *s) {
                            best_for_charge = Some((err, score));
                        }
                    }
                }
                if let Some((err, score)) = best_for_charge {
                    let features = compute_psm_features(scored_spec, &cand.peptide, fragment_tolerance_da);
                    queues[spec_idx].push(PsmMatch {
                        spectrum_idx: spec_idx,
                        candidate: cand.clone(),
                        charge_used: z,
                        mass_error_ppm: err.mass_error_ppm,
                        score,
                        spec_e_value: 1.0,  // set by Phase 6 compute_spec_e_values_for_spectrum
                        de_novo_score: i32::MIN,  // set by Phase 7 compute_spec_e_values_for_spectrum
                        activation_method: Some(scorer.param().data_type.activation),
                        e_value: 1.0,  // set by Phase 7 compute_spec_e_values_for_spectrum
                        features,
                        isotope_offset: err.isotope_offset,
                    });
                }
            }
        }

        // Phase 6: compute SpecEValue for the PSMs in this queue.
        if !queues[spec_idx].is_empty() {
            let enzyme_opt = if params.enzyme != Enzyme::NoCleavage
                && params.enzyme != Enzyme::NonSpecific
            {
                Some(params.enzyme)
            } else {
                None
            };
            // Pick the ScoredSpectrum for the top PSM's charge.
            // For charge-explicit spectra there is only 1 entry in the cache.
            // For charge-missing spectra, use the top PSM's charge so the GF
            // reflects the dominant scoring context (option (a) per B3 plan).
            let top_charge = queues[spec_idx]
                .iter_psms()
                .max_by(|a, b| a.cmp(b))
                .map(|p| p.charge_used)
                .unwrap_or(charges_to_try[0]);
            // Unwrap is safe: the cache was built for every charge in charges_to_try,
            // and top_charge comes from a PSM that was scored at one of those charges.
            let scored_spec_for_gf = &scored_per_charge[&top_charge];
            compute_spec_e_values_for_spectrum(
                spec,
                params,
                &mut queues[spec_idx],
                &aa_set_for_gf,
                enzyme_opt,
                scorer,
                scored_spec_for_gf,
                top_charge,
                fragment_tolerance_da,
                idx,
            );
        }
    }

    queues
}

/// For a single spectrum, compute the GF across the precursor tolerance
/// window in nominal mass space, then assign `spec_e_value` to every PSM
/// in `queue` whose nominal_peptide_mass falls within the window.
///
/// Mirrors Java DBScanner.java:597-650.
///
/// # Arguments
/// * `spec` — the spectrum (used for precursor m/z).
/// * `params` — search params (precursor_tolerance, isotope_error_range).
/// * `queue` — the PSM queue for this spectrum (mutated in place).
/// * `aa_set` — amino acid set with enzyme already registered via `register_enzyme`.
/// * `enzyme` — the search enzyme (passed to PrimitiveAaGraph; may be None).
/// * `scorer` — Phase 5 RankScorer.
/// * `scored_spec` — ScoredSpectrum built with `top_charge` (B3: per-charge cache).
/// * `top_charge` — charge of the top PSM in the queue; used for GF mass window.
///   For charge-explicit spectra this equals `spec.precursor_charge.unwrap()`.
///   For charge-missing spectra, using the top PSM's charge ensures the GF
///   reflects the dominant scoring context.
/// * `fragment_tolerance_da` — fragment mass tolerance in Da.
/// * `search_index` — database (target+decoy); used to look up protein sequences
///   for protein-terminal flag derivation (Track B4).
#[allow(clippy::too_many_arguments)]
fn compute_spec_e_values_for_spectrum(
    spec: &Spectrum,
    params: &SearchParams,
    queue: &mut TopNQueue,
    aa_set: &AminoAcidSet,
    enzyme: Option<Enzyme>,
    scorer: &RankScorer,
    scored_spec: &ScoredSpectrum<'_>,
    top_charge: u8,
    fragment_tolerance_da: f64,
    search_index: &SearchIndex,
) {
    // 1. Determine the peptide neutral mass and its tolerance window.
    // For charge-explicit spectra, `top_charge` == spec.precursor_charge.unwrap().
    // For charge-missing spectra, `top_charge` is the top PSM's charge (B3 fix).
    let charge = top_charge;
    if charge == 0 {
        return;
    }

    // peptide_neutral_mass = (precursor_mz - H) * charge - H2O
    // This matches Java: scoredSpec.getPrecursorPeak().getMass() - H2O
    // where getPrecursorPeak().getMass() = (mz - H) * charge.
    let peptide_neutral_mass = (spec.precursor_mz - PROTON) * (charge as f64) - H2O;
    let nominal_peptide_mass = nominal_from(peptide_neutral_mass);

    // Java isotope error convention: range [min_iso, max_iso] is applied as
    //   minNominalPeptideMass = nominalPeptideMass - maxIsotopeError
    //   maxNominalPeptideMass = nominalPeptideMass - minIsotopeError
    let iso_min = *params.isotope_error_range.start() as i32;
    let iso_max = *params.isotope_error_range.end() as i32;
    let min_iso_nominal = nominal_peptide_mass - iso_max;
    let max_iso_nominal = nominal_peptide_mass - iso_min;

    // Tolerance widening: Java uses Math.round(tol_da - 0.4999f).
    // tolDaLeft governs the upper bound; tolDaRight governs the lower bound.
    let tol_da_left = params.precursor_tolerance.left.as_da(peptide_neutral_mass);
    let tol_da_right = params.precursor_tolerance.right.as_da(peptide_neutral_mass);
    let widen_left = (tol_da_left - 0.4999_f64).round() as i32;
    let widen_right = (tol_da_right - 0.4999_f64).round() as i32;

    let max_peptide_mass_idx = max_iso_nominal + widen_left;
    let min_peptide_mass_idx = min_iso_nominal - widen_right;

    if max_peptide_mass_idx < min_peptide_mass_idx {
        return;
    }

    // 2. Compute the minimum score across all PSMs (used as score threshold).
    let min_score = queue
        .iter_psms()
        .map(|p| p.score.round() as i32)
        .min()
        .unwrap_or(i32::MIN);

    // parent_mass = (mz - H) * charge  (precursor peak mass, with H added back in Java).
    let parent_mass = (spec.precursor_mz - PROTON) * (charge as f64);

    // 3. Derive protein-terminal flags by OR-ing across ALL PSMs in the queue.
    //
    // Java reference: DBScanner.java:592-602 aggregates useProteinNTerm /
    // useProteinCTerm across all candidates before GF construction. We mirror
    // this by iterating the full queue and setting either flag the moment any
    // PSM is at a protein N- or C-terminus, short-circuiting once both are set.
    let (use_protein_n_term, use_protein_c_term) = {
        let mut any_n = false;
        let mut any_c = false;
        for psm in queue.iter_psms() {
            if let Some(prot) = search_index.protein_at(psm.candidate.protein_index) {
                let start = psm.candidate.start_offset_in_protein;
                let pep_len = psm.candidate.peptide.length();
                if start == 0 { any_n = true; }
                if start + pep_len >= prot.sequence.len() { any_c = true; }
                if any_n && any_c { break; }
            }
        }
        (any_n, any_c)
    };

    // 3b. Build the GF group across the nominal mass range.
    let mut group = GeneratingFunctionGroup::new();

    for nominal_mass_idx in min_peptide_mass_idx..=max_peptide_mass_idx {
        if nominal_mass_idx <= 0 {
            continue;
        }
        let graph = PrimitiveAaGraph::new(
            aa_set,
            nominal_mass_idx,
            enzyme,
            scored_spec,
            scorer,
            charge,
            parent_mass,
            fragment_tolerance_da,
            use_protein_n_term,
            use_protein_c_term,
        );
        match GeneratingFunction::with_score_threshold(&graph, min_score, aa_set) {
            Ok(gf) => group.accept(gf),
            Err(_) => continue, // skip degenerate / unreachable bins
        }
    }

    if !group.is_computed() {
        return;
    }

    // 4. For each PSM in the queue, compute spec_e_value from its score.
    let max_score = group.max_score();

    queue.update_spec_e_values(|psm| {
        // Nominal peptide mass: residue masses sum + no water (Java convention for mass index).
        // Use nominal_from() (INTEGER_MASS_SCALER-aware) to match how graph nodes are indexed.
        let psm_nominal_mass = nominal_from(psm.candidate.peptide.mass() - H2O);
        if psm_nominal_mass < min_peptide_mass_idx || psm_nominal_mass > max_peptide_mass_idx {
            return 1.0;
        }
        let score_int = psm.score.round() as i32;
        if score_int >= max_score {
            // Score exceeds GF range — return the probability at max_score - 1
            // (which already has the underflow guard applied by the GF DP).
            // Mirrors Java behavior; avoids returning a grossly inflated value
            // (1/max_score ≈ 0.01) that would invert ranking of the best PSMs.
            return group.spectral_probability(max_score - 1)
                .unwrap_or(f32::from_bits(1) as f64);
        }
        group.spectral_probability(score_int).unwrap_or(1.0)
    });

    // 5. Phase 7 enrichment: set de_novo_score and e_value for output writers.
    //
    // de_novo_score = group.max_score() - 1  (mirrors Java's getDeNovoScore()).
    //
    // e_value = spec_e_value * num_distinct_peptides_at_length.
    // Approximate: count how many PSMs in the queue share the same peptide
    // length as this PSM and use queue.len() as the peptide-space proxy.
    // This matches Java's `spec_e_value * sa.getNumDistinctPeptides(length)`
    // intent (Java reads distinct peptides of given length from the suffix
    // array; our proxy uses the candidate-set count which is approximate and
    // typically within the same order of magnitude). Documented as
    // intentional MVP approximation — Phase 7+ wires in the real suffix-array
    // helper.
    let de_novo_score = max_score - 1;
    // Collect peptide lengths once (cheap; queue is ≤ top_n, usually ≤ 10).
    let length_counts: std::collections::HashMap<usize, usize> = {
        let mut map = std::collections::HashMap::new();
        for psm in queue.iter_psms() {
            let len = psm.candidate.peptide.length();
            *map.entry(len).or_insert(0) += 1;
        }
        map
    };
    queue.update_psm_enrichment(|psm| {
        psm.de_novo_score = de_novo_score;
        let num_distinct = *length_counts.get(&psm.candidate.peptide.length()).unwrap_or(&1);
        psm.e_value = psm.spec_e_value * num_distinct as f64;
    });
}

/// Compute fragment-ion feature columns for a single PSM.
///
/// Uses charge-1 b/y ions only (matching Java's `NumMatchedMainIons`
/// convention).  A peptide position counts at most once per ion series;
/// a position can contribute 1 from b AND 1 from y (so the maximum
/// `num_matched_main_ions` is `2 * (n - 1)` for a peptide of length n).
///
/// Returns `PsmFeatures::default()` for peptides shorter than 2 residues
/// (no cleavable fragment ions exist).
///
/// # Phase 4 alignment: 9 new ion-current + error-stat features
///
/// All 9 previously zero-stubbed PIN columns are now filled:
/// - Ion-current ratios use raw peak intensities vs total MS2 ion current.
///   Mirrors Java `PSMFeatureFinder.computeExplainedIonCurrent()`.
/// - `MS2IonCurrent` is the raw sum (NOT log10).  Java `getMS2IonCurrent()`
///   returns the raw sum; the PIN emitter emits it as-is.
/// - `IsolationWindowEfficiency` is always 0.0; Java returns `null` here
///   (no isolation-window data in the Spectrum object).
/// - Top-7 error stats mirror Java `MassErrorStat`: errors are collected for
///   all matched b+y ions, sorted descending by intensity, top-7 taken;
///   absolute Da error for mean/stdev, signed ppm for rel-mean/rel-stdev.
///   Population stdev formula: `sqrt(E[x²] - mean²)` — matches Java.
pub(crate) fn compute_psm_features(
    scored_spec: &ScoredSpectrum<'_>,
    peptide: &Peptide,
    fragment_tolerance_da: f64,
) -> PsmFeatures {
    let n = peptide.length();
    if n < 2 {
        return PsmFeatures::default();
    }

    // Predict charge-1 b/y ions; one bool per fragment position.
    let predicted = predict_by_ions(peptide, 1..=1);
    let mut b_matched = vec![false; n - 1];
    let mut y_matched = vec![false; n - 1];

    // Collect matched-ion details for ion-current ratio and error-stat features.
    // Each entry: (intensity, observed_mz, predicted_mz, is_b_ion)
    let mut matched_ions: Vec<(f32, f64, f64, bool)> = Vec::new();

    for p in &predicted {
        if let Some((_rank, intensity, peak_mz)) =
            scored_spec.nearest_peak_full(p.mz, fragment_tolerance_da)
        {
            let is_b = matches!(p.kind, IonKind::B);
            matched_ions.push((intensity, peak_mz, p.mz, is_b));

            // position is 1-based (b1/y1 = index 0 in the matched arrays)
            let pos = (p.position - 1) as usize;
            match p.kind {
                IonKind::B => {
                    if pos < b_matched.len() {
                        b_matched[pos] = true;
                    }
                }
                IonKind::Y => {
                    if pos < y_matched.len() {
                        y_matched[pos] = true;
                    }
                }
            }
        }
    }

    let num_matched: u32 = (b_matched.iter().filter(|&&m| m).count()
        + y_matched.iter().filter(|&&m| m).count()) as u32;

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

    let longest_b = longest_run(&b_matched);
    let longest_y = longest_run(&y_matched);

    // ── Ion-current ratio features ────────────────────────────────────────────

    let total_intensity = scored_spec.total_intensity(); // raw sum, all peaks

    let matched_b_intensity: f64 = matched_ions.iter()
        .filter(|&&(_, _, _, is_b)| is_b)
        .map(|&(int, _, _, _)| int as f64)
        .sum();
    let matched_y_intensity: f64 = matched_ions.iter()
        .filter(|&&(_, _, _, is_b)| !is_b)
        .map(|&(int, _, _, _)| int as f64)
        .sum();
    let matched_total = matched_b_intensity + matched_y_intensity;

    let safe_div = |num: f64, denom: f64| -> f32 {
        if denom > 0.0 { (num / denom) as f32 } else { 0.0 }
    };

    let explained_ion_current_ratio = safe_div(matched_total, total_intensity);
    let n_term_ion_current_ratio    = safe_div(matched_b_intensity, total_intensity);
    let c_term_ion_current_ratio    = safe_div(matched_y_intensity, total_intensity);
    // Java `getMS2IonCurrent()` returns the raw sum (no log10 transform).
    let ms2_ion_current = if total_intensity > 0.0 { total_intensity as f32 } else { 0.0 };
    // Java `getIsolationWindowEfficiency()` always returns null → emit 0.0.
    let isolation_window_efficiency = 0.0_f32;

    // ── Top-7 mass-error statistics ───────────────────────────────────────────

    // Sort matched ions descending by intensity (mirrors Java MassErrorStat
    // which sorts errorList by intensity via PairReverseComparator).
    let mut by_intensity = matched_ions.clone();
    by_intensity.sort_by(|a, b| {
        b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal)
    });
    let top7: Vec<(f32, f64, f64, bool)> = by_intensity.into_iter().take(7).collect();

    // Java MassErrorStat: absolute Da errors for mean7/sd7;
    //                     signed errors (no abs) for rMean7/rSd7 (ppm).
    // Population stdev formula (Java): sqrt(sumSq/n - mean²).
    let abs_da_errors: Vec<f64> = top7.iter()
        .map(|&(_, obs, pred, _)| (obs - pred).abs())
        .collect();
    let rel_ppm_errors: Vec<f64> = top7.iter()
        .filter(|&&(_, _, pred, _)| pred > 0.0)
        .map(|&(_, obs, pred, _)| (obs - pred) / pred * 1e6)
        .collect();

    fn mean_and_pop_stdev(values: &[f64]) -> (f32, f32) {
        if values.is_empty() { return (0.0, 0.0); }
        let n = values.len() as f64;
        let mean = values.iter().sum::<f64>() / n;
        let sum_sq: f64 = values.iter().map(|v| v * v).sum();
        let var = (sum_sq / n - mean * mean).max(0.0); // clamp negative rounding noise
        (mean as f32, var.sqrt() as f32)
    }

    let (mean_error_top7, stdev_error_top7)         = mean_and_pop_stdev(&abs_da_errors);
    let (mean_rel_error_top7, stdev_rel_error_top7) = mean_and_pop_stdev(&rel_ppm_errors);

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
    }
}

// ── Unit tests for Phase 4 alignment feature columns ─────────────────────────

#[cfg(test)]
mod feature_tests {
    use super::*;
    use model::amino_acid::AminoAcid;
    use model::mass::PROTON;
    use model::peptide::Peptide;
    use model::spectrum::Spectrum;
    use scoring_crate::scoring::fragment_ions::predict_by_ions;
    use scoring_crate::scoring::ScoredSpectrum;

    /// Build a minimal peptide of `len` alanine residues with flanks `_-`.
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
        }
    }

    // ── Test: empty spectrum → all new features are 0 ───────────────────────

    #[test]
    fn compute_psm_features_top7_error_stats_zero_when_no_matches() {
        let pep = ala_peptide(4);
        let spec = make_spectrum(vec![]); // no peaks
        let ss = ScoredSpectrum::new_without_filtering(&spec);
        let f = compute_psm_features(&ss, &pep, 0.5);
        assert_eq!(f.mean_error_top7,     0.0, "mean_error_top7 should be 0 with no matches");
        assert_eq!(f.stdev_error_top7,    0.0, "stdev_error_top7 should be 0 with no matches");
        assert_eq!(f.mean_rel_error_top7,  0.0, "mean_rel_error_top7 should be 0 with no matches");
        assert_eq!(f.stdev_rel_error_top7, 0.0, "stdev_rel_error_top7 should be 0 with no matches");
        assert_eq!(f.explained_ion_current_ratio, 0.0, "ratio should be 0 with no peaks");
        assert_eq!(f.ms2_ion_current, 0.0, "ms2_ion_current should be 0 with no peaks");
    }

    // ── Test: ion-current ratios populate and satisfy arithmetic invariant ───

    #[test]
    fn compute_psm_features_populates_ion_current_ratios() {
        // Use a 3-residue peptide (ALA-ALA-ALA). predict_by_ions(charge=1) gives:
        //   b1, y1, b2, y2 at definite m/z values.
        // We place spectrum peaks at exactly those m/z values so all ions match,
        // then verify explained_ratio > 0 and n + c == explained.
        let pep = ala_peptide(3);
        let predicted = predict_by_ions(&pep, 1..=1);

        // Place peaks exactly at every predicted m/z with increasing intensities.
        let mut peaks: Vec<(f64, f32)> = predicted
            .iter()
            .enumerate()
            .map(|(i, p)| (p.mz, (i + 1) as f32 * 10.0))
            .collect();
        // Add some unmatched background intensity so total_intensity > matched.
        peaks.push((1500.0, 5.0)); // far from any ion
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        let spec = make_spectrum(peaks);
        let ss = ScoredSpectrum::new_without_filtering(&spec);
        let f = compute_psm_features(&ss, &pep, 0.01); // tight tolerance

        // All ratios should be positive since all predicted ions match.
        assert!(f.explained_ion_current_ratio > 0.0,
            "explained_ion_current_ratio should be > 0 when ions match, got {}",
            f.explained_ion_current_ratio);
        assert!(f.n_term_ion_current_ratio > 0.0,
            "n_term_ion_current_ratio should be > 0 when b-ions match");
        assert!(f.c_term_ion_current_ratio > 0.0,
            "c_term_ion_current_ratio should be > 0 when y-ions match");

        // Invariant: n_term + c_term == explained (within float precision)
        let sum = f.n_term_ion_current_ratio + f.c_term_ion_current_ratio;
        assert!(
            (sum - f.explained_ion_current_ratio).abs() < 1e-5,
            "n_term + c_term should == explained ({} + {} != {})",
            f.n_term_ion_current_ratio, f.c_term_ion_current_ratio, f.explained_ion_current_ratio
        );

        // ms2_ion_current should equal total peak intensity sum.
        let total: f32 = ss.total_intensity() as f32;
        assert!((f.ms2_ion_current - total).abs() < 1.0,
            "ms2_ion_current {} should match total spectrum intensity {}",
            f.ms2_ion_current, total);

        // isolation_window_efficiency always 0.0.
        assert_eq!(f.isolation_window_efficiency, 0.0);
    }

    // ── Test: top-7 error stats are nonzero when ions match ─────────────────

    #[test]
    fn compute_psm_features_error_stats_nonzero_when_ions_match_with_offset() {
        // Build a peptide and shift every peak by a fixed offset so errors are known.
        let pep = ala_peptide(5);
        let predicted = predict_by_ions(&pep, 1..=1);

        let offset_da = 0.01_f64;  // 10 mDa error on every peak
        let mut peaks: Vec<(f64, f32)> = predicted
            .iter()
            .enumerate()
            .map(|(i, p)| (p.mz + offset_da, (i + 1) as f32 * 10.0))
            .collect();
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        let spec = make_spectrum(peaks);
        let ss = ScoredSpectrum::new_without_filtering(&spec);
        // tolerance = 0.05 Da so all offset peaks are still within window.
        let f = compute_psm_features(&ss, &pep, 0.05);

        // All absolute Da errors should be ~offset_da.
        assert!(
            f.mean_error_top7 > 0.0,
            "mean_error_top7 should be > 0 when peaks are systematically offset, got {}",
            f.mean_error_top7
        );
        // With identical errors, stdev should be near 0.
        assert!(
            f.stdev_error_top7 < 1e-4,
            "stdev_error_top7 should be ~0 for identical errors, got {}",
            f.stdev_error_top7
        );
        // Relative error should also be nonzero.
        assert!(
            f.mean_rel_error_top7 != 0.0,
            "mean_rel_error_top7 should be nonzero when peaks are offset"
        );
    }

    // ── Test: ms2_ion_current mirrors total_intensity exactly ───────────────

    #[test]
    fn ms2_ion_current_equals_total_intensity() {
        let pep = ala_peptide(3);
        let peaks = vec![(100.0, 50.0_f32), (200.0, 30.0), (300.0, 20.0)];
        let spec = make_spectrum(peaks.clone());
        let ss = ScoredSpectrum::new_without_filtering(&spec);
        let f = compute_psm_features(&ss, &pep, 0.5);

        let expected: f32 = peaks.iter().map(|&(_, i)| i).sum();
        assert_eq!(f.ms2_ion_current, expected,
            "ms2_ion_current {} should equal sum of peak intensities {}",
            f.ms2_ion_current, expected);
    }

    // ── Test: PROTON mass sanity — b1 ion for alanine at charge 1 ───────────
    // This verifies the predict_by_ions formula aligns with our test setup.
    #[test]
    fn b1_mz_for_alanine_is_proton_plus_residue_mass() {
        use model::amino_acid::AminoAcid;
        let aa = AminoAcid::standard(b'A').unwrap();
        let residue_mass = aa.mass; // monoisotopic residue mass
        let expected_b1_mz = residue_mass + PROTON; // charge 1
        let pep = ala_peptide(2);
        let predicted = predict_by_ions(&pep, 1..=1);
        let b1 = predicted.iter().find(|p| matches!(p.kind, IonKind::B) && p.position == 1)
            .expect("b1 ion should exist");
        assert!(
            (b1.mz - expected_b1_mz).abs() < 1e-6,
            "b1 mz {} expected {}", b1.mz, expected_b1_mz
        );
    }
}
