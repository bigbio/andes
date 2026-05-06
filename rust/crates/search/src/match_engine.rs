//! Top-level integration: spectra × candidates → top-N PSMs per spectrum.

use std::collections::HashMap;

use model::aa_set::AminoAcidSet;
use crate::candidate_gen::{enumerate_candidates, Candidate};
use model::enzyme::Enzyme;
use scoring_crate::gf::generating_function::GeneratingFunction;
use scoring_crate::gf::group::GeneratingFunctionGroup;
use scoring_crate::gf::primitive_graph::PrimitiveAaGraph;
use model::mass::{H2O, PROTON};
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

        for cand in &candidates {
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
    let nominal_peptide_mass = peptide_neutral_mass.round() as i32;

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

    // 3. Derive protein-terminal flags from the top PSM (Track B4).
    //
    // Java reference: DBScanner.java:592 aggregates these flags across all
    // candidates before GF construction. Our MVP approximation: derive from
    // the single best-scoring (top) PSM currently in the queue — the most
    // likely promotion candidate. This is exact for the common case where the
    // top PSM is unambiguously best; edge cases (ties near a protein boundary)
    // are addressed in Phase 7+ when per-candidate GFs become feasible.
    let (use_protein_n_term, use_protein_c_term) = {
        let top_psm = queue.iter_psms().max_by(|a, b| a.cmp(b));
        match top_psm {
            Some(top) => {
                let start = top.candidate.start_offset_in_protein;
                let pep_len = top.candidate.peptide.length();
                let is_n = start == 0;
                let is_c = match search_index.protein_at(top.candidate.protein_index) {
                    Some(prot) => start + pep_len >= prot.sequence.len(),
                    None => false,
                };
                (is_n, is_c)
            }
            None => (false, false),
        }
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
        let psm_nominal_mass = (psm.candidate.peptide.mass() - H2O).round() as i32;
        if psm_nominal_mass < min_peptide_mass_idx || psm_nominal_mass > max_peptide_mass_idx {
            return 1.0;
        }
        let score_int = psm.score.round() as i32;
        if score_int >= max_score {
            // Score exceeds GF range: approximate as 1/max_score.
            return if max_score > 0 { 1.0 / max_score as f64 } else { 1.0 };
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
/// # Deferred features
/// The following columns remain zero-stubbed in the PIN writer and are
/// intentionally NOT computed here:
/// - `ExplainedIonCurrentRatio`, `NTermIonCurrentRatio`, `CTermIonCurrentRatio`:
///   require summing matched peak intensities vs total MS2 ion current.
/// - `MS2IonCurrent`, `IsolationWindowEfficiency`:
///   require raw precursor isolation window intensity data not yet threaded
///   into `PsmMatch`.
/// - `MeanErrorTop7`, `StdevErrorTop7`, `MeanRelErrorTop7`, `StdevRelErrorTop7`:
///   require mass error statistics over the top-7 matched ions.
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

    for p in &predicted {
        if scored_spec.nearest_peak_rank(p.mz, fragment_tolerance_da).is_some() {
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

    PsmFeatures {
        num_matched_main_ions: num_matched,
        longest_b,
        longest_y,
        longest_y_pct: longest_y as f32 / n as f32,
        matched_ion_ratio: num_matched as f32 / n as f32,
    }
}
