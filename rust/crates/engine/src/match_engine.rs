//! Top-level integration: spectra × candidates → top-N PSMs per spectrum.

use crate::aa_set::AminoAcidSet;
use crate::candidate_gen::{enumerate_candidates, Candidate};
use crate::enzyme::Enzyme;
use crate::gf::generating_function::GeneratingFunction;
use crate::gf::group::GeneratingFunctionGroup;
use crate::gf::primitive_graph::PrimitiveAaGraph;
use crate::mass::{H2O, PROTON};
use crate::precursor_matching::{matches_precursor, MassError};
use crate::psm::{PsmMatch, TopNQueue};
use crate::search_index::SearchIndex;
use crate::search_params::SearchParams;
use crate::scoring::{score_psm, RankScorer, ScoredSpectrum};
use crate::spectrum::Spectrum;

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
        // Build ScoredSpectrum once per spectrum; reuse across all candidates.
        // Phase 5b Task 1: pass the scorer's param and the precursor charge so
        // that precursor peaks (and their charge-reduced / neutral-loss neighbours)
        // are filtered out before peak ranking.
        let precursor_z = spec.precursor_charge.unwrap_or(2) as u8;
        let scored_spec = ScoredSpectrum::new(spec, scorer.param(), precursor_z);

        let charges_to_try: Vec<u8> = match spec.precursor_charge {
            Some(z) if z > 0 => vec![z as u8],
            _ => params.charge_range.clone().collect(),
        };

        for cand in &candidates {
            for &z in &charges_to_try {
                let mut best_for_charge: Option<(MassError, f32)> = None;
                for offset in params.isotope_error_range.clone() {
                    if let Some(err) = matches_precursor(spec, &cand.peptide, z, offset, &params.precursor_tolerance) {
                        // Phase 5: use real score_psm instead of -|mass_error_ppm| placeholder.
                        let score = score_psm(&scored_spec, &cand.peptide, scorer, z, fragment_tolerance_da);
                        if best_for_charge.as_ref().map_or(true, |(_, s)| score > *s) {
                            best_for_charge = Some((err, score));
                        }
                    }
                }
                if let Some((err, score)) = best_for_charge {
                    queues[spec_idx].push(PsmMatch {
                        spectrum_idx: spec_idx,
                        candidate: cand.clone(),
                        charge_used: z,
                        mass_error_ppm: err.mass_error_ppm,
                        score,
                        spec_e_value: 1.0,  // set by Phase 6 compute_spec_e_values_for_spectrum
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
            compute_spec_e_values_for_spectrum(
                spec,
                params,
                &mut queues[spec_idx],
                &aa_set_for_gf,
                enzyme_opt,
                scorer,
                &scored_spec,
                fragment_tolerance_da,
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
/// * `spec` — the spectrum (used for precursor m/z and charge).
/// * `params` — search params (precursor_tolerance, isotope_error_range).
/// * `queue` — the PSM queue for this spectrum (mutated in place).
/// * `aa_set` — amino acid set with enzyme already registered via `register_enzyme`.
/// * `enzyme` — the search enzyme (passed to PrimitiveAaGraph; may be None).
/// * `scorer` — Phase 5 RankScorer.
/// * `scored_spec` — Phase 5 ScoredSpectrum for this spectrum.
/// * `fragment_tolerance_da` — fragment mass tolerance in Da.
#[allow(clippy::too_many_arguments)]
fn compute_spec_e_values_for_spectrum(
    spec: &Spectrum,
    params: &SearchParams,
    queue: &mut TopNQueue,
    aa_set: &AminoAcidSet,
    enzyme: Option<Enzyme>,
    scorer: &RankScorer,
    scored_spec: &ScoredSpectrum<'_>,
    fragment_tolerance_da: f64,
) {
    // 1. Determine the peptide neutral mass and its tolerance window.
    let charge = spec.precursor_charge.unwrap_or(2) as u8;
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

    // 3. Build the GF group across the nominal mass range.
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
            false, // use_protein_n_term — defer per-candidate protein-term flags to future work
            false, // use_protein_c_term
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
}
