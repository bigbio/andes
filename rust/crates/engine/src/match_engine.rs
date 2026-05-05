//! Top-level integration: spectra × candidates → top-N PSMs per spectrum.

use crate::candidate_gen::{enumerate_candidates, Candidate};
use crate::precursor_matching::{matches_precursor, MassError};
use crate::psm::{PsmMatch, TopNQueue};
use crate::search_index::SearchIndex;
use crate::search_params::SearchParams;
use crate::spectrum::Spectrum;

/// Match every spectrum against every candidate from the SearchIndex.
/// Returns one top-N PSM queue per spectrum, in input order.
///
/// Phase 4e MVP: O(spectra × candidates) — acceptable for small DBs.
/// Phase 4f optimization: bucket candidates by mass for sub-linear lookup.
pub fn match_spectra(
    spectra: &[Spectrum],
    idx: &SearchIndex,
    params: &SearchParams,
    decoy_prefix: &str,
) -> Vec<TopNQueue> {
    let mut queues: Vec<TopNQueue> = (0..spectra.len())
        .map(|_| TopNQueue::new(params.top_n_psms_per_spectrum))
        .collect();

    let candidates: Vec<Candidate> = enumerate_candidates(idx, params, decoy_prefix).collect();

    for (spec_idx, spec) in spectra.iter().enumerate() {
        let charges_to_try: Vec<u8> = match spec.precursor_charge {
            Some(z) if z > 0 => vec![z as u8],
            _ => params.charge_range.clone().collect(),
        };

        for cand in &candidates {
            for &z in &charges_to_try {
                let mut best_for_charge: Option<(MassError, f32)> = None;
                for offset in params.isotope_error_range.clone() {
                    if let Some(err) = matches_precursor(spec, &cand.peptide, z, offset, &params.precursor_tolerance) {
                        let score = -(err.mass_error_ppm.abs() as f32);
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
                    });
                }
            }
        }
    }

    queues
}
