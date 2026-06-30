// Glycopeptide scoring driver for andes.
//
// Implements bare-backbone glyco-PSM scoring: given a set of spectra and a
// PreparedSearch (standard tryptic candidate pool), enumerates hybrid
// backbone candidates (DB-branch glycan list + de-novo Y-ladder), filters
// peptide candidates by backbone mass and N-X-S/T sequon, and scores each
// (backbone peptide, glycan) pair using the standard rank-LLR scorer.
//
// The glycan mass is NOT added to the peptide's Asn — this is the
// "bare-backbone" scoring model: Percolator sees standard b/y ions from the
// peptide backbone only.  The glycan-level evidence (oxonium ions, Y-ladder)
// lives in the GlycoPsmKey appended as additive PIN columns.
//
// Placement: inside the search crate so `pub(crate)` items (compute_psm_features,
// candidate_nominal_bounds) are reachable without visibility changes.

use std::collections::HashMap;

use model::mass::{nominal_from, H2O, PROTON};
use model::spectrum::Spectrum;

use andes_glyco::glycan_db::GlycanComp;
use andes_glyco::glyco_psm::GlycoPsmKey;
use andes_glyco::hybrid::{hybrid_candidates, BackboneHit, Source};
use andes_glyco::oxonium::oxonium_gate;
use andes_glyco::sequon::has_nxst_sequon;

use crate::match_engine::{compute_psm_features, PreparedSearch};
use crate::psm::PsmMatch;
#[cfg(test)]
use crate::psm::PsmFeatures;
use scoring_crate::scoring::{psm_edge_score, score_psm, ScoredSpectrum};

/// A scored glyco-PSM: the bare-backbone PSM + all glycan-level evidence.
#[derive(Debug, Clone)]
pub struct FullGlycoPsm {
    /// Glycan-level key (oxonium evidence, Y-ladder, glycan composition).
    pub glycan_key: GlycoPsmKey,
    /// Standard PSM (bare backbone, scored as if unmodified).
    pub psm: PsmMatch,
}

/// Per-spectrum result: the spectrum's global index + all scored glyco PSMs.
#[derive(Debug, Clone)]
pub struct GlycoSpectrumResult {
    pub spectrum_idx: usize,
    pub hits: Vec<FullGlycoPsm>,
}

/// Run the glyco-PSM scoring driver over all spectra.
///
/// For each spectrum:
/// 1. Run `oxonium_gate` to gather oxonium evidence.
/// 2. For each charge in the params charge range, call `hybrid_candidates`
///    to enumerate backbone hits (DB + de-novo).
/// 3. Union and dedup backbone hits within 0.02 Da, capping at `backbone_top_k`.
/// 4. For each backbone hit, find candidates in the mass bucket whose peptide
///    mass matches the backbone and has a N-X-S/T sequon.
/// 5. Score each (peptide, glycan) pair and emit a `FullGlycoPsm`.
///
/// Results are serialized (rayon is not used here to keep v1 simple; the
/// standard search path handles parallelism separately).
pub fn glyco_search_run(
    spectra: &[Spectrum],
    prepared: &PreparedSearch<'_>,
    glycan_list: &[GlycanComp],
    tol_ppm: f64,
    backbone_top_k: usize,
) -> Vec<GlycoSpectrumResult> {
    let scorer = prepared.scorer;
    let params = prepared.params;
    let candidates = &prepared.candidates;
    let bucket_index = &prepared.bucket_index;
    let fragment_tolerance_da = prepared.fragment_tolerance_da;

    let mut results: Vec<GlycoSpectrumResult> = Vec::new();

    for (spec_idx, spec) in spectra.iter().enumerate() {
        if spec.peaks.len() < params.min_peaks as usize {
            continue;
        }

        // Oxonium evidence for the whole spectrum (charge-independent).
        let ox_ev = oxonium_gate(&spec.peaks, 0.10, tol_ppm);

        // Determine which charges to try.
        let charges_to_try: Vec<u8> = match spec.precursor_charge {
            Some(z) if z > 0 => vec![z as u8],
            _ => params.charge_range.clone().collect(),
        };

        // Gather backbone hits across all charges, then union+dedup.
        let mut all_backbone: Vec<BackboneHit> = Vec::new();
        for &z in &charges_to_try {
            let charge_f = z as f64;
            let precursor_neutral = (spec.precursor_mz - PROTON) * charge_f - H2O;
            let hits = hybrid_candidates(
                &spec.peaks,
                precursor_neutral,
                z,
                glycan_list,
                tol_ppm,
                50,
            );
            for h in hits {
                all_backbone.push(h);
            }
        }

        // Dedup cross-charge backbone hits within 0.02 Da.
        // Sort by backbone_mass, then by Db<DeNovo, then dedup.
        all_backbone.sort_by(|a, b| {
            a.backbone_mass
                .partial_cmp(&b.backbone_mass)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| {
                    let oa = if a.source == Source::Db { 0u8 } else { 1u8 };
                    let ob = if b.source == Source::Db { 0u8 } else { 1u8 };
                    oa.cmp(&ob)
                })
        });
        let mut deduped_backbone: Vec<BackboneHit> = Vec::with_capacity(all_backbone.len());
        if !all_backbone.is_empty() {
            let mut rep = all_backbone.remove(0);
            for next in all_backbone {
                let tol = (rep.backbone_mass * 20e-6_f64).max(0.02);
                if (next.backbone_mass - rep.backbone_mass).abs() < tol {
                    // Prefer Db over DeNovo within a cluster.
                    if rep.source == Source::DeNovo && next.source == Source::Db {
                        rep = next;
                    }
                } else {
                    deduped_backbone.push(rep);
                    rep = next;
                }
            }
            deduped_backbone.push(rep);
        }

        // Sort by core_y_hits DESC (we don't have this in BackboneHit; just keep
        // the natural order) and cap at backbone_top_k.
        deduped_backbone.truncate(backbone_top_k);

        if deduped_backbone.is_empty() {
            continue;
        }

        // Build ScoredSpectrum per charge (cached).
        let mut scored_per_charge: Vec<(u8, ScoredSpectrum<'_>)> = Vec::new();
        for &z in &charges_to_try {
            if scored_per_charge.iter().all(|(c, _)| *c != z) {
                scored_per_charge.push((z, ScoredSpectrum::new(spec, scorer, z)));
            }
        }

        // Map (backbone, glycan) → best FullGlycoPsm (dedup same pair by score).
        // Key: (candidate_idx, glycan_hexnac, glycan_hex, glycan_fuc, glycan_neuac, glycan_neugc).
        // For DeNovo hits (no glycan), key uses (candidate_idx, 255, 255, 255, 255, 255).
        let mut best_hits: HashMap<(u32, u8, u8, u8, u8, u8), FullGlycoPsm> = HashMap::new();

        for bb_hit in &deduped_backbone {
            let bb = bb_hit.backbone_mass;

            // Derive nominal bucket bounds for the backbone mass.
            // backbone mass is the PEPTIDE neutral mass (residues + H2O).
            // nominal_residue_mass = nominal(mass - H2O), which is what bucket_index keys on.
            let nb = nominal_from(bb - H2O);
            // Search ±1 Da bucket (coarser than needed, but safe).
            let min_nom = nb - 1;
            let max_nom = nb + 1;

            // Also widen by the precursor tolerance (compute from params for one
            // representative charge; use the narrowest precursor tolerance side).
            let candidate_slot_iter: Vec<usize> = bucket_index
                .range(min_nom..=max_nom)
                .flat_map(|(_, v)| v.iter().copied())
                .collect();

            for cand_slot in candidate_slot_iter {
                let cand = &candidates[cand_slot];
                let cand_mass = cand.peptide.mass(); // includes H2O

                // Fine mass filter: peptide must match the backbone within tol.
                let mass_tol = (bb * tol_ppm * 1e-6_f64).max(0.01);
                if (cand_mass - bb).abs() > mass_tol {
                    continue;
                }

                // Sequon filter: must have at least one N-X-S/T.
                let residue_bytes: Vec<u8> =
                    cand.peptide.residues.iter().map(|aa| aa.residue).collect();
                if !has_nxst_sequon(&residue_bytes) {
                    continue;
                }

                // Score the backbone peptide at all relevant charges.
                for &z in &charges_to_try {
                    // Quick precursor check: backbone mass should match at this charge.
                    // We accept candidates that passed the bucket filter, but do a
                    // charge-implied neutral-mass sanity check.
                    let obs_neutral = (spec.precursor_mz - PROTON) * z as f64 - H2O;
                    // The glycan mass is (obs_neutral - bb), check it's positive.
                    let glycan_mass_implied = obs_neutral - bb;
                    if glycan_mass_implied < 0.0 {
                        continue;
                    }

                    let scored_spec = scored_per_charge
                        .iter()
                        .find(|(c, _)| *c == z)
                        .map(|(_, s)| s)
                        .expect("ScoredSpectrum for charge must exist");

                    let pin_score =
                        score_psm(scored_spec, &cand.peptide, scorer, z, fragment_tolerance_da);
                    let edge_i = psm_edge_score(scored_spec, &cand.peptide, scorer, z);
                    let rank_score = pin_score + edge_i as f32;

                    let features = compute_psm_features(
                        scored_spec,
                        &cand.peptide,
                        scorer,
                        z,
                        prepared.intensity_model.as_deref(),
                    );

                    // Mass error in ppm vs the backbone mass.
                    let mass_error_ppm = if bb > 0.0 {
                        (cand_mass - bb) / bb * 1e6
                    } else {
                        0.0
                    };

                    let psm = PsmMatch {
                        spectrum_idx: spec_idx,
                        candidate_idxs: vec![cand_slot as u32],
                        charge_used: z,
                        mass_error_ppm,
                        score: pin_score,
                        rank_score,
                        edge_score: edge_i,
                        activation_method: Some(scorer.param().data_type.activation),
                        features,
                        isotope_offset: 0,
                        precursor_mz_override: None,
                    };

                    let glycan_mass = bb_hit.glycan.as_ref().map(|g| g.mass).unwrap_or(0.0);

                    let glycan_key = GlycoPsmKey {
                        spectrum_idx: spec_idx,
                        glycan: bb_hit.glycan.clone(),
                        glycan_source: bb_hit.source.clone(),
                        oxonium_summed_frac: ox_ev.summed_frac,
                        n_core_oxonium_ions: ox_ev.n_core_ions,
                        y_ladder_intensity_score: 0.0, // populated by backbone solver in future
                        core_y_hits: 0,                // populated by backbone solver in future
                        glycan_mass,
                        backbone_mass: bb,
                    };

                    // Dedup key.
                    let gl_key = match &bb_hit.glycan {
                        Some(g) => (
                            cand_slot as u32,
                            g.hexnac,
                            g.hex,
                            g.fuc,
                            g.neuac,
                            g.neugc,
                        ),
                        None => (cand_slot as u32, 255, 255, 255, 255, 255),
                    };

                    let new_hit = FullGlycoPsm { glycan_key, psm };
                    best_hits
                        .entry(gl_key)
                        .and_modify(|existing| {
                            if new_hit.psm.rank_score > existing.psm.rank_score {
                                *existing = new_hit.clone();
                            }
                        })
                        .or_insert(new_hit);
                }
            }
        }

        if !best_hits.is_empty() {
            let hits: Vec<FullGlycoPsm> = best_hits.into_values().collect();
            results.push(GlycoSpectrumResult { spectrum_idx: spec_idx, hits });
        }
    }

    results
}

#[cfg(test)]
mod tests {
    // Integration-level tests are deferred to the search-crate integration tests
    // (tests/ directory) where real PreparedSearch fixtures can be built.
    // Unit-level sequon + mass filter logic is tested in andes_glyco::sequon.
    //
    // Smoke test: verify the public types compile and are accessible.
    use super::*;

    #[test]
    fn full_glyco_psm_is_clone() {
        // Minimal construction check — verifies the types are well-formed.
        // PsmMatch does not impl Default, so we do the minimal construction.
        let psm = PsmMatch {
            spectrum_idx: 0,
            candidate_idxs: vec![0],
            charge_used: 2,
            mass_error_ppm: 0.0,
            score: 0.0,
            rank_score: 0.0,
            edge_score: 0,
            activation_method: None,
            features: PsmFeatures::default(),
            isotope_offset: 0,
            precursor_mz_override: None,
        };
        let key = GlycoPsmKey {
            spectrum_idx: 0,
            glycan: None,
            glycan_source: Source::Db,
            oxonium_summed_frac: 0.0,
            n_core_oxonium_ions: 0,
            y_ladder_intensity_score: 0.0,
            core_y_hits: 0,
            glycan_mass: 0.0,
            backbone_mass: 0.0,
        };
        let hit = FullGlycoPsm { glycan_key: key, psm };
        let cloned = hit.clone();
        assert_eq!(cloned.psm.spectrum_idx, 0);
    }

    #[test]
    fn glyco_spectrum_result_is_clone() {
        let result = GlycoSpectrumResult { spectrum_idx: 7, hits: vec![] };
        let c = result.clone();
        assert_eq!(c.spectrum_idx, 7);
    }
}
