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
use rayon::prelude::*;

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

    // Process spectra in parallel; filter_map returns None for spectra with no
    // glyco hits. Order within the output Vec is not guaranteed (rayon chunks),
    // but the PIN writer only needs the hits themselves (spec_idx is in each row).
    let results: Vec<GlycoSpectrumResult> = spectra
        .par_iter()
        .enumerate()
        .filter_map(|(spec_idx, spec)| {
            if spec.peaks.len() < params.min_peaks as usize {
                return None;
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

            if all_backbone.is_empty() {
                return None;
            }

            // Dedup cross-charge backbone hits within 0.02 Da.
            // Sort ascending by backbone_mass, prefer Db over DeNovo in ties.
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
            {
                let mut rep = all_backbone.remove(0);
                for next in all_backbone {
                    let tol = (rep.backbone_mass * 20e-6_f64).max(0.02);
                    if (next.backbone_mass - rep.backbone_mass).abs() < tol {
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

            // Sort DESC by backbone_mass (largest backbone = smallest glycan first),
            // then cap at backbone_top_k.  Larger backbone masses are more likely to
            // have matching candidates in the tryptic database (smaller glycan,
            // bigger peptide).
            deduped_backbone.sort_by(|a, b| {
                b.backbone_mass
                    .partial_cmp(&a.backbone_mass)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            deduped_backbone.truncate(backbone_top_k);

            // Build ScoredSpectrum per unique charge (cached, cheap amortized).
            let mut scored_per_charge: Vec<(u8, ScoredSpectrum<'_>)> = Vec::new();
            for &z in &charges_to_try {
                if scored_per_charge.iter().all(|(c, _)| *c != z) {
                    scored_per_charge.push((z, ScoredSpectrum::new(spec, scorer, z)));
                }
            }

            // Two-phase scoring:
            // Phase 1: cheap score_psm+edge for all (backbone, candidate, charge) triples.
            //          Track best (backbone_hit_idx, cand_slot, z) per dedup key.
            // Phase 2: compute_psm_features only for the unique winners.
            //
            // Dedup key: (cand_slot, glycan composition).
            // For DeNovo (no glycan): uses sentinel (255, 255, 255, 255, 255).

            // Phase 1: cheap scoring pass.
            // Accumulate winners: key → (bb_hit_idx, cand_slot, z, rank, score, edge).
            #[derive(Clone, Copy)]
            struct CheapWinner {
                bb_hit_idx: usize,
                cand_slot: usize,
                z: u8,
                rank: f32,
                score: f32,
                edge: i32,
                cand_residue_mass: f64,
            }
            let mut cheap_winners: HashMap<(u32, u8, u8, u8, u8, u8), CheapWinner> =
                HashMap::new();

            for (bb_idx, bb_hit) in deduped_backbone.iter().enumerate() {
                let bb_residue = bb_hit.backbone_mass;

                // Tight nominal bounds.
                let nb = nominal_from(bb_residue);
                let tol_da = (bb_residue * tol_ppm * 1e-6_f64).max(0.01);
                let widen = (tol_da - 0.4999_f64).max(0.0_f64).round() as i32;

                let candidate_slots: Vec<usize> = bucket_index
                    .range((nb - widen)..=(nb + widen))
                    .flat_map(|(_, v)| v.iter().copied())
                    .collect();

                for cand_slot in candidate_slots {
                    let cand = &candidates[cand_slot];
                    let cand_residue_mass = cand.peptide.mass() - H2O;

                    if (cand_residue_mass - bb_residue).abs() > tol_da {
                        continue;
                    }

                    let residue_bytes: Vec<u8> =
                        cand.peptide.residues.iter().map(|aa| aa.residue).collect();
                    if !has_nxst_sequon(&residue_bytes) {
                        continue;
                    }

                    let gl_key = match &bb_hit.glycan {
                        Some(g) => (cand_slot as u32, g.hexnac, g.hex, g.fuc, g.neuac, g.neugc),
                        None => (cand_slot as u32, 255, 255, 255, 255, 255),
                    };

                    // Pick best charge cheaply.
                    let mut best_z: Option<u8> = None;
                    let mut best_rank: f32 = f32::NEG_INFINITY;
                    let mut best_score: f32 = 0.0;
                    let mut best_edge: i32 = 0;
                    for &z in &charges_to_try {
                        let obs_rn = (spec.precursor_mz - PROTON) * z as f64 - H2O;
                        if obs_rn - bb_residue < 0.0 {
                            continue;
                        }
                        let ss = scored_per_charge
                            .iter()
                            .find(|(c, _)| *c == z)
                            .map(|(_, s)| s)
                            .expect("ScoredSpectrum must exist for this charge");
                        let sc = score_psm(ss, &cand.peptide, scorer, z, fragment_tolerance_da);
                        let ei = psm_edge_score(ss, &cand.peptide, scorer, z);
                        let rk = sc + ei as f32;
                        if rk > best_rank {
                            best_rank = rk;
                            best_z = Some(z);
                            best_score = sc;
                            best_edge = ei;
                        }
                    }
                    let z = match best_z { Some(z) => z, None => continue };

                    let w = CheapWinner {
                        bb_hit_idx: bb_idx,
                        cand_slot,
                        z,
                        rank: best_rank,
                        score: best_score,
                        edge: best_edge,
                        cand_residue_mass,
                    };
                    cheap_winners
                        .entry(gl_key)
                        .and_modify(|existing| {
                            if w.rank > existing.rank {
                                *existing = w;
                            }
                        })
                        .or_insert(w);
                }
            }

            // Phase 2: expensive feature extraction for top-K winners only.
            // Cap at backbone_top_k × 2 to bound compute_psm_features calls.
            // After Phase 1 we have potentially hundreds of unique (cand, glycan)
            // pairs; only the top-ranked ones are worth expensive feature extraction.
            let max_features = backbone_top_k * 2;
            let winners_for_features: Vec<((u32, u8, u8, u8, u8, u8), CheapWinner)> = {
                let mut v: Vec<_> = cheap_winners.into_iter().collect();
                v.sort_by(|a, b| {
                    b.1.rank
                        .partial_cmp(&a.1.rank)
                        .unwrap_or(std::cmp::Ordering::Equal)
                });
                v.truncate(max_features);
                v
            };

            let mut best_hits: HashMap<(u32, u8, u8, u8, u8, u8), FullGlycoPsm> =
                HashMap::with_capacity(winners_for_features.len());

            for (gl_key, w) in winners_for_features {
                let bb_hit = &deduped_backbone[w.bb_hit_idx];
                let bb_residue = bb_hit.backbone_mass;
                let bb_neutral = bb_residue + H2O;

                let ss = scored_per_charge
                    .iter()
                    .find(|(c, _)| *c == w.z)
                    .map(|(_, s)| s)
                    .expect("ScoredSpectrum must exist for winning charge");
                let cand = &candidates[w.cand_slot];
                let features = compute_psm_features(
                    ss,
                    &cand.peptide,
                    scorer,
                    w.z,
                    prepared.intensity_model.as_deref(),
                );

                let mass_error_ppm = if bb_residue > 0.0 {
                    (w.cand_residue_mass - bb_residue) / bb_residue * 1e6
                } else {
                    0.0
                };
                let psm = PsmMatch {
                    spectrum_idx: spec_idx,
                    candidate_idxs: vec![w.cand_slot as u32],
                    charge_used: w.z,
                    mass_error_ppm,
                    score: w.score,
                    rank_score: w.rank,
                    edge_score: w.edge,
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
                    y_ladder_intensity_score: 0.0,
                    core_y_hits: 0,
                    glycan_mass,
                    backbone_mass: bb_neutral,
                };
                best_hits.insert(gl_key, FullGlycoPsm { glycan_key, psm });
            }

            if best_hits.is_empty() {
                None
            } else {
                Some(GlycoSpectrumResult {
                    spectrum_idx: spec_idx,
                    hits: best_hits.into_values().collect(),
                })
            }
        })
        .collect();

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
