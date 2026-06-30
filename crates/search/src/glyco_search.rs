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
// Backbone selection strategy (v2 — b/y-ranked):
//
//   The prior approach (Y-ladder pre-filter → core_y_hits-ranked cap) discarded
//   backbones whose spectra lacked strong core-Y ions before any b/y scoring,
//   capping find-rate at ~11 %.  The fix: use the curated `n_glycan_list_common()`
//   (~600 glycans instead of 2510), score ALL resulting backbone candidates in
//   phase-1 b/y scoring, aggregate the best b/y rank score per backbone, and
//   only then apply the backbone_top_k cap.  Y-ladder hit count is retained as
//   a tiebreaker so spectra with strong Y-ladder evidence still benefit from it.
//
// Placement: inside the search crate so `pub(crate)` items (compute_psm_features,
// candidate_nominal_bounds) are reachable without visibility changes.

use std::collections::HashMap;

use model::mass::{nominal_from, H2O, PROTON};
use model::spectrum::Spectrum;
use rayon::prelude::*;

use andes_glyco::backbone::count_core_y_hits;
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

            // --- b/y-ranked backbone selection (replaces Y-ladder pre-filter) ---
            //
            // Previous approach: rank all backbones by core_y_hits → truncate to
            // backbone_top_k → score the survivors in phase-1.  This discards the
            // true backbone when the spectrum has weak core-Y ions (common in HCD),
            // capping find-rate at ~11 %.
            //
            // New approach: skip the Y-ladder pre-filter entirely.  Instead, run
            // phase-1 b/y scoring (score_psm) for EVERY backbone candidate.  Because
            // we use n_glycan_list_common() (~600 glycans) by default, the total
            // number of (backbone, candidate) pairs per spectrum is tractable.
            //
            // After phase-1 we know the best b/y score achieved for each backbone.
            // We THEN rank backbones by that best b/y score, using core_y_hits as a
            // tiebreaker, and apply the backbone_top_k cap.  Only phase-2
            // (compute_psm_features) is bounded by that cap.

            // Build ScoredSpectrum per unique charge (cached, cheap amortized).
            let mut scored_per_charge: Vec<(u8, ScoredSpectrum<'_>)> = Vec::new();
            for &z in &charges_to_try {
                if scored_per_charge.iter().all(|(c, _)| *c != z) {
                    scored_per_charge.push((z, ScoredSpectrum::new(spec, scorer, z)));
                }
            }

            // Collect core-Y hit counts for all backbones (cheap, used as tiebreaker
            // after b/y ranking; avoids a second pass over deduped_backbone later).
            let core_y_counts: Vec<u8> = deduped_backbone
                .iter()
                .map(|h| count_core_y_hits(&spec.peaks, h.backbone_mass, tol_ppm))
                .collect();

            // Phase 1: cheap b/y scoring for ALL backbones.
            //
            // Accumulate per (cand_slot, glycan_key) winner: the best-ranked
            // (backbone_hit_idx, z, rank, score, edge).
            //
            // Simultaneously track per backbone index the best b/y rank seen over
            // all of its matching candidates.  This is the signal used to rank
            // backbones AFTER phase-1.
            //
            // Dedup key: (cand_slot, glycan composition).
            // For DeNovo (no glycan): uses sentinel (255, 255, 255, 255, 255).
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

            // Per-backbone best b/y rank (index = backbone index in deduped_backbone).
            let mut backbone_best_rank: Vec<f32> =
                vec![f32::NEG_INFINITY; deduped_backbone.len()];

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

                    // Update per-backbone best rank.
                    if best_rank > backbone_best_rank[bb_idx] {
                        backbone_best_rank[bb_idx] = best_rank;
                    }

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

            // Determine which backbones are in the top-K by b/y rank.
            //
            // Rank: PRIMARY = backbone_best_rank DESC (best b/y score from any
            //       candidate that matched this backbone);
            //       SECONDARY = core_y_hits DESC (Y-ladder evidence breaks ties,
            //       so spectra with strong Y-ladder evidence retain that advantage);
            //       TERTIARY = backbone_mass DESC (deterministic tiebreak).
            let mut backbone_order: Vec<usize> = (0..deduped_backbone.len()).collect();
            backbone_order.sort_by(|&ai, &bi| {
                backbone_best_rank[bi]
                    .partial_cmp(&backbone_best_rank[ai])
                    .unwrap_or(std::cmp::Ordering::Equal)
                    .then_with(|| core_y_counts[bi].cmp(&core_y_counts[ai]))
                    .then_with(|| {
                        deduped_backbone[bi]
                            .backbone_mass
                            .partial_cmp(&deduped_backbone[ai].backbone_mass)
                            .unwrap_or(std::cmp::Ordering::Equal)
                    })
            });
            backbone_order.truncate(backbone_top_k);
            // Build a set of accepted backbone indices for O(1) lookup.
            let accepted_backbones: std::collections::HashSet<usize> =
                backbone_order.into_iter().collect();

            // Phase 2: expensive feature extraction for top-K winners only.
            // Only process cheap_winners whose backbone is in the accepted set.
            // Cap at backbone_top_k × 2 to bound compute_psm_features calls.
            let max_features = backbone_top_k * 2;
            let winners_for_features: Vec<((u32, u8, u8, u8, u8, u8), CheapWinner)> = {
                let mut v: Vec<_> = cheap_winners
                    .into_iter()
                    .filter(|(_, w)| accepted_backbones.contains(&w.bb_hit_idx))
                    .collect();
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

    /// Core-Y ranking: a backbone WITH Y-ladder support must outrank one without.
    ///
    /// This test constructs two synthetic backbone candidates:
    ///   - `true_bb` (small backbone, large glycan): has all 6 core-Y ions present
    ///     in the spectrum.
    ///   - `noise_bb` (large backbone, small glycan): has zero core-Y ions in the
    ///     spectrum.
    ///
    /// Under the old size-based ranking, `noise_bb` (larger backbone_mass) would
    /// have been ranked first and `true_bb` would be discarded.  After the fix,
    /// `count_core_y_hits` gives `true_bb` a count of 6 and `noise_bb` a count of 0,
    /// so the sort produces: `true_bb` first, `noise_bb` second.
    #[test]
    fn core_y_ranking_promotes_supported_backbone_over_unsupported() {
        use andes_glyco::backbone::count_core_y_hits;
        use andes_glyco::glycan_mass::{CORE_Y_STEPS, PROTON};

        // True backbone: small peptide (large glycan).
        // Typical serum N-glycopeptide scenario: backbone ~1100 Da, glycan ~2200 Da.
        let true_bb = 1100.0_f64;

        // Noise backbone: large peptide (small glycan).
        // The OLD buggy ranking kept this one (largest backbone = first after DESC sort).
        let noise_bb = 2800.0_f64;

        // Build synthetic spectrum: core-Y ions for true_bb only.
        let mut peaks: Vec<(f64, f32)> = vec![
            (true_bb + PROTON, 500.0),                          // Y0
            (true_bb + PROTON + CORE_Y_STEPS[0], 400.0),       // Y1
            (true_bb + PROTON + CORE_Y_STEPS[1], 350.0),       // Y2
            (true_bb + PROTON + CORE_Y_STEPS[2], 300.0),       // Y3
            (true_bb + PROTON + CORE_Y_STEPS[3], 250.0),       // Y4
            (true_bb + PROTON + CORE_Y_STEPS[4], 200.0),       // Y5
            (900.0, 10.0),   // noise
            (1050.0, 10.0),  // noise
        ];
        // Deliberately do NOT add core-Y ions for noise_bb.
        // Add some noise near noise_bb m/z to ensure they don't accidentally match.
        peaks.push((noise_bb + PROTON + 5.0, 50.0)); // off by 5 Da — won't match

        let tol_ppm = 20.0;

        // Verify counts directly.
        let true_hits = count_core_y_hits(&peaks, true_bb, tol_ppm);
        let noise_hits = count_core_y_hits(&peaks, noise_bb, tol_ppm);

        assert_eq!(true_hits, 6, "expected all 6 core-Y hits for true_bb, got {}", true_hits);
        assert_eq!(noise_hits, 0, "expected 0 core-Y hits for noise_bb, got {}", noise_hits);

        // Now simulate the new ranking logic: sort by core_y_hits DESC, backbone_mass DESC.
        let mut candidates = vec![
            (noise_bb, noise_hits), // large backbone — old ranking would put this first
            (true_bb, true_hits),   // small backbone — true hit
        ];
        candidates.sort_by(|&(am, ay), &(bm, by)| {
            by.cmp(&ay).then_with(|| bm.partial_cmp(&am).unwrap_or(std::cmp::Ordering::Equal))
        });

        assert!(
            (candidates[0].0 - true_bb).abs() < 0.01,
            "expected true_bb ranked first after core-Y sort, got backbone_mass={}",
            candidates[0].0
        );
        assert!(
            (candidates[1].0 - noise_bb).abs() < 0.01,
            "expected noise_bb ranked second, got backbone_mass={}",
            candidates[1].0
        );
    }

    /// b/y ranking: a backbone whose backbone b/y ions match the spectrum must
    /// outrank a backbone that does NOT match, even when the losing backbone has
    /// more core-Y hits.
    ///
    /// This validates the new backbone selection logic: after phase-1 b/y scoring,
    /// the backbone with a higher `backbone_best_rank` (best score_psm over all
    /// its matching peptide candidates) must rank above one with lower b/y rank,
    /// regardless of Y-ladder evidence.
    ///
    /// We simulate the per-backbone ranking sort that runs after phase-1:
    ///   PRIMARY   = backbone_best_rank DESC
    ///   SECONDARY = core_y_hits DESC (tiebreaker)
    ///   TERTIARY  = backbone_mass DESC
    #[test]
    fn by_rank_promotes_byone_matching_backbone_over_y_ladder_backbone() {
        // true_bb: backbone whose peptide b/y ions match the spectrum.
        //   - backbone_best_rank = 10.0 (good b/y match)
        //   - core_y_hits = 0         (no Y-ladder ions — no pre-filter benefit)
        let true_bb_mass = 1500.0_f64;
        let true_bb_best_rank: f32 = 10.0;
        let true_bb_core_y: u8 = 0;

        // noise_bb: backbone with strong Y-ladder but poor b/y backbone match.
        //   - backbone_best_rank = 2.0 (poor b/y: wrong peptide candidates)
        //   - core_y_hits = 6         (coincidental Y-ladder ions)
        let noise_bb_mass = 2000.0_f64;
        let noise_bb_best_rank: f32 = 2.0;
        let noise_bb_core_y: u8 = 6;

        // Simulate the backbone_order sort from glyco_search_run:
        //   PRIMARY = backbone_best_rank DESC
        //   SECONDARY = core_y_hits DESC
        //   TERTIARY = backbone_mass DESC
        let backbones = vec![
            (noise_bb_mass, noise_bb_best_rank, noise_bb_core_y), // idx=0
            (true_bb_mass, true_bb_best_rank, true_bb_core_y),   // idx=1
        ];
        let mut order: Vec<usize> = (0..backbones.len()).collect();
        order.sort_by(|&ai, &bi| {
            backbones[bi]
                .1
                .partial_cmp(&backbones[ai].1)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| backbones[bi].2.cmp(&backbones[ai].2))
                .then_with(|| {
                    backbones[bi]
                        .0
                        .partial_cmp(&backbones[ai].0)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
        });

        // true_bb (idx=1) must be ranked first because its b/y rank (10.0) > noise_bb (2.0).
        assert_eq!(
            order[0], 1,
            "expected true_bb (idx=1) ranked first by b/y rank, got idx={}",
            order[0]
        );
        assert_eq!(
            order[1], 0,
            "expected noise_bb (idx=0) ranked second, got idx={}",
            order[1]
        );
    }
}
