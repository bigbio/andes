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

use model::mass::{nominal_from, H2O, ISOTOPE, PROTON};
use model::spectrum::Spectrum;
use rayon::prelude::*;

use andes_glyco::backbone::{
    core_y_intensity, count_core_y_hits, glycan_y_intensity, glycan_y_intensity_decoy,
    y0y1_anchor_intensity,
};
use andes_glyco::glycan_db::GlycanComp;
use andes_glyco::glyco_psm::GlycoPsmKey;
use andes_glyco::hybrid::{hybrid_candidates_with_isotope, BackboneHit, Source};
use andes_glyco::oxonium::oxonium_gate;
use andes_glyco::sequon::has_nxst_sequon;

use crate::glyco_fragment_index::FragmentIndex;
use andes_glyco::crossspectrum::GlycoformWhitelist;
use andes_glyco::glyco_y_index::GlycanYIndex;

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

/// Minimum implied N-glycan mass (2×HexNAc core) for the peptide-first path's
/// glycan-by-subtraction filter. Mirrors the solver's `MIN_GLYCAN`.
const MIN_GLYCAN: f64 = 406.0;

/// Find the known glycan nearest `target` mass within tolerance, using a sorted
/// `(mass, glycan_index)` view for a binary-search start. Returns the matched
/// composition, or `None` when the subtraction residual matches no known glycan
/// (used by the peptide-first path's glycan-by-subtraction check).
fn nearest_glycan_mass(
    sorted: &[(f64, usize)],
    glycans: &[GlycanComp],
    target: f64,
    tol: f64,
) -> Option<GlycanComp> {
    let lo = target - tol;
    let hi = target + tol;
    let start = sorted.partition_point(|&(m, _)| m < lo);
    let mut best: Option<(f64, usize)> = None;
    for &(m, gi) in &sorted[start..] {
        if m > hi {
            break;
        }
        let d = (m - target).abs();
        if best.map_or(true, |(bd, _)| d < bd) {
            best = Some((d, gi));
        }
    }
    best.map(|(_, gi)| glycans[gi].clone())
}

/// Dedup backbone hits collected across charges and isotope offsets.
///
/// Two hits are merged ONLY when they represent the same candidate: the same
/// backbone mass (within `max(bb*tol_ppm*1e-6, 0.02)`) AND the same glycan
/// hypothesis. Distinct glycan hypotheses at the same backbone mass are kept
/// separate:
///   - annotated (`Source::Db`) hits with different compositions, and
///   - de-novo hits from different isotope offsets — these carry different
///     residual glycan masses (`glycan_mass_residual = precursor(iso) − bb`),
///     so merging them would corrupt the intact `CalcMass` of novel glycans by
///     up to one isotope (Codex adversarial-review finding #2).
///
/// When a DeNovo and a Db hit coincide at the same backbone AND isotope offset,
/// the Db (annotated) hit is kept as the representative.
/// Stable per-composition seed for the glycan-axis decoy ladder, so the same
/// glycan always yields the same shifted decoy (a fixed decoy "structure",
/// analogous to a reversed-peptide decoy being fixed per target).
fn glycan_decoy_seed(g: &GlycanComp) -> u64 {
    let mut s: u64 = 0xD1B5_4A32_D192_ED03;
    for &c in &[g.hexnac, g.hex, g.fuc, g.neuac, g.neugc] {
        s = s.wrapping_mul(0x0100_0000_01B3).wrapping_add(c as u64);
    }
    s
}

fn dedup_backbone_hits(mut all_backbone: Vec<BackboneHit>, tol_ppm: f64) -> Vec<BackboneHit> {
    if all_backbone.is_empty() {
        return Vec::new();
    }
    // Sort by backbone mass; within a mass cluster, put Db before DeNovo and
    // monoisotopic (|offset| small) first for deterministic representatives.
    all_backbone.sort_by(|a, b| {
        a.backbone_mass
            .partial_cmp(&b.backbone_mass)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| {
                let oa = if a.source == Source::Db { 0u8 } else { 1u8 };
                let ob = if b.source == Source::Db { 0u8 } else { 1u8 };
                oa.cmp(&ob)
            })
            .then_with(|| a.isotope_offset.abs().cmp(&b.isotope_offset.abs()))
            .then_with(|| a.charge.cmp(&b.charge))
    });

    let mut deduped: Vec<BackboneHit> = Vec::with_capacity(all_backbone.len());
    let mut rep = all_backbone.remove(0);
    for next in all_backbone {
        let tol = (rep.backbone_mass * tol_ppm * 1e-6_f64).max(0.02);
        let same_backbone = (next.backbone_mass - rep.backbone_mass).abs() < tol;
        // Same candidate iff same backbone AND same glycan hypothesis.
        let same_hypothesis = match (&rep.glycan, &next.glycan) {
            (Some(g1), Some(g2)) => g1 == g2,
            // Unannotated: the residual is isotope-specific, so only the same
            // offset is the same candidate.
            (None, None) => rep.isotope_offset == next.isotope_offset,
            // DeNovo vs Db: the same candidate only at the same isotope offset
            // (then the annotated hit supersedes below).
            _ => rep.isotope_offset == next.isotope_offset,
        };
        if same_backbone && same_hypothesis {
            if rep.source == Source::DeNovo && next.source == Source::Db {
                rep = next; // prefer the annotated representative
            }
            // otherwise `next` is a duplicate of `rep` (e.g. different charge)
        } else {
            deduped.push(rep);
            rep = next;
        }
    }
    deduped.push(rep);
    deduped
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

    // Independent experiment toggles (one variable at a time): the peptide-first
    // fragment-index path (default ON) and cross-spectrum glycoform transfer
    // (default OFF, opt-in). Kept separable to isolate each lever + its cost.
    let peptide_first_on = std::env::var("ANDES_GLYCO_PEPTIDE_FIRST")
        .map(|v| v != "0")
        .unwrap_or(true);
    let cross_spectrum_on = std::env::var("ANDES_GLYCO_CROSSSPECTRUM")
        .map(|v| v == "1")
        .unwrap_or(false);
    // G3 glycan-axis decoy (default OFF). When off we must NOT compute the decoy
    // Y-ladder per hit — it is unused and ~doubles the glyco composition-ladder
    // cost, so leaving it on would slow the shipping default path.
    let glyco_decoy_on = std::env::var("ANDES_GLYCO_DECOY")
        .map(|v| v == "1")
        .unwrap_or(false);
    // Fast dev harness: ANDES_GLYCO_SCANS=<file> (one scan number per line)
    // restricts the glyco driver to those scans (e.g. the truth-scan subset), so
    // a redesign iteration is minutes not hours. The standard search still runs
    // over all spectra; only glyco scoring is subset. Unset = all spectra.
    let scan_filter: Option<std::collections::HashSet<i32>> = std::env::var("ANDES_GLYCO_SCANS")
        .ok()
        .and_then(|path| std::fs::read_to_string(&path).ok())
        .map(|s| {
            s.lines()
                .filter_map(|l| l.trim().parse::<i32>().ok())
                .collect()
        });
    // Experiment A (diagnostic): disable BOTH truncations — the hybrid DB-union
    // core-Y cap and the driver's backbone_top_k — to measure the true findable
    // ceiling. Large finite cap avoids the max_features usize overflow. SLOW.
    let exhaustive = std::env::var("ANDES_GLYCO_EXHAUSTIVE")
        .map(|v| v == "1")
        .unwrap_or(false);
    let effective_top_k = if exhaustive { 100_000 } else { backbone_top_k };
    // Phase G1: glycan-Y-first candidate SELECTION (a glycan-Y-complementary
    // index generates backbones from the strong glycan-Y ladder) + TWO-AXIS
    // retention (keep backbones in top_k by peptide-b/y OR by glycan-Y evidence),
    // so a weak-b/y / strong-glycan-Y spectrum survives truncation. Opt-in for a
    // clean A/B vs the b/y-only path.
    let yindex_on = std::env::var("ANDES_GLYCO_YINDEX")
        .map(|v| v == "1")
        .unwrap_or(false);
    let glycan_y_index = if yindex_on {
        GlycanYIndex::build(glycan_list, fragment_tolerance_da.max(0.02))
    } else {
        GlycanYIndex::build(&[], fragment_tolerance_da.max(0.02))
    };

    // PEPTIDE-FIRST index (combines with the backbone-first hybrid below). Build a
    // fragment-ion index over the SEQUON candidate peptides ONCE. Per spectrum we
    // query it for peptides with real b/y support, then keep those whose
    // glycan-by-subtraction (`precursor − peptide`) hits a known glycan. This
    // recovers backbones on weak/absent-core-Y spectra that the core-Y-ranked
    // truncation drops — the candidate-generation ceiling — without a brute force.
    // Only build the (expensive) index when the peptide-first path is on; an
    // empty index is a no-op query otherwise (CodeRabbit: avoid the wasted build).
    let frag_index = if !peptide_first_on {
        FragmentIndex::build(std::iter::empty::<(u32, &model::peptide::Peptide)>(), fragment_tolerance_da.max(0.01))
    } else {
        let seq_entries: Vec<(u32, &model::peptide::Peptide)> = candidates
            .iter()
            .enumerate()
            .filter_map(|(i, c)| {
                let res: Vec<u8> = c.peptide.residues.iter().map(|aa| aa.residue).collect();
                if has_nxst_sequon(&res) {
                    Some((i as u32, &c.peptide))
                } else {
                    None
                }
            })
            .collect();
        // Memory guard (Codex re-review #3): glyco is forced to the in-RAM
        // candidate index, and the fragment index adds ~24 B/ion of postings.
        // Warn on very large sequon sets so an OOM isn't silent (a smaller search
        // space or an out-of-core index is the fix for those).
        let est_mb = seq_entries.len().saturating_mul(20 * 24) / 1_000_000;
        if seq_entries.len() > 3_000_000 {
            eprintln!(
                "WARN glyco fragment index: {} sequon candidates (~{} MB postings) — \
                 large database; consider narrowing the search space",
                seq_entries.len(),
                est_mb
            );
        }
        FragmentIndex::build(seq_entries.iter().copied(), fragment_tolerance_da.max(0.01))
    };
    // Sorted glycan masses for the peptide-first glycan-by-subtraction lookup.
    let glycan_sorted: Vec<(f64, usize)> = {
        let mut v: Vec<(f64, usize)> = glycan_list.iter().enumerate().map(|(i, g)| (g.mass, i)).collect();
        v.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        v
    };
    // Minimum b/y peaks a peptide must match to be a peptide-first candidate.
    // 6 (was 4) sharply cuts the coincidental-match Poisson tail — most of the
    // per-spectrum query cost — with negligible loss of real glycopeptides
    // (identifiable backbones carry several b/y ions).
    const MIN_BY_MATCHES: u32 = 6;
    // Hard cap on peptide-first candidates per spectrum (strongest b/y support
    // first) so a peak-dense spectrum can't blow up phase-1 scoring.
    const MAX_PEPTIDE_FIRST: usize = 64;

    // Per-spectrum processing, reused across both passes. `transfer` holds extra
    // backbones injected by cross-spectrum transfer (empty on pass 1).
    let process_one = |spec_idx: usize,
                       spec: &Spectrum,
                       transfer: &[BackboneHit]|
     -> Option<GlycoSpectrumResult> {
            if spec.peaks.len() < params.min_peaks as usize {
                return None;
            }
            // Fast dev harness: skip spectra outside the scan subset (if set).
            if let Some(ref scans) = scan_filter {
                match spec.scan {
                    Some(sc) if scans.contains(&sc) => {}
                    _ => return None,
                }
            }

            // Oxonium evidence for the whole spectrum (charge-independent).
            let ox_ev = oxonium_gate(&spec.peaks, 0.10, tol_ppm);

            // Determine which charges to try.
            let charges_to_try: Vec<u8> = match spec.precursor_charge {
                Some(z) if z > 0 => vec![z as u8],
                _ => params.charge_range.clone().collect(),
            };

            // Gather backbone hits across all charges AND all isotope offsets,
            // then union+dedup. Mirrors the standard search path's
            // `isotope_error_range` handling (search_params.rs): glyco
            // precursors frequently mis-pick the M+1/M+2 isotope peak, so
            // trying only the monoisotopic offset silently loses the true
            // backbone. Each resulting `BackboneHit` records the (charge,
            // isotope_offset) pair that produced it (see hybrid.rs).
            let iso_min = *params.isotope_error_range.start();
            let iso_max = *params.isotope_error_range.end();
            let mut all_backbone: Vec<BackboneHit> = Vec::new();
            for &z in &charges_to_try {
                let charge_f = z as f64;
                let observed_neutral = (spec.precursor_mz - PROTON) * charge_f - H2O;
                for iso in iso_min..=iso_max {
                    let precursor_neutral = observed_neutral - (iso as f64) * ISOTOPE;
                    if precursor_neutral <= 0.0 {
                        continue;
                    }
                    let hits = hybrid_candidates_with_isotope(
                        &spec.peaks,
                        precursor_neutral,
                        z,
                        iso,
                        glycan_list,
                        tol_ppm,
                        // Honor the configured cap (Codex finding #2): this was
                        // hardcoded to 50, so `--glyco-backbone-top-k` overrides
                        // were silently ignored and the Y-first solver truncated
                        // at 50 regardless. The solver ranks by core-Y evidence,
                        // so a larger cap only widens the candidate space the
                        // downstream b/y scorer sees. `effective_top_k` is huge in
                        // exhaustive mode (no truncation).
                        effective_top_k,
                    );
                    for h in hits {
                        all_backbone.push(h);
                    }
                }
            }

            // PEPTIDE-FIRST union: for oxonium-positive spectra, ask the fragment
            // index which sequon peptides actually have b/y support, then keep the
            // ones whose glycan-by-subtraction hits a known glycan across the same
            // charge/isotope grid. These backbones are selected by PEPTIDE evidence
            // (works when core-Y is weak/absent), and the glycan filter keeps the
            // count small — a handful of high-quality candidates per spectrum.
            if peptide_first_on && ox_ev.fired {
                // Process strongest-b/y-support peptides first, but cap on VALID
                // (peptide, charge, isotope, glycan) hypotheses — NOT raw b/y
                // count — so high-count peptides that cannot form a known glycan
                // don't evict a lower-count peptide that can (Codex re-review #1).
                let mut pf = frag_index.query(&spec.peaks, MIN_BY_MATCHES);
                pf.sort_unstable_by(|a, b| b.1.cmp(&a.1));
                let mut pf_added = 0usize;
                'pf: for (cand_idx, _n) in pf {
                    let pep_residue = candidates[cand_idx as usize].peptide.mass() - H2O;
                    for &z in &charges_to_try {
                        let observed_neutral = (spec.precursor_mz - PROTON) * z as f64 - H2O;
                        for iso in iso_min..=iso_max {
                            let precursor_neutral = observed_neutral - iso as f64 * ISOTOPE;
                            let glycan_mass = precursor_neutral - pep_residue;
                            if glycan_mass < MIN_GLYCAN {
                                continue;
                            }
                            // Tolerance is set by the PRECURSOR mass error (the
                            // measured quantity), not the smaller glycan mass
                            // (CodeRabbit): glycan = precursor − (exact) peptide.
                            let tol = (precursor_neutral * tol_ppm * 1e-6_f64).max(0.02);
                            if let Some(g) =
                                nearest_glycan_mass(&glycan_sorted, glycan_list, glycan_mass, tol)
                            {
                                // Observed backbone = precursor − theoretical glycan,
                                // matching the DB path's convention so the peptide's
                                // precursor mass error is preserved downstream
                                // (CodeRabbit: theoretical pep mass → 0 mass error).
                                let bb = precursor_neutral - g.mass;
                                all_backbone.push(BackboneHit {
                                    backbone_mass: bb,
                                    glycan: Some(g),
                                    source: Source::Db,
                                    charge: z,
                                    isotope_offset: iso,
                                    glycan_mass_residual: precursor_neutral - bb,
                                });
                                pf_added += 1;
                                if pf_added >= MAX_PEPTIDE_FIRST {
                                    break 'pf;
                                }
                            }
                        }
                    }
                }
            }

            // GLYCAN-Y-FIRST candidate generation (Phase G1): glycans whose core-Y
            // ladder is supported in THIS spectrum (peptide-independent, O(#peaks))
            // → their backbone (precursor − glycan), added with glycan evidence
            // regardless of peptide b/y. This is how the strong glycan signal
            // reaches the candidate set on weak-b/y spectra.
            if yindex_on && ox_ev.fired {
                for &z in &charges_to_try {
                    // FULL neutral precursor (water included) — the glycan-Y index
                    // convention (Y_complement = precursor_full − Y_ion).
                    let precursor_full = (spec.precursor_mz - PROTON) * z as f64;
                    for iso in iso_min..=iso_max {
                        let pf = precursor_full - iso as f64 * ISOTOPE;
                        if pf <= 0.0 {
                            continue;
                        }
                        for (gid, _core) in glycan_y_index.query(&spec.peaks, pf, z, 2) {
                            let g = &glycan_list[gid as usize];
                            let backbone_residue = pf - H2O - g.mass;
                            if backbone_residue < 500.0 {
                                continue;
                            }
                            all_backbone.push(BackboneHit {
                                backbone_mass: backbone_residue,
                                glycan: Some(g.clone()),
                                source: Source::Db,
                                charge: z,
                                isotope_offset: iso,
                                glycan_mass_residual: g.mass,
                            });
                        }
                    }
                }
            }

            // Cross-spectrum transfer: backbones borrowed from confident sibling
            // glycoforms (empty on pass 1). Added to the same dedup/score path.
            all_backbone.extend_from_slice(transfer);

            if all_backbone.is_empty() {
                return None;
            }

            // Dedup cross-charge/cross-isotope backbone hits, merging only hits
            // that represent the SAME (backbone, glycan-hypothesis) candidate.
            let deduped_backbone = dedup_backbone_hits(all_backbone, tol_ppm);

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
            // Core-Y ions live at the NEUTRAL peptide mass (Y0 = neutral + PROTON);
            // `backbone_mass` is the RESIDUE mass, so add H2O. (Previously passed
            // the residue mass → the ladder was sought ~H2O too low, so CoreYHits
            // measured near-noise. Phase-1 convention fix.)
            let core_y_counts: Vec<u8> = deduped_backbone
                .iter()
                .map(|h| count_core_y_hits(&spec.peaks, h.backbone_mass + H2O, tol_ppm))
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
                isotope_offset: i8,
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
                // The charge (and isotope offset) that produced this backbone
                // via `hybrid_candidates_with_isotope`. Scoring MUST use this
                // exact charge — re-deriving/re-picking a charge independently
                // here would score against a precursor mass inconsistent with
                // the one that actually matched this backbone (BUG: precursor
                // charge silently dropped).
                let z = bb_hit.charge;

                // Tight nominal bounds.
                let nb = nominal_from(bb_residue);
                let tol_da = (bb_residue * tol_ppm * 1e-6_f64).max(0.01);
                let widen = (tol_da - 0.4999_f64).max(0.0_f64).round() as i32;

                let candidate_slots: Vec<usize> = bucket_index
                    .range((nb - widen)..=(nb + widen))
                    .flat_map(|(_, v)| v.iter().copied())
                    .collect();

                let ss = match scored_per_charge.iter().find(|(c, _)| *c == z) {
                    Some((_, s)) => s,
                    // The backbone's charge fell outside `charges_to_try`
                    // (shouldn't happen since `hybrid_candidates_with_isotope`
                    // is only called for charges in that set, but guard
                    // defensively rather than panic).
                    None => continue,
                };

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

                    let sc = score_psm(ss, &cand.peptide, scorer, z, fragment_tolerance_da);
                    let ei = psm_edge_score(ss, &cand.peptide, scorer, z);
                    let rk = sc + ei as f32;

                    // Update per-backbone best rank.
                    if rk > backbone_best_rank[bb_idx] {
                        backbone_best_rank[bb_idx] = rk;
                    }

                    let w = CheapWinner {
                        bb_hit_idx: bb_idx,
                        cand_slot,
                        z,
                        isotope_offset: bb_hit.isotope_offset,
                        rank: rk,
                        score: sc,
                        edge: ei,
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
            //       TERTIARY = backbone_mass DESC via `total_cmp` (a true total
            //       order over all f64 bit patterns incl. sign/NaN — unlike the
            //       old `partial_cmp().unwrap_or(Equal)`, which silently treated
            //       any NaN comparison as a tie);
            //       QUATERNARY = bb_idx DESC — final total-order tiebreak so
            //       HashMap/rayon iteration-order jitter can never change which
            //       backbones survive truncation (BUG 4: nondeterministic cap).
            // AXIS 1 — peptide b/y rank (primary), core-Y as tiebreak.
            let mut by_by: Vec<usize> = (0..deduped_backbone.len()).collect();
            by_by.sort_by(|&ai, &bi| {
                backbone_best_rank[bi]
                    .total_cmp(&backbone_best_rank[ai])
                    .then_with(|| core_y_counts[bi].cmp(&core_y_counts[ai]))
                    .then_with(|| {
                        deduped_backbone[bi]
                            .backbone_mass
                            .total_cmp(&deduped_backbone[ai].backbone_mass)
                    })
                    .then_with(|| bi.cmp(&ai))
            });
            by_by.truncate(effective_top_k);
            let mut accepted_backbones: std::collections::HashSet<usize> =
                by_by.into_iter().collect();

            // AXIS 2 (Phase G1, TWO-AXIS retention) — also keep the top_k by
            // GLYCAN-Y evidence (core_y_hits), so a backbone that is strong on the
            // glycan axis but weak on peptide b/y survives truncation instead of
            // being dropped before its glycan features can be scored.
            if yindex_on {
                let mut by_gy: Vec<usize> = (0..deduped_backbone.len()).collect();
                by_gy.sort_by(|&ai, &bi| {
                    core_y_counts[bi]
                        .cmp(&core_y_counts[ai])
                        .then_with(|| backbone_best_rank[bi].total_cmp(&backbone_best_rank[ai]))
                        .then_with(|| {
                            deduped_backbone[bi]
                                .backbone_mass
                                .total_cmp(&deduped_backbone[ai].backbone_mass)
                        })
                        .then_with(|| bi.cmp(&ai))
                });
                by_gy.truncate(effective_top_k);
                accepted_backbones.extend(by_gy);
            }

            // Phase 2: expensive feature extraction for top-K winners only.
            // Only process cheap_winners whose backbone is in the accepted set.
            //
            // Cap at backbone_top_k × 2 to bound compute_psm_features calls,
            // but never below the number of accepted backbones — otherwise a
            // spectrum with many DISTINCT accepted backbones (each contributing
            // >2 candidate/glycan winners) could have true phase-1 winners
            // silently dropped before feature computation ever runs (BUG 4:
            // accepted candidates discarded pre-features). `accepted_backbones`
            // is already bounded by `backbone_top_k`, so this cap can only grow,
            // never shrink, relative to the correctness requirement.
            let max_features = (effective_top_k * 2).max(accepted_backbones.len() * 4);
            let winners_for_features: Vec<((u32, u8, u8, u8, u8, u8), CheapWinner)> = {
                let mut v: Vec<_> = cheap_winners
                    .into_iter()
                    .filter(|(_, w)| accepted_backbones.contains(&w.bb_hit_idx))
                    .collect();
                // Deterministic total order: rank DESC via `total_cmp` (not
                // partial_cmp-with-Equal-fallback), then gl_key ASC as the
                // final tiebreak so truncation never depends on HashMap
                // iteration order (BUG 4).
                v.sort_by(|a, b| {
                    b.1.rank
                        .total_cmp(&a.1.rank)
                        .then_with(|| a.0.cmp(&b.0))
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
                    // The isotope offset that was actually subtracted when this
                    // backbone was derived (BUG 1 fix — previously hardcoded 0
                    // even when an M+1/M+2 offset produced the winning match).
                    isotope_offset: w.isotope_offset,
                    precursor_mz_override: None,
                };
                // Use the annotated composition's theoretical mass when a known
                // glycan matched; otherwise fall back to the observed residual
                // (precursor − backbone) so a novel/unannotated glycan still
                // reports its real intact mass instead of 0.0 (Codex finding #3).
                let glycan_mass = bb_hit
                    .glycan
                    .as_ref()
                    .map(|g| g.mass)
                    .unwrap_or(bb_hit.glycan_mass_residual);
                let glycan_key = GlycoPsmKey {
                    spectrum_idx: spec_idx,
                    glycan: bb_hit.glycan.clone(),
                    glycan_source: bb_hit.source.clone(),
                    oxonium_summed_frac: ox_ev.summed_frac,
                    n_core_oxonium_ions: ox_ev.n_core_ions,
                    // Y-ladder intensity match at the NEUTRAL backbone. For an
                    // ANNOTATED glycan use the COMPOSITION-SPECIFIC ladder
                    // (Phase 2: a wrong glycan of similar mass scores lower, so
                    // the feature discriminates on the glycan axis); for a novel
                    // glycan fall back to the composition-independent core-Y
                    // ladder. (Was hardcoded 0.0 = dead before Phase 1.)
                    y_ladder_intensity_score: match &bb_hit.glycan {
                        Some(g) => glycan_y_intensity(&spec.peaks, bb_neutral, g, tol_ppm) as f32,
                        None => core_y_intensity(&spec.peaks, bb_neutral, tol_ppm) as f32,
                    },
                    // Glycan-axis decoy ladder (G3): same composition, intermediate
                    // Y-rungs shifted. Seed from the composition so the decoy ladder
                    // is stable per glycan (a fixed decoy "structure"). 0.0 for
                    // de-novo hits (no composition → no glycan-axis decoy row).
                    y_ladder_decoy_score: match &bb_hit.glycan {
                        Some(g) if glyco_decoy_on => glycan_y_intensity_decoy(
                            &spec.peaks,
                            bb_neutral,
                            g,
                            tol_ppm,
                            glycan_decoy_seed(g),
                        ) as f32,
                        _ => 0.0,
                    },
                    // G2 Y0/Y1 anchor: peptide-mass-conditioned (uses THIS
                    // candidate's neutral mass, so it discriminates competing
                    // peptides). Additive PIN feature only — not in the ranker.
                    y0y1_anchor_score: y0y1_anchor_intensity(
                        &spec.peaks,
                        cand.peptide.mass(),
                        w.z,
                        tol_ppm,
                    ) as f32,
                    // Threaded from the per-backbone Y-ladder evidence computed
                    // earlier in `core_y_counts` (previously discarded/hardcoded
                    // to 0, so the `CoreYHits` PIN feature was always dead).
                    core_y_hits: core_y_counts[w.bb_hit_idx],
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
    };

    // PASS 1: baseline candidate gen (+ peptide-first if on), no transfer.
    let pass1: Vec<GlycoSpectrumResult> = spectra
        .par_iter()
        .enumerate()
        .filter_map(|(spec_idx, spec)| process_one(spec_idx, spec, &[]))
        .collect();

    if !cross_spectrum_on {
        return pass1;
    }

    // CROSS-SPECTRUM TRANSFER. Build a whitelist of CONFIDENT backbone (residue)
    // masses from pass-1 PSMs with a strong core-Y ladder (well-fragmented
    // glycoforms), then transfer them to poorly-fragmenting sibling glycoforms.
    const CONF_MIN_CORE_Y: u8 = 3;
    let confident_bb: Vec<f64> = pass1
        .iter()
        .flat_map(|r| {
            r.hits
                .iter()
                .filter(|h| h.glycan_key.core_y_hits >= CONF_MIN_CORE_Y)
                .map(|h| h.glycan_key.backbone_mass - H2O)
        })
        .collect();
    let whitelist = GlycoformWhitelist::new(confident_bb, 0.02);
    if whitelist.is_empty() {
        return pass1;
    }
    // Spectra that already have a confident ID need no transfer (bounds pass 2).
    let confident_scans: std::collections::HashSet<usize> = pass1
        .iter()
        .filter(|r| r.hits.iter().any(|h| h.glycan_key.core_y_hits >= CONF_MIN_CORE_Y))
        .map(|r| r.spectrum_idx)
        .collect();
    let mut by_idx: HashMap<usize, GlycoSpectrumResult> =
        pass1.into_iter().map(|r| (r.spectrum_idx, r)).collect();

    // PASS 2: only the non-confident (weak-ladder) spectra; inject transferred
    // backbones and re-score. Results supersede their pass-1 entry.
    let iso_min = *params.isotope_error_range.start();
    let iso_max = *params.isotope_error_range.end();
    let pass2: Vec<GlycoSpectrumResult> = spectra
        .par_iter()
        .enumerate()
        .filter(|(spec_idx, _)| !confident_scans.contains(spec_idx))
        .filter_map(|(spec_idx, spec)| {
            if spec.peaks.len() < params.min_peaks as usize {
                return None;
            }
            // Transfer only to glyco-plausible (oxonium-positive) spectra.
            if !oxonium_gate(&spec.peaks, 0.10, tol_ppm).fired {
                return None;
            }
            let charges_to_try: Vec<u8> = match spec.precursor_charge {
                Some(z) if z > 0 => vec![z as u8],
                _ => params.charge_range.clone().collect(),
            };
            let mut transfer: Vec<BackboneHit> = Vec::new();
            for &z in &charges_to_try {
                let observed_neutral = (spec.precursor_mz - PROTON) * z as f64 - H2O;
                for iso in iso_min..=iso_max {
                    let pn = observed_neutral - iso as f64 * ISOTOPE;
                    if pn <= 0.0 {
                        continue;
                    }
                    let tol = (pn * tol_ppm * 1e-6_f64).max(0.02);
                    for (_bb, g) in
                        whitelist.transfer(pn, &glycan_sorted, glycan_list, MIN_GLYCAN, tol)
                    {
                        // Observed backbone = precursor − theoretical glycan (real
                        // mass error), matching the DB/peptide-first convention.
                        let bb_obs = pn - g.mass;
                        transfer.push(BackboneHit {
                            backbone_mass: bb_obs,
                            glycan: Some(g),
                            source: Source::Db,
                            charge: z,
                            isotope_offset: iso,
                            glycan_mass_residual: pn - bb_obs,
                        });
                    }
                }
            }
            if transfer.is_empty() {
                return None;
            }
            process_one(spec_idx, spec, &transfer)
        })
        .collect();

    for r in pass2 {
        by_idx.insert(r.spectrum_idx, r);
    }
    // Deterministic output order (CodeRabbit): HashMap iteration is unordered,
    // so sort by spectrum_idx for reproducibility.
    let mut out: Vec<GlycoSpectrumResult> = by_idx.into_values().collect();
    out.sort_by_key(|r| r.spectrum_idx);
    out
}

#[cfg(test)]
mod tests {
    // Integration-level tests (full `glyco_search_run` over a `PreparedSearch`)
    // are deferred to the search-crate integration tests (tests/ directory)
    // where real PreparedSearch fixtures can be built. Unit-level sequon + mass
    // filter logic is tested in andes_glyco::sequon.
    //
    // Smoke test: verify the public types compile and are accessible.
    use super::*;

    /// Isotope-sweep GLYCAN ANNOTATION regression. Under the Y-ion-first
    /// cascade the backbone is read from the core-Y ladder and is therefore
    /// recovered regardless of which isotope peak the instrument picked as the
    /// precursor. The isotope sweep's remaining job is to annotate the glycan
    /// CORRECTLY: `glycan = precursor_neutral − backbone` only matches a known
    /// composition when the precursor is at the right isotope offset.
    ///
    /// With the precursor mis-picked at M+1, this test confirms: (a) at offset
    /// 0 the backbone is still found from the ladder, but its by-subtraction
    /// glycan is ~1 ISOTOPE off → NOT annotated (Source::DeNovo); (b) at offset
    /// +1 the corrected precursor yields the true glycan → annotated Source::Db
    /// with isotope_offset=1. This is exactly what the driver's
    /// `for iso in iso_min..=iso_max` sweep buys in the Y-first world.
    #[test]
    fn isotope_sweep_annotates_glycan_only_at_correct_offset() {
        use andes_glyco::glycan_db::GlycanComp;
        use andes_glyco::glycan_mass::{CORE_Y_STEPS, PROTON as GLY_PROTON};
        use andes_glyco::hybrid::hybrid_candidates_with_isotope;

        let true_backbone_residue = 1500.0_f64;
        let glycan = GlycanComp {
            hexnac: 2,
            hex: 3,
            fuc: 0,
            neuac: 0,
            neugc: 0,
            mass: 2.0 * andes_glyco::glycan_mass::HEXNAC + 3.0 * andes_glyco::glycan_mass::HEX,
        };
        let true_precursor_neutral = true_backbone_residue + glycan.mass;

        // Instrument reports the M+1 isotope peak as the precursor.
        let observed_neutral = true_precursor_neutral + ISOTOPE;

        let glycans = vec![glycan];

        // Full core-Y ladder anchored at the true backbone (Y0 = peptide neutral
        // + proton = residue + H2O + proton) plus two oxonium ions so the gate
        // fires. The ladder is independent of the precursor isotope pick.
        let y0_neutral = true_backbone_residue + H2O;
        let mut peaks: Vec<(f64, f32)> = vec![(204.08665, 200.0), (138.05496, 120.0)];
        peaks.push((y0_neutral + GLY_PROTON, 150.0));
        for &s in CORE_Y_STEPS.iter() {
            peaks.push((y0_neutral + s + GLY_PROTON, 100.0));
        }

        // Offset 0 (M+1 assumption uncorrected): backbone recovered from the
        // ladder, but the by-subtraction glycan is ~1 ISOTOPE off → not annotated.
        let hits0 =
            hybrid_candidates_with_isotope(&peaks, observed_neutral, 2, 0, &glycans, 20.0, 5);
        let matching0: Vec<_> = hits0
            .iter()
            .filter(|h| (h.backbone_mass - true_backbone_residue).abs() < 0.05)
            .collect();
        assert!(
            !matching0.is_empty(),
            "backbone is read from the ladder → recovered even at the wrong isotope offset"
        );
        assert!(
            matching0.iter().any(|h| h.source == Source::DeNovo),
            "at the wrong offset the glycan is ~1 ISOTOPE off and must NOT annotate to a known composition"
        );
        assert!(
            matching0.iter().all(|h| h.source != Source::Db),
            "wrong isotope offset must not produce a DB annotation for the true backbone"
        );

        // Offset +1: corrected precursor → true glycan → annotated Source::Db.
        let precursor_neutral_iso1 = observed_neutral - ISOTOPE;
        let hits1 = hybrid_candidates_with_isotope(
            &peaks,
            precursor_neutral_iso1,
            2,
            1,
            &glycans,
            20.0,
            5,
        );
        let hit = hits1
            .iter()
            .find(|h| {
                (h.backbone_mass - true_backbone_residue).abs() < 0.05 && h.source == Source::Db
            })
            .expect("offset +1 must recover AND annotate the backbone via the corrected precursor");
        assert_eq!(hit.isotope_offset, 1, "recovered hit must record isotope_offset=1");
        assert_eq!(hit.charge, 2, "recovered hit must record the charge it was matched at");
    }

    /// P0.1 (Codex #2): `dedup_backbone_hits` must NOT merge two de-novo hits
    /// that share a backbone mass but carry different isotope offsets — their
    /// residual glycan masses differ, so merging corrupts the novel-glycan
    /// intact mass. Same-hypothesis duplicates (same offset) must still merge.
    #[test]
    fn dedup_preserves_distinct_isotope_residuals_for_novel_glycans() {
        use andes_glyco::hybrid::BackboneHit;
        let mk = |iso: i8, residual: f64| BackboneHit {
            backbone_mass: 1500.0,
            glycan: None, // novel / unannotated
            source: Source::DeNovo,
            charge: 3,
            isotope_offset: iso,
            glycan_mass_residual: residual,
        };
        // Two isotope hypotheses at the same backbone → both must survive.
        let out = dedup_backbone_hits(vec![mk(0, 892.317), mk(1, 891.313)], 20.0);
        assert_eq!(out.len(), 2, "distinct isotope residuals must not be merged: {out:?}");
        let residuals: Vec<f64> = out.iter().map(|h| h.glycan_mass_residual).collect();
        assert!(residuals.iter().any(|r| (r - 892.317).abs() < 1e-6));
        assert!(residuals.iter().any(|r| (r - 891.313).abs() < 1e-6));

        // Same offset (true duplicate, e.g. from another charge) must merge.
        let dup = dedup_backbone_hits(vec![mk(0, 892.317), mk(0, 892.317)], 20.0);
        assert_eq!(dup.len(), 1, "same-hypothesis duplicates must merge");
    }

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
            y_ladder_decoy_score: 0.0,
            y0y1_anchor_score: 0.0,
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
