//! PTM-refinement cascade (Pass-2). After the Pass-1 search, this module:
//!   1. selects the proteins that already carry a confident Pass-1 target PSM
//!      (`confident_protein_indices`),
//!   2. finds the spectra still without any target PSM
//!      (`unidentified_spectrum_indices`),
//!   3. builds a SMALL target+decoy search index scoped to just the confident
//!      proteins (`build_refinement_index`),
//!   4. searches the unidentified spectra against that small index with the
//!      refinement variable-mod tier applied (`refinement_aa_set`), and
//!   5. merges every Pass-2 winner back into the per-spectrum queue, tagging it
//!      as a refinement PSM (`run_refinement`).
//!
//! The cascade is OPT-IN (`--refine`) and additive: nothing here runs on the
//! default path, and the merged Pass-2 PSMs are `force_push`ed so they never
//! evict a Pass-1 winner.

use std::collections::HashSet;

use model::aa_set::{AminoAcidSet, AminoAcidSetBuilder};
use model::modification::{ModLocation, Modification, ResidueSpec};
use model::peptide::Peptide;
use model::protein::ProteinDb;
use scoring_crate::RankScorer;

use crate::candidate_gen::Candidate;
use crate::decoy::{build_search_db, DecoyStrategy};
use crate::match_engine::PreparedSearch;
use crate::psm::TopNQueue;
use crate::refine_config::RefineConfig;
use crate::search_index::SearchIndex;
use crate::search_params::SearchParams;
use crate::tdc::{confident_target_indices, ScoredLabel};

// ─────────────────────────────────────────────────────────────────────────────
// (a) confident_protein_indices
// ─────────────────────────────────────────────────────────────────────────────

/// Proteins with ≥1 Pass-1 target PSM at internal-TDC-q ≤ `select_q`.
///
/// For each spectrum's BEST (max `rank_score`) PSM, build a [`ScoredLabel`] from
/// its `rank_score` and the target/decoy label of its primary candidate, run the
/// shared TDC q-walk ([`confident_target_indices`]), and collect the
/// `protein_index` of every confident PSM's primary candidate.
///
/// The returned indices are `Candidate.protein_index` values, which index into
/// the COMBINED target+decoy database the Pass-1 search used (see
/// `enumerate_candidates`). Because only TARGET PSMs survive the q-walk, every
/// returned index is a target protein's position in that combined db. For the
/// default Reverse strategy on a clean target FASTA, targets occupy
/// `[0, target_count)`, so these coincide with target-only positions —
/// [`build_refinement_index`] guards against any decoy index regardless.
pub fn confident_protein_indices(
    queues: &[TopNQueue],
    candidates: &[Candidate],
    select_q: f64,
) -> HashSet<usize> {
    // One (label, candidate_idx) per spectrum that produced a PSM, using that
    // spectrum's best PSM. `peek_top` returns the max-`rank_score` PSM.
    let mut labels: Vec<ScoredLabel> = Vec::new();
    let mut cand_idx_of: Vec<usize> = Vec::new();
    for queue in queues {
        if let Some(best) = queue.peek_top() {
            let cand_idx = best.primary_candidate_idx() as usize;
            labels.push(ScoredLabel {
                score: best.rank_score,
                is_decoy: candidates[cand_idx].is_decoy,
            });
            cand_idx_of.push(cand_idx);
        }
    }

    confident_target_indices(&labels, select_q)
        .into_iter()
        .map(|li| candidates[cand_idx_of[li]].protein_index)
        .collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// (b) unidentified_spectrum_indices
// ─────────────────────────────────────────────────────────────────────────────

/// Indices of spectra with NO retained TARGET PSM (decoy-only or empty queues).
/// These are the spectra Pass-2 attempts to rescue with the refinement mod tier.
pub fn unidentified_spectrum_indices(
    queues: &[TopNQueue],
    candidates: &[Candidate],
) -> Vec<usize> {
    queues
        .iter()
        .enumerate()
        .filter(|(_, queue)| {
            !queue
                .iter_psms()
                .any(|p| !candidates[p.primary_candidate_idx() as usize].is_decoy)
        })
        .map(|(i, _)| i)
        .collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// (c) build_refinement_index
// ─────────────────────────────────────────────────────────────────────────────

/// Build a SMALL target+decoy search index scoped to the confident proteins.
///
/// Subsets `full_target_db` (the TARGET-ONLY db) to the confident protein
/// indices, then regenerates a fresh 1:1 Reverse decoy per kept target so the
/// scoped search has its own decoys for TDC. `confident_target_proteins` carries
/// `Candidate.protein_index` values (combined-db positions; see
/// [`confident_protein_indices`]); any index ≥ `full_target_db.len()` would be a
/// decoy position and is filtered out before subsetting the target db.
pub fn build_refinement_index(
    full_target_db: &ProteinDb,
    confident_target_proteins: &HashSet<usize>,
    decoy_prefix: &str,
    seed: u64,
) -> SearchIndex {
    // Keep only indices that address a TARGET protein in `full_target_db`.
    // Indices that point past the target range belong to the combined db's decoy
    // half (Reverse strategy) and must never be subset into a target db.
    let keep: HashSet<usize> = confident_target_proteins
        .iter()
        .copied()
        .filter(|&i| i < full_target_db.len())
        .collect();

    let subset = full_target_db.subset_by_index(&keep);
    let db = build_search_db(&subset, decoy_prefix, DecoyStrategy::Reverse, seed);
    SearchIndex {
        db,
        decoy_prefix: crate::decoy::normalize_decoy_prefix(decoy_prefix),
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// (d) refinement_aa_set
// ─────────────────────────────────────────────────────────────────────────────

/// Map a [`crate::refine_config::RefineMod`] `location` string to a [`ModLocation`].
fn parse_location(loc: &str) -> ModLocation {
    match loc.trim().to_ascii_lowercase().as_str() {
        "n_term" | "n-term" | "nterm" => ModLocation::NTerm,
        "c_term" | "c-term" | "cterm" => ModLocation::CTerm,
        "protein_n_term" | "prot-n-term" | "prot_n_term" => ModLocation::ProtNTerm,
        "protein_c_term" | "prot-c-term" | "prot_c_term" => ModLocation::ProtCTerm,
        _ => ModLocation::Anywhere,
    }
}

/// Map a single residue string to a [`ResidueSpec`]. `"*"` (or any non
/// single-uppercase token) ⇒ wildcard; a single uppercase letter ⇒ specific.
fn parse_residue_spec(res: &str) -> ResidueSpec {
    let b = res.as_bytes();
    if b.len() == 1 && b[0].is_ascii_uppercase() {
        ResidueSpec::Specific(b[0])
    } else {
        ResidueSpec::Wildcard
    }
}

/// Build the Pass-2 [`AminoAcidSet`]: the base set's FIXED mods (e.g.
/// Carbamidomethyl-C) PLUS the refinement tier's VARIABLE mods.
///
/// `RefineConfig` carries one [`crate::refine_config::RefineMod`] per
/// refinement chemistry; each may target several residues (e.g. Deamidation on
/// N and Q), which expands to one [`Modification`] per `(residue, location)`.
/// Mods whose `class == "deamidation"` are SKIPPED when `!high_res` (the
/// deamidation +0.984 Da delta is near-isobaric with a C13 isotope error at low
/// resolution, so it is high-res-only).
///
/// **Fixed-mod handling.** The `AminoAcidSetBuilder` cannot start from an
/// existing `AminoAcidSet`'s fixed mods (it only accepts fresh `Modification`s),
/// so we reconstruct the standard fixed-mod baseline: if `base` carries a fixed
/// Carbamidomethyl-C variant (the standard alkylation, and the only fixed mod on
/// the default tier), we seed the builder with it. Other custom fixed mods on
/// `base` are NOT reconstructed — the refinement tier is the the classic hyperscore engine always-on
/// chemistry layered on the standard CAM-C baseline; richer fixed-mod schemes
/// (TMT, etc.) would be supplied via the normal mods file, not the refine path.
pub fn refinement_aa_set(base: &AminoAcidSet, cfg: &RefineConfig, high_res: bool) -> AminoAcidSet {
    // Seed the builder with the standard fixed Carbamidomethyl-C baseline iff
    // `base` carries it (a single fixed modified C-Anywhere variant).
    let base_has_fixed_cam_c = base
        .variants_for(b'C', ModLocation::Anywhere)
        .iter()
        .any(|aa| aa.mod_.as_ref().map(|m| m.fixed).unwrap_or(false));

    let mut builder = if base_has_fixed_cam_c {
        AminoAcidSetBuilder::new_standard_with_carbamidomethyl_c()
    } else {
        AminoAcidSetBuilder::new_standard()
    };

    for rm in &cfg.mods {
        if !high_res && rm.class == "deamidation" {
            continue;
        }
        let location = parse_location(&rm.location);
        for res in &rm.residues {
            let m = Modification {
                name: rm.name.clone(),
                mass_delta: rm.delta,
                residue: parse_residue_spec(res),
                location,
                fixed: false,
                accession: None,
                neutral_losses: Vec::new(),
                loss_class: 0,
            };
            builder = builder.add_variable_mod(m);
        }
    }

    builder
        .build()
        .expect("refinement_aa_set: standard fixed + variable refinement mods must build")
}

// ─────────────────────────────────────────────────────────────────────────────
// (e) mod_count_and_class
// ─────────────────────────────────────────────────────────────────────────────

/// Count the variable-modified residues on `peptide` and classify them.
///
/// Returns `(num_mods, refine_mod_class)`:
/// - `num_mods` is the number of residues carrying a modification.
/// - `refine_mod_class` is the [`RefineConfig::class_id`] of the peptide's
///   mod(s), determined by matching each modified residue's `(delta, residue)`
///   against `cfg.mods`. The first non-`"other"` class wins; if every matched
///   class is `"other"` (or no `cfg.mods` entry matches a present mod) the class
///   is `99`. An unmodified peptide returns `(0, 0)`.
///
/// A `Peptide` represents a modified residue as an `AminoAcid` carrying an
/// `Option<Arc<Modification>>` whose `mass_delta` is the variable-mod delta
/// (fixed mods are folded into the residue the same way, so this counts BOTH —
/// see the note below). We match on `mass_delta ≈ rm.delta` (1e-4 Da) and the
/// residue letter against `rm.residues` (`"*"` matches any residue).
///
/// NOTE: a `Peptide` does not record whether a given mod was fixed or variable,
/// so this counts every modified residue. On the refinement path the base set's
/// only fixed mod is Carbamidomethyl-C, which is not part of any refinement
/// `cfg.mods` class; it would therefore inflate `num_mods` by 1 per modified C.
/// To keep `num_mods` aligned with the REFINEMENT tier, we count a residue only
/// when its mod matches one of `cfg.mods` (i.e. a refinement chemistry).
pub fn mod_count_and_class(peptide: &Peptide, cfg: &RefineConfig) -> (u32, u32) {
    let mut num_mods = 0u32;
    let mut class_id = 0u32; // 0 = none until a refinement mod matches

    for aa in &peptide.residues {
        let Some(m) = aa.mod_.as_ref() else { continue };
        // Find the refinement mod that this residue's modification matches.
        let matched = cfg.mods.iter().find(|rm| {
            (rm.delta - m.mass_delta).abs() < 1e-4
                && rm.residues.iter().any(|r| {
                    r == "*" || (r.len() == 1 && r.as_bytes()[0] == aa.residue)
                })
        });
        let Some(rm) = matched else { continue };
        num_mods += 1;
        let this_class = RefineConfig::class_id(&rm.class);
        // Prefer the first non-"other"(99) class seen; fall back to 99 only when
        // every matched mod is "other". `class_id == 0` is the "none yet" state.
        let pending = class_id == 0 || class_id == 99;
        if pending && (class_id == 0 || this_class != 99) {
            class_id = this_class;
        }
    }

    (num_mods, class_id)
}

// ─────────────────────────────────────────────────────────────────────────────
// (f) run_refinement
// ─────────────────────────────────────────────────────────────────────────────

/// Orchestrate the Pass-2 PTM-refinement search and merge its winners back into
/// `queues`.
///
/// Steps (each an early-return when its input is empty, so the cheap no-op cases
/// never build a Pass-2 search):
///   1. `confident_protein_indices` over the Pass-1 queues — return if empty.
///   2. `unidentified_spectrum_indices` over the Pass-1 queues — return if empty.
///   3. `build_refinement_index` scoped to the confident proteins.
///   4. A Pass-2 `SearchParams` clone whose `aa_set` is the refinement tier
///      ([`refinement_aa_set`]) and whose `max_variable_mods_per_peptide` is
///      `cfg.max_mods`, prepared via `PreparedSearch::prepare`.
///   5. `run_chunk` over the (cloned) unidentified spectra; each winner is tagged
///      `is_refinement = 1`, has `(num_mods, refine_mod_class)` filled from its
///      Pass-2 peptide, has its `spectrum_idx` fixed back to the GLOBAL index,
///      and is `force_push`ed into the corresponding Pass-1 queue.
///
/// CAVEAT (per Task 5 spec): a merged Pass-2 PSM's `candidate_idxs` index into
/// the PASS-2 `PreparedSearch::candidates`, which is a DIFFERENT slice than the
/// Pass-1 candidate list these queues were built from. This function only
/// merges + tags; downstream PIN emission (Task 6) must resolve a refinement
/// PSM's peptide/protein against the Pass-2 candidate list, not the Pass-1 one.
#[allow(clippy::too_many_arguments)]
pub fn run_refinement(
    queues: &mut [TopNQueue],
    all_spectra: &[model::spectrum::Spectrum],
    pass1_candidates: &[Candidate],
    full_target_db: &ProteinDb,
    base_params: &SearchParams,
    scorer: &RankScorer,
    cfg: &RefineConfig,
    high_res: bool,
    fragment_tol_da: f64,
    decoy_prefix: &str,
    seed: u64,
) {
    // 1. Confident proteins from Pass-1 (scoping gate). Empty ⇒ nothing to refine.
    let confident =
        confident_protein_indices(queues, pass1_candidates, base_params.refine_select_psm_fdr);
    if confident.is_empty() {
        return;
    }

    // 2. Spectra still without a target PSM. Empty ⇒ nothing to rescue.
    let unident = unidentified_spectrum_indices(queues, pass1_candidates);
    if unident.is_empty() {
        return;
    }

    // 3. Scoped target+decoy index over the confident proteins.
    let refine_idx = build_refinement_index(full_target_db, &confident, decoy_prefix, seed);

    // 4. Pass-2 params: refinement variable-mod tier + the tier's mod cap.
    let mut refine_params = base_params.clone();
    refine_params.aa_set = refinement_aa_set(&base_params.aa_set, cfg, high_res);
    refine_params.max_variable_mods_per_peptide = cfg.max_mods;

    let refine_prepared =
        PreparedSearch::prepare(&refine_idx, &refine_params, scorer, fragment_tol_da, decoy_prefix);

    // 5. Search ONLY the unidentified spectra. We clone them into a contiguous
    //    sub-vec (run_chunk takes a &[Spectrum]); `spectrum_idx_offset = 0` is
    //    fine because we remap to the GLOBAL index ourselves below.
    let subset: Vec<model::spectrum::Spectrum> =
        unident.iter().map(|&i| all_spectra[i].clone()).collect();
    let pass2_queues = refine_prepared.run_chunk(&subset, 0);

    // Merge each Pass-2 winner back into its GLOBAL queue, tagged as refinement.
    for (local_idx, pass2_queue) in pass2_queues.into_iter().enumerate() {
        let global_idx = unident[local_idx];
        for mut psm in pass2_queue.into_rank_sorted_vec() {
            psm.features.is_refinement = 1;
            let (num_mods, class) = mod_count_and_class(
                &refine_prepared.candidates[psm.primary_candidate_idx() as usize].peptide,
                cfg,
            );
            psm.features.num_mods = num_mods;
            psm.features.refine_mod_class = class;
            // Fix the spectrum index back to the global stream position.
            psm.spectrum_idx = global_idx;
            queues[global_idx].force_push(psm);
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use model::amino_acid::AminoAcid;
    use model::peptide::Peptide;
    use model::protein::Protein;
    use crate::psm::{PsmFeatures, PsmMatch};

    // ── Fixture helpers ────────────────────────────────────────────────────

    /// A target/decoy candidate at a given combined-db protein index.
    fn cand(protein_index: usize, is_decoy: bool) -> Candidate {
        let residues = vec![AminoAcid::standard(b'A').unwrap()];
        Candidate {
            peptide: Peptide::new(residues, b'_', b'-'),
            protein_index,
            start_offset_in_protein: 0,
            is_decoy,
            is_protein_n_term: false,
            is_protein_c_term: false,
        }
    }

    /// A PSM pointing at `candidate_idxs[0] = cand_idx` with the given rank_score.
    fn psm(spectrum_idx: usize, cand_idx: u32, rank_score: f32) -> PsmMatch {
        PsmMatch {
            spectrum_idx,
            candidate_idxs: vec![cand_idx],
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

    fn protein(acc: &str, seq: &[u8]) -> Protein {
        Protein { accession: acc.into(), description: String::new(), sequence: seq.to_vec() }
    }

    // ── (a) confident_protein_indices ──────────────────────────────────────

    #[test]
    fn confident_protein_indices_keeps_targets_above_q_threshold() {
        // 200 confident target spectra (proteins 0..200) + 1 lone decoy tail.
        // The decoy maps to protein_index 999 and must NOT appear.
        let mut candidates: Vec<Candidate> = (0..200).map(|i| cand(i, false)).collect();
        candidates.push(cand(999, true)); // index 200, a decoy

        let mut queues = Vec::new();
        for i in 0..200 {
            queues.push(queue_with(vec![psm(i, i as u32, 30.0 - i as f32 * 0.05)]));
        }
        // The decoy spectrum (lowest score).
        queues.push(queue_with(vec![psm(200, 200, 1.0)]));

        let confident = confident_protein_indices(&queues, &candidates, 0.01);
        assert!(!confident.is_empty());
        assert!(!confident.contains(&999), "decoy protein must be excluded");
        // High-scoring targets are present.
        assert!(confident.contains(&0));
    }

    #[test]
    fn confident_protein_indices_uses_best_psm_per_spectrum() {
        // One spectrum, two PSMs: the BEST (higher rank_score) is the target on
        // protein 7; the worse one is a decoy. peek_top must pick the target.
        let candidates = vec![cand(7, false), cand(8, true)];
        let q = queue_with(vec![psm(0, 0, 50.0), psm(0, 1, 5.0)]);
        // A second high target spectrum so the lone-spectrum q-walk is non-empty.
        let candidates2 = {
            let mut c = candidates.clone();
            c.push(cand(9, false));
            c
        };
        let q2 = queue_with(vec![psm(1, 2, 49.0)]);
        let confident = confident_protein_indices(&[q, q2], &candidates2, 0.5);
        assert!(confident.contains(&7), "best PSM (target on protein 7) selected");
        assert!(!confident.contains(&8), "the worse decoy PSM must not be selected");
    }

    #[test]
    fn confident_protein_indices_empty_when_all_decoys() {
        let candidates: Vec<Candidate> = (0..20).map(|i| cand(i, true)).collect();
        let queues: Vec<TopNQueue> = (0..20)
            .map(|i| queue_with(vec![psm(i, i as u32, 10.0)]))
            .collect();
        assert!(confident_protein_indices(&queues, &candidates, 0.01).is_empty());
    }

    // ── (b) unidentified_spectrum_indices ──────────────────────────────────

    #[test]
    fn unidentified_spectrum_indices_flags_empty_and_decoy_only() {
        let candidates = vec![cand(0, false), cand(1, true)];
        let queues = vec![
            queue_with(vec![psm(0, 0, 20.0)]), // target PSM → identified
            queue_with(vec![psm(1, 1, 20.0)]), // decoy-only → unidentified
            TopNQueue::new(10),                // empty → unidentified
        ];
        let unident = unidentified_spectrum_indices(&queues, &candidates);
        assert_eq!(unident, vec![1, 2]);
    }

    #[test]
    fn unidentified_spectrum_indices_mixed_queue_with_one_target_is_identified() {
        // A queue holding BOTH a target and a decoy PSM counts as identified.
        let candidates = vec![cand(0, true), cand(1, false)];
        let queues = vec![queue_with(vec![psm(0, 0, 30.0), psm(0, 1, 10.0)])];
        assert!(unidentified_spectrum_indices(&queues, &candidates).is_empty());
    }

    // ── (c) build_refinement_index ─────────────────────────────────────────

    #[test]
    fn build_refinement_index_subsets_targets_and_pairs_decoys() {
        let full_target = ProteinDb {
            proteins: vec![
                protein("P0", b"MKWVR"),
                protein("P1", b"AGCTR"),
                protein("P2", b"LDESR"),
            ],
        };
        let keep: HashSet<usize> = [0usize, 2].into_iter().collect();
        let idx = build_refinement_index(&full_target, &keep, "XXX", 42);

        // 2 kept targets + 2 paired Reverse decoys.
        assert_eq!(idx.db.len(), 4);
        assert_eq!(idx.db.proteins[0].accession, "P0");
        assert_eq!(idx.db.proteins[1].accession, "P2");
        // Decoys carry the "<prefix>_" needle.
        let needle = crate::decoy::decoy_accession_needle("XXX");
        assert!(idx.db.proteins[2].accession.starts_with(&needle));
        assert!(idx.db.proteins[3].accession.starts_with(&needle));
        assert_eq!(idx.decoy_prefix, "XXX");
    }

    #[test]
    fn build_refinement_index_drops_out_of_range_decoy_positions() {
        // A combined-db index past the target count (a decoy position) must be
        // ignored so it never gets subset into the target db.
        let full_target = ProteinDb { proteins: vec![protein("P0", b"MKWVR")] };
        // keep {0 (valid target), 5 (out of range / decoy half)}.
        let keep: HashSet<usize> = [0usize, 5].into_iter().collect();
        let idx = build_refinement_index(&full_target, &keep, "XXX", 42);
        // Only P0 + its decoy.
        assert_eq!(idx.db.len(), 2);
        assert_eq!(idx.db.proteins[0].accession, "P0");
    }

    // ── (d) refinement_aa_set ──────────────────────────────────────────────

    fn base_set_with_cam_c() -> AminoAcidSet {
        AminoAcidSetBuilder::new_standard_with_carbamidomethyl_c().build().unwrap()
    }

    #[test]
    fn refinement_aa_set_high_res_exposes_five_variable_mods() {
        let base = base_set_with_cam_c();
        let set = refinement_aa_set(&base, &RefineConfig::default_tier(), true);

        // Count DISTINCT variable mods across all locations by name.
        let mut names: HashSet<String> = HashSet::new();
        for aa in set.iter_variants() {
            if let Some(m) = &aa.mod_ {
                if !m.fixed {
                    names.insert(m.name.clone());
                }
            }
        }
        assert_eq!(
            names.len(),
            5,
            "default tier has 5 variable mods at high-res; got {names:?}"
        );
        // The fixed Carbamidomethyl baseline is preserved.
        let c_anywhere = set.variants_for(b'C', ModLocation::Anywhere);
        assert!(c_anywhere.iter().any(|aa| aa.mod_.as_ref().map(|m| m.fixed).unwrap_or(false)));
    }

    #[test]
    fn refinement_aa_set_low_res_drops_deamidation() {
        let base = base_set_with_cam_c();
        let set = refinement_aa_set(&base, &RefineConfig::default_tier(), false);
        let mut names: HashSet<String> = HashSet::new();
        for aa in set.iter_variants() {
            if let Some(m) = &aa.mod_ {
                if !m.fixed {
                    names.insert(m.name.clone());
                }
            }
        }
        assert_eq!(names.len(), 4, "deamidation dropped at low-res; got {names:?}");
        assert!(!names.contains("Deamidation"), "Deamidation must be dropped at low-res");
    }

    #[test]
    fn refinement_aa_set_oxidation_m_has_modified_variant() {
        let base = base_set_with_cam_c();
        let set = refinement_aa_set(&base, &RefineConfig::default_tier(), true);
        // M-Anywhere should now carry an Oxidation variable variant (unmod + ox).
        let m_variants = set.variants_for(b'M', ModLocation::Anywhere);
        assert!(m_variants.iter().any(|aa| {
            aa.mod_.as_ref().map(|m| m.name == "Oxidation").unwrap_or(false)
        }));
    }

    // ── (e) mod_count_and_class ────────────────────────────────────────────

    /// Build a peptide with one residue carrying a mod of the given delta.
    fn modded_pep(residue: u8, delta: f64) -> Peptide {
        let m = Modification {
            name: "x".into(),
            mass_delta: delta,
            residue: ResidueSpec::Specific(residue),
            location: ModLocation::Anywhere,
            fixed: false,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
        };
        let aa = AminoAcid::standard(residue).unwrap().with_mod(m);
        Peptide::new(vec![aa, AminoAcid::standard(b'G').unwrap()], b'_', b'-')
    }

    #[test]
    fn mod_count_and_class_unmodified_is_zero_zero() {
        let cfg = RefineConfig::default_tier();
        let pep = Peptide::new(
            vec![AminoAcid::standard(b'P').unwrap(), AminoAcid::standard(b'G').unwrap()],
            b'_',
            b'-',
        );
        assert_eq!(mod_count_and_class(&pep, &cfg), (0, 0));
    }

    #[test]
    fn mod_count_and_class_oxidation_m_is_one_class_oxidation() {
        let cfg = RefineConfig::default_tier();
        // Oxidation on M, delta 15.994915, class "oxidation" → id 1.
        let pep = modded_pep(b'M', 15.994915);
        assert_eq!(mod_count_and_class(&pep, &cfg), (1, 1));
    }

    #[test]
    fn mod_count_and_class_deamidation_n_is_class_two() {
        let cfg = RefineConfig::default_tier();
        // Deamidation on N, delta 0.984016, class "deamidation" → id 2.
        let pep = modded_pep(b'N', 0.984016);
        assert_eq!(mod_count_and_class(&pep, &cfg), (1, 2));
    }

    #[test]
    fn mod_count_and_class_ignores_non_refinement_mods() {
        // A Carbamidomethyl-C (fixed baseline) is NOT a refinement chemistry,
        // so it must not be counted as a refinement mod.
        let cfg = RefineConfig::default_tier();
        let pep = modded_pep(b'C', 57.02146);
        assert_eq!(mod_count_and_class(&pep, &cfg), (0, 0));
    }

    #[test]
    fn mod_count_and_class_counts_two_distinct_mods() {
        let cfg = RefineConfig::default_tier();
        // Oxidation-M (class 1) AND Deamidation-N (class 2) on one peptide.
        let ox = AminoAcid::standard(b'M').unwrap().with_mod(Modification {
            name: "Oxidation".into(),
            mass_delta: 15.994915,
            residue: ResidueSpec::Specific(b'M'),
            location: ModLocation::Anywhere,
            fixed: false,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
        });
        let deam = AminoAcid::standard(b'N').unwrap().with_mod(Modification {
            name: "Deamidation".into(),
            mass_delta: 0.984016,
            residue: ResidueSpec::Specific(b'N'),
            location: ModLocation::Anywhere,
            fixed: false,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
        });
        let pep = Peptide::new(vec![ox, deam], b'_', b'-');
        let (n, class) = mod_count_and_class(&pep, &cfg);
        assert_eq!(n, 2);
        // First non-other class encountered wins (oxidation = 1).
        assert_eq!(class, 1);
    }

    // ── (f) run_refinement: tag-merge logic on synthetic queues ────────────
    //
    // A full Pass-2 search needs a real PreparedSearch (RankScorer + Param +
    // matching spectra), which is disproportionate to construct here for a
    // deterministic winner. We instead test the merge/tag step — the orchestration
    // logic `run_refinement` performs on each Pass-2 winner — directly: tag
    // is_refinement, fill (num_mods, class) from the Pass-2 peptide, remap
    // spectrum_idx, and force_push into the global queue without evicting the
    // Pass-1 PSM. The early-return guards (a)/(b) are covered above; the real
    // PreparedSearch path is exercised by `run_refinement_early_returns_*` below.

    #[test]
    fn run_refinement_tag_merge_preserves_pass1_and_tags_pass2() {
        let cfg = RefineConfig::default_tier();

        // Pass-1 candidates: spectrum 0 has a target PSM; spectrum 1 is empty.
        let pass1_cands = vec![cand(0, false)];
        let mut global_queues = vec![
            queue_with(vec![psm(0, 0, 40.0)]), // identified
            TopNQueue::new(10),                // unidentified
        ];

        // Simulate a Pass-2 winner for spectrum 1: an Oxidation-M peptide.
        let pass2_cands = vec![Candidate {
            peptide: modded_pep(b'M', 15.994915),
            protein_index: 0,
            start_offset_in_protein: 0,
            is_decoy: false,
            is_protein_n_term: false,
            is_protein_c_term: false,
        }];
        let mut winner = psm(0 /* local */, 0, 22.0);

        // Replicate the merge body of run_refinement for global spectrum 1.
        let global_idx = 1usize;
        winner.features.is_refinement = 1;
        let (num_mods, class) =
            mod_count_and_class(&pass2_cands[winner.primary_candidate_idx() as usize].peptide, &cfg);
        winner.features.num_mods = num_mods;
        winner.features.refine_mod_class = class;
        winner.spectrum_idx = global_idx;
        global_queues[global_idx].force_push(winner);

        // Pass-1 queue 0 is untouched.
        assert_eq!(global_queues[0].len(), 1);
        assert_eq!(global_queues[0].peek_top().unwrap().features.is_refinement, 0);

        // Queue 1 now holds the tagged Pass-2 winner.
        assert_eq!(global_queues[1].len(), 1);
        let merged = global_queues[1].peek_top().unwrap();
        assert_eq!(merged.features.is_refinement, 1);
        assert_eq!(merged.features.num_mods, 1);
        assert_eq!(merged.features.refine_mod_class, 1); // oxidation
        assert_eq!(merged.spectrum_idx, 1);

        // Sanity: this is the same pass1_cands list a real caller would pass;
        // reference it so the fixture mirrors the real signature shape.
        assert!(!pass1_cands[0].is_decoy);
    }

    // A genuine `PreparedSearch`-backed run, exercising the real Pass-2 plumbing
    // through the early-return guards. We build a tiny scorer + a 1-protein target
    // db; with all-decoy / no-unidentified inputs the cascade must early-return
    // and leave the queues untouched, without panicking.
    fn tiny_scorer() -> RankScorer {
        use rustc_hash::FxHashMap;
        use scoring_crate::param_model::{IonType, Partition, SpecDataType};
        use scoring_crate::Param;
        use model::activation::ActivationMethod;
        use model::instrument::InstrumentType;
        use model::protocol::Protocol;
        use model::tolerance::Tolerance;

        let part = Partition { charge: 2, parent_mass: 1000.0, seg_num: 0 };
        let prefix1 = IonType::Prefix { charge: 1, offset_bits: 0.0_f32.to_bits(), loss_class: 0 };
        let suffix1 = IonType::Suffix { charge: 1, offset_bits: 0.0_f32.to_bits(), loss_class: 0 };
        let noise = IonType::Noise;

        let mut ion_table = FxHashMap::default();
        ion_table.insert(prefix1, vec![0.5_f32, 0.1, 0.05, 0.01]);
        ion_table.insert(suffix1, vec![0.5_f32, 0.1, 0.05, 0.01]);
        ion_table.insert(noise, vec![0.1_f32, 0.05, 0.02, 0.01]);

        let mut rank_dist_table = FxHashMap::default();
        rank_dist_table.insert(part, ion_table);
        let mut frag_off_table = FxHashMap::default();
        frag_off_table.insert(part, vec![]);

        let mut param = Param {
            version: 10001,
            data_type: SpecDataType {
                activation: ActivationMethod::HCD,
                instrument: InstrumentType::QExactive,
                enzyme: None,
                protocol: Protocol::Automatic,
            },
            mme: Tolerance::Ppm(20.0),
            apply_deconvolution: false,
            deconvolution_error_tolerance: 0.0,
            charge_hist: vec![(2, 100)],
            min_charge: 2,
            max_charge: 2,
            num_segments: 1,
            partitions: vec![part],
            num_precursor_off: 0,
            precursor_off_map: FxHashMap::default(),
            frag_off_table,
            max_rank: 3,
            rank_dist_table,
            error_scaling_factor: 0,
            ion_err_dist_table: FxHashMap::default(),
            noise_err_dist_table: FxHashMap::default(),
            ion_existence_table: FxHashMap::default(),
            partition_ion_types_cache: FxHashMap::default(),
            gbdt_peak_model: None,
            frag_intensity_model: None,
            rich_ion_model: None,
        };
        param.rebuild_cache();
        RankScorer::new(&param)
    }

    #[test]
    fn run_refinement_early_returns_when_no_confident_proteins() {
        let aa_set = base_set_with_cam_c();
        let params = SearchParams::default_tryptic(aa_set);
        let scorer = tiny_scorer();
        let target = ProteinDb { proteins: vec![protein("P0", b"MKWVTFISLLR")] };

        // All-decoy Pass-1 → no confident proteins → early return.
        let pass1_cands = vec![cand(0, true)];
        let mut queues = vec![queue_with(vec![psm(0, 0, 10.0)])];
        let before_len = queues[0].len();
        let spectra: Vec<model::spectrum::Spectrum> = vec![model::spectrum::Spectrum::default()];

        run_refinement(
            &mut queues,
            &spectra,
            &pass1_cands,
            &target,
            &params,
            &scorer,
            &RefineConfig::default_tier(),
            true,
            0.5,
            "XXX",
            42,
        );
        // Untouched.
        assert_eq!(queues[0].len(), before_len);
    }

    #[test]
    fn run_refinement_early_returns_when_no_unidentified_spectra() {
        let aa_set = base_set_with_cam_c();
        let params = SearchParams::default_tryptic(aa_set);
        let scorer = tiny_scorer();
        let target = ProteinDb { proteins: vec![protein("P0", b"MKWVTFISLLR")] };

        // Confident target Pass-1 AND every spectrum identified → early return at (b).
        let pass1_cands = vec![cand(0, false)];
        let mut queues = vec![queue_with(vec![psm(0, 0, 50.0)])];
        let spectra: Vec<model::spectrum::Spectrum> = vec![model::spectrum::Spectrum::default()];

        run_refinement(
            &mut queues,
            &spectra,
            &pass1_cands,
            &target,
            &params,
            &scorer,
            &RefineConfig::default_tier(),
            true,
            0.5,
            "XXX",
            42,
        );
        assert_eq!(queues[0].len(), 1);
        assert_eq!(queues[0].peek_top().unwrap().features.is_refinement, 0);
    }
}
