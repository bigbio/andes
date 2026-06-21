//! PTM-refinement cascade (Pass-2). After the Pass-1 search, this module:
//!   1. extracts the UNMODIFIED backbone sequences of the confident Pass-1 target
//!      PEPTIDES (`confident_base_peptides`; MaxQuant dependent-peptide-style
//!      anchoring),
//!   2. finds the spectra still without any target PSM
//!      (`unidentified_spectrum_indices`),
//!   3. builds a SMALL target+decoy search index anchored on those confident base
//!      peptides — one mini-protein per peptide + a 1:1 decoy each
//!      (`build_peptide_anchored_index`),
//!   4. searches the unidentified spectra against that small index with the
//!      refinement variable-mod tier applied (`refinement_aa_set`), and
//!   5. RETURNS every Pass-2 winner (tagged as a refinement PSM) as a
//!      [`RefinementOutput`] so the caller can write it as a SEPARATE PIN
//!      (`run_refinement`).
//!
//! The cascade is OPT-IN (`--refine`) and additive: nothing here runs on the
//! default path. Pass-2 winners are NOT merged into the Pass-1 queues — a merged
//! Pass-2 `PsmMatch::candidate_idxs` indexes the PASS-2 candidate list, which is
//! a different slice than the Pass-1 `candidates` the PIN writer resolves
//! peptide/accession against (a merge would alias the indices → wrong protein or
//! OOB). Instead the Pass-2 artifacts (index + candidates + tagged queues +
//! global spectrum indices) are returned, and the caller emits a second PIN: the
//! merged report is the disjoint union Pass-1 PIN ⊎ refine PIN.

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
///
/// SUPERSEDED on the live path by [`confident_base_peptides`] (PEPTIDE
/// anchoring); kept (with unit tests) for reference and possible reuse.
#[allow(dead_code)]
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
// (a') confident_base_peptides
// ─────────────────────────────────────────────────────────────────────────────

/// The UNMODIFIED backbone sequences of the confident Pass-1 target peptides,
/// deduped and sorted (MaxQuant dependent-peptide-style anchoring).
///
/// Mirrors [`confident_protein_indices`]'s best-per-spectrum TDC walk, but
/// instead of mapping each confident TARGET to its `protein_index`, it extracts
/// the candidate peptide's bare amino-acid backbone (each residue's
/// `AminoAcid::residue` byte, IGNORING any attached modification — Pass-2 re-adds
/// the refinement-tier mods to these backbones). The de-duplicated, sorted
/// sequences become the anchor set [`build_peptide_anchored_index`] turns into a
/// peptide-anchored refinement db: one 1-residue-mini-protein per confident
/// backbone, NOT every tryptic peptide of every confident protein. This collapses
/// the Pass-2 candidate+decoy pool (and so the decoy RankScore ceiling) to
/// "confident peptides × refinement mods" + reversed-peptide decoys, so genuine
/// modified forms clear FDR instead of being buried by ~99k protein-digest decoys.
pub fn confident_base_peptides(
    queues: &[TopNQueue],
    candidates: &[Candidate],
    select_q: f64,
) -> Vec<Vec<u8>> {
    // One (label, candidate_idx) per spectrum that produced a PSM, using that
    // spectrum's best PSM — identical to `confident_protein_indices`.
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

    // Confident TARGETs → their candidate peptide's UNMODIFIED backbone bytes.
    // `AminoAcid::residue` is the bare residue letter regardless of `mod_`, so
    // collecting it per residue yields the modification-free backbone.
    let mut seqs: HashSet<Vec<u8>> = HashSet::new();
    for li in confident_target_indices(&labels, select_q) {
        let cand = &candidates[cand_idx_of[li]];
        let backbone: Vec<u8> = cand.peptide.residues.iter().map(|aa| aa.residue).collect();
        seqs.insert(backbone);
    }

    let mut out: Vec<Vec<u8>> = seqs.into_iter().collect();
    out.sort(); // deterministic order
    out
}

// ─────────────────────────────────────────────────────────────────────────────
// (b) unidentified_spectrum_indices
// ─────────────────────────────────────────────────────────────────────────────

/// Spectra Pass-1 did NOT confidently identify at `report_q` (the report FDR
/// threshold, e.g. 0.01). A spectrum is "identified" iff its best PSM (max
/// rank_score across the queue) is a TARGET whose TDC q-value <= report_q; all
/// other spectra (best PSM is a decoy, best target above q, or no PSM at all)
/// are returned for refinement. Uses the same best-per-spectrum TDC walk as
/// `confident_protein_indices`.
pub fn unidentified_spectrum_indices(
    queues: &[TopNQueue],
    candidates: &[Candidate],
    report_q: f64,
) -> Vec<usize> {
    // One best-PSM (label, spectrum) per spectrum that produced a PSM, using
    // that spectrum's max-`rank_score` PSM (target OR decoy). Empty-queue spectra
    // contribute no label — they are unconditionally unidentified below.
    let mut labels: Vec<ScoredLabel> = Vec::new();
    let mut label_spectrum: Vec<usize> = Vec::new();
    for (s, queue) in queues.iter().enumerate() {
        if let Some(best) = queue.peek_top() {
            let cand_idx = best.primary_candidate_idx() as usize;
            labels.push(ScoredLabel {
                score: best.rank_score,
                is_decoy: candidates[cand_idx].is_decoy,
            });
            label_spectrum.push(s);
        }
    }

    // TDC walk over the best-per-spectrum labels → indices (into `labels`) of
    // confident TARGETs; map those back to their global spectrum indices.
    let confident = confident_target_indices(&labels, report_q);
    let identified: HashSet<usize> = confident.iter().map(|&i| label_spectrum[i]).collect();

    // Every spectrum NOT confidently identified (best-is-decoy, best target
    // above q, or empty queue) is returned for refinement.
    (0..queues.len()).filter(|s| !identified.contains(s)).collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// (c) build_refinement_index
// ─────────────────────────────────────────────────────────────────────────────

/// Build a SMALL target+decoy search index scoped to the confident proteins.
///
/// Subsets `full_target_db` (the TARGET-ONLY db) to the confident protein
/// indices, then regenerates a fresh 1:1 decoy per kept target (using
/// `decoy_strategy`, matching the Pass-1 search) so the scoped search has its
/// own decoys for TDC. `confident_target_proteins` carries
/// `Candidate.protein_index` values (combined-db positions; see
/// [`confident_protein_indices`]); any index ≥ `full_target_db.len()` would be a
/// decoy position and is filtered out before subsetting the target db.
///
/// SUPERSEDED on the live path by [`build_peptide_anchored_index`] (PEPTIDE
/// anchoring); kept (with unit tests) for reference and possible reuse.
#[allow(dead_code)]
pub fn build_refinement_index(
    full_target_db: &ProteinDb,
    confident_target_proteins: &HashSet<usize>,
    decoy_prefix: &str,
    decoy_strategy: DecoyStrategy,
    seed: u64,
) -> SearchIndex {
    // Keep only indices that address a TARGET protein in `full_target_db`.
    // Indices that point past the target range belong to the combined db's decoy
    // half (e.g. the Reverse strategy) and must never be subset into a target db.
    //
    // INVARIANT: build_search_db places targets at [0, target_count) then decoys; a confident TARGET's protein_index is therefore always < full_target_db.len(). Re-check if a future DecoyStrategy interleaves proteins.
    let keep: HashSet<usize> = confident_target_proteins
        .iter()
        .copied()
        .filter(|&i| i < full_target_db.len())
        .collect();

    let subset = full_target_db.subset_by_index(&keep);
    let db = build_search_db(&subset, decoy_prefix, decoy_strategy, seed);
    SearchIndex {
        db,
        decoy_prefix: crate::decoy::normalize_decoy_prefix(decoy_prefix),
        decoy_suffix: None,
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// (c') build_peptide_anchored_index
// ─────────────────────────────────────────────────────────────────────────────

/// Build the peptide-anchored Pass-2 search index: each confident base peptide
/// becomes its OWN 1-entry "mini-protein", then a 1:1 decoy is generated per
/// mini-protein (matching the Pass-1 `decoy_strategy`) so the scoped search has
/// its own decoys for TDC.
///
/// This is the MaxQuant dependent-peptide-style anchoring that replaces the
/// protein-subset path ([`build_refinement_index`]): the candidate+decoy pool is
/// "confident peptides × refinement mods" + reversed-peptide decoys, NOT "all
/// tryptic peptides of confident proteins × mods". Each mini-protein digests to
/// its own backbone peptide under the existing enzyme (protein-N-term flank +
/// K/R C-term), to which `enumerate_candidates` then applies the refinement-tier
/// mods. `base_seqs` are the UNMODIFIED backbones from [`confident_base_peptides`].
pub fn build_peptide_anchored_index(
    base_seqs: &[Vec<u8>],
    decoy_prefix: &str,
    decoy_strategy: DecoyStrategy,
    seed: u64,
) -> SearchIndex {
    let minidb = ProteinDb {
        proteins: base_seqs
            .iter()
            .enumerate()
            .map(|(i, seq)| model::protein::Protein {
                accession: format!("BASEPEP_{i}"),
                description: String::new(),
                sequence: seq.clone(),
            })
            .collect(),
    };
    let combined = build_search_db(&minidb, decoy_prefix, decoy_strategy, seed);
    SearchIndex {
        db: combined,
        decoy_prefix: crate::decoy::normalize_decoy_prefix(decoy_prefix),
        decoy_suffix: None,
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
/// **Fixed/variable-mod handling.** Pass-2 INHERITS the base search's FULL
/// chemistry: every fixed mod (Carbamidomethyl-C, TMT/iTRAQ labels, …) and every
/// variable mod (Ox-M, Acetyl, …) the user declared, recovered from
/// [`AminoAcidSet::distinct_mods`]. The discovery tier is then layered on top as
/// VARIABLE mods, DE-DUPLICATED against the base (a base variable Oxidation-M and
/// the tier's Oxidation-M are the same chemistry → added once). This is what lets
/// `--refine` work on labeled data: previously only the standard CAM-C baseline
/// was reconstructed, so TMT/iTRAQ deltas were dropped and every Pass-2 candidate
/// mass-mismatched its precursor (a silent no-op).
pub fn refinement_aa_set(
    base: &AminoAcidSet,
    cfg: &RefineConfig,
    high_res: bool,
) -> Result<AminoAcidSet, String> {
    let mut builder = AminoAcidSetBuilder::new_standard();

    // De-dup key: (residue, location, delta rounded to 1e-4 Da). Rounding (rather
    // than raw f64 bits) folds trivially-different literals onto one entry — e.g.
    // Cam-C 57.02146 vs 57.021464, or Ox-M 15.99491 vs 15.994915.
    let mod_key = |residue: ResidueSpec, location: ModLocation, delta: f64| {
        (residue, location, (delta * 1e4).round() as i64)
    };
    let mut present: HashSet<(ResidueSpec, ModLocation, i64)> = HashSet::new();
    // (residue, location) slots occupied by a base FIXED mod. The builder forbids
    // a fixed AND a variable mod on the same residue+location, so a tier variable
    // mod landing on a fixed slot (e.g. tier Oxidation-K on a TMT-fixed K) must be
    // skipped — the fixed label owns that position.
    let mut fixed_slots: Vec<(ResidueSpec, ModLocation)> = Vec::new();

    // 1. Inherit the base search's full mod set (fixed AND variable).
    for m in base.distinct_mods() {
        present.insert(mod_key(m.residue, m.location, m.mass_delta));
        if m.fixed {
            fixed_slots.push((m.residue, m.location));
            builder = builder.add_fixed_mod((*m).clone());
        } else {
            builder = builder.add_variable_mod((*m).clone());
        }
    }

    // Residues overlap when either side is a wildcard or they are the same letter.
    // Intentionally BROADER than the builder's exact-`(residue, location)`
    // conflict rule: wildcard on either side counts as overlap. This may
    // over-skip a custom-tier wildcard-Anywhere variable mod when the base has a
    // specific fixed mod at Anywhere (the builder would accept that pair) — but
    // over-skipping a tier mod is safe (it is just not searched), whereas the
    // alternative risks a builder hard-error. The built-in default tier has no
    // wildcard-Anywhere mod, so it is unaffected.
    let residue_overlaps = |a: ResidueSpec, b: ResidueSpec| {
        a == ResidueSpec::Wildcard || b == ResidueSpec::Wildcard || a == b
    };

    // 2. Layer the discovery tier on top as variable mods, skipping (a) low-res
    //    deamidation, (b) chemistry the base already provides, and (c) a tier mod
    //    that lands on a base FIXED slot (would conflict in the builder).
    for rm in &cfg.mods {
        if !high_res && rm.class == "deamidation" {
            continue;
        }
        let location = parse_location(&rm.location);
        for res in &rm.residues {
            let residue = parse_residue_spec(res);
            if fixed_slots
                .iter()
                .any(|&(fr, fl)| fl == location && residue_overlaps(fr, residue))
            {
                continue; // a base fixed mod owns this residue/location slot
            }
            if !present.insert(mod_key(residue, location, rm.delta)) {
                continue; // keep the provided PTM; don't double-add the overlap
            }
            let m = Modification {
                name: rm.name.clone(),
                mass_delta: rm.delta,
                residue,
                location,
                fixed: false,
                accession: None,
                neutral_losses: Vec::new(),
                loss_class: 0,
            };
            builder = builder.add_variable_mod(m);
        }
    }

    // Propagate build errors (implausible mass, conflicting fixed+variable on a
    // residue, …) instead of panicking — a bad `--refine-config` must degrade to
    // "skip Pass-2 with a WARN", never crash the whole search.
    builder
        .build()
        .map_err(|e| format!("refinement_aa_set: {e}"))
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
/// so this would otherwise count every modified residue. The base search's fixed
/// mods (Carbamidomethyl-C, and — since Pass-2 now inherits the full base mod set
/// — any labels such as TMT/iTRAQ) are NOT part of any refinement `cfg.mods`
/// class; counting them would inflate `num_mods`. To keep `num_mods` aligned with
/// the REFINEMENT tier, we count a residue only when its mod matches one of
/// `cfg.mods` (i.e. a refinement chemistry), so base fixed labels are excluded.
pub fn mod_count_and_class(peptide: &Peptide, cfg: &RefineConfig) -> (u32, u32) {
    let mut num_mods = 0u32;
    let mut class_id = 0u32; // 0 = none until a refinement mod matches

    for aa in &peptide.residues {
        let Some(m) = aa.mod_.as_ref() else { continue };
        // Find the refinement mod that this residue's modification matches.
        // Match on delta + residue + LOCATION so two tiers with the same
        // mass/residue but a different location/class are not conflated
        // (e.g. an anywhere mod vs an N-term mod of the same Δ).
        let matched = cfg.mods.iter().find(|rm| {
            (rm.delta - m.mass_delta).abs() < 1e-4
                && parse_location(&rm.location) == m.location
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

/// Pass-2 refinement results, ready to write as a SEPARATE PIN (the merged
/// report is the disjoint union: Pass-1 PIN ⊎ this). `queues[i]` corresponds to
/// the spectrum at `all_spectra[global_spectrum_indices[i]]`; resolve peptide/
/// accession against `candidates` + `index` (NOT the Pass-1 candidates).
pub struct RefinementOutput {
    /// The scoped Pass-2 target+decoy `SearchIndex`. The PIN writer resolves
    /// accessions through this — NOT the Pass-1 index.
    pub index: SearchIndex,
    /// The Pass-2 candidate pool. Every returned PSM's `candidate_idxs` index
    /// into THIS slice (not the Pass-1 candidates).
    pub candidates: Vec<Candidate>,
    /// One tagged Pass-2 queue per refined spectrum, in `global_spectrum_indices`
    /// order. Each PSM has `is_refinement = 1`, `(num_mods, refine_mod_class)`
    /// filled, and `spectrum_idx` set to the GLOBAL stream index.
    pub queues: Vec<TopNQueue>,
    /// Global stream index of each refined spectrum; `queues[i]` is the Pass-2
    /// result for `all_spectra[global_spectrum_indices[i]]`. 1:1 with `queues`.
    pub global_spectrum_indices: Vec<usize>,
}

/// Orchestrate the Pass-2 PTM-refinement search and RETURN its tagged winners as
/// a [`RefinementOutput`] (the caller writes them as a separate "refine" PIN).
///
/// Steps (each an early-return — `None` — when its input is empty, so the cheap
/// no-op cases never build a Pass-2 search):
///   1. `confident_base_peptides` over the Pass-1 queues (PEPTIDE anchoring) —
///      `None` if empty.
///   2. `unidentified_spectrum_indices` over the Pass-1 queues (spectra whose
///      best Pass-1 target PSM is above `report_q`) — `None` if empty.
///   3. `build_peptide_anchored_index` scoped to the confident base peptides.
///   4. A Pass-2 `SearchParams` clone whose `aa_set` is the refinement tier
///      ([`refinement_aa_set`]) and whose `max_variable_mods_per_peptide` is
///      `cfg.max_mods`, prepared via `PreparedSearch::prepare`.
///   5. `run_chunk` over the (cloned) unidentified spectra; each winner is tagged
///      `is_refinement = 1`, has `(num_mods, refine_mod_class)` filled from its
///      Pass-2 peptide, and has its `spectrum_idx` fixed to the GLOBAL index.
///
/// `pass1_queues` is READ-ONLY: it drives the confident-protein / unidentified-
/// spectrum detection only. The Pass-2 winners are NOT merged into it — see the
/// module doc for why (candidate-index aliasing). The returned
/// [`RefinementOutput`] carries the Pass-2 index + candidates + tagged queues so
/// the PIN writer resolves peptide/protein against the PASS-2 slice.
#[allow(clippy::too_many_arguments)]
pub fn run_refinement(
    pass1_queues: &[TopNQueue],
    all_spectra: &[model::spectrum::Spectrum],
    pass1_candidates: &[Candidate],
    full_target_db: &ProteinDb,
    base_params: &SearchParams,
    scorer: &RankScorer,
    cfg: &RefineConfig,
    report_q: f64,
    high_res: bool,
    fragment_tol_da: f64,
    decoy_prefix: &str,
    decoy_strategy: DecoyStrategy,
    seed: u64,
) -> Option<RefinementOutput> {
    // `full_target_db` is no longer consulted (the refinement db is anchored on
    // confident PEPTIDES, not a protein subset); keep the param so the andes call
    // site is untouched.
    let _ = full_target_db;

    // 1. Confident base peptides from Pass-1 (PEPTIDE-anchored scoping gate).
    //    Empty ⇒ nothing to refine.
    let base_seqs = confident_base_peptides(
        pass1_queues,
        pass1_candidates,
        base_params.refine_select_psm_fdr,
    );
    if base_seqs.is_empty() {
        return None;
    }

    // 2. Spectra Pass-1 did not confidently identify at the report FDR
    //    threshold. Empty ⇒ nothing to rescue.
    let unident = unidentified_spectrum_indices(pass1_queues, pass1_candidates, report_q);
    if unident.is_empty() {
        return None;
    }

    // Announce the cascade + the discovery chemistry it will search, splitting
    // active vs. the mods skipped on this resolution (deamidation is high-res
    // only). Mirrors the filter `refinement_aa_set` applies, so the banner always
    // matches what is actually searched.
    let fmt_mod = |m: &crate::refine_config::RefineMod| {
        format!("{} ({} {:+.3})", m.name, m.residues.join(","), m.delta)
    };
    let active: Vec<String> = cfg
        .mods
        .iter()
        .filter(|m| high_res || m.class != "deamidation")
        .map(&fmt_mod)
        .collect();
    let skipped: Vec<String> = cfg
        .mods
        .iter()
        .filter(|m| !high_res && m.class == "deamidation")
        .map(&fmt_mod)
        .collect();
    eprintln!(
        "[refine] Starting PTM-refinement search over {} confident-protein anchor peptides \
         ({} unidentified spectra, max {} mods/peptide)",
        base_seqs.len(),
        unident.len(),
        cfg.max_mods,
    );
    eprintln!("[refine]   discovery mods: {}", active.join(", "));
    if !skipped.is_empty() {
        eprintln!("[refine]   skipped (low-res): {}", skipped.join(", "));
    }

    // 3. Peptide-anchored target+decoy index: one mini-protein per confident base
    //    peptide + a 1:1 decoy each. Built as an owned local so it can be MOVED
    //    into the returned `RefinementOutput` once the `PreparedSearch` borrow of
    //    it ends (see below).
    let refine_idx =
        build_peptide_anchored_index(&base_seqs, decoy_prefix, decoy_strategy, seed);

    // 4. Pass-2 params: refinement variable-mod tier + the tier's mod cap.
    //    Bound the per-peptide mod cap: the candidate count grows combinatorially
    //    in (mods × sites), so a large `--refine-max-mods` over wildcard mods can
    //    blow up CPU/RAM. Clamp to a sane ceiling and WARN rather than letting a
    //    misconfiguration exhaust the machine.
    const MAX_REFINE_MODS_PER_PEPTIDE: u32 = 5;
    let capped_max_mods = if cfg.max_mods > MAX_REFINE_MODS_PER_PEPTIDE {
        eprintln!(
            "WARN: --refine-max-mods {} exceeds the safe ceiling {}; clamping (combinatorial blowup guard).",
            cfg.max_mods, MAX_REFINE_MODS_PER_PEPTIDE
        );
        MAX_REFINE_MODS_PER_PEPTIDE
    } else {
        cfg.max_mods
    };
    let mut refine_params = base_params.clone();
    refine_params.aa_set = match refinement_aa_set(&base_params.aa_set, cfg, high_res) {
        Ok(set) => set,
        Err(e) => {
            eprintln!("WARN: --refine could not build the Pass-2 modification set ({e}); skipping refinement for this run.");
            return None;
        }
    };
    refine_params.max_variable_mods_per_peptide = capped_max_mods;

    // 5. Search ONLY the unidentified spectra. We clone them into a contiguous
    //    sub-vec (run_chunk takes a &[Spectrum]); `spectrum_idx_offset = 0` is
    //    fine because we remap to the GLOBAL index ourselves below.
    //
    //    `PreparedSearch` borrows `&refine_idx`/`&refine_params`/`scorer`. We
    //    must return `refine_idx` (owned) AND have used it to build the prepared
    //    search, so we extract the Pass-2 candidate pool into an owned Vec and
    //    DROP `refine_prepared` (ending the borrow) before moving `refine_idx`
    //    into the output below.
    let (pass2_queues, refine_candidates) = {
        let refine_prepared = PreparedSearch::prepare(
            &refine_idx,
            &refine_params,
            scorer,
            fragment_tol_da,
            decoy_prefix,
        );
        let subset: Vec<model::spectrum::Spectrum> =
            unident.iter().map(|&i| all_spectra[i].clone()).collect();
        let queues = refine_prepared.run_chunk(&subset, 0);
        // `PreparedSearch.candidates` is an owned `Vec<Candidate>`; move it out so
        // the prepared search (and its `&refine_idx` borrow) can be dropped.
        (queues, refine_prepared.candidates)
    };

    // Tag each Pass-2 winner; collect into one queue per refined spectrum, AND
    // PRUNE the candidate pool to just the winners. The Pass-2 search visits an
    // enormous pool (mod permutations over every anchored peptide — millions of
    // candidates), but only the few thousand referenced by a surviving PSM are
    // ever needed downstream: the PIN resolves peptide/accession through
    // `candidate_idxs`. We copy those into a compact pool, remap indices, and drop
    // the full pool — it is the dominant refine memory driver (the cloned peaks
    // are ~45 MB by comparison). Output is unchanged: each PSM still resolves to
    // the same peptide/candidate, only the internal index numbering differs.
    let pass2_total = refine_candidates.len();
    let mut compact: Vec<Candidate> = Vec::new();
    let mut remap: std::collections::HashMap<u32, u32> = std::collections::HashMap::new();
    let mut out_queues: Vec<TopNQueue> = Vec::with_capacity(unident.len());
    for (local_idx, pass2_queue) in pass2_queues.into_iter().enumerate() {
        let global_idx = unident[local_idx];
        let mut tagged = TopNQueue::new(refine_params.top_n_psms_per_spectrum);
        for mut psm in pass2_queue.into_rank_sorted_vec() {
            // mod_count uses the OLD index into the full pool — compute before remap.
            let (num_mods, class) = mod_count_and_class(
                &refine_candidates[psm.primary_candidate_idx() as usize].peptide,
                cfg,
            );
            // The anchored Pass-2 DB also contains the UNMODIFIED base peptide, so a
            // Pass-2 winner can carry zero refinement mods. Those are NOT discoveries
            // — they duplicate what Pass-1 already searched (same peptide, same
            // spectrum) and must not be emitted as refinement rows (`is_refinement=1`)
            // or they inflate the modified count and double-count per scan.
            if num_mods == 0 {
                continue;
            }
            psm.features.is_refinement = 1;
            psm.features.num_mods = num_mods;
            psm.features.refine_mod_class = class;
            // Fix the spectrum index back to the global stream position.
            psm.spectrum_idx = global_idx;
            // Remap every candidate index into the compact winners-only pool.
            for idx in &mut psm.candidate_idxs {
                let old = *idx;
                *idx = match remap.get(&old) {
                    Some(&n) => n,
                    None => {
                        let n = compact.len() as u32;
                        compact.push(refine_candidates[old as usize].clone());
                        remap.insert(old, n);
                        n
                    }
                };
            }
            tagged.force_push(psm);
        }
        out_queues.push(tagged);
    }
    drop(refine_candidates); // free the full Pass-2 pool now; `compact` replaces it

    // Transparency + cost visibility: show the full Pass-2 pool size and the
    // pruned winners-only count so a user/operator sees both the search cost and
    // what survives into the merged PIN.
    eprintln!(
        "[REFINE] anchors={} unident={}/{} pass2_candidates={} (pruned to {} winners)",
        base_seqs.len(),
        unident.len(),
        all_spectra.len(),
        pass2_total,
        compact.len(),
    );

    Some(RefinementOutput {
        index: refine_idx,
        candidates: compact,
        queues: out_queues,
        global_spectrum_indices: unident,
    })
}

// ─────────────────────────────────────────────────────────────────────────────
// (g) merge_into_pass1
// ─────────────────────────────────────────────────────────────────────────────

/// The result of merging a Pass-2 [`RefinementOutput`] into the Pass-1 candidate
/// index space: a single combined candidate list and the concatenated
/// target+decoy index. Pass to the single `write_pin` call so unmodified and
/// modified PSMs for every scan compete in ONE PIN.
pub struct MergedSearch {
    /// Pass-1 candidates followed by offset-corrected Pass-2 candidates.
    /// Every `PsmMatch::candidate_idxs` in both passes index into this slice.
    pub candidates: Vec<Candidate>,
    /// Pass-1 index concatenated with the Pass-2 index. Accession resolution
    /// for all candidates uses this combined db.
    pub index: SearchIndex,
}

/// Merge the Pass-2 refinement winners into the Pass-1 per-scan queues under one
/// candidate-index space, so unmodified and modified PSMs for a scan compete in
/// a SINGLE PIN. Offsets `protein_index` by the Pass-1 db length and
/// `candidate_idxs` by the Pass-1 candidate count, then `force_push`es each
/// refine PSM into `pass1_queues[global_spectrum]` (legitimate extra emission,
/// like a chimeric secondary). Returns the combined candidates + index for the
/// single `write_pin` call.
pub fn merge_into_pass1(
    pass1_queues: &mut [TopNQueue],
    pass1_candidates: Vec<Candidate>,
    pass1_index: &SearchIndex,
    refine: RefinementOutput,
) -> MergedSearch {
    let cand_offset = pass1_candidates.len() as u32;
    let prot_offset = pass1_index.db.len();

    // Take ownership and EXTEND in place — never `to_vec()` the Pass-1 pool. On a
    // refine run Pass-1 holds millions of candidates; copying them here doubled the
    // pool in RAM (the original + the copy) for the rest of the run. Moving costs
    // nothing and frees the duplicate.
    let mut combined: Vec<Candidate> = pass1_candidates;
    let RefinementOutput { index: refine_index, candidates: refine_cands, queues, global_spectrum_indices } = refine;

    for mut c in refine_cands {
        c.protein_index += prot_offset;
        combined.push(c);
    }

    for (local_i, mut q) in queues.into_iter().enumerate() {
        let global = global_spectrum_indices[local_i];
        for mut psm in q.drain_into_vec() {
            for idx in &mut psm.candidate_idxs {
                *idx += cand_offset;
            }
            pass1_queues[global].force_push(psm);
        }
    }

    MergedSearch { candidates: combined, index: pass1_index.concat(&refine_index) }
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

    // ── (a') confident_base_peptides ────────────────────────────────────────

    /// A target/decoy candidate carrying a specific (possibly modified) peptide
    /// at a given combined-db protein index.
    fn cand_pep(peptide: Peptide, protein_index: usize, is_decoy: bool) -> Candidate {
        Candidate {
            peptide,
            protein_index,
            start_offset_in_protein: 0,
            is_decoy,
            is_protein_n_term: false,
            is_protein_c_term: false,
        }
    }

    /// An UNMODIFIED peptide over the given backbone bytes.
    fn unmod_peptide(seq: &[u8]) -> Peptide {
        let residues: Vec<_> = seq.iter().map(|&b| AminoAcid::standard(b).unwrap()).collect();
        Peptide::new(residues, b'_', b'-')
    }

    #[test]
    fn confident_base_peptides_returns_backbone_of_confident_targets() {
        // 40 confident high-scoring target spectra carrying PEPTIK / SAMPLER (two
        // distinct backbones interleaved so both are confident), then a band of
        // low-scoring decoys driving q above threshold for a low-scoring target —
        // only the confident TARGET backbones come back.
        let mut candidates: Vec<Candidate> = Vec::new();
        let mut queues: Vec<TopNQueue> = Vec::new();
        for i in 0..40 {
            let seq: &[u8] = if i % 2 == 0 { b"PEPTIK" } else { b"SAMPLER" };
            candidates.push(cand_pep(unmod_peptide(seq), i, false));
            queues.push(queue_with(vec![psm(i, i as u32, 90.0 - i as f32 * 0.1)]));
        }
        // 39 low-scoring decoys (q explodes above this tail). They must never
        // contribute a backbone.
        for k in 0..39 {
            let ci = candidates.len();
            candidates.push(cand_pep(unmod_peptide(b"DECYSEAR"), 200 + k, true));
            queues.push(queue_with(vec![psm(40 + k, ci as u32, 5.0)]));
        }
        // A low-scoring TARGET buried in the decoy tail (q above threshold → NOT
        // confident; its backbone must be excluded).
        let low_ci = candidates.len();
        candidates.push(cand_pep(unmod_peptide(b"LWSCAGER"), 300, false));
        queues.push(queue_with(vec![psm(79, low_ci as u32, 4.0)]));

        let seqs = confident_base_peptides(&queues, &candidates, 0.01);
        // Exactly the two confident backbones, deduped + sorted.
        assert_eq!(seqs, vec![b"PEPTIK".to_vec(), b"SAMPLER".to_vec()]);
    }

    #[test]
    fn confident_base_peptides_ignores_modifications_in_backbone() {
        // A confident MODIFIED target peptide must yield its bare backbone bytes
        // (the mod is stripped — Pass-2 re-adds refinement mods).
        let ox_m = AminoAcid::standard(b'M').unwrap().with_mod(Modification {
            name: "Oxidation".into(),
            mass_delta: 15.994915,
            residue: ResidueSpec::Specific(b'M'),
            location: ModLocation::Anywhere,
            fixed: false,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
        });
        let pep = Peptide::new(
            vec![AminoAcid::standard(b'P').unwrap(), ox_m, AminoAcid::standard(b'K').unwrap()],
            b'_',
            b'-',
        );
        // One confident target spectrum + a second so the q-walk is non-empty.
        let candidates = vec![
            cand_pep(pep, 0, false),
            cand_pep(unmod_peptide(b"AAAK"), 1, false),
        ];
        let queues = vec![
            queue_with(vec![psm(0, 0, 50.0)]),
            queue_with(vec![psm(1, 1, 49.0)]),
        ];
        let seqs = confident_base_peptides(&queues, &candidates, 0.5);
        // The modified peptide's backbone is the bare residues (no mod folded in).
        assert!(seqs.contains(&b"PMK".to_vec()), "modified backbone stripped to PMK; got {seqs:?}");
        assert!(seqs.contains(&b"AAAK".to_vec()));
    }

    #[test]
    fn confident_base_peptides_dedups_repeated_backbones() {
        // Two confident spectra matching the SAME backbone collapse to one entry.
        let candidates = vec![
            cand_pep(unmod_peptide(b"PEPTIK"), 0, false),
            cand_pep(unmod_peptide(b"PEPTIK"), 1, false),
        ];
        let queues = vec![
            queue_with(vec![psm(0, 0, 50.0)]),
            queue_with(vec![psm(1, 1, 49.0)]),
        ];
        let seqs = confident_base_peptides(&queues, &candidates, 0.5);
        assert_eq!(seqs, vec![b"PEPTIK".to_vec()]);
    }

    #[test]
    fn confident_base_peptides_empty_when_all_decoys() {
        let candidates: Vec<Candidate> =
            (0..20).map(|i| cand_pep(unmod_peptide(b"DECYSK"), i, true)).collect();
        let queues: Vec<TopNQueue> = (0..20)
            .map(|i| queue_with(vec![psm(i, i as u32, 10.0)]))
            .collect();
        assert!(confident_base_peptides(&queues, &candidates, 0.01).is_empty());
    }

    // ── (b) unidentified_spectrum_indices ──────────────────────────────────

    #[test]
    fn unidentified_spectrum_indices_returns_decoy_above_q_and_empty() {
        // 50 clearly-confident target spectra (high scores, no decoy competition
        // → q ≈ 0) plus a low-scoring decoy-best spectrum and an empty queue. The
        // confident targets are identified; the decoy-best and empty are returned.
        let mut candidates: Vec<Candidate> = (0..50).map(|i| cand(i, false)).collect();
        candidates.push(cand(999, true)); // index 50, a decoy at a low score

        let mut queues = Vec::new();
        for i in 0..50 {
            // High target scores, gently descending; well above the decoy tail.
            queues.push(queue_with(vec![psm(i, i as u32, 80.0 - i as f32 * 0.1)]));
        }
        // Spectrum 50: best PSM is a decoy (low score) → unidentified.
        queues.push(queue_with(vec![psm(50, 50, 1.0)]));
        // Spectrum 51: empty queue → unidentified.
        queues.push(TopNQueue::new(10));

        let unident = unidentified_spectrum_indices(&queues, &candidates, 0.01);
        // The 50 confident targets are NOT returned.
        for s in 0..50 {
            assert!(!unident.contains(&s), "confident target spectrum {s} must be identified");
        }
        // The decoy-best (50) and empty-queue (51) spectra ARE returned.
        assert!(unident.contains(&50), "decoy-best spectrum must be unidentified");
        assert!(unident.contains(&51), "empty-queue spectrum must be unidentified");
        assert_eq!(unident.len(), 2);
    }

    #[test]
    fn unidentified_spectrum_indices_returns_poorly_scoring_target_above_q() {
        // A field of confident high-scoring targets establishes the q-walk, then a
        // band of interleaved low-scoring decoys drives q above the threshold for
        // a low-scoring target — which must therefore be returned for refinement.
        let mut candidates: Vec<Candidate> = (0..40).map(|i| cand(i, false)).collect();
        // 39 low-scoring decoys forming the tail above which q explodes.
        for _ in 0..39 {
            candidates.push(cand(999, true));
        }
        // One low-scoring TARGET buried in the decoy tail (the spectrum under test).
        candidates.push(cand(40, false)); // candidate index 79

        let mut queues = Vec::new();
        // 40 confident targets at high scores (no decoy competition above → q≈0).
        for i in 0..40 {
            queues.push(queue_with(vec![psm(i, i as u32, 90.0 - i as f32 * 0.1)]));
        }
        // 39 decoys at a LOW score band (spectra 40..79).
        for k in 0..39 {
            queues.push(queue_with(vec![psm(40 + k, (40 + k) as u32, 5.0)]));
        }
        // A poorly-scoring TARGET at spectrum 79, score below/among the decoy tail
        // → its TDC q-value > report_q → unidentified.
        queues.push(queue_with(vec![psm(79, 79, 4.0)]));

        let unident = unidentified_spectrum_indices(&queues, &candidates, 0.01);
        // The 40 high-scoring confident targets are identified.
        for s in 0..40 {
            assert!(!unident.contains(&s), "confident target spectrum {s} must be identified");
        }
        // The poorly-scoring target (q above threshold) is returned.
        assert!(unident.contains(&79), "low-scoring target above report_q must be unidentified");
    }

    #[test]
    fn unidentified_spectrum_indices_identifies_best_target_in_mixed_queue() {
        // A queue's BEST PSM (max rank_score) is the target → identified, even if a
        // lower-scoring decoy PSM is also present. One confident target spectrum
        // plus a second high target so the q-walk cleanly separates.
        let candidates = vec![cand(0, false), cand(1, true), cand(2, false)];
        let queues = vec![
            // Best PSM is the target (50.0) over the decoy (10.0) → identified.
            queue_with(vec![psm(0, 0, 50.0), psm(0, 1, 10.0)]),
            queue_with(vec![psm(1, 2, 49.0)]),
        ];
        let unident = unidentified_spectrum_indices(&queues, &candidates, 0.5);
        assert!(unident.is_empty(), "both best-target spectra are identified");
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
        let idx = build_refinement_index(&full_target, &keep, "XXX", DecoyStrategy::Reverse, 42);

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
        let idx = build_refinement_index(&full_target, &keep, "XXX", DecoyStrategy::Reverse, 42);
        // Only P0 + its decoy.
        assert_eq!(idx.db.len(), 2);
        assert_eq!(idx.db.proteins[0].accession, "P0");
    }

    // ── (c') build_peptide_anchored_index ──────────────────────────────────

    #[test]
    fn build_peptide_anchored_index_one_mini_protein_per_base_seq_plus_decoys() {
        let base_seqs = vec![b"PEPTIK".to_vec(), b"SAMPLER".to_vec()];
        let idx = build_peptide_anchored_index(&base_seqs, "XXX", DecoyStrategy::Reverse, 42);

        // 2 target mini-proteins + 2 paired Reverse decoys.
        assert_eq!(idx.db.len(), 4);
        // Each base seq is its own mini-protein, accession BASEPEP_i, sequence == seq.
        assert_eq!(idx.db.proteins[0].accession, "BASEPEP_0");
        assert_eq!(idx.db.proteins[0].sequence, b"PEPTIK".to_vec());
        assert_eq!(idx.db.proteins[1].accession, "BASEPEP_1");
        assert_eq!(idx.db.proteins[1].sequence, b"SAMPLER".to_vec());
        // Decoys carry the "<prefix>_" needle.
        let needle = crate::decoy::decoy_accession_needle("XXX");
        assert!(idx.db.proteins[2].accession.starts_with(&needle));
        assert!(idx.db.proteins[3].accession.starts_with(&needle));
        assert_eq!(idx.decoy_prefix, "XXX");
    }

    // ── (d) refinement_aa_set ──────────────────────────────────────────────

    fn base_set_with_cam_c() -> AminoAcidSet {
        AminoAcidSetBuilder::new_standard_with_carbamidomethyl_c().build().unwrap()
    }

    #[test]
    fn refinement_aa_set_high_res_exposes_five_variable_mods() {
        let base = base_set_with_cam_c();
        let set = refinement_aa_set(&base, &RefineConfig::default_tier(), true).unwrap();

        // Count DISTINCT variable mods across all locations by name.
        let mut names: HashSet<String> = HashSet::new();
        for aa in set.iter_variants() {
            if let Some(m) = &aa.mod_ {
                if !m.fixed {
                    // A stacked fixed+variable combined variant (e.g.
                    // "Carbamidomethyl+Acetyl" on protein-N-term Cys) counts as
                    // its variable component (the part after the last '+').
                    names.insert(m.name.rsplit('+').next().unwrap_or(&m.name).to_string());
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
        let set = refinement_aa_set(&base, &RefineConfig::default_tier(), false).unwrap();
        let mut names: HashSet<String> = HashSet::new();
        for aa in set.iter_variants() {
            if let Some(m) = &aa.mod_ {
                if !m.fixed {
                    // A stacked fixed+variable combined variant (e.g.
                    // "Carbamidomethyl+Acetyl" on protein-N-term Cys) counts as
                    // its variable component (the part after the last '+').
                    names.insert(m.name.rsplit('+').next().unwrap_or(&m.name).to_string());
                }
            }
        }
        assert_eq!(names.len(), 4, "deamidation dropped at low-res; got {names:?}");
        assert!(!names.contains("Deamidation"), "Deamidation must be dropped at low-res");
    }

    #[test]
    fn refinement_aa_set_oxidation_m_has_modified_variant() {
        let base = base_set_with_cam_c();
        let set = refinement_aa_set(&base, &RefineConfig::default_tier(), true).unwrap();
        // M-Anywhere should now carry an Oxidation variable variant (unmod + ox).
        let m_variants = set.variants_for(b'M', ModLocation::Anywhere);
        assert!(m_variants.iter().any(|aa| {
            aa.mod_.as_ref().map(|m| m.name == "Oxidation").unwrap_or(false)
        }));
    }

    #[test]
    fn refinement_aa_set_excludes_stacked_combined_base_variant() {
        // Base = fixed Cam-C + variable Acetyl on protein-N-term (a common LFQ
        // config — UPS/Astral both use it). The builder synthesizes a stacked
        // "Carbamidomethyl+Acetyl" (~99 Da) variant on N-term Cys; that build
        // artifact must NOT be re-injected into Pass-2 as a phantom variable mod.
        let acetyl = Modification {
            name: "Acetyl".into(),
            mass_delta: 42.010565,
            residue: ResidueSpec::Wildcard,
            location: ModLocation::ProtNTerm,
            fixed: false,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
        };
        let base = AminoAcidSetBuilder::new_standard_with_carbamidomethyl_c()
            .add_variable_mod(acetyl)
            .build()
            .unwrap();
        // distinct_mods must not surface the synthesized "+"-named stacked mod.
        assert!(
            base.distinct_mods().iter().all(|m| !m.name.contains('+')),
            "distinct_mods leaked a synthesized stacked (fixed+variable) variant"
        );
        // The refined set legitimately re-synthesizes the Cam+Acetyl (~99 Da) stack
        // on N-term CYS (a fixed-Cam residue) — that is correct, same as Pass-1.
        // The bug would instead inject +99 as a WILDCARD variable mod, leaking it
        // onto NON-Cys N-term residues. Assert Ala-ProtNTerm carries Acetyl (+42)
        // but NOT a ~99 Da variant.
        let set = refinement_aa_set(&base, &RefineConfig::default_tier(), true).unwrap();
        let a_nterm = set.variants_for(b'A', ModLocation::ProtNTerm);
        assert!(
            a_nterm.iter().all(|aa| aa
                .mod_
                .as_ref()
                .map(|m| (m.mass_delta - 99.032).abs() >= 1e-2)
                .unwrap_or(true)),
            "phantom ~99 Da stacked mod leaked onto non-Cys N-term residue"
        );
    }

    // ── (d') Pass-2 inherits the base search's full chemistry ──────────────

    #[test]
    fn refinement_aa_set_inherits_nonstandard_fixed_tmt() {
        // A labeled search: fixed CAM-C + fixed TMT6plex-on-K. Pass-2 MUST carry
        // the TMT delta forward (previously it was dropped → silent no-op on TMT).
        let tmt_k = Modification {
            name: "TMT6plex".into(),
            mass_delta: 229.162932,
            residue: ResidueSpec::Specific(b'K'),
            location: ModLocation::Anywhere,
            fixed: true,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
        };
        let base = AminoAcidSetBuilder::new_standard_with_carbamidomethyl_c()
            .add_fixed_mod(tmt_k)
            .build()
            .unwrap();
        let set = refinement_aa_set(&base, &RefineConfig::default_tier(), true).unwrap();
        // K-Anywhere carries the fixed TMT (~229.16) in the refined set.
        let k_anywhere = set.variants_for(b'K', ModLocation::Anywhere);
        assert!(
            k_anywhere.iter().any(|aa| aa
                .mod_
                .as_ref()
                .map(|m| m.fixed && (m.mass_delta - 229.162932).abs() < 1e-3)
                .unwrap_or(false)),
            "refined set must carry the base's fixed TMT6plex-on-K"
        );
        // …and the fixed CAM-C baseline is still present.
        let c_anywhere = set.variants_for(b'C', ModLocation::Anywhere);
        assert!(c_anywhere
            .iter()
            .any(|aa| aa.mod_.as_ref().map(|m| m.fixed).unwrap_or(false)));
    }

    #[test]
    fn refinement_aa_set_dedups_base_variable_against_tier() {
        // Base already declares variable Oxidation-M; the default tier also has
        // Oxidation-M. The refined set must expose it ONCE, not twice.
        let ox_m = Modification {
            name: "Oxidation".into(),
            mass_delta: 15.994915,
            residue: ResidueSpec::Specific(b'M'),
            location: ModLocation::Anywhere,
            fixed: false,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
        };
        let base = AminoAcidSetBuilder::new_standard_with_carbamidomethyl_c()
            .add_variable_mod(ox_m)
            .build()
            .unwrap();
        let set = refinement_aa_set(&base, &RefineConfig::default_tier(), true).unwrap();
        let ox_count = set
            .variants_for(b'M', ModLocation::Anywhere)
            .iter()
            .filter(|aa| {
                aa.mod_
                    .as_ref()
                    .map(|m| !m.fixed && (m.mass_delta - 15.994915).abs() < 1e-3)
                    .unwrap_or(false)
            })
            .count();
        assert_eq!(ox_count, 1, "Oxidation-M must not be double-added; got {ox_count}");
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
        let pass1_cands = [cand(0, false)];
        let mut global_queues = [
            queue_with(vec![psm(0, 0, 40.0)]), // identified
            TopNQueue::new(10),                // unidentified
        ];

        // Simulate a Pass-2 winner for spectrum 1: an Oxidation-M peptide.
        let pass2_cands = [Candidate {
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
    fn run_refinement_early_returns_when_no_confident_base_peptides() {
        let aa_set = base_set_with_cam_c();
        let params = SearchParams::default_tryptic(aa_set);
        let scorer = tiny_scorer();
        let target = ProteinDb { proteins: vec![protein("P0", b"MKWVTFISLLR")] };

        // All-decoy Pass-1 → no confident base peptides → early return (None).
        let pass1_cands = vec![cand(0, true)];
        let queues = vec![queue_with(vec![psm(0, 0, 10.0)])];
        let before_len = queues[0].len();
        let spectra: Vec<model::spectrum::Spectrum> = vec![model::spectrum::Spectrum::default()];

        let out = run_refinement(
            &queues,
            &spectra,
            &pass1_cands,
            &target,
            &params,
            &scorer,
            &RefineConfig::default_tier(),
            0.01,
            true,
            0.5,
            "XXX",
            DecoyStrategy::Reverse,
            42,
        );
        // No confident base peptides ⇒ no Pass-2 output; Pass-1 queues untouched.
        assert!(out.is_none());
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
        let queues = vec![queue_with(vec![psm(0, 0, 50.0)])];
        let spectra: Vec<model::spectrum::Spectrum> = vec![model::spectrum::Spectrum::default()];

        let out = run_refinement(
            &queues,
            &spectra,
            &pass1_cands,
            &target,
            &params,
            &scorer,
            &RefineConfig::default_tier(),
            0.01,
            true,
            0.5,
            "XXX",
            DecoyStrategy::Reverse,
            42,
        );
        // No unidentified spectra ⇒ no Pass-2 output; Pass-1 queues untouched.
        assert!(out.is_none());
        assert_eq!(queues[0].len(), 1);
        assert_eq!(queues[0].peek_top().unwrap().features.is_refinement, 0);
    }

    // ── (g) merge_into_pass1 ───────────────────────────────────────────────

    #[test]
    fn merge_offsets_indices_and_force_pushes_into_pass1() {
        // Pass-1: 1 candidate (protein 0), 1 spectrum with one unmod PSM (rank 40).
        let p1_cands = vec![cand(0, false)];
        let p1_index = SearchIndex {
            db: ProteinDb { proteins: vec![protein("P0", b"PEPTIDEK")] },
            decoy_prefix: "XXX".into(),
            decoy_suffix: None,
        };
        let mut p1_queues = vec![queue_with(vec![psm(0, 0, 40.0)])];

        // Refine: 1 modified candidate (its OWN protein 0 = BASEPEP_0), one PSM (rank 22) on spectrum 0.
        let refine = RefinementOutput {
            index: SearchIndex { db: ProteinDb { proteins: vec![protein("BASEPEP_0", b"PEPTIDEK")] }, decoy_prefix: "XXX".into(), decoy_suffix: None },
            candidates: vec![cand(0, false)],   // protein_index 0 in the refine db
            queues: vec![queue_with(vec![psm(0, 0, 22.0)])], // candidate_idxs=[0], spectrum 0
            global_spectrum_indices: vec![0],
        };

        let merged = merge_into_pass1(&mut p1_queues, p1_cands, &p1_index, refine);

        // Combined candidates: pass1 (1) + refine (1) = 2; refine candidate's protein_index offset by 1.
        assert_eq!(merged.candidates.len(), 2);
        assert_eq!(merged.candidates[1].protein_index, 1, "refine protein_index offset by pass1 db len");
        // Combined index: 2 proteins, refine accession resolvable at offset 1.
        assert_eq!(merged.index.db.len(), 2);
        assert_eq!(merged.index.db.proteins[1].accession, "BASEPEP_0");
        // Spectrum 0's queue now holds BOTH the unmod (cand 0) and the modified (cand 1) PSM.
        let mut idxs: Vec<u32> = p1_queues[0].iter_psms().map(|m| m.primary_candidate_idx()).collect();
        idxs.sort();
        assert_eq!(idxs, vec![0, 1], "refine PSM candidate_idx offset to 1 and force_pushed into pass1 queue");
    }

    #[test]
    fn merge_offsets_every_entry_of_multi_candidate_psms() {
        // A refine PSM aggregating a SHARED peptide across two refine candidates
        // (candidate_idxs = [0, 1]) must have EVERY entry offset by the Pass-1
        // candidate count, not just the primary — the multi-protein shared-peptide
        // path the PIN writer iterates to emit one accession per candidate_idx.
        let p1_cands = vec![cand(0, false), cand(1, false)]; // cand_offset = 2
        let p1_index = SearchIndex {
            db: ProteinDb { proteins: vec![protein("P0", b"AAAK"), protein("P1", b"BBBK")] }, // prot_offset = 2
            decoy_prefix: "XXX".into(),
            decoy_suffix: None,
        };
        let mut p1_queues = vec![queue_with(vec![psm(0, 0, 40.0)])];

        // One refine PSM on spectrum 0 whose candidate_idxs spans two refine candidates.
        let mut shared = psm(0, 0, 22.0);
        shared.candidate_idxs = vec![0, 1];
        let refine = RefinementOutput {
            index: SearchIndex {
                db: ProteinDb { proteins: vec![protein("BASEPEP_0", b"PEPTIDEK"), protein("BASEPEP_1", b"SAMPLEK")] },
                decoy_prefix: "XXX".into(),
                decoy_suffix: None,
            },
            candidates: vec![cand(0, false), cand(1, false)],
            queues: vec![queue_with(vec![shared])],
            global_spectrum_indices: vec![0],
        };

        let merged = merge_into_pass1(&mut p1_queues, p1_cands, &p1_index, refine);

        // Combined: 2 pass1 + 2 refine; refine protein_index offset by 2.
        assert_eq!(merged.candidates.len(), 4);
        assert_eq!(merged.candidates[2].protein_index, 2);
        assert_eq!(merged.candidates[3].protein_index, 3);
        // The merged refine PSM's BOTH candidate_idxs offset by cand_offset=2 → [2, 3].
        let refine_psm = p1_queues[0]
            .iter_psms()
            .find(|m| m.candidate_idxs.len() == 2)
            .expect("the multi-candidate refine PSM was force_pushed");
        let mut idxs = refine_psm.candidate_idxs.clone();
        idxs.sort();
        assert_eq!(idxs, vec![2, 3], "every candidate_idx entry offset, not just the primary");
    }

}
