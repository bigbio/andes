// Hybrid backbone candidate generation: DB-branch ∪ de-novo Y-ladder.
//
// DB branch: for each known glycan composition, `backbone = precursor_neutral − glycan.mass`.
// De-novo branch: the existing Y-ladder complement-pair solver (backbone.rs).
// Union: merge candidates within 0.02 Da, preferring DB source over de-novo.
//
// Mass convention: `BackboneHit::backbone_mass` is always the peptide RESIDUE
// mass (water NOT included), matching the DB branch's
// `precursor_neutral − glycan.mass` (where `precursor_neutral` is itself
// computed with H2O already subtracted, see glyco_search.rs). `solve_backbone`
// derives its candidate mass from the Y0 ion (`Y0 = bare_peptide_NEUTRAL_mass +
// PROTON`), i.e. the peptide NEUTRAL mass (water included). We subtract H2O
// here so both branches agree on one convention before union/dedup/filter.

use crate::backbone::{solve_backbone_min, H2O};
use crate::glycan_db::GlycanComp;
use crate::oxonium::oxonium_gate;

/// Source of a backbone hit.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Source {
    /// Backbone computed as precursor_neutral − known glycan mass.
    Db,
    /// Backbone proposed by the de-novo Y-ladder solver.
    DeNovo,
}

/// A single backbone candidate from either the DB branch or the de-novo solver.
#[derive(Debug, Clone)]
pub struct BackboneHit {
    pub backbone_mass: f64,
    /// The glycan composition that produced this backbone (Db source only).
    pub glycan: Option<GlycanComp>,
    pub source: Source,
    /// Precursor charge state used to derive `precursor_neutral` (and thus
    /// this `backbone_mass`). Scoring/emission MUST use this charge, not an
    /// independently re-picked one, or the reported mass/backbone pairing is
    /// inconsistent with what was actually matched.
    pub charge: u8,
    /// Isotope offset (in units of `model::mass::ISOTOPE`) that was subtracted
    /// from the observed precursor mass before deriving `backbone_mass`. `0` =
    /// monoisotopic. Mirrors the standard search path's `isotope_error_range`
    /// handling (see `search_params.rs`), which glyco previously ignored.
    pub isotope_offset: i8,
    /// Observed residual glycan mass = `precursor_neutral − backbone_mass`, the
    /// intact glycan actually implied by this backbone. Populated for EVERY hit,
    /// including novel/unannotated ones where `glycan` is `None` (a de-novo
    /// glycan that matched no known composition). Downstream code must use this
    /// for the intact `CalcMass` rather than `glycan.map(|g| g.mass).unwrap_or(0.0)`,
    /// which silently reported the bare peptide for novel glycans.
    pub glycan_mass_residual: f64,
}

/// DB-branch backbone enumeration.
///
/// For each glycan in `glycans`, compute `bb = precursor_neutral − glycan.mass`.
/// Keep if `bb ≥ min_backbone`. Returns candidates sorted by backbone_mass ascending.
pub fn db_branch(
    precursor_neutral: f64,
    glycans: &[GlycanComp],
    min_backbone: f64,
    charge: u8,
    isotope_offset: i8,
) -> Vec<BackboneHit> {
    let mut out: Vec<BackboneHit> = glycans
        .iter()
        .filter_map(|g| {
            let bb = precursor_neutral - g.mass;
            if bb >= min_backbone {
                Some(BackboneHit {
                    backbone_mass: bb,
                    glycan: Some(g.clone()),
                    source: Source::Db,
                    charge,
                    isotope_offset,
                    // By construction bb = precursor − g.mass, so the observed
                    // residual is exactly the glycan mass.
                    glycan_mass_residual: precursor_neutral - bb,
                })
            } else {
                None
            }
        })
        .collect();

    // Total-order sort: backbone_mass asc, tiebreak by glycan fields for determinism.
    out.sort_by(|a, b| {
        a.backbone_mass
            .to_bits()
            .cmp(&b.backbone_mass.to_bits())
            .then_with(|| {
                // tiebreak by glycan composition (always Some for Db branch)
                let ga = a.glycan.as_ref();
                let gb = b.glycan.as_ref();
                match (ga, gb) {
                    (Some(x), Some(y)) => x
                        .hexnac
                        .cmp(&y.hexnac)
                        .then(x.hex.cmp(&y.hex))
                        .then(x.fuc.cmp(&y.fuc))
                        .then(x.neuac.cmp(&y.neuac))
                        .then(x.neugc.cmp(&y.neugc)),
                    _ => std::cmp::Ordering::Equal,
                }
            })
    });

    out
}

/// Find the known glycan composition nearest to `target_mass` within tolerance.
///
/// Used to ANNOTATE a Y-ladder-derived backbone with its glycan (glycan =
/// precursor_neutral − backbone). Returns the closest composition within
/// `max(target_mass·tol_ppm·1e-6, 0.02)` Da, or `None` when the residual matches
/// no known glycan (a novel/unexpected glycan — kept, not discarded). Linear
/// scan (glyco lists are ~hundreds of entries, a handful of backbones/spectrum).
fn nearest_glycan(glycans: &[GlycanComp], target_mass: f64, tol_ppm: f64) -> Option<GlycanComp> {
    if target_mass <= 0.0 {
        return None;
    }
    let tol = (target_mass * tol_ppm * 1e-6_f64).max(0.02);
    let mut best: Option<(f64, &GlycanComp)> = None;
    for g in glycans {
        let d = (g.mass - target_mass).abs();
        if d <= tol && best.map_or(true, |(bd, _)| d < bd) {
            best = Some((d, g));
        }
    }
    best.map(|(_, g)| g.clone())
}

/// Hybrid backbone candidate union: DB-branch ∪ de-novo Y-ladder.
///
/// 1. Run oxonium gate (min_frac=0.10, tol=20 ppm). If it doesn't fire,
///    still run DB branch (DB doesn't need oxonium evidence) but skip de-novo.
/// 2. Collect DB-branch candidates.
/// 3. Run de-novo solver (top_k = top_k); append as Source::DeNovo.
/// 4. Dedup: merge any two candidates within 0.02 Da, preferring Db source.
///    When two Db candidates cluster, keep the first (by backbone_mass order).
///
/// Returns all candidates sorted by backbone_mass ascending (deterministic).
pub fn hybrid_candidates(
    peaks: &[(f64, f32)],
    precursor_neutral: f64,
    precursor_z: u8,
    glycans: &[GlycanComp],
    tol_ppm: f64,
    top_k: usize,
) -> Vec<BackboneHit> {
    hybrid_candidates_with_isotope(peaks, precursor_neutral, precursor_z, 0, glycans, tol_ppm, top_k)
}

/// Same as [`hybrid_candidates`] but records the isotope offset that produced
/// `precursor_neutral` on every returned [`BackboneHit`] (see
/// [`BackboneHit::isotope_offset`]). Callers that try multiple isotope offsets
/// (mirroring the standard search's `isotope_error_range`) should call this
/// once per offset and union the results.
pub fn hybrid_candidates_with_isotope(
    peaks: &[(f64, f32)],
    precursor_neutral: f64,
    precursor_z: u8,
    isotope_offset: i8,
    glycans: &[GlycanComp],
    tol_ppm: f64,
    top_k: usize,
) -> Vec<BackboneHit> {
    const MIN_BACKBONE: f64 = 500.0;

    // === Y-ION-FIRST candidate generation ===
    // The de-novo Y-ladder solver is the PRIMARY (and only) backbone source: we
    // read the glycan-independent trimannosyl-core Y-ion ladder to determine the
    // backbone mass from the spectrum, then annotate the glycan by subtraction
    // (glycan = precursor_neutral − backbone). This is the architecture of
    // a glyco search engine / a cross-spectrum glyco engine / StrucGP. The previous design enumerated
    // `precursor − every glycan` (~600 brute-force backbones) ranked by a
    // glycan-blind peptide b/y score, which forced top-k ≈ 600 and ~19%
    // find-rate. `db_branch` is retained (tests / annotation) but is no longer
    // used to GENERATE backbones.

    // Oxonium gate: glyco-spectrum filter (every mature engine gates on oxonium).
    let ox = oxonium_gate(peaks, 0.10, tol_ppm);
    if !ox.fired {
        return Vec::new();
    }

    // Y-ION-FIRST confidence cascade (replaces both the deleted DB brute force
    // AND the Codex-flagged silent drop). The core-Y ladder is the primary
    // evidence, but we degrade the quorum gracefully instead of giving up:
    //   1. PRIMARY (quorum 2): ≥2 distinct core-Y rungs → high-confidence,
    //      few candidates. Handles the Y-ladder-rich majority.
    //   2. RESCUE (quorum 1): if the primary yields nothing, relax to a single
    //      core-Y rung (Y0/Y1 — the most diagnostic). This recovers the
    //      weak-ladder tail (Codex finding #1) while staying evidence-driven
    //      and bounded to `top_k` — unlike a blind `precursor − every glycan`
    //      enumeration (~600 backbones/spectrum), which is both O(glycans) slow
    //      and precision-poor. Spectra with NO core-Y ladder at all remain
    //      unsearchable by this method (a documented limitation; they need an
    //      oxonium-composition / cross-spectrum lever, not brute force).
    let mut dn = solve_backbone_min(peaks, precursor_neutral, precursor_z, tol_ppm, top_k, 2);
    if dn.is_empty() {
        dn = solve_backbone_min(peaks, precursor_neutral, precursor_z, tol_ppm, top_k, 1);
    }
    let mut combined: Vec<BackboneHit> = Vec::with_capacity(dn.len());
    for c in dn {
        // `solve_backbone` returns the Y0-derived peptide NEUTRAL mass (water
        // included); convert to the RESIDUE convention used everywhere else.
        let bb = c.backbone_mass - H2O;
        if bb < MIN_BACKBONE {
            continue;
        }
        // Annotate the glycan by subtraction. `precursor_neutral` and `bb` are
        // both residue-convention, so `precursor_neutral − bb` = glycan mass.
        let residual = precursor_neutral - bb;
        let glycan = nearest_glycan(glycans, residual, tol_ppm);
        let source = if glycan.is_some() {
            Source::Db
        } else {
            Source::DeNovo
        };
        combined.push(BackboneHit {
            backbone_mass: bb,
            glycan,
            source,
            charge: precursor_z,
            isotope_offset,
            // Keep the observed residual even when annotation returned None so a
            // novel glycan's intact mass is not lost (Codex finding #3).
            glycan_mass_residual: residual,
        });
    }

    // BOUNDED DB FALLBACK (Codex re-review #1): if even the relaxed quorum-1
    // solver found nothing — a truly 0-core-Y oxonium-positive spectrum (no
    // detectable Y-ladder at all) — fall back to `precursor − known glycan` so
    // these real glyco spectra are not silently dropped. This fires ONLY for
    // the rare 0-rung set (NOT the weak-ladder tail, which the quorum-1 pass
    // already handles), so it avoids the O(glycans) cost blowing up broadly;
    // the downstream b/y ranking + `backbone_top_k` in `glyco_search` bound the
    // per-spectrum candidate count.
    if combined.is_empty() {
        combined = db_branch(precursor_neutral, glycans, MIN_BACKBONE, precursor_z, isotope_offset);
    }

    // --- Sort all candidates by backbone_mass for dedup pass ---
    // Total-order: backbone_mass bits asc, then Source::Db before DeNovo.
    combined.sort_by(|a, b| {
        a.backbone_mass
            .to_bits()
            .cmp(&b.backbone_mass.to_bits())
            .then_with(|| {
                // Db < DeNovo so Db wins in the dedup step below
                let ord_a = if a.source == Source::Db { 0u8 } else { 1u8 };
                let ord_b = if b.source == Source::Db { 0u8 } else { 1u8 };
                ord_a.cmp(&ord_b)
            })
    });

    // --- Dedup: merge candidates within max(bb*tol_ppm*1e-6, 0.01) Da, prefer Db source ---
    // Window derives from the caller's `tol_ppm` (previously hardcoded to a
    // 20 ppm assumption) so dedup always agrees with the gate/searchable
    // window actually configured for this run.
    // Strategy: single-pass, keep a running "cluster representative".
    // Because Db is sorted before DeNovo at equal mass, when a cluster contains
    // both, the first element is always Db → the representative is Db.
    if combined.is_empty() {
        return combined;
    }

    let mut deduped: Vec<BackboneHit> = Vec::with_capacity(combined.len());
    let mut rep = combined.remove(0);

    for next in combined {
        let tol = (rep.backbone_mass * tol_ppm * 1e-6_f64).max(0.01);
        if (next.backbone_mass - rep.backbone_mass).abs() < tol {
            // Same cluster: prefer Db source. Since Db is sorted first, rep is
            // already Db if any Db candidate exists in this cluster.
            // If rep is DeNovo and next is Db (shouldn't happen due to sort order,
            // but guard defensively):
            if rep.source == Source::DeNovo && next.source == Source::Db {
                rep = next;
            }
            // else keep rep as-is
        } else {
            deduped.push(rep);
            rep = next;
        }
    }
    deduped.push(rep);

    deduped
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::glycan_db::n_glycan_list;
    use crate::glycan_mass::{HEX, HEXNAC, PROTON};

    /// DB branch must return the correct backbone when precursor = backbone + known glycan.
    #[test]
    fn db_branch_recovers_backbone_for_known_glycan() {
        let glycans = n_glycan_list();
        // Use HexNAc2Hex3 (trimannosyl core) mass as the glycan.
        let glycan_mass = 2.0 * HEXNAC + 3.0 * HEX; // ~892.317 Da
        let true_backbone = 1500.0_f64;
        let precursor = true_backbone + glycan_mass;

        let hits = db_branch(precursor, &glycans, 500.0, 2, 0);
        assert!(!hits.is_empty(), "expected DB branch hits");

        // Must include a hit within ±0.01 Da of true_backbone.
        let found = hits
            .iter()
            .any(|h| (h.backbone_mass - true_backbone).abs() < 0.01 && h.source == Source::Db);
        assert!(found, "did not find backbone at {:.4} in DB branch", true_backbone);
    }

    /// DB branch must filter out backbones below min_backbone.
    #[test]
    fn db_branch_filters_below_min() {
        let glycans = n_glycan_list();
        // Very small precursor so backbone would be < 500 Da.
        let precursor = 600.0; // glycan of ~100 Da not in list; backbone ~100 Da
        let hits = db_branch(precursor, &glycans, 500.0, 2, 0);
        for h in &hits {
            assert!(h.backbone_mass >= 500.0, "backbone below min: {}", h.backbone_mass);
        }
    }

    /// DB branch is sorted by backbone_mass ascending.
    #[test]
    fn db_branch_is_sorted() {
        let glycans = n_glycan_list();
        let precursor = 4000.0;
        let hits = db_branch(precursor, &glycans, 500.0, 2, 0);
        for w in hits.windows(2) {
            assert!(
                w[0].backbone_mass <= w[1].backbone_mass + 1e-9,
                "not sorted: {} > {}",
                w[0].backbone_mass,
                w[1].backbone_mass
            );
        }
    }

    /// hybrid_candidates (now Y-ion-first: solve_backbone is the sole
    /// generator) must recover the true backbone AND annotate it against the
    /// glycan list when the implied glycan mass matches a known composition.
    ///
    /// NOTE: `solve_backbone` anchors its Y0 rung at the peptide NEUTRAL mass
    /// (`neutral_backbone + PROTON`, water INCLUDED) — peaks must be built at
    /// `true_backbone_residue + H2O + PROTON`, not `true_backbone_residue +
    /// PROTON`, or the recovered residue mass (and thus the by-subtraction
    /// glycan mass) is off by H2O and no longer annotates (see the 4d1362be
    /// de-novo/DB H2O-convention fix).
    #[test]
    fn hybrid_union_contains_both_sources() {
        let glycans = n_glycan_list();
        let proton = PROTON;
        let steps = crate::glycan_mass::CORE_Y_STEPS;

        let glycan_mass = 2.0 * HEXNAC + 3.0 * HEX; // HexNAc2Hex3 ~892.317
        let true_backbone_residue = 1500.0_f64;
        let precursor = true_backbone_residue + glycan_mass;
        let y0_neutral = true_backbone_residue + H2O;

        // Build a synthetic spectrum: oxonium peaks + full core-Y ladder.
        let mut peaks: Vec<(f64, f32)> = vec![
            (204.08665, 200.0), // HexNAc oxonium
            (138.05496, 150.0), // HexNAc fragment
            (186.07608, 80.0),  // HexNAc ring-open
            (y0_neutral + proton, 100.0), // Y0
        ];
        for &s in steps.iter() {
            peaks.push((y0_neutral + s + proton, 90.0));
        }

        let hits = hybrid_candidates(&peaks, precursor, 2, &glycans, 20.0, 10);
        assert!(!hits.is_empty(), "expected hybrid hits");

        // The recovered backbone must annotate to the known HexNAc2Hex3
        // composition (Source::Db), since the implied glycan mass matches.
        let has_db = hits.iter().any(|h| h.source == Source::Db);
        assert!(has_db, "expected at least one annotated (Db) hit, got {:?}", hits);
    }

    /// Dedup: near-duplicate solver candidates within 0.02 Da collapse to one
    /// representative, and it must be annotated (Source::Db) when the implied
    /// glycan mass matches a known composition.
    #[test]
    fn hybrid_dedup_keeps_db_over_denovo() {
        let glycans = n_glycan_list();
        let proton = PROTON;
        let steps = crate::glycan_mass::CORE_Y_STEPS;

        // True backbone matches a DB glycan exactly.
        let glycan_mass = 2.0 * HEXNAC + 3.0 * HEX; // ~892.317
        let true_backbone_residue = 1500.0_f64;
        let precursor = true_backbone_residue + glycan_mass;
        let y0_neutral = true_backbone_residue + H2O;

        // Full core-Y ladder anchored at the correct NEUTRAL mass.
        let mut peaks: Vec<(f64, f32)> = vec![
            (204.08665, 200.0),
            (138.05496, 150.0),
            (186.07608, 80.0),
            (y0_neutral + proton, 100.0),
        ];
        for &s in steps.iter() {
            peaks.push((y0_neutral + s + proton, 90.0));
        }

        let hits = hybrid_candidates(&peaks, precursor, 2, &glycans, 20.0, 10);

        // Find candidates near true_backbone; after dedup there must be exactly one
        // in the ±0.02 Da window, and it must be Source::Db.
        let near: Vec<&BackboneHit> = hits
            .iter()
            .filter(|h| (h.backbone_mass - true_backbone_residue).abs() < 0.02)
            .collect();
        assert_eq!(near.len(), 1, "expected exactly one candidate near true backbone after dedup, got {}", near.len());
        assert_eq!(
            near[0].source,
            Source::Db,
            "expected Db source to win dedup, got {:?}",
            near[0].source
        );
    }

    /// BUG 2 regression: the DB branch and the de-novo branch must agree on
    /// the RESIDUE-mass convention (water NOT included) for `backbone_mass`.
    ///
    /// Before the fix, `solve_backbone` returned the peptide NEUTRAL mass
    /// (derived from the Y0 ion `= bare_peptide_neutral + PROTON`), while the
    /// DB branch computed `precursor_neutral - glycan.mass` where
    /// `precursor_neutral` already has H2O subtracted (residue-mass
    /// convention). The two branches were therefore ~18 Da (H2O) apart on a
    /// backbone that should be identical, silently breaking dedup/union and
    /// the driver's exact-mass candidate filter for de-novo hits.
    ///
    /// This test isolates `solve_backbone`'s raw output against the known
    /// residue mass and asserts the H2O gap directly, then verifies
    /// `hybrid_candidates` (which applies the fix) reports the de-novo
    /// backbone within tight tolerance of the residue-mass convention.
    #[test]
    fn denovo_and_db_branches_agree_on_residue_mass_convention() {
        let proton = PROTON;
        let steps = crate::glycan_mass::CORE_Y_STEPS;

        // Residue-mass convention backbone (matches the DB branch/driver).
        let true_backbone_residue = 1500.0_f64;
        let glycan_mass = 2.0 * HEXNAC + 3.0 * HEX; // HexNAc2Hex3, no DB entry needed here
        let precursor_neutral = true_backbone_residue + glycan_mass;

        // Build a synthetic core-Y ladder anchored at the RESIDUE-mass
        // backbone (Y0 = bare_peptide_NEUTRAL + PROTON = (residue+H2O)+PROTON),
        // plus oxonium peaks so `hybrid_candidates`'s de-novo branch fires
        // (it is gated behind `oxonium_gate`).
        let neutral_backbone = true_backbone_residue + crate::backbone::H2O;
        let mut peaks: Vec<(f64, f32)> = vec![
            (204.08665, 200.0), // HexNAc oxonium
            (138.05496, 150.0), // HexNAc fragment
            (186.07608, 80.0),  // HexNAc ring-open
            (neutral_backbone + proton, 100.0),
        ];
        for &s in steps.iter() {
            peaks.push((neutral_backbone + s + proton, 90.0));
        }

        // 1. Raw solve_backbone output is in NEUTRAL-mass convention: it must
        //    be ~H2O (18.0106 Da) ABOVE the residue-mass backbone.
        let raw = crate::backbone::solve_backbone(&peaks, precursor_neutral, 2, 20.0, 5);
        assert!(!raw.is_empty(), "expected solve_backbone candidates");
        let raw_gap = raw[0].backbone_mass - true_backbone_residue;
        assert!(
            (raw_gap - crate::backbone::H2O).abs() < 0.01,
            "expected solve_backbone's raw output to be ~H2O ({:.4}) above the \
             residue-mass backbone, got gap={:.4} (raw={:.4}, residue={:.4})",
            crate::backbone::H2O,
            raw_gap,
            raw[0].backbone_mass,
            true_backbone_residue
        );

        // 2. After the hybrid_candidates fix, the CLOSEST DeNovo BackboneHit
        //    to the true residue-mass backbone must land within tight
        //    tolerance (gap ~0), not ~H2O away (the pre-fix neutral-mass
        //    convention). `solve_backbone` can return multiple candidate
        //    clusters (top_k=5); we pick the nearest one to the known true
        //    value rather than assuming index 0, since cluster ORDER is not
        //    what this test is verifying — the MASS CONVENTION is.
        let hits = hybrid_candidates(&peaks, precursor_neutral, 2, &[], 20.0, 5);
        let dn_hit = hits
            .iter()
            .filter(|h| h.source == Source::DeNovo)
            .min_by(|a, b| {
                (a.backbone_mass - true_backbone_residue)
                    .abs()
                    .total_cmp(&(b.backbone_mass - true_backbone_residue).abs())
            })
            .expect("expected at least one DeNovo BackboneHit");
        assert!(
            (dn_hit.backbone_mass - true_backbone_residue).abs() < 0.01,
            "DeNovo backbone_mass must match residue-mass convention: expected \
             ~{:.4}, got {:.4} (gap={:.4}); note the pre-fix bug would leave this \
             ~H2O ({:.4}) away instead",
            true_backbone_residue,
            dn_hit.backbone_mass,
            dn_hit.backbone_mass - true_backbone_residue,
            crate::backbone::H2O
        );
    }

    /// P0.2 (Codex re-review #1): a truly 0-core-Y oxonium-positive spectrum
    /// (NO core-Y ladder at all → both the quorum-2 and quorum-1 solver passes
    /// return nothing) carrying a KNOWN glycan must still be searchable via the
    /// bounded DB fallback (`precursor − known glycan`). The Y-first-only path
    /// silently dropped these; they are a real, representable slice of the data
    /// (DB-branch ceiling 85.9% vs cascade 80.1%).
    #[test]
    fn zero_core_y_known_glycan_is_recovered_by_db_fallback() {
        let glycans = n_glycan_list();
        let glycan_mass = 2.0 * HEXNAC + 3.0 * HEX; // HexNAc2Hex3, in list
        let true_backbone_residue = 1500.0_f64;
        let precursor = true_backbone_residue + glycan_mass;

        // Oxonium only — NO core-Y ladder anywhere.
        let peaks: Vec<(f64, f32)> = vec![
            (204.08665, 200.0),
            (138.05496, 150.0),
            (366.13950, 60.0),
            (900.0, 20.0),
        ];
        // Premise: both solver quorums find nothing.
        assert!(crate::backbone::solve_backbone_min(&peaks, precursor, 2, 20.0, 50, 1).is_empty(),
            "premise: even quorum-1 solver must be empty for a 0-core-Y spectrum");

        let hits = hybrid_candidates(&peaks, precursor, 2, &glycans, 20.0, 50);
        let found = hits.iter().any(|h| {
            (h.backbone_mass - true_backbone_residue).abs() < 0.02 && h.source == Source::Db
        });
        assert!(found, "DB fallback must recover the 0-core-Y known-glycan backbone; got {hits:?}");
    }

    /// Finding #3 regression (Codex adversarial review): a Y-first backbone
    /// whose implied glycan matches NO known composition (a novel/unexpected
    /// glycan — explicitly "kept, not discarded") must still carry the observed
    /// residual glycan mass (`precursor_neutral − backbone`). Previously the hit
    /// was pushed with `glycan: None` and the residual was thrown away, so
    /// downstream `glycan_mass` fell back to 0.0 → the intact PIN `CalcMass`
    /// represented the bare peptide, not the glycopeptide.
    #[test]
    fn novel_glycan_hit_preserves_residual_mass() {
        // Empty glycan list makes the "no annotation" premise deterministic and
        // independent of the shipped DB contents (CodeRabbit): nothing can ever
        // annotate, so every hit is DeNovo and must carry its residual.
        let glycans: Vec<GlycanComp> = Vec::new();
        let proton = PROTON;
        let steps = crate::glycan_mass::CORE_Y_STEPS;

        // A glycan mass off any plausible composition grid.
        let novel_glycan_mass = 1234.5678_f64;

        let true_backbone_residue = 1500.0_f64;
        let precursor = true_backbone_residue + novel_glycan_mass;
        let y0_neutral = true_backbone_residue + H2O;

        // Full core-Y ladder so the Y-first solver anchors the backbone.
        let mut peaks: Vec<(f64, f32)> = vec![
            (204.08665, 200.0),
            (138.05496, 150.0),
            (186.07608, 80.0),
            (y0_neutral + proton, 100.0),
        ];
        for &s in steps.iter() {
            peaks.push((y0_neutral + s + proton, 90.0));
        }

        let hits = hybrid_candidates(&peaks, precursor, 2, &glycans, 20.0, 10);
        let novel = hits
            .iter()
            .filter(|h| h.source == Source::DeNovo && h.glycan.is_none())
            .min_by(|a, b| {
                (a.backbone_mass - true_backbone_residue)
                    .abs()
                    .total_cmp(&(b.backbone_mass - true_backbone_residue).abs())
            })
            .expect("expected a novel (unannotated) DeNovo hit near the true backbone");

        // The residual glycan mass must equal precursor − backbone (the observed
        // intact glycan), NOT 0.0.
        assert!(
            (novel.glycan_mass_residual - (precursor - novel.backbone_mass)).abs() < 1e-6,
            "residual must be precursor − backbone; got {}",
            novel.glycan_mass_residual
        );
        assert!(
            (novel.glycan_mass_residual - novel_glycan_mass).abs() < 0.02,
            "residual must recover the novel glycan mass ~{novel_glycan_mass}, got {}",
            novel.glycan_mass_residual
        );
    }

    /// `db_branch` must record the (charge, isotope_offset) it was called
    /// with on every returned `BackboneHit` (BUG: precursor charge silently
    /// dropped / isotope offsets ignored).
    #[test]
    fn db_branch_records_charge_and_isotope_offset() {
        let glycans = n_glycan_list();
        let glycan_mass = 2.0 * HEXNAC + 3.0 * HEX;
        let true_backbone = 1500.0_f64;
        let precursor = true_backbone + glycan_mass;

        let hits = db_branch(precursor, &glycans, 500.0, 3, -1);
        assert!(!hits.is_empty());
        for h in &hits {
            assert_eq!(h.charge, 3, "charge must be threaded onto BackboneHit");
            assert_eq!(h.isotope_offset, -1, "isotope_offset must be threaded onto BackboneHit");
        }
    }

    /// Finding #1 regression (Codex adversarial review of the Y-ion-first
    /// rewrite): an oxonium-positive spectrum whose core-Y ladder is too WEAK
    /// for the primary quorum (only a single core-Y rung present, below the ≥2
    /// primary quorum — common in HCD, the sparse stratum) must STILL yield the
    /// true backbone via the quorum-1 confidence-cascade rescue. The pure
    /// Y-ion-first path silently returned an empty Vec here, dropping the whole
    /// weak-ladder tail.
    ///
    /// (Spectra with NO core-Y ladder at all remain unsearchable by this
    /// evidence-driven method — a documented limitation, not a silent drop; a
    /// blind precursor−glycan brute force is O(glycans)-slow and precision-poor.)
    #[test]
    fn weak_ladder_spectrum_is_rescued_by_quorum1_cascade() {
        let glycans = n_glycan_list();
        // HexNAc2Hex3 (trimannosyl core) — present in the list.
        let glycan_mass = 2.0 * HEXNAC + 3.0 * HEX;
        let true_backbone_residue = 1500.0_f64;
        let precursor = true_backbone_residue + glycan_mass;
        let y0_neutral = true_backbone_residue + H2O;
        let proton = PROTON;

        // Oxonium peaks so the glyco gate fires + a SINGLE core-Y rung (Y0 only,
        // no Y1..Y5) → the primary quorum-2 solver finds nothing, but the
        // quorum-1 rescue anchors on the lone Y0.
        let peaks: Vec<(f64, f32)> = vec![
            (204.08665, 200.0),           // HexNAc oxonium
            (138.05496, 150.0),           // HexNAc fragment
            (186.07608, 80.0),            // HexNAc ring-open
            (y0_neutral + proton, 120.0), // Y0 (the single core-Y rung)
        ];

        // Premise: the primary quorum-2 solver alone recovers nothing (1 rung < 2).
        let primary = crate::backbone::solve_backbone(&peaks, precursor, 2, 20.0, 50);
        assert!(
            primary.is_empty(),
            "premise: quorum-2 solver should find nothing with a single core-Y rung, got {:?}",
            primary
        );

        let hits = hybrid_candidates(&peaks, precursor, 2, &glycans, 20.0, 50);
        let found = hits.iter().any(|h| {
            (h.backbone_mass - true_backbone_residue).abs() < 0.05 && h.source == Source::Db
        });
        assert!(
            found,
            "quorum-1 cascade must recover the true backbone for a weak-ladder \
             known-glycan spectrum; got {:?}",
            hits
        );
    }

    /// `hybrid_candidates_with_isotope` must thread a non-zero isotope offset
    /// onto both DB and DeNovo hits (BUG 1: isotope offsets ignored).
    #[test]
    fn hybrid_candidates_with_isotope_threads_offset_onto_all_hits() {
        let glycans = n_glycan_list();
        let proton = PROTON;
        let steps = crate::glycan_mass::CORE_Y_STEPS;

        let glycan_mass = 2.0 * HEXNAC + 3.0 * HEX;
        let true_backbone = 1500.0_f64;
        let precursor = true_backbone + glycan_mass;

        let mut peaks: Vec<(f64, f32)> = vec![
            (204.08665, 200.0),
            (138.05496, 150.0),
            (186.07608, 80.0),
            (true_backbone + proton, 100.0),
        ];
        for &s in steps.iter() {
            peaks.push((true_backbone + s + proton, 90.0));
        }

        let hits = hybrid_candidates_with_isotope(&peaks, precursor, 2, 2, &glycans, 20.0, 10);
        assert!(!hits.is_empty(), "expected hybrid hits");
        for h in &hits {
            assert_eq!(h.isotope_offset, 2, "every hit must carry the caller's isotope offset");
            assert_eq!(h.charge, 2, "every hit must carry the caller's charge");
        }
    }
}
