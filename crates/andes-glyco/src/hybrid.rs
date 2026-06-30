// Hybrid backbone candidate generation: DB-branch ∪ de-novo Y-ladder.
//
// DB branch: for each known glycan composition, `backbone = precursor_neutral − glycan.mass`.
// De-novo branch: the existing Y-ladder complement-pair solver (backbone.rs).
// Union: merge candidates within 0.02 Da, preferring DB source over de-novo.

use crate::backbone::solve_backbone;
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
}

/// DB-branch backbone enumeration.
///
/// For each glycan in `glycans`, compute `bb = precursor_neutral − glycan.mass`.
/// Keep if `bb ≥ min_backbone`. Returns candidates sorted by backbone_mass ascending.
pub fn db_branch(
    precursor_neutral: f64,
    glycans: &[GlycanComp],
    min_backbone: f64,
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
    const MIN_BACKBONE: f64 = 500.0;

    // --- DB branch (always; doesn't require oxonium evidence) ---
    let mut combined: Vec<BackboneHit> = db_branch(precursor_neutral, glycans, MIN_BACKBONE);

    // --- De-novo branch (requires oxonium gate to have fired) ---
    let ox = oxonium_gate(peaks, 0.10, tol_ppm);
    if ox.fired {
        let dn = solve_backbone(peaks, precursor_neutral, precursor_z, tol_ppm, top_k);
        for c in dn {
            combined.push(BackboneHit {
                backbone_mass: c.backbone_mass,
                glycan: None,
                source: Source::DeNovo,
            });
        }
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

    // --- Dedup: merge candidates within max(bb*20e-6, 0.01) Da, prefer Db source ---
    // Window matches the gate/searchable window so dedup and gate agree.
    // Strategy: single-pass, keep a running "cluster representative".
    // Because Db is sorted before DeNovo at equal mass, when a cluster contains
    // both, the first element is always Db → the representative is Db.
    if combined.is_empty() {
        return combined;
    }

    let mut deduped: Vec<BackboneHit> = Vec::with_capacity(combined.len());
    let mut rep = combined.remove(0);

    for next in combined {
        let tol = (rep.backbone_mass * 20e-6_f64).max(0.01);
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

        let hits = db_branch(precursor, &glycans, 500.0);
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
        let hits = db_branch(precursor, &glycans, 500.0);
        for h in &hits {
            assert!(h.backbone_mass >= 500.0, "backbone below min: {}", h.backbone_mass);
        }
    }

    /// DB branch is sorted by backbone_mass ascending.
    #[test]
    fn db_branch_is_sorted() {
        let glycans = n_glycan_list();
        let precursor = 4000.0;
        let hits = db_branch(precursor, &glycans, 500.0);
        for w in hits.windows(2) {
            assert!(
                w[0].backbone_mass <= w[1].backbone_mass + 1e-9,
                "not sorted: {} > {}",
                w[0].backbone_mass,
                w[1].backbone_mass
            );
        }
    }

    /// hybrid_candidates must return hits from both DB and de-novo sources when
    /// the spectrum has a real core-Y ladder (de-novo fires) and the precursor
    /// matches a known glycan (DB fires).
    #[test]
    fn hybrid_union_contains_both_sources() {
        let glycans = n_glycan_list();
        let proton = PROTON;
        let steps = crate::glycan_mass::CORE_Y_STEPS;

        let glycan_mass = 2.0 * HEXNAC + 3.0 * HEX; // HexNAc2Hex3 ~892.317
        let true_backbone = 1500.0_f64;
        let precursor = true_backbone + glycan_mass;

        // Build a synthetic spectrum: oxonium peaks + full core-Y ladder.
        let mut peaks: Vec<(f64, f32)> = vec![
            (204.08665, 200.0), // HexNAc oxonium
            (138.05496, 150.0), // HexNAc fragment
            (186.07608, 80.0),  // HexNAc ring-open
            (true_backbone + proton, 100.0), // Y0
        ];
        for &s in steps.iter() {
            peaks.push((true_backbone + s + proton, 90.0));
        }

        let hits = hybrid_candidates(&peaks, precursor, 2, &glycans, 20.0, 10);
        assert!(!hits.is_empty(), "expected hybrid hits");

        let has_db = hits.iter().any(|h| h.source == Source::Db);
        let has_dn = hits.iter().any(|h| h.source == Source::DeNovo);
        assert!(has_db, "expected at least one DB hit");
        assert!(has_dn, "expected at least one DeNovo hit");
    }

    /// Dedup: when DB and de-novo candidates cluster within 0.02 Da, the Db source wins.
    #[test]
    fn hybrid_dedup_keeps_db_over_denovo() {
        let glycans = n_glycan_list();
        let proton = PROTON;
        let steps = crate::glycan_mass::CORE_Y_STEPS;

        // True backbone matches a DB glycan exactly.
        let glycan_mass = 2.0 * HEXNAC + 3.0 * HEX; // ~892.317
        let true_backbone = 1500.0_f64;
        let precursor = true_backbone + glycan_mass;

        // Full core-Y ladder so de-novo also proposes ~1500 Da.
        let mut peaks: Vec<(f64, f32)> = vec![
            (204.08665, 200.0),
            (138.05496, 150.0),
            (186.07608, 80.0),
            (true_backbone + proton, 100.0),
        ];
        for &s in steps.iter() {
            peaks.push((true_backbone + s + proton, 90.0));
        }

        let hits = hybrid_candidates(&peaks, precursor, 2, &glycans, 20.0, 10);

        // Find candidates near true_backbone; after dedup there must be exactly one
        // in the ±0.02 Da window, and it must be Source::Db.
        let near: Vec<&BackboneHit> = hits
            .iter()
            .filter(|h| (h.backbone_mass - true_backbone).abs() < 0.02)
            .collect();
        assert_eq!(near.len(), 1, "expected exactly one candidate near true backbone after dedup, got {}", near.len());
        assert_eq!(
            near[0].source,
            Source::Db,
            "expected Db source to win dedup, got {:?}",
            near[0].source
        );
    }
}
