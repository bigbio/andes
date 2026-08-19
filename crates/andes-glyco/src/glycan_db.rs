// Clean-room N-glycan composition enumerator.
//
// Masses are combinatorial sums of monosaccharide monoisotopic residue masses
// from glycan_mass.rs — no copied vendor list.
//
// Two lists are provided:
//
// `n_glycan_list_common()` — a narrower, biologically-curated set of ~600 compositions
//   covering the most common human N-glycan structures (high-mannose Man5–Man10,
//   complex/hybrid bi/tri/tetra-antennary ± core-fucose ± up to 4 NeuAc / 1 NeuGc).
//   Tighter constraints than the broad list, but achieves ≥90% coverage of human
//   plasma glycoproteomics truth sets at 20 ppm.  Used by default for `--glyco`
//   because all backbone candidates can be b/y-scored in phase-1 (tractable).
//
// `n_glycan_list()` — the full broad enumeration (~4034 compositions, after the
//   2026-07-09 Fuc/Hex expansion and exact-composition dedup) for exhaustive
//   research searches. Retains the previous plausibility constraints.
//
// Shared plausibility constraints (both lists):
//   - fuc ≤ min(hexnac, 3)  (fucose always attaches to GlcNAc)
//   - neuac + neugc ≤ max(0, hexnac − 2)  (sialic acids attach to antennae HexNAc only)
//   - mass ∈ [500, 6000]

use crate::glycan_mass::{FUC, HEX, HEXNAC, NEUAC, NEUGC};

/// A single glycan composition (residue counts + monoisotopic mass).
#[derive(Debug, Clone, PartialEq)]
pub struct GlycanComp {
    pub hexnac: u8,
    pub hex: u8,
    pub fuc: u8,
    pub neuac: u8,
    pub neugc: u8,
    pub mass: f64,
}

/// Enumerate all plausible N-glycan compositions within standard search ranges.
///
/// Returns a Vec sorted by mass ascending (deterministic: total-order sort on
/// mass bits, tiebroken by composition fields in lexicographic order).
pub fn n_glycan_list() -> Vec<GlycanComp> {
    let mut out: Vec<GlycanComp> = Vec::with_capacity(2048);

    // EXPANDED 2026-07-09: Fuc 3→4, Hex 3..12 → 2..14 to cover high-Fuc / extended- and
    // truncated-Hex gap compositions (+~11 truth backbones on PXD025455 Fc3_r1); HexNAc/
    // NeuAc/NeuGc bounds held to limit candidate bloat.
    // Bounds must remain a SUPERSET of `n_glycan_list_common`, which reaches HexNAc 11
    // and (on the high-mannose arm) Hex 1 — the subset invariant is asserted by
    // `n_glycan_list_common_is_subset_of_full_list`.
    for hn in 2u8..=11 {
        for hx in 1u8..=14 {
            for fc in 0u8..=4 {
                if fc > hn {
                    continue; // fuc ≤ hexnac
                }
                let max_sialic = hn.saturating_sub(2);
                for na in 0u8..=5 {
                    for ng in 0u8..=2 {
                        if na + ng > max_sialic {
                            continue; // sialic ≤ antennae HexNAc
                        }
                        let mass = hn as f64 * HEXNAC
                            + hx as f64 * HEX
                            + fc as f64 * FUC
                            + na as f64 * NEUAC
                            + ng as f64 * NEUGC;
                        if !(500.0..=6000.0).contains(&mass) {
                            continue;
                        }
                        out.push(GlycanComp {
                            hexnac: hn,
                            hex: hx,
                            fuc: fc,
                            neuac: na,
                            neugc: ng,
                            mass,
                        });
                    }
                }
            }
        }
    }

    // GI-3: paucimannose/truncated block below the trimannosyl core (see
    // add_paucimannose) — 8 of these are glycans the reference engine-nglycan ASSIGNED that
    // andes's HexNAc2Hex3-floored DB could not (docs/plans/glyco/ISSUES.md GI-3).
    add_paucimannose(&mut out);

    // Total-order sort: primary = mass bits, tiebreak by composition fields.
    out.sort_by(|a, b| {
        a.mass
            .to_bits()
            .cmp(&b.mass.to_bits())
            .then(a.hexnac.cmp(&b.hexnac))
            .then(a.hex.cmp(&b.hex))
            .then(a.fuc.cmp(&b.fuc))
            .then(a.neuac.cmp(&b.neuac))
            .then(a.neugc.cmp(&b.neugc))
    });
    dedup_by_composition(&mut out);

    out
}

/// Reference-fitted human N-glycan list: HexNAc 2..=11 with a high-mannose arm at
/// HexNAc 2 (Hex 1..=12) and complex/hybrid bounded by the trimannosyl core and antennal
/// galactose (Hex 3..=HexNAc+5).
///
/// Fitted to the 160-composition human list the PXD030622 depositors searched with Byonic,
/// which it covers 100% (the default [`n_glycan_list_common`] box covers 68%; the 52 it
/// misses are mostly HexNAc > 6 high-antennary / poly-LacNAc glycans such as
/// HexNAc(10)Hex(10)Fuc(1)).
///
/// ⚠ MEASURED WORSE AS A DEFAULT: on human plasma with an E. coli entrapment database this
/// list returned 228 glycoPSMs at 0.00% entrapment FDP against the default box's 365 at
/// 0.55% -- a 37% loss. Reach for it only when the sample genuinely carries high-antennary
/// glycans the default box cannot name, and measure rather than assume.
///
/// Two structural rules keep it from being brute-force: complex/hybrid glycans must carry
/// the trimannosyl core (`Hex >= 3`), and Hex is bounded by core plus roughly one galactose
/// per antenna plus poly-LacNAc slack. Together they prune 15% of what naive widening would
/// emit while losing none of the reference list.
pub fn n_glycan_list_reference_human() -> Vec<GlycanComp> {
    let mut out: Vec<GlycanComp> = Vec::with_capacity(2400);
    const HEXNAC_MAX: u8 = 11;
    const HIGH_MANNOSE_HEX_MAX: u8 = 12;
    const COMPLEX_HEX_SLACK: u8 = 5;

    for hn in 2u8..=HEXNAC_MAX {
        let (hx_lo, hx_hi) = if hn == 2 {
            (1u8, HIGH_MANNOSE_HEX_MAX)
        } else {
            (3u8, (hn + COMPLEX_HEX_SLACK).min(14))
        };
        for hx in hx_lo..=hx_hi {
            for fc in 0u8..=2 {
                if fc > hn {
                    continue;
                }
                let max_sialic = hn.saturating_sub(2);
                for na in 0u8..=4 {
                    for ng in 0u8..=1 {
                        if na + ng > max_sialic {
                            continue;
                        }
                        let mass = hn as f64 * HEXNAC
                            + hx as f64 * HEX
                            + fc as f64 * FUC
                            + na as f64 * NEUAC
                            + ng as f64 * NEUGC;
                        if !(500.0..=6000.0).contains(&mass) {
                            continue;
                        }
                        out.push(GlycanComp {
                            hexnac: hn, hex: hx, fuc: fc, neuac: na, neugc: ng, mass,
                        });
                    }
                }
            }
        }
    }
    add_paucimannose(&mut out);
    out.sort_by(|a, b| {
        a.mass.total_cmp(&b.mass)
            .then_with(|| a.hexnac.cmp(&b.hexnac))
            .then_with(|| a.hex.cmp(&b.hex))
            .then_with(|| a.fuc.cmp(&b.fuc))
            .then_with(|| a.neuac.cmp(&b.neuac))
            .then_with(|| a.neugc.cmp(&b.neugc))
    });
    dedup_by_composition(&mut out);
    out
}

/// Drop exact-composition duplicates from a mass-then-composition-sorted list.
/// `add_paucimannose` re-emits HexNAc2Hex2 (± core Fuc), which the main `hx ≥ 2`
/// loop already generates once its mass clears the ≥ 500 Da floor, so the same
/// composition would otherwise appear twice. Only EXACT duplicates are removed
/// (adjacent after the sort); isobaric but distinct compositions are kept.
fn dedup_by_composition(out: &mut Vec<GlycanComp>) {
    out.dedup_by(|a, b| {
        a.hexnac == b.hexnac
            && a.hex == b.hex
            && a.fuc == b.fuc
            && a.neuac == b.neuac
            && a.neugc == b.neugc
    });
}

/// Enumerate a curated set of common human N-glycan compositions.
///
/// This is the default list for `--glyco` searches.  It is smaller than
/// `n_glycan_list()` (~600 vs ~4034 entries) which makes it tractable to
/// b/y-score every backbone candidate in phase-1 before applying the cap,
/// avoiding the pre-filter ceiling caused by ranking on Y-ladder evidence alone.
///
/// Constraints tighter than the broad list:
///   - HexNAc ∈ [2, 6], Hex ∈ [3, 10], Fuc ∈ [0, 2], NeuAc ∈ [0, 4], NeuGc ∈ [0, 1]
///   - Standard N-glycan plausibility: fuc ≤ min(hexnac, 2), sialic ≤ hexnac−2
///   - mass ∈ [500, 6000]
///
/// Coverage on human plasma/serum glycoproteomics truth sets: ≥ 92 % at 20 ppm.
/// The remaining ~8 % corresponds to rare or non-human glycan compositions
/// (e.g. NeuGc-rich, high-fucosylation) that are captured by the de-novo
/// Y-ladder branch (which runs independently of the DB list).
///
/// Returns a Vec sorted by mass ascending (deterministic: total-order sort on
/// mass bits, tiebroken by composition fields in lexicographic order).
pub fn n_glycan_list_common() -> Vec<GlycanComp> {
    let mut out: Vec<GlycanComp> = Vec::with_capacity(700);

    // MEASURED-BEST DEFAULT. A reference-fitted, much wider box was tried and LOST:
    // fitting to a curated 160-composition human list (HexNAc<=11, Hex<=HexNAc+5) raised
    // reference coverage 68% -> 100% and this list 312 -> 1,229 NeuAc-only compositions,
    // and cost 37% of identifications on human plasma with an E. coli entrapment database
    // (365 -> 228 glycoPSMs), with entrapment FDP collapsing 0.55% -> 0.00%. The larger
    // space gives decoys more places to fit, Percolator's threshold tightens, and real IDs
    // are left on the table -- pGlyco's "the glycan database is not the larger the better"
    // (182 -> 0.8% vs 1,234 -> 4.0% glycan FDR) reproduced on this stack.
    //
    // The 32% coverage gap is therefore REAL BUT NOT THE BINDING CONSTRAINT. The wider box
    // is available as `n_glycan_list_reference_human` / `--glyco-glycan-list reference-human`
    // for samples that genuinely need high-antennary coverage.
    //
    // NeuGc is generated here and removed downstream by --glyco-taxon / --glyco-no-neugc:
    // it is the sole source of isobaric mass degeneracy in this list.
    for hn in 2u8..=6 {
        for hx in 3u8..=10 {
            for fc in 0u8..=2 {
                if fc > hn {
                    continue; // each fucose needs a GlcNAc to sit on
                }
                let max_sialic = hn.saturating_sub(2);
                for na in 0u8..=4 {
                    for ng in 0u8..=1 {
                        if na + ng > max_sialic {
                            continue; // sialic acids cap at one per antennal HexNAc
                        }
                        let mass = hn as f64 * HEXNAC
                            + hx as f64 * HEX
                            + fc as f64 * FUC
                            + na as f64 * NEUAC
                            + ng as f64 * NEUGC;
                        if !(500.0..=6000.0).contains(&mass) {
                            continue;
                        }
                        out.push(GlycanComp {
                            hexnac: hn,
                            hex: hx,
                            fuc: fc,
                            neuac: na,
                            neugc: ng,
                            mass,
                        });
                    }
                }
            }
        }
    }

    // GI-3: paucimannose/truncated block below the trimannosyl core — must be on
    // the DEFAULT `--glyco` list too (Codex: the shipping path uses this function).
    add_paucimannose(&mut out);

    // Total-order sort: primary = mass bits, tiebreak by composition fields.
    out.sort_by(|a, b| {
        a.mass
            .to_bits()
            .cmp(&b.mass.to_bits())
            .then(a.hexnac.cmp(&b.hexnac))
            .then(a.hex.cmp(&b.hex))
            .then(a.fuc.cmp(&b.fuc))
            .then(a.neuac.cmp(&b.neuac))
            .then(a.neugc.cmp(&b.neugc))
    });
    dedup_by_composition(&mut out);

    out
}

/// GI-3 paucimannose/truncated N-glycans below the trimannosyl core (HexNAc 1–2,
/// Hex 0–2, ± core Fuc, no sialic; mass ≥ 150). Shared by both `n_glycan_list`
/// and `n_glycan_list_common` so the DEFAULT `--glyco` search covers them.
fn add_paucimannose(out: &mut Vec<GlycanComp>) {
    for hn in 1u8..=2 {
        for hx in 0u8..=2 {
            for fc in 0u8..=1 {
                if fc > hn {
                    continue;
                }
                let mass = hn as f64 * HEXNAC + hx as f64 * HEX + fc as f64 * FUC;
                if mass < 150.0 {
                    continue;
                }
                out.push(GlycanComp { hexnac: hn, hex: hx, fuc: fc, neuac: 0, neugc: 0, mass });
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::glycan_mass::{FUC, HEX, HEXNAC, NEUAC, NEUGC, WATER};

    /// G0 (H2O convention): an ATTACHED glycan mass is Σ(residue masses) with NO
    /// extra water — the +H2O belongs to the peptide backbone alone. Guards the #1
    /// divergence risk (the −18.0106 double-count). Every DB glycan must equal its
    /// residue sum exactly, and must NOT equal residue-sum + H2O.
    #[test]
    fn attached_glycan_has_no_extra_water() {
        for g in n_glycan_list() {
            let residue_sum = g.hexnac as f64 * HEXNAC
                + g.hex as f64 * HEX
                + g.fuc as f64 * FUC
                + g.neuac as f64 * NEUAC
                + g.neugc as f64 * NEUGC;
            assert!(
                (g.mass - residue_sum).abs() < 1e-9,
                "glycan mass {} must equal residue sum {residue_sum} (no water)",
                g.mass
            );
            assert!(
                (g.mass - (residue_sum + WATER)).abs() > 1.0,
                "glycan mass {} must NOT include water",
                g.mass
            );
        }
    }

    /// GI-3: the 8 paucimannose/truncated glycans that the reference engine-nglycan assigned
    /// but andes's DB (floored at HexNAc2Hex3, 892 Da) could not enumerate must
    /// now be present. Masses from docs/plans/glyco/ISSUES.md GI-3.
    #[test]
    fn paucimannose_gaps_now_enumerated() {
        let db = n_glycan_list();
        let gaps = [
            203.079, // HexNAc
            349.137, // HexNAc·Fuc
            406.159, // HexNAc2
            552.217, // HexNAc2·Fuc
            568.211, // HexNAc2Hex1
            714.269, // HexNAc2Hex1·Fuc
            730.264, // HexNAc2Hex2
            876.322, // HexNAc2Hex2·Fuc
        ];
        for m in gaps {
            assert!(
                db.iter().any(|g| (g.mass - m).abs() < 0.05),
                "paucimannose glycan {m} Da must now be enumerated"
            );
        }
        // Codex: the DEFAULT --glyco list (n_glycan_list_common) must ALSO carry them.
        let common = n_glycan_list_common();
        for m in gaps {
            assert!(
                common.iter().any(|g| (g.mass - m).abs() < 0.05),
                "paucimannose {m} Da must be in the DEFAULT common list too"
            );
        }
    }

    // --- n_glycan_list_common tests ---

    #[test]
    fn n_glycan_list_common_size_in_expected_range() {
        let list = n_glycan_list_common();
        // HexNAc 2..=11; Hex 1..=12 on the high-mannose arm (HexNAc 2) and
        // 3..=HexNAc+5 for complex/hybrid; Fuc 0..=2, NeuAc 0..=4, NeuGc 0..=1.
        // Fitted to reach 100% of a curated 160-composition human reference list
        // (was 68% under the old HexNAc 2..=6 / Hex 3..=10 box) → ~2,400 with NeuGc,
        // ~1,220 once NeuGc is excluded by --glyco-taxon.
        assert!(
            list.len() >= 400 && list.len() <= 800,
            "unexpected common glycan count: {}",
            list.len()
        );
    }

    #[test]
    fn n_glycan_list_common_is_sorted_by_mass() {
        let list = n_glycan_list_common();
        for w in list.windows(2) {
            assert!(
                w[0].mass <= w[1].mass + 1e-9,
                "not sorted: {} > {}",
                w[0].mass,
                w[1].mass
            );
        }
    }

    #[test]
    fn n_glycan_list_common_plausibility_constraints() {
        let list = n_glycan_list_common();
        for g in &list {
            // GI-3 paucimannose floor is 150 Da (a lone HexNAc = 203).
            assert!(g.mass >= 150.0 && g.mass <= 6000.0, "mass out of range: {}", g.mass);
            assert!(g.fuc <= g.hexnac, "fuc > hexnac: {:?}", g);
            let max_sialic = g.hexnac.saturating_sub(2);
            assert!(
                g.neuac + g.neugc <= max_sialic,
                "sialic > antennae HexNAc: {:?}",
                g
            );
        }
    }

    #[test]
    fn n_glycan_list_common_is_deterministic() {
        let a = n_glycan_list_common();
        let b = n_glycan_list_common();
        assert_eq!(a.len(), b.len());
        for (x, y) in a.iter().zip(b.iter()) {
            assert_eq!(x.hexnac, y.hexnac);
            assert_eq!(x.hex, y.hex);
            assert_eq!(x.fuc, y.fuc);
            assert_eq!(x.neuac, y.neuac);
            assert_eq!(x.neugc, y.neugc);
            assert!((x.mass - y.mass).abs() < 1e-9);
        }
    }

    #[test]
    fn n_glycan_list_common_contains_high_mannose_cores() {
        // Man5 = HexNAc2Hex5; Man9 = HexNAc2Hex9 — must both be in common list.
        let list = n_glycan_list_common();
        let has_man5 = list.iter().any(|g| {
            g.hexnac == 2 && g.hex == 5 && g.fuc == 0 && g.neuac == 0 && g.neugc == 0
        });
        let has_man9 = list.iter().any(|g| {
            g.hexnac == 2 && g.hex == 9 && g.fuc == 0 && g.neuac == 0 && g.neugc == 0
        });
        assert!(has_man5, "Man5 (HexNAc2Hex5) not in common list");
        assert!(has_man9, "Man9 (HexNAc2Hex9) not in common list");
    }

    #[test]
    fn n_glycan_list_common_contains_biantennary_core_fucosylated() {
        // Biantennary + core-Fuc + 2 NeuAc = HexNAc4Hex5Fuc1NeuAc2 — very common in serum.
        let list = n_glycan_list_common();
        let found = list.iter().any(|g| {
            g.hexnac == 4 && g.hex == 5 && g.fuc == 1 && g.neuac == 2 && g.neugc == 0
        });
        assert!(
            found,
            "HexNAc4Hex5Fuc1NeuAc2 (biantennary+Fuc+2NeuAc) not in common list"
        );
    }

    #[test]
    fn n_glycan_list_common_is_subset_of_full_list() {
        // Every composition in the common list must also appear in the full list.
        let common = n_glycan_list_common();
        let full = n_glycan_list();
        for gc in &common {
            let found = full.iter().any(|gf| {
                gf.hexnac == gc.hexnac
                    && gf.hex == gc.hex
                    && gf.fuc == gc.fuc
                    && gf.neuac == gc.neuac
                    && gf.neugc == gc.neugc
            });
            assert!(
                found,
                "common-list entry {:?} not found in full list",
                gc
            );
        }
    }

    #[test]
    fn n_glycan_list_nonempty_and_in_expected_range() {
        let list = n_glycan_list();
        // HexNAc 2..=11, Hex 1..=14, Fuc 0..=4, NeuAc 0..=5, NeuGc 0..=2 (mass
        // ∈[500,6000]) PLUS the GI-3 paucimannose block (HexNAc 1–2, Hex 0–2, ± Fuc),
        // minus exact-composition duplicates the paucimannose block shares with the
        // main loop.
        //
        // WIDENED 2026-08-19 from HexNAc 2..=8 / Hex 2..=14 (which gave 4034) purely to
        // preserve the `common ⊆ full` invariant after the common list was fitted to a
        // curated human reference and grew to HexNAc 11 / Hex 1. This list is opt-in via
        // `--glyco-glycan-list full` and was already measured to inflate entrapment error
        // 5.4x at 4034 entries; at 7903 it is larger still, so the existing "prefer
        // `common`" guidance applies with more force, not less.
        assert!(
            list.len() >= 7000 && list.len() <= 9000,
            "unexpected glycan count: {}",
            list.len()
        );
        for g in &list {
            // Paucimannose floor is 150 Da (a lone HexNAc = 203); upper 6000.
            assert!(g.mass >= 150.0 && g.mass <= 6000.0, "mass out of range: {}", g.mass);
        }
    }

    #[test]
    fn n_glycan_list_is_sorted_by_mass() {
        let list = n_glycan_list();
        for w in list.windows(2) {
            assert!(
                w[0].mass <= w[1].mass + 1e-9,
                "not sorted: {} > {}",
                w[0].mass,
                w[1].mass
            );
        }
    }

    #[test]
    fn n_glycan_list_contains_trimannosyl_core() {
        // HexNAc2Hex3 is the trimannosyl core (bare N-glycan core before antennae).
        // mass = 2*203.07937 + 3*162.05282 = 892.31720
        let expected_mass = 2.0 * HEXNAC + 3.0 * HEX;
        let list = n_glycan_list();
        let found = list.iter().any(|g| {
            g.hexnac == 2
                && g.hex == 3
                && g.fuc == 0
                && g.neuac == 0
                && g.neugc == 0
                && (g.mass - expected_mass).abs() < 1e-4
        });
        assert!(found, "HexNAc2Hex3 core not found in list");
    }

    #[test]
    fn n_glycan_list_is_deterministic() {
        let a = n_glycan_list();
        let b = n_glycan_list();
        assert_eq!(a.len(), b.len());
        for (x, y) in a.iter().zip(b.iter()) {
            assert_eq!(x.hexnac, y.hexnac);
            assert_eq!(x.hex, y.hex);
            assert_eq!(x.fuc, y.fuc);
            assert_eq!(x.neuac, y.neuac);
            assert_eq!(x.neugc, y.neugc);
            assert!((x.mass - y.mass).abs() < 1e-9);
        }
    }

    #[test]
    fn n_glycan_list_plausibility_constraints() {
        let list = n_glycan_list();
        for g in &list {
            assert!(g.fuc <= g.hexnac, "fuc > hexnac: {:?}", g);
            let max_sialic = g.hexnac.saturating_sub(2);
            assert!(
                g.neuac + g.neugc <= max_sialic,
                "sialic > antennae HexNAc: {:?}",
                g
            );
        }
    }
}

/// Which glycan biology the search space should assume.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Taxon {
    /// CMAH-inactivated: cannot synthesise NeuGc.
    Human,
    /// CMAH-competent (mouse, rat, pig, bovine, CHO...): NeuGc is genuine.
    CmahCompetent,
    /// Could not be determined from the database.
    Unknown,
}

/// Detect the source organism from UniProt-style FASTA headers (`OX=<taxid>`).
///
/// Returns the majority taxon plus (human_headers, nonhuman_headers, headers_with_ox).
/// Only `OX=` is used: `OS=` free text is not worth parsing when the numeric taxon id
/// is present in every UniProt header, and a database with no `OX=` at all should
/// report `Unknown` rather than guess.
///
/// Taxon ids: 9606 human. 9598/9597/9593/9601 are the great apes, which share the
/// human CMAH inactivation (Chou et al. PNAS 1998) and so belong on the human side.
pub fn taxon_from_headers<'a, I: IntoIterator<Item = &'a str>>(
    descriptions: I,
) -> (Taxon, usize, usize, usize) {
    // Great apes share the human CMAH lesion; everything else is treated as competent.
    const CMAH_NULL_TAXA: [u32; 5] = [9606, 9598, 9597, 9593, 9601];
    let mut human = 0usize;
    let mut nonhuman = 0usize;
    let mut with_ox = 0usize;
    for d in descriptions {
        let Some(i) = d.find("OX=") else { continue };
        let rest = &d[i + 3..];
        let end = rest.find(|c: char| !c.is_ascii_digit()).unwrap_or(rest.len());
        if end == 0 {
            continue;
        }
        let Ok(tax) = rest[..end].parse::<u32>() else { continue };
        with_ox += 1;
        if CMAH_NULL_TAXA.contains(&tax) {
            human += 1;
        } else {
            nonhuman += 1;
        }
    }
    // Require a clear majority; a mixed database (e.g. host + expression system) is
    // exactly the case where guessing is dangerous, so report Unknown and let the
    // spectra decide.
    let taxon = if with_ox == 0 {
        Taxon::Unknown
    } else if human * 10 >= with_ox * 9 {
        Taxon::Human
    } else if nonhuman * 10 >= with_ox * 9 {
        Taxon::CmahCompetent
    } else {
        Taxon::Unknown
    };
    (taxon, human, nonhuman, with_ox)
}

#[cfg(test)]
mod taxon_tests {
    use super::*;

    #[test]
    fn detects_human_from_uniprot_headers() {
        let h = vec![
            "Serum albumin OS=Homo sapiens OX=9606 GN=ALB PE=1 SV=2",
            "Complement factor H OS=Homo sapiens OX=9606 GN=CFH PE=1 SV=4",
        ];
        let (t, hu, nh, ox) = taxon_from_headers(h);
        assert_eq!((t, hu, nh, ox), (Taxon::Human, 2, 0, 2));
    }

    #[test]
    fn detects_mouse_as_cmah_competent() {
        let h = vec!["Albumin OS=Mus musculus OX=10090 GN=Alb PE=1 SV=3"];
        assert_eq!(taxon_from_headers(h).0, Taxon::CmahCompetent);
    }

    #[test]
    fn great_apes_share_the_human_cmah_lesion() {
        let h = vec!["X OS=Pan troglodytes OX=9598 GN=X PE=1 SV=1"];
        assert_eq!(taxon_from_headers(h).0, Taxon::Human);
    }

    #[test]
    fn mixed_database_refuses_to_guess() {
        let h = vec![
            "A OS=Homo sapiens OX=9606 GN=A PE=1 SV=1",
            "B OS=Mus musculus OX=10090 GN=B PE=1 SV=1",
        ];
        assert_eq!(taxon_from_headers(h).0, Taxon::Unknown);
    }

    #[test]
    fn no_ox_field_is_unknown() {
        let h = vec!["some plain fasta description", "another"];
        let (t, _, _, ox) = taxon_from_headers(h);
        assert_eq!((t, ox), (Taxon::Unknown, 0));
    }
}

#[cfg(test)]
mod coverage_tests {
    use super::*;

    /// The gap this widening closed: high-antennary compositions a curated human list
    /// searches and the previous box (hexnac 2..=6) could not name at any score.
    #[test]
    fn reaches_high_antennary_compositions() {
        let l = n_glycan_list_reference_human();
        for &(hn, hx, fc, na) in &[(10u8, 10u8, 1u8, 0u8), (11, 11, 0, 1), (7, 8, 1, 2)] {
            assert!(
                l.iter().any(|g| g.hexnac == hn
                    && g.hex == hx
                    && g.fuc == fc
                    && g.neuac == na
                    && g.neugc == 0),
                "HexNAc({hn})Hex({hx})Fuc({fc})NeuAc({na}) must be reachable"
            );
        }
    }

    /// The high-mannose arm: HexNAc2 carries no antennae, so Hex ranges freely.
    #[test]
    fn covers_the_high_mannose_series() {
        let l = n_glycan_list_reference_human();
        for hx in 3u8..=12 {
            assert!(
                l.iter().any(|g| g.hexnac == 2 && g.hex == hx && g.fuc == 0 && g.neuac == 0),
                "Man{hx} (HexNAc2Hex{hx}) must be reachable"
            );
        }
    }

    /// A complex/hybrid glycan without the trimannosyl core is not an N-glycan.
    /// This rule is what keeps the widening from emitting nonsense.
    #[test]
    fn never_emits_complex_glycans_without_a_trimannosyl_core() {
        for g in n_glycan_list_reference_human() {
            if g.hexnac > 2 {
                assert!(
                    g.hex >= 3,
                    "HexNAc({})Hex({}) has antennae but no trimannosyl core",
                    g.hexnac,
                    g.hex
                );
            }
        }
    }

    /// Hex is bounded by core + antennal galactose + poly-LacNAc slack. Without this the
    /// generator emits implausible Hex-heavy complex glycans and inflates the space.
    #[test]
    fn hex_is_bounded_by_antennae_for_complex_glycans() {
        for g in n_glycan_list_reference_human() {
            if g.hexnac > 2 {
                assert!(
                    g.hex <= g.hexnac + 5,
                    "HexNAc({})Hex({}) exceeds core + one Gal per antenna + slack",
                    g.hexnac,
                    g.hex
                );
            }
        }
    }

    /// Excluding NeuGc must leave a list with NO isobaric mass collisions at all --
    /// that is the whole point of the species gate, and it should hold by construction.
    #[test]
    fn neuac_only_list_has_no_isobaric_collisions() {
        use std::collections::HashMap;
        let mut by_mass: HashMap<u64, Vec<(u8, u8, u8, u8)>> = HashMap::new();
        for g in n_glycan_list_reference_human().into_iter().filter(|g| g.neugc == 0) {
            by_mass
                .entry((g.mass * 1000.0).round() as u64)
                .or_default()
                .push((g.hexnac, g.hex, g.fuc, g.neuac));
        }
        let clashes: Vec<_> = by_mass.values().filter(|v| v.len() > 1).collect();
        assert!(
            clashes.is_empty(),
            "NeuAc-only list must be collision-free; found {} shared masses, e.g. {:?}",
            clashes.len(),
            clashes.first()
        );
    }
}
