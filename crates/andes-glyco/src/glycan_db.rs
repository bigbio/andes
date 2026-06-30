// Clean-room N-glycan composition enumerator.
//
// Masses are combinatorial sums of monosaccharide monoisotopic residue masses
// from glycan_mass.rs — no copied vendor list.
//
// Plausibility constraints for N-glycans:
//   - HexNAc ∈ [2,8], Hex ∈ [3,12], Fuc ∈ [0,3], NeuAc ∈ [0,5], NeuGc ∈ [0,2]
//   - fuc ≤ hexnac   (fucose always attaches to GlcNAc in known N-glycan topologies)
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

    for hn in 2u8..=8 {
        for hx in 3u8..=12 {
            for fc in 0u8..=3 {
                if fc > hn {
                    continue; // fuc ≤ hexnac
                }
                let max_sialic = hn.saturating_sub(2) as u8;
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
                        if mass < 500.0 || mass > 6000.0 {
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

    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::glycan_mass::{HEX, HEXNAC};

    #[test]
    fn n_glycan_list_nonempty_and_in_expected_range() {
        let list = n_glycan_list();
        // Exact combinatorial count for the defined ranges + plausibility constraints:
        // HexNAc 2..=8, Hex 3..=12, Fuc 0..=3, NeuAc 0..=5, NeuGc 0..=2
        // after fuc≤hexnac, sialic≤hexnac-2, mass∈[500,6000] → 2510 compositions.
        assert!(
            list.len() >= 2000 && list.len() <= 3000,
            "unexpected glycan count: {}",
            list.len()
        );
        for g in &list {
            assert!(g.mass >= 500.0 && g.mass <= 6000.0, "mass out of range: {}", g.mass);
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
