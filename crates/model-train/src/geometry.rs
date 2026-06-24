//! Data-derived partition/segment geometry.
//!
//! Builds the scoring [`Param`]'s structural geometry (charge span, per-charge
//! `parent_mass` tier boundaries, segment count, rank cap, and per-partition ion
//! membership) from andes's own labeled corpus, so a trained model owns its
//! geometry rather than inheriting it from a seed template.

use scoring_crate::param_model::Partition;

/// Equal-occupancy `parent_mass` tier lower-bounds for one charge.
///
/// Returns up to `n_tiers` ascending lower-bounds derived as the masses at the
/// `k/n_tiers` quantiles of `masses` (`k = 0..n_tiers`), so each tier holds an
/// approximately equal number of training PSMs. The first bound is the minimum
/// mass (tier 0 catches everything at or above it under the floor lookup).
/// Consecutive duplicate bounds are collapsed, so a highly-degenerate mass
/// distribution can yield fewer than `n_tiers` tiers.
pub fn derive_mass_tiers(masses: &[f32], n_tiers: usize) -> Vec<f32> {
    if masses.is_empty() || n_tiers == 0 {
        return Vec::new();
    }
    let mut sorted: Vec<f32> = masses.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).expect("parent masses must be finite"));
    let n = sorted.len();
    let mut tiers: Vec<f32> = Vec::with_capacity(n_tiers);
    for k in 0..n_tiers {
        let idx = ((k * n) / n_tiers).min(n - 1);
        let bound = sorted[idx];
        // Collapse consecutive duplicate bounds: a zero-width tier would create
        // an unreachable partition under the floor lookup.
        if tiers.last() != Some(&bound) {
            tiers.push(bound);
        }
    }
    tiers
}

/// Build the partition skeleton: the Cartesian product of charges, their
/// `parent_mass` tier lower-bounds, and the `0..num_segments` segments. Returned
/// sorted by the `Partition` lex order (the loader / `find_partition` invariant).
/// The returned partitions have empty learned tables — they define geometry only.
pub fn build_partition_skeleton(
    tiers_by_charge: &[(i32, Vec<f32>)],
    num_segments: i32,
) -> Vec<Partition> {
    let mut parts: Vec<Partition> = Vec::new();
    for &(charge, ref tiers) in tiers_by_charge {
        for &parent_mass in tiers {
            for seg_num in 0..num_segments.max(1) {
                parts.push(Partition { charge, parent_mass, seg_num });
            }
        }
    }
    parts.sort();
    parts
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn equal_occupancy_tiers_split_uniform_distribution() {
        // 8 distinct masses, 4 tiers -> boundaries at indices 0,2,4,6.
        let masses = vec![
            100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0,
        ];
        let tiers = derive_mass_tiers(&masses, 4);
        assert_eq!(tiers, vec![100.0, 300.0, 500.0, 700.0]);
    }

    #[test]
    fn skeleton_is_cartesian_product_sorted_lex() {
        let tiers_by_charge = vec![
            (2, vec![100.0, 500.0]),
            (3, vec![200.0, 600.0]),
        ];
        let parts = build_partition_skeleton(&tiers_by_charge, 2);
        // Canonical Partition Ord is charge -> seg_num -> parent_mass
        // (param_model.rs: load-bearing for find_partition's floor lookup).
        let expect = vec![
            Partition { charge: 2, parent_mass: 100.0, seg_num: 0 },
            Partition { charge: 2, parent_mass: 500.0, seg_num: 0 },
            Partition { charge: 2, parent_mass: 100.0, seg_num: 1 },
            Partition { charge: 2, parent_mass: 500.0, seg_num: 1 },
            Partition { charge: 3, parent_mass: 200.0, seg_num: 0 },
            Partition { charge: 3, parent_mass: 600.0, seg_num: 0 },
            Partition { charge: 3, parent_mass: 200.0, seg_num: 1 },
            Partition { charge: 3, parent_mass: 600.0, seg_num: 1 },
        ];
        assert_eq!(parts, expect);
    }
}
