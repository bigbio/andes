//! Data-derived partition/segment geometry.
//!
//! Builds the scoring [`Param`]'s structural geometry (charge span, per-charge
//! `parent_mass` tier boundaries, segment count, rank cap, and per-partition ion
//! membership) from andes's own labeled corpus, so a trained model owns its
//! geometry rather than inheriting it from a seed template.

use rustc_hash::FxHashMap;
use scoring_crate::param_model::{FragmentOffsetFrequency, IonType, Partition};

/// b-ion m/z offset (proton). y-ion offset adds a water. Both CODATA-sourced
/// in `model::mass` — chemistry, not seed data.
fn b_offset_bits() -> u32 {
    (model::mass::PROTON as f32).to_bits()
}
fn y_offset_bits() -> u32 {
    ((model::mass::H2O + model::mass::PROTON) as f32).to_bits()
}

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

/// Group corpus `(charge, parent_mass)` PSMs by charge and derive each charge's
/// equal-occupancy `parent_mass` tier bounds, returned ascending by charge — the
/// `tiers_by_charge` input to [`build_partition_skeleton`].
pub fn derive_tiers_by_charge(
    charge_masses: &[(i32, f32)],
    n_tiers: usize,
) -> Vec<(i32, Vec<f32>)> {
    let mut by_charge: std::collections::BTreeMap<i32, Vec<f32>> =
        std::collections::BTreeMap::new();
    for &(charge, mass) in charge_masses {
        by_charge.entry(charge).or_default().push(mass);
    }
    by_charge
        .into_iter()
        .map(|(charge, masses)| (charge, derive_mass_tiers(&masses, n_tiers)))
        .collect()
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

/// Build the per-partition ion membership (`frag_off_table`, G6) by chemistry,
/// not by copying a seed: for a partition of precursor charge `C`, emit b
/// (`Prefix`) and y (`Suffix`) ions at fragment charges `1..=max(1, C-1)`,
/// `loss_class 0`, plus a `Noise` entry (RankScorer requires one per populated
/// partition). Frequencies are placeholders here (refined by a later count
/// pass); the ion *set* is what makes the geometry own-derived.
pub fn build_frag_off_table(
    partitions: &[Partition],
) -> FxHashMap<Partition, Vec<FragmentOffsetFrequency>> {
    let mut table: FxHashMap<Partition, Vec<FragmentOffsetFrequency>> = FxHashMap::default();
    for &part in partitions {
        let max_frag_charge = (part.charge - 1).max(1);
        let mut frags: Vec<FragmentOffsetFrequency> = Vec::new();
        for fc in 1..=max_frag_charge {
            frags.push(FragmentOffsetFrequency {
                ion_type: IonType::Prefix { charge: fc, offset_bits: b_offset_bits(), loss_class: 0 },
                frequency: 0.0,
            });
            frags.push(FragmentOffsetFrequency {
                ion_type: IonType::Suffix { charge: fc, offset_bits: y_offset_bits(), loss_class: 0 },
                frequency: 0.0,
            });
        }
        frags.push(FragmentOffsetFrequency { ion_type: IonType::Noise, frequency: 0.0 });
        table.insert(part, frags);
    }
    table
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ion_charges(frags: &[FragmentOffsetFrequency], want_prefix: bool) -> Vec<i32> {
        let mut cs: Vec<i32> = frags
            .iter()
            .filter_map(|f| match f.ion_type {
                IonType::Prefix { charge, .. } if want_prefix => Some(charge),
                IonType::Suffix { charge, .. } if !want_prefix => Some(charge),
                _ => None,
            })
            .collect();
        cs.sort();
        cs
    }

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
    fn tiers_by_charge_groups_and_quantiles_per_charge() {
        // charge 2: 4 masses, 2 tiers -> [100,300]; charge 3: 2 masses -> [500,600].
        let pairs = vec![
            (2, 100.0), (2, 200.0), (2, 300.0), (2, 400.0),
            (3, 500.0), (3, 600.0),
        ];
        let got = derive_tiers_by_charge(&pairs, 2);
        assert_eq!(got, vec![(2, vec![100.0, 300.0]), (3, vec![500.0, 600.0])]);
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

    #[test]
    fn frag_off_table_ion_set_scales_with_precursor_charge() {
        let p2 = Partition { charge: 2, parent_mass: 500.0, seg_num: 0 };
        let p3 = Partition { charge: 3, parent_mass: 800.0, seg_num: 0 };
        let table = build_frag_off_table(&[p2, p3]);

        // charge 2 -> fragment charges {1}: b1, y1, + Noise = 3 entries.
        let f2 = table.get(&p2).expect("partition present");
        assert_eq!(ion_charges(f2, true), vec![1], "b-ion charges");
        assert_eq!(ion_charges(f2, false), vec![1], "y-ion charges");
        assert_eq!(f2.iter().filter(|f| f.ion_type.is_noise()).count(), 1, "one Noise");
        assert_eq!(f2.len(), 3);

        // charge 3 -> fragment charges {1,2}: b1,b2,y1,y2,+Noise = 5 entries.
        let f3 = table.get(&p3).expect("partition present");
        assert_eq!(ion_charges(f3, true), vec![1, 2]);
        assert_eq!(ion_charges(f3, false), vec![1, 2]);
        assert_eq!(f3.len(), 5);
    }
}
