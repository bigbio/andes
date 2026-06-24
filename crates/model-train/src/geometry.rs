//! Data-derived partition/segment geometry.
//!
//! Builds the scoring [`Param`]'s structural geometry (charge span, per-charge
//! `parent_mass` tier boundaries, segment count, rank cap, and per-partition ion
//! membership) from andes's own labeled corpus, so a trained model owns its
//! geometry rather than inheriting it from a seed template.

use rustc_hash::FxHashMap;
use scoring_crate::param_model::{FragmentOffsetFrequency, IonType, Param, Partition};

use crate::labeled::LabeledMatch;

/// Knobs for [`derive_geometry`] — the structural choices that, in a seed model,
/// were inherited verbatim. Swept by the benchmark harness; collapse to fixed
/// defaults once the optimum is validated.
#[derive(Debug, Clone)]
pub struct GeometryConfig {
    pub num_segments: i32,
    pub max_rank: i32,
    /// Target training PSMs per mass tier. The tier count is derived PER CHARGE
    /// as `n_psms_for_charge / mass_tier_occupancy` (clamped to
    /// `1..=max_mass_tiers`), so a data-rich charge gets many tiers and a sparse
    /// charge few — matching how a seed geometry's resolution scales with data
    /// (the dominant charges carry ~33 tiers, sparse charges ~4). A fixed small
    /// tier count instead under-partitions the data-rich charges, leaving the
    /// learned tables badly undertrained.
    pub mass_tier_occupancy: usize,
    /// Upper bound on mass tiers per charge.
    pub max_mass_tiers: usize,
}

/// Adaptive per-charge mass-tier count: `n_psms / occupancy`, clamped to
/// `1..=max_tiers`. An `occupancy` of 0 falls back to `max_tiers` (treat every
/// charge as data-rich). Guarantees at least one tier so a charge with any data
/// is never dropped.
fn adaptive_n_tiers(n_psms: usize, occupancy: usize, max_tiers: usize) -> usize {
    let cap = max_tiers.max(1);
    if occupancy == 0 {
        return cap;
    }
    (n_psms / occupancy).clamp(1, cap)
}

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
    occupancy: usize,
    max_tiers: usize,
) -> Vec<(i32, Vec<f32>)> {
    let mut by_charge: std::collections::BTreeMap<i32, Vec<f32>> =
        std::collections::BTreeMap::new();
    for &(charge, mass) in charge_masses {
        by_charge.entry(charge).or_default().push(mass);
    }
    by_charge
        .into_iter()
        .map(|(charge, masses)| {
            // Adaptive: tier count scales with this charge's data volume.
            let n_tiers = adaptive_n_tiers(masses.len(), occupancy, max_tiers);
            (charge, derive_mass_tiers(&masses, n_tiers))
        })
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

/// Assemble a geometry-only [`Param`] from corpus `(charge, parent_mass)` PSMs:
/// derive per-charge tiers → partition skeleton → chemistry ion table, take the
/// segment count / rank cap from `cfg`, and clone the **non-geometry** metadata
/// (data_type, tolerance, deconvolution, version, precursor offsets) from `base`.
/// Learned tables are left empty for [`Estimator::estimate`] to fill.
pub fn derive_geometry(charge_masses: &[(i32, f32)], base: &Param, cfg: &GeometryConfig) -> Param {
    let tiers_by_charge =
        derive_tiers_by_charge(charge_masses, cfg.mass_tier_occupancy, cfg.max_mass_tiers);
    let partitions = build_partition_skeleton(&tiers_by_charge, cfg.num_segments);
    let frag_off_table = build_frag_off_table(&partitions);

    // Charge histogram + span from the corpus.
    let mut charge_counts: std::collections::BTreeMap<i32, i32> = std::collections::BTreeMap::new();
    for &(charge, _) in charge_masses {
        *charge_counts.entry(charge).or_insert(0) += 1;
    }
    let charge_hist: Vec<(i32, i32)> = charge_counts.iter().map(|(&c, &n)| (c, n)).collect();
    let min_charge = charge_counts.keys().next().copied().unwrap_or(base.min_charge);
    let max_charge = charge_counts.keys().next_back().copied().unwrap_or(base.max_charge);

    let mut param = Param {
        // Non-geometry metadata cloned from the base/seed.
        version: base.version,
        data_type: base.data_type.clone(),
        mme: base.mme,
        apply_deconvolution: base.apply_deconvolution,
        deconvolution_error_tolerance: base.deconvolution_error_tolerance,
        num_precursor_off: base.num_precursor_off,
        precursor_off_map: base.precursor_off_map.clone(),
        error_scaling_factor: base.error_scaling_factor,
        // Geometry derived from config + corpus.
        num_segments: cfg.num_segments,
        max_rank: cfg.max_rank,
        partitions,
        frag_off_table,
        charge_hist,
        min_charge,
        max_charge,
        // Learned tables left empty for `Estimator::estimate` to fill.
        rank_dist_table: FxHashMap::default(),
        ion_err_dist_table: FxHashMap::default(),
        noise_err_dist_table: FxHashMap::default(),
        ion_existence_table: FxHashMap::default(),
        partition_ion_types_cache: FxHashMap::default(),
        gbdt_peak_model: None,
        frag_intensity_model: None,
        rich_ion_model: None,
    };
    param.rebuild_cache();
    param
}

/// Extract the `(charge, parent_mass)` corpus from confident labels — the input
/// to [`derive_geometry`]. `parent_mass = peptide.mass()` (the neutral mass the
/// partitioner keys on).
pub fn corpus_charge_masses(labels: &[LabeledMatch]) -> Vec<(i32, f32)> {
    labels
        .iter()
        .map(|l| (l.charge as i32, l.peptide.mass() as f32))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::activation::ActivationMethod;
    use model::instrument::InstrumentType;
    use model::protocol::Protocol;
    use model::tolerance::Tolerance;
    use scoring_crate::param_model::SpecDataType;

    /// Minimal non-geometry base `Param` (the metadata `derive_geometry` clones).
    fn base_param() -> Param {
        let part = Partition { charge: 2, parent_mass: 1500.0, seg_num: 0 };
        let mut frag_off_table: FxHashMap<Partition, Vec<FragmentOffsetFrequency>> =
            FxHashMap::default();
        frag_off_table.insert(part, vec![]);
        let mut p = Param {
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
            rank_dist_table: FxHashMap::default(),
            error_scaling_factor: 0,
            ion_err_dist_table: FxHashMap::default(),
            noise_err_dist_table: FxHashMap::default(),
            ion_existence_table: FxHashMap::default(),
            partition_ion_types_cache: FxHashMap::default(),
            gbdt_peak_model: None,
            frag_intensity_model: None,
            rich_ion_model: None,
        };
        p.rebuild_cache();
        p
    }

    #[test]
    fn corpus_extraction_maps_charge_and_peptide_mass() {
        use model::amino_acid::AminoAcid;
        use model::peptide::Peptide;
        fn pep(seq: &[u8]) -> Peptide {
            let residues = seq.iter().map(|&r| AminoAcid::standard(r).unwrap()).collect();
            Peptide::new(residues, b'_', b'-')
        }
        let p1 = pep(b"PEPTIDE");
        let p2 = pep(b"PEPTIDER");
        let (m1, m2) = (p1.mass() as f32, p2.mass() as f32);
        let labels = vec![
            LabeledMatch { spectrum_index: 0, peptide: p1, charge: 2, confidence: 0.001 },
            LabeledMatch { spectrum_index: 1, peptide: p2, charge: 3, confidence: 0.001 },
        ];
        assert_eq!(corpus_charge_masses(&labels), vec![(2, m1), (3, m2)]);
    }

    #[test]
    fn derive_geometry_assembles_geometry_only_param() {
        let base = base_param();
        let charge_masses = vec![
            (2, 1000.0), (2, 1200.0), (2, 1400.0), (2, 1600.0),
            (3, 2000.0), (3, 2400.0),
        ];
        let cfg = GeometryConfig {
            num_segments: 2,
            max_rank: 150,
            mass_tier_occupancy: 1,
            max_mass_tiers: 2,
        };
        let p = derive_geometry(&charge_masses, &base, &cfg);

        // Geometry comes from config + data.
        assert_eq!(p.num_segments, 2);
        assert_eq!(p.max_rank, 150);
        // 2 charges × 2 tiers × 2 segments = 8 partitions.
        assert_eq!(p.partitions.len(), 8);
        assert_eq!(p.min_charge, 2);
        assert_eq!(p.max_charge, 3);
        // Every partition carries a Noise entry (RankScorer requires it).
        for part in &p.partitions {
            let frags = p.frag_off_table.get(part).expect("frag entry per partition");
            assert!(frags.iter().any(|f| f.ion_type.is_noise()));
        }
        // Non-geometry metadata cloned from base.
        assert_eq!(p.mme, base.mme);
        assert_eq!(p.version, base.version);
        // Learned tables start empty (estimate fills them).
        assert!(p.rank_dist_table.is_empty());
    }

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
        // occupancy=1, max=2 → each charge gets min(n_masses, 2) tiers = 2.
        let got = derive_tiers_by_charge(&pairs, 1, 2);
        assert_eq!(got, vec![(2, vec![100.0, 300.0]), (3, vec![500.0, 600.0])]);
    }

    #[test]
    fn adaptive_tiers_scale_with_per_charge_data() {
        assert_eq!(adaptive_n_tiers(0, 100, 33), 1); // never drop a charge with data
        assert_eq!(adaptive_n_tiers(50, 100, 33), 1); // below one occupancy → 1
        assert_eq!(adaptive_n_tiers(100, 100, 33), 1);
        assert_eq!(adaptive_n_tiers(350, 100, 33), 3);
        assert_eq!(adaptive_n_tiers(100_000, 2500, 33), 33); // data-rich → capped
        assert_eq!(adaptive_n_tiers(10, 0, 33), 33); // occupancy 0 → cap
    }

    #[test]
    fn tiers_by_charge_data_rich_gets_more_tiers_than_sparse() {
        // charge 2: 40 masses, charge 4: 4 masses; occupancy=10, cap=8.
        let mut pairs: Vec<(i32, f32)> = Vec::new();
        for i in 0..40 {
            pairs.push((2, 500.0 + i as f32));
        }
        for i in 0..4 {
            pairs.push((4, 800.0 + i as f32));
        }
        let got = derive_tiers_by_charge(&pairs, 10, 8);
        let t2 = &got.iter().find(|(c, _)| *c == 2).unwrap().1;
        let t4 = &got.iter().find(|(c, _)| *c == 4).unwrap().1;
        assert_eq!(t2.len(), 4); // 40/10 = 4 tiers
        assert_eq!(t4.len(), 1); // 4/10 < 1 → clamped to 1
        assert!(t2.len() > t4.len());
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
