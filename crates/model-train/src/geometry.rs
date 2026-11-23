use rustc_hash::FxHashMap;
use scoring_crate::params::{FragmentOffsetFrequency, IonType, Param, Partition};

use crate::labeled::LabeledMatch;

#[derive(Debug, Clone)]
pub struct GeometryConfig {
    pub num_segments: i32,
    pub max_rank: i32,
    pub mass_tier_occupancy: usize,
    pub max_mass_tiers: usize,
    pub max_fragment_charge: i32,
}

fn adaptive_n_tiers(n_psms: usize, occupancy: usize, max_tiers: usize) -> usize {
    let cap = max_tiers.max(1);
    if occupancy == 0 {
        return cap;
    }
    (n_psms / occupancy).clamp(1, cap)
}

fn b_offset_bits() -> u32 {
    (model::mass::PROTON as f32).to_bits()
}
fn y_offset_bits() -> u32 {
    ((model::mass::H2O + model::mass::PROTON) as f32).to_bits()
}
fn c_offset_bits() -> u32 {
    ((model::mass::PROTON + model::mass::NH3) as f32).to_bits()
}
fn z_offset_bits() -> u32 {
    ((model::mass::H2O + model::mass::PROTON + model::mass::Z_DOT_OFFSET) as f32).to_bits()
}

pub fn derive_mass_tiers(masses: &[f32], n_tiers: usize) -> Vec<f32> {
    let mut sorted: Vec<f32> = masses.iter().copied().filter(|m| m.is_finite()).collect();
    if sorted.is_empty() || n_tiers == 0 {
        return Vec::new();
    }
    sorted.sort_by(f32::total_cmp);
    let n = sorted.len();
    let mut tiers: Vec<f32> = Vec::with_capacity(n_tiers);
    for k in 0..n_tiers {
        let idx = ((k * n) / n_tiers).min(n - 1);
        let bound = sorted[idx];
        if tiers.last() != Some(&bound) {
            tiers.push(bound);
        }
    }
    tiers
}

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
            let n_tiers = adaptive_n_tiers(masses.len(), occupancy, max_tiers);
            (charge, derive_mass_tiers(&masses, n_tiers))
        })
        .collect()
}

pub fn build_partition_skeleton(
    tiers_by_charge: &[(i32, Vec<f32>)],
    num_segments: i32,
) -> Vec<Partition> {
    let mut parts: Vec<Partition> = Vec::new();
    for &(charge, ref tiers) in tiers_by_charge {
        for &parent_mass in tiers {
            for seg_num in 0..num_segments.max(1) {
                parts.push(Partition {
                    charge,
                    parent_mass,
                    seg_num,
                });
            }
        }
    }
    parts.sort();
    parts
}

pub fn build_frag_off_table(
    partitions: &[Partition],
    max_fragment_charge: i32,
    activation: model::activation::ActivationMethod,
) -> FxHashMap<Partition, Vec<FragmentOffsetFrequency>> {
    let mut table: FxHashMap<Partition, Vec<FragmentOffsetFrequency>> = FxHashMap::default();
    let cap = max_fragment_charge.max(1);
    use model::activation::ActivationMethod;
    let offsets: &[(u32, u32)] = match activation {
        ActivationMethod::ETD => &[(c_offset_bits(), z_offset_bits())],
        ActivationMethod::EThcD => &[
            (b_offset_bits(), y_offset_bits()),
            (c_offset_bits(), z_offset_bits()),
        ],
        _ => &[(b_offset_bits(), y_offset_bits())],
    };
    for &part in partitions {
        let max_frag_charge = (part.charge - 1).clamp(1, cap);
        let mut frags: Vec<FragmentOffsetFrequency> = Vec::new();
        for fc in 1..=max_frag_charge {
            for &(prefix_off, suffix_off) in offsets {
                frags.push(FragmentOffsetFrequency {
                    ion_type: IonType::Prefix {
                        charge: fc,
                        offset_bits: prefix_off,
                        loss_class: 0,
                    },
                    frequency: 0.0,
                });
                frags.push(FragmentOffsetFrequency {
                    ion_type: IonType::Suffix {
                        charge: fc,
                        offset_bits: suffix_off,
                        loss_class: 0,
                    },
                    frequency: 0.0,
                });
            }
        }
        frags.push(FragmentOffsetFrequency {
            ion_type: IonType::Noise,
            frequency: 0.0,
        });
        table.insert(part, frags);
    }
    table
}

pub fn derive_geometry(charge_masses: &[(i32, f32)], base: &Param, cfg: &GeometryConfig) -> Param {
    let tiers_by_charge =
        derive_tiers_by_charge(charge_masses, cfg.mass_tier_occupancy, cfg.max_mass_tiers);
    let num_segments = cfg.num_segments.max(1);
    let partitions = build_partition_skeleton(&tiers_by_charge, num_segments);
    let is_etd = base.data_type.activation == model::activation::ActivationMethod::ETD;
    let effective_frag_charge = if base.apply_deconvolution && !is_etd {
        1
    } else {
        cfg.max_fragment_charge
    };
    let frag_off_table = build_frag_off_table(
        &partitions,
        effective_frag_charge,
        base.data_type.activation,
    );

    let mut charge_counts: std::collections::BTreeMap<i32, i32> = std::collections::BTreeMap::new();
    for &(charge, _) in charge_masses {
        *charge_counts.entry(charge).or_insert(0) += 1;
    }
    let charge_hist: Vec<(i32, i32)> = charge_counts.iter().map(|(&c, &n)| (c, n)).collect();
    let min_charge = charge_counts
        .keys()
        .next()
        .copied()
        .unwrap_or(base.min_charge);
    let max_charge = charge_counts
        .keys()
        .next_back()
        .copied()
        .unwrap_or(base.max_charge);

    let mut param = Param {
        version: base.version,
        data_type: base.data_type.clone(),
        mme: base.mme,
        apply_deconvolution: base.apply_deconvolution,
        deconvolution_error_tolerance: base.deconvolution_error_tolerance,
        num_precursor_off: base.num_precursor_off,
        precursor_off_map: base.precursor_off_map.clone(),
        error_scaling_factor: base.error_scaling_factor,
        num_segments,
        max_rank: cfg.max_rank,
        partitions,
        frag_off_table,
        charge_hist,
        min_charge,
        max_charge,
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
    use scoring_crate::params::SpecDataType;

    fn base_param() -> Param {
        let part = Partition {
            charge: 2,
            parent_mass: 1500.0,
            seg_num: 0,
        };
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
            let residues = seq
                .iter()
                .map(|&r| AminoAcid::standard(r).unwrap())
                .collect();
            Peptide::new(residues, b'_', b'-')
        }
        let p1 = pep(b"PEPTIDE");
        let p2 = pep(b"PEPTIDER");
        let (m1, m2) = (p1.mass() as f32, p2.mass() as f32);
        let labels = vec![
            LabeledMatch {
                spectrum_index: 0,
                peptide: p1,
                charge: 2,
                confidence: 0.001,
            },
            LabeledMatch {
                spectrum_index: 1,
                peptide: p2,
                charge: 3,
                confidence: 0.001,
            },
        ];
        assert_eq!(corpus_charge_masses(&labels), vec![(2, m1), (3, m2)]);
    }

    #[test]
    fn derive_geometry_deconv_forces_fragment_charge_one() {
        let mut base = base_param();
        base.apply_deconvolution = true;
        let charge_masses = vec![(2, 1000.0), (4, 2500.0)];
        let cfg = GeometryConfig {
            num_segments: 1,
            max_rank: 150,
            mass_tier_occupancy: 1,
            max_mass_tiers: 1,
            max_fragment_charge: 3,
        };
        let p = derive_geometry(&charge_masses, &base, &cfg);
        for frags in p.frag_off_table.values() {
            assert!(
                ion_charges(frags, true).iter().all(|&c| c == 1),
                "deconv -> b charge 1 only"
            );
            assert!(
                ion_charges(frags, false).iter().all(|&c| c == 1),
                "deconv -> y charge 1 only"
            );
        }
    }

    #[test]
    fn derive_geometry_ethcd_models_by_and_cz_at_charge_one() {
        let mut base = base_param();
        base.apply_deconvolution = true;
        base.data_type.activation = ActivationMethod::EThcD;
        let charge_masses = vec![(2, 1000.0), (4, 2500.0)];
        let cfg = GeometryConfig {
            num_segments: 1,
            max_rank: 150,
            mass_tier_occupancy: 1,
            max_mass_tiers: 1,
            max_fragment_charge: 3,
        };
        let p = derive_geometry(&charge_masses, &base, &cfg);
        let want = [
            (true, b_offset_bits()),
            (false, y_offset_bits()),
            (true, c_offset_bits()),
            (false, z_offset_bits()),
        ];
        for frags in p.frag_off_table.values() {
            let mut got: Vec<(bool, u32)> = frags
                .iter()
                .filter_map(|f| match f.ion_type {
                    IonType::Prefix {
                        charge,
                        offset_bits,
                        ..
                    } => Some((true, offset_bits, charge)),
                    IonType::Suffix {
                        charge,
                        offset_bits,
                        ..
                    } => Some((false, offset_bits, charge)),
                    IonType::Noise => None,
                })
                .map(|(prefix, off, charge)| {
                    assert_eq!(charge, 1, "deconvolved EThcD models charge-1 ions");
                    (prefix, off)
                })
                .collect();
            got.sort_unstable();
            let mut want = want.to_vec();
            want.sort_unstable();
            assert_eq!(got, want);
            assert!(frags.iter().any(|f| f.ion_type.is_noise()));
        }
    }

    #[test]
    fn derive_geometry_etd_models_multicharge_cz() {
        let mut base = base_param();
        base.apply_deconvolution = true;
        base.data_type.activation = ActivationMethod::ETD;
        let charge_masses = vec![(2, 1000.0), (5, 3000.0)];
        let cfg = GeometryConfig {
            num_segments: 1,
            max_rank: 150,
            mass_tier_occupancy: 1,
            max_mass_tiers: 1,
            max_fragment_charge: 5,
        };
        let p = derive_geometry(&charge_masses, &base, &cfg);
        let c_off = ((model::mass::PROTON + model::mass::NH3) as f32).to_bits();
        let z_off =
            ((model::mass::H2O + model::mass::PROTON + model::mass::Z_DOT_OFFSET) as f32).to_bits();
        let mut saw_multicharge = false;
        for (part, frags) in &p.frag_off_table {
            for f in frags {
                match f.ion_type {
                    IonType::Prefix {
                        charge,
                        offset_bits,
                        ..
                    } => {
                        assert_eq!(offset_bits, c_off, "ETD prefix must be c-ion");
                        if charge > 1 {
                            saw_multicharge = true;
                        }
                    }
                    IonType::Suffix { offset_bits, .. } => {
                        assert_eq!(offset_bits, z_off, "ETD suffix must be z-ion");
                    }
                    IonType::Noise => {}
                }
            }
            if part.charge == 5 {
                assert!(
                    ion_charges(frags, true).contains(&4),
                    "z5 ETD partition models charge-4 c ions"
                );
            }
        }
        assert!(saw_multicharge, "ETD must model multi-charge c/z ions");
    }

    #[test]
    fn derive_geometry_low_res_keeps_multicharge_fragments() {
        let mut base = base_param();
        base.apply_deconvolution = false;
        let charge_masses = vec![(3, 2000.0)];
        let cfg = GeometryConfig {
            num_segments: 1,
            max_rank: 150,
            mass_tier_occupancy: 1,
            max_mass_tiers: 1,
            max_fragment_charge: 3,
        };
        let p = derive_geometry(&charge_masses, &base, &cfg);
        let frags = p.frag_off_table.values().next().expect("a partition");
        assert_eq!(
            ion_charges(frags, true),
            vec![1, 2],
            "no deconv, charge-3 -> b at 1,2"
        );
    }

    #[test]
    fn derive_geometry_assembles_geometry_only_param() {
        let base = base_param();
        let charge_masses = vec![
            (2, 1000.0),
            (2, 1200.0),
            (2, 1400.0),
            (2, 1600.0),
            (3, 2000.0),
            (3, 2400.0),
        ];
        let cfg = GeometryConfig {
            num_segments: 2,
            max_rank: 150,
            mass_tier_occupancy: 1,
            max_mass_tiers: 2,
            max_fragment_charge: 2,
        };
        let p = derive_geometry(&charge_masses, &base, &cfg);

        assert_eq!(p.num_segments, 2);
        assert_eq!(p.max_rank, 150);
        assert_eq!(p.partitions.len(), 8);
        assert_eq!(p.min_charge, 2);
        assert_eq!(p.max_charge, 3);
        for part in &p.partitions {
            let frags = p
                .frag_off_table
                .get(part)
                .expect("frag entry per partition");
            assert!(frags.iter().any(|f| f.ion_type.is_noise()));
        }
        assert_eq!(p.mme, base.mme);
        assert_eq!(p.version, base.version);
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
        let masses = vec![100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0];
        let tiers = derive_mass_tiers(&masses, 4);
        assert_eq!(tiers, vec![100.0, 300.0, 500.0, 700.0]);
    }

    #[test]
    fn non_finite_parent_masses_are_ignored() {
        let masses = vec![300.0, f32::NAN, 100.0, f32::INFINITY, 200.0, 400.0];
        assert_eq!(derive_mass_tiers(&masses, 2), vec![100.0, 300.0]);
        assert!(derive_mass_tiers(&[f32::NAN], 2).is_empty());
    }

    #[test]
    fn tiers_by_charge_groups_and_quantiles_per_charge() {
        let pairs = vec![
            (2, 100.0),
            (2, 200.0),
            (2, 300.0),
            (2, 400.0),
            (3, 500.0),
            (3, 600.0),
        ];
        let got = derive_tiers_by_charge(&pairs, 1, 2);
        assert_eq!(got, vec![(2, vec![100.0, 300.0]), (3, vec![500.0, 600.0])]);
    }

    #[test]
    fn adaptive_tiers_scale_with_per_charge_data() {
        assert_eq!(adaptive_n_tiers(0, 100, 33), 1);
        assert_eq!(adaptive_n_tiers(50, 100, 33), 1);
        assert_eq!(adaptive_n_tiers(100, 100, 33), 1);
        assert_eq!(adaptive_n_tiers(350, 100, 33), 3);
        assert_eq!(adaptive_n_tiers(100_000, 2500, 33), 33);
        assert_eq!(adaptive_n_tiers(10, 0, 33), 33);
    }

    #[test]
    fn tiers_by_charge_data_rich_gets_more_tiers_than_sparse() {
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
        assert_eq!(t2.len(), 4);
        assert_eq!(t4.len(), 1);
        assert!(t2.len() > t4.len());
    }

    #[test]
    fn skeleton_is_cartesian_product_sorted_lex() {
        let tiers_by_charge = vec![(2, vec![100.0, 500.0]), (3, vec![200.0, 600.0])];
        let parts = build_partition_skeleton(&tiers_by_charge, 2);
        let expect = vec![
            Partition {
                charge: 2,
                parent_mass: 100.0,
                seg_num: 0,
            },
            Partition {
                charge: 2,
                parent_mass: 500.0,
                seg_num: 0,
            },
            Partition {
                charge: 2,
                parent_mass: 100.0,
                seg_num: 1,
            },
            Partition {
                charge: 2,
                parent_mass: 500.0,
                seg_num: 1,
            },
            Partition {
                charge: 3,
                parent_mass: 200.0,
                seg_num: 0,
            },
            Partition {
                charge: 3,
                parent_mass: 600.0,
                seg_num: 0,
            },
            Partition {
                charge: 3,
                parent_mass: 200.0,
                seg_num: 1,
            },
            Partition {
                charge: 3,
                parent_mass: 600.0,
                seg_num: 1,
            },
        ];
        assert_eq!(parts, expect);
    }

    #[test]
    fn frag_off_table_caps_fragment_charge() {
        let p2 = Partition {
            charge: 2,
            parent_mass: 500.0,
            seg_num: 0,
        };
        let p3 = Partition {
            charge: 3,
            parent_mass: 800.0,
            seg_num: 0,
        };

        let t1 = build_frag_off_table(&[p2, p3], 1, ActivationMethod::HCD);
        assert_eq!(ion_charges(t1.get(&p2).unwrap(), true), vec![1]);
        assert_eq!(
            ion_charges(t1.get(&p3).unwrap(), true),
            vec![1],
            "charge-3 capped to frag 1"
        );
        assert_eq!(t1.get(&p3).unwrap().len(), 3, "b1,y1,Noise");

        let t2 = build_frag_off_table(&[p2, p3], 2, ActivationMethod::HCD);
        assert_eq!(
            ion_charges(t2.get(&p2).unwrap(), true),
            vec![1],
            "charge-2 bounded by C-1"
        );
        assert_eq!(ion_charges(t2.get(&p3).unwrap(), true), vec![1, 2]);
        assert_eq!(t2.get(&p3).unwrap().len(), 5, "b1,b2,y1,y2,Noise");
    }
}
