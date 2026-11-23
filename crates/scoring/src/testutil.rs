use rustc_hash::FxHashMap;

use crate::params::{FragmentOffsetFrequency, IonType, Param, Partition, SpecDataType};
use model::activation::ActivationMethod;
use model::instrument::InstrumentType;
use model::protocol::Protocol;
use model::tolerance::Tolerance;

pub fn tiny_param() -> Param {
    let part = Partition {
        charge: 2,
        parent_mass: 1500.0,
        seg_num: 0,
    };
    let prefix_ion = IonType::Prefix {
        charge: 1,
        offset_bits: 0.0_f32.to_bits(),
        loss_class: 0,
    };
    let noise_ion = IonType::Noise;

    let max_rank = 3;
    let ion_freqs = vec![0.6_f32, 0.3, 0.05, 0.001];
    let noise_freqs = vec![0.1_f32, 0.2, 0.3, 0.4];

    let mut ion_table_inner: FxHashMap<IonType, Vec<f32>> = FxHashMap::default();
    ion_table_inner.insert(prefix_ion, ion_freqs);
    ion_table_inner.insert(noise_ion, noise_freqs);

    let mut rank_dist_table: FxHashMap<Partition, FxHashMap<IonType, Vec<f32>>> =
        FxHashMap::default();
    rank_dist_table.insert(part, ion_table_inner);

    let mut frag_off_table = FxHashMap::default();
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
        max_rank,
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
    p.rebuild_cache();
    p
}

pub fn tiny_param_with_ions() -> Param {
    let part = Partition {
        charge: 2,
        parent_mass: 1000.0,
        seg_num: 0,
    };
    let prefix1 = IonType::Prefix {
        charge: 1,
        offset_bits: 0.0_f32.to_bits(),
        loss_class: 0,
    };
    let noise = IonType::Noise;

    let ion_freqs = vec![0.6_f32, 0.3, 0.05, 0.001];
    let noise_freqs = vec![0.1_f32, 0.2, 0.3, 0.4];

    let mut ion_table: FxHashMap<IonType, Vec<f32>> = FxHashMap::default();
    ion_table.insert(prefix1, ion_freqs);
    ion_table.insert(noise, noise_freqs);

    let mut rank_dist_table: FxHashMap<Partition, FxHashMap<IonType, Vec<f32>>> =
        FxHashMap::default();
    rank_dist_table.insert(part, ion_table);

    let mut frag_off_table = FxHashMap::default();
    frag_off_table.insert(
        part,
        vec![FragmentOffsetFrequency {
            ion_type: prefix1,
            frequency: 0.7,
        }],
    );

    let mut p = Param {
        version: 10001,
        data_type: SpecDataType {
            activation: ActivationMethod::HCD,
            instrument: InstrumentType::QExactive,
            enzyme: None,
            protocol: Protocol::Automatic,
        },
        mme: Tolerance::Da(0.5),
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
    p.rebuild_cache();
    p
}
