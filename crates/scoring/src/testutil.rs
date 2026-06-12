//! Test fixtures shared across engine module tests.
//!
//! `cfg(test)` only — does not appear in release builds.

use rustc_hash::FxHashMap;

use model::activation::ActivationMethod;
use model::instrument::InstrumentType;
use crate::param_model::{FragmentOffsetFrequency, IonType, Param, Partition, SpecDataType};
use model::protocol::Protocol;
use model::tolerance::Tolerance;

/// Minimal `Param` for testing: 1 partition (charge=2, parent_mass=1500.0,
/// seg_num=0), 1 prefix ion (charge=1, offset=0) + Noise, max_rank=3, empty
/// frag_off_table, Ppm(20.0) tolerance.
///
/// This is the canonical fixture from `scoring/rank_scorer.rs:185`, promoted
/// to a shared helper so every duplicate site can import it instead of
/// rebuilding 50 lines of boilerplate.
pub fn tiny_param() -> Param {
    let part = Partition { charge: 2, parent_mass: 1500.0, seg_num: 0 };
    let prefix_ion = IonType::Prefix { charge: 1, offset_bits: 0.0_f32.to_bits() };
    let noise_ion = IonType::Noise;

    // max_rank = 3 means each rank-distribution array has length 4
    // (indices 0..2 for ranks 1..3, index 3 for "missing ion" slot).
    let max_rank = 3;
    // ion_freqs[i] / noise_freqs[i] computed manually:
    //   index 0: 0.6 / 0.1 = 6.0
    //   index 1: 0.3 / 0.2 = 1.5
    //   index 2: 0.05 / 0.3 = 0.166...
    //   index 3 (missing): 0.001 / 0.4 = 0.0025
    let ion_freqs = vec![0.6_f32, 0.3, 0.05, 0.001];
    let noise_freqs = vec![0.1_f32, 0.2, 0.3, 0.4];

    let mut ion_table_inner: FxHashMap<IonType, Vec<f32>> = FxHashMap::default();
    ion_table_inner.insert(prefix_ion, ion_freqs);
    ion_table_inner.insert(noise_ion, noise_freqs);

    let mut rank_dist_table: FxHashMap<Partition, FxHashMap<IonType, Vec<f32>>> = FxHashMap::default();
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
    };
    p.rebuild_cache();
    p
}

/// Richer `Param` for testing the GF / ScoredSpectrum scoring paths.
///
/// Differs from `tiny_param()` in three ways that matter for the GF tests:
/// - `parent_mass = 1000.0` (smaller, so GF DP exercises fewer nodes)
/// - `mme = Tolerance::Da(0.5)` (simpler tolerance arithmetic in fragment lookup)
/// - `frag_off_table` seeded with one `FragmentOffsetFrequency` entry for the
///   prefix ion, so `ion_types_for_segment(0)` returns a non-empty list and
///   `node_score` / `edge_score` can exercise the live scoring paths.
///
/// Used by tests in `scoring/scored_spectrum.rs`, `gf/group.rs`, and
/// `gf/primitive_graph.rs`.
pub fn tiny_param_with_ions() -> Param {
    let part = Partition { charge: 2, parent_mass: 1000.0, seg_num: 0 };
    let prefix1 = IonType::Prefix { charge: 1, offset_bits: 0.0_f32.to_bits() };
    let noise = IonType::Noise;

    // max_rank=3 → 4 slots. Ion has higher freq at rank 1.
    let ion_freqs = vec![0.6_f32, 0.3, 0.05, 0.001];
    let noise_freqs = vec![0.1_f32, 0.2, 0.3, 0.4];

    let mut ion_table: FxHashMap<IonType, Vec<f32>> = FxHashMap::default();
    ion_table.insert(prefix1, ion_freqs);
    ion_table.insert(noise, noise_freqs);

    let mut rank_dist_table: FxHashMap<Partition, FxHashMap<IonType, Vec<f32>>> = FxHashMap::default();
    rank_dist_table.insert(part, ion_table);

    // frag_off_table: one prefix ion entry so ion_types_for_segment returns it.
    let mut frag_off_table = FxHashMap::default();
    frag_off_table.insert(part, vec![FragmentOffsetFrequency {
        ion_type: prefix1,
        frequency: 0.7,
    }]);

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
    };
    p.rebuild_cache();
    p
}
