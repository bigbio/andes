//! Verify pooled and non-pooled PrimitiveAaGraph construction produce
//! bit-identical output for the same inputs across multiple fixtures.
//!
//! Task 1 of `docs/superpowers/plans/2026-05-11-suffix-array-refactor-plan.md`:
//! thread-local arena pool for `PrimitiveAaGraph::new`'s 11 per-call Vec
//! allocations. Bit-identical output required.

use rustc_hash::FxHashMap;

use model::{AminoAcidSetBuilder, Spectrum, Tolerance};
use model::activation::ActivationMethod;
use model::instrument::InstrumentType;
use model::protocol::Protocol;
use scoring::gf::PrimitiveAaGraph;
use scoring::param_model::{FragmentOffsetFrequency, IonType, Partition, SpecDataType};
use scoring::{Param, RankScorer, ScoredSpectrum};

/// Local mirror of `tiny_param_with_ions`. testutil is `pub(crate) cfg(test)`
/// so integration tests can't import it directly. Matches the fixture used in
/// `gf_graph_dp.rs`.
fn tiny_param() -> Param {
    let part = Partition { charge: 2, parent_mass: 1000.0, seg_num: 0 };
    let prefix1 = IonType::Prefix { charge: 1, offset_bits: 0.0_f32.to_bits() };
    let noise = IonType::Noise;

    let mut ion_table: FxHashMap<IonType, Vec<f32>> = FxHashMap::default();
    ion_table.insert(prefix1, vec![0.6_f32, 0.3, 0.05, 0.001]);
    ion_table.insert(noise, vec![0.1_f32, 0.2, 0.3, 0.4]);

    let mut rank_dist_table: FxHashMap<Partition, FxHashMap<IonType, Vec<f32>>> = FxHashMap::default();
    rank_dist_table.insert(part, ion_table);

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

fn empty_spec() -> Spectrum {
    Spectrum {
        title: "parity_test".into(),
        precursor_mz: 500.0,
        precursor_intensity: None,
        precursor_charge: Some(2),
        rt_seconds: None,
        scan: None,
        peaks: vec![],
        activation_method: None,
        isolation_lower_offset: None,
        isolation_upper_offset: None,
    }
}

/// Assert all observable fields of two `PrimitiveAaGraph` are bit-identical.
///
/// Fields compared:
/// - Scalars: `peptide_mass`, `direction`, `enzyme`, `min_node_mass`,
///   `mass_offset`, `node_count`, `source_node_idx`, `sink_node_idx`.
/// - Vectors: `active_nodes`, `mass_to_node_idx`, `edge_offset`,
///   `edge_prev_node`, `edge_prob` (compared as raw bit-patterns via
///   `f32::to_bits`), `edge_score`, `node_scores`.
fn assert_graphs_equal(a: &PrimitiveAaGraph, b: &PrimitiveAaGraph, label: &str) {
    assert_eq!(a.peptide_mass, b.peptide_mass, "{label}: peptide_mass");
    assert_eq!(a.direction, b.direction, "{label}: direction");
    assert_eq!(a.enzyme, b.enzyme, "{label}: enzyme");
    assert_eq!(a.min_node_mass, b.min_node_mass, "{label}: min_node_mass");
    assert_eq!(a.mass_offset, b.mass_offset, "{label}: mass_offset");
    assert_eq!(a.node_count, b.node_count, "{label}: node_count");
    assert_eq!(a.source_node_idx, b.source_node_idx, "{label}: source_node_idx");
    assert_eq!(a.sink_node_idx, b.sink_node_idx, "{label}: sink_node_idx");

    assert_eq!(a.active_nodes, b.active_nodes, "{label}: active_nodes");
    assert_eq!(a.mass_to_node_idx, b.mass_to_node_idx, "{label}: mass_to_node_idx");
    assert_eq!(a.edge_offset, b.edge_offset, "{label}: edge_offset");
    assert_eq!(a.edge_prev_node, b.edge_prev_node, "{label}: edge_prev_node");
    assert_eq!(a.edge_score, b.edge_score, "{label}: edge_score");
    assert_eq!(a.node_scores, b.node_scores, "{label}: node_scores");

    // Compare f32 vectors bit-for-bit (NaN-safe and detects any rounding drift).
    assert_eq!(a.edge_prob.len(), b.edge_prob.len(), "{label}: edge_prob len");
    for (i, (x, y)) in a.edge_prob.iter().zip(b.edge_prob.iter()).enumerate() {
        assert_eq!(
            x.to_bits(), y.to_bits(),
            "{label}: edge_prob[{i}] bit-mismatch ({x} vs {y})"
        );
    }
}

#[test]
fn pooled_graph_matches_unpooled_bit_for_bit() {
    let aa = AminoAcidSetBuilder::new_standard().build().unwrap();
    let param = tiny_param();
    let scorer = RankScorer::new(&param);
    let spec = empty_spec();
    let ss = ScoredSpectrum::new_without_filtering(&spec);

    // Six peptide masses spanning typical PXD001819 mass range.
    for &peptide_mass in &[500_i32, 800, 1200, 1800, 2400, 3000] {
        let g_unpooled = PrimitiveAaGraph::new(
            &aa, peptide_mass, None, &ss, &scorer, 2, 1000.0, 0.5, false, false,
        );
        let g_pooled = PrimitiveAaGraph::new_pooled(
            &aa, peptide_mass, None, &ss, &scorer, 2, 1000.0, 0.5, false, false,
        );
        assert_graphs_equal(&g_unpooled, &g_pooled, &format!("pep_mass={peptide_mass}"));
    }
}

#[test]
fn pooled_graph_repeated_calls_remain_correct() {
    // Calling new_pooled many times must continue to produce the same result
    // as new (catches stale-state bugs in the arena).
    let aa = AminoAcidSetBuilder::new_standard().build().unwrap();
    let param = tiny_param();
    let scorer = RankScorer::new(&param);
    let spec = empty_spec();
    let ss = ScoredSpectrum::new_without_filtering(&spec);

    let masses = [600_i32, 1500, 2200, 700, 1900, 1100];
    for &peptide_mass in &masses {
        let g_unpooled = PrimitiveAaGraph::new(
            &aa, peptide_mass, None, &ss, &scorer, 2, 1000.0, 0.5, false, false,
        );
        let g_pooled = PrimitiveAaGraph::new_pooled(
            &aa, peptide_mass, None, &ss, &scorer, 2, 1000.0, 0.5, false, false,
        );
        assert_graphs_equal(&g_unpooled, &g_pooled, &format!("repeat pep_mass={peptide_mass}"));
    }

    // And once more in reverse order for good measure.
    for &peptide_mass in masses.iter().rev() {
        let g_unpooled = PrimitiveAaGraph::new(
            &aa, peptide_mass, None, &ss, &scorer, 2, 1000.0, 0.5, false, false,
        );
        let g_pooled = PrimitiveAaGraph::new_pooled(
            &aa, peptide_mass, None, &ss, &scorer, 2, 1000.0, 0.5, false, false,
        );
        assert_graphs_equal(&g_unpooled, &g_pooled, &format!("reverse pep_mass={peptide_mass}"));
    }
}
