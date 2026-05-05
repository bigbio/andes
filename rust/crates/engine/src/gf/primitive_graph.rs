//! Primitive-array–based amino acid graph for the generating function.
//!
//! Ports Java `edu.ucsd.msjava.msgf.PrimitiveAminoAcidGraph` (~290 LOC).
//! Replaces FlexAminoAcidGraph in the DB search hot path to eliminate
//! HashMap/ArrayList/NominalMass object overhead.
//!
//! Graph topology is stored in CSR (Compressed Sparse Row) format:
//!   `edge_offset[node+1] - edge_offset[node]` = number of incoming edges for node
//!   `edge_prev_node[e]`, `edge_prob[e]`, `edge_score[e]` = edge data
//!
//! Node scores are stored in a flat `Vec<i32>` indexed by node index.
//!
//! # Construction phases (Java-faithful)
//!
//! 1. Resolve source/sink AA lists from `direction` and protein-term flags.
//! 2. Compute `min_node_mass` and `mass_offset` from minimum nominal masses.
//! 3. Reachability sweep + per-mass incoming-edge counts.
//! 4. Build `active_nodes` and `mass_to_node_idx` dense lookup.
//! 5. Build CSR `edge_offset` and fill `edge_prev_node`, `edge_prob`, `edge_score`.
//! 6. Compute edge error scores via `scored_spec.edge_score`.
//! 7. Compute node scores via `scored_spec.node_score`.

use crate::aa_set::AminoAcidSet;
use crate::amino_acid::AminoAcid;
use crate::enzyme::Enzyme;
use crate::modification::ModLocation;
use crate::scoring::rank_scorer::RankScorer;
use crate::scoring::scored_spectrum::ScoredSpectrum;

/// Primitive CSR amino-acid graph used by the generating-function DP.
///
/// All fields are `pub` so that Task 6's GF DP can read them without
/// accessor overhead. The graph is built once per (spectrum, peptide-mass)
/// pair and is never mutated after construction.
#[derive(Debug, Clone)]
pub struct PrimitiveAaGraph {
    /// Nominal peptide mass (sum of residue nominal masses).
    pub peptide_mass: i32,
    /// `true` = prefix-ion direction (b-ions dominate); derived from
    /// `scored_spec.main_ion_direction()`. Governs which end is the source.
    pub direction: bool,
    /// The smallest nominal mass that can appear as a node (may be negative
    /// for very light residues or N-terminal mods).
    pub min_node_mass: i32,
    /// `-min_node_mass`: added to a nominal mass to get its dense index.
    pub mass_offset: i32,
    /// Number of active (reachable) nodes, including source and sink.
    pub node_count: usize,
    /// Node index of the source (mass = 0).
    pub source_node_idx: usize,
    /// Node index of the sink (mass = `peptide_mass`).
    pub sink_node_idx: usize,
    /// Sorted ascending list of active nominal masses. `active_nodes[ni]` is
    /// the nominal mass of node `ni`.
    pub active_nodes: Vec<i32>,
    /// Dense array: `mass_to_node_idx[mass + mass_offset]` → node index, or
    /// `-1` if that mass is not an active node.
    pub mass_to_node_idx: Vec<i32>,
    /// CSR row offsets: incoming edges of node `ni` are stored in
    /// `edge_prev_node[edge_offset[ni]..edge_offset[ni+1]]`.
    pub edge_offset: Vec<usize>,
    /// Predecessor nominal mass for each edge.
    pub edge_prev_node: Vec<i32>,
    /// Amino-acid prior probability for each edge (default: `1/20 = 0.05`).
    pub edge_prob: Vec<f32>,
    /// Combined (cleavage + error) score for each edge.
    pub edge_score: Vec<i32>,
    /// Per-node score from the spectrum. Indexed by node index.
    /// Source (ni=0) and sink always have score 0.
    pub node_scores: Vec<i32>,
}

impl PrimitiveAaGraph {
    /// Build the graph. Mirrors Java `PrimitiveAminoAcidGraph` constructor
    /// phases 1-5 exactly.
    ///
    /// # Parameters
    ///
    /// - `aa_set`: the amino acid set (determines which AAs appear at each
    ///   position and their cleavage credits/penalties).
    /// - `peptide_mass`: nominal precursor mass in integer Da.
    /// - `enzyme`: optional enzyme for cleavage scoring at source/sink edges.
    /// - `scored_spec`: per-spectrum precomputed scoring state (node/edge scores).
    /// - `scorer`: the rank-based scoring model.
    /// - `charge`: precursor charge state.
    /// - `parent_mass`: neutral precursor mass in Da (for scoring).
    /// - `fragment_tolerance_da`: fragment mass tolerance in Da (for node scoring).
    /// - `use_protein_n_term` / `use_protein_c_term`: whether the peptide is at
    ///   the protein terminus (affects which AA list is used for source/sink).
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        aa_set: &AminoAcidSet,
        peptide_mass: i32,
        enzyme: Option<Enzyme>,
        scored_spec: &ScoredSpectrum<'_>,
        scorer: &RankScorer,
        charge: u8,
        parent_mass: f64,
        fragment_tolerance_da: f64,
        use_protein_n_term: bool,
        use_protein_c_term: bool,
    ) -> Self {
        // ---------------------------------------------------------------
        // Phase 1: Resolve source / sink AA lists.
        // ---------------------------------------------------------------
        let direction = scored_spec.main_ion_direction();

        let (source_location, sink_location) = if direction {
            // prefix direction: source = N-term, sink = C-term
            let src = if use_protein_n_term { ModLocation::ProtNTerm } else { ModLocation::NTerm };
            let snk = if use_protein_c_term { ModLocation::ProtCTerm } else { ModLocation::CTerm };
            (src, snk)
        } else {
            // suffix direction: source = C-term, sink = N-term
            let src = if use_protein_c_term { ModLocation::ProtCTerm } else { ModLocation::CTerm };
            let snk = if use_protein_n_term { ModLocation::ProtNTerm } else { ModLocation::NTerm };
            (src, snk)
        };

        let source_aas: Vec<AminoAcid> = aa_set
            .aa_list_for(source_location)
            .into_iter()
            .cloned()
            .collect();
        let anywhere_aas: Vec<AminoAcid> = aa_set
            .aa_list_for(ModLocation::Anywhere)
            .into_iter()
            .cloned()
            .collect();
        let sink_aas: Vec<AminoAcid> = aa_set
            .aa_list_for(sink_location)
            .into_iter()
            .cloned()
            .collect();

        // ---------------------------------------------------------------
        // Phase 2: Compute min_node_mass and mass_offset.
        // ---------------------------------------------------------------
        let mut min_mass: i32 = 0;
        for aa in &source_aas {
            min_mass = min_mass.min(aa.nominal_mass());
        }
        for aa in &anywhere_aas {
            min_mass = min_mass.min(1 + aa.nominal_mass());
        }
        for aa in &sink_aas {
            min_mass = min_mass.min(peptide_mass - aa.nominal_mass());
        }
        let min_node_mass = min_mass;
        let mass_offset = -min_node_mass;

        // ---------------------------------------------------------------
        // Phase 3: Reachability sweep + per-mass incoming edge counts.
        // ---------------------------------------------------------------
        let dense_len = (peptide_mass - min_node_mass + 1) as usize;
        let mut reachable = vec![false; dense_len];
        let mut in_edge_count_by_mass = vec![0_i32; dense_len];

        let to_dense = |mass: i32| -> usize { (mass + mass_offset) as usize };
        let is_representable = |mass: i32| -> bool {
            mass >= min_node_mass && mass <= peptide_mass
        };

        reachable[to_dense(0)] = true;

        // Cleavage flags (Java: addCleavageFromSource / addCleavageToSink).
        // direction == enzyme.isNTerm() → cleavage credit added at source edges.
        let add_cleavage_from_source = enzyme.is_some_and(|e| direction == e.is_n_term());
        let add_cleavage_to_sink     = enzyme.is_some_and(|e| direction != e.is_n_term());

        // Forward edges from source (mass 0).
        for aa in &source_aas {
            let next_mass = aa.nominal_mass();
            if next_mass >= peptide_mass || !is_representable(next_mass) {
                continue;
            }
            reachable[to_dense(next_mass)] = true;
            in_edge_count_by_mass[to_dense(next_mass)] += 1;
        }

        // Forward edges from intermediate nodes.
        for cur_mass in 1..peptide_mass {
            if !reachable[to_dense(cur_mass)] {
                continue;
            }
            for aa in &anywhere_aas {
                let next_mass = cur_mass + aa.nominal_mass();
                if next_mass >= peptide_mass || !is_representable(next_mass) {
                    continue;
                }
                reachable[to_dense(next_mass)] = true;
                in_edge_count_by_mass[to_dense(next_mass)] += 1;
            }
        }

        // Backward edges to sink (peptide_mass): counted in sink's in_edge_count.
        for aa in &sink_aas {
            let prev_mass = peptide_mass - aa.nominal_mass();
            if !is_representable(prev_mass) || !reachable[to_dense(prev_mass)] {
                continue;
            }
            in_edge_count_by_mass[to_dense(peptide_mass)] += 1;
        }
        reachable[to_dense(peptide_mass)] = true;

        // ---------------------------------------------------------------
        // Phase 4: Build active_nodes and mass_to_node_idx.
        // ---------------------------------------------------------------
        let count = reachable.iter().filter(|&&r| r).count();
        let node_count = count;
        let mut active_nodes = Vec::with_capacity(node_count);
        let mut mass_to_node_idx = vec![-1_i32; dense_len];

        // Source node (mass = 0) is always index 0.
        active_nodes.push(0_i32);
        mass_to_node_idx[to_dense(0)] = 0;
        let source_node_idx = 0_usize;

        for m in min_node_mass..=peptide_mass {
            if m == 0 || !reachable[to_dense(m)] {
                continue;
            }
            let idx = active_nodes.len();
            active_nodes.push(m);
            mass_to_node_idx[to_dense(m)] = idx as i32;
        }

        let sink_node_idx = mass_to_node_idx[to_dense(peptide_mass)] as usize;

        // ---------------------------------------------------------------
        // Phase 5: Build CSR edge_offset and fill edges.
        // ---------------------------------------------------------------
        let mut edge_offset = vec![0_usize; node_count + 1];
        for ni in 0..node_count {
            let mass = active_nodes[ni];
            let in_count = in_edge_count_by_mass[to_dense(mass)] as usize;
            edge_offset[ni + 1] = edge_offset[ni] + in_count;
        }
        let total_edges = edge_offset[node_count];

        let mut edge_prev_node = vec![0_i32; total_edges];
        let mut edge_prob = vec![0.0_f32; total_edges];
        let mut edge_mass_scratch = vec![0.0_f64; total_edges]; // AA accurate mass, for error score
        let mut edge_score_arr = vec![0_i32; total_edges];

        // Write cursor per node (starts at edge_offset[ni], advances as edges are written).
        let mut write_cursor = edge_offset[..node_count].to_vec();

        // Helper: write one edge into the CSR arrays.
        let get_node_idx = |mass: i32| -> i32 {
            if !is_representable(mass) {
                return -1;
            }
            mass_to_node_idx[to_dense(mass)]
        };

        // Source → intermediate edges.
        for aa in &source_aas {
            let next_mass = aa.nominal_mass();
            if next_mass >= peptide_mass || !is_representable(next_mass) {
                continue;
            }
            let target_ni = get_node_idx(next_mass);
            if target_ni < 0 {
                continue;
            }
            let target_ni = target_ni as usize;
            let cleavage_score = if add_cleavage_from_source {
                if let Some(e) = enzyme {
                    if e.is_cleavable(aa.residue) {
                        aa_set.peptide_cleavage_credit()
                    } else {
                        aa_set.peptide_cleavage_penalty()
                    }
                } else {
                    0
                }
            } else {
                0
            };
            let e_idx = write_cursor[target_ni];
            write_cursor[target_ni] += 1;
            edge_prev_node[e_idx] = 0; // prev is source (mass 0)
            edge_prob[e_idx] = aa_total_probability(aa);
            edge_mass_scratch[e_idx] = aa_total_mass(aa);
            edge_score_arr[e_idx] = cleavage_score;
        }

        // Intermediate → intermediate edges.
        for cur_mass in 1..peptide_mass {
            if !reachable[to_dense(cur_mass)] {
                continue;
            }
            for aa in &anywhere_aas {
                let next_mass = cur_mass + aa.nominal_mass();
                if next_mass >= peptide_mass || !is_representable(next_mass) {
                    continue;
                }
                let target_ni = get_node_idx(next_mass);
                if target_ni < 0 {
                    continue;
                }
                let target_ni = target_ni as usize;
                let e_idx = write_cursor[target_ni];
                write_cursor[target_ni] += 1;
                edge_prev_node[e_idx] = cur_mass;
                edge_prob[e_idx] = aa_total_probability(aa);
                edge_mass_scratch[e_idx] = aa_total_mass(aa);
                edge_score_arr[e_idx] = 0;
            }
        }

        // Backward sink edges.
        for aa in &sink_aas {
            let prev_mass = peptide_mass - aa.nominal_mass();
            if !is_representable(prev_mass) || !reachable[to_dense(prev_mass)] {
                continue;
            }
            let target_ni = get_node_idx(peptide_mass);
            if target_ni < 0 {
                continue;
            }
            let target_ni = target_ni as usize;
            let cleavage_score = if add_cleavage_to_sink {
                if let Some(e) = enzyme {
                    if e.is_cleavable(aa.residue) {
                        aa_set.peptide_cleavage_credit()
                    } else {
                        aa_set.peptide_cleavage_penalty()
                    }
                } else {
                    0
                }
            } else {
                0
            };
            let e_idx = write_cursor[target_ni];
            write_cursor[target_ni] += 1;
            edge_prev_node[e_idx] = prev_mass;
            edge_prob[e_idx] = aa_total_probability(aa);
            edge_mass_scratch[e_idx] = aa_total_mass(aa);
            edge_score_arr[e_idx] = cleavage_score;
        }

        // ---------------------------------------------------------------
        // Phase 6: Compute edge error scores.
        // ---------------------------------------------------------------
        compute_edge_error_scores(
            &active_nodes,
            &edge_offset,
            &edge_prev_node,
            &edge_mass_scratch,
            &mut edge_score_arr,
            peptide_mass,
            scored_spec,
            scorer,
            charge,
            parent_mass,
        );

        // ---------------------------------------------------------------
        // Phase 7: Compute node scores.
        // ---------------------------------------------------------------
        let node_scores = compute_node_scores(
            &active_nodes,
            peptide_mass,
            direction,
            scored_spec,
            scorer,
            charge,
            parent_mass,
            fragment_tolerance_da,
        );

        Self {
            peptide_mass,
            direction,
            min_node_mass,
            mass_offset,
            node_count,
            source_node_idx,
            sink_node_idx,
            active_nodes,
            mass_to_node_idx,
            edge_offset,
            edge_prev_node,
            edge_prob,
            edge_score: edge_score_arr,
            node_scores,
        }
    }

    // -----------------------------------------------------------------------
    // Accessors
    // -----------------------------------------------------------------------

    /// Look up the node index for a nominal mass, or `None` if the mass is
    /// not an active node. Mirrors Java `getNodeIndexForMass`.
    pub fn node_index_for_mass(&self, mass: i32) -> Option<usize> {
        if mass < self.min_node_mass || mass > self.peptide_mass {
            return None;
        }
        let idx = self.mass_to_node_idx[(mass + self.mass_offset) as usize];
        if idx < 0 { None } else { Some(idx as usize) }
    }

    /// The nominal mass of node `ni`.
    pub fn node_mass(&self, ni: usize) -> i32 {
        self.active_nodes[ni]
    }

    /// Score of node `ni`. Source and sink always have score 0.
    pub fn node_score(&self, ni: usize) -> i32 {
        self.node_scores[ni]
    }

    /// The total number of edges in the CSR graph.
    pub fn total_edges(&self) -> usize {
        self.edge_offset[self.node_count]
    }
}

// -----------------------------------------------------------------------
// Private helpers
// -----------------------------------------------------------------------

/// Standard amino acid prior probability: `1 / 20 = 0.05`.
/// Java: `AminoAcid.getProbability()`. For standard AAs the default uniform
/// prior is used. Modified AAs share the same probability as their parent.
#[inline]
fn aa_total_probability(aa: &AminoAcid) -> f32 {
    // Standard uniform prior: 1/20. Java: aa.getProbability() defaults to
    // uniform 0.05 unless a frequency model is loaded.
    const UNIFORM_PRIOR: f32 = 1.0 / 20.0;
    let _ = aa; // no per-AA prior stored yet; future: aa.probability field
    UNIFORM_PRIOR
}

/// Accurate (float) mass of the amino acid including any modification delta.
/// Java: `aa.getMass()` (which is accurate mass, not nominal).
#[inline]
fn aa_total_mass(aa: &AminoAcid) -> f64 {
    aa.mass + aa.mod_.as_ref().map_or(0.0, |m| m.mass_delta)
}

/// Phase 6 (Java `computeEdgeErrorScores`): for each intermediate node's
/// incoming edges, accumulate `scored_spec.edge_score(cur, prev, theo_aa_mass)`.
///
/// Source (mass = 0) and sink (mass = peptide_mass) nodes are skipped,
/// matching Java's `if (curMass == 0 || curMass == peptideMass) continue`.
///
/// Scores outside `[-100, 100]` are replaced with `-4` (Java behavior).
#[allow(clippy::too_many_arguments)]
fn compute_edge_error_scores(
    active_nodes: &[i32],
    edge_offset: &[usize],
    edge_prev_node: &[i32],
    edge_mass_scratch: &[f64],
    edge_score: &mut [i32],
    peptide_mass: i32,
    scored_spec: &ScoredSpectrum<'_>,
    scorer: &RankScorer,
    charge: u8,
    parent_mass: f64,
) {
    let node_count = active_nodes.len();
    for ni in 0..node_count {
        let cur_mass = active_nodes[ni];
        if cur_mass == 0 || cur_mass == peptide_mass {
            continue;
        }
        for e in edge_offset[ni]..edge_offset[ni + 1] {
            let prev_mass = edge_prev_node[e];
            let theo_aa_mass = edge_mass_scratch[e];
            let mut error_score =
                scored_spec.edge_score(cur_mass, prev_mass, theo_aa_mass, scorer, charge, parent_mass);
            if !(-100..=100).contains(&error_score) {
                eprintln!(
                    "WARN PrimitiveAaGraph: edge_score {error_score} out of range \
                     [cur={cur_mass}, prev={prev_mass}]; using -4"
                );
                error_score = -4;
            }
            edge_score[e] += error_score;
        }
    }
}

/// Phase 7 (Java `computeNodeScores`): for each intermediate node, compute
/// `scored_spec.node_score(prefix_nominal, suffix_nominal, scorer, charge,
/// parent_mass, fragment_tolerance_da)`.
///
/// - If `direction` (prefix direction): `prefix = nominal_mass`, `suffix = complement`.
/// - Else: `prefix = complement`, `suffix = nominal_mass`.
///
/// Source (ni = 0) and sink get score 0.
#[allow(clippy::too_many_arguments)]
fn compute_node_scores(
    active_nodes: &[i32],
    peptide_mass: i32,
    direction: bool,
    scored_spec: &ScoredSpectrum<'_>,
    scorer: &RankScorer,
    charge: u8,
    parent_mass: f64,
    fragment_tolerance_da: f64,
) -> Vec<i32> {
    let node_count = active_nodes.len();
    let mut node_scores = vec![0_i32; node_count];

    // ni = 0 is source; skip. Also skip sink.
    for ni in 1..node_count {
        let mass = active_nodes[ni];
        if mass == peptide_mass {
            node_scores[ni] = 0;
            continue;
        }
        let comp_mass = peptide_mass - mass;
        let (prefix_nom, suffix_nom) = if direction {
            (mass as f64, comp_mass as f64)
        } else {
            (comp_mass as f64, mass as f64)
        };
        node_scores[ni] = scored_spec.node_score(
            prefix_nom,
            suffix_nom,
            scorer,
            charge,
            parent_mass,
            fragment_tolerance_da,
        );
    }

    node_scores
}

// -----------------------------------------------------------------------
// Tests
// -----------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::aa_set::AminoAcidSetBuilder;
    use crate::amino_acid::AminoAcid;
    use crate::enzyme::Enzyme;
    use crate::param_model::{FragmentOffsetFrequency, IonType, Param, Partition, SpecDataType};
    use crate::scoring::rank_scorer::RankScorer;
    use crate::scoring::scored_spectrum::ScoredSpectrum;
    use crate::spectrum::Spectrum;
    use std::collections::HashMap;

    // -----------------------------------------------------------------------
    // Test fixtures
    // -----------------------------------------------------------------------

    fn empty_spectrum() -> Spectrum {
        Spectrum {
            title: "test".into(),
            precursor_mz: 500.0,
            precursor_intensity: None,
            precursor_charge: Some(2),
            rt_seconds: None,
            scan: None,
            peaks: vec![],
        }
    }

    fn tiny_param() -> Param {
        use crate::activation::ActivationMethod;
        use crate::instrument::InstrumentType;
        use crate::protocol::Protocol;

        let part = Partition { charge: 2, parent_mass: 1000.0, seg_num: 0 };
        let prefix1 = IonType::Prefix { charge: 1, offset_bits: 0.0_f32.to_bits() };
        let noise = IonType::Noise;

        let mut ion_table: HashMap<IonType, Vec<f32>> = HashMap::new();
        ion_table.insert(prefix1, vec![0.6_f32, 0.3, 0.05, 0.001]);
        ion_table.insert(noise, vec![0.1_f32, 0.2, 0.3, 0.4]);

        let mut rank_dist_table: HashMap<Partition, HashMap<IonType, Vec<f32>>> = HashMap::new();
        rank_dist_table.insert(part, ion_table);

        let mut frag_off_table = HashMap::new();
        frag_off_table.insert(part, vec![FragmentOffsetFrequency {
            ion_type: prefix1,
            frequency: 0.7,
        }]);

        Param {
            version: 10001,
            data_type: SpecDataType {
                activation: ActivationMethod::HCD,
                instrument: InstrumentType::QExactive,
                enzyme: None,
                protocol: Protocol::Automatic,
            },
            mme: crate::tolerance::Tolerance::Da(0.5),
            apply_deconvolution: false,
            deconvolution_error_tolerance: 0.0,
            charge_hist: vec![(2, 100)],
            min_charge: 2,
            max_charge: 2,
            num_segments: 1,
            partitions: vec![part],
            num_precursor_off: 0,
            precursor_off_map: HashMap::new(),
            frag_off_table,
            max_rank: 3,
            rank_dist_table,
            error_scaling_factor: 0,
            ion_err_dist_table: HashMap::new(),
            noise_err_dist_table: HashMap::new(),
            ion_existence_table: HashMap::new(),
        }
    }

    fn build_graph(peptide_mass: i32, enzyme: Option<Enzyme>) -> PrimitiveAaGraph {
        let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
        let spec = empty_spectrum();
        let param = tiny_param();
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new_without_filtering(&spec);
        PrimitiveAaGraph::new(
            &aa_set,
            peptide_mass,
            enzyme,
            &ss,
            &scorer,
            2,
            1000.0,
            0.5,
            false,
            false,
        )
    }

    // -----------------------------------------------------------------------
    // Required tests from the plan
    // -----------------------------------------------------------------------

    #[test]
    fn graph_for_peptide_mass_zero_has_only_source_and_sink() {
        // peptide_mass = 0: source (mass 0) == sink (mass 0) so the graph
        // degenerates to a single node.
        let g = build_graph(0, None);
        assert_eq!(g.node_count, 1, "peptide_mass=0 should yield 1 node (source=sink)");
        assert_eq!(g.source_node_idx, g.sink_node_idx);
    }

    #[test]
    fn graph_active_nodes_contain_source_and_sink() {
        // For a non-degenerate mass, source (0) and sink (peptide_mass)
        // must both be reachable.
        let g = build_graph(1000, None);
        assert!(
            g.active_nodes.contains(&0),
            "source mass 0 must be in active_nodes"
        );
        assert!(
            g.active_nodes.contains(&1000),
            "sink mass 1000 must be in active_nodes"
        );
        assert_eq!(g.active_nodes[g.source_node_idx], 0);
        assert_eq!(g.active_nodes[g.sink_node_idx], 1000);
    }

    #[test]
    fn csr_edge_offsets_are_monotonic() {
        let g = build_graph(500, None);
        for i in 0..g.node_count {
            assert!(
                g.edge_offset[i] <= g.edge_offset[i + 1],
                "edge_offset must be non-decreasing at index {i}"
            );
        }
    }

    #[test]
    fn enzyme_credit_added_to_source_edges_when_n_term_enzyme() {
        // LysN is N-terminal → direction (b-ion prefix) == enzyme.is_n_term() (true).
        // So addCleavageFromSource = true. The source edges for K should receive
        // cleavage credit, and for non-K residues should receive penalty.
        // With the default aa_set (no enzyme registered → credit=0, penalty=0),
        // the score stays 0. To test the branch we use a set with register_enzyme.
        let mut aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
        // LysN: 2 K residues, prob = 0.05 → 0.05, efficiency ≈ 0.89
        aa_set.register_enzyme(Enzyme::LysN, 0.89, 0.79);
        let spec = empty_spectrum();
        let param = tiny_param();
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new_without_filtering(&spec);
        // direction is true (prefix), LysN.is_n_term() = true → addCleavageFromSource.
        let g = PrimitiveAaGraph::new(
            &aa_set,
            1000,
            Some(Enzyme::LysN),
            &ss,
            &scorer,
            2,
            1000.0,
            0.5,
            false,
            false,
        );
        // Look at source-outgoing edges: they're stored as incoming edges of their
        // target nodes. Target nodes with edge from source have prev=0.
        let credit  = aa_set.peptide_cleavage_credit();
        let penalty = aa_set.peptide_cleavage_penalty();
        let mut found_credit  = false;
        let mut found_penalty = false;
        for ni in 0..g.node_count {
            for e in g.edge_offset[ni]..g.edge_offset[ni + 1] {
                if g.edge_prev_node[e] == 0 {
                    // This is a source edge.
                    // Target node has mass == active_nodes[ni].
                    let target_mass = g.active_nodes[ni];
                    // K residue nominal mass ≈ 128; if target == 128, it's a K edge.
                    let k_nom = AminoAcid::standard(b'K').unwrap().nominal_mass();
                    if target_mass == k_nom {
                        if g.edge_score[e] == credit { found_credit = true; }
                    } else if g.edge_score[e] == penalty {
                        found_penalty = true;
                    }
                }
            }
        }
        assert!(
            found_credit,
            "expected a source edge with K (cleavage credit {credit}) for LysN"
        );
        assert!(
            found_penalty,
            "expected a source edge with non-K residue (penalty {penalty}) for LysN"
        );
    }

    // -----------------------------------------------------------------------
    // Additional tests
    // -----------------------------------------------------------------------

    #[test]
    fn sink_node_idx_points_to_peptide_mass() {
        let pep_mass = 800_i32;
        let g = build_graph(pep_mass, None);
        assert_eq!(
            g.active_nodes[g.sink_node_idx], pep_mass,
            "sink_node_idx must point to a node with mass = peptide_mass"
        );
    }

    #[test]
    fn node_index_for_mass_returns_none_for_non_reachable() {
        let g = build_graph(500, None);
        // Mass 499 is an intermediate mass; 499 - minAA > 0 so it may or may
        // not be reachable. Mass -1 is definitely unreachable.
        assert!(
            g.node_index_for_mass(-1).is_none(),
            "negative mass is never reachable"
        );
        assert!(
            g.node_index_for_mass(g.peptide_mass + 1).is_none(),
            "mass > peptide_mass is never reachable"
        );
    }

    #[test]
    fn node_count_is_at_least_two_for_nonzero_mass() {
        // Any peptide_mass > 0 must have at least source and sink.
        let g = build_graph(200, None);
        assert!(g.node_count >= 2, "must have at least source and sink");
    }

    #[test]
    fn source_always_index_zero() {
        let g = build_graph(600, None);
        assert_eq!(g.source_node_idx, 0);
        assert_eq!(g.active_nodes[0], 0);
    }

    #[test]
    fn with_no_enzyme_no_cleavage_scores_on_intermediate_edges() {
        // Without enzyme, all cleavage scores are 0 (error score may be 0 too
        // since error_scaling_factor = 0 in tiny_param).
        let g = build_graph(300, None);
        // All edge scores should be 0 because: no enzyme → no cleavage score,
        // and tiny_param.error_scaling_factor = 0 → edge_score returns 0.
        for e in 0..g.total_edges() {
            assert_eq!(
                g.edge_score[e], 0,
                "without enzyme + zero error_scaling_factor, all edge scores must be 0"
            );
        }
    }

    #[test]
    fn node_scores_source_and_sink_are_zero() {
        let g = build_graph(400, None);
        // Source (ni = 0) must be 0.
        assert_eq!(g.node_scores[g.source_node_idx], 0);
        // Sink must be 0.
        assert_eq!(g.node_scores[g.sink_node_idx], 0);
    }

    #[test]
    fn known_peptide_node_count_peptide() {
        // PEPTIDE nominal masses: P=97, E=129, P=97, T=101, I=113, D=115, E=129.
        // Sum = 97+129+97+101+113+115+129 = 781.
        let pep_mass = 781_i32;
        let g = build_graph(pep_mass, None);
        // Source (0) and sink (781) must be present.
        assert!(g.node_index_for_mass(0).is_some());
        assert!(g.node_index_for_mass(pep_mass).is_some());
        // The graph should have intermediate nodes between 0 and 781.
        assert!(g.node_count >= 2);
    }

    #[test]
    fn trypsin_c_term_adds_cleavage_to_sink_edges() {
        // Trypsin: C-terminal enzyme → direction (true, prefix) != is_n_term (false)
        // → addCleavageToSink = true.
        // Register Trypsin with non-zero efficiencies so credit/penalty are computed.
        let mut aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
        aa_set.register_enzyme(Enzyme::Trypsin, 0.99999, 0.99999);
        let credit  = aa_set.peptide_cleavage_credit();
        let penalty = aa_set.peptide_cleavage_penalty();
        // Ensure register_enzyme produced non-trivial scores.
        assert_ne!(credit, 0, "Trypsin should produce a non-zero cleavage credit");

        let spec = empty_spectrum();
        let param = tiny_param();
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new_without_filtering(&spec);
        let g = PrimitiveAaGraph::new(
            &aa_set,
            781,
            Some(Enzyme::Trypsin),
            &ss,
            &scorer,
            2,
            1000.0,
            0.5,
            false,
            false,
        );
        // Sink edges (ni = sink_node_idx) should carry cleavage score.
        let sink_ni = g.sink_node_idx;
        let mut saw_credit  = false;
        let mut saw_penalty = false;
        for e in g.edge_offset[sink_ni]..g.edge_offset[sink_ni + 1] {
            let prev_mass = g.edge_prev_node[e];
            // The AA spanning from prev_mass to peptide_mass (781).
            let aa_nom = 781 - prev_mass;
            // K = 128, R = 156 — if aa_nom matches K or R, it should get credit.
            let k_nom = AminoAcid::standard(b'K').unwrap().nominal_mass();
            let r_nom = AminoAcid::standard(b'R').unwrap().nominal_mass();
            if aa_nom == k_nom || aa_nom == r_nom {
                if g.edge_score[e] == credit { saw_credit = true; }
            } else if g.edge_score[e] == penalty {
                saw_penalty = true;
            }
        }
        // We should observe at least one credit (K or R ending peptide) and
        // at least one penalty if the peptide has a non-KR residue at C-term.
        // Both K (128) and R (156) lead to edges if 781-128 and 781-156 are reachable.
        assert!(saw_credit, "expected at least one sink edge with cleavage credit (cleavable residue like K or R)");
        assert!(saw_penalty, "expected at least one sink edge with cleavage penalty (non-cleavable residue)");
        // Verify at least some edge has a non-zero score at the sink.
        let has_nonzero = (g.edge_offset[sink_ni]..g.edge_offset[sink_ni + 1])
            .any(|e| g.edge_score[e] != 0);
        assert!(has_nonzero, "Trypsin cleavage scoring should produce non-zero scores at sink edges");
    }

    #[test]
    fn graph_with_suffix_main_ion_swaps_node_score_arg_order() {
        // Exercise the suffix direction code path (direction = false).
        // When direction = false:
        //   - source = C-term, sink = N-term (swapped from prefix direction)
        //   - compute_node_scores swaps prefix/suffix args: (comp_mass, mass) instead of (mass, comp_mass)
        //
        // Build a ScoredSpectrum with the default prefix main ion, then mutate it to Suffix.
        let spec = empty_spectrum();
        let param = tiny_param();
        let scorer = RankScorer::new(&param);
        let mut ss = ScoredSpectrum::new_without_filtering(&spec);
        // Mutate to a Suffix ion to exercise direction = false.
        ss.set_main_ion_for_test(IonType::Suffix { charge: 1, offset_bits: 0.0_f32.to_bits() });

        // Verify main_ion_direction returns false for suffix.
        assert!(!ss.main_ion_direction(), "Suffix ion should return direction = false");

        let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
        let g = PrimitiveAaGraph::new(
            &aa_set,
            200,
            None,
            &ss,
            &scorer,
            2,
            1000.0,
            0.5,
            false,
            false,
        );

        // With direction = false:
        // - source at mass 0 becomes the C-term end (sink in prefix direction)
        // - sink at mass peptide_mass becomes the N-term end (source in prefix direction)
        assert!(!g.direction, "direction should be false for suffix ion");
        assert_eq!(g.source_node_idx, 0, "source node is always index 0");
        assert_eq!(g.active_nodes[g.source_node_idx], 0, "source node mass is always 0");
        assert_eq!(g.active_nodes[g.sink_node_idx], 200, "sink node mass is peptide_mass");
        assert!(g.node_count > 1, "graph must be non-empty (source != sink)");
    }
}
