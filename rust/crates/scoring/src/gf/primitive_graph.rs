//! Primitive-array–based amino acid graph for the generating function.
//!
//! A flat CSR replacement for the HashMap/ArrayList/NominalMass-object graph,
//! used in the DB search hot path. Graph topology is stored in CSR
//! (Compressed Sparse Row) format:
//!   `edge_offset[node+1] - edge_offset[node]` = number of incoming edges for node
//!   `edge_prev_node[e]`, `edge_prob[e]`, `edge_score[e]` = edge data
//!
//! Node scores are stored in a flat `Vec<i32>` indexed by node index.
//!
//! # Construction phases
//!
//! 1. Resolve source/sink AA lists from `direction` and protein-term flags.
//! 2. Compute `min_node_mass` and `mass_offset` from minimum nominal masses.
//! 3. Reachability sweep + per-mass incoming-edge counts.
//! 4. Build `active_nodes` and `mass_to_node_idx` dense lookup.
//! 5. Build CSR `edge_offset` and fill `edge_prev_node`, `edge_prob`, `edge_score`.
//! 6. Compute edge error scores via `scored_spec.edge_score`.
//! 7. Compute node scores via `scored_spec.node_score`.

use std::cell::RefCell;
use std::mem;

use model::aa_set::AminoAcidSet;
use model::amino_acid::AminoAcid;
use model::enzyme::Enzyme;
use model::modification::ModLocation;
use crate::scoring::rank_scorer::RankScorer;
use crate::scoring::scored_spectrum::ScoredSpectrum;

// -----------------------------------------------------------------------
// Thread-local arena pool
// -----------------------------------------------------------------------
//
// `PrimitiveAaGraph::new` allocates 11 `Vec`s per call (4 scratch +
// 7 graph-owned). On PXD001819 the graph is built ~380k times per pass,
// so we re-use the buffers across calls via a thread-local arena.
//
// Mechanism (Option B from the plan):
//   - `new_pooled` lifts the buffers out of the arena with `mem::take`
//     and `clear()`s them (length 0, capacity preserved).
//   - The buffers are populated and length-reshaped in place (`resize` /
//     `fill` / `push`) without (re)allocating provided peak capacity is
//     already sufficient — after a few hundred calls it always is.
//   - The 7 graph-owned buffers move into the returned `PrimitiveAaGraph`.
//     When the graph is dropped, `Drop` returns them to the arena
//     (`pooled = true`).
//   - The 4 scratch buffers are returned to the arena at end of
//     `new_pooled` directly.
//
// Graphs built via the legacy `new` keep `pooled = false` and skip the
// `Drop`-side roundtrip — they allocate and free as before, so existing
// callers (including tests that build many graphs without an arena) are
// unaffected.

/// Per-thread buffer pool for `PrimitiveAaGraph::new_pooled`.
///
/// Holds the 7 graph-owned buffers AND the 4 scratch buffers needed during
/// construction. When the pool is empty (first call on a thread) each Vec
/// is heap-default (no allocation); after the first build all buffers carry
/// their accumulated capacity.
#[derive(Default)]
struct PrimitiveGraphArena {
    // Graph-owned (lifted into the returned graph, returned by Drop):
    active_nodes: Vec<i32>,
    mass_to_node_idx: Vec<i32>,
    edge_offset: Vec<usize>,
    edge_prev_node: Vec<i32>,
    edge_prob: Vec<f32>,
    edge_score: Vec<i32>,
    node_scores: Vec<i32>,
    // Scratch (returned at end of new_pooled):
    reachable: Vec<bool>,
    in_edge_count_by_mass: Vec<i32>,
    edge_mass_scratch: Vec<f64>,
    write_cursor: Vec<usize>,
}

thread_local! {
    static GRAPH_ARENA: RefCell<PrimitiveGraphArena> =
        RefCell::new(PrimitiveGraphArena::default());
}

/// Take a `Vec<T>` out of the arena (length 0, capacity preserved).
#[inline]
fn take_clear<T>(slot: &mut Vec<T>) -> Vec<T> {
    let mut v = mem::take(slot);
    v.clear();
    v
}

/// Primitive CSR amino-acid graph used by the generating-function DP.
///
/// All fields are `pub` so that the GF DP can read them without accessor
/// overhead. The graph is built once per (spectrum, peptide-mass) pair and
/// is never mutated after construction.
#[derive(Debug, Clone)]
pub struct PrimitiveAaGraph {
    /// Nominal peptide mass (sum of residue nominal masses).
    pub peptide_mass: i32,
    /// `true` = prefix-ion direction (b-ions dominate); derived from
    /// `scored_spec.main_ion_direction()`. Governs which end is the source.
    pub direction: bool,
    /// Optional enzyme used during graph construction. Stored so that the
    /// GF DP can apply the neighboring-AA cleavage adjustment.
    pub enzyme: Option<Enzyme>,
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
    /// When `true`, the graph borrows its buffers from the thread-local
    /// arena and returns them on `Drop`. Set by `new_pooled`; legacy `new`
    /// keeps it `false` so existing callers behave identically.
    pooled: bool,
}

impl PrimitiveAaGraph {
    /// Build the graph by running construction phases 1-5.
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
        // Use fresh, unpooled buffers; allocates on every call. Kept for
        // tests + non-hot-path callers.
        let mut active_nodes: Vec<i32> = Vec::new();
        let mut mass_to_node_idx: Vec<i32> = Vec::new();
        let mut edge_offset: Vec<usize> = Vec::new();
        let mut edge_prev_node: Vec<i32> = Vec::new();
        let mut edge_prob: Vec<f32> = Vec::new();
        let mut edge_score: Vec<i32> = Vec::new();
        let mut node_scores: Vec<i32> = Vec::new();
        let mut reachable: Vec<bool> = Vec::new();
        let mut in_edge_count_by_mass: Vec<i32> = Vec::new();
        let mut edge_mass_scratch: Vec<f64> = Vec::new();
        let mut write_cursor: Vec<usize> = Vec::new();

        let (
            direction,
            min_node_mass,
            mass_offset,
            node_count,
            source_node_idx,
            sink_node_idx,
        ) = Self::build_in_place(
            aa_set,
            peptide_mass,
            enzyme,
            scored_spec,
            scorer,
            charge,
            parent_mass,
            fragment_tolerance_da,
            use_protein_n_term,
            use_protein_c_term,
            &mut active_nodes,
            &mut mass_to_node_idx,
            &mut edge_offset,
            &mut edge_prev_node,
            &mut edge_prob,
            &mut edge_score,
            &mut node_scores,
            &mut reachable,
            &mut in_edge_count_by_mass,
            &mut edge_mass_scratch,
            &mut write_cursor,
        );

        Self {
            peptide_mass,
            direction,
            enzyme,
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
            edge_score,
            node_scores,
            pooled: false,
        }
    }

    /// Same algorithm as `new`, but draws its 11 buffers from a thread-local
    /// arena instead of allocating fresh. The graph keeps its `pooled` flag
    /// set so `Drop` returns the 7 graph-owned buffers back to the arena.
    ///
    /// First call on a thread allocates (arena is empty); subsequent calls
    /// re-use the buffers at their accumulated peak capacity. Eliminates
    /// 11 per-call Vec allocations (~4.4M allocs per PXD001819 run).
    #[allow(clippy::too_many_arguments)]
    pub fn new_pooled(
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
        // Lift all 11 buffers out of the arena (length 0, capacity preserved).
        let (
            mut active_nodes,
            mut mass_to_node_idx,
            mut edge_offset,
            mut edge_prev_node,
            mut edge_prob,
            mut edge_score,
            mut node_scores,
            mut reachable,
            mut in_edge_count_by_mass,
            mut edge_mass_scratch,
            mut write_cursor,
        ) = GRAPH_ARENA.with(|cell| {
            let mut a = cell.borrow_mut();
            (
                take_clear(&mut a.active_nodes),
                take_clear(&mut a.mass_to_node_idx),
                take_clear(&mut a.edge_offset),
                take_clear(&mut a.edge_prev_node),
                take_clear(&mut a.edge_prob),
                take_clear(&mut a.edge_score),
                take_clear(&mut a.node_scores),
                take_clear(&mut a.reachable),
                take_clear(&mut a.in_edge_count_by_mass),
                take_clear(&mut a.edge_mass_scratch),
                take_clear(&mut a.write_cursor),
            )
        });

        let (
            direction,
            min_node_mass,
            mass_offset,
            node_count,
            source_node_idx,
            sink_node_idx,
        ) = Self::build_in_place(
            aa_set,
            peptide_mass,
            enzyme,
            scored_spec,
            scorer,
            charge,
            parent_mass,
            fragment_tolerance_da,
            use_protein_n_term,
            use_protein_c_term,
            &mut active_nodes,
            &mut mass_to_node_idx,
            &mut edge_offset,
            &mut edge_prev_node,
            &mut edge_prob,
            &mut edge_score,
            &mut node_scores,
            &mut reachable,
            &mut in_edge_count_by_mass,
            &mut edge_mass_scratch,
            &mut write_cursor,
        );

        // Return scratch buffers to the arena immediately (they outlive
        // construction but not the graph). The 7 graph-owned buffers go
        // back via Drop.
        GRAPH_ARENA.with(|cell| {
            let mut a = cell.borrow_mut();
            a.reachable = reachable;
            a.in_edge_count_by_mass = in_edge_count_by_mass;
            a.edge_mass_scratch = edge_mass_scratch;
            a.write_cursor = write_cursor;
        });

        Self {
            peptide_mass,
            direction,
            enzyme,
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
            edge_score,
            node_scores,
            pooled: true,
        }
    }

    /// Core construction algorithm. Operates in place on the 11 buffer
    /// Vecs; clears, resizes, and fills them. Returns the scalar fields
    /// that `new` / `new_pooled` need to assemble the struct.
    ///
    /// Pre-condition: all 11 buffers may be in any state (length, capacity).
    /// They will be `clear()`-ed and then resized to the lengths used by
    /// this build.
    #[allow(clippy::too_many_arguments)]
    fn build_in_place(
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
        active_nodes: &mut Vec<i32>,
        mass_to_node_idx: &mut Vec<i32>,
        edge_offset: &mut Vec<usize>,
        edge_prev_node: &mut Vec<i32>,
        edge_prob: &mut Vec<f32>,
        edge_score: &mut Vec<i32>,
        node_scores: &mut Vec<i32>,
        reachable: &mut Vec<bool>,
        in_edge_count_by_mass: &mut Vec<i32>,
        edge_mass_scratch: &mut Vec<f64>,
        write_cursor: &mut Vec<usize>,
    ) -> (bool, i32, i32, usize, usize, usize) {
        // Defensive: ensure buffers start empty (no-op when called from
        // new/new_pooled which always pass freshly-cleared Vecs).
        active_nodes.clear();
        mass_to_node_idx.clear();
        edge_offset.clear();
        edge_prev_node.clear();
        edge_prob.clear();
        edge_score.clear();
        node_scores.clear();
        reachable.clear();
        in_edge_count_by_mass.clear();
        edge_mass_scratch.clear();
        write_cursor.clear();

        // ---------------------------------------------------------------
        // Step 1: Resolve source / sink AA lists.
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

        // Borrow precomputed AA lists from the AminoAcidSet cache (populated
        // in `AminoAcidSetBuilder::build`). Avoids per-call Vec + per-AA
        // String clones; this matters because PrimitiveAaGraph::new is called
        // once per mass-bin × per spectrum (~10 × 38k = 380k calls on
        // PXD001819).
        let source_aas: &[AminoAcid] = aa_set.cached_aa_list(source_location);
        let anywhere_aas: &[AminoAcid] = aa_set.cached_aa_list(ModLocation::Anywhere);
        let sink_aas: &[AminoAcid] = aa_set.cached_aa_list(sink_location);

        // ---------------------------------------------------------------
        // Step 2: Compute min_node_mass and mass_offset.
        // ---------------------------------------------------------------
        let mut min_mass: i32 = 0;
        for aa in source_aas {
            min_mass = min_mass.min(aa.nominal_mass());
        }
        for aa in anywhere_aas {
            min_mass = min_mass.min(1 + aa.nominal_mass());
        }
        for aa in sink_aas {
            min_mass = min_mass.min(peptide_mass - aa.nominal_mass());
        }
        let min_node_mass = min_mass;
        let mass_offset = -min_node_mass;

        // ---------------------------------------------------------------
        // Step 3: Reachability sweep + per-mass incoming edge counts.
        // ---------------------------------------------------------------
        let dense_len = (peptide_mass - min_node_mass + 1) as usize;
        reachable.resize(dense_len, false);
        in_edge_count_by_mass.resize(dense_len, 0_i32);

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
        for aa in source_aas {
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
            for aa in anywhere_aas {
                let next_mass = cur_mass + aa.nominal_mass();
                if next_mass >= peptide_mass || !is_representable(next_mass) {
                    continue;
                }
                reachable[to_dense(next_mass)] = true;
                in_edge_count_by_mass[to_dense(next_mass)] += 1;
            }
        }

        // Backward edges to sink (peptide_mass): counted in sink's in_edge_count.
        for aa in sink_aas {
            let prev_mass = peptide_mass - aa.nominal_mass();
            if !is_representable(prev_mass) || !reachable[to_dense(prev_mass)] {
                continue;
            }
            in_edge_count_by_mass[to_dense(peptide_mass)] += 1;
        }
        reachable[to_dense(peptide_mass)] = true;

        // ---------------------------------------------------------------
        // Step 4: Build active_nodes and mass_to_node_idx.
        // ---------------------------------------------------------------
        let count = reachable.iter().filter(|&&r| r).count();
        let node_count = count;
        active_nodes.reserve(node_count);
        mass_to_node_idx.resize(dense_len, -1_i32);

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
        // Step 5: Build CSR edge_offset and fill edges.
        // ---------------------------------------------------------------
        edge_offset.resize(node_count + 1, 0_usize);
        // edge_offset[0] must be 0 after the resize from len=0; the loop fills
        // the rest. (resize from 0 -> node_count+1 appends node_count+1 zeros.)
        for ni in 0..node_count {
            let mass = active_nodes[ni];
            let in_count = in_edge_count_by_mass[to_dense(mass)] as usize;
            edge_offset[ni + 1] = edge_offset[ni] + in_count;
        }
        let total_edges = edge_offset[node_count];

        edge_prev_node.resize(total_edges, 0_i32);
        edge_prob.resize(total_edges, 0.0_f32);
        edge_mass_scratch.resize(total_edges, 0.0_f64); // AA accurate mass, for error score
        edge_score.resize(total_edges, 0_i32);

        // Write cursor per node (starts at edge_offset[ni], advances as edges are written).
        write_cursor.extend_from_slice(&edge_offset[..node_count]);

        // Helper: write one edge into the CSR arrays.
        let get_node_idx = |mass: i32| -> i32 {
            if !is_representable(mass) {
                return -1;
            }
            mass_to_node_idx[to_dense(mass)]
        };

        // Source → intermediate edges.
        for aa in source_aas {
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
            edge_score[e_idx] = cleavage_score;
        }

        // Intermediate → intermediate edges.
        for cur_mass in 1..peptide_mass {
            if !reachable[to_dense(cur_mass)] {
                continue;
            }
            for aa in anywhere_aas {
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
                edge_score[e_idx] = 0;
            }
        }

        // Backward sink edges.
        for aa in sink_aas {
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
            edge_score[e_idx] = cleavage_score;
        }

        // ---------------------------------------------------------------
        // Step 6: Compute edge error scores.
        // ---------------------------------------------------------------
        compute_edge_error_scores(
            active_nodes,
            edge_offset,
            edge_prev_node,
            edge_mass_scratch,
            edge_score,
            peptide_mass,
            mass_offset,
            scored_spec,
            scorer,
            charge,
            parent_mass,
        );

        // ---------------------------------------------------------------
        // Step 7: Compute node scores.
        // ---------------------------------------------------------------
        compute_node_scores_in_place(
            active_nodes,
            peptide_mass,
            direction,
            scored_spec,
            scorer,
            charge,
            parent_mass,
            fragment_tolerance_da,
            node_scores,
        );

        (direction, min_node_mass, mass_offset, node_count, source_node_idx, sink_node_idx)
    }

    // -----------------------------------------------------------------------
    // Accessors
    // -----------------------------------------------------------------------

    /// Look up the node index for a nominal mass, or `None` if the mass is
    /// not an active node.
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

impl Drop for PrimitiveAaGraph {
    fn drop(&mut self) {
        if !self.pooled {
            return;
        }
        // Return the 7 graph-owned buffers to the thread-local arena.
        // Each `mem::take` swaps in an empty Vec (capacity 0) — but that
        // empty Vec gets dropped immediately, while the populated buffer
        // (with grown capacity) goes back into the arena slot.
        //
        // If a borrow on the arena is already held (e.g. panic during
        // arena callback), we silently leak the capacity rather than
        // double-borrow-panic; the buffers themselves get freed normally.
        let _ = GRAPH_ARENA.try_with(|cell| {
            if let Ok(mut a) = cell.try_borrow_mut() {
                a.active_nodes = mem::take(&mut self.active_nodes);
                a.mass_to_node_idx = mem::take(&mut self.mass_to_node_idx);
                a.edge_offset = mem::take(&mut self.edge_offset);
                a.edge_prev_node = mem::take(&mut self.edge_prev_node);
                a.edge_prob = mem::take(&mut self.edge_prob);
                a.edge_score = mem::take(&mut self.edge_score);
                a.node_scores = mem::take(&mut self.node_scores);
            }
        });
    }
}

// -----------------------------------------------------------------------
// Private helpers
// -----------------------------------------------------------------------

/// Standard amino acid prior probability: `1 / 20 = 0.05`. Modified AAs
/// share the same probability as their parent.
#[inline]
fn aa_total_probability(aa: &AminoAcid) -> f32 {
    // Uniform prior 1/20 unless a frequency model is loaded.
    const UNIFORM_PRIOR: f32 = 1.0 / 20.0;
    let _ = aa; // no per-AA prior stored yet; future: aa.probability field
    UNIFORM_PRIOR
}

/// Accurate (float) mass of the amino acid including any modification delta.
#[inline]
fn aa_total_mass(aa: &AminoAcid) -> f64 {
    aa.mass + aa.mod_.as_ref().map_or(0.0, |m| m.mass_delta)
}

/// For each intermediate node's incoming edges, accumulate an inlined
/// edge-error score: a precomputed `ion_existence_score[idx]` plus,
/// when both endpoints have an observed peak (`idx == 3`), an additional
/// `scorer.error_score(part, delta)` term where `delta = obs(cur) -
/// obs(prev) - theo_aa_mass`.
///
/// Constants hoisted out of the per-edge inner loop (see the perf commit
/// in the git history for the rationale):
/// - the `partition_for(charge, parent_mass, last_seg)` lookup
/// - the 4-entry `ion_existence_score[0..=3]` table
/// - per-node `observed_node_mass` results (cached into a dense
///   `Vec<Option<f64>>` keyed by `mass + mass_offset`)
///
/// Source (mass = 0) and sink (mass = peptide_mass) nodes are skipped.
/// Scores outside `[-100, 100]` are replaced with `-4`.
#[allow(clippy::too_many_arguments)]
fn compute_edge_error_scores(
    active_nodes: &[i32],
    edge_offset: &[usize],
    edge_prev_node: &[i32],
    edge_mass_scratch: &[f64],
    edge_score: &mut [i32],
    peptide_mass: i32,
    mass_offset: i32,
    scored_spec: &ScoredSpectrum<'_>,
    scorer: &RankScorer,
    charge: u8,
    parent_mass: f64,
) {
    let node_count = active_nodes.len();

    // Spectrum-constant short-circuit: if either fast-out condition is true,
    // every edge gets score 0. Done once for the whole graph instead of
    // per-edge inside ScoredSpectrum::edge_score (~24k calls saved per graph
    // on PXD001819).
    if scorer.param().error_scaling_factor == 0
        || scorer.param().ion_existence_table.is_empty()
    {
        return;
    }

    // Spectrum-constant: partition for this (charge, parent_mass, last_seg).
    // Hoisted out of the per-edge inner loop — was the per-call partition_for
    // binary search inside edge_score, now done once per graph build.
    let last_seg = (scorer.param().num_segments - 1).max(0) as usize;
    let part = scorer.param().partition_for(charge, parent_mass, last_seg);
    let prob_peak = scored_spec.prob_peak;

    // Spectrum-constant: ion_existence_score for each of the 4 possible
    // ion_existence_index values (0..=3). Replaces the per-edge table lookup
    // in scorer.ion_existence_score.
    let ies = [
        scorer.ion_existence_score(part, 0, prob_peak),
        scorer.ion_existence_score(part, 1, prob_peak),
        scorer.ion_existence_score(part, 2, prob_peak),
        scorer.ion_existence_score(part, 3, prob_peak),
    ];

    // Graph-constant: observed peak mass for each node, keyed by dense mass
    // index `(mass + mass_offset)`. Each unique node mass is observed at most
    // once instead of once per outgoing edge. With ~18 edges per node on
    // PXD001819 that's an ~18× reduction in the dominant inner cost of
    // edge_score. `None` entries mark masses with no qualifying peak in
    // the tolerance window.
    let dense_len = (peptide_mass + mass_offset + 1) as usize;
    // Pre-fill with None so unreachable masses don't need explicit insertion.
    // Allocates once per graph (~1.3k entries on PXD001819); cheaper than
    // re-observing 18× per node.
    let mut observed_by_mass: Vec<Option<f64>> = vec![None; dense_len];
    for &m in active_nodes {
        let idx = (m + mass_offset) as usize;
        observed_by_mass[idx] = scored_spec.observed_node_mass(m, scorer, charge, parent_mass);
    }

    let mut clamp_count: u32 = 0;
    for ni in 0..node_count {
        let cur_mass = active_nodes[ni];
        if cur_mass == 0 || cur_mass == peptide_mass {
            continue;
        }
        let cur_obs = observed_by_mass[(cur_mass + mass_offset) as usize];
        for e in edge_offset[ni]..edge_offset[ni + 1] {
            let prev_mass = edge_prev_node[e];
            // prev_mass should always be a valid representable mass for an
            // edge written by build_in_place — fall through to None for
            // safety if it somehow isn't.
            let prev_obs = if prev_mass + mass_offset >= 0
                && (prev_mass + mass_offset) as usize <= peptide_mass as usize + mass_offset as usize
            {
                observed_by_mass
                    .get((prev_mass + mass_offset) as usize)
                    .copied()
                    .flatten()
            } else {
                None
            };

            // ion_existence_index: 1 if cur observed, +2 if prev observed.
            let mut idx = 0usize;
            if cur_obs.is_some() { idx += 1; }
            if prev_obs.is_some() { idx += 2; }

            let mut s = ies[idx];
            if idx == 3 {
                let delta = cur_obs.unwrap() - prev_obs.unwrap() - edge_mass_scratch[e];
                s += scorer.error_score(part, delta as f32);
            }
            let mut error_score = s.round() as i32;
            if !(-100..=100).contains(&error_score) {
                clamp_count += 1;
                error_score = -4;
            }
            edge_score[e] += error_score;
        }
    }
    // Emit a single aggregated warning rather than one line per offending edge
    // (this loop is hot — per-edge stderr output can spam millions of lines).
    if clamp_count > 0 {
        eprintln!(
            "WARN: PrimitiveAaGraph: {} edge score(s) clamped (out of [-100, 100] range)",
            clamp_count
        );
    }
}

/// For each intermediate node, compute
/// `scored_spec.node_score(prefix_nominal, suffix_nominal, scorer, charge,
/// parent_mass, fragment_tolerance_da)`.
///
/// - If `direction` (prefix direction): `prefix = nominal_mass`, `suffix = complement`.
/// - Else: `prefix = complement`, `suffix = nominal_mass`.
///
/// Source (ni = 0) and sink get score 0.
///
/// Writes results into `node_scores` (pre-condition: empty Vec, gets resized
/// to `active_nodes.len()`).
#[allow(clippy::too_many_arguments)]
fn compute_node_scores_in_place(
    active_nodes: &[i32],
    peptide_mass: i32,
    direction: bool,
    scored_spec: &ScoredSpectrum<'_>,
    scorer: &RankScorer,
    charge: u8,
    parent_mass: f64,
    fragment_tolerance_da: f64,
    node_scores: &mut Vec<i32>,
) {
    let node_count = active_nodes.len();
    node_scores.resize(node_count, 0_i32);

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
}

// -----------------------------------------------------------------------
// Tests
// -----------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use model::aa_set::AminoAcidSetBuilder;
    use model::amino_acid::AminoAcid;
    use model::enzyme::Enzyme;
    use crate::param_model::IonType;
    use crate::scoring::rank_scorer::RankScorer;
    use crate::scoring::scored_spectrum::ScoredSpectrum;
    use model::spectrum::Spectrum;
    use crate::testutil::tiny_param_with_ions;

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

    fn build_graph(peptide_mass: i32, enzyme: Option<Enzyme>) -> PrimitiveAaGraph {
        let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
        let spec = empty_spectrum();
        let param = tiny_param_with_ions();
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
        let param = tiny_param_with_ions();
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
        let param = tiny_param_with_ions();
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
        let param = tiny_param_with_ions();
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
