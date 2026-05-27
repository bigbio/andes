//! Generating-function DP: computes the score distribution
//! `P(score | random peptide of given nominal mass)`.
//!
//! # Uniform-prior DP
//!
//! `compute_uniform` takes a generic increment-score callback and uses a
//! uniform AA prior (`1/N`). Kept for tests and reference; not used in the
//! production search path.
//!
//! # Graph-based DP
//!
//! `compute` (and `with_score_threshold`) operate on a pre-built
//! `PrimitiveAaGraph` and produce a single final `ScoreDist` (plus enzyme
//! adjustment).

use model::aa_set::AminoAcidSet;
use crate::gf::primitive_graph::PrimitiveAaGraph;
use crate::gf::score_dist::{ScoreBound, ScoreDist};

/// Errors returned by the graph-based GF DP.
#[derive(thiserror::Error, Debug)]
pub enum GfError {
    #[error("score range is empty: min_score {min} >= max_score {max}")]
    EmptyScoreRange { min: i32, max: i32 },
    #[error("aa_masses is empty")]
    NoAminoAcids,
    #[error("sink node has no reachable distribution")]
    SinkUnreachable,
}

/// Result of the generating-function DP. Stores the final per-peptide score
/// distribution and allows querying the spectral probability.
#[derive(Debug, Clone)]
pub struct GeneratingFunction {
    /// One ScoreDist per nominal mass in 0..=max_mass (for `compute_uniform`),
    /// or exactly one element (the final adjusted dist) for `compute`.
    score_dists: Vec<ScoreDist>,
    score_bound: ScoreBound,
    /// Diagnostic-only — exposes internal DP state for tracing.
    ///
    /// Populated only when the GF is built via
    /// [`GeneratingFunction::with_score_threshold_retain_node_dists`] (the
    /// production `compute` / `with_score_threshold` paths leave this `None`
    /// so the per-node DP buffer is freed at the end of `compute_inner`).
    /// Tuples are `(node_idx, node_mass, dist)`, in node-index order.
    node_dists: Option<Vec<(usize, i32, ScoreDist)>>,
}

impl GeneratingFunction {
    // -----------------------------------------------------------------------
    // Graph-based public API
    // -----------------------------------------------------------------------

    /// Compute the GF over a precomputed primitive graph.
    pub fn compute(graph: &PrimitiveAaGraph, aa_set: &AminoAcidSet) -> Result<Self, GfError> {
        compute_inner(graph, aa_set, None, false)
    }

    /// Pre-pass: prune nodes whose maximum possible final score is below
    /// `score_threshold`. Computes `min_score_by_node` and uses it to skip
    /// irrelevant DP work.
    pub fn with_score_threshold(
        graph: &PrimitiveAaGraph,
        score_threshold: i32,
        aa_set: &AminoAcidSet,
    ) -> Result<Self, GfError> {
        let min_score_by_node = setup_score_threshold(graph, aa_set, score_threshold);
        compute_inner(graph, aa_set, Some(min_score_by_node), false)
    }

    /// Diagnostic-only — same DP as [`with_score_threshold`] but additionally
    /// retains the per-node `ScoreDist` buffer so `iter_node_dists` can dump
    /// it for tracing. Do NOT use on the production search path: it disables
    /// the per-node `.take()` cleanup, increasing peak memory by the size of
    /// the DP table.
    pub fn with_score_threshold_retain_node_dists(
        graph: &PrimitiveAaGraph,
        score_threshold: i32,
        aa_set: &AminoAcidSet,
    ) -> Result<Self, GfError> {
        let min_score_by_node = setup_score_threshold(graph, aa_set, score_threshold);
        compute_inner(graph, aa_set, Some(min_score_by_node), true)
    }

    /// The final (enzyme-adjusted) score distribution.
    pub fn score_dist(&self) -> &ScoreDist {
        // For the graph-based path this is always index 0.
        &self.score_dists[0]
    }

    /// Minimum score (inclusive) of the final distribution.
    pub fn min_score(&self) -> i32 {
        self.score_bound.min_score()
    }

    /// Maximum score (exclusive) of the final distribution.
    pub fn max_score(&self) -> i32 {
        self.score_bound.max_score()
    }

    /// Cumulative tail probability `P(random_score >= score)`.
    pub fn spectral_probability(&self, score: i32) -> f64 {
        let dist = &self.score_dists[0];
        if !dist.is_prob_set() {
            return 1.0;
        }
        dist.get_spectral_probability(score)
    }

    /// Diagnostic-only — exposes internal DP state for tracing.
    ///
    /// Yields `(node_idx, node_mass, &ScoreDist)` for every node retained by
    /// the DP. Returns an empty iterator unless the GF was built via
    /// [`Self::with_score_threshold_retain_node_dists`].
    pub fn iter_node_dists(&self) -> impl Iterator<Item = (usize, i32, &ScoreDist)> {
        self.node_dists
            .iter()
            .flat_map(|v| v.iter().map(|(ni, m, d)| (*ni, *m, d)))
    }

    // -----------------------------------------------------------------------
    // Uniform-prior DP
    // -----------------------------------------------------------------------

    /// Compute the generating function up to `max_mass`. The
    /// `increment_score` callback returns the score added when the
    /// peptide is extended by amino acid `aa_idx` (an index into
    /// `aa_masses`) at mass position `mass`.
    ///
    /// Probability prior over amino acids: uniform `1 / aa_masses.len()`.
    pub fn compute_uniform<F>(
        max_mass: i32,
        score_bound: ScoreBound,
        aa_masses: &[i32],
        increment_score: F,
    ) -> Self
    where
        F: Fn(i32, u8) -> i32,
    {
        if aa_masses.is_empty() {
            // Caller error; return an empty GF.
            return Self {
                score_dists: Vec::new(),
                score_bound,
                node_dists: None,
            };
        }
        let num_aas = aa_masses.len();
        let prior = 1.0 / num_aas as f64;

        let mut score_dists: Vec<ScoreDist> = (0..=max_mass)
            .map(|_| ScoreDist::new(score_bound.min_score(), score_bound.max_score(), false, true))
            .collect();

        // Base case: mass 0 has full probability at score 0.
        if score_bound.min_score() <= 0 && 0 < score_bound.max_score() {
            score_dists[0].set_prob(0, 1.0);
        }

        // Forward DP.
        for m in 1..=max_mass {
            let m_idx = m as usize;
            for (aa_idx, &aa_mass) in aa_masses.iter().enumerate() {
                if m - aa_mass < 0 {
                    continue;
                }
                let pred_idx = (m - aa_mass) as usize;
                let inc = increment_score(m, aa_idx as u8);

                // Iterate over the predecessor's entire score range.
                let pred_min = score_dists[pred_idx].min_score();
                let pred_max = score_dists[pred_idx].max_score();
                for s in pred_min..pred_max {
                    let p = score_dists[pred_idx].get_probability(s);
                    if p == 0.0 {
                        continue;
                    }
                    let target_s = s + inc;
                    if target_s < score_bound.min_score() || target_s >= score_bound.max_score() {
                        continue;
                    }
                    score_dists[m_idx].add_prob(target_s, p * prior);
                }
            }
        }

        Self {
            score_dists,
            score_bound,
            node_dists: None,
        }
    }

    pub fn score_bound(&self) -> ScoreBound {
        self.score_bound
    }

    pub fn score_dist_at(&self, mass: i32) -> Option<&ScoreDist> {
        if mass < 0 {
            return None;
        }
        self.score_dists.get(mass as usize)
    }

    /// Total spectral probability at the given mass and score: P(X >= score).
    /// Used by the uniform-prior path.
    pub fn spectral_probability_at(&self, mass: i32, score: i32) -> Option<f64> {
        self.score_dist_at(mass).map(|d| d.get_spectral_probability(score))
    }
}

// -----------------------------------------------------------------------
// Graph-based DP — private implementation
// -----------------------------------------------------------------------

/// Pre-pass that propagates the score threshold backward through the graph.
///
/// Returns a `min_score_by_node` array of length `graph.node_count` where
/// `min_score_by_node[ni]` is the minimum score needed at node `ni` for a
/// path from `ni` to the sink to reach >= `score_threshold`.
/// Nodes that cannot reach `score_threshold` keep `i32::MAX`.
fn setup_score_threshold(
    graph: &PrimitiveAaGraph,
    aa_set: &AminoAcidSet,
    score_threshold: i32,
) -> Vec<i32> {
    let node_count = graph.node_count;
    let source_idx = graph.source_node_idx;
    let sink_idx = graph.sink_node_idx;

    // Adjust threshold for enzyme neighboring-AA credit.
    let adjusted_score = if graph.enzyme.is_some() {
        score_threshold - aa_set.neighboring_aa_cleavage_credit()
    } else {
        score_threshold
    };

    let mut min_score_by_node = vec![i32::MAX; node_count];
    min_score_by_node[sink_idx] = adjusted_score;

    // Propagate from sink backward through sink's own incoming edges.
    for e in graph.edge_offset[sink_idx]..graph.edge_offset[sink_idx + 1] {
        let prev_mass = graph.edge_prev_node[e];
        if let Some(prev_idx) = graph.node_index_for_mass(prev_mass) {
            let new_min = adjusted_score.saturating_sub(graph.edge_score[e]);
            if new_min < min_score_by_node[prev_idx] {
                min_score_by_node[prev_idx] = new_min;
            }
        }
    }

    // Walk nodes in reverse order (from sink toward source).
    for ni in (0..node_count).rev() {
        if ni == source_idx || ni == sink_idx {
            continue;
        }
        if min_score_by_node[ni] == i32::MAX {
            continue;
        }
        let cur_mass = graph.active_nodes[ni];
        if cur_mass == graph.peptide_mass {
            continue;
        }
        let cur_node_score = graph.node_scores[ni];

        for e in graph.edge_offset[ni]..graph.edge_offset[ni + 1] {
            let prev_mass = graph.edge_prev_node[e];
            if let Some(prev_idx) = graph.node_index_for_mass(prev_mass) {
                let new_min = min_score_by_node[ni]
                    .saturating_sub(cur_node_score)
                    .saturating_sub(graph.edge_score[e]);
                if new_min < min_score_by_node[prev_idx] {
                    min_score_by_node[prev_idx] = new_min;
                }
            }
        }
    }

    min_score_by_node
}

/// Per-node header into the flat `ScoreDistArena` storage.
///
/// `start..start+len` is the half-open f64 slice for this node's
/// `prob_distribution`. `min_score` is the lowest score covered; the
/// score at storage index `start + k` is `min_score + k`. `is_set` flips
/// `false → true` the first time the node is populated by the DP, taking
/// the role of the `Option::None` sentinel in the legacy DP.
#[derive(Debug, Clone, Copy)]
struct NodeSlice {
    start: u32,
    len: u32,
    min_score: i32,
    is_set: bool,
}

impl NodeSlice {
    const UNSET: NodeSlice = NodeSlice {
        start: 0,
        len: 0,
        min_score: 0,
        is_set: false,
    };

    #[inline]
    fn range(&self) -> std::ops::Range<usize> {
        let s = self.start as usize;
        s..s + self.len as usize
    }
}

/// Flat-arena replacement for `Vec<Option<ScoreDist>>`. A single contiguous
/// `Vec<f64>` backs the probability arrays of every node; per-node headers
/// describe slice ranges. Replaces ~node_count tiny `Vec<f64>` allocations
/// (one per node, summed to ~55M per PXD001819 run) with one moderately
/// sized allocation per graph (~96 KB typical).
struct ScoreDistArena {
    storage: Vec<f64>,
    headers: Vec<NodeSlice>,
    /// Length of the next free region in `storage`; `storage[..fill]` is
    /// the populated prefix. Used by `reserve_slot` to bump-allocate
    /// per-node slices as nodes are visited.
    fill: usize,
}

impl ScoreDistArena {
    fn new(node_count: usize, initial_capacity: usize) -> Self {
        Self {
            storage: Vec::with_capacity(initial_capacity),
            headers: vec![NodeSlice::UNSET; node_count],
            fill: 0,
        }
    }

    /// Reserve a slot for node `ni` spanning scores `[min_score, max_score)`.
    /// Returns the offset of the freshly zeroed slice within `storage`.
    ///
    /// Grows `storage` if necessary. Callers must NOT hold any borrows into
    /// `storage` across a `reserve_slot` call (growth may relocate the
    /// backing buffer). The DP body honors this: it only calls
    /// `reserve_slot` once per outer-loop iteration, before any
    /// `split_at_mut` borrows are taken.
    fn reserve_slot(&mut self, ni: usize, min_score: i32, max_score: i32) -> usize {
        let len = (max_score - min_score) as usize;
        let start = self.fill;
        let needed = start + len;
        if needed > self.storage.len() {
            // Grow with zero-fill so the slice we hand out is initialized.
            self.storage.resize(needed, 0.0);
        } else {
            // Reusing existing capacity (unlikely on first pass, but the
            // resize() above might over-allocate on subsequent growth
            // cycles; either way zero the slice).
            for slot in &mut self.storage[start..start + len] {
                *slot = 0.0;
            }
        }
        self.headers[ni] = NodeSlice {
            start: start as u32,
            len: len as u32,
            min_score,
            is_set: true,
        };
        self.fill += len;
        start
    }

    /// Materialize the slice for node `ni` as an owned `ScoreDist` (used
    /// for the sink and for `retain_node_dists` snapshots).
    fn to_score_dist(&self, ni: usize) -> Option<ScoreDist> {
        let hdr = self.headers[ni];
        if !hdr.is_set {
            return None;
        }
        let mut d = ScoreDist::new(
            hdr.min_score,
            hdr.min_score + hdr.len as i32,
            false,
            true,
        );
        let slice = &self.storage[hdr.range()];
        for (i, &v) in slice.iter().enumerate() {
            // get_probability/set_prob both index from min_score, so
            // index k corresponds to score (min_score + k).
            d.set_prob(hdr.min_score + i as i32, v);
        }
        Some(d)
    }
}

/// Core DP for the graph-based generating function.
///
/// Uses a flat-arena `ScoreDistArena` for per-node probability buffers: one
/// `Vec<f64>` allocation per graph instead of `node_count` tiny allocations
/// (one per `Option<ScoreDist>::Some(_)`). Semantics are bit-identical to
/// the previous `Vec<Option<ScoreDist>>` implementation; the equivalence
/// is gated by per-peptide-mass parity fixtures.
///
/// `retain_node_dists` is a diagnostic-only flag: when `true`, each visited
/// node's probability slice is materialized into a `ScoreDist` and stashed
/// on `GeneratingFunction.node_dists` so the caller can dump it via
/// `iter_node_dists`. The production path passes `false`.
fn compute_inner(
    graph: &PrimitiveAaGraph,
    aa_set: &AminoAcidSet,
    min_score_by_node: Option<Vec<i32>>,
    retain_node_dists: bool,
) -> Result<GeneratingFunction, GfError> {
    let node_count = graph.node_count;
    let source_idx = graph.source_node_idx;
    let sink_idx = graph.sink_node_idx;

    // Estimate initial arena capacity: typical per-node score range is ~80;
    // we pick 256 to absorb deeper, higher-mass graphs without reallocating
    // mid-DP. The arena grows via `Vec::resize` if a node exceeds the
    // estimate — growth happens BEFORE any in-flight slice borrows are
    // taken, so it cannot invalidate a `split_at_mut` view.
    let initial_capacity = 1 // source slot
        + node_count.saturating_mul(256);
    let mut arena = ScoreDistArena::new(node_count, initial_capacity);

    // Debug-only counter: tracks how many nodes were skipped due to the
    // score-range guard (|score| > 10000). Fires only in debug builds;
    // release builds compile this out entirely (no perf regression).
    #[cfg(debug_assertions)]
    let mut score_range_overflow_count: u32 = 0;

    // Source has full probability at score 0.
    {
        let start = arena.reserve_slot(source_idx, 0, 1);
        arena.storage[start] = 1.0;
    }

    // Scratch buffer for valid edge indices.
    let max_edges_per_node = (0..node_count)
        .map(|ni| graph.edge_offset[ni + 1] - graph.edge_offset[ni])
        .max()
        .unwrap_or(0);
    let mut valid_edges: Vec<usize> = Vec::with_capacity(max_edges_per_node);

    // Forward DP over nodes in index order.
    for ni in 0..node_count {
        if ni == source_idx {
            continue;
        }

        let cur_node_score = graph.node_scores[ni];

        // Skip if this node is pruned by the threshold pre-pass.
        if let Some(ref msbn) = min_score_by_node {
            if msbn[ni] == i32::MAX {
                continue;
            }
        }

        // Determine initial cur_min_score.
        let mut cur_min_score: i32 = match min_score_by_node {
            Some(ref msbn) => msbn[ni],
            None => i32::MAX,
        };
        let mut cur_max_score: i32 = i32::MIN;

        valid_edges.clear();

        // Scan incoming edges.
        for e in graph.edge_offset[ni]..graph.edge_offset[ni + 1] {
            let prev_mass = graph.edge_prev_node[e];
            let prev_idx = match graph.node_index_for_mass(prev_mass) {
                Some(idx) => idx,
                None => continue,
            };
            let prev_hdr = arena.headers[prev_idx];
            if !prev_hdr.is_set {
                continue;
            }

            let combined_score = cur_node_score + graph.edge_score[e];
            let prev_max = prev_hdr.min_score + prev_hdr.len as i32;
            let possible_max = prev_max + combined_score;
            if possible_max > cur_max_score {
                cur_max_score = possible_max;
            }

            // Only update min from predecessor when NOT using threshold pre-pass.
            if min_score_by_node.is_none() {
                let possible_min = prev_hdr.min_score + combined_score;
                if possible_min < cur_min_score {
                    cur_min_score = possible_min;
                }
            }

            valid_edges.push(e);
        }

        // Skip degenerate or out-of-bound ranges.
        let valid_count = valid_edges.len();
        if cur_min_score >= cur_max_score || valid_count == 0 {
            continue;
        }
        if cur_min_score < -10000 || cur_max_score > 10000 {
            #[cfg(debug_assertions)]
            {
                score_range_overflow_count += 1;
            }
            continue;
        }

        // Reserve cur_dist slice in the arena.
        let cur_start = arena.reserve_slot(ni, cur_min_score, cur_max_score);
        let cur_len = (cur_max_score - cur_min_score) as usize;

        // Fill cur_dist by accumulating from each predecessor.
        // `split_at_mut` is required to borrow `storage` immutably (predecessor
        // slice) and mutably (cur_dist slice) simultaneously. The cur_dist
        // slice was just appended to the end of `storage`, so all predecessor
        // slices live in `storage[..cur_start]`.
        let (prev_region, cur_region) = arena.storage.split_at_mut(cur_start);
        let cur_slice = &mut cur_region[..cur_len];

        for &e in &valid_edges {
            let prev_mass = graph.edge_prev_node[e];
            // Safety: we already verified these are valid above.
            let prev_idx = graph.node_index_for_mass(prev_mass).unwrap();
            let prev_hdr = arena.headers[prev_idx];
            let prev_slice = &prev_region[prev_hdr.range()];
            let combined_score = cur_node_score + graph.edge_score[e];
            let aa_prob = graph.edge_prob[e] as f64;

            // Mirror ScoreDist::add_prob_dist:
            //   for t in max(other_min, self_min - score_diff)
            //          .. min(other_max, self_max - score_diff):
            //     self[t + score_diff - self_min] += other[t - other_min] * aa_prob
            //
            // Inner loop is split into 4-wide chunks so LLVM can auto-vectorize
            // on AVX2 / NEON. `dst_idx - src_idx = combined_score + other_min -
            // self_min` is a constant offset, so each chunk's 4 writes hit
            // distinct indices and the chunked form is bit-identical to the
            // scalar loop. Parity is gated by
            // `tests/add_prob_dist_chunked_parity.rs` (covers the standalone
            // `ScoreDist::add_prob_dist` method, which has the same structure).
            let other_min = prev_hdr.min_score;
            let other_max = prev_hdr.min_score + prev_hdr.len as i32;
            let self_min = cur_min_score;
            let self_max = cur_max_score;
            let t_start = other_min.max(self_min - combined_score);
            let t_end = other_max.min(self_max - combined_score);
            if t_end > t_start {
                let len = (t_end - t_start) as usize;
                let src_base = (t_start - other_min) as usize;
                let dst_base = (t_start + combined_score - self_min) as usize;
                let chunks = len / 4;
                for c in 0..chunks {
                    let s = src_base + c * 4;
                    let d = dst_base + c * 4;
                    cur_slice[d    ] += prev_slice[s    ] * aa_prob;
                    cur_slice[d + 1] += prev_slice[s + 1] * aa_prob;
                    cur_slice[d + 2] += prev_slice[s + 2] * aa_prob;
                    cur_slice[d + 3] += prev_slice[s + 3] * aa_prob;
                }
                let tail_start = chunks * 4;
                for r in tail_start..len {
                    cur_slice[dst_base + r] += prev_slice[src_base + r] * aa_prob;
                }
            }
        }

        // compute_inner already tightly written; further perf needs algorithmic changes
        // outside this iteration (e.g. caching prev_idx alongside valid_edges to avoid
        // the second node_index_for_mass call, or SIMD-widening the inner multiply loop).

        // Underflow guard at max_score - 1.
        // Read-then-write on the same slice; `cur_slice` is already &mut.
        let guard_idx = (cur_max_score - 1 - cur_min_score) as usize;
        if cur_slice[guard_idx] == 0.0 {
            // Use the smallest positive denormal f32 (~1.4e-45) as the
            // underflow floor — NOT `f32::MIN_POSITIVE` (smallest positive
            // normal ~1.18e-38). The denormal value matches the GF tail's
            // expected dynamic range.
            cur_slice[guard_idx] = f32::from_bits(1) as f64;
        }
    }

    // Debug-only: surface score-range overflow count before returning.
    #[cfg(debug_assertions)]
    if score_range_overflow_count > 0 {
        eprintln!(
            "[GF DP debug] score-range cutoff fired for {} node(s); \
             some nodes may not be reachable",
            score_range_overflow_count
        );
    }

    // Diagnostic-only: snapshot per-node dists. Production path leaves this
    // as `None`, identical to prior behavior.
    let node_dists_snapshot: Option<Vec<(usize, i32, ScoreDist)>> = if retain_node_dists {
        let mut snap: Vec<(usize, i32, ScoreDist)> = Vec::new();
        for ni in 0..node_count {
            if let Some(d) = arena.to_score_dist(ni) {
                snap.push((ni, graph.node_mass(ni), d));
            }
        }
        Some(snap)
    } else {
        None
    };

    // Extract sink distribution.
    let sink_dist = arena
        .to_score_dist(sink_idx)
        .ok_or(GfError::SinkUnreachable)?;

    let min_score = sink_dist.min_score();
    let max_score = sink_dist.max_score();

    if max_score <= min_score {
        return Err(GfError::EmptyScoreRange { min: min_score, max: max_score });
    }

    // Enzyme neighboring-AA adjustment.
    let final_dist: ScoreDist = if let Some(enzyme) = graph.enzyme {
        if !enzyme.residues().is_empty() {
            let credit  = aa_set.neighboring_aa_cleavage_credit();
            let penalty = aa_set.neighboring_aa_cleavage_penalty();
            let prob_clv = aa_set.prob_cleavage_sites(enzyme) as f64;

            let mut fd = ScoreDist::new(min_score + penalty, max_score + credit, false, true);
            fd.add_prob_dist(&sink_dist, credit, prob_clv);
            fd.add_prob_dist(&sink_dist, penalty, 1.0 - prob_clv);
            fd
        } else {
            sink_dist
        }
    } else {
        sink_dist
    };

    let final_min = final_dist.min_score();
    let final_max = final_dist.max_score();

    Ok(GeneratingFunction {
        score_dists: vec![final_dist],
        score_bound: ScoreBound::new(final_min, final_max),
        node_dists: node_dists_snapshot,
    })
}

// -----------------------------------------------------------------------
// Tests (uniform-prior DP — renamed from compute to compute_uniform)
// -----------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Trivial increment_score: every (mass, aa) gives score 0.
    /// Result: full probability mass at score 0 for every reachable mass.
    fn zero_inc(_mass: i32, _aa: u8) -> i32 { 0 }

    /// All amino acids have nominal mass 1, and there are 2 AAs.
    /// At mass M, the only reachable score is 0 with prob 1.0.
    fn aa_masses_uniform_one() -> Vec<i32> {
        vec![1, 1]  // 2 AAs each with nominal mass 1
    }

    #[test]
    fn empty_peptide_at_mass_zero() {
        let aa_masses = aa_masses_uniform_one();
        let gf = GeneratingFunction::compute_uniform(
            10,                    // max_mass
            ScoreBound::new(0, 5), // score range [0, 5)
            &aa_masses,
            zero_inc,
        );
        // At mass 0: only score 0 has probability, equal to 1.0.
        let d0 = gf.score_dist_at(0).expect("dist at mass 0");
        assert!((d0.get_probability(0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn dist_at_mass_one_with_zero_increment() {
        let aa_masses = aa_masses_uniform_one();
        let gf = GeneratingFunction::compute_uniform(
            5,
            ScoreBound::new(0, 5),
            &aa_masses,
            zero_inc,
        );
        // At mass 1, both AAs (each mass 1, prior 1/2) contribute. Each adds
        // (prob_at_mass_0 / 2) at score 0+0=0. So total prob at score 0 = 1.0.
        let d1 = gf.score_dist_at(1).expect("dist at mass 1");
        assert!((d1.get_probability(0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn nonzero_increment_shifts_score() {
        // Increment = 1 always. At mass 1: prob mass moves from score 0 (mass 0)
        // to score 1 (mass 1).
        let aa_masses = aa_masses_uniform_one();
        let gf = GeneratingFunction::compute_uniform(
            5,
            ScoreBound::new(0, 5),
            &aa_masses,
            |_m, _a| 1,
        );
        let d1 = gf.score_dist_at(1).expect("dist at mass 1");
        assert!((d1.get_probability(1) - 1.0).abs() < 1e-12);
        assert!(d1.get_probability(0).abs() < 1e-12);
        // At mass 2, increment +1 again: prob shifts to score 2.
        let d2 = gf.score_dist_at(2).expect("dist at mass 2");
        assert!((d2.get_probability(2) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn unreachable_mass_has_zero_prob() {
        // AA masses 2 and 3; mass 1 is unreachable.
        let gf = GeneratingFunction::compute_uniform(
            5,
            ScoreBound::new(0, 5),
            &[2, 3],
            zero_inc,
        );
        let d1 = gf.score_dist_at(1).expect("dist at mass 1 exists (zero)");
        // Total prob at mass 1 should be 0 (can't reach with AA masses 2 or 3).
        assert!(d1.get_probability(0).abs() < 1e-12);
    }

    #[test]
    fn two_aa_with_different_increments() {
        // AAs of mass 1 each. AA[0] gives +0 score, AA[1] gives +1 score.
        // At mass 1: prob 0.5 at score 0 (from AA[0]), prob 0.5 at score 1 (from AA[1]).
        let inc = |_m: i32, aa: u8| if aa == 0 { 0 } else { 1 };
        let gf = GeneratingFunction::compute_uniform(
            3,
            ScoreBound::new(0, 5),
            &[1, 1],
            inc,
        );
        let d1 = gf.score_dist_at(1).expect("dist at 1");
        assert!((d1.get_probability(0) - 0.5).abs() < 1e-12);
        assert!((d1.get_probability(1) - 0.5).abs() < 1e-12);
    }

    #[test]
    fn spectral_probability_at_target_mass() {
        // AA[0] = +1 always, AA[1] = -1 always. At mass 5, distribution
        // is binomial-like over scores -5..+5.
        let inc = |_m: i32, aa: u8| if aa == 0 { 1 } else { -1 };
        let gf = GeneratingFunction::compute_uniform(
            5,
            ScoreBound::new(-10, 10),
            &[1, 1],
            inc,
        );
        let d5 = gf.score_dist_at(5).expect("dist at 5");
        // Sum of all probabilities at this mass should be ~1.0
        let mut total = 0.0;
        for s in -10..10 {
            total += d5.get_probability(s);
        }
        assert!((total - 1.0).abs() < 1e-9, "total prob = {}", total);
    }

}
