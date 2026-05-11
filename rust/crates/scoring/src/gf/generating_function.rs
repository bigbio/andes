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
//! `compute` (and `with_score_threshold`) mirror Java's
//! `PrimitiveGeneratingFunction.computeGeneratingFunction()` and
//! `setUpScoreThreshold()`. They operate on a pre-built `PrimitiveAaGraph`
//! and produce a single final `ScoreDist` (plus enzyme adjustment).

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

    /// Compute the GF over a precomputed primitive graph. Mirrors Java
    /// `PrimitiveGeneratingFunction.computeGeneratingFunction()`.
    pub fn compute(graph: &PrimitiveAaGraph, aa_set: &AminoAcidSet) -> Result<Self, GfError> {
        compute_inner(graph, aa_set, None, false)
    }

    /// Pre-pass: prune nodes whose maximum possible final score is below
    /// `score_threshold`. Mirrors Java `setUpScoreThreshold`; computes
    /// `min_score_by_node` and uses it to skip irrelevant DP work.
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

    /// Cumulative tail probability `P(random_score >= score)`. Mirrors Java
    /// `PrimitiveGeneratingFunction.getSpectralProbability(score)`.
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

/// Translate Java `setUpScoreThreshold` (lines 36-87).
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

    // Adjust threshold for enzyme neighboring-AA credit (Java line 48-50).
    let adjusted_score = if graph.enzyme.is_some() {
        score_threshold - aa_set.neighboring_aa_cleavage_credit()
    } else {
        score_threshold
    };

    let mut min_score_by_node = vec![i32::MAX; node_count];
    min_score_by_node[sink_idx] = adjusted_score;

    // Propagate from sink backward through sink's own incoming edges
    // (Java lines 58-66).
    for e in graph.edge_offset[sink_idx]..graph.edge_offset[sink_idx + 1] {
        let prev_mass = graph.edge_prev_node[e];
        if let Some(prev_idx) = graph.node_index_for_mass(prev_mass) {
            let new_min = adjusted_score.saturating_sub(graph.edge_score[e]);
            if new_min < min_score_by_node[prev_idx] {
                min_score_by_node[prev_idx] = new_min;
            }
        }
    }

    // Walk nodes in reverse order (Java line 68: `for (int ni = nodeCount-1; ni >= 0; ni--)`).
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

/// Core DP. Translates Java `computeGeneratingFunction` (lines 89-205).
///
/// `retain_node_dists` is a diagnostic-only flag: when `true`, the per-node
/// `ScoreDist` buffer is cloned into `GeneratingFunction.node_dists` so the
/// caller can inspect intermediate distributions via `iter_node_dists`. When
/// `false` (the production path) the DP behavior is unchanged.
fn compute_inner(
    graph: &PrimitiveAaGraph,
    aa_set: &AminoAcidSet,
    min_score_by_node: Option<Vec<i32>>,
    retain_node_dists: bool,
) -> Result<GeneratingFunction, GfError> {
    let node_count = graph.node_count;
    let source_idx = graph.source_node_idx;
    let sink_idx = graph.sink_node_idx;

    // dist_by_node[ni] = Some(dist) once computed; None = not yet reachable.
    let mut dist_by_node: Vec<Option<ScoreDist>> = vec![None; node_count];

    // Debug-only counter: tracks how many nodes were skipped due to the
    // score-range guard (|score| > 10000). Fires only in debug builds;
    // release builds compile this out entirely (no perf regression).
    #[cfg(debug_assertions)]
    let mut score_range_overflow_count: u32 = 0;

    // Source has full probability at score 0 (Java lines 101-103).
    let mut source_dist = ScoreDist::new(0, 1, false, true);
    source_dist.set_prob(0, 1.0);
    dist_by_node[source_idx] = Some(source_dist);

    // Scratch buffer for valid edge indices (Java lines 106-111).
    let max_edges_per_node = (0..node_count)
        .map(|ni| graph.edge_offset[ni + 1] - graph.edge_offset[ni])
        .max()
        .unwrap_or(0);
    let mut valid_edges: Vec<usize> = Vec::with_capacity(max_edges_per_node);

    // Forward DP over nodes in index order (Java lines 114-176).
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

        // Determine initial curMinScore (Java lines 124-129).
        let mut cur_min_score: i32 = match min_score_by_node {
            Some(ref msbn) => msbn[ni],
            None => i32::MAX,
        };
        let mut cur_max_score: i32 = i32::MIN;

        valid_edges.clear();

        // Scan incoming edges (Java lines 132-149).
        for e in graph.edge_offset[ni]..graph.edge_offset[ni + 1] {
            let prev_mass = graph.edge_prev_node[e];
            let prev_idx = match graph.node_index_for_mass(prev_mass) {
                Some(idx) => idx,
                None => continue,
            };
            let prev_dist = match dist_by_node[prev_idx].as_ref() {
                Some(d) => d,
                None => continue,
            };

            let combined_score = cur_node_score + graph.edge_score[e];
            let possible_max = prev_dist.max_score() + combined_score;
            if possible_max > cur_max_score {
                cur_max_score = possible_max;
            }

            // Only update min from predecessor when NOT using threshold pre-pass.
            if min_score_by_node.is_none() {
                let possible_min = prev_dist.min_score() + combined_score;
                if possible_min < cur_min_score {
                    cur_min_score = possible_min;
                }
            }

            valid_edges.push(e);
        }

        // Skip degenerate or out-of-bound ranges (Java lines 152-158).
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

        // Allocate and fill cur_dist (Java lines 160-169).
        let mut cur_dist = ScoreDist::new(cur_min_score, cur_max_score, false, true);

        for &e in &valid_edges {
            let prev_mass = graph.edge_prev_node[e];
            // Safety: we already verified these are valid above.
            let prev_idx = graph.node_index_for_mass(prev_mass).unwrap();
            let prev_dist = dist_by_node[prev_idx].as_ref().unwrap();
            let combined_score = cur_node_score + graph.edge_score[e];
            cur_dist.add_prob_dist(prev_dist, combined_score, graph.edge_prob[e] as f64);
        }

        // Underflow guard at max_score - 1 (Java lines 171-173).
        let guard_score = cur_max_score - 1;
        if cur_dist.get_probability(guard_score) == 0.0 {
            // Mirrors Java's `Float.MIN_VALUE` (smallest positive denormal f32 ~1.4e-45),
            // NOT `f32::MIN_POSITIVE` (smallest positive normal ~1.18e-38). The Java GF
            // uses denormal min as the underflow floor; using normal min instead biases
            // SpecEValues low by ~7 OOM for short peptides.
            cur_dist.set_prob(guard_score, f32::from_bits(1) as f64);
        }

        dist_by_node[ni] = Some(cur_dist);
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

    // Diagnostic-only: if requested, snapshot per-node dists before the
    // `.take()` below would consume the sink. Production path leaves this as
    // `None`, identical to prior behavior.
    let node_dists_snapshot: Option<Vec<(usize, i32, ScoreDist)>> = if retain_node_dists {
        let mut snap: Vec<(usize, i32, ScoreDist)> = Vec::new();
        for ni in 0..node_count {
            if let Some(d) = dist_by_node[ni].as_ref() {
                snap.push((ni, graph.node_mass(ni), d.clone()));
            }
        }
        Some(snap)
    } else {
        None
    };

    // Extract sink distribution (Java lines 179-185).
    let sink_dist = dist_by_node[sink_idx]
        .take()
        .ok_or(GfError::SinkUnreachable)?;

    let min_score = sink_dist.min_score();
    let max_score = sink_dist.max_score();

    if max_score <= min_score {
        return Err(GfError::EmptyScoreRange { min: min_score, max: max_score });
    }

    // Enzyme neighboring-AA adjustment (Java lines 188-200).
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
