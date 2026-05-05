//! PSM (peptide-spectrum match) data + top-N ranking queue.

use std::cmp::Reverse;
use std::collections::BinaryHeap;

use crate::candidate_gen::Candidate;

#[derive(Debug, Clone)]
pub struct PsmMatch {
    pub spectrum_idx: usize,
    pub candidate: Candidate,
    pub charge_used: u8,
    /// Signed: positive when peptide mass exceeds spectrum's implied mass.
    pub mass_error_ppm: f64,
    /// Higher is better. Phase 5 fills with real spectral-similarity score.
    /// Phase 4e MVP uses negative |mass_error_ppm| as a placeholder.
    pub score: f32,
    /// Phase 6 SpecEValue: lower is better. Default 1.0 = "not yet computed"
    /// / "no signal". Set by `compute_spec_e_values_for_spectrum` after the
    /// per-candidate scoring loop.
    pub spec_e_value: f64,
}

impl PartialEq for PsmMatch {
    fn eq(&self, other: &Self) -> bool {
        self.spec_e_value == other.spec_e_value && self.score == other.score
    }
}

impl Eq for PsmMatch {}

impl PartialOrd for PsmMatch {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

/// Primary: `spec_e_value` ascending (lower = better).
/// Secondary: `score` descending (higher = better).
///
/// This ordering is used by `TopNQueue`'s min-heap (via `Reverse<PsmMatch>`):
/// the heap's "minimum" element is the one with the *largest* spec_e_value
/// (worst), so `push` evicts it when over capacity.
impl Ord for PsmMatch {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // spec_e_value ascending: `other` compared to `self` gives ascending sort
        // when embedded in the `Reverse` wrapper of `BinaryHeap`.
        // But here we define natural order "better = greater" so that the heap
        // min sits at the top for easy eviction.
        // "Better" = smaller spec_e_value; then larger score.
        match other.spec_e_value.partial_cmp(&self.spec_e_value) {
            Some(std::cmp::Ordering::Equal) | None => {
                // Secondary: score descending (larger score = better)
                self.score.partial_cmp(&other.score).unwrap_or(std::cmp::Ordering::Equal)
            }
            Some(ord) => ord,
        }
    }
}

#[derive(Debug, Clone)]
pub struct TopNQueue {
    capacity: u32,
    /// Min-heap (via Reverse): smallest score sits at top, easy to pop
    /// when over capacity.
    heap: BinaryHeap<Reverse<PsmMatch>>,
}

impl TopNQueue {
    pub fn new(capacity: u32) -> Self {
        Self { capacity, heap: BinaryHeap::with_capacity(capacity as usize) }
    }

    /// Insert a PSM. The queue keeps the `capacity` *best* PSMs.
    ///
    /// "Best" = smallest `spec_e_value` first (then largest `score` for ties).
    /// The min-heap (via `Reverse<PsmMatch>`) puts the *worst* PSM at the top
    /// so it can be evicted when over capacity.
    ///
    /// Before Phase 6 computes spec_e_value, all PSMs have `spec_e_value = 1.0`
    /// and the secondary `score` key governs eviction — same behaviour as before.
    pub fn push(&mut self, m: PsmMatch) {
        if self.heap.len() < self.capacity as usize {
            self.heap.push(Reverse(m));
        } else if let Some(Reverse(top)) = self.heap.peek() {
            // `m > top` in natural ordering means m is better.
            if m.cmp(top) == std::cmp::Ordering::Greater {
                self.heap.pop();
                self.heap.push(Reverse(m));
            }
        }
    }

    pub fn len(&self) -> usize { self.heap.len() }
    pub fn is_empty(&self) -> bool { self.heap.is_empty() }

    /// Iterate over all PSMs in the queue (order not guaranteed).
    pub fn iter_psms(&self) -> impl Iterator<Item = &PsmMatch> {
        self.heap.iter().map(|Reverse(m)| m)
    }

    /// Apply `f` to each PSM to compute its `spec_e_value`, then rebuild
    /// the heap so the ordering invariant holds.
    ///
    /// Draining + re-inserting is O(N log N) — cheap for small N (top-10).
    pub fn update_spec_e_values<F: Fn(&PsmMatch) -> f64>(&mut self, f: F) {
        let mut psms: Vec<PsmMatch> = self.heap.drain().map(|Reverse(m)| m).collect();
        for psm in &mut psms {
            psm.spec_e_value = f(psm);
        }
        for psm in psms {
            self.heap.push(Reverse(psm));
        }
    }

    /// Drain into a Vec sorted best-first (smallest spec_e_value, then largest score).
    pub fn into_sorted_vec(self) -> Vec<PsmMatch> {
        let mut v: Vec<PsmMatch> = self.heap.into_iter().map(|Reverse(m)| m).collect();
        v.sort_by(|a, b| b.cmp(a));
        v
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::amino_acid::AminoAcid;
    use crate::peptide::Peptide;

    fn make_match(spectrum_idx: usize, score: f32) -> PsmMatch {
        let aa = AminoAcid::standard(b'A').unwrap();
        let peptide = Peptide::new(vec![aa], b'_', b'-');
        PsmMatch {
            spectrum_idx,
            candidate: Candidate {
                peptide,
                protein_index: 0,
                start_offset_in_protein: 0,
                is_decoy: false,
            },
            charge_used: 2,
            mass_error_ppm: 0.0,
            score,
            spec_e_value: 1.0,  // default sentinel: "not yet computed"
        }
    }

    fn make_match_with_evalue(spectrum_idx: usize, score: f32, spec_e_value: f64) -> PsmMatch {
        let mut m = make_match(spectrum_idx, score);
        m.spec_e_value = spec_e_value;
        m
    }

    #[test]
    fn empty_queue() {
        let q = TopNQueue::new(5);
        assert!(q.is_empty());
        assert_eq!(q.len(), 0);
    }

    #[test]
    fn queue_below_capacity_keeps_everything() {
        let mut q = TopNQueue::new(5);
        for s in [1.0, 2.0, 3.0] { q.push(make_match(0, s)); }
        assert_eq!(q.len(), 3);
        let sorted = q.into_sorted_vec();
        // All spec_e_value = 1.0 (default) → secondary sort by score descending.
        assert_eq!(sorted.iter().map(|m| m.score).collect::<Vec<_>>(),
                   vec![3.0, 2.0, 1.0]);
    }

    #[test]
    fn queue_at_capacity_keeps_top_n_by_score() {
        let mut q = TopNQueue::new(3);
        for s in [1.0, 5.0, 2.0, 4.0, 3.0] { q.push(make_match(0, s)); }
        assert_eq!(q.len(), 3);
        let sorted = q.into_sorted_vec();
        // All spec_e_value = 1.0 → secondary score keeps top-3 by score.
        assert_eq!(sorted.iter().map(|m| m.score).collect::<Vec<_>>(),
                   vec![5.0, 4.0, 3.0]);
    }

    #[test]
    fn lower_score_dropped_when_full() {
        let mut q = TopNQueue::new(2);
        q.push(make_match(0, 5.0));
        q.push(make_match(0, 3.0));
        assert_eq!(q.len(), 2);
        q.push(make_match(0, 1.0));
        let sorted = q.into_sorted_vec();
        assert_eq!(sorted.iter().map(|m| m.score).collect::<Vec<_>>(),
                   vec![5.0, 3.0]);
    }

    #[test]
    fn psm_match_clones_correctly() {
        let m = make_match(7, 4.2);
        let cloned = m.clone();
        assert_eq!(cloned.spectrum_idx, 7);
        assert_eq!(cloned.score, 4.2);
        assert_eq!(cloned.spec_e_value, 1.0);
    }

    // -----------------------------------------------------------------------
    // Phase 6 / Task 8: SpecEValue ordering tests
    // -----------------------------------------------------------------------

    #[test]
    fn psm_match_orders_by_spec_e_value_ascending_then_score_descending() {
        // Lower spec_e_value means "better" → should sort before (greater in
        // natural Ord so the min-heap can evict the worst).
        let better = make_match_with_evalue(0, 5.0, 0.001);
        let worse  = make_match_with_evalue(0, 5.0, 0.5);
        // "better" is greater in natural order (because lower e-value wins).
        assert!(better > worse,
            "PSM with lower spec_e_value should be Ord-greater (better in the min-heap)");

        // Tie-break by score descending.
        let high_score = make_match_with_evalue(0, 10.0, 0.01);
        let low_score  = make_match_with_evalue(0, 3.0,  0.01);
        assert!(high_score > low_score,
            "when spec_e_value equal, higher score should be Ord-greater");
    }

    #[test]
    fn queue_keeps_best_spec_e_value_psms_when_full() {
        // Three PSMs with same score but different spec_e_values; capacity = 2.
        let mut q = TopNQueue::new(2);
        q.push(make_match_with_evalue(0, 5.0, 0.5));   // worst
        q.push(make_match_with_evalue(0, 5.0, 0.001)); // best
        assert_eq!(q.len(), 2);
        // Push a medium one; it should evict the worst (0.5).
        q.push(make_match_with_evalue(0, 5.0, 0.1));
        assert_eq!(q.len(), 2);
        let sorted = q.into_sorted_vec();
        // Should keep 0.001 and 0.1 (best two).
        let evalues: Vec<f64> = sorted.iter().map(|m| m.spec_e_value).collect();
        assert!(evalues.contains(&0.001), "best e-value 0.001 should be retained");
        assert!(evalues.contains(&0.1),   "medium e-value 0.1 should be retained");
        assert!(!evalues.contains(&0.5),  "worst e-value 0.5 should be evicted");
    }

    #[test]
    fn update_spec_e_values_applies_to_all_psms() {
        let mut q = TopNQueue::new(5);
        for s in [1.0_f32, 2.0, 3.0] {
            q.push(make_match(0, s));
        }
        // Set spec_e_value = 1.0 / score for each PSM.
        q.update_spec_e_values(|psm| 1.0 / psm.score as f64);
        let sorted = q.into_sorted_vec();
        // After update: score 3.0 → e=0.333, score 2.0 → e=0.5, score 1.0 → e=1.0.
        // Best e-value first.
        assert!((sorted[0].spec_e_value - 1.0 / 3.0).abs() < 1e-9);
        assert!((sorted[1].spec_e_value - 0.5).abs() < 1e-9);
        assert!((sorted[2].spec_e_value - 1.0).abs() < 1e-9);
    }

    #[test]
    fn iter_psms_yields_all_psms() {
        let mut q = TopNQueue::new(5);
        for s in [1.0_f32, 2.0, 3.0] { q.push(make_match(0, s)); }
        let scores: Vec<f32> = {
            let mut v: Vec<f32> = q.iter_psms().map(|m| m.score).collect();
            v.sort_by(|a, b| b.partial_cmp(a).unwrap());
            v
        };
        assert_eq!(scores, vec![3.0, 2.0, 1.0]);
    }
}
