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
}

impl PartialEq for PsmMatch {
    fn eq(&self, other: &Self) -> bool { self.score == other.score }
}

impl Eq for PsmMatch {}

impl PartialOrd for PsmMatch {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for PsmMatch {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.score.partial_cmp(&other.score).unwrap_or(std::cmp::Ordering::Equal)
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

    pub fn push(&mut self, m: PsmMatch) {
        if self.heap.len() < self.capacity as usize {
            self.heap.push(Reverse(m));
        } else if let Some(Reverse(top)) = self.heap.peek() {
            if m.score > top.score {
                self.heap.pop();
                self.heap.push(Reverse(m));
            }
        }
    }

    pub fn len(&self) -> usize { self.heap.len() }
    pub fn is_empty(&self) -> bool { self.heap.is_empty() }

    /// Drain into a Vec sorted by score descending.
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
        }
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
        assert_eq!(sorted.iter().map(|m| m.score).collect::<Vec<_>>(),
                   vec![3.0, 2.0, 1.0]);
    }

    #[test]
    fn queue_at_capacity_keeps_top_n_by_score() {
        let mut q = TopNQueue::new(3);
        for s in [1.0, 5.0, 2.0, 4.0, 3.0] { q.push(make_match(0, s)); }
        assert_eq!(q.len(), 3);
        let sorted = q.into_sorted_vec();
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
    }
}
