//! `ScoreBound` + `ScoreDist` data structures for the GF DP.
//! Mirrors Java `edu.ucsd.msjava.msgf.{ScoreBound, ScoreDist}`.
//!
//! `ScoreDist` stores per-score arrays of probabilities and/or counts
//! over an integer score range `[min_score, max_score)`. Index = score - min_score.
//!
//! ## Buffer pooling
//!
//! `ScoreDist::new` pulls its `Vec<f64>` storage from a thread-local pool
//! and `Drop` returns the buffer to the pool. This eliminates ~19M Vec
//! allocations per PXD001819 search (~50 active nodes × 10 mass bins ×
//! 38k spectra). The pool is per-thread (no contention) and capped at
//! `MAX_POOL_SIZE` buffers to bound memory.

use std::cell::RefCell;

const MAX_POOL_SIZE: usize = 256;

thread_local! {
    /// Free list of `Vec<f64>` buffers reusable by `ScoreDist`. Borrowed
    /// only briefly during `new` / `Drop` — no long-lived borrows.
    static F64_POOL: RefCell<Vec<Vec<f64>>> = const { RefCell::new(Vec::new()) };
}

/// Take a buffer from the pool, or allocate a fresh one. Returned Vec
/// is empty (length 0) but has at least `min_capacity` capacity.
fn pool_take(min_capacity: usize) -> Vec<f64> {
    F64_POOL.with(|p| {
        let mut p = p.borrow_mut();
        // Walk the pool from the back (most-recently-returned first, hot in
        // L1/L2 cache) for a buffer with sufficient capacity.
        if let Some(idx) = p.iter().rposition(|v| v.capacity() >= min_capacity) {
            let mut v = p.swap_remove(idx);
            v.clear();
            v
        } else {
            Vec::with_capacity(min_capacity)
        }
    })
}

/// Return a buffer to the pool. Pool is capped at `MAX_POOL_SIZE` to bound
/// per-thread memory; surplus buffers are dropped.
fn pool_return(mut v: Vec<f64>) {
    if v.capacity() == 0 {
        return;
    }
    v.clear();
    F64_POOL.with(|p| {
        let mut p = p.borrow_mut();
        if p.len() < MAX_POOL_SIZE {
            p.push(v);
        }
    });
}

#[derive(Debug, Clone, Copy)]
pub struct ScoreBound {
    /// inclusive
    min_score: i32,
    /// exclusive
    max_score: i32,
}

impl ScoreBound {
    pub fn new(min_score: i32, max_score: i32) -> Self {
        Self { min_score, max_score }
    }

    pub fn min_score(&self) -> i32 { self.min_score }
    pub fn max_score(&self) -> i32 { self.max_score }
    pub fn range(&self) -> i32 { self.max_score - self.min_score }

    pub fn set_min_score(&mut self, v: i32) { self.min_score = v; }
    pub fn set_max_score(&mut self, v: i32) { self.max_score = v; }
}

#[derive(Debug)]
pub struct ScoreDist {
    bound: ScoreBound,
    num_distribution: Option<Vec<f64>>,
    prob_distribution: Option<Vec<f64>>,
}

/// Manual `Clone` (replaces the derive that came with the original struct).
/// We can't pool through a clone because both copies would race on `Drop`;
/// cloned buffers are always freshly allocated and freshly zeroed via the
/// underlying `Vec::clone` (no pool involvement).
impl Clone for ScoreDist {
    fn clone(&self) -> Self {
        Self {
            bound: self.bound,
            num_distribution: self.num_distribution.clone(),
            prob_distribution: self.prob_distribution.clone(),
        }
    }
}

impl Drop for ScoreDist {
    fn drop(&mut self) {
        if let Some(v) = self.num_distribution.take() {
            pool_return(v);
        }
        if let Some(v) = self.prob_distribution.take() {
            pool_return(v);
        }
    }
}

impl ScoreDist {
    pub fn new(min_score: i32, max_score: i32, calc_number: bool, calc_prob: bool) -> Self {
        let range = (max_score - min_score) as usize;
        let make = || {
            let mut v = pool_take(range);
            v.resize(range, 0.0);
            v
        };
        Self {
            bound: ScoreBound::new(min_score, max_score),
            num_distribution: if calc_number { Some(make()) } else { None },
            prob_distribution: if calc_prob { Some(make()) } else { None },
        }
    }

    pub fn bound(&self) -> ScoreBound { self.bound }
    pub fn min_score(&self) -> i32 { self.bound.min_score }
    pub fn max_score(&self) -> i32 { self.bound.max_score }

    pub fn is_prob_set(&self) -> bool { self.prob_distribution.is_some() }
    pub fn is_num_set(&self) -> bool { self.num_distribution.is_some() }

    pub fn set_prob(&mut self, score: i32, prob: f64) {
        let idx = (score - self.bound.min_score) as usize;
        if let Some(p) = self.prob_distribution.as_mut() {
            p[idx] = prob;
        }
    }

    pub fn add_prob(&mut self, score: i32, prob: f64) {
        let idx = (score - self.bound.min_score) as usize;
        if let Some(p) = self.prob_distribution.as_mut() {
            p[idx] += prob;
        }
    }

    pub fn set_number(&mut self, score: i32, n: f64) {
        let idx = (score - self.bound.min_score) as usize;
        if let Some(p) = self.num_distribution.as_mut() {
            p[idx] = n;
        }
    }

    pub fn add_number(&mut self, score: i32, n: f64) {
        let idx = (score - self.bound.min_score) as usize;
        if let Some(p) = self.num_distribution.as_mut() {
            p[idx] += n;
        }
    }

    /// Java: `getProbability` returns `probDistribution[max(0, score - minScore)]`.
    /// A score below minScore returns the entry at index 0; above maxScore is
    /// undefined behavior in Java (would index out of bounds). We mirror Java
    /// for in-range and clamp-to-zero for below-range; above-range is caller's
    /// responsibility (panics if out of bounds).
    pub fn get_probability(&self, score: i32) -> f64 {
        let p = self.prob_distribution.as_ref().expect("prob distribution not allocated");
        let idx = if score >= self.bound.min_score {
            (score - self.bound.min_score) as usize
        } else {
            0
        };
        p[idx]
    }

    pub fn get_number_recs(&self, score: i32) -> f64 {
        let n = self.num_distribution.as_ref().expect("num distribution not allocated");
        let idx = if score >= self.bound.min_score {
            (score - self.bound.min_score) as usize
        } else {
            0
        };
        n[idx]
    }

    /// Cumulative tail probability `P(X >= score)`, clamped to 1.0.
    pub fn get_spectral_probability(&self, score: i32) -> f64 {
        let p = self.prob_distribution.as_ref().expect("prob distribution not allocated");
        let min_index = if score >= self.bound.min_score {
            (score - self.bound.min_score) as usize
        } else {
            0
        };
        let sum: f64 = p[min_index..].iter().sum();
        sum.min(1.0)
    }

    /// Mirror Java's `ScoreDist.addProbDist(other, scoreDiff, aaProb)`:
    /// for each `t` in `other`'s score range, accumulate
    /// `other.prob[t] * aa_prob` into `self.prob[t + score_diff]`,
    /// clipping the destination to `self`'s range.
    pub fn add_prob_dist(&mut self, other: &ScoreDist, score_diff: i32, aa_prob: f64) {
        let other_p = match other.prob_distribution.as_ref() {
            Some(p) => p,
            None => return,
        };
        let self_p = match self.prob_distribution.as_mut() {
            Some(p) => p,
            None => return,
        };
        let other_min = other.bound.min_score;
        let other_max = other.bound.max_score;
        let self_min = self.bound.min_score;
        let self_max = self.bound.max_score;
        let t_start = other_min.max(self_min - score_diff);
        let t_end = other_max.min(self_max - score_diff);
        for t in t_start..t_end {
            let src_idx = (t - other_min) as usize;
            let dst_idx = (t + score_diff - self_min) as usize;
            self_p[dst_idx] += other_p[src_idx] * aa_prob;
        }
    }

    /// Mirror Java's `ScoreDist.addNumDist(other, scoreDiff, coeff)`.
    pub fn add_num_dist(&mut self, other: &ScoreDist, score_diff: i32, coeff: f64) {
        let other_n = match other.num_distribution.as_ref() {
            Some(n) => n,
            None => return,
        };
        let self_n = match self.num_distribution.as_mut() {
            Some(n) => n,
            None => return,
        };
        let other_min = other.bound.min_score;
        let other_max = other.bound.max_score;
        let self_min = self.bound.min_score;
        let self_max = self.bound.max_score;
        let t_start = other_min.max(self_min - score_diff);
        let t_end = other_max.min(self_max - score_diff);
        for t in t_start..t_end {
            let src_idx = (t - other_min) as usize;
            let dst_idx = (t + score_diff - self_min) as usize;
            self_n[dst_idx] += other_n[src_idx] * coeff;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn score_bound_range() {
        let b = ScoreBound::new(-3, 7);
        assert_eq!(b.min_score(), -3);
        assert_eq!(b.max_score(), 7);
        assert_eq!(b.range(), 10);
    }

    #[test]
    fn score_dist_set_get_prob() {
        let mut d = ScoreDist::new(-2, 5, false, true);
        d.set_prob(0, 0.5);
        d.set_prob(-2, 0.1);
        d.set_prob(4, 0.2);
        assert_eq!(d.get_probability(0), 0.5);
        assert_eq!(d.get_probability(-2), 0.1);
        assert_eq!(d.get_probability(4), 0.2);
    }

    #[test]
    fn score_dist_add_prob_accumulates() {
        let mut d = ScoreDist::new(0, 5, false, true);
        d.set_prob(2, 0.1);
        d.add_prob(2, 0.3);
        assert!((d.get_probability(2) - 0.4).abs() < 1e-9);
    }

    #[test]
    fn score_dist_set_get_number() {
        let mut d = ScoreDist::new(0, 5, true, false);
        d.set_number(3, 100.0);
        d.add_number(3, 50.0);
        assert!((d.get_number_recs(3) - 150.0).abs() < 1e-9);
    }

    #[test]
    fn is_prob_set_and_is_num_set() {
        let only_prob = ScoreDist::new(0, 5, false, true);
        assert!(only_prob.is_prob_set());
        assert!(!only_prob.is_num_set());

        let only_num = ScoreDist::new(0, 5, true, false);
        assert!(!only_num.is_prob_set());
        assert!(only_num.is_num_set());

        let both = ScoreDist::new(0, 5, true, true);
        assert!(both.is_prob_set());
        assert!(both.is_num_set());
    }

    #[test]
    fn score_below_min_clamped_to_min_index() {
        let mut d = ScoreDist::new(0, 5, false, true);
        d.set_prob(0, 0.5);
        // Java: getProbability returns probDistribution[max(0, score - minScore)],
        // so a score below minScore returns the entry at index 0.
        assert_eq!(d.get_probability(-10), 0.5);
    }

    #[test]
    fn spectral_probability_is_cumulative_sum() {
        let mut d = ScoreDist::new(0, 5, false, true);
        d.set_prob(0, 0.1);
        d.set_prob(1, 0.2);
        d.set_prob(2, 0.3);
        d.set_prob(3, 0.05);
        d.set_prob(4, 0.05);
        // Sum from score=2 onward = 0.3 + 0.05 + 0.05 = 0.4
        assert!((d.get_spectral_probability(2) - 0.4).abs() < 1e-9);
        // Sum from score=0 onward = 0.7
        assert!((d.get_spectral_probability(0) - 0.7).abs() < 1e-9);
    }

    #[test]
    fn spectral_probability_clamped_to_one() {
        // Even if the sum exceeds 1.0 (numerical overshoot), output clamped.
        let mut d = ScoreDist::new(0, 5, false, true);
        for s in 0..5 { d.set_prob(s, 0.5); }  // sum = 2.5
        assert!((d.get_spectral_probability(0) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn spectral_probability_below_min_uses_index_zero() {
        let mut d = ScoreDist::new(2, 5, false, true);
        d.set_prob(2, 0.1);
        d.set_prob(3, 0.2);
        d.set_prob(4, 0.3);
        // score < minScore: minIndex = 0, sum from there = 0.1 + 0.2 + 0.3 = 0.6
        assert!((d.get_spectral_probability(-100) - 0.6).abs() < 1e-9);
    }

    #[test]
    fn add_prob_dist_offset_zero_scalar_one() {
        // self range [0, 5), other range [0, 5). After add_prob_dist(other, 0, 1.0)
        // each self[s] += other[s].
        let mut a = ScoreDist::new(0, 5, false, true);
        let mut b = ScoreDist::new(0, 5, false, true);
        for s in 0..5 { b.set_prob(s, 0.1 * (s + 1) as f64); }
        a.add_prob_dist(&b, 0, 1.0);
        for s in 0..5 {
            assert!((a.get_probability(s) - 0.1 * (s + 1) as f64).abs() < 1e-12);
        }
    }

    #[test]
    fn add_prob_dist_with_score_offset() {
        // self [0, 10), other [0, 5). add(other, +3, 1.0) shifts other's scores
        // by +3: self[3..8] += other[0..5].
        let mut a = ScoreDist::new(0, 10, false, true);
        let mut b = ScoreDist::new(0, 5, false, true);
        for s in 0..5 { b.set_prob(s, 0.2); }
        a.add_prob_dist(&b, 3, 1.0);
        for s in 0..3 { assert_eq!(a.get_probability(s), 0.0); }
        for s in 3..8 { assert!((a.get_probability(s) - 0.2).abs() < 1e-12); }
        for s in 8..10 { assert_eq!(a.get_probability(s), 0.0); }
    }

    #[test]
    fn add_prob_dist_with_negative_offset() {
        // self [-3, 5), other [0, 5). add(other, -2, 1.0) shifts down by 2.
        let mut a = ScoreDist::new(-3, 5, false, true);
        let mut b = ScoreDist::new(0, 5, false, true);
        for s in 0..5 { b.set_prob(s, 0.1); }
        a.add_prob_dist(&b, -2, 1.0);
        // other[0]→self[-2], other[4]→self[2]; self[-3] and self[3..5) untouched.
        assert_eq!(a.get_probability(-3), 0.0);
        for s in -2..3 { assert!((a.get_probability(s) - 0.1).abs() < 1e-12); }
        for s in 3..5 { assert_eq!(a.get_probability(s), 0.0); }
    }

    #[test]
    fn add_prob_dist_clips_to_self_range() {
        // self [0, 3), other [0, 5). add(other, 0, 1.0) only fills self[0..3].
        let mut a = ScoreDist::new(0, 3, false, true);
        let mut b = ScoreDist::new(0, 5, false, true);
        for s in 0..5 { b.set_prob(s, 0.2); }
        a.add_prob_dist(&b, 0, 1.0);
        for s in 0..3 { assert!((a.get_probability(s) - 0.2).abs() < 1e-12); }
    }

    #[test]
    fn add_prob_dist_scales_by_aa_prob() {
        let mut a = ScoreDist::new(0, 5, false, true);
        let mut b = ScoreDist::new(0, 5, false, true);
        for s in 0..5 { b.set_prob(s, 0.1); }
        a.add_prob_dist(&b, 0, 0.5);
        for s in 0..5 { assert!((a.get_probability(s) - 0.05).abs() < 1e-12); }
    }

    #[test]
    fn add_num_dist_with_coefficient() {
        let mut a = ScoreDist::new(0, 5, true, false);
        let mut b = ScoreDist::new(0, 5, true, false);
        for s in 0..5 { b.set_number(s, 2.0); }
        a.add_num_dist(&b, 0, 3.0);
        for s in 0..5 { assert!((a.get_number_recs(s) - 6.0).abs() < 1e-12); }
    }
}
