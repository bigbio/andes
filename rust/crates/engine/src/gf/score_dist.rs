//! `ScoreBound` + `ScoreDist` data structures for the GF DP.
//! Mirrors Java `edu.ucsd.msjava.msgf.{ScoreBound, ScoreDist}`.
//!
//! `ScoreDist` stores per-score arrays of probabilities and/or counts
//! over an integer score range `[min_score, max_score)`. Index = score - min_score.

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

#[derive(Debug, Clone)]
pub struct ScoreDist {
    bound: ScoreBound,
    num_distribution: Option<Vec<f64>>,
    prob_distribution: Option<Vec<f64>>,
}

impl ScoreDist {
    pub fn new(min_score: i32, max_score: i32, calc_number: bool, calc_prob: bool) -> Self {
        let range = (max_score - min_score) as usize;
        Self {
            bound: ScoreBound::new(min_score, max_score),
            num_distribution: if calc_number { Some(vec![0.0; range]) } else { None },
            prob_distribution: if calc_prob { Some(vec![0.0; range]) } else { None },
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
}
