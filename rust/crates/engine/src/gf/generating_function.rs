//! Generating-function DP: computes the score distribution
//! `P(score | random peptide of given nominal mass)`.
//!
//! Phase 6 Task 2: bounded DP machinery taking a generic increment-score
//! callback. Phase 6 Task 3+ will bind the callback to RankScorer.
//!
//! Algorithm (simplified vs Java's GeneratingFunction.java):
//!
//! For each nominal_mass M in 0..=max_mass:
//!   if M == 0:
//!     score_dist[0].add_prob(0, 1.0)
//!   else:
//!     For each amino-acid index aa_idx with nominal_mass aa_mass:
//!       if M - aa_mass < 0: continue
//!       inc = increment_score(M, aa_idx)
//!       For each score s in score_dist[M - aa_mass]:
//!         p = score_dist[M - aa_mass].get_probability(s)
//!         score_dist[M].add_prob(s + inc, p / num_aas)
//!
//! Java reference: edu.ucsd.msjava.msgf.GeneratingFunction.

use crate::gf::score_dist::{ScoreBound, ScoreDist};

#[derive(thiserror::Error, Debug)]
pub enum GfError {
    #[error("score range is empty: min_score {min} >= max_score {max}")]
    EmptyScoreRange { min: i32, max: i32 },
    #[error("aa_masses is empty")]
    NoAminoAcids,
}

#[derive(Debug, Clone)]
pub struct GeneratingFunction {
    /// One ScoreDist per nominal mass in 0..=max_mass.
    score_dists: Vec<ScoreDist>,
    score_bound: ScoreBound,
}

impl GeneratingFunction {
    /// Compute the generating function up to `max_mass`. The
    /// `increment_score` callback returns the score added when the
    /// peptide is extended by amino acid `aa_idx` (an index into
    /// `aa_masses`) at mass position `mass`.
    ///
    /// Probability prior over amino acids: uniform `1 / aa_masses.len()`.
    pub fn compute<F>(
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
    pub fn spectral_probability(&self, mass: i32, score: i32) -> Option<f64> {
        self.score_dist_at(mass).map(|d| d.get_spectral_probability(score))
    }
}

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
        let gf = GeneratingFunction::compute(
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
        let gf = GeneratingFunction::compute(
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
        let gf = GeneratingFunction::compute(
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
        let gf = GeneratingFunction::compute(
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
        let gf = GeneratingFunction::compute(
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
        let gf = GeneratingFunction::compute(
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
