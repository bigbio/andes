//! Streaming merger for `GeneratingFunction` distributions across
//! precursor-mass bins. Mirrors Java
//! `edu.ucsd.msjava.msgf.PrimitiveGeneratingFunctionGroup`.
//!
//! Math identity: `ScoreDist::add_prob_dist(other, 0, 1.0)` is a linear sum
//! over the probability arrays, so register-all-then-merge and
//! streaming-merge produce the same aggregate.

use crate::gf::generating_function::GeneratingFunction;
use crate::gf::score_dist::ScoreDist;

#[derive(Debug, Default)]
pub struct GeneratingFunctionGroup {
    min_score: i32,
    max_score: i32,
    merged: Option<ScoreDist>,
}

impl GeneratingFunctionGroup {
    pub fn new() -> Self {
        Self {
            min_score: i32::MAX,
            max_score: i32::MIN,
            merged: None,
        }
    }

    /// Merge `gf`'s score distribution into the running aggregate.
    /// Takes `gf` by value so its memory can be released after merging.
    pub fn accept(&mut self, gf: GeneratingFunction) {
        let dist = gf.score_dist();
        let gf_min = dist.min_score();
        let gf_max = dist.max_score();

        if self.merged.is_none() {
            self.min_score = gf_min;
            self.max_score = gf_max;
            let mut m = ScoreDist::new(gf_min, gf_max, false, true);
            m.add_prob_dist(dist, 0, 1.0);
            self.merged = Some(m);
            return;
        }

        let new_min = self.min_score.min(gf_min);
        let new_max = self.max_score.max(gf_max);
        if new_min != self.min_score || new_max != self.max_score {
            let mut expanded = ScoreDist::new(new_min, new_max, false, true);
            expanded.add_prob_dist(self.merged.as_ref().unwrap(), 0, 1.0);
            self.merged = Some(expanded);
            self.min_score = new_min;
            self.max_score = new_max;
        }
        self.merged.as_mut().unwrap().add_prob_dist(dist, 0, 1.0);
    }

    pub fn is_computed(&self) -> bool {
        self.merged.is_some()
    }

    pub fn min_score(&self) -> i32 {
        self.min_score
    }

    pub fn max_score(&self) -> i32 {
        self.max_score
    }

    pub fn score_dist(&self) -> Option<&ScoreDist> {
        self.merged.as_ref()
    }

    pub fn spectral_probability(&self, score: i32) -> Option<f64> {
        self.merged.as_ref().map(|d| d.get_spectral_probability(score))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::aa_set::{AminoAcidSet, AminoAcidSetBuilder};
    use crate::gf::primitive_graph::PrimitiveAaGraph;
    use crate::scoring::{RankScorer, ScoredSpectrum};
    use crate::spectrum::Spectrum;
    use crate::testutil::tiny_param_with_ions;

    fn aa() -> AminoAcidSet {
        AminoAcidSetBuilder::new_standard().build().unwrap()
    }

    fn empty_spec() -> Spectrum {
        Spectrum {
            title: "t".into(),
            precursor_mz: 500.0,
            precursor_intensity: None,
            precursor_charge: Some(2),
            rt_seconds: None,
            scan: None,
            peaks: vec![],
        }
    }

    fn build_gf(peptide_mass: i32) -> GeneratingFunction {
        let aa = aa();
        let s = empty_spec();
        let param = tiny_param_with_ions();
        let scorer = RankScorer::new(&param);
        let ss = ScoredSpectrum::new_without_filtering(&s);
        let g = PrimitiveAaGraph::new(&aa, peptide_mass, None, &ss, &scorer, 2, 1000.0, 0.5, false, false);
        GeneratingFunction::compute(&g, &aa).expect("non-empty GF")
    }

    #[test]
    fn empty_group_is_not_computed() {
        let g = GeneratingFunctionGroup::new();
        assert!(!g.is_computed());
        assert!(g.score_dist().is_none());
        assert!(g.spectral_probability(0).is_none());
    }

    #[test]
    fn single_gf_merge_preserves_distribution() {
        // Accept one GF; merged dist should equal the single GF's dist.
        let gf = build_gf(200);
        let dist_min = gf.min_score();
        let dist_max = gf.max_score();
        let p_at_min = gf.score_dist().get_probability(dist_min);
        let p_at_max_minus_1 = gf.score_dist().get_probability(dist_max - 1);
        let mut group = GeneratingFunctionGroup::new();
        group.accept(gf);
        assert!(group.is_computed());
        assert_eq!(group.min_score(), dist_min);
        assert_eq!(group.max_score(), dist_max);
        let merged = group.score_dist().unwrap();
        assert!((merged.get_probability(dist_min) - p_at_min).abs() < 1e-12);
        assert!((merged.get_probability(dist_max - 1) - p_at_max_minus_1).abs() < 1e-12);
    }

    #[test]
    fn two_gfs_merge_sum_of_probabilities() {
        let gf1 = build_gf(200);
        let gf2 = build_gf(210);
        let dist1_clone = gf1.score_dist().clone();
        let dist2_clone = gf2.score_dist().clone();

        let mut group = GeneratingFunctionGroup::new();
        group.accept(gf1);
        group.accept(gf2);
        assert!(group.is_computed());
        let merged = group.score_dist().unwrap();
        // For each score in either range, merged should equal sum of inputs.
        let test_score = merged.min_score();
        let p_merged = merged.get_probability(test_score);
        let p1 = if test_score >= dist1_clone.min_score() && test_score < dist1_clone.max_score() {
            dist1_clone.get_probability(test_score)
        } else {
            0.0
        };
        let p2 = if test_score >= dist2_clone.min_score() && test_score < dist2_clone.max_score() {
            dist2_clone.get_probability(test_score)
        } else {
            0.0
        };
        assert!(
            (p_merged - (p1 + p2)).abs() < 1e-9,
            "merged at {test_score} = {p_merged}, expected {p1} + {p2}"
        );
    }

    #[test]
    fn expanding_range_keeps_existing_mass() {
        // Accept a small-range GF first, then a wider-range GF. The merged
        // dist's min/max should expand. The sum of all merged probabilities
        // should equal sum of input probs (no probability lost in re-allocation).
        let gf_a = build_gf(200);
        let gf_b = build_gf(300); // typically wider score range due to more nodes
        let total_a: f64 = (gf_a.min_score()..gf_a.max_score())
            .map(|s| gf_a.score_dist().get_probability(s))
            .sum();
        let total_b: f64 = (gf_b.min_score()..gf_b.max_score())
            .map(|s| gf_b.score_dist().get_probability(s))
            .sum();
        let mut group = GeneratingFunctionGroup::new();
        group.accept(gf_a);
        group.accept(gf_b);
        let merged = group.score_dist().unwrap();
        let total_merged: f64 = (merged.min_score()..merged.max_score())
            .map(|s| merged.get_probability(s))
            .sum();
        assert!(
            (total_merged - (total_a + total_b)).abs() < 1e-9,
            "merged total {total_merged} != {total_a} + {total_b}"
        );
    }

    #[test]
    fn spectral_probability_after_merge_clamped_to_one() {
        // After merging multiple GFs, get_spectral_probability is clamped to 1.0
        // (Java behavior). Verify the API returns at most 1.0.
        let mut group = GeneratingFunctionGroup::new();
        for mass in [200, 210, 220, 230, 240] {
            group.accept(build_gf(mass));
        }
        let p_at_min = group.spectral_probability(group.min_score()).unwrap();
        assert!(p_at_min <= 1.0 + 1e-9, "spec prob {p_at_min} > 1.0");
    }
}
