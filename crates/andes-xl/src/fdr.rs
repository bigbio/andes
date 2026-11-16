//! Crosslink FDR from target-target, target-decoy and decoy-decoy counts.
//!
//! Above a score threshold, TD matches estimate the TT matches with one false
//! peptide plus those with two, and DD matches estimate those with two false
//! peptides counted twice over in TD, so the false TT matches number
//! `TD − DD` and `FDR = (TD − DD) / TT`. The estimate counts one TD more,
//! `(TD + 1 − DD) / TT`, which removes most of the excess the monotone pass
//! picks up from chance dips in mid-sized sets. It is made separately for
//! each link class (intra- and inter-protein) by the caller.

use crate::decoy::PairClass;

pub const MIN_FALSE: f64 = 2.0;

/// Counts of a set of CSMs by class.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ClassCounts {
    /// Target-target matches.
    pub tt: usize,
    /// Target-decoy matches.
    pub td: usize,
    /// Decoy-decoy matches.
    pub dd: usize,
}

impl ClassCounts {
    /// Add one CSM of class `c`.
    pub fn add(&mut self, c: PairClass) {
        match c {
            PairClass::TT => self.tt += 1,
            PairClass::TD => self.td += 1,
            PairClass::DD => self.dd += 1,
        }
    }

    pub fn fdr(&self) -> f64 {
        if self.tt == 0 {
            return 1.0;
        }
        let false_tt = (self.td as f64 + 1.0 - self.dd as f64).max(MIN_FALSE);
        (false_tt / self.tt as f64).min(1.0)
    }

    /// FDR of homotypic pairs, the same peptide sequence on both sides. A random match
    /// pairs a random peptide with itself, so it is TT or DD with equal odds and
    /// practically never TD: the false TT number about DD, not TD − DD.
    pub fn fdr_homotypic(&self) -> f64 {
        if self.tt == 0 {
            return 1.0;
        }
        let false_tt = (self.dd as f64 + self.td as f64 + 1.0).max(MIN_FALSE);
        (false_tt / self.tt as f64).min(1.0)
    }
}

/// q-values of CSMs with scores `scores` (higher is better) and classes
/// `classes`, by the floored `(TD − DD) / TT` estimate at each score threshold, made
/// monotone from the bottom up. Ties share one threshold.
pub fn q_values(scores: &[f64], classes: &[PairClass]) -> Vec<f64> {
    q_values_by(scores, classes, ClassCounts::fdr)
}

/// [`q_values`] for homotypic pairs ([`ClassCounts::fdr_homotypic`]).
pub fn q_values_homotypic(scores: &[f64], classes: &[PairClass]) -> Vec<f64> {
    q_values_by(scores, classes, ClassCounts::fdr_homotypic)
}

fn q_values_by(scores: &[f64], classes: &[PairClass], fdr: fn(&ClassCounts) -> f64) -> Vec<f64> {
    assert_eq!(scores.len(), classes.len());
    let mut order: Vec<usize> = (0..scores.len()).collect();
    order.sort_by(|&a, &b| scores[b].total_cmp(&scores[a]));
    let mut fdr_at = vec![0.0; scores.len()];
    let mut counts = ClassCounts::default();
    let mut k = 0;
    while k < order.len() {
        let mut end = k;
        while end < order.len() && scores[order[end]].total_cmp(&scores[order[k]]).is_eq() {
            counts.add(classes[order[end]]);
            end += 1;
        }
        let f = fdr(&counts);
        for &i in &order[k..end] {
            fdr_at[i] = f;
        }
        k = end;
    }
    // Monotone: q at a threshold is the lowest FDR at that threshold or below.
    let mut q = fdr_at.clone();
    let mut running = f64::INFINITY;
    for &i in order.iter().rev() {
        running = running.min(fdr_at[i]);
        q[i] = running;
    }
    q
}

#[cfg(test)]
mod tests {
    use super::*;
    use PairClass::*;

    #[test]
    fn q_values_terminate_on_nan_scores() {
        let q = q_values(&[2.0, f64::NAN, f64::NAN, 1.0], &[TT, TD, TT, TT]);
        assert_eq!(q.len(), 4);
    }

    #[test]
    fn fdr_subtracts_dd_and_is_floored() {
        let c = ClassCounts {
            tt: 100,
            td: 7,
            dd: 2,
        };
        assert!((c.fdr() - 0.06).abs() < 1e-12);
        let c = ClassCounts {
            tt: 10,
            td: 1,
            dd: 3,
        };
        assert!((c.fdr() - 0.2).abs() < 1e-12);
        assert_eq!(ClassCounts::default().fdr(), 1.0);
    }

    #[test]
    fn homotypic_fdr_counts_dd_as_false_targets() {
        let c = ClassCounts {
            tt: 100,
            td: 0,
            dd: 9,
        };
        assert!((c.fdr_homotypic() - 0.10).abs() < 1e-12);
        assert!((c.fdr() - 0.02).abs() < 1e-12);
        let q = q_values_homotypic(&[3.0, 2.0, 1.0], &[TT, TT, DD]);
        assert!((q[0] - 1.0).abs() < 1e-12 && (q[2] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn small_sets_cannot_reach_a_low_q() {
        let scores: Vec<f64> = (0..25).map(|i| 100.0 - i as f64).collect();
        let q = q_values(&scores, &[TT; 25]);
        assert!(q.iter().all(|&v| (v - 0.08).abs() < 1e-12));
    }

    #[test]
    fn q_values_are_monotone_and_share_ties() {
        let mut scores: Vec<f64> = (0..205).map(|i| 1000.0 - i as f64).collect();
        let mut classes = vec![TT; 200];
        classes.extend([TD, TD, TD, TD, DD]);
        scores.extend([1.0, 1.0]);
        classes.extend([TT, TD]);
        let q = q_values(&scores, &classes);
        assert!((q[0] - 0.01).abs() < 1e-12);
        assert!((q[204] - 0.02).abs() < 1e-12);
        assert_eq!(q[205], q[206]);
        assert!((q[205] - 5.0 / 201.0).abs() < 1e-12);
        assert!(q.windows(2).take(204).all(|w| w[0] <= w[1]));
    }
}
