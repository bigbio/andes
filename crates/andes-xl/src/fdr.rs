//! Crosslink FDR from target-target, target-decoy and decoy-decoy counts.
//!
//! Above a score threshold, TD matches estimate the TT matches with one false
//! peptide plus those with two, and DD matches estimate those with two false
//! peptides counted twice over in TD, so the false TT matches number
//! `TD − DD` and `FDR = (TD − DD) / TT`, clamped at 0 when `DD > TD`. The estimate is made
//! separately for each link class (intra- and inter-protein) by the caller.

use crate::decoy::PairClass;

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

    /// `(TD − DD) / TT`, clamped at 0; 0 when there are no TT matches.
    pub fn fdr(&self) -> f64 {
        if self.tt == 0 {
            return 0.0;
        }
        (self.td as f64 - self.dd as f64).max(0.0) / self.tt as f64
    }
}

/// q-values of CSMs with scores `scores` (higher is better) and classes
/// `classes`, by the `(TD − DD) / TT` estimate at each score threshold, made
/// monotone from the bottom up. Ties share one threshold.
pub fn q_values(scores: &[f64], classes: &[PairClass]) -> Vec<f64> {
    assert_eq!(scores.len(), classes.len());
    let mut order: Vec<usize> = (0..scores.len()).collect();
    order.sort_by(|&a, &b| scores[b].total_cmp(&scores[a]));
    let mut fdr_at = vec![0.0; scores.len()];
    let mut counts = ClassCounts::default();
    let mut k = 0;
    while k < order.len() {
        let mut end = k;
        while end < order.len() && scores[order[end]] == scores[order[k]] {
            counts.add(classes[order[end]]);
            end += 1;
        }
        let f = counts.fdr();
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
    fn fdr_subtracts_dd_and_clamps() {
        let c = ClassCounts {
            tt: 100,
            td: 7,
            dd: 2,
        };
        assert!((c.fdr() - 0.05).abs() < 1e-12);
        let c = ClassCounts {
            tt: 10,
            td: 1,
            dd: 3,
        };
        assert_eq!(c.fdr(), 0.0);
        assert_eq!(ClassCounts::default().fdr(), 0.0);
    }

    #[test]
    fn q_values_are_monotone_and_share_ties() {
        let scores = [10.0, 9.0, 8.0, 8.0, 7.0, 6.0];
        let classes = [TT, TT, TD, TT, TT, TD];
        let q = q_values(&scores, &classes);
        // Thresholds: 10 -> 0/1, 9 -> 0/2, 8 -> 1/3, 7 -> 1/4, 6 -> 2/4.
        assert_eq!(q[0], 0.0);
        assert_eq!(q[1], 0.0);
        assert!((q[2] - 0.25).abs() < 1e-12);
        assert_eq!(q[2], q[3]);
        assert!((q[4] - 0.25).abs() < 1e-12);
        assert!((q[5] - 0.5).abs() < 1e-12);
    }
}
