//! Crosslink FDR from target-target, target-decoy and decoy-decoy counts.
//!
//! Above a score threshold, TD matches estimate the TT matches with one false
//! peptide plus those with two, and DD matches estimate those with two false
//! peptides counted twice over in TD, so the false TT matches number
//! `TD − DD` and `FDR = (TD − DD) / TT`. The estimate is made separately for
//! each link class (intra- and inter-protein) by the caller.
//!
//! At small counts `TD − DD` fluctuates to 0 or below by chance, and a
//! monotone pass would carry that dip down to every lower threshold. The
//! estimate is therefore floored at [`MIN_FALSE`]` / TT`: a set is accepted at
//! FDR `α` only when `TT · α ≥ MIN_FALSE`, the minimum-TD validity check of
//! xiFDR (Fischer & Rappsilber 2017).

use crate::decoy::PairClass;

/// Fewest false matches an accepted set is assumed to hold.
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

    /// `max(TD − DD, MIN_FALSE) / TT`, capped at 1; 1 when there are no TT
    /// matches.
    pub fn fdr(&self) -> f64 {
        if self.tt == 0 {
            return 1.0;
        }
        let false_tt = (self.td as f64 - self.dd as f64).max(MIN_FALSE);
        (false_tt / self.tt as f64).min(1.0)
    }
}

/// q-values of CSMs with scores `scores` (higher is better) and classes
/// `classes`, by the floored `(TD − DD) / TT` estimate at each score threshold, made
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
        while end < order.len() && scores[order[end]].total_cmp(&scores[order[k]]).is_eq() {
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
        assert!((c.fdr() - 0.05).abs() < 1e-12);
        // DD above TD is a chance dip, not zero false matches.
        let c = ClassCounts {
            tt: 10,
            td: 1,
            dd: 3,
        };
        assert!((c.fdr() - 0.2).abs() < 1e-12);
        assert_eq!(ClassCounts::default().fdr(), 1.0);
    }

    #[test]
    fn small_sets_cannot_reach_a_low_q() {
        // 25 TT and no decoys: q is 2/25, never 0.
        let scores: Vec<f64> = (0..25).map(|i| 100.0 - i as f64).collect();
        let q = q_values(&scores, &[TT; 25]);
        assert!(q.iter().all(|&v| (v - 0.08).abs() < 1e-12));
    }

    #[test]
    fn q_values_are_monotone_and_share_ties() {
        // 200 TT above 4 TD and 1 DD, then a tie of one TT and one TD.
        let mut scores: Vec<f64> = (0..205).map(|i| 1000.0 - i as f64).collect();
        let mut classes = vec![TT; 200];
        classes.extend([TD, TD, TD, TD, DD]);
        scores.extend([1.0, 1.0]);
        classes.extend([TT, TD]);
        let q = q_values(&scores, &classes);
        // Above the decoys: 2/200 floor. All 205: (4-1)/200. Tie: (5-1)/201.
        assert!((q[0] - 0.01).abs() < 1e-12);
        assert!((q[204] - 0.015).abs() < 1e-12);
        assert_eq!(q[205], q[206]);
        assert!((q[205] - 4.0 / 201.0).abs() < 1e-12);
        assert!(q.windows(2).take(204).all(|w| w[0] <= w[1]));
    }
}
