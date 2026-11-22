use model::fdr::{ranked_q, Eval};

use crate::decoy::PairClass;

/// Floor of the false-target estimate in a small class (see `SMALL_CLASS_TT`).
pub const MIN_FALSE: f64 = 2.0;

/// A class with fewer target-target rows than this is estimated conservatively: a pseudo-count
/// of 1 is added to the false-target estimate, which is floored at `MIN_FALSE`. A larger class
/// uses the plain estimate (TD - DD, or TD + DD for homotypic links).
pub const SMALL_CLASS_TT: usize = 50;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ClassCounts {
    pub tt: usize,
    pub td: usize,
    pub dd: usize,
}

impl ClassCounts {
    pub fn add(&mut self, c: PairClass) {
        match c {
            PairClass::TT => self.tt += 1,
            PairClass::TD => self.td += 1,
            PairClass::DD => self.dd += 1,
        }
    }

    /// (TD - DD) / TT, at least 0 and at most 1; 1 without TT.
    pub fn fdr(&self) -> f64 {
        self.ratio((self.td as f64 - self.dd as f64).max(0.0))
    }

    /// (TD + 1 - DD) / TT with the false count floored at `MIN_FALSE`: the estimate of a small
    /// class.
    pub fn fdr_small(&self) -> f64 {
        self.ratio((self.td as f64 + 1.0 - self.dd as f64).max(MIN_FALSE))
    }

    /// (TD + DD) / TT: homotypic links have no distinct TD/DT decoy pair, so their DD rows
    /// stand for false targets and are added.
    pub fn fdr_homotypic(&self) -> f64 {
        self.ratio((self.dd + self.td) as f64)
    }

    /// The homotypic estimate of a small class: pseudo-count 1, floored at `MIN_FALSE`.
    pub fn fdr_homotypic_small(&self) -> f64 {
        self.ratio((self.dd as f64 + self.td as f64 + 1.0).max(MIN_FALSE))
    }

    fn ratio(&self, false_tt: f64) -> f64 {
        if self.tt == 0 {
            return 1.0;
        }
        (false_tt / self.tt as f64).min(1.0)
    }
}

/// q-values of one link class: the plain estimate, or the small-class one when the class has
/// fewer than `SMALL_CLASS_TT` target-target rows.
pub fn q_values(scores: &[f64], classes: &[PairClass]) -> Vec<f64> {
    let est = if is_small(classes) {
        ClassCounts::fdr_small
    } else {
        ClassCounts::fdr
    };
    q_values_by(scores, classes, est)
}

/// q-values of homotypic links, with the small-class estimate as in `q_values`.
pub fn q_values_homotypic(scores: &[f64], classes: &[PairClass]) -> Vec<f64> {
    let est = if is_small(classes) {
        ClassCounts::fdr_homotypic_small
    } else {
        ClassCounts::fdr_homotypic
    };
    q_values_by(scores, classes, est)
}

fn is_small(classes: &[PairClass]) -> bool {
    classes.iter().filter(|&&c| c == PairClass::TT).count() < SMALL_CLASS_TT
}

fn q_values_by(scores: &[f64], classes: &[PairClass], fdr: fn(&ClassCounts) -> f64) -> Vec<f64> {
    assert_eq!(scores.len(), classes.len());
    let mut order: Vec<usize> = (0..scores.len()).collect();
    order.sort_by(|&a, &b| scores[b].total_cmp(&scores[a]));
    let ranked = ranked_q::<ClassCounts>(
        order.len(),
        |c, rank| c.add(classes[order[rank]]),
        |start, rank| scores[order[rank]].total_cmp(&scores[order[start]]).is_eq(),
        Eval::GroupEnd,
        fdr,
        f64::INFINITY,
    );
    let mut q = vec![0.0; scores.len()];
    for (&i, qi) in order.iter().zip(ranked) {
        q[i] = qi;
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
    fn fdr_subtracts_dd_and_a_small_class_is_floored() {
        let c = ClassCounts {
            tt: 100,
            td: 7,
            dd: 2,
        };
        assert!((c.fdr() - 0.05).abs() < 1e-12);
        assert!((c.fdr_small() - 0.06).abs() < 1e-12);
        let c = ClassCounts {
            tt: 10,
            td: 1,
            dd: 3,
        };
        assert_eq!(c.fdr(), 0.0);
        assert!((c.fdr_small() - 0.2).abs() < 1e-12);
        assert_eq!(ClassCounts::default().fdr(), 1.0);
    }

    #[test]
    fn homotypic_fdr_counts_dd_as_false_targets() {
        let c = ClassCounts {
            tt: 100,
            td: 0,
            dd: 9,
        };
        assert!((c.fdr_homotypic() - 0.09).abs() < 1e-12);
        assert!((c.fdr_homotypic_small() - 0.10).abs() < 1e-12);
        let q = q_values_homotypic(&[3.0, 2.0, 1.0], &[TT, TT, DD]);
        assert!((q[0] - 1.0).abs() < 1e-12 && (q[2] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn a_small_class_cannot_reach_a_low_q_and_a_large_one_can() {
        let small: Vec<f64> = (0..25).map(|i| 100.0 - i as f64).collect();
        let q = q_values(&small, &[TT; 25]);
        assert!(q.iter().all(|&v| (v - 0.08).abs() < 1e-12));
        let large: Vec<f64> = (0..60).map(|i| 100.0 - i as f64).collect();
        let q = q_values(&large, &[TT; 60]);
        assert!(q.iter().all(|&v| v == 0.0));
    }

    #[test]
    fn q_values_are_monotone_and_share_ties() {
        let mut scores: Vec<f64> = (0..205).map(|i| 1000.0 - i as f64).collect();
        let mut classes = vec![TT; 200];
        classes.extend([TD, TD, TD, TD, DD]);
        scores.extend([1.0, 1.0]);
        classes.extend([TT, TD]);
        let q = q_values(&scores, &classes);
        assert_eq!(q[0], 0.0);
        assert!((q[204] - 3.0 / 200.0).abs() < 1e-12);
        assert_eq!(q[205], q[206]);
        assert!((q[205] - 4.0 / 201.0).abs() < 1e-12);
        assert!(q.windows(2).take(204).all(|w| w[0] <= w[1]));
    }
}
