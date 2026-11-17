use crate::decoy::PairClass;

pub const MIN_FALSE: f64 = 2.0;

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

    pub fn fdr(&self) -> f64 {
        if self.tt == 0 {
            return 1.0;
        }
        let false_tt = (self.td as f64 + 1.0 - self.dd as f64).max(MIN_FALSE);
        (false_tt / self.tt as f64).min(1.0)
    }

    pub fn fdr_homotypic(&self) -> f64 {
        if self.tt == 0 {
            return 1.0;
        }
        let false_tt = (self.dd as f64 + self.td as f64 + 1.0).max(MIN_FALSE);
        (false_tt / self.tt as f64).min(1.0)
    }
}

pub fn q_values(scores: &[f64], classes: &[PairClass]) -> Vec<f64> {
    q_values_by(scores, classes, ClassCounts::fdr)
}

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
