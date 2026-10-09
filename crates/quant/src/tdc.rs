//! Picked target–decoy competition for feature-level q-values.
//!
//! Every quantified target has a decoy twin extracted with the same
//! procedure at a shifted m/z. Unrelated real peptides can occur there; the
//! twin is intended to model a spurious match to the requested ion. The
//! competition keeps the better of the pair, and the usual
//! `(decoys + 1) / targets` running ratio over the winners gives each surviving
//! target an estimated q-value. This is not IonQuant's mixture-model MBR
//! estimator, nor a substitute for validating the decoy model on real data.
//! Apply any feature eligibility filters identically to both sides BEFORE
//! competition; filtering only the reported targets changes the tested set.

/// Score pair of one target and its decoy twin; `None` when no peak was found.
#[derive(Debug, Clone, Copy)]
pub struct Pair {
    pub target: Option<f32>,
    pub decoy: Option<f32>,
}

/// Per-pair q-value of the target, or `None` when the target had no peak or
/// lost to its decoy. Nonfinite scores are treated as missing. Pairwise ties
/// conservatively go to the decoy; equal winning scores share one threshold.
pub fn picked_qvalues(pairs: &[Pair]) -> Vec<Option<f64>> {
    // Winner per pair: (score, is_decoy, pair index).
    let mut winners: Vec<(f32, bool, usize)> = Vec::with_capacity(pairs.len());
    for (i, p) in pairs.iter().enumerate() {
        match (
            p.target.filter(|s| s.is_finite()),
            p.decoy.filter(|s| s.is_finite()),
        ) {
            (Some(t), Some(d)) => {
                if d >= t {
                    winners.push((d, true, i));
                } else {
                    winners.push((t, false, i));
                }
            }
            (Some(t), None) => winners.push((t, false, i)),
            (None, Some(d)) => winners.push((d, true, i)),
            (None, None) => {}
        }
    }
    winners.sort_by(|a, b| b.0.total_cmp(&a.0));
    let mut q = vec![1.0; winners.len()];
    let (mut n_t, mut n_d) = (0usize, 0usize);
    let mut start = 0;
    while start < winners.len() {
        let mut end = start + 1;
        while end < winners.len() && winners[end].0 == winners[start].0 {
            end += 1;
        }
        for &(_, is_decoy, _) in &winners[start..end] {
            if is_decoy {
                n_d += 1;
            } else {
                n_t += 1;
            }
        }
        q[start..end].fill((n_d as f64 + 1.0) / n_t.max(1) as f64);
        start = end;
    }
    // Monotone q: minimum of the raw FDR at this score and every lower score.
    let mut running = f64::INFINITY;
    for value in q.iter_mut().rev() {
        running = running.min(*value);
        *value = running.min(1.0);
    }
    let mut out = vec![None; pairs.len()];
    for (k, &(_, is_decoy, i)) in winners.iter().enumerate() {
        if !is_decoy {
            out[i] = Some(q[k]);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(t: Option<f32>, d: Option<f32>) -> Pair {
        Pair {
            target: t,
            decoy: d,
        }
    }

    #[test]
    fn clean_targets_get_small_q_and_losers_none() {
        let pairs = vec![
            p(Some(0.9), Some(0.1)),
            p(Some(0.8), None),
            p(Some(0.2), Some(0.5)), // loses to decoy
            p(None, None),
        ];
        let q = picked_qvalues(&pairs);
        assert!(q[0].unwrap() <= 0.5 && q[1].unwrap() <= 0.5);
        assert_eq!(q[2], None);
        assert_eq!(q[3], None);
        // q is monotone non-decreasing with decreasing score.
        assert!(q[0].unwrap() <= q[1].unwrap());
    }

    #[test]
    fn decoy_heavy_tail_raises_q() {
        let mut pairs: Vec<Pair> = (0..50)
            .map(|i| p(Some(1.0 - i as f32 * 0.01), None))
            .collect();
        pairs.extend((0..50).map(|i| p(Some(0.2 - i as f32 * 0.001), Some(0.3))));
        let q = picked_qvalues(&pairs);
        let strong = q[0].unwrap();
        assert!(strong < 0.05, "{strong}");
        // Everything in the tail lost to a decoy → None.
        assert!(q[50..].iter().all(|v| v.is_none()));
    }

    #[test]
    fn equal_target_and_decoy_scores_do_not_create_confident_targets() {
        let pairs = vec![p(Some(0.5), Some(0.5)); 100];
        assert!(picked_qvalues(&pairs).iter().all(Option::is_none));
    }

    #[test]
    fn equal_score_threshold_includes_all_decoy_winners() {
        let mut pairs = vec![p(Some(0.5), None); 100];
        pairs.extend(vec![p(None, Some(0.5)); 100]);
        let q = picked_qvalues(&pairs);
        assert!(q[..100].iter().all(|&q| q == Some(1.0)), "{q:?}");
        pairs.reverse();
        let mut reversed = picked_qvalues(&pairs);
        reversed.reverse();
        assert_eq!(q, reversed);
    }

    #[test]
    fn nonfinite_scores_are_missing_observations() {
        let q = picked_qvalues(&[
            p(Some(f32::NAN), Some(0.5)),
            p(Some(f32::INFINITY), None),
            p(Some(0.8), Some(f32::NAN)),
        ]);
        assert_eq!(q, vec![None, None, Some(1.0)]);
    }
}
