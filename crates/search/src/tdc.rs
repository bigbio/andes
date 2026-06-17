//! Target-decoy competition q-values over scored PSMs. Single source of truth
//! for the calibration pre-pass and the refinement confident-protein gate.

/// Per-item input to the TDC q-walk: a discriminating score (higher = better)
/// and whether the item is a decoy.
#[derive(Debug, Clone, Copy)]
pub struct ScoredLabel {
    pub score: f32,
    pub is_decoy: bool,
}

/// Return the indices (into `items`, original order) of TARGET items whose
/// monotone TDC q-value is <= `q_threshold`. The walk: sort by score desc,
/// q = decoys/max(targets,1), monotone-from-bottom, conservative ties (every
/// item in an equal-score bucket takes the worst q in the bucket). Empty input
/// -> empty output.
pub fn confident_target_indices(items: &[ScoredLabel], q_threshold: f64) -> Vec<usize> {
    if items.is_empty() {
        return Vec::new();
    }
    let n = items.len();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| {
        let sa = if items[a].score.is_nan() { f32::NEG_INFINITY } else { items[a].score };
        let sb = if items[b].score.is_nan() { f32::NEG_INFINITY } else { items[b].score };
        sb.partial_cmp(&sa).unwrap_or(std::cmp::Ordering::Equal)
    });
    let mut q = vec![1.0_f64; n];
    let (mut targets, mut decoys) = (0u64, 0u64);
    for (rank, &i) in order.iter().enumerate() {
        if items[i].is_decoy {
            decoys += 1;
        } else {
            targets += 1;
        }
        q[rank] = decoys as f64 / targets.max(1) as f64;
    }
    let mut min_q = 1.0_f64;
    for qi in q.iter_mut().rev() {
        if *qi < min_q {
            min_q = *qi;
        }
        *qi = min_q;
    }
    // Conservative ties: equal-score contiguous buckets take the worst q.
    let mut start = 0usize;
    while start < n {
        let s = items[order[start]].score;
        let mut end = start + 1;
        let tie = |x: f32, y: f32| (x.is_nan() && y.is_nan()) || x.to_bits() == y.to_bits();
        while end < n && tie(items[order[end]].score, s) {
            end += 1;
        }
        let worst = q[start..end].iter().cloned().fold(0.0_f64, f64::max);
        for qi in &mut q[start..end] {
            *qi = worst;
        }
        start = end;
    }
    order
        .iter()
        .zip(q.iter())
        .filter(|(&i, &qi)| !items[i].is_decoy && qi <= q_threshold)
        .map(|(&i, _)| i)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    fn sl(score: f32, is_decoy: bool) -> ScoredLabel {
        ScoredLabel { score, is_decoy }
    }

    #[test]
    fn decoy_only_yields_nothing() {
        let items: Vec<_> = (0..50).map(|i| sl(20.0 - i as f32 * 0.1, true)).collect();
        assert!(confident_target_indices(&items, 0.01).is_empty());
    }

    #[test]
    fn confident_targets_pass_and_indices_are_targets() {
        let mut items: Vec<_> = (0..200).map(|i| sl(30.0 - i as f32 * 0.05, false)).collect();
        items.push(sl(1.0, true)); // lone tail decoy
        let conf = confident_target_indices(&items, 0.01);
        assert!(!conf.is_empty());
        assert!(conf.iter().all(|&i| !items[i].is_decoy));
    }

    #[test]
    fn interleaved_low_confidence_yields_nothing() {
        let items: Vec<_> = (0..100).map(|i| sl(5.0, i % 2 == 0)).collect();
        assert!(confident_target_indices(&items, 0.01).is_empty());
    }
}
