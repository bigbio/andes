use model::fdr::{ranked_q, Counts, Eval};

#[derive(Debug, Clone, Copy)]
pub struct ScoredLabel {
    pub score: f32,
    pub is_decoy: bool,
}

fn ranked(items: &[ScoredLabel]) -> (Vec<usize>, Vec<f64>) {
    let n = items.len();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| {
        let sa = if items[a].score.is_nan() {
            f32::NEG_INFINITY
        } else {
            items[a].score
        };
        let sb = if items[b].score.is_nan() {
            f32::NEG_INFINITY
        } else {
            items[b].score
        };
        sb.partial_cmp(&sa).unwrap_or(std::cmp::Ordering::Equal)
    });
    let tie = |x: f32, y: f32| (x.is_nan() && y.is_nan()) || x.to_bits() == y.to_bits();
    let q = ranked_q::<Counts>(
        n,
        |c, rank| c.add(items[order[rank]].is_decoy),
        |start, rank| tie(items[order[rank]].score, items[order[start]].score),
        Eval::PerRowGroupMax,
        Counts::ratio,
        1.0,
    );
    (order, q)
}

pub fn confident_target_indices(items: &[ScoredLabel], q_threshold: f64) -> Vec<usize> {
    let (order, q) = ranked(items);
    order
        .iter()
        .zip(q.iter())
        .filter(|(&i, &qi)| !items[i].is_decoy && qi <= q_threshold)
        .map(|(&i, _)| i)
        .collect()
}

pub fn qvalues(items: &[ScoredLabel]) -> Vec<f64> {
    let (order, q) = ranked(items);
    let mut out = vec![1.0_f64; items.len()];
    for (rank, &i) in order.iter().enumerate() {
        out[i] = q[rank];
    }
    out
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
        let mut items: Vec<_> = (0..200)
            .map(|i| sl(30.0 - i as f32 * 0.05, false))
            .collect();
        items.push(sl(1.0, true));
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
