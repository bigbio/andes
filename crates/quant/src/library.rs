//! Match-between-runs library confidence.
//!
//! The functions here decide which precursors enter the transfer library. The
//! values are internal quant-quality gates, not reported identification FDR;
//! identification FDR is only ever reported by the PSM-level rescorer.

/// Precursor-level library entry gate (internal; not an identification FDR).
///
/// `rows` holds `(precursor_id, is_decoy, psm_score)` with higher scores
/// better. Each precursor keeps its best PSM score across all runs, then
/// precursors are ranked by `(score desc, id asc)` and run through
/// target/decoy competition: `(D + 1) / max(T, 1)` at each rank, a running
/// minimum from the worst rank upwards, clamped to at most 1.
///
/// A NaN score is treated as the worst possible score. A precursor id is
/// expected to be all-target or all-decoy; if it is mixed, the first row seen
/// decides. The result is indexed by precursor id, `None` for ids without rows.
pub fn precursor_qvalues(rows: &[(u32, bool, f64)]) -> Vec<Option<f64>> {
    let Some(max_id) = rows.iter().map(|r| r.0).max() else {
        return Vec::new();
    };
    let n = max_id as usize + 1;
    // (is_decoy, best score) per id; NaN is mapped to -inf so it is the worst.
    let mut best: Vec<Option<(bool, f64)>> = vec![None; n];
    for &(id, is_decoy, score) in rows {
        let score = if score.is_nan() {
            f64::NEG_INFINITY
        } else {
            score
        };
        let slot = &mut best[id as usize];
        match slot {
            Some((_, s)) => {
                if score > *s {
                    *s = score;
                }
            }
            None => *slot = Some((is_decoy, score)),
        }
    }
    let mut ranked: Vec<(u32, bool, f64)> = best
        .iter()
        .enumerate()
        .filter_map(|(id, b)| b.map(|(d, s)| (id as u32, d, s)))
        .collect();
    ranked.sort_by(|a, b| b.2.total_cmp(&a.2).then(a.0.cmp(&b.0)));

    let mut fdr = Vec::with_capacity(ranked.len());
    let (mut t, mut d) = (0u64, 0u64);
    for &(_, is_decoy, _) in &ranked {
        if is_decoy {
            d += 1;
        } else {
            t += 1;
        }
        fdr.push((d + 1) as f64 / t.max(1) as f64);
    }
    let mut running = f64::INFINITY;
    let mut out = vec![None; n];
    for (i, &(id, _, _)) in ranked.iter().enumerate().rev() {
        running = running.min(fdr[i]);
        out[id as usize] = Some(running.min(1.0));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn best_psm_per_precursor_then_tdc() {
        // p0 target scores 10 and 3 -> best 10; p1 decoy 9; p2 target 8; p3 target 1
        let rows = [
            (0, false, 10.0),
            (0, false, 3.0),
            (1, true, 9.0),
            (2, false, 8.0),
            (3, false, 1.0),
        ];
        let q = precursor_qvalues(&rows);
        // raw (D+1)/T down the list: p0 1/1, p1 2/1 -> 1 (clamp), p2 2/2, p3 2/3;
        // the running minimum from the bottom makes every one 2/3
        assert_eq!(q.len(), 4);
        for (i, qi) in q.iter().enumerate() {
            assert!((qi.unwrap() - 2.0 / 3.0).abs() < 1e-12, "{i}: {qi:?}");
        }
    }

    #[test]
    fn many_targets_reach_low_q() {
        let mut rows: Vec<(u32, bool, f64)> =
            (0..1000).map(|i| (i, false, 1000.0 - i as f64)).collect();
        rows.push((1000, true, 0.5));
        let q = precursor_qvalues(&rows);
        assert!(q[0].unwrap() <= 0.01 && q[999].unwrap() <= 0.01);
    }

    #[test]
    fn missing_ids_are_none_and_empty_input_is_empty() {
        assert!(precursor_qvalues(&[]).is_empty());
        let q = precursor_qvalues(&[(2, false, 1.0)]);
        assert_eq!(q.len(), 3);
        assert!(q[0].is_none() && q[1].is_none() && q[2].is_some());
    }

    #[test]
    fn nan_scores_rank_worst_and_do_not_panic() {
        // The NaN target ranks below the finite decoy, so it cannot beat it.
        let rows = [(0, false, f64::NAN), (1, true, 1.0), (2, false, 5.0)];
        let q = precursor_qvalues(&rows);
        assert!(q[2].unwrap() < q[0].unwrap() + 1e-12);
        assert_eq!(q[2], Some(1.0));
        // A NaN in one PSM does not hide a finite score of the same precursor.
        let q = precursor_qvalues(&[(0, false, f64::NAN), (0, false, 4.0)]);
        assert_eq!(q[0], Some(1.0 / 1.0));
    }
}
