//! Match-between-runs library confidence.
//!
//! The functions here decide which precursors enter the transfer library. The
//! values are internal quant-quality gates, not reported identification FDR;
//! identification FDR is only ever reported by the PSM-level rescorer.

/// Precursor-level library entry gate (internal; not an identification FDR).
///
/// `rows` holds `(precursor_id, is_decoy, psm_score)` with higher scores
/// better. Each precursor keeps its best PSM score across all runs, then
/// precursors are ranked by `(score desc, decoys before targets, id asc)` (ties
/// are conservative: a target tied with a decoy never gets a lower q than it
/// would with the decoy first) and run through
/// target/decoy competition: `(D + 1) / max(T, 1)` at each rank, a running
/// minimum from the worst rank upwards, clamped to at most 1.
///
/// A NaN score is treated as the worst possible score. A precursor id is
/// expected to be all-target or all-decoy (checked by `debug_assert!`); in release
/// builds a mixed id is decided by the first row seen. Ids are assumed dense:
/// the result is indexed by precursor id and has length `max id + 1`, `None`
/// for ids without rows.
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
            Some((d, s)) => {
                debug_assert_eq!(
                    *d, is_decoy,
                    "precursor id {id} mixes target and decoy rows"
                );
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
    // Score desc; on ties decoys first (conservative), then id for determinism.
    ranked.sort_by(|a, b| b.2.total_cmp(&a.2).then(b.1.cmp(&a.1)).then(a.0.cmp(&b.0)));

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
        // 100 strong targets, a decoy at 1.0, and a NaN-scored target. If NaN
        // were ranked best the NaN target would get the low q of the strong
        // targets; ranked worst it sits below the decoy.
        let mut rows: Vec<(u32, bool, f64)> =
            (0..100).map(|i| (i, false, 100.0 + i as f64)).collect();
        rows.push((100, true, 1.0));
        rows.push((101, false, f64::NAN));
        let q = precursor_qvalues(&rows);
        let (strong, decoy, nan) = (q[0].unwrap(), q[100].unwrap(), q[101].unwrap());
        assert!(nan >= decoy, "{nan} < {decoy}");
        assert!(nan > strong, "{nan} <= {strong}");
        // A NaN in one PSM does not hide a finite score of the same precursor.
        let q = precursor_qvalues(&[(0, false, f64::NAN), (0, false, 4.0)]);
        assert_eq!(q[0], Some(1.0));
    }

    #[test]
    fn target_decoy_tie_is_conservative() {
        // 9 strong targets, then a target and a decoy tied at 5.0. Decoy-first
        // gives the tied target (D+1)/T = 2/10; target-first would give 1/10.
        let strong: Vec<(u32, bool, f64)> = (0..9).map(|i| (i, false, 100.0 + i as f64)).collect();
        let mut a = strong.clone();
        a.extend([(9, false, 5.0), (10, true, 5.0)]);
        let mut b = strong;
        b.extend([(9, true, 5.0), (10, false, 5.0)]);
        let ta = precursor_qvalues(&a)[9].unwrap();
        let tb = precursor_qvalues(&b)[10].unwrap();
        assert!((ta - 0.2).abs() < 1e-12, "{ta}");
        assert_eq!(ta, tb);
    }
}
