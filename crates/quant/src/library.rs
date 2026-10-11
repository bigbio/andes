//! Match-between-runs library confidence.
//!
//! The functions here decide which precursors enter the transfer library. The
//! values are internal quant-quality gates, not reported identification FDR;
//! identification FDR is only ever reported by the PSM-level rescorer.

use crate::align::{align, AlignParams, Alignment, Observation};

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

/// One PSM of one run, as the library sees it. `precursor` is a dense id in
/// `0..meta.len()` shared by every run.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PsmIn {
    pub run: u16,
    pub precursor: u32,
    pub rt: f64,
    pub score: f64,
    pub is_decoy: bool,
    pub faims_cv: Option<f32>,
}

/// Identity of a precursor: peptidoform and charge key, monoisotopic m/z and
/// the theoretical isotope fractions.
#[derive(Debug, Clone, PartialEq)]
pub struct PrecursorMeta {
    pub key: String,
    pub mz: f64,
    pub charge: u8,
    pub envelope: [f32; 4],
}

/// A confident precursor that is measured (identified or transferred) in
/// every run.
#[derive(Debug, Clone, PartialEq)]
pub struct LibraryEntry {
    pub precursor: u32,
    /// Library confidence (internal quant gate, not an identification FDR).
    pub q: f64,
    /// Consensus RT from the joint alignment (NaN without an observation in
    /// an aligned run).
    pub consensus_rt: f64,
    /// Most frequent FAIMS compensation voltage among its PSMs; the smallest
    /// value wins a tie; `None` when no PSM has one.
    pub faims_cv: Option<f32>,
    /// Runs with any PSM of the precursor, ascending.
    pub donors: Vec<u16>,
    /// Per donor run (ascending), the RT of its best-scoring PSM (earliest RT
    /// on equal scores).
    pub best_rt: Vec<(u16, f64)>,
}

/// The transfer library: entries sorted by m/z, the joint alignment of the
/// entries' observations, and the precursor table indexed by id.
#[derive(Debug, Clone)]
pub struct Library {
    pub entries: Vec<LibraryEntry>,
    pub alignment: Alignment,
    pub meta: Vec<PrecursorMeta>,
}

/// Builds the library from the PSMs of all runs.
///
/// Precursors are independent: confidence is [`precursor_qvalues`] over the
/// best score per precursor id, with no per-scan competition. Entries are
/// target precursors with `q <= max_q`, sorted by `(mz, id)`. The alignment
/// is fitted on one observation per `(entry, donor run)`, the RT of the best
/// PSM there; decoys never enter it. PSMs with an id outside `0..meta.len()`
/// or a run outside `0..n_runs` are ignored.
pub fn build_library(
    psms: &[PsmIn],
    meta: Vec<PrecursorMeta>,
    n_runs: usize,
    max_q: f64,
    ap: &AlignParams,
) -> Library {
    let n_prec = meta.len();
    let valid: Vec<&PsmIn> = psms
        .iter()
        .filter(|p| (p.precursor as usize) < n_prec && (p.run as usize) < n_runs)
        .collect();
    let rows: Vec<(u32, bool, f64)> = valid
        .iter()
        .map(|p| (p.precursor, p.is_decoy, p.score))
        .collect();
    let qs = precursor_qvalues(&rows);

    // Group the target PSMs of confident precursors by precursor id.
    let mut by_prec: Vec<&PsmIn> = valid
        .iter()
        .copied()
        .filter(|p| {
            !p.is_decoy
                && qs
                    .get(p.precursor as usize)
                    .copied()
                    .flatten()
                    .is_some_and(|q| q <= max_q)
        })
        .collect();
    by_prec.sort_by_key(|p| (p.precursor, p.run));

    let mut entries: Vec<LibraryEntry> = Vec::new();
    let mut obs: Vec<Observation> = Vec::new();
    let mut i = 0;
    while i < by_prec.len() {
        let id = by_prec[i].precursor;
        let mut j = i;
        while j < by_prec.len() && by_prec[j].precursor == id {
            j += 1;
        }
        let group = &by_prec[i..j];
        i = j;

        let mut donors: Vec<u16> = Vec::new();
        let mut best_rt: Vec<(u16, f64)> = Vec::new();
        let mut k = 0;
        while k < group.len() {
            let run = group[k].run;
            let mut best: Option<&PsmIn> = None;
            while k < group.len() && group[k].run == run {
                let c = group[k];
                // NaN scores rank worst, as in the confidence gate.
                let cs = if c.score.is_nan() {
                    f64::NEG_INFINITY
                } else {
                    c.score
                };
                let better = match best {
                    None => true,
                    Some(b) => {
                        let bs = if b.score.is_nan() {
                            f64::NEG_INFINITY
                        } else {
                            b.score
                        };
                        cs > bs || (cs == bs && c.rt < b.rt)
                    }
                };
                if better {
                    best = Some(c);
                }
                k += 1;
            }
            let b = best.expect("run group is non-empty");
            donors.push(run);
            best_rt.push((run, b.rt));
            if b.rt.is_finite() {
                obs.push(Observation {
                    precursor: id,
                    run,
                    rt: b.rt,
                });
            }
        }

        // Majority CV over all PSMs of the precursor; smallest value on ties.
        let mut cvs: Vec<f32> = group.iter().filter_map(|p| p.faims_cv).collect();
        cvs.sort_by(f32::total_cmp);
        let mut faims_cv: Option<f32> = None;
        let mut best_n = 0usize;
        let mut a = 0;
        while a < cvs.len() {
            let mut b = a;
            while b < cvs.len() && cvs[b] == cvs[a] {
                b += 1;
            }
            if b - a > best_n {
                best_n = b - a;
                faims_cv = Some(cvs[a]);
            }
            a = b;
        }

        entries.push(LibraryEntry {
            precursor: id,
            q: qs[id as usize].expect("entry has a q-value"),
            consensus_rt: f64::NAN,
            faims_cv,
            donors,
            best_rt,
        });
    }

    let alignment = align(&obs, n_prec, n_runs, ap);
    for e in &mut entries {
        e.consensus_rt = alignment
            .consensus
            .get(e.precursor as usize)
            .copied()
            .unwrap_or(f64::NAN);
    }
    entries.sort_by(|a, b| {
        meta[a.precursor as usize]
            .mz
            .total_cmp(&meta[b.precursor as usize].mz)
            .then(a.precursor.cmp(&b.precursor))
    });
    Library {
        entries,
        alignment,
        meta,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn meta(n: usize) -> Vec<PrecursorMeta> {
        (0..n)
            .map(|i| PrecursorMeta {
                key: format!("P{i}/2"),
                mz: 500.0 + i as f64,
                charge: 2,
                envelope: [0.5, 0.3, 0.15, 0.05],
            })
            .collect()
    }

    fn psm(run: u16, precursor: u32, rt: f64, score: f64, is_decoy: bool) -> PsmIn {
        PsmIn {
            run,
            precursor,
            rt,
            score,
            is_decoy,
            faims_cv: None,
        }
    }

    /// 200 strong targets over two runs plus the supplied extra PSMs.
    fn base(extra: Vec<PsmIn>) -> Vec<PsmIn> {
        let mut v = Vec::new();
        for p in 0..200u32 {
            for r in 0..2u16 {
                v.push(psm(r, p, 10.0 + p as f64, 100.0 + p as f64, false));
            }
        }
        v.extend(extra);
        v
    }

    #[test]
    fn decoys_never_enter_library() {
        // Decoy 200 outscores every target; it must still not be an entry or
        // an alignment observation.
        let psms = base(vec![psm(0, 200, 5.0, 1000.0, true)]);
        let lib = build_library(&psms, meta(201), 2, 1.0, &AlignParams::default());
        assert!(lib.entries.iter().all(|e| e.precursor != 200));
        assert_eq!(lib.entries.len(), 200);
        assert!(lib.alignment.consensus[200].is_nan());
    }

    #[test]
    fn entries_sorted_by_mz() {
        let mut m = meta(200);
        // Reverse the m/z order relative to ids.
        for (i, x) in m.iter_mut().enumerate() {
            x.mz = 900.0 - i as f64;
        }
        let lib = build_library(&base(vec![]), m, 2, 0.01, &AlignParams::default());
        assert_eq!(lib.entries.len(), 200);
        let mzs: Vec<f64> = lib
            .entries
            .iter()
            .map(|e| lib.meta[e.precursor as usize].mz)
            .collect();
        assert!(mzs.windows(2).all(|w| w[0] <= w[1]));
        assert_eq!(lib.entries[0].precursor, 199);
    }

    #[test]
    fn faims_cv_is_majority_of_donor_ids() {
        let cv = |run, rt, c: Option<f32>| PsmIn {
            faims_cv: c,
            ..psm(run, 0, rt, 50.0, false)
        };
        let mut psms: Vec<PsmIn> = base(vec![])
            .into_iter()
            .filter(|p| p.precursor != 0)
            .collect();
        psms.extend([
            cv(0, 10.0, Some(-50.0)),
            cv(1, 10.0, Some(-65.0)),
            cv(1, 11.0, Some(-65.0)),
            cv(0, 12.0, None),
        ]);
        let lib = build_library(&psms, meta(200), 2, 1.0, &AlignParams::default());
        let e0 = lib.entries.iter().find(|e| e.precursor == 0).unwrap();
        assert_eq!(e0.faims_cv, Some(-65.0));
        // A tie goes to the smallest value; no CV gives None.
        let mut psms2: Vec<PsmIn> = psms.iter().copied().filter(|p| p.precursor != 0).collect();
        psms2.extend([cv(0, 10.0, Some(-50.0)), cv(1, 10.0, Some(-65.0))]);
        let lib = build_library(&psms2, meta(200), 2, 1.0, &AlignParams::default());
        let e0 = lib.entries.iter().find(|e| e.precursor == 0).unwrap();
        assert_eq!(e0.faims_cv, Some(-65.0));
        let e1 = lib.entries.iter().find(|e| e.precursor == 1).unwrap();
        assert_eq!(e1.faims_cv, None);
    }

    #[test]
    fn best_rt_uses_best_scoring_psm_per_run() {
        let mut psms: Vec<PsmIn> = base(vec![])
            .into_iter()
            .filter(|p| p.precursor != 0)
            .collect();
        psms.extend([
            psm(0, 0, 10.0, 40.0, false),
            psm(0, 0, 12.5, 90.0, false),
            psm(2, 0, 30.0, 60.0, false),
            psm(2, 0, 20.0, 60.0, false),
        ]);
        let lib = build_library(&psms, meta(200), 3, 1.0, &AlignParams::default());
        let e0 = lib.entries.iter().find(|e| e.precursor == 0).unwrap();
        assert_eq!(e0.donors, vec![0, 2]);
        // Later PSM scores higher in run 0; equal scores take the earliest RT.
        assert_eq!(e0.best_rt, vec![(0, 12.5), (2, 20.0)]);
    }

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
