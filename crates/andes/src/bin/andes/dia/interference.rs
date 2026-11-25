//! DIA-NN-style interference removal before the precursor-level q (`--dia-interference`).

use std::collections::HashMap;

use model::Spectrum;
use onnx_infer::predictors::form::{Mass, H2O, PROTON};

use crate::dia::rt::NONE;

/// Precursor q at or below which a precursor takes part in interference removal.
const INTERFERENCE_POOL_Q: f64 = 0.2;
/// Fragment m/z bin (Th) of interference removal.
const INTERFERENCE_BIN: f64 = 0.02;

/// b/y fragment m/z bins of a precursor: charge 1, and 2 for a precursor of charge 3 or more.
fn fragment_bins(mass: &Mass, charge: u8) -> Vec<i64> {
    let total: f64 = mass.res.iter().sum::<f64>() + mass.nterm;
    let mut prefix = mass.nterm;
    let mut out = Vec::new();
    for r in &mass.res[..mass.res.len().saturating_sub(1)] {
        prefix += r;
        for z in 1..=if charge >= 3 { 2u8 } else { 1 } {
            let zf = z as f64;
            for neutral in [prefix, total - prefix + H2O] {
                out.push(((neutral + zf * PROTON) / zf / INTERFERENCE_BIN).round() as i64);
            }
        }
    }
    out.sort_unstable();
    out.dedup();
    out
}

/// DIA-NN-style interference removal, blind to labels: a precursor is removed when a
/// higher-scoring precursor of the same run, whose best row is in the same isolation window, an
/// overlapping mobility slice and within `rule.1` seconds of RT, shares at least `rule.0` of its
/// b/y fragment m/z. Returns, per precursor key, whether it is removed.
pub(super) fn interference_removed(
    windows: &mut HashMap<(u32, i64), Vec<Coeluting>>,
    n_keys: usize,
    rule: (f64, f64),
) -> Vec<bool> {
    let mut removed = vec![false; n_keys];
    for ps in windows.values_mut() {
        let hits = interfered(ps, rule);
        for (p, hit) in ps.iter().zip(hits) {
            removed[p.key] = hit;
        }
    }
    removed
}

/// The best row of each precursor at q <= [`INTERFERENCE_POOL_Q`], grouped by run and isolation
/// window, for [`interference_removed`] and [`estimate_interference`].
pub(super) fn coeluting_windows(
    rows: &[(u32, f64, bool, String)],
    line_of: &[usize],
    q0: &[Option<(f64, bool, f64)>],
    run_of: &[u32],
    table: &crate::dia::rt::PinTable,
    spectra: &[Spectrum],
) -> HashMap<(u32, i64), Vec<Coeluting>> {
    let n = q0.len();
    let mut best_row: Vec<Option<usize>> = vec![None; n];
    for (i, (k, s, _, _)) in rows.iter().enumerate() {
        let b = &mut best_row[*k as usize];
        if b.is_none_or(|b| *s > rows[b].1) {
            *b = Some(i);
        }
    }
    let mut by_window: HashMap<(u32, i64), Vec<Coeluting>> = HashMap::new();
    for key in 0..n {
        let (Some(i), Some((score, _, q))) = (best_row[key], q0[key]) else {
            continue;
        };
        if q > INTERFERENCE_POOL_Q {
            continue;
        }
        let r = &table.rows[line_of[i]];
        if r.spec == NONE || r.form == NONE {
            continue;
        }
        let s = &spectra[r.spec as usize];
        let Some(rt) = s.rt_seconds else { continue };
        let form = &table.forms[r.form as usize];
        let frags = fragment_bins(&form.1, r.charge);
        if frags.is_empty() {
            continue;
        }
        let w = (s.isolation_target_mz.unwrap_or(s.precursor_mz) * 100.0).round() as i64;
        by_window
            .entry((run_of[key], w))
            .or_default()
            .push(Coeluting {
                key,
                score,
                rt,
                mob: s.mobility_window.unwrap_or((0, u32::MAX)),
                seq: form.0.seq.clone(),
                frags,
            });
    }
    by_window
}

/// The best row of a precursor in one isolation window of one run, for [`interfered`].
pub(super) struct Coeluting {
    key: usize,
    score: f64,
    rt: f64,
    mob: (u32, u32),
    /// Bare sequence: forms of one peptide share fragments by construction.
    seq: String,
    /// Sorted, distinct fragment m/z bins.
    frags: Vec<i64>,
}

/// Whether each precursor of one window (sorted here by RT) has a higher-scoring one within
/// `rule.1` s of RT, in an overlapping mobility slice, holding at least `rule.0` of its fragments.
/// Forms of the same sequence never remove each other, as in [`estimate_interference`].
fn interfered(ps: &mut [Coeluting], rule: (f64, f64)) -> Vec<bool> {
    ps.sort_by(|a, b| a.rt.total_cmp(&b.rt));
    (0..ps.len())
        .map(|a| {
            let p = &ps[a];
            let lo = ps.partition_point(|x| x.rt < p.rt - rule.1);
            ps[lo..]
                .iter()
                .take_while(|x| x.rt <= p.rt + rule.1)
                .any(|x| {
                    if x.score <= p.score
                        || x.mob.1 <= p.mob.0
                        || p.mob.1 <= x.mob.0
                        || x.seq == p.seq
                    {
                        return false;
                    }
                    let shared = p.frags.iter().filter(|f| x.frags.binary_search(f).is_ok());
                    shared.count() as f64 >= rule.0 * p.frags.len() as f64
                })
        })
        .collect()
}

/// --dia-interference: fixed `F:DRT`, or estimated from the run (`auto`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum Interference {
    Fixed(f64, f64),
    Auto,
}

/// Pairs further apart in RT than this (s) cannot co-elute: their fragment share is the null.
const NULL_MIN_DRT: f64 = 60.0;
/// Pairs are formed up to this RT gap (s).
const NULL_MAX_DRT: f64 = 200.0;
/// |dRT| bin (s) and span (s) over which the excess sharing is traced.
const DRT_BIN: f64 = 0.5;
const DRT_SPAN: f64 = 20.0;
/// `auto` removes nothing unless the high-share rate of pairs at the same RT exceeds the null
/// rate by this factor: below it, co-eluting precursors share fragments little more than chance.
const MIN_EXCESS_OVER_NULL: f64 = 5.0;

/// An estimated interference rule and the evidence behind it.
pub(crate) struct InterferenceEstimate {
    /// `None`: the run shows no co-elution signal to act on, and nothing is removed.
    pub(crate) rule: Option<(f64, f64)>,
    /// Pairs 60-200 s apart, the null.
    null_pairs: usize,
    /// Null 0.999 quantile of the share: a pair at or above it is a hit.
    s_ref: f64,
    /// Null rate of hits, and the hit rate above it in the first |dRT| bin.
    null_rate: f64,
    excess0: f64,
    /// Why removal is off, if it is.
    off: Option<&'static str>,
}

impl std::fmt::Display for InterferenceEstimate {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let rule = match self.rule {
            Some((fr, drt)) => format!("{fr}:{drt}"),
            None => format!("off, {}", self.off.unwrap_or("no rule")),
        };
        write!(
            f,
            "{rule} ({} null pairs, hit share {:.3}, hit rate {:.4} over a null {:.4} at 0 s)",
            self.null_pairs, self.s_ref, self.excess0, self.null_rate
        )
    }
}

/// Estimates `F:DRT` from the run, without labels, over pairs of pool precursors in one window
/// with overlapping mobility slices (share = shared fragment bins over the lower-scoring one's).
/// Pairs of one bare sequence (charge states, modified forms) are skipped: they share fragments
/// at any RT gap, so they would raise the null's top share without being chance.
///
/// - Pairs [`NULL_MIN_DRT`]-[`NULL_MAX_DRT`] s apart cannot co-elute; their 0.999 share quantile
///   is `s_ref`, and their rate of shares >= `s_ref` the null rate.
/// - DRT: the hit rate (share >= `s_ref`) in 0.5 s bins of |dRT|, less the null rate, falls as
///   pairs stop co-eluting. DRT is the start of the first two non-empty bins at or below 5% of
///   the first bin's excess, at least 1 s. Empty bins (gaps on the acquisition cycle) are skipped.
/// - F: the smallest F in 0.20, 0.25 .. 0.80 at which the pairs within DRT that chance alone
///   would put at share >= F (the null rate at F) are at most 1% of those observed.
///
/// Removal is off (no rule) with fewer than 1000 null pairs or 20 hits in the first bin, an
/// excess there under [`MIN_EXCESS_OVER_NULL`] times the null rate, an excess that does not fade
/// within 20 s, or no F that meets the 1% bound.
pub(super) fn estimate_interference(
    windows: &HashMap<(u32, i64), Vec<Coeluting>>,
) -> InterferenceEstimate {
    let mut pairs: Vec<(f64, f64)> = Vec::new();
    for ps in windows.values() {
        let mut ix: Vec<usize> = (0..ps.len()).collect();
        ix.sort_by(|&a, &b| ps[a].rt.total_cmp(&ps[b].rt));
        for (n, &a) in ix.iter().enumerate() {
            for &c in &ix[n + 1..] {
                let (p, x) = (&ps[a], &ps[c]);
                let dt = x.rt - p.rt;
                if dt > NULL_MAX_DRT {
                    break;
                }
                if x.mob.1 <= p.mob.0 || p.mob.1 <= x.mob.0 || x.seq == p.seq {
                    continue;
                }
                let low = if p.score < x.score { p } else { x };
                let shared = p
                    .frags
                    .iter()
                    .filter(|f| x.frags.binary_search(f).is_ok())
                    .count();
                pairs.push((dt, shared as f64 / low.frags.len() as f64));
            }
        }
    }
    estimate_from_pairs(&pairs)
}

/// The rule of [`estimate_interference`] from `(|dRT|, share)` pairs.
fn estimate_from_pairs(pairs: &[(f64, f64)]) -> InterferenceEstimate {
    let mut e = InterferenceEstimate {
        rule: None,
        null_pairs: 0,
        s_ref: f64::NAN,
        null_rate: f64::NAN,
        excess0: f64::NAN,
        off: None,
    };
    let mut null: Vec<f64> = pairs
        .iter()
        .filter(|p| p.0 >= NULL_MIN_DRT)
        .map(|p| p.1)
        .collect();
    e.null_pairs = null.len();
    if null.len() < 1000 {
        e.off = Some("fewer than 1000 pairs 60-200 s apart");
        return e;
    }
    null.sort_by(f64::total_cmp);
    let quantile = |q: f64| {
        let h = (null.len() - 1) as f64 * q;
        let (lo, hi) = (h.floor() as usize, h.ceil() as usize);
        null[lo] + (null[hi] - null[lo]) * (h - lo as f64)
    };
    let null_at =
        |s: f64| (null.len() - null.partition_point(|&x| x < s)) as f64 / null.len() as f64;
    e.s_ref = quantile(0.999);
    e.null_rate = null_at(e.s_ref);
    let nb = (DRT_SPAN / DRT_BIN) as usize;
    let (mut n, mut hits) = (vec![0usize; nb], vec![0usize; nb]);
    for &(dt, s) in pairs {
        let k = (dt / DRT_BIN) as usize;
        if k < nb {
            n[k] += 1;
            hits[k] += usize::from(s >= e.s_ref);
        }
    }
    let excess: Vec<f64> = (0..nb)
        .map(|k| hits[k] as f64 / n[k].max(1) as f64 - e.null_rate)
        .collect();
    e.excess0 = excess[0];
    if hits[0] < 20 {
        e.off = Some("fewer than 20 high-share pairs at the same RT");
        return e;
    }
    if excess[0] < MIN_EXCESS_OVER_NULL * e.null_rate {
        e.off = Some("co-eluting pairs share fragments little more than chance");
        return e;
    }
    let cut = 0.05 * excess[0];
    let filled: Vec<usize> = (1..nb).filter(|&k| n[k] > 0).collect();
    let Some(k) = filled
        .windows(2)
        .find(|w| excess[w[0]] <= cut && excess[w[1]] <= cut)
        .map(|w| w[0])
    else {
        e.off = Some("the excess share does not fade within 20 s");
        return e;
    };
    let drt = (k as f64 * DRT_BIN).max(1.0);
    let near: Vec<f64> = pairs.iter().filter(|p| p.0 <= drt).map(|p| p.1).collect();
    let f = (0..=12)
        .map(|i| ((0.2 + 0.05 * i as f64) * 100.0).round() / 100.0)
        .find(|&f| {
            let observed = near.iter().filter(|&&s| s >= f).count() as f64;
            observed > 0.0 && null_at(f) * near.len() as f64 <= 0.01 * observed
        });
    match f {
        Some(f) => e.rule = Some((f, drt)),
        None => e.off = Some("no F up to 0.8 keeps chance pairs under 1%"),
    }
    e
}

#[cfg(test)]
#[path = "tests/interference.rs"]
mod tests;
