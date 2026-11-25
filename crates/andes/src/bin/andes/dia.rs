pub(crate) mod blocks;
pub(crate) mod interference;
pub(crate) mod joint;
pub(crate) mod llr;
pub(crate) mod pass2;
pub(crate) mod pin;
pub(crate) mod rt;
pub(crate) mod score;
pub(crate) mod study;
pub(crate) mod transfer;

use std::collections::HashMap;
use std::path::Path;
use std::time::Instant;

use model::isotope::averagine_isotope_envelope;
use model::mass::ISOTOPE;
use model::scan::Ms1Scan;
use model::stats::MAD_TO_SD;
use model::Spectrum;
use onnx_infer::predictors::form::{Mass, H2O, PROTON};
use onnx_infer::predictors::spectral::FragmentTolerance;
use quant::Ms1RunIndex;
use rayon::prelude::*;

use crate::dia::interference::removed_keys;
pub(crate) use crate::dia::interference::Interference;
use crate::dia::rt::{read_pin, write_pin, NONE};

pub(crate) const COLUMNS: [&str; 16] = [
    "dia_frag_corr",
    "dia_frag_coelute",
    "dia_apex_offset",
    "dia_apex_frac",
    "dia_ms1_corr",
    "dia_ms1_found",
    "dia_ms1_abs_ppm",
    "dia_ms1_log_int",
    "dia_ms1_iso_cos",
    "dia_frag_log_int",
    "dia_n_frag",
    "dia_support",
    "dia_support_top1",
    "dia_sum_n2",
    "dia_sum_frac2",
    "dia_sum_gain",
];

const HALF_CYCLES: usize = 5;
const TOP_FRAGMENTS: usize = 6;
const MIN_ION_LEN: usize = 3;
const COELUTE_MIN_CORR: f64 = 0.6;
const FRAG_PPM: f64 = 20.0;
const SUM_PPM: f64 = 10.0;
const MS1_PPM: f64 = 10.0;
pub(crate) const FRAGMENT_TOL_PPM: f64 = 7.0;
/// Rows kept per scan in the search, before the pool scores re-rank them down to `--top-n`.
pub(crate) const RERANK_POOL: u32 = 100;
const RERANK_FEATURE: &str = "andesms2_spectral_angle";
/// |scan RT - calibrated andes-rt RT| (minutes) of a PIN row, and its missing flag.
const RERANK_RT_DIFF: &str = "andesrt_rt_diff";
const RERANK_RT_MISSING: &str = "andesrt_missing";
/// Rounds of 3-sigma trimming of the per-file RT spread in the re-rank.
const RERANK_RT_TRIM_ROUNDS: usize = 5;
const RERANK_RT_TRIM_SIGMAS: f64 = 3.0;

fn spectrum_key<'a>(spec_id: &'a str, scan: &str) -> &'a str {
    let p: Vec<&str> = spec_id.rsplitn(4, '_').collect();
    if p.len() == 4 && p[2] == scan {
        p[3]
    } else {
        spec_id.rsplitn(3, '_').nth(2).unwrap_or(spec_id)
    }
}

/// Keeps the best `keep` rows of every scan of the PIN (by RawScore, the MS2 spectral angle and
/// the RT term) and returns the SpecIds of the dropped rows when `want_dropped`; an empty set
/// otherwise, which saves holding every dropped id when no output needs them.
pub(crate) fn rerank_pin(
    pin: &Path,
    keep: usize,
    want_dropped: bool,
) -> Result<std::collections::HashSet<String>, String> {
    use std::io::{BufRead, BufReader, BufWriter, Write};
    let t0 = Instant::now();
    let err = |e: std::io::Error| format!("{}: {e}", pin.display());
    let header = {
        let mut lines = BufReader::new(std::fs::File::open(pin).map_err(err)?).lines();
        lines.next().transpose().map_err(err)?.unwrap_or_default()
    };
    let h: Vec<&str> = header.split('\t').collect();
    let at = |n: &str| h.iter().position(|c| *c == n);
    let (Some(c_scan), Some(c_raw), Some(c_ms2)) =
        (at("ScanNr"), at("RawScore"), at(RERANK_FEATURE))
    else {
        return Ok(Default::default());
    };
    let c_rt = at(RERANK_RT_DIFF).zip(at(RERANK_RT_MISSING));
    let c_label = at("Label");
    let c_run = at("FileName");
    let need = c_scan
        .max(c_raw)
        .max(c_ms2)
        .max(c_rt.map_or(0, |(a, b)| a.max(b)))
        .max(c_label.unwrap_or(0))
        .max(c_run.unwrap_or(0))
        + 1;
    let mut group_of: rustc_hash::FxHashMap<String, u32> = rustc_hash::FxHashMap::default();
    let mut key = String::new();
    let mut run_of: rustc_hash::FxHashMap<String, u32> = rustc_hash::FxHashMap::default();
    let mut group = Vec::new();
    let mut raw = Vec::new();
    let mut ms2 = Vec::new();
    let mut rt_diff = Vec::new();
    let mut target = Vec::new();
    let mut run = Vec::new();
    struct Row<'a> {
        spectrum: &'a str,
        scan: &'a str,
        raw: f64,
        ms2: f64,
        rt_diff: f64,
        target: bool,
        run: &'a str,
    }
    crate::dia::blocks::read_blocks(
        pin,
        |_| Ok(()),
        |_, lines| {
            let parsed: Vec<Option<Row<'_>>> = lines
                .par_iter()
                .map(|l| {
                    let x: Vec<&str> = l.splitn(need + 1, '\t').collect();
                    if x.len() < need {
                        return None;
                    }
                    let (r, m) = (x[c_raw].parse::<f64>(), x[c_ms2].parse::<f64>());
                    let (rt_diff, target, run) = match c_rt {
                        Some((cd, cm)) => {
                            let missing = x[cm].parse::<f64>() != Ok(0.0);
                            (
                                if missing {
                                    f64::NAN
                                } else {
                                    x[cd].parse::<f64>().unwrap_or(f64::NAN)
                                },
                                c_label.is_some_and(|c| x[c] == "1"),
                                c_run.map_or("", |c| x[c]),
                            )
                        }
                        None => (f64::NAN, false, ""),
                    };
                    Some(Row {
                        spectrum: spectrum_key(x[0], x[c_scan]),
                        scan: x[c_scan],
                        raw: r.unwrap_or(f64::NAN),
                        ms2: m.unwrap_or(f64::NAN),
                        rt_diff,
                        target,
                        run,
                    })
                })
                .collect();
            for r in parsed {
                let Some(r) = r else {
                    return Err(format!("{}: short row while re-ranking", pin.display()));
                };
                key.clear();
                key.push_str(r.spectrum);
                key.push('\t');
                key.push_str(r.scan);
                let g = match group_of.get(&key) {
                    Some(&g) => g,
                    None => {
                        let g = group_of.len() as u32;
                        group_of.insert(key.clone(), g);
                        g
                    }
                };
                group.push(g);
                raw.push(r.raw);
                ms2.push(r.ms2);
                if c_rt.is_some() {
                    rt_diff.push(r.rt_diff);
                    target.push(r.target);
                    let id = match run_of.get(r.run) {
                        Some(&id) => id,
                        None => {
                            let id = run_of.len() as u32;
                            run_of.insert(r.run.to_string(), id);
                            id
                        }
                    };
                    run.push(id);
                }
            }
            Ok(())
        },
    )?;
    drop(group_of);
    let n_rows = group.len();
    let z = |v: &[f64]| -> Vec<f64> {
        let f: Vec<f64> = v.iter().copied().filter(|x| x.is_finite()).collect();
        let n = f.len().max(1) as f64;
        let mean = f.iter().sum::<f64>() / n;
        let sd = (f.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n).sqrt();
        v.iter()
            .map(|x| {
                if x.is_finite() && sd > 0.0 {
                    (x - mean) / sd
                } else {
                    0.0
                }
            })
            .collect()
    };
    let (zr, zm) = (z(&raw), z(&ms2));
    drop((raw, ms2));
    let mut score: Vec<f64> = zr.iter().zip(&zm).map(|(a, b)| a + b).collect();
    drop((zr, zm));
    if !rt_diff.is_empty() {
        let mut names: Vec<(String, u32)> = run_of.into_iter().collect();
        names.sort_unstable();
        add_rt_term(&mut score, &rt_diff, &target, &run, &group, &names);
    }
    drop((rt_diff, target, run));
    let mut order: Vec<usize> = (0..n_rows).collect();
    order.par_sort_unstable_by(|&a, &b| {
        group[a]
            .cmp(&group[b])
            .then(score[b].total_cmp(&score[a]))
            .then(a.cmp(&b))
    });
    let mut keep_row = vec![false; n_rows];
    let mut i = 0;
    while i < order.len() {
        let mut j = i;
        while j < order.len() && group[order[j]] == group[order[i]] {
            if j - i < keep {
                keep_row[order[j]] = true;
            }
            j += 1;
        }
        i = j;
    }
    drop((order, score, group));
    let kept = keep_row.iter().filter(|k| **k).count();
    let mut dropped =
        std::collections::HashSet::with_capacity(if want_dropped { n_rows - kept } else { 0 });
    let tmp = pin.with_file_name(format!(
        "{}.rerank.tmp",
        pin.file_name()
            .map_or_else(|| "pin".into(), |n| n.to_string_lossy())
    ));
    let io = |e: std::io::Error| format!("{}: {e}", tmp.display());
    let mut write_kept = || -> Result<(), String> {
        let mut w = BufWriter::new(std::fs::File::create(&tmp).map_err(io)?);
        w.write_all(header.as_bytes()).map_err(io)?;
        w.write_all(b"\n").map_err(io)?;
        let mut rows = 0usize;
        crate::dia::blocks::read_blocks(
            pin,
            |_| Ok(()),
            |first, lines| {
                if first + lines.len() > n_rows {
                    return Err(format!("{}: changed while re-ranking", pin.display()));
                }
                let parts: Vec<(String, Vec<String>)> = lines
                    .par_chunks(crate::dia::blocks::LINES_PER_TASK)
                    .enumerate()
                    .map(|(c, chunk)| {
                        let mut buf = String::new();
                        let mut gone = Vec::new();
                        for (j, l) in chunk.iter().enumerate() {
                            if keep_row[first + c * crate::dia::blocks::LINES_PER_TASK + j] {
                                buf.push_str(l);
                                buf.push('\n');
                            } else if want_dropped {
                                gone.push(l.split('\t').next().unwrap_or_default().to_string());
                            }
                        }
                        (buf, gone)
                    })
                    .collect();
                for (buf, gone) in parts {
                    w.write_all(buf.as_bytes()).map_err(io)?;
                    dropped.extend(gone);
                }
                rows = first + lines.len();
                Ok(())
            },
        )?;
        if rows != n_rows {
            return Err(format!("{}: changed while re-ranking", pin.display()));
        }
        w.flush().map_err(io)?;
        drop(w);
        std::fs::rename(&tmp, pin)
            .map_err(|e| format!("{} -> {}: {e}", tmp.display(), pin.display()))
    };
    if let Err(e) = write_kept() {
        let _ = std::fs::remove_file(&tmp);
        return Err(e);
    }
    log::info!(
        "kept {kept} of {n_rows} PIN rows, {keep} per scan by RawScore + {RERANK_FEATURE} [{:.2}s]",
        t0.elapsed().as_secs_f64()
    );
    Ok(dropped)
}

/// Adds the RT term of [`search::dia`] to the re-rank score of every row, with the
/// spread fitted per run (`run`, the id of the PIN's FileName in `names`, sorted by name; one
/// run without it): the target winners (best row of a scan `group` by the re-rank score) in the
/// better half of their run by that score give `MAD_TO_SD * median |RT difference|`,
/// re-estimated on the winners within `RERANK_RT_TRIM_SIGMAS` sigmas for
/// `RERANK_RT_TRIM_ROUNDS` rounds. A run with fewer than
/// [`search::dia::MIN_ANCHORS`] winners keeps the score without the term.
fn add_rt_term(
    score: &mut [f64],
    rt_diff: &[f64],
    target: &[bool],
    run: &[u32],
    group: &[u32],
    names: &[(String, u32)],
) {
    let n_groups = group.iter().max().map_or(0, |g| *g as usize + 1);
    let mut best = vec![usize::MAX; n_groups];
    for (i, &g) in group.iter().enumerate() {
        let b = &mut best[g as usize];
        if *b == usize::MAX || score[i].total_cmp(&score[*b]).is_gt() {
            *b = i;
        }
    }
    let mut winners: Vec<Vec<usize>> = vec![Vec::new(); names.len()];
    for &i in &best {
        if i != usize::MAX && target[i] && rt_diff[i].is_finite() {
            winners[run[i] as usize].push(i);
        }
    }
    drop(best);
    let mut sigma = vec![f64::NAN; names.len()];
    for (f, id) in names {
        let w = &mut winners[*id as usize];
        if w.is_empty() {
            continue;
        }
        w.sort_by(|&a, &b| score[b].total_cmp(&score[a]).then(a.cmp(&b)));
        w.truncate(w.len().div_ceil(2));
        let d: Vec<f64> = w.iter().map(|&i| rt_diff[i].abs()).collect();
        let mut s = f64::INFINITY;
        for _ in 0..RERANK_RT_TRIM_ROUNDS {
            let mut kept: Vec<f64> = d
                .iter()
                .copied()
                .filter(|v| *v < RERANK_RT_TRIM_SIGMAS * s)
                .collect();
            if kept.len() < search::dia::MIN_ANCHORS {
                break;
            }
            s = MAD_TO_SD * median_mut(&mut kept);
        }
        if s.is_finite() && s > 0.0 {
            log::info!(
                target: "dia",
                "re-rank RT term: run {f:?}: sigma {s:.3} min from {} target winners",
                w.len()
            );
            sigma[*id as usize] = s;
        }
    }
    for i in 0..score.len() {
        let s = sigma[run[i] as usize];
        if s.is_finite() {
            score[i] += search::dia::rt_penalty(rt_diff[i], s);
        }
    }
}

fn median_mut(v: &mut [f64]) -> f64 {
    v.sort_by(f64::total_cmp);
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2]
    } else {
        0.5 * (v[n / 2 - 1] + v[n / 2])
    }
}

type WindowKey = (usize, i64, Option<(u32, u32)>);

pub(crate) fn windows(
    spectra: &[Spectrum],
    run_ends: &[usize],
) -> (Vec<(u32, u32)>, Vec<Vec<u32>>) {
    let run_of = |i: usize| run_ends.partition_point(|&e| e <= i);
    let mut by_window: std::collections::BTreeMap<WindowKey, Vec<u32>> =
        std::collections::BTreeMap::new();
    for (i, s) in spectra.iter().enumerate() {
        let t = s.isolation_target_mz.unwrap_or(s.precursor_mz);
        by_window
            .entry((run_of(i), (t * 100.0).round() as i64, s.mobility_window))
            .or_default()
            .push(i as u32);
    }
    let rt = |i: u32| spectra[i as usize].rt_seconds.unwrap_or(0.0);
    let mut slot = vec![(0u32, 0u32); spectra.len()];
    let mut windows: Vec<Vec<u32>> = Vec::with_capacity(by_window.len());
    for (_, mut v) in by_window {
        v.sort_by(|&a, &b| rt(a).total_cmp(&rt(b)));
        for (p, &i) in v.iter().enumerate() {
            slot[i as usize] = (windows.len() as u32, p as u32);
        }
        windows.push(v);
    }
    (slot, windows)
}

pub(crate) fn neighbours<'a>(
    spectra: &'a [Spectrum],
    slot: &[(u32, u32)],
    windows: &[Vec<u32>],
    spec: usize,
) -> (Vec<&'a [(f64, f32)]>, usize) {
    let (w, p) = slot[spec];
    let win = &windows[w as usize];
    let (p, lo) = (p as usize, (p as usize).saturating_sub(HALF_CYCLES));
    let hi = (p + HALF_CYCLES + 1).min(win.len());
    let scans = win[lo..hi]
        .iter()
        .map(|&i| spectra[i as usize].peaks.as_slice())
        .collect();
    (scans, p - lo)
}

/// Adds the co-elution `COLUMNS` to every row of the PIN; `ms1` is [`ms1_indexes`] of the run's
/// MS1 scans, `n_ms1` their number before indexing (for the log).
pub(crate) fn add_to_pin(
    pin: &Path,
    spectra: &[Spectrum],
    run_ends: &[usize],
    ms1: &[Option<Ms1RunIndex>],
    n_ms1: usize,
) -> Result<(), String> {
    let t0 = Instant::now();
    let cols: Vec<String> = COLUMNS.iter().map(|c| c.to_string()).collect();
    let table = read_pin(pin, spectra, &cols)?;
    let run_of = |i: usize| run_ends.partition_point(|&e| e <= i);
    let rt = |i: u32| spectra[i as usize].rt_seconds.unwrap_or(0.0);
    let (slot, windows) = windows(spectra, run_ends);

    let mut best_of_spec: rustc_hash::FxHashMap<u32, f64> = rustc_hash::FxHashMap::default();
    for r in &table.rows {
        if r.spec != NONE {
            let b = best_of_spec.entry(r.spec).or_insert(f64::NEG_INFINITY);
            *b = b.max(r.score);
        }
    }
    let mut by_precursor: rustc_hash::FxHashMap<(u32, u8, u32), Vec<(u32, bool)>> =
        rustc_hash::FxHashMap::default();
    for r in &table.rows {
        if r.spec == NONE || r.form == NONE {
            continue;
        }
        let (w, p) = slot[r.spec as usize];
        by_precursor
            .entry((r.form, r.charge, w))
            .or_default()
            .push((p, r.score >= best_of_spec[&r.spec]));
    }
    let support = |form: u32, charge: u8, w: u32, p: u32| -> (f64, f64) {
        let mut seen: Vec<(u32, bool)> = by_precursor
            .get(&(form, charge, w))
            .map(|v| {
                v.iter()
                    .copied()
                    .filter(|&(q, _)| q != p && q.abs_diff(p) as usize <= HALF_CYCLES)
                    .collect()
            })
            .unwrap_or_default();
        seen.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(b.1.cmp(&a.1)));
        seen.dedup_by_key(|s| s.0);
        (
            seen.len() as f64,
            seen.iter().filter(|s| s.1).count() as f64,
        )
    };

    let feats: Vec<[f64; COLUMNS.len()]> = table
        .rows
        .par_iter()
        .map(|r| {
            if r.spec == NONE || r.form == NONE || r.charge == 0 {
                return [0.0; COLUMNS.len()];
            }
            let (w, p) = slot[r.spec as usize];
            let win = &windows[w as usize];
            let (p, lo) = (p as usize, (p as usize).saturating_sub(HALF_CYCLES));
            let hi = (p + HALF_CYCLES + 1).min(win.len());
            let scans: Vec<&[(f64, f32)]> = win[lo..hi]
                .iter()
                .map(|&i| spectra[i as usize].peaks.as_slice())
                .collect();
            let rts: Vec<f64> = win[lo..hi].iter().map(|&i| rt(i)).collect();
            let run = run_of(r.spec as usize);
            let mut f = row_features(
                &scans,
                p - lo,
                &rts,
                &table.forms[r.form as usize].1,
                r.charge,
                ms1.get(run).and_then(Option::as_ref),
            );
            let (s_any, s_top) = support(r.form, r.charge, w, p as u32);
            f[11] = s_any;
            f[12] = s_top;
            let nb = |k: Option<usize>| k.map(|k| spectra[win[k] as usize].peaks.as_slice());
            let (n2, frac2, gain) = apex_sum(
                spectra[r.spec as usize].peaks.as_slice(),
                nb(p.checked_sub(1)),
                nb((p + 1 < win.len()).then_some(p + 1)),
                &table.forms[r.form as usize].1,
                r.charge,
            );
            f[13] = n2;
            f[14] = frac2;
            f[15] = gain;
            f
        })
        .collect();
    let values: Vec<Vec<f64>> = (0..COLUMNS.len())
        .into_par_iter()
        .map(|c| feats.iter().map(|f| f[c]).collect())
        .collect();
    write_pin(pin, &table, &cols, &values)?;
    log::info!(
        "added {} co-elution columns to {} PIN rows ({} isolation windows, {n_ms1} MS1 scans) in {:.2} s",
        COLUMNS.len(),
        table.rows.len(),
        windows.len(),
        t0.elapsed().as_secs_f64()
    );
    Ok(())
}

/// The MS1 index of each run's scans, `None` for a run without MS1 scans; built once and shared
/// by the pool scores, pass 2 and the co-elution columns.
pub(crate) fn ms1_indexes(ms1: Vec<Vec<Ms1Scan>>) -> Vec<Option<Ms1RunIndex>> {
    ms1.into_iter()
        .map(|s| (!s.is_empty()).then(|| Ms1RunIndex::new(s)))
        .collect()
}

pub(crate) const EVIDENCE_COLUMNS: [&str; 5] = [
    "dia_sum3_angle",
    "dia_sum3_explained",
    "dia_pf_corr",
    "dia_pf_clean",
    "dia_proj_angle",
];
const EVIDENCE_TOP: usize = 8;

fn angle(p: &[f64], o: &[f64]) -> f64 {
    let dot: f64 = p.iter().zip(o).map(|(a, b)| a * b).sum();
    let (np, no) = (
        p.iter().map(|v| v * v).sum::<f64>().sqrt(),
        o.iter().map(|v| v * v).sum::<f64>().sqrt(),
    );
    if np > 0.0 && no > 0.0 {
        1.0 - 2.0 * (dot / (np * no)).clamp(-1.0, 1.0).acos() / std::f64::consts::PI
    } else {
        0.0
    }
}

#[cfg(test)]
pub(crate) fn evidence(
    pred: &[f32],
    ions: &[(bool, u8)],
    mass: &Mass,
    charge: u8,
    scans: &[&[(f64, f32)]],
    center: usize,
    tol: FragmentTolerance,
) -> [f64; EVIDENCE_COLUMNS.len()] {
    let (frags, traces) = fragment_evidence(pred, ions, mass, charge, scans, tol);
    evidence_of(&frags, &traces, scans.len(), center)
}

/// The predicted fragments of a row (`(m/z, predicted intensity)`, in the order of
/// [`onnx_infer::predictors::spectral::ms2_features_with`]) and their [`fragment_traces`] over
/// `scans`.
pub(crate) fn fragment_evidence(
    pred: &[f32],
    ions: &[(bool, u8)],
    mass: &Mass,
    charge: u8,
    scans: &[&[(f64, f32)]],
    tol: FragmentTolerance,
) -> (Vec<(f64, f64)>, Vec<f64>) {
    let n = mass.res.len();
    let total: f64 = mass.res.iter().sum::<f64>() + mass.nterm;
    let mut frags: Vec<(f64, f64)> = Vec::new();
    let mut prefix = mass.nterm;
    for i in 0..n.saturating_sub(1) {
        prefix += mass.res[i];
        for (k, &(b, z)) in ions.iter().enumerate() {
            if z > charge {
                continue;
            }
            let neutral = if b { prefix } else { total - prefix + H2O };
            let p = pred.get(i * ions.len() + k).copied().unwrap_or(0.0) as f64;
            frags.push(((neutral + z as f64 * PROTON) / z as f64, p.max(0.0)));
        }
    }
    let traces = fragment_traces(&frags, scans, tol);
    (frags, traces)
}

/// The `EVIDENCE_COLUMNS` of a row from its [`fragment_evidence`]: `ns` scans per fragment
/// trace, the row's own scan at `center`.
pub(crate) fn evidence_of(
    frags: &[(f64, f64)],
    traces: &[f64],
    ns: usize,
    center: usize,
) -> [f64; EVIDENCE_COLUMNS.len()] {
    let trace = |k: usize| &traces[k * ns..(k + 1) * ns];
    let p: Vec<f64> = frags.iter().map(|f| f.1).collect();
    let lo = center.saturating_sub(1);
    let hi = (center + 2).min(ns);
    let apex = (lo..hi)
        .max_by(|&x, &y| {
            let at = |s: usize| (0..frags.len()).map(|k| trace(k)[s] * p[k]).sum::<f64>();
            at(x)
                .total_cmp(&at(y))
                .then(y.abs_diff(center).cmp(&x.abs_diff(center)))
        })
        .unwrap_or(center);
    let (a, b) = (apex.saturating_sub(1), (apex + 2).min(ns));
    let sum3: Vec<f64> = (0..frags.len())
        .map(|k| trace(k)[a..b].iter().sum())
        .collect();
    let sp: f64 = p.iter().sum();
    let sum3_angle = angle(&p, &sum3);
    let sum3_explained = if sp > 0.0 {
        p.iter()
            .zip(&sum3)
            .filter(|(_, o)| **o > 0.0)
            .map(|(v, _)| v)
            .sum::<f64>()
            / sp
    } else {
        0.0
    };

    let mut top: Vec<usize> = (0..frags.len()).collect();
    top.sort_by(|&x, &y| p[y].total_cmp(&p[x]));
    top.truncate(EVIDENCE_TOP);
    let norm: Vec<Vec<f64>> = top
        .iter()
        .map(|&k| {
            let m = trace(k).iter().cloned().fold(0.0, f64::max);
            trace(k)
                .iter()
                .map(|v| if m > 0.0 { v / m } else { 0.0 })
                .collect()
        })
        .collect();
    let mut v: Vec<f64> = Vec::with_capacity(norm.len());
    let profile: Vec<f64> = (0..ns)
        .map(|s| {
            v.clear();
            v.extend(norm.iter().map(|t| t[s]));
            v.sort_by(f64::total_cmp);
            match v.len() {
                0 => 0.0,
                n if n % 2 == 1 => v[n / 2],
                n => (v[n / 2 - 1] + v[n / 2]) / 2.0,
            }
        })
        .collect();
    let pn = profile.len().max(1) as f64;
    let pm = profile.iter().sum::<f64>() / pn;
    let pvar: f64 = profile.iter().map(|v| (v - pm) * (v - pm)).sum();
    let corrs: Vec<f64> = top.iter().map(|&k| pearson(trace(k), &profile)).collect();
    let pf_corr = corrs.iter().sum::<f64>() / corrs.len().max(1) as f64;
    let pf_clean =
        corrs.iter().filter(|&&c| c >= COELUTE_MIN_CORR).count() as f64 / top.len().max(1) as f64;
    let proj: Vec<f64> = top
        .iter()
        .map(|&k| {
            let tm = trace(k).iter().sum::<f64>() / pn;
            let cov: f64 = trace(k)
                .iter()
                .zip(&profile)
                .map(|(x, y)| (x - tm) * (y - pm))
                .sum();
            if pvar > 0.0 {
                (cov / pvar).max(0.0)
            } else {
                0.0
            }
        })
        .collect();
    let ptop: Vec<f64> = top.iter().map(|&k| p[k]).collect();
    [
        sum3_angle,
        sum3_explained,
        pf_corr,
        pf_clean,
        angle(&ptop, &proj),
    ]
}

/// The highest peak intensity within `tol` of each fragment's m/z in each scan, fragment-major
/// (`frags.len()` rows of `scans.len()` values; 0 without a peak). Each scan's peaks are walked
/// once, in the fragments' m/z order, instead of a binary search per fragment; a window whose
/// lower edge falls below the previous one is searched afresh, so every fragment sees the peaks
/// a binary search would give it.
fn fragment_traces(
    frags: &[(f64, f64)],
    scans: &[&[(f64, f32)]],
    tol: FragmentTolerance,
) -> Vec<f64> {
    let ns = scans.len();
    let mut order: Vec<usize> = (0..frags.len()).collect();
    order.sort_by(|&a, &b| frags[a].0.total_cmp(&frags[b].0));
    let mut out = vec![0.0f64; frags.len() * ns];
    for (s, peaks) in scans.iter().enumerate() {
        let mut lo = 0usize;
        let mut prev = f64::NEG_INFINITY;
        for &k in &order {
            let mz = frags[k].0;
            let w = tol.window(mz);
            let edge = mz - w;
            if edge >= prev {
                while lo < peaks.len() && peaks[lo].0 < edge {
                    lo += 1;
                }
            } else {
                lo = peaks.partition_point(|p| p.0 < edge);
            }
            prev = edge;
            out[k * ns + s] = peaks[lo..]
                .iter()
                .take_while(|p| p.0 <= mz + w)
                .fold(0.0, |m, p| f64::max(m, p.1 as f64));
        }
    }
    out
}

fn best_peak(peaks: &[(f64, f32)], mz: f64, tol_ppm: f64) -> Option<(f64, f32)> {
    let tol = mz * tol_ppm * 1e-6;
    let lo = peaks.partition_point(|p| p.0 < mz - tol);
    peaks[lo..]
        .iter()
        .take_while(|p| p.0 <= mz + tol)
        .copied()
        .max_by(|a, b| a.1.total_cmp(&b.1))
}

fn pearson(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len() as f64;
    if a.len() < 3 {
        return 0.0;
    }
    let (ma, mb) = (a.iter().sum::<f64>() / n, b.iter().sum::<f64>() / n);
    let (mut sab, mut saa, mut sbb) = (0.0, 0.0, 0.0);
    for (x, y) in a.iter().zip(b) {
        sab += (x - ma) * (y - mb);
        saa += (x - ma) * (x - ma);
        sbb += (y - mb) * (y - mb);
    }
    if saa <= 0.0 || sbb <= 0.0 {
        0.0
    } else {
        sab / (saa * sbb).sqrt()
    }
}

fn row_features(
    scans: &[&[(f64, f32)]],
    center: usize,
    rts: &[f64],
    mass: &Mass,
    charge: u8,
    ms1: Option<&Ms1RunIndex>,
) -> [f64; COLUMNS.len()] {
    let n = mass.res.len();
    let total: f64 = mass.res.iter().sum::<f64>() + mass.nterm;
    let mut ions: Vec<f64> = Vec::new();
    let mut prefix = mass.nterm;
    for i in 0..n.saturating_sub(1) {
        prefix += mass.res[i];
        let (b_len, y_len) = (i + 1, n - i - 1);
        for z in 1..=(if charge >= 3 { 2 } else { 1 }) {
            let zf = z as f64;
            if b_len >= MIN_ION_LEN {
                ions.push((prefix + zf * PROTON) / zf);
            }
            if y_len >= MIN_ION_LEN {
                ions.push((total - prefix + H2O + zf * PROTON) / zf);
            }
        }
    }
    let here = scans[center];
    let mut matched: Vec<(f64, f32)> = ions
        .iter()
        .filter_map(|&mz| best_peak(here, mz, FRAG_PPM).map(|(_, i)| (mz, i)))
        .collect();
    let n_frag = matched.len() as f64;
    matched.sort_by(|a, b| b.1.total_cmp(&a.1));
    matched.truncate(TOP_FRAGMENTS);
    let traces: Vec<Vec<f64>> = matched
        .iter()
        .map(|&(mz, _)| {
            scans
                .iter()
                .map(|s| best_peak(s, mz, FRAG_PPM).map_or(0.0, |p| p.1 as f64))
                .collect()
        })
        .collect();
    let summed: Vec<f64> = (0..scans.len())
        .map(|k| traces.iter().map(|t| t[k]).sum())
        .collect();
    let corrs: Vec<f64> = traces
        .iter()
        .map(|t| {
            let rest: Vec<f64> = summed.iter().zip(t).map(|(s, x)| s - x).collect();
            pearson(t, &rest)
        })
        .collect();
    let frag_corr = if corrs.len() >= 2 {
        corrs.iter().sum::<f64>() / corrs.len() as f64
    } else {
        0.0
    };
    let coelute = corrs.iter().filter(|&&c| c >= COELUTE_MIN_CORR).count() as f64;
    let (apex, apex_v) =
        summed
            .iter()
            .enumerate()
            .fold((center, 0.0), |b, (k, &v)| if v > b.1 { (k, v) } else { b });
    let apex_offset = apex.abs_diff(center) as f64;
    let apex_frac = if apex_v > 0.0 {
        summed[center] / apex_v
    } else {
        0.0
    };

    let z = charge as f64;
    let mono = (total + H2O) / z + PROTON;
    let (mut ms1_corr, mut found, mut abs_ppm, mut log_int, mut iso_cos) =
        (0.0, 0.0, MS1_PPM, 0.0, 0.0);
    if let Some(ms1) = ms1 {
        let trace: Vec<f64> = rts
            .iter()
            .map(|&t| {
                ms1.nearest_scan(t)
                    .and_then(|i| best_peak(&ms1.scan(i).peaks, mono, MS1_PPM))
                    .map_or(0.0, |p| p.1 as f64)
            })
            .collect();
        ms1_corr = pearson(&trace, &summed);
        if let Some(i) = ms1.nearest_scan(rts[center]) {
            let peaks = &ms1.scan(i).peaks;
            if let Some((mz, int)) = best_peak(peaks, mono, MS1_PPM) {
                found = 1.0;
                abs_ppm = ((mz - mono) / mono * 1e6).abs();
                log_int = (int as f64).max(1.0).log10();
                let obs: Vec<f64> = (0..3)
                    .map(|k| {
                        best_peak(peaks, mono + k as f64 * ISOTOPE / z, MS1_PPM)
                            .map_or(0.0, |p| p.1 as f64)
                    })
                    .collect();
                let theo = averagine_isotope_envelope(total + H2O, 3);
                let dot: f64 = obs.iter().zip(&theo).map(|(a, b)| a * b).sum();
                let (na, nb) = (
                    obs.iter().map(|v| v * v).sum::<f64>().sqrt(),
                    theo.iter().map(|v| v * v).sum::<f64>().sqrt(),
                );
                iso_cos = if na > 0.0 && nb > 0.0 {
                    dot / (na * nb)
                } else {
                    0.0
                };
            }
        }
    }
    let frag_log_int = traces
        .iter()
        .map(|t| t[center])
        .sum::<f64>()
        .max(1.0)
        .log10();
    [
        frag_corr,
        coelute,
        apex_offset,
        apex_frac,
        ms1_corr,
        found,
        abs_ppm,
        log_int,
        iso_cos,
        frag_log_int,
        n_frag,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ]
}

fn apex_sum(
    here: &[(f64, f32)],
    prev: Option<&[(f64, f32)]>,
    next: Option<&[(f64, f32)]>,
    mass: &Mass,
    charge: u8,
) -> (f64, f64, f64) {
    let n = mass.res.len();
    let total: f64 = mass.res.iter().sum::<f64>() + mass.nterm;
    let mut prefix = mass.nterm;
    let (mut ions, mut n2, mut gain) = (0usize, 0usize, 0usize);
    let has =
        |s: Option<&[(f64, f32)]>, mz: f64| s.is_some_and(|s| best_peak(s, mz, SUM_PPM).is_some());
    for i in 0..n.saturating_sub(1) {
        prefix += mass.res[i];
        let (b_len, y_len) = (i + 1, n - i - 1);
        for z in 1..=charge.min(2) {
            let zf = z as f64;
            let b = (b_len >= MIN_ION_LEN).then(|| (prefix + zf * PROTON) / zf);
            let y = (y_len >= MIN_ION_LEN).then(|| (total - prefix + H2O + zf * PROTON) / zf);
            for mz in [b, y].into_iter().flatten() {
                ions += 1;
                let (a, b, c) = (has(Some(here), mz), has(prev, mz), has(next, mz));
                if a as u8 + b as u8 + c as u8 >= 2 {
                    n2 += 1;
                }
                if !a && (b || c) {
                    gain += 1;
                }
            }
        }
    }
    let frac = if ions > 0 {
        n2 as f64 / ions as f64
    } else {
        0.0
    };
    (n2 as f64, frac, gain as f64)
}

pub(crate) const RT_REFINED: &str = "dia_rt_refined_err";
const RT_MIN_ANCHORS: usize = 300;

pub(crate) fn refine_rt(
    pin: &Path,
    spectra: &[Spectrum],
    accepted: &(dyn Fn(&str) -> bool + Sync),
) -> Result<bool, String> {
    use model_train::gbdt::train::{train_gbdt_regression, RegressionDataset, TrainParams};

    let t0 = Instant::now();
    let table = read_pin(pin, spectra, &[RT_REFINED.to_string()])?;
    let cols: std::cell::OnceCell<Option<(usize, usize, usize)>> = std::cell::OnceCell::new();
    let mut files: Vec<String> = Vec::new();
    let mut run = Vec::with_capacity(table.rows.len());
    let mut anchor = Vec::with_capacity(table.rows.len());
    let mut pred = Vec::with_capacity(table.rows.len());
    crate::dia::blocks::read_blocks(
        pin,
        |header| {
            let h: Vec<&str> = header.split('\t').collect();
            let at = |n: &str| h.iter().position(|c| *c == n);
            let _ = cols.set(
                match (at("SpecId"), at("FileName"), at("andesrt_predicted_rt")) {
                    (Some(a), Some(b), Some(c)) => Some((a, b, c)),
                    _ => None,
                },
            );
            Ok(())
        },
        |_, lines| {
            let Some(&Some((c_id, c_file, c_pred))) = cols.get() else {
                return Ok(());
            };
            let parsed: Vec<(&str, bool, f64)> = lines
                .par_iter()
                .map(|line| {
                    let x: Vec<&str> = line.split('\t').collect();
                    (
                        x[c_file],
                        accepted(x[c_id]),
                        x[c_pred].parse::<f64>().unwrap_or(f64::NAN),
                    )
                })
                .collect();
            for (f, a, p) in parsed {
                run.push(files.iter().position(|g| g == f).unwrap_or_else(|| {
                    files.push(f.to_string());
                    files.len() - 1
                }));
                anchor.push(a);
                pred.push(p);
            }
            Ok(())
        },
    )?;
    if !matches!(cols.get(), Some(Some(_))) {
        return Ok(false);
    }
    if run.len() != table.rows.len() {
        return Err(format!("{}: rows changed while refining RT", pin.display()));
    }
    let obs: Vec<f64> = table
        .rows
        .iter()
        .map(|r| {
            (r.spec != NONE)
                .then(|| spectra[r.spec as usize].rt_seconds)
                .flatten()
                .map_or(f64::NAN, |t| t / 60.0)
        })
        .collect();
    const AA: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";
    let form_x: Vec<Vec<f32>> = table
        .forms
        .iter()
        .map(|(f, _)| {
            let mut v = vec![0f32; AA.len() + 2];
            for b in f.seq.bytes() {
                if let Some(i) = AA.iter().position(|&a| a == b) {
                    v[i] += 1.0;
                }
            }
            v[AA.len()] = f.seq.len() as f32;
            v[AA.len() + 1] = f
                .mods
                .iter()
                .filter(|(_, m)| matches!(m, onnx_infer::predictors::form::Mod::Oxidation))
                .count() as f32;
            v
        })
        .collect();
    let nf = AA.len() + 3;
    let names: Vec<String> = (0..nf).map(|i| format!("f{i}")).collect();
    let row_x = |i: usize| {
        let mut v = form_x[table.rows[i].form as usize].clone();
        v.push(pred[i] as f32);
        v
    };
    let usable = |i: usize| {
        let r = &table.rows[i];
        r.form != NONE && pred[i].is_finite() && obs[i].is_finite()
    };
    let fold = |i: usize| (table.rows[i].form as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15) >> 40;

    let mut err = vec![f64::NAN; table.rows.len()];
    let mut report = Vec::new();
    let mut skipped = Vec::new();
    let mut by_run: Vec<Vec<usize>> = vec![Vec::new(); files.len()];
    for i in (0..table.rows.len()).filter(|&i| usable(i)) {
        by_run[run[i]].push(i);
    }
    for (r, rows) in by_run.into_iter().enumerate() {
        let anchors: Vec<usize> = rows
            .iter()
            .copied()
            .filter(|&i| anchor[i] && !table.rows[i].decoy)
            .collect();
        if anchors.len() < RT_MIN_ANCHORS {
            report.push(format!("{}: {} anchors, skipped", files[r], anchors.len()));
            skipped.push(files[r].to_string());
            continue;
        }
        let models: Vec<_> = (0..3u64)
            .into_par_iter()
            .map(|k| {
                let train: Vec<usize> = anchors
                    .iter()
                    .copied()
                    .filter(|&i| fold(i) % 3 != k)
                    .collect();
                let ds = RegressionDataset {
                    x: train.iter().flat_map(|&i| row_x(i)).collect(),
                    y: train.iter().map(|&i| (obs[i] - pred[i]) as f32).collect(),
                    groups: train.iter().map(|&i| table.rows[i].form).collect(),
                    n_features: nf,
                    feature_names: names.clone(),
                };
                let params = TrainParams {
                    n_rounds: 200,
                    max_depth: 3,
                    min_hessian: 20.0,
                    ..TrainParams::default()
                };
                train_gbdt_regression(&ds, &params, 42).ok()
            })
            .collect();
        let errs: Vec<f64> = rows
            .par_iter()
            .map(|&i| {
                let shift = models[(fold(i) % 3) as usize]
                    .as_ref()
                    .map_or(0.0, |m| m.predict_value(&row_x(i)) as f64);
                (obs[i] - pred[i] - shift).abs()
            })
            .collect();
        for (&i, e) in rows.iter().zip(errs) {
            err[i] = e;
        }
        let mut a: Vec<f64> = anchors.iter().map(|&i| err[i]).collect();
        let mut b: Vec<f64> = anchors.iter().map(|&i| (obs[i] - pred[i]).abs()).collect();
        a.sort_by(f64::total_cmp);
        b.sort_by(f64::total_cmp);
        report.push(format!(
            "{}: {} anchors, median |error| {:.2} -> {:.2} min",
            files[r],
            anchors.len(),
            b[b.len() / 2],
            a[a.len() / 2]
        ));
    }
    if err.iter().all(|e| !e.is_finite()) {
        log::warn!(
            target: "dia-rt",
            "RT refinement skipped, no dia_rt_refined_err column: {}",
            report.join("; ")
        );
        return Ok(false);
    }
    let mut fin: Vec<f64> = err.iter().copied().filter(|e| e.is_finite()).collect();
    fin.par_sort_unstable_by(f64::total_cmp);
    let med = fin[fin.len() / 2];
    for e in &mut err {
        if !e.is_finite() {
            *e = med;
        }
    }
    write_pin(pin, &table, &[RT_REFINED.to_string()], &[err])?;
    if !skipped.is_empty() {
        log::warn!(
            target: "dia-rt",
            "RT refinement skipped for {} of {} runs (fewer than {} anchors); their dia_rt_refined_err is the median: {}",
            skipped.len(),
            files.len(),
            RT_MIN_ANCHORS,
            skipped.join(", ")
        );
    }
    log::info!(
        target: "dia-rt",
        "{} [{:.2}s]",
        report.join("; "),
        t0.elapsed().as_secs_f64()
    );
    Ok(true)
}

pub(crate) fn default_train_fdr(args: &mut Vec<String>) {
    if !args
        .iter()
        .any(|a| a == "-F" || a == "--trainFDR" || a.starts_with("--trainFDR="))
    {
        args.extend(["--trainFDR".to_string(), "0.05".to_string()]);
    }
}

/// Precursor counts of [`precursor_q`], the best `(score, is decoy)` of every precursor, and the
/// `(FileName, Peptide, charge)` of the precursors, of both labels, that interference removal
/// dropped.
pub(crate) type PrecursorQ = (
    (usize, usize, usize),
    Vec<(f64, bool)>,
    std::collections::HashSet<(String, String, u8)>,
);

/// Reports precursor-level q-values for DIA: target-decoy competition over the best
/// Percolator score of each precursor (run, peptide form, charge), with q = (D + 1) / T
/// and monotone in the score. A row receives the q of the threshold at its own score, so
/// keeping rows with q <= x keeps the precursors whose best row passes, at precursor FDR x.
/// Rows Percolator did not report keep no entry. Returns the precursor counts (targets,
/// decoys, target precursors at q <= 0.01) and the best `(score, is decoy)` of every precursor,
/// for comparing the run's precursor-level target-decoy curve at other q-values. A precursor
/// PIN (one with the `score::PRESENT_COLUMN` column) is keyed by its Peptide string and
/// charge columns and needs no spectra. With `interference` (--dia-interference) and a PSM
/// PIN, [`interference_removed`] drops precursors from the competition first; their target rows
/// get q = 1.
pub(crate) fn precursor_q(
    pin: &Path,
    spectra: &[Spectrum],
    targets: &mut HashMap<String, output::PercolatorPsm>,
    decoys: &HashMap<String, output::PercolatorPsm>,
    interference: Option<Interference>,
) -> Result<PrecursorQ, String> {
    let header = {
        use std::io::BufRead;
        let f = std::fs::File::open(pin).map_err(|e| format!("{}: {e}", pin.display()))?;
        std::io::BufReader::new(f)
            .lines()
            .find(|l| !matches!(l, Ok(l) if l.is_empty()))
            .transpose()
            .map_err(|e| format!("{}: {e}", pin.display()))?
            .unwrap_or_default()
    };
    let h: Vec<&str> = header.split('\t').collect();
    let at = |n: &str| {
        h.iter()
            .position(|c| *c == n)
            .ok_or_else(|| format!("{}: PIN has no {n} column", pin.display()))
    };
    let (c_id, c_pep, c_label) = (at("SpecId")?, at("Peptide")?, at("Label")?);
    let c_file = h.iter().position(|c| *c == "FileName");
    let charge_cols: Vec<(usize, u8)> = h
        .iter()
        .enumerate()
        .filter_map(|(i, c)| Some((i, c.strip_prefix("charge")?.parse().ok()?)))
        .collect();
    let table = if h.contains(&crate::dia::score::PRESENT_COLUMN) {
        None
    } else {
        Some(read_pin(pin, spectra, &[])?)
    };
    let mut keys: HashMap<(String, String, u8), u32> = HashMap::new();
    let mut rows: Vec<(u32, f64, bool, String)> = Vec::new();
    let mut raw: Vec<(String, String, u8)> = Vec::new();
    let mut line_of: Vec<usize> = Vec::new();
    type Parsed<'a> = Option<(bool, String, u8, &'a str, f64, &'a str, &'a str, usize)>;
    let targets_ro: &HashMap<String, output::PercolatorPsm> = targets;
    crate::dia::blocks::read_blocks(
        pin,
        |_| Ok(()),
        |first, lines| {
            let parsed: Vec<Result<Parsed<'_>, String>> = lines
                .par_iter()
                .enumerate()
                .map(|(j, line)| {
                    let i = first + j;
                    let x: Vec<&str> = line.split('\t').collect();
                    if x.len() <= c_pep.max(c_label).max(c_id) {
                        return Err(format!("{}: short PIN row {}", pin.display(), i + 2));
                    }
                    let (decoy, peptide, charge) = match &table {
                        Some(t) => {
                            let r = t.rows.get(i).ok_or_else(|| {
                                format!("{}: row {} not read", pin.display(), i + 2)
                            })?;
                            let peptide = if r.form == NONE {
                                x[c_pep].to_string()
                            } else {
                                format!("#{}", r.form)
                            };
                            (r.decoy, peptide, r.charge)
                        }
                        None => (
                            x[c_label].trim() == "-1",
                            x[c_pep].to_string(),
                            charge_cols
                                .iter()
                                .find(|(c, _)| x.get(*c) == Some(&"1"))
                                .map_or(0, |c| c.1),
                        ),
                    };
                    let id = x[c_id];
                    let psm = if decoy {
                        decoys.get(id)
                    } else {
                        targets_ro.get(id)
                    };
                    let Some(score) = psm.and_then(|p| p.score) else {
                        return Ok(None);
                    };
                    Ok(Some((
                        decoy,
                        peptide,
                        charge,
                        id,
                        score,
                        c_file.map_or("", |i| x[i]),
                        x[c_pep],
                        i,
                    )))
                })
                .collect();
            for p in parsed {
                let Some((decoy, peptide, charge, id, score, run, raw_pep, i)) = p? else {
                    continue;
                };
                let n = keys.len() as u32;
                let k = *keys.entry((run.to_string(), peptide, charge)).or_insert(n);
                raw.push((run.to_string(), raw_pep.to_string(), charge));
                rows.push((k, score, decoy, id.to_string()));
                line_of.push(i);
            }
            Ok(())
        },
    )?;
    let removed = match (&table, interference) {
        (Some(t), Some(choice)) => {
            let q0 = precursor_q_values(keys.len(), &rows);
            removed_keys(&keys, &rows, &line_of, &q0, t, spectra, choice)
        }
        _ => vec![false; keys.len()],
    };
    for (k, _, decoy, id) in &rows {
        if removed[*k as usize] && !*decoy {
            if let Some(p) = targets.get_mut(id) {
                p.q_value = 1.0;
            }
        }
    }
    let removed_precursors = rows
        .iter()
        .zip(raw)
        .filter(|(r, _)| removed[r.0 as usize])
        .map(|(_, n)| n)
        .collect();
    rows.retain(|r| !removed[r.0 as usize]);
    let q = precursor_q_values(keys.len(), &rows);
    let n_targets = q.iter().filter(|p| p.is_some_and(|p| !p.1)).count();
    let n_decoys = q.iter().filter(|p| p.is_some_and(|p| p.1)).count();
    let passing = q.iter().flatten().filter(|p| !p.1 && p.2 <= 0.01).count();
    let curve = q_curve(&q);
    for (_, score, decoy, id) in &rows {
        if *decoy {
            continue;
        }
        if let Some(p) = targets.get_mut(id) {
            p.q_value = curve_at(&curve, *score);
        }
    }
    let best = q.iter().flatten().map(|p| (p.0, p.1)).collect();
    Ok(((n_targets, n_decoys, passing), best, removed_precursors))
}

/// Best score of each precursor and its precursor-level q: `(best score, is decoy, q)`,
/// `None` for a key no row reached. Ties in score count decoys first.
fn precursor_q_values(
    n_keys: usize,
    rows: &[(u32, f64, bool, String)],
) -> Vec<Option<(f64, bool, f64)>> {
    let mut best: Vec<Option<(f64, bool)>> = vec![None; n_keys];
    for (k, score, decoy, _) in rows {
        let b = &mut best[*k as usize];
        if b.is_none_or(|(s, _)| *score > s) {
            *b = Some((*score, *decoy));
        }
    }
    let mut order: Vec<usize> = (0..n_keys).filter(|&i| best[i].is_some()).collect();
    order.sort_by(|&a, &b| {
        let (sa, da) = best[a].unwrap_or((f64::NEG_INFINITY, true));
        let (sb, db) = best[b].unwrap_or((f64::NEG_INFINITY, true));
        sb.total_cmp(&sa).then(db.cmp(&da))
    });
    let mut q = vec![1.0f64; n_keys];
    let (mut t, mut d) = (0usize, 0usize);
    let mut i = 0;
    while i < order.len() {
        let s = best[order[i]].map_or(f64::NEG_INFINITY, |b| b.0);
        let mut j = i;
        while j < order.len() && best[order[j]].map_or(f64::NEG_INFINITY, |b| b.0) == s {
            if best[order[j]].is_some_and(|b| b.1) {
                d += 1;
            } else {
                t += 1;
            }
            j += 1;
        }
        let fdr = ((d + 1) as f64 / t.max(1) as f64).min(1.0);
        for &k in &order[i..j] {
            q[k] = fdr;
        }
        i = j;
    }
    let mut running = 1.0f64;
    for &k in order.iter().rev() {
        running = running.min(q[k]);
        q[k] = running;
    }
    best.iter()
        .zip(q)
        .map(|(b, qv)| b.map(|(s, d)| (s, d, qv)))
        .collect()
}

/// Score thresholds in descending order with the precursor q at each.
fn q_curve(q: &[Option<(f64, bool, f64)>]) -> Vec<(f64, f64)> {
    let mut c: Vec<(f64, f64)> = q.iter().flatten().map(|p| (p.0, p.2)).collect();
    c.sort_by(|a, b| b.0.total_cmp(&a.0).then(a.1.total_cmp(&b.1)));
    c
}

/// The precursor q a row with this score gets: the q of the lowest precursor best score
/// that is still at or above it; 1 below every precursor.
fn curve_at(curve: &[(f64, f64)], score: f64) -> f64 {
    let n = curve.partition_point(|p| p.0 >= score);
    if n == 0 {
        1.0
    } else {
        curve[n - 1].1
    }
}

/// Precursor-level q of every row computed within its class only: one target-decoy competition
/// per class (rows of class `false` against each other, rows of class `true` against each
/// other), with the rules of [`precursor_q`]. Returns the q of each row, 1 for a row whose
/// class holds no precursor at or above its score.
pub(crate) fn classwise_q_values(
    n_keys: usize,
    rows: &[(u32, f64, bool, String)],
    class: &[bool],
) -> Vec<f64> {
    let mut q = vec![1.0; rows.len()];
    for c in [false, true] {
        let idx: Vec<usize> = (0..rows.len()).filter(|&i| class[i] == c).collect();
        let sub: Vec<(u32, f64, bool, String)> = idx.iter().map(|&i| rows[i].clone()).collect();
        let curve = q_curve(&precursor_q_values(n_keys, &sub));
        for &i in &idx {
            q[i] = curve_at(&curve, rows[i].1);
        }
    }
    q
}

/// Replaces the targets' q-values of a precursor PIN by the class-wise q
/// ([`classwise_q_values`], classes by `class_col` 0 or 1) and counts the target precursors at
/// q <= 0.01 by class under the pooled q they had and under the class-wise q. Returns
/// `[[pooled class 0, pooled class 1], [class-wise class 0, class-wise class 1]]`.
pub(crate) fn precursor_q_by_class(
    pin: &Path,
    class_col: &str,
    targets: &mut HashMap<String, output::PercolatorPsm>,
    decoys: &HashMap<String, output::PercolatorPsm>,
) -> Result<[[usize; 2]; 2], String> {
    let err = |e: &dyn std::fmt::Display| format!("{}: {e}", pin.display());
    let f = std::fs::File::open(pin).map_err(|e| err(&e))?;
    let mut lines = std::io::BufRead::lines(std::io::BufReader::new(f));
    let header = lines
        .next()
        .ok_or_else(|| err(&"empty PIN"))?
        .map_err(|e| err(&e))?;
    let h: Vec<&str> = header.split('\t').collect();
    let at = |n: &str| {
        h.iter()
            .position(|c| *c == n)
            .ok_or_else(|| format!("{}: PIN has no {n} column", pin.display()))
    };
    let (c_id, c_pep, c_label, c_class) = (at("SpecId")?, at("Peptide")?, at("Label")?, at(class_col)?);
    let charge_cols: Vec<(usize, u8)> = h
        .iter()
        .enumerate()
        .filter_map(|(i, c)| Some((i, c.strip_prefix("charge")?.parse().ok()?)))
        .collect();
    let mut keys: HashMap<(String, u8), u32> = HashMap::new();
    let mut rows: Vec<(u32, f64, bool, String)> = Vec::new();
    let mut class: Vec<bool> = Vec::new();
    for line in lines {
        let line = line.map_err(|e| err(&e))?;
        if line.is_empty() {
            continue;
        }
        let x: Vec<&str> = line.split('\t').collect();
        if x.len() <= c_pep.max(c_label).max(c_id).max(c_class) {
            return Err(err(&"short PIN row"));
        }
        let decoy = x[c_label].trim() == "-1";
        let psm = if decoy { decoys.get(x[c_id]) } else { targets.get(x[c_id]) };
        let Some(score) = psm.and_then(|p| p.score) else {
            continue;
        };
        let charge = charge_cols
            .iter()
            .find(|(c, _)| x.get(*c) == Some(&"1"))
            .map_or(0, |c| c.1);
        let n = keys.len() as u32;
        let k = *keys.entry((x[c_pep].to_string(), charge)).or_insert(n);
        rows.push((k, score, decoy, x[c_id].to_string()));
        class.push(x[c_class].trim().parse::<f32>().unwrap_or(0.0) > 0.5);
    }
    let q = classwise_q_values(keys.len(), &rows, &class);
    let mut counts = [[0usize; 2]; 2];
    for (i, (_, _, decoy, id)) in rows.iter().enumerate() {
        if *decoy {
            continue;
        }
        let c = usize::from(class[i]);
        if let Some(p) = targets.get_mut(id) {
            counts[0][c] += usize::from(p.q_value <= 0.01);
            counts[1][c] += usize::from(q[i] <= 0.01);
            p.q_value = q[i];
        }
    }
    Ok(counts)
}

#[cfg(test)]
#[path = "dia/tests/dia.rs"]
mod tests;
