use std::collections::HashMap;
use std::path::Path;
use std::time::Instant;

use model::isotope::averagine_isotope_envelope;
use model::mass::ISOTOPE;
use model::scan::Ms1Scan;
use model::Spectrum;
use onnx_infer::predictors::form::{Mass, H2O, PROTON};
use quant::Ms1RunIndex;
use rayon::prelude::*;

use crate::predicted::{read_pin, write_pin, NONE};

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
pub(crate) const RERANK_POOL: u32 = 50;
const RERANK_FEATURE: &str = "andesms2_spectral_angle";

fn spectrum_key<'a>(spec_id: &'a str, scan: &str) -> &'a str {
    let p: Vec<&str> = spec_id.rsplitn(4, '_').collect();
    if p.len() == 4 && p[2] == scan {
        p[3]
    } else {
        spec_id.rsplitn(3, '_').nth(2).unwrap_or(spec_id)
    }
}

pub(crate) fn rerank_pin(
    pin: &Path,
    keep: usize,
) -> Result<std::collections::HashSet<String>, String> {
    let t0 = Instant::now();
    let text = std::fs::read_to_string(pin).map_err(|e| format!("{}: {e}", pin.display()))?;
    let mut lines = text.lines();
    let header = lines.next().unwrap_or_default();
    let h: Vec<&str> = header.split('\t').collect();
    let at = |n: &str| h.iter().position(|c| *c == n);
    let (Some(c_scan), Some(c_raw), Some(c_ms2)) =
        (at("ScanNr"), at("RawScore"), at(RERANK_FEATURE))
    else {
        return Ok(Default::default());
    };
    let need = c_scan.max(c_raw).max(c_ms2) + 1;
    let rows: Vec<&str> = lines.filter(|l| !l.is_empty()).collect();
    let mut file = Vec::with_capacity(rows.len());
    let mut scan = Vec::with_capacity(rows.len());
    let mut raw = Vec::with_capacity(rows.len());
    let mut ms2 = Vec::with_capacity(rows.len());
    for l in &rows {
        let x: Vec<&str> = l.splitn(need + 1, '\t').collect();
        if x.len() < need {
            return Err(format!("{}: short row while re-ranking", pin.display()));
        }
        let (r, m) = (x[c_raw].parse::<f64>(), x[c_ms2].parse::<f64>());
        file.push(spectrum_key(x[0], x[c_scan]));
        scan.push(x[c_scan]);
        raw.push(r.unwrap_or(f64::NAN));
        ms2.push(m.unwrap_or(f64::NAN));
    }
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
    let score: Vec<f64> = zr.iter().zip(&zm).map(|(a, b)| a + b).collect();
    let mut order: Vec<usize> = (0..rows.len()).collect();
    order.sort_by(|&a, &b| {
        (file[a], scan[a])
            .cmp(&(file[b], scan[b]))
            .then(score[b].total_cmp(&score[a]))
    });
    let mut keep_row = vec![false; rows.len()];
    let mut i = 0;
    while i < order.len() {
        let mut j = i;
        while j < order.len()
            && (file[order[j]], scan[order[j]]) == (file[order[i]], scan[order[i]])
        {
            if j - i < keep {
                keep_row[order[j]] = true;
            }
            j += 1;
        }
        i = j;
    }
    let kept = keep_row.iter().filter(|k| **k).count();
    let mut dropped = std::collections::HashSet::with_capacity(rows.len() - kept);
    let mut out = String::with_capacity(text.len());
    out.push_str(header);
    out.push('\n');
    for (l, k) in rows.iter().zip(&keep_row) {
        if *k {
            out.push_str(l);
            out.push('\n');
        } else {
            dropped.insert(l.split('\t').next().unwrap_or_default().to_string());
        }
    }
    std::fs::write(pin, out).map_err(|e| format!("{}: {e}", pin.display()))?;
    log::info!(
        "kept {kept} of {} PIN rows, {keep} per scan by RawScore + {RERANK_FEATURE} [{:.2}s]",
        rows.len(),
        t0.elapsed().as_secs_f64()
    );
    Ok(dropped)
}

pub(crate) fn add_to_pin(
    pin: &Path,
    spectra: &[Spectrum],
    run_ends: &[usize],
    ms1: Vec<Vec<Ms1Scan>>,
) -> Result<(), String> {
    let t0 = Instant::now();
    let cols: Vec<String> = COLUMNS.iter().map(|c| c.to_string()).collect();
    let table = read_pin(pin, spectra, &cols)?;
    let run_of = |i: usize| run_ends.partition_point(|&e| e <= i);

    let mut by_window: HashMap<(usize, i64), Vec<u32>> = HashMap::new();
    for (i, s) in spectra.iter().enumerate() {
        let t = s.isolation_target_mz.unwrap_or(s.precursor_mz);
        by_window
            .entry((run_of(i), (t * 100.0).round() as i64))
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
    let n_ms1: usize = ms1.iter().map(Vec::len).sum();
    let ms1: Vec<Option<Ms1RunIndex>> = ms1
        .into_iter()
        .map(|s| (!s.is_empty()).then(|| Ms1RunIndex::new(s)))
        .collect();

    let mut best_of_spec: HashMap<u32, f64> = HashMap::new();
    for r in &table.rows {
        if r.spec != NONE {
            let b = best_of_spec.entry(r.spec).or_insert(f64::NEG_INFINITY);
            *b = b.max(r.score);
        }
    }
    let mut by_precursor: HashMap<(u32, u8, u32), Vec<(u32, bool)>> = HashMap::new();
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
    accepted: &dyn Fn(&str) -> bool,
) -> Result<bool, String> {
    use model_train::gbdt::train::{train_gbdt_regression, RegressionDataset, TrainParams};

    let t0 = Instant::now();
    let table = read_pin(pin, spectra, &[RT_REFINED.to_string()])?;
    let text = std::fs::read_to_string(pin).map_err(|e| format!("{}: {e}", pin.display()))?;
    let mut lines = text.lines().filter(|l| !l.is_empty());
    let h: Vec<&str> = lines.next().unwrap_or_default().split('\t').collect();
    let at = |n: &str| h.iter().position(|c| *c == n);
    let (Some(c_id), Some(c_file), Some(c_pred)) =
        (at("SpecId"), at("FileName"), at("andesrt_predicted_rt"))
    else {
        return Ok(false);
    };
    let mut files: Vec<&str> = Vec::new();
    let mut run = Vec::with_capacity(table.rows.len());
    let mut anchor = Vec::with_capacity(table.rows.len());
    let mut pred = Vec::with_capacity(table.rows.len());
    for line in lines {
        let x: Vec<&str> = line.split('\t').collect();
        let f = x[c_file];
        run.push(files.iter().position(|g| *g == f).unwrap_or_else(|| {
            files.push(f);
            files.len() - 1
        }));
        anchor.push(accepted(x[c_id]));
        pred.push(x[c_pred].parse::<f64>().unwrap_or(f64::NAN));
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
        for k in 0..3u64 {
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
            let model = train_gbdt_regression(&ds, &params, 42).ok();
            for &i in rows.iter().filter(|&&i| fold(i) % 3 == k) {
                let shift = model
                    .as_ref()
                    .map_or(0.0, |m| m.predict_value(&row_x(i)) as f64);
                err[i] = (obs[i] - pred[i] - shift).abs();
            }
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
    fin.sort_by(f64::total_cmp);
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

#[cfg(test)]
mod tests {
    use super::*;
    use onnx_infer::predictors::form::parse_peptide;

    #[test]
    fn coeluting_fragments_and_precursor_score_high() {
        let (_, mass) = parse_peptide("K.PEPTIDEKLR.A").unwrap();
        let z = 2u8;
        let total: f64 = mass.res.iter().sum::<f64>() + mass.nterm;
        let mono = (total + H2O) / z as f64 + PROTON;
        let mut prefix = mass.nterm;
        let mut ys = Vec::new();
        for r in &mass.res[..mass.res.len() - 3] {
            prefix += r;
            ys.push(total - prefix + H2O + PROTON);
        }
        let shape = [0.0, 0.1, 0.4, 0.9, 1.0, 0.7, 0.3, 0.1, 0.0, 0.0, 0.0];
        let rts: Vec<f64> = (0..shape.len()).map(|k| k as f64 * 1.6).collect();
        let mut scans: Vec<Vec<(f64, f32)>> = Vec::new();
        let mut survey: Vec<Ms1Scan> = Vec::new();
        for (k, &h) in shape.iter().enumerate() {
            let mut p: Vec<(f64, f32)> = ys
                .iter()
                .enumerate()
                .map(|(j, &mz)| (mz, (1e5 * h * (1.0 + j as f64)) as f32 + 10.0))
                .collect();
            p.sort_by(|a, b| a.0.total_cmp(&b.0));
            scans.push(p);
            let env = averagine_isotope_envelope(total + H2O, 3);
            survey.push(Ms1Scan {
                rt: rts[k],
                peaks: (0..3)
                    .map(|i| {
                        (
                            mono + i as f64 * ISOTOPE / 2.0,
                            (1e6 * h * env[i]) as f32 + 1.0,
                        )
                    })
                    .collect(),
                ..Default::default()
            });
        }
        let refs: Vec<&[(f64, f32)]> = scans.iter().map(Vec::as_slice).collect();
        let ms1 = Ms1RunIndex::new(survey);
        let f = row_features(&refs, 4, &rts, &mass, z, Some(&ms1));
        assert!(f[0] > 0.95, "frag corr {}", f[0]);
        assert_eq!(f[2], 0.0, "apex offset");
        assert!(f[4] > 0.95, "ms1 corr {}", f[4]);
        assert_eq!(f[5], 1.0);
        assert!(f[8] > 0.99, "iso cos {}", f[8]);

        let (_, other) = parse_peptide("K.GGGGWWWWK.A").unwrap();
        let g = row_features(&refs, 4, &rts, &other, z, Some(&ms1));
        assert_eq!((g[0], g[5], g[10]), (0.0, 0.0, 0.0));

        let n = ys.len() as f64;
        let (n2, _, gain) = apex_sum(refs[4], Some(refs[3]), Some(refs[5]), &mass, z);
        assert_eq!((n2, gain), (n, 0.0));
        assert_eq!(apex_sum(refs[4], None, None, &mass, z).0, 0.0);
        assert_eq!(apex_sum(refs[4], Some(refs[3]), None, &mass, z).0, n);
        assert_eq!(apex_sum(refs[4], None, None, &other, z).0, 0.0);
    }

    #[test]
    fn refine_rt_learns_a_run_offset() {
        let dir = tempfile::tempdir().unwrap();
        let pin = dir.path().join("x.pin");
        let mut text = String::from(
            "SpecId\tLabel\tScanNr\tFileName\tRawScore\tcharge2\tandesrt_predicted_rt\tPeptide\tProteins\n",
        );
        let mut spectra = Vec::new();
        for i in 0..1200usize {
            let core: String = (0..4)
                .map(|k| b"AELKSVTG"[(i >> (3 * k)) & 7] as char)
                .collect();
            let w = i % 3 == 0;
            let pep = format!("K.{core}{}R.A", if w { "W" } else { "F" });
            let pred = 10.0 + (i % 100) as f64 * 0.5;
            let obs = pred + if w { 2.0 } else { 0.0 };
            let scan = i as i32 + 1;
            spectra.push(Spectrum {
                title: format!("t.{scan}"),
                precursor_mz: 500.0,
                precursor_intensity: None,
                precursor_charge: Some(2),
                rt_seconds: Some(obs * 60.0),
                scan: Some(scan),
                peaks: vec![(100.0, 1.0)],
                activation_method: None,
                isolation_lower_offset: None,
                isolation_upper_offset: None,
                isolation_target_mz: None,
            });
            text.push_str(&format!(
                "t.{scan}_{scan}_1\t1\t{scan}\tr\t{}\t1\t{pred}\t{pep}\tP1\n",
                100 - i % 50
            ));
        }
        std::fs::write(&pin, &text).unwrap();
        assert!(refine_rt(&pin, &spectra, &|_| true).unwrap());
        let out = std::fs::read_to_string(&pin).unwrap();
        let h: Vec<&str> = out.lines().next().unwrap().split('\t').collect();
        let c = h.iter().position(|c| *c == RT_REFINED).unwrap();
        assert_eq!(h[c + 1], "Peptide");
        let mut e: Vec<f64> = out
            .lines()
            .skip(1)
            .map(|l| l.split('\t').nth(c).unwrap().parse().unwrap())
            .collect();
        e.sort_by(f64::total_cmp);
        assert!(
            e[e.len() * 9 / 10] < 0.5,
            "90th percentile {}",
            e[e.len() * 9 / 10]
        );

        let bare = text.replace("\tandesrt_predicted_rt", "\tother");
        std::fs::write(&pin, &bare).unwrap();
        assert!(!refine_rt(&pin, &spectra, &|_| true).unwrap());
        assert_eq!(std::fs::read_to_string(&pin).unwrap(), bare);

        for (i, s) in spectra.iter_mut().enumerate() {
            let pred = 10.0 + (i % 100) as f64 * 0.5;
            s.rt_seconds = Some((pred + ((i * 7919) % 13) as f64 * 0.1) * 60.0);
        }
        std::fs::write(&pin, &text).unwrap();
        assert!(refine_rt(&pin, &spectra, &|_| true).unwrap());
    }

    #[test]
    fn train_fdr_defaults_unless_set() {
        let mut a = vec!["--testFDR".to_string(), "0.05".to_string()];
        default_train_fdr(&mut a);
        assert_eq!(a, ["--testFDR", "0.05", "--trainFDR", "0.05"]);
        for set in [vec!["-F", "0.02"], vec!["--trainFDR=0.1"]] {
            let mut b: Vec<String> = set.iter().map(|s| s.to_string()).collect();
            default_train_fdr(&mut b);
            assert_eq!(b, set);
        }
    }

    #[test]
    fn rerank_keeps_top_rows_per_scan_by_score_and_angle() {
        let dir = tempfile::tempdir().unwrap();
        let pin = dir.path().join("x.pin");
        let rows = [
            ("r/c=1 scan=1_1_1_0", 1, 10.0, 0.1),
            ("r/c=1 scan=1_1_2_1", 1, 9.0, 0.9),
            ("r/c=1 scan=1_1_3_2", 1, 8.0, 0.2),
            ("r/c=2 scan=1_1_1_0", 1, 1.0, 0.1),
            ("r/c=2 scan=1_1_2_1", 1, 0.5, 0.1),
            ("r/c=1 scan=2_2_1", 2, 1.0, 0.1),
        ];
        let mut s = String::from("SpecId\tLabel\tScanNr\tRawScore\tandesms2_spectral_angle\n");
        for (id, scan, raw, ms2) in rows {
            s.push_str(&format!("{id}\t1\t{scan}\t{raw}\t{ms2}\n"));
        }
        std::fs::write(&pin, s).unwrap();
        let dropped = rerank_pin(&pin, 1).unwrap();
        let want: std::collections::HashSet<String> = [
            "r/c=1 scan=1_1_1_0",
            "r/c=1 scan=1_1_3_2",
            "r/c=2 scan=1_1_2_1",
        ]
        .map(String::from)
        .into();
        assert_eq!(dropped, want);
        let ids: Vec<String> = std::fs::read_to_string(&pin)
            .unwrap()
            .lines()
            .skip(1)
            .map(|l| l.split('\t').next().unwrap().to_string())
            .collect();
        assert_eq!(
            ids,
            [
                "r/c=1 scan=1_1_2_1",
                "r/c=2 scan=1_1_1_0",
                "r/c=1 scan=2_2_1"
            ]
        );
    }
}
