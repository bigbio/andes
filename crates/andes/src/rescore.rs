#![allow(clippy::needless_range_loop)]
//! Native (Percolator-free) PSM rescorer: a cross-validated, semi-supervised
//! linear model over the PIN features, then target-decoy competition → per-PSM
//! q-value + PEP.
//!
//! ROLE: the self-contained fallback for when a Percolator backend is not
//! available. Percolator (`--rescore`) remains the recommended path and the
//! reference cross-check.
//!
//! METHOD (the Percolator algorithm, `-Y` semantics):
//!
//! * **Competition groups.** Rows sharing the input file, `ScanNr` and the
//!   `ExpMass` text form one group, the same groups Percolator's target-decoy
//!   competition (`-Y`) uses with the PIN `FileName` column. Each group
//!   reports one PSM: its best-scoring row, target or decoy.
//!   q-values are computed over these winners only, so a spectrum never
//!   contributes several rows at one precursor mass.
//! * **Cross-validation.** Rows are folded by spectrum (3 folds); every row is
//!   scored by a model trained on the other folds. Held-out scores are put on
//!   a common scale per fold (the training split's q = 0.01 threshold maps to
//!   0, its median decoy to -1) before the folds are pooled.
//! * **Semi-supervised training.** Within a training split, the initial score
//!   is the single feature that accepts the most target winners at q ≤ 0.01
//!   (over the search's rank-1 row of each group). Each round, the target
//!   winners at q ≤ 0.01 are the positives and every decoy winner a negative;
//!   an L2-regularised logistic regression on standardised features gives the
//!   next score. Training on every target as a positive instead (half of them
//!   false) teaches the model what separates false targets from decoys, which
//!   breaks the target-decoy assumption the q-values rest on.
//! * **PEP.** Isotonic fit of the decoy share among winners against the final
//!   score: PEP = decoys / targets at that score, capped at 1.
//!
//! `ExpMass`, `CalcMass` and `FileName` are not features (as in Percolator).
//! The output maps match `output::run_percolator` with `-Y`: winning targets,
//! and the winning decoys separately; rows that lost their group are in
//! neither.

use std::collections::HashMap;

use output::{spec_id_run as run_of, PercolatorPsm};
use rayon::prelude::*;
use search::tdc::{qvalues, ScoredLabel};

const N_FOLDS: usize = 3;
/// Training rounds per fold; the positives are re-selected every round.
const N_ROUNDS: usize = 10;
/// q-value at which target winners become positive training examples.
const TRAIN_Q: f64 = 0.01;
/// Inverse L2 strength of the logistic regression (loss weight; the penalty
/// is `|w|² / 2`).
const L2_C: f64 = 0.1;
/// A training split with fewer positives or negatives than this keeps the
/// score it has (the initial single feature when no round could be trained).
const MIN_TRAIN_CLASS: usize = 10;
/// PIN columns between `ScanNr` and `Peptide` that are not features.
const NON_FEATURES: [&str; 3] = ["ExpMass", "CalcMass", "FileName"];
const PROB_EPS: f64 = 1e-6;
/// Rows per partial sum in the logistic fit.
const CHUNK: usize = 4096;

/// Numeric feature matrix + labels parsed from a Percolator PIN.
struct PinData {
    spec_ids: Vec<String>,
    peptides: Vec<String>,
    proteins: Vec<String>,
    is_decoy: Vec<bool>,
    /// Hash of each row's spectrum key, for the fold assignment: the rows of
    /// one spectrum (its ranks and co-isolated peptides) share a fold.
    spectrum_hashes: Vec<u64>,
    /// Competition group of each row (input file + `ScanNr` + `ExpMass` text),
    /// densely numbered.
    groups: Vec<u32>,
    /// The first row of its group in PIN order: the search's rank-1 PSM there.
    first_in_group: Vec<bool>,
    /// row-major, `n_rows × n_features`
    x: Vec<f32>,
    n_features: usize,
}

/// Parse a Percolator PIN. Layout:
/// `SpecId \t Label \t ScanNr \t <numeric features…> \t Peptide \t Proteins…`.
/// Feature columns = everything strictly between `ScanNr` and `Peptide`,
/// except `ExpMass`, `CalcMass` and `FileName`.
fn parse_pin(text: &str) -> Result<PinData, String> {
    let mut lines = text.lines();
    let header = lines.next().ok_or("empty PIN")?;
    let cols: Vec<&str> = header.split('\t').collect();
    let find = |name: &str| cols.iter().position(|c| c.eq_ignore_ascii_case(name));
    // Require an explicit PSM identifier column — never silently fall back to
    // column 0, which would let a malformed header join on an arbitrary field.
    let id_i = find("SpecId")
        .or_else(|| find("PSMId"))
        .ok_or("PIN missing SpecId/PSMId column")?;
    let label_i = find("Label").ok_or("PIN missing Label column")?;
    let scan_i = find("ScanNr").ok_or("PIN missing ScanNr column")?;
    let pep_i = find("Peptide").ok_or("PIN missing Peptide column")?;
    let exp_mass_i = find("ExpMass");
    let feat_cols: Vec<usize> = ((scan_i + 1)..pep_i)
        .filter(|&c| !NON_FEATURES.iter().any(|n| cols[c].eq_ignore_ascii_case(n)))
        .collect();
    if feat_cols.is_empty() {
        return Err("PIN has no feature columns between ScanNr and Peptide".into());
    }
    let n_features = feat_cols.len();

    let mut d = PinData {
        spec_ids: Vec::new(),
        peptides: Vec::new(),
        proteins: Vec::new(),
        is_decoy: Vec::new(),
        spectrum_hashes: Vec::new(),
        groups: Vec::new(),
        first_in_group: Vec::new(),
        x: Vec::new(),
        n_features,
    };
    let mut group_ids: HashMap<String, u32> = HashMap::new();
    for (lineno, line) in lines.enumerate() {
        if line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        // A truncated row is corruption, not something to silently drop — a
        // dropped row would quietly undercount the rescored set. (lineno+2 = the
        // 1-based file line, accounting for the header consumed above.)
        if f.len() <= pep_i {
            return Err(format!(
                "PIN line {} is truncated: {} columns, expected > {}",
                lineno + 2,
                f.len(),
                pep_i
            ));
        }
        // Label is FDR-critical — require exactly 1 (target) or -1 (decoy); fail
        // loud on anything else rather than defaulting to target and corrupting
        // the FDR.
        let raw = f[label_i].trim();
        let label: i32 = raw
            .parse()
            .ok()
            .filter(|&l| l == 1 || l == -1)
            .ok_or_else(|| {
                format!(
                    "PIN line {} Label is {raw:?} (expected 1 or -1)",
                    lineno + 2
                )
            })?;
        d.is_decoy.push(label < 0);
        d.spec_ids.push(f[id_i].to_string());
        // ScanNr is part of the spectrum key the CV folds use and of the
        // competition group; a bad value would silently mis-fold or mis-group
        // the spectrum, so fail loud instead of defaulting.
        let scan_raw = f[scan_i].trim();
        let scan: u32 = scan_raw.parse().map_err(|_| {
            format!(
                "PIN line {} ScanNr is {scan_raw:?} (expected integer)",
                lineno + 2
            )
        })?;
        let key = spectrum_key(f[id_i], scan);
        d.spectrum_hashes.push(fnv1a(key));
        // Without an ExpMass column the spectrum is the group. With one, the
        // group is the input file, ScanNr and ExpMass: scan numbers repeat
        // across the files of one search, and spectra of different files must
        // never compete.
        let group_key = match exp_mass_i {
            Some(i) => format!("{}\t{scan}\t{}", run_of(f[id_i]), f[i].trim()),
            None => format!("{key}\t{scan}"),
        };
        let next = group_ids.len() as u32;
        let mut first = false;
        let g = *group_ids.entry(group_key).or_insert_with(|| {
            first = true;
            next
        });
        d.groups.push(g);
        d.first_in_group.push(first);
        d.peptides.push(f[pep_i].to_string());
        d.proteins
            .push(f.get(pep_i + 1..).map(|s| s.join(";")).unwrap_or_default());
        for &ci in &feat_cols {
            let raw = f[ci].trim();
            let v: f32 = raw.parse().map_err(|_| {
                format!(
                    "PIN line {} feature column {ci} is {raw:?} (expected number)",
                    lineno + 2
                )
            })?;
            if !v.is_finite() {
                return Err(format!(
                    "PIN line {} feature column {ci} is non-finite ({v})",
                    lineno + 2
                ));
            }
            d.x.push(v);
        }
    }
    Ok(d)
}

/// A row's spectrum: its SpecId (`{title}_{scan}_{rank}`, multi-row scans add
/// `_{row}`) without that suffix, the first `_{scan}_` followed by one or two
/// numbers and nothing else. The title, not ScanNr, identifies the spectrum:
/// ScanNr is 0 when a native id carries no scan number, and repeats across
/// the files of one search.
fn spectrum_key(spec_id: &str, scan: u32) -> &str {
    let sep = format!("_{scan}_");
    let numbers = |s: &str| {
        let parts: Vec<&str> = s.split('_').collect();
        parts.len() <= 2
            && parts
                .iter()
                .all(|p| !p.is_empty() && p.bytes().all(|b| b.is_ascii_digit()))
    };
    spec_id
        .match_indices(&sep)
        .map(|(i, _)| i)
        .find(|&i| numbers(&spec_id[i + sep.len()..]))
        .map_or(spec_id, |i| &spec_id[..i])
}

fn fnv1a(s: &str) -> u64 {
    s.bytes().fold(0xcbf2_9ce4_8422_2325u64, |h, b| {
        (h ^ b as u64).wrapping_mul(0x0000_0100_0000_01b3)
    })
}

/// Deterministic fold assignment by spectrum (splitmix64 of its key hash).
fn fold_of(spectrum: u64, seed: u64) -> usize {
    let mut z =
        (seed ^ spectrum.wrapping_mul(0x9e37_79b9_7f4a_7c15)).wrapping_add(0x9e37_79b9_7f4a_7c15);
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^= z >> 31;
    (z % N_FOLDS as u64) as usize
}

/// The best-scoring row of each competition group among `rows` (ties: the
/// earlier row), in `rows` order.
fn winners(d: &PinData, rows: &[usize], score: &[f64]) -> Vec<usize> {
    let mut best: HashMap<u32, usize> = HashMap::with_capacity(rows.len());
    for (k, &r) in rows.iter().enumerate() {
        best.entry(d.groups[r])
            .and_modify(|b| {
                if score[k] > score[*b] {
                    *b = k;
                }
            })
            .or_insert(k);
    }
    let mut out: Vec<usize> = best.into_values().collect();
    out.sort_unstable();
    out
}

/// Monotone TDC q-values of `idx` (positions into `score` / `rows`).
fn tdc_q(d: &PinData, rows: &[usize], score: &[f64], idx: &[usize]) -> Vec<f64> {
    let items: Vec<ScoredLabel> = idx
        .iter()
        .map(|&k| ScoredLabel {
            score: score[k] as f32,
            is_decoy: d.is_decoy[rows[k]],
        })
        .collect();
    qvalues(&items)
}

/// A linear score over standardised features: `b + Σ w_j (x_j − μ_j) / σ_j`.
struct LinearModel {
    w: Vec<f64>,
    b: f64,
    mean: Vec<f64>,
    inv_sd: Vec<f64>,
}

impl LinearModel {
    fn score(&self, x: &[f32]) -> f64 {
        let mut s = self.b;
        for j in 0..self.w.len() {
            s += self.w[j] * (x[j] as f64 - self.mean[j]) * self.inv_sd[j];
        }
        s
    }
}

/// L2-regularised logistic regression by Newton's method on standardised
/// features (`z`, row-major, `nf` columns), starting from `w`/`b`.
/// Minimises `C Σ log(1 + e^{-y f}) + |w|²/2`; the intercept is not penalised.
fn fit_logistic(z: &[f64], y: &[bool], nf: usize, w: &mut [f64], b: &mut f64) {
    let p = nf + 1;
    let n = y.len();
    // Partial sums over fixed row chunks, added in chunk order: the result
    // does not depend on the thread count or scheduling.
    let n_chunks = n.div_ceil(CHUNK);
    let objective = |w: &[f64], b: f64| -> f64 {
        let parts: Vec<f64> = (0..n_chunks)
            .into_par_iter()
            .map(|ch| {
                let mut loss = 0.0f64;
                for i in ch * CHUNK..((ch + 1) * CHUNK).min(n) {
                    let row = &z[i * nf..(i + 1) * nf];
                    let f = b + row.iter().zip(w).map(|(a, c)| a * c).sum::<f64>();
                    let m = if y[i] { f } else { -f };
                    // log(1 + e^{-m}), stable for large |m|
                    loss += if m > 0.0 {
                        (-m).exp().ln_1p()
                    } else {
                        -m + m.exp().ln_1p()
                    };
                }
                loss
            })
            .collect();
        L2_C * parts.iter().sum::<f64>() + 0.5 * w.iter().map(|v| v * v).sum::<f64>()
    };
    let mut obj = objective(w, *b);
    for _ in 0..50 {
        // Gradient and Hessian over [w, b]; partial sums per chunk.
        let (wr, br) = (&*w, *b);
        let parts: Vec<(Vec<f64>, Vec<f64>)> = (0..n_chunks)
            .into_par_iter()
            .map(|ch| {
                let mut g = vec![0.0f64; p];
                let mut h = vec![0.0f64; p * p];
                for i in ch * CHUNK..((ch + 1) * CHUNK).min(n) {
                    let row = &z[i * nf..(i + 1) * nf];
                    let f = br + row.iter().zip(wr).map(|(a, c)| a * c).sum::<f64>();
                    let pr = 1.0 / (1.0 + (-f).exp());
                    let r = pr - if y[i] { 1.0 } else { 0.0 };
                    let s = pr * (1.0 - pr);
                    for a in 0..nf {
                        g[a] += r * row[a];
                        let sa = s * row[a];
                        for c in 0..=a {
                            h[a * p + c] += sa * row[c];
                        }
                    }
                    g[nf] += r;
                    for c in 0..nf {
                        h[nf * p + c] += s * row[c];
                    }
                    h[nf * p + nf] += s;
                }
                (g, h)
            })
            .collect();
        let mut g = vec![0.0f64; p];
        let mut h = vec![0.0f64; p * p];
        for (pg, ph) in &parts {
            g.iter_mut().zip(pg).for_each(|(a, c)| *a += c);
            h.iter_mut().zip(ph).for_each(|(a, c)| *a += c);
        }
        for a in 0..p {
            g[a] *= L2_C;
            for c in 0..=a {
                h[a * p + c] *= L2_C;
                h[c * p + a] = h[a * p + c];
            }
        }
        for a in 0..nf {
            g[a] += w[a];
            h[a * p + a] += 1.0;
        }
        // The intercept is unpenalised; a tiny ridge keeps the system
        // positive definite when every row has the same label.
        h[nf * p + nf] += 1e-9;
        let Some(step) = cholesky_solve(&h, &g, p) else {
            return;
        };
        // Damped Newton: halve the step until the objective does not rise.
        let mut t = 1.0;
        loop {
            let nw: Vec<f64> = (0..nf).map(|a| w[a] - t * step[a]).collect();
            let nb = *b - t * step[nf];
            let nobj = objective(&nw, nb);
            if nobj <= obj || t < 1e-4 {
                w.copy_from_slice(&nw);
                *b = nb;
                let done = (obj - nobj).abs() <= 1e-10 * obj.abs().max(1.0);
                obj = nobj;
                if done {
                    return;
                }
                break;
            }
            t *= 0.5;
        }
        if step.iter().map(|v| (v * t).abs()).fold(0.0, f64::max) < 1e-8 {
            return;
        }
    }
}

/// Solve `H x = g` for a symmetric positive definite `H` (row-major `p × p`).
fn cholesky_solve(h: &[f64], g: &[f64], p: usize) -> Option<Vec<f64>> {
    let mut l = vec![0.0f64; p * p];
    for i in 0..p {
        for j in 0..=i {
            let mut s = h[i * p + j];
            for k in 0..j {
                s -= l[i * p + k] * l[j * p + k];
            }
            if i == j {
                if s <= 0.0 || !s.is_finite() {
                    return None;
                }
                l[i * p + i] = s.sqrt();
            } else {
                l[i * p + j] = s / l[j * p + j];
            }
        }
    }
    let mut y = vec![0.0f64; p];
    for i in 0..p {
        let mut s = g[i];
        for k in 0..i {
            s -= l[i * p + k] * y[k];
        }
        y[i] = s / l[i * p + i];
    }
    let mut x = vec![0.0f64; p];
    for i in (0..p).rev() {
        let mut s = y[i];
        for k in (i + 1)..p {
            s -= l[k * p + i] * x[k];
        }
        x[i] = s / l[i * p + i];
    }
    Some(x)
}

/// One fold's model and the affine map that puts its scores on the common
/// scale; `trained` is false when no round could be trained (the initial
/// single-feature score is used).
struct FoldModel {
    model: LinearModel,
    shift: f64,
    scale: f64,
    trained: bool,
}

/// Train on `rows` (the training split): standardise, pick the initial
/// feature, then `N_ROUNDS` semi-supervised rounds.
fn train_fold(d: &PinData, rows: &[usize]) -> FoldModel {
    let nf = d.n_features;
    let n = rows.len();
    let mut mean = vec![0.0f64; nf];
    let mut m2 = vec![0.0f64; nf];
    for &r in rows {
        for j in 0..nf {
            mean[j] += d.x[r * nf + j] as f64;
        }
    }
    mean.iter_mut().for_each(|m| *m /= n.max(1) as f64);
    for &r in rows {
        for j in 0..nf {
            let v = d.x[r * nf + j] as f64 - mean[j];
            m2[j] += v * v;
        }
    }
    let inv_sd: Vec<f64> = m2
        .iter()
        .map(|&v| {
            let sd = (v / n.max(1) as f64).sqrt();
            if sd > 1e-12 {
                1.0 / sd
            } else {
                0.0
            }
        })
        .collect();
    let z: Vec<f64> = rows
        .iter()
        .flat_map(|&r| (0..nf).map(move |j| (r, j)))
        .map(|(r, j)| (d.x[r * nf + j] as f64 - mean[j]) * inv_sd[j])
        .collect();

    // Initial score: the feature (and sign) accepting the most target
    // winners at q ≤ TRAIN_Q over the search's rank-1 rows.
    let rank1: Vec<usize> = (0..n).filter(|&k| d.first_in_group[rows[k]]).collect();
    let (mut best_j, mut best_sign, mut best_n) = (0usize, 1.0f64, 0usize);
    for j in 0..nf {
        if inv_sd[j] == 0.0 {
            continue;
        }
        for sign in [1.0f64, -1.0] {
            let s: Vec<f64> = (0..n).map(|k| sign * z[k * nf + j]).collect();
            let q = tdc_q(d, rows, &s, &rank1);
            let acc = rank1
                .iter()
                .zip(&q)
                .filter(|(&k, &qv)| !d.is_decoy[rows[k]] && qv <= TRAIN_Q)
                .count();
            if acc > best_n {
                (best_j, best_sign, best_n) = (j, sign, acc);
            }
        }
    }
    let mut w = vec![0.0f64; nf];
    w[best_j] = best_sign;
    let mut b = 0.0f64;
    let score_all = |w: &[f64], b: f64| -> Vec<f64> {
        (0..n)
            .into_par_iter()
            .map(|k| {
                b + z[k * nf..(k + 1) * nf]
                    .iter()
                    .zip(w)
                    .map(|(a, c)| a * c)
                    .sum::<f64>()
            })
            .collect()
    };
    let mut score = score_all(&w, b);
    let mut trained = false;
    for _ in 0..N_ROUNDS {
        let win = winners(d, rows, &score);
        let q = tdc_q(d, rows, &score, &win);
        let mut sel = Vec::new();
        let mut y = Vec::new();
        for (&k, &qv) in win.iter().zip(&q) {
            if d.is_decoy[rows[k]] {
                sel.push(k);
                y.push(false);
            } else if qv <= TRAIN_Q {
                sel.push(k);
                y.push(true);
            }
        }
        let n_pos = y.iter().filter(|&&v| v).count();
        if n_pos < MIN_TRAIN_CLASS || y.len() - n_pos < MIN_TRAIN_CLASS {
            break;
        }
        let zs: Vec<f64> = sel
            .iter()
            .flat_map(|&k| z[k * nf..(k + 1) * nf].iter().copied())
            .collect();
        if !trained {
            // Start the first fit from zero, not from the single-feature
            // direction (whose scale is arbitrary).
            w.iter_mut().for_each(|v| *v = 0.0);
            b = 0.0;
        }
        fit_logistic(&zs, &y, nf, &mut w, &mut b);
        trained = true;
        score = score_all(&w, b);
    }

    // Common scale: the q = TRAIN_Q threshold of the training split maps to
    // 0 and its median decoy winner to -1.
    let win = winners(d, rows, &score);
    let q = tdc_q(d, rows, &score, &win);
    let threshold = win
        .iter()
        .zip(&q)
        .filter(|(&k, &qv)| !d.is_decoy[rows[k]] && qv <= TRAIN_Q)
        .map(|(&k, _)| score[k])
        .fold(f64::INFINITY, f64::min);
    let mut decoy_scores: Vec<f64> = win
        .iter()
        .filter(|&&k| d.is_decoy[rows[k]])
        .map(|&k| score[k])
        .collect();
    decoy_scores.sort_by(f64::total_cmp);
    let (shift, scale) = match decoy_scores.get(decoy_scores.len() / 2) {
        Some(&median) if threshold.is_finite() && threshold > median => {
            (threshold, 1.0 / (threshold - median))
        }
        _ => (0.0, 1.0),
    };
    FoldModel {
        model: LinearModel { w, b, mean, inv_sd },
        shift,
        scale,
        trained,
    }
}

/// Cross-validated scores of every row on the common scale, and the number of
/// folds that kept the initial single-feature score.
fn cv_scores(d: &PinData, seed: u64) -> (Vec<f64>, usize) {
    let n = d.is_decoy.len();
    let nf = d.n_features;
    let folds: Vec<usize> = d
        .spectrum_hashes
        .iter()
        .map(|&h| fold_of(h, seed))
        .collect();
    let mut out = vec![0.0f64; n];
    let mut untrained = 0;
    for fold in 0..N_FOLDS {
        let train: Vec<usize> = (0..n).filter(|&r| folds[r] != fold).collect();
        let test: Vec<usize> = (0..n).filter(|&r| folds[r] == fold).collect();
        if test.is_empty() {
            continue;
        }
        let fm = train_fold(d, &train);
        if !fm.trained {
            untrained += 1;
        }
        for r in test {
            out[r] = (fm.model.score(&d.x[r * nf..(r + 1) * nf]) - fm.shift) * fm.scale;
        }
    }
    (out, untrained)
}

/// PEP of each winner (in `idx` order, `score` descending not required):
/// decoys / targets at that score from an isotonic fit of the decoy share,
/// capped at 1. Winners with exactly equal scores start as one block, so they
/// share one PEP whatever their order in `idx`.
fn winner_pep(d: &PinData, idx: &[usize], score: &[f64]) -> Vec<f64> {
    let mut order: Vec<usize> = (0..idx.len()).collect();
    order.sort_by(|&a, &b| score[idx[b]].total_cmp(&score[idx[a]]));
    // Pool adjacent violators: the decoy share must not fall as the score
    // falls. Blocks of (sum, count) in descending-score order; a run of tied
    // scores enters as one block.
    let mut sums: Vec<f64> = Vec::new();
    let mut counts: Vec<f64> = Vec::new();
    let mut start = 0usize;
    while start < order.len() {
        let tied = score[idx[order[start]]];
        let end = start
            + order[start..]
                .iter()
                .take_while(|&&o| score[idx[o]].total_cmp(&tied).is_eq())
                .count();
        sums.push(
            order[start..end]
                .iter()
                .filter(|&&o| d.is_decoy[idx[o]])
                .count() as f64,
        );
        counts.push((end - start) as f64);
        start = end;
        while sums.len() > 1 {
            let k = sums.len();
            if sums[k - 2] / counts[k - 2] <= sums[k - 1] / counts[k - 1] {
                break;
            }
            sums[k - 2] += sums[k - 1];
            counts[k - 2] += counts[k - 1];
            sums.pop();
            counts.pop();
        }
    }
    let mut pep = vec![1.0f64; idx.len()];
    let mut pos = 0usize;
    for (s, c) in sums.iter().zip(&counts) {
        let share = s / c;
        let v = if share >= 0.5 {
            1.0
        } else {
            (share / (1.0 - share)).clamp(PROB_EPS, 1.0)
        };
        for _ in 0..(*c as usize) {
            pep[order[pos]] = v;
            pos += 1;
        }
    }
    pep
}

/// Native rescore a PIN → `SpecId → PercolatorPsm` of the winning targets
/// (monotone TDC q-value + PEP). Same map shape as `output::run_percolator`.
#[cfg(test)]
pub fn native_rescore_pin(
    pin_text: &str,
    seed: u64,
) -> Result<HashMap<String, PercolatorPsm>, String> {
    native_rescore_pin_all(pin_text, seed).map(|r| r.targets)
}

/// The native rescorer's result: the winning target and decoy PSMs of every
/// competition group, keyed by SpecId (the decoys feed FDR within each
/// class, protein inference and the match-between-runs library), plus run
/// diagnostics.
pub struct NativeRescore {
    pub targets: HashMap<String, PercolatorPsm>,
    pub decoys: HashMap<String, PercolatorPsm>,
    /// Rows that lost target-decoy competition in their group.
    pub competed_away: usize,
    /// Cross-validation folds whose training split had too few confident
    /// targets or decoys to train on; they ranked by the best single feature.
    pub untrained_folds: usize,
}

/// Native rescore of a PIN, keeping the winning decoys; see
/// [`native_rescore_pin`].
pub fn native_rescore_pin_all(pin_text: &str, seed: u64) -> Result<NativeRescore, String> {
    let d = parse_pin(pin_text)?;
    let n = d.is_decoy.len();
    if n == 0 {
        return Ok(NativeRescore {
            targets: HashMap::new(),
            decoys: HashMap::new(),
            competed_away: 0,
            untrained_folds: 0,
        });
    }
    let (scores, untrained_folds) = cv_scores(&d, seed);
    let all: Vec<usize> = (0..n).collect();
    let win = winners(&d, &all, &scores);
    let q = tdc_q(&d, &all, &scores, &win);
    let pep = winner_pep(&d, &win, &scores);
    let mut targets = HashMap::with_capacity(win.len());
    let mut decoys = HashMap::new();
    for (k, &i) in win.iter().enumerate() {
        let psm = PercolatorPsm {
            psm_id: d.spec_ids[i].clone(),
            q_value: q[k],
            pep: pep[k],
            score: Some(scores[i]),
            peptide: d.peptides[i].clone(),
            proteins: d.proteins[i].clone(),
        };
        // Target-only result contract (as `output::run_percolator`): a decoy
        // in the target map would reach the filtered TSV / QPX as an accepted
        // identification.
        let map = if d.is_decoy[i] {
            &mut decoys
        } else {
            &mut targets
        };
        if map.insert(d.spec_ids[i].clone(), psm).is_some() {
            // SpecIds are unique per emitted PSM (multi-row scans carry a
            // per-row suffix); a collision means a corrupt PIN.
            return Err(format!(
                "duplicate PSM identifier {:?} in PIN",
                d.spec_ids[i]
            ));
        }
    }
    Ok(NativeRescore {
        targets,
        decoys,
        competed_away: n - win.len(),
        untrained_folds,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a synthetic PIN with one separating feature `f1` (targets high,
    /// decoys low) plus a pure-noise feature `f2`.
    fn synth_pin(n: usize, separable: bool) -> String {
        let mut s = String::from("SpecId\tLabel\tScanNr\tf1\tf2\tPeptide\tProteins\n");
        for i in 0..n {
            let is_decoy = i % 2 == 1;
            let label = if is_decoy { -1 } else { 1 };
            // deterministic pseudo-noise in [0,1)
            let noise = (((i as u64).wrapping_mul(2654435761) >> 8) % 1000) as f32 / 1000.0;
            let f1 = if separable {
                if is_decoy {
                    noise
                } else {
                    1.0 + noise
                }
            } else {
                noise
            };
            let f2 = ((((i as u64).wrapping_mul(40503) >> 4) % 1000) as f32) / 1000.0;
            let sid = format!("scan={}_{}", i, label);
            s.push_str(&format!(
                "{sid}\t{label}\t{i}\t{f1:.4}\t{f2:.4}\tK.PEPTIDE{i}K.R\tP{i}\n"
            ));
        }
        s
    }

    #[test]
    fn separable_targets_get_low_q() {
        let pin = synth_pin(600, true);
        let map = native_rescore_pin(&pin, 42).unwrap();
        // many target PSMs should reach q <= 0.01
        let confident = map
            .values()
            .filter(|p| !p.psm_id.ends_with("_-1") && p.q_value <= 0.01)
            .count();
        assert!(
            confident > 50,
            "expected many confident targets, got {confident}"
        );
        // a clearly-high-scoring target must have a small PEP
        let min_pep = map
            .values()
            .filter(|p| !p.psm_id.ends_with("_-1"))
            .map(|p| p.pep)
            .fold(1.0f64, f64::min);
        assert!(
            min_pep < 0.5,
            "best target PEP should be small, got {min_pep}"
        );
    }

    #[test]
    fn rows_without_scan_numbers_are_folded_by_spectrum() {
        // Every ScanNr 0 (native ids without a scan number): the folds must
        // still split the spectra, or one fold holds everything and nothing
        // is scored.
        let pin = synth_pin(600, true)
            .lines()
            .enumerate()
            .map(|(i, l)| {
                if i == 0 {
                    return l.to_string();
                }
                let mut f: Vec<&str> = l.split('\t').collect();
                f[2] = "0";
                f.join("\t")
            })
            .collect::<Vec<_>>()
            .join("\n");
        let map = native_rescore_pin(&pin, 42).unwrap();
        let confident = map.values().filter(|p| p.q_value <= 0.01).count();
        assert!(confident > 50, "got {confident}");
    }

    #[test]
    fn spectrum_key_drops_scan_rank_and_row() {
        assert_eq!(spectrum_key("run/index=101_101_1", 101), "run/index=101");
        assert_eq!(spectrum_key("run/index=101_101_2_1", 101), "run/index=101");
        assert_eq!(
            spectrum_key("run/cycle=4 experiment=2_0_1", 0),
            "run/cycle=4 experiment=2"
        );
        // Rank or row equal to the scan number: every row of the spectrum
        // keeps the same key.
        for id in ["run/x_1_1", "run/x_1_1_1", "run/x_1_2_1", "run/x_1_1_0"] {
            assert_eq!(spectrum_key(id, 1), "run/x", "{id}");
        }
        assert_eq!(spectrum_key("run/t_5_3_5_1", 5), "run/t_5_3");
    }

    #[test]
    fn pure_noise_yields_few_confident() {
        // No real signal → the CV must NOT manufacture confident IDs (no leakage).
        let pin = synth_pin(600, false);
        let map = native_rescore_pin(&pin, 7).unwrap();
        let confident = map
            .values()
            .filter(|p| !p.psm_id.ends_with("_-1") && p.q_value <= 0.01)
            .count();
        assert!(
            confident < 30,
            "pure noise should yield ~no confident IDs (leakage!), got {confident}"
        );
    }

    #[test]
    fn decoys_never_appear_in_results() {
        // Target-only contract: no decoy row may surface in the result map, even
        // a high-scoring (low-q) one — otherwise it leaks into the filtered TSV.
        let pin = synth_pin(600, true);
        let map = native_rescore_pin(&pin, 42).unwrap();
        assert!(
            map.keys().all(|k| !k.ends_with("_-1")),
            "a decoy SpecId leaked into the target-only result map"
        );
        // And the synthetic PIN had n/2 decoys, so the map must be target-only.
        assert!(
            map.len() <= 300,
            "result map must exclude decoys, got {} entries",
            map.len()
        );
    }

    #[test]
    fn all_keeps_decoy_scores_and_the_same_targets() {
        let pin = synth_pin(600, true);
        let all = native_rescore_pin_all(&pin, 42).unwrap();
        assert_eq!(all.targets, native_rescore_pin(&pin, 42).unwrap());
        assert_eq!(all.decoys.len(), 300);
        assert!(all.decoys.keys().all(|k| k.ends_with("_-1")));
        assert!(all.targets.values().all(|p| p.score.is_some()));
        // Separable data: the targets outscore the decoys on average.
        let mean = |v: Vec<f64>| v.iter().sum::<f64>() / v.len() as f64;
        let t = mean(all.targets.values().map(|p| p.score.unwrap()).collect());
        let d = mean(all.decoys.values().map(|p| p.score.unwrap()).collect());
        assert!(t > d, "{t} <= {d}");
    }

    #[test]
    fn malformed_rows_fail_loud() {
        // Truncated row, bad label, bad scan, and non-numeric feature must all error.
        let base = "SpecId\tLabel\tScanNr\tf1\tPeptide\tProteins\n";
        assert!(
            native_rescore_pin(&format!("{base}s1\t1\t10\n"), 1).is_err(),
            "truncated"
        );
        assert!(
            native_rescore_pin(&format!("{base}s1\t2\t10\t0.5\tK.P.K\tP1\n"), 1).is_err(),
            "label=2"
        );
        assert!(
            native_rescore_pin(&format!("{base}s1\t1\tNaN\t0.5\tK.P.K\tP1\n"), 1).is_err(),
            "bad scan"
        );
        assert!(
            native_rescore_pin(&format!("{base}s1\t1\t10\tabc\tK.P.K\tP1\n"), 1).is_err(),
            "bad feature"
        );
    }

    /// Deterministic pseudo-noise in [0, 1).
    /// Rows of two files with the same ScanNr and ExpMass are different
    /// spectra: each file's rows compete only among themselves.
    #[test]
    fn spectra_of_different_files_never_compete() {
        let mut s = String::from("SpecId\tLabel\tScanNr\tExpMass\tf1\tPeptide\tProteins\n");
        for i in 0..200u64 {
            let f1 = if i % 2 == 0 {
                2.0 + noise(i, 1)
            } else {
                noise(i, 1)
            };
            s.push_str(&format!(
                "a/scan={i}_{i}_1\t1\t{i}\t1000.5\t{f1}\tK.AAA{i}K.R\tP{i}\n"
            ));
            s.push_str(&format!(
                "b/scan={i}_{i}_1\t-1\t{i}\t1000.5\t{}\tK.DDD{i}K.R\tXXX_P{i}\n",
                noise(i, 3)
            ));
        }
        let r = native_rescore_pin_all(&s, 42).unwrap();
        assert_eq!(
            r.competed_away, 0,
            "no row of one file may beat a row of another"
        );
        assert_eq!(r.targets.len(), 200);
        assert_eq!(r.decoys.len(), 200);
        assert_eq!(run_of("a/scan=1_1_1"), "a");
        assert_eq!(run_of("scan=1_1_1"), "");
    }

    fn noise(i: u64, salt: u64) -> f32 {
        ((i.wrapping_add(salt).wrapping_mul(2654435761) >> 7) % 1000) as f32 / 1000.0
    }

    #[test]
    fn one_psm_per_scan_and_precursor_mass() {
        // Every spectrum has a target and a decoy row at one ExpMass and a
        // second target at another ExpMass (another charge).
        let mut s = String::from("SpecId\tLabel\tScanNr\tExpMass\tf1\tf2\tPeptide\tProteins\n");
        for i in 0..400u64 {
            let good = i % 2 == 0;
            let t1 = if good { 2.0 + noise(i, 1) } else { noise(i, 1) };
            s.push_str(&format!(
                "s{i}_{i}_1_0\t1\t{i}\t1000.5\t{t1}\t{}\tK.AAA{i}K.R\tP{i}\n",
                noise(i, 2)
            ));
            s.push_str(&format!(
                "s{i}_{i}_2_1\t-1\t{i}\t1000.5\t{}\t{}\tK.DDD{i}K.R\tXXX_P{i}\n",
                noise(i, 3),
                noise(i, 4)
            ));
            s.push_str(&format!(
                "s{i}_{i}_1_2\t1\t{i}\t1500.75\t{}\t{}\tK.CCC{i}K.R\tQ{i}\n",
                noise(i, 5),
                noise(i, 6)
            ));
        }
        let r = native_rescore_pin_all(&s, 42).unwrap();
        // 800 groups, one winner each.
        assert_eq!(r.targets.len() + r.decoys.len(), 800);
        assert_eq!(r.competed_away, 400);
        for i in 0..400 {
            let a = format!("s{i}_{i}_1_0");
            let d = format!("s{i}_{i}_2_1");
            let n = r.targets.contains_key(&a) as u8 + r.decoys.contains_key(&d) as u8;
            assert_eq!(
                n, 1,
                "scan {i}: the 1000.5 group must report exactly one PSM"
            );
            assert!(r.targets.contains_key(&format!("s{i}_{i}_1_2")));
        }
        // PEP: within [0, 1], never above 1, and falls as the score rises.
        let mut v: Vec<(f64, f64)> = r
            .targets
            .values()
            .chain(r.decoys.values())
            .map(|p| (p.score.unwrap(), p.pep))
            .collect();
        v.sort_by(|a, b| a.0.total_cmp(&b.0));
        assert!(v.iter().all(|&(_, p)| (0.0..=1.0).contains(&p)));
        assert!(v.windows(2).all(|w| w[1].1 <= w[0].1 + 1e-12));
    }

    /// The root cause of the old rescorer's excess false discoveries: training
    /// every target as a positive teaches the model what separates FALSE
    /// targets from decoys. Here false targets score higher on a feature
    /// (`art`) than decoys and correct targets do. A rescorer trained
    /// on all targets ranks them above the decoys; the semi-supervised one
    /// learns only what the confident targets show.
    #[test]
    fn false_targets_unlike_decoys_are_not_rewarded() {
        let mut s = String::from("SpecId\tLabel\tScanNr\tsig\tart\tPeptide\tProteins\n");
        let n = 3000u64;
        for i in 0..n {
            // 1/3 correct targets, 1/3 false targets, 1/3 decoys.
            let (label, sig, art, prot) = match i % 3 {
                0 => (1, 1.8 + 2.0 * noise(i, 1), 0.4 * noise(i, 2), "TRUE"),
                1 => (1, 2.0 * noise(i, 1), 0.2 + 0.4 * noise(i, 2), "FALSE"),
                _ => (-1, 2.0 * noise(i, 1), 0.4 * noise(i, 2), "DECOY"),
            };
            s.push_str(&format!(
                "s{i}_{i}_1\t{label}\t{i}\t{sig}\t{art}\tK.P{i}K.R\t{prot}{i}\n"
            ));
        }
        let map = native_rescore_pin(&s, 42).unwrap();
        for q in [0.01, 0.05] {
            let acc: Vec<&PercolatorPsm> = map.values().filter(|p| p.q_value <= q).collect();
            let false_hits = acc
                .iter()
                .filter(|p| p.proteins.starts_with("FALSE"))
                .count();
            assert!(acc.len() > 300, "q {q}: {} accepted", acc.len());
            let fdp = false_hits as f64 / acc.len() as f64;
            assert!(
                fdp <= 2.0 * q + 0.01,
                "q {q}: FDP {fdp:.3} ({false_hits}/{})",
                acc.len()
            );
        }
    }

    #[test]
    fn rescoring_is_deterministic_across_thread_counts() {
        let pin = synth_pin(20_000, true);
        let one = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .unwrap();
        let four = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .unwrap();
        let a = one.install(|| native_rescore_pin(&pin, 3).unwrap());
        let b = four.install(|| native_rescore_pin(&pin, 3).unwrap());
        assert_eq!(a, b);
    }

    #[test]
    fn mass_and_file_name_columns_are_not_features() {
        let pin = "SpecId\tLabel\tScanNr\tExpMass\tCalcMass\tFileName\tf1\tPeptide\tProteins\n\
                   a/s_1_1\t1\t1\t500.1\t500.1\ta\t0.5\tK.P.K\tP1\n";
        assert_eq!(parse_pin(pin).unwrap().n_features, 1);
        let only_masses =
            "SpecId\tLabel\tScanNr\tExpMass\tCalcMass\tFileName\tPeptide\tProteins\n";
        assert!(parse_pin(only_masses).is_err());
    }

    #[test]
    fn logistic_fit_separates_and_solver_handles_spd() {
        // x > 0 → positive.
        let z: Vec<f64> = (0..200).map(|i| (i as f64 - 99.5) / 50.0).collect();
        let y: Vec<bool> = z.iter().map(|&v| v > 0.0).collect();
        let (mut w, mut b) = (vec![0.0], 0.0);
        fit_logistic(&z, &y, 1, &mut w, &mut b);
        assert!(w[0] > 0.5, "{w:?}");
        let x = cholesky_solve(&[4.0, 2.0, 2.0, 3.0], &[2.0, 1.0], 2).unwrap();
        assert!((x[0] - 0.5).abs() < 1e-12 && x[1].abs() < 1e-12, "{x:?}");
        assert!(cholesky_solve(&[0.0, 0.0, 0.0, 1.0], &[1.0, 1.0], 2).is_none());
    }

    fn pin_data(is_decoy: &[bool]) -> PinData {
        let n = is_decoy.len();
        PinData {
            spec_ids: (0..n).map(|i| format!("r{i}")).collect(),
            peptides: vec![String::new(); n],
            proteins: vec![String::new(); n],
            is_decoy: is_decoy.to_vec(),
            spectrum_hashes: vec![0; n],
            groups: (0..n as u32).collect(),
            first_in_group: vec![true; n],
            x: Vec::new(),
            n_features: 0,
        }
    }

    /// A target and a decoy with exactly equal scores get one PEP, the same
    /// whichever comes first.
    #[test]
    fn tied_scores_share_one_pep_in_any_order() {
        let score = [3.0, 2.0, 2.0, 1.0];
        let idx = [0, 1, 2, 3];
        // Order A: the tied target comes first; order B: the tied decoy does.
        let a = winner_pep(&pin_data(&[false, false, true, false]), &idx, &score);
        let b = winner_pep(&pin_data(&[false, true, false, false]), &idx, &score);
        assert_eq!(a[1], a[2], "tied rows share a PEP");
        assert_eq!(b[1], b[2]);
        assert_eq!(a, b, "the PEPs do not depend on the order of tied rows");
        assert!(a[0] <= a[1] && a[1] <= a[3]);
    }

    #[test]
    fn empty_and_headeronly_are_safe() {
        assert!(native_rescore_pin("", 1).is_err());
        let only_header = "SpecId\tLabel\tScanNr\tf1\tPeptide\tProteins\n";
        assert!(native_rescore_pin(only_header, 1).unwrap().is_empty());
    }
}
