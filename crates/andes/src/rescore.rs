#![allow(clippy::needless_range_loop)]
//! Native (Percolator-free) PSM rescorer: a Percolator-style semi-supervised
//! **regularized linear** model over the PIN feature matrix → per-PSM q-value + PEP.
//!
//! ROLE: a self-contained, no-external-dependency rescorer. Percolator (`--rescore`)
//! remains the reference; this aims to match it standalone.
//!
//! METHOD (leakage-safe, mokapot/Percolator-style): 3-fold target-decoy CV folded
//! by spectrum (ScanNr) so a spectrum's target+decoy PSMs never split across folds.
//! Per held-out fold: z-score the features on the train rows, bootstrap a ranking
//! from the single best feature, then run the semi-supervised loop — {select
//! confident targets (q < 1%) + ALL decoys → fit an **L2-regularized logistic
//! regression** (class-balanced, warm-started) → rescore the train rows} — and the
//! final weights score the held-out fold. The L2 regularization is what lets the
//! confident-set iteration refine WITHOUT the self-confirmation overfit a GBDT
//! suffers (a tree memorizes the trivially-separable confident set; a regularized
//! linear hyperplane generalizes). This is why Percolator/Sage/mokapot use linear
//! models for rescoring. The test fold never trains, so downstream TDC q-values
//! stay honest. Output map matches `output::run_percolator` so the QPX-injection +
//! filtered-TSV path is reused unchanged.

use std::collections::HashMap;

use output::PercolatorPsm;
use search::tdc::{qvalues, ScoredLabel};

const N_FOLDS: usize = 3;
const PROB_EPS: f64 = 1e-6;

/// Numeric feature matrix + labels parsed from a Percolator PIN.
struct PinData {
    spec_ids: Vec<String>,
    peptides: Vec<String>,
    proteins: Vec<String>,
    is_decoy: Vec<bool>,
    scans: Vec<u32>,
    /// row-major, `n_rows × n_features`
    x: Vec<f32>,
    n_features: usize,
}

/// Parse a Percolator PIN. Layout:
/// `SpecId \t Label \t ScanNr \t <numeric features…> \t Peptide \t Proteins…`.
/// Feature columns = everything strictly between `ScanNr` and `Peptide`.
fn parse_pin(text: &str) -> Result<PinData, String> {
    let mut lines = text.lines();
    let header = lines.next().ok_or("empty PIN")?;
    let cols: Vec<&str> = header.split('\t').collect();
    let find = |name: &str| cols.iter().position(|c| c.eq_ignore_ascii_case(name));
    let id_i = find("SpecId").or_else(|| find("PSMId")).unwrap_or(0);
    let label_i = find("Label").ok_or("PIN missing Label column")?;
    let scan_i = find("ScanNr").ok_or("PIN missing ScanNr column")?;
    let pep_i = find("Peptide").ok_or("PIN missing Peptide column")?;
    if pep_i <= scan_i + 1 {
        return Err("PIN has no feature columns between ScanNr and Peptide".into());
    }
    let feat_cols: Vec<usize> = ((scan_i + 1)..pep_i).collect();
    let n_features = feat_cols.len();

    let mut d = PinData {
        spec_ids: Vec::new(),
        peptides: Vec::new(),
        proteins: Vec::new(),
        is_decoy: Vec::new(),
        scans: Vec::new(),
        x: Vec::new(),
        n_features,
    };
    for line in lines {
        if line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() <= pep_i {
            continue;
        }
        let label: i32 = f[label_i].trim().parse().unwrap_or(1);
        d.is_decoy.push(label < 0);
        d.spec_ids.push(f[id_i].to_string());
        d.scans
            .push(f.get(scan_i).and_then(|s| s.trim().parse().ok()).unwrap_or(0));
        d.peptides.push(f[pep_i].to_string());
        d.proteins
            .push(f.get(pep_i + 1..).map(|s| s.join(";")).unwrap_or_default());
        for &ci in &feat_cols {
            let v: f32 = f.get(ci).and_then(|s| s.trim().parse().ok()).unwrap_or(0.0);
            d.x.push(if v.is_finite() { v } else { 0.0 });
        }
    }
    Ok(d)
}

/// Deterministic fold assignment by spectrum (splitmix64 of ScanNr).
fn fold_of(scan: u32, seed: u64) -> usize {
    let mut z = (seed ^ (scan as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15))
        .wrapping_add(0x9e37_79b9_7f4a_7c15);
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^= z >> 31;
    (z % N_FOLDS as u64) as usize
}

/// Semi-supervised refinement iterations per fold (Percolator-style; converges fast).
const N_ITERS: usize = 10;
/// q-value threshold for selecting confident targets as positives.
const SEL_FDR: f64 = 0.01;
/// L2 regularization strength (on z-scored features) — the overfit guard.
const LAMBDA: f32 = 1e-3;
/// Gradient-descent learning rate + iterations (z-scored features → well-conditioned).
const LR: f32 = 0.5;
const GD_ITERS: usize = 100;

/// Build the z-scored feature matrix for ALL rows, using per-feature mean/std
/// computed over `train` (std floored). Precomputing once per fold turns every
/// later dot-product into a plain `w·z` (no re-standardizing in the GD loop).
fn zscore_matrix(d: &PinData, train: &[usize]) -> Vec<f32> {
    let nf = d.n_features;
    let n = d.is_decoy.len();
    let mut mean = vec![0.0f32; nf];
    for &r in train {
        for j in 0..nf {
            mean[j] += d.x[r * nf + j];
        }
    }
    let inv = 1.0 / train.len().max(1) as f32;
    for j in 0..nf {
        mean[j] *= inv;
    }
    let mut m2 = vec![0.0f32; nf];
    for &r in train {
        for j in 0..nf {
            let dlt = d.x[r * nf + j] - mean[j];
            m2[j] += dlt * dlt;
        }
    }
    let std: Vec<f32> = (0..nf)
        .map(|j| {
            let v = (m2[j] * inv).sqrt();
            if v > 1e-6 { v } else { 1.0 }
        })
        .collect();
    let mut zx = vec![0.0f32; n * nf];
    for r in 0..n {
        for j in 0..nf {
            zx[r * nf + j] = (d.x[r * nf + j] - mean[j]) / std[j];
        }
    }
    zx
}

/// Linear score `w·zx_r + bias` (`w` has length nf+1; bias is last).
#[inline]
fn dot(zx: &[f32], nf: usize, r: usize, w: &[f32]) -> f32 {
    let base = r * nf;
    let mut z = w[nf];
    for j in 0..nf {
        z += w[j] * zx[base + j];
    }
    z
}

/// L2-regularized logistic regression by gradient descent over `rows` (already
/// z-scored in `zx`). `y`: +1 target / -1 decoy. `wpos`/`wneg` are class weights.
/// `w` is warm-started (len nf+1; last = bias). Bias is not regularized.
fn train_logistic(
    zx: &[f32],
    nf: usize,
    rows: &[usize],
    y: &[i8],
    mut w: Vec<f32>,
    wpos: f32,
    wneg: f32,
) -> Vec<f32> {
    if w.len() != nf + 1 {
        w = vec![0.0f32; nf + 1];
    }
    let inv_n = 1.0 / rows.len().max(1) as f32;
    let mut grad = vec![0.0f32; nf + 1];
    for _ in 0..GD_ITERS {
        grad.iter_mut().for_each(|g| *g = 0.0);
        for (i, &r) in rows.iter().enumerate() {
            let base = r * nf;
            let z = dot(zx, nf, r, &w);
            let yi = y[i] as f32;
            let cw = if y[i] > 0 { wpos } else { wneg };
            // d/dz log(1+exp(-yi·z)) = -yi · sigmoid(-yi·z)
            let s = 1.0 / (1.0 + (yi * z).exp());
            let g = -yi * s * cw;
            for j in 0..nf {
                grad[j] += g * zx[base + j];
            }
            grad[nf] += g;
        }
        for j in 0..nf {
            w[j] -= LR * (grad[j] * inv_n + LAMBDA * w[j]);
        }
        w[nf] -= LR * grad[nf] * inv_n;
    }
    w
}

/// Initial per-train-row ranking: the single (feature, sign) that yields the most
/// train targets at q < SEL_FDR. Bootstraps the semi-supervised loop from a
/// known-discriminative direction.
fn init_scores(d: &PinData, train: &[usize]) -> Vec<f32> {
    let nf = d.n_features;
    let mut best = vec![0.0f32; train.len()];
    let mut best_count: i64 = -1;
    for f in 0..nf {
        for &dir in &[1.0f32, -1.0f32] {
            let items: Vec<ScoredLabel> = train
                .iter()
                .map(|&r| ScoredLabel { score: dir * d.x[r * nf + f], is_decoy: d.is_decoy[r] })
                .collect();
            let q = qvalues(&items);
            let cnt = (0..train.len())
                .filter(|&k| !d.is_decoy[train[k]] && q[k] < SEL_FDR)
                .count() as i64;
            if cnt > best_count {
                best_count = cnt;
                best = items.iter().map(|it| it.score).collect();
            }
        }
    }
    best
}

/// Semi-supervised, leakage-safe CV scores via the regularized LINEAR model.
fn cv_scores(d: &PinData, seed: u64) -> Vec<f32> {
    let n = d.is_decoy.len();
    let nf = d.n_features;
    let folds: Vec<usize> = d.scans.iter().map(|&s| fold_of(s, seed)).collect();
    let mut out = vec![0.0f32; n];
    for fold in 0..N_FOLDS {
        let train: Vec<usize> = (0..n).filter(|&r| folds[r] != fold).collect();
        let test: Vec<usize> = (0..n).filter(|&r| folds[r] == fold).collect();
        let n_dec = train.iter().filter(|&&r| d.is_decoy[r]).count();
        if train.is_empty() || n_dec == 0 || n_dec == train.len() {
            continue;
        }
        let zx = zscore_matrix(d, &train);
        let mut train_scores = init_scores(d, &train);
        let mut w: Vec<f32> = vec![0.0f32; nf + 1];
        let mut have_model = false;
        for _ in 0..N_ITERS {
            let items: Vec<ScoredLabel> = (0..train.len())
                .map(|k| ScoredLabel { score: train_scores[k], is_decoy: d.is_decoy[train[k]] })
                .collect();
            let q = qvalues(&items);
            let (mut sub, mut y) = (Vec::new(), Vec::new());
            for k in 0..train.len() {
                let r = train[k];
                // positives = confident targets; negatives = ALL decoys.
                if !(d.is_decoy[r] || q[k] < SEL_FDR) {
                    continue;
                }
                sub.push(r);
                y.push(if d.is_decoy[r] { -1i8 } else { 1i8 });
            }
            let npos = y.iter().filter(|&&v| v > 0).count();
            let nneg = y.len() - npos;
            if npos == 0 || nneg == 0 {
                break;
            }
            // Class-balanced weights (Cpos/Cneg ≈ inverse class frequency).
            let wpos = y.len() as f32 / (2.0 * npos as f32);
            let wneg = y.len() as f32 / (2.0 * nneg as f32);
            w = train_logistic(&zx, nf, &sub, &y, w, wpos, wneg);
            train_scores = train.iter().map(|&r| dot(&zx, nf, r, &w)).collect();
            have_model = true;
        }
        if have_model {
            for &r in &test {
                out[r] = dot(&zx, nf, r, &w);
            }
        } else {
            // Degenerate fold → fall back to the bootstrap ranking.
            let ti = init_scores(d, &test);
            for (k, &r) in test.iter().enumerate() {
                out[r] = ti[k];
            }
        }
    }
    out
}

/// Native rescore a PIN → `SpecId → PercolatorPsm` (monotone TDC q-value +
/// calibrated PEP). Same map shape as `output::run_percolator`, so the caller
/// reuses the QPX-injection + filtered-TSV path unchanged.
pub fn native_rescore_pin(
    pin_text: &str,
    seed: u64,
) -> Result<HashMap<String, PercolatorPsm>, String> {
    let d = parse_pin(pin_text)?;
    let n = d.is_decoy.len();
    if n == 0 {
        return Ok(HashMap::new());
    }
    let scores = cv_scores(&d, seed);
    let items: Vec<ScoredLabel> = (0..n)
        .map(|i| ScoredLabel { score: scores[i], is_decoy: d.is_decoy[i] })
        .collect();
    let q = qvalues(&items);
    let mut map = HashMap::with_capacity(n);
    for i in 0..n {
        // PEP = P(decoy-like) = sigmoid(−score); high-scoring targets → ~0, decoys → ~1.
        let pep = (1.0 / (1.0 + (scores[i] as f64).exp())).clamp(PROB_EPS, 1.0);
        map.insert(
            d.spec_ids[i].clone(),
            PercolatorPsm {
                psm_id: d.spec_ids[i].clone(),
                q_value: q[i],
                pep,
                peptide: d.peptides[i].clone(),
                proteins: d.proteins[i].clone(),
            },
        );
    }
    Ok(map)
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
                if is_decoy { noise } else { 1.0 + noise }
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
        assert!(confident > 50, "expected many confident targets, got {confident}");
        // a clearly-high-scoring target must have a small PEP
        let min_pep = map
            .values()
            .filter(|p| !p.psm_id.ends_with("_-1"))
            .map(|p| p.pep)
            .fold(1.0f64, f64::min);
        assert!(min_pep < 0.5, "best target PEP should be small, got {min_pep}");
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
    fn empty_and_headeronly_are_safe() {
        assert!(native_rescore_pin("", 1).is_err());
        let only_header = "SpecId\tLabel\tScanNr\tf1\tPeptide\tProteins\n";
        assert!(native_rescore_pin(only_header, 1).unwrap().is_empty());
    }
}
