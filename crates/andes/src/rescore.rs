#![allow(clippy::needless_range_loop)]
//! Native (Percolator-free) PSM rescorer: a cross-validated GBDT over the PIN
//! feature matrix → per-PSM q-value + PEP.
//!
//! ROLE: the self-contained FALLBACK for when a Percolator backend is not
//! available (benchmarking / offline / zero-dependency). It is NOT a
//! production-grade FDR engine — Percolator (`--rescore`) remains the
//! recommended path and the reference cross-check; this is honestly surfaced as
//! a *native* rescore.
//!
//! METHOD (mokapot-style, leakage-safe): 3-fold target-decoy cross-validation.
//! Rows are folded by spectrum (ScanNr) so a spectrum's target+decoy PSMs never
//! split across folds. For each fold a GBDT classifier (target=1 / decoy=0) is
//! trained on the OTHER folds and scores the held-out fold — so every PSM is
//! scored by a model that never trained on it (no decoy leakage → no fake-good
//! q-values). The pooled CV scores feed a monotone TDC q-value and a calibrated
//! PEP. The output map matches `output::run_percolator`'s shape so the
//! downstream QPX-injection + filtered-TSV path is reused unchanged.

use std::collections::HashMap;

use model_train::gbdt::train::{train_gbdt, Dataset, TrainParams};
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
    // Require an explicit PSM identifier column — never silently fall back to
    // column 0, which would let a malformed header join on an arbitrary field.
    let id_i = find("SpecId")
        .or_else(|| find("PSMId"))
        .ok_or("PIN missing SpecId/PSMId column")?;
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
        // ScanNr drives the CV fold assignment — a bad value would silently
        // mis-fold the spectrum, so fail loud instead of defaulting to 0.
        let scan_raw = f[scan_i].trim();
        let scan: u32 = scan_raw.parse().map_err(|_| {
            format!(
                "PIN line {} ScanNr is {scan_raw:?} (expected integer)",
                lineno + 2
            )
        })?;
        d.scans.push(scan);
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

/// Deterministic fold assignment by spectrum (splitmix64 of ScanNr).
fn fold_of(scan: u32, seed: u64) -> usize {
    let mut z = (seed ^ (scan as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15))
        .wrapping_add(0x9e37_79b9_7f4a_7c15);
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^= z >> 31;
    (z % N_FOLDS as u64) as usize
}

/// 3-fold CV GBDT scores (calibrated P(target)) per row — leakage-safe.
fn cv_scores(d: &PinData, seed: u64) -> Vec<f32> {
    let n = d.is_decoy.len();
    let nf = d.n_features;
    let folds: Vec<usize> = d.scans.iter().map(|&s| fold_of(s, seed)).collect();
    let mut out = vec![0.0f32; n];
    let params = TrainParams::default();
    for fold in 0..N_FOLDS {
        let mut x = Vec::new();
        let mut y = Vec::new();
        let mut groups = Vec::new();
        for r in 0..n {
            if folds[r] == fold {
                continue;
            }
            x.extend_from_slice(&d.x[r * nf..(r + 1) * nf]);
            y.push(if d.is_decoy[r] { 0u8 } else { 1u8 });
            groups.push(d.scans[r]);
        }
        // Degenerate fold (all one class) → leave held-out scores at 0.
        if y.is_empty() || y.iter().all(|&v| v == 0) || y.iter().all(|&v| v == 1) {
            continue;
        }
        let ds = Dataset {
            x,
            y,
            groups,
            n_features: nf,
        };
        // A fold that fails the trainer's quality gate (too few rows, no
        // held-out signal, etc.) contributes no held-out scores — same as the
        // degenerate-fold skip above. The native rescorer is an opt-in
        // benchmarking fallback, so a non-deployable fold is left at 0 rather
        // than aborting the whole rescore.
        let model = match train_gbdt(&ds, &params, seed.wrapping_add(fold as u64 + 1)) {
            Ok(m) => m,
            Err(_) => continue,
        };
        for r in 0..n {
            if folds[r] != fold {
                continue;
            }
            out[r] = model.predict_proba(&d.x[r * nf..(r + 1) * nf]);
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
        .map(|i| ScoredLabel {
            score: scores[i],
            is_decoy: d.is_decoy[i],
        })
        .collect();
    let q = qvalues(&items);
    let mut map = HashMap::with_capacity(n);
    for i in 0..n {
        // TARGET-ONLY result contract: match `output::run_percolator`, whose
        // results file contains only target PSMs. Decoys above drive the TDC
        // q-value, but emitting one here would let a low-q decoy leak into the
        // filtered TSV / QPX as an accepted identification (the downstream path
        // filters by q-value/PEP and has no decoy label).
        if d.is_decoy[i] {
            continue;
        }
        // PEP = 1 − calibrated P(target).
        let pep = (1.0 - scores[i] as f64).clamp(PROB_EPS, 1.0);
        if map
            .insert(
                d.spec_ids[i].clone(),
                PercolatorPsm {
                    psm_id: d.spec_ids[i].clone(),
                    q_value: q[i],
                    pep,
                    peptide: d.peptides[i].clone(),
                    proteins: d.proteins[i].clone(),
                },
            )
            .is_some()
        {
            // SpecIds are unique per emitted PSM (multi-row scans carry a per-row
            // suffix); a collision means a corrupt PIN whose last-wins overwrite
            // would silently drop a target — fail loud instead.
            return Err(format!(
                "duplicate target PSM identifier {:?} in PIN",
                d.spec_ids[i]
            ));
        }
    }
    Ok(map)
}

/// Native-rescore a PIN and return per-row (spec_id, is_decoy, q_value, score)
/// for EVERY row (targets AND decoys), unlike `native_rescore_pin` which is
/// target-only. Used to seed the symmetric cross-spectrum decoy-transfer graph.
pub fn native_rescore_qvalues(
    pin_text: &str,
    seed: u64,
) -> Result<Vec<(String, bool, f64, f64)>, String> {
    let d = parse_pin(pin_text)?;
    let n = d.is_decoy.len();
    if n == 0 {
        return Ok(Vec::new());
    }
    let scores = cv_scores(&d, seed);
    let items: Vec<ScoredLabel> = (0..n)
        .map(|i| ScoredLabel {
            score: scores[i],
            is_decoy: d.is_decoy[i],
        })
        .collect();
    let q = qvalues(&items);
    Ok((0..n)
        .map(|i| (d.spec_ids[i].clone(), d.is_decoy[i], q[i], scores[i] as f64))
        .collect())
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

    #[test]
    fn empty_and_headeronly_are_safe() {
        assert!(native_rescore_pin("", 1).is_err());
        let only_header = "SpecId\tLabel\tScanNr\tf1\tPeptide\tProteins\n";
        assert!(native_rescore_pin(only_header, 1).unwrap().is_empty());
    }

    #[test]
    fn qvalues_includes_targets_and_decoys_with_finite_q_in_range() {
        let pin = synth_pin(600, true);
        let rows = native_rescore_qvalues(&pin, 42).unwrap();
        // Unlike native_rescore_pin (target-only), every PIN row must come back.
        assert_eq!(
            rows.len(),
            600,
            "expected all rows (targets+decoys), got {}",
            rows.len()
        );
        let n_decoys = rows.iter().filter(|(_, is_decoy, _, _)| *is_decoy).count();
        let n_targets = rows.len() - n_decoys;
        assert_eq!(
            n_decoys, 300,
            "expected the synthetic PIN's 300 decoy rows preserved"
        );
        assert_eq!(
            n_targets, 300,
            "expected the synthetic PIN's 300 target rows preserved"
        );
        for (spec_id, _is_decoy, q, score) in &rows {
            assert!(
                q.is_finite() && *q >= 0.0 && *q <= 1.0,
                "q-value out of range for {spec_id}: {q}"
            );
            assert!(score.is_finite(), "score not finite for {spec_id}: {score}");
        }
        // Separable synthetic data: many confident targets should reach q<=0.01,
        // same signal `native_rescore_pin`'s test checks, but now decoys are visible too.
        let confident_targets = rows
            .iter()
            .filter(|(_, is_decoy, q, _)| !is_decoy && *q <= 0.01)
            .count();
        assert!(
            confident_targets > 50,
            "expected many confident targets, got {confident_targets}"
        );
    }

    #[test]
    fn qvalues_empty_pin_is_safe() {
        assert!(native_rescore_qvalues("", 1).is_err());
        let only_header = "SpecId\tLabel\tScanNr\tf1\tPeptide\tProteins\n";
        assert!(native_rescore_qvalues(only_header, 1).unwrap().is_empty());
    }
}
