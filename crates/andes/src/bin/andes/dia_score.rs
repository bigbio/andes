//! DIA precursor rows for rescoring: one row per nominated precursor with its pass-2 evidence,
//! the best pass-1 PIN row of the same precursor, and a cross-fitted learned score.
//!
//! The learned score is a gradient-boosted classifier trained out of fold: folds are keyed by
//! protein, so a decoy always falls in the fold of its source target whatever the decoy strategy,
//! and no row is scored by a model that saw it.

use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use model_train::gbdt::train::{train_gbdt_from, Dataset, TrainParams, MIN_VAL_AUC};
use output::dia_pass2_table::Pass2Table;
use output::PercolatorPsm;
use rayon::prelude::*;

pub(crate) const GB_COLUMN: &str = "dia_gb_score";
/// Column present in every precursor PIN (1 when the precursor has a pass-1 row).
pub(crate) const PRESENT_COLUMN: &str = "p1_present";
/// Precursor q-value of the reported counts and of the pass-2 size check.
const COMPARE_Q: f64 = 0.01;
/// Precursor q-values at which the precursor PIN must keep pass 1's target count.
const SELECT_Q: [f64; 2] = [0.005, 0.01];
const FOLDS: u64 = 3;
const ROUND_Q: [f64; 3] = [f64::INFINITY, 0.05, 0.01];
/// Trees a round after the first may add to its fold's model of the previous round: rounds 1
/// and 2 continue that model on their sharper positive set instead of growing up to 600 trees
/// anew (3-raw study, 5 seeds: N at 1 % entrapment FDP within seed spread, finalize -25 %).
const WARM_ROUNDS: u32 = 200;
const MISSING_DIA_LLR: f32 = -99.0;
const PASS1_SKIP: [&str; 6] = ["SpecId", "Label", "ScanNr", "Peptide", "Proteins", "FileName"];

fn params() -> TrainParams {
    TrainParams {
        n_rounds: 600,
        learning_rate: 0.08,
        max_depth: 8,
        n_bins: 256,
        lambda: 1.0,
        min_hessian: 5.0,
        neg_pos_ratio: f32::MAX,
        val_fraction: 0.1,
        early_stop_rounds: 50,
        pairwise: false,
        allow_degenerate: true,
    }
}

/// One-sided 95 % normal quantile: how many null standard errors above chance an AUC margin of
/// `MIN_VAL_AUC - 0.5` must be for a fold's held-out set to tell a model from chance.
const GATE_Z: f64 = 1.645;

/// Minimum positives in the training rows of a fold. The GBDT holds out `val_fraction` of its rows
/// to stop early. Without signal, the AUC of n positives against more negatives has standard error
/// sqrt(1 / (12 n)) (the Mann-Whitney null variance); the held-out set tells a model that beats
/// chance by the GBDT's AUC margin (`MIN_VAL_AUC - 0.5`) at one-sided 5 % only with
/// n >= (GATE_Z / (MIN_VAL_AUC - 0.5))^2 / 12 held-out positives, or n / `val_fraction` training
/// positives (910 with the parameters here). The AUC gate itself does not make a fold degenerate:
/// round 0 trains on every target, most of them false, and its held-out AUC sits near 0.52 even
/// when the model ranks well; the GBDT then keeps its trees uncalibrated, and the score uses the
/// raw margin.
fn min_fold_positives(p: &TrainParams) -> usize {
    let held_out = ((GATE_Z / (MIN_VAL_AUC - 0.5)).powi(2) / 12.0).ceil();
    (held_out / f64::from(p.val_fraction)).ceil() as usize
}

/// Accessions of a `Proteins` field (`;`- or tab-separated).
fn accessions(proteins: &str) -> impl Iterator<Item = &str> {
    proteins
        .split([';', '\t'])
        .map(str::trim)
        .filter(|a| !a.is_empty())
}

/// Fold key of a precursor row: a hash of the smallest of its protein accessions with
/// `decoy_prefix` taken off. A decoy gets the key of the target protein it was made from, under
/// every decoy strategy (reversed, shuffled or mutated), and the rows of one peptide get one key
/// even when it is in several proteins. A row without proteins falls back to its sorted residue
/// letters, equal for a target and every permutation of it.
pub(crate) fn fold_key(seq: &str, proteins: &str, decoy_prefix: &str) -> u64 {
    let smallest = accessions(proteins)
        .map(|a| a.strip_prefix(decoy_prefix).unwrap_or(a))
        .min();
    if let Some(a) = smallest {
        return model::hash::fnv1a64_extend(model::hash::FNV1A64_OFFSET, a.as_bytes().iter());
    }
    let mut b: Vec<u8> = seq.bytes().filter(u8::is_ascii_uppercase).collect();
    b.sort_unstable();
    model::hash::fnv1a64_extend(model::hash::FNV1A64_OFFSET, b.iter())
}

/// The decoy prefix of the table's accessions, read from the data: the leading part of the
/// smallest decoy accession whose removal turns the most decoy accessions into target accessions
/// of the table (the longest such part at ties). Empty when the table has no decoys.
pub(crate) fn decoy_prefix(t: &Pass2Table) -> String {
    let of = |decoy: bool| -> HashSet<&str> {
        (0..t.len())
            .filter(|&i| t.decoy[i] == decoy)
            .flat_map(|i| accessions(&t.proteins[i]))
            .collect()
    };
    let (targets, decoys) = (of(false), of(true));
    let Some(first) = decoys.iter().min() else {
        return String::new();
    };
    let hits = |k: usize| {
        let p = &first[..k];
        decoys
            .iter()
            .filter(|d| d.strip_prefix(p).is_some_and(|rest| targets.contains(rest)))
            .count()
    };
    let k = (0..=first.len())
        .filter(|&k| first.is_char_boundary(k))
        .max_by_key(|&k| (hits(k), k))
        .unwrap_or(0);
    first[..k].to_string()
}

/// Target-decoy q-values of a score ((D + 1) / T, decoys first at ties, monotone).
pub(crate) fn tdc_q(score: &[f32], target: &[bool]) -> Vec<f64> {
    let mut o: Vec<usize> = (0..score.len()).collect();
    o.sort_by(|&a, &b| {
        score[b]
            .total_cmp(&score[a])
            .then(target[a].cmp(&target[b]))
    });
    let (mut t, mut d) = (0usize, 0usize);
    let mut q = vec![1.0; score.len()];
    for &i in &o {
        if target[i] {
            t += 1;
        } else {
            d += 1;
        }
        q[i] = ((d + 1) as f64 / t.max(1) as f64).min(1.0);
    }
    let mut run = 1.0f64;
    for &i in o.iter().rev() {
        run = run.min(q[i]);
        q[i] = run;
    }
    q
}

/// The out-of-fold learned score and how its training went.
pub(crate) struct LearnedScore {
    /// The score of every row, from the last round whose folds all trained; `None` when a fold of
    /// the first round did not.
    pub(crate) score: Option<Vec<f32>>,
    /// Folds that did not train (too few positives, or a GBDT that returned no trees and so a
    /// constant score), over the rounds run.
    pub(crate) degenerate_folds: usize,
}

/// Out-of-fold learned score of every row: the raw boosted margin (log-odds), the same scale in
/// every fold. `x` is row-major with `names.len()` columns. A round with a degenerate fold ends
/// the training: the previous round's score is kept, or none after the first round.
pub(crate) fn learned_score(
    x: &[f32],
    names: &[String],
    target: &[bool],
    key: &[u64],
    seed: u64,
) -> LearnedScore {
    let n = target.len();
    let nf = names.len();
    let p = params();
    let min_pos = min_fold_positives(&p);
    let fold: Vec<u64> = key
        .iter()
        .map(|k| (k ^ seed).wrapping_mul(0x9e37_79b9_7f4a_7c15) % FOLDS)
        .collect();
    let mut positive: Vec<bool> = target.to_vec();
    let mut kept: Option<Vec<f32>> = None;
    let mut degenerate_folds = 0usize;
    let mut prev: Vec<Vec<scoring_crate::gbdt_eval::Tree>> = vec![Vec::new(); FOLDS as usize];
    for (round, &round_q) in ROUND_Q.iter().enumerate() {
        let t_round = std::time::Instant::now();
        type Fit = (Vec<(usize, f32)>, Vec<scoring_crate::gbdt_eval::Tree>);
        let fits: Vec<Result<Fit, String>> = (0..FOLDS)
            .into_par_iter()
            .map(|f| {
                let train: Vec<usize> = (0..n)
                    .filter(|&i| fold[i] != f && (!target[i] || positive[i]))
                    .collect();
                let n_pos = train.iter().filter(|&&i| target[i]).count();
                if n_pos < min_pos {
                    return Err(format!("fold {f}: {n_pos} positives < {min_pos}"));
                }
                let ds = Dataset {
                    x: train
                        .iter()
                        .flat_map(|&i| x[i * nf..(i + 1) * nf].iter().copied())
                        .collect(),
                    y: train.iter().map(|&i| u8::from(target[i])).collect(),
                    groups: train.iter().map(|&i| (key[i] >> 32) as u32).collect(),
                    n_features: nf,
                    feature_names: names.to_vec(),
                };
                let mut pr = params();
                if round > 0 {
                    pr.n_rounds = WARM_ROUNDS;
                }
                let init = prev[f as usize].clone();
                let m = train_gbdt_from(&ds, &pr, seed ^ (round as u64) << 8 ^ f, init)
                    .map_err(|e| format!("fold {f}: {e:?}"))?;
                let held: Vec<usize> = (0..n).filter(|&i| fold[i] == f).collect();
                let out: Vec<(usize, f32)> = held
                    .into_par_iter()
                    .map(|i| (i, m.predict_value(&x[i * nf..(i + 1) * nf])))
                    .collect();
                if out.len() > 1 && out.iter().all(|&(_, s)| s == out[0].1) {
                    return Err(format!("fold {f}: constant score"));
                }
                Ok((out, m.trees.clone()))
            })
            .collect();
        let bad: Vec<&String> = fits.iter().filter_map(|r| r.as_ref().err()).collect();
        if !bad.is_empty() {
            degenerate_folds += bad.len();
            log::info!(
                target: "dia",
                "learned score round {round}: {} of {FOLDS} folds degenerate ({}); {}",
                bad.len(),
                bad[0],
                if round == 0 {
                    format!("{GB_COLUMN} left out of the precursor PIN")
                } else {
                    format!("keeping the round {} score", round - 1)
                }
            );
            break;
        }
        let mut score = vec![0f32; n];
        for (f, (rows, trees)) in fits.into_iter().flatten().enumerate() {
            for (i, s) in rows {
                score[i] = s;
            }
            prev[f] = trees;
        }
        let q = tdc_q(&score, target);
        let cut = ROUND_Q.get(round + 1).copied().unwrap_or(round_q);
        for i in 0..n {
            positive[i] = target[i] && q[i] <= cut;
        }
        log::info!(
            target: "dia",
            "learned score round {round}: {} targets at q <= 0.01 (out of fold) [{:.2}s]",
            (0..n).filter(|&i| target[i] && q[i] <= 0.01).count(),
            t_round.elapsed().as_secs_f64()
        );
        kept = Some(score);
    }
    LearnedScore {
        score: kept,
        degenerate_folds,
    }
}

/// Best pass-1 row per (FileName, Peptide, charge) of the PIN: its numeric columns, `dia_llr`
/// deciding. `run_file` is the FileName of each pass-2 run, read from the PIN rows of the run's
/// own nominees (their SpecIds are pass-1 SpecIds).
struct Pass1 {
    names: Vec<String>,
    rows: HashMap<(String, String, u8), Vec<f32>>,
    run_file: HashMap<u32, String>,
}

impl Pass1 {
    /// The pass-1 row of each pass-2 row: same run, peptide and charge.
    fn of_rows(&self, t: &Pass2Table) -> Vec<Option<&Vec<f32>>> {
        (0..t.len())
            .map(|i| {
                let file = self.run_file.get(&t.run[i])?;
                self.rows
                    .get(&(file.clone(), t.peptide[i].clone(), t.charge[i]))
            })
            .collect()
    }
}

/// Reads the pass-1 PIN. `own` maps the SpecId of every pass-2 row to its run.
fn read_pass1(pin: &Path, own: &HashMap<&str, u32>) -> Result<Pass1, String> {
    let err = |e: &dyn std::fmt::Display| format!("{}: {e}", pin.display());
    let f = std::fs::File::open(pin).map_err(|e| err(&e))?;
    let mut lines = BufReader::new(f).lines();
    let header = lines.next().ok_or_else(|| err(&"empty PIN"))?.map_err(|e| err(&e))?;
    let h: Vec<&str> = header.split('\t').collect();
    let c_pep = h.iter().position(|c| *c == "Peptide").ok_or_else(|| err(&"no Peptide"))?;
    let c_file = h.iter().position(|c| *c == "FileName");
    let charge: Vec<(usize, u8)> = h
        .iter()
        .enumerate()
        .filter_map(|(i, c)| Some((i, c.strip_prefix("charge")?.parse().ok()?)))
        .collect();
    let numeric: Vec<usize> = (0..c_pep)
        .filter(|&i| !PASS1_SKIP.contains(&h[i]) && !h[i].starts_with("charge"))
        .collect();
    let c_llr = numeric.iter().position(|&i| h[i] == "dia_llr");
    let names: Vec<String> = numeric.iter().map(|&i| format!("p1_{}", h[i])).collect();
    let mut rows: HashMap<(String, String, u8), Vec<f32>> = HashMap::new();
    let mut run_file: HashMap<u32, String> = HashMap::new();
    for line in lines {
        let line = line.map_err(|e| err(&e))?;
        if line.is_empty() {
            continue;
        }
        let x: Vec<&str> = line.split('\t').collect();
        if x.len() <= c_pep {
            return Err(err(&"short PIN row"));
        }
        let file = c_file.map_or("", |c| x[c]);
        if let Some(&run) = own.get(x[0]) {
            match run_file.get(&run) {
                Some(f) if f != file => {
                    return Err(err(&format!(
                        "the nominees of pass-2 run {run} come from PIN files {f:?} and {file:?}"
                    )))
                }
                Some(_) => {}
                None => {
                    run_file.insert(run, file.to_string());
                }
            }
        }
        let z = charge.iter().find(|(i, _)| x[*i] == "1").map_or(0, |c| c.1);
        let v: Vec<f32> = numeric
            .iter()
            .map(|&i| x[i].trim().parse::<f32>().unwrap_or(f32::NAN))
            .map(|v| if v.is_finite() { v } else { 0.0 })
            .collect();
        let llr = |v: &[f32]| c_llr.map_or(MISSING_DIA_LLR, |k| v[k]);
        let k = (file.to_string(), x[c_pep].to_string(), z);
        match rows.get_mut(&k) {
            Some(old) if llr(old) >= llr(&v) => {}
            Some(old) => *old = v,
            None => {
                rows.insert(k, v);
            }
        }
    }
    Ok(Pass1 {
        names,
        rows,
        run_file,
    })
}

/// Seed of the learned score: Percolator's `--seed` when given, else 42.
pub(crate) fn percolator_seed(args: &[String]) -> u64 {
    args.iter()
        .position(|a| a == "--seed" || a == "-S")
        .and_then(|i| args.get(i + 1))
        .and_then(|v| v.parse().ok())
        .unwrap_or(42)
}

/// SpecIds of every row of a PIN.
pub(crate) fn pin_ids(pin: &Path) -> Result<std::collections::HashSet<String>, String> {
    let err = |e: &dyn std::fmt::Display| format!("{}: {e}", pin.display());
    let f = std::fs::File::open(pin).map_err(|e| err(&e))?;
    let mut ids = std::collections::HashSet::new();
    for line in BufReader::new(f).lines().skip(1) {
        let line = line.map_err(|e| err(&e))?;
        if let Some(id) = line.split('\t').next().filter(|s| !s.is_empty()) {
            ids.insert(id.to_string());
        }
    }
    Ok(ids)
}

/// The rows written to the precursor PIN and how its learned score trained.
pub(crate) struct PrecursorPin {
    pub(crate) targets: Vec<String>,
    pub(crate) decoys: Vec<String>,
    /// Folds of the learned score that did not train, over its rounds.
    pub(crate) degenerate_folds: usize,
    /// Whether the PIN has the learned-score column.
    pub(crate) learned: bool,
}

/// Writes the precursor PIN `out` from the pass-1 PIN `pin` (they may be the same file): one row
/// per pass-2 precursor with the learned score (left out when its first round is degenerate),
/// `extra` columns (e.g. cross-run), the pass-2 columns and the best pass-1 row.
pub(crate) fn write_precursor_pin(
    pin: &Path,
    out_pin: &Path,
    t: &Pass2Table,
    extra: &[(String, Vec<f32>)],
    seed: u64,
) -> Result<PrecursorPin, String> {
    let t0 = std::time::Instant::now();
    let own: HashMap<&str, u32> = t
        .spec_id
        .iter()
        .zip(&t.run)
        .map(|(s, &r)| (s.as_str(), r))
        .collect();
    let p1 = read_pass1(pin, &own)?;
    drop(own);
    log::info!(
        target: "dia",
        "precursor PIN: read the pass-1 PIN {} [{:.2}s]",
        pin.display(),
        t0.elapsed().as_secs_f64()
    );
    let p1_of = p1.of_rows(t);
    let n = t.len();
    let mut names: Vec<String> = extra.iter().map(|(n, _)| n.clone()).collect();
    names.extend(t.features.iter().map(|(n, _)| n.clone()));
    names.extend((1..=5).map(|z| format!("charge{z}")));
    names.extend(["pep_len".to_string(), "mz".to_string()]);
    let nf = names.len();
    let mut x = vec![0f32; n * nf];
    x.par_chunks_mut(nf).enumerate().for_each(|(i, row)| {
        let mut k = 0;
        for (_, v) in extra.iter().chain(t.features.iter()) {
            row[k] = if v[i].is_finite() { v[i] } else { 0.0 };
            k += 1;
        }
        for z in 1..=5u8 {
            row[k] = f32::from(u8::from(t.charge[i] == z));
            k += 1;
        }
        row[k] = t.seq[i].len() as f32;
        row[k + 1] = t.prec_mz[i] as f32;
    });
    let target: Vec<bool> = t.decoy.iter().map(|d| !d).collect();
    let prefix = decoy_prefix(t);
    let key: Vec<u64> = (0..n)
        .map(|i| fold_key(&t.seq[i], &t.proteins[i], &prefix))
        .collect();
    let learned = learned_score(&x, &names, &target, &key, seed);
    let gb = learned.score.as_deref();
    let t_write = std::time::Instant::now();

    let err = |e: &dyn std::fmt::Display| format!("{}: {e}", out_pin.display());
    let tmp = out_pin.with_extension("pin.partial");
    let mut out = BufWriter::new(std::fs::File::create(&tmp).map_err(|e| err(&e))?);
    let mut head = vec![
        "SpecId".to_string(),
        "Label".into(),
        "ScanNr".into(),
        "FileName".into(),
    ];
    if gb.is_some() {
        head.push(GB_COLUMN.into());
    }
    head.extend(names.iter().cloned());
    head.extend(p1.names.iter().cloned());
    head.extend([PRESENT_COLUMN.to_string(), "Peptide".into(), "Proteins".into()]);
    writeln!(out, "{}", head.join("\t")).map_err(|e| err(&e))?;
    let zeros = vec![0f32; p1.names.len()];
    let mut ids = PrecursorPin {
        targets: Vec::new(),
        decoys: Vec::new(),
        degenerate_folds: learned.degenerate_folds,
        learned: gb.is_some(),
    };
    let mut line = String::new();
    for i in 0..n {
        use std::fmt::Write as _;
        line.clear();
        let p = p1_of[i];
        let _ = write!(
            line,
            "{}\t{}\t{}\t{}",
            t.spec_id[i],
            if t.decoy[i] { -1 } else { 1 },
            i + 1,
            p1.run_file.get(&t.run[i]).map_or("", String::as_str),
        );
        if let Some(gb) = gb {
            let _ = write!(line, "\t{:.5}", gb[i]);
        }
        for v in &x[i * nf..(i + 1) * nf] {
            let _ = write!(line, "\t{v}");
        }
        for v in p.unwrap_or(&zeros) {
            let _ = write!(line, "\t{v}");
        }
        let _ = write!(
            line,
            "\t{}\t{}\t{}",
            u8::from(p.is_some()),
            t.peptide[i],
            t.proteins[i].replace(';', "\t")
        );
        writeln!(out, "{line}").map_err(|e| err(&e))?;
        if t.decoy[i] {
            ids.decoys.push(t.spec_id[i].clone());
        } else {
            ids.targets.push(t.spec_id[i].clone());
        }
    }
    out.flush().map_err(|e| err(&e))?;
    drop(out);
    std::fs::rename(&tmp, out_pin).map_err(|e| err(&e))?;
    log::info!(
        target: "dia",
        "precursor PIN: {n} rows ({} targets), {} learned-score features ({}, {} degenerate \
         folds), {} pass-1 columns, {} with a pass-1 row [written in {:.2}s; total {:.2}s]",
        ids.targets.len(),
        nf,
        if ids.learned { "kept" } else { "score left out" },
        ids.degenerate_folds,
        p1.names.len(),
        p1_of.iter().filter(|p| p.is_some()).count(),
        t_write.elapsed().as_secs_f64(),
        t0.elapsed().as_secs_f64()
    );
    Ok(ids)
}

/// Percolator's results on a PIN: target rows, decoy ids and decoy rows, by SpecId.
pub(crate) struct Rescored {
    pub(crate) targets: HashMap<String, PercolatorPsm>,
    pub(crate) decoy_ids: HashSet<String>,
    pub(crate) decoys: HashMap<String, PercolatorPsm>,
}

/// The rescored precursor PIN, when it is used in place of pass 1.
pub(crate) struct PrecursorRescore {
    pub(crate) pin: PathBuf,
    pub(crate) ids: PrecursorPin,
    pub(crate) rescored: Rescored,
    /// SpecIds of the pass-1 PIN, which the precursor rows replace.
    pub(crate) replaced: HashSet<String>,
    /// Target precursors at q <= `COMPARE_Q`.
    pub(crate) n_q: usize,
}

/// The precursor PIN written beside the pass-1 PIN `pin`.
pub(crate) fn precursor_pin_path(pin: &Path) -> PathBuf {
    pin.with_extension("precursor.pin")
}

/// Target and decoy precursors at precursor q <= `q` among best-per-precursor
/// `(score, is decoy)` pairs: q = (D + 1) / T over the score-ordered pairs, made monotone, with
/// tied scores counting decoys first (the q of [`crate::dia::precursor_q`]).
pub(crate) fn tdc_counts(best: &[(f64, bool)], q: f64) -> (usize, usize) {
    let mut v = best.to_vec();
    v.sort_by(|a, b| b.0.total_cmp(&a.0).then(b.1.cmp(&a.1)));
    let (mut t, mut d, mut at_q) = (0usize, 0usize, (0, 0));
    let mut i = 0;
    while i < v.len() {
        let mut j = i;
        while j < v.len() && v[j].0 == v[i].0 {
            if v[j].1 {
                d += 1;
            } else {
                t += 1;
            }
            j += 1;
        }
        if (d + 1) as f64 / t.max(1) as f64 <= q {
            at_q = (t, d);
        }
        i = j;
    }
    at_q
}

/// Why pass 1 is preferred to the precursor PIN, or `None` when the precursor PIN wins: at each
/// q of `SELECT_Q`, the precursor PIN counted at its own q tightened by one standard error of
/// its decoy count, `q (1 - 1 / sqrt(D + 1))` with D its decoys at q, must keep at least
/// pass 1's target precursors at q. The larger of two noisy counts favours the arm whose decoys
/// came out low; the Poisson error of D + 1 is that fraction of the estimated FDR, so a win
/// smaller than the decoy noise keeps pass 1, and so does a win at only one point of the curve.
fn pass1_preferred(pass1: &[(f64, bool)], precursor: &[(f64, bool)]) -> Option<String> {
    SELECT_Q.iter().find_map(|&q| {
        let n1 = tdc_counts(pass1, q).0;
        let d = tdc_counts(precursor, q).1;
        let strict = q * (1.0 - 1.0 / ((d + 1) as f64).sqrt());
        let n = tdc_counts(precursor, strict).0;
        (n < n1).then(|| {
            format!(
                "{n} target precursors at q <= {strict:.4} (q {q} less one decoy standard \
                 error, {d} decoys), fewer than pass 1's {n1} at q <= {q}"
            )
        })
    })
}

/// Writes the precursor PIN beside the pass-1 PIN and rescores it with `rescore`. Its results
/// replace pass 1 only when they are sound; otherwise the error says why pass 1 is kept:
/// - the pass-2 table has decoys and at least `pass1_q` targets (`pass1_q`: pass 1's target
///   precursors at q <= `COMPARE_Q`), since a PIN with fewer targets cannot match pass 1 and
///   target-decoy competition needs both labels;
/// - Percolator's training is not degenerate (at least two target scores, not every PEP 1);
/// - the precursor-level target-decoy curve beats pass 1's (`pass1`: the best
///   `(score, is decoy)` of every pass-1 precursor) by more than the decoy noise at q 0.005
///   and 0.01 (`pass1_preferred`).
///
/// Any error while writing or rescoring also keeps pass 1.
pub(crate) fn rescore_precursors(
    pin: &Path,
    t: &Pass2Table,
    seed: u64,
    pass1: &[(f64, bool)],
    rescore: impl FnOnce(&Path) -> Result<Rescored, String>,
) -> Result<PrecursorRescore, String> {
    let pass1_q = tdc_counts(pass1, COMPARE_Q).0;
    let n_dec = t.decoy.iter().filter(|&&d| d).count();
    let n_tgt = t.len() - n_dec;
    if n_dec == 0 || n_tgt == 0 || n_tgt < pass1_q {
        return Err(format!(
            "pass 2 has {n_tgt} target and {n_dec} decoy rows, pass 1 {pass1_q} target \
             precursors at q <= {COMPARE_Q}"
        ));
    }
    let replaced = pin_ids(pin)?;
    let out = precursor_pin_path(pin);
    let ids = write_precursor_pin(pin, &out, t, &[], seed)?;
    let rescored = rescore(&out)?;
    let mut scores: Vec<f64> = rescored.targets.values().filter_map(|p| p.score).collect();
    scores.sort_by(f64::total_cmp);
    scores.dedup();
    if scores.len() < 2 || rescored.targets.values().all(|p| p.pep >= 1.0) {
        return Err(format!(
            "Percolator on {} is degenerate ({} distinct target scores)",
            out.display(),
            scores.len()
        ));
    }
    let best: Vec<(f64, bool)> = rescored
        .targets
        .values()
        .map(|p| (p, false))
        .chain(rescored.decoys.values().map(|p| (p, true)))
        .filter_map(|(p, decoy)| p.score.map(|s| (s, decoy)))
        .collect();
    if let Some(why) = pass1_preferred(pass1, &best) {
        return Err(format!("{}: {why}", out.display()));
    }
    let n_q = tdc_counts(&best, COMPARE_Q).0;
    Ok(PrecursorRescore {
        pin: out,
        ids,
        rescored,
        replaced,
        n_q,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn two_runs_with_the_same_peptide_get_their_own_pass1_rows() {
        let dir = tempfile::tempdir().unwrap();
        let pin = dir.path().join("p.pin");
        std::fs::write(
            &pin,
            "SpecId\tLabel\tScanNr\tFileName\tdia_llr\tx\tcharge2\tPeptide\tProteins\n\
             a/s1\t1\t1\ta\t5\t1\t1\tK.PEPTIDEK.A\tP1\n\
             b/s1\t1\t1\tb\t3\t2\t1\tK.PEPTIDEK.A\tP1\n\
             b/s2\t1\t2\tb\t4\t3\t1\tK.PEPTIDEK.A\tP1\n",
        )
        .unwrap();
        let t = Pass2Table {
            spec_id: vec!["b/s1".into(), "a/s1".into(), "xfer_a_b/s9".into()],
            run: vec![1, 0, 1],
            decoy: vec![false; 3],
            charge: vec![2; 3],
            peptide: vec!["K.PEPTIDEK.A".into(); 3],
            proteins: vec!["P1".into(); 3],
            seq: vec!["PEPTIDEK".into(); 3],
            prec_mz: vec![500.0; 3],
            features: Vec::new(),
            ..Default::default()
        };
        let own: HashMap<&str, u32> =
            t.spec_id.iter().zip(&t.run).map(|(s, &r)| (s.as_str(), r)).collect();
        let p1 = read_pass1(&pin, &own).unwrap();
        assert_eq!(p1.names, ["p1_dia_llr", "p1_x"]);
        assert_eq!(p1.run_file[&0], "a");
        assert_eq!(p1.run_file[&1], "b");
        let rows = p1.of_rows(&t);
        assert_eq!(rows[0].unwrap(), &vec![4.0, 3.0], "run b's best row");
        assert_eq!(rows[1].unwrap(), &vec![5.0, 1.0], "run a's own row");
        assert_eq!(rows[2].unwrap(), &vec![4.0, 3.0], "a row scored in run b reads run b");
    }

    #[test]
    fn fold_key_keeps_decoy_with_its_target() {
        let p = "XXX_";
        let target = fold_key("PEPTIDEK", "P1", p);
        assert_eq!(target, fold_key("PEPTLDEK", "XXX_P1", p), "one-residue mutant decoy");
        assert_eq!(target, fold_key("EDITPEPK", "XXX_P1", p), "reversed decoy");
        assert_ne!(target, fold_key("PEPTIDEK", "P2", p));
        let shared = fold_key("SHAREDPEPK", "P2;P1", p);
        assert_eq!(shared, fold_key("SHAREDPEPK", "P1;P2", p));
        assert_eq!(shared, fold_key("SHAREDPEPK", "P1\tP2", p));
        assert_eq!(shared, fold_key("SHARLDPEPK", "XXX_P1;XXX_P2", p), "mutant of a shared peptide");
        assert_eq!(shared, target, "every peptide of P1 shares its fold");
        assert_eq!(fold_key("PEPTIDEK", "", p), fold_key("KEDITPEP", "", p));
    }

    #[test]
    fn decoy_prefix_is_read_from_the_accessions() {
        let t = Pass2Table {
            spec_id: vec!["a".into(), "b".into(), "c".into(), "d".into()],
            run: vec![0; 4],
            decoy: vec![false, false, true, true],
            charge: vec![2; 4],
            peptide: vec!["K.PEPTIDEK.A".into(); 4],
            proteins: vec![
                "sp|P1|A_HUMAN".into(),
                "sp|P2|B_HUMAN".into(),
                "DECOY_sp|P1|A_HUMAN".into(),
                "DECOY_sp|P2|B_HUMAN;DECOY_sp|P3|C_HUMAN".into(),
            ],
            seq: vec!["PEPTIDEK".into(); 4],
            prec_mz: vec![500.0; 4],
            features: Vec::new(),
            ..Default::default()
        };
        assert_eq!(decoy_prefix(&t), "DECOY_");
        let none = Pass2Table {
            decoy: vec![false; 4],
            ..t
        };
        assert_eq!(decoy_prefix(&none), "");
    }

    #[test]
    fn tdc_q_counts_decoys_first_at_ties() {
        let q = tdc_q(&[3.0, 2.0, 2.0, 1.0], &[true, true, false, true]);
        let q_ties_decoy_first = tdc_q(&[2.0, 2.0], &[true, false]);
        assert_eq!(q_ties_decoy_first[0], 1.0);
        for i in [0, 1, 3] {
            assert!((q[i] - 2.0 / 3.0).abs() < 1e-12, "q[{i}] = {}", q[i]);
        }
    }

    fn mix(i: u64) -> u64 {
        i.wrapping_mul(0x2545_f491_4f6c_dd1d).rotate_left(17) ^ 0x9e37_79b9
    }

    fn assert_no_constant_fold(s: &[f32], key: &[u64], seed: u64) {
        for f in 0..FOLDS {
            let v: Vec<f32> = (0..s.len())
                .filter(|&i| (key[i] ^ seed).wrapping_mul(0x9e37_79b9_7f4a_7c15) % FOLDS == f)
                .map(|i| s[i])
                .collect();
            assert!(
                v.len() < 2 || v.iter().any(|&x| x != v[0]),
                "fold {f} has a constant score"
            );
        }
    }

    #[test]
    fn learned_score_is_out_of_fold() {
        let n = 3000u64;
        let target: Vec<bool> = (0..n).map(|i| mix(i) % 2 == 0).collect();
        let key: Vec<u64> = (0..n).map(|i| mix(i + 77_777)).collect();
        let x: Vec<f32> = (0..n).map(|i| (mix(i + 5_000_000) % 10_000) as f32).collect();
        let r = learned_score(&x, &["noise".to_string()], &target, &key, 7);
        let Some(s) = r.score else {
            assert!(r.degenerate_folds > 0, "no score without a degenerate fold");
            return;
        };
        assert_no_constant_fold(&s, &key, 7);
        let (mut hit, mut tot) = (0f64, 0f64);
        for i in 0..n as usize {
            for j in (0..n as usize).step_by(13) {
                if target[i] && !target[j] {
                    tot += 1.0;
                    hit += if s[i] > s[j] { 1.0 } else if s[i] == s[j] { 0.5 } else { 0.0 };
                }
            }
        }
        let auc = hit / tot;
        assert!((auc - 0.5).abs() < 0.06, "out-of-fold AUC of a noise feature {auc}");
    }

    #[test]
    fn a_small_table_falls_back_without_a_constant_score() {
        let n = 500u64;
        let target: Vec<bool> = (0..n).map(|i| mix(i) % 2 == 0).collect();
        let key: Vec<u64> = (0..n).map(|i| mix(i + 77_777)).collect();
        let x: Vec<f32> = (0..n)
            .map(|i| {
                let shift = if target[i as usize] { 60.0 } else { 0.0 };
                (mix(i + 5_000_000) % 100) as f32 + shift
            })
            .collect();
        let r = learned_score(&x, &["signal".to_string()], &target, &key, 7);
        match r.score {
            Some(s) => assert_no_constant_fold(&s, &key, 7),
            None => assert_eq!(r.degenerate_folds, FOLDS as usize, "every fold of round 0"),
        }
    }

    #[test]
    fn a_degenerate_later_round_keeps_the_previous_score() {
        let n = 6000u64;
        let target: Vec<bool> = (0..n).map(|i| mix(i) % 2 == 0).collect();
        let strong: Vec<bool> = (0..n as usize)
            .map(|i| target[i] && mix(i as u64 + 9) % 5 < 2)
            .collect();
        let key: Vec<u64> = (0..n).map(|i| mix(i + 77_777)).collect();
        let x: Vec<f32> = (0..n as usize)
            .map(|i| {
                let shift = if strong[i] { 1000.0 } else { 0.0 };
                (mix(i as u64 + 5_000_000) % 100) as f32 + shift
            })
            .collect();
        let r = learned_score(&x, &["signal".to_string()], &target, &key, 7);
        let s = r.score.expect("round 0 trains on every target");
        assert!(r.degenerate_folds > 0, "round 1 has too few positives per fold");
        assert_no_constant_fold(&s, &key, 7);
    }

    fn pass1_pin(dir: &Path) -> PathBuf {
        let pin = dir.join("run.pin");
        std::fs::write(
            &pin,
            "SpecId\tLabel\tScanNr\tFileName\tdia_llr\tcharge2\tPeptide\tProteins\n\
             s1\t1\t1\trun\t5\t1\tK.PEPTIDEK.A\tP1\n\
             s2\t-1\t2\trun\t1\t1\tK.EDITPEPK.A\tXXX_P1\n",
        )
        .unwrap();
        pin
    }

    fn pass2(rows: usize) -> Pass2Table {
        Pass2Table {
            spec_id: (0..rows).map(|i| format!("s{}", i + 1)).collect(),
            run: vec![0; rows],
            decoy: (0..rows).map(|i| i % 2 == 1).collect(),
            charge: vec![2; rows],
            peptide: (0..rows).map(|i| format!("K.PEPT{i}K.A")).collect(),
            proteins: (0..rows)
                .map(|i| if i % 2 == 1 { "XXX_P1".into() } else { "P1".into() })
                .collect(),
            seq: vec!["PEPTIDEK".into(); rows],
            prec_mz: vec![500.0; rows],
            features: vec![("p2_llr".to_string(), (0..rows).map(|i| i as f32).collect())],
            ..Default::default()
        }
    }

    /// Percolator results with one target row per score, all at q 0.001, and one decoy row per
    /// `decoys` score.
    fn rescored_with(scores: &[f64], decoys: &[f64]) -> Rescored {
        let psm = |id: String, score: f64| PercolatorPsm {
            psm_id: id,
            q_value: 0.001,
            pooled_q_value: None,
            pep: 0.01,
            score: Some(score),
            peptide: "K.PEPTIDEK.A".into(),
            proteins: "P1".into(),
        };
        let targets = scores
            .iter()
            .enumerate()
            .map(|(i, &score)| {
                let id = format!("s{}", 2 * i + 1);
                (id.clone(), psm(id, score))
            })
            .collect();
        let decoys: HashMap<String, PercolatorPsm> = decoys
            .iter()
            .enumerate()
            .map(|(i, &score)| {
                let id = format!("s{}", 2 * i + 2);
                (id.clone(), psm(id, score))
            })
            .collect();
        Rescored {
            targets,
            decoy_ids: decoys.keys().cloned().collect(),
            decoys,
        }
    }

    fn rescored(scores: &[f64]) -> Rescored {
        rescored_with(scores, &[])
    }

    /// `n` best-per-precursor targets with `top` decoys under the best target and `c r^2` more
    /// decoys above the r-th target, so the estimated FDR at the r-th target is about
    /// `(top + 1) / r + c r`.
    fn curve(n: usize, c: f64, top: usize) -> Vec<(f64, bool)> {
        let mut v = Vec::new();
        let mut d = 0usize;
        for r in 1..=n {
            v.push((-(r as f64), false));
            let want = if r == 1 { top } else { top + (c * (r * r) as f64) as usize };
            while d < want {
                v.push((-(r as f64) - 0.5, true));
                d += 1;
            }
        }
        v
    }

    fn clean(n: usize) -> Vec<(f64, bool)> {
        (0..n).map(|i| (-(i as f64), false)).collect()
    }

    #[test]
    fn tdc_counts_use_d_plus_one_over_t_monotone_with_tied_decoys_first() {
        assert_eq!(tdc_counts(&clean(200), 0.005), (200, 0));
        assert_eq!(tdc_counts(&clean(199), 0.005), (0, 0));
        let mut v = clean(300);
        v.push((-0.5, true));
        assert_eq!(tdc_counts(&v, 0.01), (300, 1), "q 2/300 after the decoy");
        v.push((-299.0, true));
        assert_eq!(
            tdc_counts(&v, 0.009),
            (299, 1),
            "the tied decoy counts with its target, so q 3/300 holds at the last target"
        );
        assert_eq!(tdc_counts(&[], 0.01), (0, 0));
    }

    #[test]
    fn the_precursor_pin_must_win_by_more_than_the_decoy_noise_at_both_q() {
        let pass1 = curve(12_000, 1e-6, 0);
        assert!(pass1_preferred(&pass1, &pass1).is_some(), "a tie keeps pass 1");
        let close = curve(12_000, 0.95e-6, 0);
        for q in SELECT_Q {
            assert!(tdc_counts(&close, q).0 > tdc_counts(&pass1, q).0, "ahead at q {q}");
        }
        let why = pass1_preferred(&pass1, &close).expect("within the decoy noise");
        assert!(why.contains("fewer than pass 1"), "{why}");
        let clear = curve(25_000, 0.5e-6, 0);
        assert_eq!(pass1_preferred(&pass1, &clear), None);
        let top_heavy = curve(25_000, 0.5e-6, 30);
        assert!(tdc_counts(&top_heavy, 0.01).0 > 13_000);
        assert_eq!(tdc_counts(&top_heavy, 0.005).0, 0);
        let why = pass1_preferred(&pass1, &top_heavy).expect("behind at q 0.005");
        assert!(why.contains("at q <= 0.005"), "{why}");
        assert_eq!(pass1_preferred(&[], &clean(150)), None, "pass 1 has none at q 0.005");
    }

    #[test]
    fn an_empty_or_tiny_pass2_table_keeps_pass1() {
        let dir = tempfile::tempdir().unwrap();
        let pin = pass1_pin(dir.path());
        let never = |_: &Path| -> Result<Rescored, String> { panic!("no rescoring expected") };
        assert!(rescore_precursors(&pin, &pass2(0), 42, &[], never).is_err());
        assert!(rescore_precursors(&pin, &pass2(8), 42, &clean(200), never).is_err());
        assert!(!precursor_pin_path(&pin).exists());
        let before = std::fs::read_to_string(&pin).unwrap();
        let fewer = |_: &Path| -> Result<Rescored, String> {
            Ok(rescored_with(&[0.0, 1.0, 2.0], &[0.5]))
        };
        let why = rescore_precursors(&pin, &pass2(250), 42, &clean(100), fewer)
            .err()
            .unwrap();
        assert!(why.contains("fewer than pass 1"), "{why}");
        assert_eq!(std::fs::read_to_string(&pin).unwrap(), before, "pass-1 PIN untouched");
        let head = std::fs::read_to_string(precursor_pin_path(&pin)).unwrap();
        assert!(!head.starts_with("SpecId\tLabel\tScanNr\tFileName\tdia_gb_score"));
        let failing = |_: &Path| -> Result<Rescored, String> { Err("Percolator failed".into()) };
        assert!(rescore_precursors(&pin, &pass2(8), 42, &[], failing).is_err());
    }

    #[test]
    fn a_sound_precursor_pin_replaces_pass1() {
        let dir = tempfile::tempdir().unwrap();
        let pin = pass1_pin(dir.path());
        let scores: Vec<f64> = (0..150).map(f64::from).collect();
        let more = |p: &Path| -> Result<Rescored, String> {
            assert_eq!(p, precursor_pin_path(&pin));
            Ok(rescored(&scores))
        };
        let r = rescore_precursors(&pin, &pass2(8), 42, &[], more).unwrap();
        assert_eq!(r.n_q, 150);
        assert_eq!(r.ids.targets.len(), 4);
        assert!(!r.ids.learned, "8 rows are too few to train the learned score");
        assert_eq!(r.replaced.len(), 2);
        let degenerate = |_: &Path| -> Result<Rescored, String> {
            Ok(rescored(&[1.0; 4]))
        };
        let why = rescore_precursors(&pin, &pass2(8), 42, &[], degenerate).err().unwrap();
        assert!(why.contains("degenerate"), "{why}");
    }

    #[test]
    fn min_fold_positives_follows_the_auc_gate() {
        assert_eq!(min_fold_positives(&params()), 910);
    }
}
