//! DIA precursor rows for rescoring: one row per nominated precursor with its pass-2 evidence,
//! the best pass-1 PIN row of the same precursor, and a cross-fitted learned score.
//!
//! The learned score is a gradient-boosted classifier trained out of fold: folds are keyed by
//! the residue composition of the peptide, so a shuffled or reversed decoy always falls in the
//! fold of its source target, and no row is scored by a model that saw it.

use std::collections::HashMap;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use model_train::gbdt::train::{train_gbdt, Dataset, TrainParams};
use output::dia_pass2_table::Pass2Table;
use rayon::prelude::*;

pub(crate) const GB_COLUMN: &str = "dia_gb_score";
const FOLDS: u64 = 3;
const ROUND_Q: [f64; 3] = [f64::INFINITY, 0.05, 0.01];
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

/// Fold key of a peptide: a hash of its sorted residue letters, equal for a target and every
/// permutation of it (shuffled or reversed decoys).
pub(crate) fn fold_key(seq: &str) -> u64 {
    let mut b: Vec<u8> = seq.bytes().filter(u8::is_ascii_uppercase).collect();
    b.sort_unstable();
    model::hash::fnv1a64_extend(model::hash::FNV1A64_OFFSET, b.iter())
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

/// Out-of-fold learned score of every row: the raw boosted margin (log-odds), the same scale in
/// every fold. `x` is row-major with `names.len()` columns.
pub(crate) fn learned_score(
    x: &[f32],
    names: &[String],
    target: &[bool],
    key: &[u64],
    seed: u64,
) -> Result<Vec<f32>, String> {
    let n = target.len();
    let nf = names.len();
    let fold: Vec<u64> = key
        .iter()
        .map(|k| (k ^ seed).wrapping_mul(0x9e37_79b9_7f4a_7c15) % FOLDS)
        .collect();
    let mut positive: Vec<bool> = target.to_vec();
    let mut score = vec![0f32; n];
    for (round, &q_cut) in ROUND_Q.iter().enumerate() {
        let parts: Vec<Result<Vec<(usize, f32)>, String>> = (0..FOLDS)
            .into_par_iter()
            .map(|f| {
                let train: Vec<usize> = (0..n)
                    .filter(|&i| fold[i] != f && (!target[i] || positive[i]))
                    .collect();
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
                let m = train_gbdt(&ds, &params(), seed ^ (round as u64) << 8 ^ f)
                    .map_err(|e| format!("learned DIA score, round {round} fold {f}: {e:?}"))?;
                Ok((0..n)
                    .filter(|&i| fold[i] == f)
                    .map(|i| (i, m.predict_value(&x[i * nf..(i + 1) * nf])))
                    .collect())
            })
            .collect();
        for p in parts {
            for (i, s) in p? {
                score[i] = s;
            }
        }
        let q = tdc_q(&score, target);
        let cut = ROUND_Q.get(round + 1).copied().unwrap_or(q_cut);
        for i in 0..n {
            positive[i] = target[i] && q[i] <= cut;
        }
        log::info!(
            target: "dia",
            "learned score round {round}: {} targets at q <= 0.01 (out of fold)",
            (0..n).filter(|&i| target[i] && q[i] <= 0.01).count()
        );
    }
    Ok(score)
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

/// Ids of the rows written to the precursor PIN: (targets, decoys).
pub(crate) struct PrecursorPin {
    pub(crate) targets: Vec<String>,
    pub(crate) decoys: Vec<String>,
}

/// Rewrites `pin` (the pass-1 PIN) as the precursor PIN: one row per pass-2 precursor with the
/// learned score, `extra` columns (e.g. cross-run), the pass-2 columns and the best pass-1 row.
pub(crate) fn write_precursor_pin(
    pin: &Path,
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
    let key: Vec<u64> = t.seq.iter().map(|s| fold_key(s)).collect();
    let gb = learned_score(&x, &names, &target, &key, seed)?;

    let err = |e: &dyn std::fmt::Display| format!("{}: {e}", pin.display());
    let tmp = pin.with_extension("pin.partial");
    let mut out = BufWriter::new(std::fs::File::create(&tmp).map_err(|e| err(&e))?);
    let mut head = vec![
        "SpecId".to_string(),
        "Label".into(),
        "ScanNr".into(),
        "FileName".into(),
        GB_COLUMN.into(),
    ];
    head.extend(names.iter().cloned());
    head.extend(p1.names.iter().cloned());
    head.extend(["p1_present".to_string(), "Peptide".into(), "Proteins".into()]);
    writeln!(out, "{}", head.join("\t")).map_err(|e| err(&e))?;
    let zeros = vec![0f32; p1.names.len()];
    let mut ids = PrecursorPin {
        targets: Vec::new(),
        decoys: Vec::new(),
    };
    let mut line = String::new();
    for i in 0..n {
        use std::fmt::Write as _;
        line.clear();
        let p = p1_of[i];
        let _ = write!(
            line,
            "{}\t{}\t{}\t{}\t{:.5}",
            t.spec_id[i],
            if t.decoy[i] { -1 } else { 1 },
            i + 1,
            p1.run_file.get(&t.run[i]).map_or("", String::as_str),
            gb[i]
        );
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
    std::fs::rename(&tmp, pin).map_err(|e| err(&e))?;
    log::info!(
        target: "dia",
        "precursor PIN: {n} rows ({} targets), {} learned-score features, {} pass-1 columns, \
         {} with a pass-1 row [{:.2}s]",
        ids.targets.len(),
        nf,
        p1.names.len(),
        p1_of.iter().filter(|p| p.is_some()).count(),
        t0.elapsed().as_secs_f64()
    );
    Ok(ids)
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
        assert_eq!(fold_key("PEPTIDEK"), fold_key("EDITPEPK"));
        assert_eq!(fold_key("PEPTIDEK"), fold_key("KEDITPEP"));
        assert_ne!(fold_key("PEPTIDEK"), fold_key("PEPTIDER"));
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

    #[test]
    fn learned_score_is_out_of_fold() {
        let n = 3000u64;
        let mix = |i: u64| i.wrapping_mul(0x2545_f491_4f6c_dd1d).rotate_left(17) ^ 0x9e37_79b9;
        let target: Vec<bool> = (0..n).map(|i| mix(i) % 2 == 0).collect();
        let key: Vec<u64> = (0..n).map(|i| mix(i + 77_777)).collect();
        let x: Vec<f32> = (0..n).map(|i| (mix(i + 5_000_000) % 10_000) as f32).collect();
        let s = learned_score(&x, &["noise".to_string()], &target, &key, 7).unwrap();
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
}
