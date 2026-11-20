//! Study mode of `andes --dia`: cross-run columns over the pass-2 tables of every run.
//!
//! With `--quant-dir`, each per-run `andes --dia` process stores its pass-1 PIN and its pass-2
//! table in `<quant-dir>/dia_pass2/` (`<stem>.pin` and `<stem>.parquet`, the run's rescoring
//! settings in the table's metadata) instead of writing the precursor PIN, and reports its
//! pass-1 results. `andes quant finalize` on that directory gives every run's precursors the
//! cross-run columns of the other runs, then per run writes the precursor PIN (learned score,
//! cross-run and pass-2 columns, best pass-1 row), runs the native Percolator engine (with the
//! runs' storage budget when they were rescored natively) and the precursor-level q, and infers
//! the protein groups of all runs.
//!
//! Cross-run columns of a precursor (peptide string, charge, target or decoy) of run r, from
//! every other run o that nominated the same key (first row of the key in o's table):
//! - `xr_max`: the largest `p2_llr` of the other runs ([`ABSENT_LLR`] when none has the key);
//! - `xr_n`: the number of other runs with `p2_llr` > [`DETECTED_LLR`];
//! - `xr_rtd`: the smallest |`p2_rtstar` of o - `p2_rtstar` of the row|, capped at
//!   [`MAX_RTD_MIN`] minutes (also the value when no other run has the key). Retention times
//!   are compared unaligned.
//!
//! A decoy is keyed by its own (shuffled) peptide string, so it only ever meets the same decoy
//! in other runs, never its source target, and targets and decoys get the columns by the same
//! rule.

use std::collections::BTreeSet;
use std::hash::{BuildHasher, BuildHasherDefault, DefaultHasher, Hash, Hasher};
use std::path::{Path, PathBuf};

use output::dia_pass2_table::{
    read_pass2, read_pass2_keys, read_pass2_meta, write_pass2_with, Pass2Keys, Pass2Table,
};
use output::PercolatorBackend;
use rayon::prelude::*;

/// Subdirectory of the quant directory with every DIA run's pass-1 PIN and pass-2 table.
pub(crate) const PASS2_DIR: &str = "dia_pass2";
/// Metadata key of the pass-2 table that holds the run's Percolator backend and arguments.
const RESCORE_KEY: &str = "andes.dia.rescore";
/// `xr_max` of a precursor no other run nominated (a `p2_llr`, natural log units).
const ABSENT_LLR: f32 = -30.0;
/// `p2_llr` above which another run counts toward `xr_n`.
const DETECTED_LLR: f32 = 8.0;
/// Cap of `xr_rtd` and its value when no other run has the precursor, in minutes.
const MAX_RTD_MIN: f32 = 5.0;
/// Median `xr_rtd` of the detected targets above which the runs' retention times likely differ
/// by more than replicate drift, in minutes (a warning only; the columns stay unaligned).
const RTD_WARN_MIN: f32 = 1.0;
/// Names of the cross-run columns, in the order they enter the precursor PIN.
pub(crate) const XR_COLUMNS: [&str; 3] = ["xr_max", "xr_n", "xr_rtd"];
const LLR: &str = "p2_llr";
const RTSTAR: &str = "p2_rtstar";

/// Join key of a precursor across runs: a 128-bit hash of (peptide string, charge, decoy).
pub(crate) fn precursor_key(peptide: &str, charge: u8, decoy: bool) -> u128 {
    let half = |salt: u8| {
        let mut h = BuildHasherDefault::<DefaultHasher>::default().build_hasher();
        salt.hash(&mut h);
        peptide.hash(&mut h);
        charge.hash(&mut h);
        decoy.hash(&mut h);
        h.finish()
    };
    (u128::from(half(0)) << 64) | u128::from(half(1))
}

/// One run's precursors for the join: `p2_llr` and `p2_rtstar` of the first row of every key
/// (in table order), sorted by key. About 24 bytes per precursor.
pub(crate) struct RunKeys {
    key: Vec<u128>,
    llr: Vec<f32>,
    rt: Vec<f32>,
}

impl RunKeys {
    pub(crate) fn new(t: &Pass2Keys) -> Self {
        let keys: Vec<u128> = (0..t.peptide.len())
            .into_par_iter()
            .map(|i| precursor_key(&t.peptide[i], t.charge[i], t.decoy[i]))
            .collect();
        let mut order: Vec<u32> = (0..keys.len() as u32).collect();
        order.sort_by_key(|&i| (keys[i as usize], i));
        order.dedup_by_key(|i| keys[*i as usize]);
        RunKeys {
            key: order.iter().map(|&i| keys[i as usize]).collect(),
            llr: order.iter().map(|&i| t.features[0][i as usize]).collect(),
            rt: order.iter().map(|&i| t.features[1][i as usize]).collect(),
        }
    }

    pub(crate) fn len(&self) -> usize {
        self.key.len()
    }

    fn get(&self, k: u128) -> Option<(f32, f32)> {
        self.key
            .binary_search(&k)
            .ok()
            .map(|i| (self.llr[i], self.rt[i]))
    }
}

/// The cross-run columns of every row of `this` from the other runs (module doc). A
/// non-finite `p2_llr` of another run counts as absent; a non-finite `p2_rtstar` of the row or
/// of another run gives that run the cap.
pub(crate) fn cross_run_columns(
    this: &Pass2Table,
    others: &[&RunKeys],
) -> Result<Vec<(String, Vec<f32>)>, String> {
    let rt = this
        .feature(RTSTAR)
        .ok_or_else(|| format!("pass-2 table has no {RTSTAR} column"))?;
    let cols: Vec<[f32; 3]> = (0..this.len())
        .into_par_iter()
        .map(|i| {
            let k = precursor_key(&this.peptide[i], this.charge[i], this.decoy[i]);
            let (mut best, mut n, mut rtd) = (ABSENT_LLR, 0u32, MAX_RTD_MIN);
            for o in others {
                let Some((llr, rt_o)) = o.get(k) else {
                    continue;
                };
                if !llr.is_finite() {
                    continue;
                }
                best = best.max(llr);
                n += u32::from(llr > DETECTED_LLR);
                let d = (rt_o - rt[i]).abs();
                if d.is_finite() {
                    rtd = rtd.min(d);
                }
            }
            [best, n as f32, rtd]
        })
        .collect();
    Ok(XR_COLUMNS
        .iter()
        .enumerate()
        .map(|(j, name)| (name.to_string(), cols.iter().map(|c| c[j]).collect()))
        .collect())
}

/// The run's Percolator backend and arguments as stored in its pass-2 table.
fn encode_rescore(backend: &PercolatorBackend, args: &[String]) -> String {
    let head = match backend {
        PercolatorBackend::Binary(p) => format!("binary\n{}", p.display()),
        PercolatorBackend::Docker { image } => format!("docker\n{image}"),
        PercolatorBackend::Native { memory } => format!(
            "native\n{}\n{}",
            memory.bytes,
            memory.spill_dir.display()
        ),
    };
    format!("{head}\n{}", args.join(" "))
}

fn decode_rescore(s: &str) -> Result<(PercolatorBackend, Vec<String>), String> {
    let f: Vec<&str> = s.split('\n').collect();
    let bad = || format!("unreadable {RESCORE_KEY} metadata {s:?}");
    let args = |a: &str| a.split_whitespace().map(str::to_string).collect::<Vec<_>>();
    match f.as_slice() {
        ["binary", path, a] => Ok((PercolatorBackend::Binary(PathBuf::from(path)), args(a))),
        ["docker", image, a] => Ok((
            PercolatorBackend::Docker {
                image: image.to_string(),
            },
            args(a),
        )),
        ["native", bytes, spill, a] => Ok((
            PercolatorBackend::Native {
                memory: output::MemLimit {
                    bytes: bytes.parse().map_err(|_| bad())?,
                    spill_dir: PathBuf::from(spill),
                },
            },
            args(a),
        )),
        _ => Err(bad()),
    }
}

/// Stores a run of a DIA study: its pass-1 PIN as `<stem>.pin` and its pass-2 table as
/// `<stem>.parquet` in `<quant_dir>/dia_pass2/`, the table last (its presence marks a complete
/// run). `rescore` is the run's Percolator backend and arguments; `None` without rescoring.
pub(crate) fn store_run(
    quant_dir: &Path,
    spectrum: &Path,
    pin: &Path,
    table: &Pass2Table,
    rescore: Option<(&PercolatorBackend, &[String])>,
) -> Result<(), String> {
    let runs: BTreeSet<u32> = table.run.iter().copied().collect();
    if runs.len() > 1 {
        return Err(format!(
            "--dia --quant-dir stores one run per process; this pass-2 table has {} runs",
            runs.len()
        ));
    }
    let dir = quant_dir.join(PASS2_DIR);
    std::fs::create_dir_all(&dir).map_err(|e| format!("create {}: {e}", dir.display()))?;
    let stem = crate::quant_run::run_stem(spectrum);
    let pin_out = dir.join(format!("{stem}.pin"));
    let partial = dir.join(format!("{stem}.pin.partial"));
    std::fs::copy(pin, &partial)
        .map_err(|e| format!("copy {} to {}: {e}", pin.display(), partial.display()))?;
    std::fs::rename(&partial, &pin_out)
        .map_err(|e| format!("rename {}: {e}", partial.display()))?;
    let meta: Vec<(&str, String)> = rescore
        .map(|(b, a)| (RESCORE_KEY, encode_rescore(b, a)))
        .into_iter()
        .collect();
    let table_out = dir.join(format!("{stem}.parquet"));
    write_pass2_with(&table_out, table, &meta)?;
    log::info!(
        target: "dia",
        "study mode: stored {} and {} ({} pass-2 precursors); `andes quant finalize --quant-dir \
         {}` adds the cross-run columns and rescores every run",
        pin_out.display(),
        table_out.display(),
        table.len(),
        quant_dir.display()
    );
    Ok(())
}

/// Stems of the DIA runs stored in a quant directory, sorted; empty when it holds none.
pub(crate) fn study_runs(quant_dir: &Path) -> Result<Vec<String>, String> {
    let dir = quant_dir.join(PASS2_DIR);
    if !dir.is_dir() {
        return Ok(Vec::new());
    }
    let mut stems: Vec<String> = std::fs::read_dir(&dir)
        .map_err(|e| format!("{}: {e}", dir.display()))?
        .filter_map(|e| e.ok())
        .filter_map(|e| {
            let name = e.file_name().to_string_lossy().into_owned();
            name.strip_suffix(".parquet").map(str::to_string)
        })
        .collect();
    stems.sort();
    Ok(stems)
}

/// Replaces any seed in the rescoring arguments by `--seed <seed>`.
fn set_seed(args: &mut Vec<String>, seed: u64) {
    let mut i = 0;
    while i < args.len() {
        if args[i] == "--seed" || args[i] == "-S" {
            args.drain(i..(i + 2).min(args.len()));
        } else {
            i += 1;
        }
    }
    args.extend(["--seed".to_string(), seed.to_string()]);
}

fn with_suffix(base: &Path, suffix: &str) -> PathBuf {
    let mut p = base.as_os_str().to_owned();
    p.push(suffix);
    PathBuf::from(p)
}

/// Median of the finite values; `None` when there are none.
fn median(mut v: Vec<f32>) -> Option<f32> {
    v.retain(|x| x.is_finite());
    if v.is_empty() {
        return None;
    }
    let m = v.len() / 2;
    v.select_nth_unstable_by(m, f32::total_cmp);
    Some(v[m])
}

/// `andes quant finalize` on a DIA study (module doc). Writes, per run, `<base>.<stem>.pin`
/// and its Percolator results beside it, and with `--database` `<base>.pg.parquet`.
pub(crate) fn finalize(
    a: &crate::cli::QuantFinalizeArgs,
    stems: &[String],
) -> Result<(), String> {
    let t0 = std::time::Instant::now();
    if a.output_parquet.is_some() {
        return Err(
            "--output-parquet is not available for a DIA study; give --database for the \
             protein groups (<base>.pg.parquet)"
                .into(),
        );
    }
    let dir = a.quant_dir.join(PASS2_DIR);
    let table_path = |s: &str| dir.join(format!("{s}.parquet"));
    let mut rescore: Option<(String, &str)> = None;
    for s in stems {
        let p = table_path(s);
        let meta = read_pass2_meta(&p)?;
        let r = meta
            .into_iter()
            .find(|(k, _)| k == RESCORE_KEY)
            .map(|(_, v)| v)
            .ok_or_else(|| {
                format!(
                    "{}: run {s} was searched without --rescore or --rescore-native; a DIA study \
                     is rescored with the settings of its runs",
                    p.display()
                )
            })?;
        match &rescore {
            None => rescore = Some((r, s.as_str())),
            Some((r0, s0)) if *r0 != r => {
                return Err(format!(
                    "runs {s0} and {s} were searched with different Percolator settings \
                     ({r0:?} and {r:?}); search every run of a study with the same ones"
                ));
            }
            _ => {}
        }
    }
    let Some((rescore, _)) = rescore else {
        return Err(format!("{}: no DIA runs", dir.display()));
    };
    let (searched, searched_args) = decode_rescore(&rescore)?;
    let seed = a.seed;
    let (backend, mut args) = match searched {
        PercolatorBackend::Native { memory } => (PercolatorBackend::Native { memory }, searched_args),
        other => {
            let spill_dir = a
                .output
                .parent()
                .filter(|d| !d.as_os_str().is_empty())
                .map_or_else(|| PathBuf::from("."), Path::to_path_buf);
            log::info!(
                target: "dia",
                "the runs were rescored with Percolator {}; the study is rescored with the \
                 native engine ({} GB storage budget, spill to {})",
                other.describe(),
                crate::cli::DEFAULT_RESCORE_MEMORY_GB,
                spill_dir.display()
            );
            let memory = output::MemLimit {
                bytes: crate::cli::DEFAULT_RESCORE_MEMORY_GB.saturating_mul(1 << 30),
                spill_dir,
            };
            (PercolatorBackend::Native { memory }, Vec::new())
        }
    };
    crate::dia::default_train_fdr(&mut args);
    set_seed(&mut args, seed);
    if let PercolatorBackend::Native { memory } = &backend {
        std::fs::create_dir_all(&memory.spill_dir)
            .map_err(|e| format!("create {}: {e}", memory.spill_dir.display()))?;
    }
    let threads = a.threads.unwrap_or_else(rayon::current_num_threads);
    let extra = output::percolator_args_with_threads(&args, threads);
    log::info!(
        target: "dia",
        "DIA study: {} runs, Percolator {} with {:?}, learned-score seed {seed}",
        stems.len(),
        backend.describe(),
        args
    );

    let mut keys: Vec<RunKeys> = Vec::with_capacity(stems.len());
    for s in stems {
        let k = read_pass2_keys(&table_path(s), [LLR, RTSTAR])?;
        keys.push(RunKeys::new(&k));
    }
    log::info!(
        target: "dia",
        "cross-run keys: {} precursors over {} runs [{:.2}s]",
        keys.iter().map(RunKeys::len).sum::<usize>(),
        stems.len(),
        t0.elapsed().as_secs_f64()
    );

    let mut all_targets = std::collections::HashMap::new();
    let mut all_decoys = std::collections::HashMap::new();
    for (r, s) in stems.iter().enumerate() {
        let t_run = std::time::Instant::now();
        let t = read_pass2(&table_path(s), None)?;
        let others: Vec<&RunKeys> = keys
            .iter()
            .enumerate()
            .filter(|(o, _)| *o != r)
            .map(|(_, k)| k)
            .collect();
        let xr = cross_run_columns(&t, &others)?;
        let n_seen = xr[1].1.iter().filter(|&&v| v > 0.0).count();
        let own = t.feature(LLR).unwrap_or(&[]);
        let rtd_detected = median(
            (0..t.len())
                .filter(|&i| !t.decoy[i] && xr[1].1[i] > 0.0 && own.get(i).is_some_and(|&l| l > DETECTED_LLR))
                .map(|i| xr[2].1[i])
                .collect(),
        );
        log::info!(
            target: "dia",
            "run {s}: {} precursors, {n_seen} detected in another run (p2_llr > {DETECTED_LLR}), \
             median xr_rtd of the detected targets {}",
            t.len(),
            rtd_detected.map_or("n/a".to_string(), |m| format!("{m:.3} min"))
        );
        if rtd_detected.is_some_and(|m| m > RTD_WARN_MIN) {
            log::warn!(
                "run {s}: the detected targets are {:.2} min (median) from their apex in the \
                 nearest other run; the cross-run columns compare unaligned retention times",
                rtd_detected.unwrap_or(0.0)
            );
        }
        let pin_in = dir.join(format!("{s}.pin"));
        let pin = with_suffix(&a.output, &format!(".{s}.pin"));
        if let Some(parent) = pin.parent().filter(|d| !d.as_os_str().is_empty()) {
            std::fs::create_dir_all(parent)
                .map_err(|e| format!("create {}: {e}", parent.display()))?;
        }
        crate::dia_score::write_precursor_pin(&pin_in, &pin, &t, &xr, seed)?;
        drop(t);
        drop(xr);
        let (mut map, _) = output::run_percolator_with_decoys(&backend, &pin, &extra)
            .map_err(|e| format!("Percolator on {}: {e}", pin.display()))?;
        let dpath = output::decoy_results_path(&pin);
        let text =
            std::fs::read_to_string(&dpath).map_err(|e| format!("{}: {e}", dpath.display()))?;
        let decoys =
            output::parse_psm_results(&text).map_err(|e| format!("{}: {e}", dpath.display()))?;
        drop(text);
        let ((nt, nd, n01), _) = crate::dia::precursor_q(&pin, &[], &mut map, &decoys)?;
        log::info!(
            target: "rescore",
            "run {s}: precursor-level q (DIA): {nt} target and {nd} decoy precursors, {n01} at \
             q <= 0.01 [{:.2}s]",
            t_run.elapsed().as_secs_f64()
        );
        if a.database.is_some() {
            all_targets.extend(map.into_iter().map(|(id, p)| (format!("{s}\t{id}"), p)));
            all_decoys.extend(decoys.into_iter().map(|(id, p)| (format!("{s}\t{id}"), p)));
        }
    }
    drop(keys);
    if let Some(database) = &a.database {
        let inf = crate::inference_cmd::infer_from_maps(
            database,
            &a.proteins,
            all_targets,
            &all_decoys,
        )
        .map_err(|e| e.to_string())?;
        let path = with_suffix(&a.output, ".pg.parquet");
        output::write_pg_parquet(
            &path,
            &inf,
            stems,
            stems.len() == 1,
            a.proteins.protein_fdr,
            None,
        )
        .map_err(|e| format!("{}: {e}", path.display()))?;
        let at = |q: f64| inf.groups.iter().filter(|g| g.q_value <= q).count();
        log::info!(
            target: "inference",
            "wrote {} ({} of {} target groups at protein q<={} over {} runs)",
            path.display(),
            at(a.proteins.protein_fdr),
            inf.groups.len(),
            a.proteins.protein_fdr,
            stems.len()
        );
    }
    log::info!(
        target: "total",
        "DIA study finalized in {:.2} s",
        t0.elapsed().as_secs_f64()
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn table(rows: &[(&str, u8, bool, f32, f32)]) -> Pass2Table {
        Pass2Table {
            spec_id: rows.iter().enumerate().map(|(i, _)| format!("s{i}")).collect(),
            run: vec![0; rows.len()],
            decoy: rows.iter().map(|r| r.2).collect(),
            charge: rows.iter().map(|r| r.1).collect(),
            peptide: rows.iter().map(|r| r.0.to_string()).collect(),
            proteins: vec!["P".to_string(); rows.len()],
            seq: rows.iter().map(|r| r.0.to_string()).collect(),
            prec_mz: vec![500.0; rows.len()],
            features: vec![
                (LLR.to_string(), rows.iter().map(|r| r.3).collect()),
                (RTSTAR.to_string(), rows.iter().map(|r| r.4).collect()),
            ],
        }
    }

    fn keys(rows: &[(&str, u8, bool, f32, f32)]) -> RunKeys {
        let t = table(rows);
        RunKeys::new(&Pass2Keys {
            decoy: t.decoy.clone(),
            charge: t.charge.clone(),
            peptide: t.peptide.clone(),
            features: [
                t.feature(LLR).unwrap().to_vec(),
                t.feature(RTSTAR).unwrap().to_vec(),
            ],
        })
    }

    fn col<'a>(xr: &'a [(String, Vec<f32>)], name: &str) -> &'a [f32] {
        &xr.iter().find(|c| c.0 == name).unwrap().1
    }

    #[test]
    fn cross_run_columns_symmetric_for_decoys() {
        let this = table(&[("PEPTIDEK", 2, false, 10.0, 20.0), ("EPTIPEDK", 2, true, 10.0, 20.0)]);
        let o1 = keys(&[
            ("PEPTIDEK", 2, false, 12.0, 20.5),
            ("EPTIPEDK", 2, true, 3.0, 24.0),
        ]);
        let o2 = keys(&[("PEPTIDEK", 2, true, 50.0, 20.0), ("EPTIPEDK", 2, false, 50.0, 20.0)]);
        let xr = cross_run_columns(&this, &[&o1, &o2]).unwrap();
        assert_eq!(col(&xr, "xr_max"), &[12.0, 3.0]);
        assert_eq!(col(&xr, "xr_n"), &[1.0, 0.0]);
        assert_eq!(col(&xr, "xr_rtd"), &[0.5, 4.0]);
    }

    #[test]
    fn cross_run_fills_absent_runs() {
        let this = table(&[("PEPTIDEK", 2, false, 10.0, 20.0), ("PEPTIDEK", 3, false, 10.0, 20.0)]);
        let o = keys(&[("PEPTIDEK", 3, false, 9.0, 40.0), ("PEPTIDER", 2, false, 9.0, 20.0)]);
        let xr = cross_run_columns(&this, &[&o]).unwrap();
        assert_eq!(col(&xr, "xr_max"), &[ABSENT_LLR, 9.0]);
        assert_eq!(col(&xr, "xr_n"), &[0.0, 1.0]);
        assert_eq!(col(&xr, "xr_rtd"), &[MAX_RTD_MIN, MAX_RTD_MIN]);
        let none = cross_run_columns(&this, &[]).unwrap();
        assert_eq!(col(&none, "xr_max"), &[ABSENT_LLR, ABSENT_LLR]);
        assert_eq!(col(&none, "xr_n"), &[0.0, 0.0]);
        assert_eq!(col(&none, "xr_rtd"), &[MAX_RTD_MIN, MAX_RTD_MIN]);
    }

    #[test]
    fn cross_run_uses_the_first_row_of_a_key_in_another_run() {
        let this = table(&[("PEPTIDEK", 2, false, 10.0, 20.0)]);
        let o = keys(&[
            ("PEPTIDEK", 2, false, 2.0, 21.0),
            ("PEPTIDEK", 2, false, 30.0, 20.0),
        ]);
        let xr = cross_run_columns(&this, &[&o]).unwrap();
        assert_eq!(col(&xr, "xr_max"), &[2.0]);
        assert_eq!(col(&xr, "xr_n"), &[0.0]);
        assert_eq!(col(&xr, "xr_rtd"), &[1.0]);
    }

    #[test]
    fn cross_run_takes_the_best_and_nearest_over_runs() {
        let this = table(&[("PEPTIDEK", 2, false, 10.0, 20.0)]);
        let a = keys(&[("PEPTIDEK", 2, false, 9.0, 23.0)]);
        let b = keys(&[("PEPTIDEK", 2, false, 20.0, 20.25)]);
        let c = keys(&[("PEPTIDEK", 2, false, f32::NAN, 20.0)]);
        let xr = cross_run_columns(&this, &[&a, &b, &c]).unwrap();
        assert_eq!(col(&xr, "xr_max"), &[20.0]);
        assert_eq!(col(&xr, "xr_n"), &[2.0]);
        assert_eq!(col(&xr, "xr_rtd"), &[0.25]);
    }

    #[test]
    fn the_study_seed_replaces_the_searched_one() {
        let mut args = vec!["--seed".to_string(), "7".into(), "--trainFDR".into(), "0.05".into()];
        set_seed(&mut args, 3);
        assert_eq!(args, ["--trainFDR", "0.05", "--seed", "3"]);
        let mut args = Vec::new();
        set_seed(&mut args, 42);
        assert_eq!(args, ["--seed", "42"]);
    }

    #[test]
    fn rescore_settings_round_trip() {
        let args = vec!["--seed".to_string(), "7".into(), "--trainFDR".into(), "0.05".into()];
        let b = PercolatorBackend::Binary(PathBuf::from("/opt/perc/percolator"));
        let (b2, a2) = decode_rescore(&encode_rescore(&b, &args)).unwrap();
        assert_eq!(a2, args);
        assert_eq!(b2.describe(), b.describe());
        let n = PercolatorBackend::Native {
            memory: output::MemLimit {
                bytes: 1 << 33,
                spill_dir: PathBuf::from("/tmp/x"),
            },
        };
        let (n2, _) = decode_rescore(&encode_rescore(&n, &[])).unwrap();
        assert_eq!(n2.describe(), n.describe());
        assert!(decode_rescore("binary").is_err());
    }
}
