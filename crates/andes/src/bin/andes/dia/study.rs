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
//! With cross-run candidate transfer (`andes quant extract --run <stem>` for every run first,
//! see `dia_xfer`), each run's table is the one over its own and its transferred nominees,
//! the cross-run columns are computed over those tables, and the reported q is class-wise: one
//! precursor-level target-decoy competition among the run's own nominees and one among its
//! transferred nominees, so the q holds for each class a user filters (the counts under one
//! pooled competition are logged for comparison).
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
//!
//! Finalize also writes `<base>.dia_quant.parquet` (`output::dia::quantities`): every target
//! precursor of a run at reported q <= `--quant-fdr` with its pass-2 quantities, own and
//! transferred rows alike (one row per peptide and charge, the lowest q), nothing below it.
//! `quantity` is the [`DEFAULT_QUANTITY`] arm; `quantity_normalized` divides it by the run's
//! median-level factor, `2^d` with `d` the median over the precursors quantified in every
//! normalised run of the run's log2 quantity minus the precursor's mean over those runs. A run
//! that shares too few quantified precursors with the others (fewer than
//! [`MIN_SHARED_PRECURSORS`] in every remaining run) is dropped from the shared set one at a
//! time, fewest shared first, and left unscaled.

use std::collections::BTreeSet;
use std::hash::{BuildHasher, BuildHasherDefault, DefaultHasher, Hash, Hasher};
use std::path::{Path, PathBuf};

use output::dia::pass2::{
    read_pass2, read_pass2_keys, read_pass2_meta, write_pass2_with, Pass2Keys, Pass2Table,
};
use output::dia::quantities::{write_dia_quant, DiaQuantRow};
use output::PercolatorBackend;
use rayon::prelude::*;

/// Subdirectory of the quant directory with every DIA run's pass-1 PIN and pass-2 table.
pub(crate) const PASS2_DIR: &str = "dia_pass2";
/// Metadata key of the pass-2 table that holds the run's Percolator backend and arguments.
pub(crate) const RESCORE_KEY: &str = "andes.dia.rescore";
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
/// Quantity arms of the DIA quantity table.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum QuantArm {
    /// Summed area of the best co-eluting fragments (`quant_frag` of pass 2).
    Fragment,
    /// Area of the monoisotopic and M+1 MS1 traces (`quant_ms1` of pass 2).
    Ms1,
}
impl QuantArm {
    /// The arm's column in the DIA quantity table.
    pub(crate) fn column(self) -> &'static str {
        match self {
            QuantArm::Fragment => "quantity_frag",
            QuantArm::Ms1 => "quantity_ms1",
        }
    }
}
/// Every quantity arm, in table column order.
const QUANT_ARMS: [QuantArm; 2] = [QuantArm::Fragment, QuantArm::Ms1];
/// The arm reported as `quantity` (and normalised).
pub(crate) const DEFAULT_QUANTITY: QuantArm = QuantArm::Fragment;
/// Precursors quantified in every normalised run that median-level normalisation needs.
const MIN_SHARED_PRECURSORS: usize = 20;
/// Metadata key of the DIA quantity table naming the default arm.
const DEFAULT_ARM_KEY: &str = "andes.dia_quant.default";

/// A target row of a run's pass-2 table, kept for the quantity table until its q is known.
pub(crate) struct QuantCandidate {
    spec_id: String,
    peptide: String,
    charge: u8,
    proteins: String,
    transfer: bool,
    frag: Option<f64>,
    ms1: Option<f64>,
}

/// The target rows of a pass-2 table with their quantities.
pub(crate) fn quant_candidates(t: &Pass2Table) -> Vec<QuantCandidate> {
    let xfer = t.feature(crate::dia::pass2::XFER_COLUMN);
    (0..t.len())
        .filter(|&i| !t.decoy[i])
        .map(|i| QuantCandidate {
            spec_id: t.spec_id[i].clone(),
            peptide: t.peptide[i].clone(),
            charge: t.charge[i],
            proteins: t.proteins[i].clone(),
            transfer: xfer.is_some_and(|x| x[i] > 0.5),
            frag: t.quant_frag.get(i).copied().flatten(),
            ms1: t.quant_ms1.get(i).copied().flatten(),
        })
        .collect()
}

/// The rows of run `stem` at reported q <= `max_q` (`q` by PIN SpecId), one per peptide and
/// charge (lowest q, then the run's own row, then SpecId), sorted by peptide and charge.
pub(crate) fn reported_rows(
    stem: &str,
    cands: &[QuantCandidate],
    q: &std::collections::HashMap<String, output::PercolatorPsm>,
    max_q: f64,
) -> Vec<DiaQuantRow> {
    let mut best: std::collections::BTreeMap<(&str, u8), (f64, &QuantCandidate)> =
        std::collections::BTreeMap::new();
    for c in cands {
        let Some(p) = q.get(&c.spec_id) else {
            continue;
        };
        if p.q_value.is_nan() || p.q_value > max_q {
            continue;
        }
        let e = best.entry((c.peptide.as_str(), c.charge)).or_insert((p.q_value, c));
        let better = p
            .q_value
            .total_cmp(&e.0)
            .then(c.transfer.cmp(&e.1.transfer))
            .then(c.spec_id.cmp(&e.1.spec_id))
            .is_lt();
        if better {
            *e = (p.q_value, c);
        }
    }
    best.into_values()
        .map(|(q_value, c)| DiaQuantRow {
            filename: stem.to_string(),
            peptide: c.peptide.clone(),
            charge: c.charge,
            proteins: c.proteins.clone(),
            transfer: c.transfer,
            q_value,
            quantity: match DEFAULT_QUANTITY {
                QuantArm::Fragment => c.frag,
                QuantArm::Ms1 => c.ms1,
            },
            quantity_frag: c.frag,
            quantity_ms1: c.ms1,
            quantity_normalized: None,
        })
        .collect()
}

/// Sets `quantity_normalized` of every row (module doc); returns each run's factor (`None`
/// for a run left unscaled), runs by `stems` order.
pub(crate) fn normalise(rows: &mut [DiaQuantRow], stems: &[String]) -> Vec<Option<f64>> {
    use std::collections::BTreeMap;
    let run_of: BTreeMap<&str, usize> = stems.iter().enumerate().map(|(i, s)| (s.as_str(), i)).collect();
    let mut by_key: BTreeMap<(String, u8), Vec<(usize, f64)>> = BTreeMap::new();
    for r in rows.iter() {
        if let (Some(&run), Some(v)) = (run_of.get(r.filename.as_str()), r.quantity) {
            if v > 0.0 && v.is_finite() {
                by_key
                    .entry((r.peptide.clone(), r.charge))
                    .or_default()
                    .push((run, v.log2()));
            }
        }
    }
    let n_runs = stems.len();
    let mut active: BTreeSet<usize> = by_key.values().flatten().map(|c| c.0).collect();
    let shared = loop {
        if active.len() < 2 {
            break Vec::new();
        }
        let shared: Vec<&Vec<(usize, f64)>> = by_key
            .values()
            .filter(|v| active.iter().all(|r| v.iter().any(|c| c.0 == *r)))
            .collect();
        if shared.len() >= MIN_SHARED_PRECURSORS {
            break shared;
        }
        let score = |r: usize| {
            let mut overlap = 0usize;
            let mut total = 0usize;
            for v in by_key.values() {
                if v.iter().any(|c| c.0 == r) {
                    total += 1;
                    overlap += usize::from(v.iter().any(|c| c.0 != r && active.contains(&c.0)));
                }
            }
            (overlap, total, std::cmp::Reverse(r))
        };
        if let Some(drop) = active.iter().copied().min_by_key(|&r| score(r)) {
            active.remove(&drop);
        }
    };
    let mut factor: Vec<Option<f64>> = vec![None; n_runs];
    if !shared.is_empty() {
        let mut dev: Vec<Vec<f64>> = vec![Vec::new(); n_runs];
        for v in &shared {
            let cells: Vec<&(usize, f64)> = v.iter().filter(|c| active.contains(&c.0)).collect();
            let mean = cells.iter().map(|c| c.1).sum::<f64>() / cells.len() as f64;
            for c in cells {
                dev[c.0].push(c.1 - mean);
            }
        }
        for &r in &active {
            factor[r] = model::stats::median_mut(&mut dev[r], model::stats::Even::Mean).map(f64::exp2);
        }
    }
    for (r, s) in stems.iter().enumerate() {
        match factor[r] {
            Some(f) => log::info!(
                target: "dia",
                "quantities, run {s}: median-level factor {f:.4} over {} precursors quantified \
                 in every normalised run",
                shared.len()
            ),
            None => log::info!(
                target: "dia",
                "quantities, run {s}: left unscaled (fewer than {MIN_SHARED_PRECURSORS} \
                 quantified precursors shared with the other runs)"
            ),
        }
    }
    for r in rows.iter_mut() {
        let f = run_of.get(r.filename.as_str()).and_then(|&i| factor[i]).unwrap_or(1.0);
        r.quantity_normalized = r.quantity.filter(|v| *v > 0.0 && v.is_finite()).map(|v| v / f);
    }
    factor
}

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
            match memory.sizing {
                output::Sizing::Fixed => memory.bytes.to_string(),
                output::Sizing::FitInput { .. } => "auto".to_string(),
            },
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
        ["native", "auto", spill, a] => Ok((
            PercolatorBackend::Native {
                memory: crate::memory::rescore_memory(None, PathBuf::from(spill)),
            },
            args(a),
        )),
        ["native", bytes, spill, a] => Ok((
            PercolatorBackend::Native {
                memory: output::MemLimit {
                    bytes: bytes.parse().map_err(|_| bad())?,
                    spill_dir: PathBuf::from(spill),
                    sizing: output::Sizing::Fixed,
                },
            },
            args(a),
        )),
        _ => Err(bad()),
    }
}

/// The scoring part of stored rescoring settings: the engine and its arguments. The native
/// engine's storage budget and spill directory only bound memory, so runs searched with
/// per-job budgets or temporary directories still compare equal.
fn scoring_settings(stored: &str) -> Result<String, String> {
    let (backend, args) = decode_rescore(stored)?;
    let engine = match &backend {
        PercolatorBackend::Binary(p) => format!("binary {}", p.display()),
        PercolatorBackend::Docker { image } => format!("docker {image}"),
        PercolatorBackend::Native { .. } => "native".to_string(),
    };
    Ok(format!("{engine}; {}", args.join(" ")))
}

/// The study's rescoring backend and arguments from the settings stored by each run (`(stem,
/// stored)`, study order, at least one): an error when two runs differ in their scoring
/// settings. The study is always rescored by the native engine, with the first run's fixed
/// storage budget, else a budget sized to the study's PIN within the memory free now, spilling
/// to `scratch`.
fn study_rescore(
    stored: &[(&str, String)],
    scratch: &Path,
) -> Result<(PercolatorBackend, Vec<String>), String> {
    let (s0, r0) = &stored[0];
    let k0 = scoring_settings(r0)?;
    for (s, r) in &stored[1..] {
        if scoring_settings(r)? != k0 {
            return Err(format!(
                "runs {s0} and {s} were searched with different Percolator settings \
                 ({r0:?} and {r:?}); search every run of a study with the same ones"
            ));
        }
    }
    let (searched, searched_args) = decode_rescore(r0)?;
    Ok(match searched {
        PercolatorBackend::Native { memory } => (
            PercolatorBackend::Native {
                memory: match memory.sizing {
                    output::Sizing::Fixed => output::MemLimit {
                        bytes: memory.bytes,
                        spill_dir: scratch.to_path_buf(),
                        sizing: output::Sizing::Fixed,
                    },
                    output::Sizing::FitInput { .. } => {
                        crate::memory::rescore_memory(None, scratch.to_path_buf())
                    }
                },
            },
            searched_args,
        ),
        other => {
            log::info!(
                target: "dia",
                "the runs were rescored with Percolator {}; the study is rescored with the \
                 native engine (storage budget sized to the study's PIN, spill to {})",
                other.describe(),
                scratch.display()
            );
            let memory = crate::memory::rescore_memory(None, scratch.to_path_buf());
            (PercolatorBackend::Native { memory }, Vec::new())
        }
    })
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
    let stem = crate::quant::run::run_stem(spectrum);
    let pin_out = dir.join(format!("{stem}.pin"));
    let partial = dir.join(format!("{stem}.pin.partial"));
    std::fs::copy(pin, &partial)
        .map_err(|e| format!("copy {} to {}: {e}", pin.display(), partial.display()))?;
    std::fs::rename(&partial, &pin_out)
        .map_err(|e| format!("rename {}: {e}", partial.display()))?;
    let meta: Vec<(&str, String)> = rescore
        .map(|(b, a)| (RESCORE_KEY, encode_rescore(b, a)))
        .into_iter()
        .chain([(crate::dia::transfer::STAMP_KEY, crate::dia::transfer::pass2_stamp(table))])
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

/// Whether `stem` is a DIA run stored in the quant directory (its pass-2 table exists).
pub(crate) fn is_study_run(quant_dir: &Path, stem: &str) -> bool {
    quant_dir
        .join(PASS2_DIR)
        .join(format!("{stem}.parquet"))
        .is_file()
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
    let with_xfer: Vec<&String> = stems
        .iter()
        .filter(|s| crate::dia::transfer::xfer_path(&a.quant_dir, s).is_file())
        .collect();
    let transfer = !with_xfer.is_empty();
    if transfer && with_xfer.len() != stems.len() {
        let missing: Vec<&str> = stems
            .iter()
            .filter(|s| !with_xfer.contains(s))
            .map(String::as_str)
            .collect();
        return Err(format!(
            "{}: {} of {} runs have a candidate-transfer table; run `andes quant extract \
             --quant-dir {} --run <stem>` for {}",
            dir.display(),
            with_xfer.len(),
            stems.len(),
            a.quant_dir.display(),
            missing.join(", ")
        ));
    }
    let table_path = |s: &str| {
        if transfer {
            crate::dia::transfer::xfer_path(&a.quant_dir, s)
        } else {
            dir.join(format!("{s}.parquet"))
        }
    };
    if transfer {
        check_transfer_tables(a, stems, &table_path)?;
    }
    let stored = rescore_settings(stems, &table_path)?;
    if stored.is_empty() {
        return Err(format!("{}: no DIA runs", dir.display()));
    }
    let scratch = a
        .temp_dir
        .clone()
        .or_else(|| {
            a.output
                .parent()
                .filter(|d| !d.as_os_str().is_empty())
                .map(Path::to_path_buf)
        })
        .unwrap_or_else(|| PathBuf::from("."));
    let (backend, mut args) = study_rescore(&stored, &scratch)?;
    drop(stored);
    let seed = a.seed;
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
    let mut quant_rows: Vec<DiaQuantRow> = Vec::new();
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
        crate::dia::score::write_precursor_pin(&pin_in, &pin, &t, &xr, seed)?;
        let has_xfer = t.feature(crate::dia::pass2::XFER_COLUMN).is_some();
        let cands = quant_candidates(&t);
        drop(t);
        drop(xr);
        let (map, decoys) = rescore_run(s, &pin, &backend, &extra, has_xfer, t_run)?;
        let rows = reported_rows(s, &cands, &map, a.quant_fdr);
        drop(cands);
        log::info!(
            target: "dia",
            "run {s}: {} target precursors at q <= {}, {} with a fragment and {} with an MS1 \
             quantity",
            rows.len(),
            a.quant_fdr,
            rows.iter().filter(|r| r.quantity_frag.is_some()).count(),
            rows.iter().filter(|r| r.quantity_ms1.is_some()).count()
        );
        quant_rows.extend(rows);
        if a.database.is_some() {
            all_targets.extend(map.into_iter().map(|(id, p)| (format!("{s}\t{id}"), p)));
            all_decoys.extend(decoys.into_iter().map(|(id, p)| (format!("{s}\t{id}"), p)));
        }
    }
    drop(keys);
    normalise(&mut quant_rows, stems);
    let quant_path = with_suffix(&a.output, ".dia_quant.parquet");
    write_dia_quant(
        &quant_path,
        &quant_rows,
        &[(DEFAULT_ARM_KEY, DEFAULT_QUANTITY.column().to_string())],
    )?;
    log::info!(
        target: "dia",
        "wrote {} ({} run-precursor rows; quantity = {}, arms {})",
        quant_path.display(),
        quant_rows.len(),
        DEFAULT_QUANTITY.column(),
        QUANT_ARMS.map(QuantArm::column).join(", ")
    );
    drop(quant_rows);
    if let Some(database) = &a.database {
        write_protein_groups(a, database, stems, all_targets, &all_decoys)?;
    }
    log::info!(
        target: "total",
        "DIA study finalized in {:.2} s",
        t0.elapsed().as_secs_f64()
    );
    Ok(())
}

type PsmMap = std::collections::HashMap<String, output::PercolatorPsm>;

/// With cross-run transfer: every run's pass-2 table must be built from the study's runs and
/// from its current inputs.
fn check_transfer_tables(
    a: &crate::cli::QuantFinalizeArgs,
    stems: &[String],
    table_path: &dyn Fn(&str) -> PathBuf,
) -> Result<(), String> {
    let want = stems.join("\n");
    let current: Vec<String> = stems
        .iter()
        .map(|s| crate::dia::transfer::run_inputs(&a.quant_dir, s))
        .collect::<Result<_, _>>()?;
    for s in stems {
        let p = table_path(s);
        let built = read_pass2_meta(&p)?
            .into_iter()
            .find(|(k, _)| k == crate::dia::transfer::RUNS_KEY)
            .map(|(_, v)| v);
        if built.as_deref() != Some(want.as_str()) {
            return Err(format!(
                "{}: built from other runs than the study's {} (stale transfer); run \
                 `andes quant extract --quant-dir {} --run {s}` again",
                p.display(),
                stems.join(", "),
                a.quant_dir.display()
            ));
        }
        let built = read_pass2_meta(&p)?
            .into_iter()
            .find(|(k, _)| k == crate::dia::transfer::INPUTS_KEY)
            .map(|(_, v)| v);
        crate::dia::transfer::check_inputs(&a.quant_dir, &p, s, built.as_deref(), &current)?;
    }
    log::info!(
        target: "dia",
        "DIA study with cross-run candidate transfer: every run's pass-2 table holds its own \
         and its transferred nominees (column {})",
        crate::dia::pass2::XFER_COLUMN
    );
    Ok(())
}

/// Each run's stored rescoring settings, from its pass-2 table.
fn rescore_settings<'s>(
    stems: &'s [String],
    table_path: &dyn Fn(&str) -> PathBuf,
) -> Result<Vec<(&'s str, String)>, String> {
    let mut stored: Vec<(&str, String)> = Vec::with_capacity(stems.len());
    for s in stems {
        let p = table_path(s);
        let meta = read_pass2_meta(&p)?;
        let r = meta
            .into_iter()
            .find(|(k, _)| k == RESCORE_KEY)
            .map(|(_, v)| v)
            .ok_or_else(|| {
                format!(
                    "{}: run {s} was searched without --rescore; a DIA study \
                     is rescored with the settings of its runs",
                    p.display()
                )
            })?;
        stored.push((s.as_str(), r));
    }
    Ok(stored)
}

/// Percolator on run `s`'s precursor PIN, then precursor-level q-values (per transfer class
/// too when `has_xfer`). Returns the target and decoy results.
fn rescore_run(
    s: &str,
    pin: &Path,
    backend: &PercolatorBackend,
    extra: &[String],
    has_xfer: bool,
    t_run: std::time::Instant,
) -> Result<(PsmMap, PsmMap), String> {
    let (mut map, _) = output::run_percolator_with_decoys(backend, pin, extra)
        .map_err(|e| format!("Percolator on {}: {e}", pin.display()))?;
    let dpath = output::decoy_results_path(pin);
    let text = std::fs::read_to_string(&dpath).map_err(|e| format!("{}: {e}", dpath.display()))?;
    let decoys =
        output::parse_psm_results(&text).map_err(|e| format!("{}: {e}", dpath.display()))?;
    drop(text);
    let ((nt, nd, n01), _) = crate::dia::precursor_q(pin, &[], &mut map, &decoys)?;
    log::info!(
        target: "rescore",
        "run {s}: precursor-level q (DIA): {nt} target and {nd} decoy precursors, {n01} at \
         q <= 0.01 [{:.2}s]",
        t_run.elapsed().as_secs_f64()
    );
    if has_xfer {
        let c = crate::dia::precursor_q_by_class(
            pin,
            crate::dia::pass2::XFER_COLUMN,
            &mut map,
            &decoys,
        )?;
        log::info!(
            target: "rescore",
            "run {s}: target precursors at q <= 0.01 own / transferred: class-wise q {} / {} \
             (reported), one pooled competition would give {} / {}",
            c[1][0],
            c[1][1],
            c[0][0],
            c[0][1]
        );
    }
    Ok((map, decoys))
}

/// Protein inference over every run's PSMs; writes `<base>.pg.parquet`.
fn write_protein_groups(
    a: &crate::cli::QuantFinalizeArgs,
    database: &Path,
    stems: &[String],
    all_targets: PsmMap,
    all_decoys: &PsmMap,
) -> Result<(), String> {
    let inf = crate::inference::infer_from_maps(database, &a.proteins, all_targets, all_decoys)
        .map_err(|e| e.to_string())?;
    let path = with_suffix(&a.output, ".pg.parquet");
    output::write_pg_parquet(
        &path,
        &inf,
        stems,
        stems.len() == 1,
        a.proteins.protein_fdr,
        None,
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
            ..Default::default()
        }
    }

    fn psm(id: &str, q: f64) -> (String, output::PercolatorPsm) {
        (
            id.to_string(),
            output::PercolatorPsm {
                psm_id: id.to_string(),
                q_value: q,
                pooled_q_value: None,
                pep: 0.0,
                score: Some(1.0),
                peptide: String::new(),
                proteins: String::new(),
            },
        )
    }

    fn cand(id: &str, pep: &str, transfer: bool, frag: Option<f64>) -> QuantCandidate {
        QuantCandidate {
            spec_id: id.into(),
            peptide: pep.into(),
            charge: 2,
            proteins: "P1".into(),
            transfer,
            frag,
            ms1: frag.map(|v| 10.0 * v),
        }
    }

    #[test]
    fn only_rows_at_reported_q_are_quantified() {
        let mut t = table(&[
            ("AAK", 2, false, 1.0, 1.0),
            ("CCK", 2, true, 1.0, 1.0),
            ("DDK", 2, false, 1.0, 1.0),
        ]);
        t.quant_frag = vec![Some(5.0), Some(6.0), None];
        t.quant_ms1 = vec![Some(50.0), Some(60.0), Some(70.0)];
        t.quant_n_frag = vec![6, 6, 0];
        let c = quant_candidates(&t);
        assert_eq!(c.len(), 2, "decoys are never quantified");
        let cands = vec![
            cand("a", "AAK", false, Some(5.0)),
            cand("b", "BBK", false, Some(6.0)),
            cand("c", "CCK", true, Some(7.0)),
            cand("d", "DDK", false, None),
            cand("e", "EEK", false, Some(9.0)),
            cand("f", "AAK", true, Some(1.0)),
        ];
        let q: std::collections::HashMap<_, _> = [
            psm("a", 0.002),
            psm("b", 0.02),
            psm("c", 0.01),
            psm("d", 0.004),
            psm("f", 0.002),
        ]
        .into_iter()
        .collect();
        let rows = reported_rows("run1", &cands, &q, 0.01);
        let got: Vec<(&str, bool, Option<f64>)> = rows
            .iter()
            .map(|r| (r.peptide.as_str(), r.transfer, r.quantity))
            .collect();
        assert_eq!(
            got,
            [("AAK", false, Some(5.0)), ("CCK", true, Some(7.0)), ("DDK", false, None)],
            "above q (b) and without a q (e) are out; at equal q the own row wins"
        );
        assert!(rows.iter().all(|r| r.filename == "run1" && r.q_value <= 0.01));
        assert_eq!(rows[1].quantity_ms1, Some(70.0));
    }

    fn qrow(run: &str, pep: &str, v: f64) -> DiaQuantRow {
        DiaQuantRow {
            filename: run.into(),
            peptide: pep.into(),
            charge: 2,
            proteins: "P".into(),
            transfer: false,
            q_value: 0.001,
            quantity: Some(v),
            quantity_frag: Some(v),
            quantity_ms1: None,
            quantity_normalized: None,
        }
    }

    #[test]
    fn median_normalisation_over_shared_precursors() {
        let stems: Vec<String> = ["A", "B", "C"].iter().map(|s| s.to_string()).collect();
        let mut rows = Vec::new();
        for p in 0..25 {
            let base = 1000.0 * (1 + p) as f64;
            rows.push(qrow("A", &format!("P{p}"), base));
            rows.push(qrow("B", &format!("P{p}"), 2.0 * base));
            rows.push(qrow("C", &format!("P{p}"), 0.5 * base));
        }
        for p in 0..10 {
            rows.push(qrow("B", &format!("X{p}"), 1e9));
        }
        let mut missing = qrow("C", "M0", 1.0);
        missing.quantity = None;
        rows.push(missing);
        let f = normalise(&mut rows, &stems);
        for (got, want) in f.iter().zip([1.0, 2.0, 0.5]) {
            assert!((got.unwrap() - want).abs() < 1e-9, "{got:?} vs {want}");
        }
        for p in 0..25 {
            let n: Vec<f64> = rows
                .iter()
                .filter(|r| r.peptide == format!("P{p}"))
                .map(|r| r.quantity_normalized.unwrap())
                .collect();
            assert!(n.iter().all(|v| (v - n[0]).abs() < 1e-6 * n[0]), "{n:?}");
        }
        assert!(rows
            .iter()
            .filter(|r| r.peptide.starts_with('X'))
            .all(|r| r.quantity_normalized == Some(5e8)));
        assert_eq!(rows.last().unwrap().quantity_normalized, None, "missing stays missing");
        assert!(rows.iter().flat_map(|r| r.quantity_normalized).all(|v| v > 0.0 && v.is_finite()));
    }

    #[test]
    fn run_without_shared_precursors_is_left_unscaled() {
        let stems: Vec<String> = ["A", "B", "C"].iter().map(|s| s.to_string()).collect();
        let mut rows = Vec::new();
        for p in 0..30 {
            let base = 100.0 * (1 + p) as f64;
            rows.push(qrow("A", &format!("P{p}"), base));
            rows.push(qrow("B", &format!("P{p}"), 4.0 * base));
            rows.push(qrow("C", &format!("Q{p}"), 3.0 * base));
        }
        let f = normalise(&mut rows, &stems);
        assert!((f[0].unwrap() - 0.5).abs() < 1e-9);
        assert!((f[1].unwrap() - 2.0).abs() < 1e-9);
        assert_eq!(f[2], None);
        for r in &rows {
            if r.filename == "C" {
                assert_eq!(r.quantity_normalized, r.quantity);
            }
        }
        let mut alone = vec![qrow("A", "P0", 7.0)];
        assert_eq!(normalise(&mut alone, &stems), vec![None; 3]);
        assert_eq!(alone[0].quantity_normalized, Some(7.0));
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
    fn runs_differing_only_in_the_spill_dir_and_budget_share_their_settings() {
        let native = |gb: u64, dir: &str, args: &[&str]| {
            let args: Vec<String> = args.iter().map(|a| a.to_string()).collect();
            encode_rescore(
                &PercolatorBackend::Native {
                    memory: output::MemLimit {
                        bytes: gb << 30,
                        spill_dir: PathBuf::from(dir),
                        sizing: output::Sizing::Fixed,
                    },
                },
                &args,
            )
        };
        let stored = vec![
            ("run01", native(8, "/scratch/job1/run01", &["--trainFDR", "0.05"])),
            ("run02", native(16, "/scratch/job2/run02", &["--trainFDR", "0.05"])),
        ];
        let (b, args) = study_rescore(&stored, Path::new("/finalize/tmp")).unwrap();
        assert_eq!(args, ["--trainFDR", "0.05"]);
        match b {
            PercolatorBackend::Native { memory } => {
                assert_eq!(memory.spill_dir, PathBuf::from("/finalize/tmp"));
                assert_eq!(memory.bytes, 8 << 30);
            }
            other => panic!("expected the native engine, got {}", other.describe()),
        }
        let differ = vec![
            ("run01", native(8, "/a", &["--trainFDR", "0.05"])),
            ("run02", native(8, "/a", &["--trainFDR", "0.01"])),
        ];
        let e = study_rescore(&differ, Path::new(".")).unwrap_err();
        assert!(e.contains("run01") && e.contains("run02"), "{e}");
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
                sizing: output::Sizing::Fixed,
            },
        };
        let (n2, _) = decode_rescore(&encode_rescore(&n, &[])).unwrap();
        assert_eq!(n2.describe(), n.describe());
        let auto = PercolatorBackend::Native {
            memory: crate::memory::rescore_memory(None, PathBuf::from("/tmp/x")),
        };
        let stored = encode_rescore(&auto, &args);
        assert!(stored.starts_with("native\nauto\n/tmp/x\n"), "{stored:?}");
        match decode_rescore(&stored).unwrap() {
            (PercolatorBackend::Native { memory }, a) => {
                assert!(matches!(memory.sizing, output::Sizing::FitInput { .. }));
                assert_eq!(memory.spill_dir, PathBuf::from("/tmp/x"));
                assert_eq!(a, args);
            }
            (other, _) => panic!("expected the native engine, got {}", other.describe()),
        }
        let (study, _) = study_rescore(&[("run01", stored)], Path::new("/fin")).unwrap();
        match study {
            PercolatorBackend::Native { memory } => {
                assert!(matches!(memory.sizing, output::Sizing::FitInput { .. }));
                assert_eq!(memory.spill_dir, PathBuf::from("/fin"));
            }
            other => panic!("expected the native engine, got {}", other.describe()),
        }
        assert!(decode_rescore("binary").is_err());
    }
}
