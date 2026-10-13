//! PeptDeep MS2 + RT predictions for every row of a Percolator PIN, with tract.
//!
//! Reads an andes PIN and the mzML it was searched from, predicts fragment
//! intensities and retention time for each distinct peptidoform with the
//! exported PeptDeep ONNX graphs, and writes one feature row per PIN row
//! (targets and decoys alike) to a TSV keyed by `SpecId`. RT calibration and
//! PIN assembly happen downstream.
//!
//! ```text
//! cargo run --release -p onnx-infer --features tract --example rescore_pin -- \
//!   --pin search.pin --mzml run.mzML --models DIR --out features.tsv \
//!   --frag-tol 0.5 --frag-unit da --threads 8
//! ```
//!
//! `DIR` holds `peptdeep_rt_dynamic.json`, `peptdeep_ms2_dynamic.json` and
//! `peptdeep_mods.json` (written by `scripts/peptdeep_reference.py`).
//!
//! `--parity peptdeep_reference.tsv` instead compares tract predictions with
//! peptdeep's own `predict` output and exits.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

use onnx_infer::andes_peptide::split_nterm_residue_delta;
use onnx_infer::peptdeep::{PeptDeepPeptide, pack_ms2_batch, pack_rt_batch};
use onnx_infer::{BackendKind, OnnxSession};
use serde::Deserialize;

const PROTON: f64 = 1.007_276_466_8;
const H2O: f64 = 18.010_564_684;

fn residue_mass(c: char) -> Option<f64> {
    Some(match c {
        'G' => 57.021_463_72,
        'A' => 71.037_113_79,
        'S' => 87.032_028_41,
        'P' => 97.052_763_85,
        'V' => 99.068_413_91,
        'T' => 101.047_678_5,
        'C' => 103.009_184_5,
        'L' | 'I' => 113.084_064,
        'N' => 114.042_927_4,
        'D' => 115.026_943,
        'Q' => 128.058_577_5,
        'K' => 128.094_963,
        'E' => 129.042_593_1,
        'M' => 131.040_484_6,
        'H' => 137.058_911_9,
        'F' => 147.068_413_9,
        'R' => 156.101_111,
        'Y' => 163.063_328_6,
        'W' => 186.079_312_9,
        _ => return None,
    })
}

/// andes mass delta → (alphabase mod name, applies at the peptide N-terminus).
/// N-terminal deltas are printed by andes after the first residue.
const MOD_TABLE: &[(Option<char>, f64, &str, bool)] = &[
    (Some('C'), 57.02146, "Carbamidomethyl@C", false),
    (Some('M'), 15.99491, "Oxidation@M", false),
    (None, 42.01056, "Acetyl@Protein_N-term", true),
];

#[derive(Deserialize)]
struct ModsMeta {
    mod_dim: usize,
    mods: HashMap<String, Vec<f32>>,
    ms2_frag_types: Vec<String>,
    instrument_index: HashMap<String, i64>,
    unknown_instrument_index: i64,
}

/// One parsed peptidoform.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct Form {
    seq: String,
    /// (site, mod name) — site 0 = N-term, 1..=n residue.
    mods: Vec<(usize, &'static str)>,
}

struct Mass {
    /// Prefix residue+mod masses, len n.
    res: Vec<f64>,
    nterm: f64,
}

/// Parse `K.PEPM+15.99491K.D` into a peptidoform; `None` for unsupported text.
fn parse_peptide(s: &str) -> Option<(Form, Mass)> {
    let b = s.as_bytes();
    let core = if b.len() > 4 && b[1] == b'.' && b[b.len() - 2] == b'.' {
        &s[2..s.len() - 2]
    } else {
        s
    };
    let mut seq = String::new();
    let mut res = Vec::new();
    let mut mods = Vec::new();
    let mut nterm = 0.0;
    let chars: Vec<char> = core.chars().collect();
    let mut i = 0;
    while i < chars.len() {
        let c = chars[i];
        if c.is_ascii_uppercase() {
            seq.push(c);
            res.push(residue_mass(c)?);
            i += 1;
        } else if c == '+' || c == '-' {
            let start = i;
            i += 1;
            while i < chars.len() && (chars[i].is_ascii_digit() || chars[i] == '.') {
                i += 1;
            }
            let delta: f64 = chars[start..i].iter().collect::<String>().parse().ok()?;
            let pos = seq.len(); // 1-based residue the delta follows
            if pos == 0 {
                return None;
            }
            let aa = seq.as_bytes()[pos - 1] as char;
            let single = MOD_TABLE.iter().find_map(|(r, m, name, nt)| {
                let ok = (delta - m).abs() < 2e-3 && (*nt && pos == 1 || r.is_some_and(|r| r == aa));
                ok.then_some((*name, *nt))
            });
            let Some((name, nt)) = single else {
                // An N-terminal and a residue mod on residue 1 print as one summed delta.
                if pos != 1 || !mods.is_empty() {
                    return None;
                }
                let nterm_cands: Vec<(&'static str, f64)> = MOD_TABLE
                    .iter()
                    .filter(|e| e.3)
                    .map(|e| (e.2, e.1))
                    .collect();
                let residue_cands: Vec<(&'static str, f64)> = MOD_TABLE
                    .iter()
                    .filter(|e| !e.3 && e.0 == Some(aa))
                    .map(|e| (e.2, e.1))
                    .collect();
                let (n, r) = split_nterm_residue_delta(delta, 2e-3, &nterm_cands, &residue_cands)?;
                let n_mass = nterm_cands.iter().find(|e| e.0 == n).map(|e| e.1)?;
                mods.push((0, n));
                mods.push((1, r));
                nterm += n_mass;
                res[0] += delta - n_mass;
                continue;
            };
            if nt && pos == 1 && mods.iter().all(|(s, _)| *s != 0) && name.starts_with("Acetyl") {
                mods.push((0, name));
                nterm += delta;
            } else if !nt {
                mods.push((pos, name));
                res[pos - 1] += delta;
            } else {
                return None;
            }
        } else {
            return None;
        }
    }
    if seq.is_empty() {
        return None;
    }
    mods.sort();
    Some((Form { seq, mods }, Mass { res, nterm }))
}

struct Args {
    pin: PathBuf,
    mzml: PathBuf,
    models: PathBuf,
    out: PathBuf,
    frag_tol: f64,
    frag_ppm: bool,
    threads: usize,
    batch: usize,
    nce: f32,
    instrument: String,
    parity: Option<PathBuf>,
}

fn parse_args() -> Args {
    let mut a = Args {
        pin: PathBuf::new(),
        mzml: PathBuf::new(),
        models: PathBuf::new(),
        out: PathBuf::new(),
        frag_tol: 0.5,
        frag_ppm: false,
        threads: 8,
        batch: 256,
        nce: 30.0,
        instrument: "Lumos".into(),
        parity: None,
    };
    let v: Vec<String> = std::env::args().skip(1).collect();
    let mut i = 0;
    while i < v.len() {
        let val = v.get(i + 1).cloned().unwrap_or_default();
        match v[i].as_str() {
            "--pin" => a.pin = val.into(),
            "--mzml" => a.mzml = val.into(),
            "--models" => a.models = val.into(),
            "--out" => a.out = val.into(),
            "--frag-tol" => a.frag_tol = val.parse().expect("--frag-tol"),
            "--frag-unit" => a.frag_ppm = val == "ppm",
            "--threads" => a.threads = val.parse().expect("--threads"),
            "--batch" => a.batch = val.parse().expect("--batch"),
            "--nce" => a.nce = val.parse().expect("--nce"),
            "--instrument" => a.instrument = val,
            "--parity" => a.parity = Some(val.into()),
            other => panic!("unknown argument {other}"),
        }
        i += 2;
    }
    a
}

fn load(models: &Path, name: &str) -> OnnxSession {
    OnnxSession::from_manifest(models.join(name), BackendKind::Tract)
        .unwrap_or_else(|e| panic!("load {name}: {e}"))
}

struct Predictor<'m> {
    meta: &'m ModsMeta,
    nce: f32,
    instrument: i64,
    /// Column of b_z1, b_z2, y_z1, y_z2 in the MS2 output.
    cols: [usize; 4],
}

impl Predictor<'_> {
    fn peptides<'f>(&'f self, forms: &[&'f Form]) -> Vec<PeptDeepPeptide<'f>> {
        forms
            .iter()
            .map(|f| PeptDeepPeptide {
                sequence: &f.seq,
                mods: f
                    .mods
                    .iter()
                    .map(|(s, n)| (*s, self.meta.mods[*n].as_slice()))
                    .collect(),
            })
            .collect()
    }

    fn rt(&self, s: &mut OnnxSession, forms: &[&Form]) -> Vec<f32> {
        let len = forms[0].seq.len();
        let batch = pack_rt_batch(&self.peptides(forms), len, self.meta.mod_dim).unwrap();
        let out = s.run(&batch.as_inputs()).expect("rt run");
        out[0].data.clone()
    }

    /// Per form: (len-1) × [b1, b2, y1, y2] intensities, clipped at 0.
    fn ms2(&self, s: &mut OnnxSession, forms: &[(&Form, u8)]) -> Vec<Vec<f32>> {
        let len = forms[0].0.seq.len();
        let only: Vec<&Form> = forms.iter().map(|(f, _)| *f).collect();
        let charges: Vec<u8> = forms.iter().map(|(_, z)| *z).collect();
        let batch = pack_ms2_batch(
            &self.peptides(&only),
            &charges,
            len,
            self.meta.mod_dim,
            self.nce,
            self.instrument,
        )
        .unwrap();
        let out = &s.run(&batch.as_inputs()).expect("ms2 run")[0];
        let (rows, ncol) = (out.shape[1], out.shape[2]);
        assert_eq!(rows, len - 1, "MS2 output rows");
        (0..forms.len())
            .map(|b| {
                let mut v = Vec::with_capacity(rows * 4);
                for r in 0..rows {
                    for (k, &c) in self.cols.iter().enumerate() {
                        let frag_z = 1 + (k % 2) as u8;
                        let x = out.data[(b * rows + r) * ncol + c].max(0.0);
                        v.push(if frag_z > charges[b] { 0.0 } else { x });
                    }
                }
                v
            })
            .collect()
    }
}

fn pearson(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len() as f64;
    if a.len() < 2 {
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

const MS2_FEATURES: &[&str] = &[
    "spec_pearson_norm",
    "ionb_pearson_norm",
    "iony_pearson_norm",
    "spec_mse_norm",
    "dotprod_norm",
    "cos_norm",
    "spectral_angle",
    "pred_explained",
];

/// Spectral-match features between predicted and observed b/y intensities.
fn ms2_features(
    pred: &[f32],
    mass: &Mass,
    charge: u8,
    peaks: &[(f64, f32)],
    tol: f64,
    ppm: bool,
) -> [f64; 8] {
    let n = mass.res.len();
    let total: f64 = mass.res.iter().sum::<f64>() + mass.nterm;
    let mut p = Vec::new();
    let mut o = Vec::new();
    let mut is_b = Vec::new();
    let mut prefix = mass.nterm;
    for i in 0..n - 1 {
        prefix += mass.res[i];
        let b = prefix;
        let y = total - prefix + H2O;
        for (k, (neutral, bion)) in [(b, true), (b, true), (y, false), (y, false)].iter().enumerate() {
            let z = 1 + (k % 2) as u8;
            if z > charge {
                continue;
            }
            let mz = (neutral + z as f64 * PROTON) / z as f64;
            let w = if ppm { mz * tol * 1e-6 } else { tol };
            let lo = peaks.partition_point(|(m, _)| *m < mz - w);
            let obs = peaks[lo..]
                .iter()
                .take_while(|(m, _)| *m <= mz + w)
                .map(|(_, i)| *i as f64)
                .fold(0.0, f64::max);
            p.push(pred[i * 4 + k] as f64);
            o.push(obs);
            is_b.push(*bion);
        }
    }
    let (sp, so) = (p.iter().sum::<f64>(), o.iter().sum::<f64>());
    let pn: Vec<f64> = p.iter().map(|v| if sp > 0.0 { v / sp } else { 0.0 }).collect();
    let on: Vec<f64> = o.iter().map(|v| if so > 0.0 { v / so } else { 0.0 }).collect();
    let sel = |want_b: bool, v: &[f64]| -> Vec<f64> {
        v.iter().zip(&is_b).filter(|(_, b)| **b == want_b).map(|(x, _)| *x).collect()
    };
    let dot: f64 = pn.iter().zip(&on).map(|(a, b)| a * b).sum();
    let (np, no) = (
        pn.iter().map(|v| v * v).sum::<f64>().sqrt(),
        on.iter().map(|v| v * v).sum::<f64>().sqrt(),
    );
    let cos = if np > 0.0 && no > 0.0 { (dot / (np * no)).clamp(-1.0, 1.0) } else { 0.0 };
    let sa = if np > 0.0 && no > 0.0 { 1.0 - 2.0 * cos.acos() / std::f64::consts::PI } else { 0.0 };
    let mse = pn.iter().zip(&on).map(|(a, b)| (a - b) * (a - b)).sum::<f64>() / pn.len().max(1) as f64;
    let explained: f64 = pn.iter().zip(&o).filter(|(_, ob)| **ob > 0.0).map(|(v, _)| *v).sum();
    [
        pearson(&pn, &on),
        pearson(&sel(true, &pn), &sel(true, &on)),
        pearson(&sel(false, &pn), &sel(false, &on)),
        mse,
        dot,
        cos,
        sa,
        explained,
    ]
}

fn rss_mb() -> f64 {
    std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|s| {
            s.lines()
                .find(|l| l.starts_with("VmHWM:"))
                .and_then(|l| l.split_whitespace().nth(1)?.parse::<f64>().ok())
        })
        .map_or(0.0, |kb| kb / 1024.0)
}

/// Run `jobs` across `threads` workers, each with its own RT + MS2 sessions.
fn run_parallel<J: Sync, F>(args: &Args, jobs: &[J], work: F) -> (Vec<f64>, f64)
where
    F: Fn(&mut OnnxSession, &mut OnnxSession, &J) + Sync,
{
    let next = AtomicUsize::new(0);
    let loads = Mutex::new(Vec::new());
    let t0 = Instant::now();
    std::thread::scope(|sc| {
        for _ in 0..args.threads {
            sc.spawn(|| {
                let tl = Instant::now();
                let mut rt = load(&args.models, "peptdeep_rt_dynamic.json");
                let mut ms2 = load(&args.models, "peptdeep_ms2_dynamic.json");
                loads.lock().unwrap().push(tl.elapsed().as_secs_f64());
                loop {
                    let j = next.fetch_add(1, Ordering::Relaxed);
                    if j >= jobs.len() {
                        break;
                    }
                    work(&mut rt, &mut ms2, &jobs[j]);
                }
            });
        }
    });
    (loads.into_inner().unwrap(), t0.elapsed().as_secs_f64())
}

fn parity(args: &Args, meta: &ModsMeta, pred: &Predictor) {
    let mut rt = load(&args.models, "peptdeep_rt_dynamic.json");
    let mut ms2 = load(&args.models, "peptdeep_ms2_dynamic.json");
    let f = BufReader::new(File::open(args.parity.as_ref().unwrap()).unwrap());
    let mut lines = f.lines();
    let header: Vec<String> = lines.next().unwrap().unwrap().split('\t').map(String::from).collect();
    let col = |n: &str| header.iter().position(|h| h == n).unwrap();
    let (mut max_rt, mut max_ms2) = (0f64, 0f64);
    for line in lines {
        let line = line.unwrap();
        let x: Vec<&str> = line.split('\t').collect();
        let seq = x[col("sequence")];
        let names: Vec<&str> = x[col("mods")].split(';').filter(|s| !s.is_empty()).collect();
        let sites: Vec<usize> = x[col("mod_sites")].split(';').filter(|s| !s.is_empty()).map(|s| s.parse().unwrap()).collect();
        let charge: u8 = x[col("charge")].parse().unwrap();
        let mods: Vec<(usize, &'static str)> = names
            .iter()
            .zip(&sites)
            .map(|(n, s)| (*s, MOD_TABLE.iter().find(|m| m.2 == *n).expect("mod in table").2))
            .collect();
        let form = Form { seq: seq.to_string(), mods };
        let r = pred.rt(&mut rt, &[&form])[0] as f64;
        let ref_rt: f64 = x[col("rt_peptdeep")].parse().unwrap();
        max_rt = max_rt.max((r - ref_rt).abs());
        let m = &pred.ms2(&mut ms2, &[(&form, charge)])[0];
        let ref_ms2: Vec<f32> = x[col("ms2_peptdeep")].split(',').map(|v| v.parse().unwrap()).collect();
        // reference rows carry all frag-type columns; compare the four b/y columns
        let nft = meta.ms2_frag_types.len();
        for r in 0..seq.len() - 1 {
            for (k, &c) in pred.cols.iter().enumerate() {
                max_ms2 = max_ms2.max((m[r * 4 + k] - ref_ms2[r * nft + c]).abs() as f64);
            }
        }
        println!("{seq}\t{charge}\trt tract {r:.6} peptdeep {ref_rt:.6}");
    }
    println!("PARITY max|rt tract - peptdeep| = {max_rt:.3e}  max|ms2 tract - peptdeep| = {max_ms2:.3e}");
}

fn main() {
    let args = parse_args();
    let meta: ModsMeta =
        serde_json::from_reader(File::open(args.models.join("peptdeep_mods.json")).unwrap()).unwrap();
    let col = |n: &str| meta.ms2_frag_types.iter().position(|f| f == n).unwrap();
    let instrument = *meta
        .instrument_index
        .get(&args.instrument.to_uppercase())
        .unwrap_or(&meta.unknown_instrument_index);
    let pred = Predictor {
        meta: &meta,
        nce: args.nce,
        instrument,
        cols: [col("b_z1"), col("b_z2"), col("y_z1"), col("y_z2")],
    };
    if args.parity.is_some() {
        let t = Instant::now();
        load(&args.models, "peptdeep_rt_dynamic.json");
        let rt_load = t.elapsed().as_secs_f64();
        let t = Instant::now();
        load(&args.models, "peptdeep_ms2_dynamic.json");
        println!("load: rt {rt_load:.2}s ms2 {:.2}s", t.elapsed().as_secs_f64());
        parity(&args, &meta, &pred);
        return;
    }

    // PIN rows.
    let t_all = Instant::now();
    let f = BufReader::new(File::open(&args.pin).unwrap());
    let mut lines = f.lines();
    let header: Vec<String> = lines.next().unwrap().unwrap().split('\t').map(String::from).collect();
    let hc = |n: &str| header.iter().position(|h| h == n).unwrap_or_else(|| panic!("PIN column {n}"));
    let (c_id, c_scan, c_pep) = (hc("SpecId"), hc("ScanNr"), hc("Peptide"));
    let charge_cols: Vec<(usize, u8)> = header
        .iter()
        .enumerate()
        .filter_map(|(i, h)| Some((i, h.strip_prefix("charge")?.parse().ok()?)))
        .collect();
    struct Row {
        id: String,
        scan: i32,
        charge: u8,
        form: Option<usize>,
    }
    let mut forms: Vec<(Form, Mass)> = Vec::new();
    let mut form_ix: HashMap<Form, usize> = HashMap::new();
    let mut rows = Vec::new();
    for line in lines {
        let line = line.unwrap();
        let x: Vec<&str> = line.split('\t').collect();
        let charge = charge_cols.iter().find(|(i, _)| x[*i] == "1").map_or(0, |(_, z)| *z);
        let form = parse_peptide(x[c_pep]).map(|(f, m)| {
            *form_ix.entry(f.clone()).or_insert_with(|| {
                forms.push((f, m));
                forms.len() - 1
            })
        });
        rows.push(Row {
            id: x[c_id].to_string(),
            scan: x[c_scan].parse().unwrap(),
            charge,
            form,
        });
    }
    let unsupported = rows.iter().filter(|r| r.form.is_none() || r.charge == 0).count();
    eprintln!(
        "PIN rows {} distinct peptidoforms {} unsupported rows {} ({:.1}s)",
        rows.len(),
        forms.len(),
        unsupported,
        t_all.elapsed().as_secs_f64()
    );

    // Spectra.
    let t = Instant::now();
    let (spectra, errs) = input::read_spectra_parallel(&args.mzml, 2, 2, args.threads).unwrap();
    assert!(errs.is_empty(), "mzML errors: {:?}", &errs[..errs.len().min(3)]);
    let mut by_scan: HashMap<i32, usize> = HashMap::new();
    for (i, s) in spectra.iter().enumerate() {
        let scan = s.scan.expect("spectrum without scan number");
        assert!(by_scan.insert(scan, i).is_none(), "duplicate scan {scan}");
    }
    eprintln!("mzML {} MS2 spectra ({:.1}s)", spectra.len(), t.elapsed().as_secs_f64());

    // RT jobs: distinct forms grouped by length.
    let mut by_len: HashMap<usize, Vec<usize>> = HashMap::new();
    for (i, (f, _)) in forms.iter().enumerate() {
        by_len.entry(f.seq.len()).or_default().push(i);
    }
    let mut rt_jobs: Vec<Vec<usize>> = Vec::new();
    for (_, ix) in by_len {
        rt_jobs.extend(ix.chunks(args.batch).map(|c| c.to_vec()));
    }
    let rt_pred = Mutex::new(vec![f32::NAN; forms.len()]);
    let (rt_loads, rt_wall) = run_parallel(&args, &rt_jobs, |rt, _, job| {
        let fs: Vec<&Form> = job.iter().map(|&i| &forms[i].0).collect();
        let out = pred.rt(rt, &fs);
        let mut g = rt_pred.lock().unwrap();
        for (&i, v) in job.iter().zip(out) {
            g[i] = v;
        }
    });
    let rt_pred = rt_pred.into_inner().unwrap();
    eprintln!(
        "RT: {} forms in {} batches, wall {:.1}s, per-thread load {:.2}s",
        forms.len(),
        rt_jobs.len(),
        rt_wall,
        rt_loads.iter().sum::<f64>() / rt_loads.len() as f64
    );

    // MS2 jobs: distinct (form, charge) grouped by length.
    let mut fz_ix: HashMap<(usize, u8), usize> = HashMap::new();
    let mut fz: Vec<(usize, u8)> = Vec::new();
    for r in &rows {
        if let (Some(f), z) = (r.form, r.charge) {
            if z > 0 && forms[f].0.seq.len() >= 2 {
                fz_ix.entry((f, z)).or_insert_with(|| {
                    fz.push((f, z));
                    fz.len() - 1
                });
            }
        }
    }
    let mut by_len: HashMap<usize, Vec<usize>> = HashMap::new();
    for (i, (f, _)) in fz.iter().enumerate() {
        by_len.entry(forms[*f].0.seq.len()).or_default().push(i);
    }
    let mut ms2_jobs: Vec<Vec<usize>> = Vec::new();
    for (_, ix) in by_len {
        ms2_jobs.extend(ix.chunks(args.batch).map(|c| c.to_vec()));
    }
    let ms2_pred: Mutex<Vec<Vec<f32>>> = Mutex::new(vec![Vec::new(); fz.len()]);
    let (ms2_loads, ms2_wall) = run_parallel(&args, &ms2_jobs, |_, ms2, job| {
        let items: Vec<(&Form, u8)> = job.iter().map(|&i| (&forms[fz[i].0].0, fz[i].1)).collect();
        let out = pred.ms2(ms2, &items);
        let mut g = ms2_pred.lock().unwrap();
        for (&i, v) in job.iter().zip(out) {
            g[i] = v;
        }
    });
    let ms2_pred = ms2_pred.into_inner().unwrap();
    eprintln!(
        "MS2: {} (form, charge) in {} batches, wall {:.1}s, per-thread load {:.2}s",
        fz.len(),
        ms2_jobs.len(),
        ms2_wall,
        ms2_loads.iter().sum::<f64>() / ms2_loads.len() as f64
    );

    // Features per PIN row.
    let mut w = BufWriter::new(File::create(&args.out).unwrap());
    write!(w, "SpecId\tobserved_retention_time\tpdeep_irt\tpdeep_missing").unwrap();
    for n in MS2_FEATURES {
        write!(w, "\t{n}").unwrap();
    }
    writeln!(w).unwrap();
    let mut n_missing = 0usize;
    for r in &rows {
        let spec = by_scan.get(&r.scan).map(|&i| &spectra[i]);
        let obs_rt = spec.and_then(|s| s.rt_seconds).map_or(f64::NAN, |s| s / 60.0);
        let feats = match (r.form, spec, fz_ix.get(&(r.form.unwrap_or(usize::MAX), r.charge))) {
            (Some(f), Some(s), Some(&k)) => Some((
                rt_pred[f],
                ms2_features(&ms2_pred[k], &forms[f].1, r.charge, &s.peaks, args.frag_tol, args.frag_ppm),
            )),
            _ => None,
        };
        match feats {
            Some((irt, fv)) => {
                write!(w, "{}\t{obs_rt:.4}\t{irt:.6}\t0", r.id).unwrap();
                for v in fv {
                    write!(w, "\t{v:.6}").unwrap();
                }
            }
            None => {
                n_missing += 1;
                write!(w, "{}\t{obs_rt:.4}\tNA\t1", r.id).unwrap();
                for _ in MS2_FEATURES {
                    write!(w, "\tNA").unwrap();
                }
            }
        }
        writeln!(w).unwrap();
    }
    w.flush().unwrap();
    eprintln!(
        "wrote {} rows ({} missing) to {}; total {:.1}s; peak RSS {:.0} MB; threads {}",
        rows.len(),
        n_missing,
        args.out.display(),
        t_all.elapsed().as_secs_f64(),
        rss_mb(),
        args.threads
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_andes_peptides() {
        let (f, m) = parse_peptide("K.M+15.99491EEHILR.P").unwrap();
        assert_eq!(f.seq, "MEEHILR");
        assert_eq!(f.mods, vec![(1, "Oxidation@M")]);
        assert!((m.res[0] - (131.0404846 + 15.99491)).abs() < 1e-6);
        let (f, m) = parse_peptide("M.R+42.01056AC+57.02146INNQQIR.H").unwrap();
        assert_eq!(f.seq, "RACINNQQIR");
        assert_eq!(f.mods, vec![(0, "Acetyl@Protein_N-term"), (3, "Carbamidomethyl@C")]);
        assert!((m.nterm - 42.01056).abs() < 1e-9);
        let (f, _) = parse_peptide("_.M+42.01056+15.99491LRK.A").unwrap();
        assert_eq!(f.mods, vec![(0, "Acetyl@Protein_N-term"), (1, "Oxidation@M")]);
        assert!(parse_peptide("K.PEP+79.96633TIDE.K").is_none());
        assert!(parse_peptide("K.PEPXIDE.K").is_none());
    }

    /// The PIN string andes writes for a protein N-terminal Cys under fixed
    /// carbamidomethyl-C and variable protein N-terminal acetyl.
    fn andes_cam_acetyl_nterm_cys() -> String {
        use model::{AminoAcidSetBuilder, ModLocation, Modification, Peptide};
        let aa_set = AminoAcidSetBuilder::new_standard()
            .add_fixed_mod(
                Modification::from_mods_txt_line("57.021464,C,fix,any,Carbamidomethyl").unwrap(),
            )
            .add_variable_mod(
                Modification::from_mods_txt_line("42.010565,*,opt,Prot-N-term,Acetyl").unwrap(),
            )
            .build()
            .unwrap();
        let first = aa_set
            .variants_for(b'C', ModLocation::ProtNTerm)
            .iter()
            .find(|aa| {
                aa.mod_
                    .as_ref()
                    .is_some_and(|m| (m.mass_delta - 99.032029).abs() < 1e-3)
            })
            .expect("stacked CAM + acetyl variant")
            .clone();
        let mut residues = vec![first];
        for r in b"PEPTIDEK" {
            residues.push(aa_set.standard(*r).unwrap().clone());
        }
        Peptide::new(residues, b'M', b'R').to_string()
    }

    #[test]
    fn parses_summed_nterm_and_residue_delta_from_the_andes_writer() {
        let s = andes_cam_acetyl_nterm_cys();
        assert_eq!(s, "M.C+99.03203PEPTIDEK.R");
        let (f, m) = parse_peptide(&s).expect("summed CAM + acetyl delta");
        assert_eq!(f.seq, "CPEPTIDEK");
        assert_eq!(
            f.mods,
            vec![(0, "Acetyl@Protein_N-term"), (1, "Carbamidomethyl@C")]
        );
        assert!((m.nterm - 42.01056).abs() < 1e-4);
        assert!((m.res[0] - (103.009_184_5 + 57.02146)).abs() < 1e-4);
        assert!((m.nterm + m.res[0] - (103.009_184_5 + 99.03203)).abs() < 1e-9);
    }
}
