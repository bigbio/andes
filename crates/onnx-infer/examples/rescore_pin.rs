//! Predicted-spectrum and RT features for every row of a Percolator PIN, with tract.
//!
//! Reads an andes PIN and the mzML it was searched from, predicts fragment
//! intensities and/or retention time for each distinct peptidoform with one
//! model family's exported ONNX graphs, and writes one feature row per PIN row
//! (targets and decoys alike) to a TSV keyed by `SpecId`. RT calibration and
//! PIN assembly happen downstream (`scripts/rescore_pin_arms.py`).
//!
//! ```text
//! cargo run --release -p onnx-infer --features tract --example rescore_pin -- \
//!   --model peptdeep --pin search.pin --mzml run.mzML --models DIR --out features.tsv \
//!   --frag-tol 0.5 --frag-unit da --threads 8
//! ```
//!
//! `--model` selects the family (feature-column prefix in brackets):
//!
//! | `--model` | Graphs in `DIR` | Outputs |
//! |---|---|---|
//! | `peptdeep` (`pdeep`) | `peptdeep_{rt,ms2}_dynamic.json`, `peptdeep_mods.json` | RT + b/y z1–2 |
//! | `ms2pip` (`ms2pip`) | `ms2pip_<M>_{b,y}.json`, `M` = `--ms2pip-model` | b/y z1 |
//! | `deeplc` (`deeplc`) | `deeplc_rt.json` | RT |
//! | `prosit` (`prosit`) | `prosit_irt.json`, `prosit_intensity_hcd.json` (`--nce`) | iRT + b/y z1–3 |
//!
//! Output columns: `SpecId`, `observed_retention_time`, `<prefix>_irt` (NA
//! without an RT model), `<prefix>_missing` (1 when the model cannot encode the
//! row), then the spectral-match features (NA without an MS2 model).

#[path = "common/models.rs"]
mod models;

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

use models::{DeepLc, Form, H2O, Ion, Mass, Ms2pip, PROTON, PeptDeep, Predictor, Prosit, parse_peptide};
use onnx_infer::OnnxSession;

struct Args {
    model: String,
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
    ms2pip_model: String,
}

fn parse_args() -> Args {
    let mut a = Args {
        model: "peptdeep".into(),
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
        ms2pip_model: "HCD2021".into(),
    };
    let v: Vec<String> = std::env::args().skip(1).collect();
    let mut i = 0;
    while i < v.len() {
        let val = v.get(i + 1).cloned().unwrap_or_default();
        match v[i].as_str() {
            "--model" => a.model = val,
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
            "--ms2pip-model" => a.ms2pip_model = val,
            other => panic!("unknown argument {other}"),
        }
        i += 2;
    }
    a
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

/// Spectral-match features between predicted and observed fragment intensities.
/// `pred` holds `ions.len()` values per cleavage site.
fn ms2_features(
    pred: &[f32],
    ions: &[Ion],
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
        for (k, &(bion, z)) in ions.iter().enumerate() {
            if z > charge {
                continue;
            }
            let neutral = if bion { prefix } else { total - prefix + H2O };
            let mz = (neutral + z as f64 * PROTON) / z as f64;
            let w = if ppm { mz * tol * 1e-6 } else { tol };
            let lo = peaks.partition_point(|(m, _)| *m < mz - w);
            let obs = peaks[lo..]
                .iter()
                .take_while(|(m, _)| *m <= mz + w)
                .map(|(_, i)| *i as f64)
                .fold(0.0, f64::max);
            p.push(pred[i * ions.len() + k] as f64);
            o.push(obs);
            is_b.push(bion);
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

/// Run `jobs` across worker threads, each with shared handles on `master`.
fn run_parallel<J: Sync, F>(master: &[OnnxSession], threads: usize, jobs: &[J], work: F) -> f64
where
    F: Fn(&mut [OnnxSession], &J) + Sync,
{
    let next = AtomicUsize::new(0);
    let t0 = Instant::now();
    let handles: Vec<Vec<OnnxSession>> = (0..threads)
        .map(|_| master.iter().map(|s| s.share().expect("share session")).collect())
        .collect();
    std::thread::scope(|sc| {
        for mut sessions in handles {
            let (next, work) = (&next, &work);
            sc.spawn(move || {
                loop {
                    let j = next.fetch_add(1, Ordering::Relaxed);
                    if j >= jobs.len() {
                        break;
                    }
                    work(&mut sessions, &jobs[j]);
                }
            });
        }
    });
    t0.elapsed().as_secs_f64()
}

/// Group item indices by `key` and chunk each group into batches.
fn batches(keys: impl Iterator<Item = usize>, batch: usize) -> Vec<Vec<usize>> {
    let mut by: HashMap<usize, Vec<usize>> = HashMap::new();
    for (i, k) in keys.enumerate() {
        by.entry(k).or_default().push(i);
    }
    let mut groups: Vec<_> = by.into_iter().collect();
    groups.sort_by_key(|(k, _)| *k);
    groups
        .into_iter()
        .flat_map(|(_, ix)| ix.chunks(batch).map(|c| c.to_vec()).collect::<Vec<_>>())
        .collect()
}

fn main() {
    let args = parse_args();
    let pred: Box<dyn Predictor> = match args.model.as_str() {
        "peptdeep" => Box::new(PeptDeep::new(&args.models, args.nce, &args.instrument)),
        "ms2pip" => Box::new(Ms2pip { dir: args.models.clone(), model: args.ms2pip_model.clone() }),
        "deeplc" => Box::new(DeepLc { dir: args.models.clone() }),
        "prosit" => Box::new(Prosit { dir: args.models.clone(), nce: args.nce }),
        other => panic!("unknown --model {other}"),
    };
    let pred = pred.as_ref();

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
        rows.push(Row { id: x[c_id].to_string(), scan: x[c_scan].parse().unwrap(), charge, form });
    }
    let unparsed = rows.iter().filter(|r| r.form.is_none() || r.charge == 0).count();
    eprintln!(
        "model {} | PIN rows {} distinct peptidoforms {} unparsed rows {} ({:.1}s)",
        args.model,
        rows.len(),
        forms.len(),
        unparsed,
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

    let t = Instant::now();
    let master = pred.sessions();
    eprintln!("model load {:.2}s, RSS {:.0} MB", t.elapsed().as_secs_f64(), rss_mb());

    // RT over distinct supported forms.
    let mut rt_pred = vec![f32::NAN; forms.len()];
    if pred.has_rt() {
        let rt_forms: Vec<usize> = (0..forms.len()).filter(|&i| pred.supports_rt(&forms[i].0)).collect();
        let jobs = batches(rt_forms.iter().map(|&i| pred.batch_key(&forms[i].0)), args.batch);
        let out = Mutex::new(&mut rt_pred);
        let wall = run_parallel(&master, args.threads, &jobs, |s, job| {
            let fs: Vec<&Form> = job.iter().map(|&j| &forms[rt_forms[j]].0).collect();
            let v = pred.rt(s, &fs);
            let mut g = out.lock().unwrap();
            for (&j, x) in job.iter().zip(v) {
                g[rt_forms[j]] = x;
            }
        });
        eprintln!("RT: {} forms in {} batches, wall {wall:.1}s", rt_forms.len(), jobs.len());
    }

    // MS2 over distinct supported (form, charge).
    let mut fz_ix: HashMap<(usize, u8), usize> = HashMap::new();
    let mut fz: Vec<(usize, u8)> = Vec::new();
    if !pred.ions().is_empty() {
        for r in &rows {
            if let Some(f) = r.form {
                if pred.supports_ms2(&forms[f].0, r.charge) {
                    fz_ix.entry((f, r.charge)).or_insert_with(|| {
                        fz.push((f, r.charge));
                        fz.len() - 1
                    });
                }
            }
        }
    }
    let mut ms2_pred: Vec<Vec<f32>> = vec![Vec::new(); fz.len()];
    if !fz.is_empty() {
        let jobs = batches(fz.iter().map(|(f, _)| pred.batch_key(&forms[*f].0)), args.batch);
        let out = Mutex::new(&mut ms2_pred);
        let wall = run_parallel(&master, args.threads, &jobs, |s, job| {
            let items: Vec<(&Form, u8)> = job.iter().map(|&i| (&forms[fz[i].0].0, fz[i].1)).collect();
            let v = pred.ms2(s, &items);
            let mut g = out.lock().unwrap();
            for (&i, x) in job.iter().zip(v) {
                g[i] = x;
            }
        });
        eprintln!("MS2: {} (form, charge) in {} batches, wall {wall:.1}s", fz.len(), jobs.len());
    }

    // Features per PIN row.
    let p = pred.prefix();
    let mut w = BufWriter::new(File::create(&args.out).unwrap());
    write!(w, "SpecId\tobserved_retention_time\t{p}_irt\t{p}_missing").unwrap();
    for n in MS2_FEATURES {
        write!(w, "\t{n}").unwrap();
    }
    writeln!(w).unwrap();
    let mut n_missing = 0usize;
    for r in &rows {
        let spec = by_scan.get(&r.scan).map(|&i| &spectra[i]);
        let obs_rt = spec.and_then(|s| s.rt_seconds).map_or(f64::NAN, |s| s / 60.0);
        let irt = r.form.map_or(f32::NAN, |f| rt_pred[f]);
        let ms2 = match (r.form, spec) {
            (Some(f), Some(s)) => fz_ix.get(&(f, r.charge)).map(|&k| {
                ms2_features(&ms2_pred[k], pred.ions(), &forms[f].1, r.charge, &s.peaks, args.frag_tol, args.frag_ppm)
            }),
            _ => None,
        };
        let missing = (pred.has_rt() && irt.is_nan()) || (!pred.ions().is_empty() && ms2.is_none());
        n_missing += missing as usize;
        let irt = if irt.is_nan() { "NA".to_string() } else { format!("{irt:.6}") };
        write!(w, "{}\t{obs_rt:.4}\t{irt}\t{}", r.id, missing as u8).unwrap();
        match ms2 {
            Some(fv) => fv.iter().for_each(|v| write!(w, "\t{v:.6}").unwrap()),
            None => MS2_FEATURES.iter().for_each(|_| write!(w, "\tNA").unwrap()),
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
