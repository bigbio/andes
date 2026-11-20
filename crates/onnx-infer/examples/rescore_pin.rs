#[path = "common/models.rs"]
mod models;

use std::collections::hash_map::Entry;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;
use std::time::Instant;

use models::{parse_peptide, AndesMs2, DeepLc, Form, Mass, Ms2pip, PeptDeep, Predictor, Prosit};
use onnx_infer::predictors::spectral::{ms2_features, FragmentTolerance, MS2_FEATURES};
use onnx_infer::{BackendKind, BackendOptions, OnnxSession};

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
    activation: String,
    ms2pip_model: String,
    backend: BackendKind,
    intra_threads: usize,
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
        activation: "HCD-FTMS".into(),
        ms2pip_model: "HCD2021".into(),
        backend: BackendKind::Tract,
        intra_threads: 1,
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
            "--nce" => {
                a.nce = if val == "NA" {
                    f32::NAN
                } else {
                    val.parse().expect("--nce")
                }
            }
            "--instrument" => a.instrument = val,
            "--activation" => a.activation = val,
            "--ms2pip-model" => a.ms2pip_model = val,
            "--backend" => a.backend = models::parse_backend(&val),
            "--intra-threads" => a.intra_threads = val.parse().expect("--intra-threads"),
            other => panic!("unknown argument {other}"),
        }
        i += 2;
    }
    a
}

fn cpu_s() -> f64 {
    std::fs::read_to_string("/proc/self/stat")
        .ok()
        .and_then(|s| {
            let rest = s.rsplit_once(')')?.1;
            let f: Vec<&str> = rest.split_whitespace().collect();
            let t = f.get(11)?.parse::<f64>().ok()? + f.get(12)?.parse::<f64>().ok()?;
            Some(t / 100.0)
        })
        .unwrap_or(0.0)
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

fn run_parallel<J: Sync, F>(
    master: &[OnnxSession],
    threads: usize,
    jobs: &[J],
    work: F,
) -> (f64, f64)
where
    F: Fn(&mut [OnnxSession], &J) + Sync,
{
    let next = AtomicUsize::new(0);
    let t0 = Instant::now();
    let c0 = cpu_s();
    let handles: Vec<Vec<OnnxSession>> = (0..threads)
        .map(|_| {
            master
                .iter()
                .map(|s| s.share().expect("share session"))
                .collect()
        })
        .collect();
    std::thread::scope(|sc| {
        for mut sessions in handles {
            let (next, work) = (&next, &work);
            sc.spawn(move || loop {
                let j = next.fetch_add(1, Ordering::Relaxed);
                if j >= jobs.len() {
                    break;
                }
                work(&mut sessions, &jobs[j]);
            });
        }
    });
    (t0.elapsed().as_secs_f64(), cpu_s() - c0)
}

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
    models::set_backend(
        args.backend,
        BackendOptions {
            intra_threads: args.intra_threads,
        },
    );
    let pred: Box<dyn Predictor> = match args.model.as_str() {
        "peptdeep" => Box::new(
            PeptDeep::new(
                &args.models,
                args.nce,
                &args.instrument,
                models::engine(),
                true,
                true,
            )
            .unwrap_or_else(|e| panic!("{e}")),
        ),
        "ms2pip" => Box::new(Ms2pip {
            dir: args.models.clone(),
            model: args.ms2pip_model.clone(),
        }),
        "deeplc" => {
            Box::new(DeepLc::new(&args.models, models::engine()).unwrap_or_else(|e| panic!("{e}")))
        }
        "prosit" => Box::new(Prosit {
            dir: args.models.clone(),
            nce: args.nce,
        }),
        "andes-ms2" => Box::new(AndesMs2::new(
            &args.models,
            Some(args.nce).filter(|v| v.is_finite()),
            &args.activation,
            &args.instrument,
        )),
        other => panic!("unknown --model {other}"),
    };
    let pred = pred.as_ref();

    let t_all = Instant::now();
    let f = BufReader::new(File::open(&args.pin).unwrap());
    let mut lines = f.lines();
    let header: Vec<String> = lines
        .next()
        .unwrap()
        .unwrap()
        .split('\t')
        .map(String::from)
        .collect();
    let hc = |n: &str| {
        header
            .iter()
            .position(|h| h == n)
            .unwrap_or_else(|| panic!("PIN column {n}"))
    };
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
    let n_cols = charge_cols
        .iter()
        .map(|(i, _)| *i)
        .chain([c_id, c_scan, c_pep])
        .max()
        .unwrap()
        + 1;
    let mut rows = Vec::new();
    let mut skipped_rows = 0usize;
    for line in lines {
        let line = line.unwrap();
        let x: Vec<&str> = line.split('\t').collect();
        let scan = match x.get(c_scan).map(|v| v.trim().parse::<i32>()) {
            Some(Ok(scan)) if x.len() >= n_cols => scan,
            _ => {
                skipped_rows += 1;
                continue;
            }
        };
        let charge = charge_cols
            .iter()
            .find(|(i, _)| x[*i] == "1")
            .map_or(0, |(_, z)| *z);
        let form = parse_peptide(x[c_pep]).map(|(f, m)| {
            *form_ix.entry(f.clone()).or_insert_with(|| {
                forms.push((f, m));
                forms.len() - 1
            })
        });
        rows.push(Row {
            id: x[c_id].to_string(),
            scan,
            charge,
            form,
        });
    }
    if skipped_rows > 0 {
        eprintln!("warning: skipped {skipped_rows} PIN lines (blank, short or non-numeric ScanNr)");
    }
    let unparsed = rows
        .iter()
        .filter(|r| r.form.is_none() || r.charge == 0)
        .count();
    eprintln!(
        "model {} backend {:?} | PIN rows {} distinct peptidoforms {} unparsed rows {} ({:.1}s)",
        args.model,
        args.backend,
        rows.len(),
        forms.len(),
        unparsed,
        t_all.elapsed().as_secs_f64()
    );

    let t = Instant::now();
    let (spectra, errs) = input::read_spectra_parallel(&args.mzml, 2, 2, args.threads).unwrap();
    assert!(
        errs.is_empty(),
        "mzML errors: {:?}",
        &errs[..errs.len().min(3)]
    );
    let mut by_scan: HashMap<i32, usize> = HashMap::new();
    let (mut no_scan, mut dup_scan) = (0usize, 0usize);
    for (i, s) in spectra.iter().enumerate() {
        match s.scan {
            None => no_scan += 1,
            Some(scan) => match by_scan.entry(scan) {
                Entry::Occupied(_) => dup_scan += 1,
                Entry::Vacant(e) => {
                    e.insert(i);
                }
            },
        }
    }
    if no_scan + dup_scan > 0 {
        eprintln!(
            "warning: skipped {no_scan} spectra without a scan number and {dup_scan} with a duplicate scan number"
        );
    }
    eprintln!(
        "mzML {} MS2 spectra ({:.1}s)",
        spectra.len(),
        t.elapsed().as_secs_f64()
    );

    let t = Instant::now();
    let master = pred.sessions().unwrap_or_else(|e| panic!("{e}"));
    eprintln!(
        "model load {:.2}s, RSS {:.0} MB",
        t.elapsed().as_secs_f64(),
        rss_mb()
    );

    let mut rt_pred = vec![f32::NAN; forms.len()];
    if pred.has_rt() {
        let rt_forms: Vec<usize> = (0..forms.len())
            .filter(|&i| pred.supports_rt(&forms[i].0))
            .collect();
        let jobs = batches(
            rt_forms.iter().map(|&i| pred.batch_key(&forms[i].0)),
            args.batch,
        );
        let out = Mutex::new(&mut rt_pred);
        let (wall, cpu) = run_parallel(&master, args.threads, &jobs, |s, job| {
            let fs: Vec<&Form> = job.iter().map(|&j| &forms[rt_forms[j]].0).collect();
            let v = pred.rt(s, &fs).expect("RT prediction");
            let mut g = out.lock().unwrap();
            for (&j, x) in job.iter().zip(v) {
                g[rt_forms[j]] = x;
            }
        });
        eprintln!(
            "RT: {} forms in {} batches, wall {wall:.1}s cpu {cpu:.1}s ({:.0} us/form cpu)",
            rt_forms.len(),
            jobs.len(),
            cpu * 1e6 / rt_forms.len().max(1) as f64
        );
    }

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
        let jobs = batches(
            fz.iter().map(|(f, _)| pred.batch_key(&forms[*f].0)),
            args.batch,
        );
        let out = Mutex::new(&mut ms2_pred);
        let (wall, cpu) = run_parallel(&master, args.threads, &jobs, |s, job| {
            let items: Vec<(&Form, u8)> =
                job.iter().map(|&i| (&forms[fz[i].0].0, fz[i].1)).collect();
            let v = pred.ms2(s, &items).expect("MS2 prediction");
            let mut g = out.lock().unwrap();
            for (&i, x) in job.iter().zip(v) {
                g[i] = x;
            }
        });
        eprintln!(
            "MS2: {} (form, charge) in {} batches, wall {wall:.1}s cpu {cpu:.1}s ({:.0} us/item cpu)",
            fz.len(),
            jobs.len(),
            cpu * 1e6 / fz.len() as f64
        );
    }

    let p = pred.prefix();
    let tol = if args.frag_ppm {
        FragmentTolerance::Ppm(args.frag_tol)
    } else {
        FragmentTolerance::Da(args.frag_tol)
    };
    let mut w = BufWriter::new(File::create(&args.out).unwrap());
    write!(w, "SpecId\tobserved_retention_time\t{p}_irt\t{p}_missing").unwrap();
    for n in MS2_FEATURES {
        write!(w, "\t{n}").unwrap();
    }
    writeln!(w).unwrap();
    let mut n_missing = 0usize;
    for r in &rows {
        let spec = by_scan.get(&r.scan).map(|&i| &spectra[i]);
        let obs_rt = spec
            .and_then(|s| s.rt_seconds)
            .map_or(f64::NAN, |s| s / 60.0);
        let irt = r.form.map_or(f32::NAN, |f| rt_pred[f]);
        let irt = if irt.is_finite() { irt } else { f32::NAN };
        let ms2 = match (r.form, spec) {
            (Some(f), Some(s)) => fz_ix.get(&(f, r.charge)).map(|&k| {
                ms2_features(
                    &ms2_pred[k],
                    pred.ions(),
                    None,
                    &forms[f].1,
                    r.charge,
                    &s.peaks,
                    tol,
                )
            }),
            _ => None,
        };
        let missing = (pred.has_rt() && irt.is_nan()) || (!pred.ions().is_empty() && ms2.is_none());
        n_missing += missing as usize;
        let irt = if irt.is_nan() {
            "NA".to_string()
        } else {
            format!("{irt:.6}")
        };
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
