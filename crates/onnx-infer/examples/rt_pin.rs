//! andes-rt retention index for every row of a Percolator PIN.
//!
//! ```text
//! cargo run --release -p onnx-infer --example rt_pin -- \
//!   --pin search.pin --mzml run.mzML --weights andes_rt.artn --out rt.tsv [--threads 8]
//! ```
//!
//! Writes `SpecId, observed_retention_time (min), andes_rt, rt_missing, proforma`
//! per PIN row (targets and decoys alike). `proforma` (Unimod accessions) lets
//! other predictors score the same peptidoforms. Calibration to minutes happens
//! downstream (`scripts/rt_calibration.py`).

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::time::Instant;

use onnx_infer::rt_native::{RtModel, RtPeptide, Scratch};

fn main() {
    let v: Vec<String> = std::env::args().skip(1).collect();
    let get = |k: &str| v.iter().position(|x| x == k).map(|i| v[i + 1].clone());
    let need = |k: &str| get(k).unwrap_or_else(|| panic!("missing {k}"));
    let threads: usize = get("--threads").map_or(8, |t| t.parse().unwrap());
    let model = RtModel::load(need("--weights")).expect("weights");

    let t0 = Instant::now();
    let f = BufReader::new(File::open(need("--pin")).unwrap());
    let mut lines = f.lines();
    let h: Vec<String> = lines
        .next()
        .unwrap()
        .unwrap()
        .split('\t')
        .map(String::from)
        .collect();
    let hc = |n: &str| {
        h.iter()
            .position(|x| x == n)
            .unwrap_or_else(|| panic!("PIN column {n}"))
    };
    let (c_id, c_scan, c_pep) = (hc("SpecId"), hc("ScanNr"), hc("Peptide"));
    let mut forms: Vec<RtPeptide> = Vec::new();
    let mut ix: HashMap<RtPeptide, usize> = HashMap::new();
    let mut rows: Vec<(String, i32, Option<usize>)> = Vec::new();
    for line in lines {
        let line = line.unwrap();
        let x: Vec<&str> = line.split('\t').collect();
        let form = RtPeptide::from_andes(x[c_pep]).map(|p| {
            *ix.entry(p.clone()).or_insert_with(|| {
                forms.push(p);
                forms.len() - 1
            })
        });
        rows.push((x[c_id].to_string(), x[c_scan].parse().unwrap(), form));
    }

    let (spectra, errs) =
        input::read_spectra_parallel(std::path::Path::new(&need("--mzml")), 2, 2, threads).unwrap();
    assert!(errs.is_empty(), "mzML errors");
    let rt_by_scan: HashMap<i32, f64> = spectra
        .iter()
        .map(|s| {
            (
                s.scan.expect("scan"),
                s.rt_seconds.map_or(f64::NAN, |t| t / 60.0),
            )
        })
        .collect();

    let tp = Instant::now();
    let mut pred = vec![f32::NAN; forms.len()];
    std::thread::scope(|sc| {
        let chunk = forms.len().div_ceil(threads).max(1);
        for (fs, out) in forms.chunks(chunk).zip(pred.chunks_mut(chunk)) {
            let model = &model;
            sc.spawn(move || {
                let mut s = Scratch::default();
                for (p, o) in fs.iter().zip(out.iter_mut()) {
                    *o = model.predict_with(p, &mut s).unwrap_or(f32::NAN);
                }
            });
        }
    });
    let pred_s = tp.elapsed().as_secs_f64();

    let mut w = BufWriter::new(File::create(need("--out")).unwrap());
    writeln!(
        w,
        "SpecId\tobserved_retention_time\tandes_rt\trt_missing\tproforma"
    )
    .unwrap();
    let mut miss = 0;
    for (id, scan, form) in &rows {
        let obs = rt_by_scan.get(scan).copied().unwrap_or(f64::NAN);
        match form.filter(|&f| pred[f].is_finite()) {
            Some(f) => writeln!(
                w,
                "{id}\t{obs:.4}\t{:.6}\t0\t{}",
                pred[f],
                forms[f].to_proforma()
            )
            .unwrap(),
            None => {
                miss += 1;
                writeln!(w, "{id}\t{obs:.4}\tNA\t1\tNA").unwrap()
            }
        }
    }
    w.flush().unwrap();
    eprintln!(
        "rows {} forms {} missing {miss}; prediction {pred_s:.2}s on {threads} threads ({:.1} us/form/thread); total {:.1}s",
        rows.len(),
        forms.len(),
        pred_s * threads as f64 / forms.len().max(1) as f64 * 1e6,
        t0.elapsed().as_secs_f64()
    );
}
