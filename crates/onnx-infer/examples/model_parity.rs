//! Parity of onnx-infer predictions (tract or ONNX Runtime, `--backend`) with a
//! model's reference Python implementation, from the TSV its
//! `scripts/<model>_reference.py` writes.
//!
//! ```text
//! cargo run --release -p onnx-infer --features tract --example model_parity -- \
//!   --model ms2pip --reference ms2pip_reference.tsv --models DIR [--backend ort]
//! ```
//!
//! Prints the largest absolute and relative deviation per model variant.
//! `--model peptdeep` reads `peptdeep_reference.tsv` from
//! `scripts/peptdeep_reference.py`.

#[path = "common/models.rs"]
mod models;

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

use onnx_infer::{BackendOptions, OnnxSession, ms2pip};

struct Args {
    model: String,
    reference: PathBuf,
    models: PathBuf,
}

fn parse_args() -> Args {
    let v: Vec<String> = std::env::args().skip(1).collect();
    let opt = |k: &str| v.iter().position(|a| a == k).and_then(|i| v.get(i + 1)).cloned();
    let get = |k: &str| opt(k).unwrap_or_else(|| panic!("missing {k}"));
    let backend = models::parse_backend(opt("--backend").as_deref().unwrap_or("tract"));
    models::set_backend(backend, BackendOptions::default());
    println!("backend {backend:?}");
    Args {
        model: get("--model"),
        reference: get("--reference").into(),
        models: get("--models").into(),
    }
}

fn floats(s: &str) -> Vec<f32> {
    s.split(',').filter(|x| !x.is_empty()).map(|x| x.parse().expect("float")).collect()
}

#[derive(Default)]
struct Dev {
    abs: f64,
    rel: f64,
    n: usize,
}

impl Dev {
    fn add(&mut self, got: f32, want: f32) {
        let d = (got as f64 - want as f64).abs();
        self.abs = self.abs.max(d);
        self.rel = self.rel.max(d / (want as f64).abs().max(1e-6));
        self.n += 1;
    }
}

fn load(dir: &Path, stem: &str) -> OnnxSession {
    let t = Instant::now();
    let s = models::load(dir, &format!("{stem}.json"));
    println!("load {stem}: {:.2}s", t.elapsed().as_secs_f64());
    s
}

fn ms2pip_parity(a: &Args) {
    let text = fs::read_to_string(&a.reference).expect("reference TSV");
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let col = |n: &str| header.iter().position(|h| *h == n).unwrap();
    let mut by_model: BTreeMap<String, Vec<Vec<&str>>> = BTreeMap::new();
    for l in lines {
        let x: Vec<&str> = l.split('\t').collect();
        by_model.entry(x[col("model")].to_string()).or_default().push(x);
    }
    for (model, rows) in by_model {
        if !a.models.join(format!("ms2pip_{model}_b.json")).is_file() {
            println!("SKIP ms2pip {model}: no ms2pip_{model}_b.json in {}", a.models.display());
            continue;
        }
        let mut sb = load(&a.models, &format!("ms2pip_{model}_b"));
        let mut sy = load(&a.models, &format!("ms2pip_{model}_y"));
        let (mut feat, mut pb, mut py) = (Dev::default(), Dev::default(), Dev::default());
        for x in &rows {
            let seq = x[col("sequence")];
            let z: u8 = x[col("charge")].parse().unwrap();
            let batch = ms2pip::pack_batch(&[(seq, z)]).unwrap();
            for (g, w) in batch.features.data.iter().zip(floats(x[col("features")])) {
                feat.add(*g, w);
            }
            let input = std::slice::from_ref(&batch.features);
            let b = ms2pip::log2_predictions(&sb.run(input).unwrap()[0]).unwrap();
            let y = ms2pip::log2_predictions(&sy.run(input).unwrap()[0]).unwrap();
            let (rb, ry) = (floats(x[col("b_log2")]), floats(x[col("y_log2")]));
            assert_eq!((b.len(), y.len()), (rb.len(), ry.len()), "{seq} row count");
            b.iter().zip(&rb).for_each(|(g, w)| pb.add(*g, *w));
            y.iter().zip(&ry).for_each(|(g, w)| py.add(*g, *w));
        }
        println!(
            "PARITY ms2pip {model}: {} peptides | features max|Δ| {:.3e} ({} values) | \
             b log2 max|Δ| {:.3e} rel {:.3e} | y log2 max|Δ| {:.3e} rel {:.3e}",
            rows.len(),
            feat.abs,
            feat.n,
            pb.abs,
            pb.rel,
            py.abs,
            py.rel
        );
    }
}

fn peptdeep_parity(a: &Args) {
    use models::{PeptDeep, Predictor, parse_peptide};
    let pred = PeptDeep::new(&a.models, 30.0, "Lumos", models::engine(), true, true).expect("PeptDeep");
    let mut s = pred.sessions().expect("PeptDeep sessions");
    let text = fs::read_to_string(&a.reference).expect("reference TSV");
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let col = |n: &str| header.iter().position(|h| *h == n).unwrap();
    let (mut rt, mut ms2) = (Dev::default(), Dev::default());
    let nft = pred.ms2_frag_types().len();
    for l in lines {
        let x: Vec<&str> = l.split('\t').collect();
        let seq = x[col("sequence")];
        // The andes spelling must parse to peptdeep's own mod names and sites.
        let andes = x[col("andes")];
        let (form, _) = parse_peptide(andes).unwrap_or_else(|| panic!("parse {andes}"));
        let got: Vec<(String, String)> = form
            .mods
            .iter()
            .map(|(site, m)| (m.alphabase(*site, form.residue(*site)), site.to_string()))
            .collect();
        let names = x[col("mods")].split(';').filter(|s| !s.is_empty());
        let sites = x[col("mod_sites")].split(';').filter(|s| !s.is_empty());
        let mut want: Vec<(String, String)> = names.zip(sites).map(|(n, s)| (n.into(), s.into())).collect();
        let mut got_sorted = got.clone();
        got_sorted.sort();
        want.sort();
        assert_eq!(form.seq, seq, "{andes}");
        assert_eq!(got_sorted, want, "{andes}: parsed mods differ from peptdeep's");
        let z: u8 = x[col("charge")].parse().unwrap();
        rt.add(pred.rt(&mut s, &[&form]).unwrap()[0], x[col("rt_peptdeep")].parse().unwrap());
        let m = &pred.ms2(&mut s, &[(&form, z)]).unwrap()[0];
        let r = floats(x[col("ms2_peptdeep")]);
        for site in 0..seq.len() - 1 {
            for (k, &c) in pred.ms2_columns().iter().enumerate() {
                ms2.add(m[site * 4 + k], r[site * nft + c]);
            }
        }
    }
    println!(
        "PARITY peptdeep: rt max|Δ| {:.3e} rel {:.3e} | ms2 max|Δ| {:.3e} ({} values)",
        rt.abs, rt.rel, ms2.abs, ms2.n
    );
}

fn deeplc_parity(a: &Args) {
    use models::{DeepLc, Predictor, parse_peptide};
    let pred = DeepLc::new(&a.models, models::engine()).expect("DeepLC");
    let mut s = pred.sessions().expect("DeepLC session");
    let text = fs::read_to_string(&a.reference).expect("reference TSV");
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let col = |n: &str| header.iter().position(|h| *h == n).unwrap();
    let (mut enc, mut rt) = (Dev::default(), Dev::default());
    let mut n = 0;
    for l in lines {
        let x: Vec<&str> = l.split('\t').collect();
        let (form, _) = parse_peptide(x[col("andes")]).unwrap_or_else(|| panic!("parse {}", x[col("andes")]));
        let f = onnx_infer::deeplc::encode(&DeepLc::peptide(&form)).unwrap();
        for (k, got) in [("x_atom", &f.x_atom), ("x_global", &f.x_global), ("x_one_hot", &f.x_one_hot)] {
            let want = floats(x[col(k)]);
            assert_eq!(got.len(), want.len(), "{k} length");
            for (i, (g, w)) in got.iter().zip(&want).enumerate() {
                if g != w {
                    println!("{} {k}[{i}]: rust {g} deeplc {w}", x[col("proforma")]);
                }
                enc.add(*g, *w);
            }
        }
        let r = pred.rt(&mut s, &[&form]).unwrap()[0];
        let w: f32 = x[col("rt_deeplc")].parse().unwrap();
        rt.add(r, w);
        println!("{:60} onnx-infer {r:12.4} deeplc {w:12.4}", x[col("proforma")]);
        n += 1;
    }
    println!(
        "PARITY deeplc: {n} peptides | encoding max|Δ| {:.3e} ({} values) | rt max|Δ| {:.3e} min, max rel {:.3e}",
        enc.abs, enc.n, rt.abs, rt.rel
    );
}

fn prosit_parity(a: &Args) {
    use models::{Predictor, Prosit, parse_peptide};
    let text = fs::read_to_string(&a.reference).expect("reference TSV");
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let col = |n: &str| header.iter().position(|h| *h == n).unwrap();
    let mut s = Prosit { dir: a.models.clone(), nce: 30.0 }.sessions().unwrap();
    let (mut int, mut irt) = (Dev::default(), Dev::default());
    let mut n = 0;
    for l in lines {
        let x: Vec<&str> = l.split('\t').collect();
        let (form, _) = parse_peptide(x[col("andes")]).unwrap();
        let z: u8 = x[col("charge")].parse().unwrap();
        let pred = Prosit { dir: a.models.clone(), nce: x[col("nce")].parse().unwrap() };
        let got = &pred.ms2(&mut s, &[(&form, z)]).unwrap()[0];
        // Koina: 29 positions x (y+1..3, b+1..3), max-normalised, -1 = impossible.
        let koina = floats(x[col("intensities_koina")]);
        let len = form.seq.len();
        for site in 0..len - 1 {
            for k in 0..3 {
                let kb = koina[site * 6 + 3 + k];
                if kb >= 0.0 {
                    int.add(got[site * 6 + k], kb);
                }
                let ky = koina[(len - 2 - site) * 6 + k];
                if ky >= 0.0 {
                    int.add(got[site * 6 + 3 + k], ky);
                }
            }
        }
        let r = pred.rt(&mut s, &[&form]).unwrap()[0] * 43.39373 + 56.353_634;
        irt.add(r, x[col("irt_koina")].parse().unwrap());
        n += 1;
    }
    println!(
        "PARITY prosit (vs Koina): {n} peptides | normalised intensity max|Δ| {:.3e} ({} values) | \
         iRT max|Δ| {:.3e} rel {:.3e}",
        int.abs, int.n, irt.abs, irt.rel
    );
}

fn main() {
    let a = parse_args();
    match a.model.as_str() {
        "prosit" => prosit_parity(&a),
        "deeplc" => deeplc_parity(&a),
        "ms2pip" => ms2pip_parity(&a),
        "peptdeep" => peptdeep_parity(&a),
        other => panic!("unknown --model {other}"),
    }
}
