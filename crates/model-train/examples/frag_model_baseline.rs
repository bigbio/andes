//! How good is the in-house fragment-intensity GBDT as a spectrum predictor?
//!
//! For a list of identified PSMs (`SpecId`, andes peptide string) and the mzML
//! they came from, predicts b/y 1+/2+ intensities with the store's
//! frag-intensity GBDT and compares them with the observed spectrum using the
//! same feature definitions as `onnx-infer/examples/rescore_pin.rs` (spectral
//! angle, cosine, Pearson on sum-normalised b/y vectors; fragment charge above
//! the precursor charge skipped), so the numbers line up column for column with
//! the PeptDeep features that example writes. A flat prediction (every ion
//! equal) is reported as a floor. The andes RT index features (`rt_features`)
//! and seed index are written as well for an RT baseline.
//!
//! ```text
//! cargo run --release -p model-train --example frag_model_baseline -- \
//!   --psms conf.tsv --mzml run.mzML --store resources/models --model-id hcd_astral_tryp \
//!   --frag-tol 20 --frag-unit ppm --out gbdt.tsv
//! ```
//!
//! `--psms` is a TSV with a header containing `SpecId` and `peptide`
//! (`K.PEPM+15.99491K.D`); the scan number is parsed from `scan=<n>` in SpecId.
//! The precursor charge is the charge whose theoretical m/z matches the
//! spectrum's precursor within 30 ppm (isotope offsets 0..=3).

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;

use model::amino_acid::AminoAcid;
use model::modification::{leak_mod, ModLocation, Modification, ResidueSpec};
use model::peptide::Peptide;
use model_train::store::ModelStore;
use scoring_crate::rt_model::{rt_features, RtIndexModel};
use scoring_crate::scoring::fragment_ions::IonKind;
use scoring_crate::scoring::{predict_by_ions, predict_frag_intensities};

const PROTON: f64 = 1.007_276_466_8;
const C13: f64 = 1.003_354_835;

/// andes mass delta → (residues it may sit on, empty = N-term only; mass; name).
const MODS: &[(&[u8], f64, &str)] = &[
    (b"C", 57.02146, "Carbamidomethyl"),
    (b"M", 15.99491, "Oxidation"),
    (b"", 42.01056, "Acetyl"),
    (b"STY", 79.96633, "Phospho"),
    (b"K", 229.16293, "TMT6plex"),
    (b"", 229.16293, "TMT6plex"),
    (b"NQ", 0.98402, "Deamidated"),
    (b"", -17.02655, "Gln->pyro-Glu"),
    (b"", -18.01056, "Glu->pyro-Glu"),
];

/// Parse `K.PEPM+15.99491K.D` → (residue bytes, per-residue mod delta).
/// N-terminal deltas follow the first residue and are folded into it.
fn parse(s: &str) -> Option<(Vec<u8>, Vec<f64>)> {
    let b = s.as_bytes();
    let core = if b.len() > 4 && b[1] == b'.' && b[b.len() - 2] == b'.' { &s[2..s.len() - 2] } else { s };
    let c: Vec<char> = core.chars().collect();
    let (mut seq, mut dm) = (Vec::new(), Vec::new());
    let mut i = 0;
    while i < c.len() {
        if c[i].is_ascii_uppercase() {
            seq.push(c[i] as u8);
            dm.push(0.0);
            i += 1;
        } else if c[i] == '+' || c[i] == '-' {
            let st = i;
            i += 1;
            while i < c.len() && (c[i].is_ascii_digit() || c[i] == '.') {
                i += 1;
            }
            let d: f64 = c[st..i].iter().collect::<String>().parse().ok()?;
            let pos = seq.len();
            if pos == 0 {
                return None;
            }
            let aa = seq[pos - 1];
            MODS.iter().find(|(r, m, _)| (d - m).abs() < 2e-3 && if r.is_empty() { pos == 1 } else { r.contains(&aa) })?;
            dm[pos - 1] += d;
        } else {
            return None;
        }
    }
    (!seq.is_empty()).then_some((seq, dm))
}

fn peptide(seq: &[u8], dm: &[f64], cache: &mut HashMap<(u8, u64), &'static Modification>) -> Option<Peptide> {
    let mut res = Vec::with_capacity(seq.len());
    for (&r, &d) in seq.iter().zip(dm) {
        let mut aa = AminoAcid::standard(r)?;
        if d != 0.0 {
            let m = *cache.entry((r, d.to_bits())).or_insert_with(|| {
                leak_mod(Modification {
                    name: format!("{}{d:+.5}", r as char),
                    mass_delta: d,
                    residue: ResidueSpec::Specific(r),
                    location: ModLocation::Anywhere,
                    fixed: false,
                    accession: None,
                    neutral_losses: Vec::new(),
                    loss_class: 0,
                })
            });
            aa = aa.with_mod(m);
        }
        res.push(aa);
    }
    Some(Peptide::new(res, b'_', b'_'))
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
    if saa <= 0.0 || sbb <= 0.0 { 0.0 } else { sab / (saa * sbb).sqrt() }
}

/// (spectral angle, cosine, Pearson) on sum-normalised vectors; identical to
/// `rescore_pin.rs::ms2_features` for these three columns.
fn sim(p: &[f64], o: &[f64]) -> [f64; 3] {
    let (sp, so) = (p.iter().sum::<f64>(), o.iter().sum::<f64>());
    let pn: Vec<f64> = p.iter().map(|v| if sp > 0.0 { v / sp } else { 0.0 }).collect();
    let on: Vec<f64> = o.iter().map(|v| if so > 0.0 { v / so } else { 0.0 }).collect();
    let dot: f64 = pn.iter().zip(&on).map(|(a, b)| a * b).sum();
    let np = pn.iter().map(|v| v * v).sum::<f64>().sqrt();
    let no = on.iter().map(|v| v * v).sum::<f64>().sqrt();
    let cos = if np > 0.0 && no > 0.0 { (dot / (np * no)).clamp(-1.0, 1.0) } else { 0.0 };
    let sa = if np > 0.0 && no > 0.0 { 1.0 - 2.0 * cos.acos() / std::f64::consts::PI } else { 0.0 };
    [sa, cos, pearson(&pn, &on)]
}

fn main() {
    let v: Vec<String> = std::env::args().skip(1).collect();
    let get = |k: &str| v.iter().position(|x| x == k).map(|i| v[i + 1].clone());
    let need = |k: &str| get(k).unwrap_or_else(|| panic!("missing {k}"));
    let psms = PathBuf::from(need("--psms"));
    let mzml = PathBuf::from(need("--mzml"));
    let store = PathBuf::from(need("--store"));
    let model_id = need("--model-id");
    let out = PathBuf::from(need("--out"));
    let tol: f64 = need("--frag-tol").parse().unwrap();
    let ppm = get("--frag-unit").as_deref() == Some("ppm");
    let threads: usize = get("--threads").map_or(8, |t| t.parse().unwrap());

    let st = ModelStore::open(&store).expect("open store");
    let param = st.load_param(&model_id).expect("load param");
    let frag = param.frag_intensity_model.clone().expect("model has no frag-intensity GBDT");
    eprintln!("model {model_id}: frag GBDT {} trees", frag.trees.len());
    let rt_seed = RtIndexModel::seed();

    let (spectra, errs) = input::read_spectra_parallel(&mzml, 2, 2, threads).unwrap();
    assert!(errs.is_empty(), "mzML errors");
    let by_scan: HashMap<i32, usize> =
        spectra.iter().enumerate().map(|(i, s)| (s.scan.expect("scan"), i)).collect();
    eprintln!("mzML {} MS2 spectra", spectra.len());

    let f = BufReader::new(File::open(&psms).unwrap());
    let mut lines = f.lines();
    let h: Vec<String> = lines.next().unwrap().unwrap().split('\t').map(String::from).collect();
    let (ci, cp) = (h.iter().position(|x| x == "SpecId").unwrap(), h.iter().position(|x| x == "peptide").unwrap());
    let mut w = BufWriter::new(File::create(&out).unwrap());
    writeln!(w, "SpecId\tcharge\tobs_rt_min\tn_ions\tgbdt_sa\tgbdt_cos\tgbdt_pearson\tflat_sa\tflat_cos\tflat_pearson\trt_seed_index\trt_features").unwrap();
    let mut cache = HashMap::new();
    let (mut n_ok, mut n_skip) = (0usize, 0usize);
    for line in lines {
        let line = line.unwrap();
        let x: Vec<&str> = line.split('\t').collect();
        let id = x[ci];
        let scan: i32 = id
            .split("scan=")
            .nth(1)
            .and_then(|r| r.split(|c: char| !c.is_ascii_digit()).next())
            .and_then(|s| s.parse().ok())
            .expect("scan in SpecId");
        let (Some((seq, dm)), Some(&si)) = (parse(x[cp]), by_scan.get(&scan)) else {
            n_skip += 1;
            continue;
        };
        let Some(pep) = peptide(&seq, &dm, &mut cache) else {
            n_skip += 1;
            continue;
        };
        let s = &spectra[si];
        let m = pep.mass();
        let charge = (1u8..=6).find(|&z| {
            (0..=3).any(|iso| {
                let mz = (m + iso as f64 * C13 + z as f64 * PROTON) / z as f64;
                ((mz - s.precursor_mz) / mz).abs() < 30e-6
            })
        });
        let Some(z) = charge else {
            n_skip += 1;
            continue;
        };
        // Observed b/y by max peak within tolerance; predicted from the GBDT.
        let ions = predict_by_ions(&pep, 1..=2);
        let pred = predict_frag_intensities(&frag, &pep, z);
        let (mut p, mut o) = (Vec::new(), Vec::new());
        for (ion, &pv) in ions.iter().zip(&pred) {
            if ion.charge > z || !matches!(ion.kind, IonKind::B | IonKind::Y) {
                continue;
            }
            let wd = if ppm { ion.mz * tol * 1e-6 } else { tol };
            let lo = s.peaks.partition_point(|(mz, _)| *mz < ion.mz - wd);
            let ob = s.peaks[lo..].iter().take_while(|(mz, _)| *mz <= ion.mz + wd).map(|(_, i)| *i as f64).fold(0.0, f64::max);
            p.push(pv);
            o.push(ob);
        }
        let g = sim(&p, &o);
        let fl = sim(&vec![1.0; p.len()], &o);
        let rf = rt_features(&pep);
        let rts = s.rt_seconds.map_or(f64::NAN, |t| t / 60.0);
        let rfs: Vec<String> = rf.iter().map(|v| format!("{v}")).collect();
        writeln!(
            w,
            "{id}\t{z}\t{rts:.4}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{}",
            p.len(),
            g[0], g[1], g[2], fl[0], fl[1], fl[2],
            rt_seed.predict_index(&rf),
            rfs.join(",")
        )
        .unwrap();
        n_ok += 1;
    }
    w.flush().unwrap();
    eprintln!("wrote {n_ok} rows, skipped {n_skip}");
}
