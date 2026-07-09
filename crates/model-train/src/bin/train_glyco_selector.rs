//! Train the learned glyco backbone SELECTOR (native GBDT) from an all-hits glyco
//! PIN + the reference engine truth, with a leakage-free held-out split.
//!
//! Usage:
//!   train_glyco_selector <all_hits.pin> <truth.tsv> <out_model.gbdt> [train_frac=0.6] [seed=42]
//!
//! The PIN must be produced with `ANDES_GLYCO_ALL_HITS=1` (features for every
//! candidate). Columns are read by the SAME names/order as
//! [`search::glyco_selector::glyco_selector_feature_names`] so the training vector is
//! positionally identical to what the collapse computes at inference (feature parity).
//!
//! Label = backbone-correct: |CalcMass − GlycanMass − (truth_backbone + H2O)| ≤ 0.05.
//! Split = by scan (scan % 5 < 3 → train 60%, else held-out 40%), so no scan appears in
//! both sets. Reports held-out top-1 recovery (the honest selector-quality number) and
//! writes the trained model bytes for the engine to load.

use std::collections::HashMap;
use std::io::{BufRead, BufReader, Write};

use model_train::gbdt::train::{train_gbdt, Dataset, TrainParams};
use search::glyco_selector::glyco_selector_feature_names;

const H2O: f64 = 18.010565;
const TOL: f64 = 0.05;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 4 {
        eprintln!("usage: train_glyco_selector <all_hits.pin> <truth.tsv> <out_model.gbdt> [train_frac=0.6] [seed=42]");
        std::process::exit(2);
    }
    let pin_path = &args[1];
    let truth_path = &args[2];
    let out_path = &args[3];
    let seed: u64 = args.get(5).and_then(|s| s.parse().ok()).unwrap_or(42);

    // ── truth: scan → neutral backbone mass ──
    let mut truth: HashMap<u32, f64> = HashMap::new();
    let tf = BufReader::new(std::fs::File::open(truth_path).expect("open truth"));
    for (i, line) in tf.lines().enumerate() {
        let line = line.unwrap();
        if i == 0 {
            continue;
        }
        let c: Vec<&str> = line.split('\t').collect();
        if c.len() < 2 {
            continue;
        }
        if let (Ok(scan), Ok(bb)) = (c[0].parse::<f64>(), c[1].parse::<f64>()) {
            truth.insert(scan as u32, bb + H2O);
        }
    }
    eprintln!("truth scans: {}", truth.len());

    // ── PIN header → column indices ──
    let feat_names = glyco_selector_feature_names();
    let pf = BufReader::new(std::fs::File::open(pin_path).expect("open pin"));
    let mut lines = pf.lines();
    let header = lines.next().unwrap().unwrap();
    let cols: Vec<&str> = header.split('\t').collect();
    let idx = |name: &str| -> usize {
        cols.iter().position(|c| *c == name).unwrap_or_else(|| panic!("PIN missing column {name}"))
    };
    let feat_idx: Vec<usize> = feat_names.iter().map(|n| idx(n)).collect();
    let (i_spec, i_calc, i_gly, i_label) = (idx("SpecId"), idx("CalcMass"), idx("GlycanMass"), idx("Label"));

    // ── read candidate rows (targets only) ──
    let nf = feat_names.len();
    // per split: x (row-major), y, groups(scan)
    let mut tr_x: Vec<f32> = Vec::new();
    let mut tr_y: Vec<u8> = Vec::new();
    let mut tr_g: Vec<u32> = Vec::new();
    // held-out: keep per-scan (features, label) to compute top-1
    let mut ho: HashMap<u32, Vec<(Vec<f32>, u8)>> = HashMap::new();
    let mut n_rows = 0u64;
    for line in lines {
        let line = line.unwrap();
        let p: Vec<&str> = line.split('\t').collect();
        if p.len() <= *feat_idx.iter().max().unwrap() {
            continue;
        }
        let label_raw: f64 = p[i_label].parse().unwrap_or(0.0);
        if label_raw < 0.0 {
            continue; // targets only
        }
        // scan from SpecId `..._glyco_<scan>_<row>`
        let scan = match p[i_spec].split("_glyco_").nth(1).and_then(|s| s.split('_').next()).and_then(|s| s.parse::<u32>().ok()) {
            Some(s) if truth.contains_key(&s) => s,
            _ => continue,
        };
        let calc: f64 = p[i_calc].parse().unwrap_or(0.0);
        let gly: f64 = p[i_gly].parse().unwrap_or(0.0);
        let y: u8 = if (calc - gly - truth[&scan]).abs() <= TOL { 1 } else { 0 };
        let fv: Vec<f32> = feat_idx.iter().map(|&i| p[i].parse::<f32>().unwrap_or(0.0)).collect();
        n_rows += 1;
        if scan % 5 < 3 {
            tr_x.extend_from_slice(&fv);
            tr_y.push(y);
            tr_g.push(scan);
        } else {
            ho.entry(scan).or_default().push((fv, y));
        }
    }
    let n_pos: u64 = tr_y.iter().map(|&v| v as u64).sum();
    eprintln!("rows={n_rows} train_rows={} train_pos={n_pos} heldout_scans={}", tr_y.len(), ho.len());

    // ── train native GBDT ──
    let ds = Dataset { x: tr_x, y: tr_y, groups: tr_g, n_features: nf };
    let params = TrainParams::default();
    let model = match train_gbdt(&ds, &params, seed) {
        Ok(m) => m,
        Err(e) => {
            eprintln!("train_gbdt failed: {e}");
            std::process::exit(1);
        }
    };

    // ── held-out top-1 recovery (honest selector quality) ──
    let (mut top1, mut tot) = (0u32, 0u32);
    for (_scan, cands) in &ho {
        if cands.is_empty() {
            continue;
        }
        tot += 1;
        let mut best_i = 0usize;
        let mut best_p = f32::NEG_INFINITY;
        for (i, (fv, _)) in cands.iter().enumerate() {
            let pr = model.predict_proba(fv);
            if pr > best_p {
                best_p = pr;
                best_i = i;
            }
        }
        if cands[best_i].1 == 1 {
            top1 += 1;
        }
    }
    eprintln!("HELD-OUT top-1 (native GBDT): {top1}/{tot} scans (40% split; gp fusion on the same held-out is the baseline to beat)");

    // ── write model bytes ──
    let bytes = model.to_bytes();
    let mut out = std::fs::File::create(out_path).expect("create out model");
    out.write_all(&bytes).expect("write model");
    eprintln!("wrote model: {out_path} ({} bytes, {nf} features)", bytes.len());
}
