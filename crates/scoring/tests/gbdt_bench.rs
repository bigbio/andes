//! Micro-benchmark of the shipped GBDT evaluator for the compiled-ensemble comparison
//! (competitive plan, standard-search stop condition: "no hand-written SIMD/bitvector
//! tree traversal before the exact compiled ensemble benchmark").
//!
//! Ignored by default. Run with
//!   ANDES_GBDT_DIR=<dir from scripts/gbdt_bench/export_gbdt.py> \
//!   cargo test --release -p scoring --test gbdt_bench -- --ignored --nocapture
//! It loads `frag.agbd` and `rows.f32`, times `predict_value_batch` over all rows in
//! per-PSM-sized batches, prints ns/row and rows/s, and writes `preds_rust.f32` so the
//! compiled evaluators can be checked against it.
use scoring::gbdt_eval::GbdtPeakModel;
use std::time::Instant;

#[test]
#[ignore = "benchmark: needs ANDES_GBDT_DIR"]
fn bench_predict_value_batch() {
    let dir = match std::env::var("ANDES_GBDT_DIR") {
        Ok(d) => d,
        Err(_) => {
            eprintln!("ANDES_GBDT_DIR not set; skipping");
            return;
        }
    };
    let blob = std::fs::read(format!("{dir}/frag.agbd")).expect("frag.agbd");
    let model = GbdtPeakModel::from_bytes(&blob).expect("parse blob");
    let nf = model.n_features as usize;
    let raw = std::fs::read(format!("{dir}/rows.f32")).expect("rows.f32");
    let rows_flat: Vec<f32> = raw
        .chunks_exact(4)
        .map(|b| f32::from_le_bytes([b[0], b[1], b[2], b[3]]))
        .collect();
    let n = rows_flat.len() / nf;
    let rows: Vec<&[f32]> = (0..n).map(|i| &rows_flat[i * nf..(i + 1) * nf]).collect();
    eprintln!(
        "model: {} trees, {} features; rows: {n}",
        model.trees.len(),
        nf
    );

    let mut out = vec![0.0f32; n];
    // Per-PSM batch size ~ 4*(n-1) ions for a 13-mer at two charges: 48 rows.
    for &batch in &[48usize, 512, n] {
        let mut best = f64::MAX;
        for _ in 0..5 {
            let t = Instant::now();
            for (chunk_rows, chunk_out) in rows.chunks(batch).zip(out.chunks_mut(batch)) {
                model.predict_value_batch(chunk_rows, chunk_out);
            }
            best = best.min(t.elapsed().as_secs_f64());
        }
        let ns_row = best * 1e9 / n as f64;
        eprintln!(
            "predict_value_batch batch={batch:>5}: {:.1} ns/row  {:.2} Mrows/s  ({:.2} ns/row/tree)",
            ns_row,
            1e-6 / (best / n as f64),
            ns_row / model.trees.len() as f64
        );
    }
    // Reference predictions for the exactness check of the compiled evaluators.
    model.predict_value_batch(&rows, &mut out);
    let bytes: Vec<u8> = out.iter().flat_map(|v| v.to_le_bytes()).collect();
    std::fs::write(format!("{dir}/preds_rust.f32"), bytes).expect("write preds");
    let checksum: f64 = out.iter().map(|&v| v as f64).sum();
    eprintln!("preds_rust.f32 written; sum={checksum:.6}");
}
