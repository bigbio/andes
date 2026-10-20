//! `andes inference` over two runs' Percolator results: one dataset-level set
//! of protein groups in the QPX pg view.

use std::process::Command;

use arrow::array::{Array, BooleanArray, Float64Array, ListArray, StringArray};
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

const HDR: &str = "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds\n";

fn write(path: &std::path::Path, rows: &[(&str, f64, f64, &str, &str)]) {
    let mut text = HDR.to_string();
    for (id, score, q, pep, prots) in rows {
        text.push_str(&format!("{id}\t{score}\t{q}\t0.001\t{pep}\t{prots}\n"));
    }
    std::fs::write(path, text).unwrap();
}

#[test]
fn inference_writes_one_dataset_level_pg_view() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    std::fs::write(
        d.join("db.fasta"),
        ">sp|P1|A\nMAAAKCCCKR\n>sp|P2|B\nMDDDKEEEKR\n>sp|P3|C\nMGGGKR\n",
    )
    .unwrap();
    // Run 1 sees P1, run 2 sees P1 and P2; P3 only through a decoy-beaten PSM.
    // Both runs share SpecIds, as per-run files do.
    write(
        &d.join("r1.percolator.target.psms.txt"),
        &[("s_1_1", 5.0, 0.001, "K.CCCK.R", "sp|P1|A")],
    );
    write(
        &d.join("r1.percolator.decoy.psms.txt"),
        &[("s_2_1", -1.0, 0.5, "K.KGGG.M", "XXX_sp|P3|C")],
    );
    write(
        &d.join("r2.percolator.target.psms.txt"),
        &[
            ("s_1_1", 4.0, 0.001, "M.AAAK.C", "sp|P1|A"),
            ("s_2_1", 3.0, 0.001, "K.EEEK.R", "sp|P2|B"),
            ("s_3_1", -2.0, 0.6, "M.GGGK.R", "sp|P3|C"),
        ],
    );
    write(
        &d.join("r2.percolator.decoy.psms.txt"),
        &[("s_9_1", -3.0, 0.7, "K.KDDD.M", "XXX_sp|P2|B")],
    );
    let out = Command::new(env!("CARGO_BIN_EXE_andes"))
        .arg("inference")
        .arg("--psms")
        .arg(d.join("r1.percolator.target.psms.txt"))
        .arg(d.join("r2.percolator.target.psms.txt"))
        .arg("--decoy-psms")
        .arg(d.join("r1.percolator.decoy.psms.txt"))
        .arg(d.join("r2.percolator.decoy.psms.txt"))
        .arg("--database")
        .arg(d.join("db.fasta"))
        .arg("-o")
        .arg(d.join("qpx/exp"))
        .output()
        .expect("run andes inference");
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );

    let file = std::fs::File::open(d.join("qpx/exp.pg.parquet")).unwrap();
    let builder = ParquetRecordBatchReaderBuilder::try_new(file).unwrap();
    assert_eq!(builder.schema().metadata()["file_type"], "pg_file");
    let batch = builder.build().unwrap().next().unwrap().unwrap();
    let col = |n: &str| batch.column(batch.schema().index_of(n).unwrap()).clone();
    let anchor = col("anchor_protein");
    let anchor = anchor.as_any().downcast_ref::<StringArray>().unwrap();
    let decoy = col("is_decoy");
    let decoy = decoy.as_any().downcast_ref::<BooleanArray>().unwrap();
    let q = col("global_qvalue");
    let q = q.as_any().downcast_ref::<Float64Array>().unwrap();
    let rows: Vec<(String, bool, f64)> = (0..batch.num_rows())
        .map(|i| (anchor.value(i).to_string(), decoy.value(i), q.value(i)))
        .collect();
    // P1 has two peptides from two runs; P3 lost to its decoy.
    assert_eq!(rows[0].0, "sp|P1|A");
    assert!(!rows[0].1 && rows[0].2 == 0.0);
    assert!(rows.iter().any(|r| r.0 == "sp|P2|B" && !r.1));
    assert!(rows.iter().any(|r| r.0 == "sp|P3|C" && !r.1 && r.2 == 1.0));
    assert!(rows.iter().any(|r| r.0 == "XXX_sp|P3|C" && r.1));
    // Groups span both runs; two runs means no run-level q.
    let runs = col("grouped_runs");
    let runs = runs.as_any().downcast_ref::<ListArray>().unwrap();
    let names = runs.value(0);
    let names = names.as_any().downcast_ref::<StringArray>().unwrap();
    assert_eq!((names.value(0), names.value(1)), ("r1", "r2"));
    assert_eq!(col("pg_qvalue").null_count(), batch.num_rows());
}
