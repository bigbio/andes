use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::Command;

use arrow::util::display::{ArrayFormatter, FormatOptions};

fn fixture(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(rel)
        .canonicalize()
        .unwrap_or_else(|e| panic!("canonicalize {rel}: {e}"))
}

fn rows_of(t: &output::ReadTable) -> Vec<HashMap<String, String>> {
    let opts = FormatOptions::default().with_null("");
    let mut out = Vec::new();
    for b in &t.batches {
        let names: Vec<String> = b
            .schema()
            .fields()
            .iter()
            .map(|f| f.name().clone())
            .collect();
        let fmts: Vec<ArrayFormatter> = b
            .columns()
            .iter()
            .map(|c| ArrayFormatter::try_new(c.as_ref(), &opts).unwrap())
            .collect();
        for i in 0..b.num_rows() {
            out.push(
                names
                    .iter()
                    .zip(&fmts)
                    .map(|(n, f)| (n.clone(), f.value(i).to_string()))
                    .collect(),
            );
        }
    }
    out
}

fn andes(args: &[&str]) -> String {
    let out = Command::new(env!("CARGO_BIN_EXE_andes"))
        .args(args)
        .output()
        .expect("run andes");
    let stderr = String::from_utf8_lossy(&out.stderr).into_owned();
    assert!(out.status.success(), "andes {args:?} failed:\n{stderr}");
    stderr
}

const MODS: &str = "Carbamidomethyl (C), TMT10plex (K), TMT10plex (N-term)";

fn search(spectrum: &Path, pin: &Path, qpx: &Path) {
    andes(&[
        "--spectrum",
        spectrum.to_str().unwrap(),
        "--database",
        fixture("test-fixtures/BSA.fasta").to_str().unwrap(),
        "--fixed-mods",
        MODS,
        "--variable-mods",
        "",
        "--max-mods",
        "1",
        "--output-pin",
        pin.to_str().unwrap(),
        "--output-parquet",
        qpx.to_str().unwrap(),
        "--tmt",
        "tmt10",
        "--tmt-min-purity",
        "0.9",
        "--threads",
        "2",
    ]);
}

const CHANNELS: [&str; 10] = [
    "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131",
];

fn sdrf(path: &Path, runs: &[(&str, &str)]) {
    let mut s = String::from(
        "source name\tcharacteristics[organism part]\tcharacteristics[biological replicate]\t\
         comment[instrument]\tcomment[label]\tcomment[fraction identifier]\tcomment[data file]\n",
    );
    for (file, fraction) in runs {
        for (i, c) in CHANNELS.iter().enumerate() {
            let (name, part, bio) = match i {
                0 | 9 => ("GIS", "pool", "pooled"),
                1 => ("x-blank", "blank", "1"),
                _ => ("", "liver", "1"),
            };
            let name = if name.is_empty() {
                format!("S{i}")
            } else {
                name.to_string()
            };
            s.push_str(&format!(
                "{name}\t{part}\t{bio}\tNT=Orbitrap Fusion;AC=MS:1002416\tTMT{c}\t{fraction}\t{file}\n"
            ));
        }
    }
    std::fs::write(path, s).unwrap();
}

#[test]
fn tmt_study_joins_pooled_scores_gates_groups_and_maps_channels() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let tables = d.join("tmt");
    std::fs::create_dir_all(&tables).unwrap();
    let mut run_tables = Vec::new();
    for run in ["runA", "runB"] {
        let spectrum = d.join(format!("{run}.mzML.gz"));
        std::fs::copy(fixture("test-fixtures/tmt10_synthetic.mzML.gz"), &spectrum).unwrap();
        search(
            &spectrum,
            &tables.join(format!("{run}.pin")),
            &d.join(format!("{run}_qpx")),
        );
        run_tables.push(tables.join(format!("{run}.tmt.parquet")));
    }
    let before: Vec<Vec<u8>> = run_tables
        .iter()
        .map(|p| std::fs::read(p).unwrap())
        .collect();

    let header = "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds\n";
    let mut targets = header.to_string();
    let mut expected: HashMap<String, (f64, f64)> = HashMap::new();
    let mut n_run_rows = 0;
    for (r, path) in run_tables.iter().enumerate() {
        let rows = rows_of(&output::read_tmt_parquet(path).unwrap());
        assert!(rows.len() >= 4, "{} rows in {}", rows.len(), path.display());
        let high_pep = rows
            .iter()
            .position(|row| row["quant_gate_passed"] == "true")
            .expect("a row passing the run's gate");
        for (i, row) in rows.iter().enumerate() {
            n_run_rows += 1;
            assert!(row["psm_q_value"].is_empty(), "searched without rescoring");
            let q = if r == 1 && i == 0 { 0.5 } else { 0.001 };
            let pep = if r == 0 && i == high_pep { 0.2 } else { 0.01 };
            let seq: String = {
                let mut out = String::new();
                let mut depth = 0;
                for c in row["peptide"].chars() {
                    match c {
                        '[' => depth += 1,
                        ']' => depth -= 1,
                        _ if depth == 0 => out.push(c),
                        _ => {}
                    }
                }
                out
            };
            targets.push_str(&format!(
                "{}\t{}\t{q}\t{pep}\tK.{seq}.A\t{}\n",
                row["spec_id"],
                5.0 - i as f64,
                row["proteins"].replace(';', "\t")
            ));
            if q <= 0.01 {
                expected.insert(row["spec_id"].clone(), (q, pep));
            }
        }
    }
    let t_path = d.join("pooled.t.psms");
    let d_path = d.join("pooled.d.psms");
    std::fs::write(&t_path, targets).unwrap();
    std::fs::write(
        &d_path,
        format!("{header}runA/decoy_1_1\t-3\t0.9\t1\tK.KAAAK.A\tXXX_sp|P02769|ALBU_BOVIN\n"),
    )
    .unwrap();
    let sdrf_path = d.join("study.sdrf.tsv");
    sdrf(
        &sdrf_path,
        &[
            ("runA.mzML.gz", "1"),
            ("runB.mzML.gz", "2"),
            ("runC.raw", "3"),
        ],
    );

    let qpx = d.join("study_qpx");
    let log = andes(&[
        "quant",
        "finalize",
        "--quant-dir",
        tables.to_str().unwrap(),
        "-o",
        d.join("out/study").to_str().unwrap(),
        "--database",
        fixture("test-fixtures/BSA.fasta").to_str().unwrap(),
        "--fixed-mods",
        MODS,
        "--variable-mods",
        "",
        "--psms",
        t_path.to_str().unwrap(),
        "--decoy-psms",
        d_path.to_str().unwrap(),
        "--protein-fdr",
        "1",
        "--sdrf",
        sdrf_path.to_str().unwrap(),
        "--output-parquet",
        qpx.to_str().unwrap(),
    ]);
    assert!(log.contains("TMT study"), "{log}");
    assert!(log.contains("1 of its 3 runs have no table"), "{log}");

    let after: Vec<Vec<u8>> = run_tables
        .iter()
        .map(|p| std::fs::read(p).unwrap())
        .collect();
    assert_eq!(before, after, "the per-run tables are not touched");

    let study = output::tmt_study_table::read_tmt_study(&d.join("out/study.tmt.parquet")).unwrap();
    assert_eq!(study.meta["andes.tmt.plex"], "tmt10plex");
    let rows = rows_of(&study);
    assert_eq!(rows.len(), expected.len());
    assert_eq!(rows.len(), n_run_rows - 1);
    let run_rows: HashMap<String, HashMap<String, String>> = run_tables
        .iter()
        .flat_map(|p| rows_of(&output::read_tmt_parquet(p).unwrap()))
        .map(|r| (r["spec_id"].clone(), r))
        .collect();
    let mut grouped = 0;
    for r in &rows {
        let (q, pep) = expected[&r["spec_id"]];
        let run = &run_rows[&r["spec_id"]];
        assert_eq!(r["psm_q_value"].parse::<f64>().unwrap(), q);
        assert_eq!(r["psm_pep"].parse::<f64>().unwrap(), pep);
        for c in CHANNELS {
            let k = format!("tmt_{c}");
            assert_eq!(r[&k], run[&k], "{k} is the run's value");
        }
        assert_eq!(r["run_gate_passed"], run["quant_gate_passed"]);
        let run_passed = run["quant_gate_passed"] == "true";
        let passed = run_passed && pep <= 0.05;
        assert_eq!(
            r["quant_gate_passed"],
            passed.to_string(),
            "{}",
            r["spec_id"]
        );
        let purity: f64 = r["purity"].parse().unwrap_or(1.0);
        let w: f64 = r["quant_weight"].parse().unwrap();
        assert!((w - if passed { purity * purity } else { 0.0 }).abs() < 1e-4);
        assert_eq!(r["plex"], "1");
        let fraction = if r["filename"] == "runA" { "1" } else { "2" };
        assert_eq!(r["fraction"], fraction);
        assert_eq!(r["instrument"], "Orbitrap Fusion");
        assert_eq!(r["sample_tmt_126"], "GIS");
        assert_eq!(r["role_tmt_126"], "reference");
        assert_eq!(r["role_tmt_131"], "reference");
        assert_eq!(r["role_tmt_127N"], "blank");
        assert_eq!(r["sample_tmt_128C"], "S4");
        assert_eq!(r["role_tmt_128C"], "sample");
        if !r["pg_accessions"].is_empty() {
            grouped += 1;
            assert!(
                r["pg_accessions"].contains("ALBU_BOVIN"),
                "{}",
                r["pg_accessions"]
            );
        }
    }
    assert!(grouped > 0);
    assert!(
        rows.iter().any(|r| r["run_gate_passed"] == "false"),
        "the interfered peptide fails the run's purity gate"
    );
    assert!(
        rows.iter()
            .any(|r| r["run_gate_passed"] == "true" && r["quant_gate_passed"] == "false"),
        "the PEP 0.2 row fails the study gate"
    );

    let read_feature = |p: &Path| -> Vec<HashMap<String, String>> {
        let f = std::fs::File::open(p).unwrap();
        let b = parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder::try_new(f).unwrap();
        let t = output::ReadTable {
            batches: b.build().unwrap().map(|b| b.unwrap()).collect(),
            meta: HashMap::new(),
        };
        rows_of(&t)
    };
    let features = read_feature(&qpx.join("quantms.feature.parquet"));
    assert_eq!(features.len(), grouped);
    let run_mz: HashMap<(String, String), f64> = ["runA", "runB"]
        .iter()
        .flat_map(|run| read_feature(&d.join(format!("{run}_qpx/quantms.feature.parquet"))))
        .map(|f| {
            (
                (f["peptidoform"].clone(), f["charge"].clone()),
                f["calculated_mz"].parse().unwrap(),
            )
        })
        .collect();
    for f in &features {
        let want = run_mz[&(f["peptidoform"].clone(), f["charge"].clone())];
        let got: f64 = f["calculated_mz"].parse().unwrap();
        assert!(
            (got - want).abs() < 1e-3,
            "{}: {got} vs {want}",
            f["peptidoform"]
        );
    }
    assert!(qpx.join("quantms.pg.parquet").is_file());
}

#[test]
fn a_rescored_run_table_is_refused() {
    let dir = tempfile::tempdir().unwrap();
    let d = dir.path();
    let tables = d.join("tmt");
    std::fs::create_dir_all(&tables).unwrap();
    let row = output::TmtRow {
        id: output::QuantId {
            file_idx: 0,
            run: "r".into(),
            spec_id: "r/s_1_1".into(),
            scan: 1,
            rt_seconds: None,
            sequence: "PEPTIDEK".into(),
            peptidoform: "PEPTIDEK".into(),
            modifications: Vec::new(),
            charge: 2,
            proteins: vec!["P1".into()],
            is_decoy: false,
            calculated_mz: 500.0,
            observed_mz: 500.0,
            q_value: Some(0.001),
            pep: Some(0.01),
            missed_cleavages: None,
            psm_index: 0,
        },
        quant_scan_id: "s".into(),
        quant_ms_level: 2,
        purity: None,
        raw: vec![1.0; 6],
        corrected: vec![1.0; 6],
        sn: None,
        quant_weight: 1.0,
        quant_gate_passed: true,
        interference: None,
        unmixed: None,
        transfer: None,
    };
    output::write_tmt_parquet(
        &tables.join("r.tmt.parquet"),
        quant::isobaric::Plex::Tmt6,
        &[row],
        false,
        None,
    )
    .unwrap();
    std::fs::write(d.join("db.fasta"), ">sp|P1|A\nMPEPTIDEKR\n").unwrap();
    let hdr = "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds\n";
    std::fs::write(
        d.join("t.psms"),
        format!("{hdr}r/s_1_1\t1\t0.001\t0.01\tK.PEPTIDEK.R\tsp|P1|A\n"),
    )
    .unwrap();
    std::fs::write(d.join("d.psms"), hdr).unwrap();
    let out = Command::new(env!("CARGO_BIN_EXE_andes"))
        .args(["quant", "finalize", "--quant-dir"])
        .arg(&tables)
        .arg("-o")
        .arg(d.join("study"))
        .arg("--database")
        .arg(d.join("db.fasta"))
        .arg("--psms")
        .arg(d.join("t.psms"))
        .arg("--decoy-psms")
        .arg(d.join("d.psms"))
        .output()
        .unwrap();
    assert!(!out.status.success());
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("carries a per-run q-value"), "{stderr}");
}
