use output::percolator::{parse_psm_results, PsmResultReader};
use rescore::{read_pin, run, write_results, MemLimit, Results};
use std::path::{Path, PathBuf};

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../rescore/tests/data")
        .join(name)
}

fn limit(dir: &Path) -> MemLimit {
    MemLimit {
        bytes: 1 << 30,
        spill_dir: dir.to_path_buf(),
    }
}

fn p(path: &Path) -> String {
    path.to_string_lossy().into_owned()
}

fn strs(v: &[&str]) -> Vec<String> {
    v.iter().map(|s| s.to_string()).collect()
}

struct Out {
    target: PathBuf,
    decoy: PathBuf,
    weights: PathBuf,
}

fn outs(dir: &Path, stem: &str) -> Out {
    Out {
        target: dir.join(format!("{stem}.target.psms")),
        decoy: dir.join(format!("{stem}.decoy.psms")),
        weights: dir.join(format!("{stem}.weights")),
    }
}

fn training_args(o: &Out, threads: &str, pin: &Path) -> Vec<String> {
    strs(&[
        "--seed",
        "42",
        "--results-psms",
        &p(&o.target),
        "--decoy-results-psms",
        &p(&o.decoy),
        "--weights",
        &p(&o.weights),
        "--only-psms",
        "-Y",
        "--num-threads",
        threads,
        &p(pin),
    ])
}

fn static_args(o: &Out, weights: &Path, tdc: bool, pin: &Path) -> Vec<String> {
    let mut a = strs(&[
        "--seed",
        "42",
        "--results-psms",
        &p(&o.target),
        "--decoy-results-psms",
        &p(&o.decoy),
        "--only-psms",
        "--static",
    ]);
    if tdc {
        a.push("--override".into());
    }
    a.push("--init-weights".into());
    a.push(p(weights));
    if tdc {
        a.push("-Y".into());
    }
    a.extend(strs(&["--num-threads", "2", &p(pin)]));
    a
}

fn row(id: &str, label: i8, file: &str, scan: u64, proteins: &str) -> String {
    format!(
        "{id}\t{label}\t{scan}\t{}\t{file}\t0\tK.PEP{scan}K.R\t{proteins}\n",
        1000.0 + scan as f64
    )
}

#[test]
fn results_parse_with_output_crate() {
    let triples = [
        ("a_1", 1, "run_a.mzML", 1, "P1\tP2"),
        ("b_2", 1, "run_b.mzML", 2, "P3\t\tP4"),
        ("a_3", -1, "run_a.mzML", 3, "decoy_P5"),
        ("b_4", 1, "run_b.mzML", 4, ""),
        ("a_5", -1, "run_a.mzML", 5, "decoy_P6\tdecoy_P7"),
    ];
    let dir = tempfile::tempdir().unwrap();
    let pin_path = dir.path().join("rows.pin");
    let mut pin_text =
        String::from("SpecId\tLabel\tScanNr\tExpMass\tFileName\tf1\tPeptide\tProteins\n");
    for &(id, label, file, scan, proteins) in &triples {
        pin_text.push_str(&row(id, label, file, scan, proteins));
    }
    std::fs::write(&pin_path, pin_text).unwrap();
    let pin = read_pin(&pin_path, limit(dir.path())).unwrap();
    let (meta, files) = (&pin.meta, &pin.file_names);
    assert_eq!(files, &["run_a.mzML", "run_b.mzML"]);
    let scores = [1.2345678, -0.000123456789, -2.5, 0.0, -1.0];
    let q = [0.0098039215686, 1.234567e-7, 0.5, 0.25, 1.0];
    let pep = [6.305116760146989e-16, 1.0, 0.75, 0.031, 1.0];
    let rows: Vec<usize> = (0..meta.len()).collect();
    let (t, d) = (dir.path().join("t.psms"), dir.path().join("d.psms"));
    let res = Results {
        rows: &rows,
        scores: &scores,
        q: &q,
        pep: &pep,
    };
    write_results(&t, &d, meta, &pin.texts, files, &res).unwrap();
    let text = std::fs::read_to_string(&t).unwrap();
    assert_eq!(
        text.lines().next().unwrap(),
        "PSMId\tfilename\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds"
    );
    assert_eq!(
        text.lines().nth(1).unwrap(),
        "a_1\trun_a.mzML\t1.23457\t0.00980392\t6.30512e-16\tK.PEP1K.R\tP1\tP2"
    );
    let ids: Vec<&str> = text
        .lines()
        .skip(1)
        .map(|l| l.split('\t').next().unwrap())
        .collect();
    assert_eq!(ids, ["a_1", "b_4", "b_2"]);
    let parsed = parse_psm_results(&text).unwrap();
    assert_eq!(parsed.len(), 3);
    for (i, _) in meta.iter().enumerate().filter(|(_, m)| m.label == 1) {
        let id = triples[i].0;
        let p = &parsed[id];
        let close = |a: f64, b: f64| (a - b).abs() <= 1e-5 * b.abs().max(1e-300);
        assert!(close(p.q_value, q[i]), "{} q {} vs {}", id, p.q_value, q[i]);
        assert!(close(p.pep, pep[i]), "{} pep {} vs {}", id, p.pep, pep[i]);
        let sc = p.score.unwrap();
        assert!(close(sc, scores[i]) || sc == scores[i], "{} score {sc}", id);
        assert_eq!(p.peptide, format!("K.PEP{}K.R", triples[i].3));
    }
    assert_eq!(parsed["a_1"].proteins, "P1;P2");
    assert_eq!(parsed["b_2"].proteins, "P3;P4");
    assert_eq!(parsed["b_4"].proteins, "");
    let decoys = std::fs::read_to_string(&d).unwrap();
    let mut reader = PsmResultReader::new(decoys.as_bytes()).unwrap();
    let mut seen = Vec::new();
    while let Some(r) = reader.next_row().unwrap() {
        seen.push((r.psm_id.to_string(), r.proteins(), r.pep));
    }
    assert_eq!(
        seen,
        [
            ("a_5".to_string(), "decoy_P6;decoy_P7".to_string(), 1.0),
            ("a_3".to_string(), "decoy_P5".to_string(), 0.75)
        ]
    );
}

#[test]
fn andes_argument_vectors_run() {
    let dir = tempfile::tempdir().unwrap();
    let pin = data("svm_ref.pin");
    let main = outs(dir.path(), "main");
    run(&training_args(&main, "2", &pin), limit(dir.path())).unwrap();
    for f in [&main.target, &main.decoy, &main.weights] {
        assert!(f.exists(), "{}", f.display());
    }
    let t = parse_psm_results(&std::fs::read_to_string(&main.target).unwrap()).unwrap();
    assert!(t.values().any(|psm| psm.q_value < 0.01));
    for tdc in [true, false] {
        let o = outs(dir.path(), if tdc { "run0" } else { "pep" });
        run(
            &static_args(&o, &main.weights, tdc, &pin),
            limit(dir.path()),
        )
        .unwrap();
        parse_psm_results(&std::fs::read_to_string(&o.target).unwrap()).unwrap();
        parse_psm_results(&std::fs::read_to_string(&o.decoy).unwrap()).unwrap();
        assert!(!o.weights.exists());
    }
}
