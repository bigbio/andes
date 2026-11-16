use output::percolator::{parse_psm_results, PsmResultReader};
use percolator_native::{
    fmt_fixed, fmt_g, merge_folds, normalize_fold_scores, read_pin, write_results, write_weights,
    CvResult, MemLimit, Results, RowMeta,
};
use std::path::PathBuf;

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

fn fold_data() -> (Vec<f64>, Vec<bool>) {
    let mut s = Vec::new();
    let mut d = Vec::new();
    for i in 0..200 {
        s.push(5.0 + i as f64 * 0.01);
        d.push(false);
    }
    for i in 0..100 {
        s.push(i as f64 * 0.03);
        d.push(true);
    }
    s.push(5.505);
    d.push(true);
    (s, d)
}

#[test]
fn fold_scaling_threshold_zero_median_decoy_minus_one() {
    let (s, d) = fold_data();
    let (sub, scale) = normalize_fold_scores(&s, &d, 0.01);
    let threshold = 5.0 + 51.0 * 0.01;
    let median_decoy = 51.0 * 0.03;
    assert_eq!(sub, threshold);
    assert_eq!(scale, threshold - median_decoy);
    assert_eq!((threshold - sub) / scale, 0.0);
    assert!(((median_decoy - sub) / scale + 1.0).abs() < 1e-12);
    let accepted_below = s
        .iter()
        .zip(&d)
        .filter(|&(&v, &dec)| !dec && (v - sub) / scale < 0.0)
        .count();
    assert_eq!(accepted_below, 51);
}

fn cv_of(copies: usize) -> (CvResult, Vec<bool>) {
    let (s, d) = fold_data();
    let mut test_scores = Vec::new();
    let mut fold = Vec::new();
    let mut decoy = Vec::new();
    for f in 0..copies {
        for (&v, &dec) in s.iter().zip(&d) {
            test_scores.push(v * (f + 1) as f64 + f as f64);
            fold.push(f as u8);
            decoy.push(dec);
        }
    }
    let w = |f: usize| vec![2.0 + f as f64, -1.0, 0.5 - f as f64];
    (
        CvResult {
            w_norm: [w(0), w(1), w(2)],
            w_raw: [
                vec![0.3, 0.7, -2.0],
                vec![0.2, 0.1, -1.0],
                vec![1.5, -0.5, 4.0],
            ],
            test_scores,
            fold,
        },
        decoy,
    )
}

#[test]
fn merge_rescales_scores_and_weights_per_fold() {
    let (mut cv, decoy) = cv_of(3);
    let before = cv.test_scores.clone();
    let (wn, wr) = (cv.w_norm.clone(), cv.w_raw.clone());
    merge_folds(&mut cv, &decoy, 0.01).unwrap();
    let n = before.len() / 3;
    for f in 0..3 {
        let rows: Vec<usize> = (f * n..(f + 1) * n).collect();
        let s: Vec<f64> = rows.iter().map(|&i| before[i]).collect();
        let d: Vec<bool> = rows.iter().map(|&i| decoy[i]).collect();
        let (sub, scale) = normalize_fold_scores(&s, &d, 0.01);
        for (&i, &v) in rows.iter().zip(&s) {
            assert_eq!(cv.test_scores[i], (v - sub) / scale);
        }
        for (got, orig) in [(&cv.w_norm[f], &wn[f]), (&cv.w_raw[f], &wr[f])] {
            assert_eq!(got[0], orig[0] / scale);
            assert_eq!(got[1], orig[1] / scale);
            assert_eq!(got[2], (orig[2] - sub) / scale);
        }
        let threshold = rows[51];
        assert_eq!(cv.test_scores[threshold], 0.0);
    }
}

#[test]
fn merge_aborts_when_median_decoy_is_accepted() {
    let mut test_scores: Vec<f64> = (0..5000).map(|i| i as f64).collect();
    let mut decoy = vec![false; 5000];
    for v in [2000.5, 1000.5, 0.5, -1.0] {
        test_scores.push(v);
        decoy.push(true);
    }
    let n = test_scores.len();
    let mut cv = CvResult {
        w_norm: [vec![1.0, 0.0], vec![1.0, 0.0], vec![1.0, 0.0]],
        w_raw: [vec![1.0, 0.0], vec![1.0, 0.0], vec![1.0, 0.0]],
        test_scores,
        fold: vec![0; n],
    };
    let (sub, scale) = normalize_fold_scores(&cv.test_scores, &decoy, 0.01);
    assert_eq!((sub, scale), (1000.5, 0.0));
    let err = merge_folds(&mut cv, &decoy, 0.01).unwrap_err();
    assert!(
        err.starts_with("Error: median decoy score <= score at 1% FDR."),
        "{err}"
    );
}

fn svm_ref_rows() -> (Vec<RowMeta>, Vec<String>) {
    let dir = tempfile::tempdir().unwrap();
    let pin = read_pin(
        &data("svm_ref.pin"),
        MemLimit {
            bytes: 1 << 30,
            spill_dir: dir.path().to_path_buf(),
        },
    )
    .unwrap();
    (pin.meta, pin.file_names)
}

#[test]
fn results_byte_identical_to_percolator() {
    let (meta, file_names) = svm_ref_rows();
    let mut by_id = std::collections::HashMap::new();
    for name in ["svm_ref_perc.target.psms", "svm_ref_perc.decoy.psms"] {
        for l in std::fs::read_to_string(data(name)).unwrap().lines().skip(1) {
            let c: Vec<&str> = l.split('\t').collect();
            let v: Vec<f64> = c[1..4].iter().map(|x| x.parse().unwrap()).collect();
            by_id.insert(c[0].to_string(), v);
        }
    }
    let rows: Vec<usize> = (0..meta.len()).collect();
    let col = |k: usize| -> Vec<f64> { meta.iter().map(|m| by_id[&m.psm_id][k]).collect() };
    let (scores, q, pep) = (col(0), col(1), col(2));
    let dir = tempfile::tempdir().unwrap();
    let (t, d) = (dir.path().join("t.psms"), dir.path().join("d.psms"));
    let res = Results {
        rows: &rows,
        scores: &scores,
        q: &q,
        pep: &pep,
    };
    write_results(&t, &d, &meta, &file_names, &res).unwrap();
    for (ours, theirs) in [
        (&t, "svm_ref_perc.target.psms"),
        (&d, "svm_ref_perc.decoy.psms"),
    ] {
        assert_eq!(
            std::fs::read_to_string(ours).unwrap(),
            std::fs::read_to_string(data(theirs)).unwrap(),
            "{theirs}"
        );
    }
}

fn row(id: &str, label: i8, file: u32, scan: u64, proteins: &str) -> RowMeta {
    RowMeta {
        psm_id: id.to_string(),
        label,
        file,
        scan,
        exp_mass: 1000.0 + scan as f64,
        peptide: format!("K.PEP{scan}K.R"),
        proteins: proteins.to_string(),
    }
}

#[test]
fn results_parse_with_output_crate() {
    let meta = [
        row("a_1", 1, 0, 1, "P1\tP2"),
        row("b_2", 1, 1, 2, "P3\t\tP4"),
        row("a_3", -1, 0, 3, "decoy_P5"),
        row("b_4", 1, 1, 4, ""),
        row("a_5", -1, 0, 5, "decoy_P6\tdecoy_P7"),
    ];
    let files = ["run_a.mzML".to_string(), "run_b.mzML".to_string()];
    let scores = [1.2345678, -0.000123456789, -2.5, 0.0, -1.0];
    let q = [0.0098039215686, 1.234567e-7, 0.5, 0.25, 1.0];
    let pep = [6.305116760146989e-16, 1.0, 0.75, 0.031, 1.0];
    let rows: Vec<usize> = (0..meta.len()).collect();
    let dir = tempfile::tempdir().unwrap();
    let (t, d) = (dir.path().join("t.psms"), dir.path().join("d.psms"));
    let res = Results {
        rows: &rows,
        scores: &scores,
        q: &q,
        pep: &pep,
    };
    write_results(&t, &d, &meta, &files, &res).unwrap();
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
    for (i, m) in meta.iter().enumerate().filter(|(_, m)| m.label == 1) {
        let p = &parsed[&m.psm_id];
        let close = |a: f64, b: f64| (a - b).abs() <= 1e-5 * b.abs().max(1e-300);
        assert!(
            close(p.q_value, q[i]),
            "{} q {} vs {}",
            m.psm_id,
            p.q_value,
            q[i]
        );
        assert!(
            close(p.pep, pep[i]),
            "{} pep {} vs {}",
            m.psm_id,
            p.pep,
            pep[i]
        );
        let sc = p.score.unwrap();
        assert!(
            close(sc, scores[i]) || sc == scores[i],
            "{} score {sc}",
            m.psm_id
        );
        assert_eq!(p.peptide, m.peptide);
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
fn number_format_matches_cpp_streams() {
    let g = [
        (1.0 / 102.0, "0.00980392"),
        (6.305116760146989e-16, "6.30512e-16"),
        (1.0, "1"),
        (-0.5647412, "-0.564741"),
        (0.0, "0"),
        (-0.0, "-0"),
        (1e-5, "1e-05"),
        (0.0001, "0.0001"),
        (123456789.0, "1.23457e+08"),
        (100000.0, "100000"),
        (1e6, "1e+06"),
        (9.999996, "10"),
        (2.5e-100, "2.5e-100"),
    ];
    for (v, want) in g {
        assert_eq!(fmt_g(v, 6), want, "{v:e}");
    }
    assert_eq!(fmt_g(0.39012, 3), "0.39");
    assert_eq!(fmt_g(-0.000123456, 3), "-0.000123");
    assert_eq!(fmt_fixed(0.35612, 3), "0.356");
    assert_eq!(fmt_fixed(-0.011849, 4), "-0.0118");
    assert_eq!(fmt_fixed(-0.00001, 4), "-0.0000");
}

fn percolator_read_weights(text: &str, n: usize) -> Vec<f64> {
    let mut lines = text.split_inclusive('\n');
    for line in lines.by_ref() {
        let first = line.as_bytes().first().copied().unwrap_or(b' ');
        if first == b'-' || first.is_ascii_digit() {
            break;
        }
    }
    let rest: String = lines.collect();
    rest.split_whitespace()
        .take(n)
        .map(|v| v.parse().unwrap())
        .collect()
}

#[test]
fn weights_file_readable_by_init_weights_reader() {
    let golden = std::fs::read_to_string(data("svm_ref_perc.weights")).unwrap();
    let lines: Vec<&str> = golden.lines().collect();
    let names: Vec<String> = lines[3]
        .split('\t')
        .take_while(|c| *c != "m0")
        .map(str::to_string)
        .collect();
    let nums = |l: &str| -> Vec<f64> { l.split('\t').map(|v| v.parse().unwrap()).collect() };
    let w_norm: Vec<Vec<f64>> = (0..3).map(|f| nums(lines[4 + 3 * f])).collect();
    let w_raw: Vec<Vec<f64>> = (0..3).map(|f| nums(lines[5 + 3 * f])).collect();
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("w.txt");
    write_weights(&path, &names, &w_norm, &w_raw).unwrap();
    let ours = std::fs::read_to_string(&path).unwrap();
    assert_eq!(ours, golden);

    let w_norm = [
        vec![2.1452191, -0.0648600, 0.2210942, 0.0150415],
        vec![-1.8136169, 0.0416757, -0.2749317, -0.0546565],
        vec![1.9548566, 0.0134639, 0.2468673, 0.0449189],
    ];
    let w_raw = [
        vec![1.5115123, -0.0555, 0.21184, -2.61768],
        vec![-1.2778, -0.03571, 0.26339, -2.41301],
        vec![1.37731, 0.011535, 0.23651, -2.38479],
    ];
    let names: Vec<String> = ["sep", "f2", "f4"].iter().map(|s| s.to_string()).collect();
    write_weights(&path, &names, &w_norm, &w_raw).unwrap();
    let text = std::fs::read_to_string(&path).unwrap();
    let read = percolator_read_weights(&text, names.len() + 1);
    let want: Vec<f64> = w_raw[0]
        .iter()
        .map(|&v| fmt_fixed(v, 4).parse().unwrap())
        .collect();
    assert_eq!(read, want);
    assert_eq!(
        text.lines().nth(4).unwrap(),
        "2.15\t-0.0649\t0.2211\t0.0150"
    );
    assert_eq!(
        text.lines().nth(7).unwrap(),
        "-1.814\t0.0417\t-0.2749\t-0.0547"
    );
}
