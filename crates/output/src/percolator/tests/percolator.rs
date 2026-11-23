use super::*;

fn psm(id: &str, q: f64, pep: f64) -> PercolatorPsm {
    PercolatorPsm {
        pooled_q_value: None,
        psm_id: id.to_string(),
        q_value: q,
        pep,
        score: None,
        peptide: "K.AAK.R".to_string(),
        proteins: "P".to_string(),
    }
}

#[test]
fn report_order_is_q_then_pep_then_id_whatever_the_map_order() {
    let rows: Vec<PercolatorPsm> = (0..300)
        .map(|i| {
            psm(
                &format!("run/scan={i}"),
                (i % 3) as f64 * 0.001,
                (i % 5) as f64 * 0.01,
            )
        })
        .collect();
    let sorted_ids = |map: &HashMap<String, PercolatorPsm>| -> Vec<String> {
        let mut v: Vec<&PercolatorPsm> = map.values().collect();
        v.sort_by(|a, b| a.report_cmp(b));
        v.iter().map(|p| p.psm_id.clone()).collect()
    };
    let first: HashMap<String, PercolatorPsm> =
        rows.iter().map(|p| (p.psm_id.clone(), p.clone())).collect();
    let want = sorted_ids(&first);
    for _ in 0..8 {
        let m: HashMap<String, PercolatorPsm> = rows
            .iter()
            .rev()
            .map(|p| (p.psm_id.clone(), p.clone()))
            .collect();
        assert_eq!(sorted_ids(&m), want);
    }
    let by_key: Vec<_> = want.iter().map(|id| &first[id]).collect();
    assert!(by_key.windows(2).all(|w| {
        (w[0].q_value, w[0].pep, w[0].psm_id.as_str())
            < (w[1].q_value, w[1].pep, w[1].psm_id.as_str())
    }));
}

#[test]
fn percolator_gets_the_thread_count_unless_the_user_set_it() {
    let user = vec!["--testFDR".to_string(), "0.05".to_string()];
    assert_eq!(
        percolator_args_with_threads(&user, 16),
        ["--num-threads", "16", "--testFDR", "0.05"]
    );
    let set = vec!["--num-threads".to_string(), "4".to_string()];
    assert_eq!(percolator_args_with_threads(&set, 16), set);
    let set_eq = vec!["--num-threads=4".to_string()];
    assert_eq!(percolator_args_with_threads(&set_eq, 16), set_eq);
    assert_eq!(percolator_args_with_threads(&[], 0), ["--num-threads", "1"]);
}

#[test]
fn percolator_threads_clamped_to_percolator_range() {
    assert_eq!(
        percolator_args_with_threads(&[], 192),
        ["--num-threads", "128"]
    );
    assert_eq!(
        percolator_args_with_threads(&[], 128),
        ["--num-threads", "128"]
    );
    let dir = tempfile::tempdir().unwrap();
    let native = PercolatorBackend::Native {
        memory: rescore::MemLimit {
            bytes: 1 << 30,
            spill_dir: dir.path().to_path_buf(),
            sizing: rescore::Sizing::Fixed,
        },
    };
    let pin = native_fixture_pin(dir.path());
    let targets =
        run_percolator_pooled(&native, &pin, &percolator_args_with_threads(&[], 500)).unwrap();
    assert!(!targets.is_empty());
}

#[test]
fn parse_pin_header() {
    let text = "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds\n\
                title_100_1\t2.5\t0.001\t0.0005\tR.PEPTIDEK.A\tsp|P1\tsp|P2\n\
                title_101_1\t1.0\t0.02\t0.03\tK.ANOTHERK.R\tsp|P3\n";
    let m = parse_psm_results(text).unwrap();
    assert_eq!(m.len(), 2);
    let a = &m["title_100_1"];
    assert_eq!(a.q_value, 0.001);
    assert_eq!(a.pep, 0.0005);
    assert_eq!(a.peptide, "R.PEPTIDEK.A");
    assert_eq!(a.proteins, "sp|P1;sp|P2");
    assert_eq!(m["title_101_1"].proteins, "sp|P3");
    assert_eq!(a.score, Some(2.5));
    assert_eq!(m["title_101_1"].score, Some(1.0));
}

#[test]
fn reader_splits_lines_as_str_lines() {
    let text = "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds\r\n\
                a_1_1\t2.5\t0.001\t0.0005\tR.PEPTIDEK.A\tsp|P1\tsp|P2\r\n\
                \r\n\
                b_1_1\t1.0\t0.02\t0.03\tK.ANOTHERK.R";
    let mut r = PsmResultReader::new(text.as_bytes()).unwrap();
    let a = r.next_row().unwrap().unwrap();
    assert_eq!(
        (a.psm_id, a.score, a.peptide),
        ("a_1_1", Some(2.5), "R.PEPTIDEK.A")
    );
    assert_eq!(a.proteins(), "sp|P1;sp|P2");
    assert_eq!(a.run, None);
    let b = r.next_row().unwrap().unwrap();
    assert_eq!((b.psm_id, b.q_value, b.pep), ("b_1_1", 0.02, 0.03));
    assert_eq!(b.proteins(), "");
    assert!(r.next_row().unwrap().is_none());
    let m = parse_psm_results(text).unwrap();
    assert_eq!(m.len(), 2);
    assert_eq!(m["a_1_1"].proteins, "sp|P1;sp|P2");
    assert!(matches!(
        PsmResultReader::new(&b""[..]),
        Err(PercolatorError::MissingColumn { .. })
    ));
}

#[test]
fn score_column_is_optional_and_lenient() {
    let text = "PSMId\tq-value\tposterior_error_prob\tpeptide\tproteinIds\n\
                a_1_1\t0.001\t0.0005\tR.PEPTIDEK.A\tsp|P1\n";
    assert_eq!(parse_psm_results(text).unwrap()["a_1_1"].score, None);
    let text = "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds\n\
                a_1_1\tnan\t0.001\t0.0005\tR.PEPTIDEK.A\tsp|P1\n\
                b_1_1\tinf\t0.002\t0.0005\tR.PEPTIDEK.A\tsp|P1\n\
                c_1_1\tx\t0.003\t0.0005\tR.PEPTIDEK.A\tsp|P1\n";
    let m = parse_psm_results(text).unwrap();
    assert_eq!(m.len(), 3);
    assert!(m.values().all(|p| p.score.is_none()));
    assert_eq!(m["a_1_1"].q_value, 0.001);
    let text = "PSMId\tq-value\tposterior_error_prob\tpeptide\tproteinIds\tscore\n\
                a_1_1\t0.001\t0.0005\tR.PEPTIDEK.A\n";
    assert_eq!(parse_psm_results(text).unwrap()["a_1_1"].score, None);
}

#[test]
fn pep_of_one_for_every_accepted_psm_is_degenerate() {
    let map = |rows: &[(f64, f64)]| -> HashMap<String, PercolatorPsm> {
        rows.iter()
            .enumerate()
            .map(|(i, (q, pep))| (format!("r/{i}"), psm(&format!("r/{i}"), *q, *pep)))
            .collect()
    };
    assert!(pep_is_degenerate(&map(&[
        (0.0002, 1.0),
        (0.004, 1.0),
        (0.2, 0.3)
    ])));
    assert!(!pep_is_degenerate(&map(&[(0.0002, 1e-20), (0.004, 1.0)])));
    assert!(!pep_is_degenerate(&map(&[(0.001, 0.0003), (0.3, 1.0)])));
    assert!(!pep_is_degenerate(&map(&[(0.3, 1.0)])));
    assert!(!pep_is_degenerate(&HashMap::new()));
}

#[test]
fn pep_pass_is_static_mix_max_without_tdc() {
    let args = static_pass_args(
        "/p/x.pin".into(),
        "/p/t".into(),
        "/p/d".into(),
        "/p/w".into(),
        &[
            "--num-threads".into(),
            "4".into(),
            "-Y".into(),
            "--post-processing-tdc".into(),
        ],
        false,
    );
    assert!(args.windows(2).any(|w| w == ["--init-weights", "/p/w"]));
    assert!(args.contains(&"--static".to_string()));
    assert!(!args
        .iter()
        .any(|a| a == "-Y" || a == "--post-processing-tdc"));
    assert!(args.windows(2).any(|w| w == ["--num-threads", "4"]));
    assert_eq!(args.last().unwrap(), "/p/x.pin");
}

#[test]
fn per_run_pass_is_static_with_tdc_and_keeps_the_weights() {
    let args = static_pass_args(
        "/p/x.pin".into(),
        "/p/t".into(),
        "/p/d".into(),
        "/p/w".into(),
        &["-Y".into(), "--num-threads".into(), "4".into()],
        true,
    );
    assert_eq!(args.iter().filter(|a| *a == "-Y").count(), 1);
    for flag in ["--static", "--override"] {
        assert!(args.contains(&flag.to_string()), "{flag}");
    }
    assert!(args.windows(2).any(|w| w == ["--init-weights", "/p/w"]));
    assert_eq!(args.last().unwrap(), "/p/x.pin");
}

#[test]
fn decoy_results_sit_next_to_the_pin() {
    assert_eq!(
        decoy_results_path(Path::new("/x/run.pin")),
        Path::new("/x/run.percolator.decoy.psms.txt")
    );
}

#[test]
fn parse_leading_filename_column() {
    let text =
        "PSMId\tFileName\tscore\tq-value\tposterior_error_probability\tpeptide\tproteinIds\n\
                sc_5_1\trun.mzML\t9.1\t0.004\t0.002\tK.SAGTPEPK.L\tDECOY_x\n";
    let m = parse_psm_results(text).unwrap();
    let r = &m["sc_5_1"];
    let mut rows = PsmResultReader::new(text.as_bytes()).unwrap();
    assert_eq!(rows.next_row().unwrap().unwrap().run, Some("run.mzML"));
    let lower = text.replacen("FileName", "filename", 1);
    let mut rows = PsmResultReader::new(lower.as_bytes()).unwrap();
    assert_eq!(rows.next_row().unwrap().unwrap().run, Some("run.mzML"));
    assert_eq!(r.q_value, 0.004);
    assert_eq!(r.pep, 0.002);
    assert_eq!(r.peptide, "K.SAGTPEPK.L");
    assert_eq!(r.proteins, "DECOY_x");
}

#[test]
fn clear_stale_removes_prior_outputs() {
    let dir = std::env::temp_dir().join(format!("andes_perc_stale_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let out = ResultPaths::new(&dir, "x");
    std::fs::write(&out.target, b"stale").unwrap();
    std::fs::write(&out.weights, b"stale").unwrap();
    assert!(out.target.exists());
    out.clear_stale().unwrap();
    assert!(!out.target.exists(), "stale target must be removed");
    assert!(!out.weights.exists(), "stale weights must be removed");
    out.clear_stale().unwrap();
    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn parse_errors_on_non_finite_qvalue() {
    let text = "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds\n\
                x_1_1\t2.0\tnan\t0.01\tK.PK.A\tsp|P1\n";
    let err = parse_psm_results(text).unwrap_err();
    assert!(
        matches!(err, PercolatorError::NonFiniteScore { column, .. } if column == "q-value"),
        "got {err:?}"
    );
}

#[test]
fn parse_errors_on_unparseable_pep() {
    let text = "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds\n\
                x_1_1\t2.0\t0.01\tnotanumber\tK.PK.A\tsp|P1\n";
    let err = parse_psm_results(text).unwrap_err();
    assert!(
        matches!(err, PercolatorError::NonFiniteScore { .. }),
        "got {err:?}"
    );
}

#[test]
fn parse_errors_on_short_row() {
    let text = "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds\n\
                x_1_1\t2.0\t0.01\n";
    let err = parse_psm_results(text).unwrap_err();
    assert!(
        matches!(err, PercolatorError::MalformedRow { .. }),
        "got {err:?}"
    );
}

#[test]
fn parse_errors_on_duplicate_psmid() {
    let text = "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds\n\
                dup_1_1\t2.0\t0.01\t0.005\tK.PK.A\tsp|P1\n\
                dup_1_1\t1.0\t0.02\t0.006\tK.QK.A\tsp|P2\n";
    let err = parse_psm_results(text).unwrap_err();
    assert!(
        matches!(err, PercolatorError::DuplicatePsmId { ref psm_id } if psm_id == "dup_1_1"),
        "got {err:?}"
    );
}

#[test]
fn parse_missing_qvalue_column_errors() {
    let text = "PSMId\tscore\tpeptide\tproteinIds\nx_1_1\t2.0\tK.PK.A\tsp|P1\n";
    let err = parse_psm_results(text).unwrap_err();
    match err {
        PercolatorError::MissingColumn { column, .. } => assert_eq!(column, "q-value"),
        other => panic!("expected MissingColumn, got {other:?}"),
    }
}

#[test]
fn resolve_explicit_missing_binary_errors() {
    let err = resolve_backend(
        Some(Path::new("/no/such/percolator")),
        false,
        DEFAULT_PERCOLATOR_IMAGE,
    )
    .unwrap_err();
    assert!(matches!(err, PercolatorError::BinaryNotFound(_)));
}

#[test]
fn resolve_explicit_existing_binary() {
    let me = std::env::current_exe().unwrap();
    let b = resolve_backend(Some(&me), false, DEFAULT_PERCOLATOR_IMAGE).unwrap();
    assert_eq!(b, PercolatorBackend::Binary(me));
}

#[test]
fn native_args_shape() {
    let out = ResultPaths::new(Path::new("/tmp/d"), "x");
    let args = native_args(
        Path::new("/tmp/d/x.pin"),
        &out,
        &["--testFDR".into(), "0.05".into()],
    );
    assert_eq!(args[0], "--seed");
    assert_eq!(args[1], "42");
    assert!(args.contains(&"--only-psms".to_string()));
    assert_eq!(args.iter().filter(|a| *a == "-Y").count(), 1);
    for user in ["-Y", "--post-processing-tdc"] {
        let args = native_args(Path::new("/tmp/d/x.pin"), &out, &[user.into()]);
        assert!(!args.iter().any(|a| a == "-Y" && user != "-Y"));
        assert_eq!(args.iter().filter(|a| *a == user).count(), 1);
    }
    assert_eq!(args.last().unwrap(), "/tmp/d/x.pin");
    let fdr_pos = args.iter().position(|a| a == "--testFDR").unwrap();
    assert!(fdr_pos < args.len() - 1);
    let seeded = native_args(Path::new("/tmp/d/x.pin"), &out, &["--seed".into(), "7".into()]);
    assert_eq!(seeded.iter().filter(|a| *a == "--seed").count(), 1);
    assert!(seeded.windows(2).any(|w| w == ["--seed", "7"]));
    let seeded = static_pass_args(
        "/p/x.pin".into(),
        "/p/t".into(),
        "/p/d".into(),
        "/p/w".into(),
        &["--seed".into(), "7".into()],
        true,
    );
    assert_eq!(seeded.iter().filter(|a| *a == "--seed").count(), 1);
    assert!(seeded.windows(2).any(|w| w == ["--seed", "7"]));
}

#[test]
fn docker_args_mounts_and_image() {
    let out = ResultPaths::new(Path::new("/data"), "run");
    let args = docker_args(
        DEFAULT_PERCOLATOR_IMAGE,
        Path::new("/data/run.pin"),
        Path::new("/data"),
        &out,
        &[],
    )
    .unwrap();
    assert_eq!(args[0], "run");
    assert!(args.contains(&"--platform".to_string()));
    let pin_mount = format!(
        "{}:/pin:ro",
        std::path::absolute("/data").unwrap().display()
    );
    assert!(args.contains(&pin_mount));
    assert!(args.contains(&DEFAULT_PERCOLATOR_IMAGE.to_string()));
    assert_eq!(args.last().unwrap(), "/pin/run.pin");
    assert!(args.contains(&"-Y".to_string()));
    assert!(args
        .iter()
        .any(|a| a == "/out/run.percolator.target.psms.txt"));
}

#[test]
fn docker_args_relative_pin_mounts_absolute_dirs() {
    let cwd = std::env::current_dir().unwrap();
    for (pin, dir) in [("out.pin", cwd.clone()), ("res/out.pin", cwd.join("res"))] {
        let pin = Path::new(pin);
        let out_dir = pin.parent().unwrap();
        let out = result_paths(pin);
        let args = docker_args(DEFAULT_PERCOLATOR_IMAGE, pin, out_dir, &out, &[]).unwrap();
        let mounts: Vec<&String> = args
            .iter()
            .zip(args.iter().skip(1))
            .filter(|(flag, _)| *flag == "-v")
            .map(|(_, m)| m)
            .collect();
        assert_eq!(
            mounts,
            [
                &format!("{}:/pin:ro", dir.display()),
                &format!("{}:/out", dir.display())
            ]
        );
    }
}

#[test]
fn unreadable_decoy_file_leaves_the_target_file_pooled() {
    let dir = tempfile::tempdir().unwrap();
    let paths = ResultPaths::new(dir.path(), "s");
    let target = "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds\n\
                  a\t1.0\t0.01\t0.001\tK.AAK.R\tP\n";
    std::fs::write(&paths.target, target).unwrap();
    std::fs::write(
        &paths.decoy,
        "PSMId\tscore\tpeptide\tproteinIds\nd\t0.1\tK.KAA.R\tXXX_P\n",
    )
    .unwrap();
    let target_q = HashMap::from([("a".to_string(), 0.5)]);
    let decoy_q = HashMap::from([("d".to_string(), 0.5)]);
    assert!(rewrite_result_q_values(&paths, &target_q, &decoy_q).is_err());
    assert_eq!(std::fs::read_to_string(&paths.target).unwrap(), target);

    std::fs::write(
        &paths.decoy,
        "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds\nd\t0.1\t0.2\t0.9\tK.KAA.R\tXXX_P\n",
    )
    .unwrap();
    rewrite_result_q_values(&paths, &target_q, &decoy_q).unwrap();
    for (path, pooled) in [(&paths.target, "0.01"), (&paths.decoy, "0.2")] {
        let text = std::fs::read_to_string(path).unwrap();
        let mut lines = text.lines();
        let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
        let row: Vec<&str> = lines.next().unwrap().split('\t').collect();
        let col = |name: &str| row[header.iter().position(|c| *c == name).unwrap()];
        assert_eq!(col("q-value").parse::<f64>().unwrap(), 0.5);
        assert_eq!(col(POOLED_Q_COLUMN), pooled);
    }
}

fn native_fixture_pin(dir: &Path) -> PathBuf {
    let src =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../rescore/tests/data/svm_ref.pin");
    let text = std::fs::read_to_string(src).unwrap();
    let mut lines = text.lines();
    let (id, rest) = lines.next().unwrap().split_once('\t').unwrap();
    let mut out = format!("{id}\t{RUN_COLUMN}\t{rest}\n");
    let rows: Vec<&str> = lines.collect();
    for run in 0..2u64 {
        for line in &rows {
            let c: Vec<&str> = line.split('\t').collect();
            let scan: u64 = c[2].parse().unwrap();
            out.push_str(&format!(
                "r{run}_{}\trun{run}.mzML\t{}\t{}\t{}\n",
                c[0],
                c[1],
                scan + 1000 * run,
                c[3..].join("\t")
            ));
        }
    }
    let pin = dir.join("fixture.pin");
    std::fs::write(&pin, out).unwrap();
    pin
}

#[test]
fn native_backend_runs_main_and_static_passes() {
    let dir = tempfile::tempdir().unwrap();
    let pin = native_fixture_pin(dir.path());
    let backend = PercolatorBackend::Native {
        memory: rescore::MemLimit {
            bytes: 1 << 30,
            spill_dir: dir.path().to_path_buf(),
            sizing: rescore::Sizing::Fixed,
        },
    };
    assert!(backend.describe().contains("native"));
    let extra = percolator_args_with_threads(&[], 2);
    let mut targets = run_percolator_pooled(&backend, &pin, &extra).unwrap();
    assert!(targets.values().any(|p| p.q_value < 0.01));
    let summary = apply_per_run_q(&backend, &pin, &extra, &mut targets)
        .unwrap()
        .unwrap();
    assert!(summary.contains("of 2 runs"), "{summary}");
    assert!(!summary.contains("kept the pooled"), "{summary}");
    assert!(targets.values().any(|p| p.pooled_q_value.is_some()));
    let out = result_paths(&pin);
    parse_psm_results(&std::fs::read_to_string(&out.target).unwrap()).unwrap();
    parse_psm_results(&std::fs::read_to_string(&out.decoy).unwrap()).unwrap();
    let (t, d) = score_all_rows(&backend, &pin, &extra).unwrap();
    assert_eq!(t.len() + d.len(), 402);
    let (map, decoys) = run_percolator_with_decoys(&backend, &pin, &extra).unwrap();
    assert!(!map.is_empty() && !decoys.is_empty());
}
