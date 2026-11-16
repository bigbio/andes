use percolator_native::{fmt_g, mix_max_q, run, tdc_q, MemLimit};
use std::collections::HashMap;
use std::path::{Path, PathBuf};

fn data(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
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

struct PinRow {
    id: String,
    label: i8,
    scan: u64,
    exp_mass: f64,
    features: Vec<f32>,
}

fn pin_rows(path: &Path) -> Vec<PinRow> {
    let text = std::fs::read_to_string(path).unwrap();
    text.lines()
        .skip(1)
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            PinRow {
                id: c[0].to_string(),
                label: c[1].parse().unwrap(),
                scan: c[2].parse().unwrap(),
                exp_mass: c[3].parse().unwrap(),
                features: c[4..9]
                    .iter()
                    .map(|v| v.parse::<f64>().unwrap() as f32)
                    .collect(),
            }
        })
        .collect()
}

fn subset_pin(dir: &Path, name: &str, keep: impl Fn(i8) -> bool) -> PathBuf {
    let text = std::fs::read_to_string(data("svm_ref.pin")).unwrap();
    let mut lines = text.lines();
    let mut out = String::from(lines.next().unwrap());
    out.push('\n');
    for l in lines {
        let label: i8 = l.split('\t').nth(1).unwrap().parse().unwrap();
        if keep(label) {
            out.push_str(l);
            out.push('\n');
        }
    }
    let path = dir.join(name);
    std::fs::write(&path, out).unwrap();
    path
}

fn shown(v: f64) -> f64 {
    fmt_g(v, 6).parse().unwrap()
}

fn raw_bin1() -> Vec<f64> {
    vec![0.2747, -0.0101, -0.0074, 0.0385, 0.0111, -1.1262]
}

fn dot(features: &[f32], w: &[f64]) -> f64 {
    let nf = features.len();
    let mut s = w[nf];
    for (&x, &wj) in features.iter().zip(&w[..nf]).rev() {
        s += f64::from(x) * wj;
    }
    s
}

fn scores(path: &Path) -> HashMap<String, f64> {
    let text = std::fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let col = header.iter().position(|c| *c == "score").unwrap();
    lines
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            (c[0].to_string(), c[col].parse().unwrap())
        })
        .collect()
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

#[test]
fn static_pass_equals_dot_product() {
    let dir = tempfile::tempdir().unwrap();
    let pin = data("svm_ref.pin");
    let o = outs(dir.path(), "all");
    let mut a = static_args(&o, &data("svm_ref_perc.weights"), false, &pin);
    a.insert(0, "--override".into());
    run(&a, limit(dir.path())).unwrap();
    let mut got = scores(&o.target);
    got.extend(scores(&o.decoy));
    let rows = pin_rows(&pin);
    assert_eq!(got.len(), rows.len());
    let w = raw_bin1();
    for r in &rows {
        let want = dot(&r.features, &w);
        let have = got[&r.id];
        assert!(
            (have - want).abs() <= 1e-5 * want.abs().max(1.0),
            "{}: {have} vs {want}",
            r.id
        );
    }

    let o = outs(dir.path(), "tdc");
    run(
        &static_args(&o, &data("svm_ref_perc.weights"), true, &pin),
        limit(dir.path()),
    )
    .unwrap();
    let mut got = scores(&o.target);
    got.extend(scores(&o.decoy));
    let mut best: HashMap<(u64, u64), f64> = HashMap::new();
    for r in &rows {
        let s = dot(&r.features, &w);
        let e = best
            .entry((r.scan, r.exp_mass.to_bits()))
            .or_insert(f64::NEG_INFINITY);
        *e = e.max(s);
    }
    assert_eq!(got.len(), best.len());
    for r in &rows {
        if let Some(&have) = got.get(&r.id) {
            let want = best[&(r.scan, r.exp_mass.to_bits())];
            assert!((have - want).abs() <= 1e-5 * want.abs().max(1.0));
        }
    }
    let labels: HashMap<&str, i8> = rows.iter().map(|r| (r.id.as_str(), r.label)).collect();
    assert!(scores(&o.target).keys().all(|id| labels[id.as_str()] == 1));
    assert!(scores(&o.decoy).keys().all(|id| labels[id.as_str()] == -1));
}

#[test]
fn static_tdc_run_without_decoys_aborts_as_percolator() {
    let dir = tempfile::tempdir().unwrap();
    let pin = subset_pin(dir.path(), "targets.pin", |l| l == 1);
    let o = outs(dir.path(), "t");
    let err = run(
        &static_args(&o, &data("svm_ref_perc.weights"), true, &pin),
        limit(dir.path()),
    )
    .unwrap_err();
    assert!(err.contains("Error: no decoy PSMs were provided."), "{err}");
    assert!(!o.target.exists());
}

#[test]
fn static_run_without_targets_aborts_as_percolator() {
    let dir = tempfile::tempdir().unwrap();
    let pin = subset_pin(dir.path(), "decoys.pin", |l| l == -1);
    for tdc in [true, false] {
        let o = outs(dir.path(), if tdc { "tdc" } else { "mm" });
        let err = run(
            &static_args(&o, &data("svm_ref_perc.weights"), tdc, &pin),
            limit(dir.path()),
        )
        .unwrap_err();
        assert!(
            err.contains("Error: no target PSMs were provided."),
            "{err}"
        );
        assert!(!o.target.exists());
    }
}

#[test]
fn unknown_argument_is_an_error() {
    let dir = tempfile::tempdir().unwrap();
    let o = outs(dir.path(), "u");
    let pin = p(&data("svm_ref.pin"));
    for bad in [
        "--trainFDR",
        "-y",
        "--subset-max-train",
        "-v",
        "--post-processing-mix-max",
    ] {
        let mut a = training_args(&o, "1", &data("svm_ref.pin"));
        a.insert(0, bad.to_string());
        let err = run(&a, limit(dir.path())).unwrap_err();
        assert!(err.contains(bad), "{bad}: {err}");
    }
    let mut a = training_args(&o, "1", &data("svm_ref.pin"));
    a.insert(0, "extra.pin".to_string());
    let err = run(&a, limit(dir.path())).unwrap_err();
    assert!(err.contains("extra.pin"), "{err}");
    for (flag, value) in [("--seed", "0"), ("--num-threads", "0"), ("--seed", "x")] {
        let a = strs(&[
            flag,
            value,
            "--results-psms",
            "t",
            "--decoy-results-psms",
            "d",
            &pin,
        ]);
        let err = run(&a, limit(dir.path())).unwrap_err();
        assert!(err.contains(flag), "{flag}: {err}");
    }
    let a = strs(&["--results-psms", "t", "--decoy-results-psms", "d", "--seed"]);
    assert!(run(&a, limit(dir.path())).unwrap_err().contains("--seed"));
    let a = strs(&[
        "--results-psms",
        "t",
        "--decoy-results-psms",
        "d",
        "--static",
        &pin,
    ]);
    assert!(run(&a, limit(dir.path()))
        .unwrap_err()
        .contains("--init-weights"));
    assert!(!o.target.exists());
}

#[test]
fn mixmax_q_monotone() {
    let mut s = Vec::new();
    let mut d = Vec::new();
    let mut x: u64 = 7;
    for i in 0..600 {
        x = x
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let u = (x >> 11) as f64 / (1u64 << 53) as f64;
        let decoy = i % 2 == 1;
        let shift = if decoy || i % 4 == 0 { 0.0 } else { 2.0 };
        s.push(((u * 4.0 + shift) * 8.0).round() / 8.0);
        d.push(decoy);
    }
    for pi0 in [0.3, 0.6, 0.95] {
        let q = mix_max_q(&s, &d, pi0);
        let mut order: Vec<usize> = (0..s.len()).collect();
        order.sort_by(|&a, &b| s[b].total_cmp(&s[a]));
        for w in order.windows(2) {
            assert!(q[w[0]] <= q[w[1]], "pi0 {pi0}");
            if s[w[0]] == s[w[1]] {
                assert_eq!(q[w[0]], q[w[1]]);
            }
        }
        assert!(q.iter().all(|&v| (0.0..=1.0).contains(&v)));
        let plain = tdc_q(&s, &d, false);
        assert!(q.iter().zip(&plain).all(|(a, b)| a <= b));
    }
    assert_eq!(mix_max_q(&s, &d, 1.0), tdc_q(&s, &d, false));
}

#[test]
fn training_output_independent_of_threads() {
    let dir = tempfile::tempdir().unwrap();
    let pin = data("svm_ref.pin");
    let a = outs(dir.path(), "a");
    let b = outs(dir.path(), "b");
    run(&training_args(&a, "1", &pin), limit(dir.path())).unwrap();
    run(&training_args(&b, "4", &pin), limit(dir.path())).unwrap();
    for (x, y) in [
        (&a.target, &b.target),
        (&a.decoy, &b.decoy),
        (&a.weights, &b.weights),
    ] {
        assert_eq!(
            std::fs::read_to_string(x).unwrap(),
            std::fs::read_to_string(y).unwrap()
        );
    }
}

#[test]
fn static_override_keeps_weights_that_separate_nothing() {
    let dir = tempfile::tempdir().unwrap();
    let pin = data("svm_ref.pin");
    let w = dir.path().join("flat.weights");
    std::fs::write(
        &w,
        "# header\nsep\tf2\tf3\tf4\tf5\tm0\n0\t0\t0\t0\t0\t0\n0\t0\t0\t0\t0\t0.5\n",
    )
    .unwrap();
    let o = outs(dir.path(), "keep");
    run(&static_args(&o, &w, true, &pin), limit(dir.path())).unwrap();
    let mut got = scores(&o.target);
    got.extend(scores(&o.decoy));
    assert!(got.values().all(|&s| s == 0.5));

    let o = outs(dir.path(), "reset");
    let mut a = static_args(&o, &w, true, &pin);
    a.retain(|x| x != "--override");
    run(&a, limit(dir.path())).unwrap();
    let mut got = scores(&o.target);
    got.extend(scores(&o.decoy));
    let rows = pin_rows(&pin);
    let by_id: HashMap<&str, &PinRow> = rows.iter().map(|r| (r.id.as_str(), r)).collect();
    let single = (0..5).any(|j| {
        [1.0f64, -1.0].iter().any(|&sign| {
            got.iter()
                .all(|(id, &s)| s == shown(sign * f64::from(by_id[id.as_str()].features[j]) + 0.5))
        })
    });
    assert!(
        single,
        "reset scores are not one raw feature plus the kept bias"
    );
}

#[test]
fn missing_weights_values_read_as_zero() {
    let dir = tempfile::tempdir().unwrap();
    let pin = data("svm_ref.pin");
    let w = dir.path().join("short.weights");
    std::fs::write(&w, "names\n1\t2\n0\t0\t0\n").unwrap();
    let o = outs(dir.path(), "short");
    let mut a = static_args(&o, &w, false, &pin);
    a.insert(0, "--override".into());
    run(&a, limit(dir.path())).unwrap();
    let mut got = scores(&o.target);
    got.extend(scores(&o.decoy));
    assert!(got.values().all(|&s| s == 0.0));
}

#[test]
fn training_without_tdc_is_an_error_on_concatenated_input() {
    let dir = tempfile::tempdir().unwrap();
    let o = outs(dir.path(), "plain");
    let mut a = training_args(&o, "1", &data("svm_ref.pin"));
    a.retain(|x| x != "-Y");
    let err = run(&a, limit(dir.path())).unwrap_err();
    assert!(err.contains("-Y"), "{err}");
    assert!(!o.target.exists());
}

#[test]
fn unreadable_init_weights_without_direction_aborts_even_with_override() {
    let dir = tempfile::tempdir().unwrap();
    let pin = dir.path().join("flat.pin");
    let mut text = String::from("SpecId\tLabel\tScanNr\tExpMass\tf1\tf2\tPeptide\tProteins\n");
    for i in 0..60 {
        let label = if i % 2 == 0 { 1 } else { -1 };
        text.push_str(&format!(
            "r{i}\t{label}\t{i}\t{}\t1\t2\tK.A.R\tP\n",
            1000 + i
        ));
    }
    std::fs::write(&pin, text).unwrap();
    let missing = dir.path().join("missing.weights");
    let o = outs(dir.path(), "miss");
    for tdc in [true, false] {
        let mut a = static_args(&o, &missing, tdc, &pin);
        if !tdc {
            a.insert(0, "--override".into());
        }
        let err = run(&a, limit(dir.path())).unwrap_err();
        assert!(err.contains("cannot find an initial direction"), "{err}");
        assert!(!o.target.exists());
    }
}
