use std::path::PathBuf;
use std::process::Command;

fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .expect("workspace root")
        .to_path_buf()
}

fn read_sorted_rows(p: &std::path::Path) -> (String, Vec<String>) {
    let text = std::fs::read_to_string(p).unwrap_or_else(|e| panic!("read {}: {e}", p.display()));
    let mut lines = text.lines();
    let header = lines.next().unwrap_or_default().to_string();
    let mut rows: Vec<String> = lines
        .filter(|l| !l.trim().is_empty())
        .map(str::to_string)
        .collect();
    rows.sort_by(|x, y| spec_id(x).cmp(spec_id(y)));
    (header, rows)
}

fn row_diff(golden: &str, actual: &str) -> Option<String> {
    let g: Vec<&str> = golden.split('\t').collect();
    let a: Vec<&str> = actual.split('\t').collect();
    if g.len() != a.len() {
        return Some(format!("field count {} vs {}", g.len(), a.len()));
    }
    for (i, (gf, af)) in g.iter().zip(a.iter()).enumerate() {
        if gf == af {
            continue;
        }
        match (gf.parse::<f64>(), af.parse::<f64>()) {
            (Ok(gv), Ok(av)) => {
                if gv.is_nan() || av.is_nan() {
                    if gv.is_nan() != av.is_nan() {
                        return Some(format!("field {i}: {gf} vs {af} (NaN mismatch)"));
                    }
                    continue;
                }
                if gv.is_infinite() || av.is_infinite() {
                    if gv != av {
                        return Some(format!("field {i}: {gf} vs {af} (infinity mismatch)"));
                    }
                    continue;
                }
                let tol = 1e-5 * gv.abs().max(av.abs()).max(1.0);
                if (gv - av).abs() > tol {
                    return Some(format!("field {i}: {gf} vs {af}"));
                }
            }
            _ => return Some(format!("field {i}: {gf} vs {af}")),
        }
    }
    None
}

fn spec_id(row: &str) -> &str {
    row.split('\t').next().unwrap_or(row)
}

#[test]
fn glyco_pin_matches_golden_after_sort() {
    let root = workspace_root();
    let binary = PathBuf::from(env!("CARGO_BIN_EXE_andes"));
    let spectra = root.join("test-fixtures/glyco_fixture.mgf.gz");
    let fasta = root.join("test-fixtures/glyco_fixture.fasta");
    let gdb = root.join("test-fixtures/glyco_fixture.gdb");
    let golden = root.join("test-fixtures/parity/goldens/glyco.pin");
    for f in [&spectra, &fasta, &gdb, &golden] {
        assert!(f.exists(), "fixture missing: {}", f.display());
    }

    let outdir = tempfile::tempdir().expect("tempdir");
    let out_pin = outdir.path().join("actual.pin");
    let status = Command::new(&binary)
        .arg("--spectrum")
        .arg(&spectra)
        .arg("--database")
        .arg(&fasta)
        .arg("--glyco")
        .arg("--glyco-tol-ppm")
        .arg("20")
        .arg("--fragmentation")
        .arg("HCD")
        .arg("--fragment-tol-da")
        .arg("0.5")
        .arg("--model")
        .arg("cid_lowres_tryp")
        .arg("--glyco-glycan-gdb")
        .arg(&gdb)
        .arg("--output-pin")
        .arg(&out_pin)
        .status()
        .expect("run andes");
    assert!(status.success(), "andes exited {status}");

    let actual = outdir.path().join("actual.glyco.pin");
    assert!(
        actual.exists(),
        "glyco PIN not written: {}",
        actual.display()
    );

    let (g_hdr, g_rows) = read_sorted_rows(&golden);
    let (a_hdr, a_rows) = read_sorted_rows(&actual);

    assert_eq!(g_hdr, a_hdr, "glyco PIN column set changed");
    assert_eq!(
        g_rows.len(),
        a_rows.len(),
        "glyco PIN row count changed: golden {} vs actual {}",
        g_rows.len(),
        a_rows.len()
    );
    for (i, (g, a)) in g_rows.iter().zip(a_rows.iter()).enumerate() {
        if let Some(d) = row_diff(g, a) {
            panic!("glyco PIN row mismatch at sorted index {i} ({d})\n golden: {g}\n actual: {a}");
        }
    }
}

#[test]
fn glyco_precursor_mono_on_mgf_is_identical_to_off() {
    let root = workspace_root();
    let binary = PathBuf::from(env!("CARGO_BIN_EXE_andes"));
    let spectra = root.join("test-fixtures/glyco_fixture.mgf.gz");
    let fasta = root.join("test-fixtures/glyco_fixture.fasta");
    let gdb = root.join("test-fixtures/glyco_fixture.gdb");
    let outdir = tempfile::tempdir().expect("tempdir");
    let run = |name: &str, extra: &[&str]| -> (String, String) {
        let out_pin = outdir.path().join(format!("{name}.pin"));
        let mut cmd = Command::new(&binary);
        cmd.arg("--spectrum")
            .arg(&spectra)
            .arg("--database")
            .arg(&fasta)
            .args(["--glyco", "--glyco-tol-ppm", "20", "--fragmentation", "HCD"])
            .arg("--glyco-glycan-gdb")
            .arg(&gdb)
            .args(extra)
            .arg("--output-pin")
            .arg(&out_pin);
        let out = cmd.output().expect("run andes");
        assert!(out.status.success(), "andes exited {}", out.status);
        let glyco = outdir.path().join(format!("{name}.glyco.pin"));
        let pin = std::fs::read_to_string(&glyco)
            .unwrap_or_else(|e| panic!("read {}: {e}", glyco.display()));
        (pin, String::from_utf8_lossy(&out.stderr).into_owned())
    };
    let (off, _) = run("off", &[]);
    let (auto, log) = run("auto", &["--precursor-mono", "auto"]);
    assert!(off.lines().count() > 1, "fixture produced no glyco rows");
    assert_eq!(
        off, auto,
        "--precursor-mono auto on MGF must be byte-identical to off"
    );
    assert!(
        log.contains("--precursor-mono auto needs an input with MS1"),
        "expected the no-MS1 warning, got:\n{log}"
    );
}
