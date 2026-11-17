use std::path::PathBuf;
use std::process::Command;

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .expect("repo root")
        .to_path_buf()
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
fn glyco_highres_pin_matches_golden() {
    let root = repo_root();
    let spectra = root.join("test-fixtures/orbitrap_lumos_120.mzML.gz");
    let fasta = root.join("test-fixtures/glyco_fixture.fasta");
    let gdb = root.join("test-fixtures/glyco_fixture.gdb");
    let golden = root.join("test-fixtures/parity/goldens/glyco_highres.pin");
    for f in [&spectra, &fasta, &gdb, &golden] {
        assert!(f.exists(), "fixture missing: {}", f.display());
    }

    let tmpdir = tempfile::tempdir().expect("tempdir");
    let tmp = tmpdir.path();
    let out = tmp.join("out.pin");

    let binary = PathBuf::from(env!("CARGO_BIN_EXE_andes"));
    let result = Command::new(&binary)
        .arg("--spectrum")
        .arg(&spectra)
        .arg("--database")
        .arg(&fasta)
        .arg("--glyco")
        .arg("--glyco-tol-ppm")
        .arg("20")
        .arg("--glyco-glycan-gdb")
        .arg(&gdb)
        .arg("--output-pin")
        .arg(&out)
        .output()
        .expect("run andes");
    assert!(
        result.status.success(),
        "andes exited {}\n{}",
        result.status,
        String::from_utf8_lossy(&result.stderr)
    );

    let log = format!(
        "{}{}",
        String::from_utf8_lossy(&result.stdout),
        String::from_utf8_lossy(&result.stderr)
    );
    assert!(
        log.contains("hcd_qexactive_tryp"),
        "expected the high-res model to be auto-detected from the gzipped mzML, \
         but the run did not select hcd_qexactive_tryp. Log:\n{log}"
    );

    let actual = tmp.join("out.glyco.pin");
    assert!(
        actual.exists(),
        "glyco PIN not written: {}",
        actual.display()
    );

    let g_txt = std::fs::read_to_string(&golden).expect("read golden");
    let a_txt = std::fs::read_to_string(&actual).expect("read actual");
    let mut g_lines: Vec<&str> = g_txt.lines().collect();
    let mut a_lines: Vec<&str> = a_txt.lines().collect();

    let g_hdr = g_lines.remove(0);
    let a_hdr = a_lines.remove(0);
    assert_eq!(g_hdr, a_hdr, "glyco high-res PIN column set changed");
    assert_eq!(
        g_lines.len(),
        a_lines.len(),
        "glyco high-res PIN row count changed: golden {} vs actual {}",
        g_lines.len(),
        a_lines.len()
    );

    g_lines.sort_unstable_by_key(|r| spec_id(r));
    a_lines.sort_unstable_by_key(|r| spec_id(r));
    for (i, (g, a)) in g_lines.iter().zip(a_lines.iter()).enumerate() {
        if let Some(d) = row_diff(g, a) {
            panic!(
                "glyco PIN row mismatch at sorted index {i} ({d})
 golden: {g}
 actual: {a}"
            );
        }
    }
}

#[test]
fn glyco_highres_precursor_mono_without_ms1_is_identical_to_golden() {
    let root = repo_root();
    let spectra = root.join("test-fixtures/orbitrap_lumos_120.mzML.gz");
    let fasta = root.join("test-fixtures/glyco_fixture.fasta");
    let gdb = root.join("test-fixtures/glyco_fixture.gdb");
    let golden = root.join("test-fixtures/parity/goldens/glyco_highres.pin");
    let tmpdir = tempfile::tempdir().expect("tempdir");
    let out = tmpdir.path().join("out.pin");
    let dump = tmpdir.path().join("mono.tsv");
    let binary = PathBuf::from(env!("CARGO_BIN_EXE_andes"));
    let result = Command::new(&binary)
        .arg("--spectrum")
        .arg(&spectra)
        .arg("--database")
        .arg(&fasta)
        .args(["--glyco", "--glyco-tol-ppm", "20"])
        .arg("--glyco-glycan-gdb")
        .arg(&gdb)
        .args(["--precursor-mono", "auto", "--precursor-mono-dump"])
        .arg(&dump)
        .arg("--output-pin")
        .arg(&out)
        .output()
        .expect("run andes");
    let log = String::from_utf8_lossy(&result.stderr).into_owned();
    assert!(
        result.status.success(),
        "andes exited {}\n{log}",
        result.status
    );
    assert!(
        log.contains("precursor-mono: streamed 120 MS2 spectra (0 MS1 scans linked)"),
        "expected the MS1-linked stream summary, got:\n{log}"
    );
    assert!(
        log.contains("no MS2 could be fitted"),
        "expected the nothing-fitted warning, got:\n{log}"
    );

    let actual = std::fs::read_to_string(tmpdir.path().join("out.glyco.pin")).expect("read pin");
    let golden_txt = std::fs::read_to_string(&golden).expect("read golden");
    let a_hdr = actual.lines().next().unwrap();
    assert!(
        a_hdr.contains("MonoShift"),
        "Mono columns are written without a fitted MS2: {a_hdr}"
    );
    let hdr: Vec<&str> = a_hdr.split('\t').collect();
    let first = hdr.iter().position(|c| *c == "MonoShift").unwrap();
    for row in actual.lines().skip(1) {
        let v: Vec<&str> = row.split('\t').collect();
        assert_eq!(
            &v[first..first + 4],
            ["0", "0", "0", "0"],
            "unfitted Mono values"
        );
    }
    let mut a_rows: Vec<&str> = actual.lines().skip(1).collect();
    let mut g_rows: Vec<&str> = golden_txt.lines().skip(1).collect();
    a_rows.sort_by_key(|r| spec_id(r).to_string());
    g_rows.sort_by_key(|r| spec_id(r).to_string());
    assert_eq!(a_hdr, golden_txt.lines().next().unwrap());
    assert_eq!(
        a_rows.len(),
        g_rows.len(),
        "row count differs from the golden"
    );
    for (i, (g, a)) in g_rows.iter().zip(a_rows.iter()).enumerate() {
        if let Some(d) = row_diff(g, a) {
            panic!("row {i} differs from the golden ({d})\n golden: {g}\n actual: {a}");
        }
    }

    let dump_txt = std::fs::read_to_string(&dump).expect("read dump");
    assert_eq!(dump_txt.lines().count(), 121, "header + one row per MS2");
    assert!(
        dump_txt.lines().skip(1).all(|l| l.ends_with("\tNA")),
        "no MS2 has a linked MS1, so every row is NA"
    );
}
