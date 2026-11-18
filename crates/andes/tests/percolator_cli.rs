use std::collections::BTreeSet;
use std::path::{Path, PathBuf};
use std::process::Command;

fn pin() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../percolator-native/tests/data/svm_ref.pin")
        .canonicalize()
        .expect("svm_ref.pin")
}

fn ids(path: &Path, column: usize) -> BTreeSet<String> {
    std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("{}: {e}", path.display()))
        .lines()
        .skip(1)
        .map(|l| l.split('\t').nth(column).unwrap_or("").to_string())
        .collect()
}

#[test]
fn percolator_subcommand_runs_the_native_engine() {
    let dir = tempfile::tempdir().expect("tempdir");
    let d = dir.path();
    let (t, dc, w) = (d.join("t.psms"), d.join("d.psms"), d.join("w.txt"));
    let out = Command::new(env!("CARGO_BIN_EXE_andes"))
        .arg("percolator")
        .args(["--rescore-memory-gb", "1", "--temp-dir"])
        .arg(d.join("scratch"))
        .args(["--seed", "42", "--only-psms", "-Y", "--num-threads", "2"])
        .arg("--results-psms")
        .arg(&t)
        .arg("--decoy-results-psms")
        .arg(&dc)
        .arg("--weights")
        .arg(&w)
        .arg(pin())
        .output()
        .expect("run andes percolator");
    assert!(
        out.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let header = std::fs::read_to_string(&t).unwrap();
    assert!(
        header.starts_with("PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds"),
        "{}",
        header.lines().next().unwrap_or("")
    );
    assert!(w.exists());
    assert!(d.join("scratch").is_dir());
    let mut written = ids(&t, 0);
    written.extend(ids(&dc, 0));
    assert_eq!(written, ids(&pin(), 0));
}

#[test]
fn percolator_subcommand_reports_engine_errors() {
    let dir = tempfile::tempdir().expect("tempdir");
    let d = dir.path();
    let out = Command::new(env!("CARGO_BIN_EXE_andes"))
        .args(["percolator", "--temp-dir"])
        .arg(d)
        .arg("--results-psms")
        .arg(d.join("t.psms"))
        .arg("--decoy-results-psms")
        .arg(d.join("d.psms"))
        .arg(pin())
        .output()
        .expect("run andes percolator");
    assert!(!out.status.success());
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("native training needs target-decoy competition: add -Y"),
        "stderr: {stderr}"
    );
    assert!(!d.join("t.psms").exists());
}

#[test]
fn percolator_subcommand_is_hidden_from_help() {
    for flag in ["-h", "--help"] {
        let out = Command::new(env!("CARGO_BIN_EXE_andes"))
            .arg(flag)
            .output()
            .expect("run andes --help");
        let help = String::from_utf8_lossy(&out.stdout);
        assert!(help.contains("inference"), "{help}");
        assert!(
            !help
                .lines()
                .any(|l| l.trim_start().starts_with("percolator")),
            "{help}"
        );
    }
}
