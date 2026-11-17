use std::collections::BTreeSet;
use std::process::Command;

fn repo_root() -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .expect("repo root")
        .to_path_buf()
}

#[test]
fn standard_pin_constant_column_set_is_pinned() {
    let root = repo_root();
    let dir = tempfile::tempdir().expect("tempdir");
    let pin = dir.path().join("out.pin");

    let status = Command::new(env!("CARGO_BIN_EXE_andes"))
        .current_dir(&root)
        .arg("--spectrum")
        .arg("test-fixtures/test.mgf.gz")
        .arg("--database")
        .arg("test-fixtures/BSA.fasta")
        .arg("--output-pin")
        .arg(&pin)
        .status()
        .expect("run andes");
    assert!(status.success(), "fixture search should succeed: {status}");

    let text = std::fs::read_to_string(&pin).expect("read pin");
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().expect("header").split('\t').collect();
    let rows: Vec<Vec<&str>> = lines
        .map(|l| l.split('\t').collect::<Vec<_>>())
        .filter(|r| r.len() == header.len())
        .collect();

    assert!(
        rows.len() > 100,
        "fixture must produce many rows for constancy to mean anything, got {}",
        rows.len()
    );

    let constant: BTreeSet<&str> = header
        .iter()
        .enumerate()
        .filter(|(i, _)| {
            let first = rows[0][*i];
            rows.iter().all(|r| r[*i] == first)
        })
        .map(|(_, name)| *name)
        .collect();

    let expected: BTreeSet<&str> = [
        "IsolationWindowEfficiency",
        "PrecursorIsotopeKL",
        "PrecursorSNR",
        "EdgeScore",
        "IsRefinement",
        "NumMods",
        "RefinementModClass",
        "DeltaRT",
        "AbsDeltaRT",
        "DeltaRTNorm",
        "andesrt_rt_diff",
        "andesrt_predicted_rt",
        "andesrt_missing",
        "FileName",
        "IsVariant",
        "IsIsoform",
        "IsChain",
        "RefCompetitor",
        "DeltaRefScore",
        "IsotopeMimic",
    ]
    .into_iter()
    .collect();

    let unexpected: Vec<&&str> = constant.difference(&expected).collect();
    assert!(
        unexpected.is_empty(),
        "new always-constant PIN column(s) {unexpected:?} — a feature that never \
         varies cannot help Percolator separate anything. Either populate it or \
         drop it; if it is genuinely fixture-dependent, add it here with a reason."
    );
    let now_varying: Vec<&&str> = expected.difference(&constant).collect();
    assert!(
        now_varying.is_empty(),
        "column(s) {now_varying:?} are no longer constant — good, but remove them \
         from this guard's list so it keeps its teeth"
    );

    assert!(
        header.len() - constant.len() > 40,
        "expected >40 varying columns, got {}",
        header.len() - constant.len()
    );
}
