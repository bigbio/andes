//! Guard: the STANDARD PIN's always-zero column set must not grow silently.
//!
//! Unlike the glyco PIN, this one does not prune its dead columns — downstream
//! consumers index PIN columns positionally, and the QPX writer shares
//! `psm_feature_values` with it, so removing a column is a breaking change rather
//! than hygiene. What is cheap, and what this test does, is pin the set so a NEW
//! always-zero feature cannot be added without someone noticing.
//!
//! Worth recording why pruning is not proposed here: the glyco prune was worth
//! 256.8 -> 384.6 glycoPSMs, but the benchmark record attributes that to a
//! data-starved classifier at a few hundred PSMs and states that at ~7,250 PSMs
//! "the pruning buys nothing". The standard search is in the second regime, so this
//! is a defect-class guard, not an expected identification win.
//!
//! THREE of these are structural — hardcoded 0.0 for every row in every mode,
//! including under `--chimeric`:
//!   IsolationWindowEfficiency, PrecursorIsotopeKL, PrecursorSNR
//! The rest are constant only because of what this fixture is: an MGF (so no
//! retention time and no linked MS1), searched with a low-resolution model (where
//! `EdgeScore` is documented as returning 0), with `--refine` off and no charge-5
//! precursors. Those are expected and are listed so the structural three stand out.
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

    // Anti-vacuity: without rows every column is trivially "constant" and this
    // test would pass while asserting nothing.
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
        // Structural: hardcoded 0.0 in every mode.
        "IsolationWindowEfficiency",
        "PrecursorIsotopeKL",
        "PrecursorSNR",
        // Fixture-dependent, listed so the three above stand out.
        "ScanNr",
        "charge5",
        "EdgeScore",
        "IsRefinement",
        "NumMods",
        "RefinementModClass",
        "DeltaRT",
        "AbsDeltaRT",
        "DeltaRTNorm",
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

    // Also assert the fixture exercises enough live features for the comparison
    // above to be meaningful.
    assert!(
        header.len() - constant.len() > 40,
        "expected >40 varying columns, got {}",
        header.len() - constant.len()
    );
}
