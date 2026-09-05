//! Flag-ON guard for the glyco scoring-redesign columns.
//!
//! WHY THIS EXISTS. The goldens prove that every redesign flag is inert when OFF.
//! Nothing proved the flag-ON paths produce a column that VARIES: a flag-gated column
//! was once a constant and no test noticed, because the only automated check was
//! the OFF-path golden. This runs the fixture with every
//! redesign flag on and requires each new column to take at least two distinct values
//! across the rows. A column that is constant with its flag on is the silent-defect
//! shape this repo's path-parity guards exist for.
//!
//! It is NOT a golden: values are not pinned, only their non-constancy.

use std::collections::HashSet;
use std::path::PathBuf;
use std::process::Command;

fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .expect("workspace root")
        .to_path_buf()
}

#[test]
fn redesign_columns_vary_when_their_flags_are_on() {
    let root = workspace_root();
    let binary = PathBuf::from(env!("CARGO_BIN_EXE_andes"));
    let spectra = root.join("test-fixtures/glyco_fixture.mgf.gz");
    let fasta = root.join("test-fixtures/glyco_fixture.fasta");
    for f in [&spectra, &fasta] {
        assert!(f.exists(), "fixture missing: {}", f.display());
    }

    let outdir = tempfile::tempdir().expect("tempdir");
    let out_pin = outdir.path().join("flags.pin");
    let status = Command::new(&binary)
        .arg("--spectrum").arg(&spectra)
        .arg("--database").arg(&fasta)
        .arg("--glyco")
        .arg("--glyco-tol-ppm").arg("20")
        .arg("--fragmentation").arg("HCD")
        .arg("--glyco-taxon").arg("human")
        .arg("--glyco-y-tree")
        .arg("--glyco-oxonium-llr")
        .arg("--glyco-rank-masked")
        .arg("--glyco-chance-llr-masked")
        .arg("--output-pin").arg(&out_pin)
        .status()
        .expect("run andes");
    assert!(status.success(), "andes exited {status}");

    let actual = outdir.path().join("flags.glyco.pin");
    let text = std::fs::read_to_string(&actual)
        .unwrap_or_else(|e| panic!("read {}: {e}", actual.display()));
    let mut lines = text.lines().filter(|l| !l.trim().is_empty());
    let header: Vec<&str> = lines.next().expect("header").split('\t').collect();
    let rows: Vec<Vec<&str>> = lines.map(|l| l.split('\t').collect()).collect();
    assert!(rows.len() >= 20, "fixture produced only {} rows", rows.len());

    // Columns that must vary across rows with their flag on. `YTreeHighPriorMissing`
    // and `MaskedPeakCount` are integers and may legitimately be constant on a small
    // fixture, so they are checked for non-zero instead.
    let must_vary = [
        "YTreeLLR", "YTreeHitFrac", "YTreeDecoyGap", "OxoniumCompLLR",
        "RankScoreMasked", "ChanceLlrMasked", "ExplainedMasked",
    ];
    let mut failures = Vec::new();
    for col in must_vary {
        let Some(i) = header.iter().position(|h| *h == col) else {
            failures.push(format!("{col}: column missing"));
            continue;
        };
        let distinct: HashSet<&str> = rows.iter().map(|r| r[i]).collect();
        if distinct.len() < 2 {
            failures.push(format!("{col}: constant ({:?}) with its flag on", distinct));
        }
    }
    for col in ["MaskedPeakCount", "YTreeHighPriorMissing"] {
        let Some(i) = header.iter().position(|h| *h == col) else {
            failures.push(format!("{col}: column missing"));
            continue;
        };
        if rows.iter().all(|r| r[i] == "0" || r[i] == "0.0") {
            failures.push(format!("{col}: all zero with its flag on"));
        }
    }
    assert!(failures.is_empty(), "flag-on columns not live:\n  {}", failures.join("\n  "));
}
