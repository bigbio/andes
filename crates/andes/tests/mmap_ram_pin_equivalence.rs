//! Regression gate for GitHub issue #65: `--candidate-index mmap` and
//! `--candidate-index ram` must retrieve the SAME candidate multiset, in the
//! SAME order, and therefore emit the same PIN.
//!
//! Why this lives here and not next to the unit tests in
//! `crates/search/tests/mmap_candidate_backing.rs`: all four tests there passed
//! while the bug was live. Both defects it guards need REAL scale to surface —
//! a fixed modification on a real proteome (the `base_record_multiplicity`
//! fixed-mod mass offset) and proteins with an N-terminal Met whose two
//! enumeration passes both reach the precursor window (the candidate ORDER).
//! So the gate drives the real CLI over the in-repo `test.mgf.gz` +
//! `ecoli.fasta` fixtures.
//!
//! Three cases run with `--precursor-cal off` (before issue #70 was fixed the
//! calibration pre-pass silently downgraded mmap to in-RAM, which would have
//! made them vacuous); the fourth runs with calibration ON and is the gate for
//! #70: the pre-pass must run on the out-of-core index and the two backings
//! must still agree.
//!
//! Marked `#[ignore]` — it runs two full E. coli searches (~40 s release,
//! several minutes debug). Run it with:
//!
//! ```text
//! cargo test -p andes --release --test mmap_ram_pin_equivalence -- --ignored
//! ```

use std::path::PathBuf;
use std::process::Command;

fn workspace_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("canonicalize workspace root")
}

/// Run one search and return the PIN as (header, sorted rows).
///
/// Rows are sorted before comparison: the rayon pipeline has a known
/// row-ORDER nondeterminism across runs (see `precursor_cal_bit_identical`),
/// which is orthogonal to the candidate-set divergence under test. Both
/// defects this gate covers change row CONTENT (the retrieved multiset, and
/// the order-dependent `RawScoreCal` accumulator), so sorted comparison is
/// sufficient — and it was verified to FAIL on the unfixed code.
fn run_search(
    backing: &str,
    mods: Option<&str>,
    precursor_cal: &str,
    out: &PathBuf,
) -> (String, Vec<String>) {
    let root = workspace_root();
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_andes"));
    cmd.current_dir(&root)
        .arg("--spectrum")
        .arg("test-fixtures/test.mgf.gz")
        .arg("--database")
        .arg("test-fixtures/ecoli.fasta")
        .arg("--fragmentation")
        .arg("HCD")
        .arg("--fragment-tol-ppm")
        .arg("20")
        .arg("--threads")
        .arg("1")
        .arg("--precursor-cal")
        .arg(precursor_cal)
        .arg("--candidate-index")
        .arg(backing)
        .arg("--output-pin")
        .arg(out);
    if let Some(m) = mods {
        cmd.arg("--mods").arg(m);
    }
    let output = cmd.output().expect("spawn andes");
    assert!(
        output.status.success(),
        "andes --candidate-index {backing} failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    if backing == "mmap" {
        assert!(
            stderr.contains("out-of-core candidate-index: mmap"),
            "mmap run did not report the out-of-core backing (issue #70 downgrade?):\n{stderr}"
        );
        assert!(
            !stderr.contains("running in-RAM for this search"),
            "mmap was downgraded to in-RAM:\n{stderr}"
        );
    }

    let text = std::fs::read_to_string(out).expect("read pin");
    let mut lines = text.lines();
    let header = lines.next().expect("pin header").to_string();
    let mut rows: Vec<String> = lines.map(str::to_string).collect();
    rows.sort();
    (header, rows)
}

fn assert_backings_agree(label: &str, mods: Option<&str>, precursor_cal: &str) {
    let dir = tempfile::tempdir().expect("tempdir");
    let ram_pin = dir.path().join("ram.pin");
    let mmap_pin = dir.path().join("mmap.pin");

    let (ram_header, ram_rows) = run_search("ram", mods, precursor_cal, &ram_pin);
    let (mmap_header, mmap_rows) = run_search("mmap", mods, precursor_cal, &mmap_pin);

    assert_eq!(ram_header, mmap_header, "[{label}] PIN header differs");
    assert!(
        ram_rows.len() > 10_000,
        "[{label}] fixture produced only {} rows — the search did not run at \
         the intended scale, so this gate would be vacuous",
        ram_rows.len()
    );
    assert_eq!(
        ram_rows.len(),
        mmap_rows.len(),
        "[{label}] row COUNT differs: ram={} mmap={} — the two backings \
         retrieved different candidate sets (issue #65)",
        ram_rows.len(),
        mmap_rows.len()
    );
    if ram_rows != mmap_rows {
        let first = ram_rows
            .iter()
            .zip(&mmap_rows)
            .position(|(a, b)| a != b)
            .expect("row counts already asserted equal");
        panic!(
            "[{label}] PIN rows differ at sorted index {first}\n  ram : {}\n  mmap: {}",
            ram_rows[first], mmap_rows[first]
        );
    }
}

#[test]
#[ignore = "runs two full E. coli searches (~40 s release)"]
fn mmap_matches_ram_pin_with_default_mods() {
    assert_backings_agree("default mods", None, "off");
}

#[test]
#[ignore = "runs two full E. coli searches (~40 s release)"]
fn mmap_matches_ram_pin_with_fixed_mods_only() {
    // NumMods=0 with a fixed Cam-C: isolates the `base_record_multiplicity`
    // fixed-mod defect from any variable-mod expansion.
    let dir = tempfile::tempdir().expect("tempdir");
    let mods = dir.path().join("mods_fixed.txt");
    std::fs::write(&mods, "NumMods=0\n57.021464,C,fix,any,Carbamidomethyl\n").expect("write mods");
    assert_backings_agree("fixed mods only", Some(mods.to_str().unwrap()), "off");
}

#[test]
#[ignore = "runs two full E. coli searches (~40 s release)"]
fn mmap_matches_ram_pin_with_one_variable_mod() {
    let dir = tempfile::tempdir().expect("tempdir");
    let mods = dir.path().join("mods_var1.txt");
    std::fs::write(
        &mods,
        "NumMods=1\n57.021464,C,fix,any,Carbamidomethyl\n15.994915,M,opt,any,Oxidation\n",
    )
    .expect("write mods");
    assert_backings_agree("one variable mod", Some(mods.to_str().unwrap()), "off");
}

#[test]
#[ignore = "runs two full E. coli searches (~40 s release)"]
fn mmap_matches_ram_pin_with_precursor_calibration_on() {
    // Issue #70: with calibration on, the pre-pass used to build the in-RAM
    // index regardless of the backing (and OOM-killed searches whose index did
    // not fit). Now it runs on the chosen backing; both must still agree.
    assert_backings_agree("precursor-cal auto", None, "auto");
}
