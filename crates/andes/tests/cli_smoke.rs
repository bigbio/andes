//! End-to-end smoke tests: invoke andes on various fixtures and verify
//! the PIN and TSV outputs exist with sensible content.

use std::path::PathBuf;
use std::process::Command;

/// Resolve a path relative to the workspace root (three levels above the
/// cli crate's manifest directory: cli → crates → rust → astral-speed).
fn fixture(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(rel)
        .canonicalize()
        .unwrap_or_else(|e| panic!("canonicalize {rel}: {e}"))
}

/// Build a base Command with the mandatory arguments that every test requires.
fn base_cmd(spectrum: &str, database: &str, pin: &std::path::Path) -> Command {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_andes"));
    cmd.arg("--spectrum")
        .arg(fixture(spectrum))
        .arg("--database")
        .arg(fixture(database))
        .arg("--output-pin")
        .arg(pin);
    cmd
}

// ── BSA / MGF end-to-end test (original smoke test) ─────────────────────────

#[test]
fn cli_runs_end_to_end_on_bsa_test_mgf() {
    let dir = tempfile::tempdir().expect("tempdir");
    let pin_path = dir.path().join("rust.pin");
    let tsv_path = dir.path().join("rust.tsv");

    let status = base_cmd(
        "test-fixtures/test.mgf.gz",
        "test-fixtures/BSA.fasta",
        &pin_path,
    )
    .arg("--output-tsv")
    .arg(&tsv_path)
    .arg("--decoy-prefix")
    .arg("XXX_")
    .status()
    .expect("run andes");

    assert!(status.success(), "andes exit code: {status}");
    assert!(pin_path.exists(), "PIN output not written");
    assert!(tsv_path.exists(), "TSV output not written");

    // Validate PIN header and content.
    let pin_content = std::fs::read_to_string(&pin_path).unwrap();
    assert!(
        pin_content.lines().count() > 1,
        "PIN should have header + at least 1 row"
    );
    let pin_header = pin_content.lines().next().unwrap();
    assert!(
        pin_header.starts_with("SpecId\tLabel\tScanNr"),
        "unexpected PIN header: {pin_header}"
    );

    // Assert that at least one data row carries a real BSA accession (P02769)
    // in the Proteins column — confirms real accessions are threaded through.
    let pin_has_bsa_accession = pin_content
        .lines()
        .skip(1) // skip header
        .any(|line| line.contains("P02769"));
    assert!(
        pin_has_bsa_accession,
        "PIN should contain at least one row with BSA accession 'P02769' \
         in the Proteins column (got PROT_N placeholder instead?)"
    );

    // Validate TSV header and content.
    let tsv_content = std::fs::read_to_string(&tsv_path).unwrap();
    assert!(
        tsv_content.lines().count() > 1,
        "TSV should have header + at least 1 row"
    );
    let tsv_header = tsv_content.lines().next().unwrap();
    assert!(
        tsv_header.starts_with("#SpecFile\tSpecID\tScanNum"),
        "unexpected TSV header: {tsv_header}"
    );

    // Assert TSV also has a real BSA accession.
    let tsv_has_bsa_accession = tsv_content
        .lines()
        .skip(1)
        .any(|line| line.contains("P02769"));
    assert!(
        tsv_has_bsa_accession,
        "TSV should contain at least one row with BSA accession 'P02769' \
         in the Protein column (got PROT_N placeholder instead?)"
    );
}

// ── New flag smoke tests: verify the flags parse and the binary exits 0 ──────

#[test]
fn cli_accepts_max_missed_cleavages_flag() {
    let dir = tempfile::tempdir().expect("tempdir");
    let pin_path = dir.path().join("out.pin");

    let status = base_cmd(
        "test-fixtures/test.mgf.gz",
        "test-fixtures/BSA.fasta",
        &pin_path,
    )
    .arg("--max-missed-cleavages")
    .arg("2")
    .status()
    .expect("run andes");

    assert!(status.success(), "--max-missed-cleavages 2 should exit 0, got: {status}");
}

#[test]
fn cli_accepts_min_peaks_flag() {
    let dir = tempfile::tempdir().expect("tempdir");
    let pin_path = dir.path().join("out.pin");

    let status = base_cmd(
        "test-fixtures/test.mgf.gz",
        "test-fixtures/BSA.fasta",
        &pin_path,
    )
    .arg("--min-peaks")
    .arg("5")
    .status()
    .expect("run andes");

    assert!(status.success(), "--min-peaks 5 should exit 0, got: {status}");
}

#[test]
fn cli_accepts_min_length_max_length_flags() {
    let dir = tempfile::tempdir().expect("tempdir");
    let pin_path = dir.path().join("out.pin");

    let status = base_cmd(
        "test-fixtures/test.mgf.gz",
        "test-fixtures/BSA.fasta",
        &pin_path,
    )
    .arg("--min-length")
    .arg("7")
    .arg("--max-length")
    .arg("35")
    .status()
    .expect("run andes");

    assert!(status.success(), "--min-length 7 --max-length 35 should exit 0, got: {status}");
}

// ── mzML integration smoke test: format dispatch + non-empty PIN ─────────────

// ── New flag smoke tests: --mod, --fragmentation, --protocol ──────────────────

#[test]
fn cli_accepts_mod_fragmentation_protocol_flags() {
    // Verify the TMT-CLI flags parse and the param resolver picks up a real
    // bundled model. We use the existing BSA fixture (no actual TMT spectra)
    // and pass a tiny TMT-style mods file — the binary should exit 0 because
    // all flags are valid. (--instrument was removed: analyzer resolution is
    // metadata-detected for mzML/.raw/.d and `--fragment-tol-*` for MGF.)
    let dir = tempfile::tempdir().expect("tempdir");
    let pin_path = dir.path().join("out.pin");
    let mods_path = dir.path().join("mods.txt");
    std::fs::write(
        &mods_path,
        "NumMods=2\n\
         229.162932,K,fix,any,TMT6plex\n\
         229.162932,*,fix,N-term,TMT6plex\n\
         57.021464,C,fix,any,Carbamidomethyl\n\
         15.994915,M,opt,any,Oxidation\n",
    ).unwrap();

    let status = base_cmd(
        "test-fixtures/test.mgf.gz",
        "test-fixtures/BSA.fasta",
        &pin_path,
    )
    .arg("--mods").arg(&mods_path)
    .arg("--fragmentation").arg("HCD")
    .arg("--protocol").arg("TMT")
    // Allow a wider tolerance — the TMT-labelled candidates differ in mass
    // and we just want to confirm the binary exits cleanly, not assert
    // recall on a non-TMT fixture.
    .arg("--precursor-tol").arg("100ppm")
    .status()
    .expect("run andes with TMT flags");

    assert!(
        status.success(),
        "andes should exit 0 with --mods + TMT flags, got: {status}"
    );
    assert!(pin_path.exists(), "PIN output should still be written");
}

#[test]
fn cli_rejects_invalid_protocol_index() {
    // Out-of-range --protocol must produce a non-zero exit with the
    // helpful error message from `parse_protocol`.
    let dir = tempfile::tempdir().expect("tempdir");
    let pin_path = dir.path().join("out.pin");

    let status = base_cmd(
        "test-fixtures/test.mgf.gz",
        "test-fixtures/BSA.fasta",
        &pin_path,
    )
    .arg("--protocol").arg("42")
    .status()
    .expect("run andes with bad protocol");

    assert!(!status.success(), "out-of-range --protocol must fail");
}

#[test]
fn cli_runs_end_to_end_on_tiny_mzml() {
    // tiny.pwiz.mzML is the standard fixture used by the mzML reader unit tests.
    // It is a real mzML file with MS2 spectra.  Because there is no matched FASTA,
    // we expect few or zero PSMs — but the binary must exit 0 and the PIN must be
    // written (even if it contains only the header row).
    //
    // We use BSA.fasta as the target database: it is the only fixture available.
    // The point of this test is NOT PSM recall but that the mzML code path runs
    // end-to-end without a crash or panic.
    let dir = tempfile::tempdir().expect("tempdir");
    let pin_path = dir.path().join("mzml_out.pin");

    let status = base_cmd(
        "test-fixtures/tiny.pwiz.mzML",
        "test-fixtures/BSA.fasta",
        &pin_path,
    )
    // Lower min-peaks so we don't filter out the tiny fixture's sparse spectra.
    .arg("--min-peaks")
    .arg("1")
    .status()
    .expect("run andes on mzML");

    assert!(
        status.success(),
        "andes should exit 0 on mzML input, got: {status}"
    );
    assert!(pin_path.exists(), "PIN output should be written for mzML input");

    // The PIN must at least contain a header row.
    let pin_content = std::fs::read_to_string(&pin_path).unwrap();
    let first_line = pin_content.lines().next().unwrap_or("");
    assert!(
        first_line.starts_with("SpecId\tLabel\tScanNr"),
        "PIN header should be present for mzML output; got: {first_line}"
    );
}

#[test]
fn bench_mode_max_spectra_produces_nonempty_pin() {
    // Regression for send_chunks bench-cap bug: --max-spectra 100 must not
    // drop the entire final partial chunk (which used to truncate to zero).
    let dir = tempfile::tempdir().expect("tempdir");
    let pin_path = dir.path().join("bench.pin");

    let status = base_cmd(
        "test-fixtures/test.mgf.gz",
        "test-fixtures/BSA.fasta",
        &pin_path,
    )
    .arg("--max-spectra")
    .arg("100")
    .status()
    .expect("run andes bench mode");

    assert!(status.success(), "bench mode should exit 0, got: {status}");
    assert!(pin_path.exists(), "PIN should be written in bench mode");

    let content = std::fs::read_to_string(&pin_path).unwrap();
    assert!(
        content.lines().count() > 1,
        "bench mode with --max-spectra 100 should produce header + data rows"
    );
}

#[test]
fn cli_rejects_inverted_charge_range() {
    let dir = tempfile::tempdir().expect("tempdir");
    let pin_path = dir.path().join("out.pin");

    let status = base_cmd(
        "test-fixtures/test.mgf.gz",
        "test-fixtures/BSA.fasta",
        &pin_path,
    )
    .arg("--charge")
    .arg("4..2")
    .status()
    .expect("run andes with inverted charge range");

    assert!(!status.success(), "inverted charge range must fail");
}

#[test]
fn cli_rejects_inverted_isotope_error_range() {
    let dir = tempfile::tempdir().expect("tempdir");
    let pin_path = dir.path().join("out.pin");

    let status = base_cmd(
        "test-fixtures/test.mgf.gz",
        "test-fixtures/BSA.fasta",
        &pin_path,
    )
    .arg("--isotope-error")
    .arg("3..-1")
    .status()
    .expect("run andes with inverted isotope range");

    assert!(!status.success(), "inverted isotope error range must fail");
}

#[test]
fn cli_accepts_isotope_error_min_negative_one() {
    let dir = tempfile::tempdir().expect("tempdir");
    let pin_path = dir.path().join("out.pin");

    let status = base_cmd(
        "test-fixtures/test.mgf.gz",
        "test-fixtures/BSA.fasta",
        &pin_path,
    )
    .arg("--isotope-error")
    .arg("-1..2")
    .arg("--max-spectra")
    .arg("10")
    .status()
    .expect("run andes with isotope-error -1..2");

    assert!(status.success(), "negative isotope-error MIN must parse");
    assert!(pin_path.exists());
}

#[test]
fn cli_accepts_precursor_cal_off() {
    let dir = tempfile::tempdir().expect("tempdir");
    let pin_path = dir.path().join("out.pin");

    let status = base_cmd(
        "test-fixtures/test.mgf.gz",
        "test-fixtures/BSA.fasta",
        &pin_path,
    )
    .arg("--precursor-cal")
    .arg("off")
    .arg("--max-spectra")
    .arg("10")
    .status()
    .expect("run andes with precursor-cal off");

    assert!(status.success());
    assert!(pin_path.exists());
}

/// Smoke guard for the canonical (named + merged-range) CLI surface that the
/// quantms andes module passes: named `--fragmentation/--protocol`, the
/// `--enzyme-specificity` name, the merged `--charge`/`--isotope-error` ranges,
/// and the unit-bearing `--precursor-tol`. Legacy numeric forms and split-range
/// flags were removed (andes is pre-release, no back-compat).
#[test]
fn cli_accepts_canonical_named_param_values() {
    let bsa_fasta = fixture("test-fixtures/BSA.fasta");
    let test_mgf = fixture("test-fixtures/test.mgf.gz");

    let dir = tempfile::tempdir().expect("tempdir");
    let mods_path = dir.path().join("mods.txt");
    std::fs::write(
        &mods_path,
        "NumMods=2\n\
         229.162932,K,fix,any,TMT6plex\n\
         229.162932,*,fix,N-term,TMT6plex\n\
         57.021464,C,fix,any,Carbamidomethyl\n\
         15.994915,M,opt,any,Oxidation\n",
    ).unwrap();

    let tmp = tempfile::tempdir().expect("tmpdir");
    let pin = tmp.path().join("named.pin");

    let status = base_cmd(test_mgf.to_str().unwrap(),
                          bsa_fasta.to_str().unwrap(),
                          &pin)
        .arg("--mods").arg(&mods_path)
        .arg("--fragmentation").arg("HCD")
        .arg("--protocol").arg("TMT")
        .arg("--enzyme-specificity").arg("fully")
        .arg("--charge").arg("2..5")
        .arg("--isotope-error").arg("-1..2")
        .arg("--precursor-tol").arg("100ppm")
        .status()
        .expect("named form exit");
    assert!(status.success(), "canonical named CLI form failed");

    let pin_content = std::fs::read_to_string(&pin).expect("read pin");
    assert!(!pin_content.lines().next().unwrap_or("").is_empty(), "PIN must have a header");
}

// ── MGF metadata-less model-selection routing tests ──────────────────────────

#[test]
fn mgf_no_flags_defaults_to_cid_lowres_with_warning() {
    // MGF carries no analyzer metadata. With no --fragmentation/--fragment-tol
    // flags, decision E applies: assume CID / low-res / 0.5 Da (cid_lowres_tryp)
    // and emit a warning so the user knows a default was chosen.
    let dir = tempfile::tempdir().expect("tempdir");
    let pin_path = dir.path().join("out.pin");

    let output = base_cmd(
        "test-fixtures/test.mgf.gz",
        "test-fixtures/BSA.fasta",
        &pin_path,
    )
    .output()
    .unwrap();

    let err = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "andes exited non-zero; stderr: {err}");
    assert!(err.to_lowercase().contains("cid_lowres"), "stderr: {err}");
    assert!(
        err.to_lowercase().contains("assuming") || err.to_lowercase().contains("warn"),
        "expected metadata-less default warning; stderr: {err}"
    );
}

#[test]
fn mgf_fragment_tol_ppm_selects_high_res_model() {
    // --fragment-tol-ppm on MGF input declares high-resolution MS/MS, so the
    // resolver selects a QExactive (high-res) model rather than the low-res
    // default.
    let dir = tempfile::tempdir().expect("tempdir");
    let pin_path = dir.path().join("out.pin");

    let output = base_cmd(
        "test-fixtures/test.mgf.gz",
        "test-fixtures/BSA.fasta",
        &pin_path,
    )
    .arg("--fragment-tol-ppm")
    .arg("20")
    .output()
    .unwrap();

    let err = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "andes exited non-zero; stderr: {err}");
    assert!(err.to_lowercase().contains("qexactive"), "stderr: {err}");
}

// ── Fragment-tolerance CLI flag tests ─────────────────────────────────────────

#[test]
fn fragment_tol_flags_are_mutually_exclusive() {
    let dir = tempfile::tempdir().expect("tempdir");
    let pin_path = dir.path().join("out.pin");

    let output = base_cmd(
        "test-fixtures/test.mgf.gz",
        "test-fixtures/BSA.fasta",
        &pin_path,
    )
    .arg("--fragment-tol-ppm")
    .arg("20")
    .arg("--fragment-tol-da")
    .arg("0.5")
    .output()
    .expect("run andes");

    assert!(!output.status.success(), "both flags together must fail");
    let err = String::from_utf8_lossy(&output.stderr);
    assert!(
        err.contains("cannot be used with") || err.contains("conflicts"),
        "expected conflict error in stderr, got: {err}"
    );
}
