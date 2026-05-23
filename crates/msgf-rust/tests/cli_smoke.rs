//! End-to-end smoke tests: invoke msgf-rust on various fixtures and verify
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
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_msgf-rust"));
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
        "test-fixtures/test.mgf",
        "test-fixtures/BSA.fasta",
        &pin_path,
    )
    .arg("--output-tsv")
    .arg(&tsv_path)
    .arg("--decoy-prefix")
    .arg("XXX_")
    .status()
    .expect("run msgf-rust");

    assert!(status.success(), "msgf-rust exit code: {status}");
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
        "test-fixtures/test.mgf",
        "test-fixtures/BSA.fasta",
        &pin_path,
    )
    .arg("--max-missed-cleavages")
    .arg("2")
    .status()
    .expect("run msgf-rust");

    assert!(status.success(), "--max-missed-cleavages 2 should exit 0, got: {status}");
}

#[test]
fn cli_accepts_min_peaks_flag() {
    let dir = tempfile::tempdir().expect("tempdir");
    let pin_path = dir.path().join("out.pin");

    let status = base_cmd(
        "test-fixtures/test.mgf",
        "test-fixtures/BSA.fasta",
        &pin_path,
    )
    .arg("--min-peaks")
    .arg("5")
    .status()
    .expect("run msgf-rust");

    assert!(status.success(), "--min-peaks 5 should exit 0, got: {status}");
}

#[test]
fn cli_accepts_min_length_max_length_flags() {
    let dir = tempfile::tempdir().expect("tempdir");
    let pin_path = dir.path().join("out.pin");

    let status = base_cmd(
        "test-fixtures/test.mgf",
        "test-fixtures/BSA.fasta",
        &pin_path,
    )
    .arg("--min-length")
    .arg("7")
    .arg("--max-length")
    .arg("35")
    .status()
    .expect("run msgf-rust");

    assert!(status.success(), "--min-length 7 --max-length 35 should exit 0, got: {status}");
}

// ── mzML integration smoke test: format dispatch + non-empty PIN ─────────────

// ── New flag smoke tests: --mod, --fragmentation, --instrument, --protocol ────

#[test]
fn cli_accepts_mod_fragmentation_instrument_protocol_flags() {
    // Verify the new TMT-CLI flags parse and the param resolver picks up a
    // real bundled .param file. We use the existing BSA fixture (no actual
    // TMT spectra) and pass a tiny TMT-style mods file — the binary should
    // exit 0 because all flags are valid and the resolver finds
    // HCD_QExactive_Tryp_TMT.param.
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
        "test-fixtures/test.mgf",
        "test-fixtures/BSA.fasta",
        &pin_path,
    )
    .arg("--mod").arg(&mods_path)
    .arg("--fragmentation").arg("3")
    .arg("--instrument").arg("3")
    .arg("--protocol").arg("4")
    // Allow a wider tolerance — the TMT-labelled candidates differ in mass
    // and we just want to confirm the binary exits cleanly, not assert
    // recall on a non-TMT fixture.
    .arg("--precursor-tol-ppm").arg("100")
    .status()
    .expect("run msgf-rust with TMT flags");

    assert!(
        status.success(),
        "msgf-rust should exit 0 with --mod + TMT flags, got: {status}"
    );
    assert!(pin_path.exists(), "PIN output should still be written");
}

#[test]
fn cli_rejects_invalid_protocol_index() {
    // Out-of-range --protocol must produce a non-zero exit with the
    // helpful error message from `resolve_bundled_param`.
    let dir = tempfile::tempdir().expect("tempdir");
    let pin_path = dir.path().join("out.pin");

    let status = base_cmd(
        "test-fixtures/test.mgf",
        "test-fixtures/BSA.fasta",
        &pin_path,
    )
    .arg("--protocol").arg("42")
    .status()
    .expect("run msgf-rust with bad protocol");

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
    .expect("run msgf-rust on mzML");

    assert!(
        status.success(),
        "msgf-rust should exit 0 on mzML input, got: {status}"
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
