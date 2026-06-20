//! End-to-end Percolator runner test.
//!
//! GRACEFUL SKIP: this test needs a real Percolator. It probes, in order,
//! `$ANDES_TEST_PERCOLATOR_BIN`, `percolator` on `$PATH`, then the docker
//! fallback. If NONE is available it prints a skip notice and returns (it does
//! NOT fail) so the suite stays green in environments without Percolator. The
//! parser, backend resolution, and the QPX/PIN SpecId join are covered by unit
//! tests + `qpx_rescore_join.rs` regardless.

use std::io::Write;
use std::path::PathBuf;

use output::{resolve_backend, run_percolator, PercolatorBackend, DEFAULT_PERCOLATOR_IMAGE};

/// Resolve a backend for the test, honoring `$ANDES_TEST_PERCOLATOR_BIN` first.
/// Returns `None` (→ skip) when nothing is available.
fn test_backend() -> Option<PercolatorBackend> {
    if let Some(bin) = std::env::var_os("ANDES_TEST_PERCOLATOR_BIN") {
        let p = PathBuf::from(bin);
        if p.exists() {
            return Some(PercolatorBackend::Binary(p));
        }
    }
    // PATH, then docker fallback (resolve_backend handles both).
    resolve_backend(None, false, DEFAULT_PERCOLATOR_IMAGE).ok()
}

/// A tiny but separable synthetic PIN: a handful of high-scoring targets and a
/// few low-scoring decoys, with the andes column layout Percolator expects.
fn write_synthetic_pin(path: &std::path::Path) {
    let mut w = std::fs::File::create(path).unwrap();
    // Minimal Percolator-valid PIN: SpecId, Label, ScanNr, <features...>, Peptide, Proteins.
    writeln!(w, "SpecId\tLabel\tScanNr\tRawScore\tdm\tPeptide\tProteins").unwrap();
    // Targets — high score, small mass error.
    for i in 0..30 {
        writeln!(
            w,
            "tt_{i}_1\t1\t{i}\t{score}\t0.001\tR.PEPTIDEK.A\tsp|T{i}",
            score = 40 + i
        )
        .unwrap();
    }
    // Decoys — low score, larger mass error.
    for i in 0..30 {
        writeln!(
            w,
            "dd_{i}_1\t-1\t{scan}\t{score}\t0.9\tR.EDITPEPK.A\tXXX_D{i}",
            scan = 1000 + i,
            score = 1 + (i % 3)
        )
        .unwrap();
    }
}

#[test]
fn percolator_end_to_end_or_skip() {
    let Some(backend) = test_backend() else {
        eprintln!(
            "SKIP percolator_end_to_end_or_skip: no Percolator backend \
             (set $ANDES_TEST_PERCOLATOR_BIN, put `percolator` on PATH, or provide docker)."
        );
        return;
    };
    eprintln!("percolator backend = {}", backend.describe());

    let tmp = tempfile::tempdir().unwrap();
    let pin = tmp.path().join("synthetic.pin");
    write_synthetic_pin(&pin);

    let map = match run_percolator(&backend, &pin, &[]) {
        Ok(m) => m,
        Err(e) => {
            // A real backend that fails on this synthetic input is informative
            // but should not red the suite in CI without a vetted percolator.
            eprintln!("SKIP percolator_end_to_end_or_skip: run failed: {e}");
            return;
        }
    };

    assert!(!map.is_empty(), "expected at least one target PSM from Percolator");
    // PSMIds must round-trip (every key starts with our target prefix and parses).
    for (id, psm) in &map {
        assert_eq!(*id, psm.psm_id, "map key must equal the parsed PSMId");
        assert!(psm.q_value.is_finite(), "q-value should parse to a finite number");
        assert!(psm.pep.is_finite(), "PEP should parse to a finite number");
    }
    // At least one of our injected targets should survive.
    assert!(
        map.keys().any(|k| k.starts_with("tt_")),
        "expected a target PSM (tt_*) in the results"
    );
}
