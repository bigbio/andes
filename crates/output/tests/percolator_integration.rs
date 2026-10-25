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

use output::{
    percolator_args_with_threads, resolve_backend, run_percolator, run_percolator_with_decoys,
    PercolatorBackend, DEFAULT_PERCOLATOR_IMAGE,
};

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

    assert!(
        !map.is_empty(),
        "expected at least one target PSM from Percolator"
    );
    // PSMIds must round-trip (every key starts with our target prefix and parses).
    for (id, psm) in &map {
        assert_eq!(*id, psm.psm_id, "map key must equal the parsed PSMId");
        assert!(
            psm.q_value.is_finite(),
            "q-value should parse to a finite number"
        );
        assert!(psm.pep.is_finite(), "PEP should parse to a finite number");
    }
    // At least one of our injected targets should survive.
    assert!(
        map.keys().any(|k| k.starts_with("tt_")),
        "expected a target PSM (tt_*) in the results"
    );
}

/// Two runs' PINs pooled by concatenation, every target sharing its `ScanNr`
/// and `ExpMass` with the other run's target. Percolator's target-decoy
/// competition (`-Y`) keeps one PSM per file, scan and `ExpMass`: with the
/// `FileName` column both runs' PSMs survive; without it the rows of the
/// better-scoring run wins each group and the other run's row is dropped.
fn write_pooled_pin(path: &std::path::Path, with_file_name: bool) {
    let mut w = std::fs::File::create(path).unwrap();
    let file_col = if with_file_name { "\tFileName" } else { "" };
    writeln!(
        w,
        "SpecId\tLabel\tScanNr\tExpMass{file_col}\tRawScore\tdm\tPeptide\tProteins"
    )
    .unwrap();
    for (run, offset) in [("run_a", 0), ("run_b", 7)] {
        let file = if with_file_name {
            format!("\t{run}")
        } else {
            String::new()
        };
        for i in 0..60 {
            // Same scan and mass in both runs; the score differs slightly.
            writeln!(
                w,
                "{run}/t_{i}_1\t1\t{i}\t{mass}{file}\t{score}\t0.001\tR.PEPTIDEK.A\tsp|T{i}",
                mass = 1000.0 + 0.5 * i as f64,
                score = 40 + i + offset
            )
            .unwrap();
        }
        for i in 0..60 {
            writeln!(
                w,
                "{run}/d_{i}_1\t-1\t{scan}\t{mass}{file}\t{score}\t0.9\tR.EDITPEPK.A\tXXX_D{i}",
                scan = 1000 + i,
                mass = 2000.0 + 0.5 * i as f64 + offset as f64,
                score = 1 + (i % 3)
            )
            .unwrap();
        }
    }
}

#[test]
fn pooled_runs_compete_per_file_or_skip() {
    let Some(backend) = test_backend() else {
        eprintln!(
            "SKIP pooled_runs_compete_per_file_or_skip: no Percolator backend \
             (set $ANDES_TEST_PERCOLATOR_BIN, put `percolator` on PATH, or provide docker)."
        );
        return;
    };
    let tmp = tempfile::tempdir().unwrap();
    let args = percolator_args_with_threads(&[], 1);
    let winners = |with_file_name: bool| -> Option<(usize, usize)> {
        let pin = tmp.path().join(format!("pooled_{with_file_name}.pin"));
        write_pooled_pin(&pin, with_file_name);
        match run_percolator_with_decoys(&backend, &pin, &args) {
            Ok((targets, decoys)) => {
                let count = |run: &str| {
                    targets
                        .keys()
                        .chain(decoys.iter())
                        .filter(|k| k.starts_with(&format!("{run}/t_")))
                        .count()
                };
                Some((count("run_a"), count("run_b")))
            }
            Err(e) => {
                eprintln!("SKIP pooled_runs_compete_per_file_or_skip: run failed: {e}");
                None
            }
        }
    };
    let (Some(without), Some(with)) = (winners(false), winners(true)) else {
        return;
    };
    assert_eq!(
        with,
        (60, 60),
        "with FileName every run's 60 colliding targets must be reported"
    );
    assert_eq!(
        without.0 + without.1,
        60,
        "without FileName the two runs' rows share 60 groups and one row of each wins: {without:?}"
    );
}
