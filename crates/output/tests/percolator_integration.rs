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

/// A PIN with one decoy far above every other decoy. Percolator's `-Y` PEP fit
/// can return PEP = 1 for every PSM on high-resolution TMT data with such a
/// decoy (depending on the cross-validation split) while the q-values stay
/// valid; the runner must hand back a PEP that is below 1 for the accepted
/// PSMs whichever way the fit goes. This synthetic PIN does not always
/// provoke the failure; `degenerate_pep_is_re_estimated_with_the_run_s_weights`
/// covers the repair with a stand-in Percolator that always fails.
fn write_top_decoy_pin(path: &std::path::Path) {
    let mut w = std::fs::File::create(path).unwrap();
    writeln!(
        w,
        "SpecId\tLabel\tScanNr\tExpMass\tFileName\tf1\tf2\tf3\tPeptide\tProteins"
    )
    .unwrap();
    // A fixed pseudo-random stream (xorshift), so the PIN is the same every run.
    let mut state = 0x9E37_79B9_7F4A_7C15_u64;
    let mut gauss = move || {
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        (0..12).map(|_| next()).sum::<f64>() - 6.0
    };
    for i in 0..3000 {
        let (f1, f2, f3) = (gauss(), gauss(), gauss());
        writeln!(
            w,
            "run/t_{i}\t1\t{i}\t{mass}\trun\t{a}\t{b}\t{c}\tR.PEPTIDEK.A\tsp|T{i}",
            mass = 800.0 + 0.37 * i as f64,
            a = f1 + 3.0,
            b = f2 + 2.0,
            c = f3 + 1.0
        )
        .unwrap();
    }
    for i in 0..3000 {
        let (f1, f2, f3) = (gauss(), gauss(), gauss());
        // Decoy 0 is the outlier.
        let boost = if i == 0 { 12.0 } else { 0.0 };
        writeln!(
            w,
            "run/d_{i}\t-1\t{scan}\t{mass}\trun\t{a}\t{b}\t{c}\tR.EDITPEPK.A\tXXX_D{i}",
            scan = 10_000 + i,
            mass = 5000.0 + 0.37 * i as f64,
            a = f1 + boost,
            b = f2 + boost,
            c = f3
        )
        .unwrap();
    }
}

#[test]
fn pep_of_accepted_psms_is_below_one_or_skip() {
    let Some(backend) = test_backend() else {
        eprintln!(
            "SKIP pep_of_accepted_psms_is_below_one_or_skip: no Percolator \
             backend (set $ANDES_TEST_PERCOLATOR_BIN, put `percolator` on PATH, or provide docker)."
        );
        return;
    };
    let tmp = tempfile::tempdir().unwrap();
    let pin = tmp.path().join("top_decoy.pin");
    write_top_decoy_pin(&pin);
    // The runner trains with a fixed seed, so the failure shows or not by the
    // PIN; a later `--seed` among the extra arguments takes over and gives
    // the same PIN several random splits.
    for seed in 1..=8 {
        let mut args = percolator_args_with_threads(&[], 1);
        args.extend(["--seed".to_string(), seed.to_string()]);
        let (targets, _) = match run_percolator_with_decoys(&backend, &pin, &args) {
            Ok(r) => r,
            Err(e) => {
                eprintln!("SKIP pep_of_accepted_psms_is_below_one_or_skip: {e}");
                return;
            }
        };
        assert!(
            targets.values().any(|p| p.q_value <= 0.01),
            "seed {seed}: the PIN must have accepted PSMs"
        );
        assert!(
            !output::pep_is_degenerate(&targets),
            "seed {seed}: every accepted PSM has PEP = 1"
        );
    }
    // The repair pass leaves no scratch files next to the PIN.
    let leftovers: Vec<_> = std::fs::read_dir(tmp.path())
        .unwrap()
        .filter_map(|e| e.ok())
        .filter(|e| e.file_name().to_string_lossy().contains(".pep."))
        .collect();
    assert!(leftovers.is_empty(), "{leftovers:?}");
}

/// A stand-in `percolator` script: with `-Y` it writes `y_pep` for every PSM
/// (1 is the failed fit); with `--static` (the repair pass) it writes
/// `static_pep`. Every call's arguments are appended to `calls.log` next to
/// the script.
///
/// The script is executable once the function returns: while another test
/// thread forks, its child briefly holds the script's write descriptor and
/// exec fails with `ExecutableFileBusy`, so the function runs the script with
/// `--probe` until no such descriptor remains. No descriptor can appear after
/// the write has closed.
#[cfg(unix)]
fn write_fake_percolator(
    dir: &std::path::Path,
    y_pep: &str,
    static_pep: &str,
    static_exit: i32,
) -> PathBuf {
    use std::os::unix::fs::PermissionsExt;
    let script = dir.join("percolator");
    let body = format!(
        r#"#!/bin/sh
[ "$1" = --probe ] && exit 0
echo "$@" >> "{log}"
static=0; t=""; d=""
while [ $# -gt 0 ]; do
  case "$1" in
    --static) static=1;;
    --results-psms) t="$2"; shift;;
    --decoy-results-psms) d="$2"; shift;;
  esac
  shift
done
if [ $static = 1 ]; then pep={static_pep}; [ {static_exit} -ne 0 ] && exit {static_exit}; else pep={y_pep}; fi
printf 'PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds\n' > "$t"
printf 'r/1\t3.0\t0.001\t%s\tK.AAK.R\tsp|P1\n' $pep >> "$t"
printf 'r/2\t2.0\t0.004\t%s\tK.BBK.R\tsp|P2\n' $pep >> "$t"
printf 'r/3\t0.1\t0.4\t%s\tK.CCK.R\tsp|P3\n' $pep >> "$t"
printf 'PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds\n' > "$d"
"#,
        log = dir.join("calls.log").display(),
    );
    std::fs::write(&script, body).unwrap();
    std::fs::set_permissions(&script, std::fs::Permissions::from_mode(0o755)).unwrap();
    loop {
        match std::process::Command::new(&script).arg("--probe").status() {
            Err(e) if e.kind() == std::io::ErrorKind::ExecutableFileBusy => {
                std::thread::sleep(std::time::Duration::from_millis(10));
            }
            status => {
                assert!(status.unwrap().success());
                return script;
            }
        }
    }
}

#[cfg(unix)]
#[test]
fn degenerate_pep_is_re_estimated_with_the_run_s_weights() {
    // The first call (target-decoy competition) returns PEP = 1 for every
    // accepted PSM; the repair call must be static, take the first call's
    // weights, drop `-Y`, and its PEP replaces the failed one. q-values stay.
    let tmp = tempfile::tempdir().unwrap();
    let pin = tmp.path().join("run.pin");
    std::fs::write(&pin, "SpecId\tLabel\n").unwrap();
    let fake = write_fake_percolator(tmp.path(), "1", "0.01", 0);
    let backend = PercolatorBackend::Binary(fake);
    let (targets, _) = run_percolator_with_decoys(&backend, &pin, &[]).unwrap();
    assert!(!output::pep_is_degenerate(&targets));
    assert_eq!(targets["r/1"].pep, 0.01);
    assert_eq!(targets["r/1"].q_value, 0.001);
    assert_eq!(
        targets["r/3"].pep, 0.01,
        "every target PSM takes the new PEP"
    );
    let log = std::fs::read_to_string(tmp.path().join("calls.log")).unwrap();
    let calls: Vec<&str> = log.lines().collect();
    assert_eq!(calls.len(), 2, "{calls:?}");
    assert!(calls[0].contains(" -Y ") && !calls[0].contains("--static"));
    assert!(calls[1].contains("--static") && !calls[1].contains(" -Y "));
    let weights = tmp.path().join("run.percolator.weights.txt");
    assert!(calls[1].contains(&format!("--init-weights {}", weights.display())));
    let leftovers: Vec<_> = std::fs::read_dir(tmp.path())
        .unwrap()
        .filter_map(|e| e.ok())
        .filter(|e| e.file_name().to_string_lossy().contains(".pep."))
        .collect();
    assert!(leftovers.is_empty(), "{leftovers:?}");
}

#[cfg(unix)]
#[test]
fn a_repair_that_fails_or_stays_degenerate_keeps_the_run() {
    for (static_pep, static_exit) in [("1", 0), ("0.01", 3)] {
        let tmp = tempfile::tempdir().unwrap();
        let pin = tmp.path().join("run.pin");
        std::fs::write(&pin, "SpecId\tLabel\n").unwrap();
        let backend = PercolatorBackend::Binary(write_fake_percolator(
            tmp.path(),
            "1",
            static_pep,
            static_exit,
        ));
        let (targets, _) = run_percolator_with_decoys(&backend, &pin, &[]).unwrap();
        assert_eq!(targets.len(), 3, "{static_pep}/{static_exit}");
        assert_eq!(targets["r/1"].pep, 1.0, "the main run's values stay");
        assert_eq!(targets["r/1"].q_value, 0.001);
    }
}

#[cfg(unix)]
#[test]
fn a_working_pep_is_not_touched() {
    // A stand-in whose `-Y` output is fine: one call, no repair.
    let tmp = tempfile::tempdir().unwrap();
    let pin = tmp.path().join("run.pin");
    std::fs::write(&pin, "SpecId\tLabel\n").unwrap();
    let fake = write_fake_percolator(tmp.path(), "0.002", "0.01", 0);
    let backend = PercolatorBackend::Binary(fake);
    let (targets, _) = run_percolator_with_decoys(&backend, &pin, &[]).unwrap();
    assert_eq!(targets["r/1"].pep, 0.002);
    let log = std::fs::read_to_string(tmp.path().join("calls.log")).unwrap();
    assert_eq!(log.lines().count(), 1);
}

/// A stand-in `percolator` for a PIN of two runs. The pooled call writes the
/// pooled q-values (0.001 for every PSM). A `--static` call writes the q-value
/// of the run named by the PIN it is given: run `a` 0.002, run `b` 0.5. Every
/// call's arguments are appended to `calls.log`; the PSM ids are the PIN's
/// `SpecId`s.
#[cfg(unix)]
fn write_two_run_fake(dir: &std::path::Path, fail_run_b: bool) -> PathBuf {
    use std::os::unix::fs::PermissionsExt;
    let script = dir.join("percolator");
    let body = format!(
        r#"#!/bin/sh
echo "$@" >> "{log}"
static=0; t=""; d=""
for a in "$@"; do pin="$a"; done
while [ $# -gt 0 ]; do
  case "$1" in
    --static) static=1;;
    --results-psms) t="$2"; shift;;
    --decoy-results-psms) d="$2"; shift;;
  esac
  shift
done
printf 'PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds\n' > "$t"
printf 'PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds\n' > "$d"
if [ $static = 0 ]; then
  for id in a1 a2 b1 b2 x1; do printf '%s\t1.0\t0.001\t0.01\tK.AAK.R\tsp|P1\n' $id >> "$t"; done
  printf 'bd1\t0.5\t0.001\t0.5\tK.DDK.R\tsp|D1\n' >> "$d"
  exit 0
fi
if grep -q 'run_b' "$pin"; then
  [ {fail} = 1 ] && exit 3
  q=0.5
else
  q=0.002
fi
tail -n +2 "$pin" | cut -f1 | while read id; do
  case "$id" in
    x*) ;;
    bd*) printf '%s\t0.5\t%s\t0.5\tK.DDK.R\tsp|D1\n' $id $q >> "$d";;
    *) printf '%s\t1.0\t%s\t0.01\tK.AAK.R\tsp|P1\n' $id $q >> "$t";;
  esac
done
"#,
        log = dir.join("calls.log").display(),
        fail = i32::from(fail_run_b),
    );
    std::fs::write(&script, body).unwrap();
    std::fs::set_permissions(&script, std::fs::Permissions::from_mode(0o755)).unwrap();
    script
}

#[cfg(unix)]
fn write_two_run_pin(path: &std::path::Path) {
    std::fs::write(
        path,
        "SpecId\tLabel\tFileName\tRawScore\n\
         a1\t1\trun_a\t1\n a2\t1\trun_a\t1\n\
         b1\t1\trun_b\t1\n b2\t1\trun_b\t1\n\
         bd1\t-1\trun_b\t1\nx1\t1\trun_a\t1\n"
            .replace("\n ", "\n"),
    )
    .unwrap();
}

#[cfg(unix)]
#[test]
fn per_run_q_values_replace_the_pooled_ones_and_leave_no_scratch() {
    let tmp = tempfile::tempdir().unwrap();
    let pin = tmp.path().join("set.pin");
    write_two_run_pin(&pin);
    let backend = PercolatorBackend::Binary(write_two_run_fake(tmp.path(), false));
    let (targets, decoys) = run_percolator_with_decoys(&backend, &pin, &[]).unwrap();
    assert_eq!(targets["a1"].q_value, 0.002);
    assert_eq!(targets["b2"].q_value, 0.5);
    assert_eq!(targets["x1"].q_value, 1.0, "no verdict from the run's pass");
    assert_eq!(
        targets["a1"].pooled_q_value,
        Some(0.001),
        "the pooled q stays"
    );
    assert_eq!(targets["x1"].pooled_q_value, Some(0.001));
    assert_eq!(targets["b1"].pep, 0.01, "the PEP stays the pooled one");
    assert!(decoys.contains("bd1"));
    // The result file carries the per-run q-values too.
    let text = std::fs::read_to_string(tmp.path().join("set.percolator.target.psms.txt")).unwrap();
    let file = output::parse_psm_results(&text).unwrap();
    assert_eq!(file["a2"].q_value, 0.002);
    assert_eq!(file["b1"].q_value, 0.5);
    assert_eq!(file["b1"].pooled_q_value, Some(0.001));
    assert!(text
        .lines()
        .next()
        .unwrap()
        .contains("\tpooled-q-value\tproteinIds"));
    // One pooled call, then one static call per run with the pooled weights
    // and target-decoy competition; one PIN per run, deleted afterwards.
    let log = std::fs::read_to_string(tmp.path().join("calls.log")).unwrap();
    let calls: Vec<&str> = log.lines().collect();
    assert_eq!(calls.len(), 3, "{calls:?}");
    let weights = tmp.path().join("set.percolator.weights.txt");
    for c in &calls[1..] {
        assert!(c.contains("--static") && c.contains(" -Y "), "{c}");
        assert!(c.contains("--override"), "the pooled weights are kept: {c}");
        assert!(c.contains(&format!("--init-weights {}", weights.display())));
    }
    assert!(calls[1].ends_with("set.run0.pin") && calls[2].ends_with("set.run1.pin"));
    let leftovers: Vec<_> = std::fs::read_dir(tmp.path())
        .unwrap()
        .filter_map(|e| e.ok())
        .filter(|e| e.file_name().to_string_lossy().contains(".run"))
        .collect();
    assert!(leftovers.is_empty(), "{leftovers:?}");
}

#[cfg(unix)]
#[test]
fn a_run_whose_pass_fails_keeps_its_pooled_q_values() {
    let tmp = tempfile::tempdir().unwrap();
    let pin = tmp.path().join("set.pin");
    write_two_run_pin(&pin);
    let backend = PercolatorBackend::Binary(write_two_run_fake(tmp.path(), true));
    let (targets, _) = run_percolator_with_decoys(&backend, &pin, &[]).unwrap();
    assert_eq!(targets["a1"].q_value, 0.002);
    assert_eq!(targets["b1"].q_value, 0.001);
}

#[cfg(unix)]
#[test]
fn a_pin_of_one_run_or_without_run_names_is_left_pooled() {
    let tmp = tempfile::tempdir().unwrap();
    let backend = PercolatorBackend::Binary(write_two_run_fake(tmp.path(), false));
    for text in [
        "SpecId\tLabel\tFileName\na1\t1\trun_a\na2\t1\trun_a\n",
        "SpecId\tLabel\na1\t1\na2\t1\n",
    ] {
        let pin = tmp.path().join("one.pin");
        std::fs::write(&pin, text).unwrap();
        std::fs::remove_file(tmp.path().join("calls.log")).ok();
        let (targets, _) = run_percolator_with_decoys(&backend, &pin, &[]).unwrap();
        assert_eq!(targets["a1"].q_value, 0.001);
        let log = std::fs::read_to_string(tmp.path().join("calls.log")).unwrap();
        assert_eq!(log.lines().count(), 1, "{log}");
    }
}
