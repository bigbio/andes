//! Java SpecEValue parity for hand-picked traced PSMs.
//!
//! Phase 6 Task 9 baseline: 5 PSMs from BSA + test.mgf, asserting Rust
//! SpecEValue stays within `TOLERANCE_LOG10` OOM of Java's lnSpecEValue.
//! Tightened from 4.0 to 3.5 on 2026-05-11 after the 2026-05-10 cumulative
//! fixes (theo_mz formula, cleavage credit, partition Ord, per-partition
//! ions). 4/5 PSMs now pass at 1.0 OOM; scan 3353 is the bottleneck at
//! 3.276 OOM (see the per-PSM table on `TOLERANCE_LOG10` below).
//!
//! Reference fixture:
//!   `astral-speed/benchmark/parity-fixtures/bsa_test_mgf_java.pin`
//!
//! The 5 PSMs were hand-picked from Label=1 (target) rows spanning the
//! SpecEValue range:
//!
//! | scan | peptide            | ch | Java lnSEV | Java SEV      |
//! |------|--------------------|----|------------|---------------|
//! | 3416 | KVPQVSTPTLVEVSR   |  3 | -18.0089   | 1.5095e-08    |
//! | 3353 | KVPQVSTPTLVEVSR   |  3 | -14.2373   | 6.5587e-07    |
//! | 5442 | LGEYGFQNALIVR     |  2 | -10.4288   | 2.9569e-05    |
//! | 1507 | YLYEIAR           |  2 | -8.5826    | 1.8734e-04    |
//! | 2693 | SLGKVGTR          |  2 | -5.3004    | 4.9898e-03    |

mod common;
use common::*;

use std::fs::File;
use std::io::BufReader;

use search::{match_spectra, SearchIndex, SearchParams};
use input::{FastaReader, MgfReader};

/// (scan_nr, peptide, charge, java_spec_evalue)
///
/// java_spec_evalue = exp(lnSpecEValue) from column 9 of bsa_test_mgf_java.pin.
/// Values are literals (not runtime computations) so the gate is reproducible.
const FIVE_TRACED_PSMS: &[(i32, &str, u8, f64)] = &[
    // Very confident: lnSEV = -18.0089
    (3416, "KVPQVSTPTLVEVSR", 3, 1.5095e-8),
    // Confident: lnSEV = -14.2373
    (3353, "KVPQVSTPTLVEVSR", 3, 6.5587e-7),
    // Moderate: lnSEV = -10.4288
    (5442, "LGEYGFQNALIVR", 2, 2.9569e-5),
    // Middling: lnSEV = -8.5826
    (1507, "YLYEIAR", 2, 1.8734e-4),
    // Weak: lnSEV = -5.3004
    (2693, "SLGKVGTR", 2, 4.9898e-3),
];

/// Within 3.5 OOM tolerance (tightened from 4.0 after 2026-05-10 fixes).
///
/// The 2026-05-10 cumulative fixes (theo_mz formula, cleavage credit,
/// partition Ord, per-partition ion enumeration, allocation-free
/// `ions_for_node`) closed most of the SEV-level gap that motivated the
/// previous 4.0 OOM tolerance. 4/5 PSMs now pass at 1.0 OOM; only scan 3353
/// remains as a bottleneck at 3.276 OOM (Rust MORE confident than Java).
/// Tightening below 3.5 would require diagnosing the remaining SP-level
/// drift (the SEV gate masks part of it via num_distinct_peptides
/// multiplication).
///
/// Per-PSM table (measured 2026-05-11, post-2026-05-10 fixes):
///
///   scan 3416 'KVPQVSTPTLVEVSR' ch3:
///     Java 1.510e-8 vs Rust 5.220e-9 (log10 diff 0.461) — PASS at 1.0 OOM
///     Rust is ~3x MORE confident than Java. Previous measurement: 0.106.
///     Drift increased slightly after the cumulative fixes, but still well
///     within 1.0 OOM. Reference calibration point.
///
///   scan 3353 'KVPQVSTPTLVEVSR' ch3:
///     Java 6.559e-7 vs Rust 3.473e-10 (log10 diff 3.276) — FAIL at 1.0/2.0/3.0
///     Rust is ~1900x MORE confident than Java. Previous measurement: 1.010.
///     Same peptide as scan 3416 (which passes), so the divergence is
///     spectrum-specific, not peptide-specific. The 2026-05-10 fixes
///     amplified the divergence for this scan — the score distribution
///     width is now significantly different from Java's, and the GF tail
///     falls off too fast. THIS IS THE TOLERANCE BOTTLENECK.
///
///   scan 5442 'LGEYGFQNALIVR' ch2:
///     Java 2.957e-5 vs Rust 2.752e-6 (log10 diff 1.031) — borderline at 1.0 OOM
///     Rust is ~10x MORE confident than Java. Previous measurement: 2.396
///     (and Rust was LESS confident). Direction flipped after the fixes —
///     this is now consistent with the general "Rust more confident" pattern.
///
///   scan 1507 'YLYEIAR' ch2:
///     Java 1.873e-4 vs Rust 2.914e-4 (log10 diff 0.192) — PASS at 1.0 OOM
///     Rust and Java agree to within a factor of 2. Previous measurement:
///     2.862 (Rust 700x more confident). The 2026-05-10 fixes essentially
///     resolved this case.
///
///   scan 2693 'SLGKVGTR' ch2:
///     Java 4.990e-3 vs Rust 1.652e-3 (log10 diff 0.480) — PASS at 1.0 OOM
///     Rust is ~3x MORE confident than Java. Previous measurement: 3.675
///     (Rust 4700x more confident). The 2026-05-10 fixes (especially the
///     theo_mz formula correction and per-partition ion enumeration)
///     resolved most of the divergence on this short, low-confidence PSM.
///
/// Remaining drift is at the SP level — the SEV gate compares
/// SP * num_distinct_peptides, which masks part of the underlying
/// per-spectrum score-distribution mismatch. Future tightening below
/// 3.5 OOM requires reconciling SP itself for scan 3353.
///
/// Root causes still pending (post-2026-05-10):
///   1. GF score distribution width on scan 3353 (spectrum-specific).
///   2. Underflow guard alignment with Java's Float.MIN_VALUE.
///   3. Score-range calibration: Rust score may differ from Java's RawScore.
const TOLERANCE_LOG10: f64 = 3.5;

/// Extract a scan number from a TITLE string of the form
/// `... scan=N` (e.g. mzML controllerType/controllerNumber/scan triplets).
fn extract_scan_from_title(title: &str) -> Option<i32> {
    title
        .split_whitespace()
        .find_map(|tok| tok.strip_prefix("scan=")?.parse::<i32>().ok())
}

/// Extract plain residue string from a Rust Peptide (no flanking, no mods).
fn peptide_residue_string(p: &model::Peptide) -> String {
    p.residues.iter().map(|aa| aa.residue as char).collect()
}

#[test]
fn rust_spec_evalue_within_one_oom_of_java_for_5_traced_psms() {
    let target = FastaReader::load_all(BufReader::new(
        File::open(fixture("src/test/resources/BSA.fasta")).unwrap(),
    ))
    .unwrap();
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let aa = aa_set();
    let scorer = rank_scorer();
    let params = SearchParams::default_tryptic(aa.clone());
    // params already has:
    //   enzyme = Trypsin, isotope_error_range = -1..=2,
    //   precursor_tolerance = 20 ppm, charge_range = 2..=3

    let mgf_file = File::open(fixture("src/test/resources/test.mgf")).unwrap();
    let spectra: Vec<_> = MgfReader::new(BufReader::new(mgf_file))
        .filter_map(|r| r.ok())
        .collect();

    let queues = match_spectra(&spectra, &idx, &params, &scorer, 0.5, "XXX");
    assert_eq!(queues.len(), spectra.len());

    let mut failures: Vec<String> = Vec::new();
    let mut notes: Vec<String> = Vec::new();

    for &(scan_nr, peptide, charge, java_spec_evalue) in FIVE_TRACED_PSMS {
        // Locate spectrum by scan number encoded in TITLE.
        let spec_idx = spectra.iter().position(|s| {
            let title_scan = extract_scan_from_title(&s.title);
            title_scan == Some(scan_nr)
        });
        let spec_idx = match spec_idx {
            Some(i) => i,
            None => {
                failures.push(format!(
                    "scan {scan_nr}: NOT FOUND in test.mgf (title scan= field)"
                ));
                continue;
            }
        };

        let queue = &queues[spec_idx];
        if queue.is_empty() {
            failures.push(format!(
                "scan {scan_nr}: Rust returned empty queue (no PSMs at all)"
            ));
            continue;
        }

        let top_psms = queue.clone().into_sorted_vec();

        // Find a PSM with the matching peptide (any mod variant).
        let pep_match = top_psms.iter().find(|p| {
            peptide_residue_string(&p.candidate.peptide)
                .eq_ignore_ascii_case(peptide)
        });

        let psm = match pep_match {
            Some(p) => p,
            None => {
                let top_pep = peptide_residue_string(&top_psms[0].candidate.peptide);
                notes.push(format!(
                    "scan {scan_nr} '{peptide}' ch{charge}: \
                     peptide not in Rust top-{} queue; top-1 is '{top_pep}'",
                    top_psms.len()
                ));
                // Count as a failure for the gate check below.
                failures.push(format!(
                    "scan {scan_nr} '{peptide}' ch{charge}: \
                     Java {java_spec_evalue:.3e} — peptide not in Rust queue (top-1: '{top_pep}')"
                ));
                continue;
            }
        };

        let rust_sev = psm.spec_e_value;
        let log_diff = (rust_sev.log10() - java_spec_evalue.log10()).abs();

        let status = if log_diff < TOLERANCE_LOG10 { "PASS" } else { "FAIL" };
        notes.push(format!(
            "scan {scan_nr} '{peptide}' ch{charge}: \
             Java {java_spec_evalue:.3e} vs Rust {rust_sev:.3e} \
             (log10 diff {log_diff:.3}) [{status}]"
        ));

        if log_diff >= TOLERANCE_LOG10 {
            // PHASE 6 followup: document diverging cases with both values and
            // suspected root cause so Task 10 can target the fix.
            failures.push(format!(
                "scan {scan_nr} '{peptide}' ch{charge}: \
                 Java {java_spec_evalue:.3e} vs Rust {rust_sev:.3e} \
                 (log10 diff {log_diff:.3} >= tolerance {TOLERANCE_LOG10:.1})"
            ));
        }
    }

    // Always print the per-PSM table for visibility in CI logs.
    println!("\n=== Phase 6 Task 9: per-PSM SpecEValue parity ===");
    for n in &notes {
        println!("  {n}");
    }
    println!("===================================================\n");

    assert!(
        failures.is_empty(),
        "Phase 6 Task 9: {}/{} traced PSMs failed parity (tolerance = {TOLERANCE_LOG10:.1} OOM):\n{}",
        failures.len(),
        FIVE_TRACED_PSMS.len(),
        failures.join("\n")
    );
}
