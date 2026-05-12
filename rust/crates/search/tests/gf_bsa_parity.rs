//! Bulk SpecEValue Java parity histogram.
//!
//! For all 217 Java-identified PSMs from BSA + test.mgf:
//!   - Compute abs(log10(rust_spec_evalue) - log10(java_spec_evalue))
//!   - Bucket by tolerance: ≤1 OOM, ≤2 OOM, ≤3 OOM, ≤4 OOM, >4 OOM
//!   - Print the histogram and summary stats (median, max diff)
//!   - SOFT gate: ≥50% within 4 OOM (not the aspirational 95% gate)
//!
//! Reference fixture:
//!   `astral-speed/benchmark/parity-fixtures/bsa_test_mgf_java.pin`

mod common;
use common::*;

use std::fs::File;
use std::io::{BufRead, BufReader};

use search::{match_spectra, SearchIndex, SearchParams};
use input::{FastaReader, MgfReader};

/// Extract a scan number from a TITLE string of the form `... scan=N`.
fn extract_scan_from_title(title: &str) -> Option<i32> {
    title
        .split_whitespace()
        .find_map(|tok| tok.strip_prefix("scan=")?.parse::<i32>().ok())
}

/// Extract plain residue string from a Rust Peptide (no flanking, no mods).
fn peptide_residue_string(p: &model::Peptide) -> String {
    p.residues.iter().map(|aa| aa.residue as char).collect()
}

#[derive(Debug, Clone)]
struct JavaRef {
    scan_nr: i32,
    peptide: String,
    charge: u8,
    spec_evalue: f64,
}

fn load_java_reference() -> Vec<JavaRef> {
    let path = fixture("benchmark/parity-fixtures/bsa_test_mgf_java.pin");
    let f = File::open(&path).unwrap_or_else(|e| panic!("open fixture: {e}"));
    let r = BufReader::new(f);
    let mut lines = r.lines();
    let header = lines
        .next()
        .expect("header line missing")
        .expect("header read error");
    let cols: Vec<&str> = header.split('\t').collect();
    let scan_idx = cols.iter().position(|c| *c == "ScanNr").expect("ScanNr");
    let label_idx = cols.iter().position(|c| *c == "Label").expect("Label");
    let lnsev_idx = cols
        .iter()
        .position(|c| *c == "lnSpecEValue")
        .expect("lnSpecEValue");
    let pep_idx = cols.iter().position(|c| *c == "Peptide").expect("Peptide");
    let charge2_idx = cols
        .iter()
        .position(|c| *c == "charge2")
        .expect("charge2");
    let charge3_idx = cols
        .iter()
        .position(|c| *c == "charge3")
        .expect("charge3");

    let mut out = Vec::new();
    for line in lines {
        let line = line.unwrap();
        let fields: Vec<&str> = line.split('\t').collect();
        let max_idx = [scan_idx, label_idx, lnsev_idx, pep_idx, charge2_idx, charge3_idx]
            .iter()
            .copied()
            .max()
            .unwrap_or(0);
        if fields.len() <= max_idx {
            continue;
        }
        // Target PSMs only (Label = 1).
        if fields[label_idx] != "1" {
            continue;
        }
        let scan: i32 = match fields[scan_idx].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let lnsev: f64 = match fields[lnsev_idx].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let spec_evalue = lnsev.exp();

        // Strip flanking + mod-mass tokens via the shared correct parser.
        // Earlier inline `split('.').nth(1)` was buggy for peptides with mods
        // (e.g. `K.GAC+57.021LLPK.E` parsed to `"GAC"`), wildly understating
        // the population of comparable PSMs.
        let peptide = strip_flanking_and_mods(fields[pep_idx]);

        let charge = if fields[charge2_idx] == "1" {
            2
        } else if fields[charge3_idx] == "1" {
            3
        } else {
            0
        };

        out.push(JavaRef {
            scan_nr: scan,
            peptide,
            charge,
            spec_evalue,
        });
    }
    out
}

#[test]
fn phase6_task10_bsa_specevalue_parity_histogram() {
    let java_refs = load_java_reference();
    eprintln!("Loaded {} Java reference PSMs", java_refs.len());

    let target = FastaReader::load_all(BufReader::new(
        File::open(fixture("src/test/resources/BSA.fasta")).unwrap(),
    ))
    .unwrap();
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let aa = aa_set();
    let scorer = rank_scorer();
    let params = SearchParams::default_tryptic(aa.clone());
    // default_tryptic already sets: enzyme=Trypsin, isotope_error_range=-1..=2,
    // precursor_tolerance=20ppm, charge_range=2..=3.

    let mgf_file = File::open(fixture("src/test/resources/test.mgf")).unwrap();
    let spectra: Vec<_> = MgfReader::new(BufReader::new(mgf_file))
        .filter_map(|r| r.ok())
        .collect();

    // Use a broad decoy fraction (0.5) so we get a large top-N queue to search
    // for matching peptides, consistent with gf_java_parity.rs.
    let queues = match_spectra(&spectra, &idx, &params, &scorer, 0.5, "XXX");

    // Track per-PSM outcomes.
    #[derive(Debug)]
    struct MeasuredPsm {
        scan_nr: i32,
        peptide: String,
        charge: u8,
        java_sev: f64,
        rust_sev: f64,
        log_diff: f64,
    }

    let mut measured: Vec<MeasuredPsm> = Vec::new();
    let mut peptide_mismatches = 0usize;
    let mut spec_not_found = 0usize;
    let mut empty_queues = 0usize;

    for jref in &java_refs {
        // Locate the spectrum by scan number (try .scan field first, fall back to title parse).
        let spec_idx = spectra.iter().position(|s| {
            let scan_from_field = s.scan;
            let scan_from_title = extract_scan_from_title(&s.title);
            scan_from_field == Some(jref.scan_nr) || scan_from_title == Some(jref.scan_nr)
        });
        let spec_idx = match spec_idx {
            Some(i) => i,
            None => {
                spec_not_found += 1;
                continue;
            }
        };

        let queue = &queues[spec_idx];
        if queue.is_empty() {
            empty_queues += 1;
            continue;
        }

        // Search all PSMs in the queue for one whose plain residues match Java's reference.
        let top_psms = queue.clone().into_sorted_vec();
        let matched = top_psms.iter().find(|p| {
            peptide_residue_string(&p.candidate.peptide)
                .eq_ignore_ascii_case(&jref.peptide)
        });

        let psm = match matched {
            Some(p) => p,
            None => {
                peptide_mismatches += 1;
                continue;
            }
        };

        let rust_sev = psm.spec_e_value;
        // Guard against zero/negative values that would make log10 undefined.
        if rust_sev <= 0.0 || jref.spec_evalue <= 0.0 {
            peptide_mismatches += 1;
            continue;
        }
        let log_diff = (rust_sev.log10() - jref.spec_evalue.log10()).abs();

        measured.push(MeasuredPsm {
            scan_nr: jref.scan_nr,
            peptide: jref.peptide.clone(),
            charge: jref.charge,
            java_sev: jref.spec_evalue,
            rust_sev,
            log_diff,
        });
    }

    // Bucket the log10 differences: [<=1, <=2, <=3, <=4, >4].
    let mut buckets = [0_usize; 5];
    for m in &measured {
        if m.log_diff <= 1.0 {
            buckets[0] += 1;
        } else if m.log_diff <= 2.0 {
            buckets[1] += 1;
        } else if m.log_diff <= 3.0 {
            buckets[2] += 1;
        } else if m.log_diff <= 4.0 {
            buckets[3] += 1;
        } else {
            buckets[4] += 1;
        }
    }

    let total = measured.len();
    let mut sorted_diffs: Vec<f64> = measured.iter().map(|m| m.log_diff).collect();
    sorted_diffs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let median = if total > 0 {
        sorted_diffs[total / 2]
    } else {
        0.0
    };
    let max = sorted_diffs.last().copied().unwrap_or(0.0);

    // Cumulative percentage within k OOM (k = 0..=3 → buckets <=1,<=2,<=3,<=4).
    let cumulative_pct = |max_bucket: usize| -> f64 {
        if total == 0 {
            return 0.0;
        }
        let cum: usize = buckets[..=max_bucket.min(4)].iter().sum();
        cum as f64 / total as f64 * 100.0
    };

    // Identify the top 3 outliers (largest log_diff) for the commit body.
    let mut outliers: Vec<&MeasuredPsm> = measured.iter().collect();
    outliers.sort_by(|a, b| b.log_diff.partial_cmp(&a.log_diff).unwrap());
    let top_outliers: Vec<&MeasuredPsm> = outliers.into_iter().take(3).collect();

    // Print the full histogram to stderr (visible with --nocapture or in CI logs).
    eprintln!();
    eprintln!("BSA SpecEValue parity histogram");
    eprintln!("  Java reference PSMs:  {}", java_refs.len());
    eprintln!("  Spectra not found:    {}", spec_not_found);
    eprintln!("  Empty Rust queues:    {}", empty_queues);
    eprintln!("  Peptide mismatches:   {}", peptide_mismatches);
    eprintln!("  PSMs measured:        {}", total);
    eprintln!();
    eprintln!("  log10 diff buckets (per-bucket):");
    eprintln!(
        "    <=1 OOM:   {:>4}  ({:.1}%)",
        buckets[0],
        buckets[0] as f64 / total.max(1) as f64 * 100.0
    );
    eprintln!(
        "    <=2 OOM:   {:>4}  ({:.1}%)",
        buckets[1],
        buckets[1] as f64 / total.max(1) as f64 * 100.0
    );
    eprintln!(
        "    <=3 OOM:   {:>4}  ({:.1}%)",
        buckets[2],
        buckets[2] as f64 / total.max(1) as f64 * 100.0
    );
    eprintln!(
        "    <=4 OOM:   {:>4}  ({:.1}%)",
        buckets[3],
        buckets[3] as f64 / total.max(1) as f64 * 100.0
    );
    eprintln!(
        "     >4 OOM:   {:>4}  ({:.1}%)",
        buckets[4],
        buckets[4] as f64 / total.max(1) as f64 * 100.0
    );
    eprintln!();
    eprintln!("  cumulative within:");
    eprintln!("    1 OOM: {:.1}%", cumulative_pct(0));
    eprintln!("    2 OOM: {:.1}%", cumulative_pct(1));
    eprintln!("    3 OOM: {:.1}%", cumulative_pct(2));
    eprintln!("    4 OOM: {:.1}%", cumulative_pct(3));
    eprintln!();
    eprintln!("  median log10 diff: {:.3}", median);
    eprintln!("  max log10 diff:    {:.3}", max);
    eprintln!();
    eprintln!("  Top 3 outliers (largest log10 diff):");
    for (i, m) in top_outliers.iter().enumerate() {
        eprintln!(
            "    [{}] scan {:>5}  '{}'  ch{}  Java {:.3e}  Rust {:.3e}  diff {:.3}",
            i + 1,
            m.scan_nr,
            m.peptide,
            m.charge,
            m.java_sev,
            m.rust_sev,
            m.log_diff
        );
    }
    eprintln!();

    // SOFT gate: at least 50% of measured PSMs must be within 4 OOM.
    // A failure here indicates a structural bug, not just calibration drift.
    let pct_within_4 = cumulative_pct(3);
    assert!(
        total > 0,
        "no PSMs were measured (all spectra missing or queues empty)"
    );
    assert!(
        pct_within_4 >= 50.0,
        "SOFT GATE FAILED: only {:.1}% of {} measured PSMs within 4 OOM \
         (gate is 50%). This indicates a structural scoring bug worth \
         investigating.",
        pct_within_4,
        total
    );
}
