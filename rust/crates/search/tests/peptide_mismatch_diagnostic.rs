//! One-shot diagnostic: split BSA peptide mismatches into enumerator-gap vs
//! scoring-gap. Picks up to 10 mismatching scans where Rust's top-1 target
//! peptide differs from Java's; for each, checks whether Java's peptide appears
//! anywhere in Rust's global candidate set (enumerator gap) or in the top-N
//! queue for that spectrum (scoring gap).
//!
//! Run with:
//!   cargo test --release -p search --test peptide_mismatch_diagnostic \
//!     -- --ignored --nocapture
//!
//! Output:
//!   scan 3416 ch3: Java pep "KVPQVSTPTLVEVSR" — RUST_NOT_GENERATED (enumerator gap)
//! or
//!   scan 5442 ch2: Java pep "LGEYGFQNALIVR" — generated, ranked 4 (top-1 was "GEYGFQNALIVRR")

mod common;
use common::*;

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};

use input::{FastaReader, MgfReader};
use search::{enumerate_candidates, match_spectra, SearchIndex, SearchParams};

// ── helpers ─────────────────────────────────────────────────────────────────

/// Strip flanking residues `X.PEPTIDE.Y` → `PEPTIDE`, then remove any
/// modification annotations (e.g. `+57.021`, `C+57.021`) — keeps only
/// ASCII uppercase letters.
fn strip_flanking_and_mods(pin_pep: &str) -> String {
    let core = if let Some(mid) = pin_pep.split('.').nth(1) {
        mid
    } else {
        pin_pep
    };
    core.chars().filter(|c| c.is_ascii_uppercase()).collect()
}

/// Extract a scan number from a TITLE string of the form `... scan=N`.
fn extract_scan_from_title(title: &str) -> Option<i32> {
    title
        .split_whitespace()
        .find_map(|tok| tok.strip_prefix("scan=")?.parse::<i32>().ok())
}

/// Residue-only string from a Rust Peptide (no flanking, no mod masses).
fn peptide_residue_string(p: &model::Peptide) -> String {
    p.residues.iter().map(|aa| aa.residue as char).collect()
}

// ── Java reference fixture ───────────────────────────────────────────────────

#[derive(Debug, Clone)]
struct JavaRef {
    scan_nr: i32,
    peptide: String,   // bare residues, uppercase, no mods, no flanking
    charge:  u8,
}

fn load_java_reference() -> Vec<JavaRef> {
    let path = fixture("benchmark/parity-fixtures/bsa_test_mgf_java.pin");
    let f = File::open(&path)
        .unwrap_or_else(|e| panic!("open {:?}: {}", path, e));
    let r = BufReader::new(f);
    let mut lines = r.lines();

    let header = lines.next().unwrap().unwrap();
    let cols: Vec<&str> = header.split('\t').collect();

    let scan_idx    = cols.iter().position(|&c| c == "ScanNr").expect("ScanNr");
    let label_idx   = cols.iter().position(|&c| c == "Label").expect("Label");
    let pep_idx     = cols.iter().position(|&c| c == "Peptide").expect("Peptide");
    let charge2_idx = cols.iter().position(|&c| c == "charge2").expect("charge2");
    let charge3_idx = cols.iter().position(|&c| c == "charge3").expect("charge3");

    let mut out: HashMap<i32, JavaRef> = HashMap::new();
    for line in lines {
        let line = line.unwrap();
        let fields: Vec<&str> = line.split('\t').collect();
        let max_idx = [scan_idx, label_idx, pep_idx, charge2_idx, charge3_idx]
            .iter()
            .copied()
            .max()
            .unwrap();
        if fields.len() <= max_idx {
            continue;
        }
        let label: i32 = fields[label_idx].parse().unwrap_or(0);
        if label != 1 {
            continue; // targets only
        }
        let scan: i32 = match fields[scan_idx].parse() {
            Ok(s) => s,
            Err(_) => continue,
        };
        // Keep only the first entry per scan (top-1).
        if out.contains_key(&scan) {
            continue;
        }
        let peptide = strip_flanking_and_mods(fields[pep_idx]);
        let charge = if fields[charge2_idx] == "1" {
            2u8
        } else if fields[charge3_idx] == "1" {
            3u8
        } else {
            0u8
        };
        out.insert(scan, JavaRef { scan_nr: scan, peptide, charge });
    }
    out.into_values().collect()
}

// ── diagnostic test ──────────────────────────────────────────────────────────

#[test]
#[ignore]
fn diagnose_peptide_mismatches() {
    // ── 1. Load Java reference ───────────────────────────────────────────────
    let java_refs = load_java_reference();
    eprintln!("Loaded {} Java reference PSMs", java_refs.len());

    // ── 2. Build search index + params (same as match_engine_java_parity) ───
    let target = FastaReader::load_all(BufReader::new(
        File::open(fixture("src/test/resources/BSA.fasta")).unwrap(),
    ))
    .unwrap();
    let idx    = SearchIndex::from_target_db(&target, "XXX");
    let params = SearchParams::default_tryptic(aa_set());
    let scorer = rank_scorer();

    // ── 3. Load spectra ──────────────────────────────────────────────────────
    let mgf_file = File::open(fixture("src/test/resources/test.mgf")).unwrap();
    let spectra: Vec<_> = MgfReader::new(BufReader::new(mgf_file))
        .filter_map(|r| r.ok())
        .collect();
    eprintln!("Loaded {} spectra from test.mgf", spectra.len());

    // ── 4. Run full search ───────────────────────────────────────────────────
    let queues = match_spectra(&spectra, &idx, &params, &scorer, 0.05, "XXX");

    // ── 5. Build global enumerator peptide set ───────────────────────────────
    // Collect every residue-only string that Rust's enumerator can generate
    // for BSA (target side only — Java's references are target peptides).
    let all_pep_strings: HashSet<String> = enumerate_candidates(&idx, &params, "XXX")
        .filter(|c| !c.is_decoy)
        .map(|c| peptide_residue_string(&c.peptide))
        .collect();
    eprintln!(
        "Enumerator produced {} distinct target peptide strings",
        all_pep_strings.len()
    );

    // ── 6. Build scan → spectrum index ───────────────────────────────────────
    let scan_to_spec_idx: HashMap<i32, usize> = spectra
        .iter()
        .enumerate()
        .filter_map(|(i, s)| {
            let scan = s.scan.or_else(|| extract_scan_from_title(&s.title))?;
            Some((scan, i))
        })
        .collect();

    // ── 7. Classify mismatches ───────────────────────────────────────────────
    let mut enumerator_gap_count = 0usize;
    let mut scoring_gap_count    = 0usize;
    let mut total_mismatches     = 0usize;
    let mut classify_log: Vec<String> = Vec::new();
    let mut report_remaining = 10usize;

    for jref in &java_refs {
        let spec_idx = match scan_to_spec_idx.get(&jref.scan_nr) {
            Some(&i) => i,
            None => continue, // scan not in MGF
        };
        let queue = &queues[spec_idx];
        if queue.is_empty() {
            continue;
        }

        let sorted = queue.clone().into_sorted_vec();

        // Top-1 TARGET PSM (skip decoys to match the parity test convention).
        let top_target = match sorted.iter().find(|m| !m.candidate.is_decoy) {
            Some(t) => t,
            None => continue,
        };
        let rust_top_pep = peptide_residue_string(&top_target.candidate.peptide);

        if rust_top_pep == jref.peptide {
            continue; // top-1 match — not a mismatch
        }

        // ── Mismatch: classify ───────────────────────────────────────────────
        total_mismatches += 1;

        let in_enumerator = all_pep_strings.contains(&jref.peptide);

        // Find Java's peptide's rank in this spectrum's top-N queue (if present).
        let rank_in_queue: Option<usize> = sorted
            .iter()
            .position(|m| !m.candidate.is_decoy && peptide_residue_string(&m.candidate.peptide) == jref.peptide);

        let classification = if !in_enumerator {
            enumerator_gap_count += 1;
            format!("RUST_NOT_GENERATED (enumerator gap)")
        } else {
            scoring_gap_count += 1;
            match rank_in_queue {
                Some(rank) => format!(
                    "generated, ranked {} in queue (top-1 target was '{}', spec_e_value {:.2e})",
                    rank + 1,
                    rust_top_pep,
                    top_target.spec_e_value
                ),
                None => format!(
                    "generated globally but NOT in top-N for this spectrum \
                     (evicted or precursor-filtered; top-1 target was '{}')",
                    rust_top_pep
                ),
            }
        };

        if report_remaining > 0 {
            classify_log.push(format!(
                "  scan {} ch{}: Java pep '{}' — {}",
                jref.scan_nr, jref.charge, jref.peptide, classification
            ));
            report_remaining -= 1;
        }
    }

    // ── 8. Print report ───────────────────────────────────────────────────────
    eprintln!();
    eprintln!("=== PEPTIDE MISMATCH DIAGNOSTIC ===");
    eprintln!("Java reference PSMs (target):         {}", java_refs.len());
    eprintln!("Total mismatches classified:          {}", total_mismatches);
    eprintln!(
        "  Enumerator gap (RUST_NOT_GENERATED): {} ({:.1}%)",
        enumerator_gap_count,
        if total_mismatches > 0 {
            100.0 * enumerator_gap_count as f64 / total_mismatches as f64
        } else {
            0.0
        }
    );
    eprintln!(
        "  Scoring/ranking gap:                 {} ({:.1}%)",
        scoring_gap_count,
        if total_mismatches > 0 {
            100.0 * scoring_gap_count as f64 / total_mismatches as f64
        } else {
            0.0
        }
    );
    eprintln!();
    eprintln!(
        "=== Sample of {} mismatches (first {} chronologically): ===",
        classify_log.len(),
        classify_log.len()
    );
    for line in &classify_log {
        eprintln!("{}", line);
    }
    eprintln!();
    eprintln!("Verdict: {} dominates.",
        if enumerator_gap_count >= scoring_gap_count { "ENUMERATOR GAP" } else { "SCORING/RANKING GAP" }
    );

    // Sanity check: the diagnostic found at least one mismatch.
    assert!(
        total_mismatches > 0,
        "no mismatches detected — either parity is fully closed or the diagnostic is broken"
    );
}
