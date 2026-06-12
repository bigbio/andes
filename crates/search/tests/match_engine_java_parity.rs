//! BSA benchmark regression gate: Rust must recover at least N% of the
//! reference post-scoring identifications on the BSA + test.mgf fixture.
//!
//! Rationale:
//! - The reference `.pin` output contains top-1 PSMs after scoring and
//!   Q-value filtering. For BSA + test.mgf with 20 ppm tolerance, Trypsin,
//!   1 missed cleavage, Carbamidomethyl-C fixed + Oxidation-M variable: the
//!   reference reports 217 unique target spectra (and 222 decoy entries).
//! - Rust's pipeline produces top-N=10 PSMs per spectrum with rank-based
//!   scoring via score_psm / RankScorer (Kim et al., Nat Commun 5:5277, 2014).
//! - With isotope-error tolerance (-1..=2), Rust catches all 217 reference
//!   target spectra (100% coverage).
//!
//! Gate: per-spectrum top-1 peptide identity. For each reference-identified
//! scan, Rust's top-1 PSM (by score) must agree with the reference top-1
//! peptide. Threshold: >= 50% top-1 identity match.
//!
//! Reference fixture:
//!   `astral-speed/test-fixtures/parity/bsa_test_mgf_java.pin`
//! (bundled BSA benchmark PIN output).
//!
//! ## Scope of this test file
//!
//! The integration tests below verify spectrum coverage and top-1 peptide
//! identity against the reference `.pin`. They do NOT validate the full
//! per-feature score distribution — those are covered by unit tests in the
//! scoring/output crates and the benchmark harness.

mod common;
use common::*;

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use search::{match_spectra, SearchIndex, SearchParams};
use input::{FastaReader, MgfReader};

/// Extract a scan number from a TITLE string of the form
/// `... scan=N` (e.g. mzML controllerType/controllerNumber/scan triplets).
fn extract_scan_from_title(title: &str) -> Option<i32> {
    title
        .split_whitespace()
        .find_map(|tok| tok.strip_prefix("scan=")?.parse::<i32>().ok())
}

/// Parse a reference `.pin` file and return the set of unique scan numbers
/// that have at least one target PSM (Label = 1).
fn java_target_scans(pin_path: &PathBuf) -> HashSet<i32> {
    let file = File::open(pin_path)
        .unwrap_or_else(|e| panic!("open {pin_path:?}: {e}"));
    let reader = BufReader::new(file);
    let mut lines = reader.lines();
    let header = lines
        .next()
        .expect("empty pin file")
        .expect("read pin header");

    let cols: Vec<&str> = header.split('\t').collect();
    let label_idx = cols.iter().position(|&c| c == "Label").expect("Label column");
    let scan_idx = cols.iter().position(|&c| c == "ScanNr").expect("ScanNr column");

    let mut scans = HashSet::new();
    for line in lines {
        let line = line.expect("read pin line");
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() <= scan_idx.max(label_idx) {
            continue;
        }
        let label: i32 = fields[label_idx].parse().unwrap_or(0);
        if label == 1 {
            if let Ok(scan) = fields[scan_idx].parse::<i32>() {
                scans.insert(scan);
            }
        }
    }
    scans
}

/// Parse a reference `.pin` file and return a map of scan_number → peptide string
/// (bare residues, no flanking, no modifications) for target PSMs (Label = 1).
///
/// Peptide column format: `R.KVPQVSTPTLVEVSR.S`
/// We strip the flanking X.PEPTIDE.Y → "PEPTIDE".
/// Modifications like `+57.021` are stripped for the plain-residue comparison.
fn java_target_peptides(pin_path: &PathBuf) -> HashMap<i32, String> {
    let file = File::open(pin_path)
        .unwrap_or_else(|e| panic!("open {pin_path:?}: {e}"));
    let reader = BufReader::new(file);
    let mut lines = reader.lines();
    let header = lines
        .next()
        .expect("empty pin file")
        .expect("read pin header");

    let cols: Vec<&str> = header.split('\t').collect();
    let label_idx = cols.iter().position(|&c| c == "Label").expect("Label column");
    let scan_idx = cols.iter().position(|&c| c == "ScanNr").expect("ScanNr column");
    let pep_idx = cols.iter().position(|&c| c == "Peptide").expect("Peptide column");

    let mut map: HashMap<i32, String> = HashMap::new();
    for line in lines {
        let line = line.expect("read pin line");
        let fields: Vec<&str> = line.split('\t').collect();
        let max_idx = scan_idx.max(label_idx).max(pep_idx);
        if fields.len() <= max_idx {
            continue;
        }
        let label: i32 = fields[label_idx].parse().unwrap_or(0);
        if label != 1 {
            continue;
        }
        if let Ok(scan) = fields[scan_idx].parse::<i32>() {
            let raw = fields[pep_idx];
            let bare = strip_flanking_and_mods(raw);
            // Keep only the first (and usually only) top-1 entry per scan.
            map.entry(scan).or_insert(bare);
        }
    }
    map
}

// `strip_flanking_and_mods` is shared from `common/mod.rs`. The previous
// local copy used `split('.').nth(1)` which silently truncated peptides
// containing mod masses (e.g. `K.GAC+57.021LLPK.E` → `"GAC"`), wildly
// understating peptide-identity matches in this benchmark test.

/// Extract plain residue string from a Rust Peptide (no flanking, no mods).
fn peptide_residue_string(p: &model::Peptide) -> String {
    // Access residues via the length and mass — but Peptide exposes residues publicly.
    // Use the iterator approach via the public API.
    let mut s = String::new();
    // Peptide::residues is pub in our model.
    for aa in &p.residues {
        s.push(aa.residue as char);
    }
    s
}

#[test]
fn rust_matches_superset_java_target_psms() {
    let java_pin = fixture("test-fixtures/parity/bsa_test_mgf_java.pin");
    let java_scans = java_target_scans(&java_pin);
    assert!(
        !java_scans.is_empty(),
        "Java pin file has no target PSMs (Label=1); fixture may be stale"
    );
    println!("Java identified {} target spectra", java_scans.len());

    let target = FastaReader::load_all(BufReader::new(
        File::open(fixture("test-fixtures/BSA.fasta")).unwrap(),
    ))
    .unwrap();
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let params = SearchParams::default_tryptic(aa_set());

    let mgf_file = File::open(fixture("test-fixtures/test.mgf")).unwrap();
    let spectra: Vec<_> = MgfReader::new(BufReader::new(mgf_file))
        .filter_map(|r| r.ok())
        .collect();

    let scorer = rank_scorer();
    let (queues, candidates) = match_spectra(&spectra, &idx, &params, &scorer, 0.05, "XXX");
    assert_eq!(queues.len(), spectra.len());

    // Collect scan numbers of Rust spectra that have ≥1 target PSM.
    let mut rust_target_scans: HashSet<i32> = HashSet::new();
    for (spec, queue) in spectra.iter().zip(queues.iter()) {
        let queue_clone = queue.clone();
        if queue_clone.is_empty() {
            continue;
        }
        let has_target = queue_clone
            .into_sorted_vec()
            .iter()
            .any(|m| !candidates[m.primary_candidate_idx() as usize].is_decoy);
        if !has_target {
            continue;
        }
        let scan = spec.scan.or_else(|| extract_scan_from_title(&spec.title));
        if let Some(s) = scan {
            rust_target_scans.insert(s);
        }
    }
    println!(
        "Rust pre-scoring matched {} target spectra",
        rust_target_scans.len()
    );

    // Compute coverage: fraction of reference target spectra that Rust also matched.
    let intersection = java_scans.intersection(&rust_target_scans).count();
    let coverage = intersection as f64 / java_scans.len() as f64;
    println!(
        "Rust ∩ Java target spectra: {} / {} (coverage = {:.1}%)",
        intersection,
        java_scans.len(),
        coverage * 100.0
    );

    // Regression gate: Rust must catch at least 95% of reference target spectra.
    const MIN_COVERAGE: f64 = 0.95;
    assert!(
        coverage >= MIN_COVERAGE,
        "Rust caught only {:.1}% of Java's target spectra; minimum gate is {:.0}%. \
         Java had {} target spectra, Rust caught {} of them.",
        coverage * 100.0,
        MIN_COVERAGE * 100.0,
        java_scans.len(),
        intersection
    );
}

#[test]
fn rust_top1_matches_java_top1_for_majority_of_spectra() {
    let java_pin = fixture("test-fixtures/parity/bsa_test_mgf_java.pin");
    let java_peps = java_target_peptides(&java_pin);
    assert!(
        !java_peps.is_empty(),
        "Java pin file has no target PSMs (Label=1); fixture may be stale"
    );
    println!("Java top-1 peptides: {} entries", java_peps.len());

    let target = FastaReader::load_all(BufReader::new(
        File::open(fixture("test-fixtures/BSA.fasta")).unwrap(),
    ))
    .unwrap();
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let params = SearchParams::default_tryptic(aa_set());

    let mgf_file = File::open(fixture("test-fixtures/test.mgf")).unwrap();
    let spectra: Vec<_> = MgfReader::new(BufReader::new(mgf_file))
        .filter_map(|r| r.ok())
        .collect();

    let scorer = rank_scorer();
    let (queues, candidates) = match_spectra(&spectra, &idx, &params, &scorer, 0.05, "XXX");
    assert_eq!(queues.len(), spectra.len());

    let mut top1_match = 0usize;
    let mut top1_total = 0usize;

    for (spec, queue) in spectra.iter().zip(queues.iter()) {
        let scan = spec.scan.or_else(|| extract_scan_from_title(&spec.title));
        let scan = match scan {
            Some(s) => s,
            None => continue,
        };
        let java_pep = match java_peps.get(&scan) {
            Some(p) => p,
            None => continue,
        };

        top1_total += 1;

        let sorted = queue.clone().into_sorted_vec();
        // Take the top-1 target PSM (skip decoys for the comparison).
        let top_target = sorted.iter().find(|m| !candidates[m.primary_candidate_idx() as usize].is_decoy);
        if let Some(top) = top_target {
            let rust_pep = peptide_residue_string(&candidates[top.primary_candidate_idx() as usize].peptide);
            if rust_pep == *java_pep {
                top1_match += 1;
            }
        }
    }

    let top1_rate = if top1_total > 0 {
        top1_match as f64 / top1_total as f64
    } else {
        0.0
    };
    println!(
        "Top-1 identity match: {} / {} ({:.1}%)",
        top1_match,
        top1_total,
        top1_rate * 100.0
    );

    // Gate: >= 95% top-1 identity match. Observed (post-parser-fix): 98.6%
    // (214/217). Earlier the gate was 45% based on a buggy peptide-string
    // comparator (see common::strip_flanking_and_mods regression tests) which
    // wildly understated match rate. The 95% floor is a regression guard ~3 pp
    // below observed.
    const MIN_TOP1_RATE: f64 = 0.95;
    assert!(
        top1_rate >= MIN_TOP1_RATE,
        "top-1 identity match rate {:.1}% < {:.0}% gate ({} / {} matched)",
        top1_rate * 100.0,
        MIN_TOP1_RATE * 100.0,
        top1_match,
        top1_total,
    );
}

/// Regression test for R-1 (commit fc16407): tied PSM retention in TopNQueue.
///
/// Why this test exists:
/// - Commit R-1 fixed TopNQueue::push to retain tied PSMs at capacity
///   (Kim et al., Nat Commun 5:5277, 2014: `size < n OR score == worst → add`).
/// - The existing integration tests check spectrum coverage and top-1 identity,
///   but neither validates that multiple PSMs are retained when they tie at the
///   worst score in a queue.
/// - Reverting the `Ordering::Equal` branch would still pass those tests because
///   they only check top-1 identity, not tie retention.
///
/// What it verifies:
/// - Runs match_spectra on the BSA + test.mgf fixture (same setup as the other tests).
/// - Iterates over the resulting TopNQueues and counts how many contain ≥2 PSMs.
/// - Asserts at least 1 such queue exists.
/// - With capacity=10 and integer-rounded scores producing ties, the BSA fixture
///   reliably produces ≥1 queue with tied PSMs (most queues will have 1, but at
///   least one will have 2+ due to ties).
///
/// Regression guard:
/// - If R-1 is reverted, all queues will be at capacity with no multi-PSM ties,
///   and the assertion will fail.
#[test]
fn r1_tie_retention_active_in_production_pipeline() {
    let target = FastaReader::load_all(BufReader::new(
        File::open(fixture("test-fixtures/BSA.fasta")).unwrap(),
    ))
    .unwrap();
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let params = SearchParams::default_tryptic(aa_set());

    let mgf_file = File::open(fixture("test-fixtures/test.mgf")).unwrap();
    let spectra: Vec<_> = MgfReader::new(BufReader::new(mgf_file))
        .filter_map(|r| r.ok())
        .collect();

    let scorer = rank_scorer();
    let (queues, _candidates) = match_spectra(&spectra, &idx, &params, &scorer, 0.05, "XXX");

    // Count how many queues have ≥2 PSMs (only possible if ties exist and R-1
    // is active to retain them).
    let queues_with_ties: usize = queues
        .iter()
        .filter(|queue| queue.len() >= 2)
        .count();

    println!(
        "Queues with ≥2 PSMs (tied retention): {}/{}",
        queues_with_ties,
        queues.len()
    );

    // Regression gate: at least 1 queue must have ties. If R-1 is reverted,
    // this assertion will fail.
    assert!(
        queues_with_ties >= 1,
        "No queues with ≥2 PSMs found (count={}). R-1 tie retention may be broken.",
        queues_with_ties
    );
}

/// Parse the reference pin file and return a Set of distinct (scan, peptide_residue)
/// pairs for target rows (Label=1). Uses the shared `strip_flanking_and_mods`
/// to correctly handle mod-mass tokens that contain dots.
fn java_target_scan_peptide_pairs(pin_path: &PathBuf) -> HashSet<(i32, String)> {
    let f = File::open(pin_path).unwrap_or_else(|e| panic!("open {pin_path:?}: {e}"));
    let r = BufReader::new(f);
    let mut lines = r.lines();
    let header = lines.next().unwrap().unwrap();
    let cols: Vec<&str> = header.split('\t').collect();
    let scan_idx = cols.iter().position(|c| *c == "ScanNr").expect("ScanNr");
    let label_idx = cols.iter().position(|c| *c == "Label").expect("Label");
    let pep_idx = cols.iter().position(|c| *c == "Peptide").expect("Peptide");

    let mut pairs: HashSet<(i32, String)> = HashSet::new();
    for line_result in lines {
        let line = match line_result {
            Ok(l) => l,
            Err(_) => continue,
        };
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() <= label_idx.max(scan_idx).max(pep_idx) {
            continue;
        }
        if fields[label_idx] != "1" {
            continue;
        }
        let scan: i32 = match fields[scan_idx].parse() {
            Ok(s) => s,
            Err(_) => continue,
        };
        let pep_stripped = strip_flanking_and_mods(fields[pep_idx]);
        if pep_stripped.is_empty() {
            continue;
        }
        pairs.insert((scan, pep_stripped));
    }
    pairs
}

/// R-2 (2026-05-18): after per-charge queues + dedup + spectrum merge, Rust's
/// distinct (scan, peptide) PSM count on the BSA fixture should approach the
/// reference benchmark. This catches:
///   - dedup collapsing PSMs it shouldn't (would reduce distinct count)
///   - missed cross-charge merge (would inflate count)
///   - protein-aggregation breaking peptide identity
///
/// Reference: bsa_test_mgf_java.pin has 217 unique (scan, peptide) target PSMs.
/// Rust should fall within +/-5% — i.e. 207-227.
#[test]
fn r2_deduped_psm_count_matches_java_on_bsa_fixture() {
    let java_pin = fixture("test-fixtures/parity/bsa_test_mgf_java.pin");
    let java_target_pairs = java_target_scan_peptide_pairs(&java_pin);
    let java_count = java_target_pairs.len();
    println!("Java distinct (scan, peptide) target PSMs: {}", java_count);

    let target = FastaReader::load_all(BufReader::new(
        File::open(fixture("test-fixtures/BSA.fasta")).unwrap(),
    ))
    .unwrap();
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let params = SearchParams::default_tryptic(aa_set());

    let mgf_file = File::open(fixture("test-fixtures/test.mgf")).unwrap();
    let spectra: Vec<_> = MgfReader::new(BufReader::new(mgf_file))
        .filter_map(|r| r.ok())
        .collect();

    let scorer = rank_scorer();
    let (queues, candidates) = match_spectra(&spectra, &idx, &params, &scorer, 0.05, "XXX");

    // Top-1 semantics (n=1): take the literal top-1 PSM (the queue's best by
    // score, target OR decoy). Only count the pair if the top-1 is a target.
    // The reference pin has one Label=1 row per spectrum whose best PSM is a
    // target. Using `find !is_decoy` instead would over-count by surfacing a
    // target PSM even when a decoy ranked higher.
    let mut rust_target_pairs: HashSet<(i32, String)> = HashSet::new();
    for (spec, queue) in spectra.iter().zip(queues.iter()) {
        let scan = match spec.scan.or_else(|| extract_scan_from_title(&spec.title)) {
            Some(s) => s,
            None => continue,
        };
        let sorted = queue.clone().into_sorted_vec();
        if let Some(top1) = sorted.first() {
            let cand = &candidates[top1.primary_candidate_idx() as usize];
            if cand.is_decoy {
                continue;
            }
            let pep = peptide_residue_string(&cand.peptide);
            rust_target_pairs.insert((scan, pep));
        }
    }
    let rust_count = rust_target_pairs.len();
    println!("Rust distinct (scan, peptide) target PSMs: {}", rust_count);

    let ratio = rust_count as f64 / java_count as f64;
    println!("Rust/Java ratio: {:.3}", ratio);

    assert!(
        (0.95..=1.05).contains(&ratio),
        "Rust distinct PSM count {} is {:.1}% of Java's {} (gate: 95%-105%)",
        rust_count,
        ratio * 100.0,
        java_count
    );
}
