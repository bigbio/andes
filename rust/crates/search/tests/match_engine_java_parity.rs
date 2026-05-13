//! Java parity regression gate: Rust must catch at least N% of Java's
//! post-scoring identifications.
//!
//! Rationale:
//! - Java MS-GF+'s `.pin` output contains top-1 PSMs after scoring + Q-value
//!   filtering. For BSA + test.mgf with 20 ppm tolerance, Trypsin, 1 missed
//!   cleavage, Carbamidomethyl-C fixed + Oxidation-M variable: Java reports
//!   217 unique target spectra (and 222 decoy entries).
//! - Rust's Phase 5 pipeline produces top-N=10 PSMs per spectrum with real
//!   rank-based scoring via score_psm / RankScorer.
//! - With isotope-error tolerance (`-ti -1..=2` matching Java's default),
//!   Rust catches ALL 217 of Java's target spectra (100% coverage).
//!
//! Gate: per-spectrum top-1 peptide identity. For each Java-identified scan,
//! Rust's top-1 PSM (by score) must agree with Java's top-1 peptide.
//! Threshold: >= 50% top-1 identity match.
//!
//! Reference fixture:
//!   `astral-speed/benchmark/parity-fixtures/bsa_test_mgf_java.pin`
//! generated via:
//!   java -Xmx4g -jar target/MSGFPlus.jar \
//!     -s src/test/resources/test.mgf \
//!     -d src/test/resources/BSA.fasta \
//!     -mod benchmark/parity-fixtures/bsa_test_mgf_mods.txt \
//!     -o /tmp/bsa.pin -tda 1 -t 20ppm -ti -1,2 -m 3 -inst 0 -e 1 -ntt 2 \
//!     -minLength 6 -maxLength 40 -minCharge 2 -maxCharge 3 \
//!     -maxMissedCleavages 1 -n 1 -addFeatures 1 -msLevel 2

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

/// Parse a Java `.pin` file and return the set of unique scan numbers
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

/// Parse a Java `.pin` file and return a map of scan_number → peptide string
/// (bare residues, no flanking, no modifications) for target PSMs (Label = 1).
///
/// Java's Peptide column format: `R.KVPQVSTPTLVEVSR.S`
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
// understating peptide-identity matches in this parity test.

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
    let java_pin = fixture("benchmark/parity-fixtures/bsa_test_mgf_java.pin");
    let java_scans = java_target_scans(&java_pin);
    assert!(
        !java_scans.is_empty(),
        "Java pin file has no target PSMs (Label=1); fixture may be stale"
    );
    println!("Java identified {} target spectra", java_scans.len());

    let target = FastaReader::load_all(BufReader::new(
        File::open(fixture("src/test/resources/BSA.fasta")).unwrap(),
    ))
    .unwrap();
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let params = SearchParams::default_tryptic(aa_set());

    let mgf_file = File::open(fixture("src/test/resources/test.mgf")).unwrap();
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
            .any(|m| !candidates[m.candidate_idx as usize].is_decoy);
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

    // Compute coverage: fraction of Java's target spectra that Rust also matched.
    let intersection = java_scans.intersection(&rust_target_scans).count();
    let coverage = intersection as f64 / java_scans.len() as f64;
    println!(
        "Rust ∩ Java target spectra: {} / {} (coverage = {:.1}%)",
        intersection,
        java_scans.len(),
        coverage * 100.0
    );

    // Regression gate: Rust must catch at least 95% of Java's target spectra.
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
    let java_pin = fixture("benchmark/parity-fixtures/bsa_test_mgf_java.pin");
    let java_peps = java_target_peptides(&java_pin);
    assert!(
        !java_peps.is_empty(),
        "Java pin file has no target PSMs (Label=1); fixture may be stale"
    );
    println!("Java top-1 peptides: {} entries", java_peps.len());

    let target = FastaReader::load_all(BufReader::new(
        File::open(fixture("src/test/resources/BSA.fasta")).unwrap(),
    ))
    .unwrap();
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let params = SearchParams::default_tryptic(aa_set());

    let mgf_file = File::open(fixture("src/test/resources/test.mgf")).unwrap();
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
        let top_target = sorted.iter().find(|m| !candidates[m.candidate_idx as usize].is_decoy);
        if let Some(top) = top_target {
            let rust_pep = peptide_residue_string(&candidates[top.candidate_idx as usize].peptide);
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
    // wildly understated parity. The 95% floor is a regression guard ~3 pp
    // below observed — tighten further once any further parity improvements
    // land.
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
