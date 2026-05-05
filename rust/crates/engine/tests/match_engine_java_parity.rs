//! Java parity regression gate: Rust must catch at least N% of Java's
//! post-scoring identifications.
//!
//! Rationale:
//! - Java MS-GF+'s `.pin` output contains top-1 PSMs after scoring + Q-value
//!   filtering. For BSA + test.mgf with 20 ppm tolerance, Trypsin, 1 missed
//!   cleavage, Carbamidomethyl-C fixed + Oxidation-M variable: Java reports
//!   217 unique target spectra (and 222 decoy entries).
//! - Rust's Phase 4e pipeline produces top-N=10 PSMs per spectrum AFTER
//!   precursor-mass filtering only (no scoring yet).
//! - With isotope-error tolerance (`-ti -1..=2` matching Java's default),
//!   Rust catches ALL 217 of Java's target spectra (100% coverage). Rust's
//!   total of 308 target spectra exceeds Java's 217 because Rust hasn't
//!   yet applied scoring/Q-value filtering (Phase 5+ work).
//!
//! When Phase 5 scoring lands, this gate evolves to per-spectrum top-1
//! peptide-identity parity (Rust's top-1 should equal Java's top-1).
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

use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use engine::{
    match_spectra, AminoAcidSetBuilder, ModLocation, Modification, ResidueSpec,
    SearchIndex, SearchParams,
};
use input::{FastaReader, MgfReader};

fn fixture(path: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .join(path)
        .canonicalize()
        .unwrap_or_else(|e| panic!("canonicalize {path}: {e}"))
}

fn aa_set() -> engine::AminoAcidSet {
    let cam = Modification {
        name: "Carbamidomethyl".into(),
        mass_delta: 57.02146,
        residue: ResidueSpec::Specific(b'C'),
        location: ModLocation::Anywhere,
        fixed: true,
        accession: None,
    };
    let ox = Modification {
        name: "Oxidation".into(),
        mass_delta: 15.99491,
        residue: ResidueSpec::Specific(b'M'),
        location: ModLocation::Anywhere,
        fixed: false,
        accession: None,
    };
    AminoAcidSetBuilder::new_standard()
        .add_fixed_mod(cam)
        .add_variable_mod(ox)
        .build()
        .unwrap()
}

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

    let queues = match_spectra(&spectra, &idx, &params, "XXX");
    assert_eq!(queues.len(), spectra.len());

    // Collect scan numbers of Rust spectra that have ≥1 target PSM.
    // test.mgf has no SCANS= lines; scan numbers are embedded in
    // TITLE as `scan=N`. We extract them inline (the Phase 3a MGF
    // reader keeps scan=None for these spectra; a follow-up could
    // teach the reader to parse TITLE for `scan=N` as a fallback).
    let mut rust_target_scans: HashSet<i32> = HashSet::new();
    for (spec, queue) in spectra.iter().zip(queues.iter()) {
        let queue_clone = queue.clone();
        if queue_clone.is_empty() {
            continue;
        }
        let has_target = queue_clone
            .into_sorted_vec()
            .iter()
            .any(|m| !m.candidate.is_decoy);
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
    // Current baseline: 100% (217/217). The 95% floor catches accidental
    // regressions while leaving a tiny buffer for spectrum-level edge cases.
    // When Phase 5 scoring lands, this gate evolves to per-spectrum top-1
    // peptide-identity match (Rust's top-1 == Java's top-1) at the same 95%
    // threshold.
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
