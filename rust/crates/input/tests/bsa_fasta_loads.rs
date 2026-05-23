//! Load `astral-speed/rust/test-fixtures/BSA.fasta` (1 protein, ~607
//! residues) and assert basic invariants.

use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

use input::FastaReader;

fn fixture_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .join("rust/test-fixtures/BSA.fasta")
        .canonicalize()
        .expect("canonicalize BSA.fasta path")
}

#[test]
fn bsa_loads_exactly_one_protein() {
    let path = fixture_path();
    let file = File::open(&path).unwrap_or_else(|e| panic!("open {path:?}: {e}"));
    let db = FastaReader::load_all(BufReader::new(file)).unwrap();
    assert_eq!(db.len(), 1, "expected 1 protein in BSA.fasta");
}

#[test]
fn bsa_protein_has_expected_accession_and_length() {
    let path = fixture_path();
    let file = File::open(&path).unwrap();
    let db = FastaReader::load_all(BufReader::new(file)).unwrap();
    let p = &db.proteins[0];
    assert_eq!(p.accession, "sp|P02769|ALBU_BOVIN");
    assert!(p.sequence.len() >= 500, "BSA sequence too short: {}", p.sequence.len());
    assert!(p.sequence.iter().all(|&b| b.is_ascii_uppercase() && b.is_ascii_alphabetic()),
        "non-uppercase or non-alpha residue found");
}
