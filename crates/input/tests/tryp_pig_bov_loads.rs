//! Load `astral-speed/test-fixtures/Tryp_Pig_Bov.fasta` (16
//! proteins) and assert count + per-protein invariants.

use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

use input::FastaReader;

fn fixture_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join("test-fixtures/Tryp_Pig_Bov.fasta")
        .canonicalize()
        .expect("canonicalize Tryp_Pig_Bov.fasta path")
}

#[test]
fn tryp_pig_bov_loads_16_proteins() {
    let path = fixture_path();
    let file = File::open(&path).unwrap_or_else(|e| panic!("open {path:?}: {e}"));
    let db = FastaReader::load_all(BufReader::new(file)).unwrap();
    assert_eq!(db.len(), 16, "expected 16 proteins, got {}", db.len());
}

#[test]
fn each_protein_well_formed() {
    let path = fixture_path();
    let file = File::open(&path).unwrap();
    let db = FastaReader::load_all(BufReader::new(file)).unwrap();
    for (i, p) in db.iter().enumerate() {
        assert!(!p.accession.is_empty(), "protein {} has empty accession", i);
        assert!(!p.sequence.is_empty(), "protein {} ({}) has empty sequence", i, p.accession);
        assert!(p.sequence.iter().all(|&b| b.is_ascii_uppercase() && b.is_ascii_alphabetic()),
            "protein {} ({}) has non-uppercase or non-alpha residue", i, p.accession);
    }
}
