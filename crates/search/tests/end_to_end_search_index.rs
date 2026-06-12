//! End-to-end: load FASTA → build SearchIndex → assert decoy-generation
//! shape invariants. Exercises the production pipeline (FASTA reader →
//! decoy gen) on real fixtures.

use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

use search::SearchIndex;
use input::FastaReader;

fn fasta(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join("test-fixtures")
        .join(name)
        .canonicalize()
        .unwrap_or_else(|e| panic!("canonicalize {name}: {e}"))
}

#[test]
fn bsa_end_to_end() {
    let target = FastaReader::load_all(BufReader::new(File::open(fasta("BSA.fasta")).unwrap())).unwrap();
    let idx = SearchIndex::from_target_db(&target, "XXX");
    assert_eq!(idx.db.len(), 2);  // 1 target + 1 decoy
}

#[test]
fn tryp_pig_bov_end_to_end() {
    let target = FastaReader::load_all(BufReader::new(File::open(fasta("Tryp_Pig_Bov.fasta")).unwrap())).unwrap();
    let idx = SearchIndex::from_target_db(&target, "XXX");
    assert_eq!(idx.db.len(), 32);
}
