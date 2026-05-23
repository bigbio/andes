//! End-to-end Phase 4b+4c: load FASTA → build SearchIndex → assert
//! shape invariants. Exercises the full pipeline (FASTA reader →
//! decoy gen → CompactFastaSequence → SA build) on real fixtures.

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
    assert!(idx.compact.size > 1000);  // BSA ~607 residues × 2 + sentinels
    assert_eq!(idx.sa.indices.len(), idx.compact.size as usize);
}

#[test]
fn tryp_pig_bov_end_to_end() {
    let target = FastaReader::load_all(BufReader::new(File::open(fasta("Tryp_Pig_Bov.fasta")).unwrap())).unwrap();
    let idx = SearchIndex::from_target_db(&target, "XXX");
    assert_eq!(idx.db.len(), 32);
    assert_eq!(idx.sa.indices.len(), idx.compact.size as usize);
}
