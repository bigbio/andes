//! Decoy generation parity test against Tryp_Pig_Bov.fasta.

use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

use engine::{reverse_db, target_plus_decoy};
use input::FastaReader;

fn fixture_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .join("src/test/resources/Tryp_Pig_Bov.fasta")
        .canonicalize()
        .expect("canonicalize Tryp_Pig_Bov.fasta path")
}

#[test]
fn tryp_pig_bov_reverses_to_16_decoys() {
    let path = fixture_path();
    let target = FastaReader::load_all(BufReader::new(File::open(&path).unwrap())).unwrap();
    let decoy = reverse_db(&target, "XXX");
    assert_eq!(decoy.len(), 16);
    for (t, d) in target.iter().zip(decoy.iter()) {
        assert_eq!(d.accession, format!("XXX_{}", t.accession));
        assert_eq!(d.description, t.description);
        let reversed: Vec<u8> = t.sequence.iter().rev().copied().collect();
        assert_eq!(d.sequence, reversed);
        assert_eq!(d.sequence.len(), t.sequence.len());
    }
}

#[test]
fn tryp_pig_bov_target_plus_decoy_has_32_proteins() {
    let path = fixture_path();
    let target = FastaReader::load_all(BufReader::new(File::open(&path).unwrap())).unwrap();
    let combined = target_plus_decoy(&target, "XXX");
    assert_eq!(combined.len(), 32);
    for i in 0..16 {
        assert_eq!(combined.proteins[i].accession, target.proteins[i].accession);
        assert!(combined.proteins[16 + i].accession.starts_with("XXX_"));
    }
}
