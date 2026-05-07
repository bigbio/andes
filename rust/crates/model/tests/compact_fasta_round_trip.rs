//! Round-trip + Java fixture parity tests for CompactFastaSequence I/O.

use std::io::Cursor;
use std::path::PathBuf;

use model::{CompactFastaSequence, Protein, ProteinDb};

fn small_db() -> ProteinDb {
    ProteinDb {
        proteins: vec![
            Protein {
                accession: "P1".into(),
                description: "first".into(),
                sequence: b"MKWVTFISLL".to_vec(),
            },
            Protein {
                accession: "P2".into(),
                description: "second".into(),
                sequence: b"AGCTAGCTAGCT".to_vec(),
            },
        ],
    }
}

#[test]
fn cseq_canno_round_trip_preserves_structure() {
    let db = small_db();
    let cf = CompactFastaSequence::from_protein_db(&db);

    let mut cseq_bytes = Vec::new();
    let mut canno_bytes = Vec::new();
    cf.write_to(&mut cseq_bytes, &mut canno_bytes).unwrap();

    let parsed = CompactFastaSequence::read_from(
        &mut Cursor::new(&cseq_bytes),
        &mut Cursor::new(&canno_bytes),
    )
    .unwrap();

    assert_eq!(parsed.size, cf.size);
    assert_eq!(parsed.sequence, cf.sequence);
    assert_eq!(parsed.annotations.len(), cf.annotations.len());
    for (a, b) in parsed.annotations.iter().zip(cf.annotations.iter()) {
        assert_eq!(a.start, b.start);
        assert_eq!(a.accession, b.accession);
        assert_eq!(a.description, b.description);
    }
}

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../../target/test-classes")
        .join(name)
        .canonicalize()
        .unwrap_or_else(|e| panic!("canonicalize {name}: {e}"))
}

#[test]
fn read_bsa_canno_text_format() {
    let cseq_bytes = std::fs::read(fixture("BSA.cseq")).unwrap();
    let canno_bytes = std::fs::read(fixture("BSA.canno")).unwrap();
    let cf = CompactFastaSequence::read_from(
        &mut Cursor::new(&cseq_bytes),
        &mut Cursor::new(&canno_bytes),
    )
    .unwrap();
    assert_eq!(cf.protein_count(), 1);
    assert_eq!(cf.annotations[0].accession, "sp|P02769|ALBU_BOVIN");
    assert!(cf.size > 500);
}
