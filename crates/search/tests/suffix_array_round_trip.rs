//! Round-trip + Java fixture parity tests for SuffixArray I/O.

use std::io::Cursor;
use std::path::PathBuf;

use model::{CompactFastaSequence, Protein, ProteinDb};
use search::SuffixArray;

#[test]
fn sa_round_trip_preserves_arrays() {
    let db = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "".into(),
            sequence: b"MKWVTFISLLLLFSSAYSRGV".to_vec(),
        }],
    };
    let cf = CompactFastaSequence::from_protein_db(&db);
    let sa = SuffixArray::build(&cf);

    let mut csarr_bytes = Vec::new();
    let mut cnlcp_bytes = Vec::new();
    sa.write_to(&mut csarr_bytes, &mut cnlcp_bytes).unwrap();

    let parsed = SuffixArray::read_from(
        &mut Cursor::new(&csarr_bytes),
        &mut Cursor::new(&cnlcp_bytes),
    )
    .unwrap();

    assert_eq!(parsed.indices, sa.indices);
    assert_eq!(parsed.nlcps, sa.nlcps);
}

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../target/test-classes")
        .join(name)
        .canonicalize()
        .unwrap_or_else(|e| panic!("canonicalize {name}: {e}"))
}

#[test]
fn read_tryp_pig_bov_revcat_csarr_cnlcp() {
    let csarr_bytes = std::fs::read(fixture("Tryp_Pig_Bov.revCat.csarr")).unwrap();
    let cnlcp_bytes = std::fs::read(fixture("Tryp_Pig_Bov.revCat.cnlcp")).unwrap();
    let sa = SuffixArray::read_from(
        &mut Cursor::new(&csarr_bytes),
        &mut Cursor::new(&cnlcp_bytes),
    )
    .unwrap();
    assert!(!sa.indices.is_empty());
    assert_eq!(sa.indices.len(), sa.nlcps.len());
    // Tryp_Pig_Bov.revCat has ~32 proteins ~5K residues; SA has ~9565 entries.
    assert!(sa.indices.len() > 1000);
}
