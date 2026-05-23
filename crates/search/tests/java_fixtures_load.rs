//! Cross-file Java fixture parity: load Tryp_Pig_Bov.revCat.{cseq,canno,csarr,cnlcp}
//! and verify SA size matches CompactFastaSequence size.

use std::io::Cursor;
use std::path::PathBuf;

use model::CompactFastaSequence;
use search::SuffixArray;

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../target/test-classes")
        .join(name)
        .canonicalize()
        .unwrap_or_else(|e| panic!("canonicalize {name}: {e}"))
}

#[test]
fn tryp_pig_bov_revcat_full_set_loads() {
    let cseq = std::fs::read(fixture("Tryp_Pig_Bov.revCat.cseq")).unwrap();
    let canno = std::fs::read(fixture("Tryp_Pig_Bov.revCat.canno")).unwrap();
    let cf = CompactFastaSequence::read_from(
        &mut Cursor::new(&cseq),
        &mut Cursor::new(&canno),
    ).unwrap();

    let csarr = std::fs::read(fixture("Tryp_Pig_Bov.revCat.csarr")).unwrap();
    let cnlcp = std::fs::read(fixture("Tryp_Pig_Bov.revCat.cnlcp")).unwrap();
    let sa = SuffixArray::read_from(
        &mut Cursor::new(&csarr),
        &mut Cursor::new(&cnlcp),
    ).unwrap();

    // 32 = 16 target + 16 decoy.
    assert_eq!(cf.protein_count(), 32);

    // SA length must match CompactFastaSequence size.
    assert_eq!(sa.indices.len() as u64, cf.size,
        "SA indices length {} != .cseq size {}", sa.indices.len(), cf.size);
}
