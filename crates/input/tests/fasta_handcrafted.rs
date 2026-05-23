//! Handcrafted FASTA strings exercising parser edge cases.

use std::io::Cursor;
use input::{FastaParseError, FastaReader, Protein};

fn parse_all(s: &str) -> Vec<Result<Protein, FastaParseError>> {
    FastaReader::new(Cursor::new(s)).collect()
}

fn parse_ok(s: &str) -> Vec<Protein> {
    parse_all(s).into_iter().map(|r| r.unwrap()).collect()
}

#[test]
fn empty_input_emits_nothing() {
    let v = parse_ok("");
    assert!(v.is_empty());
}

#[test]
fn single_protein_single_sequence_line() {
    let fa = ">P1 description here\nMKWVTFISLL\n";
    let v = parse_ok(fa);
    assert_eq!(v.len(), 1);
    assert_eq!(v[0].accession, "P1");
    assert_eq!(v[0].description, "description here");
    assert_eq!(v[0].sequence, b"MKWVTFISLL");
}

#[test]
fn single_protein_multi_line_sequence() {
    let fa = ">P1\n\
              MKWVTFISLL\n\
              LFSSAYSRGV\n";
    let v = parse_ok(fa);
    assert_eq!(v.len(), 1);
    assert_eq!(v[0].sequence, b"MKWVTFISLLLFSSAYSRGV");
}

#[test]
fn multiple_proteins() {
    let fa = ">P1 first\n\
              MKWV\n\
              >P2 second\n\
              TFIS\n\
              >P3 third\n\
              LLLF\n";
    let v = parse_ok(fa);
    assert_eq!(v.len(), 3);
    assert_eq!(v[0].accession, "P1");
    assert_eq!(v[1].accession, "P2");
    assert_eq!(v[2].accession, "P3");
}

#[test]
fn semicolon_comments_skipped() {
    let fa = "; this is a comment\n\
              >P1\n\
              MKWV\n\
              ; another comment\n\
              TFIS\n";
    let v = parse_ok(fa);
    assert_eq!(v.len(), 1);
    assert_eq!(v[0].sequence, b"MKWVTFIS");
}

#[test]
fn blank_lines_tolerated() {
    let fa = "\n\
              >P1\n\
              \n\
              MKWV\n\
              \n\
              \n\
              TFIS\n";
    let v = parse_ok(fa);
    assert_eq!(v.len(), 1);
    assert_eq!(v[0].sequence, b"MKWVTFIS");
}

#[test]
fn lowercase_residues_uppercased() {
    let fa = ">P1\nmKwVtFiSlL\n";
    let v = parse_ok(fa);
    assert_eq!(v[0].sequence, b"MKWVTFISLL");
}

#[test]
fn whitespace_inside_sequence_stripped() {
    let fa = ">P1\nM K W V\nT F I S\n";
    let v = parse_ok(fa);
    assert_eq!(v[0].sequence, b"MKWVTFIS");
}

#[test]
fn header_no_description() {
    let fa = ">P1\nMKWV\n";
    let v = parse_ok(fa);
    assert_eq!(v[0].accession, "P1");
    assert_eq!(v[0].description, "");
}

#[test]
fn header_multi_word_description() {
    let fa = ">sp|P02769|ALBU_BOVIN Serum albumin OS=Bos taurus\nMKWV\n";
    let v = parse_ok(fa);
    assert_eq!(v[0].accession, "sp|P02769|ALBU_BOVIN");
    assert_eq!(v[0].description, "Serum albumin OS=Bos taurus");
}

#[test]
fn empty_accession_errors() {
    let fa = ">\nMKWV\n";
    let err = parse_all(fa).into_iter().next().unwrap().unwrap_err();
    assert!(matches!(err, FastaParseError::EmptyAccession { .. }));
}

#[test]
fn orphan_sequence_errors() {
    let fa = "MKWV\n>P1\nTFIS\n";
    let err = parse_all(fa).into_iter().next().unwrap().unwrap_err();
    assert!(matches!(err, FastaParseError::OrphanSequence { .. }));
}

#[test]
fn last_protein_terminated_by_eof() {
    let fa = ">P1\nMKWV\n>P2\nTFIS";  // no trailing newline
    let v = parse_ok(fa);
    assert_eq!(v.len(), 2);
    assert_eq!(v[1].sequence, b"TFIS");
}
