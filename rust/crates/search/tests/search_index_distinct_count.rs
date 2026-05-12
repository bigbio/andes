//! Verifies `SearchIndex::num_distinct_peptides_at_length` returns the count
//! of distinct residue sequences (no mods, no flanking, target+decoy combined)
//! enumerated by `enumerate_candidates` for each peptide length.
//!
//! Test fixture: 3 synthetic proteins with controlled overlap to exercise
//! per-length deduplication across both target and decoy proteomes.
//!
//! NOTE: The plan's draft fixture used non-standard residues (B, X, Y, Z) and
//! only counted target peptides. We use a fully-standard-AA fixture and
//! account for decoy contributions in the expected counts.

mod common;
#[allow(unused_imports)]
use common::*;

use model::{AminoAcidSetBuilder, Protein, ProteinDb};
use search::{SearchIndex, SearchParams};

/// Build a fixture with 3 proteins designed to share specific tryptic
/// peptides at known lengths. All sequences use only standard residues.
///
/// Target tryptic peptides (Trypsin, missed=0):
///   prot1 = AGTLPDQVIK + LMNPQR        → "AGTLPDQVIK" (10), "LMNPQR" (6)
///   prot2 = AGTLPDQVIK + STVCYHK       → "AGTLPDQVIK" (10), "STVCYHK" (7)
///   prot3 = LMNPQR     + WWWK          → "LMNPQR" (6),     "WWWK" (4)
///
/// Decoy tryptic peptides (reversed sequences):
///   prot1 decoy "RQPNMLKIVQDPLTGA" → "QPNMLK" (6), "IVQDPLTGA" (9)
///   prot2 decoy "KHYCVTSKIVQDPLTGA" → "HYCVTSK" (7), "IVQDPLTGA" (9)
///   prot3 decoy "KWWWRQPNML"          → "WWWR" (4), "QPNML" (5)
///
/// Distinct counts per length (target ∪ decoy, deduplicated):
///   len  4: {WWWK, WWWR}                    → 2
///   len  5: {QPNML}                         → 1
///   len  6: {LMNPQR, QPNMLK}                → 2  (LMNPQR shared p1+p3 → counted once)
///   len  7: {STVCYHK, HYCVTSK}              → 2
///   len  9: {IVQDPLTGA}                     → 1  (shared by both decoys → counted once)
///   len 10: {AGTLPDQVIK}                    → 1  (shared p1+p2 → counted once)
fn build_fixture() -> (SearchIndex, SearchParams) {
    let target = ProteinDb {
        proteins: vec![
            Protein {
                accession: "prot1".into(),
                description: "".into(),
                sequence: b"AGTLPDQVIKLMNPQR".to_vec(),
            },
            Protein {
                accession: "prot2".into(),
                description: "".into(),
                sequence: b"AGTLPDQVIKSTVCYHK".to_vec(),
            },
            Protein {
                accession: "prot3".into(),
                description: "".into(),
                sequence: b"LMNPQRWWWK".to_vec(),
            },
        ],
    };
    let idx = SearchIndex::from_target_db(&target, "XXX");

    let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
    let mut params = SearchParams::default_tryptic(aa_set);
    params.min_length = 4;
    params.max_length = 12;
    params.max_missed_cleavages = 0;
    params.max_variable_mods_per_peptide = 0;
    params.num_tolerable_termini = 2;

    let idx = idx.with_distinct_peptide_counts(&params, "XXX");
    (idx, params)
}

#[test]
fn distinct_count_at_length_10_dedups_shared_target_peptide() {
    let (idx, _) = build_fixture();
    // "AGTLPDQVIK" appears in prot1 + prot2 targets; counted once.
    assert_eq!(idx.num_distinct_peptides_at_length(10), 1);
}

#[test]
fn distinct_count_at_length_6_includes_decoy() {
    let (idx, _) = build_fixture();
    // Targets: "LMNPQR" (shared p1+p3, 1 distinct).
    // Decoys: "QPNMLK" (prot1 decoy).
    // Total distinct: 2.
    assert_eq!(idx.num_distinct_peptides_at_length(6), 2);
}

#[test]
fn distinct_count_at_length_7_includes_decoy() {
    let (idx, _) = build_fixture();
    // Targets: "STVCYHK" (prot2). Decoys: "HYCVTSK" (prot2 decoy). Distinct: 2.
    assert_eq!(idx.num_distinct_peptides_at_length(7), 2);
}

#[test]
fn distinct_count_at_length_4_includes_decoy() {
    let (idx, _) = build_fixture();
    // Targets: "WWWK" (prot3). Decoys: "WWWR" (prot3 decoy). Distinct: 2.
    assert_eq!(idx.num_distinct_peptides_at_length(4), 2);
}

#[test]
fn distinct_count_at_length_9_dedups_shared_decoy_peptide() {
    let (idx, _) = build_fixture();
    // Decoys: "IVQDPLTGA" appears in both prot1 + prot2 decoys; counted once.
    assert_eq!(idx.num_distinct_peptides_at_length(9), 1);
}

#[test]
fn distinct_count_at_unseen_length_is_zero() {
    let (idx, _) = build_fixture();
    // No peptide in the fixture has length 99.
    assert_eq!(idx.num_distinct_peptides_at_length(99), 0);
}

#[test]
fn distinct_count_at_length_below_min_length_is_zero() {
    let (idx, _) = build_fixture();
    // min_length=4, so length=1 is excluded from enumeration.
    assert_eq!(idx.num_distinct_peptides_at_length(1), 0);
}
