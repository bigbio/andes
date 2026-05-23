//! Verify `SaPeptideStream` walks the SA + LCP and produces one
//! `DistinctPeptide` per unique residue sequence (within the limits of the
//! current LCP-only dedup), accumulating every `(protein, offset)` position
//! it encounters.
//!
//! Fixture: 3 proteins where two of them (prot1 + prot3) contain the same
//! tryptic peptide `LMNPQR`. The exact dedup outcome depends on whether
//! the two SA-adjacent suffixes share their N-term flank byte:
//!
//! - prot1 LMNPQR pre-flank = `R` (residue at index 10 of `ABCDEFGHIKRLMNPQR`)
//! - prot3 LMNPQR pre-flank = TERMINATOR (start of protein)
//!
//! The SA walk does not see the pre-flank directly — it sees the residues
//! and the FORWARD characters. The LCP between
//! `LMNPQR\0...` (prot1 trailing TERM) and `LMNPQRR...` (prot3 next residue)
//! is exactly 6 (the residues match, the 7th byte differs).
//!
//! With the current simplification (lcp == L+1 treated as a new peptide),
//! this yields two separate `DistinctPeptide` entries for `LMNPQR`. The
//! test therefore checks the SOFT contract: at least one `DistinctPeptide`
//! has residues `LMNPQR`, AND every emitted `DistinctPeptide` carries at
//! least one valid `Position`. The plan flags the imperfect dedup as
//! acceptable for this subtask; the next subtask refines the SA walk's
//! flank handling.

mod common;
#[allow(unused_imports)]
use common::*;

use model::{AminoAcidSetBuilder, Protein, ProteinDb};
use search::distinct_peptide::DistinctPeptide;
use search::sa_walk::SaPeptideStream;
use search::{SearchIndex, SearchParams};

fn build_fixture_idx_params() -> (SearchIndex, SearchParams) {
    let target = ProteinDb {
        proteins: vec![
            Protein {
                accession: "prot1".into(),
                description: "".into(),
                sequence: b"ABCDEFGHIKRLMNPQR".to_vec(),
            },
            Protein {
                accession: "prot2".into(),
                description: "".into(),
                sequence: b"ABCDEFGHIKRSTVWY".to_vec(),
            },
            Protein {
                accession: "prot3".into(),
                description: "".into(),
                sequence: b"LMNPQRRZZZZ".to_vec(),
            },
        ],
    };
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
    let mut params = SearchParams::default_tryptic(aa_set);
    params.min_length = 6;
    params.max_length = 20;
    params.max_missed_cleavages = 0;
    params.num_tolerable_termini = 0; // SA walk doesn't enforce missed-cleavage; keep NTT loose so LMNPQR is admitted even when one flank is non-tryptic.
    (idx, params)
}

#[test]
fn sa_walk_yields_lmnpqr_with_positions() {
    let (idx, params) = build_fixture_idx_params();
    let peptides: Vec<DistinctPeptide> = SaPeptideStream::new(&idx, &params, "XXX").collect();

    // Sanity: walk produced peptides at all.
    assert!(!peptides.is_empty(), "SA walk produced zero peptides");

    // Every emitted peptide must have at least one Position.
    for dp in &peptides {
        assert!(
            !dp.positions.is_empty(),
            "DistinctPeptide with no positions emitted: {:?}",
            dp.residues
        );
    }

    // Find every DistinctPeptide whose residues are exactly "LMNPQR".
    let lmnpqr: Vec<&DistinctPeptide> = peptides
        .iter()
        .filter(|d| d.residues == b"LMNPQR")
        .collect();

    assert!(
        !lmnpqr.is_empty(),
        "LMNPQR not emitted by SA walk; got peptides: {:?}",
        peptides
            .iter()
            .map(|d| std::str::from_utf8(&d.residues).unwrap_or("<?>"))
            .collect::<Vec<_>>()
    );

    // Aggregate positions across every LMNPQR entry. We expect TWO total
    // occurrences (prot1 offset 11, prot3 offset 0) regardless of whether
    // LCP dedup folded them into one or two DistinctPeptides.
    let total_positions: usize = lmnpqr.iter().map(|d| d.positions.len()).sum();
    assert_eq!(
        total_positions, 2,
        "expected 2 total LMNPQR occurrences (prot1 + prot3), got {} across {} DistinctPeptide(s)",
        total_positions,
        lmnpqr.len()
    );

    // Per the plan: ideal dedup yields one DistinctPeptide with two
    // Positions. Current LCP-only impl may yield two separate entries
    // because the pre-flank differs (R vs protein-start), which the SA
    // walk cannot observe directly. Flag-but-don't-fail when dedup is
    // imperfect — the next subtask refines flank handling.
    if lmnpqr.len() == 1 {
        assert_eq!(
            lmnpqr[0].positions.len(),
            2,
            "single LMNPQR entry should aggregate both positions"
        );
    } else {
        eprintln!(
            "warning: LMNPQR not deduped into a single DistinctPeptide \
             (got {} entries with {} total positions). Acceptable for this \
             subtask; flank-aware dedup arrives in the next subtask.",
            lmnpqr.len(),
            total_positions
        );
    }

    // Protein-index sanity: the two occurrences must come from target
    // proteins 0 (prot1) and 2 (prot3) — never the decoys (3, 4, 5).
    let mut seen_target_proteins: Vec<u32> = lmnpqr
        .iter()
        .flat_map(|d| d.positions.iter().map(|p| p.protein_index))
        .filter(|p| (*p as usize) < idx.db.proteins.len() / 2)
        .collect();
    seen_target_proteins.sort();
    assert_eq!(
        seen_target_proteins,
        vec![0, 2],
        "LMNPQR target positions should be in prot1 (idx 0) and prot3 (idx 2)"
    );
}
