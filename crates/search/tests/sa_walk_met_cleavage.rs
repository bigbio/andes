//! Verify Met-cleaved peptides yield a SEPARATE `DistinctPeptide`
//! (distinguished by `is_protein_n_term`) when their residues happen to
//! match a non-cleaved peptide elsewhere in the database.
//!
//! Fixture: two proteins both contain the tryptic peptide `SAMPLEPEPTIDEK`.
//! - prot1 is M-prefixed: `MSAMPLEPEPTIDEKAGCDR` — Met-cleavage emits
//!   SAMPLEPEPTIDEK at offset 1 with `is_protein_n_term = true` (post-Met
//!   biological N-terminus).
//! - prot2 is `LLSAMPLEPEPTIDEKAGCDR` — SAMPLEPEPTIDEK appears at offset 2
//!   with `is_protein_n_term = false` (interior tryptic peptide).
//!
//! All residues used in the fixture are standard amino acids (no B/J/O/U/X/Z),
//! so the residue-validity gate inside the SA walk admits every length-6+
//! span. NTT is loosened to 0 so SAMPLEPEPTIDEK is admitted from prot2
//! regardless of its non-tryptic pre-flank (L).
//!
//! Contract: residues alone are NOT a sufficient dedup key. The
//! `(residues, is_protein_n_term)` pair must distinguish the two
//! variants, otherwise terminal-mod search space differs between
//! Java and Rust.

mod common;
#[allow(unused_imports)]
use common::*;

use model::{AminoAcidSetBuilder, Protein, ProteinDb};
use search::distinct_peptide::DistinctPeptide;
use search::sa_walk::SaPeptideStream;
use search::{SearchIndex, SearchParams};

fn build_fixture() -> (SearchIndex, SearchParams) {
    let target = ProteinDb {
        proteins: vec![
            Protein {
                accession: "prot1".into(),
                description: "".into(),
                sequence: b"MSAMPLEPEPTIDEKAGCDR".to_vec(),
            },
            Protein {
                accession: "prot2".into(),
                description: "".into(),
                sequence: b"LLSAMPLEPEPTIDEKAGCDR".to_vec(),
            },
        ],
    };
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
    let mut params = SearchParams::default_tryptic(aa_set);
    params.min_length = 6;
    params.max_length = 20;
    params.max_missed_cleavages = 0;
    // Loosen NTT so SAMPLEPEPTIDEK is admitted from prot2 regardless of the
    // pre-flank being X (non-tryptic). This is the SA walk's NTT, not the
    // candidate_gen pass.
    params.num_tolerable_termini = 0;
    (idx, params)
}

#[test]
fn met_cleavage_produces_separate_distinct_peptide() {
    let (idx, params) = build_fixture();
    let peptides: Vec<DistinctPeptide> =
        SaPeptideStream::new(&idx, &params, "XXX").collect();

    // Every emitted peptide should have at least one Position.
    for dp in &peptides {
        assert!(
            !dp.positions.is_empty(),
            "DistinctPeptide with no positions emitted: {:?}",
            std::str::from_utf8(&dp.residues).unwrap_or("<?>")
        );
    }

    let sek: Vec<&DistinctPeptide> = peptides
        .iter()
        .filter(|d| d.residues == b"SAMPLEPEPTIDEK")
        .collect();

    assert!(
        !sek.is_empty(),
        "SAMPLEPEPTIDEK not emitted at all; got {} peptides: {:?}",
        peptides.len(),
        peptides
            .iter()
            .map(|d| (
                std::str::from_utf8(&d.residues).unwrap_or("<?>").to_string(),
                d.positions
                    .iter()
                    .map(|p| (p.protein_index, p.offset, p.is_protein_n_term))
                    .collect::<Vec<_>>()
            ))
            .collect::<Vec<_>>()
    );

    let has_n_term = sek
        .iter()
        .any(|d| d.positions.iter().any(|p| p.is_protein_n_term));
    let has_non_n_term = sek
        .iter()
        .any(|d| d.positions.iter().any(|p| !p.is_protein_n_term));

    assert!(
        has_n_term,
        "Met-cleaved SAMPLEPEPTIDEK (is_protein_n_term=true) must be present; \
         got {} SAMPLEPEPTIDEK entries: {:?}",
        sek.len(),
        sek.iter()
            .map(|d| d
                .positions
                .iter()
                .map(|p| (p.protein_index, p.offset, p.is_protein_n_term))
                .collect::<Vec<_>>())
            .collect::<Vec<_>>()
    );
    assert!(
        has_non_n_term,
        "non-cleaved SAMPLEPEPTIDEK from prot2 (is_protein_n_term=false) must be present"
    );

    // The two variants must NOT collapse into a single DistinctPeptide
    // whose positions vector contains both `is_protein_n_term` values.
    // Either there are >= 2 entries (separate by the n-term axis), or a
    // single entry whose positions all share the same is_protein_n_term.
    let collapsed_into_one_with_mixed = sek.len() == 1
        && sek[0]
            .positions
            .iter()
            .any(|p| p.is_protein_n_term)
        && sek[0]
            .positions
            .iter()
            .any(|p| !p.is_protein_n_term);
    assert!(
        !collapsed_into_one_with_mixed,
        "Met-cleaved + non-cleaved SAMPLEPEPTIDEK were merged into ONE DistinctPeptide; \
         dedup key must include is_protein_n_term"
    );

    // Met-cleaved variant: should be at prot1 (target idx 0), offset 1.
    let met_cleaved_position = sek
        .iter()
        .flat_map(|d| d.positions.iter())
        .find(|p| p.is_protein_n_term);
    let mc = met_cleaved_position.expect("Met-cleaved Position present");
    assert_eq!(mc.offset, 1, "Met-cleaved SAMPLEPEPTIDEK must have offset=1");
}
