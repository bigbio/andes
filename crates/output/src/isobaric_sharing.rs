//! Variant peptides that a reference protein holds once I and L are taken as
//! equal. I and L have the same mass and fragment the same, so such a peptide
//! is as much an identification of the reference protein as of the variant,
//! yet the search groups proteins by exact sequence and reports it as
//! variant-specific (a variant V→L can spell a paralog's `…VIDEVR` as
//! `…VLDEVR`). [`share_isobaric_references`] adds those references to the
//! row, so it counts as the identification it is.

use std::collections::HashMap;

use model::EntryKind;
use rayon::prelude::*;
use search::SearchIndex;

use crate::percolator::PercolatorPsm;
use crate::proteoform_table::bare_sequence;

/// References a variant-class row names at most, besides its variants.
const MAX_ADDED: usize = 8;

/// For every target row at q ≤ `q_max` whose proteins are all variant entries
/// (real or knockoff, so both are treated alike), add the target Reference
/// proteins whose sequence holds its peptide with I and L taken as equal.
/// Rows above `q_max` are left alone: they are never reported. Returns how
/// many rows gained a reference.
pub fn share_isobaric_references(
    targets: &mut HashMap<String, PercolatorPsm>,
    idx: &SearchIndex,
    q_max: f64,
) -> usize {
    let kinds: HashMap<&str, EntryKind> = idx
        .db
        .proteins
        .iter()
        .map(|p| (p.accession.as_str(), p.entry_kind()))
        .collect();
    let variant_only = |p: &PercolatorPsm| {
        p.proteins
            .split(';')
            .filter(|a| !a.is_empty())
            .all(|a| kinds.get(a) == Some(&EntryKind::Variant))
    };
    let il = |s: &[u8]| -> Vec<u8> {
        s.iter()
            .map(|&c| if c == b'I' { b'L' } else { c })
            .collect()
    };
    let mut peptides: Vec<String> = targets
        .values()
        .filter(|p| p.q_value <= q_max && variant_only(p))
        .map(|p| bare_sequence(&p.peptide))
        .collect();
    peptides.sort();
    peptides.dedup();
    if peptides.is_empty() {
        return 0;
    }
    let references: Vec<(&str, Vec<u8>)> = idx
        .iter_target_proteins()
        .filter(|p| p.entry_kind() == EntryKind::Reference)
        .map(|p| (p.accession.as_str(), il(&p.sequence)))
        .collect();
    let holders: HashMap<String, Vec<&str>> = peptides
        .par_iter()
        .filter_map(|pep| {
            let needle = il(pep.as_bytes());
            let found: Vec<&str> = references
                .iter()
                .filter(|(_, seq)| contains(seq, &needle))
                .map(|(acc, _)| *acc)
                .take(MAX_ADDED)
                .collect();
            (!found.is_empty()).then(|| (pep.clone(), found))
        })
        .collect();
    let mut changed = 0;
    for p in targets.values_mut() {
        if p.q_value > q_max || !variant_only(p) {
            continue;
        }
        if let Some(found) = holders.get(&bare_sequence(&p.peptide)) {
            for acc in found {
                p.proteins.push(';');
                p.proteins.push_str(acc);
            }
            changed += 1;
        }
    }
    changed
}

fn contains(haystack: &[u8], needle: &[u8]) -> bool {
    !needle.is_empty() && memchr::memmem::find(haystack, needle).is_some()
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::protein::{EntryOrigin, Protein, ProteinAnnotation, ProteinDb};
    use std::sync::Arc;

    fn protein(acc: &str, seq: &[u8], kind: EntryKind) -> Protein {
        let annotation = (kind != EntryKind::Reference).then(|| {
            Arc::new(
                ProteinAnnotation::new(Vec::new(), None).with_origin(Some(EntryOrigin {
                    kind,
                    family: Arc::from("P1"),
                    feature: Arc::from("V"),
                    variant_delta_da: 0.0,
                    knockoff_of: None,
                    source: None,
                    tier: model::PriorTier::Catalogue,
                })),
            )
        });
        Protein {
            accession: acc.into(),
            description: String::new(),
            sequence: seq.to_vec(),
            annotation,
        }
    }

    fn psm(id: &str, q: f64, peptide: &str, proteins: &str) -> (String, PercolatorPsm) {
        (
            id.to_string(),
            PercolatorPsm {
                psm_id: id.to_string(),
                q_value: q,
                pep: 0.0,
                score: Some(1.0),
                peptide: peptide.into(),
                proteins: proteins.into(),
            },
        )
    }

    #[test]
    fn a_variant_peptide_a_reference_holds_up_to_i_and_l_is_shared() {
        let target = ProteinDb {
            proteins: vec![
                protein("sp|P1|A", b"MKAVFVDLEPTVVDEVRK", EntryKind::Reference),
                protein("sp|P2|B", b"MKAVFVDLEPTVIDEVRK", EntryKind::Reference),
                protein("sp|P1_V|A", b"MKAVFVDLEPTVLDEVRK", EntryKind::Variant),
                protein("sp|P1_W|A", b"MKAVFVDLEPTVWDEVRK", EntryKind::Variant),
            ],
        };
        let idx = SearchIndex::from_target_db(&target, "XXX");
        let mut t: HashMap<_, _> = [
            psm("v", 0.001, "K.AVFVDLEPTVLDEVR.K", "sp|P1_V|A"),
            psm("w", 0.001, "K.AVFVDLEPTVWDEVR.K", "sp|P1_W|A"),
            psm("late", 0.5, "K.AVFVDLEPTVLDEVR.K", "sp|P1_V|A"),
        ]
        .into_iter()
        .collect();
        assert_eq!(share_isobaric_references(&mut t, &idx, 0.05), 1);
        assert_eq!(t["v"].proteins, "sp|P1_V|A;sp|P2|B");
        assert_eq!(
            t["w"].proteins, "sp|P1_W|A",
            "a truly variant peptide is left alone"
        );
        assert_eq!(
            t["late"].proteins, "sp|P1_V|A",
            "rows never reported are left alone"
        );
    }
}
