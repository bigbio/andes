//! Decoy database generation. Mirrors Java
//! `edu.ucsd.msjava.msdbsearch.ReverseDB`.

use crate::protein::{Protein, ProteinDb};

/// Default decoy prefix matching Java MSGFPlus.DEFAULT_DECOY_PROTEIN_PREFIX.
pub const DEFAULT_DECOY_PREFIX: &str = "XXX";

/// Reverse each protein's sequence and prepend `<prefix>_` to its
/// accession. `prefix` is normalized: trailing `_`s stripped; empty
/// prefix → `DEFAULT_DECOY_PREFIX`.
pub fn reverse_db(db: &ProteinDb, prefix: &str) -> ProteinDb {
    let normalized = normalize_prefix(prefix);
    let proteins = db.proteins.iter().map(|p| Protein {
        accession: format!("{}_{}", normalized, p.accession),
        description: p.description.clone(),
        sequence: p.sequence.iter().rev().copied().collect(),
    }).collect();
    ProteinDb { proteins }
}

/// Concatenate target + decoy. Equivalent to Java's `concat=true` mode.
pub fn target_plus_decoy(target: &ProteinDb, prefix: &str) -> ProteinDb {
    let decoy = reverse_db(target, prefix);
    let mut proteins = target.proteins.clone();
    proteins.extend(decoy.proteins);
    ProteinDb { proteins }
}

fn normalize_prefix(prefix: &str) -> String {
    let trimmed = prefix.trim().trim_end_matches('_');
    if trimmed.is_empty() {
        DEFAULT_DECOY_PREFIX.to_string()
    } else {
        trimmed.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_db(proteins: &[(&str, &[u8])]) -> ProteinDb {
        ProteinDb {
            proteins: proteins.iter().map(|(acc, seq)| Protein {
                accession: acc.to_string(),
                description: String::new(),
                sequence: seq.to_vec(),
            }).collect(),
        }
    }

    #[test]
    fn reverse_db_reverses_sequences() {
        let db = make_db(&[("P1", b"MKWV"), ("P2", b"AGCT")]);
        let decoy = reverse_db(&db, "XXX");
        assert_eq!(decoy.len(), 2);
        assert_eq!(decoy.proteins[0].sequence, b"VWKM");
        assert_eq!(decoy.proteins[1].sequence, b"TCGA");
    }

    #[test]
    fn reverse_db_prepends_prefix() {
        let db = make_db(&[("P1", b"AB")]);
        let decoy = reverse_db(&db, "XXX");
        assert_eq!(decoy.proteins[0].accession, "XXX_P1");
    }

    #[test]
    fn reverse_db_strips_trailing_underscores_in_prefix() {
        let db = make_db(&[("P1", b"AB")]);
        let decoy = reverse_db(&db, "XXX_");
        assert_eq!(decoy.proteins[0].accession, "XXX_P1");
    }

    #[test]
    fn reverse_db_empty_prefix_uses_default() {
        let db = make_db(&[("P1", b"AB")]);
        let decoy = reverse_db(&db, "");
        assert_eq!(decoy.proteins[0].accession, "XXX_P1");
    }

    #[test]
    fn reverse_db_preserves_description() {
        let mut db = make_db(&[("P1", b"AB")]);
        db.proteins[0].description = "Some description".into();
        let decoy = reverse_db(&db, "XXX");
        assert_eq!(decoy.proteins[0].description, "Some description");
    }

    #[test]
    fn target_plus_decoy_concats() {
        let target = make_db(&[("P1", b"AB"), ("P2", b"CD")]);
        let combined = target_plus_decoy(&target, "XXX");
        assert_eq!(combined.len(), 4);
        assert_eq!(combined.proteins[0].accession, "P1");
        assert_eq!(combined.proteins[2].accession, "XXX_P1");
        assert_eq!(combined.proteins[2].sequence, b"BA");
    }
}
