//! Protein records loaded from a FASTA database.

#[derive(Debug, Clone)]
pub struct Protein {
    /// First whitespace-delimited token after the leading `>` on the
    /// header line.
    pub accession: String,
    /// Remainder of the header line (after the first whitespace),
    /// trimmed. Empty string if absent.
    pub description: String,
    /// Concatenated sequence lines, uppercase ASCII, whitespace stripped.
    pub sequence: Vec<u8>,
}

impl Protein {
    pub fn len(&self) -> usize { self.sequence.len() }
    pub fn is_empty(&self) -> bool { self.sequence.is_empty() }
}

#[derive(Debug, Clone, Default)]
pub struct ProteinDb {
    pub proteins: Vec<Protein>,
}

impl ProteinDb {
    pub fn new() -> Self { Self::default() }
    pub fn len(&self) -> usize { self.proteins.len() }
    pub fn is_empty(&self) -> bool { self.proteins.is_empty() }
    pub fn iter(&self) -> std::slice::Iter<'_, Protein> { self.proteins.iter() }

    /// Return a new `ProteinDb` containing only the proteins whose positional
    /// index is in `keep`, preserving the original order. Indices in `keep`
    /// that are out of range are ignored. Used by the PTM-refinement cascade to
    /// scope the Pass-2 search to the subset of confidently-identified proteins.
    pub fn subset_by_index(&self, keep: &std::collections::HashSet<usize>) -> ProteinDb {
        let proteins = self
            .proteins
            .iter()
            .enumerate()
            .filter(|(i, _)| keep.contains(i))
            .map(|(_, p)| p.clone())
            .collect();
        ProteinDb { proteins }
    }

    /// Append `other`'s proteins after `self`'s, preserving order. `other`'s
    /// local protein index `j` maps to `self.len() + j` in the result — callers
    /// offset `protein_index` by `self.len()` when merging `other`'s candidates.
    pub fn concat(&self, other: &ProteinDb) -> ProteinDb {
        let mut proteins = self.proteins.clone();
        proteins.extend(other.proteins.iter().cloned());
        ProteinDb { proteins }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_protein() -> Protein {
        Protein {
            accession: "sp|P02769|ALBU_BOVIN".to_string(),
            description: "Serum albumin".to_string(),
            sequence: b"MKWVTFISLL".to_vec(),
        }
    }

    #[test]
    fn protein_len_returns_sequence_length() {
        let p = make_protein();
        assert_eq!(p.len(), 10);
    }

    #[test]
    fn protein_is_empty_false_with_sequence() {
        let p = make_protein();
        assert!(!p.is_empty());
    }

    #[test]
    fn protein_is_empty_true_no_sequence() {
        let p = Protein {
            accession: "x".into(),
            description: "".into(),
            sequence: vec![],
        };
        assert!(p.is_empty());
        assert_eq!(p.len(), 0);
    }

    #[test]
    fn protein_db_default_is_empty() {
        let db = ProteinDb::new();
        assert!(db.is_empty());
        assert_eq!(db.len(), 0);
    }

    #[test]
    fn protein_db_iter() {
        let db = ProteinDb {
            proteins: vec![make_protein(), make_protein()],
        };
        assert_eq!(db.len(), 2);
        let count = db.iter().count();
        assert_eq!(count, 2);
    }

    #[test]
    fn subset_by_index_keeps_selected_in_order() {
        use std::collections::HashSet;
        let p = |acc: &str, seq: &[u8]| Protein {
            accession: acc.into(),
            description: String::new(),
            sequence: seq.to_vec(),
        };
        let db = ProteinDb {
            proteins: vec![p("A", b"AA"), p("B", b"BB"), p("C", b"CC")],
        };
        let keep: HashSet<usize> = [0usize, 2].into_iter().collect();
        let sub = db.subset_by_index(&keep);
        assert_eq!(sub.len(), 2);
        assert_eq!(sub.proteins[0].accession, "A");
        assert_eq!(sub.proteins[1].accession, "C");
    }

    #[test]
    fn subset_by_index_ignores_out_of_range_and_empty() {
        use std::collections::HashSet;
        let db = ProteinDb { proteins: vec![make_protein()] };
        // Index 5 is out of range → ignored, leaving only the in-range pick.
        let keep: HashSet<usize> = [0usize, 5].into_iter().collect();
        assert_eq!(db.subset_by_index(&keep).len(), 1);
        // Empty keep → empty db.
        assert!(db.subset_by_index(&HashSet::new()).is_empty());
    }

    #[test]
    fn concat_appends_proteins_preserving_order() {
        let a = ProteinDb { proteins: vec![
            Protein { accession: "A0".into(), description: String::new(), sequence: b"AK".to_vec() },
            Protein { accession: "A1".into(), description: String::new(), sequence: b"BK".to_vec() },
        ]};
        let b = ProteinDb { proteins: vec![
            Protein { accession: "B0".into(), description: String::new(), sequence: b"CK".to_vec() },
        ]};
        let c = a.concat(&b);
        assert_eq!(c.len(), 3);
        assert_eq!(c.proteins[0].accession, "A0");
        assert_eq!(c.proteins[2].accession, "B0"); // b's index 0 → a.len()+0 = 2
    }
}
