//! Protein records loaded from a FASTA database. Foundation data
//! consumed by Phase 4b-e search code.

#[derive(Debug, Clone)]
pub struct Protein {
    /// First whitespace-delimited token after the leading `>` on the
    /// header line. Java equivalent: the first column of
    /// `FastaSequence.getAnnotation()`.
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
}
