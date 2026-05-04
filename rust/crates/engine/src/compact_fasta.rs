//! Concatenated-byte representation of a ProteinDb with separator
//! bytes. Mirrors Java
//! `edu.ucsd.msjava.msdbsearch.CompactFastaSequence`. Used as input
//! to suffix-array construction. Phase 4c Task 2 ships only the
//! in-memory builder; file I/O lands in Task 4.

use crate::protein::ProteinDb;

/// Java's protein delimiter byte. Java uses `_` (0x5F) per
/// `Constants.PROTEIN_DELIMITER`. Phase 4c Task 4 will verify this
/// against the Java source.
pub const SEPARATOR: u8 = b'_';

/// End-of-sequence terminator byte. Java uses 0 per
/// `Constants.SEQUENCE_TERMINATOR`. Verified against existing fixtures
/// in Task 4.
pub const TERMINATOR: u8 = 0;

#[derive(Debug, Clone)]
pub struct CompactFastaSequence {
    /// `[SEP] <protein0> [SEP] <protein1> [SEP] ... [TERM]`
    pub sequence: Vec<u8>,
    pub annotations: Vec<ProteinAnnotation>,
    pub size: u64,
}

#[derive(Debug, Clone)]
pub struct ProteinAnnotation {
    /// Offset into `sequence` of this protein's first residue byte.
    pub start: u64,
    pub accession: String,
    pub description: String,
}

impl CompactFastaSequence {
    pub fn from_protein_db(db: &ProteinDb) -> Self {
        if db.proteins.is_empty() {
            // Empty DB: no sequence content, no annotations.
            return Self {
                sequence: Vec::new(),
                annotations: Vec::new(),
                size: 0,
            };
        }

        let mut sequence = Vec::with_capacity(
            db.proteins.iter().map(|p| p.sequence.len() + 1).sum::<usize>() + 2,
        );
        let mut annotations = Vec::with_capacity(db.proteins.len());

        // Lead with a separator byte (matches Java's CompactFastaSequence layout).
        sequence.push(SEPARATOR);
        for p in &db.proteins {
            let start = sequence.len() as u64;
            sequence.extend_from_slice(&p.sequence);
            sequence.push(SEPARATOR);
            annotations.push(ProteinAnnotation {
                start,
                accession: p.accession.clone(),
                description: p.description.clone(),
            });
        }
        // Replace the final SEPARATOR with TERMINATOR (matches Java end-of-stream).
        if let Some(last) = sequence.last_mut() {
            *last = TERMINATOR;
        }

        let size = sequence.len() as u64;
        Self {
            sequence,
            annotations,
            size,
        }
    }

    pub fn protein_count(&self) -> usize {
        self.annotations.len()
    }

    /// Binary-search the annotation array for the protein containing
    /// position `pos`. Returns `None` for positions before the first
    /// protein.
    pub fn protein_index_at(&self, pos: u64) -> Option<usize> {
        if self.annotations.is_empty() {
            return None;
        }
        match self.annotations.binary_search_by(|a| a.start.cmp(&pos)) {
            Ok(idx) => Some(idx),
            Err(0) => None,
            Err(idx) => Some(idx - 1),
        }
    }
}

#[derive(thiserror::Error, Debug)]
pub enum CompactFastaError {
    #[error("I/O error: {source}")]
    Io {
        #[from]
        source: std::io::Error,
    },
    #[error("file I/O not yet implemented (Phase 4c/Task 4 stub)")]
    NotYetImplemented,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::protein::{Protein, ProteinDb};

    fn make_db(proteins: &[(&str, &[u8])]) -> ProteinDb {
        ProteinDb {
            proteins: proteins
                .iter()
                .map(|(acc, seq)| Protein {
                    accession: acc.to_string(),
                    description: String::new(),
                    sequence: seq.to_vec(),
                })
                .collect(),
        }
    }

    #[test]
    fn empty_db_produces_zero_proteins() {
        let db = ProteinDb::new();
        let cf = CompactFastaSequence::from_protein_db(&db);
        assert_eq!(cf.protein_count(), 0);
        assert_eq!(cf.annotations.len(), 0);
    }

    #[test]
    fn single_protein_sequence_is_preserved() {
        let db = make_db(&[("P1", b"MKWV")]);
        let cf = CompactFastaSequence::from_protein_db(&db);
        assert_eq!(cf.protein_count(), 1);
        assert_eq!(cf.annotations[0].accession, "P1");
        let start = cf.annotations[0].start as usize;
        assert_eq!(&cf.sequence[start..start + 4], b"MKWV");
    }

    #[test]
    fn two_proteins_have_separator_between() {
        let db = make_db(&[("P1", b"AB"), ("P2", b"CD")]);
        let cf = CompactFastaSequence::from_protein_db(&db);
        assert_eq!(cf.protein_count(), 2);
        let start1 = cf.annotations[0].start as usize;
        let start2 = cf.annotations[1].start as usize;
        // Each protein 2 bytes; at least one separator byte between them.
        assert!(
            start2 > start1 + 2,
            "expected separator between proteins; start1={start1}, start2={start2}"
        );
        // The byte between protein 1's end and protein 2's start should be SEPARATOR.
        assert_eq!(cf.sequence[start1 + 2], SEPARATOR);
    }

    #[test]
    fn protein_index_at_returns_correct_index() {
        let db = make_db(&[("P1", b"ABC"), ("P2", b"DEF"), ("P3", b"GHI")]);
        let cf = CompactFastaSequence::from_protein_db(&db);
        let p1_start = cf.annotations[0].start;
        assert_eq!(cf.protein_index_at(p1_start), Some(0));
        let p2_start = cf.annotations[1].start;
        assert_eq!(cf.protein_index_at(p2_start), Some(1));
        let p3_start = cf.annotations[2].start;
        assert_eq!(cf.protein_index_at(p3_start), Some(2));
    }

    #[test]
    fn description_preserved() {
        let mut db = make_db(&[("P1", b"AB")]);
        db.proteins[0].description = "test description".into();
        let cf = CompactFastaSequence::from_protein_db(&db);
        assert_eq!(cf.annotations[0].description, "test description");
    }

    #[test]
    fn size_matches_sequence_length() {
        let db = make_db(&[("P1", b"AB"), ("P2", b"CD")]);
        let cf = CompactFastaSequence::from_protein_db(&db);
        assert_eq!(cf.size, cf.sequence.len() as u64);
    }
}
