//! Bundled search database: target+decoy ProteinDb, CompactFastaSequence,
//! and SuffixArray. Output of Phase 4b+4c, input of Phase 4d.

use crate::compact_fasta::{CompactFastaError, CompactFastaSequence};
use crate::decoy::target_plus_decoy;
use crate::protein::ProteinDb;
use crate::suffix_array::{SuffixArray, SuffixArrayError};

#[derive(Debug, Clone)]
pub struct SearchIndex {
    pub db: ProteinDb,
    pub compact: CompactFastaSequence,
    pub sa: SuffixArray,
}

impl SearchIndex {
    /// Pipeline: target ProteinDb → reverse for decoys → concat target+decoy
    /// → CompactFastaSequence → SA + LCP.
    pub fn from_target_db(target: &ProteinDb, decoy_prefix: &str) -> Self {
        let db = target_plus_decoy(target, decoy_prefix);
        let compact = CompactFastaSequence::from_protein_db(&db);
        let sa = SuffixArray::build(&compact);
        Self { db, compact, sa }
    }

    /// Look up the `Protein` at the given index in the combined target+decoy
    /// database.
    ///
    /// Target proteins occupy `[0, target_count)` and their accessions are the
    /// raw FASTA accessions.  Decoy proteins occupy `[target_count, 2 *
    /// target_count)` and their accessions already carry the decoy prefix (set
    /// by [`target_plus_decoy`]).  Returns `None` when `idx` is out of range.
    pub fn protein_at(&self, idx: usize) -> Option<&crate::protein::Protein> {
        self.db.proteins.get(idx)
    }
}

#[derive(thiserror::Error, Debug)]
pub enum SearchIndexError {
    #[error("compact fasta error: {0}")]
    CompactFasta(#[from] CompactFastaError),
    #[error("suffix array error: {0}")]
    SuffixArray(#[from] SuffixArrayError),
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::protein::Protein;

    #[test]
    fn from_target_db_doubles_protein_count() {
        let target = ProteinDb {
            proteins: vec![
                Protein { accession: "P1".into(), description: "".into(), sequence: b"MKWV".to_vec() },
                Protein { accession: "P2".into(), description: "".into(), sequence: b"AGCT".to_vec() },
            ],
        };
        let idx = SearchIndex::from_target_db(&target, "XXX");
        assert_eq!(idx.db.len(), 4);
        assert_eq!(idx.sa.indices.len(), idx.compact.size as usize);
    }

    #[test]
    fn from_target_db_first_half_is_target_second_half_is_decoy() {
        let target = ProteinDb {
            proteins: vec![
                Protein { accession: "P1".into(), description: "".into(), sequence: b"AB".to_vec() },
            ],
        };
        let idx = SearchIndex::from_target_db(&target, "XXX");
        assert_eq!(idx.db.proteins[0].accession, "P1");
        assert_eq!(idx.db.proteins[1].accession, "XXX_P1");
        assert_eq!(idx.db.proteins[1].sequence, b"BA");
    }
}
