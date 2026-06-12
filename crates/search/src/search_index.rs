//! Bundled search database: target+decoy ProteinDb. Consumed by candidate
//! generation.

use crate::decoy::target_plus_decoy;
use model::protein::ProteinDb;

#[derive(Debug, Clone)]
pub struct SearchIndex {
    pub db: ProteinDb,
}

impl SearchIndex {
    /// Pipeline: target ProteinDb → reverse for decoys → concat target+decoy.
    pub fn from_target_db(target: &ProteinDb, decoy_prefix: &str) -> Self {
        let db = target_plus_decoy(target, decoy_prefix);
        Self { db }
    }

    /// Look up the `Protein` at the given index in the combined target+decoy
    /// database.
    ///
    /// Target proteins occupy `[0, target_count)` and their accessions are the
    /// raw FASTA accessions.  Decoy proteins occupy `[target_count, 2 *
    /// target_count)` and their accessions already carry the decoy prefix (set
    /// by [`target_plus_decoy`]).  Returns `None` when `idx` is out of range.
    pub fn protein_at(&self, idx: usize) -> Option<&model::protein::Protein> {
        self.db.proteins.get(idx)
    }

    /// Iterate over target proteins only (the first half of the combined db).
    ///
    /// `target_plus_decoy` always appends decoys after targets, so target
    /// proteins occupy `[0, total/2)` in `self.db.proteins`.
    pub fn iter_target_proteins(&self) -> impl Iterator<Item = &model::protein::Protein> {
        let target_count = self.db.proteins.len() / 2;
        self.db.proteins[..target_count].iter()
    }

    /// Returns `true` iff `residues` (peptide sequence, no flanking) appears as
    /// a substring in ANY target protein. Used by the PIN writer to compute
    /// Label semantics: Label=-1 only when ALL explaining proteins are decoy.
    ///
    /// Naive scan: O(target_count × len). Acceptable at BSA scale; for real
    /// databases a substring index could accelerate — deferred to a perf pass.
    pub fn peptide_has_target_match(&self, residues: &[u8]) -> bool {
        for prot in self.iter_target_proteins() {
            if Self::contains_subsequence(prot.sequence.as_slice(), residues) {
                return true;
            }
        }
        false
    }

    fn contains_subsequence(haystack: &[u8], needle: &[u8]) -> bool {
        if needle.is_empty() { return true; }
        if needle.len() > haystack.len() { return false; }
        haystack.windows(needle.len()).any(|w| w == needle)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::protein::Protein;

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

    // -----------------------------------------------------------------------
    // peptide_has_target_match (all-decoy Label rule)
    // -----------------------------------------------------------------------

    #[test]
    fn peptide_has_target_match_finds_substring() {
        // Target protein: MABCDEFGHIK (as bytes: M=77, A=65, B=66, ...)
        // Use a realistic amino acid sequence the model will accept.
        let target = ProteinDb {
            proteins: vec![
                Protein {
                    accession: "P1".into(),
                    description: "".into(),
                    sequence: b"MABCDEFGHIK".to_vec(),
                },
            ],
        };
        let idx = SearchIndex::from_target_db(&target, "XXX");
        assert!(
            idx.peptide_has_target_match(b"BCDEF"),
            "BCDEF should be found as a substring of the target protein"
        );
    }

    #[test]
    fn peptide_has_target_match_misses_when_only_in_decoy() {
        // The decoy of MABCDEFGHIK is KIHLGFEDCBAM (reversed).
        // A peptide in the decoy but not the target should return false.
        let target = ProteinDb {
            proteins: vec![
                Protein {
                    accession: "P1".into(),
                    description: "".into(),
                    sequence: b"MABCDEFGHIK".to_vec(),
                },
            ],
        };
        let idx = SearchIndex::from_target_db(&target, "XXX");
        // "KIHLG" appears only in the reversed (decoy) sequence, not in the target.
        assert!(
            !idx.peptide_has_target_match(b"KIHLG"),
            "KIHLG is only in the decoy sequence and should not match any target protein"
        );
    }

    #[test]
    fn peptide_has_target_match_empty_peptide_matches_any_target_protein() {
        // An empty peptide is trivially a substring of any non-empty protein.
        let target = ProteinDb {
            proteins: vec![
                Protein {
                    accession: "P1".into(),
                    description: "".into(),
                    sequence: b"MABCDEFGHIK".to_vec(),
                },
            ],
        };
        let idx = SearchIndex::from_target_db(&target, "XXX");
        assert!(
            idx.peptide_has_target_match(b""),
            "empty peptide is trivially a substring of any target protein"
        );
    }
}
