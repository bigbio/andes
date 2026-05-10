//! Bundled search database: target+decoy ProteinDb, CompactFastaSequence,
//! and SuffixArray. Consumed by candidate generation.

use std::collections::{HashMap, HashSet};

use model::compact_fasta::{CompactFastaError, CompactFastaSequence};
use crate::candidate_gen::enumerate_candidates;
use crate::decoy::target_plus_decoy;
use model::protein::ProteinDb;
use crate::search_params::SearchParams;
use crate::suffix_array::{SuffixArray, SuffixArrayError};

#[derive(Debug, Clone)]
pub struct SearchIndex {
    pub db: ProteinDb,
    pub compact: CompactFastaSequence,
    pub sa: SuffixArray,
    /// Number of distinct residue sequences enumerated per length.
    /// Populated by [`SearchIndex::with_distinct_peptide_counts`]; mirrors Java's
    /// `CompactSuffixArray.getNumDistinctPeptides(int length)`. Consumed by the
    /// `e_value` computation in `match_engine` to replace the queue-derived
    /// proxy with a per-length distinct-peptide count over the full target+decoy
    /// candidate set.
    distinct_peptide_counts: HashMap<usize, usize>,
}

impl SearchIndex {
    /// Pipeline: target ProteinDb → reverse for decoys → concat target+decoy
    /// → CompactFastaSequence → SA + LCP.
    ///
    /// `distinct_peptide_counts` is left empty; call
    /// [`SearchIndex::with_distinct_peptide_counts`] to populate it once
    /// `SearchParams` are known.
    pub fn from_target_db(target: &ProteinDb, decoy_prefix: &str) -> Self {
        let db = target_plus_decoy(target, decoy_prefix);
        let compact = CompactFastaSequence::from_protein_db(&db);
        let sa = SuffixArray::build(&compact);
        Self {
            db,
            compact,
            sa,
            distinct_peptide_counts: HashMap::new(),
        }
    }

    /// Walk every candidate emitted by [`enumerate_candidates`] for `params`
    /// and `decoy_prefix`, then store the count of distinct residue sequences
    /// per peptide length. Returns the index with the populated map.
    ///
    /// Mirrors Java `CompactSuffixArray.computeNumDistinctPeptides`
    /// (`src/main/java/edu/ucsd/msjava/msdbsearch/CompactSuffixArray.java:198`):
    /// counts distinct prefixes of length `l` across the entire suffix array
    /// (target + decoy combined, modulo the still-open mod-context divergence
    /// tracked in `docs/parity-analysis/known-divergences.md`).
    ///
    /// Distinct identity is the residue byte sequence with no mods and no
    /// flanking residues — matching Java's prefix comparison over residue bytes.
    /// Two candidates with identical residues but different mod variants count
    /// as one (residues only); two candidates that differ in flanking context
    /// also count as one.
    ///
    /// Cost is bounded by the candidate-set walk (a HashSet of residue
    /// `Vec<u8>` per length); the HashSets are dropped after the count is
    /// frozen. See `docs/superpowers/specs/2026-05-10-evalue-search-index-design.md`
    /// for the memory analysis at PXD001819 scale.
    pub fn with_distinct_peptide_counts(
        mut self,
        params: &SearchParams,
        decoy_prefix: &str,
    ) -> Self {
        let mut seen_per_length: HashMap<usize, HashSet<Vec<u8>>> = HashMap::new();
        for cand in enumerate_candidates(&self, params, decoy_prefix) {
            let residues: Vec<u8> = cand.peptide.residues.iter().map(|aa| aa.residue).collect();
            seen_per_length
                .entry(residues.len())
                .or_default()
                .insert(residues);
        }
        self.distinct_peptide_counts = seen_per_length
            .into_iter()
            .map(|(len, set)| (len, set.len()))
            .collect();
        self
    }

    /// Number of distinct residue sequences (no mods, no flanking) of length
    /// `len` enumerated during candidate generation. Returns `0` for unseen
    /// lengths (including any length never populated because
    /// [`SearchIndex::with_distinct_peptide_counts`] was not called).
    ///
    /// Mirrors Java `CompactSuffixArray.getNumDistinctPeptides(int length)`.
    pub fn num_distinct_peptides_at_length(&self, len: usize) -> usize {
        self.distinct_peptide_counts.get(&len).copied().unwrap_or(0)
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
    /// Java-faithful Label semantics (Java: DirectPinWriter.java:188-191 —
    /// Label=-1 only when ALL explaining proteins are decoy).
    ///
    /// Naive scan: O(target_count × len). Acceptable at BSA scale; for real
    /// databases the suffix array could accelerate — deferred to a perf pass.
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
