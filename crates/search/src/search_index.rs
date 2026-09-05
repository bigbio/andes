//! Bundled search database: target+decoy ProteinDb. Consumed by candidate
//! generation.

use crate::decoy::{
    build_search_db, is_decoy_accession_affix, normalize_decoy_prefix, DecoyStrategy,
    DEFAULT_DECOY_SEED,
};
use model::protein::ProteinDb;

#[derive(Debug, Clone)]
pub struct SearchIndex {
    pub db: ProteinDb,
    /// Normalized decoy-accession prefix (no trailing `_`). A protein is a decoy
    /// iff its accession is `"<decoy_prefix>_…"`, independent of FASTA layout
    /// (supports target-only and shuffle modes). The single source of truth for
    /// target/decoy membership, shared with candidate generation.
    pub decoy_prefix: String,
    /// Optional decoy-accession SUFFIX (e.g. `"rev"` for quantms/OpenMS
    /// `<orig>_rev` decoys). When set, a protein is a decoy iff its accession
    /// starts with `"<decoy_prefix>_"` OR ends with this suffix. `None` =
    /// andes-native prefix-only decoys. Set with `--decoy-suffix` (typically
    /// alongside `--decoy-strategy none`) to consume an externally-built
    /// target+decoy database without generating extra decoys.
    pub decoy_suffix: Option<String>,
}

impl SearchIndex {
    /// Target/decoy membership for `accession` under this index's affix config
    /// (prefix needle OR optional suffix). Single entry point used by candidate
    /// generation and the target-protein iterators so RAM and mmap paths agree.
    pub fn is_decoy_accession(&self, accession: &str) -> bool {
        is_decoy_accession_affix(accession, &self.decoy_prefix, self.decoy_suffix.as_deref())
    }

    /// Pipeline: target ProteinDb → reversed decoys → concat target+decoy.
    /// Convenience wrapper for the default (Reverse) strategy.
    pub fn from_target_db(target: &ProteinDb, decoy_prefix: &str) -> Self {
        Self::from_target_db_with_strategy(
            target,
            decoy_prefix,
            None,
            DecoyStrategy::Reverse,
            DEFAULT_DECOY_SEED,
        )
    }

    /// Build the search index for a given decoy `strategy`. `seed` is only used
    /// by [`DecoyStrategy::Shuffle`]. The combined db is target+decoy for
    /// Reverse/Shuffle and the input verbatim for None.
    pub fn from_target_db_with_strategy(
        target: &ProteinDb,
        decoy_prefix: &str,
        decoy_suffix: Option<&str>,
        strategy: DecoyStrategy,
        seed: u64,
    ) -> Self {
        let norm = normalize_decoy_prefix(decoy_prefix);
        let pre = target
            .proteins
            .iter()
            .filter(|p| is_decoy_accession_affix(&p.accession, decoy_prefix, decoy_suffix))
            .count();
        match strategy {
            DecoyStrategy::None => {
                // Target-only / externally-decoyed mode: we generate nothing.
                if pre == 0 {
                    eprintln!(
                        "WARN: --decoy-strategy none and no input proteins carry the '{norm}_' \
                         prefix — the search will have NO decoys, so target/decoy FDR cannot be \
                         estimated. Supply a FASTA that already contains decoys, or use \
                         --decoy-strategy reverse/shuffle."
                    );
                } else {
                    eprintln!(
                        "INFO: --decoy-strategy none — using the input FASTA verbatim; \
                         {pre}/{} proteins carry the '{norm}_' decoy prefix.",
                        target.proteins.len()
                    );
                }
            }
            _ => {
                // Reverse/Shuffle GENERATE a 1:1 decoy per target. If the FASTA
                // ALREADY contains decoys, this doubles them and biases FDR.
                if pre > 0 {
                    eprintln!(
                        "WARN: {pre}/{} input proteins already carry the '{norm}_' decoy prefix — \
                         the FASTA appears to already contain decoys. andes generates its own \
                         decoys, so this DOUBLES decoys and biases FDR. Use --decoy-strategy none, \
                         a target-only FASTA, or a different --decoy-prefix.",
                        target.proteins.len()
                    );
                }
            }
        }
        let db = build_search_db(target, decoy_prefix, strategy, seed);
        Self {
            db,
            decoy_prefix: norm,
            decoy_suffix: decoy_suffix.map(str::to_string),
        }
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

    /// Iterate over target proteins only — every protein whose accession does
    /// NOT carry the decoy prefix. This matches the target/decoy test used in
    /// candidate generation exactly, and (unlike the old positional `len()/2`
    /// split) stays correct under shuffle decoys and target-only/pre-decoyed
    /// inputs. For the default Reverse strategy on a clean target FASTA it is
    /// equivalent to the first half of the combined db.
    pub fn iter_target_proteins(&self) -> impl Iterator<Item = &model::protein::Protein> {
        self.db
            .proteins
            .iter()
            .filter(move |p| !self.is_decoy_accession(&p.accession))
    }

    /// Clone the TARGET-only proteins into a fresh `ProteinDb` — every protein
    /// whose accession does NOT carry the decoy prefix (the same membership test
    /// as [`Self::iter_target_proteins`]). The PTM-refinement cascade uses this
    /// to recover the target-only db from the combined Pass-1 index, then
    /// generates fresh decoys for the scoped Pass-2 search.
    pub fn target_db(&self) -> ProteinDb {
        ProteinDb {
            proteins: self.iter_target_proteins().cloned().collect(),
        }
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
        if needle.is_empty() {
            return true;
        }
        if needle.len() > haystack.len() {
            return false;
        }
        haystack.windows(needle.len()).any(|w| w == needle)
    }

    /// Append `other`'s proteins after `self`'s into one index space, keeping
    /// `self`'s `decoy_prefix` (`other`'s is dropped). Mirrors
    /// [`ProteinDb::concat`]: `other`'s protein at local index `j` lands at
    /// `self.db.len() + j`, so callers must offset `other`'s `protein_index`
    /// values by `self.db.len()` when merging its candidates.
    pub fn concat(&self, other: &SearchIndex) -> SearchIndex {
        SearchIndex {
            db: self.db.concat(&other.db),
            decoy_prefix: self.decoy_prefix.clone(),
            decoy_suffix: self.decoy_suffix.clone(),
        }
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
                Protein {
                    accession: "P1".into(),
                    description: "".into(),
                    sequence: b"MKWV".to_vec(),
                },
                Protein {
                    accession: "P2".into(),
                    description: "".into(),
                    sequence: b"AGCT".to_vec(),
                },
            ],
        };
        let idx = SearchIndex::from_target_db(&target, "XXX");
        assert_eq!(idx.db.len(), 4);
    }

    #[test]
    fn from_target_db_first_half_is_target_second_half_is_decoy() {
        let target = ProteinDb {
            proteins: vec![Protein {
                accession: "P1".into(),
                description: "".into(),
                sequence: b"AB".to_vec(),
            }],
        };
        let idx = SearchIndex::from_target_db(&target, "XXX");
        assert_eq!(idx.db.proteins[0].accession, "P1");
        assert_eq!(idx.db.proteins[1].accession, "XXX_P1");
        assert_eq!(idx.db.proteins[1].sequence, b"BA");
    }

    #[test]
    fn target_db_returns_only_non_decoy_proteins() {
        // Combined db = 2 targets + 2 reverse decoys; target_db() recovers the
        // 2 targets with their original accessions, dropping the prefixed decoys.
        let target = ProteinDb {
            proteins: vec![
                Protein {
                    accession: "P1".into(),
                    description: "".into(),
                    sequence: b"MKWV".to_vec(),
                },
                Protein {
                    accession: "P2".into(),
                    description: "".into(),
                    sequence: b"AGCT".to_vec(),
                },
            ],
        };
        let idx = SearchIndex::from_target_db(&target, "XXX");
        assert_eq!(idx.db.len(), 4, "2 targets + 2 decoys in the combined db");

        let recovered = idx.target_db();
        assert_eq!(recovered.len(), 2, "only the 2 target proteins survive");
        let accs: Vec<&str> = recovered
            .proteins
            .iter()
            .map(|p| p.accession.as_str())
            .collect();
        assert_eq!(accs, vec!["P1", "P2"]);
        // Sequences are preserved verbatim (not reversed).
        assert_eq!(recovered.proteins[0].sequence, b"MKWV");
    }

    // -----------------------------------------------------------------------
    // peptide_has_target_match (all-decoy Label rule)
    // -----------------------------------------------------------------------

    #[test]
    fn peptide_has_target_match_finds_substring() {
        // Target protein: MABCDEFGHIK (as bytes: M=77, A=65, B=66, ...)
        // Use a realistic amino acid sequence the model will accept.
        let target = ProteinDb {
            proteins: vec![Protein {
                accession: "P1".into(),
                description: "".into(),
                sequence: b"MABCDEFGHIK".to_vec(),
            }],
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
            proteins: vec![Protein {
                accession: "P1".into(),
                description: "".into(),
                sequence: b"MABCDEFGHIK".to_vec(),
            }],
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
            proteins: vec![Protein {
                accession: "P1".into(),
                description: "".into(),
                sequence: b"MABCDEFGHIK".to_vec(),
            }],
        };
        let idx = SearchIndex::from_target_db(&target, "XXX");
        assert!(
            idx.peptide_has_target_match(b""),
            "empty peptide is trivially a substring of any target protein"
        );
    }

    #[test]
    fn searchindex_concat_keeps_self_decoy_prefix_and_appends_db() {
        let a = SearchIndex {
            db: ProteinDb { proteins: vec![] },
            decoy_prefix: "XXX".into(),
            decoy_suffix: None,
        };
        let b = SearchIndex {
            db: ProteinDb {
                proteins: vec![model::protein::Protein {
                    accession: "P".into(),
                    description: String::new(),
                    sequence: b"AK".to_vec(),
                }],
            },
            decoy_prefix: "ZZZ".into(),
            decoy_suffix: None,
        };
        let c = a.concat(&b);
        assert_eq!(c.db.len(), 1);
        assert_eq!(c.decoy_prefix, "XXX");
    }
}
