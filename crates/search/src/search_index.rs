use crate::decoy::{
    build_search_db_for, is_decoy_accession_affix, normalize_decoy_prefix, DecoyStrategy,
    Proteases, DEFAULT_DECOY_SEED,
};
use model::protein::ProteinDb;

#[derive(Debug, Clone)]
pub struct SearchIndex {
    pub db: ProteinDb,
    pub has_alternative_entries: bool,
    pub decoy_prefix: String,
    pub decoy_suffix: Option<String>,
}

impl SearchIndex {
    pub fn is_decoy_accession(&self, accession: &str) -> bool {
        is_decoy_accession_affix(accession, &self.decoy_prefix, self.decoy_suffix.as_deref())
    }

    pub fn from_target_db(target: &ProteinDb, decoy_prefix: &str) -> Self {
        Self::from_target_db_with_strategy(
            target,
            decoy_prefix,
            None,
            DecoyStrategy::default(),
            DEFAULT_DECOY_SEED,
        )
    }

    pub fn from_target_db_with_strategy(
        target: &ProteinDb,
        decoy_prefix: &str,
        decoy_suffix: Option<&str>,
        strategy: DecoyStrategy,
        seed: u64,
    ) -> Self {
        Self::from_target_db_for(
            target,
            decoy_prefix,
            decoy_suffix,
            strategy,
            seed,
            &Proteases::default(),
        )
    }

    pub fn from_target_db_for(
        target: &ProteinDb,
        decoy_prefix: &str,
        decoy_suffix: Option<&str>,
        strategy: DecoyStrategy,
        seed: u64,
        proteases: &Proteases,
    ) -> Self {
        let norm = normalize_decoy_prefix(decoy_prefix);
        let pre = target
            .proteins
            .iter()
            .filter(|p| is_decoy_accession_affix(&p.accession, decoy_prefix, decoy_suffix))
            .count();
        match strategy {
            DecoyStrategy::None => {
                if pre == 0 {
                    log::warn!(
                        target: "database",
                        "--decoy-strategy none and no input proteins carry the '{norm}_' \
                         prefix — the search will have NO decoys, so target/decoy FDR cannot be \
                         estimated. Supply a FASTA that already contains decoys, or use \
                         --decoy-strategy reverse/shuffle."
                    );
                } else {
                    log::info!(
                        target: "database",
                        "--decoy-strategy none — using the input FASTA verbatim; \
                         {pre}/{} proteins carry the '{norm}_' decoy prefix.",
                        target.proteins.len()
                    );
                }
            }
            _ => {
                if pre > 0 {
                    log::warn!(
                        target: "database",
                        "{pre}/{} input proteins already carry the '{norm}_' decoy prefix — \
                         the FASTA appears to already contain decoys. andes generates its own \
                         decoys, so this DOUBLES decoys and biases FDR. Use --decoy-strategy none, \
                         a target-only FASTA, or a different --decoy-prefix.",
                        target.proteins.len()
                    );
                }
            }
        }
        let db = build_search_db_for(target, decoy_prefix, strategy, seed, proteases);
        Self {
            has_alternative_entries: holds_alternative_entries(&db),
            db,
            decoy_prefix: norm,
            decoy_suffix: decoy_suffix.map(str::to_string),
        }
    }

    pub fn protein_at(&self, idx: usize) -> Option<&model::protein::Protein> {
        self.db.proteins.get(idx)
    }

    pub fn iter_target_proteins(&self) -> impl Iterator<Item = &model::protein::Protein> {
        self.db
            .proteins
            .iter()
            .filter(move |p| !self.is_decoy_accession(&p.accession))
    }

    pub fn target_db(&self) -> ProteinDb {
        ProteinDb {
            proteins: self.iter_target_proteins().cloned().collect(),
        }
    }

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

    pub fn concat(&self, other: &SearchIndex) -> SearchIndex {
        SearchIndex {
            has_alternative_entries: self.has_alternative_entries || other.has_alternative_entries,
            db: self.db.concat(&other.db),
            decoy_prefix: self.decoy_prefix.clone(),
            decoy_suffix: self.decoy_suffix.clone(),
        }
    }
}

pub fn holds_alternative_entries(db: &ProteinDb) -> bool {
    db.proteins
        .iter()
        .any(|p| p.entry_kind() != model::EntryKind::Reference)
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
                    annotation: None,
                },
                Protein {
                    accession: "P2".into(),
                    description: "".into(),
                    sequence: b"AGCT".to_vec(),
                    annotation: None,
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
                sequence: b"PEPTIDEK".to_vec(),
                annotation: None,
            }],
        };
        let idx = SearchIndex::from_target_db_with_strategy(
            &target,
            "XXX",
            None,
            DecoyStrategy::Reverse,
            DEFAULT_DECOY_SEED,
        );
        assert_eq!(idx.db.proteins[0].accession, "P1");
        assert_eq!(idx.db.proteins[1].accession, "XXX_P1");
        assert_eq!(idx.db.proteins[1].sequence, b"EDITPEPK");
    }

    #[test]
    fn target_db_returns_only_non_decoy_proteins() {
        let target = ProteinDb {
            proteins: vec![
                Protein {
                    accession: "P1".into(),
                    description: "".into(),
                    sequence: b"MKWV".to_vec(),
                    annotation: None,
                },
                Protein {
                    accession: "P2".into(),
                    description: "".into(),
                    sequence: b"AGCT".to_vec(),
                    annotation: None,
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
        assert_eq!(recovered.proteins[0].sequence, b"MKWV");
    }

    #[test]
    fn peptide_has_target_match_finds_substring() {
        let target = ProteinDb {
            proteins: vec![Protein {
                accession: "P1".into(),
                description: "".into(),
                sequence: b"MABCDEFGHIK".to_vec(),
                annotation: None,
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
        let target = ProteinDb {
            proteins: vec![Protein {
                accession: "P1".into(),
                description: "".into(),
                sequence: b"MABCDEFGHIK".to_vec(),
                annotation: None,
            }],
        };
        let idx = SearchIndex::from_target_db_with_strategy(
            &target,
            "XXX",
            None,
            DecoyStrategy::Reverse,
            DEFAULT_DECOY_SEED,
        );
        assert_eq!(idx.db.proteins[1].sequence, b"MIHGFEDCBAK");
        assert!(
            !idx.peptide_has_target_match(b"IHGFE"),
            "IHGFE is only in the decoy sequence and should not match any target protein"
        );
    }

    #[test]
    fn peptide_has_target_match_empty_peptide_matches_any_target_protein() {
        let target = ProteinDb {
            proteins: vec![Protein {
                accession: "P1".into(),
                description: "".into(),
                sequence: b"MABCDEFGHIK".to_vec(),
                annotation: None,
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
            has_alternative_entries: false,
            decoy_prefix: "XXX".into(),
            decoy_suffix: None,
        };
        let b = SearchIndex {
            db: ProteinDb {
                proteins: vec![model::protein::Protein {
                    accession: "P".into(),
                    description: String::new(),
                    sequence: b"AK".to_vec(),
                    annotation: None,
                }],
            },
            has_alternative_entries: false,
            decoy_prefix: "ZZZ".into(),
            decoy_suffix: None,
        };
        let c = a.concat(&b);
        assert_eq!(c.db.len(), 1);
        assert_eq!(c.decoy_prefix, "XXX");
    }
}
