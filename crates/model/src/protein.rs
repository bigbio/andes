//! Protein records loaded from a FASTA or UniProt XML database.

use std::sync::Arc;

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
    /// Database annotation of this sequence: annotated modification sites and,
    /// for a sequence-variant entry, the residues the variant changed. Only a
    /// UniProt XML database carries one; `None` for FASTA proteins. Shared, so
    /// cloning a protein does not copy it.
    pub annotation: Option<Arc<ProteinAnnotation>>,
}

impl Protein {
    pub fn len(&self) -> usize {
        self.sequence.len()
    }
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    /// Annotated modification sites, sorted by position (empty without an
    /// annotation).
    pub fn mod_sites(&self) -> &[ModSite] {
        self.annotation
            .as_deref()
            .map_or(&[], |a| a.mod_sites.as_slice())
    }

    /// The residues a sequence-variant entry changed, `None` for an ordinary
    /// protein.
    pub fn variant_region(&self) -> Option<VariantRegion> {
        self.annotation.as_deref().and_then(|a| a.variant)
    }
}

/// Tolerance (Da) within which an annotated site's mass and a modification's
/// mass are taken to be the same modification. Tight enough to keep phospho
/// (79.9663) and sulfo (79.9568) apart, loose enough for a mass typed to two
/// decimals.
pub const SITE_MASS_TOLERANCE: f64 = 0.005;

/// Annotation a database supplies for one protein sequence.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ProteinAnnotation {
    /// Annotated modification sites, sorted by `position`. A variable
    /// modification declared `sites=annotated` is only placed on these.
    pub mod_sites: Vec<ModSite>,
    /// Set when this protein is a sequence variant of another entry: only
    /// candidate peptides that carry the variant are enumerated from it, so
    /// the peptides it shares with the canonical sequence are not searched
    /// twice and every hit on it is a hit on the variant.
    pub variant: Option<VariantRegion>,
}

impl ProteinAnnotation {
    /// Whether a site at 0-based `position` is annotated with a modification
    /// of mass `mass_delta` (within [`SITE_MASS_TOLERANCE`]).
    pub fn has_site(&self, position: u32, mass_delta: f64) -> bool {
        has_site(&self.mod_sites, position, mass_delta)
    }
}

/// [`ProteinAnnotation::has_site`] on a sorted site list.
pub fn has_site(sites: &[ModSite], position: u32, mass_delta: f64) -> bool {
    let first = sites.partition_point(|s| s.position < position);
    sites[first..]
        .iter()
        .take_while(|s| s.position == position)
        .any(|s| (s.mass_delta - mass_delta).abs() <= SITE_MASS_TOLERANCE)
}

/// One annotated modification site.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ModSite {
    /// 0-based residue position in the protein sequence.
    pub position: u32,
    /// Monoisotopic mass change of the annotated modification (Da).
    pub mass_delta: f64,
}

/// The residues a sequence variant changed, in the variant sequence's own
/// 0-based coordinates: `[start, end)` holds the replacement residues. A
/// deletion replaces with nothing, so `start == end` and the variant is the
/// junction between residues `start - 1` and `start`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct VariantRegion {
    pub start: u32,
    pub end: u32,
}

impl VariantRegion {
    /// Whether the residue span `[span_start, span_end)` carries the variant:
    /// it overlaps the replacement residues or, for a deletion, holds both
    /// residues on either side of the junction. Any other span reads the same
    /// as the canonical sequence.
    pub fn covered_by(&self, span_start: u32, span_end: u32) -> bool {
        if self.start < self.end {
            span_start < self.end && span_end > self.start
        } else {
            span_start < self.start && span_end > self.start
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct ProteinDb {
    pub proteins: Vec<Protein>,
}

impl ProteinDb {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn len(&self) -> usize {
        self.proteins.len()
    }
    pub fn is_empty(&self) -> bool {
        self.proteins.is_empty()
    }
    pub fn iter(&self) -> std::slice::Iter<'_, Protein> {
        self.proteins.iter()
    }

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
            annotation: None,
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
            annotation: None,
        };
        assert!(p.is_empty());
        assert_eq!(p.len(), 0);
    }

    #[test]
    fn variant_region_covers_spans_carrying_the_change() {
        // Substitution of residue 5.
        let sub = VariantRegion { start: 5, end: 6 };
        assert!(sub.covered_by(5, 6));
        assert!(sub.covered_by(0, 9));
        assert!(sub.covered_by(3, 6));
        assert!(!sub.covered_by(0, 5), "ends just before the variant");
        assert!(!sub.covered_by(6, 9), "starts just after the variant");
        // Deletion: the junction between residues 4 and 5.
        let del = VariantRegion { start: 5, end: 5 };
        assert!(del.covered_by(4, 6));
        assert!(del.covered_by(0, 9));
        assert!(
            !del.covered_by(0, 5),
            "only the residue before the junction"
        );
        assert!(!del.covered_by(5, 9), "only the residue after the junction");
    }

    #[test]
    fn annotation_accessors_and_site_lookup() {
        let mut p = make_protein();
        assert!(p.mod_sites().is_empty());
        assert_eq!(p.variant_region(), None);
        p.annotation = Some(Arc::new(ProteinAnnotation {
            mod_sites: vec![
                ModSite {
                    position: 2,
                    mass_delta: 15.994915,
                },
                ModSite {
                    position: 2,
                    mass_delta: 79.966331,
                },
            ],
            variant: Some(VariantRegion { start: 1, end: 2 }),
        }));
        assert_eq!(p.mod_sites().len(), 2);
        assert_eq!(p.variant_region(), Some(VariantRegion { start: 1, end: 2 }));
        let ann = p.annotation.as_deref().unwrap();
        assert!(ann.has_site(2, 79.966));
        assert!(ann.has_site(2, 15.9949));
        assert!(!ann.has_site(2, 42.010565));
        assert!(!ann.has_site(3, 79.966331));
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
            annotation: None,
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
        let db = ProteinDb {
            proteins: vec![make_protein()],
        };
        // Index 5 is out of range → ignored, leaving only the in-range pick.
        let keep: HashSet<usize> = [0usize, 5].into_iter().collect();
        assert_eq!(db.subset_by_index(&keep).len(), 1);
        // Empty keep → empty db.
        assert!(db.subset_by_index(&HashSet::new()).is_empty());
    }

    #[test]
    fn concat_appends_proteins_preserving_order() {
        let a = ProteinDb {
            proteins: vec![
                Protein {
                    accession: "A0".into(),
                    description: String::new(),
                    sequence: b"AK".to_vec(),
                    annotation: None,
                },
                Protein {
                    accession: "A1".into(),
                    description: String::new(),
                    sequence: b"BK".to_vec(),
                    annotation: None,
                },
            ],
        };
        let b = ProteinDb {
            proteins: vec![Protein {
                accession: "B0".into(),
                description: String::new(),
                sequence: b"CK".to_vec(),
                annotation: None,
            }],
        };
        let c = a.concat(&b);
        assert_eq!(c.len(), 3);
        assert_eq!(c.proteins[0].accession, "A0");
        assert_eq!(c.proteins[2].accession, "B0"); // b's index 0 → a.len()+0 = 2
    }
}
