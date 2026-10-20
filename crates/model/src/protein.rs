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

    /// The region of a sequence-variant entry that holds all of its changed
    /// parts (for most variants the one part), `None` for an ordinary
    /// protein. A span carrying the variant lies within `max_length` of it.
    pub fn variant_region(&self) -> Option<VariantRegion> {
        let parts = self.variant_parts();
        Some(VariantRegion {
            start: parts.iter().map(|v| v.start).min()?,
            end: parts.iter().map(|v| v.end).max()?,
        })
    }

    /// The parts of a sequence-variant entry that read differently from its
    /// canonical sequence; a span carries the variant iff it covers one of
    /// them. Empty for an ordinary protein.
    pub fn variant_parts(&self) -> &[VariantRegion] {
        self.annotation
            .as_deref()
            .map_or(&[], |a| a.variant.as_slice())
    }

    /// The new termini of a mature-chain entry, `None` for an ordinary
    /// protein.
    pub fn chain_ends(&self) -> Option<ChainEnds> {
        self.annotation.as_deref().and_then(|a| a.chain)
    }

    /// The origin of an entry a UniProt XML database added, `None` otherwise.
    pub fn origin(&self) -> Option<&EntryOrigin> {
        self.annotation.as_deref().and_then(|a| a.origin.as_ref())
    }

    /// The single-residue substitutions of a variant entry (empty otherwise).
    pub fn substitutions(&self) -> &[Substitution] {
        self.annotation
            .as_deref()
            .map_or(&[], |a| a.substitutions.as_slice())
    }

    /// The entry's kind; `Reference` for any protein without an origin.
    pub fn entry_kind(&self) -> EntryKind {
        self.origin().map_or(EntryKind::Reference, |o| o.kind)
    }

    /// Whether this is a knockoff entry (see [`EntryOrigin::knockoff_of`]).
    pub fn is_knockoff(&self) -> bool {
        self.origin().is_some_and(|o| o.knockoff_of.is_some())
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
    /// Set when this protein is a sequence variant of another entry: the
    /// parts of it that read differently from that entry, usually one (see
    /// [`VariantRegion::of_edit`]). Only candidate peptides covering one of
    /// them are enumerated from it, so the peptides it shares with the
    /// canonical sequence are not searched twice and every hit on it is a
    /// hit on the variant. Empty for any other protein.
    pub variant: Vec<VariantRegion>,
    /// Set when this protein is a mature chain cut out of a precursor entry:
    /// only the candidate peptides at its new termini are enumerated, which
    /// the precursor does not have as protein termini.
    pub chain: Option<ChainEnds>,
    /// Set on the entries a UniProt XML database adds (isoforms, variants,
    /// chains).
    pub origin: Option<EntryOrigin>,
    /// The single-residue substitutions of a variant entry (equal-length
    /// changes only), sorted by position: where a residue differs from the
    /// sequence the variant was built from, and which residue it replaced.
    pub substitutions: Vec<Substitution>,
    /// The annotated sites of a reference protein as claims (each also in
    /// `mod_sites`), and their knockoff sites, sorted by position. Empty for
    /// alternative entries and decoys, which make no site claims.
    pub site_claims: Vec<SiteClaim>,
}

/// A modification-site claim of a reference protein: modification
/// `mass_delta` at `position`, with the tier of its source. A knockoff site
/// (`knockoff_of`: the real site's position) carries the same modification on
/// a matched residue no annotation names, an almost certainly unmodified
/// position that a placement or localization error fits as well as the real
/// site.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SiteClaim {
    pub position: u32,
    pub mass_delta: f64,
    pub tier: PriorTier,
    pub knockoff_of: Option<u32>,
}

impl ProteinAnnotation {
    /// An annotation with `mod_sites` sorted by position (then mass) and
    /// same-position sites of the same mass (within [`SITE_MASS_TOLERANCE`])
    /// merged, the order [`has_site`] relies on.
    pub fn new(mut mod_sites: Vec<ModSite>, variant: Option<VariantRegion>) -> Self {
        mod_sites.sort_by(|a, b| {
            a.position
                .cmp(&b.position)
                .then(a.mass_delta.total_cmp(&b.mass_delta))
        });
        mod_sites.dedup_by(|b, a| {
            a.position == b.position && (a.mass_delta - b.mass_delta).abs() <= SITE_MASS_TOLERANCE
        });
        Self {
            mod_sites,
            variant: variant.into_iter().collect(),
            chain: None,
            origin: None,
            substitutions: Vec::new(),
            site_claims: Vec::new(),
        }
    }

    /// This annotation with `claims` as its site claims, sorted by position.
    pub fn with_site_claims(mut self, mut claims: Vec<SiteClaim>) -> Self {
        claims.sort_by(|a, b| {
            a.position
                .cmp(&b.position)
                .then(a.mass_delta.total_cmp(&b.mass_delta))
        });
        self.site_claims = claims;
        self
    }

    /// This annotation for a sequence-variant entry whose differing parts are
    /// `parts`.
    pub fn with_variant(mut self, parts: Vec<VariantRegion>) -> Self {
        self.variant = parts;
        self
    }

    /// This annotation for a variant entry with `substitutions`, sorted here.
    pub fn with_substitutions(mut self, mut substitutions: Vec<Substitution>) -> Self {
        substitutions.sort_by_key(|s| s.position);
        self.substitutions = substitutions;
        self
    }

    /// This annotation for an entry of origin `origin`.
    pub fn with_origin(mut self, origin: Option<EntryOrigin>) -> Self {
        self.origin = origin;
        self
    }

    /// This annotation for a mature-chain entry with new termini `chain`.
    pub fn with_chain(mut self, chain: Option<ChainEnds>) -> Self {
        self.chain = chain;
        self
    }
}

/// A residue of a variant entry that replaced `original` at `position` (0-based,
/// in the variant's own coordinates).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Substitution {
    pub position: u32,
    pub original: u8,
}

/// What kind of database entry a protein is. Every FASTA protein, and every
/// UniProt entry's displayed sequence, is `Reference`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum EntryKind {
    #[default]
    Reference,
    Isoform,
    Variant,
    Chain,
}

impl EntryKind {
    pub fn name(self) -> &'static str {
        match self {
            EntryKind::Reference => "reference",
            EntryKind::Isoform => "isoform",
            EntryKind::Variant => "variant",
            EntryKind::Chain => "chain",
        }
    }
}

/// Where a non-reference entry comes from: its kind, the accession of the
/// family it belongs to (`P04637`), the feature that made it (`VAR_005851_D7H`,
/// `P04637-2`, `PRO_0000001058`) and, for a variant, the mass of its replacement
/// residues minus the mass of the residues they replace. A knockoff entry (the
/// same substitution moved to a position no known variant has, an almost
/// certainly false target) names the feature it copies in `knockoff_of`.
/// `source` names the database the claim comes from, and `tier` how likely
/// claims of its kind are to be real, which keeps claims of different
/// reliability from sharing one error estimate.
#[derive(Debug, Clone, PartialEq)]
pub struct EntryOrigin {
    pub kind: EntryKind,
    pub family: Arc<str>,
    pub feature: Arc<str>,
    pub variant_delta_da: f64,
    pub knockoff_of: Option<Arc<str>>,
    pub source: Option<Arc<str>>,
    pub tier: PriorTier,
}

/// How likely a claim is to be real, from what its source says about it.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Default)]
pub enum PriorTier {
    /// Called in the sample itself (its own genotype).
    Sample,
    /// Common in a reference population (allele frequency at least 1%).
    Common,
    /// Curated or experimental evidence, or classed benign.
    Curated,
    /// Anything else: a catalogue entry, weak evidence, no prior given.
    #[default]
    Catalogue,
}

impl PriorTier {
    pub fn name(self) -> &'static str {
        match self {
            PriorTier::Sample => "sample",
            PriorTier::Common => "common",
            PriorTier::Curated => "curated",
            PriorTier::Catalogue => "catalogue",
        }
    }
}

/// Which ends of a mature-chain entry are new termini, i.e. not termini of
/// the precursor (or the residue after an initiator Met) it was cut from.
/// The entry yields only the spans that start at a new N-terminus or end at
/// a new C-terminus; they are protein-terminal there, so protein-terminal
/// mods apply. Every other span reads like the precursor's.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ChainEnds {
    pub n_term: bool,
    pub c_term: bool,
}

impl ChainEnds {
    /// Whether the residue span `[span_start, span_end)` of a chain of
    /// length `len` starts at a new N-terminus or ends at a new C-terminus.
    pub fn covered_by(&self, span_start: u32, span_end: u32, len: u32) -> bool {
        (self.n_term && span_start == 0) || (self.c_term && span_end == len)
    }
}

/// Whether `sites` (sorted by position, see [`ProteinAnnotation::new`]) hold
/// a site at 0-based `position` with a modification of mass `mass_delta`
/// (within [`SITE_MASS_TOLERANCE`]).
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
    /// What the modification is tied to besides its residue.
    pub anchor: SiteAnchor,
}

/// What an annotated modification is tied to: its residue wherever it sits,
/// or the terminus of the chain it modifies (UniProt's "position in the
/// polypeptide", e.g. N-terminal for N-acetylserine, C-terminal for an
/// amidation). Decoy generation keeps the residue of a terminal site at the
/// decoy's own terminus, so the site stays both terminal and on its residue.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SiteAnchor {
    #[default]
    Residue,
    NTerminus,
    CTerminus,
}

/// Where a sequence variant differs from its canonical sequence, in the
/// variant sequence's own 0-based coordinates. `start` is the length of the
/// prefix the two sequences share, and `end` is the variant's length minus
/// the length of the suffix they share (I and L count as the same residue:
/// no mass spectrum tells them apart). Every span ending at or before
/// `start`, or starting at or after `end`, reads like the canonical sequence
/// at the aligned position. A variant that changes residues on both sides of
/// an unchanged run has one region per changed part
/// ([`VariantRegion::of_edit`]).
///
/// - A substitution or insertion has `start < end`: `[start, end)` holds the
///   changed residues.
/// - A deletion has `start == end`: the junction between residues
///   `start - 1` and `start`.
/// - An indel inside a repeat (`AQR` → `AQQR`, `KQQQR` → `KQQR`) can slide
///   along the repeat, which makes `start > end`. Only a span holding all of
///   `[end - 1, start]` then differs from the canonical sequence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct VariantRegion {
    pub start: u32,
    pub end: u32,
}

impl VariantRegion {
    /// The region of the variant sequence `variant` of `canonical`, from the
    /// prefix and suffix the two share.
    pub fn between(canonical: &[u8], variant: &[u8]) -> Self {
        let prefix = canonical
            .iter()
            .zip(variant)
            .take_while(|&(&a, &b)| same_residue(a, b))
            .count();
        let suffix = canonical
            .iter()
            .rev()
            .zip(variant.iter().rev())
            .take_while(|&(&a, &b)| same_residue(a, b))
            .count();
        Self {
            start: prefix as u32,
            end: (variant.len() - suffix) as u32,
        }
    }

    /// The parts of the variant that replaces `canonical[start..end]` by
    /// `replacement` which read differently from `canonical`, in the variant
    /// sequence's coordinates. Usually one, [`VariantRegion::between`] the
    /// two sequences. A replacement that keeps a run of the original
    /// residues between two changes (`KRAPEPTIDEKR` → `ERAPEPTIDEKE`) gets one
    /// region per change, so a peptide inside the kept run (`APEPTIDEK`), which
    /// reads like the canonical one, is not taken for a variant peptide. The
    /// changes come from a longest-common-subsequence alignment of the
    /// original and replacement residues; when a change could slide along a
    /// repeat into its neighbour, the variant keeps the single region instead.
    pub fn of_edit(canonical: &[u8], start: usize, end: usize, replacement: &[u8]) -> Vec<Self> {
        let edited = |from: usize, to: usize, with: &[u8]| -> Vec<u8> {
            let mut s = Vec::with_capacity(canonical.len() + with.len());
            s.extend_from_slice(&canonical[..from]);
            s.extend_from_slice(with);
            s.extend_from_slice(&canonical[to..]);
            s
        };
        let whole = vec![Self::between(canonical, &edited(start, end, replacement))];
        let original = &canonical[start..end];
        let (m, k) = (original.len(), replacement.len());
        // Alignment table cost; UniProt variants are a few residues long.
        if m == 0 || k == 0 || m * k > 1 << 16 {
            return whole;
        }
        // lcs[i][j]: longest common subsequence of original[i..], replacement[j..].
        let mut lcs = vec![vec![0u16; k + 1]; m + 1];
        for i in (0..m).rev() {
            for j in (0..k).rev() {
                lcs[i][j] = if same_residue(original[i], replacement[j]) {
                    lcs[i + 1][j + 1] + 1
                } else {
                    lcs[i + 1][j].max(lcs[i][j + 1])
                };
            }
        }
        // Changed blocks between aligned runs: (original range, replacement range).
        let mut blocks: Vec<(usize, usize, usize, usize)> = Vec::new();
        let (mut i, mut j) = (0, 0);
        let mut open: Option<(usize, usize)> = None;
        while i < m || j < k {
            if i < m
                && j < k
                && same_residue(original[i], replacement[j])
                && lcs[i][j] == lcs[i + 1][j + 1] + 1
            {
                if let Some((oi, oj)) = open.take() {
                    blocks.push((oi, i, oj, j));
                }
                i += 1;
                j += 1;
                continue;
            }
            open.get_or_insert((i, j));
            if j == k || (i < m && lcs[i + 1][j] >= lcs[i][j + 1]) {
                i += 1;
            } else {
                j += 1;
            }
        }
        if let Some((oi, oj)) = open {
            blocks.push((oi, m, oj, k));
        }
        if blocks.len() < 2 {
            return whole;
        }
        let mut parts = Vec::with_capacity(blocks.len());
        for (oi, oe, ri, re) in blocks {
            // This change alone: its region must be exactly its own residues
            // (or its junction), i.e. it must not slide along a repeat.
            let alone = Self::between(
                canonical,
                &edited(start + oi, start + oe, &replacement[ri..re]),
            );
            let own = Self {
                start: (start + oi) as u32,
                end: (start + oi + re - ri) as u32,
            };
            if alone != own {
                return whole;
            }
            // Earlier changes shift it by their length difference (`ri - oi`).
            parts.push(Self {
                start: (start + ri) as u32,
                end: (start + re) as u32,
            });
        }
        parts
    }

    /// Whether the residue span `[span_start, span_end)` carries the variant,
    /// i.e. reads differently from the canonical sequence at the aligned
    /// position: it starts before `end` and ends after `start`. That covers
    /// overlapping a substitution or insertion, holding both residues of a
    /// deletion junction, and holding a whole repeat stretch. Any other span
    /// has a canonical sequence. Its cleavage context can still differ (a
    /// variant to K or R makes the next peptide tryptic); such peptides are
    /// not variant peptides and are left out.
    pub fn covered_by(&self, span_start: u32, span_end: u32) -> bool {
        span_start < self.end && span_end > self.start
    }
}

/// Whether residues `a` and `b` read the same to a mass spectrum: equal, or
/// I and L.
fn same_residue(a: u8, b: u8) -> bool {
    a == b || matches!((a, b), (b'I', b'L') | (b'L', b'I'))
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

    #[test]
    fn entry_kind_defaults_to_reference() {
        let plain = Protein {
            accession: "P1".into(),
            description: String::new(),
            sequence: b"PEPTIDEK".to_vec(),
            annotation: None,
        };
        assert_eq!(plain.entry_kind(), EntryKind::Reference);
        assert!(plain.origin().is_none());
        let origin = EntryOrigin {
            kind: EntryKind::Variant,
            family: Arc::from("P04637"),
            feature: Arc::from("VAR_005851_D7H"),
            variant_delta_da: 0.984016,
            knockoff_of: None,
            source: None,
            tier: PriorTier::Catalogue,
        };
        let var = Protein {
            annotation: Some(Arc::new(
                ProteinAnnotation::new(Vec::new(), None).with_origin(Some(origin.clone())),
            )),
            ..plain.clone()
        };
        assert_eq!(var.entry_kind(), EntryKind::Variant);
        assert_eq!(var.origin(), Some(&origin));
        assert_eq!(EntryKind::Chain.name(), "chain");
    }

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
    fn variant_region_in_a_repeat_keeps_only_differing_spans() {
        // Duplication AQR -> AQQR: only spans holding both Qs differ.
        let dup = VariantRegion::between(b"AQR", b"AQQR");
        assert_eq!(dup, VariantRegion { start: 2, end: 2 });
        assert!(dup.covered_by(1, 3), "QQ");
        assert!(!dup.covered_by(2, 4), "QR is canonical");
        assert!(!dup.covered_by(0, 2), "AQ is canonical");
        // Deletion KQQQR -> KQQR: only the whole KQQR differs.
        let del = VariantRegion::between(b"KQQQR", b"KQQR");
        assert_eq!(del, VariantRegion { start: 3, end: 1 });
        assert!(del.covered_by(0, 4));
        assert!(!del.covered_by(1, 4), "QQR is canonical");
        assert!(!del.covered_by(0, 3), "KQQ is canonical");
        // Outside repeats the region is the replaced residues, as UniProt
        // gives them.
        assert_eq!(
            VariantRegion::between(b"MKDEF", b"MKHEF"),
            VariantRegion { start: 2, end: 3 }
        );
        assert_eq!(
            VariantRegion::between(b"MKDEF", b"MKEF"),
            VariantRegion { start: 2, end: 2 }
        );
        // An insertion whose first residue repeats the original is narrowed
        // to what is new: P -> PRAA.
        assert_eq!(
            VariantRegion::between(b"MPEK", b"MPRAAEK"),
            VariantRegion { start: 2, end: 5 }
        );
        // I and L read the same: only K -> E differs in IK -> LE.
        assert_eq!(
            VariantRegion::between(b"AIKDE", b"ALEDE"),
            VariantRegion { start: 2, end: 3 }
        );
    }

    #[test]
    fn edits_keeping_a_run_get_one_region_per_change() {
        // KRAPEPTIDEKR -> ERAPEPTIDEKE keeps RAPEPTIDEK: K2E and R13E.
        let canonical = b"MKRAPEPTIDEKRGH";
        let parts = VariantRegion::of_edit(canonical, 1, 13, b"ERAPEPTIDEKE");
        assert_eq!(
            parts,
            [
                VariantRegion { start: 1, end: 2 },
                VariantRegion { start: 12, end: 13 }
            ]
        );
        // APEPTIDEK (3..12) reads like the canonical peptide.
        assert!(!parts.iter().any(|v| v.covered_by(3, 12)));
        assert!(parts.iter().any(|v| v.covered_by(1, 12)));
        // Kept runs between an insertion and a deletion: KRAPEPTIDEKR ->
        // KRWAPEPTIDEK (W inserted, final R deleted).
        assert_eq!(
            VariantRegion::of_edit(canonical, 1, 13, b"KRWAPEPTIDEK"),
            [
                VariantRegion { start: 3, end: 4 },
                VariantRegion { start: 13, end: 13 }
            ]
        );
        // A single change is the region between the sequences.
        assert_eq!(
            VariantRegion::of_edit(b"MKDEF", 2, 3, b"H"),
            [VariantRegion { start: 2, end: 3 }]
        );
        assert_eq!(
            VariantRegion::of_edit(b"MKDEF", 2, 3, b""),
            [VariantRegion { start: 2, end: 2 }]
        );
        // A change that could slide along a repeat (one Q of QQ deleted)
        // keeps the single region holding everything.
        let canonical = b"KQQRAPEPTIDEKRG";
        let variant = b"KQRAPEPTIDEKEG";
        assert_eq!(
            VariantRegion::of_edit(canonical, 1, 14, &variant[1..13]),
            [VariantRegion::between(canonical, variant)]
        );
    }

    #[test]
    fn variant_region_holds_all_parts() {
        let p = Protein {
            accession: "v".to_string(),
            description: String::new(),
            sequence: b"MERAPEPTIDEKEGH".to_vec(),
            annotation: Some(Arc::new(
                ProteinAnnotation::new(Vec::new(), None).with_variant(vec![
                    VariantRegion { start: 1, end: 2 },
                    VariantRegion { start: 12, end: 13 },
                ]),
            )),
        };
        assert_eq!(p.variant_parts().len(), 2);
        assert_eq!(
            p.variant_region(),
            Some(VariantRegion { start: 1, end: 13 })
        );
    }

    #[test]
    fn annotation_accessors_and_site_lookup() {
        let mut p = make_protein();
        assert!(p.mod_sites().is_empty());
        assert_eq!(p.variant_region(), None);
        p.annotation = Some(Arc::new(ProteinAnnotation::new(
            vec![
                ModSite {
                    position: 2,
                    mass_delta: 79.966331,
                    anchor: SiteAnchor::Residue,
                },
                ModSite {
                    position: 2,
                    mass_delta: 15.994915,
                    anchor: SiteAnchor::Residue,
                },
                ModSite {
                    position: 2,
                    mass_delta: 79.9663,
                    anchor: SiteAnchor::Residue,
                },
            ],
            Some(VariantRegion { start: 1, end: 2 }),
        )));
        // Sorted by mass, the near-duplicate phospho merged.
        assert_eq!(p.mod_sites().len(), 2);
        assert!(p.mod_sites()[0].mass_delta < p.mod_sites()[1].mass_delta);
        assert_eq!(p.variant_region(), Some(VariantRegion { start: 1, end: 2 }));
        let ann = p.annotation.as_deref().unwrap();
        assert!(has_site(&ann.mod_sites, 2, 79.966));
        assert!(has_site(&ann.mod_sites, 2, 15.9949));
        assert!(!has_site(&ann.mod_sites, 2, 42.010565));
        assert!(!has_site(&ann.mod_sites, 3, 79.966331));
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
