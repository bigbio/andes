//! Streaming UniProtKB XML reader (`uniprot_sprot.xml[.gz]`, per-proteome
//! downloads, `rest.uniprot.org/uniprotkb/<acc>.xml`).
//!
//! A UniProt XML database searches like the equivalent UniProt FASTA: every
//! entry becomes one [`Protein`] whose accession and description are the
//! UniProt FASTA header (`sp|P02769|ALBU_BOVIN`, `Albumin OS=Bos taurus
//! OX=9913 GN=ALB PE=1 SV=4`), so a sequences-only search is the FASTA search.
//! On top of that the XML carries annotation FASTA drops, which this reader
//! keeps:
//!
//! * **Annotated modification sites** (`modified residue`, `lipid
//!   moiety-binding region`, `glycosylation site`, and ubiquitin-like `cross-link`s
//!   that leave a Gly-Gly remnant) become [`ModSite`]s on the protein's
//!   [`ProteinAnnotation`]. They only change a search when a variable mod is
//!   declared `sites=annotated` (see `model::Modification::annotated_site_mass`),
//!   which then applies at the annotated positions alone.
//! * **Natural sequence variants** (`sequence variant`), with
//!   [`UniProtOptions::variants`], become extra target proteins, one per
//!   variant, whose accession names the change (`sp|P04637_VAR_005851_D7H|P53_HUMAN`).
//!   Candidate generation enumerates only the peptides of such a protein that
//!   carry the variant ([`VariantRegion`]), so every hit on it is a hit on
//!   the variant and the canonical peptides are not searched twice.
//!
//! Both are filtered by UniProt's evidence codes ([`EvidenceLevel`]): expanding
//! the search on weak annotation adds search space, and so error, without
//! adding truth.
//!
//! Features located on another isoform (`<location sequence="P04637-2">`) are
//! in that isoform's coordinates and are not applied to the canonical sequence;
//! isoforms themselves are not reconstructed (as in OpenMS' `UniProtXMLHandler`,
//! which this parser follows for element scoping).

use std::collections::BTreeMap;
use std::io::BufRead;
use std::sync::Arc;

use model::{ModSite, Protein, ProteinAnnotation, ProteinDb, SiteAnchor, VariantRegion};
use quick_xml::events::{BytesStart, Event};
use quick_xml::Reader;

use crate::uniprot_ptm::{
    clean_description, ptm_mass, GLYGLY_CROSSLINK_PREFIX, GLYGLY_MASS, GLYGLY_MODIFIERS,
};

/// How well supported a UniProt feature must be to expand the search.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Default)]
pub enum EvidenceLevel {
    /// Experimental evidence: ECO:0000269 (experimental evidence used in
    /// manual assertion) and ECO:0007744 / ECO:0007829 (combined experimental
    /// and computational evidence, e.g. large-scale proteomics studies).
    #[default]
    Experimental,
    /// Also manually curated assertions without direct experimental support:
    /// combinatorial evidence (ECO:0000244), by similarity (ECO:0000250),
    /// sequence model (ECO:0000255), author statement (ECO:0000303), curator
    /// inference (ECO:0000305) and imported information (ECO:0000312), each
    /// used in manual assertion.
    Curated,
    /// Every feature, including automatic annotation and features with no
    /// evidence at all (e.g. most dbSNP-only variants).
    Any,
}

impl EvidenceLevel {
    /// Parse `experimental`, `curated` or `any` (case-insensitive).
    pub fn from_name(s: &str) -> Option<Self> {
        match s.trim().to_ascii_lowercase().as_str() {
            "experimental" => Some(Self::Experimental),
            "curated" => Some(Self::Curated),
            "any" => Some(Self::Any),
            _ => None,
        }
    }

    /// The name [`EvidenceLevel::from_name`] parses.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Experimental => "experimental",
            Self::Curated => "curated",
            Self::Any => "any",
        }
    }

    /// The level an evidence code (`ECO:…`) satisfies, `None` for codes only
    /// [`EvidenceLevel::Any`] admits (automatic annotation).
    fn of_code(code: &str) -> Option<Self> {
        match code {
            "ECO:0000269" | "ECO:0007744" | "ECO:0007829" => Some(Self::Experimental),
            "ECO:0000244" | "ECO:0000250" | "ECO:0000255" | "ECO:0000303" | "ECO:0000305"
            | "ECO:0000312" => Some(Self::Curated),
            _ => None,
        }
    }

    /// Whether a feature with these evidence codes reaches this level: one
    /// code at or above it suffices. A feature without evidence reaches only
    /// [`EvidenceLevel::Any`].
    pub fn admits<'a>(self, mut codes: impl Iterator<Item = &'a str>) -> bool {
        self == Self::Any || codes.any(|c| Self::of_code(c).is_some_and(|l| l <= self))
    }
}

/// What a UniProt XML load keeps beyond the canonical sequences.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct UniProtOptions {
    /// Add one target protein per natural sequence variant (off by default:
    /// it grows the search space, so it must be measured before being used).
    pub variants: bool,
    /// Minimum evidence for a variant or a modification site to be used.
    pub evidence: EvidenceLevel,
}

/// Counts from a UniProt XML load, for the run log.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct UniProtLoadStats {
    /// `<entry>` elements read (= canonical proteins).
    pub entries: usize,
    /// Variant proteins added.
    pub variant_proteins: usize,
    /// Annotated modification sites kept (over canonical proteins).
    pub mod_sites: usize,
    /// Variant and modification features left out for weak evidence.
    pub below_evidence: usize,
    /// Modification features whose name has no single known mass, by name.
    pub unknown_modifications: BTreeMap<String, usize>,
    /// Variant or modification features that do not fit the canonical
    /// sequence (unknown position, original residues or modified residue not
    /// matching, replacement that is not a residue string), and variants a
    /// search cannot tell from the canonical sequence (I↔L only).
    pub unusable_features: usize,
    /// Features located on a non-canonical isoform (not applied).
    pub isoform_features: usize,
}

/// One `<feature>` of an entry.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct UniProtFeature {
    /// `@type`, e.g. `modified residue`, `sequence variant`.
    pub kind: String,
    /// `@id`, e.g. `VAR_005851` (absent on many features).
    pub id: String,
    /// `@description`, raw.
    pub description: String,
    /// `@evidence` keys, resolved against the entry's `<evidence>` elements.
    pub evidence_keys: Vec<u32>,
    /// `<location @sequence>`: the isoform the coordinates refer to; empty
    /// for the canonical sequence.
    pub location_sequence: String,
    /// 1-based inclusive residue range (`<position>` gives a one-residue
    /// range); `None` when either end is unknown or absent.
    pub location: Option<(u32, u32)>,
    /// `<original>` residues (sequence variants), empty when absent.
    pub original: String,
    /// `<variation>` residues; several for a multi-allelic variant, none for
    /// a deletion.
    pub variations: Vec<String>,
}

/// One `<entry>` of a UniProtKB XML file.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct UniProtEntry {
    /// `@dataset`: `Swiss-Prot` or `TrEMBL`.
    pub dataset: String,
    /// First `<accession>` (the primary accession).
    pub accession: String,
    /// Entry `<name>` (mnemonic, e.g. `ALBU_BOVIN`).
    pub name: String,
    /// Recommended (else first submitted) protein full name.
    pub full_name: String,
    /// Primary gene name of the first `<gene>`.
    pub gene: String,
    /// Organism scientific name.
    pub organism: String,
    /// NCBI taxonomy id of the organism.
    pub taxon_id: String,
    /// `<proteinExistence @type>`, e.g. `evidence at protein level`.
    pub protein_existence: String,
    /// `<sequence @version>`.
    pub sequence_version: String,
    /// Canonical sequence, uppercase, whitespace stripped.
    pub sequence: Vec<u8>,
    /// Features in document order.
    pub features: Vec<UniProtFeature>,
    /// `(key, ECO code)` of the entry's `<evidence>` elements.
    pub evidence: Vec<(u32, String)>,
    /// Ids of the isoforms whose sequence is the displayed (canonical) one.
    pub displayed_isoforms: Vec<String>,
}

/// A sequence variant applied to its entry's canonical sequence.
#[derive(Debug, Clone, PartialEq, Eq)]
struct SequenceVariant {
    /// 0-based half-open canonical range the variant replaces.
    start: usize,
    end: usize,
    /// Replacement residues; empty for a deletion.
    replacement: Vec<u8>,
    /// Accession tag: `<feature id>_<change>` or `VAR_<change>`.
    tag: String,
}

impl UniProtEntry {
    /// UniProt FASTA accession token: `sp|P02769|ALBU_BOVIN`.
    pub fn fasta_accession(&self) -> String {
        format!("{}|{}|{}", self.db_code(), self.accession, self.name)
    }

    fn db_code(&self) -> &'static str {
        if self.dataset.eq_ignore_ascii_case("TrEMBL") {
            "tr"
        } else {
            "sp"
        }
    }

    /// UniProt FASTA description: `Albumin OS=Bos taurus OX=9913 GN=ALB PE=1 SV=4`
    /// (fields UniProt has no value for are left out).
    pub fn fasta_description(&self) -> String {
        let mut d = self.full_name.clone();
        let mut field = |key: &str, value: &str| {
            if !value.is_empty() {
                if !d.is_empty() {
                    d.push(' ');
                }
                d.push_str(key);
                d.push('=');
                d.push_str(value);
            }
        };
        field("OS", &self.organism);
        field("OX", &self.taxon_id);
        field("GN", &self.gene);
        field("PE", protein_existence_code(&self.protein_existence));
        field("SV", &self.sequence_version);
        d
    }

    /// ECO codes of `feature`'s evidence keys.
    fn evidence_codes<'a>(&'a self, feature: &'a UniProtFeature) -> impl Iterator<Item = &'a str> {
        feature.evidence_keys.iter().filter_map(|k| {
            self.evidence
                .iter()
                .find(|(key, _)| key == k)
                .map(|(_, code)| code.as_str())
        })
    }

    /// Whether `feature`'s coordinates are the canonical sequence's.
    fn on_canonical(&self, feature: &UniProtFeature) -> bool {
        let s = &feature.location_sequence;
        s.is_empty() || *s == self.accession || self.displayed_isoforms.iter().any(|d| d == s)
    }

    /// The annotated modification sites of the canonical sequence that pass
    /// `opts.evidence`, sorted by position; tallies the rest in `stats`.
    fn mod_sites(&self, opts: &UniProtOptions, stats: &mut UniProtLoadStats) -> Vec<ModSite> {
        let mut sites = Vec::new();
        for f in &self.features {
            let resolved = match f.kind.as_str() {
                "modified residue" | "lipid moiety-binding region" | "glycosylation site" => {
                    let name = clean_description(&f.description);
                    match ptm_mass(name) {
                        Some(m) => m,
                        None => {
                            *stats
                                .unknown_modifications
                                .entry(name.to_string())
                                .or_insert(0) += 1;
                            continue;
                        }
                    }
                }
                "cross-link" => {
                    let name = clean_description(&f.description);
                    match name.strip_prefix(GLYGLY_CROSSLINK_PREFIX) {
                        Some(rest) if GLYGLY_MODIFIERS.contains(&rest) => {
                            (b'K', GLYGLY_MASS, SiteAnchor::Residue)
                        }
                        _ => {
                            *stats
                                .unknown_modifications
                                .entry(name.to_string())
                                .or_insert(0) += 1;
                            continue;
                        }
                    }
                }
                _ => continue,
            };
            if !self.on_canonical(f) {
                stats.isoform_features += 1;
                continue;
            }
            let (target, mass_delta, anchor) = resolved;
            let position = match f.location {
                Some((b, e)) if b == e && b >= 1 && (b as usize) <= self.sequence.len() => b - 1,
                _ => {
                    stats.unusable_features += 1;
                    continue;
                }
            };
            if self.sequence[position as usize] != target {
                stats.unusable_features += 1;
                continue;
            }
            if !opts.evidence.admits(self.evidence_codes(f)) {
                stats.below_evidence += 1;
                continue;
            }
            sites.push(ModSite {
                position,
                mass_delta,
                anchor,
            });
        }
        ProteinAnnotation::new(sites, None).mod_sites
    }

    /// The natural sequence variants that pass `opts.evidence` and fit the
    /// canonical sequence, in document order with duplicates dropped.
    fn sequence_variants(
        &self,
        opts: &UniProtOptions,
        stats: &mut UniProtLoadStats,
    ) -> Vec<SequenceVariant> {
        let mut out: Vec<SequenceVariant> = Vec::new();
        for f in self
            .features
            .iter()
            .filter(|f| f.kind == "sequence variant")
        {
            if !self.on_canonical(f) {
                stats.isoform_features += 1;
                continue;
            }
            let Some((begin, end)) = f
                .location
                .filter(|&(b, e)| b >= 1 && b <= e && (e as usize) <= self.sequence.len())
            else {
                stats.unusable_features += 1;
                continue;
            };
            let (start, end) = (begin as usize - 1, end as usize);
            let canonical = &self.sequence[start..end];
            // A deletion ("Missing") carries neither <original> nor <variation>:
            // the location is the deleted residues.
            let replacements: Vec<&[u8]> = if f.original.is_empty() && f.variations.is_empty() {
                vec![&[]]
            } else if f.original.as_bytes() == canonical {
                f.variations.iter().map(|v| v.as_bytes()).collect()
            } else {
                Vec::new()
            };
            if replacements.is_empty() {
                stats.unusable_features += 1;
                continue;
            }
            if !opts.evidence.admits(self.evidence_codes(f)) {
                stats.below_evidence += 1;
                continue;
            }
            for replacement in replacements {
                // A deletion at either end only moves a protein terminus: its
                // peptides read like canonical ones, so it has nothing to add.
                let terminal_deletion =
                    replacement.is_empty() && (start == 0 || end == self.sequence.len());
                // I and L have the same mass: an I↔L change is the canonical
                // peptide again as far as a mass spectrum can tell.
                let il = |c: &u8| if *c == b'I' { b'L' } else { *c };
                let isobaric = replacement.len() == canonical.len()
                    && replacement.iter().map(il).eq(canonical.iter().map(il));
                if isobaric || terminal_deletion || !replacement.iter().all(u8::is_ascii_uppercase)
                {
                    stats.unusable_features += 1;
                    continue;
                }
                let change = change_label(canonical, begin, end as u32, replacement);
                let tag = if f.id.is_empty() {
                    format!("VAR_{change}")
                } else {
                    format!("{}_{change}", f.id)
                };
                let v = SequenceVariant {
                    start,
                    end,
                    replacement: replacement.to_vec(),
                    tag,
                };
                if !out
                    .iter()
                    .any(|o| o.start == v.start && o.end == v.end && o.replacement == v.replacement)
                {
                    out.push(v);
                }
            }
        }
        out
    }

    /// The proteins this entry contributes: the canonical one, then (with
    /// `opts.variants`) one per sequence variant.
    pub fn into_proteins(
        self,
        opts: &UniProtOptions,
        stats: &mut UniProtLoadStats,
    ) -> Vec<Protein> {
        let accession = self.fasta_accession();
        let description = self.fasta_description();
        let sites = self.mod_sites(opts, stats);
        stats.mod_sites += sites.len();
        let variants = if opts.variants {
            self.sequence_variants(opts, stats)
        } else {
            Vec::new()
        };
        let mut out = Vec::with_capacity(1 + variants.len());
        for v in variants {
            let mut sequence = Vec::with_capacity(self.sequence.len() + v.replacement.len());
            sequence.extend_from_slice(&self.sequence[..v.start]);
            sequence.extend_from_slice(&v.replacement);
            sequence.extend_from_slice(&self.sequence[v.end..]);
            // Sites keep their residue: before the variant unchanged, after it
            // shifted by the length change, inside it dropped.
            let shift = v.replacement.len() as i64 - (v.end - v.start) as i64;
            let mod_sites = sites
                .iter()
                .filter(|s| (s.position as usize) < v.start || (s.position as usize) >= v.end)
                .map(|s| ModSite {
                    position: if (s.position as usize) < v.start {
                        s.position
                    } else {
                        (s.position as i64 + shift) as u32
                    },
                    ..*s
                })
                .collect();
            out.push(Protein {
                accession: format!(
                    "{}|{}_{}|{}",
                    self.db_code(),
                    self.accession,
                    v.tag,
                    self.name
                ),
                description: format!("{description} VARIANT={}", v.tag),
                sequence,
                annotation: Some(Arc::new(ProteinAnnotation::new(
                    mod_sites,
                    Some(VariantRegion {
                        start: v.start as u32,
                        end: (v.start + v.replacement.len()) as u32,
                    }),
                ))),
            });
        }
        stats.variant_proteins += out.len();
        let canonical = Protein {
            accession,
            description,
            sequence: self.sequence,
            annotation: (!sites.is_empty()).then(|| Arc::new(ProteinAnnotation::new(sites, None))),
        };
        out.insert(0, canonical);
        out
    }
}

/// UniProt FASTA `PE=` code of a `<proteinExistence @type>`.
fn protein_existence_code(kind: &str) -> &'static str {
    match kind {
        "evidence at protein level" => "1",
        "evidence at transcript level" => "2",
        "inferred from homology" => "3",
        "predicted" => "4",
        "uncertain" => "5",
        _ => "",
    }
}

/// The change a variant makes, 1-based: `D7H`, `KR120-121EE`, `F508del`.
fn change_label(original: &[u8], begin: u32, end: u32, replacement: &[u8]) -> String {
    let original = String::from_utf8_lossy(original);
    let range = if begin == end {
        begin.to_string()
    } else {
        format!("{begin}-{end}")
    };
    if replacement.is_empty() {
        format!("{original}{range}del")
    } else {
        format!("{original}{range}{}", String::from_utf8_lossy(replacement))
    }
}

/// Load every entry of a UniProt XML document into a target [`ProteinDb`]:
/// canonical proteins in document order, each followed by its variant
/// proteins when `opts.variants` is set.
pub fn load_uniprot_xml<R: BufRead>(
    reader: R,
    opts: &UniProtOptions,
) -> Result<(ProteinDb, UniProtLoadStats), UniProtParseError> {
    let mut stats = UniProtLoadStats::default();
    let mut proteins = Vec::new();
    for entry in UniProtXmlReader::new(reader) {
        let entry = entry?;
        stats.entries += 1;
        proteins.extend(entry.into_proteins(opts, &mut stats));
    }
    if stats.entries == 0 {
        return Err(UniProtParseError::NoEntries);
    }
    Ok((ProteinDb { proteins }, stats))
}

#[derive(thiserror::Error, Debug)]
pub enum UniProtParseError {
    #[error("UniProt XML parse error at byte {position}: {source}")]
    Xml {
        position: usize,
        #[source]
        source: quick_xml::Error,
    },
    #[error("UniProt XML entry {accession:?} has no sequence")]
    MissingSequence { accession: String },
    #[error("UniProt XML ends inside <{element}> (truncated file?)")]
    Truncated { element: String },
    #[error("no UniProt <entry> elements found (not a UniProtKB XML file?)")]
    NoEntries,
}

/// Elements the reader tracks; everything else is [`Tag::Other`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Tag {
    Entry,
    Accession,
    Name,
    Protein,
    RecommendedName,
    SubmittedName,
    FullName,
    Gene,
    Organism,
    DbReference,
    ProteinExistence,
    Comment,
    Isoform,
    Id,
    Feature,
    Original,
    Variation,
    Location,
    Position,
    Begin,
    End,
    Evidence,
    Sequence,
    Other,
}

impl Tag {
    /// The element name, `"…"` for an untracked one.
    fn name(self) -> &'static str {
        match self {
            Self::Entry => "entry",
            Self::Accession => "accession",
            Self::Name => "name",
            Self::Protein => "protein",
            Self::RecommendedName => "recommendedName",
            Self::SubmittedName => "submittedName",
            Self::FullName => "fullName",
            Self::Gene => "gene",
            Self::Organism => "organism",
            Self::DbReference => "dbReference",
            Self::ProteinExistence => "proteinExistence",
            Self::Comment => "comment",
            Self::Isoform => "isoform",
            Self::Id => "id",
            Self::Feature => "feature",
            Self::Original => "original",
            Self::Variation => "variation",
            Self::Location => "location",
            Self::Position => "position",
            Self::Begin => "begin",
            Self::End => "end",
            Self::Evidence => "evidence",
            Self::Sequence => "sequence",
            Self::Other => "…",
        }
    }

    fn of(local_name: &[u8]) -> Self {
        match local_name {
            b"entry" => Self::Entry,
            b"accession" => Self::Accession,
            b"name" => Self::Name,
            b"protein" => Self::Protein,
            b"recommendedName" => Self::RecommendedName,
            b"submittedName" => Self::SubmittedName,
            b"fullName" => Self::FullName,
            b"gene" => Self::Gene,
            b"organism" => Self::Organism,
            b"dbReference" => Self::DbReference,
            b"proteinExistence" => Self::ProteinExistence,
            b"comment" => Self::Comment,
            b"isoform" => Self::Isoform,
            b"id" => Self::Id,
            b"feature" => Self::Feature,
            b"original" => Self::Original,
            b"variation" => Self::Variation,
            b"location" => Self::Location,
            b"position" => Self::Position,
            b"begin" => Self::Begin,
            b"end" => Self::End,
            b"evidence" => Self::Evidence,
            b"sequence" => Self::Sequence,
            _ => Self::Other,
        }
    }
}

/// Where the text of the element being read goes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Capture {
    Accession,
    EntryName,
    RecommendedFullName,
    SubmittedFullName,
    Gene,
    Organism,
    IsoformId,
    Original,
    Variation,
    Sequence,
}

/// Iterator over the `<entry>` elements of a UniProtKB XML document, one
/// [`UniProtEntry`] at a time (constant memory in the number of entries).
pub struct UniProtXmlReader<R: BufRead> {
    xml: Reader<R>,
    buf: Vec<u8>,
    /// Open elements, outermost first.
    stack: Vec<Tag>,
    entry: UniProtEntry,
    feature: UniProtFeature,
    /// Feature `<begin>`/`<end>` positions while the feature is open.
    begin: Option<u32>,
    end: Option<u32>,
    /// Inside `<comment type="alternative products">`.
    in_alt_products: bool,
    /// Current isoform: first `<id>` and whether its sequence is displayed.
    isoform_id: Option<String>,
    isoform_displayed: bool,
    genes_seen: usize,
    submitted_name: String,
    capture: Option<(Capture, usize)>,
    text: String,
    done: bool,
}

impl<R: BufRead> UniProtXmlReader<R> {
    pub fn new(reader: R) -> Self {
        let mut xml = Reader::from_reader(reader);
        xml.trim_text(true);
        Self {
            xml,
            buf: Vec::new(),
            stack: Vec::new(),
            entry: UniProtEntry::default(),
            feature: UniProtFeature::default(),
            begin: None,
            end: None,
            in_alt_products: false,
            isoform_id: None,
            isoform_displayed: false,
            genes_seen: 0,
            submitted_name: String::new(),
            capture: None,
            text: String::new(),
            done: false,
        }
    }

    fn err(&self, source: quick_xml::Error) -> UniProtParseError {
        UniProtParseError::Xml {
            position: self.xml.buffer_position(),
            source,
        }
    }

    /// The innermost open element the reader tracks (`entry` inside an entry
    /// whose children it does not track).
    fn open_element(&self) -> &'static str {
        self.stack
            .iter()
            .rev()
            .find(|t| **t != Tag::Other)
            .map_or("…", |t| t.name())
    }

    /// `stack[len - 1 - up]`, i.e. `up = 0` is the element itself.
    fn ancestor(&self, up: usize) -> Option<Tag> {
        self.stack.len().checked_sub(1 + up).map(|i| self.stack[i])
    }

    /// Handle an element start (an empty element is a start plus an end).
    fn on_start(&mut self, e: &BytesStart<'_>) -> Result<(), quick_xml::Error> {
        let tag = Tag::of(e.local_name().as_ref());
        self.stack.push(tag);
        self.capture = None;
        self.text.clear();
        let parent = self.ancestor(1);
        let in_entry = self.stack.contains(&Tag::Entry);
        if !in_entry {
            return Ok(());
        }
        let depth = self.stack.len();
        let capture = |c| Some((c, depth));
        match (tag, parent) {
            (Tag::Entry, _) => {
                self.entry = UniProtEntry {
                    dataset: attr(e, b"dataset")?.unwrap_or_default(),
                    ..UniProtEntry::default()
                };
                self.genes_seen = 0;
                self.submitted_name.clear();
                self.in_alt_products = false;
            }
            (Tag::Accession, Some(Tag::Entry)) => self.capture = capture(Capture::Accession),
            (Tag::Name, Some(Tag::Entry)) => self.capture = capture(Capture::EntryName),
            (Tag::FullName, Some(Tag::RecommendedName | Tag::SubmittedName))
                if self.ancestor(2) == Some(Tag::Protein)
                    && self.ancestor(3) == Some(Tag::Entry) =>
            {
                self.capture = if parent == Some(Tag::RecommendedName) {
                    capture(Capture::RecommendedFullName)
                } else {
                    capture(Capture::SubmittedFullName)
                };
            }
            (Tag::Gene, Some(Tag::Entry)) => self.genes_seen += 1,
            (Tag::Name, Some(Tag::Gene))
                if self.genes_seen == 1
                    && self.entry.gene.is_empty()
                    && attr(e, b"type")?.as_deref() == Some("primary") =>
            {
                self.capture = capture(Capture::Gene);
            }
            (Tag::Name, Some(Tag::Organism))
                if self.entry.organism.is_empty()
                    && attr(e, b"type")?.as_deref() == Some("scientific") =>
            {
                self.capture = capture(Capture::Organism);
            }
            (Tag::DbReference, Some(Tag::Organism))
                if self.entry.taxon_id.is_empty()
                    && attr(e, b"type")?.as_deref() == Some("NCBI Taxonomy") =>
            {
                self.entry.taxon_id = attr(e, b"id")?.unwrap_or_default();
            }
            (Tag::ProteinExistence, Some(Tag::Entry)) => {
                self.entry.protein_existence = attr(e, b"type")?.unwrap_or_default();
            }
            (Tag::Comment, Some(Tag::Entry)) => {
                self.in_alt_products = attr(e, b"type")?.as_deref() == Some("alternative products");
            }
            (Tag::Isoform, Some(Tag::Comment)) if self.in_alt_products => {
                self.isoform_id = None;
                self.isoform_displayed = false;
            }
            (Tag::Id, Some(Tag::Isoform)) if self.isoform_id.is_none() => {
                self.capture = capture(Capture::IsoformId);
            }
            (Tag::Sequence, Some(Tag::Isoform)) => {
                self.isoform_displayed = attr(e, b"type")?.as_deref() == Some("displayed");
            }
            (Tag::Feature, Some(Tag::Entry)) => {
                self.feature = UniProtFeature {
                    kind: attr(e, b"type")?.unwrap_or_default(),
                    id: attr(e, b"id")?.unwrap_or_default(),
                    description: attr(e, b"description")?.unwrap_or_default(),
                    evidence_keys: attr(e, b"evidence")?
                        .map(|s| {
                            s.split_whitespace()
                                .filter_map(|k| k.parse().ok())
                                .collect()
                        })
                        .unwrap_or_default(),
                    ..UniProtFeature::default()
                };
                self.begin = None;
                self.end = None;
            }
            (Tag::Original, Some(Tag::Feature)) => self.capture = capture(Capture::Original),
            (Tag::Variation, Some(Tag::Feature)) => self.capture = capture(Capture::Variation),
            (Tag::Location, Some(Tag::Feature)) => {
                self.feature.location_sequence = attr(e, b"sequence")?.unwrap_or_default();
            }
            (Tag::Position | Tag::Begin | Tag::End, Some(Tag::Location))
                if self.ancestor(2) == Some(Tag::Feature) =>
            {
                // `status="unknown"` coordinates carry no `position`.
                let p = attr(e, b"position")?.and_then(|s| s.trim().parse::<u32>().ok());
                match tag {
                    Tag::Position => {
                        self.begin = p;
                        self.end = p;
                    }
                    Tag::Begin => self.begin = p,
                    _ => self.end = p,
                }
            }
            (Tag::Evidence, Some(Tag::Entry)) => {
                if let (Some(code), Some(key)) = (attr(e, b"type")?, attr(e, b"key")?) {
                    if let Ok(key) = key.trim().parse() {
                        self.entry.evidence.push((key, code));
                    }
                }
            }
            (Tag::Sequence, Some(Tag::Entry)) => {
                self.entry.sequence_version = attr(e, b"version")?.unwrap_or_default();
                self.capture = capture(Capture::Sequence);
            }
            _ => {}
        }
        Ok(())
    }

    /// Handle an element end; returns the finished entry at `</entry>`.
    fn on_end(&mut self) -> Result<Option<UniProtEntry>, UniProtParseError> {
        let depth = self.stack.len();
        if let Some((target, at)) = self.capture.take() {
            if at == depth {
                let text = std::mem::take(&mut self.text);
                let e = &mut self.entry;
                match target {
                    Capture::Accession if e.accession.is_empty() => e.accession = text,
                    Capture::EntryName if e.name.is_empty() => e.name = text,
                    Capture::RecommendedFullName if e.full_name.is_empty() => e.full_name = text,
                    Capture::SubmittedFullName if self.submitted_name.is_empty() => {
                        self.submitted_name = text
                    }
                    Capture::Gene => e.gene = text,
                    Capture::Organism => e.organism = text,
                    Capture::IsoformId => self.isoform_id = Some(text),
                    Capture::Original => self.feature.original = text,
                    Capture::Variation => self.feature.variations.push(text),
                    Capture::Sequence => {
                        e.sequence = text
                            .bytes()
                            .filter(|b| !b.is_ascii_whitespace())
                            .map(|b| b.to_ascii_uppercase())
                            .collect();
                    }
                    _ => {}
                }
            }
        }
        let tag = self.stack.pop();
        let parent = self.ancestor(0);
        match (tag, parent) {
            (Some(Tag::Isoform), Some(Tag::Comment)) if self.in_alt_products => {
                if self.isoform_displayed {
                    if let Some(id) = self.isoform_id.take() {
                        self.entry.displayed_isoforms.push(id);
                    }
                }
            }
            (Some(Tag::Comment), Some(Tag::Entry)) => self.in_alt_products = false,
            (Some(Tag::Feature), Some(Tag::Entry)) => {
                let mut f = std::mem::take(&mut self.feature);
                f.location = self.begin.zip(self.end);
                self.entry.features.push(f);
            }
            (Some(Tag::Entry), _) => {
                let mut entry = std::mem::take(&mut self.entry);
                if entry.full_name.is_empty() {
                    entry.full_name = std::mem::take(&mut self.submitted_name);
                }
                if entry.sequence.is_empty() {
                    return Err(UniProtParseError::MissingSequence {
                        accession: entry.accession,
                    });
                }
                return Ok(Some(entry));
            }
            _ => {}
        }
        Ok(None)
    }
}

impl<R: BufRead> Iterator for UniProtXmlReader<R> {
    type Item = Result<UniProtEntry, UniProtParseError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        // The event borrows the read buffer, not `self`, so the handlers can
        // take `&mut self` without copying every event.
        let mut buf = std::mem::take(&mut self.buf);
        let result = loop {
            buf.clear();
            let event = match self.xml.read_event_into(&mut buf) {
                Ok(ev) => ev,
                Err(e) => break Some(Err(self.err(e))),
            };
            let finished = match event {
                Event::Start(ref e) => self.on_start(e).map(|()| None).map_err(|x| self.err(x)),
                Event::Empty(ref e) => match self.on_start(e) {
                    Ok(()) => self.on_end(),
                    Err(x) => Err(self.err(x)),
                },
                Event::End(_) => self.on_end(),
                Event::Text(ref t) if self.capture.is_some() => match t.unescape() {
                    Ok(s) => {
                        self.text.push_str(&s);
                        Ok(None)
                    }
                    Err(x) => Err(self.err(x)),
                },
                Event::CData(ref t) if self.capture.is_some() => {
                    self.text.push_str(&String::from_utf8_lossy(t));
                    Ok(None)
                }
                Event::Eof if self.stack.is_empty() => break None,
                Event::Eof => {
                    break Some(Err(UniProtParseError::Truncated {
                        element: self.open_element().to_string(),
                    }))
                }
                _ => Ok(None),
            };
            match finished {
                Ok(Some(entry)) => break Some(Ok(entry)),
                Ok(None) => {}
                Err(e) => break Some(Err(e)),
            }
        };
        self.buf = buf;
        if !matches!(result, Some(Ok(_))) {
            self.done = true;
        }
        result
    }
}

/// Unescaped value of attribute `name` (matched on its local name).
fn attr(e: &BytesStart<'_>, name: &[u8]) -> Result<Option<String>, quick_xml::Error> {
    for a in e.attributes() {
        let a = a?;
        if a.key.local_name().as_ref() == name {
            return Ok(Some(a.unescape_value()?.into_owned()));
        }
    }
    Ok(None)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Two entries exercising every branch the reader takes. Sequence of
    /// P00001 (1-based): M1 S2 E3 K4 R5 S6 P7 T8 K9 A10 G11 N12 L13 S14 T15
    /// Y16 K17 D18 E19 R20.
    const XML: &str = r#"<?xml version="1.0" encoding="UTF-8"?>
<uniprot xmlns="http://uniprot.org/uniprot">
<entry dataset="Swiss-Prot" version="3">
  <accession>P00001</accession>
  <accession>Q99999</accession>
  <name>TEST1_HUMAN</name>
  <protein>
    <recommendedName><fullName evidence="1">Test protein &amp; friend</fullName></recommendedName>
    <alternativeName><fullName>Other name</fullName></alternativeName>
    <component><recommendedName><fullName>Chain A</fullName></recommendedName></component>
  </protein>
  <gene><name type="synonym">SYN1</name><name type="primary">TST1</name></gene>
  <gene><name type="primary">TST2</name></gene>
  <organism>
    <name type="scientific">Homo sapiens</name>
    <name type="common">Human</name>
    <dbReference type="NCBI Taxonomy" id="9606"/>
  </organism>
  <organismHost><name type="scientific">Mus musculus</name></organismHost>
  <comment type="alternative products">
    <isoform><id>P00001-1</id><name>1</name><sequence type="displayed"/></isoform>
    <isoform><id>P00001-2</id><name>2</name><sequence type="described" ref="VSP_1"/></isoform>
  </comment>
  <comment type="function"><text>Some function.</text></comment>
  <dbReference type="PDB" id="1ABC"/>
  <proteinExistence type="evidence at protein level"/>
  <feature type="modified residue" description="N-acetylserine" evidence="1">
    <location><position position="2"/></location>
  </feature>
  <feature type="modified residue" description="Phosphoserine; by PKA" evidence="1 2">
    <location><position position="6"/></location>
  </feature>
  <feature type="modified residue" description="Phosphothreonine" evidence="2">
    <location><position position="8"/></location>
  </feature>
  <feature type="modified residue" description="Phosphoserine" evidence="1">
    <location sequence="P00001-2"><position position="14"/></location>
  </feature>
  <feature type="modified residue" description="Phosphoserine" evidence="1">
    <location sequence="P00001-1"><position position="6"/></location>
  </feature>
  <feature type="modified residue" description="N6-lactoyllysine" evidence="1">
    <location><position position="9"/></location>
  </feature>
  <feature type="modified residue" description="Phosphoserine" evidence="1">
    <location><position position="10"/></location>
  </feature>
  <feature type="cross-link" description="Glycyl lysine isopeptide (Lys-Gly) (interchain with G-Cter in ubiquitin); alternate" evidence="1">
    <location><position position="17"/></location>
  </feature>
  <feature type="cross-link" description="Glycyl lysine isopeptide (Lys-Gly) (interchain with G-Cter in SUMO2)" evidence="1">
    <location><position position="9"/></location>
  </feature>
  <feature type="glycosylation site" description="N-linked (GlcNAc...) asparagine" evidence="1">
    <location><position position="12"/></location>
  </feature>
  <feature type="modified residue" description="Phosphotyrosine" evidence="3">
    <location><position position="16"/></location>
  </feature>
  <feature type="sequence variant" id="VAR_000001" description="In dbSNP:rs1." evidence="1">
    <original>K</original><variation>E</variation>
    <location><position position="4"/></location>
  </feature>
  <feature type="sequence variant" id="VAR_000002" evidence="1">
    <original>S</original><variation>L</variation><variation>P</variation>
    <location><position position="14"/></location>
  </feature>
  <feature type="sequence variant" evidence="1">
    <original>AG</original><variation>VVV</variation>
    <location><begin position="10"/><end position="11"/></location>
  </feature>
  <feature type="sequence variant" id="VAR_000003" description="Missing." evidence="1">
    <location><position position="13"/></location>
  </feature>
  <feature type="sequence variant" id="VAR_000004" evidence="1">
    <original>Q</original><variation>H</variation>
    <location><position position="3"/></location>
  </feature>
  <feature type="sequence variant" id="VAR_000005" description="In dbSNP:rs5.">
    <original>E</original><variation>K</variation>
    <location><position position="19"/></location>
  </feature>
  <feature type="sequence variant" id="VAR_000006" evidence="1">
    <original>R</original><variation>W</variation>
    <location><begin status="unknown"/><end position="5"/></location>
  </feature>
  <feature type="sequence variant" id="VAR_000007" evidence="1">
    <location><position position="1"/></location>
  </feature>
  <feature type="sequence variant" id="VAR_000008" evidence="1">
    <original>S</original><variation>A</variation>
    <location sequence="P00001-2"><position position="2"/></location>
  </feature>
  <feature type="sequence variant" id="VAR_000009" evidence="1">
    <original>L</original><variation>I</variation>
    <location><position position="13"/></location>
  </feature>
  <feature type="mutagenesis site" description="Loss of function." evidence="1">
    <original>S</original><variation>A</variation>
    <location><position position="6"/></location>
  </feature>
  <evidence type="ECO:0000269" key="1"><source><dbReference type="PubMed" id="1"/></source></evidence>
  <evidence type="ECO:0000250" key="2"/>
  <evidence type="ECO:0000256" key="3"/>
  <sequence length="20" mass="2300" checksum="0" modified="2000-01-01" version="2">MSEKRSPTKAGNLSTYKDER</sequence>
</entry>
<entry dataset="TrEMBL" version="1">
  <accession>A0A000</accession>
  <name>A0A000_YEAST</name>
  <protein><submittedName><fullName>Uncharacterized protein</fullName></submittedName></protein>
  <organism>
    <name type="scientific">Saccharomyces cerevisiae</name>
    <dbReference type="NCBI Taxonomy" id="4932"/>
  </organism>
  <proteinExistence type="predicted"/>
  <sequence length="8" version="1">
    MKVLL
    AGR
  </sequence>
</entry>
</uniprot>
"#;

    const CANONICAL: &[u8] = b"MSEKRSPTKAGNLSTYKDER";

    fn load(opts: UniProtOptions) -> (ProteinDb, UniProtLoadStats) {
        load_uniprot_xml(XML.as_bytes(), &opts).unwrap()
    }

    fn sites(p: &Protein) -> Vec<(u32, f64)> {
        p.mod_sites()
            .iter()
            .map(|s| (s.position, s.mass_delta))
            .collect()
    }

    #[test]
    fn entries_parse_with_uniprot_fasta_headers() {
        let entries: Vec<UniProtEntry> = UniProtXmlReader::new(XML.as_bytes())
            .collect::<Result<_, _>>()
            .unwrap();
        assert_eq!(entries.len(), 2);
        let e = &entries[0];
        assert_eq!(e.accession, "P00001");
        assert_eq!(e.name, "TEST1_HUMAN");
        assert_eq!(e.full_name, "Test protein & friend");
        assert_eq!(e.gene, "TST1", "primary name of the first gene only");
        assert_eq!(e.organism, "Homo sapiens", "not the organism host");
        assert_eq!(e.taxon_id, "9606");
        assert_eq!(e.sequence, CANONICAL);
        assert_eq!(e.displayed_isoforms, vec!["P00001-1".to_string()]);
        assert_eq!(e.features.len(), 22);
        assert_eq!(e.evidence.len(), 3);
        assert_eq!(e.fasta_accession(), "sp|P00001|TEST1_HUMAN");
        assert_eq!(
            e.fasta_description(),
            "Test protein & friend OS=Homo sapiens OX=9606 GN=TST1 PE=1 SV=2"
        );
        let range = &e.features[17];
        assert_eq!(range.id, "VAR_000006");
        assert_eq!(range.location, None, "an unknown begin leaves no location");
        let multi = &e.features[12];
        assert_eq!(multi.variations, vec!["L".to_string(), "P".to_string()]);

        let t = &entries[1];
        assert_eq!(t.fasta_accession(), "tr|A0A000|A0A000_YEAST");
        assert_eq!(
            t.fasta_description(),
            "Uncharacterized protein OS=Saccharomyces cerevisiae OX=4932 PE=4 SV=1"
        );
        assert_eq!(t.sequence, b"MKVLLAGR");
    }

    #[test]
    fn sequences_only_by_default() {
        let (db, stats) = load(UniProtOptions::default());
        let accs: Vec<&str> = db.proteins.iter().map(|p| p.accession.as_str()).collect();
        assert_eq!(accs, ["sp|P00001|TEST1_HUMAN", "tr|A0A000|A0A000_YEAST"]);
        assert_eq!(db.proteins[0].sequence, CANONICAL);
        assert_eq!(stats.entries, 2);
        assert_eq!(stats.variant_proteins, 0);
        // An entry without usable sites carries no annotation at all.
        assert!(db.proteins[1].annotation.is_none());
    }

    #[test]
    fn modification_sites_follow_the_evidence_level() {
        const ACETYL: f64 = 42.010565;
        const PHOSPHO: f64 = 79.966331;
        let (db, stats) = load(UniProtOptions::default());
        assert_eq!(
            sites(&db.proteins[0]),
            vec![(1, ACETYL), (5, PHOSPHO), (16, GLYGLY_MASS)],
            "experimental sites; the displayed-isoform duplicate of S6 is merged"
        );
        assert_eq!(stats.mod_sites, 3);
        assert_eq!(
            stats.below_evidence, 2,
            "T8 (by similarity), Y16 (automatic)"
        );
        assert_eq!(stats.unusable_features, 1, "phospho on A10");
        assert_eq!(stats.isoform_features, 1, "S14 of isoform 2");
        assert_eq!(
            stats.unknown_modifications.keys().collect::<Vec<_>>(),
            [
                "Glycyl lysine isopeptide (Lys-Gly) (interchain with G-Cter in SUMO2)",
                "N-linked (GlcNAc...) asparagine",
                "N6-lactoyllysine"
            ]
        );
        // N-acetylserine is tied to the N-terminus, the rest to their residue.
        let anchors: Vec<SiteAnchor> = db.proteins[0]
            .mod_sites()
            .iter()
            .map(|s| s.anchor)
            .collect();
        assert_eq!(
            anchors,
            [
                SiteAnchor::NTerminus,
                SiteAnchor::Residue,
                SiteAnchor::Residue
            ]
        );

        let curated = UniProtOptions {
            evidence: EvidenceLevel::Curated,
            ..UniProtOptions::default()
        };
        let (db, _) = load(curated);
        assert_eq!(
            sites(&db.proteins[0]),
            vec![(1, ACETYL), (5, PHOSPHO), (7, PHOSPHO), (16, GLYGLY_MASS)]
        );
        let any = UniProtOptions {
            evidence: EvidenceLevel::Any,
            ..UniProtOptions::default()
        };
        let (db, _) = load(any);
        assert_eq!(
            sites(&db.proteins[0]),
            vec![
                (1, ACETYL),
                (5, PHOSPHO),
                (7, PHOSPHO),
                (15, PHOSPHO),
                (16, GLYGLY_MASS)
            ]
        );
    }

    #[test]
    fn variants_become_proteins_carrying_only_the_change() {
        let opts = UniProtOptions {
            variants: true,
            ..UniProtOptions::default()
        };
        let (db, stats) = load(opts);
        let accs: Vec<&str> = db.proteins.iter().map(|p| p.accession.as_str()).collect();
        assert_eq!(
            accs,
            [
                "sp|P00001|TEST1_HUMAN",
                "sp|P00001_VAR_000001_K4E|TEST1_HUMAN",
                "sp|P00001_VAR_000002_S14L|TEST1_HUMAN",
                "sp|P00001_VAR_000002_S14P|TEST1_HUMAN",
                "sp|P00001_VAR_AG10-11VVV|TEST1_HUMAN",
                "sp|P00001_VAR_000003_L13del|TEST1_HUMAN",
                "tr|A0A000|A0A000_YEAST",
            ]
        );
        assert_eq!(stats.variant_proteins, 5);
        // + E19K without evidence; T8 and Y16 as before.
        assert_eq!(stats.below_evidence, 3);
        // + Q3H (original mismatch), unknown range, M1 deletion at the terminus,
        // and L13I, which no mass spectrum tells from the canonical L13.
        assert_eq!(stats.unusable_features, 5);
        assert_eq!(stats.isoform_features, 2);

        let k4e = &db.proteins[1];
        assert_eq!(k4e.sequence, b"MSEERSPTKAGNLSTYKDER");
        assert_eq!(
            k4e.variant_region(),
            Some(VariantRegion { start: 3, end: 4 })
        );
        assert!(k4e.description.ends_with("SV=2 VARIANT=VAR_000001_K4E"));
        assert_eq!(sites(k4e), sites(&db.proteins[0]), "same-length change");

        let ins = &db.proteins[4];
        assert_eq!(ins.sequence, b"MSEKRSPTKVVVNLSTYKDER");
        assert_eq!(
            ins.variant_region(),
            Some(VariantRegion { start: 9, end: 12 })
        );
        // The GlyGly K17 moves one residue right and still sits on its K.
        assert_eq!(
            sites(ins).iter().map(|s| s.0).collect::<Vec<_>>(),
            [1, 5, 17]
        );
        assert_eq!(ins.sequence[17], b'K');

        let del = &db.proteins[5];
        assert_eq!(del.sequence, b"MSEKRSPTKAGNSTYKDER");
        assert_eq!(
            del.variant_region(),
            Some(VariantRegion { start: 12, end: 12 })
        );
        assert_eq!(
            sites(del).iter().map(|s| s.0).collect::<Vec<_>>(),
            [1, 5, 15]
        );
        // The canonical protein is not a variant.
        assert_eq!(db.proteins[0].variant_region(), None);
    }

    #[test]
    fn evidence_levels_nest() {
        let codes = |c: &'static [&'static str]| c.iter().copied();
        const EXP: &[&str] = &["ECO:0000269"];
        const SIM: &[&str] = &["ECO:0000250"];
        const AUTO: &[&str] = &["ECO:0000256"];
        for (level, admitted) in [
            (EvidenceLevel::Experimental, [true, false, false, false]),
            (EvidenceLevel::Curated, [true, true, false, false]),
            (EvidenceLevel::Any, [true, true, true, true]),
        ] {
            assert_eq!(level.admits(codes(EXP)), admitted[0], "{level:?}");
            assert_eq!(level.admits(codes(SIM)), admitted[1], "{level:?}");
            assert_eq!(level.admits(codes(AUTO)), admitted[2], "{level:?}");
            assert_eq!(level.admits(codes(&[])), admitted[3], "{level:?}");
        }
        assert_eq!(
            EvidenceLevel::from_name(" Curated "),
            Some(EvidenceLevel::Curated)
        );
        assert_eq!(EvidenceLevel::from_name("strong"), None);
        for l in [
            EvidenceLevel::Experimental,
            EvidenceLevel::Curated,
            EvidenceLevel::Any,
        ] {
            assert_eq!(EvidenceLevel::from_name(l.as_str()), Some(l));
        }
    }

    #[test]
    fn malformed_documents_are_errors() {
        // Cut inside the second entry: the first is read, then the error.
        let cut = XML.rfind("<sequence length=\"8\"").unwrap();
        let r: Vec<_> = UniProtXmlReader::new(&XML.as_bytes()[..cut]).collect();
        assert_eq!(r.len(), 2);
        assert!(r[0].is_ok());
        assert!(matches!(
            &r[1],
            Err(UniProtParseError::Truncated { element }) if element == "entry"
        ));
        let opts = UniProtOptions::default();
        assert!(load_uniprot_xml(&XML.as_bytes()[..cut], &opts).is_err());
        // Well-formed XML that is not UniProt holds no entries.
        let html = "<html><body><p>not UniProt</p></body></html>";
        assert!(matches!(
            load_uniprot_xml(html.as_bytes(), &opts),
            Err(UniProtParseError::NoEntries)
        ));
        let mismatched = "<uniprot><entry><accession>X</name></entry></uniprot>";
        let r: Result<Vec<_>, _> = UniProtXmlReader::new(mismatched.as_bytes()).collect();
        assert!(matches!(r, Err(UniProtParseError::Xml { .. })));
        let no_seq =
            r#"<uniprot><entry dataset="Swiss-Prot"><accession>X1</accession></entry></uniprot>"#;
        let r: Result<Vec<_>, _> = UniProtXmlReader::new(no_seq.as_bytes()).collect();
        assert!(matches!(
            r,
            Err(UniProtParseError::MissingSequence { accession }) if accession == "X1"
        ));
    }

    #[test]
    fn change_labels() {
        assert_eq!(change_label(b"D", 7, 7, b"H"), "D7H");
        assert_eq!(change_label(b"KR", 120, 121, b"EE"), "KR120-121EE");
        assert_eq!(change_label(b"F", 508, 508, b""), "F508del");
    }
}
