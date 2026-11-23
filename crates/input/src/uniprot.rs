use std::collections::BTreeMap;
use std::io::BufRead;
use std::sync::Arc;

use model::{
    Annotation, ChainEnds, Edit, EntryKind, EntryOrigin, ModSite, PriorTier, Processing, Protein,
    ProteinAnnotation, ProteinDb, SiteAnchor, SiteClaim,
};

use crate::annotation::{build_entry, change_label, knockoff_annotations, KnockoffOptions};
use quick_xml::events::{BytesStart, Event};
use quick_xml::Reader;

use crate::uniprot_ptm::{
    clean_description, ptm_mass, GLYGLY_CROSSLINK_PREFIX, GLYGLY_MASS, GLYGLY_MODIFIERS,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Default)]
pub enum EvidenceLevel {
    #[default]
    Experimental,
    Curated,
    Any,
}

impl EvidenceLevel {
    pub fn from_name(s: &str) -> Option<Self> {
        match s.trim().to_ascii_lowercase().as_str() {
            "experimental" => Some(Self::Experimental),
            "curated" => Some(Self::Curated),
            "any" => Some(Self::Any),
            _ => None,
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Self::Experimental => "experimental",
            Self::Curated => "curated",
            Self::Any => "any",
        }
    }

    fn of_code(code: &str) -> Option<Self> {
        match code {
            "ECO:0000269" | "ECO:0007744" | "ECO:0007829" => Some(Self::Experimental),
            "ECO:0000244" | "ECO:0000250" | "ECO:0000255" | "ECO:0000303" | "ECO:0000305"
            | "ECO:0000312" => Some(Self::Curated),
            _ => None,
        }
    }

    pub fn admits<'a>(self, mut codes: impl Iterator<Item = &'a str>) -> bool {
        self == Self::Any || codes.any(|c| Self::of_code(c).is_some_and(|l| l <= self))
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct UniProtOptions {
    pub variants: bool,
    pub isoforms: bool,
    pub chains: bool,
    pub evidence: EvidenceLevel,
    pub knockoffs: Option<KnockoffOptions>,
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct UniProtLoadStats {
    pub entries: usize,
    pub variant_proteins: usize,
    pub mod_sites: usize,
    pub below_evidence: usize,
    pub unknown_modifications: BTreeMap<String, usize>,
    pub unusable_features: usize,
    pub isoform_features: usize,
    pub isoform_proteins: usize,
    pub unbuilt_isoforms: usize,
    pub chain_proteins: usize,
    pub knockoff_proteins: usize,
    pub variants_without_knockoff: usize,
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct UniProtFeature {
    pub kind: String,
    pub id: String,
    pub description: String,
    pub evidence_keys: Vec<u32>,
    pub location_sequence: String,
    pub location: Option<(u32, u32)>,
    pub original: String,
    pub variations: Vec<String>,
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct UniProtIsoform {
    pub id: String,
    pub name: String,
    pub sequence_type: String,
    pub splice_refs: Vec<String>,
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct UniProtEntry {
    pub dataset: String,
    pub accession: String,
    pub name: String,
    pub full_name: String,
    pub gene: String,
    pub organism: String,
    pub taxon_id: String,
    pub protein_existence: String,
    pub sequence_version: String,
    pub sequence: Vec<u8>,
    pub features: Vec<UniProtFeature>,
    pub evidence: Vec<(u32, String)>,
    pub isoforms: Vec<UniProtIsoform>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct MatureChain {
    start: usize,
    end: usize,
    ends: ChainEnds,
    tag: String,
}

const UNIPROT: &str = "UniProt";

const PROCESSING_FEATURES: [&str; 4] = [
    "signal peptide",
    "transit peptide",
    "propeptide",
    "initiator methionine",
];

impl UniProtEntry {
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

    pub fn fasta_description(&self) -> String {
        self.description_with(&self.full_name, true)
    }

    fn isoform_description(&self, iso: &UniProtIsoform) -> String {
        self.description_with(
            &format!("Isoform {} of {}", iso.name, self.full_name),
            false,
        )
    }

    fn description_with(&self, title: &str, existence_and_version: bool) -> String {
        let mut d = title.to_string();
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
        if existence_and_version {
            field("PE", protein_existence_code(&self.protein_existence));
            field("SV", &self.sequence_version);
        }
        d
    }

    fn evidence_codes<'a>(&'a self, feature: &'a UniProtFeature) -> impl Iterator<Item = &'a str> {
        feature.evidence_keys.iter().filter_map(|k| {
            self.evidence
                .iter()
                .find(|(key, _)| key == k)
                .map(|(_, code)| code.as_str())
        })
    }

    pub fn displayed_isoforms(&self) -> impl Iterator<Item = &str> {
        self.isoforms
            .iter()
            .filter(|i| i.sequence_type == "displayed")
            .map(|i| i.id.as_str())
    }

    fn on_canonical(&self, feature: &UniProtFeature) -> bool {
        let s = feature.location_sequence.as_str();
        s.is_empty() || s == self.accession || self.displayed_isoforms().any(|d| d == s)
    }

    fn located_on(&self, feature: &UniProtFeature, isoform: Option<&str>) -> bool {
        match isoform {
            None => self.on_canonical(feature),
            Some(id) => feature.location_sequence == id,
        }
    }

    fn mod_sites(
        &self,
        isoform: Option<&str>,
        seq: &[u8],
        opts: &UniProtOptions,
        stats: &mut UniProtLoadStats,
    ) -> Vec<ModSite> {
        let mut sites = Vec::new();
        for f in &self.features {
            let site_kind = matches!(
                f.kind.as_str(),
                "modified residue"
                    | "lipid moiety-binding region"
                    | "glycosylation site"
                    | "cross-link"
            );
            if !site_kind || !self.located_on(f, isoform) {
                continue;
            }
            let name = clean_description(&f.description);
            let resolved = if f.kind == "cross-link" {
                match name.strip_prefix(GLYGLY_CROSSLINK_PREFIX) {
                    Some(rest) if GLYGLY_MODIFIERS.contains(&rest) => {
                        Some((b'K', GLYGLY_MASS, SiteAnchor::Residue))
                    }
                    _ => None,
                }
            } else {
                ptm_mass(name)
            };
            let Some((target, mass_delta, anchor)) = resolved else {
                *stats
                    .unknown_modifications
                    .entry(name.to_string())
                    .or_insert(0) += 1;
                continue;
            };
            let position = match f.location {
                Some((b, e)) if b == e && b >= 1 && (b as usize) <= seq.len() => b - 1,
                _ => {
                    stats.unusable_features += 1;
                    continue;
                }
            };
            if seq[position as usize] != target {
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

    fn sequence_variants(
        &self,
        isoform: Option<&str>,
        seq: &[u8],
        opts: &UniProtOptions,
        stats: &mut UniProtLoadStats,
    ) -> Vec<Annotation> {
        let mut out: Vec<Annotation> = Vec::new();
        for f in &self.features {
            if f.kind != "sequence variant" || !self.located_on(f, isoform) {
                continue;
            }
            let Some((begin, end)) = f
                .location
                .filter(|&(b, e)| b >= 1 && b <= e && (e as usize) <= seq.len())
            else {
                stats.unusable_features += 1;
                continue;
            };
            let (start, end) = (begin as usize - 1, end as usize);
            let canonical = &seq[start..end];
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
                let terminal_deletion = replacement.is_empty() && (start == 0 || end == seq.len());
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
                let mut v = Annotation::new(&tag, EntryKind::Variant);
                v.source = Some(Arc::from(UNIPROT));
                v.tier = if EvidenceLevel::Curated.admits(self.evidence_codes(f)) {
                    PriorTier::Curated
                } else {
                    PriorTier::Catalogue
                };
                v.edits = vec![Edit {
                    start: start as u32,
                    end: end as u32,
                    alt: replacement.to_vec(),
                }];
                if !out.iter().any(|o| o.edits == v.edits) {
                    out.push(v);
                }
            }
        }
        out
    }

    fn knockoff_variants(
        &self,
        real: &[Annotation],
        ko: KnockoffOptions,
        stats: &mut UniProtLoadStats,
    ) -> Vec<Annotation> {
        let listed: Vec<usize> = self
            .features
            .iter()
            .filter(|f| f.kind == "sequence variant" && self.on_canonical(f))
            .filter_map(|f| f.location)
            .flat_map(|(b, e)| (b as usize).saturating_sub(1)..e as usize)
            .collect();
        let (out, untested) =
            knockoff_annotations(&self.sequence, &self.accession, real, &listed, ko);
        stats.variants_without_knockoff += untested;
        out
    }

    fn origin(&self, kind: EntryKind, feature: &str) -> EntryOrigin {
        EntryOrigin {
            kind,
            family: Arc::from(self.accession.as_str()),
            feature: Arc::from(feature),
            variant_delta_da: 0.0,
            knockoff_of: None,
            source: Some(Arc::from(UNIPROT)),
            tier: PriorTier::Curated,
        }
    }

    fn processing_admitted(&self, f: &UniProtFeature, opts: &UniProtOptions) -> bool {
        if f.evidence_keys.is_empty() && self.dataset.eq_ignore_ascii_case("Swiss-Prot") {
            opts.evidence >= EvidenceLevel::Curated
        } else {
            opts.evidence.admits(self.evidence_codes(f))
        }
    }

    fn mature_chains(
        &self,
        opts: &UniProtOptions,
        stats: &mut UniProtLoadStats,
    ) -> Vec<MatureChain> {
        let len = self.sequence.len();
        let mut out: Vec<MatureChain> = Vec::new();
        for f in &self.features {
            if !matches!(f.kind.as_str(), "chain" | "peptide") || !self.on_canonical(f) {
                continue;
            }
            let Some((begin, end)) = f
                .location
                .filter(|&(b, e)| b >= 1 && b <= e && (e as usize) <= len)
            else {
                stats.unusable_features += 1;
                continue;
            };
            let (start, end_excl) = (begin as usize - 1, end as usize);
            let ends = ChainEnds {
                n_term: start > 0 && !(start == 1 && self.sequence[0] == b'M'),
                c_term: end_excl < len,
            };
            if !ends.n_term && !ends.c_term {
                continue;
            }
            let cut_here = |g: &UniProtFeature| {
                PROCESSING_FEATURES.contains(&g.kind.as_str())
                    && self.on_canonical(g)
                    && g.location.is_some_and(|(gb, ge)| {
                        (ends.n_term && ge + 1 == begin) || (ends.c_term && gb == end + 1)
                    })
            };
            let admitted = self.processing_admitted(f, opts)
                || self
                    .features
                    .iter()
                    .any(|g| cut_here(g) && self.processing_admitted(g, opts));
            if !admitted {
                stats.below_evidence += 1;
                continue;
            }
            let tag = if f.id.is_empty() {
                format!("CHAIN_{begin}-{end}")
            } else {
                f.id.clone()
            };
            if !out.iter().any(|c| c.start == start && c.end == end_excl) {
                out.push(MatureChain {
                    start,
                    end: end_excl,
                    ends,
                    tag,
                });
            }
        }
        out
    }

    fn isoform_edits(&self, iso: &UniProtIsoform) -> Option<Vec<Edit>> {
        if iso.sequence_type != "described" || iso.splice_refs.is_empty() {
            return None;
        }
        let mut edits: Vec<(usize, usize, &[u8])> = Vec::new();
        for r in &iso.splice_refs {
            let f = self
                .features
                .iter()
                .find(|f| f.kind == "splice variant" && f.id == *r && self.on_canonical(f))?;
            let (begin, end) = f
                .location
                .filter(|&(b, e)| b >= 1 && b <= e && (e as usize) <= self.sequence.len())?;
            let (start, end) = (begin as usize - 1, end as usize);
            if !f.original.is_empty() && f.original.as_bytes() != &self.sequence[start..end] {
                return None;
            }
            let replacement = f.variations.first().map_or(&[][..], |v| v.as_bytes());
            if !replacement.iter().all(u8::is_ascii_uppercase) {
                return None;
            }
            edits.push((start, end, replacement));
        }
        edits.sort_by_key(|e| e.0);
        if edits.windows(2).any(|w| w[1].0 < w[0].1) {
            return None;
        }
        let mut sequence = Vec::with_capacity(self.sequence.len());
        let mut at = 0;
        for &(start, end, replacement) in &edits {
            sequence.extend_from_slice(&self.sequence[at..start]);
            sequence.extend_from_slice(replacement);
            at = end;
        }
        sequence.extend_from_slice(&self.sequence[at..]);
        if sequence.is_empty() || sequence == self.sequence {
            return None;
        }
        Some(
            edits
                .iter()
                .map(|&(start, end, alt)| Edit {
                    start: start as u32,
                    end: end as u32,
                    alt: alt.to_vec(),
                })
                .collect(),
        )
    }

    pub fn into_proteins(
        self,
        opts: &UniProtOptions,
        stats: &mut UniProtLoadStats,
    ) -> Vec<Protein> {
        let accession = self.fasta_accession();
        let description = self.fasta_description();
        let sites = self.mod_sites(None, &self.sequence, opts, stats);
        stats.mod_sites += sites.len();
        let mut out = vec![Protein {
            accession,
            description: description.clone(),
            sequence: self.sequence.clone(),
            annotation: Some(Arc::new(
                ProteinAnnotation::new(sites.clone(), None)
                    .with_origin(Some(self.origin(EntryKind::Reference, "")))
                    .with_site_claims(
                        sites
                            .iter()
                            .map(|m| SiteClaim {
                                position: m.position,
                                mass_delta: m.mass_delta,
                                tier: PriorTier::Curated,
                                knockoff_of: None,
                            })
                            .collect(),
                    ),
            )),
        }];

        let canonical = out[0].clone();
        let build = |reference: &Protein, a: &Annotation| {
            build_entry(reference, a).expect("an alternative kind makes an entry")
        };
        if opts.variants {
            let variants = self.sequence_variants(None, &self.sequence, opts, stats);
            let knockoffs = match opts.knockoffs.filter(|k| k.per_variant > 0) {
                Some(ko) => self.knockoff_variants(&variants, ko, stats),
                None => Vec::new(),
            };
            stats.variant_proteins += variants.len();
            out.extend(variants.iter().map(|a| build(&canonical, a)));
            stats.knockoff_proteins += knockoffs.len();
            out.extend(knockoffs.iter().map(|a| build(&canonical, a)));
        }

        if opts.chains {
            for c in self.mature_chains(opts, stats) {
                let mut a = Annotation::new(&c.tag, EntryKind::Chain);
                a.source = Some(Arc::from(UNIPROT));
                a.tier = PriorTier::Curated;
                a.processing = Some(Processing {
                    start: c.start as u32,
                    end: c.end as u32,
                    ends: c.ends,
                });
                out.push(build(&canonical, &a));
                stats.chain_proteins += 1;
            }
        }

        let mut built: Vec<&str> = Vec::new();
        if opts.isoforms {
            for iso in self
                .isoforms
                .iter()
                .filter(|i| i.sequence_type != "displayed")
            {
                let Some(edits) = self.isoform_edits(iso) else {
                    stats.unbuilt_isoforms += 1;
                    continue;
                };
                built.push(&iso.id);
                let mut a = Annotation::new(&iso.id, EntryKind::Isoform);
                a.source = Some(Arc::from(UNIPROT));
                a.tier = PriorTier::Curated;
                a.edits = edits;
                a.accession = Some(format!("{}|{}|{}", self.db_code(), iso.id, self.name));
                a.description = Some(self.isoform_description(iso));
                let isoform = {
                    let sequence = build(&canonical, &a).sequence;
                    let own = self.mod_sites(Some(&iso.id), &sequence, opts, stats);
                    stats.mod_sites += own.len();
                    a.mods = own;
                    build(&canonical, &a)
                };
                let variants = if opts.variants {
                    self.sequence_variants(Some(&iso.id), &isoform.sequence, opts, stats)
                } else {
                    Vec::new()
                };
                stats.variant_proteins += variants.len();
                let proteins: Vec<Protein> = variants.iter().map(|v| build(&isoform, v)).collect();
                out.push(isoform);
                out.extend(proteins);
                stats.isoform_proteins += 1;
            }
        }

        stats.isoform_features += self
            .features
            .iter()
            .filter(|f| {
                let used_kind = matches!(
                    f.kind.as_str(),
                    "modified residue"
                        | "lipid moiety-binding region"
                        | "glycosylation site"
                        | "cross-link"
                ) || (opts.variants && f.kind == "sequence variant");
                used_kind && !self.on_canonical(f) && !built.contains(&f.location_sequence.as_str())
            })
            .count();
        out
    }
}

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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Capture {
    Accession,
    EntryName,
    RecommendedFullName,
    SubmittedFullName,
    Gene,
    Organism,
    IsoformId,
    IsoformName,
    Original,
    Variation,
    Sequence,
}

pub struct UniProtXmlReader<R: BufRead> {
    xml: Reader<R>,
    buf: Vec<u8>,
    stack: Vec<Tag>,
    entry: UniProtEntry,
    feature: UniProtFeature,
    begin: Option<u32>,
    end: Option<u32>,
    in_alt_products: bool,
    isoform: UniProtIsoform,
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
            isoform: UniProtIsoform::default(),
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

    fn open_element(&self) -> &'static str {
        self.stack
            .iter()
            .rev()
            .find(|t| **t != Tag::Other)
            .map_or("…", |t| t.name())
    }

    fn ancestor(&self, up: usize) -> Option<Tag> {
        self.stack.len().checked_sub(1 + up).map(|i| self.stack[i])
    }

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
                self.isoform = UniProtIsoform::default();
            }
            (Tag::Id, Some(Tag::Isoform)) if self.isoform.id.is_empty() => {
                self.capture = capture(Capture::IsoformId);
            }
            (Tag::Name, Some(Tag::Isoform)) if self.isoform.name.is_empty() => {
                self.capture = capture(Capture::IsoformName);
            }
            (Tag::Sequence, Some(Tag::Isoform)) => {
                self.isoform.sequence_type = attr(e, b"type")?.unwrap_or_default();
                self.isoform.splice_refs = attr(e, b"ref")?
                    .map(|r| r.split_whitespace().map(str::to_string).collect())
                    .unwrap_or_default();
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
                    Capture::IsoformId => self.isoform.id = text,
                    Capture::IsoformName => self.isoform.name = text,
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
                let iso = std::mem::take(&mut self.isoform);
                if !iso.id.is_empty() {
                    self.entry.isoforms.push(iso);
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
#[path = "uniprot/tests/uniprot.rs"]
mod tests;
