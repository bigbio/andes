//! Annotations of a reference sequence: the claims an alternative entry is
//! built from, whatever database they come from.

use std::sync::Arc;

use crate::protein::{ChainEnds, EntryKind, ModSite, PriorTier};

/// One claim about a reference protein: sequence edits, modification sites
/// and processing that occur together on one molecule, or any one of them
/// alone. The entry it makes is built by the database reader from the
/// reference protein and this annotation.
#[derive(Debug, Clone, PartialEq)]
pub struct Annotation {
    /// Stable id, used in the accession and as the entry's feature
    /// (`VAR_005851_D7H`, `P04637-2`, `PRO_0000001058`).
    pub id: Arc<str>,
    /// The kind of entry it makes: `Variant`, `Isoform` or `Chain`.
    pub kind: EntryKind,
    /// Sequence changes on the reference, sorted and disjoint.
    pub edits: Vec<Edit>,
    /// Modification sites of the entry itself, in the coordinates after the
    /// edits, besides the reference's sites carried through them.
    pub mods: Vec<ModSite>,
    /// For a mature chain, the span cut out of the reference.
    pub processing: Option<Processing>,
    /// For a knockoff, the id of the annotation it copies.
    pub knockoff_of: Option<Arc<str>>,
    /// The entry's accession when the default (the reference accession with
    /// `_<id>` added to its middle field) does not fit.
    pub accession: Option<String>,
    /// The entry's description when the default does not fit.
    pub description: Option<String>,
    /// The database the claim comes from.
    pub source: Option<Arc<str>>,
    /// How likely claims of its kind are to be real.
    pub tier: PriorTier,
}

/// Reference residues `[start, end)` (0-based) replaced by `alt`; an empty
/// `alt` is a deletion.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Edit {
    pub start: u32,
    pub end: u32,
    pub alt: Vec<u8>,
}

/// A mature chain: reference residues `[start, end)`, with the ends that are
/// new termini.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Processing {
    pub start: u32,
    pub end: u32,
    pub ends: ChainEnds,
}

impl Annotation {
    /// An annotation of `kind` named `id`, with nothing else set.
    pub fn new(id: &str, kind: EntryKind) -> Self {
        Annotation {
            id: Arc::from(id),
            kind,
            edits: Vec::new(),
            mods: Vec::new(),
            processing: None,
            knockoff_of: None,
            accession: None,
            description: None,
            source: None,
            tier: PriorTier::Catalogue,
        }
    }

    /// Whether this is a single-residue substitution.
    pub fn is_single_substitution(&self) -> bool {
        matches!(self.edits.as_slice(), [e] if e.end == e.start + 1 && e.alt.len() == 1)
    }
}
