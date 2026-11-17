use std::sync::Arc;

use crate::protein::{ChainEnds, EntryKind, ModSite, PriorTier};

#[derive(Debug, Clone, PartialEq)]
pub struct Annotation {
    pub id: Arc<str>,
    pub kind: EntryKind,
    pub edits: Vec<Edit>,
    pub mods: Vec<ModSite>,
    pub processing: Option<Processing>,
    pub knockoff_of: Option<Arc<str>>,
    pub accession: Option<String>,
    pub description: Option<String>,
    pub source: Option<Arc<str>>,
    pub tier: PriorTier,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Edit {
    pub start: u32,
    pub end: u32,
    pub alt: Vec<u8>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Processing {
    pub start: u32,
    pub end: u32,
    pub ends: ChainEnds,
}

impl Annotation {
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

    pub fn is_single_substitution(&self) -> bool {
        matches!(self.edits.as_slice(), [e] if e.end == e.start + 1 && e.alt.len() == 1)
    }
}
