//! Modifications. The lists users write (`Name (sites)`) are parsed in
//! [`crate::mod_spec`].

/// Where a modification can attach within (or at the ends of) a peptide.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ModLocation {
    /// Any internal or terminal position. Subsumes the four terminal
    /// locations for matching purposes.
    Anywhere,
    /// Peptide N-terminus (any residue), but not protein N-terminus.
    NTerm,
    /// Peptide C-terminus (any residue), but not protein C-terminus.
    CTerm,
    /// Protein N-terminus (only when the residue is the protein's first AA).
    ProtNTerm,
    /// Protein C-terminus (only when the residue is the protein's last AA).
    ProtCTerm,
}

/// Which residues a modification is allowed to target.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ResidueSpec {
    /// Exactly one residue (e.g. `b'C'` for Carbamidomethyl).
    Specific(u8),
    /// Any residue (e.g. terminal-only mods like protein-N-term Acetyl).
    Wildcard,
}

#[derive(Debug, Clone)]
pub struct Modification {
    pub name: String,
    pub mass_delta: f64,
    pub residue: ResidueSpec,
    pub location: ModLocation,
    pub fixed: bool,
    pub accession: Option<String>,
    /// User-declared neutral-loss masses (Da) for this mod's fragment ions.
    /// Empty ⇒ no loss ions predicted (default; byte-identical to pre-feature).
    pub neutral_losses: Vec<f64>,
    /// Loss-class pool id for this mod's neutral-loss ions: 0 = none/intact,
    /// 1.. = a per-mod-class pool (glyco=1, phospho=2, sulfo=3, generic=255).
    /// Loss ions of the same class share one trained distribution.
    pub loss_class: u8,
    /// `Some(delta)` restricts a variable mod to the positions the protein
    /// database annotates with a modification of mass `delta` (config
    /// `restrict_to: annotated`; the annotations come from a UniProt XML database, see
    /// [`crate::protein::ProteinAnnotation`]). `None` (the default) places the
    /// mod wherever its residue and location allow. For a fixed-anywhere mod
    /// stacked with a restricted terminal mod, `delta` is the restricted mod's
    /// own mass, not the stacked sum.
    pub annotated_site_mass: Option<f64>,
    /// The enzyme does not cut at a residue carrying this mod (config
    /// `cleavage: blocked`), so the mod never sits on the residue a peptide was
    /// cleaved at: its C-terminal residue for an enzyme that cuts after it (its
    /// N-terminal residue for one that cuts before), unless that residue is the
    /// protein's own terminus. True of most lysine acylations, GlyGly and di-
    /// and trimethylation under trypsin; not of labels added after digestion.
    pub blocks_cleavage: bool,
    /// Singly charged ions this mod gives off when fragmented (config
    /// `diagnostic_ions`), such as the immonium ion of the modified residue.
    /// Evidence for which modification a site carries; not used in scoring.
    pub diagnostic_ions: Vec<f64>,
    /// For a record that stacks several declared modifications on one
    /// residue (a terminal and a side-chain mod, see
    /// [`crate::aa_set::AminoAcidSetBuilder::build`]), those modifications in
    /// stacking order: the record itself has no accession, its parts do.
    /// Empty for a single modification.
    pub parts: Vec<&'static Modification>,
}

impl Modification {
    /// A modification of `mass_delta` Da on `residue` at `location`, with no
    /// accession, neutral losses, diagnostic ions, site restriction or
    /// cleavage blocking.
    pub fn new(
        name: impl Into<String>,
        mass_delta: f64,
        residue: ResidueSpec,
        location: ModLocation,
        fixed: bool,
    ) -> Self {
        Modification {
            name: name.into(),
            mass_delta,
            residue,
            location,
            fixed,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
            annotated_site_mass: None,
            blocks_cleavage: false,
            diagnostic_ions: Vec::new(),
            parts: Vec::new(),
        }
    }

    /// The single modifications this record stands for: its `parts` when it
    /// is stacked, else itself.
    pub fn components(&self) -> Vec<&Modification> {
        if self.parts.is_empty() {
            vec![self]
        } else {
            self.parts.iter().map(|m| &**m).collect()
        }
    }
}

impl Modification {
    /// Test whether this mod is allowed on `residue` at the given
    /// `location`. `Anywhere`-targeting mods match any of the four
    /// non-Anywhere locations; otherwise the mod's `location` must equal
    /// the queried location exactly.
    pub fn applies_to(&self, residue: u8, location: ModLocation) -> bool {
        let residue_ok = match self.residue {
            ResidueSpec::Specific(r) => r == residue,
            ResidueSpec::Wildcard => true,
        };
        let location_ok = match (self.location, location) {
            (ModLocation::Anywhere, _) => true,
            (a, b) => a == b,
        };
        residue_ok && location_ok
    }

    /// Whether this mod may sit at a position annotated with `sites`: always
    /// for an unrestricted mod; for a `restrict_to: annotated` mod only when one of
    /// `sites` carries a modification of the same mass.
    pub fn allowed_at(&self, sites: &[crate::protein::ModSite], position: u32) -> bool {
        match self.annotated_site_mass {
            None => true,
            Some(mass) => crate::protein::has_site(sites, position, mass),
        }
    }
}

/// Return the stable loss-class id for a recognised class name, or `None`.
///
/// Loss ions are pooled within a class so that a shared trained distribution
/// can be learned across mods of the same biochemical type.
///
/// | name      | id  |
/// |-----------|-----|
/// | `glyco`   |   1 |
/// | `phospho` |   2 |
/// | `sulfo`   |   3 |
/// | `generic` | 255 |
pub fn loss_class_id(name: &str) -> Option<u8> {
    match name.trim().to_ascii_lowercase().as_str() {
        "glyco" => Some(1),
        "phospho" => Some(2),
        "sulfo" => Some(3),
        "generic" => Some(255),
        _ => None,
    }
}

/// Hand a `Modification` out as a `'static` reference by leaking it. The
/// per-search modification table is a few dozen small records built once, so
/// the leak is negligible; in exchange every `AminoAcid` clone during candidate
/// enumeration copies a pointer instead of bumping a shared atomic refcount
/// (which all worker threads contended on: 16 threads burned 4.7x the CPU of
/// one thread for identical output on a phospho search).
pub fn leak_mod(m: Modification) -> &'static Modification {
    Box::leak(Box::new(m))
}

impl From<Modification> for &'static Modification {
    fn from(m: Modification) -> Self {
        leak_mod(m)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn carbamidomethyl_c() -> Modification {
        Modification {
            name: "Carbamidomethyl".to_string(),
            mass_delta: 57.02146,
            residue: ResidueSpec::Specific(b'C'),
            location: ModLocation::Anywhere,
            fixed: true,
            accession: Some("UNIMOD:4".to_string()),
            neutral_losses: Vec::new(),
            loss_class: 0,
            annotated_site_mass: None,
            blocks_cleavage: false,
            parts: Vec::new(),
            diagnostic_ions: Vec::new(),
        }
    }

    #[test]
    fn applies_to_matching_residue_anywhere() {
        let m = carbamidomethyl_c();
        assert!(m.applies_to(b'C', ModLocation::Anywhere));
        assert!(m.applies_to(b'C', ModLocation::NTerm)); // Anywhere subsumes
        assert!(m.applies_to(b'C', ModLocation::CTerm));
    }

    #[test]
    fn applies_to_wrong_residue() {
        let m = carbamidomethyl_c();
        assert!(!m.applies_to(b'A', ModLocation::Anywhere));
    }

    #[test]
    fn applies_to_wildcard_residue() {
        let m = Modification {
            name: "Acetyl".to_string(),
            mass_delta: 42.01057,
            residue: ResidueSpec::Wildcard,
            location: ModLocation::ProtNTerm,
            fixed: false,
            accession: Some("UNIMOD:1".to_string()),
            neutral_losses: Vec::new(),
            loss_class: 0,
            annotated_site_mass: None,
            blocks_cleavage: false,
            parts: Vec::new(),
            diagnostic_ions: Vec::new(),
        };
        // Wildcard matches any residue at the specified location only.
        assert!(m.applies_to(b'A', ModLocation::ProtNTerm));
        assert!(m.applies_to(b'M', ModLocation::ProtNTerm));
        // ...but not at other locations.
        assert!(!m.applies_to(b'A', ModLocation::Anywhere));
        assert!(!m.applies_to(b'A', ModLocation::NTerm));
    }

    #[test]
    fn applies_to_specific_location() {
        let m = Modification {
            name: "TMT6plex".to_string(),
            mass_delta: 229.16293,
            residue: ResidueSpec::Specific(b'K'),
            location: ModLocation::Anywhere,
            fixed: true,
            accession: Some("UNIMOD:737".to_string()),
            neutral_losses: Vec::new(),
            loss_class: 0,
            annotated_site_mass: None,
            blocks_cleavage: false,
            parts: Vec::new(),
            diagnostic_ions: Vec::new(),
        };
        assert!(m.applies_to(b'K', ModLocation::Anywhere));
        assert!(!m.applies_to(b'R', ModLocation::Anywhere));
    }

    #[test]
    fn applies_to_nterm_only() {
        let m = Modification {
            name: "TMT6plex_NT".to_string(),
            mass_delta: 229.16293,
            residue: ResidueSpec::Wildcard,
            location: ModLocation::NTerm,
            fixed: true,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
            annotated_site_mass: None,
            blocks_cleavage: false,
            parts: Vec::new(),
            diagnostic_ions: Vec::new(),
        };
        assert!(m.applies_to(b'A', ModLocation::NTerm));
        assert!(!m.applies_to(b'A', ModLocation::Anywhere));
        assert!(!m.applies_to(b'A', ModLocation::CTerm));
    }

    #[test]
    fn allowed_at_needs_an_annotated_site_of_the_same_mass() {
        use crate::protein::{ModSite, SiteAnchor};
        let free = Modification::new(
            "Phospho",
            79.966331,
            ResidueSpec::Specific(b'S'),
            ModLocation::Anywhere,
            false,
        );
        let phospho = Modification {
            annotated_site_mass: Some(79.966331),
            ..free.clone()
        };
        let sites = [
            ModSite {
                position: 3,
                mass_delta: 79.966331,
                anchor: SiteAnchor::Residue,
            },
            ModSite {
                position: 7,
                mass_delta: 79.956815, // sulfation: 9.5 mDa lighter
                anchor: SiteAnchor::Residue,
            },
        ];
        assert!(phospho.allowed_at(&sites, 3));
        assert!(!phospho.allowed_at(&sites, 4), "no site at 4");
        assert!(
            !phospho.allowed_at(&sites, 7),
            "a sulfo site is not a phospho site"
        );
        assert!(!phospho.allowed_at(&[], 3), "no annotation, no site");
        // An unrestricted mod is allowed anywhere.
        assert!(free.allowed_at(&[], 4));
    }
}
