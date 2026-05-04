//! Modifications. Mirrors Java
//! `edu.ucsd.msjava.msutil.Modification` and the Mods.txt parser in
//! `edu.ucsd.msjava.params.ParamObject`.

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
    pub name:       String,
    pub mass_delta: f64,
    pub residue:    ResidueSpec,
    pub location:   ModLocation,
    pub fixed:      bool,
    pub accession:  Option<String>,
}

impl Modification {
    /// Test whether this mod is allowed on `residue` at the given
    /// `location`. `Anywhere`-targeting mods match any of the four
    /// non-Anywhere locations; otherwise the mod's `location` must equal
    /// the queried location exactly.
    pub fn applies_to(&self, residue: u8, location: ModLocation) -> bool {
        let residue_ok = match self.residue {
            ResidueSpec::Specific(r) => r == residue,
            ResidueSpec::Wildcard    => true,
        };
        let location_ok = match (self.location, location) {
            (ModLocation::Anywhere, _) => true,
            (a, b) => a == b,
        };
        residue_ok && location_ok
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
        }
    }

    fn oxidation_m() -> Modification {
        Modification {
            name: "Oxidation".to_string(),
            mass_delta: 15.99491,
            residue: ResidueSpec::Specific(b'M'),
            location: ModLocation::Anywhere,
            fixed: false,
            accession: Some("UNIMOD:35".to_string()),
        }
    }

    #[test]
    fn applies_to_matching_residue_anywhere() {
        let m = carbamidomethyl_c();
        assert!(m.applies_to(b'C', ModLocation::Anywhere));
        assert!(m.applies_to(b'C', ModLocation::NTerm));   // Anywhere subsumes
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
        };
        assert!(m.applies_to(b'A', ModLocation::NTerm));
        assert!(!m.applies_to(b'A', ModLocation::Anywhere));
        assert!(!m.applies_to(b'A', ModLocation::CTerm));
    }
}
