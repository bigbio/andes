#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ModLocation {
    Anywhere,
    NTerm,
    CTerm,
    ProtNTerm,
    ProtCTerm,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ResidueSpec {
    Specific(u8),
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
    pub neutral_losses: Vec<f64>,
    pub loss_class: u8,
    pub annotated_site_mass: Option<f64>,
    pub blocks_cleavage: bool,
    pub diagnostic_ions: Vec<f64>,
    pub parts: Vec<&'static Modification>,
}

impl Modification {
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

    pub fn components(&self) -> Vec<&Modification> {
        if self.parts.is_empty() {
            vec![self]
        } else {
            self.parts.iter().map(|m| &**m).collect()
        }
    }
}

impl Modification {
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

    pub fn allowed_at(&self, sites: &[crate::protein::ModSite], position: u32) -> bool {
        match self.annotated_site_mass {
            None => true,
            Some(mass) => crate::protein::has_site(sites, position, mass),
        }
    }
}

pub fn loss_class_id(name: &str) -> Option<u8> {
    match name.trim().to_ascii_lowercase().as_str() {
        "glyco" => Some(1),
        "phospho" => Some(2),
        "sulfo" => Some(3),
        "generic" => Some(255),
        _ => None,
    }
}

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
        assert!(m.applies_to(b'C', ModLocation::NTerm));
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
        assert!(m.applies_to(b'A', ModLocation::ProtNTerm));
        assert!(m.applies_to(b'M', ModLocation::ProtNTerm));
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
                mass_delta: 79.956815,
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
        assert!(free.allowed_at(&[], 4));
    }
}
