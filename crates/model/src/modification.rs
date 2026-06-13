//! Modifications and the Mods.txt parser.

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
    /// User-declared neutral-loss masses (Da) for this mod's fragment ions.
    /// Empty ⇒ no loss ions predicted (default; byte-identical to pre-feature).
    pub neutral_losses: Vec<f64>,
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

#[derive(thiserror::Error, Debug)]
pub enum ModParseError {
    #[error("expected 5 comma-separated fields, got {got}")]
    WrongFieldCount { got: usize },
    #[error("invalid mass delta {field:?}: {source}")]
    BadMass { field: String, #[source] source: std::num::ParseFloatError },
    #[error("invalid residue spec {field:?} (expected single ASCII upper char or `*`)")]
    BadResidue { field: String },
    #[error("invalid location {field:?} (expected `any|N-term|C-term|Prot-N-term|Prot-C-term`)")]
    BadLocation { field: String },
    #[error("invalid fixed/variable flag {field:?} (expected `fix|opt`)")]
    BadFixedFlag { field: String },
    #[error("unknown mod attribute key {key:?} (expected loss|accession)")]
    UnknownModAttr { key: String },
    #[error("malformed mod attribute {field:?} (expected key=value)")]
    BadModAttr { field: String },
    #[error("invalid neutral-loss value {value:?} (expected positive number < 2000)")]
    BadNeutralLoss { value: String },
}

impl Modification {
    /// Parse a single non-empty, non-comment line from a Mods.txt file.
    /// Empty lines and `# ...` comment lines should be filtered by the
    /// caller (see `aa_set::AminoAcidSetBuilder::add_mods_from_file`).
    pub fn from_mods_txt_line(line: &str) -> Result<Self, ModParseError> {
        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() < 5 {
            return Err(ModParseError::WrongFieldCount { got: parts.len() });
        }
        let (mass_s, residues_s, fixity_s, location_s, name_s) = (
            parts[0].trim(), parts[1].trim(), parts[2].trim(),
            parts[3].trim(), parts[4].trim(),
        );

        let mass_delta: f64 = mass_s.parse()
            .map_err(|source| ModParseError::BadMass { field: mass_s.to_string(), source })?;

        let residue = match residues_s {
            "*" => ResidueSpec::Wildcard,
            s if s.len() == 1 && s.as_bytes()[0].is_ascii_uppercase() => {
                ResidueSpec::Specific(s.as_bytes()[0])
            }
            _ => return Err(ModParseError::BadResidue { field: residues_s.to_string() }),
        };

        let fixed = match fixity_s.to_ascii_lowercase().as_str() {
            "fix" => true,
            "opt" => false,
            _ => return Err(ModParseError::BadFixedFlag { field: fixity_s.to_string() }),
        };

        let location = match location_s.to_ascii_lowercase().as_str() {
            "any"          => ModLocation::Anywhere,
            "n-term"       => ModLocation::NTerm,
            "c-term"       => ModLocation::CTerm,
            "prot-n-term"  => ModLocation::ProtNTerm,
            "prot-c-term"  => ModLocation::ProtCTerm,
            _ => return Err(ModParseError::BadLocation { field: location_s.to_string() }),
        };

        let mut neutral_losses: Vec<f64> = Vec::new();
        let mut accession: Option<String> = None;
        for attr in &parts[5..] {
            let attr = attr.trim();
            if attr.is_empty() { continue; }
            let (key, value) = attr.split_once('=')
                .ok_or_else(|| ModParseError::BadModAttr { field: attr.to_string() })?;
            match key.trim().to_ascii_lowercase().as_str() {
                "loss" => {
                    for tok in value.split(';') {
                        let tok = tok.trim();
                        if tok.is_empty() { continue; }
                        let v: f64 = tok.parse()
                            .map_err(|_| ModParseError::BadNeutralLoss { value: tok.to_string() })?;
                        if !(v > 0.0 && v < 2000.0) {
                            return Err(ModParseError::BadNeutralLoss { value: tok.to_string() });
                        }
                        neutral_losses.push(v);
                    }
                }
                "accession" => accession = Some(value.trim().to_string()),
                other => return Err(ModParseError::UnknownModAttr { key: other.to_string() }),
            }
        }

        Ok(Modification {
            name: name_s.to_string(),
            mass_delta,
            residue,
            location,
            fixed,
            accession,
            neutral_losses,
        })
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
        }
    }

    #[allow(dead_code)]
    fn oxidation_m() -> Modification {
        Modification {
            name: "Oxidation".to_string(),
            mass_delta: 15.99491,
            residue: ResidueSpec::Specific(b'M'),
            location: ModLocation::Anywhere,
            fixed: false,
            accession: Some("UNIMOD:35".to_string()),
            neutral_losses: Vec::new(),
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
            neutral_losses: Vec::new(),
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
        };
        assert!(m.applies_to(b'A', ModLocation::NTerm));
        assert!(!m.applies_to(b'A', ModLocation::Anywhere));
        assert!(!m.applies_to(b'A', ModLocation::CTerm));
    }

    #[test]
    fn parse_carbamidomethyl_c() {
        let line = "57.021464,C,fix,any,Carbamidomethyl";
        let m = Modification::from_mods_txt_line(line).unwrap();
        assert_eq!(m.name, "Carbamidomethyl");
        assert_eq!(m.mass_delta, 57.021464);
        assert_eq!(m.residue, ResidueSpec::Specific(b'C'));
        assert_eq!(m.location, ModLocation::Anywhere);
        assert!(m.fixed);
    }

    #[test]
    fn parse_oxidation_m_variable() {
        let line = "15.994915,M,opt,any,Oxidation";
        let m = Modification::from_mods_txt_line(line).unwrap();
        assert!(!m.fixed);
        assert_eq!(m.mass_delta, 15.994915);
    }

    #[test]
    fn parse_wildcard_nterm() {
        let line = "229.162932,*,fix,N-term,TMT6plex";
        let m = Modification::from_mods_txt_line(line).unwrap();
        assert_eq!(m.residue, ResidueSpec::Wildcard);
        assert_eq!(m.location, ModLocation::NTerm);
    }

    #[test]
    fn parse_protein_nterm_acetyl() {
        let line = "42.010565,*,opt,Prot-N-term,Acetyl";
        let m = Modification::from_mods_txt_line(line).unwrap();
        assert_eq!(m.location, ModLocation::ProtNTerm);
    }

    #[test]
    fn parse_negative_mass_delta() {
        let line = "-17.026549,Q,opt,N-term,Pyro-glu";
        let m = Modification::from_mods_txt_line(line).unwrap();
        assert_eq!(m.mass_delta, -17.026549);
    }

    #[test]
    fn parse_wrong_field_count() {
        let line = "57.021464,C,fix,any";  // 4 fields
        let err = Modification::from_mods_txt_line(line).unwrap_err();
        assert!(matches!(err, ModParseError::WrongFieldCount { got: 4 }));
    }

    #[test]
    fn parse_bad_mass() {
        let line = "abc,C,fix,any,Bad";
        let err = Modification::from_mods_txt_line(line).unwrap_err();
        assert!(matches!(err, ModParseError::BadMass { .. }));
    }

    #[test]
    fn parse_bad_residue() {
        let line = "57.0,CC,fix,any,Bad";
        let err = Modification::from_mods_txt_line(line).unwrap_err();
        assert!(matches!(err, ModParseError::BadResidue { .. }));
    }

    #[test]
    fn parse_bad_location() {
        let line = "57.0,C,fix,middle,Bad";
        let err = Modification::from_mods_txt_line(line).unwrap_err();
        assert!(matches!(err, ModParseError::BadLocation { .. }));
    }

    #[test]
    fn parse_bad_fixity() {
        let line = "57.0,C,maybe,any,Bad";
        let err = Modification::from_mods_txt_line(line).unwrap_err();
        assert!(matches!(err, ModParseError::BadFixedFlag { .. }));
    }

    #[test]
    fn parse_location_case_insensitive() {
        let line = "229.162932,*,fix,n-term,TMT";
        let m = Modification::from_mods_txt_line(line).unwrap();
        assert_eq!(m.location, ModLocation::NTerm);
    }

    #[test]
    fn parses_loss_and_accession_attributes() {
        let m = Modification::from_mods_txt_line(
            "340.100562,K,opt,any,Glucosylgalactosyl,loss=162.0528;324.1056,accession=UNIMOD:393"
        ).unwrap();
        assert_eq!(m.residue, ResidueSpec::Specific(b'K'));
        assert!(!m.fixed);
        assert_eq!(m.neutral_losses, vec![162.0528, 324.1056]);
        assert_eq!(m.accession.as_deref(), Some("UNIMOD:393"));
    }

    #[test]
    fn five_field_line_has_no_losses_or_accession() {
        let m = Modification::from_mods_txt_line("57.02146,C,fix,any,Carbamidomethyl").unwrap();
        assert!(m.neutral_losses.is_empty());
        assert_eq!(m.accession, None);
    }

    #[test]
    fn rejects_unknown_attr_and_bad_loss() {
        assert!(matches!(Modification::from_mods_txt_line("1.0,K,opt,any,X,frobnicate=7"), Err(ModParseError::UnknownModAttr { .. })));
        assert!(matches!(Modification::from_mods_txt_line("1.0,K,opt,any,X,loss=abc"), Err(ModParseError::BadNeutralLoss { .. })));
        assert!(matches!(Modification::from_mods_txt_line("1.0,K,opt,any,X,nokey"), Err(ModParseError::BadModAttr { .. })));
    }
}
