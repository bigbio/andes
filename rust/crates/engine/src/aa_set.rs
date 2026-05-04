//! Heavyweight residue-and-modification set. Mirrors Java
//! `edu.ucsd.msjava.msutil.AminoAcidSet`. Built via
//! `AminoAcidSetBuilder`; queried by Phase 4's candidate generator.

use std::collections::HashMap;
use std::fs;
use std::path::Path;

use crate::amino_acid::AminoAcid;
use crate::modification::{ModLocation, ModParseError, Modification, ResidueSpec};

const STANDARD_RESIDUES: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
const IMPLAUSIBLE_MASS_THRESHOLD: f64 = 1000.0;

#[derive(Debug, Clone)]
pub struct AminoAcidSet {
    /// (residue, location) → all variants (unmodified + modified) at that position.
    table: HashMap<(u8, ModLocation), Vec<AminoAcid>>,
    has_cterm_mods: bool,
    min_aa_mass: f64,
    max_aa_mass: f64,
    max_residue_mod_mass: f64,
    max_fixed_term_mod_mass: f64,
}

impl AminoAcidSet {
    /// All variants of `residue` valid at the given `location`.
    pub fn variants_for(&self, residue: u8, location: ModLocation) -> &[AminoAcid] {
        self.table
            .get(&(residue, location))
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    pub fn standard(&self, residue: u8) -> Option<&AminoAcid> {
        self.variants_for(residue, ModLocation::Anywhere)
            .iter()
            .find(|aa| !aa.is_modified())
    }

    pub fn contains_cterm_mods(&self) -> bool { self.has_cterm_mods }
    pub fn min_aa_mass(&self) -> f64           { self.min_aa_mass }
    pub fn max_aa_mass(&self) -> f64           { self.max_aa_mass }
    pub fn max_residue_mod_mass(&self) -> f64  { self.max_residue_mod_mass }
    pub fn max_fixed_term_mod_mass(&self) -> f64 { self.max_fixed_term_mod_mass }

    pub fn iter_variants(&self) -> impl Iterator<Item = &AminoAcid> {
        self.table.values().flat_map(|v| v.iter())
    }
}

/// Accumulator. Each `add_*` call validates lazily; `build()` does final
/// checks and produces the immutable `AminoAcidSet`.
#[derive(Debug, Clone)]
pub struct AminoAcidSetBuilder {
    fixed_mods:    Vec<Modification>,
    variable_mods: Vec<Modification>,
}

impl AminoAcidSetBuilder {
    pub fn new_standard() -> Self {
        Self { fixed_mods: vec![], variable_mods: vec![] }
    }

    pub fn new_standard_with_carbamidomethyl_c() -> Self {
        let cam = Modification {
            name: "Carbamidomethyl".to_string(),
            mass_delta: 57.02146,
            residue: ResidueSpec::Specific(b'C'),
            location: ModLocation::Anywhere,
            fixed: true,
            accession: Some("UNIMOD:4".to_string()),
        };
        Self {
            fixed_mods: vec![cam],
            variable_mods: vec![],
        }
    }

    pub fn add_fixed_mod(mut self, m: Modification) -> Self {
        self.fixed_mods.push(m);
        self
    }

    pub fn add_variable_mod(mut self, m: Modification) -> Self {
        self.variable_mods.push(m);
        self
    }

    pub fn add_mods_from_file(mut self, path: &Path) -> Result<Self, AaSetError> {
        let text = fs::read_to_string(path)?;
        for (line_no, raw) in text.lines().enumerate() {
            let line = raw.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let m = Modification::from_mods_txt_line(line)
                .map_err(|source| AaSetError::ModsTxtParse { line_no: line_no + 1, source })?;
            if m.fixed {
                self.fixed_mods.push(m);
            } else {
                self.variable_mods.push(m);
            }
        }
        Ok(self)
    }

    pub fn build(self) -> Result<AminoAcidSet, AaSetError> {
        // 1. Reject implausible mod masses.
        for m in self.fixed_mods.iter().chain(self.variable_mods.iter()) {
            if m.mass_delta.abs() > IMPLAUSIBLE_MASS_THRESHOLD {
                return Err(AaSetError::ImplausibleMassDelta {
                    name: m.name.clone(),
                    delta: m.mass_delta,
                });
            }
        }

        // 2. Detect (residue, location) overlap between fixed and variable.
        for fm in &self.fixed_mods {
            for vm in &self.variable_mods {
                if mods_target_same_slot(fm, vm) {
                    let res_char = match fm.residue {
                        ResidueSpec::Specific(r) => r as char,
                        ResidueSpec::Wildcard    => '*',
                    };
                    return Err(AaSetError::ConflictingMods {
                        residue: res_char,
                        location: fm.location,
                    });
                }
            }
        }

        // 3. Build the table.
        let mut table: HashMap<(u8, ModLocation), Vec<AminoAcid>> = HashMap::new();
        let locations = [
            ModLocation::Anywhere, ModLocation::NTerm, ModLocation::CTerm,
            ModLocation::ProtNTerm, ModLocation::ProtCTerm,
        ];

        for &r in STANDARD_RESIDUES {
            let std_aa = AminoAcid::standard(r).expect("STANDARD_RESIDUES has only valid residues");

            for &loc in &locations {
                let fixed_match = self.fixed_mods
                    .iter()
                    .find(|m| m.applies_to(r, loc))
                    .cloned();

                let variable_matches: Vec<_> = self.variable_mods
                    .iter()
                    .filter(|m| m.applies_to(r, loc))
                    .cloned()
                    .collect();

                let mut variants = Vec::new();
                if loc == ModLocation::Anywhere {
                    if let Some(fm) = &fixed_match {
                        variants.push(std_aa.clone().with_mod(fm.clone()));
                    } else {
                        variants.push(std_aa.clone());
                    }
                    for vm in &variable_matches {
                        variants.push(std_aa.clone().with_mod(vm.clone()));
                    }
                } else {
                    if let Some(fm) = &fixed_match {
                        if fm.location == loc {
                            variants.push(std_aa.clone().with_mod(fm.clone()));
                        }
                    }
                    for vm in &variable_matches {
                        if vm.location == loc {
                            variants.push(std_aa.clone().with_mod(vm.clone()));
                        }
                    }
                }

                if !variants.is_empty() {
                    table.insert((r, loc), variants);
                }
            }
        }

        // 4. Aggregates.
        let standard_masses: Vec<f64> = STANDARD_RESIDUES.iter()
            .filter_map(|&r| AminoAcid::standard(r).map(|aa| aa.mass))
            .collect();
        let min_aa_mass = standard_masses.iter().copied().fold(f64::INFINITY, f64::min);
        let max_aa_mass = standard_masses.iter().copied().fold(f64::NEG_INFINITY, f64::max);

        let mut max_mod_delta = 0.0_f64;
        for m in self.fixed_mods.iter().chain(self.variable_mods.iter()) {
            if m.mass_delta > max_mod_delta {
                max_mod_delta = m.mass_delta;
            }
        }
        let max_residue_mod_mass = max_aa_mass + max_mod_delta;

        let max_fixed_term_mod_mass = self.fixed_mods
            .iter()
            .filter(|m| matches!(m.location,
                ModLocation::NTerm | ModLocation::CTerm |
                ModLocation::ProtNTerm | ModLocation::ProtCTerm))
            .map(|m| m.mass_delta)
            .fold(0.0_f64, f64::max);

        let has_cterm_mods = self.fixed_mods.iter().chain(self.variable_mods.iter())
            .any(|m| matches!(m.location, ModLocation::CTerm | ModLocation::ProtCTerm));

        Ok(AminoAcidSet {
            table,
            has_cterm_mods,
            min_aa_mass,
            max_aa_mass,
            max_residue_mod_mass,
            max_fixed_term_mod_mass,
        })
    }
}

/// Two mods target the same slot iff residue and location overlap after
/// wildcard expansion.
fn mods_target_same_slot(a: &Modification, b: &Modification) -> bool {
    let residue_overlap = match (a.residue, b.residue) {
        (ResidueSpec::Specific(x), ResidueSpec::Specific(y)) => x == y,
        (ResidueSpec::Wildcard, _) | (_, ResidueSpec::Wildcard) => true,
    };
    let location_overlap = match (a.location, b.location) {
        (ModLocation::Anywhere, _) | (_, ModLocation::Anywhere) => true,
        (x, y) => x == y,
    };
    residue_overlap && location_overlap
}

#[derive(thiserror::Error, Debug)]
pub enum AaSetError {
    #[error("conflicting fixed and variable mod for residue {residue:?} at {location:?}")]
    ConflictingMods { residue: char, location: ModLocation },
    #[error("mod {name:?} mass delta {delta} is implausible (>1000 Da)")]
    ImplausibleMassDelta { name: String, delta: f64 },
    #[error("malformed Mods.txt line {line_no}: {source}")]
    ModsTxtParse { line_no: usize, #[source] source: ModParseError },
    #[error("Mods.txt I/O error: {source}")]
    Io { #[from] source: std::io::Error },
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::amino_acid::AminoAcid;
    use crate::modification::{Modification, ModLocation, ResidueSpec};

    fn carbamidomethyl_c() -> Modification {
        Modification {
            name: "Carbamidomethyl".to_string(),
            mass_delta: 57.02146,
            residue: ResidueSpec::Specific(b'C'),
            location: ModLocation::Anywhere,
            fixed: true,
            accession: None,
        }
    }

    fn oxidation_m() -> Modification {
        Modification {
            name: "Oxidation".to_string(),
            mass_delta: 15.99491,
            residue: ResidueSpec::Specific(b'M'),
            location: ModLocation::Anywhere,
            fixed: false,
            accession: None,
        }
    }

    #[test]
    fn standard_set_has_20_residues() {
        let set = AminoAcidSetBuilder::new_standard().build().unwrap();
        let mut seen = std::collections::HashSet::new();
        for aa in set.iter_variants() {
            seen.insert(aa.residue);
        }
        assert_eq!(seen.len(), 20);
    }

    #[test]
    fn standard_set_no_mods() {
        let set = AminoAcidSetBuilder::new_standard().build().unwrap();
        for aa in set.iter_variants() {
            assert!(!aa.is_modified());
        }
    }

    #[test]
    fn fixed_mod_replaces_residue() {
        let set = AminoAcidSetBuilder::new_standard()
            .add_fixed_mod(carbamidomethyl_c())
            .build().unwrap();
        let c_variants = set.variants_for(b'C', ModLocation::Anywhere);
        assert_eq!(c_variants.len(), 1);
        assert!(c_variants[0].is_modified());
    }

    #[test]
    fn variable_mod_adds_residue_variant() {
        let set = AminoAcidSetBuilder::new_standard()
            .add_variable_mod(oxidation_m())
            .build().unwrap();
        let m_variants = set.variants_for(b'M', ModLocation::Anywhere);
        assert_eq!(m_variants.len(), 2);
        assert!(m_variants.iter().any(|aa| !aa.is_modified()));
        assert!(m_variants.iter().any(|aa| aa.is_modified()));
    }

    #[test]
    fn conflicting_fixed_and_variable_errors() {
        let cam_fixed = carbamidomethyl_c();
        let mut cam_variable = carbamidomethyl_c();
        cam_variable.fixed = false;

        let err = AminoAcidSetBuilder::new_standard()
            .add_fixed_mod(cam_fixed)
            .add_variable_mod(cam_variable)
            .build()
            .unwrap_err();
        assert!(matches!(err, AaSetError::ConflictingMods { residue: 'C', location: ModLocation::Anywhere }));
    }

    #[test]
    fn implausible_mass_errors() {
        let bad = Modification {
            name: "Bad".to_string(),
            mass_delta: 1500.0,
            residue: ResidueSpec::Specific(b'C'),
            location: ModLocation::Anywhere,
            fixed: true,
            accession: None,
        };
        let err = AminoAcidSetBuilder::new_standard()
            .add_fixed_mod(bad)
            .build().unwrap_err();
        assert!(matches!(err, AaSetError::ImplausibleMassDelta { .. }));
    }

    #[test]
    fn standard_lookup() {
        let set = AminoAcidSetBuilder::new_standard().build().unwrap();
        let g = set.standard(b'G').unwrap();
        assert_eq!(g.residue, b'G');
        assert!(set.standard(b'!').is_none());
    }

    #[test]
    fn min_max_aa_mass() {
        let set = AminoAcidSetBuilder::new_standard().build().unwrap();
        // Min: G ≈ 57.02, Max: W ≈ 186.08
        let g = AminoAcid::standard(b'G').unwrap().mass;
        let w = AminoAcid::standard(b'W').unwrap().mass;
        assert_eq!(set.min_aa_mass(), g);
        assert_eq!(set.max_aa_mass(), w);
    }

    #[test]
    fn max_residue_mod_mass_includes_mods() {
        let set = AminoAcidSetBuilder::new_standard()
            .add_variable_mod(oxidation_m())
            .build().unwrap();
        let w = AminoAcid::standard(b'W').unwrap().mass;
        let expected = w + 15.99491;
        assert!((set.max_residue_mod_mass() - expected).abs() < 1e-9);
    }

    #[test]
    fn contains_cterm_mods_default_false() {
        let set = AminoAcidSetBuilder::new_standard().build().unwrap();
        assert!(!set.contains_cterm_mods());
    }

    #[test]
    fn contains_cterm_mods_when_added() {
        let cterm_mod = Modification {
            name: "Amide".to_string(),
            mass_delta: -0.984016,
            residue: ResidueSpec::Wildcard,
            location: ModLocation::CTerm,
            fixed: false,
            accession: None,
        };
        let set = AminoAcidSetBuilder::new_standard()
            .add_variable_mod(cterm_mod)
            .build().unwrap();
        assert!(set.contains_cterm_mods());
    }

    #[test]
    fn standard_with_carbamidomethyl_c_convenience() {
        let set = AminoAcidSetBuilder::new_standard_with_carbamidomethyl_c().build().unwrap();
        let c_variants = set.variants_for(b'C', ModLocation::Anywhere);
        assert_eq!(c_variants.len(), 1);
        assert!(c_variants[0].is_modified());
    }

    #[test]
    fn add_mods_from_file_parses_real_format() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(tmp.path(),
            "# comment line\n\
             \n\
             57.021464,C,fix,any,Carbamidomethyl\n\
             15.994915,M,opt,any,Oxidation\n").unwrap();

        let set = AminoAcidSetBuilder::new_standard()
            .add_mods_from_file(tmp.path()).unwrap()
            .build().unwrap();

        assert_eq!(set.variants_for(b'C', ModLocation::Anywhere).len(), 1);
        assert!(set.variants_for(b'C', ModLocation::Anywhere)[0].is_modified());
        assert_eq!(set.variants_for(b'M', ModLocation::Anywhere).len(), 2);
    }

    #[test]
    fn add_mods_from_file_reports_line_number() {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        std::fs::write(tmp.path(),
            "57.021464,C,fix,any,Carbamidomethyl\n\
             garbage_line\n").unwrap();

        let err = AminoAcidSetBuilder::new_standard()
            .add_mods_from_file(tmp.path()).unwrap_err();
        match err {
            AaSetError::ModsTxtParse { line_no, .. } => assert_eq!(line_no, 2),
            other => panic!("expected ModsTxtParse, got {:?}", other),
        }
    }
}
