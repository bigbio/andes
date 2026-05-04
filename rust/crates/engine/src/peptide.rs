//! Peptide. Mirrors Java `edu.ucsd.msjava.msutil.Peptide`. The
//! `Display` impl (Task 9) is byte-parity-gated by
//! `tests/peptide_display_parity.rs` (Task 12).

use std::hash::{Hash, Hasher};

use crate::amino_acid::AminoAcid;
use crate::mass::{nominal_from, H2O};

#[derive(Debug, Clone)]
pub struct Peptide {
    pub residues: Vec<AminoAcid>,
    /// Flanking residue at the N-terminus (the AA *before* this peptide
    /// in its source protein). `_` for protein N-term, `-` for protein
    /// C-term. Matches Java's `Constants.PROTEIN_N_TERM` / `_C_TERM`.
    pub pre:  u8,
    pub post: u8,
    pub charge: Option<u8>,
}

impl Peptide {
    pub fn new(residues: Vec<AminoAcid>, pre: u8, post: u8) -> Self {
        Self { residues, pre, post, charge: None }
    }

    pub fn with_charge(mut self, charge: u8) -> Self {
        self.charge = Some(charge);
        self
    }

    pub fn length(&self) -> usize {
        self.residues.len()
    }

    /// Total monoisotopic mass: sum of residue masses + sum of mod deltas
    /// + `H2O`. Matches Java's `Peptide.getMass()`.
    pub fn mass(&self) -> f64 {
        let residue_sum: f64 = self.residues.iter().map(|aa| aa.mass).sum();
        let mod_sum: f64 = self.residues
            .iter()
            .filter_map(|aa| aa.mod_.as_ref().map(|m| m.mass_delta))
            .sum();
        residue_sum + mod_sum + H2O
    }

    pub fn nominal_mass(&self) -> i32 {
        nominal_from(self.mass())
    }
}

// Custom Eq/Hash: relies on AminoAcid's custom impls (which route f64
// through to_bits). Same rationale as AminoAcid: f64 doesn't impl Eq/Hash.
impl PartialEq for Peptide {
    fn eq(&self, other: &Self) -> bool {
        self.pre == other.pre
            && self.post == other.post
            && self.charge == other.charge
            && self.residues == other.residues
    }
}

impl Eq for Peptide {}

impl Hash for Peptide {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.pre.hash(state);
        self.post.hash(state);
        self.charge.hash(state);
        self.residues.hash(state);
    }
}

impl std::fmt::Display for Peptide {
    /// Phase 1 canonical text form: `pre.SEQ_WITH_MODS.post`.
    /// Mod deltas render as `{:+.5}` (signed, 5 decimals) after each
    /// modified residue. Charge is not rendered. This format is the
    /// inverse of `Peptide::from_str` (Task 10) and is NOT a byte-parity
    /// match against Java's `Peptide.toString()` — the PIN/TSV output
    /// formats live in the Phase 7 output crate.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}.", self.pre as char)?;
        for aa in &self.residues {
            write!(f, "{}", aa.residue as char)?;
            if let Some(m) = &aa.mod_ {
                write!(f, "{:+.5}", m.mass_delta)?;
            }
        }
        write!(f, ".{}", self.post as char)
    }
}

use crate::aa_set::AminoAcidSet;

#[derive(thiserror::Error, Debug)]
pub enum PeptideParseError {
    #[error("empty peptide string")]
    Empty,
    #[error("malformed flanking residue pattern: expected `X.SEQ.Y`, got {got:?}")]
    BadFlanking { got: String },
    #[error("unknown residue {residue:?} at position {position}")]
    UnknownResidue { residue: char, position: usize },
    #[error("malformed mod-mass token {token:?} at position {position}: {source}")]
    BadModMass { token: String, position: usize, #[source] source: std::num::ParseFloatError },
    #[error("mod {token:?} at position {position} does not match any variant in AminoAcidSet")]
    UnknownMod { token: String, position: usize },
}

impl Peptide {
    /// Parse `pre.SEQ.post` form. `aa_set` provides the variant lookup
    /// for modified residues (matches mass deltas to known
    /// `(residue, mass_delta)` pairs).
    pub fn from_str(s: &str, aa_set: &AminoAcidSet) -> Result<Self, PeptideParseError> {
        if s.is_empty() {
            return Err(PeptideParseError::Empty);
        }
        let bytes = s.as_bytes();
        let first_dot = bytes.iter().position(|&b| b == b'.')
            .ok_or_else(|| PeptideParseError::BadFlanking { got: s.to_string() })?;
        let last_dot = bytes.iter().rposition(|&b| b == b'.')
            .ok_or_else(|| PeptideParseError::BadFlanking { got: s.to_string() })?;
        if first_dot == last_dot || first_dot != 1 || last_dot != bytes.len() - 2 {
            return Err(PeptideParseError::BadFlanking { got: s.to_string() });
        }
        let pre = bytes[0];
        let post = bytes[bytes.len() - 1];
        let middle = &s[first_dot + 1..last_dot];

        let residues = parse_middle(middle, aa_set)?;
        Ok(Peptide::new(residues, pre, post))
    }
}

fn parse_middle(s: &str, aa_set: &AminoAcidSet) -> Result<Vec<AminoAcid>, PeptideParseError> {
    let bytes = s.as_bytes();
    let mut out = Vec::with_capacity(bytes.len());
    let mut i = 0;
    while i < bytes.len() {
        let r = bytes[i];
        if !r.is_ascii_uppercase() {
            return Err(PeptideParseError::UnknownResidue { residue: r as char, position: i });
        }
        i += 1;

        if i < bytes.len() && (bytes[i] == b'+' || bytes[i] == b'-') {
            let start = i;
            i += 1;
            while i < bytes.len() && (bytes[i].is_ascii_digit() || bytes[i] == b'.') {
                i += 1;
            }
            let token = &s[start..i];
            let delta: f64 = token.parse().map_err(|source| {
                PeptideParseError::BadModMass { token: token.to_string(), position: start, source }
            })?;

            let variant = aa_set
                .variants_for(r, crate::modification::ModLocation::Anywhere)
                .iter()
                .find(|aa| aa.mod_.as_ref()
                    .map(|m| m.mass_delta.to_bits() == delta.to_bits())
                    .unwrap_or(false))
                .cloned()
                .ok_or_else(|| PeptideParseError::UnknownMod {
                    token: format!("{}{}", r as char, token), position: start - 1
                })?;
            out.push(variant);
        } else {
            let aa = AminoAcid::standard(r)
                .ok_or_else(|| PeptideParseError::UnknownResidue {
                    residue: r as char, position: i - 1
                })?;
            out.push(aa);
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::amino_acid::AminoAcid;
    use crate::mass::H2O;
    use crate::modification::{Modification, ModLocation, ResidueSpec};

    fn unmod_pep(seq: &[u8]) -> Peptide {
        let residues: Vec<_> = seq.iter().map(|&r| AminoAcid::standard(r).unwrap()).collect();
        Peptide::new(residues, b'_', b'-')
    }

    #[test]
    fn length_counts_residues() {
        let p = unmod_pep(b"PEPTIDE");
        assert_eq!(p.length(), 7);
    }

    #[test]
    fn mass_is_sum_plus_h2o() {
        let p = unmod_pep(b"GA");  // G + A masses
        let g = AminoAcid::standard(b'G').unwrap().mass;
        let a = AminoAcid::standard(b'A').unwrap().mass;
        let expected = g + a + H2O;
        assert_eq!(p.mass().to_bits(), expected.to_bits());
    }

    #[test]
    fn mass_includes_mod_deltas() {
        let oxidation = Modification {
            name: "Oxidation".to_string(),
            mass_delta: 15.99491,
            residue: ResidueSpec::Specific(b'M'),
            location: ModLocation::Anywhere,
            fixed: false,
            accession: None,
        };
        let m = AminoAcid::standard(b'M').unwrap().with_mod(oxidation);
        let g = AminoAcid::standard(b'G').unwrap();
        let m_mass = AminoAcid::standard(b'M').unwrap().mass;
        let p = Peptide::new(vec![m, g.clone()], b'_', b'-');
        let expected = m_mass + 15.99491 + g.mass + H2O;
        assert_eq!(p.mass().to_bits(), expected.to_bits());
    }

    #[test]
    fn nominal_mass_for_GA() {
        let p = unmod_pep(b"GA");
        // G + A + H2O ≈ 146.069 → nominal 146
        assert_eq!(p.nominal_mass(), 146);
    }

    #[test]
    fn with_charge_attaches_charge() {
        let p = unmod_pep(b"PEPTIDE").with_charge(2);
        assert_eq!(p.charge, Some(2));
    }

    #[test]
    fn flanking_bytes_preserved() {
        let p = unmod_pep(b"PEPTIDE");
        assert_eq!(p.pre, b'_');
        assert_eq!(p.post, b'-');
    }

    #[test]
    fn eq_compares_structurally() {
        let p1 = unmod_pep(b"PEPTIDE");
        let p2 = unmod_pep(b"PEPTIDE");
        assert_eq!(p1, p2);

        let p3 = unmod_pep(b"PEPTIDQ");
        assert_ne!(p1, p3);
    }

    #[test]
    fn hash_consistent_with_eq() {
        use std::collections::HashSet;
        let p1 = unmod_pep(b"PEPTIDE");
        let p2 = unmod_pep(b"PEPTIDE");
        let set: HashSet<_> = [p1, p2].into_iter().collect();
        assert_eq!(set.len(), 1);
    }

    fn modded(residue: u8, mod_name: &str, delta: f64) -> AminoAcid {
        let aa = AminoAcid::standard(residue).unwrap();
        let m = Modification {
            name: mod_name.to_string(),
            mass_delta: delta,
            residue: ResidueSpec::Specific(residue),
            location: ModLocation::Anywhere,
            fixed: false,
            accession: None,
        };
        aa.with_mod(m)
    }

    #[test]
    fn display_unmodified() {
        let p = unmod_pep(b"PEPTIDE");
        assert_eq!(p.to_string(), "_.PEPTIDE.-");
    }

    #[test]
    fn display_real_flanking() {
        let mut p = unmod_pep(b"PEPTIDE");
        p.pre = b'K';
        p.post = b'R';
        assert_eq!(p.to_string(), "K.PEPTIDE.R");
    }

    #[test]
    fn display_single_mod() {
        let residues = vec![
            AminoAcid::standard(b'P').unwrap(),
            AminoAcid::standard(b'E').unwrap(),
            modded(b'C', "Carbamidomethyl", 57.02146),
            AminoAcid::standard(b'I').unwrap(),
            AminoAcid::standard(b'D').unwrap(),
            AminoAcid::standard(b'E').unwrap(),
        ];
        let p = Peptide::new(residues, b'_', b'-');
        assert_eq!(p.to_string(), "_.PEC+57.02146IDE.-");
    }

    #[test]
    fn display_oxidation_m() {
        let residues = vec![
            AminoAcid::standard(b'M').unwrap(),
            AminoAcid::standard(b'E').unwrap(),
            modded(b'M', "Oxidation", 15.99491),
            AminoAcid::standard(b'D').unwrap(),
            AminoAcid::standard(b'E').unwrap(),
        ];
        let p = Peptide::new(residues, b'_', b'-');
        assert_eq!(p.to_string(), "_.MEM+15.99491DE.-");
    }

    #[test]
    fn display_negative_mass_mod() {
        let residues = vec![
            modded(b'K', "Pyro-glu", -17.02655),
            AminoAcid::standard(b'R').unwrap(),
            AminoAcid::standard(b'I').unwrap(),
            AminoAcid::standard(b'P').unwrap(),
            modded(b'M', "Oxidation", 15.99491),
        ];
        let p = Peptide::new(residues, b'_', b'-');
        assert_eq!(p.to_string(), "_.K-17.02655RIPM+15.99491.-");
    }

    #[test]
    fn display_multi_mod() {
        let residues = vec![
            AminoAcid::standard(b'P').unwrap(),
            modded(b'C', "Carbamidomethyl", 57.02146),
            AminoAcid::standard(b'P').unwrap(),
            modded(b'M', "Oxidation", 15.99491),
            AminoAcid::standard(b'D').unwrap(),
            AminoAcid::standard(b'E').unwrap(),
        ];
        let p = Peptide::new(residues, b'_', b'-');
        assert_eq!(p.to_string(), "_.PC+57.02146PM+15.99491DE.-");
    }

    #[test]
    fn display_charge_not_rendered() {
        let p = unmod_pep(b"AG").with_charge(2);
        assert_eq!(p.to_string(), "_.AG.-");
        assert_eq!(p.charge, Some(2));
    }

    use crate::aa_set::AminoAcidSetBuilder;

    fn aa_set_with_carbamidomethyl_and_oxidation() -> crate::aa_set::AminoAcidSet {
        let cam = Modification {
            name: "Carbamidomethyl".to_string(),
            mass_delta: 57.02146,
            residue: ResidueSpec::Specific(b'C'),
            location: ModLocation::Anywhere,
            fixed: true,
            accession: None,
        };
        let ox = Modification {
            name: "Oxidation".to_string(),
            mass_delta: 15.99491,
            residue: ResidueSpec::Specific(b'M'),
            location: ModLocation::Anywhere,
            fixed: false,
            accession: None,
        };
        AminoAcidSetBuilder::new_standard()
            .add_fixed_mod(cam)
            .add_variable_mod(ox)
            .build()
            .unwrap()
    }

    #[test]
    fn from_str_unmodified() {
        let aa_set = aa_set_with_carbamidomethyl_and_oxidation();
        let p = Peptide::from_str("_.PEPTIDE.-", &aa_set).unwrap();
        assert_eq!(p.length(), 7);
        assert_eq!(p.pre, b'_');
        assert_eq!(p.post, b'-');
        assert_eq!(p.residues[0].residue, b'P');
    }

    #[test]
    fn from_str_with_carbamidomethyl() {
        let aa_set = aa_set_with_carbamidomethyl_and_oxidation();
        let p = Peptide::from_str("K.PEC+57.02146IDE.R", &aa_set).unwrap();
        assert_eq!(p.length(), 6);
        assert!(p.residues[2].is_modified());
        assert_eq!(p.residues[2].mod_.as_ref().unwrap().name, "Carbamidomethyl");
    }

    #[test]
    fn from_str_with_oxidation_m() {
        let aa_set = aa_set_with_carbamidomethyl_and_oxidation();
        let p = Peptide::from_str("_.MEM+15.99491DE.-", &aa_set).unwrap();
        assert!(!p.residues[0].is_modified());
        assert!(p.residues[2].is_modified());
    }

    #[test]
    fn from_str_round_trip_unmodified() {
        let aa_set = aa_set_with_carbamidomethyl_and_oxidation();
        let s = "_.PEPTIDE.-";
        let p = Peptide::from_str(s, &aa_set).unwrap();
        assert_eq!(p.to_string(), s);
    }

    #[test]
    fn from_str_round_trip_with_mods() {
        let aa_set = aa_set_with_carbamidomethyl_and_oxidation();
        let s = "K.PEC+57.02146PM+15.99491DE.R";
        let p = Peptide::from_str(s, &aa_set).unwrap();
        assert_eq!(p.to_string(), s);
    }

    #[test]
    fn from_str_empty() {
        let aa_set = aa_set_with_carbamidomethyl_and_oxidation();
        let err = Peptide::from_str("", &aa_set).unwrap_err();
        assert!(matches!(err, PeptideParseError::Empty));
    }

    #[test]
    fn from_str_bad_flanking() {
        let aa_set = aa_set_with_carbamidomethyl_and_oxidation();
        let err = Peptide::from_str("PEPTIDE", &aa_set).unwrap_err();
        assert!(matches!(err, PeptideParseError::BadFlanking { .. }));
    }

    #[test]
    fn from_str_unknown_residue() {
        let aa_set = aa_set_with_carbamidomethyl_and_oxidation();
        let err = Peptide::from_str("_.PEPxIDE.-", &aa_set).unwrap_err();
        assert!(matches!(err, PeptideParseError::UnknownResidue { .. }));
    }

    #[test]
    fn from_str_unknown_mod() {
        let aa_set = aa_set_with_carbamidomethyl_and_oxidation();
        let err = Peptide::from_str("_.PEC+99.99999IDE.-", &aa_set).unwrap_err();
        assert!(matches!(err, PeptideParseError::UnknownMod { .. }));
    }
}
