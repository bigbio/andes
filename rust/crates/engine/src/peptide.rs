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
}
