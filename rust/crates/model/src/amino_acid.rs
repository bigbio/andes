//! Amino acid residue with optional modification. Standard residue masses
//! are computed from atomic composition (C/H/N/O/S counts) so they are
//! bit-equal to the canonical composition-based mass. Pinned by
//! `tests/standard_aa_masses_match_java.rs`.

use std::hash::{Hash, Hasher};

use crate::mass::{nominal_from, C, H, N, O, S};
use crate::modification::Modification;

#[derive(Debug, Clone)]
pub struct AminoAcid {
    pub residue: u8,
    pub mass:    f64,
    pub mod_:    Option<Modification>,
}

impl AminoAcid {
    /// Look up the standard (unmodified) residue table. Returns `None`
    /// for any byte not in the 20-residue standard set.
    pub fn standard(residue: u8) -> Option<Self> {
        let (c, h, n, o, s) = standard_composition(residue)?;
        let mass = c as f64 * C + h as f64 * H + n as f64 * N
                 + o as f64 * O + s as f64 * S;
        Some(AminoAcid { residue, mass, mod_: None })
    }

    /// Attach a modification, returning the modified residue. The `mass`
    /// field is unchanged; consumers compute total mass as `aa.mass +
    /// mod_.mass_delta` separately (see `Peptide::mass`).
    pub fn with_mod(mut self, m: Modification) -> Self {
        self.mod_ = Some(m);
        self
    }

    pub fn nominal_mass(&self) -> i32 {
        let total = self.mass + self.mod_.as_ref().map_or(0.0, |m| m.mass_delta);
        nominal_from(total)
    }

    pub fn is_modified(&self) -> bool {
        self.mod_.is_some()
    }
}

// Custom Eq/Hash via to_bits() — bit-exact comparison (NOT IEEE 754).
// Needed because AminoAcid contains f64, which doesn't implement Eq/Hash
// directly.
impl PartialEq for AminoAcid {
    fn eq(&self, other: &Self) -> bool {
        self.residue == other.residue
            && self.mass.to_bits() == other.mass.to_bits()
            && mods_eq(&self.mod_, &other.mod_)
    }
}

impl Eq for AminoAcid {}

impl Hash for AminoAcid {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.residue.hash(state);
        self.mass.to_bits().hash(state);
        match &self.mod_ {
            None => 0u8.hash(state),
            Some(m) => {
                1u8.hash(state);
                m.name.hash(state);
                m.mass_delta.to_bits().hash(state);
            }
        }
    }
}

fn mods_eq(a: &Option<Modification>, b: &Option<Modification>) -> bool {
    match (a, b) {
        (None, None) => true,
        (Some(x), Some(y)) => {
            x.name == y.name && x.mass_delta.to_bits() == y.mass_delta.to_bits()
        }
        _ => false,
    }
}

/// 20 standard AA atomic compositions (C, H, N, O, S). Computing mass
/// from these integer counts at runtime guarantees bit-equal parity with
/// a canonical composition-based mass.
fn standard_composition(residue: u8) -> Option<(u32, u32, u32, u32, u32)> {
    Some(match residue {
        b'G' => (2,  3, 1, 1, 0),
        b'A' => (3,  5, 1, 1, 0),
        b'S' => (3,  5, 1, 2, 0),
        b'P' => (5,  7, 1, 1, 0),
        b'V' => (5,  9, 1, 1, 0),
        b'T' => (4,  7, 1, 2, 0),
        b'C' => (3,  5, 1, 1, 1),
        b'L' => (6, 11, 1, 1, 0),
        b'I' => (6, 11, 1, 1, 0),
        b'N' => (4,  6, 2, 2, 0),
        b'D' => (4,  5, 1, 3, 0),
        b'Q' => (5,  8, 2, 2, 0),
        b'K' => (6, 12, 2, 1, 0),
        b'E' => (5,  7, 1, 3, 0),
        b'M' => (5,  9, 1, 1, 1),
        b'H' => (6,  7, 3, 1, 0),
        b'F' => (9,  9, 1, 1, 0),
        b'R' => (6, 12, 4, 1, 0),
        b'Y' => (9,  9, 1, 2, 0),
        b'W' => (11, 10, 2, 1, 0),
        _ => return None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::modification::{Modification, ModLocation, ResidueSpec};

    #[test]
    fn standard_g_mass_matches_composition() {
        let g = AminoAcid::standard(b'G').unwrap();
        assert_eq!(g.residue, b'G');
        // Glycine = C2H3NO = 2*12 + 3*1.007825035 + 1*14.003074 + 1*15.99491463
        let expected = 2.0 * crate::mass::C + 3.0 * crate::mass::H
                     + 1.0 * crate::mass::N + 1.0 * crate::mass::O;
        assert_eq!(g.mass.to_bits(), expected.to_bits());
        assert!(g.mod_.is_none());
    }

    #[test]
    fn standard_unknown_residue_is_none() {
        assert!(AminoAcid::standard(b'X').is_none());
        assert!(AminoAcid::standard(b'!').is_none());
    }

    #[test]
    fn nominal_mass_for_glycine() {
        // Gly mass ≈ 57.02146 → nominal 57
        let g = AminoAcid::standard(b'G').unwrap();
        assert_eq!(g.nominal_mass(), 57);
    }

    #[test]
    fn nominal_mass_for_tryptophan() {
        let w = AminoAcid::standard(b'W').unwrap();
        assert_eq!(w.nominal_mass(), 186);
    }

    #[test]
    fn with_mod_attaches_modification() {
        let oxidation = Modification {
            name: "Oxidation".to_string(),
            mass_delta: 15.99491,
            residue: ResidueSpec::Specific(b'M'),
            location: ModLocation::Anywhere,
            fixed: false,
            accession: None,
        };
        let m = AminoAcid::standard(b'M').unwrap().with_mod(oxidation.clone());
        assert!(m.is_modified());
        assert_eq!(m.mod_.as_ref().unwrap().mass_delta, 15.99491);
    }

    #[test]
    fn nominal_mass_includes_mod_delta() {
        let oxidation = Modification {
            name: "Oxidation".to_string(),
            mass_delta: 15.99491,
            residue: ResidueSpec::Specific(b'M'),
            location: ModLocation::Anywhere,
            fixed: false,
            accession: None,
        };
        let m = AminoAcid::standard(b'M').unwrap().with_mod(oxidation);
        // M (131) + Ox (16) = 147 nominal
        assert_eq!(m.nominal_mass(), 147);
    }

    #[test]
    fn eq_compares_by_to_bits() {
        let a = AminoAcid::standard(b'G').unwrap();
        let b = AminoAcid::standard(b'G').unwrap();
        assert_eq!(a, b);

        // Two AAs with the same residue but different mass are NOT equal.
        let mut c = a.clone();
        c.mass = 57.0214637_f64;  // slightly off
        assert_ne!(a, c);
    }

    #[test]
    fn hash_consistent_with_eq() {
        use std::collections::HashSet;
        let a = AminoAcid::standard(b'G').unwrap();
        let b = AminoAcid::standard(b'G').unwrap();
        let set: HashSet<_> = [a, b].into_iter().collect();
        assert_eq!(set.len(), 1);
    }
}
