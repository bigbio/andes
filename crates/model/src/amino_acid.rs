pub mod set;
use crate::mass::{nominal_from, C, H, N, O, S};
use crate::modification::Modification;
use std::hash::{Hash, Hasher};

#[derive(Debug, Clone)]
pub struct AminoAcid {
    pub residue: u8,
    pub mass: f64,
    pub mod_: Option<&'static Modification>,
}

impl AminoAcid {
    pub fn standard(residue: u8) -> Option<Self> {
        if residue == b'U' {
            return Some(AminoAcid {
                residue: b'U',
                mass: 150.953636,
                mod_: None,
            });
        }
        let (c, h, n, o, s) = standard_composition(residue)?;
        let mass = c as f64 * C + h as f64 * H + n as f64 * N + o as f64 * O + s as f64 * S;
        Some(AminoAcid {
            residue,
            mass,
            mod_: None,
        })
    }

    pub fn with_mod<M: Into<&'static Modification>>(mut self, m: M) -> Self {
        self.mod_ = Some(m.into());
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

fn mods_eq(a: &Option<&'static Modification>, b: &Option<&'static Modification>) -> bool {
    match (a, b) {
        (None, None) => true,
        (Some(x), Some(y)) => {
            if std::ptr::eq(*x, *y) {
                return true;
            }
            x.name == y.name && x.mass_delta.to_bits() == y.mass_delta.to_bits()
        }
        _ => false,
    }
}

pub fn standard_composition(residue: u8) -> Option<(u32, u32, u32, u32, u32)> {
    Some(match residue {
        b'G' => (2, 3, 1, 1, 0),
        b'A' => (3, 5, 1, 1, 0),
        b'S' => (3, 5, 1, 2, 0),
        b'P' => (5, 7, 1, 1, 0),
        b'V' => (5, 9, 1, 1, 0),
        b'T' => (4, 7, 1, 2, 0),
        b'C' => (3, 5, 1, 1, 1),
        b'L' => (6, 11, 1, 1, 0),
        b'I' => (6, 11, 1, 1, 0),
        b'N' => (4, 6, 2, 2, 0),
        b'D' => (4, 5, 1, 3, 0),
        b'Q' => (5, 8, 2, 2, 0),
        b'K' => (6, 12, 2, 1, 0),
        b'E' => (5, 7, 1, 3, 0),
        b'M' => (5, 9, 1, 1, 1),
        b'H' => (6, 7, 3, 1, 0),
        b'F' => (9, 9, 1, 1, 0),
        b'R' => (6, 12, 4, 1, 0),
        b'Y' => (9, 9, 1, 2, 0),
        b'W' => (11, 10, 2, 1, 0),
        _ => return None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::modification::{ModLocation, Modification, ResidueSpec};

    #[test]
    fn standard_g_mass_matches_composition() {
        let g = AminoAcid::standard(b'G').unwrap();
        assert_eq!(g.residue, b'G');
        let expected = 2.0 * crate::mass::C
            + 3.0 * crate::mass::H
            + 1.0 * crate::mass::N
            + 1.0 * crate::mass::O;
        assert_eq!(g.mass.to_bits(), expected.to_bits());
        assert!(g.mod_.is_none());
    }

    #[test]
    fn standard_unknown_residue_is_none() {
        assert!(AminoAcid::standard(b'X').is_none());
        assert!(AminoAcid::standard(b'!').is_none());
    }

    #[test]
    fn selenocysteine_u_is_supported() {
        let u = AminoAcid::standard(b'U').expect("U (selenocysteine) must be recognized");
        assert!(
            (u.mass - 150.953636).abs() < 1e-4,
            "Sec mass {} off",
            u.mass
        );
        assert!(u.mod_.is_none());
        assert!(AminoAcid::standard(b'B').is_none());
    }

    #[test]
    fn nominal_mass_for_glycine() {
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
            neutral_losses: Vec::new(),
            loss_class: 0,
            annotated_site_mass: None,
            blocks_cleavage: false,
            parts: Vec::new(),
            diagnostic_ions: Vec::new(),
        };
        let m = AminoAcid::standard(b'M')
            .unwrap()
            .with_mod(oxidation.clone());
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
            neutral_losses: Vec::new(),
            loss_class: 0,
            annotated_site_mass: None,
            blocks_cleavage: false,
            parts: Vec::new(),
            diagnostic_ions: Vec::new(),
        };
        let m = AminoAcid::standard(b'M').unwrap().with_mod(oxidation);
        assert_eq!(m.nominal_mass(), 147);
    }

    #[test]
    fn eq_compares_by_to_bits() {
        let a = AminoAcid::standard(b'G').unwrap();
        let b = AminoAcid::standard(b'G').unwrap();
        assert_eq!(a, b);

        let mut c = a.clone();
        c.mass = 57.0214637_f64;
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
