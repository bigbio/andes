//! Positional-isomer class of a peptidoform: its residues and the multiset of
//! its modification mass deltas, without the sites that carry them.
//!
//! Two peptidoforms of one class are the same peptide with the same
//! modifications, possibly on other residues. The class is read from a
//! [`Peptide`] or from peptide text in any of the notations andes writes or
//! reads:
//!
//! - Percolator / PIN: `K.GAC+57.02146LLPK.I` (flanks, signed mass after the
//!   residue),
//! - ProForma-style accessions: `GAC[UNIMOD:4]LLPK`, `[UNIMOD:1]-PEPTIDEK`,
//! - bracketed masses: `PEPM[+15.99491]K`.
//!
//! Modification names and accessions are turned into masses by a resolver
//! (see [`bundled_unimod_mass`]).

use crate::peptide::Peptide;

/// A modification mass delta in units of 1e-5 Da, the precision of the PIN
/// `Peptide` column.
#[inline]
pub fn mod_mass_units(mass_delta: f64) -> i64 {
    (mass_delta * 100_000.0).round() as i64
}

/// Largest difference between two modification masses that [`IsomerClass::matches`]
/// counts as one modification: a mass written to two decimals still matches
/// its full-precision value.
pub const MOD_MASS_TOL_DA: f64 = 0.005;

/// The mass of a modification name or `UNIMOD:NN` accession in the bundled
/// Unimod table.
pub fn bundled_unimod_mass(token: &str) -> Option<f64> {
    crate::unimod::resolve(token)
        .ok()
        .flatten()
        .map(|m| m.mono_mass)
}

/// Residues and sorted modification mass deltas (in [`mod_mass_units`]) of a
/// peptidoform. Terminal modifications count with the others. Equality is
/// exact: two peptidoforms whose masses come from the same values (one
/// search's output) compare equal when they are positional isomers. Across
/// notations, whose masses are rounded differently, compare with
/// [`Self::matches`].
#[derive(Debug, Clone, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct IsomerClass {
    residues: Vec<u8>,
    mod_units: Vec<i64>,
}

impl IsomerClass {
    /// The class of `peptide`.
    pub fn of_peptide(peptide: &Peptide) -> Self {
        let residues = peptide.residues.iter().map(|aa| aa.residue).collect();
        let mut mod_units: Vec<i64> = peptide
            .residues
            .iter()
            .filter_map(|aa| aa.mod_.map(|m| mod_mass_units(m.mass_delta)))
            .collect();
        mod_units.sort_unstable();
        Self {
            residues,
            mod_units,
        }
    }

    /// The class of peptide text; `resolve` gives the mass of a bracketed
    /// modification that is not a signed number (a name or accession).
    /// Errors: a modification without a mass, an unreadable mass, or no
    /// residues.
    pub fn parse(text: &str, resolve: impl Fn(&str) -> Option<f64>) -> Result<Self, String> {
        let mut c = Self::default();
        c.parse_into(text, resolve)?;
        Ok(c)
    }

    /// [`Self::parse`] into `self`, reusing its buffers.
    pub fn parse_into(
        &mut self,
        text: &str,
        resolve: impl Fn(&str) -> Option<f64>,
    ) -> Result<(), String> {
        self.residues.clear();
        self.mod_units.clear();
        let b = text.as_bytes();
        // Percolator-style flanks: `K.` and `.R`.
        let core = if b.len() >= 4 && b[1] == b'.' && b[b.len() - 2] == b'.' {
            &text[2..text.len() - 2]
        } else {
            text
        };
        let bytes = core.as_bytes();
        let mut i = 0;
        while i < bytes.len() {
            let c = bytes[i];
            if c.is_ascii_uppercase() {
                self.residues.push(c);
                i += 1;
            } else if c == b'[' {
                let end = core[i..]
                    .find(']')
                    .map(|j| i + j)
                    .ok_or_else(|| format!("{text}: unclosed modification bracket"))?;
                let tag = &core[i + 1..end];
                let mass = if tag.starts_with(['+', '-'])
                    || tag.starts_with(|c: char| c.is_ascii_digit())
                {
                    tag.parse::<f64>()
                        .map_err(|_| format!("{text}: unreadable mass [{tag}]"))?
                } else {
                    resolve(tag)
                        .ok_or_else(|| format!("{text}: no mass for modification [{tag}]"))?
                };
                self.mod_units.push(mod_mass_units(mass));
                i = end + 1;
            } else if (c == b'+' || c == b'-')
                && bytes
                    .get(i + 1)
                    .is_some_and(|d| d.is_ascii_digit() || *d == b'.')
            {
                let start = i;
                i += 1;
                while i < bytes.len() && (bytes[i].is_ascii_digit() || bytes[i] == b'.') {
                    i += 1;
                }
                let mass = core[start..i]
                    .parse::<f64>()
                    .map_err(|_| format!("{text}: unreadable mass {}", &core[start..i]))?;
                self.mod_units.push(mod_mass_units(mass));
            } else {
                // ProForma terminal separators (`]-`), lowercase terminal
                // markers and other punctuation carry no residue or mass.
                i += 1;
            }
        }
        if self.residues.is_empty() {
            return Err(format!("{text}: no residues"));
        }
        self.mod_units.sort_unstable();
        Ok(())
    }

    /// The residues, in order.
    pub fn residues(&self) -> &[u8] {
        &self.residues
    }

    /// The modification mass deltas in [`mod_mass_units`], sorted.
    pub fn mod_units(&self) -> &[i64] {
        &self.mod_units
    }

    /// Whether `self` and `other` are one class up to the rounding of their
    /// masses: the same residues and, pairing the sorted masses, every pair
    /// within [`MOD_MASS_TOL_DA`].
    pub fn matches(&self, other: &Self) -> bool {
        self.residues == other.residues && mod_units_match(&self.mod_units, &other.mod_units)
    }
}

/// Whether two sorted modification mass multisets (in [`mod_mass_units`]) are
/// equal within [`MOD_MASS_TOL_DA`]: same size, and the sorted masses pair up
/// within the tolerance.
pub fn mod_units_match(a: &[i64], b: &[i64]) -> bool {
    let tol = mod_mass_units(MOD_MASS_TOL_DA);
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| (x - y).abs() <= tol)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{AminoAcid, ModLocation, Modification, ResidueSpec};

    fn resolve(tag: &str) -> Option<f64> {
        match tag {
            "UNIMOD:4" => Some(57.021464),
            "UNIMOD:35" => Some(15.994915),
            "UNIMOD:21" => Some(79.966331),
            "UNIMOD:1" => Some(42.010565),
            _ => None,
        }
    }

    fn class(text: &str) -> IsomerClass {
        IsomerClass::parse(text, resolve).unwrap()
    }

    #[test]
    fn notations_give_one_class() {
        let pin = class("K.GAC+57.02146LLPK.I");
        let proforma = class("GAC[UNIMOD:4]LLPK");
        let bracketed = class("GAC[+57.02146]LLPK");
        assert_eq!(pin.residues(), b"GACLLPK");
        assert_eq!(pin.mod_units(), &[5_702_146]);
        assert!(pin.matches(&proforma) && proforma.matches(&bracketed));
        assert_eq!(pin, bracketed, "same masses, same rounding");
        let flanked = class("-.PEPM[+15.99491]K.-");
        assert_eq!(flanked.residues(), b"PEPMK");
        assert!(flanked.matches(&class("PEPM[UNIMOD:35]K")));
        assert!(flanked.matches(&class("PEPM+15.99K")), "two decimals");
        assert!(!flanked.matches(&class("PEPM+15.98K")));
    }

    #[test]
    fn terminal_modifications_count() {
        let n = class("[UNIMOD:1]-PEPTIDEK");
        assert_eq!(n.residues(), b"PEPTIDEK");
        assert!(n.matches(&class("K.P+42.01057EPTIDEK.A")));
        assert!(!n.matches(&class("PEPTIDEK")));
        assert!(class("PEPTIDEK-[+0.98402]").matches(&class("PEPTIDEK+0.98402")));
        assert!(
            class("n[42.0106]PEPTIDEK").matches(&n),
            "unsigned bracketed mass"
        );
    }

    #[test]
    fn site_is_ignored_and_composition_is_not() {
        let s = class("PEPS[UNIMOD:21]TK");
        let t = class("K.PEPST+79.96633K.A");
        assert!(s.matches(&t), "positional isomers");
        assert!(!s.matches(&class("PEPSTK")), "unmodified");
        assert!(
            !s.matches(&class("PEPS[UNIMOD:21]T[UNIMOD:21]K")),
            "two sites"
        );
        assert!(
            !s.matches(&class("PEPST[+79.95682]K")),
            "sulfation is 0.0095 Da from phosphorylation"
        );
        assert!(!s.matches(&class("PEPTSK[UNIMOD:21]")), "other residues");
        assert!(class("M[UNIMOD:35]PEM[+15.99492]K").matches(&class("M[+15.99492]PEM[UNIMOD:35]K")));
    }

    #[test]
    fn unknown_and_unreadable_modifications_fail() {
        let e = IsomerClass::parse("PEPC[UNIMOD:999]K", resolve).unwrap_err();
        assert!(e.contains("UNIMOD:999"), "{e}");
        assert!(IsomerClass::parse("PEPC[+5.7.0]K", resolve).is_err());
        assert!(IsomerClass::parse("PEPC[UNIMOD:4K", resolve).is_err());
        assert!(IsomerClass::parse("K.+57.0.A", resolve).is_err());
    }

    #[test]
    fn bundled_unimod_names_and_accessions_resolve() {
        let ox = bundled_unimod_mass("UNIMOD:35").unwrap();
        assert!((ox - 15.994915).abs() < 1e-6);
        assert_eq!(bundled_unimod_mass("Oxidation"), Some(ox));
        let a = IsomerClass::parse("PEPM[Oxidation]K", bundled_unimod_mass).unwrap();
        assert!(a.matches(&class("PEPM+15.99491K")));
    }

    #[test]
    fn a_peptide_and_its_text_agree() {
        let mut m = AminoAcid::standard(b'S').unwrap();
        m.mod_ = Some(crate::modification::leak_mod(Modification {
            name: "Phospho".into(),
            mass_delta: 79.966331,
            residue: ResidueSpec::Specific(b'S'),
            location: ModLocation::Anywhere,
            fixed: false,
            accession: Some("UNIMOD:21".into()),
            neutral_losses: Vec::new(),
            loss_class: 0,
            annotated_site_mass: None,
            blocks_cleavage: false,
            diagnostic_ions: Vec::new(),
        }));
        let p = Peptide::new(
            vec![
                AminoAcid::standard(b'P').unwrap(),
                m,
                AminoAcid::standard(b'K').unwrap(),
            ],
            b'R',
            b'A',
        );
        let c = IsomerClass::of_peptide(&p);
        assert_eq!(c, class(&p.to_string()), "the PIN text");
        assert_eq!(c, class("PS[UNIMOD:21]K"));
    }
}
