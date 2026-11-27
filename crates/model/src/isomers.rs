use crate::peptide::Peptide;

#[inline]
pub fn mod_mass_units(mass_delta: f64) -> i64 {
    (mass_delta * 100_000.0).round() as i64
}

pub const MOD_MASS_TOL_DA: f64 = 0.005;

pub fn bundled_unimod_mass(token: &str) -> Option<f64> {
    crate::unimod::resolve(token)
        .ok()
        .flatten()
        .map(|m| m.mono_mass)
}

#[derive(Debug, Clone, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct IsomerClass {
    residues: Vec<u8>,
    mod_units: Vec<i64>,
}

impl IsomerClass {
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

    pub fn parse(text: &str, resolve: impl Fn(&str) -> Option<f64>) -> Result<Self, String> {
        let mut c = Self::default();
        c.parse_into(text, resolve)?;
        Ok(c)
    }

    /// Parses a peptidoform (`K.PEPT+79.97IDE.R`, `[UNIMOD:1]-PEPM[UNIMOD:35]K`). Modification
    /// masses are summed per residue, a terminal modification counting on the first or last
    /// residue, so a terminal and a first-residue modification written apart give the same class
    /// as the two written as one summed mass, as `of_peptide` gives it.
    pub fn parse_into(
        &mut self,
        text: &str,
        resolve: impl Fn(&str) -> Option<f64>,
    ) -> Result<(), String> {
        self.residues.clear();
        self.mod_units.clear();
        let b = text.as_bytes();
        let core = if b.len() >= 4 && b[1] == b'.' && b[b.len() - 2] == b'.' {
            &text[2..text.len() - 2]
        } else {
            text
        };
        let bytes = core.as_bytes();
        let mut pending: Option<f64> = None;
        let mut on_residue: Option<f64> = None;
        let add = |mass: f64, residues: &[u8], pending: &mut Option<f64>, on: &mut Option<f64>| {
            let slot = if residues.is_empty() { pending } else { on };
            *slot = Some(slot.unwrap_or(0.0) + mass);
        };
        let mut i = 0;
        while i < bytes.len() {
            let c = bytes[i];
            if c.is_ascii_uppercase() {
                if let Some(m) = on_residue.take() {
                    self.mod_units.push(mod_mass_units(m));
                }
                on_residue = pending.take();
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
                add(mass, &self.residues, &mut pending, &mut on_residue);
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
                add(mass, &self.residues, &mut pending, &mut on_residue);
            } else {
                i += 1;
            }
        }
        if self.residues.is_empty() {
            return Err(format!("{text}: no residues"));
        }
        if let Some(m) = on_residue {
            self.mod_units.push(mod_mass_units(m));
        }
        self.mod_units.sort_unstable();
        Ok(())
    }

    pub fn residues(&self) -> &[u8] {
        &self.residues
    }

    pub fn mod_units(&self) -> &[i64] {
        &self.mod_units
    }

    pub fn matches(&self, other: &Self) -> bool {
        self.residues == other.residues && mod_units_match(&self.mod_units, &other.mod_units)
    }
}

pub fn mod_units_match(a: &[i64], b: &[i64]) -> bool {
    let tol = mod_mass_units(MOD_MASS_TOL_DA);
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| (x - y).abs() <= tol)
}

/// The residue with I read as L. I and L have the same residue mass, so two peptides that differ
/// only by I and L ("I/L twins") have the same precursor and fragment masses and no MS2 spectrum
/// tells them apart. Every key of peptide identity reads residues through this function, so the
/// twins are one peptide whose proteins are the union of both; the sequence a row reports stays
/// one of the database sequences.
#[inline]
pub fn fold_isoleucine(residue: u8) -> u8 {
    if residue == b'I' {
        b'L'
    } else {
        residue
    }
}

/// `residues` with every I read as L (see [`fold_isoleucine`]).
pub fn fold_isoleucine_bytes(residues: &[u8]) -> Vec<u8> {
    residues.iter().map(|&r| fold_isoleucine(r)).collect()
}

/// A peptide text (a bare sequence, a PIN `Peptide` cell or a ProForma string) with every residue
/// I read as L; text inside `[...]` or `(...)`, such as a modification name or accession, is kept
/// as written. The identity key of a peptide or peptidoform text (see [`fold_isoleucine`]).
pub fn fold_isoleucine_text(peptide: &str) -> String {
    let mut depth = 0u32;
    peptide
        .chars()
        .map(|c| match c {
            '[' | '(' => {
                depth += 1;
                c
            }
            ']' | ')' => {
                depth = depth.saturating_sub(1);
                c
            }
            'I' if depth == 0 => 'L',
            _ => c,
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn modifications_are_summed_per_residue_with_terminal_ones_on_the_end_residue() {
        let r = |t: &str| match t {
            "UNIMOD:5" => Some(43.005814),
            "UNIMOD:35" => Some(15.994915),
            "UNIMOD:4" => Some(57.021464),
            _ => None,
        };
        let pin = IsomerClass::parse("K.M+59.00073QWNVPR.T", r).unwrap();
        let table = IsomerClass::parse("[UNIMOD:5]-M[UNIMOD:35]QWNVPR", r).unwrap();
        assert!(pin.matches(&table));
        let cam = IsomerClass::parse("K.C+100.02728AC+57.02146TVVNR.D", r).unwrap();
        let cam_table = IsomerClass::parse("[UNIMOD:5]-C[UNIMOD:4]AC[UNIMOD:4]TVVNR", r).unwrap();
        assert!(cam.matches(&cam_table));
        let carboxymethyl = IsomerClass::parse("AC[+58.005479]MK", r).unwrap();
        let acetyl_oxidation = IsomerClass::parse("[+42.010565]-ACM[+15.994915]K", r).unwrap();
        assert!(!carboxymethyl.matches(&acetyl_oxidation));
        let c_term = IsomerClass::parse("PEPTIDEK-[+0.984016]", r).unwrap();
        assert_eq!(c_term.mod_units(), &[mod_mass_units(0.984016)]);
    }
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
            parts: Vec::new(),
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

    #[test]
    fn i_l_twins_fold_to_one_key_and_modification_text_is_kept() {
        assert_eq!(fold_isoleucine(b'I'), b'L');
        assert_eq!(fold_isoleucine(b'K'), b'K');
        assert_eq!(fold_isoleucine_bytes(b"LDIK"), b"LDLK".to_vec());
        assert_eq!(fold_isoleucine_text("K.PEPTIDEK.R"), "K.PEPTLDEK.R");
        assert_eq!(
            fold_isoleucine_text("PEPTIDEM[UNIMOD:35]K(Ile->Leu)"),
            "PEPTLDEM[UNIMOD:35]K(Ile->Leu)"
        );
        assert_eq!(fold_isoleucine_text("ELIVK"), fold_isoleucine_text("EILVK"));
    }
}
