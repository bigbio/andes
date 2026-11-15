//! Crosslinker chemistry: intact mass, MS-cleavage stubs and reactive sites.

/// Monoisotopic mass of the 13C-12C isotope spacing.
pub const ISOTOPE_SPACING: f64 = 1.003_354_835;
/// Proton mass.
pub const PROTON: f64 = 1.007_276_49;

/// Stub masses left on a peptide when an MS-cleavable linker breaks in MS2.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cleavage {
    /// Alkene stub (DSSO: C3H2O).
    pub alkene: f64,
    /// Unsaturated thiol stub (DSSO: C3H2OS), the sulfenic acid minus water.
    pub thiol: f64,
    /// Sulfenic acid stub (DSSO: C3H4O2S); alkene + sulfenic = intact linker.
    pub sulfenic: f64,
}

impl Cleavage {
    /// The reporter-doublet spacing: thiol minus alkene stub (neutral mass).
    pub fn doublet_delta(&self) -> f64 {
        self.thiol - self.alkene
    }
}

/// A crosslinking reagent.
#[derive(Debug, Clone, PartialEq)]
pub struct Crosslinker {
    /// Lower-case name used on the command line.
    pub name: &'static str,
    /// Mass added to the two peptides by the intact linker.
    pub intact: f64,
    /// Stub masses when the linker is MS-cleavable.
    pub cleavage: Option<Cleavage>,
    /// Residues the linker reacts with (one-letter codes).
    pub residues: &'static [u8],
    /// Whether the linker also reacts with the protein N-terminal amine.
    pub protein_n_term: bool,
}

impl Crosslinker {
    /// DSSO (disuccinimidyl sulfoxide): C6H6O3S2 on two primary amines.
    pub fn dsso() -> Self {
        Crosslinker {
            name: "dsso",
            intact: 158.003_765,
            cleavage: Some(Cleavage {
                alkene: 54.010_565,
                thiol: 85.982_635,
                sulfenic: 103.993_200,
            }),
            residues: b"K",
            protein_n_term: true,
        }
    }

    /// The linker named `name` (case-insensitive), if andes knows it.
    pub fn by_name(name: &str) -> Option<Self> {
        match name.to_ascii_lowercase().as_str() {
            "dsso" => Some(Self::dsso()),
            _ => None,
        }
    }

    /// Names accepted by [`Self::by_name`].
    pub fn known_names() -> &'static [&'static str] {
        &["dsso"]
    }

    /// Whether the linker reacts with residue `r`.
    pub fn reacts_with(&self, r: u8) -> bool {
        self.residues.contains(&r)
    }
}

/// What the link residue of a peptide carries in a fragment spectrum.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum LinkForm {
    /// Alkene stub after MS2 cleavage.
    Alkene,
    /// Thiol stub after MS2 cleavage.
    Thiol,
    /// Sulfenic acid stub after MS2 cleavage.
    Sulfenic,
    /// The uncleaved linker with the partner peptide.
    Intact,
}

impl LinkForm {
    /// Every form, in a fixed order.
    pub const ALL: [LinkForm; 4] = [
        LinkForm::Alkene,
        LinkForm::Thiol,
        LinkForm::Sulfenic,
        LinkForm::Intact,
    ];

    /// The stub forms (cleaved linker).
    pub const STUBS: [LinkForm; 3] = [LinkForm::Alkene, LinkForm::Thiol, LinkForm::Sulfenic];

    /// Short name used in output columns.
    pub fn name(self) -> &'static str {
        match self {
            LinkForm::Alkene => "alkene",
            LinkForm::Thiol => "thiol",
            LinkForm::Sulfenic => "sulfenic",
            LinkForm::Intact => "intact",
        }
    }

    /// Index into [`Self::ALL`].
    pub fn index(self) -> usize {
        self as usize
    }

    /// Mass this form puts on the link residue. `partner_mass` is the neutral
    /// mass of the other peptide (used by [`LinkForm::Intact`] only). `None`
    /// when the linker is not cleavable and the form is a stub.
    pub fn delta(self, linker: &Crosslinker, partner_mass: f64) -> Option<f64> {
        match self {
            LinkForm::Intact => Some(partner_mass + linker.intact),
            LinkForm::Alkene => linker.cleavage.map(|c| c.alkene),
            LinkForm::Thiol => linker.cleavage.map(|c| c.thiol),
            LinkForm::Sulfenic => linker.cleavage.map(|c| c.sulfenic),
        }
    }
}

/// Neutral mass of a crosslink of peptides `a` and `b` (neutral masses).
pub fn crosslink_mass(linker: &Crosslinker, a: f64, b: f64) -> f64 {
    a + b + linker.intact
}

#[cfg(test)]
mod tests {
    use super::*;

    const H2O: f64 = 18.010_565;

    #[test]
    fn dsso_stubs_add_up_to_the_intact_linker() {
        let x = Crosslinker::dsso();
        let c = x.cleavage.unwrap();
        // alkene + sulfenic acid = intact linker; thiol = sulfenic - water.
        assert!((c.alkene + c.sulfenic - x.intact).abs() < 1e-5);
        assert!((c.sulfenic - H2O - c.thiol).abs() < 1e-5);
        assert!((c.doublet_delta() - 31.972_07).abs() < 1e-4);
    }

    #[test]
    fn linker_lookup_is_case_insensitive() {
        assert_eq!(Crosslinker::by_name("DSSO"), Some(Crosslinker::dsso()));
        assert_eq!(Crosslinker::by_name("bs3"), None);
    }

    #[test]
    fn intact_form_carries_the_partner() {
        let x = Crosslinker::dsso();
        assert!((LinkForm::Intact.delta(&x, 1000.0).unwrap() - 1158.003_765).abs() < 1e-9);
        assert!((LinkForm::Alkene.delta(&x, 1000.0).unwrap() - 54.010_565).abs() < 1e-9);
        assert!((crosslink_mass(&x, 1000.0, 2000.0) - 3158.003_765).abs() < 1e-9);
    }
}
