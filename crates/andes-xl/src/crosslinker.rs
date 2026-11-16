//! Crosslinker chemistry: intact mass, MS-cleavage stubs, reactive sites and
//! monolinks. A [`Crosslinker`] comes from a spec (`--xl "DSSO (K, Protein
//! N-term)"`, see [`crate::spec`]) resolved against the reagent table
//! ([`crate::reagents`]) and UNIMOD.

/// Monoisotopic mass of the 13C-12C isotope spacing.
pub const ISOTOPE_SPACING: f64 = 1.003_354_835;
/// Proton mass.
pub const PROTON: f64 = 1.007_276_49;

/// Most stub forms a cleavable linker can have.
pub const MAX_STUBS: usize = model::xl_reagents::MAX_STUBS;
/// Number of [`LinkForm`] slots: the stubs, then the intact form.
pub const N_FORMS: usize = MAX_STUBS + 1;

/// A mass the cleaved linker leaves on a peptide in MS2.
#[derive(Debug, Clone, PartialEq)]
pub struct Stub {
    /// Short name used in output columns (`alkene`, `thiol`, `Bu`, ...).
    pub name: String,
    /// Mass added to the link residue.
    pub mass: f64,
}

/// How an MS-cleavable linker breaks in MS2.
#[derive(Debug, Clone, PartialEq)]
pub struct Cleavage {
    /// The stub masses a released peptide can carry, at most [`MAX_STUBS`].
    pub stubs: Vec<Stub>,
    /// Reporter doublets as `(light, heavy)` indices into `stubs`, the light
    /// stub lighter than the heavy one: one peptide released with either stub
    /// gives two ions `(heavy - light) / z` apart.
    pub doublets: Vec<(usize, usize)>,
}

impl Cleavage {
    /// The reporter-doublet spacing of doublet `k` (neutral mass).
    pub fn doublet_delta(&self, k: usize) -> f64 {
        let (l, h) = self.doublets[k];
        self.stubs[h].mass - self.stubs[l].mass
    }

    /// The light stub mass of doublet `k`.
    pub fn doublet_light(&self, k: usize) -> f64 {
        self.stubs[self.doublets[k].0].mass
    }

    /// The stub masses, in order.
    pub fn stub_masses(&self) -> impl Iterator<Item = f64> + '_ {
        self.stubs.iter().map(|s| s.mass)
    }

    /// A copy whose heavy doublet stubs are moved by `shift` Da (a null
    /// spacing for the doublet probe).
    pub fn with_heavy_shift(&self, shift: f64) -> Cleavage {
        let mut c = self.clone();
        let mut heavy: Vec<usize> = c.doublets.iter().map(|d| d.1).collect();
        heavy.sort_unstable();
        heavy.dedup();
        for h in heavy {
            c.stubs[h].mass += shift;
        }
        c
    }
}

/// The residues and termini one end of a linker reacts with.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct LinkEnd {
    /// Reactive residues (one-letter codes), in written order.
    pub residues: Vec<u8>,
    /// The protein N-terminal amine.
    pub protein_n_term: bool,
    /// Every peptide N-terminal amine: for synthetic peptide libraries, or a
    /// sample linked after digestion.
    pub peptide_n_term: bool,
    /// The protein C-terminal carboxyl.
    pub protein_c_term: bool,
}

impl LinkEnd {
    /// Whether this end reacts with residue `r`.
    pub fn reacts_with(&self, r: u8) -> bool {
        self.residues.contains(&r)
    }

    /// The end as a site list, e.g. `K, S, Protein N-term`.
    pub fn describe(&self) -> String {
        let mut parts: Vec<String> = self
            .residues
            .iter()
            .map(|&r| (r as char).to_string())
            .collect();
        if self.peptide_n_term {
            parts.push("N-term".into());
        }
        if self.protein_n_term {
            parts.push("Protein N-term".into());
        }
        if self.protein_c_term {
            parts.push("Protein C-term".into());
        }
        parts.join(", ")
    }
}

/// A quenched linker hanging from one residue (a dead end).
#[derive(Debug, Clone, PartialEq)]
pub struct Monolink {
    /// UNIMOD name, or a descriptive name for a reagent UNIMOD lacks.
    pub name: String,
    /// UNIMOD accession, when UNIMOD has the record.
    pub accession: Option<String>,
    /// Mass added to the residue.
    pub mass: f64,
}

/// A crosslinking reagent.
#[derive(Debug, Clone, PartialEq)]
pub struct Crosslinker {
    /// Reagent name (`DSSO`, `DSS`, or a custom label).
    pub name: String,
    /// UNIMOD accession of the intact cross-link, when UNIMOD has it.
    pub accession: Option<String>,
    /// Mass added to the two peptides by the intact linker.
    pub intact: f64,
    /// The MS2 cleavage, for an MS-cleavable linker.
    pub cleavage: Option<Cleavage>,
    /// What each end reacts with; equal for a homobifunctional linker.
    pub ends: [LinkEnd; 2],
    /// Quenched monolinks, searchable as variable modifications.
    pub monolinks: Vec<Monolink>,
}

impl Crosslinker {
    /// DSSO with its default sites (K, protein N-terminus), from the reagent
    /// table and UNIMOD.
    pub fn dsso() -> Self {
        crate::spec::XlDecl::parse("DSSO")
            .and_then(|d| d.resolve(model::unimod::bundled()))
            .expect("DSSO resolves from the bundled UNIMOD table")
    }

    /// Whether both ends react with the same sites.
    pub fn is_homobifunctional(&self) -> bool {
        self.ends[0] == self.ends[1]
    }

    /// Whether either end reacts with residue `r`.
    pub fn reacts_with(&self, r: u8) -> bool {
        self.ends.iter().any(|e| e.reacts_with(r))
    }

    /// Whether either end can sit on a peptide's N-terminal amine (any
    /// peptide, or a protein N-terminal one).
    pub fn n_term(&self) -> bool {
        self.ends
            .iter()
            .any(|e| e.protein_n_term || e.peptide_n_term)
    }

    /// The forms a link residue can carry: each stub, then the intact form.
    pub fn forms(&self) -> impl Iterator<Item = LinkForm> {
        let n = self.cleavage.as_ref().map_or(0, |c| c.stubs.len());
        (0..n as u8)
            .map(LinkForm::Stub)
            .chain(std::iter::once(LinkForm::Intact))
    }

    /// Short name of `form` used in output columns.
    pub fn form_name(&self, form: LinkForm) -> &str {
        match form {
            LinkForm::Intact => "intact",
            LinkForm::Stub(k) => self
                .cleavage
                .as_ref()
                .and_then(|c| c.stubs.get(k as usize))
                .map_or("stub", |s| s.name.as_str()),
        }
    }

    /// The sites as written in a spec: one end, or `end / end`.
    pub fn describe_sites(&self) -> String {
        if self.is_homobifunctional() {
            self.ends[0].describe()
        } else {
            format!("{} / {}", self.ends[0].describe(), self.ends[1].describe())
        }
    }
}

/// What the link residue of a peptide carries in a fragment spectrum.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum LinkForm {
    /// Stub `k` of the linker's [`Cleavage`] after MS2 cleavage.
    Stub(u8),
    /// The uncleaved linker with the partner peptide.
    Intact,
}

impl LinkForm {
    /// Slot of the form in per-form arrays of length [`N_FORMS`]: the stub
    /// index, or [`MAX_STUBS`] for the intact form.
    pub fn index(self) -> usize {
        match self {
            LinkForm::Stub(k) => k as usize,
            LinkForm::Intact => MAX_STUBS,
        }
    }

    /// Mass this form puts on the link residue. `partner_mass` is the neutral
    /// mass of the other peptide (used by [`LinkForm::Intact`] only). `None`
    /// for a stub the linker does not have.
    pub fn delta(self, linker: &Crosslinker, partner_mass: f64) -> Option<f64> {
        match self {
            LinkForm::Intact => Some(partner_mass + linker.intact),
            LinkForm::Stub(k) => linker
                .cleavage
                .as_ref()
                .and_then(|c| c.stubs.get(k as usize))
                .map(|s| s.mass),
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
        let c = x.cleavage.as_ref().unwrap();
        let m: Vec<f64> = c.stub_masses().collect();
        // alkene + sulfenic acid = intact linker; thiol = sulfenic - water.
        assert!((m[0] + m[2] - x.intact).abs() < 1e-5);
        assert!((m[2] - H2O - m[1]).abs() < 1e-5);
        assert_eq!(c.doublets, vec![(0, 1)]);
        assert!((c.doublet_delta(0) - 31.972_07).abs() < 1e-4);
        assert!((c.doublet_light(0) - 54.010_565).abs() < 1e-9);
    }

    #[test]
    fn dsso_keeps_the_constants_of_the_first_search() {
        // The values the DSSO search used before the reagent table.
        let x = Crosslinker::dsso();
        let c = x.cleavage.as_ref().unwrap();
        assert_eq!(x.intact, 158.003_765);
        let m: Vec<f64> = c.stub_masses().collect();
        assert_eq!(m, vec![54.010_565, 85.982_635, 103.993_200]);
        let names: Vec<&str> = c.stubs.iter().map(|s| s.name.as_str()).collect();
        assert_eq!(names, ["alkene", "thiol", "sulfenic"]);
        assert_eq!(x.ends[0].residues, b"K");
        assert!(x.ends[0].protein_n_term && !x.ends[0].protein_c_term);
        assert!(!x.ends[0].peptide_n_term);
        assert!(x.is_homobifunctional());
    }

    #[test]
    fn forms_are_the_stubs_then_intact() {
        let x = Crosslinker::dsso();
        let f: Vec<LinkForm> = x.forms().collect();
        assert_eq!(
            f,
            vec![
                LinkForm::Stub(0),
                LinkForm::Stub(1),
                LinkForm::Stub(2),
                LinkForm::Intact
            ]
        );
        assert_eq!(x.form_name(LinkForm::Stub(1)), "thiol");
        assert_eq!(LinkForm::Intact.index(), MAX_STUBS);
        let mut y = x.clone();
        y.cleavage = None;
        assert_eq!(y.forms().collect::<Vec<_>>(), vec![LinkForm::Intact]);
        assert_eq!(LinkForm::Stub(0).delta(&y, 0.0), None);
    }

    #[test]
    fn intact_form_carries_the_partner() {
        let x = Crosslinker::dsso();
        assert!((LinkForm::Intact.delta(&x, 1000.0).unwrap() - 1_158.003_765).abs() < 1e-9);
        assert!((LinkForm::Stub(0).delta(&x, 1000.0).unwrap() - 54.010_565).abs() < 1e-9);
        assert!((crosslink_mass(&x, 1000.0, 2000.0) - 3_158.003_765).abs() < 1e-9);
    }

    #[test]
    fn heavy_shift_moves_only_the_heavy_doublet_stub() {
        let c = Crosslinker::dsso().cleavage.unwrap();
        let s = c.with_heavy_shift(5.0);
        assert_eq!(s.stubs[0].mass, c.stubs[0].mass);
        assert!((s.stubs[1].mass - c.stubs[1].mass - 5.0).abs() < 1e-12);
        assert_eq!(s.stubs[2].mass, c.stubs[2].mass);
    }
}
