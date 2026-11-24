pub const ISOTOPE_SPACING: f64 = model::mass::ISOTOPE;
pub const PROTON: f64 = model::mass::PROTON;

pub const MAX_STUBS: usize = model::crosslinkers::MAX_STUBS;
pub const N_FORMS: usize = MAX_STUBS + 1;

#[derive(Debug, Clone, PartialEq)]
pub struct Stub {
    pub name: String,
    pub mass: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Cleavage {
    pub stubs: Vec<Stub>,
    pub doublets: Vec<(usize, usize)>,
}

impl Cleavage {
    pub fn doublet_delta(&self, k: usize) -> f64 {
        let (l, h) = self.doublets[k];
        self.stubs[h].mass - self.stubs[l].mass
    }

    pub fn doublet_light(&self, k: usize) -> f64 {
        self.stubs[self.doublets[k].0].mass
    }

    pub fn stub_masses(&self) -> impl Iterator<Item = f64> + '_ {
        self.stubs.iter().map(|s| s.mass)
    }

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

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct LinkEnd {
    pub residues: Vec<u8>,
    pub protein_n_term: bool,
    pub peptide_n_term: bool,
    pub protein_c_term: bool,
}

impl LinkEnd {
    pub fn reacts_with(&self, r: u8) -> bool {
        self.residues.contains(&r)
    }

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

#[derive(Debug, Clone, PartialEq)]
pub struct Monolink {
    pub name: String,
    pub accession: Option<String>,
    pub mass: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Crosslinker {
    pub name: String,
    pub accession: Option<String>,
    pub intact: f64,
    pub cleavage: Option<Cleavage>,
    pub ends: [LinkEnd; 2],
    pub monolinks: Vec<Monolink>,
}

impl Crosslinker {
    pub fn dsso() -> Self {
        crate::spec::XlDecl::parse("DSSO")
            .and_then(|d| d.resolve(model::unimod::bundled()))
            .expect("DSSO resolves from the bundled UNIMOD table")
    }

    pub fn is_homobifunctional(&self) -> bool {
        self.ends[0] == self.ends[1]
    }

    pub fn reacts_with(&self, r: u8) -> bool {
        self.ends.iter().any(|e| e.reacts_with(r))
    }

    pub fn n_term(&self) -> bool {
        self.ends
            .iter()
            .any(|e| e.protein_n_term || e.peptide_n_term)
    }

    pub fn forms(&self) -> impl Iterator<Item = LinkForm> {
        let n = self.cleavage.as_ref().map_or(0, |c| c.stubs.len());
        (0..n as u8)
            .map(LinkForm::Stub)
            .chain(std::iter::once(LinkForm::Intact))
    }

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

    pub fn describe_sites(&self) -> String {
        if self.is_homobifunctional() {
            self.ends[0].describe()
        } else {
            format!("{} / {}", self.ends[0].describe(), self.ends[1].describe())
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum LinkForm {
    Stub(u8),
    Intact,
}

impl LinkForm {
    pub fn index(self) -> usize {
        match self {
            LinkForm::Stub(k) => k as usize,
            LinkForm::Intact => MAX_STUBS,
        }
    }

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
        assert!((m[0] + m[2] - x.intact).abs() < 1e-5);
        assert!((m[2] - H2O - m[1]).abs() < 1e-5);
        assert_eq!(c.doublets, vec![(0, 1)]);
        assert!((c.doublet_delta(0) - 31.972_07).abs() < 1e-4);
        assert!((c.doublet_light(0) - 54.010_565).abs() < 1e-9);
    }

    #[test]
    fn dsso_keeps_the_constants_of_the_first_search() {
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
