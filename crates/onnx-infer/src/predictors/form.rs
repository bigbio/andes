//! Peptidoforms parsed from andes PIN peptide strings, with the fragment
//! masses the spectral-match features need.

use crate::rt_native::{self, RtPeptide};

/// Proton mass (Da).
pub const PROTON: f64 = 1.007_276_466_8;
/// Water mass (Da).
pub const H2O: f64 = 18.010_564_684;

/// Monoisotopic residue mass of a one-letter amino-acid code.
pub fn residue_mass(c: char) -> Option<f64> {
    Some(match c {
        'G' => 57.021_463_72,
        'A' => 71.037_113_79,
        'S' => 87.032_028_41,
        'P' => 97.052_763_85,
        'V' => 99.068_413_91,
        'T' => 101.047_678_5,
        'C' => 103.009_184_5,
        'L' | 'I' => 113.084_064,
        'N' => 114.042_927_4,
        'D' => 115.026_943,
        'Q' => 128.058_577_5,
        'K' => 128.094_963,
        'E' => 129.042_593_1,
        'M' => 131.040_484_6,
        'H' => 137.058_911_9,
        'F' => 147.068_413_9,
        'R' => 156.101_111,
        'Y' => 163.063_328_6,
        'W' => 186.079_312_9,
        _ => return None,
    })
}

/// Modifications the predictors recognise in andes peptide strings.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Mod {
    /// Carbamidomethyl (Unimod 4) on C.
    Carbamidomethyl,
    /// Oxidation (Unimod 35) on M.
    Oxidation,
    /// Acetyl (Unimod 1) on the N-terminus.
    Acetyl,
    /// TMT6plex (Unimod 737) on K or the N-terminus.
    Tmt6plex,
    /// Phospho (Unimod 21) on S, T or Y.
    Phospho,
    /// Deamidated (Unimod 7) on N or Q.
    Deamidated,
    /// Gln->pyro-Glu (Unimod 28) on an N-terminal Q.
    GlnPyroGlu,
    /// Glu->pyro-Glu (Unimod 27) on an N-terminal E.
    GluPyroGlu,
}

/// Where a modification may sit.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Placement {
    /// On one of the listed residues.
    Residue,
    /// Peptide N-terminus only (site 0).
    NTerm,
    /// On one of the listed residues, or the peptide N-terminus.
    ResidueOrNTerm,
    /// Peptide N-terminus (site 0) of a peptide whose first residue is one of the
    /// listed residues (alphabase `X^Any_N-term`).
    NTermOf,
}

/// (mod, Unimod accession, residues, placement)
const MODS: &[(Mod, u32, &str, Placement)] = &[
    (Mod::Carbamidomethyl, 4, "C", Placement::Residue),
    (Mod::Oxidation, 35, "M", Placement::Residue),
    (Mod::Acetyl, 1, "", Placement::NTerm),
    (Mod::Tmt6plex, 737, "K", Placement::ResidueOrNTerm),
    (Mod::Phospho, 21, "STY", Placement::Residue),
    (Mod::Deamidated, 7, "NQ", Placement::Residue),
    (Mod::GlnPyroGlu, 28, "Q", Placement::NTermOf),
    (Mod::GluPyroGlu, 27, "E", Placement::NTermOf),
];

impl Mod {
    /// alphabase / peptdeep name at a site (0 = N-terminus).
    pub fn alphabase(self, site: usize, aa: char) -> String {
        match self {
            Mod::Carbamidomethyl => "Carbamidomethyl@C".into(),
            Mod::Oxidation => "Oxidation@M".into(),
            Mod::Acetyl => "Acetyl@Protein_N-term".into(),
            Mod::Tmt6plex if site == 0 => "TMT6plex@Any_N-term".into(),
            Mod::Tmt6plex => "TMT6plex@K".into(),
            Mod::Phospho => format!("Phospho@{aa}"),
            Mod::Deamidated => format!("Deamidated@{aa}"),
            Mod::GlnPyroGlu => "Gln->pyro-Glu@Q^Any_N-term".into(),
            Mod::GluPyroGlu => "Glu->pyro-Glu@E^Any_N-term".into(),
        }
    }

    /// Unimod accession.
    pub fn unimod(self) -> u32 {
        MODS.iter()
            .find(|m| m.0 == self)
            .expect("every Mod is in MODS")
            .1
    }

    /// The mod with Unimod accession `id`, if it sits at `site` of `seq` (0 =
    /// N-terminus, 1..=n residue) by its [`Placement`].
    fn at(id: u32, site: usize, seq: &[u8]) -> Option<Self> {
        let &(m, _, residues, placement) = MODS.iter().find(|m| m.1 == id)?;
        let on_residue =
            site >= 1 && site <= seq.len() && residues.as_bytes().contains(&seq[site - 1]);
        let ok = match placement {
            Placement::Residue => on_residue,
            Placement::NTerm => site == 0,
            Placement::ResidueOrNTerm => site == 0 || on_residue,
            Placement::NTermOf => site == 0 && residues.as_bytes().contains(seq.first()?),
        };
        ok.then_some(m)
    }

    /// Elemental composition delta (C, H, N, O, S, P), Unimod; heavy
    /// isotopes counted as their element.
    pub fn composition(self) -> [i32; 6] {
        match self {
            Mod::Carbamidomethyl => [2, 3, 1, 1, 0, 0],
            Mod::Oxidation => [0, 0, 0, 1, 0, 0],
            Mod::Acetyl => [2, 2, 0, 1, 0, 0],
            Mod::Tmt6plex => [12, 20, 2, 2, 0, 0],
            Mod::Phospho => [0, 1, 0, 3, 0, 1],
            Mod::Deamidated => [0, -1, -1, 1, 0, 0],
            Mod::GlnPyroGlu => [0, -3, -1, 0, 0, 0],
            Mod::GluPyroGlu => [0, -2, 0, -1, 0, 0],
        }
    }
}

/// One parsed peptidoform.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Form {
    /// Unmodified residues.
    pub seq: String,
    /// (site, mod) sorted — site 0 = N-term, 1..=n residue.
    pub mods: Vec<(usize, Mod)>,
}

impl Form {
    /// The residue at `site` (the first residue for site 0).
    pub fn residue(&self, site: usize) -> char {
        self.seq.as_bytes()[site.max(1) - 1] as char
    }

    /// The same peptidoform as an [`RtPeptide`] (Unimod accessions by site).
    pub fn rt_peptide(&self) -> RtPeptide {
        let mut mods: Vec<(usize, u32)> = self.mods.iter().map(|(s, m)| (*s, m.unimod())).collect();
        mods.sort_unstable();
        RtPeptide {
            sequence: self.seq.clone(),
            mods,
        }
    }
}

/// Fragment-mass inputs for one peptidoform.
#[derive(Clone, Debug)]
pub struct Mass {
    /// Residue + residue-mod masses, len n.
    pub res: Vec<f64>,
    /// N-terminal mod mass.
    pub nterm: f64,
}

/// Parse an andes PIN peptide (`K.PEPM+15.99491K.D`) with
/// [`RtPeptide::from_andes`], the parser andes-rt is trained and run with, so
/// every predictor sees the same sites: a delta printed after the first residue
/// goes to the N-terminus (site 0) when an N-terminal Unimod entry matches it
/// (acetyl, TMT, pyro-Glu), as in alphabase, where `Gln->pyro-Glu@Q^Any_N-term`
/// sits at site 0. A summed N-terminal + residue delta on residue 1
/// (`C+99.03203`, acetyl + carbamidomethyl) is split by `from_andes` into both
/// mods. `None` when `from_andes` rejects the string or a modification is
/// outside [`Mod`] at its site. Masses are the Unimod monoisotopic deltas.
pub fn parse_peptide(s: &str) -> Option<(Form, Mass)> {
    let p = RtPeptide::from_andes(s)?;
    let seq = p.sequence.as_bytes();
    let mut res = seq
        .iter()
        .map(|&c| residue_mass(c as char))
        .collect::<Option<Vec<f64>>>()?;
    let mut nterm = 0.0;
    let mut mods = Vec::with_capacity(p.mods.len());
    for &(site, id) in &p.mods {
        let m = Mod::at(id, site, seq)?;
        let delta = rt_native::mods::by_id(id)?.mono;
        if site == 0 {
            nterm += delta;
        } else {
            res[site - 1] += delta;
        }
        mods.push((site, m));
    }
    mods.sort();
    Some((
        Form {
            seq: p.sequence,
            mods,
        },
        Mass { res, nterm },
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_andes_peptides() {
        let (f, m) = parse_peptide("K.M+15.99491EEHILR.P").unwrap();
        assert_eq!(f.seq, "MEEHILR");
        assert_eq!(f.mods, vec![(1, Mod::Oxidation)]);
        assert!((m.res[0] - (131.0404846 + 15.994915)).abs() < 1e-6);
        let (f, m) = parse_peptide("M.R+42.01056AC+57.02146INNQQIR.H").unwrap();
        assert_eq!(f.seq, "RACINNQQIR");
        assert_eq!(f.mods, vec![(0, Mod::Acetyl), (3, Mod::Carbamidomethyl)]);
        assert!((m.nterm - 42.010565).abs() < 1e-9);
        let (f, _) = parse_peptide("_.M+42.01056+15.99491LRK.A").unwrap();
        assert_eq!(f.mods, vec![(0, Mod::Acetyl), (1, Mod::Oxidation)]);
        assert!(parse_peptide("K.PEPXIDE.K").is_none());
        assert!(parse_peptide("K.PEP+12.34TIDE.K").is_none());
    }

    #[test]
    fn places_tmt_and_phospho() {
        let (f, m) = parse_peptide("K.K+229.16293+229.16293PEPTK+229.16293.D").unwrap();
        assert_eq!(
            f.mods,
            vec![(0, Mod::Tmt6plex), (1, Mod::Tmt6plex), (6, Mod::Tmt6plex)]
        );
        assert!((m.nterm - 229.162932).abs() < 1e-9);
        let (f, _) = parse_peptide("K.A+229.16293PES+79.96633K+229.16293.D").unwrap();
        assert_eq!(
            f.mods,
            vec![(0, Mod::Tmt6plex), (4, Mod::Phospho), (5, Mod::Tmt6plex)]
        );
    }

    /// N-terminal mods sit at site 0, as in alphabase and in andes-rt training:
    /// pyro-Glu, and a lone TMT delta printed after a leading K.
    #[test]
    fn places_n_terminal_mods_at_site_zero() {
        let (f, m) = parse_peptide("K.Q-17.02655PEPTIDEK.D").unwrap();
        assert_eq!(f.mods, vec![(0, Mod::GlnPyroGlu)]);
        assert_eq!(
            f.mods[0].1.alphabase(0, f.residue(0)),
            "Gln->pyro-Glu@Q^Any_N-term"
        );
        assert!((m.nterm + 17.026549).abs() < 1e-9);
        let (f, _) = parse_peptide("K.E-18.01056AGHY+79.96633LLR.A").unwrap();
        assert_eq!(f.mods, vec![(0, Mod::GluPyroGlu), (5, Mod::Phospho)]);
        assert!(parse_peptide("K.AQ-17.02655PEPTIDEK.D").is_none());
        assert!(parse_peptide("K.E-17.02655PEPTIDEK.D").is_none());
        let (f, _) = parse_peptide("K.K+229.16293PEPTIDE.D").unwrap();
        assert_eq!(f.mods, vec![(0, Mod::Tmt6plex)]);
        assert_eq!(f.mods[0].1.alphabase(0, 'K'), "TMT6plex@Any_N-term");
        // Acetyl on an internal K is not an N-terminal acetyl.
        assert!(parse_peptide("K.PEPK+42.01056TIDE.D").is_none());
    }

    /// The form and andes-rt read a PIN peptide identically.
    #[test]
    fn matches_from_andes() {
        for s in [
            "K.Q-17.02655PEPTIDEK.D",
            "K.K+229.16293PEPTIDE.D",
            "K.K+229.16293+229.16293PEPTK+229.16293.D",
            "_.M+42.01056+15.99491LRK.A",
            "K.A+229.16293PES+79.96633K+229.16293.D",
            "K.YCC+57.02146QWN+0.98402K.R",
        ] {
            let (f, _) = parse_peptide(s).unwrap();
            assert_eq!(f.rt_peptide(), RtPeptide::from_andes(s).unwrap(), "{s}");
        }
    }

    /// The PIN string andes writes for a protein N-terminal Cys under fixed
    /// carbamidomethyl-C and variable protein N-terminal acetyl.
    fn andes_cam_acetyl_nterm_cys() -> String {
        use model::{AminoAcidSetBuilder, ModLocation, Modification, Peptide};
        let aa_set = AminoAcidSetBuilder::new_standard()
            .add_fixed_mod(
                Modification::from_mods_txt_line("57.021464,C,fix,any,Carbamidomethyl").unwrap(),
            )
            .add_variable_mod(
                Modification::from_mods_txt_line("42.010565,*,opt,Prot-N-term,Acetyl").unwrap(),
            )
            .build()
            .unwrap();
        let first = aa_set
            .variants_for(b'C', ModLocation::ProtNTerm)
            .iter()
            .find(|aa| {
                aa.mod_
                    .as_ref()
                    .is_some_and(|m| (m.mass_delta - 99.032029).abs() < 1e-3)
            })
            .expect("stacked CAM + acetyl variant")
            .clone();
        let mut residues = vec![first];
        for r in b"PEPTIDEK" {
            residues.push(aa_set.standard(*r).unwrap().clone());
        }
        Peptide::new(residues, b'M', b'R').to_string()
    }

    #[test]
    fn parses_summed_nterm_and_residue_delta_from_the_andes_writer() {
        let s = andes_cam_acetyl_nterm_cys();
        assert_eq!(s, "M.C+99.03203PEPTIDEK.R");
        let (f, m) = parse_peptide(&s).expect("summed CAM + acetyl delta");
        assert_eq!(f.seq, "CPEPTIDEK");
        assert_eq!(f.mods, vec![(0, Mod::Acetyl), (1, Mod::Carbamidomethyl)]);
        assert!((m.nterm - 42.010565).abs() < 1e-9);
        assert!((m.res[0] - (103.009_184_5 + 57.021464)).abs() < 1e-6);
        // fixed TMT6plex on K under variable N-terminal acetyl: K+271.17350
        let tmt = "_.K+271.17350PEPTIDEK+229.16293.R";
        let (f, _) = parse_peptide(tmt).unwrap();
        assert_eq!(
            f.mods,
            vec![(0, Mod::Acetyl), (1, Mod::Tmt6plex), (9, Mod::Tmt6plex)]
        );
        assert_eq!(f.rt_peptide(), RtPeptide::from_andes(tmt).unwrap());
    }
}
