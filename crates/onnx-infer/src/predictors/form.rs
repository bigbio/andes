use crate::rt_native::{self, RtPeptide};

pub const PROTON: f64 = model::mass::PROTON;
pub const H2O: f64 = model::mass::H2O;

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

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Mod {
    Carbamidomethyl,
    Oxidation,
    Acetyl,
    Tmt6plex,
    Phospho,
    Deamidated,
    GlnPyroGlu,
    GluPyroGlu,
    Label13C6N2,
    Label13C6N4,
    Label13C6,
    Label2H4,
    Label2H3,
    Dimethyl,
    DimethylD4,
    DimethylD4C13,
    DimethylD6C13,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Placement {
    Residue,
    NTerm,
    ResidueOrNTerm,
    NTermOf,
}

const MODS: &[(Mod, u32, &str, Placement)] = &[
    (Mod::Carbamidomethyl, 4, "C", Placement::Residue),
    (Mod::Oxidation, 35, "M", Placement::Residue),
    (Mod::Acetyl, 1, "", Placement::NTerm),
    (Mod::Tmt6plex, 737, "K", Placement::ResidueOrNTerm),
    (Mod::Phospho, 21, "STY", Placement::Residue),
    (Mod::Deamidated, 7, "NQ", Placement::Residue),
    (Mod::GlnPyroGlu, 28, "Q", Placement::NTermOf),
    (Mod::GluPyroGlu, 27, "E", Placement::NTermOf),
    (Mod::Label13C6N2, 259, "K", Placement::Residue),
    (Mod::Label13C6N4, 267, "R", Placement::Residue),
    (Mod::Label13C6, 188, "KR", Placement::Residue),
    (Mod::Label2H4, 481, "K", Placement::Residue),
    (Mod::Label2H3, 262, "L", Placement::Residue),
    (Mod::Dimethyl, 36, "KR", Placement::ResidueOrNTerm),
    (Mod::DimethylD4, 199, "KR", Placement::ResidueOrNTerm),
    (Mod::DimethylD4C13, 510, "KR", Placement::ResidueOrNTerm),
    (Mod::DimethylD6C13, 330, "KR", Placement::ResidueOrNTerm),
];

impl Mod {
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
            Mod::Label13C6N2 => "Label:13C(6)15N(2)@K".into(),
            Mod::Label13C6N4 => "Label:13C(6)15N(4)@R".into(),
            Mod::Label13C6 => format!("Label:13C(6)@{aa}"),
            Mod::Label2H4 => "Label:2H(4)@K".into(),
            Mod::Label2H3 => "Label:2H(3)@L".into(),
            Mod::Dimethyl | Mod::DimethylD4 | Mod::DimethylD4C13 | Mod::DimethylD6C13
                if site == 0 =>
            {
                "Dimethyl@Any_N-term".into()
            }
            Mod::Dimethyl | Mod::DimethylD4 | Mod::DimethylD4C13 | Mod::DimethylD6C13 => {
                format!("Dimethyl@{aa}")
            }
        }
    }

    pub fn light_equivalent(self) -> Option<Self> {
        match self {
            m if m.is_isotope_label() => None,
            Mod::DimethylD4 | Mod::DimethylD4C13 | Mod::DimethylD6C13 => Some(Mod::Dimethyl),
            m => Some(m),
        }
    }

    pub fn is_isotope_label(self) -> bool {
        matches!(
            self,
            Mod::Label13C6N2 | Mod::Label13C6N4 | Mod::Label13C6 | Mod::Label2H4 | Mod::Label2H3
        )
    }

    pub fn unimod(self) -> u32 {
        MODS.iter()
            .find(|m| m.0 == self)
            .expect("every Mod is in MODS")
            .1
    }

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
            Mod::Label13C6N2
            | Mod::Label13C6N4
            | Mod::Label13C6
            | Mod::Label2H4
            | Mod::Label2H3 => [0, 0, 0, 0, 0, 0],
            Mod::Dimethyl | Mod::DimethylD4 | Mod::DimethylD4C13 | Mod::DimethylD6C13 => {
                [2, 4, 0, 0, 0, 0]
            }
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Form {
    pub seq: String,
    pub mods: Vec<(usize, Mod)>,
}

impl Form {
    pub fn residue(&self, site: usize) -> char {
        self.seq.as_bytes()[site.max(1) - 1] as char
    }

    pub fn rt_peptide(&self) -> RtPeptide {
        let mut mods: Vec<(usize, u32)> = self
            .mods
            .iter()
            .filter_map(|(s, m)| m.light_equivalent().map(|l| (*s, l.unimod())))
            .collect();
        mods.sort_unstable();
        RtPeptide {
            sequence: self.seq.clone(),
            mods,
        }
    }
}

#[derive(Clone, Debug)]
pub struct Mass {
    pub res: Vec<f64>,
    pub nterm: f64,
}

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
        assert!(parse_peptide("K.PEPK+42.01056TIDE.D").is_none());
    }

    #[test]
    fn label_variants_score_as_their_light_form() {
        let (f, m) = parse_peptide("K.PEPTK+4.02511.D").expect("Lys4");
        assert_eq!(f.mods, vec![(5, Mod::Label2H4)]);
        assert!((m.res[4] - (128.094963 + 4.025107)).abs() < 1e-6);
        assert!(f.rt_peptide().mods.is_empty());
        let (f, _) = parse_peptide("K.PEPTR+6.02013.D").expect("Arg6");
        assert_eq!(f.mods, vec![(5, Mod::Label13C6)]);
        assert!(f.rt_peptide().mods.is_empty());
        let (f, m) = parse_peptide("K.P+36.07567EPTK+36.07567.D").expect("Dimethyl8");
        assert_eq!(
            f.mods,
            vec![(0, Mod::DimethylD6C13), (5, Mod::DimethylD6C13)]
        );
        assert!((m.nterm - 36.07567).abs() < 1e-9);
        assert_eq!(f.rt_peptide().mods, vec![(0, 36), (5, 36)]);
        let (light, _) = parse_peptide("K.P+28.03130EPTK+28.03130.D").expect("Dimethyl0");
        assert_eq!(light.rt_peptide(), f.rt_peptide());
        let (f, _) = parse_peptide("K.P+32.05641EPTK+32.05641.D").expect("Dimethyl4");
        assert_eq!(f.rt_peptide(), light.rt_peptide());
        assert_eq!(f.mods[0].1.alphabase(0, 'P'), "Dimethyl@Any_N-term");
    }

    #[test]
    fn silac_labels_parse_and_score_as_the_light_form() {
        let heavy = "K.ELC+57.02146VAK+8.01420.R";
        let (f, m) = parse_peptide(heavy).expect("heavy K parses");
        assert_eq!(f.seq, "ELCVAK");
        assert_eq!(
            f.mods,
            vec![(3, Mod::Carbamidomethyl), (6, Mod::Label13C6N2)]
        );
        assert!((m.res[5] - (128.094_963 + 8.014199)).abs() < 1e-6);
        let (light, _) = parse_peptide("K.ELC+57.02146VAK.R").unwrap();
        assert_eq!(f.rt_peptide(), light.rt_peptide());
        assert!(f.rt_peptide().encode().is_some());
        let (f, m) = parse_peptide("K.M+15.99491PEPTIDER+10.00827.A").expect("heavy R parses");
        assert_eq!(f.mods, vec![(1, Mod::Oxidation), (9, Mod::Label13C6N4)]);
        assert!((m.res[8] - (156.101_111 + 10.008269)).abs() < 1e-6);
        let (light, _) = parse_peptide("K.M+15.99491PEPTIDER.A").unwrap();
        assert_eq!(f.rt_peptide(), light.rt_peptide());
        assert!(parse_peptide("K.PEPTIDEK+10.00827.A").is_none());
    }

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

    fn andes_cam_acetyl_nterm_cys() -> String {
        use model::mod_spec::parse_mod_list;
        use model::{AminoAcidSetBuilder, ModLocation, Peptide};
        let aa_set = AminoAcidSetBuilder::new_standard()
            .add_mods(parse_mod_list("Carbamidomethyl (C)", true).unwrap())
            .add_mods(parse_mod_list("Acetyl (Protein N-term)", false).unwrap())
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
        let tmt = "_.K+271.17350PEPTIDEK+229.16293.R";
        let (f, _) = parse_peptide(tmt).unwrap();
        assert_eq!(
            f.mods,
            vec![(0, Mod::Acetyl), (1, Mod::Tmt6plex), (9, Mod::Tmt6plex)]
        );
        assert_eq!(f.rt_peptide(), RtPeptide::from_andes(tmt).unwrap());
    }
}
