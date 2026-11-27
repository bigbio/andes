use crate::rt::{self, RtPeptide};

pub const PROTON: f64 = model::mass::PROTON;
pub const H2O: f64 = model::mass::H2O;

pub fn residue_mass(c: char) -> Option<f64> {
    let residue = u8::try_from(c).ok().filter(|&r| r != b'U')?;
    model::amino_acid::AminoAcid::standard(residue).map(|a| a.mass)
}

/// A modification of a PIN peptidoform, by Unimod accession. Any entry of the Unimod table
/// [`rt::mods::UNIMOD`] is one: its name, elemental composition, mass and sites come from that
/// table, so a mod the table carries needs no code here.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Mod(u32);

/// Unimod accessions in the order [`Form::mods`] lists mods that share a site; any other
/// accession sorts after these, by accession.
const SITE_ORDER: [u32; 17] = [
    4, 35, 1, 737, 21, 7, 28, 27, 259, 267, 188, 481, 262, 36, 199, 510, 330,
];

/// N-terminal mods of one residue, by Unimod accession, with the N-terminal residues each
/// needs: cyclisations (pyro-Glu from Q or E, pyro-carbamidomethyl and Ammonia-loss from C)
/// and the loss of the initiator Met. Their table entries list the N-terminus (`^`) only.
const NTERM_RESIDUE_MODS: [(u32, &str); 6] = [
    (26, "C"),
    (27, "E"),
    (28, "Q"),
    (385, "C"),
    (765, "M"),
    (766, "M"),
];

/// Heavy dimethyl labels (Unimod accessions), scored as light dimethyl.
const HEAVY_DIMETHYL: [u32; 3] = [199, 330, 510];

impl Ord for Mod {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.site_rank().cmp(&other.site_rank())
    }
}

impl PartialOrd for Mod {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Mod {
    pub const ACETYL: Mod = Mod(1);
    pub const CARBAMIDOMETHYL: Mod = Mod(4);
    pub const DEAMIDATED: Mod = Mod(7);
    pub const PHOSPHO: Mod = Mod(21);
    pub const GLU_PYRO_GLU: Mod = Mod(27);
    pub const GLN_PYRO_GLU: Mod = Mod(28);
    pub const OXIDATION: Mod = Mod(35);
    pub const DIMETHYL: Mod = Mod(36);
    pub const LABEL_13C6: Mod = Mod(188);
    pub const DIMETHYL_2H4: Mod = Mod(199);
    pub const LABEL_13C6_15N2: Mod = Mod(259);
    pub const LABEL_2H3: Mod = Mod(262);
    pub const LABEL_13C6_15N4: Mod = Mod(267);
    pub const DIMETHYL_2H6_13C2: Mod = Mod(330);
    pub const LABEL_2H4: Mod = Mod(481);
    pub const DIMETHYL_2H4_13C2: Mod = Mod(510);
    pub const TMT6PLEX: Mod = Mod(737);

    /// The mod of Unimod accession `id`, if the table carries it with an encodable composition
    /// (C, H, N, O, S, P and their isotopes; not the Na of Cation:Na).
    pub fn from_unimod(id: u32) -> Option<Self> {
        rt::mods::by_id(id)
            .filter(|e| e.elements[6] == 0.0)
            .map(|e| Mod(e.id))
    }

    fn site_rank(self) -> (usize, u32) {
        let listed = SITE_ORDER.iter().position(|&i| i == self.0);
        (listed.unwrap_or(SITE_ORDER.len()), self.0)
    }

    fn entry(self) -> &'static rt::mods::UnimodEntry {
        rt::mods::by_id(self.0).expect("a Mod is built only from an entry of the Unimod table")
    }

    /// The N-terminal residues this mod needs when it is an N-terminal mod of one residue
    /// (see [`NTERM_RESIDUE_MODS`]).
    pub fn nterm_residues(self) -> Option<&'static str> {
        NTERM_RESIDUE_MODS
            .iter()
            .find(|c| c.0 == self.0)
            .map(|c| c.1)
    }

    /// The AlphaBase name of the mod at `site` (0 = N-terminus) on residue `aa`, the key of
    /// `peptdeep_mods.json`: `Name@Residue`, `Name@Any_N-term`, `Name@Q^Any_N-term` for an
    /// N-terminal mod of one residue, `Acetyl@Protein_N-term`, and heavy dimethyl under its light name.
    pub fn alphabase(self, site: usize, aa: char) -> String {
        let named = if HEAVY_DIMETHYL.contains(&self.0) {
            Mod::DIMETHYL
        } else {
            self
        };
        let name = named.entry().name;
        match self {
            Mod::ACETYL if site == 0 => "Acetyl@Protein_N-term".into(),
            m if m.nterm_residues().is_some() => format!("{name}@{aa}^Any_N-term"),
            _ if site == 0 => format!("{name}@Any_N-term"),
            _ => format!("{name}@{aa}"),
        }
    }

    pub fn light_equivalent(self) -> Option<Self> {
        match self {
            m if m.is_isotope_label() => None,
            m if HEAVY_DIMETHYL.contains(&m.0) => Some(Mod::DIMETHYL),
            m => Some(m),
        }
    }

    /// A stable-isotope label: a mass shift with no change in elemental composition.
    pub fn is_isotope_label(self) -> bool {
        self.entry().elements.iter().all(|&k| k == 0.0)
    }

    pub fn unimod(self) -> u32 {
        self.0
    }

    fn at(id: u32, site: usize, seq: &[u8]) -> Option<Self> {
        let m = Self::from_unimod(id)?;
        let sites = m.entry().sites.as_bytes();
        let ok = match m.nterm_residues() {
            Some(residues) => site == 0 && residues.as_bytes().contains(seq.first()?),
            None if site == 0 => sites.contains(&b'^'),
            None => site <= seq.len() && sites.contains(&seq[site - 1]),
        };
        ok.then_some(m)
    }

    /// The mod that Unimod `id` at `site` (0 = N-terminus) is on this peptide, with the site
    /// it goes to. An N-terminal mass that the entry `id` cannot carry here is read as an entry
    /// of the same mass that can, at the N-terminus or on the first residue: Ammonia-loss on an
    /// N-terminal C for Gln->pyro-Glu, Dehydrated S for Glu->pyro-Glu.
    fn place(id: u32, site: usize, seq: &[u8]) -> Option<(usize, Self)> {
        let mono = rt::mods::by_id(id)?.mono.to_bits();
        let ids: Vec<u32> = std::iter::once(id)
            .chain(
                rt::mods::UNIMOD
                    .iter()
                    .filter(|e| e.id != id && e.mono.to_bits() == mono)
                    .map(|e| e.id),
            )
            .collect();
        let sites: &[usize] = if site == 0 { &[0, 1] } else { &[site] };
        sites.iter().find_map(|&s| {
            ids.iter()
                .find_map(|&k| Self::at(k, s, seq).map(|m| (s, m)))
        })
    }

    pub fn composition(self) -> [i32; 6] {
        let e = self.entry().elements;
        std::array::from_fn(|i| e[i] as i32)
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
        let (site, m) = Mod::place(id, site, seq)?;
        let delta = m.entry().mono;
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
        assert_eq!(f.mods, vec![(1, Mod::OXIDATION)]);
        assert!((m.res[0] - (131.0404846 + 15.994915)).abs() < 1e-6);
        let (f, m) = parse_peptide("M.R+42.01056AC+57.02146INNQQIR.H").unwrap();
        assert_eq!(f.seq, "RACINNQQIR");
        assert_eq!(f.mods, vec![(0, Mod::ACETYL), (3, Mod::CARBAMIDOMETHYL)]);
        assert!((m.nterm - 42.010565).abs() < 1e-9);
        let (f, _) = parse_peptide("_.M+42.01056+15.99491LRK.A").unwrap();
        assert_eq!(f.mods, vec![(0, Mod::ACETYL), (1, Mod::OXIDATION)]);
        assert!(parse_peptide("K.PEPXIDE.K").is_none());
        assert!(parse_peptide("K.PEP+12.34TIDE.K").is_none());
    }

    #[test]
    fn places_tmt_and_phospho() {
        let (f, m) = parse_peptide("K.K+229.16293+229.16293PEPTK+229.16293.D").unwrap();
        assert_eq!(
            f.mods,
            vec![(0, Mod::TMT6PLEX), (1, Mod::TMT6PLEX), (6, Mod::TMT6PLEX)]
        );
        assert!((m.nterm - 229.162932).abs() < 1e-9);
        let (f, _) = parse_peptide("K.A+229.16293PES+79.96633K+229.16293.D").unwrap();
        assert_eq!(
            f.mods,
            vec![(0, Mod::TMT6PLEX), (4, Mod::PHOSPHO), (5, Mod::TMT6PLEX)]
        );
    }

    #[test]
    fn places_n_terminal_mods_at_site_zero() {
        let (f, m) = parse_peptide("K.Q-17.02655PEPTIDEK.D").unwrap();
        assert_eq!(f.mods, vec![(0, Mod::GLN_PYRO_GLU)]);
        assert_eq!(
            f.mods[0].1.alphabase(0, f.residue(0)),
            "Gln->pyro-Glu@Q^Any_N-term"
        );
        assert!((m.nterm + 17.026549).abs() < 1e-9);
        let (f, _) = parse_peptide("K.E-18.01056AGHY+79.96633LLR.A").unwrap();
        assert_eq!(f.mods, vec![(0, Mod::GLU_PYRO_GLU), (5, Mod::PHOSPHO)]);
        assert!(parse_peptide("K.AQ-17.02655PEPTIDEK.D").is_none());
        assert!(parse_peptide("K.E-17.02655PEPTIDEK.D").is_none());
        let (f, _) = parse_peptide("K.K+229.16293PEPTIDE.D").unwrap();
        assert_eq!(f.mods, vec![(0, Mod::TMT6PLEX)]);
        assert_eq!(f.mods[0].1.alphabase(0, 'K'), "TMT6plex@Any_N-term");
        let (f, _) = parse_peptide("K.PEPK+42.01056TIDE.D").unwrap();
        assert_eq!(f.mods, vec![(4, Mod::ACETYL)]);
        assert_eq!(f.mods[0].1.alphabase(4, 'K'), "Acetyl@K");
    }

    #[test]
    fn label_variants_score_as_their_light_form() {
        let (f, m) = parse_peptide("K.PEPTK+4.02511.D").expect("Lys4");
        assert_eq!(f.mods, vec![(5, Mod::LABEL_2H4)]);
        assert!((m.res[4] - (128.094963 + 4.025107)).abs() < 1e-6);
        assert!(f.rt_peptide().mods.is_empty());
        let (f, _) = parse_peptide("K.PEPTR+6.02013.D").expect("Arg6");
        assert_eq!(f.mods, vec![(5, Mod::LABEL_13C6)]);
        assert!(f.rt_peptide().mods.is_empty());
        let (f, m) = parse_peptide("K.P+36.07567EPTK+36.07567.D").expect("Dimethyl8");
        assert_eq!(
            f.mods,
            vec![(0, Mod::DIMETHYL_2H6_13C2), (5, Mod::DIMETHYL_2H6_13C2)]
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
            vec![(3, Mod::CARBAMIDOMETHYL), (6, Mod::LABEL_13C6_15N2)]
        );
        assert!((m.res[5] - (128.094_963 + 8.014199)).abs() < 1e-6);
        let (light, _) = parse_peptide("K.ELC+57.02146VAK.R").unwrap();
        assert_eq!(f.rt_peptide(), light.rt_peptide());
        assert!(f.rt_peptide().encode().is_some());
        let (f, m) = parse_peptide("K.M+15.99491PEPTIDER+10.00827.A").expect("heavy R parses");
        assert_eq!(f.mods, vec![(1, Mod::OXIDATION), (9, Mod::LABEL_13C6_15N4)]);
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
        use model::modification::spec::parse_mod_list;
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
        assert_eq!(f.mods, vec![(0, Mod::ACETYL), (1, Mod::CARBAMIDOMETHYL)]);
        assert!((m.nterm - 42.010565).abs() < 1e-9);
        assert!((m.res[0] - (103.009_184_5 + 57.021464)).abs() < 1e-6);
        let tmt = "_.K+271.17350PEPTIDEK+229.16293.R";
        let (f, _) = parse_peptide(tmt).unwrap();
        assert_eq!(
            f.mods,
            vec![(0, Mod::ACETYL), (1, Mod::TMT6PLEX), (9, Mod::TMT6PLEX)]
        );
        assert_eq!(f.rt_peptide(), RtPeptide::from_andes(tmt).unwrap());
    }

    #[test]
    fn mods_listed_before_the_unimod_table_keep_their_names_and_compositions() {
        let cases: [(Mod, usize, char, &str, [i32; 6]); 17] = [
            (
                Mod::CARBAMIDOMETHYL,
                1,
                'C',
                "Carbamidomethyl@C",
                [2, 3, 1, 1, 0, 0],
            ),
            (Mod::OXIDATION, 1, 'M', "Oxidation@M", [0, 0, 0, 1, 0, 0]),
            (
                Mod::ACETYL,
                0,
                'A',
                "Acetyl@Protein_N-term",
                [2, 2, 0, 1, 0, 0],
            ),
            (
                Mod::TMT6PLEX,
                0,
                'A',
                "TMT6plex@Any_N-term",
                [12, 20, 2, 2, 0, 0],
            ),
            (Mod::PHOSPHO, 2, 'T', "Phospho@T", [0, 1, 0, 3, 0, 1]),
            (
                Mod::DEAMIDATED,
                2,
                'N',
                "Deamidated@N",
                [0, -1, -1, 1, 0, 0],
            ),
            (
                Mod::GLN_PYRO_GLU,
                0,
                'Q',
                "Gln->pyro-Glu@Q^Any_N-term",
                [0, -3, -1, 0, 0, 0],
            ),
            (
                Mod::GLU_PYRO_GLU,
                0,
                'E',
                "Glu->pyro-Glu@E^Any_N-term",
                [0, -2, 0, -1, 0, 0],
            ),
            (Mod::LABEL_13C6_15N2, 3, 'K', "Label:13C(6)15N(2)@K", [0; 6]),
            (Mod::LABEL_13C6_15N4, 3, 'R', "Label:13C(6)15N(4)@R", [0; 6]),
            (Mod::LABEL_13C6, 3, 'R', "Label:13C(6)@R", [0; 6]),
            (Mod::LABEL_2H4, 3, 'K', "Label:2H(4)@K", [0; 6]),
            (Mod::LABEL_2H3, 3, 'L', "Label:2H(3)@L", [0; 6]),
            (Mod::DIMETHYL, 3, 'K', "Dimethyl@K", [2, 4, 0, 0, 0, 0]),
            (
                Mod::DIMETHYL_2H4,
                0,
                'A',
                "Dimethyl@Any_N-term",
                [2, 4, 0, 0, 0, 0],
            ),
            (
                Mod::DIMETHYL_2H4_13C2,
                3,
                'K',
                "Dimethyl@K",
                [2, 4, 0, 0, 0, 0],
            ),
            (
                Mod::DIMETHYL_2H6_13C2,
                0,
                'A',
                "Dimethyl@Any_N-term",
                [2, 4, 0, 0, 0, 0],
            ),
        ];
        for (m, site, aa, name, comp) in cases {
            assert_eq!(m.alphabase(site, aa), name);
            assert_eq!(m.composition(), comp, "{name}");
        }
        let labels: Vec<u32> = rt::mods::UNIMOD
            .iter()
            .filter_map(|e| Mod::from_unimod(e.id))
            .filter(|m| m.is_isotope_label())
            .map(Mod::unimod)
            .collect();
        assert_eq!(labels, vec![188, 259, 262, 267, 481]);
        assert!(Mod::from_unimod(30).is_none());
        let mut order = vec![
            Mod::from_unimod(5).unwrap(),
            Mod::TMT6PLEX,
            Mod::CARBAMIDOMETHYL,
        ];
        order.sort();
        assert_eq!(order, vec![Mod::CARBAMIDOMETHYL, Mod::TMT6PLEX, Mod(5)]);
    }

    #[test]
    fn carbamylated_peptides_of_an_urea_digest_are_encoded() {
        let (f, m) = parse_peptide("K.A+43.00581AQVFK.L").unwrap();
        assert_eq!(f.mods, vec![(0, Mod(5))]);
        assert_eq!(f.mods[0].1.alphabase(0, 'A'), "Carbamyl@Any_N-term");
        assert!((m.nterm - 43.005814).abs() < 1e-9);
        assert!(f.rt_peptide().encode().is_some());
        let (f, m) = parse_peptide("R.M+59.00073LPLPR.L").unwrap();
        assert_eq!(f.mods, vec![(0, Mod(5)), (1, Mod::OXIDATION)]);
        assert!((m.nterm - 43.005814).abs() < 1e-9);
        assert!((m.res[0] - (131.040_485 + 15.994915)).abs() < 1e-5);
        let (f, _) = parse_peptide("K.C+100.02728PEPTIDEK.R").unwrap();
        assert_eq!(f.mods, vec![(0, Mod(5)), (1, Mod::CARBAMIDOMETHYL)]);
        let (f, _) = parse_peptide("K.PEPR+43.00581TIDEK.R").unwrap();
        assert_eq!(f.mods, vec![(4, Mod(5))]);
    }

    #[test]
    fn an_n_terminal_mass_goes_to_the_isobaric_mod_its_residue_carries() {
        let (f, m) = parse_peptide("K.C+39.99491PEPTIDEK.R").unwrap();
        assert_eq!(f.mods, vec![(0, Mod(26))]);
        assert_eq!(
            f.mods[0].1.alphabase(0, 'C'),
            "Pyro-carbamidomethyl@C^Any_N-term"
        );
        assert!((m.nterm - 39.994915).abs() < 1e-9);
        let (f, _) = parse_peptide("K.C+57.02146-17.02655PEPTIDEK.R").unwrap();
        assert_eq!(f.mods, vec![(0, Mod(385)), (1, Mod::CARBAMIDOMETHYL)]);
        let (f, m) = parse_peptide("K.S-18.01056PEPTIDEK.R").unwrap();
        assert_eq!(f.mods, vec![(1, Mod(23))]);
        assert!((m.nterm).abs() < 1e-12);
        assert!((m.res[0] - (87.032_028 - 18.010565)).abs() < 1e-5);
        assert!(parse_peptide("K.E-17.02655PEPTIDEK.D").is_none());
    }

    fn andes_peptide(
        variable: &str,
        first: u8,
        location: model::ModLocation,
        delta: f64,
    ) -> String {
        use model::modification::spec::parse_mod_list;
        use model::{AminoAcidSetBuilder, ModLocation, Peptide};
        let aa_set = AminoAcidSetBuilder::new_standard()
            .add_mods(parse_mod_list(variable, false).unwrap())
            .build()
            .unwrap();
        let (site, rest): (usize, &[u8]) = match location {
            ModLocation::NTerm | ModLocation::ProtNTerm => (0, b"PEPTIDEK"),
            _ => (2, b"PEPTIDEK"),
        };
        let modified = aa_set
            .variants_for(first, location)
            .iter()
            .find(|aa| {
                aa.mod_
                    .as_ref()
                    .is_some_and(|m| (m.mass_delta - delta).abs() < 1e-3)
            })
            .unwrap_or_else(|| panic!("{variable}: no variant of {} with {delta}", first as char))
            .clone();
        let mut residues: Vec<_> = rest
            .iter()
            .map(|r| aa_set.standard(*r).unwrap().clone())
            .collect();
        residues.insert(site, modified);
        Peptide::new(residues, b'K', b'R').to_string()
    }

    #[test]
    fn common_sdrf_modifications_written_by_andes_are_encoded_by_the_predictors() {
        use model::ModLocation::{Anywhere, NTerm, ProtNTerm};
        let cases = [
            ("Carbamyl (N-term)", b'A', NTerm, 5),
            ("Carbamyl (K)", b'K', Anywhere, 5),
            ("Carbamyl (R)", b'R', Anywhere, 5),
            ("Acetyl (Protein N-term)", b'A', ProtNTerm, 1),
            ("Acetyl (K)", b'K', Anywhere, 1),
            ("Gln->pyro-Glu (N-term Q)", b'Q', NTerm, 28),
            ("Glu->pyro-Glu (N-term E)", b'E', NTerm, 27),
            ("Ammonia-loss (N-term C)", b'C', NTerm, 385),
            ("Oxidation (M)", b'M', Anywhere, 35),
            ("Oxidation (P)", b'P', Anywhere, 35),
            ("Oxidation (W)", b'W', Anywhere, 35),
            ("Dioxidation (M)", b'M', Anywhere, 425),
            ("Deamidated (N)", b'N', Anywhere, 7),
            ("Deamidated (R)", b'R', Anywhere, 7),
            ("Phospho (S)", b'S', Anywhere, 21),
            ("Phospho (Y)", b'Y', Anywhere, 21),
            ("GlyGly (K)", b'K', Anywhere, 121),
            ("Methyl (K)", b'K', Anywhere, 34),
            ("Methyl (R)", b'R', Anywhere, 34),
            ("Dimethyl (K)", b'K', Anywhere, 36),
            ("Dimethyl (N-term)", b'A', NTerm, 36),
            ("Trimethyl (K)", b'K', Anywhere, 37),
            ("Formyl (N-term)", b'A', NTerm, 122),
            ("Propionyl (K)", b'K', Anywhere, 58),
            ("Succinyl (K)", b'K', Anywhere, 64),
            ("Crotonyl (K)", b'K', Anywhere, 1363),
            ("Malonyl (K)", b'K', Anywhere, 747),
            ("Nitro (Y)", b'Y', Anywhere, 354),
            ("Carboxymethyl (C)", b'C', Anywhere, 6),
            ("Methylthio (C)", b'C', Anywhere, 39),
            ("Nethylmaleimide (C)", b'C', Anywhere, 108),
            ("Propionamide (C)", b'C', Anywhere, 24),
            ("ADP-Ribosyl (S)", b'S', Anywhere, 213),
            ("TMT6plex (K)", b'K', Anywhere, 737),
            ("TMTpro (K)", b'K', Anywhere, 2016),
            ("TMTpro (N-term)", b'A', NTerm, 2016),
            ("iTRAQ4plex (K)", b'K', Anywhere, 214),
            ("iTRAQ8plex (N-term)", b'A', NTerm, 730),
            ("Met-loss (Protein N-term M)", b'M', ProtNTerm, 765),
            ("Met-loss+Acetyl (Protein N-term M)", b'M', ProtNTerm, 766),
            ("Trp->Kynurenin (W)", b'W', Anywhere, 351),
            ("Trioxidation (C)", b'C', Anywhere, 345),
            ("HexNAc (S)", b'S', Anywhere, 43),
            ("Hex (T)", b'T', Anywhere, 41),
            ("Biotin (K)", b'K', Anywhere, 3),
            ("Carbamidomethyl (K)", b'K', Anywhere, 4),
            ("Carbamidomethyl (H)", b'H', Anywhere, 4),
        ];
        for (spec, first, location, id) in cases {
            let delta = rt::mods::by_id(id).unwrap().mono;
            let s = andes_peptide(spec, first, location, delta);
            let (f, m) = parse_peptide(&s).unwrap_or_else(|| panic!("{spec}: {s} not parsed"));
            assert!(
                f.mods.iter().any(|(_, m)| m.unimod() == id),
                "{spec}: {s} {:?}",
                f.mods
            );
            let total = m.nterm + m.res.iter().sum::<f64>();
            let unmodified: f64 = f.seq.chars().filter_map(residue_mass).sum();
            let expect = unmodified + delta;
            assert!((total - expect).abs() < 1e-4, "{spec}: {total} vs {expect}");
            assert!(f.rt_peptide().encode().is_some(), "{spec}: {s} not encoded");
        }
    }
}
