pub const N_COMP: usize = 8;

#[derive(Debug, Clone, Copy)]
pub struct UnimodEntry {
    pub id: u32,
    pub name: &'static str,
    pub sites: &'static str,
    pub elements: [f32; 7],
    pub mono: f64,
}

macro_rules! table {
    ($(($id:expr, $n:expr, $s:expr, [$($e:expr),*], $m:expr)),* $(,)?) => {
        &[$(UnimodEntry { id: $id, name: $n, sites: $s, elements: [$($e),*], mono: $m }),*]
    };
}

pub const UNIMOD: &[UnimodEntry] = table![
    (
        1,
        "Acetyl",
        "^K",
        [2.0, 2.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        42.010565
    ),
    (
        2,
        "Amidated",
        "$",
        [0.0, 1.0, 1.0, -1.0, 0.0, 0.0, 0.0],
        -0.984016
    ),
    (
        3,
        "Biotin",
        "K^",
        [10.0, 14.0, 2.0, 2.0, 1.0, 0.0, 0.0],
        226.077598
    ),
    (
        4,
        "Carbamidomethyl",
        "CKH",
        [2.0, 3.0, 1.0, 1.0, 0.0, 0.0, 0.0],
        57.021464
    ),
    (
        5,
        "Carbamyl",
        "^KR",
        [1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0],
        43.005814
    ),
    (
        6,
        "Carboxymethyl",
        "C",
        [2.0, 2.0, 0.0, 2.0, 0.0, 0.0, 0.0],
        58.005479
    ),
    (
        7,
        "Deamidated",
        "NQR",
        [0.0, -1.0, -1.0, 1.0, 0.0, 0.0, 0.0],
        0.984016
    ),
    (
        21,
        "Phospho",
        "STYH",
        [0.0, 1.0, 0.0, 3.0, 0.0, 1.0, 0.0],
        79.966331
    ),
    (
        23,
        "Dehydrated",
        "STDY",
        [0.0, -2.0, 0.0, -1.0, 0.0, 0.0, 0.0],
        -18.010565
    ),
    (
        24,
        "Propionamide",
        "C",
        [3.0, 5.0, 1.0, 1.0, 0.0, 0.0, 0.0],
        71.037114
    ),
    (
        26,
        "Pyro-carbamidomethyl",
        "^",
        [2.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        39.994915
    ),
    (
        27,
        "Glu->pyro-Glu",
        "^",
        [0.0, -2.0, 0.0, -1.0, 0.0, 0.0, 0.0],
        -18.010565
    ),
    (
        28,
        "Gln->pyro-Glu",
        "^",
        [0.0, -3.0, -1.0, 0.0, 0.0, 0.0, 0.0],
        -17.026549
    ),
    (
        30,
        "Cation:Na",
        "DE$",
        [0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0],
        21.981944
    ),
    (
        34,
        "Methyl",
        "KRHDE^",
        [1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        14.015650
    ),
    (
        35,
        "Oxidation",
        "MWHCP",
        [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        15.994915
    ),
    (
        36,
        "Dimethyl",
        "KR^",
        [2.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        28.031300
    ),
    (
        37,
        "Trimethyl",
        "K",
        [3.0, 6.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        42.046950
    ),
    (
        39,
        "Methylthio",
        "C",
        [1.0, 2.0, 0.0, 0.0, 1.0, 0.0, 0.0],
        45.987721
    ),
    (
        40,
        "Sulfo",
        "STY",
        [0.0, 0.0, 0.0, 3.0, 1.0, 0.0, 0.0],
        79.956815
    ),
    (
        41,
        "Hex",
        "NST",
        [6.0, 10.0, 0.0, 5.0, 0.0, 0.0, 0.0],
        162.052824
    ),
    (
        43,
        "HexNAc",
        "NST",
        [8.0, 13.0, 1.0, 5.0, 0.0, 0.0, 0.0],
        203.079373
    ),
    (
        58,
        "Propionyl",
        "K^",
        [3.0, 4.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        56.026215
    ),
    (
        64,
        "Succinyl",
        "K^",
        [4.0, 4.0, 0.0, 3.0, 0.0, 0.0, 0.0],
        100.016044
    ),
    (
        108,
        "Nethylmaleimide",
        "C",
        [6.0, 7.0, 1.0, 2.0, 0.0, 0.0, 0.0],
        125.047679
    ),
    (
        121,
        "GlyGly",
        "K",
        [4.0, 6.0, 2.0, 2.0, 0.0, 0.0, 0.0],
        114.042927
    ),
    (
        122,
        "Formyl",
        "K^",
        [1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        27.994915
    ),
    (
        188,
        "Label:13C(6)",
        "KR",
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        6.020129
    ),
    (
        199,
        "Dimethyl:2H(4)",
        "KR^",
        [2.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        32.056407
    ),
    (
        213,
        "ADP-Ribosyl",
        "CDEKNRS",
        [15.0, 21.0, 5.0, 13.0, 0.0, 2.0, 0.0],
        541.061110
    ),
    (
        214,
        "iTRAQ4plex",
        "K^Y",
        [7.0, 12.0, 2.0, 1.0, 0.0, 0.0, 0.0],
        144.102062
    ),
    (
        259,
        "Label:13C(6)15N(2)",
        "K",
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        8.014199
    ),
    (
        262,
        "Label:2H(3)",
        "L",
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        3.01883
    ),
    (
        267,
        "Label:13C(6)15N(4)",
        "R",
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        10.008269
    ),
    (
        312,
        "Cysteinyl",
        "C",
        [3.0, 5.0, 1.0, 2.0, 1.0, 0.0, 0.0],
        119.004099
    ),
    (
        330,
        "Dimethyl:2H(6)13C(2)",
        "KR^",
        [2.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        36.07567
    ),
    (
        345,
        "Trioxidation",
        "CW",
        [0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0],
        47.984744
    ),
    (
        351,
        "Trp->Kynurenin",
        "W",
        [-1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        3.994915
    ),
    (
        354,
        "Nitro",
        "YW",
        [0.0, -1.0, 1.0, 2.0, 0.0, 0.0, 0.0],
        44.985078
    ),
    (
        385,
        "Ammonia-loss",
        "^",
        [0.0, -3.0, -1.0, 0.0, 0.0, 0.0, 0.0],
        -17.026549
    ),
    (
        425,
        "Dioxidation",
        "MWC",
        [0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0],
        31.989829
    ),
    (
        481,
        "Label:2H(4)",
        "K",
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        4.025107
    ),
    (
        510,
        "Dimethyl:2H(4)13C(2)",
        "KR^",
        [2.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        34.063117
    ),
    (
        730,
        "iTRAQ8plex",
        "K^Y",
        [14.0, 24.0, 4.0, 3.0, 0.0, 0.0, 0.0],
        304.205359
    ),
    (
        737,
        "TMT6plex",
        "K^",
        [12.0, 20.0, 2.0, 2.0, 0.0, 0.0, 0.0],
        229.162932
    ),
    (
        738,
        "TMT2plex",
        "K^",
        [12.0, 20.0, 2.0, 2.0, 0.0, 0.0, 0.0],
        225.155833
    ),
    (
        739,
        "TMT",
        "K^",
        [12.0, 20.0, 2.0, 2.0, 0.0, 0.0, 0.0],
        224.152478
    ),
    (
        747,
        "Malonyl",
        "K",
        [3.0, 2.0, 0.0, 3.0, 0.0, 0.0, 0.0],
        86.000394
    ),
    (
        765,
        "Met-loss",
        "^",
        [-5.0, -9.0, -1.0, -1.0, -1.0, 0.0, 0.0],
        -131.040485
    ),
    (
        766,
        "Met-loss+Acetyl",
        "^",
        [-3.0, -7.0, -1.0, 0.0, -1.0, 0.0, 0.0],
        -89.029920
    ),
    (
        1363,
        "Crotonyl",
        "K",
        [4.0, 4.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        68.026215
    ),
    (
        2016,
        "TMTpro",
        "K^",
        [15.0, 25.0, 3.0, 3.0, 0.0, 0.0, 0.0],
        304.207145
    ),
];

pub fn by_id(id: u32) -> Option<&'static UnimodEntry> {
    UNIMOD
        .binary_search_by_key(&id, |e| e.id)
        .ok()
        .map(|i| &UNIMOD[i])
}

pub fn match_mass(delta: f64, residue: Option<u8>, nterm: bool, tol: f64) -> Option<u32> {
    UNIMOD
        .iter()
        .filter(|e| (e.mono - delta).abs() <= tol)
        .filter(|e| {
            (nterm && e.sites.contains('^'))
                || residue.is_some_and(|r| e.sites.as_bytes().contains(&r))
        })
        .min_by(|a, b| {
            (a.mono - delta)
                .abs()
                .total_cmp(&(b.mono - delta).abs())
                .then(a.id.cmp(&b.id))
        })
        .map(|e| e.id)
}

pub fn split_nterm_residue(delta: f64, residue: u8, tol: f64) -> Option<(u32, u32)> {
    let nterm: Vec<(u32, f64)> = UNIMOD
        .iter()
        .filter(|e| e.sites.contains('^'))
        .map(|e| (e.id, e.mono))
        .collect();
    let on_residue: Vec<(u32, f64)> = UNIMOD
        .iter()
        .filter(|e| e.sites.as_bytes().contains(&residue))
        .map(|e| (e.id, e.mono))
        .collect();
    crate::peptide::split_nterm_residue_delta(delta, tol, &nterm, &on_residue)
}

pub fn features(e: &UnimodEntry) -> [f32; N_COMP] {
    let mut v = [0.0f32; N_COMP];
    v[..7].copy_from_slice(&e.elements);
    v[7] = (e.mono / 100.0) as f32;
    v
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn table_sorted_and_masses() {
        assert!(UNIMOD.windows(2).all(|w| w[0].id < w[1].id));
        assert!((by_id(21).unwrap().mono - 79.966331).abs() < 1e-5);
        assert!((by_id(737).unwrap().mono - 229.162932).abs() < 1e-5);
        assert_eq!(match_mass(79.96633, Some(b'S'), false, 2e-3), Some(21));
        assert_eq!(match_mass(42.01057, Some(b'M'), true, 2e-3), Some(1));
        assert_eq!(match_mass(42.01057, Some(b'M'), false, 2e-3), None);
        assert_eq!(match_mass(-17.02655, Some(b'Q'), true, 2e-3), Some(28));
    }

    /// Unimod delta compositions as counts of [C, H, N, O, S, P, 13C, 15N, 2H] (unimod.org),
    /// for every entry but Cation:Na.
    const UNIMOD_COMPOSITIONS: &[(u32, [i32; 9])] = &[
        (1, [2, 2, 0, 1, 0, 0, 0, 0, 0]),
        (2, [0, 1, 1, -1, 0, 0, 0, 0, 0]),
        (3, [10, 14, 2, 2, 1, 0, 0, 0, 0]),
        (4, [2, 3, 1, 1, 0, 0, 0, 0, 0]),
        (5, [1, 1, 1, 1, 0, 0, 0, 0, 0]),
        (6, [2, 2, 0, 2, 0, 0, 0, 0, 0]),
        (7, [0, -1, -1, 1, 0, 0, 0, 0, 0]),
        (21, [0, 1, 0, 3, 0, 1, 0, 0, 0]),
        (23, [0, -2, 0, -1, 0, 0, 0, 0, 0]),
        (24, [3, 5, 1, 1, 0, 0, 0, 0, 0]),
        (26, [2, 0, 0, 1, 0, 0, 0, 0, 0]),
        (27, [0, -2, 0, -1, 0, 0, 0, 0, 0]),
        (28, [0, -3, -1, 0, 0, 0, 0, 0, 0]),
        (34, [1, 2, 0, 0, 0, 0, 0, 0, 0]),
        (35, [0, 0, 0, 1, 0, 0, 0, 0, 0]),
        (36, [2, 4, 0, 0, 0, 0, 0, 0, 0]),
        (37, [3, 6, 0, 0, 0, 0, 0, 0, 0]),
        (39, [1, 2, 0, 0, 1, 0, 0, 0, 0]),
        (40, [0, 0, 0, 3, 1, 0, 0, 0, 0]),
        (41, [6, 10, 0, 5, 0, 0, 0, 0, 0]),
        (43, [8, 13, 1, 5, 0, 0, 0, 0, 0]),
        (58, [3, 4, 0, 1, 0, 0, 0, 0, 0]),
        (64, [4, 4, 0, 3, 0, 0, 0, 0, 0]),
        (108, [6, 7, 1, 2, 0, 0, 0, 0, 0]),
        (121, [4, 6, 2, 2, 0, 0, 0, 0, 0]),
        (122, [1, 0, 0, 1, 0, 0, 0, 0, 0]),
        (188, [-6, 0, 0, 0, 0, 0, 6, 0, 0]),
        (199, [2, 0, 0, 0, 0, 0, 0, 0, 4]),
        (213, [15, 21, 5, 13, 0, 2, 0, 0, 0]),
        (214, [4, 12, 1, 1, 0, 0, 3, 1, 0]),
        (259, [-6, 0, -2, 0, 0, 0, 6, 2, 0]),
        (262, [0, -3, 0, 0, 0, 0, 0, 0, 3]),
        (267, [-6, 0, -4, 0, 0, 0, 6, 4, 0]),
        (312, [3, 5, 1, 2, 1, 0, 0, 0, 0]),
        (330, [0, -2, 0, 0, 0, 0, 2, 0, 6]),
        (345, [0, 0, 0, 3, 0, 0, 0, 0, 0]),
        (351, [-1, 0, 0, 1, 0, 0, 0, 0, 0]),
        (354, [0, -1, 1, 2, 0, 0, 0, 0, 0]),
        (385, [0, -3, -1, 0, 0, 0, 0, 0, 0]),
        (425, [0, 0, 0, 2, 0, 0, 0, 0, 0]),
        (481, [0, -4, 0, 0, 0, 0, 0, 0, 4]),
        (510, [0, 0, 0, 0, 0, 0, 2, 0, 4]),
        (730, [7, 24, 3, 3, 0, 0, 7, 1, 0]),
        (737, [8, 20, 1, 2, 0, 0, 4, 1, 0]),
        (738, [11, 20, 2, 2, 0, 0, 1, 0, 0]),
        (739, [12, 20, 2, 2, 0, 0, 0, 0, 0]),
        (747, [3, 2, 0, 3, 0, 0, 0, 0, 0]),
        (765, [-5, -9, -1, -1, -1, 0, 0, 0, 0]),
        (766, [-3, -7, -1, 0, -1, 0, 0, 0, 0]),
        (1363, [4, 4, 0, 1, 0, 0, 0, 0, 0]),
        (2016, [8, 25, 1, 3, 0, 0, 7, 2, 0]),
    ];

    #[test]
    fn every_entry_mass_and_elements_follow_from_its_unimod_composition() {
        use model::mass;
        let atom = [
            mass::C,
            mass::H,
            mass::N,
            mass::O,
            mass::S,
            mass::P,
            mass::C13,
            mass::N15,
            mass::H2,
        ];
        for &(id, n) in UNIMOD_COMPOSITIONS {
            let e = by_id(id).unwrap_or_else(|| panic!("UNIMOD:{id} not in the table"));
            let mono: f64 = n.iter().zip(atom).map(|(&k, m)| k as f64 * m).sum();
            assert!(
                (e.mono - mono).abs() < 1e-5,
                "{} {} vs {mono}",
                e.name,
                e.mono
            );
            let folded = [n[0] + n[6], n[1] + n[8], n[2] + n[7], n[3], n[4], n[5], 0];
            assert_eq!(e.elements, folded.map(|k| k as f32), "{}", e.name);
        }
        let unchecked: Vec<u32> = UNIMOD
            .iter()
            .map(|e| e.id)
            .filter(|id| UNIMOD_COMPOSITIONS.iter().all(|c| c.0 != *id))
            .collect();
        assert_eq!(unchecked, vec![30]);
    }

    #[test]
    fn added_sdrf_mods_resolve_at_their_sites() {
        assert_eq!(match_mass(43.005814, None, true, 2e-3), Some(5));
        assert_eq!(match_mass(43.005814, Some(b'R'), false, 2e-3), Some(5));
        assert_eq!(match_mass(58.005479, Some(b'C'), false, 2e-3), Some(6));
        assert_eq!(match_mass(45.987721, Some(b'C'), false, 2e-3), Some(39));
        assert_eq!(match_mass(125.047679, Some(b'C'), false, 2e-3), Some(108));
        assert_eq!(match_mass(15.994915, Some(b'P'), false, 2e-3), Some(35));
        for site in *b"CDEKNRS" {
            assert_eq!(match_mass(541.06111, Some(site), false, 2e-3), Some(213));
        }
        assert_eq!(match_mass(541.06111, Some(b'Y'), false, 2e-3), None);
        assert_eq!(match_mass(45.987721, Some(b'K'), false, 2e-3), None);
    }

    #[test]
    fn closest_mass_wins_and_id_breaks_exact_ties() {
        assert_eq!(match_mass(304.20715, Some(b'K'), false, 2e-3), Some(2016));
        assert_eq!(match_mass(304.20715, None, true, 2e-3), Some(2016));
        assert_eq!(match_mass(304.20536, Some(b'K'), false, 2e-3), Some(730));
        assert_eq!(match_mass(304.20715, Some(b'Y'), false, 2e-3), Some(730));
        assert_eq!(match_mass(-17.02655, None, true, 2e-3), Some(28));
        assert_eq!(match_mass(-18.01056, Some(b'S'), true, 2e-3), Some(23));
    }
}
