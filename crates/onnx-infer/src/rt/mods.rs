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
        4,
        "Carbamidomethyl",
        "C",
        [2.0, 3.0, 1.0, 1.0, 0.0, 0.0, 0.0],
        57.021464
    ),
    (
        5,
        "Carbamyl",
        "^K",
        [1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0],
        43.005814
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
        "MWHC",
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
        40,
        "Sulfo",
        "STY",
        [0.0, 0.0, 0.0, 3.0, 1.0, 0.0, 0.0],
        79.956815
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
