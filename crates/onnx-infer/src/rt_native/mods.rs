//! Built-in Unimod elemental compositions (no database lookup).
//!
//! Generated from `scripts/student/mod_composition.py`; the RT parity check
//! fails if the two tables disagree. Heavy isotopes are folded into their
//! element for the feature vector; `mono` is the exact monoisotopic delta.

/// Composition features per position: C, H, N, O, S, P, other, mass / 100.
pub const N_COMP: usize = 8;

/// One Unimod entry.
#[derive(Debug, Clone, Copy)]
pub struct UnimodEntry {
    /// Unimod accession number.
    pub id: u32,
    /// Unimod short name.
    pub name: &'static str,
    /// Allowed sites: residues, `^` = peptide N-terminus, `$` = C-terminus.
    pub sites: &'static str,
    /// Element-count delta: C, H, N, O, S, P, other.
    pub elements: [f32; 7],
    /// Monoisotopic mass delta (Da).
    pub mono: f64,
}

macro_rules! table {
    ($(($id:expr, $n:expr, $s:expr, [$($e:expr),*], $m:expr)),* $(,)?) => {
        &[$(UnimodEntry { id: $id, name: $n, sites: $s, elements: [$($e),*], mono: $m }),*]
    };
}

/// The table, sorted by id.
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

/// Entry by Unimod id.
pub fn by_id(id: u32) -> Option<&'static UnimodEntry> {
    UNIMOD
        .binary_search_by_key(&id, |e| e.id)
        .ok()
        .map(|i| &UNIMOD[i])
}

/// Unimod id for a mass delta at `residue` (or the N-terminus when `nterm`).
///
/// Candidates within `tol` Da must allow the site; ties go to the lowest id.
pub fn match_mass(delta: f64, residue: Option<u8>, nterm: bool, tol: f64) -> Option<u32> {
    UNIMOD
        .iter()
        .filter(|e| (e.mono - delta).abs() <= tol)
        .filter(|e| {
            (nterm && e.sites.contains('^'))
                || residue.is_some_and(|r| e.sites.as_bytes().contains(&r))
        })
        .map(|e| e.id)
        .min()
}

/// Feature vector of one entry (unscaled).
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
}
