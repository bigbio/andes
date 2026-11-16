use std::cmp::Ordering;
use std::collections::HashMap;

use crate::glycan_db::GlycanComp;
use crate::glycan_mass::{FUC, HEX, HEXNAC, NEUAC, NEUGC, PROTON};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Monosaccharide {
    Hex,
    HexNAc,
    Fuc,
    NeuAc,
    NeuGc,
}

impl Monosaccharide {
    pub fn residue_mass(self) -> f64 {
        match self {
            Monosaccharide::Hex => HEX,
            Monosaccharide::HexNAc => HEXNAC,
            Monosaccharide::Fuc => FUC,
            Monosaccharide::NeuAc => NEUAC,
            Monosaccharide::NeuGc => NEUGC,
        }
    }

    pub fn neutral_mass(self) -> f64 {
        self.residue_mass() + crate::glycan_mass::WATER
    }

    pub fn symbol(self) -> &'static str {
        match self {
            Monosaccharide::Hex => "Hex",
            Monosaccharide::HexNAc => "HexNAc",
            Monosaccharide::Fuc => "Fuc",
            Monosaccharide::NeuAc => "NeuAc",
            Monosaccharide::NeuGc => "NeuGc",
        }
    }
}

#[derive(Debug, Clone)]
pub struct Glycan {
    pub id: u32,
    pub composition: GlycanComp,
}

#[derive(Debug, Clone)]
struct GlycanIon {
    glycan_id: u32,
    mass: f64,
    is_core: bool,
}

#[derive(Debug, Clone, Copy)]
pub struct DiagnosticIon {
    pub mz: f64,
    pub hex: u8,
    pub hexnac: u8,
    pub fuc: u8,
    pub neuac: u8,
    pub neugc: u8,
}

impl DiagnosticIon {
    pub fn consistent_with(&self, comp: &GlycanComp) -> bool {
        comp.hex >= self.hex
            && comp.hexnac >= self.hexnac
            && comp.fuc >= self.fuc
            && comp.neuac >= self.neuac
            && comp.neugc >= self.neugc
    }
}

pub fn default_diagnostic_ions() -> Vec<DiagnosticIon> {
    use crate::glycan_mass::WATER;
    const IONS: [(f64, u8, u8, u8, u8, u8); 28] = [
        (84.04439, 0, 1, 0, 0, 0),
        (85.02841, 1, 0, 0, 0, 0),
        (97.02841, 1, 0, 0, 0, 0),
        (HEX + PROTON - 3.0 * WATER, 1, 0, 0, 0, 0),
        (115.03897, 1, 0, 0, 0, 0),
        (126.05495, 0, 1, 0, 0, 0),
        (HEX + PROTON - 2.0 * WATER, 1, 0, 0, 0, 0),
        (138.05495, 0, 1, 0, 0, 0),
        (144.06552, 0, 1, 0, 0, 0),
        (HEX + PROTON - WATER, 1, 0, 0, 0, 0),
        (HEX + PROTON, 1, 0, 0, 0, 0),
        (HEXNAC + PROTON - 2.0 * WATER, 0, 1, 0, 0, 0),
        (HEXNAC + PROTON - WATER, 0, 1, 0, 0, 0),
        (HEXNAC + PROTON, 0, 1, 0, 0, 0),
        (NEUAC + PROTON - WATER, 0, 0, 0, 1, 0),
        (NEUGC + PROTON - WATER, 0, 0, 0, 0, 1),
        (NEUAC + PROTON, 0, 0, 0, 1, 0),
        (NEUGC + PROTON, 0, 0, 0, 0, 1),
        (2.0 * HEX + PROTON, 2, 0, 0, 0, 0),
        (HEX + HEXNAC + PROTON, 1, 1, 0, 0, 0),
        (HEX + NEUAC + PROTON, 1, 0, 0, 1, 0),
        (HEX + HEXNAC + FUC + PROTON, 1, 1, 1, 0, 0),
        (2.0 * HEX + HEXNAC + PROTON, 2, 1, 0, 0, 0),
        (HEX + HEXNAC + NEUAC + PROTON, 1, 1, 0, 1, 0),
        (3.0 * HEX + HEXNAC + PROTON, 3, 1, 0, 0, 0),
        (HEX + HEXNAC + NEUAC + FUC + PROTON, 1, 1, 1, 1, 0),
        (3.0 * HEX + 2.0 * HEXNAC + PROTON, 3, 2, 0, 0, 0),
        (3.0 * HEX + 2.0 * HEXNAC + FUC + PROTON, 3, 2, 1, 0, 0),
    ];
    IONS.iter()
        .map(|&(mz, hex, hexnac, fuc, neuac, neugc)| DiagnosticIon {
            mz,
            hex,
            hexnac,
            fuc,
            neuac,
            neugc,
        })
        .collect()
}

pub fn o_glycan_diagnostic_ions() -> Vec<DiagnosticIon> {
    const IONS: [(f64, u8, u8, u8, u8, u8); 7] = [
        (138.0550, 0, 1, 0, 0, 0),
        (168.0660, 0, 1, 0, 0, 0),
        (186.0760, 0, 1, 0, 0, 0),
        (204.0867, 0, 1, 0, 0, 0),
        (366.1395, 1, 1, 0, 0, 0),
        (495.1821, 0, 1, 0, 1, 0),
        (657.2349, 1, 1, 0, 1, 0),
    ];
    IONS.iter()
        .map(|&(mz, hex, hexnac, fuc, neuac, neugc)| DiagnosticIon {
            mz,
            hex,
            hexnac,
            fuc,
            neuac,
            neugc,
        })
        .collect()
}

pub const N_CORE_Y_LADDER: [(u8, u8); 5] = [(1, 0), (2, 0), (2, 1), (2, 2), (2, 3)];

pub const O_CORE_Y_LADDER: [(u8, u8); 1] = [(1, 0)];

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum GlycanCore {
    #[default]
    NGlycan,
    OGlycan,
}

impl GlycanCore {
    pub fn core_composition(self) -> (u8, u8) {
        match self {
            GlycanCore::NGlycan => (2, 3),
            GlycanCore::OGlycan => (1, 0),
        }
    }

    pub fn core_y_ladder(self) -> &'static [(u8, u8)] {
        match self {
            GlycanCore::NGlycan => &N_CORE_Y_LADDER,
            GlycanCore::OGlycan => &O_CORE_Y_LADDER,
        }
    }

    pub fn core_fucose(self) -> bool {
        matches!(self, GlycanCore::NGlycan)
    }

    pub fn default_min_core_ions(self) -> u32 {
        match self {
            GlycanCore::NGlycan => 2,
            GlycanCore::OGlycan => 1,
        }
    }
}

fn generate_y_ions(glycan: &Glycan, core: GlycanCore) -> Vec<GlycanIon> {
    let c = &glycan.composition;
    let (core_hn, core_h) = core.core_composition();
    let core_fuc: u8 = if core.core_fucose() {
        c.core_fuc.min(1)
    } else {
        0
    };
    let core_mass = core_hn as f64 * HEXNAC + core_h as f64 * HEX + core_fuc as f64 * FUC;

    let retained = c.tree.y_ion_retained_masses();
    let mut ions = Vec::with_capacity(retained.len() + 1);
    for r in retained {
        let lost = c.mass - r;
        if lost <= 1e-6 {
            continue;
        }
        ions.push(GlycanIon {
            glycan_id: glycan.id,
            mass: lost,
            is_core: r <= core_mass + 1e-6,
        });
    }

    ions.push(GlycanIon {
        glycan_id: glycan.id,
        mass: c.mass,
        is_core: false,
    });

    ions
}

#[inline]
pub fn encode_glycan_record(glycan_id: u32, is_core: bool) -> u32 {
    debug_assert!(
        glycan_id < (1 << 31),
        "glycan_id {glycan_id} exceeds 31 bits"
    );
    (glycan_id << 1) | (is_core as u32)
}

#[inline]
pub fn decode_glycan_record(encoded: u32) -> (u32, bool) {
    (encoded >> 1, (encoded & 1) != 0)
}

#[inline]
pub fn mass_to_bin(mass: f64, bin_width: f64) -> i64 {
    (mass / bin_width).round() as i64
}

pub struct GlycanIonIndex {
    bin_width: f64,
    tol_ppm: Option<f64>,
    core: GlycanCore,
    bins: HashMap<i64, Vec<(f64, u32)>>,
}

pub const PPM_FLOOR_DA: f64 = 0.005;

impl GlycanIonIndex {
    pub fn new(bin_width: f64, tol_ppm: Option<f64>, core: GlycanCore) -> Self {
        GlycanIonIndex {
            bin_width: bin_width.max(1e-4),
            tol_ppm,
            core,
            bins: HashMap::new(),
        }
    }

    pub fn len(&self) -> usize {
        self.bins.len()
    }

    pub fn ion_count(&self) -> usize {
        self.bins.values().map(|v| v.len()).sum()
    }

    pub fn is_empty(&self) -> bool {
        self.bins.is_empty()
    }

    pub fn add_glycan(&mut self, glycan: &Glycan) {
        assert!(
            glycan.id < (1 << 31),
            "glycan id {} exceeds 31 bits",
            glycan.id
        );
        for ion in generate_y_ions(glycan, self.core) {
            let key = mass_to_bin(ion.mass, self.bin_width);
            self.bins
                .entry(key)
                .or_default()
                .push((ion.mass, encode_glycan_record(ion.glycan_id, ion.is_core)));
        }
    }

    pub fn build(
        glycans: &[Glycan],
        bin_width: f64,
        tol_ppm: Option<f64>,
        core: GlycanCore,
    ) -> Self {
        let mut idx = Self::new(bin_width, tol_ppm, core);
        for g in glycans {
            idx.add_glycan(g);
        }
        idx
    }

    fn window_da(&self, mass: f64) -> f64 {
        match self.tol_ppm {
            Some(ppm) => (mass * ppm * 1e-6).max(PPM_FLOOR_DA),
            None => self.bin_width,
        }
    }

    pub fn query(&self, query_mass: f64) -> Vec<u32> {
        let tol = self.window_da(query_mass);
        let lo = mass_to_bin(query_mass - tol, self.bin_width);
        let hi = mass_to_bin(query_mass + tol, self.bin_width);
        let mut out = Vec::new();
        for b in lo..=hi {
            if let Some(records) = self.bins.get(&b) {
                for &(mass, enc) in records {
                    if (mass - query_mass).abs() <= tol {
                        out.push(enc);
                    }
                }
            }
        }
        out
    }
}

#[derive(Debug, Clone)]
pub struct GlycanCandidate {
    pub glycan_id: u32,
    pub score: f32,
    pub matched_y_ions: u32,
    pub matched_core_ions: u32,
    pub matched_diagnostic_ions: u32,
    pub mass_errors: Vec<f64>,
}

#[derive(Debug, Clone)]
pub struct GlycanConfig {
    pub tol_ppm: f64,
    pub min_diagnostic_ions: u32,
    pub min_core_ions: u32,
    pub top_k: usize,
    pub diagnostic_ions: Vec<DiagnosticIon>,
    pub diagnostic_filter: bool,
    pub score_weights: (f32, f32, f32),
}

impl Default for GlycanConfig {
    fn default() -> Self {
        GlycanConfig {
            tol_ppm: 20.0,
            min_diagnostic_ions: 2,
            min_core_ions: 2,
            top_k: 100,
            diagnostic_ions: default_diagnostic_ions(),
            diagnostic_filter: true,
            score_weights: (1.0, 1.0, 0.5),
        }
    }
}

#[derive(Debug, Clone)]
pub struct GlycanSearchResult {
    pub is_glycopeptide: bool,
    pub candidates: Vec<GlycanCandidate>,
}

pub fn search_glycans(
    index: &GlycanIonIndex,
    config: &GlycanConfig,
    glycans: &[GlycanComp],
    precursor_mass: f64,
    precursor_charge: u8,
    peaks: &[(f64, f32)],
) -> GlycanSearchResult {
    let matched_diagnostic: Vec<&DiagnosticIon> = config
        .diagnostic_ions
        .iter()
        .filter(|di| has_peak(peaks, di.mz, config.tol_ppm))
        .collect();
    if matched_diagnostic.len() < config.min_diagnostic_ions as usize {
        return GlycanSearchResult {
            is_glycopeptide: false,
            candidates: Vec::new(),
        };
    }

    let mut y_count: HashMap<u32, u32> = HashMap::new();
    let mut core_count: HashMap<u32, u32> = HashMap::new();
    let zmax = precursor_charge.max(1);
    for &(mz, _intensity) in peaks {
        for z in 1..=zmax {
            let neutral = (mz - PROTON) * z as f64;
            let q = precursor_mass - neutral;
            if q <= 0.0 {
                continue;
            }
            for enc in index.query(q) {
                let (gid, is_core) = decode_glycan_record(enc);
                *y_count.entry(gid).or_insert(0) += 1;
                if is_core {
                    *core_count.entry(gid).or_insert(0) += 1;
                }
            }
        }
    }

    let mut gids: Vec<u32> = y_count.keys().copied().collect();
    gids.sort_unstable();
    let mut candidates = Vec::with_capacity(gids.len().min(config.top_k));
    for gid in gids {
        let y = y_count[&gid];
        let c = core_count.get(&gid).copied().unwrap_or(0);
        if c < config.min_core_ions {
            continue;
        }
        let comp = match glycans.get(gid as usize) {
            Some(comp) => comp,
            None => continue,
        };
        if config.diagnostic_filter && !diagnostic_supported(comp, &matched_diagnostic) {
            continue;
        }
        let d = matched_diagnostic
            .iter()
            .filter(|di| di.consistent_with(comp))
            .count() as u32;
        let (w_y, w_c, w_d) = config.score_weights;
        let score = w_y * y as f32 + w_c * c as f32 + w_d * d as f32;
        candidates.push(GlycanCandidate {
            glycan_id: gid,
            score,
            matched_y_ions: y,
            matched_core_ions: c,
            matched_diagnostic_ions: d,
            mass_errors: Vec::new(),
        });
    }
    candidates.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap_or(Ordering::Equal));
    candidates.truncate(config.top_k);

    GlycanSearchResult {
        is_glycopeptide: true,
        candidates,
    }
}

fn diagnostic_supported(comp: &GlycanComp, matched: &[&DiagnosticIon]) -> bool {
    if comp.fuc > 0 && !matched.iter().any(|di| di.fuc > 0) {
        return false;
    }
    if comp.neuac > 0 && !matched.iter().any(|di| di.neuac > 0) {
        return false;
    }
    if comp.neugc > 0 && !matched.iter().any(|di| di.neugc > 0) {
        return false;
    }
    true
}

fn has_peak(peaks: &[(f64, f32)], target_mz: f64, tol_ppm: f64) -> bool {
    let tol = (target_mz * tol_ppm * 1e-6).max(PPM_FLOOR_DA);
    let idx = peaks.partition_point(|&(mz, _)| mz < target_mz - tol);
    idx < peaks.len() && peaks[idx].0 <= target_mz + tol
}

#[cfg(test)]
mod tests {
    use super::*;

    fn glycan(id: u32, hex: u8, hexnac: u8, fuc: u8, neuac: u8, neugc: u8) -> Glycan {
        glycan_with_core_fuc(id, hex, hexnac, fuc, neuac, neugc, 0)
    }

    fn glycan_core_fuc(id: u32, hex: u8, hexnac: u8, fuc: u8, neuac: u8, neugc: u8) -> Glycan {
        glycan_with_core_fuc(id, hex, hexnac, fuc, neuac, neugc, 1)
    }

    fn glycan_with_core_fuc(
        id: u32,
        hex: u8,
        hexnac: u8,
        fuc: u8,
        neuac: u8,
        neugc: u8,
        core_fuc: u8,
    ) -> Glycan {
        Glycan {
            id,
            composition: GlycanComp::from_counts(hexnac, hex, fuc, neuac, neugc, core_fuc),
        }
    }

    #[test]
    fn residue_masses_match_table() {
        assert_eq!(Monosaccharide::Hex.residue_mass(), HEX);
        assert_eq!(Monosaccharide::HexNAc.residue_mass(), HEXNAC);
        assert_eq!(Monosaccharide::Fuc.residue_mass(), FUC);
        assert_eq!(Monosaccharide::NeuAc.residue_mass(), NEUAC);
    }

    #[test]
    fn diagnostic_panel_matches_cation_formulas() {
        const C: f64 = 12.0;
        const H: f64 = 1.007_825_032_07;
        const N: f64 = 14.003_074_004_8;
        const O: f64 = 15.994_914_619_56;
        const E: f64 = 0.000_548_579_909_46;
        let mz = |c: u32, h: u32, n: u32, o: u32| {
            c as f64 * C + h as f64 * H + n as f64 * N + o as f64 * O - E
        };
        let formulas: [(u32, u32, u32, u32); 28] = [
            (4, 6, 1, 1),
            (4, 5, 0, 2),
            (5, 5, 0, 2),
            (6, 5, 0, 2),
            (5, 7, 0, 3),
            (6, 8, 1, 2),
            (6, 7, 0, 3),
            (7, 8, 1, 2),
            (6, 10, 1, 3),
            (6, 9, 0, 4),
            (6, 11, 0, 5),
            (8, 10, 1, 3),
            (8, 12, 1, 4),
            (8, 14, 1, 5),
            (11, 16, 1, 7),
            (11, 16, 1, 8),
            (11, 18, 1, 8),
            (11, 18, 1, 9),
            (12, 21, 0, 10),
            (14, 24, 1, 10),
            (17, 28, 1, 13),
            (20, 34, 1, 14),
            (20, 34, 1, 15),
            (25, 41, 2, 18),
            (26, 44, 1, 20),
            (31, 51, 2, 22),
            (34, 57, 2, 25),
            (40, 67, 2, 29),
        ];
        let panel = default_diagnostic_ions();
        assert_eq!(panel.len(), formulas.len());
        for (ion, &(c, h, n, o)) in panel.iter().zip(formulas.iter()) {
            let want = mz(c, h, n, o);
            let ppm = (ion.mz - want) / want * 1e6;
            assert!(
                ppm.abs() < 1.0,
                "diagnostic ion {:.5} vs formula C{c}H{h}N{n}O{o} {want:.5}: {ppm:.2} ppm",
                ion.mz
            );
        }
    }

    #[test]
    fn encode_decode_roundtrip() {
        for (id, core) in [(0, false), (1, true), (12345, true), ((1 << 31) - 1, false)] {
            let enc = encode_glycan_record(id, core);
            assert_eq!(decode_glycan_record(enc), (id, core));
        }
    }

    #[test]
    fn mass_binning_is_deterministic() {
        assert_eq!(mass_to_bin(1000.0, 0.1), mass_to_bin(1000.0, 0.1));
        assert_eq!(mass_to_bin(203.08, 0.02), 10154);
    }

    #[test]
    fn generate_y_ions_respects_core_ladder() {
        let g = glycan(0, 3, 2, 0, 0, 0);
        let ions = generate_y_ions(&g, GlycanCore::NGlycan);
        assert_eq!(ions.len(), 5);
        let y1_comp = 3.0 * HEX + HEXNAC;
        let y1 = ions
            .iter()
            .find(|i| (i.mass - y1_comp).abs() < 1e-6)
            .expect("Y1 complement present");
        assert!(y1.is_core);
        let y0 = ions
            .iter()
            .find(|i| (i.mass - g.composition.mass).abs() < 1e-6)
            .expect("full-glycan loss present");
        assert!(!y0.is_core);
    }

    #[test]
    fn generate_y_ions_is_structure_aware_not_full_power_set() {
        let g = glycan(0, 5, 4, 0, 0, 0);
        let ions = generate_y_ions(&g, GlycanCore::NGlycan);
        assert_eq!(ions.len(), 9);
        for ion in &ions {
            if !ion.is_core {
                let remaining = g.composition.mass - ion.mass;
                let is_y0 = remaining < 1e-6;
                let keeps_core = remaining >= 3.0 * HEX + 2.0 * HEXNAC - 1e-6;
                assert!(is_y0 || keeps_core);
            }
        }
    }

    #[test]
    fn o_glycan_core_ion_is_single_galnac() {
        let g = glycan(0, 1, 1, 0, 0, 0);
        let ions = generate_y_ions(&g, GlycanCore::OGlycan);
        assert_eq!(ions.len(), 2);
        let y1 = ions
            .iter()
            .find(|i| (i.mass - HEX).abs() < 1e-6)
            .expect("Y1 complement (Gal) present");
        assert!(y1.is_core);
        let y0 = ions
            .iter()
            .find(|i| (i.mass - g.composition.mass).abs() < 1e-6)
            .expect("Y0 present");
        assert!(!y0.is_core);
    }

    #[test]
    fn n_glycan_core_fucose_is_a_separate_branch_of_the_reducing_end() {
        let g = glycan_core_fuc(0, 3, 2, 1, 0, 0);
        let ions = generate_y_ions(&g, GlycanCore::NGlycan);
        assert_eq!(ions.len(), 6);
        assert!(
            ions.iter().any(|i| (i.mass - FUC).abs() < 1e-6),
            "core Fuc must be losable alone (separate branch)"
        );
        let y1 = ions
            .iter()
            .find(|i| (i.mass - (3.0 * HEX + HEXNAC)).abs() < 1e-6)
            .expect("Y1+Fuc complement present");
        assert!(y1.is_core);
        for ion in &ions {
            if !ion.is_core {
                assert!(
                    (ion.mass - g.composition.mass).abs() < 1e-6,
                    "only Y0 is non-core"
                );
            }
        }
    }

    #[test]
    fn antenna_fucose_stays_antenna_not_pinned_to_core() {
        let g = glycan(0, 3, 2, 1, 0, 0);
        let ions = generate_y_ions(&g, GlycanCore::NGlycan);
        assert!(
            ions.iter().any(|i| (i.mass - FUC).abs() < 1e-6),
            "antenna Fuc must be lost alone"
        );
        assert!(
            !ions
                .iter()
                .any(|i| (i.mass - (3.0 * HEX + HEXNAC)).abs() < 1e-6),
            "no core-fucose Y1 signature expected for antenna fucose"
        );
    }

    #[test]
    fn index_retrieves_glycan_by_complementary_ion() {
        let glycans = vec![glycan(7, 3, 2, 0, 0, 0)];
        let idx = GlycanIonIndex::build(&glycans, 0.02, Some(20.0), GlycanCore::NGlycan);
        let q = 3.0 * HEX + HEXNAC;
        let hits = idx.query(q);
        assert!(
            hits.iter()
                .any(|&enc| decode_glycan_record(enc) == (7, true)),
            "Y1 complement should retrieve glycan 7 as core"
        );
    }

    #[test]
    fn search_returns_candidates_for_glycopeptide_spectrum() {
        let glycans = vec![glycan(0, 5, 2, 0, 0, 0), glycan(1, 3, 2, 0, 0, 0)];
        let idx = GlycanIonIndex::build(&glycans, 0.02, Some(20.0), GlycanCore::NGlycan);
        let cfg = GlycanConfig::default();
        let peptide = 1000.0;
        let precursor_mass = peptide + 3.0 * HEX + 2.0 * HEXNAC;
        let mut peaks = vec![
            (204.0867, 100.0),
            (366.1395, 80.0),
            (138.0550, 50.0),
            (peptide + HEXNAC + PROTON, 60.0),
            (peptide + 2.0 * HEXNAC + PROTON, 50.0),
        ];
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let comps: Vec<GlycanComp> = glycans.iter().map(|g| g.composition.clone()).collect();
        let result = search_glycans(&idx, &cfg, &comps, precursor_mass, 3, &peaks);
        assert!(result.is_glycopeptide);
        assert!(result.candidates.iter().any(|c| c.glycan_id == 1));
        assert!(result.candidates.iter().all(|c| c.glycan_id != 0));
    }

    #[test]
    fn search_rejects_non_glycopeptide() {
        let glycans = vec![glycan(1, 3, 2, 0, 0, 0)];
        let idx = GlycanIonIndex::build(&glycans, 0.02, Some(20.0), GlycanCore::NGlycan);
        let cfg = GlycanConfig::default();
        let peaks = vec![(100.0, 10.0), (300.0, 10.0), (500.0, 10.0)];
        let comps: Vec<GlycanComp> = glycans.iter().map(|g| g.composition.clone()).collect();
        let result = search_glycans(&idx, &cfg, &comps, 1500.0, 2, &peaks);
        assert!(!result.is_glycopeptide);
    }

    #[test]
    fn diagnostic_filter_rejects_unsupported_sialic_acid() {
        let glycans = vec![glycan(0, 3, 2, 0, 1, 0)];
        let idx = GlycanIonIndex::build(&glycans, 0.02, Some(20.0), GlycanCore::NGlycan);
        let cfg = GlycanConfig {
            diagnostic_filter: true,
            ..GlycanConfig::default()
        };
        let peptide = 1000.0;
        let precursor_mass = peptide + 3.0 * HEX + 2.0 * HEXNAC + NEUAC;
        let mut peaks = vec![
            (204.0867, 100.0),
            (366.1395, 80.0),
            (peptide + HEXNAC + PROTON, 60.0),
            (peptide + 2.0 * HEXNAC + PROTON, 50.0),
        ];
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let comps: Vec<GlycanComp> = glycans.iter().map(|g| g.composition.clone()).collect();
        let result = search_glycans(&idx, &cfg, &comps, precursor_mass, 3, &peaks);
        assert!(
            result.candidates.is_empty(),
            "unsupported NeuAc must be filtered"
        );
    }

    #[test]
    fn benchmark_index_construction_and_query() {
        let comps = match std::env::var("ANDES_GLYCO_GDB") {
            Ok(path) => {
                let content =
                    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("read {path}: {e}"));
                crate::glycan_db::load_glycan_gdb(&content)
                    .unwrap_or_else(|e| panic!("parse {path}: {e}"))
            }
            Err(_) => {
                eprintln!("skipping: ANDES_GLYCO_GDB not set");
                return;
            }
        };
        let glycans: Vec<Glycan> = comps
            .iter()
            .enumerate()
            .map(|(i, c)| Glycan {
                id: i as u32,
                composition: c.clone(),
            })
            .collect();

        let t0 = std::time::Instant::now();
        let idx = GlycanIonIndex::build(&glycans, 0.02, Some(20.0), GlycanCore::NGlycan);
        let build_time = t0.elapsed();

        let peaks: Vec<(f64, f32)> = (0..200).map(|i| (200.0 + i as f64 * 5.0, 100.0)).collect();
        let t1 = std::time::Instant::now();
        let mut total_records = 0usize;
        for &(mz, _) in &peaks {
            total_records += idx.query(2500.0 - (mz - PROTON)).len();
        }
        let query_time = t1.elapsed();

        println!(
            "glycan-first index: {} glycans -> {} ions / {} bins; build {:?}, 200-peak query {:?} ({} records)",
            glycans.len(),
            idx.ion_count(),
            idx.len(),
            build_time,
            query_time,
            total_records
        );
        assert!(idx.ion_count() > 0, "index must contain ions");
    }
}
