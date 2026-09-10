//! Glycan-first candidate retrieval (pGlyco3-style "glycan ion indexing").
//!
//! This is the search-strategy counterpart to the backbone-first driver: instead
//! of enumerating `precursor − glycan` for every glycan and then scoring the
//! backbone, we index each glycan's *Y-complementary fragment ions* once and,
//! per spectrum, use the MS/MS peaks via `precursor − peak` to retrieve candidate
//! glycans in O(peaks × charge × neighbourhood bins), independent of the
//! glycan-database size.
//!
//! Pipeline (matches the spec's "main part"):
//!
//! ```text
//! glycan database → Y-complementary ion generation → ion index → peak query
//!     → candidate retrieval → core-ion filter → diagnostic-ion filter → top-K
//! ```
//!
//! # Mass conventions (the one place to get right)
//!
//! Everything here is in **residue convention** for the glycan (Σ of residue
//! masses, no extra +H₂O), but the precursor passed to [`search_glycans`] is the
//! **neutral monoisotopic mass** `M = peptide_neutral + Σ(glycan residues)`, i.e.
//! `(m/z − H⁺)·z` *without* the −H₂O that the backbone driver subtracts. This is
//! what makes the complementary query exact:
//!
//! ```text
//! Y-ion peak (charge z):   m/z = (peptide_neutral + y_comp + z·H⁺) / z
//! neutral at charge z:     (m/z − H⁺)·z = peptide_neutral + y_comp
//! query mass:              M − (peptide_neutral + y_comp) = glycan − y_comp
//! ```
//!
//! `glycan − y_comp` is exactly the Y-complementary ion mass the index stores.
//! Subtracting H₂O from the precursor (residue convention) would shift every
//! query by −18.0106 Da and make the index miss everything.
//!
//! Masses and thresholds that the reference figure does not pin down (exact
//! Y-ion formula, core definitions, diagnostic-ion list, scoring weights) are
//! exposed in [`GlycanConfig`], [`GlycanCore`] and [`DiagnosticIon`] rather than
//! silently assumed.

use std::cmp::Ordering;
use std::collections::HashMap;

use crate::glycan_db::GlycanComp;
use crate::glycan_mass::{FUC, HEX, HEXNAC, NEUAC, NEUGC, PROTON};

/// A monosaccharide counted in a glycan composition. Extensible: add a variant
/// here, its residue mass in [`Monosaccharide::residue_mass`], and (if it has a
/// diagnostic ion) an entry in the configured [`DiagnosticIon`] table.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Monosaccharide {
    Hex,
    HexNAc,
    Fuc,
    NeuAc,
    NeuGc,
}

impl Monosaccharide {
    /// Residue mass (free monosaccharide − H₂O): the mass contributed when the
    /// monosaccharide is polymerised into a glycan. Single source of truth is
    /// [`crate::glycan_mass`]; this is the enum-indexed view.
    pub fn residue_mass(self) -> f64 {
        match self {
            Monosaccharide::Hex => HEX,
            Monosaccharide::HexNAc => HEXNAC,
            Monosaccharide::Fuc => FUC,
            Monosaccharide::NeuAc => NEUAC,
            Monosaccharide::NeuGc => NEUGC,
        }
    }

    /// Neutral (free, reducing-end) monosaccharide mass = residue + H₂O.
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

/// A glycan: a unique id plus its monosaccharide composition.
#[derive(Debug, Clone)]
pub struct Glycan {
    pub id: u32,
    pub composition: GlycanComp,
}

/// One Y-complementary fragment ion of a glycan: the glycan part *lost* when a
/// Y ion fragments, i.e. `full glycan − Y ion`. Indexed at `mass`; `is_core`
/// records whether the *corresponding Y ion* was a trimannosyl-core Y ion.
#[derive(Debug, Clone)]
struct GlycanIon {
    glycan_id: u32,
    /// Y-complementary mass (residue convention, peptide-independent).
    mass: f64,
    is_core: bool,
}

/// A diagnostic (oxonium / B-type) ion used to decide whether a spectrum is a
/// glycopeptide and which monosaccharides its glycan carries. `counts` is the
/// (hex, hexnac, fuc, neuac, neugc) the ion implies; an ion is "consistent" with
/// a composition only if every implied count is present in the composition.
#[derive(Debug, Clone, Copy)]
pub struct DiagnosticIon {
    /// Ion m/z (proton already added).
    pub mz: f64,
    pub hex: u8,
    pub hexnac: u8,
    pub fuc: u8,
    pub neuac: u8,
    pub neugc: u8,
}

impl DiagnosticIon {
    /// True when `comp` contains every monosaccharide this ion implies.
    pub fn consistent_with(&self, comp: &GlycanComp) -> bool {
        comp.hex >= self.hex
            && comp.hexnac >= self.hexnac
            && comp.fuc >= self.fuc
            && comp.neuac >= self.neuac
            && comp.neugc >= self.neugc
    }
}

/// The default N-glycan diagnostic-ion panel: the 28 oxonium (B-fragment) ions
/// from Glyco-Decipher's `OxoniumIonEnum`. `counts` are the monosaccharides each
/// ion implies (a fragment implies ≥ that many of each monosaccharide).
/// Configurable via [`GlycanConfig::diagnostic_ions`].
pub fn default_diagnostic_ions() -> Vec<DiagnosticIon> {
    // (mz, hex, hexnac, fuc, neuac, neugc)
    const IONS: [(f64, u8, u8, u8, u8, u8); 28] = [
        (84.0444, 0, 1, 0, 0, 0),  // HexNAc fragment
        (85.0284, 1, 0, 0, 0, 0),  // Hex fragment
        (97.0284, 1, 0, 0, 0, 0),  // Hex fragment
        (109.0270, 1, 0, 0, 0, 0), // Hex fragment
        (115.0390, 1, 0, 0, 0, 0), // Hex fragment
        (126.0550, 0, 1, 0, 0, 0), // HexNAc fragment
        (127.0390, 1, 0, 0, 0, 0), // Hex fragment
        (138.0550, 0, 1, 0, 0, 0), // HexNAc fragment
        (144.0640, 0, 1, 0, 0, 0), // HexNAc fragment
        (145.0670, 1, 0, 0, 0, 0), // Hex fragment
        (163.0600, 1, 0, 0, 0, 0), // Hex
        (168.0660, 0, 1, 0, 0, 0), // HexNAc fragment
        (186.0760, 0, 1, 0, 0, 0), // HexNAc fragment
        (204.0867, 0, 1, 0, 0, 0), // HexNAc
        (274.0874, 0, 0, 0, 1, 0), // NeuAc−H2O
        (290.0823, 0, 0, 0, 0, 1), // NeuGc−H2O
        (292.1027, 0, 0, 0, 1, 0), // NeuAc
        (308.0976, 0, 0, 0, 0, 1), // NeuGc
        (325.1129, 2, 0, 0, 0, 0), // Hex2
        (366.1395, 1, 1, 0, 0, 0), // HexHexNAc
        (454.1560, 1, 0, 0, 1, 0), // HexNeuAc
        (512.1970, 1, 1, 1, 0, 0), // HexHexNAcFuc
        (528.1917, 2, 1, 0, 0, 0), // Hex2HexNAc
        (657.2349, 1, 1, 0, 1, 0), // HexHexNAcNeuAc
        (690.2445, 3, 1, 0, 0, 0), // Hex3HexNAc
        (803.2930, 1, 1, 1, 1, 0), // HexHexNAcNeuAcFuc
        (893.3239, 3, 2, 0, 0, 0), // Hex3HexNAc2
        (1039.3823, 3, 2, 1, 0, 0), // Hex3HexNAc2Fuc
    ];
    IONS
        .iter()
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

/// Core Y-ion ladder for N-glycans: `(hexnac, hex)` compositions that stay on the
/// peptide, innermost (Y1) first. These are the trimannosyl-core ladder Y1..Y5.
pub const N_CORE_Y_LADDER: [(u8, u8); 5] = [
    (1, 0), // Y1 = HexNAc
    (2, 0), // Y2 = 2HexNAc
    (2, 1), // Y3 = 2HexNAc+Hex
    (2, 2), // Y4 = 2HexNAc+2Hex
    (2, 3), // Y5 = 2HexNAc+3Hex (full trimannosyl core)
];

/// Core Y-ion ladder for O-glycans: only the reducing-end GalNAc (Y1). O-glycan
/// cores 1–4 (and beyond) all share the GalNAc → Ser/Thr linkage, so a single
/// GalNAc is the only "core" rung; every other residue is antenna.
pub const O_CORE_Y_LADDER: [(u8, u8); 1] = [(1, 0)]; // Y1 = GalNAc (HexNAc)

/// The glycan core (reducing-end) type. Determines which Y ions count as *core*
/// ions, whether core-fucose is recognized, and the default core-ion filter
/// threshold. A single [`GlycanIonIndex`] indexes one core type at a time
/// (pGlyco3 likewise ships separate N- and O-glycan databases).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum GlycanCore {
    /// Canonical N-glycan trimannosyl core (GlcNAc₂Man₃).
    #[default]
    NGlycan,
    /// O-glycan core: the reducing-end GalNAc (HexNAc).
    OGlycan,
}

impl GlycanCore {
    /// Core composition `(hexnac, hex)` retained at the reducing end.
    pub fn core_composition(self) -> (u8, u8) {
        match self {
            GlycanCore::NGlycan => (2, 3),
            GlycanCore::OGlycan => (1, 0),
        }
    }

    /// Core Y-ion ladder, Y1 first, ending at the full core.
    pub fn core_y_ladder(self) -> &'static [(u8, u8)] {
        match self {
            GlycanCore::NGlycan => &N_CORE_Y_LADDER,
            GlycanCore::OGlycan => &O_CORE_Y_LADDER,
        }
    }

    /// Whether one Fuc may be core-fucose (α1-6 on the reducing-end GlcNAc of an
    /// N-glycan). When true and the composition carries ≥1 Fuc, one Fuc is pinned
    /// to the core and retained on every core Y ion instead of being an antenna
    /// residue. (Any Fuc beyond the first stays antenna/Lewis fucose.)
    pub fn core_fucose(self) -> bool {
        matches!(self, GlycanCore::NGlycan)
    }

    /// Default minimum matched core Y ions for the glycan-first core filter:
    /// 2 for N-glycans (Y-ladder), 1 for O-glycans (single GalNAc).
    pub fn default_min_core_ions(self) -> u32 {
        match self {
            GlycanCore::NGlycan => 2,
            GlycanCore::OGlycan => 1,
        }
    }
}

/// Generate every Y-complementary ion of one glycan, structure-aware.
///
/// A Y ion is a fragment that stays on the peptide and always contains the
/// reducing end. The naive analogue "every residue sub-composition" over-generates
/// physically impossible fragments (trimming a core mannose while keeping an
/// antenna, or splitting the core). We instead respect the core:
///
///   1. **Antenna losses** — every sub-composition of the antenna residues
///      (`full − core − core_fucose`). The remaining Y ion is `core + (antenna −
///      lost)` and is a *core* ion only when `lost` is the whole antenna (remaining
///      == the full core, the outermost rung).
///   2. **Core ladder** — the antenna fully removed, then core residues trimmed in
///      the fixed Yₘₐₓ→…→Y1 order; each rung's `full − rung` complement is flagged
///      `is_core`.
///   3. **Y0** — the bare peptide (`full` complement), never `is_core`.
///
/// The empty loss (mass 0 = the intact precursor) is skipped.
fn generate_y_ions(glycan: &Glycan, core: GlycanCore) -> Vec<GlycanIon> {
    let c = &glycan.composition;
    let (core_hn, core_h) = core.core_composition();
    // Core fucose consumes at most one Fuc; the rest is antenna fucose. Read the
    // structural fact from the glycan (recorded by the .gdb loader) instead of
    // guessing "first Fuc is core"; the guess is wrong for antenna/Lewis fucose.
    let core_fuc: u8 = if core.core_fucose() { c.core_fuc.min(1) } else { 0 };

    let ant_hn = c.hexnac.saturating_sub(core_hn);
    let ant_h = c.hex.saturating_sub(core_h);
    let ant_f = c.fuc - core_fuc;
    let ant_a = c.neuac;
    let ant_g = c.neugc;

    let full_mass = c.mass;
    let gid = glycan.id;

    let cap = (ant_hn as usize + 1)
        * (ant_h as usize + 1)
        * (ant_f as usize + 1)
        * (ant_a as usize + 1)
        * (ant_g as usize + 1);
    let mut ions = Vec::with_capacity(cap + core.core_y_ladder().len());

    // 1. Antenna losses (non-empty sub-compositions of the antenna).
    for lh in 0..=ant_h {
        for lhn in 0..=ant_hn {
            for lf in 0..=ant_f {
                for la in 0..=ant_a {
                    for lg in 0..=ant_g {
                        if lh == 0 && lhn == 0 && lf == 0 && la == 0 && lg == 0 {
                            continue; // empty loss = intact precursor
                        }
                        let lost = lh as f64 * HEX
                            + lhn as f64 * HEXNAC
                            + lf as f64 * FUC
                            + la as f64 * NEUAC
                            + lg as f64 * NEUGC;
                        // Remaining = core + core_fuc + (antenna − lost). It is a core
                        // ion only when the whole antenna is lost (remaining == full core).
                        let is_core = lh == ant_h
                            && lhn == ant_hn
                            && lf == ant_f
                            && la == ant_a
                            && lg == ant_g;
                        ions.push(GlycanIon {
                            glycan_id: gid,
                            mass: lost,
                            is_core,
                        });
                    }
                }
            }
        }
    }

    // 2. Core ladder below the full core (Yₘₐₓ−1 … Y1), antenna already removed.
    for &(hn, h) in core.core_y_ladder() {
        if (hn, h) == (core_hn, core_h) {
            continue; // full core already emitted as "whole antenna lost" above
        }
        if hn > c.hexnac || h > c.hex {
            continue; // rung unreachable from a truncated composition
        }
        let lost = full_mass - (hn as f64 * HEXNAC + h as f64 * HEX + core_fuc as f64 * FUC);
        ions.push(GlycanIon {
            glycan_id: gid,
            mass: lost,
            is_core: true,
        });
    }

    // 3. Y0: the bare peptide (everything lost).
    ions.push(GlycanIon {
        glycan_id: gid,
        mass: full_mass,
        is_core: false,
    });

    ions
}

/// Compact index record: `glycan_id` in the high 31 bits, `is_core` in the low
/// bit. Glycan ids must fit in 31 bits (asserted at index build).
#[inline]
pub fn encode_glycan_record(glycan_id: u32, is_core: bool) -> u32 {
    debug_assert!(glycan_id < (1 << 31), "glycan_id {glycan_id} exceeds 31 bits");
    (glycan_id << 1) | (is_core as u32)
}

/// Inverse of [`encode_glycan_record`].
#[inline]
pub fn decode_glycan_record(encoded: u32) -> (u32, bool) {
    (encoded >> 1, (encoded & 1) != 0)
}

/// Mass binning: integer key for an approximate mass lookup.
#[inline]
pub fn mass_to_bin(mass: f64, bin_width: f64) -> i64 {
    (mass / bin_width).round() as i64
}

/// Hash index from (binned) Y-complementary ion mass to compact glycan records.
///
/// Construction is O(number of glycan ions); a query is O(peaks × neighbourhood
/// bins), independent of the number of glycans. Each record stores the exact ion
/// mass so a query can apply the true tolerance window (not just the bin).
pub struct GlycanIonIndex {
    bin_width: f64,
    /// ppm tolerance, when set; the acceptance window is `mass * ppm * 1e-6`
    /// (floored at [`PPM_FLOOR_DA`]). `None` means a fixed-Da window `bin_width`.
    tol_ppm: Option<f64>,
    /// Core (reducing-end) type the indexed glycans share; drives Y-ion generation.
    core: GlycanCore,
    /// bin → `(exact_ion_mass, encoded_record)`.
    bins: HashMap<i64, Vec<(f64, u32)>>,
}

/// Absolute floor (Da) on a ppm acceptance window, so a light fragment is not
/// held to an unrealistically tight tolerance.
pub const PPM_FLOOR_DA: f64 = 0.005;

impl GlycanIonIndex {
    /// Build an empty index.
    ///
    /// `bin_width` is a fixed Da bin for the integer keys (kept fine enough that
    /// the neighbourhood walk covers every admissible match; exact tolerance is
    /// checked per record at query time via `tol_ppm` or `bin_width`).
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

    /// Total number of ion records across all bins.
    pub fn ion_count(&self) -> usize {
        self.bins.values().map(|v| v.len()).sum()
    }

    pub fn is_empty(&self) -> bool {
        self.bins.is_empty()
    }

    /// Add one glycan's Y-complementary ions to the index.
    pub fn add_glycan(&mut self, glycan: &Glycan) {
        assert!(glycan.id < (1 << 31), "glycan id {} exceeds 31 bits", glycan.id);
        for ion in generate_y_ions(glycan, self.core) {
            let key = mass_to_bin(ion.mass, self.bin_width);
            self.bins
                .entry(key)
                .or_default()
                .push((ion.mass, encode_glycan_record(ion.glycan_id, ion.is_core)));
        }
    }

    /// Build an index over a whole database of glycans of one core type.
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

    /// Half-window (Da) for a query mass.
    fn window_da(&self, mass: f64) -> f64 {
        match self.tol_ppm {
            Some(ppm) => (mass * ppm * 1e-6).max(PPM_FLOOR_DA),
            None => self.bin_width,
        }
    }

    /// Query the index for `query_mass`, returning the encoded records whose ion
    /// mass is within tolerance. Each returned `u32` is a
    /// [`encode_glycan_record`] value; decode with [`decode_glycan_record`].
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

/// A candidate glycan with the evidence accumulated from a spectrum.
#[derive(Debug, Clone)]
pub struct GlycanCandidate {
    pub glycan_id: u32,
    pub score: f32,
    pub matched_y_ions: u32,
    pub matched_core_ions: u32,
    pub matched_diagnostic_ions: u32,
    pub mass_errors: Vec<f64>,
}

/// Tunable search parameters. Fields that the reference figure leaves open are
/// explicit here (and in [`DiagnosticIon`]) rather than hard-coded.
#[derive(Debug, Clone)]
pub struct GlycanConfig {
    /// Mass tolerance in ppm (used for both the index window and scoring).
    pub tol_ppm: f64,
    /// Minimum diagnostic (oxonium) ions required for a spectrum to be treated
    /// as a glycopeptide at all.
    pub min_diagnostic_ions: u32,
    /// Minimum matched core ions a candidate must have to survive core filtering.
    /// pGlyco3 uses 2 for N-glycans and 1 for O-glycans.
    pub min_core_ions: u32,
    /// Number of top candidates to return per spectrum. pGlyco3 keeps the top 100
    /// candidate glycan compositions for the peptide search.
    pub top_k: usize,
    /// Diagnostic-ion panel (defaults to [`default_diagnostic_ions`]).
    pub diagnostic_ions: Vec<DiagnosticIon>,
    /// When true, discard a candidate whose variable monosaccharides (Fuc, NeuAc,
    /// NeuGc) are not supported by a fired diagnostic ion. Hex/HexNAc are present
    /// in every N-glycan and are therefore not gated.
    pub diagnostic_filter: bool,
    /// Weights for the score: (y-ion, core-ion, diagnostic). Configurable so the
    /// scoring formula is not silently fixed.
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

/// The result of a glycan-first search over one spectrum.
#[derive(Debug, Clone)]
pub struct GlycanSearchResult {
    pub is_glycopeptide: bool,
    pub candidates: Vec<GlycanCandidate>,
}

/// Search the glycan ion index for a spectrum.
///
/// `peaks` is `(m/z, intensity)`, sorted by m/z ascending. `glycans` is the
/// glycan list the index was built from, indexed by `glycan_id`. `precursor_mass`
/// is the **neutral monoisotopic mass** `M = peptide_neutral + Σ(glycan residues)`
/// (`(m/z − H⁺)·z`, *without* the −H₂O the backbone driver applies); see the
/// module docs. `precursor_charge` bounds the fragment-charge ladder.
///
/// Pipeline:
///
/// 1. **Diagnostic precheck** — enough oxonium ions to call the spectrum a
///    glycopeptide.
/// 2. **Y-complementary retrieval** — for each peak and each charge,
///    `q = M − (m/z − H⁺)·z` is looked up in the index; matches accumulate
///    per-glycan total and core (via the index's `is_core` bit) counts.
/// 3. **Core-ion filter** — keep a glycan only with ≥ `min_core_ions` core ions.
/// 4. **Diagnostic-ion filter** — drop glycans whose variable monosaccharides
///    lack diagnostic-ion support.
/// 5. **Score + top-K** — `w_y·y + w_c·core + w_d·diag`, top `top_k` returned.
pub fn search_glycans(
    index: &GlycanIonIndex,
    config: &GlycanConfig,
    glycans: &[GlycanComp],
    precursor_mass: f64,
    precursor_charge: u8,
    peaks: &[(f64, f32)],
) -> GlycanSearchResult {
    // 1. Diagnostic precheck.
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

    // 2. Y-complementary retrieval. Each peak is interpreted at every charge
    // from 1 up to the precursor charge, so multiply-charged Y ions (frequent in
    // glycopeptide MS/MS) are not lost to a charge-1-only read.
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

    // 3–5. Core filter, diagnostic filter, score, top-K.
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
        is_glycopeptide: !candidates.is_empty(),
        candidates,
    }
}

/// Diagnostic-ion support for a candidate's variable monosaccharides: a glycan
/// claiming Fuc/NeuAc/NeuGc must have at least one fired diagnostic ion implying
/// that residue. Hex/HexNAc are in every N-glycan and are not gated (they would
/// be for O-glycans or other core definitions — a configuration point).
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

/// Binary search for a peak within a ppm window of `target_mz`.
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

    /// Like [`glycan`], but with `core_fuc` set on the composition (a Fuc that is
    /// a direct child of the reducing-end GlcNAc).
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
        let mass = hex as f64 * HEX
            + hexnac as f64 * HEXNAC
            + fuc as f64 * FUC
            + neuac as f64 * NEUAC
            + neugc as f64 * NEUGC;
        Glycan {
            id,
            composition: GlycanComp {
                hexnac,
                hex,
                fuc,
                neuac,
                neugc,
                core_fuc,
                mass,
            },
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
    fn encode_decode_roundtrip() {
        for (id, core) in [(0, false), (1, true), (12345, true), ((1 << 31) - 1, false)] {
            let enc = encode_glycan_record(id, core);
            assert_eq!(decode_glycan_record(enc), (id, core));
        }
    }

    #[test]
    fn mass_binning_is_deterministic() {
        assert_eq!(mass_to_bin(1000.0, 0.1), mass_to_bin(1000.0, 0.1));
        assert_eq!(mass_to_bin(203.08, 0.02), 10154); // round(203.08/0.02)
    }

    #[test]
    fn generate_y_ions_respects_core_ladder() {
        // Man3GlcNAc2 (hex=3, hexnac=2): no antenna, so the only Y ions are the
        // core ladder Y4..Y1 (Y5 == intact precursor, skipped) plus Y0 → 5 ions.
        let g = glycan(0, 3, 2, 0, 0, 0);
        let ions = generate_y_ions(&g, GlycanCore::NGlycan);
        assert_eq!(ions.len(), 5);
        // The Y1 core complement: lost = full − HexNAc = 3Hex+1HexNAc, is_core.
        let y1_comp = 3.0 * HEX + HEXNAC;
        let y1 = ions
            .iter()
            .find(|i| (i.mass - y1_comp).abs() < 1e-6)
            .expect("Y1 complement present");
        assert!(y1.is_core);
        // The full-glycan loss (Y0, bare peptide) is present and NOT core.
        let y0 = ions
            .iter()
            .find(|i| (i.mass - g.composition.mass).abs() < 1e-6)
            .expect("full-glycan loss present");
        assert!(!y0.is_core);
    }

    #[test]
    fn generate_y_ions_is_structure_aware_not_full_power_set() {
        // Hex5HexNAc4: the full power set would be 6·5 = 30 sub-compositions; the
        // structure-aware ladder emits 13 (8 antenna + 4 core rungs + Y0), and no
        // ion trims the core while keeping an antenna residue.
        let g = glycan(0, 5, 4, 0, 0, 0);
        let ions = generate_y_ions(&g, GlycanCore::NGlycan);
        assert_eq!(ions.len(), 13);
        // Every non-core ion (except Y0, the full loss) must leave the whole core
        // intact on the peptide.
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
        // Core-1 O-glycan Galβ1-3GalNAc = Hex1HexNAc1. Core = GalNAc(HexNAc), so
        // the only core Y ion is Y1 (remaining == GalNAc), whose complement is the
        // Gal residue (Hex); plus Y0.
        let g = glycan(0, 1, 1, 0, 0, 0);
        let ions = generate_y_ions(&g, GlycanCore::OGlycan);
        assert_eq!(ions.len(), 2); // {Gal lost → Y1 core}, {everything lost → Y0}
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
    fn n_glycan_core_fucose_pins_fuc_to_core_ions() {
        // Man3Fuc (Hex3HexNAc2Fuc1): the single Fuc is core-fucose, retained on
        // every core rung, so there are no antenna ions — only Y4..Y1 (+Fuc) and Y0.
        let g = glycan_core_fuc(0, 3, 2, 1, 0, 0);
        let ions = generate_y_ions(&g, GlycanCore::NGlycan);
        assert_eq!(ions.len(), 5);
        // Y1 complement (remaining = HexNAc + Fuc) = 3Hex + 1HexNAc, is_core.
        let y1 = ions
            .iter()
            .find(|i| (i.mass - (3.0 * HEX + HEXNAC)).abs() < 1e-6)
            .expect("Y1+Fuc complement present");
        assert!(y1.is_core);
        // The only non-core ion is Y0 (everything lost); the fucose is core-fucose
        // and never appears as an antenna loss.
        for ion in &ions {
            if !ion.is_core {
                assert!((ion.mass - g.composition.mass).abs() < 1e-6, "only Y0 is non-core");
            }
        }
    }

    #[test]
    fn antenna_fucose_stays_antenna_not_pinned_to_core() {
        // Same composition (Hex3HexNAc2Fuc1) as the core-fucose test, but with
        // core_fuc = 0 (a Fuc on an outer mannose, not the reducing-end GlcNAc).
        // The Fuc is an antenna residue: it is lost alone (mass == FUC), and the
        // core ladder does NOT carry +Fuc on every rung.
        let g = glycan(0, 3, 2, 1, 0, 0);
        let ions = generate_y_ions(&g, GlycanCore::NGlycan);
        // The antenna Fuc is lost as a standalone ion (its complement is the core).
        assert!(
            ions.iter().any(|i| (i.mass - FUC).abs() < 1e-6),
            "antenna Fuc must be lost alone"
        );
        // Unlike the core-fucose case, there is no ion at the "3Hex + HexNAc"
        // complement (that is the core-fucose Y1 rung signature, present only
        // when the Fuc rides the core).
        assert!(
            !ions.iter().any(|i| (i.mass - (3.0 * HEX + HEXNAC)).abs() < 1e-6),
            "no core-fucose Y1 signature expected for antenna fucose"
        );
    }

    #[test]
    fn index_retrieves_glycan_by_complementary_ion() {
        let glycans = vec![glycan(7, 3, 2, 0, 0, 0)];
        let idx = GlycanIonIndex::build(&glycans, 0.02, Some(20.0), GlycanCore::NGlycan);
        // Query the Y1 complement (full − HexNAc = 3Hex + 1HexNAc), which must
        // retrieve glycan 7 with the core flag set.
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
        // glycan ids must equal their index in `glycans` (the driver assigns ids
        // via enumerate), since `search_glycans` looks up the composition by id.
        let glycans = vec![
            glycan(0, 5, 2, 0, 0, 0), // Man5 (distractor)
            glycan(1, 3, 2, 0, 0, 0), // Man3 (truth)
        ];
        let idx = GlycanIonIndex::build(&glycans, 0.02, Some(20.0), GlycanCore::NGlycan);
        let cfg = GlycanConfig::default();
        // Man3 (Hex3HexNAc2) glycopeptide: peptide 1000, glycan 892.32, charge 3.
        // precursor_mass is the NEUTRAL mass M = peptide_neutral + glycan.
        let peptide = 1000.0;
        let precursor_mass = peptide + 3.0 * HEX + 2.0 * HEXNAC;
        let mut peaks = vec![
            (204.0867, 100.0), // HexNAc oxonium
            (366.1395, 80.0),  // HexHexNAc oxonium
            (138.0550, 50.0),  // HexNAc fragment
            // Core Y ions (peptide + core fragments), m/z at charge 1.
            (peptide + HEXNAC + PROTON, 60.0),              // Y1
            (peptide + 2.0 * HEXNAC + PROTON, 50.0),        // Y2
        ];
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let comps: Vec<GlycanComp> = glycans.iter().map(|g| g.composition.clone()).collect();
        let result = search_glycans(&idx, &cfg, &comps, precursor_mass, 3, &peaks);
        assert!(result.is_glycopeptide);
        // The true Man3 glycan (id 1) must be retrieved; Man5 (id 0) must be
        // filtered out (its Y1/Y2 complements are non-core).
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
        // A sialylated glycan (NeuAc present) with no NeuAc diagnostic ion fired
        // must be filtered out when diagnostic_filter is on.
        let glycans = vec![glycan(0, 3, 2, 0, 1, 0)]; // Man3 + NeuAc
        let idx = GlycanIonIndex::build(&glycans, 0.02, Some(20.0), GlycanCore::NGlycan);
        let cfg = GlycanConfig {
            diagnostic_filter: true,
            ..GlycanConfig::default()
        };
        // Peaks: core Y ions + HexNAc/HexHexNAc diagnostics, but NO sialic ion.
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
        assert!(result.candidates.is_empty(), "unsupported NeuAc must be filtered");
    }

    /// Phase 8 benchmark: index construction + query cost. Query time must scale
    /// with peaks × bins, not with the number of glycans (the whole point of the
    /// index). Prints timing for eyeballing; assertions keep it regression-safe.
    ///
    /// The glycan DB is read from the `.gdb` at `$ANDES_GLYCO_GDB` when set (so
    /// the benchmark measures the *input* database, e.g. pGlyco-N-Mouse.gdb's
    /// 1234 compositions); otherwise it falls back to the built-in common list.
    #[test]
    fn benchmark_index_construction_and_query() {
        let comps = match std::env::var("ANDES_GLYCO_GDB") {
            Ok(path) => {
                let content = std::fs::read_to_string(&path)
                    .unwrap_or_else(|e| panic!("read {path}: {e}"));
                crate::glycan_db::load_glycan_gdb(&content)
                    .unwrap_or_else(|e| panic!("parse {path}: {e}"))
            }
            // No built-in list anymore (the enumerator was removed); this benchmark
            // only runs against an explicit input .gdb.
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

        // 200-peak synthetic spectrum; query each peak's complementary mass.
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
