//! Composition-specific glycan Y-ion NODE SET and a log-likelihood-ratio score.
//!
//! WHY this exists rather than the linear ladder in `backbone::glycan_cumulative_adds`:
//! that ladder is a single CHAIN of cumulative additions which appends Fuc AFTER
//! every antenna monosaccharide. An N-glycan is not a chain — it is a branched
//! tree, and core fucose sits on the FIRST (reducing-end) HexNAc. The consequence
//! is structural, not cosmetic: the two most diagnostic core-fucose Y ions,
//! peptide+HexNAc+Fuc and peptide+2HexNAc+Fuc, are never predicted by the chain,
//! so a core-fucosylated glycopeptide can never be distinguished from its
//! non-fucosylated isobaric alternatives by the ladder term.
//!
//! The node set below is therefore built on the canonical N-glycan topology
//! (chitobiose Y1/Y2 → trimannosyl core → antennae), with an explicit
//! core-fucosylated twin of every core node, explicit sialic-loss shadow nodes
//! (these dominate HCD, where sialic acids are the first thing lost), and the
//! intact node.
//!
//! The score is in the form used by PTM-Shepherd's glycan assignment
//! (Polasky et al., Mol Cell Proteomics 2022, doi:10.1016/j.mcpro.2022.100205):
//! per-ion-class log-likelihood ratios of observing vs not observing an ion,
//! summed. pGlyco (Liu et al., Nat Commun 2017; Zeng et al., Nat Commun 2021)
//! uses the same "which Y ions does THIS composition predict" idea to make the
//! glycan axis discriminative. Counts enter as SQUARE ROOTS so that a larger
//! composition, which simply predicts more nodes, is not rewarded for size —
//! the same size bias that was measured on the linear ladder (K·ladder inversion,
//! 76.5 on wrong winners vs 53.1 on correct; see `backbone::ladder_norm_scale`).
//!
//! The decoy scorer is the mass-shifted-Y construction that five independent
//! engines converged on (pGlyco, GlycReSoft, Glyco-Decipher, GlycanFinder,
//! PTM-Shepherd): identical node set, interior nodes displaced by a random mass,
//! Y0/Y1/intact held fixed because those are fixed by the precursor and the
//! peptide and a decoy that moved them would be discriminable for the wrong reason.
//!
//! Nothing here is wired into the selector or the PIN; that is done elsewhere.

use crate::backbone::{best_frag_intensity, SpectrumStats};
use crate::glycan_db::GlycanComp;
use crate::glycan_mass::{FUC, HEX, HEXNAC, NEUAC, NEUGC};

/// Structural class of a Y node. The class, not the individual node, carries the
/// observation prior: per-node priors would need per-composition training data we
/// do not have, whereas class-level rates are estimable from any labelled corpus.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum YClass {
    /// Bare peptide backbone, no glycan retained.
    Y0,
    /// Chitobiose + trimannosyl core (Y1, Y2, N2H1..N2H3).
    Core,
    /// Core node carrying the core fucose — the class the linear ladder omits.
    CoreFuc,
    /// Partial antenna extension above the trimannosyl core.
    Antenna,
    /// Intact glycan minus one or two sialic acids (the dominant HCD losses).
    Sialylated,
    /// Intact glycopeptide (whole glycan retained).
    Intact,
}

/// Number of variants of [`YClass`]; used to size the fixed per-class accumulators
/// so the hot path needs no map and no allocation.
const N_CLASSES: usize = 6;

impl YClass {
    #[inline]
    fn idx(self) -> usize {
        match self {
            YClass::Y0 => 0,
            YClass::Core => 1,
            YClass::CoreFuc => 2,
            YClass::Antenna => 3,
            YClass::Sialylated => 4,
            YClass::Intact => 5,
        }
    }

    /// Probability that a node of this class is observed in a TRUE match.
    ///
    /// PLACEHOLDER CONSTANTS. They encode only the qualitative ordering that is
    /// uncontroversial in the literature (core > Y0 > core-fucose > intact >
    /// sialic-loss > antenna, because HCD strips the glycan from the outside in
    /// and the chitobiose bond is the last to break). They are NOT fitted: the
    /// intended replacement is a hit-rate table measured on our own labelled
    /// glyco corpus, exactly as PTM-Shepherd fits its Y-ion probabilities.
    #[inline]
    pub fn prior(self) -> f32 {
        match self {
            YClass::Y0 => 0.55,
            YClass::Core => 0.65,
            YClass::CoreFuc => 0.45,
            YClass::Antenna => 0.25,
            YClass::Sialylated => 0.30,
            YClass::Intact => 0.35,
        }
    }
}

/// One predicted Y ion: a glycan mass ADDED to the neutral peptide backbone.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct YNode {
    /// Glycan mass added to the NEUTRAL peptide backbone. 0.0 for Y0.
    pub neutral_add: f64,
    pub class: YClass,
    /// Copy of `class.prior()`, carried on the node so a caller that filters or
    /// re-weights nodes does not have to re-derive it.
    pub prior: f32,
}

/// Chance probability that an arbitrary m/z window contains a peak, i.e. the
/// probability of a "hit" under the null. 0.05 is the value PTM-Shepherd uses as
/// its default non-specific match rate; it makes `a_hit` positive for every class
/// above (all priors > 0.05) and `a_miss` negative for every class.
pub const DEFAULT_CHANCE: f32 = 0.05;

/// Minimum base-peak fraction for a matched peak to count as a HIT. Below this the
/// node is a MISS. Without a floor a noise-level peak inside the tolerance window
/// (a few 1e-4 of base) adds ~0 to the hit sum but removes a full unit from the
/// miss penalty, so on a dense HCD spectrum with a 1..5 charge sweep the null hit
/// rate exceeds `DEFAULT_CHANCE` and large compositions gain from noise.
pub const MIN_HIT_FRAC: f64 = 0.01;

/// A decoy node shifted to within this distance (Da) of any TARGET node mass is
/// re-drawn: the decoy would otherwise score the target's own peak.
const DECOY_COLLISION_DA: f64 = 0.5;

/// Hard cap on the returned node count.
///
/// WHY a cap at all: the antenna enumeration is a cross-product over the
/// HexNAc and Hex budgets, so a large composition (HexNAc6Hex10Fuc4NeuAc4)
/// generates ~39 antenna nodes on its own and the per-call cost is linear in the
/// node count (one `best_frag_intensity` binary search per node per charge).
/// WHY 40: the structural nodes (Y0 + 5 core + 5 core-fuc + up to 4 sialic-loss +
/// intact) can never exceed 16, so 40 leaves at least 24 antenna slots — enough
/// to cover every antenna node of a bi/tri-antennary glycan, which is the regime
/// essentially all real compositions live in. Only the antenna class is ever
/// trimmed, and it is trimmed from the HIGH-mass end: low-mass antenna nodes sit
/// closer to the core and are the ones actually observed.
pub const MAX_NODES: usize = 40;

/// Two node masses closer than this are the same ion for matching purposes
/// (0.1 mDa is far below any instrument tolerance we search at).
const MASS_DEDUP_TOL: f64 = 1e-4;

/// Composition-specific Y-ion node set on the canonical N-glycan topology.
///
/// Nodes are returned sorted by `neutral_add` ascending, deduplicated to
/// [`MASS_DEDUP_TOL`] keeping the highest-prior class at each mass, and capped at
/// [`MAX_NODES`]. No node ever exceeds the parent composition in any
/// monosaccharide, so the set is always a subgraph of the assigned glycan.
pub fn y_node_set(comp: &GlycanComp) -> Vec<YNode> {
    y_node_set_topology(comp, true)
}

/// The fucose placements a composition admits. `GlycanComp` carries COUNTS, not
/// topology: one fucose may sit on the core HexNAc (then it rides on every Y ion) or
/// on an antenna (then core and antenna Y ions are bare and only the intact node and
/// the sialic-loss shadows carry it). Scoring only the core placement penalised every
/// antenna-fucosylated glycopeptide by a full antenna's worth of misses, so both are
/// scored and the better one is taken (a max over latent trees). Afucosylated
/// compositions have one topology.
pub fn y_node_topologies(comp: &GlycanComp) -> Vec<Vec<YNode>> {
    if comp.fuc >= 1 {
        vec![y_node_set_topology(comp, true), y_node_set_topology(comp, false)]
    } else {
        vec![y_node_set_topology(comp, true)]
    }
}

/// [`y_node_set`] for one fucose placement: `core_fuc` puts the first fucose on the
/// reducing-end HexNAc (core-fucose twins of the core nodes, antennae in +Fuc form);
/// `false` places it on an antenna (core and antenna nodes bare).
pub fn y_node_set_topology(comp: &GlycanComp, core_fuc: bool) -> Vec<YNode> {
    // Structural (never-trimmed) nodes and antenna (trimmable) nodes are collected
    // separately so the cap can only ever bite the antennae.
    let mut fixed: Vec<YNode> = Vec::with_capacity(16);
    let mut antenna: Vec<YNode> = Vec::new();

    let push = |v: &mut Vec<YNode>, add: f64, class: YClass| {
        v.push(YNode {
            neutral_add: add,
            class,
            prior: class.prior(),
        });
    };

    // Y0: bare backbone.
    push(&mut fixed, 0.0, YClass::Y0);

    // Chitobiose + trimannosyl core, in biosynthetic order. Each rung is emitted
    // only if the composition actually contains the monosaccharides it needs, so
    // truncated (paucimannose / Man-less) compositions do not get phantom nodes.
    let hexnac = comp.hexnac as u32;
    let hex = comp.hex as u32;
    let mut core_adds: Vec<f64> = Vec::with_capacity(5);
    if hexnac >= 1 {
        core_adds.push(HEXNAC); // Y1
    }
    if hexnac >= 2 {
        core_adds.push(2.0 * HEXNAC); // Y2
        let core_hex = hex.min(3);
        for h in 1..=core_hex {
            core_adds.push(2.0 * HEXNAC + h as f64 * HEX); // N2H1..N2H3
        }
    }
    for &a in &core_adds {
        push(&mut fixed, a, YClass::Core);
    }

    // Core fucose: a twin of EVERY core node from Y1 onward. This is precisely the
    // set the linear cumulative-add ladder cannot produce, because it appends Fuc
    // last; peptide+HexNAc+Fuc and peptide+2HexNAc+Fuc are the two ions used in
    // practice to call core fucosylation.
    if core_fuc && comp.fuc >= 1 {
        for &a in &core_adds {
            push(&mut fixed, a + FUC, YClass::CoreFuc);
        }
    }

    // Antennae above the trimannosyl core: cumulative combinations of the HexNAc
    // budget left after the chitobiose (hexnac - 2) and the Hex budget left after
    // the trimannose (hex - 3). Enumerated as counts, not as an ordered walk,
    // because the branches are parallel and a real spectrum may lose them in any
    // order. Emitted only for compositions that HAVE a full trimannosyl core.
    //
    // Core fucose rides on the reducing-end HexNAc and so survives on essentially
    // every Y ion, so for a Fuc-containing composition every antenna node is
    // emitted in its +Fuc form (N3H3F, N3H4F, ...). Emitting them bare instead
    // predicted up to ~24 nodes a core-fucosylated glycopeptide can never show and
    // penalised fucosylation by exactly that miss count. A second (antenna) fucose
    // is not enumerated per node; it is carried only by the intact node and the
    // sialic-loss shadows, which use the full composition mass.
    if hexnac >= 2 && hex >= 3 {
        let base = 2.0 * HEXNAC + 3.0 * HEX + if core_fuc && comp.fuc >= 1 { FUC } else { 0.0 };
        let extra_hexnac = hexnac - 2;
        let extra_hex = hex - 3;
        for i in 0..=extra_hexnac {
            for j in 0..=extra_hex {
                if i == 0 && j == 0 {
                    continue; // that is N2H3, already a core node
                }
                let add = base + i as f64 * HEXNAC + j as f64 * HEX;
                if add >= comp.mass - MASS_DEDUP_TOL {
                    continue; // that is the intact node, emitted below
                }
                push(&mut antenna, add, YClass::Antenna);
            }
        }
    }

    // Sialic-loss shadows of the intact glycan. Under HCD the sialic acids are the
    // first residues lost, so [intact − NeuAc] and [intact − 2 NeuAc] are usually
    // more intense than the intact node itself on a sialylated glycopeptide.
    for k in 1..=2u32 {
        if comp.neuac as u32 >= k {
            push(&mut fixed, comp.mass - k as f64 * NEUAC, YClass::Sialylated);
        }
        if comp.neugc as u32 >= k {
            push(&mut fixed, comp.mass - k as f64 * NEUGC, YClass::Sialylated);
        }
    }

    // Intact glycopeptide.
    push(&mut fixed, comp.mass, YClass::Intact);

    // Trim the antennae from the high-mass end if the budget is exceeded. Sort is
    // by mass with a stable total order (no NaN can reach here: all masses are
    // finite sums of residue constants).
    antenna.sort_by(|a, b| a.neutral_add.total_cmp(&b.neutral_add));
    let budget = MAX_NODES.saturating_sub(fixed.len());
    antenna.truncate(budget);

    let mut all = fixed;
    all.append(&mut antenna);
    all.sort_by(|a, b| a.neutral_add.total_cmp(&b.neutral_add));

    // Deduplicate by mass, keeping the highest prior at each mass. A mass collision
    // is real (e.g. an antenna combination that lands on a core-fucose twin), and
    // scoring it twice would double-count one peak.
    let mut out: Vec<YNode> = Vec::with_capacity(all.len());
    for n in all {
        match out.last_mut() {
            Some(prev) if (n.neutral_add - prev.neutral_add).abs() <= MASS_DEDUP_TOL => {
                if n.prior > prev.prior {
                    prev.class = n.class;
                    prev.prior = n.prior;
                }
            }
            _ => out.push(n),
        }
    }
    out
}

/// Result of scoring a composition's Y node set against one spectrum.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct YTreeScore {
    /// Summed per-class log-likelihood ratio (see [`score_y_tree`]).
    pub llr: f32,
    /// `n_hits / n_nodes`.
    pub hit_frac: f32,
    pub n_hits: u16,
    pub n_nodes: u16,
    /// Missed nodes whose class prior is >= 0.5, i.e. ions a true match should
    /// almost certainly have shown. A large value with a positive `llr` is the
    /// signature of a composition that fits by antenna coincidence only.
    pub high_prior_missing: u16,
    /// Sum of matched base-peak-normalised intensities over all nodes.
    pub explained: f32,
}

/// Score a composition's Y node set against a spectrum.
///
/// `bb_neutral` is the NEUTRAL peptide backbone mass; node `i` is matched at
/// `bb_neutral + node.neutral_add` via [`best_frag_intensity`], which already
/// handles charge 1..=`max_charge` and isotope confirmation.
///
/// The LLR is
/// `sum over classes t of 2 * ( a_hit[t] * sqrt(U_t) + a_miss[t] * sqrt(V_t) )`
/// with `U_t` the intensity-weighted hit mass of class t (each hit contributes
/// `min(1, I/base)`, so one enormous peak cannot stand in for a whole class),
/// `V_t` the miss COUNT, `a_hit[t] = ln(prior_t / chance)` and
/// `a_miss[t] = ln((1 - prior_t) / (1 - chance))`. Both coefficients are derived
/// from the single `prior`/`chance` pair — there is deliberately no second table
/// to drift out of sync. The square roots are what keep a 40-node composition
/// from out-scoring a 10-node one on count alone.
///
/// Cost: one [`best_frag_intensity`] call per node (each a binary search over the
/// peak list per charge), i.e. at most [`MAX_NODES`] calls. One `Vec<YNode>`
/// allocation per call from [`y_node_set`]; the scoring loop itself allocates
/// nothing and uses fixed-size per-class accumulators.
pub fn score_y_tree(
    peaks: &[(f64, f32)],
    stats: &SpectrumStats,
    bb_neutral: f64,
    comp: &GlycanComp,
    tol_ppm: f64,
    max_charge: u8,
) -> YTreeScore {
    // Max over the fucose placements the composition admits (see
    // `y_node_topologies`); ties resolve to the core placement, listed first.
    y_node_topologies(comp)
        .iter()
        .map(|nodes| score_nodes(peaks, stats, bb_neutral, nodes, tol_ppm, max_charge))
        .fold(None, |best: Option<YTreeScore>, s| match best {
            Some(b) if b.llr >= s.llr => Some(b),
            _ => Some(s),
        })
        .unwrap_or_default()
}

/// Mass-shifted-Y decoy twin of [`score_y_tree`].
///
/// The node set is IDENTICAL, but every node except Y0, Y1 and the intact node is
/// displaced by an independent pseudo-random 1.0..30.0 Da offset derived from
/// `seed` and the node index. Y0, Y1 and the intact node are held fixed because
/// they are pinned by the peptide and the precursor: a decoy that moved them would
/// be rejected for the wrong reason and would inflate the target/decoy gap. This
/// is the construction pGlyco, GlycReSoft, Glyco-Decipher, GlycanFinder and
/// PTM-Shepherd all converged on (a reversal or shuffle of the composition was
/// tested by the field and rejected — a permuted composition is often another
/// REAL glycan).
pub fn score_y_tree_decoy(
    peaks: &[(f64, f32)],
    stats: &SpectrumStats,
    bb_neutral: f64,
    comp: &GlycanComp,
    tol_ppm: f64,
    max_charge: u8,
    seed: u64,
) -> YTreeScore {
    let mut topologies = y_node_topologies(comp);
    // Collision avoidance must cover EVERY placement's node masses: a shifted node
    // of the antenna-fucose tree may otherwise land on a core-fucose twin (e.g.
    // N2H1 − 16 Da = N2+Fuc) and score a real target ion.
    let avoid: Vec<f64> = topologies.iter().flatten().map(|n| n.neutral_add).collect();
    // The twin mirrors the topology the TARGET chose (same argmax, same tie rule as
    // `score_y_tree`). Letting the decoy take its own max would let it pick the
    // smaller tree purely for its fewer miss penalties, which is not a like-for-like
    // control of the target's node set.
    let mut chosen = 0usize;
    let mut best = f32::NEG_INFINITY;
    for (k, nodes) in topologies.iter().enumerate() {
        let llr = score_nodes(peaks, stats, bb_neutral, nodes, tol_ppm, max_charge).llr;
        if llr > best {
            best = llr;
            chosen = k;
        }
    }
    let shifted = shift_decoy_nodes(topologies.swap_remove(chosen), seed, &avoid);
    score_nodes(peaks, stats, bb_neutral, &shifted, tol_ppm, max_charge)
}

/// Mass-shift every non-anchor node of one topology's node set (the decoy twin of
/// exactly that node set; the target and the twin therefore differ ONLY in the
/// shifts).
fn shift_decoy_nodes(mut nodes: Vec<YNode>, seed: u64, target_masses: &[f64]) -> Vec<YNode> {
    for (i, n) in nodes.iter_mut().enumerate() {
        if is_anchor(n) {
            continue;
        }
        // Mix the node INDEX with the seed so two nodes of the same class in the
        // same spectrum get independent shifts, and so the whole decoy is
        // reproducible from (seed, node order) alone. Up to 8 draws: a shift that
        // lands on another TARGET node (adjacent nodes can sit 16 Da apart, e.g.
        // N2+Fuc vs N2H1) is re-drawn so the decoy never scores a real target ion.
        let original = n.neutral_add;
        // Best-of-N by collision MARGIN. The previous form assigned `shifted` before
        // testing `collides` and never restored it, so when all 8 draws collided the
        // node was left at the last colliding mass -- inside DECOY_COLLISION_DA of a
        // real target node, or below HEXNAC. Such a decoy node scores a genuine target
        // Y ion, inflating the decoy LLR and shrinking YTreeDecoyGap exactly for the
        // dense, large compositions (the plasma regime) where re-draws exhaust most
        // often. Restoring `original` on exhaustion would be WORSE, not better: the
        // original IS the target mass, a perfect collision rather than a near one. So
        // keep the draw that lands FURTHEST from any target node and commit that.
        let mut committed = false;
        let mut best: Option<(f64, f64)> = None; // (margin, shifted)
        for attempt in 0..8u64 {
            let r = splitmix64(seed ^ splitmix64(i as u64 + 1 + attempt * 0x9E37));
            // 53-bit mantissa fraction -> uniform in [0,1), then 1.0..30.0 Da, sign
            // from the low bit. The 1 Da floor keeps the shifted node clear of the
            // target's isotope envelope; 30 Da keeps it inside the same region of the
            // spectrum, so the decoy samples the same peak density as the target.
            let frac = (r >> 11) as f64 / (1u64 << 53) as f64;
            let mag = 1.0 + frac * 29.0;
            let shifted = if r & 1 == 0 { original + mag } else { original - mag };
            // Distance to the nearest target node; below HEXNAC is disqualifying, so
            // score it as the worst possible margin rather than a real distance.
            let margin = if shifted < HEXNAC {
                f64::NEG_INFINITY
            } else {
                target_masses
                    .iter()
                    .map(|&m| (m - shifted).abs())
                    .fold(f64::INFINITY, f64::min)
            };
            if margin >= DECOY_COLLISION_DA {
                n.neutral_add = shifted;
                committed = true;
                break;
            }
            if best.is_none_or(|(bm, _)| margin > bm) {
                best = Some((margin, shifted));
            }
        }
        if !committed {
            if let Some((_, shifted)) = best {
                n.neutral_add = shifted;
            }
        }
    }
    nodes
}

/// Nodes the decoy must not move: Y0, Y1 (+HexNAc) and the intact node.
#[inline]
fn is_anchor(n: &YNode) -> bool {
    match n.class {
        YClass::Y0 | YClass::Intact => true,
        YClass::Core => (n.neutral_add - HEXNAC).abs() <= MASS_DEDUP_TOL,
        _ => false,
    }
}

/// Shared scoring core for the target and decoy node sets, so the only permitted
/// difference between them is the node masses.
fn score_nodes(
    peaks: &[(f64, f32)],
    stats: &SpectrumStats,
    bb_neutral: f64,
    nodes: &[YNode],
    tol_ppm: f64,
    max_charge: u8,
) -> YTreeScore {
    // Fixed-size per-class accumulators: no map, no allocation in the hot path.
    let mut u = [0.0f64; N_CLASSES]; // intensity-weighted hit mass
    let mut v = [0.0f64; N_CLASSES]; // miss count
    let mut n_hits = 0u16;
    let mut high_prior_missing = 0u16;
    let mut explained = 0.0f64;

    for n in nodes {
        let raw = best_frag_intensity(
            peaks,
            stats.sorted,
            bb_neutral + n.neutral_add,
            tol_ppm,
            max_charge,
        );
        let k = n.class.idx();
        let frac = if raw > 0.0 { raw as f64 / stats.base } else { 0.0 };
        if frac >= MIN_HIT_FRAC {
            explained += frac;
            u[k] += frac.min(1.0);
            n_hits += 1;
        } else {
            v[k] += 1.0;
            if n.prior >= 0.5 {
                high_prior_missing += 1;
            }
        }
    }

    let chance = DEFAULT_CHANCE as f64;
    let mut llr = 0.0f64;
    for (k, class) in [
        YClass::Y0,
        YClass::Core,
        YClass::CoreFuc,
        YClass::Antenna,
        YClass::Sialylated,
        YClass::Intact,
    ]
    .iter()
    .enumerate()
    {
        if u[k] == 0.0 && v[k] == 0.0 {
            continue;
        }
        let p = class.prior() as f64;
        let a_hit = (p / chance).ln();
        let a_miss = ((1.0 - p) / (1.0 - chance)).ln();
        llr += 2.0 * (a_hit * u[k].sqrt() + a_miss * v[k].sqrt());
    }

    let n_nodes = nodes.len() as u16;
    YTreeScore {
        llr: llr as f32,
        hit_frac: if n_nodes > 0 {
            n_hits as f32 / n_nodes as f32
        } else {
            0.0
        },
        n_hits,
        n_nodes,
        high_prior_missing,
        explained: explained as f32,
    }
}

/// Deterministic mixing (splitmix64). Written locally rather than shared with the
/// private copy in `backbone.rs` so this module can be reviewed and moved on its own.
#[inline]
fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = x;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

#[cfg(test)]
mod shift_decoy_tests {
    use super::*;

    fn node(mass: f64) -> YNode {
        YNode { neutral_add: mass, class: YClass::Antenna, prior: YClass::Antenna.prior() }
    }

    fn nearest(masses: &[f64], v: f64) -> f64 {
        masses.iter().map(|&m| (m - v).abs()).fold(f64::INFINITY, f64::min)
    }

    /// Happy path: with sparse targets a collision-free shift is found and committed.
    #[test]
    fn shifted_decoy_node_clears_every_target_when_a_free_slot_exists() {
        let original = 1000.0;
        let targets = vec![original];
        let out = shift_decoy_nodes(vec![node(original)], 0xDEAD_BEEF, &targets);
        let got = out[0].neutral_add;
        assert!(got >= HEXNAC, "decoy node fell below HEXNAC: {got}");
        assert!(
            nearest(&targets, got) >= DECOY_COLLISION_DA,
            "decoy node {got} sits within {DECOY_COLLISION_DA} Da of a target",
        );
    }

    /// Replicate the 8 deterministic draws the shifter makes for node index 0, so the
    /// test can name the candidate it must choose. Mirrors the draw arithmetic exactly.
    fn draws(seed: u64, original: f64) -> Vec<f64> {
        (0..8u64)
            .map(|attempt| {
                let r = splitmix64(seed ^ splitmix64(1 + attempt * 0x9E37));
                let frac = (r >> 11) as f64 / (1u64 << 53) as f64;
                let mag = 1.0 + frac * 29.0;
                if r & 1 == 0 { original + mag } else { original - mag }
            })
            .collect()
    }

    /// Regression: when EVERY draw collides, the committed node must be the draw that
    /// lands FURTHEST from any target -- not the last draw, which is what the original
    /// code left behind because it assigned before testing and never restored.
    ///
    /// This asserts the exact chosen value, and first asserts that best != last for
    /// this seed, so the test genuinely discriminates the two behaviours rather than
    /// passing under both.
    #[test]
    fn exhausted_draws_commit_the_furthest_candidate_not_the_last() {
        let original = 5000.0;
        // Dense grid across the whole +/-30 Da window: every draw collides.
        let targets: Vec<f64> = (-700..=700).map(|k| original + k as f64 * 0.05).collect();

        let cand = draws(0x1234_5678, original);
        let margins: Vec<f64> = cand.iter().map(|&v| nearest(&targets, v)).collect();
        assert!(
            margins.iter().all(|&m| m < DECOY_COLLISION_DA),
            "test setup is wrong: some draw was collision-free, so exhaustion is untested",
        );
        let (best_i, _) = margins
            .iter()
            .enumerate()
            .fold((0usize, f64::NEG_INFINITY), |(bi, bm), (i, &m)| {
                if m > bm { (i, m) } else { (bi, bm) }
            });
        assert_ne!(
            best_i, 7,
            "seed chosen badly: the best draw IS the last one, so this test cannot \
             distinguish the fixed behaviour from the bug it guards",
        );

        let out = shift_decoy_nodes(vec![node(original)], 0x1234_5678, &targets);
        let got = out[0].neutral_add;
        assert_eq!(
            got, cand[best_i],
            "expected the furthest candidate {} (draw {best_i}), got {got}; the last \
             draw was {} -- committing that is the bug this test guards",
            cand[best_i], cand[7],
        );
    }

    /// Anchors (Y0, Y1, intact) are never moved, collision or not.
    #[test]
    fn anchor_nodes_are_never_shifted() {
        let anchors = vec![
            YNode { neutral_add: 0.0, class: YClass::Y0, prior: YClass::Y0.prior() },
            YNode { neutral_add: HEXNAC, class: YClass::Core, prior: YClass::Core.prior() },
        ];
        let before: Vec<f64> = anchors.iter().map(|n| n.neutral_add).collect();
        let out = shift_decoy_nodes(anchors, 42, &[0.0, HEXNAC]);
        let after: Vec<f64> = out.iter().map(|n| n.neutral_add).collect();
        assert_eq!(before, after, "an anchor node was shifted");
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::glycan_mass::PROTON;

    fn comp(hexnac: u8, hex: u8, fuc: u8, neuac: u8, neugc: u8) -> GlycanComp {
        let mass = hexnac as f64 * HEXNAC
            + hex as f64 * HEX
            + fuc as f64 * FUC
            + neuac as f64 * NEUAC
            + neugc as f64 * NEUGC;
        GlycanComp {
            hexnac,
            hex,
            fuc,
            neuac,
            neugc,
            mass,
        }
    }

    /// Singly-protonated peak list for a set of neutral masses, m/z-ascending.
    fn peaks_for(neutrals: &[f64], intensity: f32) -> Vec<(f64, f32)> {
        let mut p: Vec<(f64, f32)> = neutrals.iter().map(|&m| (m + PROTON, intensity)).collect();
        p.sort_by(|a, b| a.0.total_cmp(&b.0));
        p
    }

    fn has_mass(nodes: &[YNode], m: f64) -> bool {
        nodes.iter().any(|n| (n.neutral_add - m).abs() < 1e-6)
    }

    #[test]
    fn core_fucose_nodes_present_and_absent_from_linear_order() {
        // HexNAc4Hex5Fuc1: a plain core-fucosylated biantennary complex glycan.
        let c = comp(4, 5, 1, 0, 0);
        let nodes = y_node_set(&c);
        // The two diagnostic core-fucose ions.
        assert!(
            has_mass(&nodes, HEXNAC + FUC),
            "peptide+HexNAc+Fuc missing from node set"
        );
        assert!(
            has_mass(&nodes, 2.0 * HEXNAC + FUC),
            "peptide+2HexNAc+Fuc missing from node set"
        );

        // The linear cumulative-add ladder (chitobiose -> trimannose -> antennae ->
        // fucose -> sialics) reproduced here: Fuc enters only after everything else,
        // so neither diagnostic mass can appear.
        let mut linear: Vec<f64> = Vec::new();
        let mut cum = 0.0;
        let order = [HEXNAC, HEXNAC, HEX, HEX, HEX, HEXNAC, HEXNAC, HEX, HEX, FUC];
        for m in order {
            cum += m;
            linear.push(cum);
        }
        assert!(
            !linear.iter().any(|&m| (m - (HEXNAC + FUC)).abs() < 1e-6),
            "linear ladder unexpectedly contains peptide+HexNAc+Fuc"
        );
        assert!(
            !linear
                .iter()
                .any(|&m| (m - (2.0 * HEXNAC + FUC)).abs() < 1e-6),
            "linear ladder unexpectedly contains peptide+2HexNAc+Fuc"
        );
    }

    #[test]
    fn antenna_nodes_carry_the_core_fucose() {
        // HexNAc4Hex5Fuc1: every antenna node must be the +Fuc form; the bare
        // antenna masses (N3H3, N3H4, ...) are ions this glycopeptide cannot yield.
        let c = comp(4, 5, 1, 0, 0);
        let nodes = y_node_set(&c);
        let antenna: Vec<&YNode> = nodes.iter().filter(|n| n.class == YClass::Antenna).collect();
        assert!(!antenna.is_empty());
        let n3h3 = 3.0 * HEXNAC + 3.0 * HEX;
        assert!(has_mass(&nodes, n3h3 + FUC), "N3H3F missing");
        assert!(!has_mass(&nodes, n3h3), "bare N3H3 must not be predicted for a Fuc glycan");
        // And the afucosylated composition predicts the bare form.
        let bare = y_node_set(&comp(4, 5, 0, 0, 0));
        assert!(has_mass(&bare, n3h3));
        assert!(!has_mass(&bare, n3h3 + FUC));
    }

    #[test]
    fn antenna_fucose_placement_is_not_penalised() {
        // Spectrum of an ANTENNA-fucosylated HexNAc4Hex5Fuc1: bare core, bare antennae,
        // intact node. Scoring only the core-fucose tree would miss every +Fuc antenna
        // node; the topology max must score it as well as its afucosylated twin
        // composition scores its own (identical) core + antenna ions.
        let c = comp(4, 5, 1, 0, 0);
        let bb = 1500.0;
        let bare = y_node_set_topology(&c, false);
        let masses: Vec<f64> = bare.iter().map(|n| bb + n.neutral_add).collect();
        let peaks = peaks_for(&masses, 1000.0);
        let stats = SpectrumStats::new(&peaks);
        let with_max = score_y_tree(&peaks, &stats, bb, &c, 20.0, 2);
        let core_only = score_nodes(&peaks, &stats, bb, &y_node_set_topology(&c, true), 20.0, 2);
        assert!(with_max.llr > core_only.llr, "{} vs {}", with_max.llr, core_only.llr);
        assert_eq!(with_max.n_hits as usize, bare.len(), "{with_max:?}");
        // And a core-fucosylated spectrum still prefers the core tree.
        let core = y_node_set_topology(&c, true);
        let peaks2 = peaks_for(&core.iter().map(|n| bb + n.neutral_add).collect::<Vec<_>>(), 1000.0);
        let stats2 = SpectrumStats::new(&peaks2);
        let s2 = score_y_tree(&peaks2, &stats2, bb, &c, 20.0, 2);
        assert_eq!(s2.n_hits as usize, core.len(), "{s2:?}");
        // Decoy twin follows the same max and stays below the target on both.
        for (pk, st, tgt) in [(&peaks, &stats, with_max), (&peaks2, &stats2, s2)] {
            let d = score_y_tree_decoy(pk, st, bb, &c, 20.0, 2, 7);
            assert!(d.llr < tgt.llr, "decoy {} >= target {}", d.llr, tgt.llr);
        }
    }

    #[test]
    fn noise_level_peak_is_a_miss_not_a_hit() {
        let c = comp(4, 5, 1, 0, 0);
        let bb = 1500.0;
        let nodes = y_node_set(&c);
        let core: Vec<f64> = nodes
            .iter()
            .filter(|n| matches!(n.class, YClass::Y0 | YClass::Core | YClass::CoreFuc))
            .map(|n| bb + n.neutral_add)
            .collect();
        // Core nodes at full intensity, plus EVERY antenna node at 1e-4 of base.
        let mut peaks = peaks_for(&core, 1000.0);
        let faint: Vec<f64> = nodes
            .iter()
            .filter(|n| n.class == YClass::Antenna)
            .map(|n| bb + n.neutral_add)
            .collect();
        peaks.extend(peaks_for(&faint, 0.1));
        peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
        let stats = SpectrumStats::new(&peaks);
        let s = score_y_tree(&peaks, &stats, bb, &c, 20.0, 2);
        assert_eq!(
            s.n_hits as usize,
            core.len(),
            "faint antenna peaks must not count as hits: {s:?}"
        );
        // Removing the faint peaks entirely must give the SAME score.
        let clean = peaks_for(&core, 1000.0);
        let cstats = SpectrumStats::new(&clean);
        let s2 = score_y_tree(&clean, &cstats, bb, &c, 20.0, 2);
        assert!((s.llr - s2.llr).abs() < 1e-6, "{} vs {}", s.llr, s2.llr);
    }

    #[test]
    fn decoy_nodes_never_land_on_a_target_node() {
        for (hn, hx, f, na) in [(4u8, 5u8, 1u8, 0u8), (5, 6, 2, 2), (2, 5, 0, 0), (6, 7, 1, 3)] {
            let c = comp(hn, hx, f, na, 0);
            let target = y_node_set(&c);
            let bb = 1200.0;
            let peaks = peaks_for(
                &target.iter().map(|n| bb + n.neutral_add).collect::<Vec<_>>(),
                1000.0,
            );
            let stats = SpectrumStats::new(&peaks);
            for seed in 1..=20u64 {
                let d = score_y_tree_decoy(&peaks, &stats, bb, &c, 20.0, 2, seed);
                // Only the anchors (Y0, Y1, intact) may hit on a target-only spectrum.
                assert!(
                    d.n_hits <= 3,
                    "decoy hit {} target nodes for {:?} seed {seed}",
                    d.n_hits,
                    (hn, hx, f, na)
                );
            }
        }
    }

    #[test]
    fn core_hits_score_positive_and_empty_scores_negative() {
        let c = comp(4, 5, 1, 0, 0);
        let bb = 1500.0;
        let nodes = y_node_set(&c);

        // Spectrum holding exactly the Y0 + core + core-fucose nodes.
        let core_masses: Vec<f64> = nodes
            .iter()
            .filter(|n| matches!(n.class, YClass::Y0 | YClass::Core | YClass::CoreFuc))
            .map(|n| bb + n.neutral_add)
            .collect();
        let peaks = peaks_for(&core_masses, 1000.0);
        let stats = SpectrumStats::new(&peaks);
        let hit = score_y_tree(&peaks, &stats, bb, &c, 20.0, 2);
        assert!(hit.llr > 0.0, "core-only spectrum scored {}", hit.llr);
        assert_eq!(hit.n_hits as usize, core_masses.len());
        assert_eq!(hit.high_prior_missing, 0);

        // A spectrum with nothing at any node mass.
        let junk = peaks_for(&[bb + 7777.7, bb + 8888.8], 1000.0);
        let jstats = SpectrumStats::new(&junk);
        let miss = score_y_tree(&junk, &jstats, bb, &c, 20.0, 2);
        assert!(miss.llr < 0.0, "empty spectrum scored {}", miss.llr);
        assert_eq!(miss.n_hits, 0);
        assert!(miss.high_prior_missing > 0);
        assert!(hit.llr > miss.llr);
    }

    #[test]
    fn near_isobaric_compositions_get_different_node_sets_and_scores() {
        // 2 Fuc (292.116) vs 1 NeuAc (291.095): ~1.02 Da apart in total mass but
        // structurally unrelated node sets.
        let a = comp(4, 5, 2, 0, 0);
        let b = comp(4, 5, 0, 1, 0);
        assert!(
            (a.mass - b.mass).abs() < 1.1,
            "compositions not near-isobaric"
        );

        let na = y_node_set(&a);
        let nb = y_node_set(&b);
        let masses_a: Vec<f64> = na.iter().map(|n| n.neutral_add).collect();
        let masses_b: Vec<f64> = nb.iter().map(|n| n.neutral_add).collect();
        assert_ne!(masses_a, masses_b, "node sets identical for A and B");
        // A has core-fucose nodes; B has a sialic-loss node.
        assert!(na.iter().any(|n| n.class == YClass::CoreFuc));
        assert!(!nb.iter().any(|n| n.class == YClass::CoreFuc));
        assert!(nb.iter().any(|n| n.class == YClass::Sialylated));

        // Build the spectrum from A's nodes: A must win.
        let bb = 1500.0;
        let peaks = peaks_for(
            &na.iter().map(|n| bb + n.neutral_add).collect::<Vec<_>>(),
            1000.0,
        );
        let stats = SpectrumStats::new(&peaks);
        let sa = score_y_tree(&peaks, &stats, bb, &a, 20.0, 2);
        let sb = score_y_tree(&peaks, &stats, bb, &b, 20.0, 2);
        assert!(
            sa.llr > sb.llr,
            "spectrum built for A scored A={} B={}",
            sa.llr,
            sb.llr
        );
    }

    #[test]
    fn node_count_capped_for_large_composition() {
        let big = comp(6, 10, 4, 4, 0);
        let nodes = y_node_set(&big);
        assert!(
            nodes.len() <= MAX_NODES,
            "node set of {} exceeds cap",
            nodes.len()
        );
        // The cap must actually bind for this composition, otherwise the test is vacuous.
        assert_eq!(nodes.len(), MAX_NODES);
        // Structural nodes survive the trim.
        assert!(has_mass(&nodes, 0.0));
        assert!(has_mass(&nodes, HEXNAC));
        assert!(has_mass(&nodes, big.mass));
        assert!(nodes.iter().any(|n| n.class == YClass::CoreFuc));
        assert!(nodes.iter().any(|n| n.class == YClass::Sialylated));
    }

    #[test]
    fn decoy_scores_below_target_on_target_spectrum() {
        let c = comp(4, 5, 1, 0, 0);
        let bb = 1500.0;
        let nodes = y_node_set(&c);
        let peaks = peaks_for(
            &nodes.iter().map(|n| bb + n.neutral_add).collect::<Vec<_>>(),
            1000.0,
        );
        let stats = SpectrumStats::new(&peaks);
        let target = score_y_tree(&peaks, &stats, bb, &c, 20.0, 2);
        for seed in 0..16u64 {
            let decoy = score_y_tree_decoy(&peaks, &stats, bb, &c, 20.0, 2, seed);
            assert!(
                decoy.llr < target.llr,
                "seed {seed}: decoy {} >= target {}",
                decoy.llr,
                target.llr
            );
            assert_eq!(decoy.n_nodes, target.n_nodes);
        }
    }

    #[test]
    fn deterministic() {
        let c = comp(5, 6, 1, 2, 0);
        let bb = 1234.5678;
        let nodes = y_node_set(&c);
        assert_eq!(nodes, y_node_set(&c));
        let peaks = peaks_for(
            &nodes.iter().map(|n| bb + n.neutral_add).collect::<Vec<_>>(),
            777.0,
        );
        let stats = SpectrumStats::new(&peaks);
        assert_eq!(
            score_y_tree(&peaks, &stats, bb, &c, 20.0, 3),
            score_y_tree(&peaks, &stats, bb, &c, 20.0, 3)
        );
        assert_eq!(
            score_y_tree_decoy(&peaks, &stats, bb, &c, 20.0, 3, 42),
            score_y_tree_decoy(&peaks, &stats, bb, &c, 20.0, 3, 42)
        );
    }

    #[test]
    fn no_node_exceeds_parent_composition() {
        for &(hn, hx, fc, na) in &[
            (2u8, 2u8, 0u8, 0u8),
            (3, 3, 1, 1),
            (2, 9, 0, 0),
            (4, 5, 3, 2),
        ] {
            let c = comp(hn, hx, fc, na, 0);
            for n in y_node_set(&c) {
                assert!(
                    n.neutral_add <= c.mass + MASS_DEDUP_TOL,
                    "node {} exceeds glycan mass {}",
                    n.neutral_add,
                    c.mass
                );
                assert!(n.neutral_add >= 0.0);
            }
        }
    }
}
