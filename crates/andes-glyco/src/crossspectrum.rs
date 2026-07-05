//! Cross-spectrum glycoform transfer (the andes-unique glyco lever).
//!
//! A single peptide is typically observed as MANY glycoforms (in serum N-glyco
//! data, ~6–7 per peptide): same peptide backbone, different glycan. Some
//! glycoforms fragment well (a strong trimannosyl-core Y-ladder → confidently
//! identified); their poorly-fragmenting siblings do not, and per-spectrum
//! candidate generation misses them because there is no core-Y ladder to anchor
//! the backbone.
//!
//! This module transfers a CONFIDENT backbone mass — learned from the well-
//! fragmented glycoforms in a first pass — to a second-pass spectrum whenever
//! `precursor − backbone` is a known glycan. The target spectrum needs NO core-Y
//! ladder of its own: the backbone is borrowed from its siblings. This is the
//! a cross-spectrum glyco engine cross-spectrum idea, and it is something per-spectrum engines
//! (the reference glyco engine) do not do.
//!
//! It is cheap: a small sorted whitelist + a binary-search glycan lookup per
//! precursor, NOT a brute force over the glycan list or a fragment-index scan.

use crate::glycan_db::GlycanComp;

/// A glyco-candidate spectrum as a graph node.
#[derive(Debug, Clone)]
pub struct GlycoNode {
    pub scan: u32,
    pub precursor_neutral: f64,
    pub rt_seconds: Option<f64>,
}

/// A confident Pass-1 seed backbone to propagate. `peptide_idx` indexes the
/// driver's candidate array so Pass-2 can re-score the exact peptide+mods.
#[derive(Debug, Clone)]
pub struct Seed {
    pub scan: u32,
    pub peptide_idx: u32,
    pub backbone_mass: f64,
    pub rt_seconds: Option<f64>,
    pub seed_score: f64,
    pub is_decoy: bool,
}

/// A backbone transferred onto an acceptor spectrum, plus its graph evidence.
#[derive(Debug, Clone)]
pub struct TransferredCandidate {
    pub acceptor_scan: u32,
    pub peptide_idx: u32,
    pub backbone_mass: f64,
    pub glycan: GlycanComp,
    pub graph_support: u32,
    pub seed_score: f64,
    pub rt_delta: f64,
    pub ungated: bool,
    pub is_decoy: bool,
}

/// One confident donor backbone plus the RETENTION-TIME window over which it (and
/// its glycoform siblings) were observed. RT gating is what makes cross-spectrum
/// transfer safe: a backbone may only be transferred to an acceptor spectrum that
/// CO-ELUTES with the donor, otherwise we would propose a backbone from an
/// unrelated peptide that merely shares a precursor-mass coincidence (this is why
/// the un-gated NULL-v1 scaffold was gated OFF — see
/// docs/plans/glyco/20-theory/why-andes-fails-and-succeed.md §4).
#[derive(Debug, Clone)]
struct BackboneEntry {
    backbone: f64,
    rt_min: f32,
    rt_max: f32,
}

/// A whitelist of confidently-identified backbone (peptide residue) masses from a
/// first pass, each carrying the donor RT window, sorted+deduplicated by backbone.
#[derive(Debug, Clone, Default)]
pub struct GlycoformWhitelist {
    entries: Vec<BackboneEntry>,
}

impl GlycoformWhitelist {
    /// Build from confident `(backbone_residue_mass, donor_rt_seconds)` donor
    /// observations, sorting by backbone and collapsing near-duplicates within
    /// `dedup_tol` Da (many glycoforms of one peptide contribute the same
    /// backbone). Collapsed donors extend the entry's `[rt_min, rt_max]` window,
    /// so the whitelist records the full elution span of each confident backbone.
    pub fn new(mut donors: Vec<(f64, f32)>, dedup_tol: f64) -> Self {
        donors.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        let mut out: Vec<BackboneEntry> = Vec::with_capacity(donors.len());
        for (bb, rt) in donors {
            match out.last_mut() {
                Some(e) if (bb - e.backbone).abs() <= dedup_tol => {
                    e.rt_min = e.rt_min.min(rt);
                    e.rt_max = e.rt_max.max(rt);
                }
                _ => out.push(BackboneEntry { backbone: bb, rt_min: rt, rt_max: rt }),
            }
        }
        GlycoformWhitelist { entries: out }
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Propose confident backbones consistent with `precursor_neutral` for an
    /// acceptor spectrum eluting at `acceptor_rt` (seconds). A backbone is
    /// proposed only if the acceptor co-elutes with the donor — `acceptor_rt`
    /// falls within `[rt_min − rt_window, rt_max + rt_window]` — AND
    /// `precursor_neutral − backbone` is within `tol` of a known glycan and
    /// ≥ `min_glycan`. Emits `(backbone, glycan)`.
    ///
    /// `rt_window` is the co-elution half-width in seconds (a cross-spectrum glyco engine uses
    /// ~±1800 s). `tol` is set from the PRECURSOR mass error, as in the
    /// peptide-first path.
    pub fn transfer(
        &self,
        precursor_neutral: f64,
        acceptor_rt: f32,
        rt_window: f32,
        glycan_sorted: &[(f64, usize)],
        glycans: &[GlycanComp],
        min_glycan: f64,
        tol: f64,
    ) -> Vec<(f64, GlycanComp)> {
        let mut out: Vec<(f64, GlycanComp)> = Vec::new();
        for e in &self.entries {
            let glycan_mass = precursor_neutral - e.backbone;
            if glycan_mass < min_glycan {
                // entries sorted by backbone ascending → larger backbone only
                // shrinks the implied glycan, so nothing further qualifies.
                break;
            }
            // RT gate: acceptor must co-elute with the donor window.
            if acceptor_rt < e.rt_min - rt_window || acceptor_rt > e.rt_max + rt_window {
                continue;
            }
            if let Some(g) = nearest_glycan(glycan_sorted, glycans, glycan_mass, tol) {
                out.push((e.backbone, g));
            }
        }
        out
    }
}

/// Nearest known glycan to `target` within `tol`, via a binary-search start on
/// the sorted `(mass, index)` view. (Local copy so the module is self-contained;
/// mirrors the driver's `nearest_glycan_mass`.)
fn nearest_glycan(
    sorted: &[(f64, usize)],
    glycans: &[GlycanComp],
    target: f64,
    tol: f64,
) -> Option<GlycanComp> {
    let lo = target - tol;
    let hi = target + tol;
    let start = sorted.partition_point(|&(m, _)| m < lo);
    let mut best: Option<(f64, usize)> = None;
    for &(m, gi) in &sorted[start..] {
        if m > hi {
            break;
        }
        let d = (m - target).abs();
        if best.map_or(true, |(bd, _)| d < bd) {
            best = Some((d, gi));
        }
    }
    best.map(|(_, gi)| glycans[gi].clone())
}

/// Propagate confident seed backbones to co-eluting, glycan-consistent acceptor
/// nodes. See module docs. Deterministic: output sorted by
/// (acceptor_scan, backbone_mass, glycan mass) with an index tiebreak.
pub fn propagate_transfers(
    seeds: &[Seed],
    nodes: &[GlycoNode],
    glycan_sorted: &[(f64, usize)],
    glycans: &[GlycanComp],
    rt_window: f32,
    min_glycan: f64,
    tol: f64,
) -> Vec<TransferredCandidate> {
    let co_elutes = |seed_rt: Option<f64>, node_rt: Option<f64>| -> (bool, bool) {
        // returns (passes_gate, ungated)
        match (seed_rt, node_rt) {
            (Some(s), Some(n)) => (((n - s).abs() as f32) <= rt_window, false),
            _ => (true, true), // missing RT on either side: accept, flag ungated
        }
    };
    let mut out: Vec<TransferredCandidate> = Vec::new();
    for seed in seeds {
        // First pass over nodes: collect accepting (non-self) acceptors so we
        // can attach the family-size support to each emitted candidate.
        let mut accepted: Vec<(u32, GlycanComp, f64, bool)> = Vec::new();
        for node in nodes {
            if node.scan == seed.scan {
                continue;
            }
            let (gate, ungated) = co_elutes(seed.rt_seconds, node.rt_seconds);
            if !gate {
                continue;
            }
            let glycan_mass = node.precursor_neutral - seed.backbone_mass;
            if glycan_mass < min_glycan {
                continue;
            }
            if let Some(g) = nearest_glycan(glycan_sorted, glycans, glycan_mass, tol) {
                let rt_delta = match (seed.rt_seconds, node.rt_seconds) {
                    (Some(s), Some(n)) => (n - s).abs(),
                    _ => 0.0,
                };
                accepted.push((node.scan, g, rt_delta, ungated));
            }
        }
        let support = accepted.len() as u32;
        for (acceptor_scan, glycan, rt_delta, ungated) in accepted {
            out.push(TransferredCandidate {
                acceptor_scan,
                peptide_idx: seed.peptide_idx,
                backbone_mass: seed.backbone_mass,
                glycan,
                graph_support: support,
                seed_score: seed.seed_score,
                rt_delta,
                ungated,
                is_decoy: seed.is_decoy,
            });
        }
    }
    // Deterministic total order (no HashMap anywhere above).
    out.sort_by(|a, b| {
        a.acceptor_scan
            .cmp(&b.acceptor_scan)
            .then(a.backbone_mass.total_cmp(&b.backbone_mass))
            .then(a.glycan.mass.total_cmp(&b.glycan.mass))
            .then(a.peptide_idx.cmp(&b.peptide_idx))
    });
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::glycan_db::n_glycan_list;
    use crate::glycan_mass::{HEX, HEXNAC};

    fn sorted_view(glycans: &[GlycanComp]) -> Vec<(f64, usize)> {
        let mut v: Vec<(f64, usize)> = glycans.iter().enumerate().map(|(i, g)| (g.mass, i)).collect();
        v.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        v
    }

    // Wide RT window so the non-RT tests below are unaffected by gating.
    const WIDE_RT: f32 = 1e9;

    #[test]
    fn whitelist_dedups_sibling_backbones() {
        let wl = GlycoformWhitelist::new(
            vec![(1500.0, 100.0), (1500.005, 110.0), (1500.01, 120.0), (2000.0, 100.0)],
            0.02,
        );
        // The three ~1500 backbones (sibling glycoforms of one peptide) collapse to one.
        assert_eq!(wl.len(), 2);
    }

    #[test]
    fn transfer_recovers_sibling_glycoform_without_a_ladder() {
        let glycans = n_glycan_list();
        let sorted = sorted_view(&glycans);
        // A peptide backbone confidently seen on a well-fragmented glycoform at RT 900 s.
        let backbone = 1500.0_f64;
        let wl = GlycoformWhitelist::new(vec![(backbone, 900.0)], 0.02);

        // A DIFFERENT glycoform of the same peptide: backbone + HexNAc2Hex3, co-eluting.
        let glycan_mass = 2.0 * HEXNAC + 3.0 * HEX;
        let precursor = backbone + glycan_mass;

        let hits = wl.transfer(precursor, 900.0, WIDE_RT, &sorted, &glycans, 406.0, 0.05);
        assert!(
            hits.iter().any(|(bb, _)| (bb - backbone).abs() < 0.02),
            "transfer must propose the confident sibling backbone, got {hits:?}"
        );
    }

    #[test]
    fn transfer_rejects_precursor_with_no_valid_glycan() {
        let glycans = n_glycan_list();
        let sorted = sorted_view(&glycans);
        let wl = GlycoformWhitelist::new(vec![(1500.0, 900.0)], 0.02);
        // Precursor implies a glycan of ~50 Da (below the N-glycan core) → nothing.
        let hits = wl.transfer(1550.0, 900.0, WIDE_RT, &sorted, &glycans, 406.0, 0.05);
        assert!(hits.is_empty(), "sub-core glycan must not transfer, got {hits:?}");
    }

    /// G4 RT gate: a backbone may transfer to a co-eluting acceptor but NOT to one
    /// outside the donor RT window — otherwise transfer propagates an unrelated
    /// peptide that only shares a precursor-mass coincidence.
    #[test]
    fn transfer_respects_rt_window() {
        let glycans = n_glycan_list();
        let sorted = sorted_view(&glycans);
        let backbone = 1500.0_f64;
        let wl = GlycoformWhitelist::new(vec![(backbone, 1000.0)], 0.02); // donor at RT 1000 s
        let precursor = backbone + 2.0 * HEXNAC + 3.0 * HEX;
        let window = 300.0_f32;

        // Acceptor co-eluting (RT 1100, |Δ| = 100 < 300) → transfer fires.
        let inside = wl.transfer(precursor, 1100.0, window, &sorted, &glycans, 406.0, 0.05);
        assert!(
            inside.iter().any(|(bb, _)| (bb - backbone).abs() < 0.02),
            "co-eluting acceptor must receive the transfer, got {inside:?}"
        );

        // Acceptor far away (RT 5000, |Δ| = 4000 > 300) → no transfer.
        let outside = wl.transfer(precursor, 5000.0, window, &sorted, &glycans, 406.0, 0.05);
        assert!(
            outside.is_empty(),
            "acceptor outside the RT window must NOT receive the transfer, got {outside:?}"
        );
    }

    #[test]
    fn propagate_transfers_recovers_sibling_and_counts_support() {
        let glycans = n_glycan_list();
        let sorted = sorted_view(&glycans);
        let backbone = 1500.0_f64;
        let g_a = 2.0 * HEXNAC + 3.0 * HEX;      // core
        let g_b = 2.0 * HEXNAC + 4.0 * HEX;      // +1 Hex sibling
        // seed on scan 1 (a well-fragmented glycoform), two co-eluting acceptors.
        let seeds = vec![Seed { scan: 1, peptide_idx: 9, backbone_mass: backbone,
            rt_seconds: Some(900.0), seed_score: 3.0, is_decoy: false }];
        let nodes = vec![
            GlycoNode { scan: 1, precursor_neutral: backbone + g_a, rt_seconds: Some(900.0) }, // self: skip
            GlycoNode { scan: 2, precursor_neutral: backbone + g_a, rt_seconds: Some(905.0) }, // sibling
            GlycoNode { scan: 3, precursor_neutral: backbone + g_b, rt_seconds: Some(910.0) }, // sibling
        ];
        let out = propagate_transfers(&seeds, &nodes, &sorted, &glycans, 300.0, 406.0, 0.05);
        // Two transfers (scan 2 and 3), NOT onto the seed's own scan 1.
        assert_eq!(out.len(), 2, "got {out:?}");
        assert!(out.iter().all(|t| t.acceptor_scan != 1));
        assert!(out.iter().all(|t| t.peptide_idx == 9 && (t.backbone_mass - backbone).abs() < 1e-6));
        // graph_support = family size (2 accepting non-self nodes).
        assert!(out.iter().all(|t| t.graph_support == 2), "support {out:?}");
    }

    #[test]
    fn propagate_transfers_respects_rt_window_and_marks_ungated() {
        let glycans = n_glycan_list();
        let sorted = sorted_view(&glycans);
        let backbone = 1500.0_f64;
        let g = 2.0 * HEXNAC + 3.0 * HEX;
        let seeds = vec![Seed { scan: 1, peptide_idx: 0, backbone_mass: backbone,
            rt_seconds: Some(1000.0), seed_score: 1.0, is_decoy: false }];
        // one co-eluting (RT 1100), one far (RT 5000), one with no RT (ungated).
        let nodes = vec![
            GlycoNode { scan: 2, precursor_neutral: backbone + g, rt_seconds: Some(1100.0) },
            GlycoNode { scan: 3, precursor_neutral: backbone + g, rt_seconds: Some(5000.0) },
            GlycoNode { scan: 4, precursor_neutral: backbone + g, rt_seconds: None },
        ];
        let out = propagate_transfers(&seeds, &nodes, &sorted, &glycans, 300.0, 406.0, 0.05);
        let scans: Vec<u32> = out.iter().map(|t| t.acceptor_scan).collect();
        assert!(scans.contains(&2), "co-eluting must transfer: {out:?}");
        assert!(!scans.contains(&3), "far RT must NOT transfer: {out:?}");
        // ungated node (no RT) still receives, flagged ungated.
        let u = out.iter().find(|t| t.acceptor_scan == 4).expect("ungated transfer");
        assert!(u.ungated, "no-RT acceptor must be marked ungated");
    }

    #[test]
    fn transfer_types_construct_and_clone() {
        let g = crate::glycan_db::n_glycan_list()[0].clone();
        let n = GlycoNode { scan: 5, precursor_neutral: 2000.0, rt_seconds: Some(900.0) };
        let s = Seed { scan: 5, peptide_idx: 3, backbone_mass: 1500.0, rt_seconds: Some(900.0), seed_score: 2.5, is_decoy: false };
        let t = TransferredCandidate { acceptor_scan: 7, peptide_idx: 3, backbone_mass: 1500.0,
            glycan: g, graph_support: 4, seed_score: 2.5, rt_delta: 12.0, ungated: false, is_decoy: false };
        assert_eq!(n.clone().scan, 5);
        assert_eq!(s.clone().peptide_idx, 3);
        assert_eq!(t.clone().graph_support, 4);
    }
}
