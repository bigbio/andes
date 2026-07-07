// Carrier struct for glyco-aware PSM scoring features.
//
// `GlycoPsmKey` bundles all glycan-level evidence gathered for a single PSM
// into one value that can be stored, cloned, and passed to downstream
// re-scorers or PIN writers without re-computing spectra.

use crate::glycan_db::GlycanComp;
use crate::hybrid::Source;
use std::cmp::Ordering;

/// Collapse-selection mode, read ONCE from the environment (process-constant).
///
/// The top-1-per-scan collapse (required for honest per-scan TDC FDR) keeps a
/// single PSM per spectrum, so *which* backbone/glycan it keeps determines the
/// recovered ID. The DEFAULT selects by the core-Y ladder first (the reliable
/// glyco evidence in HCD); `ANDES_GLYCO_SELECT=rank` forces the old peptide b/y
/// `rank_score`-first rule (kept as an escape hatch for A/B and regression).
///
/// DETERMINISTIC RESULT (2026-07-05, PXD025455 Fc3_r1; clean1==clean2
/// byte-identical): core-Y-ladder primary = 218 @1% FDR / 90 backbone-correct;
/// b/y-rank primary = 128 / 53 — a +70% win. In stepped-HCD the peptide b/y
/// ions are sparse and rank_score is noisy, so it lifts wrong backbones into the
/// per-scan winner slot; the core-Y ladder is the stronger, more reliable
/// evidence. The earlier "ladder is WORSE (260→197)" verdict was a
/// NON-DETERMINISM artifact — it was measured on the pre-fix pipeline whose ±90
/// run-to-run jitter swamped the real +90 signal (see order_peptide_first fix).
pub fn y_primary_selection() -> bool {
    y_primary_from_env(std::env::var("ANDES_GLYCO_SELECT").ok().as_deref())
}

/// Pure decision for [`y_primary_selection`], split out so the default/override
/// logic is testable without racing on a process-global env var.
fn y_primary_from_env(v: Option<&str>) -> bool {
    match v {
        // Explicit escape hatch to the OLD b/y-rank-primary collapse.
        Some(s) if s.eq_ignore_ascii_case("rank") || s.eq_ignore_ascii_case("byrank") => false,
        // Unset, "yladder", or any other value → core-Y ladder primary (default).
        _ => true,
    }
}

/// Selector mode, read ONCE from the environment (process-constant).
///
/// DEFAULT (`legacy`): the per-scan winner is chosen by the bare b/y rank-LLR
/// (`collapse_cmp`), byte-identical to the shipping 253/97 baseline. When set to
/// `combined`, a PEPTIDE-CONDITIONED re-rank of the accepted-backbone shortlist
/// is applied (P1b + L4): the winner is chosen by [`glyco_combined_selector_score`],
/// which folds the peptide-mass-conditioned Y0/Y1 anchor and (optionally) the
/// frag/rich-ion intensity-model LLR into the b/y rank so a true backbone with a
/// strong Y0/Y1 anchor but a weaker b/y rank than a competitor is selected #1.
///
/// P0b proved retention (YINDEX) surfaces the true backbone into the pool but the
/// bare-rank selector still out-ranks it; the combined selector is the lever that
/// converts a recovered candidate into a top-1 correct id.
pub fn combined_selector_on() -> bool {
    combined_selector_from_env(std::env::var("ANDES_GLYCO_SELECTOR").ok().as_deref())
}

/// Pure decision for [`combined_selector_on`], split out so the default/override
/// logic is testable without racing on a process-global env var.
fn combined_selector_from_env(v: Option<&str>) -> bool {
    matches!(v, Some(s) if s.eq_ignore_ascii_case("combined"))
}

/// Seed weights for the combined peptide-conditioned selector (P1b/L4). These
/// are SEED values only — the roadmap calls for learning them (a glyco search engine uses a
/// RankSVM w=0.35 glycan / 0.65 peptide fusion). They are exposed as constants
/// so the b/y-vs-anchor-vs-intensity balance is a single, tunable knob group.
///
/// - `W_BY` scales the bare b/y rank-LLR (the backbone-SEQUENCE axis; already the
///   sole legacy selector). Kept at 1.0 so the combined score is an ADDITIVE
///   composition on top of the existing rank, never a rewrite of `score_psm`.
/// - `W_ANCHOR` scales the peptide-mass-conditioned Y0/Y1 anchor intensity
///   (range ~0..4: 2 ions × ≤2 charges, each a base-peak fraction). Seeded high
///   enough that a strong anchor (~2-4) overcomes the R7 median rank gap (~7) in
///   favour of the true backbone — the exact failure P0b diagnosed.
/// - `W_INTENSITY` scales the frag/rich-ion intensity-model LLR (0.0 when no
///   frag/rich-ion model is loaded, so this term is inert on models without it).
pub const GLYCO_SEL_W_BY: f32 = 1.0;
pub const GLYCO_SEL_W_ANCHOR: f32 = 8.0;
pub const GLYCO_SEL_W_INTENSITY: f32 = 4.0;

/// Combined peptide-conditioned selector score (P1b + L4).
///
/// Additive composition of the three peptide-conditioned axes andes already
/// computes: the bare b/y rank-LLR (`rank`), the peptide-mass-conditioned Y0/Y1
/// anchor intensity (`y0y1_anchor`), and the frag/rich-ion intensity-model LLR
/// (`intensity_llr`; pass 0.0 to disable that axis, e.g. when no model is loaded
/// or to ship the anchor-only variant). Higher is better. This does NOT modify
/// `score_psm`; it composes on top of its output.
pub fn glyco_combined_selector_score(rank: f32, y0y1_anchor: f32, intensity_llr: f32) -> f32 {
    GLYCO_SEL_W_BY * rank
        + GLYCO_SEL_W_ANCHOR * y0y1_anchor
        + GLYCO_SEL_W_INTENSITY * intensity_llr
}

/// Total order for the top-1-per-scan collapse: `max_by(collapse_cmp(...))`
/// yields the emitted winner. This ordering is the SINGLE SOURCE OF TRUTH shared
/// by the driver's pre-feature reduction (glyco_search) and the PIN writer's
/// `select_emitted_hits` — they MUST agree, or a scan's driver-emitted winner
/// and PIN-written winner diverge (a real past bug; Codex finding). Callers
/// append their own deterministic final tiebreak (gl_key / hit index) for the
/// astronomically rare exact `(rank, ladder)` tie.
///
/// - `y_primary=false`: `rank_score` DESC, then `y_ladder` DESC.
/// - `y_primary=true` (default): `y_ladder` DESC, then `rank_score` DESC.
pub fn collapse_cmp(a_rank: f32, a_ladder: f32, b_rank: f32, b_ladder: f32, y_primary: bool) -> Ordering {
    if y_primary {
        a_ladder.total_cmp(&b_ladder).then(a_rank.total_cmp(&b_rank))
    } else {
        a_rank.total_cmp(&b_rank).then(a_ladder.total_cmp(&b_ladder))
    }
}

/// All glycan-level features attached to a single PSM.
///
/// `glycan_mass` and `backbone_mass` are stored as pre-computed `f64` values
/// so callers do not need to keep a reference to the glycan database.  The
/// canonical way to populate them is:
///
/// ```rust
/// # use andes_glyco::glycan_db::GlycanComp;
/// # use andes_glyco::hybrid::Source;
/// # use andes_glyco::glyco_psm::GlycoPsmKey;
/// let glycan: Option<GlycanComp> = None;
/// let key = GlycoPsmKey {
///     spectrum_idx: 0,
///     glycan_mass: glycan.as_ref().map(|g| g.mass).unwrap_or(0.0),
///     glycan,
///     glycan_source: Source::Db,
///     oxonium_summed_frac: 0.0,
///     n_core_oxonium_ions: 0,
///     y_ladder_intensity_score: 0.0,
///     y_ladder_decoy_score: 0.0,
///     y0y1_anchor_score: 0.0,
///     sialic_consistency: 0.0,
///     core_y_hits: 0,
///     backbone_mass: 0.0,
///     is_transferred: false,
///     transfer_graph_support: 0,
///     transfer_seed_score: 0.0,
///     transfer_rt_delta: 0.0,
///     transfer_ungated: false,
/// };
/// assert_eq!(key.glycan_mass, 0.0);
/// ```
#[derive(Debug, Clone)]
pub struct GlycoPsmKey {
    /// Index of the spectrum this PSM was scored against.
    pub spectrum_idx: usize,
    /// The glycan composition assigned to this PSM, if any.
    pub glycan: Option<GlycanComp>,
    /// Whether the glycan came from the database branch or the de-novo solver.
    pub glycan_source: Source,
    /// Sum of oxonium-ion intensities as a fraction of TIC.
    pub oxonium_summed_frac: f32,
    /// Number of distinct core oxonium ions detected (≤ total in the panel).
    pub n_core_oxonium_ions: u8,
    /// Intensity-weighted score from the core-Y ladder match.
    pub y_ladder_intensity_score: f32,
    /// Glycan-AXIS decoy of `y_ladder_intensity_score`: the same composition's
    /// ladder with intermediate Y-rungs shifted (Y0/Y1 kept). On a true-glycan
    /// spectrum this scores below the target; the gap is what a glycan-decoy PIN
    /// row exposes to Percolator for 2D FDR. 0.0 when no glycan is resolved.
    pub y_ladder_decoy_score: f32,
    /// G2 Y0/Y1 peptide-mass ANCHOR (additive PIN feature): matched intensity of
    /// Y0 (bare peptide) + Y1 (peptide+HexNAc), conditioned on the PEPTIDE mass —
    /// the one glyco feature that discriminates competing peptides at a shared
    /// backbone window. Never folded into the ranking score.
    pub y0y1_anchor_score: f32,
    /// GI-2 composition-conditioned SIALIC consistency (additive PIN feature):
    /// ±NeuAc/NeuGc oxonium signed by whether this glycan claims that sialic —
    /// the one oxonium-derived feature that discriminates glycans of different
    /// sialic content on one spectrum. 0.0 when no glycan is resolved.
    pub sialic_consistency: f32,
    /// Number of core-Y ions matched in the spectrum.
    pub core_y_hits: u8,
    /// Pre-computed monoisotopic mass of the glycan (0.0 when `glycan` is None).
    pub glycan_mass: f64,
    /// Pre-computed monoisotopic mass of the peptide backbone.
    pub backbone_mass: f64,
    /// Cross-spectrum transfer provenance + evidence (additive PIN features).
    /// All inert (false/0) for natively-generated candidates.
    pub is_transferred: bool,
    /// # co-eluting, glycan-delta-linked sibling spectra corroborating this
    /// backbone (the discriminative transfer signal).
    pub transfer_graph_support: u32,
    /// Pass-1 discriminant of the donor seed.
    pub transfer_seed_score: f32,
    /// |RT(acceptor) − RT(seed)| seconds; 0 = perfect co-elution.
    pub transfer_rt_delta: f32,
    /// RT unavailable ⇒ co-elution gate skipped (distrust signal).
    pub transfer_ungated: bool,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::glycan_mass::{HEX, HEXNAC};

    #[test]
    fn y_primary_default_is_yladder_and_rank_forces_off() {
        // DETERMINISTIC RESULT (2026-07-05, PXD025455 Fc3_r1): core-Y-ladder
        // primary collapse = 218 @1% / 90 bb-correct vs b/y-rank primary = 128 /
        // 53 — a +70% win. The earlier "ladder is worse (260→197)" verdict was a
        // non-determinism artifact. So the DEFAULT (env unset) is now yladder.
        assert!(y_primary_from_env(None), "unset -> core-Y ladder primary (new default)");
        assert!(y_primary_from_env(Some("yladder")), "explicit yladder -> true");
        assert!(y_primary_from_env(Some("YLadder")), "case-insensitive");
        // Escape hatch to force the OLD b/y-rank primary (for A/B / regression).
        assert!(!y_primary_from_env(Some("rank")), "rank -> b/y-rank primary");
        assert!(!y_primary_from_env(Some("byrank")), "byrank -> b/y-rank primary");
        // Unknown values fall back to the default (yladder), not a silent error.
        assert!(y_primary_from_env(Some("")), "empty -> default");
        assert!(y_primary_from_env(Some("garbage")), "unknown -> default");
    }

    #[test]
    fn combined_selector_off_by_default_and_only_combined_enables() {
        assert!(!combined_selector_from_env(None), "unset -> legacy bare-rank selector");
        assert!(!combined_selector_from_env(Some("")), "empty -> legacy");
        assert!(!combined_selector_from_env(Some("legacy")), "legacy -> legacy");
        assert!(!combined_selector_from_env(Some("rank")), "rank -> legacy");
        assert!(combined_selector_from_env(Some("combined")), "combined -> on");
        assert!(combined_selector_from_env(Some("Combined")), "case-insensitive");
    }

    /// The core P1b/L4 property: a true backbone whose PEPTIDE-CONDITIONED Y0/Y1
    /// anchor is strong but whose bare b/y rank is BELOW a wrong competitor's must
    /// be selected as top-1 by the combined score — the exact P0b failure where
    /// the recovered true backbone lands in truth_OUTRANKED under the bare-rank
    /// selector.
    #[test]
    fn combined_selector_recovers_true_backbone_that_loses_bare_rank() {
        // Competitor: higher b/y rank, essentially no Y0/Y1 anchor (wrong
        // short-backbone/big-glycan alternative that wins the noisy rank tiebreak).
        let competitor_rank = 22.0_f32;
        let competitor_anchor = 0.05_f32;
        // Truth: LOWER b/y rank (loses the bare-rank selector, R7 median gap ~7)
        // but a strong peptide-mass-conditioned Y0/Y1 anchor.
        let truth_rank = 15.0_f32;
        let truth_anchor = 2.4_f32;

        // Bare-rank selector (the legacy path) picks the WRONG competitor.
        assert!(
            competitor_rank > truth_rank,
            "precondition: truth loses the bare b/y rank (RED without the re-rank)"
        );

        // Combined selector (anchor-only, intensity_llr = 0.0) must flip it to truth.
        let s_comp = glyco_combined_selector_score(competitor_rank, competitor_anchor, 0.0);
        let s_truth = glyco_combined_selector_score(truth_rank, truth_anchor, 0.0);
        assert!(
            s_truth > s_comp,
            "combined selector must rank the strong-anchor true backbone above the \
             high-bare-rank competitor (GREEN): truth={s_truth} comp={s_comp}"
        );
    }

    /// The intensity-model LLR axis also contributes: with equal rank + anchor, a
    /// higher frag/rich-ion LLR wins; passing 0.0 for the LLR is inert (anchor-only
    /// variant), so the two ship stages are separable.
    #[test]
    fn combined_selector_intensity_llr_breaks_ties_and_is_inert_at_zero() {
        let base = glyco_combined_selector_score(10.0, 1.0, 0.0);
        let with_llr = glyco_combined_selector_score(10.0, 1.0, 0.5);
        assert!(with_llr > base, "positive intensity LLR must increase the combined score");
        // Inert at 0.0: anchor-only ship stage is a strict special case.
        assert_eq!(
            glyco_combined_selector_score(7.0, 2.0, 0.0),
            GLYCO_SEL_W_BY * 7.0 + GLYCO_SEL_W_ANCHOR * 2.0,
            "intensity_llr=0.0 must reduce to the anchor-only combined score"
        );
    }

    #[test]
    fn collapse_cmp_default_ranks_by_rank_score_then_ladder() {
        // default: higher rank_score wins even with a lower ladder.
        assert_eq!(collapse_cmp(5.0, 1.0, 3.0, 9.0, false), Ordering::Greater);
        // rank tie → higher ladder wins.
        assert_eq!(collapse_cmp(5.0, 9.0, 5.0, 1.0, false), Ordering::Greater);
    }

    #[test]
    fn collapse_cmp_yprimary_ranks_by_ladder_then_rank_score() {
        // y_primary: higher ladder wins even with a lower rank_score — this is
        // the whole point (recover the correct backbone the noisy b/y rank loses).
        assert_eq!(collapse_cmp(3.0, 9.0, 5.0, 1.0, true), Ordering::Greater);
        // ladder tie → higher rank_score wins.
        assert_eq!(collapse_cmp(5.0, 9.0, 3.0, 9.0, true), Ordering::Greater);
    }

    #[test]
    fn glyco_psm_key_none_glycan_has_zero_glycan_mass() {
        let key = GlycoPsmKey {
            spectrum_idx: 42,
            glycan: None,
            glycan_source: Source::Db,
            oxonium_summed_frac: 0.15,
            n_core_oxonium_ions: 2,
            y_ladder_intensity_score: 0.88,
            y_ladder_decoy_score: 0.2,
            y0y1_anchor_score: 0.4,
            sialic_consistency: 0.1,
            core_y_hits: 4,
            glycan_mass: None::<GlycanComp>.as_ref().map(|g| g.mass).unwrap_or(0.0),
            backbone_mass: 1200.5,
            is_transferred: false,
            transfer_graph_support: 0,
            transfer_seed_score: 0.0,
            transfer_rt_delta: 0.0,
            transfer_ungated: false,
        };
        assert_eq!(key.glycan_mass, 0.0);
        assert!(key.glycan.is_none());
    }

    #[test]
    fn glyco_psm_key_with_real_glycan_has_correct_mass() {
        let glycan = GlycanComp {
            hexnac: 2,
            hex: 3,
            fuc: 0,
            neuac: 0,
            neugc: 0,
            mass: 2.0 * HEXNAC + 3.0 * HEX,
        };
        let expected_mass = glycan.mass;
        let key = GlycoPsmKey {
            spectrum_idx: 7,
            glycan_mass: glycan.mass,
            glycan: Some(glycan),
            glycan_source: Source::DeNovo,
            oxonium_summed_frac: 0.30,
            n_core_oxonium_ions: 3,
            y_ladder_intensity_score: 1.5,
            y_ladder_decoy_score: 0.5,
            y0y1_anchor_score: 0.7,
            sialic_consistency: 0.2,
            core_y_hits: 5,
            backbone_mass: 1500.0,
            is_transferred: false,
            transfer_graph_support: 0,
            transfer_seed_score: 0.0,
            transfer_rt_delta: 0.0,
            transfer_ungated: false,
        };
        assert!((key.glycan_mass - expected_mass).abs() < 1e-6);
        assert!(key.glycan.is_some());
    }

    #[test]
    fn glyco_psm_key_is_clone() {
        let key = GlycoPsmKey {
            spectrum_idx: 1,
            glycan: None,
            glycan_source: Source::Db,
            oxonium_summed_frac: 0.0,
            n_core_oxonium_ions: 0,
            y_ladder_intensity_score: 0.0,
            y_ladder_decoy_score: 0.0,
            y0y1_anchor_score: 0.0,
            sialic_consistency: 0.0,
            core_y_hits: 0,
            glycan_mass: 0.0,
            backbone_mass: 0.0,
            is_transferred: false,
            transfer_graph_support: 0,
            transfer_seed_score: 0.0,
            transfer_rt_delta: 0.0,
            transfer_ungated: false,
        };
        let cloned = key.clone();
        assert_eq!(cloned.spectrum_idx, key.spectrum_idx);
    }

    #[test]
    fn glyco_psm_key_defaults_to_non_transferred() {
        let key = GlycoPsmKey {
            spectrum_idx: 0, glycan: None, glycan_source: Source::Db,
            oxonium_summed_frac: 0.0, n_core_oxonium_ions: 0,
            y_ladder_intensity_score: 0.0, y_ladder_decoy_score: 0.0,
            y0y1_anchor_score: 0.0, sialic_consistency: 0.0, core_y_hits: 0,
            glycan_mass: 0.0, backbone_mass: 0.0,
            is_transferred: false, transfer_graph_support: 0,
            transfer_seed_score: 0.0, transfer_rt_delta: 0.0, transfer_ungated: false,
        };
        assert!(!key.is_transferred);
        assert_eq!(key.transfer_graph_support, 0);
    }
}
