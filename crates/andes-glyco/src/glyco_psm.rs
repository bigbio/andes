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
/// recovered ID. Default (`false`) selects by the peptide b/y `rank_score`
/// first. `ANDES_GLYCO_SELECT=yladder` selects by the core-Y ladder first.
///
/// EMPIRICAL RESULT (2026-07-04, PXD025455 Fc3_r1): ladder-primary is WORSE —
/// 260→197 @1% FDR, 101→88 backbone-correct, PIN 1928→1284 rows. It promotes
/// de-novo / mono-offset backbones that carry strong RAW core-Y but no
/// enumerated glycan into the per-scan winner slot, where the enumerated-only
/// filter then drops them, losing the scan. So b/y `rank_score` is NOT just
/// noise; the default (rank primary, ladder tiebreak) is the better rule. Kept
/// OFF-by-default as the scaffold for a future LEARNED b/y+Y combination (SP-B)
/// — the real fix for the ranking bottleneck, not a hard primary swap.
pub fn y_primary_selection() -> bool {
    std::env::var("ANDES_GLYCO_SELECT")
        .map(|v| v.eq_ignore_ascii_case("yladder"))
        .unwrap_or(false)
}

/// Total order for the top-1-per-scan collapse: `max_by(collapse_cmp(...))`
/// yields the emitted winner. This ordering is the SINGLE SOURCE OF TRUTH shared
/// by the driver's pre-feature reduction (glyco_search) and the PIN writer's
/// `select_emitted_hits` — they MUST agree, or a scan's driver-emitted winner
/// and PIN-written winner diverge (a real past bug; Codex finding). Callers
/// append their own deterministic final tiebreak (gl_key / hit index) for the
/// astronomically rare exact `(rank, ladder)` tie.
///
/// - `y_primary=false` (default): `rank_score` DESC, then `y_ladder` DESC.
/// - `y_primary=true`: `y_ladder` DESC, then `rank_score` DESC.
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
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::glycan_mass::{HEX, HEXNAC};

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
        };
        let cloned = key.clone();
        assert_eq!(cloned.spectrum_idx, key.spectrum_idx);
    }
}
