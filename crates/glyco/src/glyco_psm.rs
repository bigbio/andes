use crate::glycan_db::GlycanComp;
use crate::hybrid::Source;
use std::cmp::Ordering;

pub const GLYCO_GP_K_DEFAULT: f32 = 10.0;

pub const GLYCO_GP_J_DEFAULT: f32 = 5.0;

pub const GLYCO_GP_H_DEFAULT: f32 = 1.0;

pub const GLYCO_GP_ISO_DEFAULT: f32 = 1.0;

pub const GLYCO_GP_CZ_DEFAULT: f32 = 15.0;

pub fn glyco_gp_fused_score(
    rank: f32,
    ladder: f32,
    core_y_hits: f32,
    hyperscore: f32,
    k: f32,
    j: f32,
    h: f32,
) -> f32 {
    glyco_gp_fused_score_iso(rank, ladder, core_y_hits, hyperscore, k, j, h, 0, 0.0)
}

#[allow(clippy::too_many_arguments)]
pub fn glyco_gp_fused_score_iso(
    rank: f32,
    ladder: f32,
    core_y_hits: f32,
    hyperscore: f32,
    k: f32,
    j: f32,
    h: f32,
    isotope_offset: i8,
    iso: f32,
) -> f32 {
    rank + k * ladder + j * core_y_hits + h * hyperscore
        - iso * (isotope_offset.unsigned_abs() as f32)
}

pub fn collapse_cmp(
    a_rank: f32,
    a_ladder: f32,
    b_rank: f32,
    b_ladder: f32,
    y_primary: bool,
) -> Ordering {
    if y_primary {
        a_ladder
            .total_cmp(&b_ladder)
            .then(a_rank.total_cmp(&b_rank))
    } else {
        a_rank
            .total_cmp(&b_rank)
            .then(a_ladder.total_cmp(&b_ladder))
    }
}

/// All glycan-level features attached to a single PSM.
///
/// `glycan_mass` and `backbone_mass` are stored as pre-computed `f64` values
/// so callers do not need to keep a reference to the glycan database.  The
/// canonical way to populate them is:
///
/// ```rust
/// # use glyco::glycan_db::GlycanComp;
/// # use glyco::hybrid::Source;
/// # use glyco::glyco_psm::GlycoPsmKey;
/// let glycan: Option<GlycanComp> = None;
/// let key = GlycoPsmKey {
///     spectrum_idx: 0,
///     glycan_mass: glycan.as_ref().map(|g| g.mass).unwrap_or(0.0),
///     glycan,
///     glycan_source: Source::Db,
///     oxonium_summed_frac: 0.0,
///     n_core_oxonium_ions: 0,
///     y_ladder_intensity_score: 0.0,
///     y_hit_frac: 0.0,
///     partial_glycan_by: 0.0,
///     y0y1_anchor_score: 0.0,
///     sialic_consistency: 0.0,
///     core_y_hits: 0,
///     backbone_mass: 0.0,
///     cz_hyperscore: 0.0,
///     cz_intensity: 0.0,
///     cz_explained: 0.0,
///     cz_chance_llr: 0.0,
/// };
/// assert_eq!(key.glycan_mass, 0.0);
/// ```
#[derive(Debug, Clone)]
pub struct GlycoPsmKey {
    pub spectrum_idx: usize,
    pub glycan: Option<GlycanComp>,
    pub glycan_source: Source,
    pub oxonium_summed_frac: f32,
    pub n_core_oxonium_ions: u8,
    pub y_ladder_intensity_score: f32,
    pub y_hit_frac: f32,
    pub partial_glycan_by: f32,
    pub y0y1_anchor_score: f32,
    pub sialic_consistency: f32,
    pub core_y_hits: u8,
    pub glycan_mass: f64,
    pub backbone_mass: f64,
    pub cz_hyperscore: f32,
    pub cz_intensity: f32,
    pub cz_explained: f32,
    pub cz_chance_llr: f32,
}

#[cfg(test)]
mod isotope_penalty_tests {
    use super::*;

    #[test]
    fn the_monoisotopic_hypothesis_wins_a_near_tie() {
        let k = GLYCO_GP_K_DEFAULT;
        let truth =
            glyco_gp_fused_score_iso(40.0, 12.00, 3.0, 0.0, k, 0.0, 0.0, 0, GLYCO_GP_ISO_DEFAULT);
        let swap =
            glyco_gp_fused_score_iso(40.0, 12.02, 3.0, 0.0, k, 0.0, 0.0, 1, GLYCO_GP_ISO_DEFAULT);
        assert!(
            truth > swap,
            "the M+1 candidate still wins a 0.2-point ladder edge: truth {truth}, swap {swap}"
        );
    }

    #[test]
    fn real_evidence_still_beats_the_penalty() {
        let (k, j) = (GLYCO_GP_K_DEFAULT, GLYCO_GP_J_DEFAULT);
        let weak_mono =
            glyco_gp_fused_score_iso(40.0, 12.0, 2.0, 0.0, k, j, 0.0, 0, GLYCO_GP_ISO_DEFAULT);
        let strong_m1 =
            glyco_gp_fused_score_iso(40.0, 12.0, 3.0, 0.0, k, j, 0.0, 1, GLYCO_GP_ISO_DEFAULT);
        assert!(
            strong_m1 > weak_mono,
            "an extra core-Y hit must outweigh one isotope step"
        );
    }

    #[test]
    fn zero_weight_reproduces_the_old_ordering() {
        let k = GLYCO_GP_K_DEFAULT;
        let a = glyco_gp_fused_score_iso(40.0, 12.00, 3.0, 0.0, k, 0.0, 0.0, 0, 0.0);
        let b = glyco_gp_fused_score_iso(40.0, 12.02, 3.0, 0.0, k, 0.0, 0.0, 1, 0.0);
        assert!(b > a, "at weight 0 the ladder edge must still decide");
        assert_eq!(
            a,
            glyco_gp_fused_score(40.0, 12.00, 3.0, 0.0, k, 0.0, 0.0),
            "the seven-argument form must equal the penalty-free score"
        );
    }

    #[test]
    fn the_penalty_scales_with_the_offset() {
        let k = GLYCO_GP_K_DEFAULT;
        let one =
            glyco_gp_fused_score_iso(40.0, 12.0, 3.0, 0.0, k, 0.0, 0.0, 1, GLYCO_GP_ISO_DEFAULT);
        let two =
            glyco_gp_fused_score_iso(40.0, 12.0, 3.0, 0.0, k, 0.0, 0.0, 2, GLYCO_GP_ISO_DEFAULT);
        let neg =
            glyco_gp_fused_score_iso(40.0, 12.0, 3.0, 0.0, k, 0.0, 0.0, -1, GLYCO_GP_ISO_DEFAULT);
        assert!(one > two, "M+2 must cost more than M+1");
        assert_eq!(one, neg, "a -1 offset is the same size of correction as +1");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn collapse_cmp_default_ranks_by_rank_score_then_ladder() {
        assert_eq!(collapse_cmp(5.0, 1.0, 3.0, 9.0, false), Ordering::Greater);
        assert_eq!(collapse_cmp(5.0, 9.0, 5.0, 1.0, false), Ordering::Greater);
    }

    #[test]
    fn collapse_cmp_yprimary_ranks_by_ladder_then_rank_score() {
        assert_eq!(collapse_cmp(3.0, 9.0, 5.0, 1.0, true), Ordering::Greater);
        assert_eq!(collapse_cmp(5.0, 9.0, 3.0, 9.0, true), Ordering::Greater);
    }

    #[test]
    fn gp_fused_score_lets_rank_rescue_truth_but_a_large_ladder_still_wins() {
        let (k, j, h) = (GLYCO_GP_K_DEFAULT, 0.0, 0.0);
        let truth = glyco_gp_fused_score(15.0, 0.05, 0.0, 0.0, k, j, h);
        let wrong = glyco_gp_fused_score(2.0, 0.06, 0.0, 0.0, k, j, h);
        assert!(collapse_cmp(15.0, 0.05, 2.0, 0.06, true) == Ordering::Less);
        assert!(
            truth > wrong,
            "a real b/y-rank advantage rescues truth under gp"
        );
        let strong_glycan = glyco_gp_fused_score(2.0, 2.0, 0.0, 0.0, k, j, h);
        assert!(
            strong_glycan > truth,
            "a large ladder difference still wins under gp"
        );
    }

    #[test]
    fn gp2_core_y_hit_count_term_rescues_coverage_strong_truth() {
        let (k, j, h) = (GLYCO_GP_K_DEFAULT, GLYCO_GP_J_DEFAULT, 0.0);
        let truth = glyco_gp_fused_score(8.0, 0.1, 6.0, 0.0, k, j, h);
        let winner = glyco_gp_fused_score(12.0, 0.1, 2.0, 0.0, k, j, h);
        assert!(
            truth > winner,
            "core-Y hit count rescues the coverage-strong truth"
        );
        assert!(
            glyco_gp_fused_score(8.0, 0.1, 6.0, 0.0, k, 0.0, h)
                < glyco_gp_fused_score(12.0, 0.1, 2.0, 0.0, k, 0.0, h)
        );
    }

    #[test]
    fn glyco_gp_weights_default_and_fusion_is_monotone() {
        assert_eq!(GLYCO_GP_K_DEFAULT, 10.0);
        assert_eq!(GLYCO_GP_J_DEFAULT, 5.0);
        assert_eq!(GLYCO_GP_H_DEFAULT, 1.0);
        assert_eq!(GLYCO_GP_CZ_DEFAULT, 15.0);
        let f = |rk, yl, cy, hs| glyco_gp_fused_score(rk, yl, cy, hs, 50.0, 5.0, 1.0);
        assert!(f(10.0, 0.2, 0.0, 0.0) > f(10.0, 0.1, 0.0, 0.0));
        assert!(f(11.0, 0.1, 0.0, 0.0) > f(10.0, 0.1, 0.0, 0.0));
        assert!(f(10.0, 0.1, 3.0, 0.0) > f(10.0, 0.1, 1.0, 0.0));
        assert!(f(10.0, 0.1, 0.0, 5.0) > f(10.0, 0.1, 0.0, 2.0));
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
            y_hit_frac: 0.0,
            partial_glycan_by: 0.0,
            y0y1_anchor_score: 0.4,
            sialic_consistency: 0.1,
            core_y_hits: 4,
            glycan_mass: None::<GlycanComp>.as_ref().map(|g| g.mass).unwrap_or(0.0),
            backbone_mass: 1200.5,
            cz_hyperscore: 0.0,
            cz_intensity: 0.0,
            cz_explained: 0.0,
            cz_chance_llr: 0.0,
        };
        assert_eq!(key.glycan_mass, 0.0);
        assert!(key.glycan.is_none());
    }

    #[test]
    fn glyco_psm_key_with_real_glycan_has_correct_mass() {
        let glycan = GlycanComp::from_counts(2, 3, 0, 0, 0, 0);
        let expected_mass = glycan.mass;
        let key = GlycoPsmKey {
            spectrum_idx: 7,
            glycan_mass: glycan.mass,
            glycan: Some(glycan),
            glycan_source: Source::DeNovo,
            oxonium_summed_frac: 0.30,
            n_core_oxonium_ions: 3,
            y_ladder_intensity_score: 1.5,
            y_hit_frac: 0.0,
            partial_glycan_by: 0.0,
            y0y1_anchor_score: 0.7,
            sialic_consistency: 0.2,
            core_y_hits: 5,
            backbone_mass: 1500.0,
            cz_hyperscore: 0.0,
            cz_intensity: 0.0,
            cz_explained: 0.0,
            cz_chance_llr: 0.0,
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
            y_hit_frac: 0.0,
            partial_glycan_by: 0.0,
            y0y1_anchor_score: 0.0,
            sialic_consistency: 0.0,
            core_y_hits: 0,
            glycan_mass: 0.0,
            backbone_mass: 0.0,
            cz_hyperscore: 0.0,
            cz_intensity: 0.0,
            cz_explained: 0.0,
            cz_chance_llr: 0.0,
        };
        let cloned = key.clone();
        assert_eq!(cloned.spectrum_idx, key.spectrum_idx);
    }
}
