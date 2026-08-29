// Carrier struct for glyco-aware PSM scoring features.
//
// `GlycoPsmKey` bundles all glycan-level evidence gathered for a single PSM
// into one value that can be stored, cloned, and passed to downstream
// re-scorers or PIN writers without re-computing spectra.

use crate::glycan_db::GlycanComp;
use crate::hybrid::Source;
use std::cmp::Ordering;

/// Default balance constant for the `gp` fused selector (leg 2). K scales the
/// glycan core-Y ladder term against the peptide b/y rank-LLR. The `--glyco-gp-k`
/// flag overrides it.
///
/// LOWERED 50 → 10 (2026-07-17 round-2). A data-backed collapse audit showed the
/// `K·ladder` term is computed PER-BACKBONE (identical for every isobaric peptide
/// at a given backbone mass), so at K=50 it dominated the fused score while
/// providing ZERO discrimination between the competing peptides — and its value
/// was empirically ANTI-CORRELATED with correctness (higher on wrong winners). A
/// 2×2 factorial on PXD011533 AI-ETD (6-frac pooled, seed-42 Percolator vs 5088
/// truth): lowering K to 10 while raising c/z to 15 took backbone-correct @1% from
/// 1912 (38%) to 2207 (43%), and combined with the c/z truncation gate to 2453
/// (48%) — z4 32→45%, z5 16→23%, decoy-safe.
pub const GLYCO_GP_K_DEFAULT: f32 = 10.0;

/// Default weight for the core-Y HIT-COUNT term (gp2, leg 2b). Scales the integer
/// `core_y_hits` against rank/ladder. The the reference engine-gap audit found the residual
/// gp-outranked truth backbones have MORE matched-fragment coverage (core-Y hits,
/// matched b/y count) than the wrong winner but a LOWER intensity-weighted
/// RankScore — i.e. andes' rank under-rewards COUNT, exactly what the reference engine's
/// hyperscore (∝ Nb!·Ny!) rewards. Adding this count term recovers them.
/// The `--glyco-gp-j` flag overrides it. Offline top1-by-mass: 318 -> 334 at J=5.
pub const GLYCO_GP_J_DEFAULT: f32 = 5.0;

/// Default weight for the count-rewarding hyperscore term (andes-glyco 2.0 peptide
/// channel `P`). Scales `ln(N_matched!)` against the rank-LLR. Offline blend
/// (rank + hyperscore + ladder + core-Y) top1-by-mass 334 → 344. The `--glyco-gp-h` flag
/// overrides it; 0.0 disables the hyperscore term.
pub const GLYCO_GP_H_DEFAULT: f32 = 1.0;

/// Default weight for the ETD c/z backbone hyperscore term in the gp collapse
/// selector (`--glyco-gp-cz`). Added to the fused score ONLY on ETD/AI-ETD
/// spectra (the per-candidate c/z hyperscore is 0.0 on HCD/CID, so this term is
/// inert on the closed-HCD path → byte-identical). On electron-transfer spectra
/// the intact-glycan c/z ladder is the primary backbone evidence (the labile b/y
/// ladder is sparse for high-charge glycopeptides), so the selector must weight
/// it to pick the true backbone. Offline on PXD011533 Frac1 AI-ETD (467 truth,
/// ceiling 356): gp alone top1-correct 218; `gp + 5·cz` = 250 (+32, z3/z4).
///
/// RAISED 5 → 15 (2026-07-17 round-2). The collapse audit measured c/z as the
/// ONLY per-candidate term that discriminates the true backbone (~12× separation
/// of correct vs wrong winners when it fires), yet at weight 5 it was dominated by
/// the non-discriminating `K·ladder`. Raising c/z to 15 (with K lowered to 10)
/// lets c/z decide the winner on ETD/AI-ETD. Inert on HCD/CID (the per-candidate
/// c/z hyperscore is 0.0 there), so this is byte-identical on the closed-HCD path.
pub const GLYCO_GP_CZ_DEFAULT: f32 = 15.0;

/// The `gp` fused selector score (leg 2): `rank + k·ladder + j·core_y_hits`.
/// Higher is better.
///
/// ADDITIVE fusion of the peptide b/y rank-LLR (`rank`), the glycan core-Y ladder
/// intensity (`ladder`), and the core-Y HIT COUNT (`core_y_hits`, gp2). Unlike the
/// legacy lexicographic [`collapse_cmp`] (ladder primary, so a wrong mass-split's
/// spurious *tiny* ladder edge always overrides truth and the rank tiebreak never
/// fires), this lets a real b/y-rank advantage rescue the true backbone while a
/// *large* ladder difference still wins; the `j·core_y_hits` count term additionally
/// rescues the residual coverage-strong / rank-weak truth backbones the intensity-
/// weighted rank misses (the axis the reference engine's hyperscore rewards). Offline on
/// PXD025455 Fc3_r1 top1-by-mass (present ceiling 375): rank+50·ladder = 318;
/// +5·core_y_hits = 334. Deterministic — fixed K/J, no per-scan normalization,
/// uses only values already available at the collapse. Both collapse sites
/// (glyco_search driver + glyco_pin `select_emitted_hits`) MUST call this with the
/// SAME `rank`/`ladder`/`core_y_hits`/`k`/`j` or their winners diverge (the
/// collapse-parity bug).
pub fn glyco_gp_fused_score(
    rank: f32,
    ladder: f32,
    core_y_hits: f32,
    hyperscore: f32,
    k: f32,
    j: f32,
    h: f32,
) -> f32 {
    rank + k * ladder + j * core_y_hits + h * hyperscore
}

/// Default weight for the matched-b/y-ion term, `M` in
/// [`glyco_gp_fused_score_with_matches`].
///
/// The collapse runs BEFORE feature extraction, so it cannot see the strong score — the
/// best discriminator the engine computes. Measured over a benchmark's reference
/// identifications, the terms it does see rank the correct candidate at median 15
/// (`rank`) and median 44 (`ladder`, the heaviest weight), while the raw count of
/// matched b/y ions ranks it at median 1-2. That count falls out of the hyperscore the
/// selector already evaluates per candidate, so including it costs nothing.
pub const GLYCO_GP_M_DEFAULT: f32 = 0.0;

/// [`glyco_gp_fused_score`] plus `M * matched_b_y_ions`.
///
/// Additive: `m = 0.0` reproduces the previous score exactly, so the term can be
/// switched on by measurement rather than by assumption.
#[allow(clippy::too_many_arguments)]
pub fn glyco_gp_fused_score_with_matches(
    rank: f32,
    ladder: f32,
    core_y_hits: f32,
    hyperscore: f32,
    matched_ions: f32,
    k: f32,
    j: f32,
    h: f32,
    m: f32,
) -> f32 {
    glyco_gp_fused_score(rank, ladder, core_y_hits, hyperscore, k, j, h) + m * matched_ions
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
///     y_hit_frac: 0.0,
///     y_hit_frac_decoy: 0.0,
///     y_ladder_decoy_score: 0.0,
///     partial_glycan_by: 0.0,
///     y0y1_anchor_score: 0.0,
///     sialic_consistency: 0.0,
///     core_y_hits: 0,
///     backbone_mass: 0.0,
///     is_transferred: false,
///     transfer_graph_support: 0,
///     transfer_seed_score: 0.0,
///     transfer_rt_delta: 0.0,
///     transfer_ungated: false,
///     cz_hyperscore: 0.0,
///     cz_intensity: 0.0,
///     cz_explained: 0.0,
///     cz_chance_llr: 0.0,
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
    /// COMPLETENESS of this composition's own Y ladder: the fraction of the Y rungs
    /// the assigned composition PREDICTS that were actually matched, in [0,1].
    ///
    /// Distinct from `y_ladder_intensity_score`, which is an unnormalised sum of
    /// observed intensity and therefore grows with glycan size: a wrong, larger
    /// composition can out-score a right, smaller one by predicting more rungs and
    /// matching a few extra. A fraction cannot be inflated that way — predicting
    /// rungs you cannot support LOWERS it. This is the form the field converged on
    /// (pGlyco's coverage ratios, StrucGP's `matched_branchY_ratio`,
    /// PTM-Shepherd's explicit miss penalty).
    ///
    /// Measured motivation: 96.9% of our decoy winners sit at a DIFFERENT backbone
    /// mass than the truth, i.e. the failing decision is which mass split /
    /// composition to believe — exactly what a completeness fraction scores.
    pub y_hit_frac: f32,
    /// Glycan-AXIS decoy of `y_hit_frac`: same composition, interior rungs shifted.
    /// A correct composition should hold its completeness while its shifted twin
    /// collapses; the gap is what a glycan-decoy PIN row exposes to Percolator.
    pub y_hit_frac_decoy: f32,
    /// PARTIAL-GLYCAN b/y evidence (idea B): matched intensity of peptide b/y fragments
    /// bearing the innermost core glycan (b_i/y_i + {HexNAc, 2HexNAc, ...}). Unlike the
    /// mass-based Y-ladder, this is SEQUENCE-specific → discriminates the true backbone
    /// from mass-preserving decoys, the evidence weak large/high-charge glycopeptides
    /// lack. Additive PIN feature. 0.0 when no backbone is resolved.
    pub partial_glycan_by: f32,
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
    /// ETD c/z backbone hyperscore (additive PIN feature `CzHyperscore`), computed
    /// only on ETD/AI-ETD spectra: `ln(N_c!) + ln(N_z!)` over distinct matched
    /// c/z ions of the glycopeptide backbone (glycan on glycosite-spanning
    /// fragments). 0.0 on collisional (HCD/CID) spectra — the orthogonal
    /// electron-transfer evidence that recovers the high-charge glycopeptides the
    /// labile-glycan b/y ladder misses. On ETD spectra this same c/z hyperscore
    /// ALSO contributes to the per-scan collapse selector (weighted by
    /// `--glyco-gp-cz`); the peptide `rank_score`/`RawScore` are unchanged.
    pub cz_hyperscore: f32,
    /// Fraction of base-peak intensity captured by matched glycopeptide-aware c/z
    /// ions (additive PIN feature `CzIntensity`), ETD/AI-ETD only (0.0 on HCD/CID).
    /// The INTENSITY companion to `cz_hyperscore` (which is count-only) — reference
    /// glyco scorers weight the matched-fragment intensity that the count-only
    /// hyperscore discards (round-4 "intensity blindness" audit). PIN feature only;
    /// not (yet) in the collapse selector.
    pub cz_intensity: f32,
    /// c/z prior-weighted EXPLAINED fraction (additive PIN `CzExplained`; ETD only,
    /// else 0). See `cz_structure_features`.
    pub cz_explained: f32,
    /// c/z local-noise chance LLR (additive PIN `CzChanceLlr`; ETD only, else 0).
    pub cz_chance_llr: f32,
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
    fn gp_fused_score_lets_rank_rescue_truth_but_a_large_ladder_still_wins() {
        // Leg-2 mechanism (the exact P0 failure): a wrong mass-split with a
        // SPURIOUS TINY ladder edge beats truth under the legacy ladder-primary
        // collapse, even though truth has the stronger b/y rank.
        let (k, j, h) = (GLYCO_GP_K_DEFAULT, 0.0, 0.0); // isolate rank/ladder (K=10)
        let truth = glyco_gp_fused_score(15.0, 0.05, 0.0, 0.0, k, j, h); // 15 + 0.5 = 15.5
        let wrong = glyco_gp_fused_score(2.0, 0.06, 0.0, 0.0, k, j, h); //  2 + 0.6 =  2.6
        // Legacy y_primary would pick `wrong` (0.06 > 0.05); gp fusion rescues truth.
        assert!(collapse_cmp(15.0, 0.05, 2.0, 0.06, true) == Ordering::Less);
        assert!(truth > wrong, "a real b/y-rank advantage rescues truth under gp");
        // A GENUINELY large ladder difference (strong glycan-Y evidence) still wins —
        // at the round-2 K=10 the ladder no longer dominates a small edge (by design),
        // but a real ladder gap does. (2 + 10·2.0 = 22 > 15.5.)
        let strong_glycan = glyco_gp_fused_score(2.0, 2.0, 0.0, 0.0, k, j, h); // 2 + 20 = 22
        assert!(strong_glycan > truth, "a large ladder difference still wins under gp");
    }

    #[test]
    fn gp2_core_y_hit_count_term_rescues_coverage_strong_truth() {
        // gp2 (leg 2b): truth loses on rank+ladder alone but has MORE core-Y hits;
        // the j·core_y_hits count term flips it (the the reference engine-hyperscore axis).
        let (k, j, h) = (GLYCO_GP_K_DEFAULT, GLYCO_GP_J_DEFAULT, 0.0); // 50, 5, no hyperscore
        // truth: rank 8, ladder 0.1, core-Y 6 ; winner: rank 12, ladder 0.1, core-Y 2
        let truth = glyco_gp_fused_score(8.0, 0.1, 6.0, 0.0, k, j, h); // 8 + 5 + 30 = 43
        let winner = glyco_gp_fused_score(12.0, 0.1, 2.0, 0.0, k, j, h); // 12 + 5 + 10 = 27
        assert!(truth > winner, "core-Y hit count rescues the coverage-strong truth");
        // Without the count term (j=0) the higher-rank winner would win.
        assert!(glyco_gp_fused_score(8.0, 0.1, 6.0, 0.0, k, 0.0, h)
            < glyco_gp_fused_score(12.0, 0.1, 2.0, 0.0, k, 0.0, h));
    }

    #[test]
    fn glyco_gp_weights_default_and_fusion_is_monotone() {
        assert_eq!(GLYCO_GP_K_DEFAULT, 10.0);
        assert_eq!(GLYCO_GP_J_DEFAULT, 5.0);
        assert_eq!(GLYCO_GP_H_DEFAULT, 1.0);
        assert_eq!(GLYCO_GP_CZ_DEFAULT, 15.0);
        // Monotone in each term (rank, ladder, core-Y, hyperscore).
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
            y_ladder_decoy_score: 0.2,
            y_hit_frac: 0.0,
            y_hit_frac_decoy: 0.0,
            partial_glycan_by: 0.0,
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
            y_hit_frac: 0.0,
            y_hit_frac_decoy: 0.0,
            partial_glycan_by: 0.0,
            y0y1_anchor_score: 0.7,
            sialic_consistency: 0.2,
            core_y_hits: 5,
            backbone_mass: 1500.0,
            is_transferred: false,
            transfer_graph_support: 0,
            transfer_seed_score: 0.0,
            transfer_rt_delta: 0.0,
            transfer_ungated: false,
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
            y_hit_frac_decoy: 0.0,
            y_ladder_decoy_score: 0.0,
            partial_glycan_by: 0.0,
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
            cz_hyperscore: 0.0,
            cz_intensity: 0.0,
            cz_explained: 0.0,
            cz_chance_llr: 0.0,
        };
        let cloned = key.clone();
        assert_eq!(cloned.spectrum_idx, key.spectrum_idx);
    }

    #[test]
    fn glyco_psm_key_defaults_to_non_transferred() {
        let key = GlycoPsmKey {
            spectrum_idx: 0, glycan: None, glycan_source: Source::Db,
            oxonium_summed_frac: 0.0, n_core_oxonium_ions: 0,
            y_ladder_intensity_score: 0.0,
            y_hit_frac: 0.0,
            y_hit_frac_decoy: 0.0, y_ladder_decoy_score: 0.0, partial_glycan_by: 0.0,
            y0y1_anchor_score: 0.0, sialic_consistency: 0.0, core_y_hits: 0,
            glycan_mass: 0.0, backbone_mass: 0.0,
            is_transferred: false, transfer_graph_support: 0,
            transfer_seed_score: 0.0, transfer_rt_delta: 0.0, transfer_ungated: false,
            cz_hyperscore: 0.0,
            cz_intensity: 0.0,
            cz_explained: 0.0,
            cz_chance_llr: 0.0,
        };
        assert!(!key.is_transferred);
        assert_eq!(key.transfer_graph_support, 0);
    }
}
