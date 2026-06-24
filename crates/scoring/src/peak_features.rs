//! Peptide-AGNOSTIC per-peak signal/noise features (the GBDT input contract).
//!
//! THIS LIST IS THE SINGLE SOURCE OF TRUTH. The Python training extractor
//! (`training/gbdt/feature_spec.py`) must mirror `FEATURE_NAMES` in the same
//! order; a cross-language parity test enforces it. Every feature is computable
//! once per spectrum (no
//! candidate peptide is consulted), so the GBDT is evaluated once per peak at
//! spectrum-prep, never in the inner candidate loop.

/// Ordered feature names. Index in this array == feature index used by the
/// GBDT tree splits. DO NOT reorder without retraining + retranscoding.
pub const FEATURE_NAMES: [&str; 18] = [
    "log_intensity",           // 0  ln(intensity)
    "intensity_over_basepeak", // 1  intensity / max intensity in scan
    "intensity_over_tic",      // 2  intensity / summed kept intensity
    "global_rank_frac",        // 3  (rank-1) / kept_count
    "local_rank_frac",         // 4  (rank within ±window) / count in window
    "is_top1_in_window",       // 5  1.0 if most intense in ±window else 0.0
    "is_top3_in_window",       // 6  1.0 if among top-3 in ±window else 0.0
    "mz",                      // 7  observed m/z
    "mz_frac_of_precursor",    // 8  mz / precursor_mz
    "local_peak_density",      // 9  peaks per Da in ±window
    "spacing_left",            // 10 mz - previous peak mz (SENTINEL if none)
    "spacing_right",           // 11 next peak mz - mz (SENTINEL if none)
    "mass_defect",             // 12 mz - floor(mz)
    "has_isotope_plus1",       // 13 peak at mz + 1.00235/charge within tol
    "has_isotope_minus1",      // 14 peak at mz - 1.00235/charge within tol
    "has_complement",          // 15 peak at (M + 2*PROTON - mz) within tol
    "has_h2o_loss_partner",    // 16 peak at mz - 18.010565/charge within tol
    "has_nh3_loss_partner",    // 17 peak at mz - 17.026549/charge within tol
];

pub const N_FEATURES: usize = FEATURE_NAMES.len();

/// Sentinel for spacing when there is no neighbor on that side (Da). Large so
/// trees can isolate edge peaks.
pub const SPACING_SENTINEL: f32 = 1000.0;

/// Half-width (Da) of the local window for local-rank / density / top-K
/// features. MUST equal `training/gbdt/feature_spec.py::FEATURE_WINDOW_DA` and
/// the value passed by `build_dataset.py` at training time — diverging would
/// silently break feature parity (A2.10 guards it).
pub const FEATURE_WINDOW_DA: f64 = 50.0;

use model::mass::{ISOTOPE, PROTON};
use model::tolerance::Tolerance;

// H2O kept as a local literal (NOT model::mass::H2O): the model crate computes
// H2O as H*2+O which is NOT bit-equal to this literal, and the Python feature
// extractor (feature_spec.py) uses the literal 18.010565 — keeping the literal
// here preserves exact cross-language feature parity (A2.10).
const H2O: f64 = 18.010565;
const NH3: f64 = 17.026549; // no model::mass equivalent

/// Per-spectrum constants needed to compute the features.
pub struct PeakFeatureCtx {
    pub precursor_mz: f64,
    pub charge: u8,
    /// (precursor_mz - PROTON) * charge — the observed neutral parent mass.
    pub parent_neutral_mass: f64,
    /// Summed intensity of the active (kept) peaks.
    pub total_intensity: f64,
    /// Maximum intensity among the active peaks.
    pub base_peak_intensity: f32,
    /// Half-width (Da) of the local window for local-rank/density/top-K feats.
    pub window_da: f64,
    /// Tolerance (Da) for partner-peak existence flags (isotope/complement/loss).
    pub match_tol_da: f64,
}

impl PeakFeatureCtx {
    /// THE single source of truth for the per-spectrum feature context. Both the
    /// scoring path (`ScoredSpectrum::new`) and the training dataset builder MUST
    /// build the context through this constructor, so train-time and inference-time
    /// features are identical by construction.
    ///
    /// `active_peaks` is the peak list features are extracted from (post precursor
    /// filter / post deconvolution). `base_peak_intensity` and `total_intensity`
    /// are BOTH measured on this same list so the intensity-normalized features
    /// share one population. `match_tol_da` is a per-spectrum SCALAR tolerance
    /// (not per-peak adaptive): `mme.as_da(parent_neutral_mass)`, evaluated once at
    /// the precursor mass — the partner-peak existence flags (isotope/complement/
    /// loss) use this single tolerance for every peak. This is a deliberate
    /// definition of the feature contract; if it is ever changed, it changes for
    /// BOTH training and inference because they share this function.
    pub fn for_spectrum(
        precursor_mz: f64,
        charge: u8,
        parent_neutral_mass: f64,
        active_peaks: &[(f64, f32)],
        mme: &Tolerance,
    ) -> Self {
        let base_peak_intensity = active_peaks.iter().map(|&(_, i)| i).fold(0.0_f32, f32::max);
        let total_intensity: f64 = active_peaks.iter().map(|&(_, i)| i as f64).sum();
        PeakFeatureCtx {
            precursor_mz,
            charge,
            parent_neutral_mass,
            total_intensity,
            base_peak_intensity,
            window_da: FEATURE_WINDOW_DA,
            match_tol_da: mme.as_da(parent_neutral_mass.max(1.0)),
        }
    }
}

/// Compute one feature vector per peak. `peaks` MUST be ascending by m/z and
/// aligned 1:1 with `ranks` (rank 1 = most intense; `u32::MAX` = filtered out).
/// Filtered-out peaks still get a row (so indices stay aligned) but their
/// rank-based features use the kept-count denominator.
pub fn extract_peak_features(
    peaks: &[(f64, f32)],
    ranks: &[u32],
    ctx: &PeakFeatureCtx,
) -> Vec<[f32; N_FEATURES]> {
    debug_assert_eq!(peaks.len(), ranks.len(), "peaks and ranks must be 1:1 aligned");
    let n = peaks.len();
    let kept_count = ranks.iter().filter(|&&r| r != u32::MAX).count().max(1);
    let tol = ctx.match_tol_da;
    let mut out = Vec::with_capacity(n);

    for i in 0..n {
        let (mz, intensity) = peaks[i];
        let rank = ranks[i];
        let mut f = [0.0_f32; N_FEATURES];

        f[0] = intensity.max(1e-6).ln();
        f[1] = if ctx.base_peak_intensity > 0.0 {
            intensity / ctx.base_peak_intensity
        } else {
            0.0
        };
        f[2] = if ctx.total_intensity > 0.0 {
            (intensity as f64 / ctx.total_intensity) as f32
        } else {
            0.0
        };
        f[3] = if rank == u32::MAX {
            1.0
        } else {
            rank.saturating_sub(1) as f32 / kept_count as f32
        };

        // Local window: peaks with m/z in [mz - window, mz + window].
        let lo = peaks.partition_point(|&(m, _)| m < mz - ctx.window_da);
        let hi = peaks.partition_point(|&(m, _)| m <= mz + ctx.window_da);
        let win = &peaks[lo..hi];
        let win_count = win.len().max(1);
        // local rank = #peaks in window strictly more intense than this one.
        let more_intense = win.iter().filter(|&&(_, pint)| pint > intensity).count();
        f[4] = more_intense as f32 / win_count as f32;
        // Ties (equal intensity) are counted with strict >, so tied peaks all share the top flag; the Python extractor uses the same strict >.
        f[5] = if more_intense == 0 { 1.0 } else { 0.0 };
        f[6] = if more_intense < 3 { 1.0 } else { 0.0 };

        f[7] = mz as f32;
        f[8] = if ctx.precursor_mz > 0.0 {
            (mz / ctx.precursor_mz) as f32
        } else {
            0.0
        };
        f[9] = (win_count as f64 / (2.0 * ctx.window_da)) as f32;

        f[10] = if i > 0 {
            (mz - peaks[i - 1].0) as f32
        } else {
            SPACING_SENTINEL
        };
        f[11] = if i + 1 < n {
            (peaks[i + 1].0 - mz) as f32
        } else {
            SPACING_SENTINEL
        };
        f[12] = (mz - mz.floor()) as f32;

        let z = ctx.charge.max(1) as f64;
        f[13] = has_peak(peaks, mz + ISOTOPE / z, tol);
        f[14] = has_peak(peaks, mz - ISOTOPE / z, tol);
        // Complement of a singly-charged fragment: b_i + y_(n-i) = M + 2*PROTON.
        f[15] = has_peak(peaks, ctx.parent_neutral_mass + 2.0 * PROTON - mz, tol);
        f[16] = has_peak(peaks, mz - H2O / z, tol);
        f[17] = has_peak(peaks, mz - NH3 / z, tol);

        out.push(f);
    }
    out
}

/// 1.0 if any peak lies within `tol` Da of `target`, else 0.0. `peaks` ascending.
fn has_peak(peaks: &[(f64, f32)], target: f64, tol: f64) -> f32 {
    if target <= 0.0 {
        return 0.0;
    }
    let lo = peaks.partition_point(|&(m, _)| m < target - tol);
    if lo < peaks.len() && peaks[lo].0 <= target + tol {
        1.0
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn feature_names_count_matches() {
        assert_eq!(N_FEATURES, 18, "feature set changed — update the Python contract (feature_spec.py) and A2.10 fixtures");
    }

    #[test]
    fn extracts_expected_features_on_tiny_scan() {
        // Three peaks; precursor m/z 500.0, charge 2 → neutral ~997.985.
        // peaks ascending by m/z: (100.0, 10.0), (200.05, 100.0), (300.0, 50.0)
        let peaks = vec![(100.0_f64, 10.0_f32), (200.05, 100.0), (300.0, 50.0)];
        // ranks by intensity desc: idx1 (100)→1, idx2 (50)→2, idx0 (10)→3
        let ranks = vec![3u32, 1, 2];
        let ctx = PeakFeatureCtx {
            precursor_mz: 500.0,
            charge: 2,
            parent_neutral_mass: (500.0 - PROTON) * 2.0,
            total_intensity: 160.0,
            base_peak_intensity: 100.0,
            window_da: 50.0,
            match_tol_da: 0.5,
        };
        let f = extract_peak_features(&peaks, &ranks, &ctx);
        assert_eq!(f.len(), 3);
        // peak 1 (mz 200.05, intensity 100, global rank 1):
        let p1 = &f[1];
        // intensity_over_basepeak = 100/100 = 1.0
        assert!((p1[idx("intensity_over_basepeak")] - 1.0).abs() < 1e-6);
        // intensity_over_tic = 100/160
        assert!((p1[idx("intensity_over_tic")] - (100.0 / 160.0)).abs() < 1e-6);
        // global_rank_frac = (1-1)/3 = 0.0
        assert!((p1[idx("global_rank_frac")] - 0.0).abs() < 1e-6);
        // mass_defect = 200.05 - 200.0 = 0.05
        assert!((p1[idx("mass_defect")] - 0.05).abs() < 1e-4);
        // mz_frac_of_precursor = 200.05/500.0
        assert!((p1[idx("mz_frac_of_precursor")] - (200.05 / 500.0)).abs() < 1e-6);
        // is_top1_in_window: within ±50 Da of 200.05 the peaks are {200.05}
        // (100.0 is >50 away, 300.0 is >50 away) → it is the top → 1.0
        assert!((p1[idx("is_top1_in_window")] - 1.0).abs() < 1e-6);

        // spacing sentinels at the ends:
        assert!((f[0][idx("spacing_left")] - SPACING_SENTINEL).abs() < 1e-6);  // first peak: no left neighbor
        assert!((f[2][idx("spacing_right")] - SPACING_SENTINEL).abs() < 1e-6); // last peak: no right neighbor
        // spacing_right of peak 0 = 200.05 - 100.0 = 100.05
        assert!((f[0][idx("spacing_right")] - 100.05).abs() < 1e-3);
        // a binary partner flag: no isotope/complement/loss partner exists in this sparse scan → 0.0
        assert!((f[1][idx("has_isotope_plus1")] - 0.0).abs() < 1e-6);
        assert!((f[1][idx("has_complement")] - 0.0).abs() < 1e-6);
    }

    #[test]
    fn filtered_peak_gets_max_rank_frac() {
        let peaks = vec![(150.0_f64, 5.0_f32), (250.0, 80.0)];
        let ranks = vec![u32::MAX, 1]; // first peak filtered out
        let ctx = PeakFeatureCtx {
            precursor_mz: 400.0, charge: 2,
            parent_neutral_mass: (400.0 - PROTON) * 2.0,
            total_intensity: 80.0, base_peak_intensity: 80.0,
            window_da: 50.0, match_tol_da: 0.5,
        };
        let f = extract_peak_features(&peaks, &ranks, &ctx);
        let gi = FEATURE_NAMES.iter().position(|&n| n == "global_rank_frac").unwrap();
        assert!((f[0][gi] - 1.0).abs() < 1e-6, "filtered peak must get rank_frac 1.0");
    }

    #[test]
    fn for_spectrum_uses_active_list_for_both_denominators() {
        use model::tolerance::Tolerance;
        let peaks = vec![(100.0_f64, 10.0_f32), (200.0, 40.0), (300.0, 50.0)];
        let ctx = PeakFeatureCtx::for_spectrum(
            500.0, 2, (500.0 - 1.00727649) * 2.0, &peaks, &Tolerance::Da(0.5),
        );
        assert!((ctx.base_peak_intensity - 50.0).abs() < 1e-6);
        assert!((ctx.total_intensity - 100.0).abs() < 1e-9);
        assert!((ctx.window_da - FEATURE_WINDOW_DA).abs() < 1e-9);
        assert!((ctx.match_tol_da - 0.5).abs() < 1e-9); // Da tolerance is m/z-independent
    }

    // Test-only helper to look up a feature column by name.
    fn idx(name: &str) -> usize {
        FEATURE_NAMES.iter().position(|&n| n == name).expect("feature exists")
    }
}
