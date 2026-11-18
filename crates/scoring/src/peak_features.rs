pub const FEATURE_NAMES: [&str; 18] = [
    "log_intensity",
    "intensity_over_basepeak",
    "intensity_over_tic",
    "global_rank_frac",
    "local_rank_frac",
    "is_top1_in_window",
    "is_top3_in_window",
    "mz",
    "mz_frac_of_precursor",
    "local_peak_density",
    "spacing_left",
    "spacing_right",
    "mass_defect",
    "has_isotope_plus1",
    "has_isotope_minus1",
    "has_complement",
    "has_h2o_loss_partner",
    "has_nh3_loss_partner",
];

pub const N_FEATURES: usize = FEATURE_NAMES.len();

pub const SPACING_SENTINEL: f32 = 1000.0;

pub const FEATURE_WINDOW_DA: f64 = 50.0;

use model::mass::{ISOTOPE, PROTON};
use model::tolerance::Tolerance;

use model::mass::{H2O, NH3};

pub struct PeakFeatureCtx {
    pub precursor_mz: f64,
    pub charge: u8,
    pub parent_neutral_mass: f64,
    pub total_intensity: f64,
    pub base_peak_intensity: f32,
    pub window_da: f64,
    pub match_tol_da: f64,
}

impl PeakFeatureCtx {
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

pub fn extract_peak_features(
    peaks: &[(f64, f32)],
    ranks: &[u32],
    ctx: &PeakFeatureCtx,
) -> Vec<[f32; N_FEATURES]> {
    debug_assert_eq!(
        peaks.len(),
        ranks.len(),
        "peaks and ranks must be 1:1 aligned"
    );
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

        let lo = peaks.partition_point(|&(m, _)| m < mz - ctx.window_da);
        let hi = peaks.partition_point(|&(m, _)| m <= mz + ctx.window_da);
        let win = &peaks[lo..hi];
        let win_count = win.len().max(1);
        let more_intense = win.iter().filter(|&&(_, pint)| pint > intensity).count();
        f[4] = more_intense as f32 / win_count as f32;
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
        f[15] = has_peak(peaks, ctx.parent_neutral_mass + 2.0 * PROTON - mz, tol);
        f[16] = has_peak(peaks, mz - H2O / z, tol);
        f[17] = has_peak(peaks, mz - NH3 / z, tol);

        out.push(f);
    }
    out
}

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
        assert_eq!(
            N_FEATURES, 18,
            "feature set changed — update the Python contract (feature_spec.py) and A2.10 fixtures"
        );
    }

    #[test]
    fn extracts_expected_features_on_tiny_scan() {
        let peaks = vec![(100.0_f64, 10.0_f32), (200.05, 100.0), (300.0, 50.0)];
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
        let p1 = &f[1];
        assert!((p1[idx("intensity_over_basepeak")] - 1.0).abs() < 1e-6);
        assert!((p1[idx("intensity_over_tic")] - (100.0 / 160.0)).abs() < 1e-6);
        assert!((p1[idx("global_rank_frac")] - 0.0).abs() < 1e-6);
        assert!((p1[idx("mass_defect")] - 0.05).abs() < 1e-4);
        assert!((p1[idx("mz_frac_of_precursor")] - (200.05 / 500.0)).abs() < 1e-6);
        assert!((p1[idx("is_top1_in_window")] - 1.0).abs() < 1e-6);

        assert!((f[0][idx("spacing_left")] - SPACING_SENTINEL).abs() < 1e-6);
        assert!((f[2][idx("spacing_right")] - SPACING_SENTINEL).abs() < 1e-6);
        assert!((f[0][idx("spacing_right")] - 100.05).abs() < 1e-3);
        assert!((f[1][idx("has_isotope_plus1")] - 0.0).abs() < 1e-6);
        assert!((f[1][idx("has_complement")] - 0.0).abs() < 1e-6);
    }

    #[test]
    fn filtered_peak_gets_max_rank_frac() {
        let peaks = vec![(150.0_f64, 5.0_f32), (250.0, 80.0)];
        let ranks = vec![u32::MAX, 1];
        let ctx = PeakFeatureCtx {
            precursor_mz: 400.0,
            charge: 2,
            parent_neutral_mass: (400.0 - PROTON) * 2.0,
            total_intensity: 80.0,
            base_peak_intensity: 80.0,
            window_da: 50.0,
            match_tol_da: 0.5,
        };
        let f = extract_peak_features(&peaks, &ranks, &ctx);
        let gi = FEATURE_NAMES
            .iter()
            .position(|&n| n == "global_rank_frac")
            .unwrap();
        assert!(
            (f[0][gi] - 1.0).abs() < 1e-6,
            "filtered peak must get rank_frac 1.0"
        );
    }

    #[test]
    fn for_spectrum_uses_active_list_for_both_denominators() {
        use model::tolerance::Tolerance;
        let peaks = vec![(100.0_f64, 10.0_f32), (200.0, 40.0), (300.0, 50.0)];
        let ctx = PeakFeatureCtx::for_spectrum(
            500.0,
            2,
            (500.0 - 1.00727649) * 2.0,
            &peaks,
            &Tolerance::Da(0.5),
        );
        assert!((ctx.base_peak_intensity - 50.0).abs() < 1e-6);
        assert!((ctx.total_intensity - 100.0).abs() < 1e-9);
        assert!((ctx.window_da - FEATURE_WINDOW_DA).abs() < 1e-9);
        assert!((ctx.match_tol_da - 0.5).abs() < 1e-9);
    }

    fn idx(name: &str) -> usize {
        FEATURE_NAMES
            .iter()
            .position(|&n| n == name)
            .expect("feature exists")
    }
}
