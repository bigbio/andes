//! Engine-wide retention-time (RT) prediction: engineered features, a
//! retention-coefficient RT-INDEX predictor, and per-run linear calibration.
//!
//! This is a general peptide/backbone capability (it benefits ordinary
//! peptide search, not just glyco). The glycan per-monosaccharide RT offset is
//! a *glyco extension* that lives in the `andes-glyco` crate and simply adds to
//! the backbone index produced here — this module stays glyco-free.
//!
//! Design: `docs/plans/glyco/50-roadmap/rt-prediction-design.md`.
//!
//! # Pipeline
//! 1. [`rt_features`] turns a [`Peptide`] into a fixed-length, deterministic
//!    feature vector (AA composition, length, mass, Kyte-Doolittle hydropathy
//!    summed/mean, position-weighted N-/C-terminal hydropathy, strong-retention
//!    vs polar/charged counts, mod count + summed mod mass-delta).
//! 2. [`RtIndexModel`] is a linear retention-coefficient model producing a
//!    dimensionless RT *index* (`Σ coef·feature`). It ships with a
//!    Kyte-Doolittle-derived seed ([`RtIndexModel::seed`]) and a
//!    least-squares [`RtIndexModel::fit`]. The linear model is intentionally
//!    model-agnostic in interface so a GBDT can replace it later.
//! 3. [`calibrate`] fits a robust per-run `rt_obs = slope·index + intercept`
//!    (the iRT mechanism) from confident-PSM anchors; [`predicted_rt`] applies
//!    it. This decouples transferable peptide chemistry from per-run LC
//!    geometry.
//!
//! All ordered outputs use `Vec`/slices and total-order comparisons — there is
//! no `HashMap` on any path that feeds a deterministic result.

use model::peptide::Peptide;

use crate::frag_features::hydropathy;

/// The 20 standard amino acids, in a fixed canonical order. This order defines
/// feature indices `0..20`, so it MUST never be reordered once coefficients are
/// trained against it.
pub const AA_ORDER: [u8; 20] = *b"ARNDCQEGHILKMFPSTWYV";

/// "Strong-retention" residues on reversed phase (bulky hydrophobic): the set
/// SSRCalc/ELUDE weight most heavily for late elution.
const STRONG_RETENTION: [u8; 6] = *b"WFLIYM";

/// Polar/charged residues (net-hydrophilic), summed as an early-elution
/// counterweight.
const POLAR_CHARGED: [u8; 8] = *b"DEKRHNQS";

// ---- Feature layout (fixed; indices are load-bearing for trained coeffs) ----

/// Number of AA composition-count features (one per standard residue).
pub const N_AA_FEATURES: usize = 20;

/// Feature index: peptide length.
pub const FEAT_LENGTH: usize = 20;
/// Feature index: peptide neutral mass (kDa, scaled so magnitude ~ O(1)).
pub const FEAT_MASS_KDA: usize = 21;
/// Feature index: summed Kyte-Doolittle hydropathy over all residues.
pub const FEAT_HYDRO_SUM: usize = 22;
/// Feature index: mean Kyte-Doolittle hydropathy (sum / length).
pub const FEAT_HYDRO_MEAN: usize = 23;
/// Feature index: position-weighted N-terminal residue hydropathy
/// (SSRCalc's single biggest lever — terminal residues retain differently).
pub const FEAT_HYDRO_NTERM: usize = 24;
/// Feature index: position-weighted C-terminal residue hydropathy.
pub const FEAT_HYDRO_CTERM: usize = 25;
/// Feature index: count of strong-retention residues (W,F,L,I,Y,M).
pub const FEAT_STRONG_COUNT: usize = 26;
/// Feature index: count of polar/charged residues (D,E,K,R,H,N,Q,S).
pub const FEAT_POLAR_COUNT: usize = 27;
/// Feature index: number of modified residues.
pub const FEAT_MOD_COUNT: usize = 28;
/// Feature index: summed modification mass-delta (Da).
pub const FEAT_MOD_MASS_SUM: usize = 29;
/// Feature index: the constant/bias term (always 1.0) — lets the linear model
/// carry an intercept without special-casing.
pub const FEAT_BIAS: usize = 30;

/// Total engineered RT feature-vector length.
pub const N_RT_FEATURES: usize = 31;

/// Weight applied to terminal-residue hydropathy. Terminal residues have an
/// outsized, well-documented effect on RP retention; a >1 weight lets the model
/// express that without a separate learned term for every residue identity.
/// This is a modelling prior, not a fitted constant — flagged for re-fit.
const TERMINAL_HYDRO_WEIGHT: f32 = 1.0;

/// Extract the deterministic engineered RT feature vector for `peptide`.
///
/// The vector length is [`N_RT_FEATURES`]. Index semantics are the `FEAT_*`
/// constants above; entries `0..20` are AA composition counts in [`AA_ORDER`].
pub fn rt_features(peptide: &Peptide) -> Vec<f32> {
    let mut f = vec![0.0f32; N_RT_FEATURES];
    let residues = &peptide.residues;
    let n = residues.len();

    // AA composition counts (indices 0..20).
    for aa in residues {
        if let Some(i) = AA_ORDER.iter().position(|&r| r == aa.residue) {
            f[i] += 1.0;
        }
    }

    // Length + mass.
    f[FEAT_LENGTH] = n as f32;
    f[FEAT_MASS_KDA] = (peptide.mass() / 1000.0) as f32;

    // Hydropathy sum + mean (reuse the single Kyte-Doolittle source).
    let hydro_sum: f32 = residues.iter().map(|aa| hydropathy(aa.residue)).sum();
    f[FEAT_HYDRO_SUM] = hydro_sum;
    f[FEAT_HYDRO_MEAN] = if n > 0 { hydro_sum / n as f32 } else { 0.0 };

    // Terminal residue hydropathy (position-weighted).
    if n > 0 {
        f[FEAT_HYDRO_NTERM] = TERMINAL_HYDRO_WEIGHT * hydropathy(residues[0].residue);
        f[FEAT_HYDRO_CTERM] = TERMINAL_HYDRO_WEIGHT * hydropathy(residues[n - 1].residue);
    }

    // Strong-retention vs polar/charged counts.
    f[FEAT_STRONG_COUNT] = residues
        .iter()
        .filter(|aa| STRONG_RETENTION.contains(&aa.residue))
        .count() as f32;
    f[FEAT_POLAR_COUNT] = residues
        .iter()
        .filter(|aa| POLAR_CHARGED.contains(&aa.residue))
        .count() as f32;

    // Modifications: count + summed mass-delta.
    let mut mod_count = 0.0f32;
    let mut mod_mass_sum = 0.0f32;
    for aa in residues {
        if let Some(m) = &aa.mod_ {
            mod_count += 1.0;
            mod_mass_sum += m.mass_delta as f32;
        }
    }
    f[FEAT_MOD_COUNT] = mod_count;
    f[FEAT_MOD_MASS_SUM] = mod_mass_sum;

    // Bias term.
    f[FEAT_BIAS] = 1.0;

    f
}

/// A linear retention-coefficient RT-INDEX model: `index = Σ coefᵢ·featureᵢ`.
///
/// The predicted value is a **dimensionless index**, not minutes — per-run
/// [`calibrate`] maps it to observed RT. The interface is deliberately kept
/// model-agnostic (feature-vec in, scalar index out) so a GBDT regressor can
/// drop in later without touching callers.
#[derive(Debug, Clone, PartialEq)]
pub struct RtIndexModel {
    /// One coefficient per feature; length is [`N_RT_FEATURES`].
    pub coef: Vec<f64>,
}

impl RtIndexModel {
    /// Build from an explicit coefficient vector. Returns `None` if the length
    /// is not [`N_RT_FEATURES`].
    pub fn new(coef: Vec<f64>) -> Option<Self> {
        if coef.len() == N_RT_FEATURES {
            Some(Self { coef })
        } else {
            None
        }
    }

    /// A sane default/seed model with per-residue coefficients set to the
    /// **Kyte-Doolittle hydropathy scale** (Kyte J. & Doolittle R.F.,
    /// *J. Mol. Biol.* 157(1):105–132, 1982). More hydrophobic residues get a
    /// larger positive coefficient, so more-hydrophobic peptides get a larger
    /// index — matching reversed-phase elution order. This is a *prior*: it
    /// gives a usable ranking with zero training data, but should be replaced
    /// by [`fit`](Self::fit) on andes data. Non-composition features seed at 0
    /// so the seed index reduces to summed hydropathy (equivalent to
    /// `FEAT_HYDRO_SUM` with unit weight, expressed per-residue here so `fit`
    /// can refine each residue independently).
    pub fn seed() -> Self {
        let mut coef = vec![0.0f64; N_RT_FEATURES];
        for (i, &r) in AA_ORDER.iter().enumerate() {
            coef[i] = hydropathy(r) as f64;
        }
        Self { coef }
    }

    /// Predict the dimensionless RT index for a feature vector.
    ///
    /// `features.len()` must equal `self.coef.len()`; otherwise the shorter
    /// length governs the dot product (defensive, never truncates a
    /// correctly-built [`rt_features`] vector).
    pub fn predict_index(&self, features: &[f32]) -> f64 {
        self.coef
            .iter()
            .zip(features.iter())
            .map(|(c, x)| c * *x as f64)
            .sum()
    }

    /// Convenience: extract features from a peptide and predict.
    pub fn predict_peptide(&self, peptide: &Peptide) -> f64 {
        self.predict_index(&rt_features(peptide))
    }

    /// Learn coefficients by ordinary least squares (normal equations) from
    /// `(features, observed_rt)` pairs. `observed_rt` is used directly as the
    /// regression target, so the learned index is already in the target's
    /// units when trained this way (calibration then remains identity); when
    /// trained against an arbitrary index scale it stays dimensionless.
    ///
    /// Returns `None` if there are fewer samples than features (under-
    /// determined) or if the normal-equations system is singular.
    pub fn fit(samples: &[(Vec<f32>, f64)]) -> Option<Self> {
        let p = N_RT_FEATURES;
        if samples.len() < p {
            return None;
        }
        // Normal equations: (XᵀX) β = Xᵀy.
        let mut xtx = vec![vec![0.0f64; p]; p];
        let mut xty = vec![0.0f64; p];
        for (feat, y) in samples {
            if feat.len() != p {
                return None;
            }
            for i in 0..p {
                let fi = feat[i] as f64;
                xty[i] += fi * y;
                for j in 0..p {
                    xtx[i][j] += fi * feat[j] as f64;
                }
            }
        }
        let coef = solve_linear_system(xtx, xty)?;
        Some(Self { coef })
    }
}

/// Minimum anchors required for a per-run calibration fit. Below this we cannot
/// robustly estimate a slope/intercept (and the outlier trim needs headroom).
pub const MIN_CALIBRATION_ANCHORS: usize = 5;

/// Per-run linear calibration `rt_obs = slope·index + intercept` (the iRT
/// mechanism). Robust to outliers: fits once, then trims the worst-residual
/// fraction and refits (a one-pass trimmed least squares).
///
/// `anchors` are `(index_pred, rt_obs)`; the fit is unit-agnostic (whatever unit
/// the `rt_obs` values carry is what the returned line predicts — the callers use
/// minutes). Non-finite anchors (NaN/inf in either coordinate) are dropped so a
/// single bad float cannot poison the fit into a NaN slope (Codex review).
/// Returns `None` if fewer than [`MIN_CALIBRATION_ANCHORS`] finite anchors remain
/// or the indices are degenerate (all equal ⇒ no slope is identifiable).
pub fn calibrate(anchors: &[(f64, f64)]) -> Option<(f64, f64)> {
    // Drop non-finite anchors up front (an inf/NaN would make the OLS means and
    // sums non-finite, and the zero-variance guard does not reject NaN).
    let finite: Vec<(f64, f64)> = anchors
        .iter()
        .copied()
        .filter(|&(x, y)| x.is_finite() && y.is_finite())
        .collect();
    if finite.len() < MIN_CALIBRATION_ANCHORS {
        return None;
    }
    let anchors = &finite[..];
    // First pass over all anchors.
    let (slope0, intercept0) = ols_line(anchors)?;

    // Compute residuals, trim the worst ~20%, refit. Deterministic: sort by
    // total-ordered residual magnitude.
    let mut resid: Vec<(f64, usize)> = anchors
        .iter()
        .enumerate()
        .map(|(i, &(x, y))| ((y - (slope0 * x + intercept0)).abs(), i))
        .collect();
    resid.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));

    let keep = ((anchors.len() as f64) * 0.8).ceil() as usize;
    let keep = keep.max(MIN_CALIBRATION_ANCHORS).min(anchors.len());
    let trimmed: Vec<(f64, f64)> = resid[..keep].iter().map(|&(_, i)| anchors[i]).collect();

    ols_line(&trimmed).or(Some((slope0, intercept0)))
}

/// Apply a calibration fit **in seconds** to a predicted index, returning RT in
/// minutes (divides by 60).
///
/// ⚠️ UNIT CONTRACT: this is only correct when [`calibrate`] was fit on
/// seconds-valued anchors. The production wiring (`output::rt_wiring` /
/// `output::glyco_rt`) fits `calibrate` on **minutes** anchors and applies the
/// line directly (`slope·index + intercept`), so it does NOT use this helper —
/// calling it on a minutes-fit line would divide by 60 twice. Kept only for a
/// seconds-domain caller; do not wire it into the minutes-based path.
pub fn predicted_rt(index: f64, slope: f64, intercept: f64) -> f64 {
    (slope * index + intercept) / 60.0
}

/// Ordinary least-squares fit of a single line `y = slope·x + intercept`.
/// Returns `None` if `x` has zero variance (slope not identifiable).
fn ols_line(points: &[(f64, f64)]) -> Option<(f64, f64)> {
    let n = points.len() as f64;
    if n < 2.0 {
        return None;
    }
    let mean_x = points.iter().map(|p| p.0).sum::<f64>() / n;
    let mean_y = points.iter().map(|p| p.1).sum::<f64>() / n;
    let mut sxx = 0.0;
    let mut sxy = 0.0;
    for &(x, y) in points {
        let dx = x - mean_x;
        sxx += dx * dx;
        sxy += dx * (y - mean_y);
    }
    if !sxx.is_finite() || sxx.abs() < 1e-12 {
        return None;
    }
    let slope = sxy / sxx;
    let intercept = mean_y - slope * mean_x;
    if !slope.is_finite() || !intercept.is_finite() {
        return None;
    }
    Some((slope, intercept))
}

/// Solve `A x = b` for a square symmetric-positive-(semi)definite `A` via
/// Gaussian elimination with partial pivoting. Returns `None` if singular.
/// Deterministic (no floating-point-keyed containers).
fn solve_linear_system(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let n = b.len();
    for col in 0..n {
        // Partial pivot: pick the row with the largest magnitude in this col.
        let mut pivot = col;
        let mut best = a[col][col].abs();
        for row in (col + 1)..n {
            let v = a[row][col].abs();
            if v > best {
                best = v;
                pivot = row;
            }
        }
        if best < 1e-12 {
            return None; // singular
        }
        a.swap(col, pivot);
        b.swap(col, pivot);

        // Eliminate below.
        for row in (col + 1)..n {
            let factor = a[row][col] / a[col][col];
            if factor != 0.0 {
                for k in col..n {
                    a[row][k] -= factor * a[col][k];
                }
                b[row] -= factor * b[col];
            }
        }
    }
    // Back-substitution.
    let mut x = vec![0.0f64; n];
    for row in (0..n).rev() {
        let mut s = b[row];
        for k in (row + 1)..n {
            s -= a[row][k] * x[k];
        }
        x[row] = s / a[row][row];
    }
    Some(x)
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::amino_acid::AminoAcid;
    use model::modification::{ModLocation, Modification, ResidueSpec};

    fn unmod_pep(seq: &[u8]) -> Peptide {
        let residues: Vec<_> = seq.iter().map(|&r| AminoAcid::standard(r).unwrap()).collect();
        Peptide::new(residues, b'_', b'-')
    }

    #[test]
    fn features_are_deterministic_and_correct() {
        // PEPTIDE = P,E,P,T,I,D,E
        let p = unmod_pep(b"PEPTIDE");
        let f = rt_features(&p);
        assert_eq!(f.len(), N_RT_FEATURES);

        // Composition counts. AA_ORDER = ARNDCQEGHILKMFPSTWYV
        let idx = |r: u8| AA_ORDER.iter().position(|&x| x == r).unwrap();
        assert_eq!(f[idx(b'P')], 2.0);
        assert_eq!(f[idx(b'E')], 2.0);
        assert_eq!(f[idx(b'T')], 1.0);
        assert_eq!(f[idx(b'I')], 1.0);
        assert_eq!(f[idx(b'D')], 1.0);
        assert_eq!(f[idx(b'A')], 0.0);

        assert_eq!(f[FEAT_LENGTH], 7.0);
        assert_eq!(f[FEAT_BIAS], 1.0);

        // Hydropathy sum: P(-1.6)+E(-3.5)+P(-1.6)+T(-0.7)+I(4.5)+D(-3.5)+E(-3.5)
        let expected_sum = -1.6 - 3.5 - 1.6 - 0.7 + 4.5 - 3.5 - 3.5;
        assert!((f[FEAT_HYDRO_SUM] - expected_sum).abs() < 1e-4);
        assert!((f[FEAT_HYDRO_MEAN] - expected_sum / 7.0).abs() < 1e-4);

        // Terminal hydropathy: N-term P = -1.6, C-term E = -3.5.
        assert!((f[FEAT_HYDRO_NTERM] - (-1.6)).abs() < 1e-4);
        assert!((f[FEAT_HYDRO_CTERM] - (-3.5)).abs() < 1e-4);

        // Strong-retention (W,F,L,I,Y,M): only I → 1.
        assert_eq!(f[FEAT_STRONG_COUNT], 1.0);
        // Polar/charged (D,E,K,R,H,N,Q,S): E,D,E → 3.
        assert_eq!(f[FEAT_POLAR_COUNT], 3.0);

        // No mods.
        assert_eq!(f[FEAT_MOD_COUNT], 0.0);
        assert_eq!(f[FEAT_MOD_MASS_SUM], 0.0);

        // Determinism: same input, same output.
        assert_eq!(rt_features(&p), f);
    }

    #[test]
    fn features_capture_modifications() {
        let ox = Modification {
            name: "Oxidation".to_string(),
            mass_delta: 15.99491,
            residue: ResidueSpec::Specific(b'M'),
            location: ModLocation::Anywhere,
            fixed: false,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
        };
        let residues = vec![
            AminoAcid::standard(b'M').unwrap().with_mod(ox),
            AminoAcid::standard(b'K').unwrap(),
        ];
        let p = Peptide::new(residues, b'_', b'-');
        let f = rt_features(&p);
        assert_eq!(f[FEAT_MOD_COUNT], 1.0);
        assert!((f[FEAT_MOD_MASS_SUM] - 15.99491).abs() < 1e-4);
    }

    #[test]
    fn seed_index_ranks_hydrophobic_higher() {
        let model = RtIndexModel::seed();
        // A hydrophobic peptide should out-index a hydrophilic one.
        let hydrophobic = model.predict_peptide(&unmod_pep(b"LLIIVV"));
        let hydrophilic = model.predict_peptide(&unmod_pep(b"DEKRDE"));
        assert!(
            hydrophobic > hydrophilic,
            "hydrophobic {hydrophobic} should exceed hydrophilic {hydrophilic}"
        );
    }

    #[test]
    fn fit_recovers_known_coefficients() {
        // Ground-truth coefficient vector (arbitrary but fixed).
        let mut truth = vec![0.0f64; N_RT_FEATURES];
        for (i, c) in truth.iter_mut().enumerate() {
            *c = 0.1 * (i as f64 + 1.0);
        }
        let truth_model = RtIndexModel::new(truth.clone()).unwrap();

        // Generate deterministic pseudo-random feature vectors and exact targets.
        let mut samples = Vec::new();
        let mut state: u64 = 0x1234_5678;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 40) as f32 / 16_777_216.0 // in [0,1)
        };
        for _ in 0..200 {
            let mut feat = vec![0.0f32; N_RT_FEATURES];
            for v in feat.iter_mut() {
                *v = next() * 3.0;
            }
            feat[FEAT_BIAS] = 1.0;
            let y = truth_model.predict_index(&feat);
            samples.push((feat, y));
        }

        let fitted = RtIndexModel::fit(&samples).expect("fit should succeed");
        for (i, (a, b)) in fitted.coef.iter().zip(truth.iter()).enumerate() {
            assert!((a - b).abs() < 1e-6, "coef[{i}] fitted {a} != truth {b}");
        }
    }

    #[test]
    fn fit_returns_none_when_underdetermined() {
        let samples: Vec<(Vec<f32>, f64)> = (0..3)
            .map(|_| (vec![1.0f32; N_RT_FEATURES], 1.0))
            .collect();
        assert!(RtIndexModel::fit(&samples).is_none());
    }

    #[test]
    fn calibrate_recovers_known_slope_intercept() {
        // rt = 2.5 * index + 100 (seconds).
        let anchors: Vec<(f64, f64)> = (0..50)
            .map(|i| {
                let index = i as f64 * 0.5;
                (index, 2.5 * index + 100.0)
            })
            .collect();
        let (slope, intercept) = calibrate(&anchors).expect("calibrate should succeed");
        assert!((slope - 2.5).abs() < 1e-6, "slope {slope}");
        assert!((intercept - 100.0).abs() < 1e-6, "intercept {intercept}");
    }

    #[test]
    fn calibrate_is_robust_to_outliers() {
        let mut anchors: Vec<(f64, f64)> = (0..50)
            .map(|i| {
                let index = i as f64 * 0.5;
                (index, 2.5 * index + 100.0)
            })
            .collect();
        // Inject a few gross outliers; trimmed fit should still recover the line.
        anchors[3].1 += 5000.0;
        anchors[10].1 -= 8000.0;
        anchors[40].1 += 12000.0;
        let (slope, intercept) = calibrate(&anchors).expect("calibrate should succeed");
        assert!((slope - 2.5).abs() < 0.2, "slope {slope} not robust");
        assert!((intercept - 100.0).abs() < 20.0, "intercept {intercept} not robust");
    }

    #[test]
    fn calibrate_too_few_anchors_is_none() {
        let anchors: Vec<(f64, f64)> = (0..MIN_CALIBRATION_ANCHORS - 1)
            .map(|i| (i as f64, i as f64))
            .collect();
        assert!(calibrate(&anchors).is_none());
    }

    #[test]
    fn calibrate_degenerate_indices_is_none() {
        // All indices identical ⇒ slope not identifiable.
        let anchors: Vec<(f64, f64)> = (0..10).map(|i| (3.0, i as f64)).collect();
        assert!(calibrate(&anchors).is_none());
    }

    #[test]
    fn calibrate_drops_non_finite_anchors_and_stays_finite() {
        // A clean line y = 2x + 1 with NaN/inf anchors injected. The bad anchors
        // must be dropped (not poison the fit into a NaN slope), leaving the
        // recovered line finite and correct.
        let mut anchors: Vec<(f64, f64)> = (0..8).map(|i| (i as f64, 2.0 * i as f64 + 1.0)).collect();
        anchors.push((f64::NAN, 5.0));
        anchors.push((3.0, f64::INFINITY));
        anchors.push((f64::NEG_INFINITY, f64::NAN));
        let (slope, intercept) = calibrate(&anchors).expect("finite anchors remain");
        assert!(slope.is_finite() && intercept.is_finite());
        assert!((slope - 2.0).abs() < 1e-9, "slope {slope}");
        assert!((intercept - 1.0).abs() < 1e-9, "intercept {intercept}");
    }

    #[test]
    fn calibrate_all_non_finite_is_none() {
        // Fewer than MIN_CALIBRATION_ANCHORS finite anchors ⇒ None (neutral path).
        let anchors: Vec<(f64, f64)> = (0..10).map(|_| (f64::NAN, f64::INFINITY)).collect();
        assert!(calibrate(&anchors).is_none());
    }

    #[test]
    fn predicted_rt_converts_seconds_to_minutes() {
        // slope/intercept in seconds; index 10 → 600 s → 10 min.
        let rt = predicted_rt(10.0, 50.0, 100.0);
        assert!((rt - 10.0).abs() < 1e-9, "rt {rt}");
    }
}
