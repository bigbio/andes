//! Spectral-match features between predicted and observed fragment intensities.

use super::form::{Mass, H2O, PROTON};

/// One fragment ion type in a predictor's per-site output: (is b ion, charge).
pub type Ion = (bool, u8);

/// Names of the values [`ms2_features`] returns, in order.
pub const MS2_FEATURES: [&str; 8] = [
    "spec_pearson_norm",
    "ionb_pearson_norm",
    "iony_pearson_norm",
    "spec_mse_norm",
    "dotprod_norm",
    "cos_norm",
    "spectral_angle",
    "pred_explained",
];

/// Fragment m/z tolerance for matching predicted ions to observed peaks.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum FragmentTolerance {
    /// Absolute window in Da.
    Da(f64),
    /// Relative window in ppm of the fragment m/z.
    Ppm(f64),
}

impl FragmentTolerance {
    fn window(self, mz: f64) -> f64 {
        match self {
            FragmentTolerance::Da(d) => d,
            FragmentTolerance::Ppm(p) => mz * p * 1e-6,
        }
    }
}

fn pearson(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len() as f64;
    if a.len() < 2 {
        return 0.0;
    }
    let (ma, mb) = (a.iter().sum::<f64>() / n, b.iter().sum::<f64>() / n);
    let (mut sab, mut saa, mut sbb) = (0.0, 0.0, 0.0);
    for (x, y) in a.iter().zip(b) {
        sab += (x - ma) * (y - mb);
        saa += (x - ma) * (x - ma);
        sbb += (y - mb) * (y - mb);
    }
    if saa <= 0.0 || sbb <= 0.0 {
        0.0
    } else {
        sab / (saa * sbb).sqrt()
    }
}

/// Spectral-match features ([`MS2_FEATURES`]) between `pred` and the observed
/// `peaks` (sorted by m/z).
///
/// `pred` holds `ions.len()` linear intensities per cleavage site (site `i` =
/// cleavage after residue `i + 1`). Ions above the precursor `charge` are
/// skipped. The observed intensity of an ion is the most intense peak within
/// `tol` of its m/z (0 when none). Predicted and observed vectors are each
/// normalised to unit sum before the comparisons; `pred_explained` is the
/// predicted share of ions that matched a peak.
pub fn ms2_features(
    pred: &[f32],
    ions: &[Ion],
    mass: &Mass,
    charge: u8,
    peaks: &[(f64, f32)],
    tol: FragmentTolerance,
) -> [f64; 8] {
    let n = mass.res.len();
    let total: f64 = mass.res.iter().sum::<f64>() + mass.nterm;
    let mut p = Vec::new();
    let mut o = Vec::new();
    let mut is_b = Vec::new();
    let mut prefix = mass.nterm;
    for i in 0..n.saturating_sub(1) {
        prefix += mass.res[i];
        for (k, &(bion, z)) in ions.iter().enumerate() {
            if z > charge {
                continue;
            }
            let neutral = if bion { prefix } else { total - prefix + H2O };
            let mz = (neutral + z as f64 * PROTON) / z as f64;
            let w = tol.window(mz);
            let lo = peaks.partition_point(|(m, _)| *m < mz - w);
            let obs = peaks[lo..]
                .iter()
                .take_while(|(m, _)| *m <= mz + w)
                .map(|(_, i)| *i as f64)
                .fold(0.0, f64::max);
            p.push(pred.get(i * ions.len() + k).copied().unwrap_or(0.0) as f64);
            o.push(obs);
            is_b.push(bion);
        }
    }
    let (sp, so) = (p.iter().sum::<f64>(), o.iter().sum::<f64>());
    let pn: Vec<f64> = p
        .iter()
        .map(|v| if sp > 0.0 { v / sp } else { 0.0 })
        .collect();
    let on: Vec<f64> = o
        .iter()
        .map(|v| if so > 0.0 { v / so } else { 0.0 })
        .collect();
    let sel = |want_b: bool, v: &[f64]| -> Vec<f64> {
        v.iter()
            .zip(&is_b)
            .filter(|(_, b)| **b == want_b)
            .map(|(x, _)| *x)
            .collect()
    };
    let dot: f64 = pn.iter().zip(&on).map(|(a, b)| a * b).sum();
    let (np, no) = (
        pn.iter().map(|v| v * v).sum::<f64>().sqrt(),
        on.iter().map(|v| v * v).sum::<f64>().sqrt(),
    );
    let cos = if np > 0.0 && no > 0.0 {
        (dot / (np * no)).clamp(-1.0, 1.0)
    } else {
        0.0
    };
    let sa = if np > 0.0 && no > 0.0 {
        1.0 - 2.0 * cos.acos() / std::f64::consts::PI
    } else {
        0.0
    };
    let mse = pn
        .iter()
        .zip(&on)
        .map(|(a, b)| (a - b) * (a - b))
        .sum::<f64>()
        / pn.len().max(1) as f64;
    let explained: f64 = pn
        .iter()
        .zip(&o)
        .filter(|(_, ob)| **ob > 0.0)
        .map(|(v, _)| *v)
        .sum();
    [
        pearson(&pn, &on),
        pearson(&sel(true, &pn), &sel(true, &on)),
        pearson(&sel(false, &pn), &sel(false, &on)),
        mse,
        dot,
        cos,
        sa,
        explained,
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::predictors::form::parse_peptide;

    const BY1: [Ion; 2] = [(true, 1), (false, 1)];

    /// Peaks at every predicted b/y z1 m/z with the predicted intensities.
    fn ideal_peaks(pred: &[f32], mass: &Mass) -> Vec<(f64, f32)> {
        let total: f64 = mass.res.iter().sum::<f64>() + mass.nterm;
        let mut prefix = mass.nterm;
        let mut peaks = Vec::new();
        for i in 0..mass.res.len() - 1 {
            prefix += mass.res[i];
            peaks.push((prefix + PROTON, pred[2 * i]));
            peaks.push((total - prefix + H2O + PROTON, pred[2 * i + 1]));
        }
        peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
        peaks
    }

    #[test]
    fn a_perfect_match_scores_one() {
        let (_, mass) = parse_peptide("K.PEPTIDEK.R").unwrap();
        let pred: Vec<f32> = (0..14).map(|i| 1.0 + (i % 5) as f32).collect();
        let peaks = ideal_peaks(&pred, &mass);
        let f = ms2_features(&pred, &BY1, &mass, 2, &peaks, FragmentTolerance::Ppm(20.0));
        assert!((f[0] - 1.0).abs() < 1e-9, "pearson {}", f[0]);
        assert!((f[5] - 1.0).abs() < 1e-9, "cosine {}", f[5]);
        assert!((f[6] - 1.0).abs() < 1e-6, "spectral angle {}", f[6]);
        assert!((f[7] - 1.0).abs() < 1e-9, "explained {}", f[7]);
        assert!(f[3] < 1e-12, "mse {}", f[3]);
    }

    #[test]
    fn no_peaks_scores_zero() {
        let (_, mass) = parse_peptide("K.PEPTIDEK.R").unwrap();
        let pred = vec![1.0f32; 14];
        let f = ms2_features(&pred, &BY1, &mass, 2, &[], FragmentTolerance::Da(0.5));
        assert_eq!(f[5], 0.0);
        assert_eq!(f[6], 0.0);
        assert_eq!(f[7], 0.0);
    }

    #[test]
    fn ions_above_the_precursor_charge_are_skipped() {
        let (_, mass) = parse_peptide("K.PEPTIDEK.R").unwrap();
        let ions: [Ion; 2] = [(true, 1), (true, 3)];
        let pred = vec![1.0f32; 14];
        let peaks = ideal_peaks(&[1.0; 14], &mass);
        let f = ms2_features(&pred, &ions, &mass, 2, &peaks, FragmentTolerance::Ppm(20.0));
        // Only b z1 is compared, and every b z1 ion is present.
        assert!((f[7] - 1.0).abs() < 1e-9);
    }
}
