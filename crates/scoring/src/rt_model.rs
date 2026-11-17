use model::peptide::Peptide;

use crate::frag_features::hydropathy;

pub const AA_ORDER: [u8; 20] = *b"ARNDCQEGHILKMFPSTWYV";

const STRONG_RETENTION: [u8; 6] = *b"WFLIYM";

const POLAR_CHARGED: [u8; 8] = *b"DEKRHNQS";

pub const N_AA_FEATURES: usize = 20;

pub const FEAT_LENGTH: usize = 20;
pub const FEAT_MASS_KDA: usize = 21;
pub const FEAT_HYDRO_SUM: usize = 22;
pub const FEAT_HYDRO_MEAN: usize = 23;
pub const FEAT_HYDRO_NTERM: usize = 24;
pub const FEAT_HYDRO_CTERM: usize = 25;
pub const FEAT_STRONG_COUNT: usize = 26;
pub const FEAT_POLAR_COUNT: usize = 27;
pub const FEAT_MOD_COUNT: usize = 28;
pub const FEAT_MOD_MASS_SUM: usize = 29;
pub const FEAT_BIAS: usize = 30;

pub const N_RT_FEATURES: usize = 31;

const TERMINAL_HYDRO_WEIGHT: f32 = 1.0;

pub fn rt_features(peptide: &Peptide) -> Vec<f32> {
    let mut f = vec![0.0f32; N_RT_FEATURES];
    let residues = &peptide.residues;
    let n = residues.len();

    for aa in residues {
        if let Some(i) = AA_ORDER.iter().position(|&r| r == aa.residue) {
            f[i] += 1.0;
        }
    }

    f[FEAT_LENGTH] = n as f32;
    f[FEAT_MASS_KDA] = (peptide.mass() / 1000.0) as f32;

    let hydro_sum: f32 = residues.iter().map(|aa| hydropathy(aa.residue)).sum();
    f[FEAT_HYDRO_SUM] = hydro_sum;
    f[FEAT_HYDRO_MEAN] = if n > 0 { hydro_sum / n as f32 } else { 0.0 };

    if n > 0 {
        f[FEAT_HYDRO_NTERM] = TERMINAL_HYDRO_WEIGHT * hydropathy(residues[0].residue);
        f[FEAT_HYDRO_CTERM] = TERMINAL_HYDRO_WEIGHT * hydropathy(residues[n - 1].residue);
    }

    f[FEAT_STRONG_COUNT] = residues
        .iter()
        .filter(|aa| STRONG_RETENTION.contains(&aa.residue))
        .count() as f32;
    f[FEAT_POLAR_COUNT] = residues
        .iter()
        .filter(|aa| POLAR_CHARGED.contains(&aa.residue))
        .count() as f32;

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

    f[FEAT_BIAS] = 1.0;

    f
}

#[derive(Debug, Clone, PartialEq)]
pub struct RtIndexModel {
    pub coef: Vec<f64>,
}

impl RtIndexModel {
    pub fn new(coef: Vec<f64>) -> Option<Self> {
        if coef.len() == N_RT_FEATURES {
            Some(Self { coef })
        } else {
            None
        }
    }

    pub fn seed() -> Self {
        let mut coef = vec![0.0f64; N_RT_FEATURES];
        for (i, &r) in AA_ORDER.iter().enumerate() {
            coef[i] = hydropathy(r) as f64;
        }
        Self { coef }
    }

    pub fn predict_index(&self, features: &[f32]) -> f64 {
        self.coef
            .iter()
            .zip(features.iter())
            .map(|(c, x)| c * *x as f64)
            .sum()
    }

    pub fn predict_peptide(&self, peptide: &Peptide) -> f64 {
        self.predict_index(&rt_features(peptide))
    }

    pub fn fit(samples: &[(Vec<f32>, f64)]) -> Option<Self> {
        let p = N_RT_FEATURES;
        if samples.len() < p {
            return None;
        }
        let mut xtx = vec![vec![0.0f64; p]; p];
        let mut xty = vec![0.0f64; p];
        for (feat, y) in samples {
            if feat.len() != p {
                return None;
            }
            for i in 0..p {
                let fi = feat[i] as f64;
                xty[i] += fi * y;
                for (j, &fj) in feat.iter().enumerate().take(p) {
                    xtx[i][j] += fi * fj as f64;
                }
            }
        }
        let coef = solve_linear_system(xtx, xty)?;
        Some(Self { coef })
    }
}

pub const MIN_CALIBRATION_ANCHORS: usize = 5;

pub fn calibrate(anchors: &[(f64, f64)]) -> Option<(f64, f64)> {
    let finite: Vec<(f64, f64)> = anchors
        .iter()
        .copied()
        .filter(|&(x, y)| x.is_finite() && y.is_finite())
        .collect();
    if finite.len() < MIN_CALIBRATION_ANCHORS {
        return None;
    }
    let anchors = &finite[..];
    let (slope0, intercept0) = ols_line(anchors)?;

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

pub fn predicted_rt(index: f64, slope: f64, intercept: f64) -> f64 {
    (slope * index + intercept) / 60.0
}

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

fn solve_linear_system(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let n = b.len();
    for col in 0..n {
        let mut pivot = col;
        let mut best = a[col][col].abs();
        for (row, a_row) in a.iter().enumerate().take(n).skip(col + 1) {
            let v = a_row[col].abs();
            if v > best {
                best = v;
                pivot = row;
            }
        }
        if best < 1e-12 {
            return None;
        }
        a.swap(col, pivot);
        b.swap(col, pivot);

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
    let mut x = vec![0.0f64; n];
    for row in (0..n).rev() {
        let mut s = b[row];
        for (k, &xk) in x.iter().enumerate().take(n).skip(row + 1) {
            s -= a[row][k] * xk;
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
        let residues: Vec<_> = seq
            .iter()
            .map(|&r| AminoAcid::standard(r).unwrap())
            .collect();
        Peptide::new(residues, b'_', b'-')
    }

    #[test]
    fn features_are_deterministic_and_correct() {
        let p = unmod_pep(b"PEPTIDE");
        let f = rt_features(&p);
        assert_eq!(f.len(), N_RT_FEATURES);

        let idx = |r: u8| AA_ORDER.iter().position(|&x| x == r).unwrap();
        assert_eq!(f[idx(b'P')], 2.0);
        assert_eq!(f[idx(b'E')], 2.0);
        assert_eq!(f[idx(b'T')], 1.0);
        assert_eq!(f[idx(b'I')], 1.0);
        assert_eq!(f[idx(b'D')], 1.0);
        assert_eq!(f[idx(b'A')], 0.0);

        assert_eq!(f[FEAT_LENGTH], 7.0);
        assert_eq!(f[FEAT_BIAS], 1.0);

        let expected_sum = -1.6 - 3.5 - 1.6 - 0.7 + 4.5 - 3.5 - 3.5;
        assert!((f[FEAT_HYDRO_SUM] - expected_sum).abs() < 1e-4);
        assert!((f[FEAT_HYDRO_MEAN] - expected_sum / 7.0).abs() < 1e-4);

        assert!((f[FEAT_HYDRO_NTERM] - (-1.6)).abs() < 1e-4);
        assert!((f[FEAT_HYDRO_CTERM] - (-3.5)).abs() < 1e-4);

        assert_eq!(f[FEAT_STRONG_COUNT], 1.0);
        assert_eq!(f[FEAT_POLAR_COUNT], 3.0);

        assert_eq!(f[FEAT_MOD_COUNT], 0.0);
        assert_eq!(f[FEAT_MOD_MASS_SUM], 0.0);

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
            annotated_site_mass: None,
            blocks_cleavage: false,
            parts: Vec::new(),
            diagnostic_ions: Vec::new(),
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
        let hydrophobic = model.predict_peptide(&unmod_pep(b"LLIIVV"));
        let hydrophilic = model.predict_peptide(&unmod_pep(b"DEKRDE"));
        assert!(
            hydrophobic > hydrophilic,
            "hydrophobic {hydrophobic} should exceed hydrophilic {hydrophilic}"
        );
    }

    #[test]
    fn fit_recovers_known_coefficients() {
        let mut truth = vec![0.0f64; N_RT_FEATURES];
        for (i, c) in truth.iter_mut().enumerate() {
            *c = 0.1 * (i as f64 + 1.0);
        }
        let truth_model = RtIndexModel::new(truth.clone()).unwrap();

        let mut samples = Vec::new();
        let mut state: u64 = 0x1234_5678;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 40) as f32 / 16_777_216.0
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
        let samples: Vec<(Vec<f32>, f64)> =
            (0..3).map(|_| (vec![1.0f32; N_RT_FEATURES], 1.0)).collect();
        assert!(RtIndexModel::fit(&samples).is_none());
    }

    #[test]
    fn calibrate_recovers_known_slope_intercept() {
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
        anchors[3].1 += 5000.0;
        anchors[10].1 -= 8000.0;
        anchors[40].1 += 12000.0;
        let (slope, intercept) = calibrate(&anchors).expect("calibrate should succeed");
        assert!((slope - 2.5).abs() < 0.2, "slope {slope} not robust");
        assert!(
            (intercept - 100.0).abs() < 20.0,
            "intercept {intercept} not robust"
        );
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
        let anchors: Vec<(f64, f64)> = (0..10).map(|i| (3.0, i as f64)).collect();
        assert!(calibrate(&anchors).is_none());
    }

    #[test]
    fn calibrate_drops_non_finite_anchors_and_stays_finite() {
        let mut anchors: Vec<(f64, f64)> =
            (0..8).map(|i| (i as f64, 2.0 * i as f64 + 1.0)).collect();
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
        let anchors: Vec<(f64, f64)> = (0..10).map(|_| (f64::NAN, f64::INFINITY)).collect();
        assert!(calibrate(&anchors).is_none());
    }

    #[test]
    fn predicted_rt_converts_seconds_to_minutes() {
        let rt = predicted_rt(10.0, 50.0, 100.0);
        assert!((rt - 10.0).abs() < 1e-9, "rt {rt}");
    }
}
