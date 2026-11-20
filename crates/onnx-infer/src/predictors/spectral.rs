use model::mass::{CH4SO, CO, H3PO4, NH3, Z_DOT_OFFSET};

use super::form::{Mass, H2O, PROTON};
use crate::ms2_native::{IonType, NeutralLoss, Series};

/// One predicted channel: ion series, fragment charge and neutral loss.
pub type Ion = IonType;

/// The m/z of a fragment of `series` at `frag_charge`.
///
/// `prefix_or_suffix_mass` is the neutral b-ion mass (N-terminal residues) for the a, b and
/// c series, and the neutral y-ion mass (C-terminal residues plus water) for the y and z•
/// series, after any neutral loss. a = b − CO, c = b + NH3, z• = y − NH3 + H.
pub fn channel_mz(series: Series, frag_charge: u8, prefix_or_suffix_mass: f64) -> f64 {
    let neutral = match series {
        Series::B | Series::Y => prefix_or_suffix_mass,
        Series::A => prefix_or_suffix_mass - CO,
        Series::C => prefix_or_suffix_mass + NH3,
        Series::Z => prefix_or_suffix_mass + Z_DOT_OFFSET,
    };
    (neutral + frag_charge as f64 * PROTON) / frag_charge as f64
}

/// The m/z of channel `ion` at the cleavage whose N-terminal residues weigh `prefix`, for a
/// peptide whose residues (with the N-terminal group) weigh `total`.
pub fn ion_mz(ion: &Ion, prefix: f64, total: f64) -> f64 {
    let base = if ion.is_prefix() {
        prefix
    } else {
        total - prefix + H2O
    };
    let neutral = match ion.loss {
        NeutralLoss::None => base,
        NeutralLoss::H3po4 => base - H3PO4,
        NeutralLoss::Ch4so => base - CH4SO,
    };
    channel_mz(ion.series, ion.charge, neutral)
}

/// One theoretical fragment of a predicted spectrum.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Fragment {
    /// Cleavage site, 0 for the bond after the first residue.
    pub site: usize,
    /// Index of the channel in the predictor's ions.
    pub channel: usize,
    pub mz: f64,
    /// Predicted intensity as the predictor returned it.
    pub pred: f64,
    pub is_prefix: bool,
}

/// The fragments of `pred` (site-major, `ions.len()` channels per site), skipping channels
/// whose fragment charge exceeds `charge` and, when `valid` is given, channels it marks
/// invalid (same layout as `pred`). A missing prediction counts as 0.
pub fn fragments(
    pred: &[f32],
    ions: &[Ion],
    valid: Option<&[bool]>,
    mass: &Mass,
    charge: u8,
) -> Vec<Fragment> {
    let n = mass.res.len();
    let total: f64 = mass.res.iter().sum::<f64>() + mass.nterm;
    let mut out = Vec::with_capacity(n.saturating_sub(1) * ions.len());
    let mut prefix = mass.nterm;
    for i in 0..n.saturating_sub(1) {
        prefix += mass.res[i];
        for (k, ion) in ions.iter().enumerate() {
            if ion.charge > charge {
                continue;
            }
            let at = i * ions.len() + k;
            if valid.is_some_and(|v| !v.get(at).copied().unwrap_or(false)) {
                continue;
            }
            out.push(Fragment {
                site: i,
                channel: k,
                mz: ion_mz(ion, prefix, total),
                pred: pred.get(at).copied().unwrap_or(0.0) as f64,
                is_prefix: ion.is_prefix(),
            });
        }
    }
    out
}

/// Indices of the `k` largest `weights`, largest first; ties keep their order.
pub fn top_by_weight(weights: &[f64], k: usize) -> Vec<usize> {
    let mut top: Vec<usize> = (0..weights.len()).collect();
    top.sort_by(|&x, &y| weights[y].total_cmp(&weights[x]));
    top.truncate(k);
    top
}

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

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum FragmentTolerance {
    Da(f64),
    Ppm(f64),
}

impl FragmentTolerance {
    pub fn window(self, mz: f64) -> f64 {
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

/// Similarity features of the predicted spectrum `pred` against `peaks`; see [`fragments`]
/// for `ions` and `valid`.
pub fn ms2_features(
    pred: &[f32],
    ions: &[Ion],
    valid: Option<&[bool]>,
    mass: &Mass,
    charge: u8,
    peaks: &[(f64, f32)],
    tol: FragmentTolerance,
) -> [f64; 8] {
    let frags = fragments(pred, ions, valid, mass, charge);
    let mut p = Vec::with_capacity(frags.len());
    let mut o = Vec::with_capacity(frags.len());
    let mut is_b = Vec::with_capacity(frags.len());
    for f in &frags {
        let mz = f.mz;
        let w = tol.window(mz);
        let lo = peaks.partition_point(|(m, _)| *m < mz - w);
        let obs = peaks[lo..]
            .iter()
            .take_while(|(m, _)| *m <= mz + w)
            .map(|(_, i)| *i as f64)
            .fold(0.0, f64::max);
        p.push(f.pred);
        o.push(obs);
        is_b.push(f.is_prefix);
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
    use crate::ms2_native::{activation_index, channel_mask, ION_TYPES, ION_TYPES_V2, N_ION_V2};
    use crate::predictors::form::parse_peptide;
    use model::mass::{C, H, N, O};

    const BY1: [Ion; 2] = [ION_TYPES[0], ION_TYPES[2]];

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
        let f = ms2_features(&pred, &BY1, None, &mass, 2, &peaks, FragmentTolerance::Ppm(20.0));
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
        let f = ms2_features(&pred, &BY1, None, &mass, 2, &[], FragmentTolerance::Da(0.5));
        assert_eq!(f[5], 0.0);
        assert_eq!(f[6], 0.0);
        assert_eq!(f[7], 0.0);
    }

    #[test]
    fn ions_above_the_precursor_charge_are_skipped() {
        let (_, mass) = parse_peptide("K.PEPTIDEK.R").unwrap();
        let ions: [Ion; 2] = [ION_TYPES[0], ION_TYPES_V2[12]];
        let pred = vec![1.0f32; 14];
        let peaks = ideal_peaks(&[1.0; 14], &mass);
        let f = ms2_features(&pred, &ions, None, &mass, 2, &peaks, FragmentTolerance::Ppm(20.0));
        assert!((f[7] - 1.0).abs() < 1e-9);
    }

    #[test]
    fn b_and_y_mz_keep_the_v1_arithmetic() {
        let (_, mass) = parse_peptide("K.PEPTIDEK.R").unwrap();
        let total: f64 = mass.res.iter().sum::<f64>() + mass.nterm;
        let frags = fragments(&[], &ION_TYPES[..4], None, &mass, 2);
        let mut prefix = mass.nterm;
        let mut want = Vec::new();
        for i in 0..mass.res.len() - 1 {
            prefix += mass.res[i];
            for (bion, z) in [(true, 1u8), (true, 2), (false, 1), (false, 2)] {
                let neutral = if bion { prefix } else { total - prefix + H2O };
                want.push((neutral + z as f64 * PROTON) / z as f64);
            }
        }
        let got: Vec<f64> = frags.iter().map(|f| f.mz).collect();
        assert_eq!(
            got.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
            want.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
        );
    }

    #[test]
    fn a_c_and_z_dot_mz_match_model_mass() {
        let (_, mass) = parse_peptide("K.PEPTIDE.R").unwrap();
        let pro = C * 5.0 + H * 7.0 + N + O;
        let glu = C * 5.0 + H * 7.0 + N + O * 3.0;
        let thr = C * 4.0 + H * 7.0 + N + O * 2.0;
        let ile = C * 6.0 + H * 11.0 + N + O;
        let asp = C * 4.0 + H * 5.0 + N + O * 3.0;
        let b2 = pro + glu;
        let b3 = b2 + pro;
        let b4 = b3 + thr;
        let y1 = glu + H2O;
        let y3 = ile + asp + glu + H2O;
        let close = |got: f64, want: f64| (got - want).abs() < 1e-6;

        assert!(close(channel_mz(Series::B, 1, b2), b2 + PROTON));
        assert!(close(channel_mz(Series::A, 1, b2), b2 - (C + O) + PROTON));
        assert!(close(channel_mz(Series::C, 1, b2), b2 + N + H * 3.0 + PROTON));
        assert!(close(channel_mz(Series::C, 2, b4), (b4 + N + H * 3.0 + 2.0 * PROTON) / 2.0));
        assert!(close(channel_mz(Series::Y, 1, y1), y1 + PROTON));
        assert!(close(channel_mz(Series::Z, 1, y1), y1 - (N + H * 2.0) + PROTON));
        assert!(close(channel_mz(Series::Z, 2, y3), (y3 - (N + H * 2.0) + 2.0 * PROTON) / 2.0));
        assert!(close(channel_mz(Series::B, 3, b3), (b3 + 3.0 * PROTON) / 3.0));

        assert!((channel_mz(Series::B, 1, b2) - 227.102633).abs() < 1e-5);
        assert!((channel_mz(Series::A, 1, b2) - 199.107719).abs() < 1e-5);
        assert!((channel_mz(Series::C, 1, b2) - 244.129183).abs() < 1e-5);
        assert!((channel_mz(Series::Y, 1, y1) - 148.060434).abs() < 1e-5);
        assert!((channel_mz(Series::Z, 1, y1) - 132.041710).abs() < 1e-5);

        let frags = fragments(&[], &ION_TYPES_V2, None, &mass, 3);
        let at = |site: usize, k: usize| {
            frags
                .iter()
                .find(|f| f.site == site && f.channel == k)
                .unwrap()
                .mz
        };
        let n = mass.res.len();
        assert!(close(at(1, 0), channel_mz(Series::B, 1, b2)));
        assert!(close(at(1, 14), channel_mz(Series::A, 1, b2)));
        assert!(close(at(1, 15), channel_mz(Series::C, 1, b2)));
        assert!(close(at(3, 16), channel_mz(Series::C, 2, b4)));
        assert!(close(at(2, 12), channel_mz(Series::B, 3, b3)));
        assert!(close(at(n - 2, 17), channel_mz(Series::Z, 1, y1)));
        assert!(close(at(n - 4, 18), channel_mz(Series::Z, 2, y3)));
        assert!(close(at(n - 4, 13), channel_mz(Series::Y, 3, y3)));
    }

    #[test]
    fn impossible_channels_masked() {
        let (form, mass) = parse_peptide("K.PEPTIDEKSAMPLER.R").unwrap();
        let pep = form.rt_peptide();
        let charge = 2;
        let hcd = activation_index("HCD-FTMS");
        let valid = channel_mask(&pep, charge, hcd, &ION_TYPES_V2);
        let sites = form.seq.len() - 1;
        assert_eq!(valid.len(), sites * N_ION_V2);
        let impossible = |k: usize| {
            let ion = ION_TYPES_V2[k];
            ion.charge > charge || ion.series.is_electron() || ion.loss == NeutralLoss::H3po4
        };
        for site in 0..sites {
            for (k, ion) in ION_TYPES_V2.iter().enumerate() {
                if impossible(k) {
                    assert!(!valid[site * N_ION_V2 + k], "site {site} channel {ion:?}");
                }
            }
        }
        let pred: Vec<f32> = (0..sites * N_ION_V2)
            .map(|i| if impossible(i % N_ION_V2) { 100.0 } else { 1.0 + (i % 7) as f32 })
            .collect();
        let frags = fragments(&pred, &ION_TYPES_V2, Some(valid.as_slice()), &mass, charge);
        assert!(!frags.is_empty());
        assert!(frags.iter().all(|f| !impossible(f.channel)));
        let weights: Vec<f64> = frags.iter().map(|f| f.pred).collect();
        let top = top_by_weight(&weights, 12);
        assert_eq!(top.len(), 12);
        assert!(top.iter().all(|&t| !impossible(frags[t].channel)));
        assert!(top.iter().all(|&t| frags[t].pred < 100.0));

        let etd = activation_index("ETD");
        let electron = channel_mask(&pep, charge, etd, &ION_TYPES_V2);
        let c1 = 15;
        assert!(electron[c1]);
        let before_proline = 1;
        assert!(!electron[before_proline * N_ION_V2 + c1]);
    }
}
