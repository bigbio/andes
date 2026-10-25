//! MS²PIP feature encoding for the exported XGBoost → ONNX tree ensembles.
//!
//! MS²PIP predicts each ion series with one XGBoost regressor applied to a
//! 139-value feature vector per cleavage site. This module computes those
//! vectors exactly as `ms2rescore_rs.ms2pip_compute_features` (ms2pip ≥ 4.2)
//! does, including the quirks the published models were trained with:
//!
//! - Modifications are ignored: a modified residue encodes as its base residue,
//!   and I / L share one index. Only fragment *masses* see the modifications;
//!   the intensity models never do.
//! - Positional feature `pi` at the first cleavage site reads the properties of
//!   residue index 0 (Ala), not of the first residue.
//! - Quartile indices use `floor(q * denom)` with `denom` = n − 1 for the
//!   peptide and the N-terminal ion but n for the C-terminal ion.
//!
//! Graph contract (`scripts/export_ms2pip_to_onnx.py`): input `features`
//! [rows, 139] f32, output `log2_intensity` [rows] = `log2(intensity + 0.001)`. One graph per ion series (`b`, `y`); row `i` of a
//! peptide is the cleavage after residue `i + 1`, i.e. b(i+1) and y(n−1−i).

use crate::error::{Error, Result};
use crate::tensor::NamedTensor;

/// Features per cleavage site.
pub const N_FEATURES: usize = 139;

/// Residue alphabet size (I and L share one index).
const N_AA: usize = 19;

/// Shortest peptide ms2pip predicts (`ms2pip._utils.psm_input`).
pub const MIN_LENGTH: usize = 4;

/// Longest peptide ms2pip predicts.
pub const MAX_LENGTH: usize = 100;

/// Lowest log2 prediction ms2pip keeps (`log2(0.001)`).
pub const LOG2_FLOOR: f32 = -9.965_784;

/// Basicity, helicity, hydrophobicity and pI per residue index.
const AA_PROPERTIES: [[u32; N_AA]; 4] = [
    [
        37, 35, 59, 129, 94, 0, 210, 81, 191, 106, 101, 117, 115, 343, 49, 90, 60, 134, 104,
    ],
    [
        68, 23, 33, 29, 70, 58, 41, 73, 32, 66, 38, 0, 40, 39, 44, 53, 71, 51, 55,
    ],
    [
        51, 75, 25, 35, 100, 16, 3, 94, 0, 82, 12, 0, 22, 22, 21, 39, 80, 98, 70,
    ],
    [
        32, 23, 0, 4, 27, 32, 48, 32, 69, 29, 26, 35, 28, 79, 29, 28, 31, 31, 28,
    ],
];

/// MS²PIP residue index (`ACDEFGHIKMNPQRSTVWY`, L → I); `None` outside it.
pub fn aa_index(c: char) -> Option<usize> {
    Some(match c.to_ascii_uppercase() {
        'A' => 0,
        'C' => 1,
        'D' => 2,
        'E' => 3,
        'F' => 4,
        'G' => 5,
        'H' => 6,
        'I' | 'L' => 7,
        'K' => 8,
        'M' => 9,
        'N' => 10,
        'P' => 11,
        'Q' => 12,
        'R' => 13,
        'S' => 14,
        'T' => 15,
        'V' => 16,
        'W' => 17,
        'Y' => 18,
        _ => return None,
    })
}

fn quartiles(sorted: &[u32], denom: usize) -> [u32; 5] {
    let n = sorted.len();
    if n == 0 {
        return [0; 5];
    }
    let at = |q: f64| sorted[((denom as f64 * q).floor() as usize).min(n - 1)];
    [sorted[0], at(0.25), at(0.5), at(0.75), sorted[n - 1]]
}

/// Feature rows for one peptide: `(len − 1) × 139` values, row-major.
///
/// `sequence` is the plain residue string (no modification text). A residue
/// outside the MS²PIP alphabet, a length outside
/// [`MIN_LENGTH`]`..=`[`MAX_LENGTH`] (ms2pip skips those peptides), or a
/// charge of 0 is an error.
pub fn compute_features(sequence: &str, charge: u8) -> Result<Vec<f32>> {
    let aa: Vec<usize> = sequence
        .chars()
        .map(|c| {
            aa_index(c).ok_or_else(|| {
                Error::InvalidArg(format!(
                    "residue `{c}` in peptide {sequence:?} is not in the MS²PIP alphabet"
                ))
            })
        })
        .collect::<Result<_>>()?;
    let n = aa.len();
    if !(MIN_LENGTH..=MAX_LENGTH).contains(&n) {
        return Err(Error::InvalidArg(format!(
            "peptide {sequence:?} has {n} residues; MS²PIP predicts {MIN_LENGTH}..={MAX_LENGTH}"
        )));
    }
    if charge == 0 {
        return Err(Error::InvalidArg("precursor charge must be ≥ 1".into()));
    }
    let z = charge as usize;
    let props: Vec<[u32; 4]> = aa
        .iter()
        .map(|&a| std::array::from_fn(|p| AA_PROPERTIES[p][a]))
        .collect();
    // 1-based buffer whose slot 0 holds residue index 0's properties.
    let buf: Vec<[u32; 4]> = std::iter::once(std::array::from_fn(|p| AA_PROPERTIES[p][0]))
        .chain(props.iter().copied())
        .collect();
    let pep_q: [[u32; 5]; 4] = std::array::from_fn(|p| {
        let mut v: Vec<u32> = props.iter().map(|x| x[p]).collect();
        v.sort_unstable();
        quartiles(&v, n - 1)
    });

    let mut out = vec![0f32; (n - 1) * N_FEATURES];
    let mut count_n = [0u32; N_AA];
    let mut count_c = [0u32; N_AA];
    for &a in &aa {
        count_c[a] += 1;
    }
    for i in 0..n - 1 {
        let f = &mut out[i * N_FEATURES..(i + 1) * N_FEATURES];
        f[0] = n as f32;
        f[1] = z as f32;
        for k in 0..5 {
            f[2 + k] = if k == 4 {
                (z >= 5) as u8 as f32
            } else {
                (z == k + 1) as u8 as f32
            };
        }
        for p in 0..4 {
            for q in 0..5 {
                f[7 + p * 5 + q] = pep_q[p][q] as f32;
            }
        }
        f[27] = (i + 1) as f32;
        f[28] = (n - i) as f32;
        count_n[aa[i]] += 1;
        count_c[aa[n - 1 - i]] -= 1;
        for a in 0..N_AA {
            f[29 + 2 * a] = count_n[a] as f32;
            f[30 + 2 * a] = count_c[a] as f32;
        }
        for p in 0..4 {
            let o = 67 + p * 18;
            f[o] = buf[1][p] as f32;
            f[o + 1] = buf[n][p] as f32;
            f[o + 2] = if i == 0 { 0.0 } else { buf[i - 1][p] as f32 };
            f[o + 3] = buf[i][p] as f32;
            f[o + 4] = buf[i + 1][p] as f32;
            f[o + 5] = buf[i + 2][p] as f32;

            let mut nv: Vec<u32> = props[..=i].iter().map(|x| x[p]).collect();
            f[o + 6] = nv.iter().sum::<u32>() as f32;
            nv.sort_unstable();
            for (k, v) in quartiles(&nv, i).into_iter().enumerate() {
                f[o + 7 + k] = v as f32;
            }
            let mut cv: Vec<u32> = props[i + 1..].iter().map(|x| x[p]).collect();
            f[o + 12] = cv.iter().sum::<u32>() as f32;
            cv.sort_unstable();
            for (k, v) in quartiles(&cv, n - i - 1).into_iter().enumerate() {
                f[o + 13 + k] = v as f32;
            }
        }
    }
    Ok(out)
}

/// A packed MS²PIP batch: all cleavage sites of all peptides as one matrix.
#[derive(Debug, Clone)]
pub struct Ms2pipBatch {
    /// `features` — [rows, 139] f32.
    pub features: NamedTensor,
    /// Row offset of each peptide; peptide `k` owns rows
    /// `offsets[k]..offsets[k + 1]` (`len − 1` rows).
    pub offsets: Vec<usize>,
}

/// Pack `(sequence, charge)` pairs into one feature matrix.
pub fn pack_batch(peptides: &[(&str, u8)]) -> Result<Ms2pipBatch> {
    if peptides.is_empty() {
        return Err(Error::InvalidArg("empty peptide batch".into()));
    }
    let mut data = Vec::new();
    let mut offsets = vec![0];
    for (seq, z) in peptides {
        data.extend(compute_features(seq, *z)?);
        offsets.push(data.len() / N_FEATURES);
    }
    let rows = data.len() / N_FEATURES;
    Ok(Ms2pipBatch {
        features: NamedTensor::new("features", vec![rows, N_FEATURES], data)?,
        offsets,
    })
}

/// The graph's `log2_intensity` output (`[rows]`), clipped at
/// [`LOG2_FLOOR`] as ms2pip does.
pub fn log2_predictions(out: &NamedTensor) -> Result<Vec<f32>> {
    if out.shape.len() != 1 {
        return Err(Error::InvalidArg(format!(
            "expected [rows] log2_intensity, got {:?}",
            out.shape
        )));
    }
    Ok(out.data.iter().map(|v| v.max(LOG2_FLOOR)).collect())
}

/// ms2pip's linear intensity from a log2 prediction: `2^x − 0.001`, ≥ 0.
pub fn to_linear(log2: f32) -> f32 {
    (log2.exp2() - 0.001).max(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// First cleavage site of `ACDE/2`, from the ms2rescore-rs C-reference test.
    #[rustfmt::skip]
    const ACDE_ION0: [u32; N_FEATURES] = [
        4, 2, 0, 1, 0, 0, 0, 35, 35, 37, 59, 129, 23, 23, 29, 33, 68, 25, 25,
        35, 51, 75, 0, 0, 4, 23, 32, 1, 4, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 37, 129, 0, 37, 37, 35, 37, 37, 37, 37, 37, 37, 223, 35, 35,
        59, 129, 129, 68, 29, 0, 68, 68, 23, 68, 68, 68, 68, 68, 68, 85, 23, 23,
        29, 33, 33, 51, 35, 0, 51, 51, 75, 51, 51, 51, 51, 51, 51, 135, 25, 25,
        35, 75, 75, 32, 4, 0, 32, 32, 23, 32, 32, 32, 32, 32, 32, 27, 0, 0, 4,
        23, 23,
    ];

    #[test]
    fn matches_c_reference_vector() {
        let f = compute_features("ACDE", 2).unwrap();
        assert_eq!(f.len(), 3 * N_FEATURES);
        for (k, (&got, &want)) in f[..N_FEATURES].iter().zip(&ACDE_ION0).enumerate() {
            assert_eq!(got, want as f32, "feature {k}");
        }
    }

    #[test]
    fn ion_lengths_and_shared_features() {
        let f = compute_features("ACDEF", 6).unwrap();
        assert_eq!((f[27], f[28]), (1.0, 5.0));
        assert_eq!((f[3 * N_FEATURES + 27], f[3 * N_FEATURES + 28]), (4.0, 2.0));
        assert_eq!(f[N_FEATURES + 1], 6.0);
        assert_eq!(f[6], 1.0, "charge ≥ 5 one-hot");
    }

    #[test]
    fn leucine_encodes_as_isoleucine() {
        assert_eq!(
            compute_features("PEPTLDE", 2).unwrap(),
            compute_features("PEPTIDE", 2).unwrap()
        );
    }

    #[test]
    fn rejects_unsupported_input() {
        assert!(compute_features("PEPXIDE", 2).is_err());
        assert!(compute_features("PEK", 2).is_err());
        assert!(compute_features("PEPTIDE", 0).is_err());
    }

    #[test]
    fn packs_offsets() {
        let b = pack_batch(&[("PEPTIDE", 2), ("ACDK", 1)]).unwrap();
        assert_eq!(b.offsets, vec![0, 6, 9]);
        assert_eq!(b.features.shape, vec![9, N_FEATURES]);
    }

    #[test]
    fn linear_intensity() {
        assert!(to_linear(LOG2_FLOOR) < 1e-9);
        assert!((to_linear(0.0) - 0.999).abs() < 1e-6);
    }
}
