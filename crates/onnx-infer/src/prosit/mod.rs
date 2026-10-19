//! Prosit input encoding and output decoding for the exported Prosit 2020
//! intensity (HCD) and Prosit 2019 iRT graphs.
//!
//! Mirrors Koina's `Prosit_Preprocess_peptide` / `_charge` /
//! `_collision_energy` and `Prosit_2019_intensity_postprocess`:
//!
//! - `peptides_in` [B, 30]: residue codes `A=1 C=2 D=3 … Y=20`, oxidised
//!   methionine `M[UNIMOD:35]` = 21, zero padding. Code 2 means
//!   carbamidomethylated cysteine: the models were not trained on free
//!   cysteine, so an unmodified C is rejected here.
//! - `precursor_charge_in` [B, 6]: one-hot charge 1..=6.
//! - `collision_energy_in` [B, 1]: NCE / 100.
//! - intensity output [B, 174]: 29 positions × (y+1, y+2, y+3, b+1, b+2, b+3),
//!   position `k` holding y(k+1) and b(k+1).
//!
//! Only carbamidomethyl-C and oxidised M are supported (Prosit 2020's
//! modification set); any other modification, a residue outside the
//! alphabet, a length over 30 or a charge outside 1..=6 is an error.

use crate::error::{Error, Result};
use crate::tensor::NamedTensor;

/// Encoded sequence length.
pub const SEQ_LEN: usize = 30;
/// Intensity output values per peptide.
pub const N_OUT: usize = 174;
/// Highest precursor charge the one-hot covers.
pub const MAX_CHARGE: u8 = 6;

/// Residue modification states Prosit 2020 can encode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ResidueState {
    /// Unmodified residue.
    Plain,
    /// Carbamidomethylated cysteine (`C[UNIMOD:4]`).
    Carbamidomethyl,
    /// Oxidised methionine (`M[UNIMOD:35]`).
    Oxidation,
}

fn code(c: char, state: ResidueState) -> Option<u8> {
    use ResidueState::*;
    Some(match (c, state) {
        ('C', Carbamidomethyl) => 2,
        ('M', Oxidation) => 21,
        (_, Plain) => match c {
            'A' => 1,
            'D' => 3,
            'E' => 4,
            'F' => 5,
            'G' => 6,
            'H' => 7,
            'I' => 8,
            'K' => 9,
            'L' => 10,
            'M' => 11,
            'N' => 12,
            'P' => 13,
            'Q' => 14,
            'R' => 15,
            'S' => 16,
            'T' => 17,
            'V' => 18,
            'W' => 19,
            'Y' => 20,
            _ => return None,
        },
        _ => return None,
    })
}

/// Encode residues with their per-residue modification state (same length).
pub fn encode_sequence(sequence: &str, states: &[ResidueState]) -> Result<[f32; SEQ_LEN]> {
    let n = sequence.chars().count();
    if n == 0 || n > SEQ_LEN || states.len() != n {
        return Err(Error::InvalidArg(format!(
            "peptide {sequence:?}: length {n} (states {}) outside Prosit's 1..={SEQ_LEN}",
            states.len()
        )));
    }
    let mut out = [0f32; SEQ_LEN];
    for (i, (c, s)) in sequence.chars().zip(states).enumerate() {
        out[i] = code(c, *s).ok_or_else(|| {
            Error::InvalidArg(format!(
                "residue `{c}` ({s:?}) at {} of {sequence:?} is not in the Prosit 2020 alphabet \
                 (free cysteine and modifications other than C[UNIMOD:4] / M[UNIMOD:35] are unsupported)",
                i + 1
            ))
        })? as f32;
    }
    Ok(out)
}

/// One Prosit peptide: residues and one state per residue.
#[derive(Debug, Clone)]
pub struct PrositPeptide<'a> {
    /// Residues.
    pub sequence: &'a str,
    /// Modification state per residue.
    pub states: Vec<ResidueState>,
}

fn sequences(peptides: &[PrositPeptide<'_>]) -> Result<Vec<f32>> {
    if peptides.is_empty() {
        return Err(Error::InvalidArg("empty peptide batch".into()));
    }
    let mut v = Vec::with_capacity(peptides.len() * SEQ_LEN);
    for p in peptides {
        v.extend(encode_sequence(p.sequence, &p.states)?);
    }
    Ok(v)
}

/// iRT graph input `peptides_in`.
pub fn pack_irt_batch(peptides: &[PrositPeptide<'_>]) -> Result<NamedTensor> {
    NamedTensor::new("peptides_in", vec![peptides.len(), SEQ_LEN], sequences(peptides)?)
}

/// Intensity graph inputs (`peptides_in`, `precursor_charge_in`, `collision_energy_in`).
pub fn pack_intensity_batch(
    peptides: &[PrositPeptide<'_>],
    charges: &[u8],
    nce: f32,
) -> Result<[NamedTensor; 3]> {
    let b = peptides.len();
    if charges.len() != b {
        return Err(Error::InvalidArg(format!("charges len {} != peptides len {b}", charges.len())));
    }
    let mut onehot = vec![0f32; b * MAX_CHARGE as usize];
    for (i, &z) in charges.iter().enumerate() {
        if !(1..=MAX_CHARGE).contains(&z) {
            return Err(Error::InvalidArg(format!("charge {z} outside Prosit's 1..={MAX_CHARGE}")));
        }
        onehot[i * MAX_CHARGE as usize + z as usize - 1] = 1.0;
    }
    Ok([
        NamedTensor::new("peptides_in", vec![b, SEQ_LEN], sequences(peptides)?)?,
        NamedTensor::new("precursor_charge_in", vec![b, MAX_CHARGE as usize], onehot)?,
        NamedTensor::new("collision_energy_in", vec![b, 1], vec![nce / 100.0; b])?,
    ])
}

/// Per cleavage site `i` (b(i+1), y(n−1−i)) of an `n`-residue peptide: the
/// six intensities `[b+1, b+2, b+3, y+1, y+2, y+3]`, as Koina post-processes
/// them: negatives → 0, fragment charges above the precursor charge → 0, then
/// scaled to the largest value.
pub fn decode_intensity(raw: &[f32], n: usize, charge: u8) -> Result<Vec<f32>> {
    if raw.len() != N_OUT || !(2..=SEQ_LEN).contains(&n) {
        return Err(Error::InvalidArg(format!(
            "decode: {} outputs for an {n}-residue peptide",
            raw.len()
        )));
    }
    let v = |pos: usize, k: usize, z: usize| -> f32 {
        if z as u8 > charge { 0.0 } else { raw[pos * 6 + k * 3 + z - 1].max(0.0) }
    };
    let mut out = Vec::with_capacity((n - 1) * 6);
    for i in 0..n - 1 {
        for z in 1..=3 {
            out.push(v(i, 1, z));
        }
        for z in 1..=3 {
            out.push(v(n - 2 - i, 0, z));
        }
    }
    let max = out.iter().copied().fold(0f32, f32::max);
    if max > 0.0 {
        out.iter_mut().for_each(|x| *x /= max);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ResidueState::*;

    #[test]
    fn encodes_alphabet_and_mods() {
        let s = encode_sequence("MCK", &[Oxidation, Carbamidomethyl, Plain]).unwrap();
        assert_eq!(&s[..4], &[21.0, 2.0, 9.0, 0.0]);
        assert!(encode_sequence("CK", &[Plain, Plain]).is_err(), "free cysteine");
        assert!(encode_sequence("SK", &[Oxidation, Plain]).is_err());
        assert!(encode_sequence(&"A".repeat(31), &[Plain; 31]).is_err());
    }

    #[test]
    fn decodes_sites_and_masks_charge() {
        let mut raw = vec![0f32; N_OUT];
        // position 0: y1+1 = 1, b1+1 = 2, b1+2 = 4 (masked at charge 1)
        raw[0] = 1.0;
        raw[3] = 2.0;
        raw[4] = 4.0;
        let out = decode_intensity(&raw, 3, 1).unwrap();
        // site 0 = b1 / y2, site 1 = b2 / y1
        assert_eq!(&out[0..6], &[1.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        assert_eq!(&out[6..12], &[0.0, 0.0, 0.0, 0.5, 0.0, 0.0]);
        let out = decode_intensity(&raw, 3, 2).unwrap();
        assert_eq!(out[1], 1.0);
    }

    #[test]
    fn packs_charge_onehot_and_nce() {
        let p = PrositPeptide { sequence: "PEPTIDEK", states: vec![Plain; 8] };
        let [s, z, ce] = pack_intensity_batch(&[p], &[3], 30.0).unwrap();
        assert_eq!(s.shape, vec![1, SEQ_LEN]);
        assert_eq!(z.data, vec![0.0, 0.0, 1.0, 0.0, 0.0, 0.0]);
        assert!((ce.data[0] - 0.3).abs() < 1e-7);
        assert!(pack_intensity_batch(&[PrositPeptide { sequence: "PK", states: vec![Plain; 2] }], &[7], 30.0).is_err());
    }
}
