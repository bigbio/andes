use crate::error::{Error, Result};
use crate::peptdeep::encode::{encode_aa_sequence, MOD_ELEMENT_COUNT};
use crate::tensor::{checked_element_count, NamedTensor};

pub const NCE_DEFAULT: f32 = 30.0;

pub const INSTRUMENT_DEFAULT: i64 = 0;

pub const CHARGE_FACTOR: f32 = 0.1;

pub const NCE_FACTOR: f32 = 0.01;

#[derive(Debug, Clone, Default)]
pub struct PeptDeepPeptide<'a> {
    pub sequence: &'a str,
    pub mods: Vec<(usize, &'a [f32])>,
}

#[derive(Debug, Clone)]
pub struct PeptDeepRtBatch {
    pub input_sequences: NamedTensor,
    pub mod_x: NamedTensor,
}

#[derive(Debug, Clone)]
pub struct PeptDeepMs2Batch {
    pub aa_indices: NamedTensor,
    pub mod_x: NamedTensor,
    pub charges: NamedTensor,
    pub nce: NamedTensor,
    pub instrument_indices: NamedTensor,
}

impl PeptDeepRtBatch {
    pub fn as_inputs(&self) -> [NamedTensor; 2] {
        [self.input_sequences.clone(), self.mod_x.clone()]
    }
}

impl PeptDeepMs2Batch {
    pub fn as_inputs(&self) -> [NamedTensor; 5] {
        [
            self.aa_indices.clone(),
            self.mod_x.clone(),
            self.charges.clone(),
            self.nce.clone(),
            self.instrument_indices.clone(),
        ]
    }
}

fn pack_sequences(
    peptides: &[PeptDeepPeptide<'_>],
    seq_len: usize,
    mod_dim: usize,
) -> Result<(Vec<f32>, Vec<f32>)> {
    if peptides.is_empty() {
        return Err(Error::InvalidArg("empty peptide batch".into()));
    }
    let width = seq_len
        .checked_add(2)
        .ok_or_else(|| Error::InvalidArg(format!("seq_len {seq_len} overflows usize")))?;
    let aa_len = checked_element_count(&[peptides.len(), width], "peptdeep aa_indices")?;
    let mod_len = checked_element_count(&[peptides.len(), width, mod_dim], "peptdeep mod_x")?;
    let mut aa = Vec::with_capacity(aa_len);
    let mut mod_x = vec![0f32; mod_len];
    for (b, p) in peptides.iter().enumerate() {
        let idx = encode_aa_sequence(p.sequence, seq_len)?;
        aa.extend(idx.into_iter().map(|v| v as f32));
        let n = p.sequence.chars().count();
        for &(site, feat) in &p.mods {
            if site > n + 1 {
                return Err(Error::InvalidArg(format!(
                    "mod site {site} outside peptide {:?} (n = {n})",
                    p.sequence
                )));
            }
            if feat.len() != mod_dim {
                return Err(Error::InvalidArg(format!(
                    "mod feature length {} != mod dim {mod_dim}",
                    feat.len()
                )));
            }
            let base = (b * width + site) * mod_dim;
            for (dst, src) in mod_x[base..base + mod_dim].iter_mut().zip(feat) {
                *dst += *src;
            }
        }
    }
    Ok((aa, mod_x))
}

pub fn pack_rt_inputs(peptides: &[&str], seq_len: usize) -> Result<PeptDeepRtBatch> {
    pack_rt_inputs_with_mod_dim(peptides, seq_len, MOD_ELEMENT_COUNT)
}

pub fn pack_rt_inputs_with_mod_dim(
    peptides: &[&str],
    seq_len: usize,
    mod_dim: usize,
) -> Result<PeptDeepRtBatch> {
    let peps: Vec<PeptDeepPeptide<'_>> = peptides
        .iter()
        .map(|s| PeptDeepPeptide {
            sequence: s,
            mods: Vec::new(),
        })
        .collect();
    pack_rt_batch(&peps, seq_len, mod_dim)
}

pub fn pack_rt_batch(
    peptides: &[PeptDeepPeptide<'_>],
    seq_len: usize,
    mod_dim: usize,
) -> Result<PeptDeepRtBatch> {
    let b = peptides.len();
    let (aa, mod_x) = pack_sequences(peptides, seq_len, mod_dim)?;
    Ok(PeptDeepRtBatch {
        input_sequences: NamedTensor::new("input_sequences", vec![b, seq_len + 2], aa)?,
        mod_x: NamedTensor::new("mod_x", vec![b, seq_len + 2, mod_dim], mod_x)?,
    })
}

pub fn pack_ms2_inputs(
    peptides: &[&str],
    charges: &[u8],
    seq_len: usize,
    nce: f32,
    instrument: i64,
) -> Result<PeptDeepMs2Batch> {
    pack_ms2_inputs_with_mod_dim(
        peptides,
        charges,
        seq_len,
        MOD_ELEMENT_COUNT,
        nce,
        instrument,
    )
}

pub fn pack_ms2_inputs_with_mod_dim(
    peptides: &[&str],
    charges: &[u8],
    seq_len: usize,
    mod_dim: usize,
    nce: f32,
    instrument: i64,
) -> Result<PeptDeepMs2Batch> {
    let peps: Vec<PeptDeepPeptide<'_>> = peptides
        .iter()
        .map(|s| PeptDeepPeptide {
            sequence: s,
            mods: Vec::new(),
        })
        .collect();
    pack_ms2_batch(&peps, charges, seq_len, mod_dim, nce, instrument)
}

pub fn pack_ms2_batch(
    peptides: &[PeptDeepPeptide<'_>],
    charges: &[u8],
    seq_len: usize,
    mod_dim: usize,
    nce: f32,
    instrument: i64,
) -> Result<PeptDeepMs2Batch> {
    let b = peptides.len();
    if charges.len() != b {
        return Err(Error::InvalidArg(format!(
            "charges len {} != peptides len {b}",
            charges.len()
        )));
    }
    let (aa, mod_x) = pack_sequences(peptides, seq_len, mod_dim)?;
    let charges_f: Vec<f32> = charges.iter().map(|&c| c as f32 * CHARGE_FACTOR).collect();
    let nce_f = vec![nce * NCE_FACTOR; b];
    let inst: Vec<f32> = vec![instrument as f32; b];
    Ok(PeptDeepMs2Batch {
        aa_indices: NamedTensor::new("aa_indices", vec![b, seq_len + 2], aa)?,
        mod_x: NamedTensor::new("mod_x", vec![b, seq_len + 2, mod_dim], mod_x)?,
        charges: NamedTensor::new("charges", vec![b, 1], charges_f)?,
        nce: NamedTensor::new("nce", vec![b, 1], nce_f)?,
        instrument_indices: NamedTensor::new("instrument_indices", vec![b], inst)?,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn packs_rt_batch_shapes() {
        let batch = pack_rt_inputs(&["PEPTIDE", "SAMPLER"], 7).unwrap();
        assert_eq!(batch.input_sequences.shape, vec![2, 9]);
        assert_eq!(batch.mod_x.shape, vec![2, 9, MOD_ELEMENT_COUNT]);
    }

    #[test]
    fn packs_ms2_batch_with_peptdeep_scaling() {
        let batch =
            pack_ms2_inputs(&["PEPTIDE"], &[2], 7, NCE_DEFAULT, INSTRUMENT_DEFAULT).unwrap();
        assert_eq!(batch.aa_indices.shape, vec![1, 9]);
        assert!((batch.charges.data[0] - 0.2).abs() < 1e-6);
        assert!((batch.nce.data[0] - 0.3).abs() < 1e-6);
    }

    #[test]
    fn mod_features_land_on_their_site_and_sum() {
        let f = [1.0f32, 2.0, 0.0];
        let pep = PeptDeepPeptide {
            sequence: "ACK",
            mods: vec![(0, &f[..]), (2, &f[..]), (2, &f[..])],
        };
        let batch = pack_rt_batch(&[pep], 3, 3).unwrap();
        let m = &batch.mod_x.data;
        assert_eq!(&m[0..3], &[1.0, 2.0, 0.0]);
        assert_eq!(&m[6..9], &[2.0, 4.0, 0.0]);
        assert!(m[3..6].iter().chain(&m[9..]).all(|&v| v == 0.0));
    }

    #[test]
    fn rejects_mod_site_past_c_terminus() {
        let f = [1.0f32];
        let pep = PeptDeepPeptide {
            sequence: "AK",
            mods: vec![(4, &f[..])],
        };
        assert!(pack_rt_batch(&[pep], 2, 1).is_err());
    }
}
