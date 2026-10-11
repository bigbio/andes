//! Pack PeptDeep RT / MS2 input tensors for [`crate::OnnxSession`].

use crate::error::{Error, Result};
use crate::peptdeep::encode::{MOD_ELEMENT_COUNT, encode_aa_sequence};
use crate::tensor::NamedTensor;

/// Default NCE when the search engine does not supply one (PeptDeep generic).
pub const NCE_DEFAULT: f32 = 30.0;

/// Default instrument index (0 = first bucket in PeptDeep instrument_group).
pub const INSTRUMENT_DEFAULT: i64 = 0;

/// One RT prediction batch (OpenMS export names).
#[derive(Debug, Clone)]
pub struct PeptDeepRtBatch {
    /// `input_sequences` — [B, L] as f32 for tract-friendly feeds (cast from i64).
    pub input_sequences: NamedTensor,
    /// `mod_x` — [B, L, M] f32.
    pub mod_x: NamedTensor,
}

/// One MS2 intensity batch (OpenMS export names).
#[derive(Debug, Clone)]
pub struct PeptDeepMs2Batch {
    /// `aa_indices` — [B, L] f32 (integer values).
    pub aa_indices: NamedTensor,
    /// `mod_x` — [B, L, M] f32.
    pub mod_x: NamedTensor,
    /// `charges` — [B, 1] f32.
    pub charges: NamedTensor,
    /// `nce` — [B, 1] f32.
    pub nce: NamedTensor,
    /// `instrument_indices` — [B] f32 (integer values).
    pub instrument_indices: NamedTensor,
}

impl PeptDeepRtBatch {
    /// Flatten to the slice expected by [`crate::OnnxSession::run`].
    pub fn as_inputs(&self) -> [NamedTensor; 2] {
        [self.input_sequences.clone(), self.mod_x.clone()]
    }
}

impl PeptDeepMs2Batch {
    /// Flatten to session inputs (owned clones).
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

/// Pack plain peptide sequences into an RT batch (zero modifications).
pub fn pack_rt_inputs(peptides: &[&str], seq_len: usize) -> Result<PeptDeepRtBatch> {
    pack_rt_inputs_with_mod_dim(peptides, seq_len, MOD_ELEMENT_COUNT)
}

/// Pack RT inputs with an explicit mod channel count (from manifest extras).
pub fn pack_rt_inputs_with_mod_dim(
    peptides: &[&str],
    seq_len: usize,
    mod_dim: usize,
) -> Result<PeptDeepRtBatch> {
    let b = peptides.len();
    if b == 0 {
        return Err(Error::InvalidArg("empty peptide batch".into()));
    }
    let mut seq = Vec::with_capacity(b * seq_len);
    for p in peptides {
        let idx = encode_aa_sequence(p, seq_len)?;
        seq.extend(idx.into_iter().map(|v| v as f32));
    }
    let mod_x = vec![0f32; b * seq_len * mod_dim];
    Ok(PeptDeepRtBatch {
        input_sequences: NamedTensor::new("input_sequences", vec![b, seq_len], seq)?,
        mod_x: NamedTensor::new("mod_x", vec![b, seq_len, mod_dim], mod_x)?,
    })
}

/// Pack MS2 inputs (zero mods, shared NCE / instrument).
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

/// Pack MS2 inputs with explicit mod channel count.
pub fn pack_ms2_inputs_with_mod_dim(
    peptides: &[&str],
    charges: &[u8],
    seq_len: usize,
    mod_dim: usize,
    nce: f32,
    instrument: i64,
) -> Result<PeptDeepMs2Batch> {
    let b = peptides.len();
    if b == 0 {
        return Err(Error::InvalidArg("empty peptide batch".into()));
    }
    if charges.len() != b {
        return Err(Error::InvalidArg(format!(
            "charges len {} != peptides len {b}",
            charges.len()
        )));
    }
    let mut aa = Vec::with_capacity(b * seq_len);
    for p in peptides {
        let idx = encode_aa_sequence(p, seq_len)?;
        aa.extend(idx.into_iter().map(|v| v as f32));
    }
    let mod_x = vec![0f32; b * seq_len * mod_dim];
    let charges_f: Vec<f32> = charges.iter().map(|&c| c as f32).collect();
    let nce_f = vec![nce; b];
    let inst: Vec<f32> = vec![instrument as f32; b];
    Ok(PeptDeepMs2Batch {
        aa_indices: NamedTensor::new("aa_indices", vec![b, seq_len], aa)?,
        mod_x: NamedTensor::new("mod_x", vec![b, seq_len, mod_dim], mod_x)?,
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
        let batch = pack_rt_inputs(&["PEPTIDE", "SAMPLER"], 12).unwrap();
        assert_eq!(batch.input_sequences.shape, vec![2, 12]);
        assert_eq!(batch.mod_x.shape, vec![2, 12, MOD_ELEMENT_COUNT]);
    }

    #[test]
    fn packs_ms2_batch_shapes() {
        let batch = pack_ms2_inputs(&["PEPTIDE"], &[2], 10, NCE_DEFAULT, INSTRUMENT_DEFAULT)
            .unwrap();
        assert_eq!(batch.aa_indices.shape, vec![1, 10]);
        assert_eq!(batch.charges.data, vec![2.0]);
        assert_eq!(batch.nce.data, vec![NCE_DEFAULT]);
    }
}
