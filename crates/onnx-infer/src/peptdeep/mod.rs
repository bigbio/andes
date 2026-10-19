//! AlphaPeptDeep / PeptDeep ONNX I/O contracts.
//!
//! Tensor names and shapes match the reference PeptDeep ONNX export, adapted
//! locally as `scripts/export_peptdeep_to_onnx.py`:
//!
//! | Task | Inputs | Output |
//! |---|---|---|
//! | RT | `input_sequences` [B,L+2] i64, `mod_x` [B,L+2,M] f32 | `rt_pred` [B] |
//! | MS2 | `aa_indices`, `mod_x`, `charges`, `nce`, `instrument_indices` | `ms2_intensities` [B,L-1,8] |
//! | CCS | `aa_indices`, `mod_x`, `charges` | `ccs_pred` |
//!
//! Encoding helpers here prepare tensors for [`crate::OnnxSession`]; they do
//! not download weights. Full PeptDeep AA/mod vocab tables live in model
//! `extras` / sidecar JSON once an export is run.

mod encode;
mod tensors;

pub use encode::{aa_index, encode_aa_sequence, MOD_ELEMENT_COUNT};
pub use tensors::{
    pack_ms2_batch, pack_ms2_inputs, pack_ms2_inputs_with_mod_dim, pack_rt_batch, pack_rt_inputs,
    pack_rt_inputs_with_mod_dim, PeptDeepMs2Batch, PeptDeepPeptide, PeptDeepRtBatch, CHARGE_FACTOR,
    INSTRUMENT_DEFAULT, NCE_DEFAULT, NCE_FACTOR,
};
