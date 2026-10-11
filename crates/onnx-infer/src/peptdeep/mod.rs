//! AlphaPeptDeep / PeptDeep ONNX I/O contracts.
//!
//! Tensor names and shapes match OpenMS
//! [`export_peptdeep_models_to_onnx.py`](https://github.com/OpenMS/OpenMS/blob/develop/tools/scripts/export_peptdeep_models_to_onnx.py)
//! and the local copy under `scripts/export_peptdeep_to_onnx.py`:
//!
//! | Task | Inputs | Output |
//! |---|---|---|
//! | RT | `input_sequences` [B,L] i64, `mod_x` [B,L,M] f32 | `rt_pred` [B] |
//! | MS2 | `aa_indices`, `mod_x`, `charges`, `nce`, `instrument_indices` | `ms2_intensities` |
//! | CCS | `aa_indices`, `mod_x`, `charges` | `ccs_pred` |
//!
//! Encoding helpers here prepare tensors for [`crate::OnnxSession`]; they do
//! not download weights. Full PeptDeep AA/mod vocab tables live in model
//! `extras` / sidecar JSON once an export is run.

mod encode;
mod tensors;

pub use encode::{AA_VOCAB_1BASED, MOD_ELEMENT_COUNT, encode_aa_sequence};
pub use tensors::{
    INSTRUMENT_DEFAULT, NCE_DEFAULT, PeptDeepMs2Batch, PeptDeepRtBatch, pack_ms2_inputs,
    pack_ms2_inputs_with_mod_dim, pack_rt_inputs, pack_rt_inputs_with_mod_dim,
};
