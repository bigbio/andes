//! ONNX inference substrate: loads model manifests and runs sessions for the fragment-intensity
//! and retention-time predictors (Prosit, MS2PIP, DeepLC, AlphaPeptDeep and the native
//! predictors). Search and rescoring call it to obtain predicted spectra and retention times.

mod error;
mod manifest;
mod session;
mod tensor;

pub mod andes_peptide;
pub mod backend;
pub mod deeplc;
pub mod ms2_native;
pub mod ms2pip;
pub mod peptdeep;
pub mod predictors;
pub mod prosit;
pub mod rescoring;
pub mod rt_native;

pub use error::{BackendKind, Error, Result};
pub use manifest::{ModelManifest, ModelTask, Provenance, TensorSpec};
pub use session::{BackendOptions, InferenceBackend, OnnxSession};
pub use tensor::NamedTensor;
