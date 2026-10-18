//! # onnx-infer
//!
//! Tool-agnostic ONNX inference substrate for proteomics pipelines:
//! database search rescoring, de novo sequencing, intensity/RT predictors, etc.
//!
//! ## Role for quantms
//!
//! This crate is the intended Rust home for
//! [quantms-rescoring](https://github.com/bigbio/quantms-rescoring) inference
//! paths (AlphaPeptDeep / DeepLC / spectrum features → PIN / idparquet). See
//! [`rescoring`] for feature-name contracts and top-N batching, and
//! [`peptdeep`] for PeptDeep ONNX I/O aligned with the quantms-rescoring exports.
//!
//! ## Design
//!
//! - **Backend-agnostic API** ([`OnnxSession`], [`InferenceBackend`]) so callers do
//!   not depend on tract vs ONNX Runtime.
//! - **Sidecar manifests** ([`ModelManifest`]) carry provenance, I/O names, task
//!   hints, and tool-specific `extras` (NCE, vocab path, …).
//! - **Feature-gated engines**: `--features tract` for the pure-Rust
//!   Apache-2.0 backend; `--features ort` for ONNX Runtime (MIT), loaded as a
//!   shared library at run time (see [`backend::ort_backend`]). ONNX Runtime
//!   covers `ai.onnx.ml` (TreeEnsemble) and is several times faster on the
//!   bundled predictor graphs.
//!
//! ## Example (tract)
//!
//! ```ignore
//! use onnx_infer::{BackendKind, NamedTensor, OnnxSession};
//!
//! let mut session = OnnxSession::from_onnx("model.onnx", BackendKind::Tract)?;
//! let input = NamedTensor::new("x", vec![2, 4], vec![1.; 8])?;
//! let outputs = session.run(&[input])?;
//! ```
//!
//! Callers (search engines, de novo tools, quantms-rescoring) should run
//! inference on **retained top-N candidates only**, never inside the full
//! candidate loop. See [`rescoring::TopNPolicy`].

#![warn(missing_docs)]

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
pub mod prosit;
pub mod rescoring;
pub mod rt_native;

pub use error::{BackendKind, Error, Result};
pub use manifest::{ModelManifest, ModelTask, Provenance, TensorSpec};
pub use session::{BackendOptions, InferenceBackend, OnnxSession};
pub use tensor::NamedTensor;
