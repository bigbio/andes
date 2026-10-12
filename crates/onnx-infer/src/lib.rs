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
//! [`peptdeep`] for PeptDeep ONNX I/O aligned with OpenMS exports.
//!
//! ## Design
//!
//! - **Backend-agnostic API** ([`OnnxSession`], [`InferenceBackend`]) so callers do
//!   not depend on tract vs ONNX Runtime.
//! - **Sidecar manifests** ([`ModelManifest`]) carry provenance, I/O names, task
//!   hints, and tool-specific `extras` (NCE, vocab path, …).
//! - **Feature-gated engines**: enable `--features tract` for the pure-Rust
//!   Apache-2.0 backend. An `ort` feature is reserved for broader op coverage
//!   (required for `ai.onnx.ml` TreeEnsemble and many Transformer graphs).
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

pub mod backend;
pub mod ms2pip;
pub mod peptdeep;
pub mod rescoring;

pub use error::{BackendKind, Error, Result};
pub use manifest::{ModelManifest, ModelTask, Provenance, TensorSpec};
pub use session::{InferenceBackend, OnnxSession};
pub use tensor::NamedTensor;
