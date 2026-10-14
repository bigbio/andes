//! Inference backends (feature-gated).
//!
//! Enable with cargo features `tract` (pure Rust) and/or `ort` (ONNX Runtime).

#[cfg(feature = "ort")]
pub mod ort_backend;
#[cfg(feature = "tract")]
pub mod tract_backend;
