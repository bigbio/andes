//! Inference backends (feature-gated).
//!
//! Enable with cargo features (`tract`, later `ort`).

#[cfg(feature = "tract")]
pub mod tract_backend;
