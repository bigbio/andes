//! Errors for the ONNX inference substrate.

#![allow(missing_docs)]

use std::path::PathBuf;

use thiserror::Error;

/// Fallible operations in this crate.
#[derive(Debug, Error)]
pub enum Error {
    #[error("IO error for {path}: {source}")]
    Io {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },

    #[error("failed to parse model manifest {path}: {source}")]
    ManifestParse {
        path: PathBuf,
        #[source]
        source: serde_json::Error,
    },

    #[error("manifest validation failed: {0}")]
    ManifestInvalid(String),

    #[error("invalid argument: {0}")]
    InvalidArg(String),

    #[error("ONNX backend {backend:?} is not enabled; rebuild with `--features {feature}`")]
    BackendDisabled {
        backend: BackendKind,
        feature: &'static str,
    },

    #[error("ONNX backend {0:?} is not implemented yet")]
    BackendUnimplemented(BackendKind),

    #[error("ONNX load/run failed: {0}")]
    Inference(String),

    #[error("unknown input name `{name}` (model inputs: {available:?})")]
    UnknownInput {
        name: String,
        available: Vec<String>,
    },

    #[error("missing required input `{name}`")]
    MissingInput { name: String },

    #[error(
        "tensor `{name}` has shape {got:?} (nelem={got_nelem}) but expected {expected_nelem} elements"
    )]
    ShapeMismatch {
        name: String,
        got: Vec<usize>,
        got_nelem: usize,
        expected_nelem: usize,
    },
}

/// Which inference engine to use.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BackendKind {
    /// Pure-Rust [`tract`](https://github.com/sonos/tract) (`--features tract`).
    Tract,
    /// ONNX Runtime via `ort` (planned; broader operator coverage).
    Ort,
}

impl BackendKind {
    /// Cargo feature that enables this backend.
    pub fn feature_name(self) -> &'static str {
        match self {
            Self::Tract => "tract",
            Self::Ort => "ort",
        }
    }
}

/// Result alias for this crate.
pub type Result<T> = std::result::Result<T, Error>;
