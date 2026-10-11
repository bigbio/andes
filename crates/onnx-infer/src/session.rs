//! Backend-agnostic ONNX session API.

use std::path::Path;

use crate::error::{BackendKind, Error, Result};
use crate::manifest::ModelManifest;
use crate::tensor::NamedTensor;

/// Runs one ONNX model. Implementations are behind cargo features (`tract`, …).
pub trait InferenceBackend: Send {
    /// Engine kind.
    fn kind(&self) -> BackendKind;

    /// Graph input names in model order.
    fn input_names(&self) -> &[String];

    /// Graph output names in model order.
    fn output_names(&self) -> &[String];

    /// Execute the model. `inputs` are matched by [`NamedTensor::name`].
    fn run(&mut self, inputs: &[NamedTensor]) -> Result<Vec<NamedTensor>>;
}

/// A loaded model: optional sidecar manifest + live backend.
pub struct OnnxSession {
    manifest: Option<ModelManifest>,
    backend: Box<dyn InferenceBackend>,
}

impl OnnxSession {
    /// Load weights from `onnx_path` with the given backend.
    pub fn from_onnx(onnx_path: impl AsRef<Path>, backend: BackendKind) -> Result<Self> {
        let backend = load_backend(onnx_path.as_ref(), backend)?;
        Ok(Self {
            manifest: None,
            backend,
        })
    }

    /// Load from a JSON manifest (resolves the ONNX path, keeps provenance).
    pub fn from_manifest(manifest_path: impl AsRef<Path>, backend: BackendKind) -> Result<Self> {
        let manifest = ModelManifest::from_path(manifest_path)?;
        let backend = load_backend(manifest.onnx_path(), backend)?;
        Ok(Self {
            manifest: Some(manifest),
            backend,
        })
    }

    /// Sidecar metadata, if loaded via [`Self::from_manifest`].
    pub fn manifest(&self) -> Option<&ModelManifest> {
        self.manifest.as_ref()
    }

    /// Which engine is running.
    pub fn backend_kind(&self) -> BackendKind {
        self.backend.kind()
    }

    /// Graph input names.
    pub fn input_names(&self) -> &[String] {
        self.backend.input_names()
    }

    /// Graph output names.
    pub fn output_names(&self) -> &[String] {
        self.backend.output_names()
    }

    /// Run inference (batched when the leading axis is batch).
    pub fn run(&mut self, inputs: &[NamedTensor]) -> Result<Vec<NamedTensor>> {
        self.backend.run(inputs)
    }

    /// Convenience: single named input → first output tensor.
    pub fn run_one(&mut self, input: NamedTensor) -> Result<NamedTensor> {
        let outs = self.run(std::slice::from_ref(&input))?;
        outs.into_iter().next().ok_or_else(|| {
            Error::Inference("model produced zero outputs".into())
        })
    }
}

fn load_backend(onnx_path: &Path, backend: BackendKind) -> Result<Box<dyn InferenceBackend>> {
    match backend {
        BackendKind::Tract => {
            #[cfg(feature = "tract")]
            {
                Ok(Box::new(crate::backend::tract_backend::TractSession::load(
                    onnx_path,
                )?))
            }
            #[cfg(not(feature = "tract"))]
            {
                let _ = onnx_path;
                Err(Error::BackendDisabled {
                    backend: BackendKind::Tract,
                    feature: "tract",
                })
            }
        }
        BackendKind::Ort => {
            let _ = onnx_path;
            #[cfg(feature = "ort")]
            {
                Err(Error::BackendUnimplemented(BackendKind::Ort))
            }
            #[cfg(not(feature = "ort"))]
            {
                Err(Error::BackendDisabled {
                    backend: BackendKind::Ort,
                    feature: "ort",
                })
            }
        }
    }
}
