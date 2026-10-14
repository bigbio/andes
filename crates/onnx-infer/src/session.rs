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

    /// A second handle on the same loaded model (weights shared, not copied),
    /// for running one model from several threads. `None` when unsupported.
    fn share(&self) -> Option<Box<dyn InferenceBackend>> {
        None
    }

    /// Execute the model. `inputs` are matched by [`NamedTensor::name`].
    fn run(&mut self, inputs: &[NamedTensor]) -> Result<Vec<NamedTensor>>;
}

/// Engine settings applied when a model is loaded.
#[derive(Clone, Debug)]
pub struct BackendOptions {
    /// Operator threads inside one `run` call. tract always runs on the
    /// calling thread and ignores this; ONNX Runtime uses it as the session's
    /// intra-op pool size (inter-op is fixed at 1). The default, 1, suits
    /// callers that parallelise across batches with [`OnnxSession::share`].
    pub intra_threads: usize,
}

impl Default for BackendOptions {
    fn default() -> Self {
        Self { intra_threads: 1 }
    }
}

/// A loaded model: optional sidecar manifest + live backend.
pub struct OnnxSession {
    manifest: Option<ModelManifest>,
    backend: Box<dyn InferenceBackend>,
}

impl OnnxSession {
    /// Load weights from `onnx_path` with the given backend.
    pub fn from_onnx(onnx_path: impl AsRef<Path>, backend: BackendKind) -> Result<Self> {
        Self::from_onnx_with(onnx_path, backend, &BackendOptions::default())
    }

    /// [`Self::from_onnx`] with explicit engine settings.
    pub fn from_onnx_with(
        onnx_path: impl AsRef<Path>,
        backend: BackendKind,
        opts: &BackendOptions,
    ) -> Result<Self> {
        let backend = load_backend(onnx_path.as_ref(), backend, opts)?;
        Ok(Self {
            manifest: None,
            backend,
        })
    }

    /// Load from a JSON manifest (resolves the ONNX path, keeps provenance).
    pub fn from_manifest(manifest_path: impl AsRef<Path>, backend: BackendKind) -> Result<Self> {
        Self::from_manifest_with(manifest_path, backend, &BackendOptions::default())
    }

    /// [`Self::from_manifest`] with explicit engine settings.
    pub fn from_manifest_with(
        manifest_path: impl AsRef<Path>,
        backend: BackendKind,
        opts: &BackendOptions,
    ) -> Result<Self> {
        let manifest = ModelManifest::from_path(manifest_path)?;
        let backend = load_backend(manifest.onnx_path(), backend, opts)?;
        Ok(Self {
            manifest: Some(manifest),
            backend,
        })
    }

    /// Another session on the same loaded weights (no reload, no copy), so
    /// worker threads can share one large model.
    pub fn share(&self) -> Result<Self> {
        let backend = self.backend.share().ok_or_else(|| {
            Error::Inference(format!(
                "{:?} backend cannot share a loaded model",
                self.backend.kind()
            ))
        })?;
        Ok(Self {
            manifest: self.manifest.clone(),
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
        outs.into_iter()
            .next()
            .ok_or_else(|| Error::Inference("model produced zero outputs".into()))
    }
}

fn load_backend(
    onnx_path: &Path,
    backend: BackendKind,
    opts: &BackendOptions,
) -> Result<Box<dyn InferenceBackend>> {
    let _ = (onnx_path, opts);
    match backend {
        #[cfg(feature = "tract")]
        BackendKind::Tract => Ok(Box::new(crate::backend::tract_backend::TractSession::load(
            onnx_path,
        )?)),
        #[cfg(feature = "ort")]
        BackendKind::Ort => Ok(Box::new(crate::backend::ort_backend::OrtSession::load(
            onnx_path, opts,
        )?)),
        #[allow(unreachable_patterns)]
        other => Err(Error::BackendDisabled {
            backend: other,
            feature: other.feature_name(),
        }),
    }
}
