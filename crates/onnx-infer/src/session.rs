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

    /// Declared shapes of the graph inputs, parallel to [`Self::input_names`]:
    /// `None` when the rank is unknown, else one entry per axis with `None`
    /// for a dynamic axis. Empty when the engine does not expose shapes.
    fn input_shapes(&self) -> Vec<Option<Vec<Option<usize>>>> {
        Vec::new()
    }

    /// Declared shapes of the graph outputs, as [`Self::input_shapes`].
    fn output_shapes(&self) -> Vec<Option<Vec<Option<usize>>>> {
        Vec::new()
    }

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
        check_manifest_io(
            &manifest,
            Side::Input,
            backend.input_names(),
            &backend.input_shapes(),
        )?;
        check_manifest_io(
            &manifest,
            Side::Output,
            backend.output_names(),
            &backend.output_shapes(),
        )?;
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

#[derive(Clone, Copy, Debug)]
enum Side {
    Input,
    Output,
}

/// Check one side of a manifest against the loaded graph. An empty list in the
/// manifest declares nothing and is not checked. Otherwise the declared names
/// must be exactly the graph's names (order free), and where both the manifest
/// and the graph state a shape, the ranks and every axis both fix must agree.
fn check_manifest_io(
    manifest: &ModelManifest,
    side: Side,
    graph_names: &[String],
    graph_shapes: &[Option<Vec<Option<usize>>>],
) -> Result<()> {
    let (specs, what) = match side {
        Side::Input => (&manifest.inputs, "input"),
        Side::Output => (&manifest.outputs, "output"),
    };
    if specs.is_empty() {
        return Ok(());
    }
    let missing: Vec<&str> = specs
        .iter()
        .map(|s| s.name.as_str())
        .filter(|n| !graph_names.iter().any(|g| g == n))
        .collect();
    let undeclared: Vec<&str> = graph_names
        .iter()
        .map(String::as_str)
        .filter(|g| !specs.iter().any(|s| s.name == *g))
        .collect();
    if !missing.is_empty() || !undeclared.is_empty() {
        return Err(Error::ManifestInvalid(format!(
            "manifest `{}` {what}s do not match {}: declared but not in the graph {missing:?}, \
             in the graph but not declared {undeclared:?}",
            manifest.id,
            manifest.onnx_path().display()
        )));
    }
    for spec in specs {
        let Some(declared) = &spec.shape else {
            continue;
        };
        let Some(graph) = graph_names
            .iter()
            .position(|g| *g == spec.name)
            .and_then(|i| graph_shapes.get(i))
            .and_then(Option::as_ref)
        else {
            continue;
        };
        let axes_agree = declared.iter().zip(graph).all(|(d, g)| match (d, g) {
            (Some(d), Some(g)) => d == g,
            _ => true,
        });
        if declared.len() != graph.len() || !axes_agree {
            return Err(Error::ManifestInvalid(format!(
                "manifest `{}` {what} `{}` declares shape {declared:?}, the graph has {graph:?}",
                manifest.id, spec.name
            )));
        }
    }
    Ok(())
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::manifest::{ModelTask, Provenance, TensorSpec};

    fn manifest(inputs: Vec<TensorSpec>) -> ModelManifest {
        ModelManifest {
            id: "m".into(),
            description: None,
            onnx: "m.onnx".into(),
            inputs,
            outputs: Vec::new(),
            task: ModelTask::Generic,
            provenance: Provenance::default(),
            extras: Default::default(),
        }
    }

    fn spec(name: &str, shape: Option<Vec<Option<usize>>>) -> TensorSpec {
        TensorSpec {
            name: name.into(),
            description: None,
            shape,
        }
    }

    #[test]
    fn manifest_names_must_match_the_graph() {
        let graph = ["x".to_string(), "nce".to_string()];
        let ok = manifest(vec![spec("nce", None), spec("x", None)]);
        assert!(check_manifest_io(&ok, Side::Input, &graph, &[]).is_ok());
        let typo = manifest(vec![spec("x", None), spec("nce_", None)]);
        let err = check_manifest_io(&typo, Side::Input, &graph, &[]).unwrap_err();
        assert!(err.to_string().contains("nce_"), "{err}");
        let partial = manifest(vec![spec("x", None)]);
        assert!(check_manifest_io(&partial, Side::Input, &graph, &[]).is_err());
        // Nothing declared: nothing to check.
        assert!(check_manifest_io(&manifest(Vec::new()), Side::Input, &graph, &[]).is_ok());
        // The output side reads the manifest's outputs, which are empty here.
        assert!(check_manifest_io(&typo, Side::Output, &graph, &[]).is_ok());
    }

    #[test]
    fn manifest_shapes_must_match_where_both_are_known() {
        let graph = ["x".to_string()];
        let shapes = [Some(vec![None, Some(4)])];
        let m = |s| manifest(vec![spec("x", s)]);
        assert!(
            check_manifest_io(&m(Some(vec![None, Some(4)])), Side::Input, &graph, &shapes).is_ok()
        );
        assert!(
            check_manifest_io(&m(Some(vec![Some(8), None])), Side::Input, &graph, &shapes).is_ok()
        );
        assert!(check_manifest_io(&m(None), Side::Input, &graph, &shapes).is_ok());
        assert!(
            check_manifest_io(&m(Some(vec![None, Some(5)])), Side::Input, &graph, &shapes).is_err()
        );
        assert!(check_manifest_io(&m(Some(vec![None])), Side::Input, &graph, &shapes).is_err());
        // Unknown graph rank: names are still checked, shapes are not.
        assert!(check_manifest_io(&m(Some(vec![None])), Side::Input, &graph, &[None]).is_ok());
    }
}
