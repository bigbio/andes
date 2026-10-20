//! Sidecar JSON describing an ONNX model for tooling (provenance, I/O names, extras).

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use crate::error::{Error, Result};

/// Declares one model input or output tensor.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct TensorSpec {
    /// ONNX graph input/output name.
    pub name: String,
    /// Optional human description (e.g. "peptide embedding", "NCE").
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    /// Optional static shape hint. Use `null` / omit for fully dynamic axes.
    /// A `null` entry inside the vec means a dynamic dimension (usually batch).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub shape: Option<Vec<Option<usize>>>,
}

/// Where the weights came from (licensing / independence paperwork).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Default)]
pub struct Provenance {
    /// Short origin label (e.g. "AlphaPeptDeep 1.2.0", "own-train 2026-10").
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub origin: Option<String>,
    /// SPDX license of the *weights* (may differ from this crate).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub license: Option<String>,
    /// URL or path to NOTICE / training corpus statement.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub notice: Option<String>,
    /// Free-form notes (training accessions, exclusions, …).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub notes: Option<String>,
}

/// Task hint so callers can pick the right post-processing.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum ModelTask {
    /// Generic / unknown.
    #[default]
    Generic,
    /// Fragment-intensity predictor (e.g. AlphaPeptDeep MS2).
    Ms2Intensity,
    /// Retention-time predictor.
    RetentionTime,
    /// Ion-mobility / CCS predictor.
    Ccs,
    /// Gradient-boosted tree ensemble (andes AGBD / ONNX TreeEnsemble).
    TreeEnsemble,
    /// De novo sequence model (logits / CTC / token probs).
    DeNovo,
    /// Other rescoring feature model.
    Rescoring,
}

/// Sidecar manifest living next to (or pointing at) an `.onnx` file.
///
/// Designed so search engines, de novo tools, and rescoring pipelines share one
/// loader: tool-specific knobs go in [`ModelManifest::extras`].
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct ModelManifest {
    /// Stable id for logs and PIN provenance (e.g. `alphapept-ms2-v1`).
    pub id: String,
    /// Optional human description.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    /// Path to the `.onnx` file. Relative paths are resolved against the
    /// directory containing the manifest JSON.
    pub onnx: PathBuf,
    /// Declared inputs (order is the default feed order when names match).
    #[serde(default)]
    pub inputs: Vec<TensorSpec>,
    /// Declared outputs.
    #[serde(default)]
    pub outputs: Vec<TensorSpec>,
    /// What the model is for.
    #[serde(default)]
    pub task: ModelTask,
    /// Licensing / training provenance for the weights.
    #[serde(default)]
    pub provenance: Provenance,
    /// Escape hatch for NCE calibration keys, RT units, vocab paths, etc.
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub extras: BTreeMap<String, serde_json::Value>,
}

impl ModelManifest {
    /// Read a manifest from JSON on disk.
    pub fn from_path(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let text = fs::read_to_string(path).map_err(|source| Error::Io {
            path: path.to_path_buf(),
            source,
        })?;
        let mut m: Self = serde_json::from_str(&text).map_err(|source| Error::ManifestParse {
            path: path.to_path_buf(),
            source,
        })?;
        m.resolve_onnx_path(path)?;
        m.validate()?;
        Ok(m)
    }

    /// Resolve a relative `onnx` path against the manifest file's parent.
    fn resolve_onnx_path(&mut self, manifest_path: &Path) -> Result<()> {
        if self.onnx.is_absolute() {
            return Ok(());
        }
        let parent = manifest_path.parent().unwrap_or_else(|| Path::new("."));
        self.onnx = parent.join(&self.onnx);
        Ok(())
    }

    /// Basic structural checks (does not open the ONNX file).
    pub fn validate(&self) -> Result<()> {
        if self.id.trim().is_empty() {
            return Err(Error::ManifestInvalid("id must be non-empty".into()));
        }
        if self.onnx.as_os_str().is_empty() {
            return Err(Error::ManifestInvalid("onnx path must be non-empty".into()));
        }
        if let Some(spec) = self
            .inputs
            .iter()
            .chain(self.outputs.iter())
            .find(|spec| spec.name.trim().is_empty())
        {
            return Err(Error::ManifestInvalid(format!(
                "tensor spec name must be non-empty (got {:?})",
                spec.name
            )));
        }
        // Names must be unique within the inputs and within the outputs; ONNX
        // allows one name on both sides (a graph input passed straight through).
        let mut in_names = std::collections::HashSet::new();
        for spec in &self.inputs {
            if !in_names.insert(&spec.name) {
                return Err(Error::ManifestInvalid(format!(
                    "duplicate input name `{}`",
                    spec.name
                )));
            }
        }
        let mut out_names = std::collections::HashSet::new();
        for spec in &self.outputs {
            if !out_names.insert(&spec.name) {
                return Err(Error::ManifestInvalid(format!(
                    "duplicate output name `{}`",
                    spec.name
                )));
            }
        }
        Ok(())
    }

    /// Absolute/resolved path to the ONNX weights.
    pub fn onnx_path(&self) -> &Path {
        &self.onnx
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn loads_and_resolves_relative_onnx() {
        let dir = tempfile::tempdir().unwrap();
        let onnx = dir.path().join("m.onnx");
        fs::write(&onnx, b"not-a-real-onnx").unwrap();
        let man_path = dir.path().join("m.json");
        let mut f = fs::File::create(&man_path).unwrap();
        write!(
            f,
            r#"{{
              "id": "demo",
              "onnx": "m.onnx",
              "inputs": [{{"name": "x"}}],
              "outputs": [{{"name": "y"}}],
              "task": "ms2_intensity",
              "provenance": {{"origin": "test", "license": "Apache-2.0"}}
            }}"#
        )
        .unwrap();
        let m = ModelManifest::from_path(&man_path).unwrap();
        assert_eq!(m.id, "demo");
        assert_eq!(m.task, ModelTask::Ms2Intensity);
        assert_eq!(m.onnx_path(), onnx.as_path());
    }

    fn manifest(inputs: &[&str], outputs: &[&str]) -> ModelManifest {
        let spec = |n: &&str| TensorSpec {
            name: (*n).to_string(),
            description: None,
            shape: None,
        };
        ModelManifest {
            id: "m".into(),
            description: None,
            onnx: PathBuf::from("m.onnx"),
            inputs: inputs.iter().map(spec).collect(),
            outputs: outputs.iter().map(spec).collect(),
            task: ModelTask::Generic,
            provenance: Provenance::default(),
            extras: BTreeMap::new(),
        }
    }

    #[test]
    fn names_are_unique_per_side_only() {
        assert!(manifest(&["x"], &["x"]).validate().is_ok());
        assert!(manifest(&["x", "x"], &["y"]).validate().is_err());
        assert!(manifest(&["x"], &["y", "y"]).validate().is_err());
        assert!(manifest(&[" "], &["y"]).validate().is_err());
    }
}
