use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use crate::error::{Error, Result};

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct TensorSpec {
    pub name: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    #[serde(
        default,
        deserialize_with = "deserialize_shape",
        skip_serializing_if = "Option::is_none"
    )]
    pub shape: Option<Vec<Option<usize>>>,
}

fn deserialize_shape<'de, D>(
    deserializer: D,
) -> std::result::Result<Option<Vec<Option<usize>>>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    use serde::de::Error as _;
    use serde_json::Value;

    let axes: Option<Vec<Value>> = Option::deserialize(deserializer)?;
    let Some(axes) = axes else {
        return Ok(None);
    };
    axes.into_iter()
        .map(|axis| match axis {
            Value::Null | Value::String(_) => Ok(None),
            Value::Number(n) => {
                if let Some(size) = n.as_u64() {
                    usize::try_from(size)
                        .map(Some)
                        .map_err(|_| D::Error::custom(format!("shape axis {size} is too large")))
                } else if n.as_i64().is_some_and(|i| i < 0) {
                    Ok(None)
                } else {
                    Err(D::Error::custom(format!(
                        "shape axis {n} is not an integer"
                    )))
                }
            }
            other => Err(D::Error::custom(format!(
                "shape axis {other} must be an integer, null or a symbol name"
            ))),
        })
        .collect::<std::result::Result<Vec<_>, _>>()
        .map(Some)
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Default)]
pub struct Provenance {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub origin: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub license: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub notice: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub notes: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum ModelTask {
    #[default]
    Generic,
    Ms2Intensity,
    RetentionTime,
    Ccs,
    TreeEnsemble,
    DeNovo,
    Rescoring,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct ModelManifest {
    pub id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    pub onnx: PathBuf,
    #[serde(default)]
    pub inputs: Vec<TensorSpec>,
    #[serde(default)]
    pub outputs: Vec<TensorSpec>,
    #[serde(default)]
    pub task: ModelTask,
    #[serde(default)]
    pub provenance: Provenance,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub extras: BTreeMap<String, serde_json::Value>,
}

impl ModelManifest {
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

    fn resolve_onnx_path(&mut self, manifest_path: &Path) -> Result<()> {
        if self.onnx.is_absolute() {
            return Ok(());
        }
        let parent = manifest_path.parent().unwrap_or_else(|| Path::new("."));
        self.onnx = parent.join(&self.onnx);
        Ok(())
    }

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

    #[test]
    fn shape_axes_accept_every_dynamic_spelling() {
        let parse = |shape: &str| -> serde_json::Result<TensorSpec> {
            serde_json::from_str(&format!(r#"{{"name": "x", "shape": {shape}}}"#))
        };
        assert_eq!(
            parse(r#"[null, -1, "seq_length", 8]"#).unwrap().shape,
            Some(vec![None, None, None, Some(8)])
        );
        assert_eq!(parse("null").unwrap().shape, None);
        let no_shape: TensorSpec = serde_json::from_str(r#"{"name": "x"}"#).unwrap();
        assert_eq!(no_shape.shape, None);
        assert!(parse("[1.5]").is_err());
        assert!(parse("[true]").is_err());
        let spec = parse(r#"[-1, 8]"#).unwrap();
        assert_eq!(
            serde_json::to_string(&spec).unwrap(),
            r#"{"name":"x","shape":[null,8]}"#
        );
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
