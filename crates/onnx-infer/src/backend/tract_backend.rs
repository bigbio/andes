//! [`tract`](https://github.com/sonos/tract) ONNX backend (pure Rust, Apache-2.0).

use std::path::Path;
use std::sync::Arc;

use tract_onnx::prelude::*;

use crate::error::{BackendKind, Error, Result};
use crate::session::InferenceBackend;
use crate::tensor::NamedTensor;

type TractRunnable = SimplePlan<TypedFact, Box<dyn TypedOp>, Graph<TypedFact, Box<dyn TypedOp>>>;

/// Tract-backed session.
pub struct TractSession {
    model: Arc<TractRunnable>,
    input_names: Vec<String>,
    output_names: Vec<String>,
}

impl TractSession {
    /// Load and optimize an ONNX file.
    pub fn load(path: &Path) -> Result<Self> {
        let inferred = tract_onnx::onnx()
            .model_for_path(path)
            .map_err(|e| Error::Inference(format!("tract load {}: {e}", path.display())))?;

        // Capture ONNX-facing names before optimization rewrites the graph.
        let input_names: Vec<String> = inferred
            .inputs
            .iter()
            .enumerate()
            .map(|(i, outlet)| name_or_fallback(&inferred.nodes[outlet.node].name, "input", i))
            .collect();
        let output_names: Vec<String> = inferred
            .outputs
            .iter()
            .enumerate()
            .map(|(i, outlet)| name_or_fallback(&inferred.nodes[outlet.node].name, "output", i))
            .collect();

        let runnable = inferred
            .into_optimized()
            .map_err(|e| Error::Inference(format!("tract optimize: {e}")))?
            .into_runnable()
            .map_err(|e| Error::Inference(format!("tract runnable: {e}")))?;

        Ok(Self {
            model: Arc::new(runnable),
            input_names,
            output_names,
        })
    }

    fn resolve_inputs<'a>(&self, inputs: &'a [NamedTensor]) -> Result<Vec<&'a NamedTensor>> {
        // Single-input convenience: accept any name.
        if self.input_names.len() == 1 && inputs.len() == 1 {
            return Ok(vec![&inputs[0]]);
        }
        if let Some(unknown) = inputs
            .iter()
            .find(|t| !self.input_names.iter().any(|n| n == &t.name))
        {
            return Err(Error::UnknownInput {
                name: unknown.name.clone(),
                available: self.input_names.clone(),
            });
        }
        let mut ordered = Vec::with_capacity(self.input_names.len());
        for required in &self.input_names {
            let nt =
                inputs
                    .iter()
                    .find(|t| &t.name == required)
                    .ok_or_else(|| Error::MissingInput {
                        name: required.clone(),
                    })?;
            ordered.push(nt);
        }
        Ok(ordered)
    }
}

fn name_or_fallback(name: &str, kind: &str, i: usize) -> String {
    if name.is_empty() {
        format!("{kind}_{i}")
    } else {
        name.to_string()
    }
}

impl InferenceBackend for TractSession {
    fn kind(&self) -> BackendKind {
        BackendKind::Tract
    }

    fn share(&self) -> Option<Box<dyn InferenceBackend>> {
        Some(Box::new(Self {
            model: Arc::clone(&self.model),
            input_names: self.input_names.clone(),
            output_names: self.output_names.clone(),
        }))
    }

    fn input_names(&self) -> &[String] {
        &self.input_names
    }

    fn output_names(&self) -> &[String] {
        &self.output_names
    }

    fn run(&mut self, inputs: &[NamedTensor]) -> Result<Vec<NamedTensor>> {
        let ordered = self.resolve_inputs(inputs)?;
        let mut tensors = TVec::new();
        for (i, nt) in ordered.into_iter().enumerate() {
            let arr = tract_ndarray::ArrayD::from_shape_vec(nt.shape.clone(), nt.data.clone())
                .map_err(|e| Error::Inference(format!("ndarray shape for `{}`: {e}", nt.name)))?;
            let tensor: Tensor = arr.into();
            // NamedTensor carries f32; graphs with integer inputs (token / index
            // tensors such as PeptDeep `aa_indices`) need a cast to the declared type.
            let dt = self
                .model
                .model()
                .input_fact(i)
                .map_err(|e| Error::Inference(format!("input fact for `{}`: {e}", nt.name)))?
                .datum_type;
            let tensor = if dt == f32::datum_type() {
                tensor
            } else {
                tensor
                    .cast_to_dt(dt)
                    .map_err(|e| Error::Inference(format!("cast `{}` to {dt:?}: {e}", nt.name)))?
                    .into_owned()
            };
            tensors.push(tensor.into());
        }

        let outputs = self
            .model
            .run(tensors)
            .map_err(|e| Error::Inference(format!("tract run: {e}")))?;

        let mut named = Vec::with_capacity(outputs.len());
        for (i, t) in outputs.into_iter().enumerate() {
            let name = self
                .output_names
                .get(i)
                .cloned()
                .unwrap_or_else(|| format!("output_{i}"));
            // Integer outputs (e.g. TreeEnsembleClassifier labels) come back as f32.
            let t = if t.datum_type() == f32::datum_type() {
                t.into_tensor()
            } else {
                t.cast_to::<f32>()
                    .map_err(|e| Error::Inference(format!("cast output `{name}` to f32: {e}")))?
                    .into_owned()
            };
            let view = t
                .to_array_view::<f32>()
                .map_err(|e| Error::Inference(format!("output `{name}` is not f32: {e}")))?;
            let shape: Vec<usize> = view.shape().to_vec();
            let data = view.iter().copied().collect::<Vec<_>>();
            named.push(NamedTensor::new(name, shape, data)?);
        }
        Ok(named)
    }
}
