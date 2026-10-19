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
    input_shapes: Vec<Option<Vec<Option<usize>>>>,
    output_shapes: Vec<Option<Vec<Option<usize>>>>,
}

/// The shape an inference fact states: `None` for an open (unknown-rank) shape,
/// else one entry per axis, `None` where the axis is not a fixed integer.
fn fact_shape(model: &InferenceModel, outlet: OutletId) -> Option<Vec<Option<usize>>> {
    use tract_onnx::tract_hir::infer::GenericFactoid;
    let fact = model.outlet_fact(outlet).ok()?;
    if fact.shape.is_open() {
        return None;
    }
    Some(
        fact.shape
            .dims()
            .map(|d| match d {
                GenericFactoid::Only(dim) => dim.to_usize().ok(),
                GenericFactoid::Any => None,
            })
            .collect(),
    )
}

impl TractSession {
    /// Load and optimize an ONNX file.
    pub fn load(path: &Path) -> Result<Self> {
        let inferred = tract_onnx::onnx()
            .model_for_path(path)
            .map_err(|e| Error::Inference(format!("tract load {}: {e}", path.display())))?;

        // Capture ONNX-facing names before optimization rewrites the graph.
        // An input is a source node named after its tensor. An output is an
        // outlet of an arbitrary node: its ONNX tensor name is the outlet
        // label, while the node name is the producing operator's name.
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
            .map(|(i, outlet)| {
                name_or_fallback(inferred.outlet_label(*outlet).unwrap_or(""), "output", i)
            })
            .collect();

        let input_shapes = inferred
            .inputs
            .iter()
            .map(|o| fact_shape(&inferred, *o))
            .collect();
        let output_shapes = inferred
            .outputs
            .iter()
            .map(|o| fact_shape(&inferred, *o))
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
            input_shapes,
            output_shapes,
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
            input_shapes: self.input_shapes.clone(),
            output_shapes: self.output_shapes.clone(),
        }))
    }

    fn input_names(&self) -> &[String] {
        &self.input_names
    }

    fn output_names(&self) -> &[String] {
        &self.output_names
    }

    fn input_shapes(&self) -> Vec<Option<Vec<Option<usize>>>> {
        self.input_shapes.clone()
    }

    fn output_shapes(&self) -> Vec<Option<Vec<Option<usize>>>> {
        self.output_shapes.clone()
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
