//! [ONNX Runtime](https://onnxruntime.ai) backend (MIT) through the
//! [`ort`](https://ort.pyke.io) crate.
//!
//! The runtime is a shared library loaded on first use (`ort`'s `load-dynamic`
//! mode): `ORT_DYLIB_PATH` names it, otherwise the platform loader searches
//! for `libonnxruntime.so` / `libonnxruntime.dylib` / `onnxruntime.dll`.
//! [`init_runtime`] sets the path from code instead. Any ONNX Runtime 1.22 or
//! newer works.
//!
//! One loaded model serves every worker thread: [`InferenceBackend::share`]
//! hands out handles on the same session, and each handle calls the C API
//! `Run` directly, which ONNX Runtime documents as safe to call concurrently
//! on one session.

use std::ffi::{c_char, CString};
use std::path::Path;
use std::ptr::{self, NonNull};
use std::sync::Arc;

use ort::session::builder::GraphOptimizationLevel;
use ort::session::{Session, SharedSessionInner};
use ort::tensor::TensorElementType;
use ort::value::{DynValue, Tensor, ValueType};
use ort::AsPointer;

use crate::error::{BackendKind, Error, Result};
use crate::session::{BackendOptions, InferenceBackend};
use crate::tensor::NamedTensor;

/// Load the ONNX Runtime shared library from `path` instead of
/// `ORT_DYLIB_PATH` / the loader search path. Call before the first model is
/// loaded; later calls have no effect. A file that cannot be loaded as the
/// runtime (not a shared library, an incompatible build, missing
/// dependencies) is an error.
pub fn init_runtime(path: impl AsRef<Path>) -> Result<()> {
    let path = path.as_ref();
    if !path.is_file() {
        return Err(Error::Inference(format!(
            "ONNX Runtime library {} not found",
            path.display()
        )));
    }
    let lib = path.display().to_string();
    catch_load(|| ort::init_from(lib).commit())
        .map_err(|msg| {
            Error::Inference(format!(
                "cannot load the ONNX Runtime library {} ({msg})",
                path.display()
            ))
        })?
        .map(|_| ())
        .map_err(|e| Error::Inference(format!("ONNX Runtime init: {e}")))
}

/// Resolve the runtime library once, turning `ort`'s load panic into an error.
fn ensure_runtime() -> Result<()> {
    catch_load(|| {
        ort::api();
    })
    .map_err(|msg| {
        Error::Inference(format!(
            "cannot load the ONNX Runtime library ({msg}); set ORT_DYLIB_PATH to libonnxruntime"
        ))
    })
}

/// Runs `f`, which may load the runtime library; `ort` panics when the
/// library cannot be loaded, and the panic message is returned as the error.
fn catch_load<T>(f: impl FnOnce() -> T) -> std::result::Result<T, String> {
    std::panic::catch_unwind(std::panic::AssertUnwindSafe(f)).map_err(|p| {
        p.downcast_ref::<String>()
            .cloned()
            .or_else(|| p.downcast_ref::<&str>().map(|s| s.to_string()))
            .unwrap_or_default()
    })
}

fn ort_err(what: &str) -> impl Fn(ort::Error) -> Error + '_ {
    move |e| Error::Inference(format!("ort {what}: {e}"))
}

/// Axis sizes as ONNX Runtime reports them: a negative size is a dynamic axis.
fn axes(shape: &[i64]) -> Vec<Option<usize>> {
    shape.iter().map(|&d| usize::try_from(d).ok()).collect()
}

/// One loaded graph, shared by every handle.
struct Model {
    inner: Arc<SharedSessionInner>,
    input_names: Vec<String>,
    input_types: Vec<TensorElementType>,
    input_shapes: Vec<Option<Vec<Option<usize>>>>,
    output_names: Vec<String>,
    output_shapes: Vec<Option<Vec<Option<usize>>>>,
    c_inputs: Vec<CString>,
    c_outputs: Vec<CString>,
}

/// ONNX Runtime session handle.
pub struct OrtSession {
    model: Arc<Model>,
}

impl OrtSession {
    /// Load an ONNX file. `opts.intra_threads` sets the operator thread pool
    /// of this session (1 = run on the calling thread); inter-op is always 1.
    pub fn load(path: &Path, opts: &BackendOptions) -> Result<Self> {
        ensure_runtime()?;
        let intra = opts.intra_threads.max(1);
        let session = Session::builder()
            .map_err(ort_err("session builder"))?
            .with_optimization_level(GraphOptimizationLevel::Level3)
            .map_err(ort_err("optimization level"))?
            .with_intra_threads(intra)
            .map_err(ort_err("intra threads"))?
            .with_inter_threads(1)
            .map_err(ort_err("inter threads"))?
            .with_intra_op_spinning(false)
            .map_err(ort_err("intra-op spinning"))?
            .commit_from_file(path)
            .map_err(|e| Error::Inference(format!("ort load {}: {e}", path.display())))?;

        let mut input_names = Vec::new();
        let mut input_types = Vec::new();
        let mut input_shapes = Vec::new();
        for input in &session.inputs {
            match &input.input_type {
                ValueType::Tensor { ty, shape, .. } => {
                    input_names.push(input.name.clone());
                    input_types.push(*ty);
                    input_shapes.push(Some(axes(shape)));
                }
                other => {
                    return Err(Error::Inference(format!(
                        "ort: input `{}` is {other:?}, only tensor inputs are supported",
                        input.name
                    )));
                }
            }
        }
        // Non-tensor outputs (ZipMap sequences of maps) are not requested.
        let (output_names, output_shapes): (Vec<String>, Vec<_>) = session
            .outputs
            .iter()
            .filter_map(|o| match &o.output_type {
                ValueType::Tensor { shape, .. } => Some((o.name.clone(), Some(axes(shape)))),
                _ => None,
            })
            .unzip();
        if output_names.is_empty() {
            return Err(Error::Inference(format!(
                "ort: {} has no tensor outputs",
                path.display()
            )));
        }
        let cstr = |n: &String| {
            CString::new(n.as_str())
                .map_err(|_| Error::Inference(format!("ort: name `{n}` contains NUL")))
        };
        let c_inputs = input_names.iter().map(cstr).collect::<Result<Vec<_>>>()?;
        let c_outputs = output_names.iter().map(cstr).collect::<Result<Vec<_>>>()?;
        let inner = session.inner();
        Ok(Self {
            model: Arc::new(Model {
                inner,
                input_names,
                input_types,
                input_shapes,
                output_names,
                output_shapes,
                c_inputs,
                c_outputs,
            }),
        })
    }

    fn resolve_inputs<'a>(&self, inputs: &'a [NamedTensor]) -> Result<Vec<&'a NamedTensor>> {
        let names = &self.model.input_names;
        // Single-input convenience: accept any name.
        if names.len() == 1 && inputs.len() == 1 {
            return Ok(vec![&inputs[0]]);
        }
        if let Some(unknown) = inputs.iter().find(|t| !names.iter().any(|n| n == &t.name)) {
            return Err(Error::UnknownInput {
                name: unknown.name.clone(),
                available: names.clone(),
            });
        }
        names
            .iter()
            .map(|required| {
                inputs
                    .iter()
                    .find(|t| &t.name == required)
                    .ok_or_else(|| Error::MissingInput {
                        name: required.clone(),
                    })
            })
            .collect()
    }
}

/// `NamedTensor` data (f32) as a tensor of the graph's declared element type.
fn to_value(nt: &NamedTensor, ty: TensorElementType) -> Result<DynValue> {
    let shape = nt.shape.clone();
    let v: DynValue = match ty {
        TensorElementType::Float32 => Tensor::from_array((shape, nt.data.clone())).map(Into::into),
        TensorElementType::Float64 => {
            Tensor::from_array((shape, nt.data.iter().map(|&x| x as f64).collect::<Vec<_>>()))
                .map(Into::into)
        }
        TensorElementType::Int64 => {
            Tensor::from_array((shape, nt.data.iter().map(|&x| x as i64).collect::<Vec<_>>()))
                .map(Into::into)
        }
        TensorElementType::Int32 => {
            Tensor::from_array((shape, nt.data.iter().map(|&x| x as i32).collect::<Vec<_>>()))
                .map(Into::into)
        }
        TensorElementType::Bool => {
            Tensor::from_array((shape, nt.data.iter().map(|&x| x != 0.0).collect::<Vec<_>>()))
                .map(Into::into)
        }
        other => {
            return Err(Error::Inference(format!(
                "ort: input `{}` has unsupported element type {other:?}",
                nt.name
            )));
        }
    }
    .map_err(|e| Error::Inference(format!("ort tensor `{}`: {e}", nt.name)))?;
    Ok(v)
}

/// A tensor output as f32, whatever its element type.
fn to_named(name: &str, v: &DynValue) -> Result<NamedTensor> {
    let err = |e: ort::Error| Error::Inference(format!("ort output `{name}`: {e}"));
    let ty = match v.dtype() {
        ValueType::Tensor { ty, .. } => *ty,
        other => {
            return Err(Error::Inference(format!(
                "ort output `{name}` is {other:?}, not a tensor"
            )))
        }
    };
    macro_rules! extract {
        ($t:ty) => {{
            let (shape, data) = v.try_extract_tensor::<$t>().map_err(err)?;
            (
                shape.iter().map(|&d| d as usize).collect::<Vec<_>>(),
                data.iter().map(|&x| x as f32).collect(),
            )
        }};
    }
    let (shape, data): (Vec<usize>, Vec<f32>) = match ty {
        TensorElementType::Float32 => {
            let (shape, data) = v.try_extract_tensor::<f32>().map_err(err)?;
            (shape.iter().map(|&d| d as usize).collect(), data.to_vec())
        }
        TensorElementType::Float64 => extract!(f64),
        TensorElementType::Int64 => extract!(i64),
        TensorElementType::Int32 => extract!(i32),
        TensorElementType::Bool => {
            let (shape, data) = v.try_extract_tensor::<bool>().map_err(err)?;
            (
                shape.iter().map(|&d| d as usize).collect(),
                data.iter().map(|&b| f32::from(u8::from(b))).collect(),
            )
        }
        other => {
            return Err(Error::Inference(format!(
                "ort output `{name}` has unsupported element type {other:?}"
            )));
        }
    };
    NamedTensor::new(name, shape, data)
}

impl InferenceBackend for OrtSession {
    fn kind(&self) -> BackendKind {
        BackendKind::Ort
    }

    fn share(&self) -> Option<Box<dyn InferenceBackend>> {
        Some(Box::new(Self {
            model: Arc::clone(&self.model),
        }))
    }

    fn input_names(&self) -> &[String] {
        &self.model.input_names
    }

    fn output_names(&self) -> &[String] {
        &self.model.output_names
    }

    fn input_shapes(&self) -> Vec<Option<Vec<Option<usize>>>> {
        self.model.input_shapes.clone()
    }

    fn output_shapes(&self) -> Vec<Option<Vec<Option<usize>>>> {
        self.model.output_shapes.clone()
    }

    fn run(&mut self, inputs: &[NamedTensor]) -> Result<Vec<NamedTensor>> {
        let m = &*self.model;
        let ordered = self.resolve_inputs(inputs)?;
        let values = ordered
            .iter()
            .zip(&m.input_types)
            .map(|(nt, ty)| to_value(nt, *ty))
            .collect::<Result<Vec<_>>>()?;

        let in_names: Vec<*const c_char> = m.c_inputs.iter().map(|c| c.as_ptr()).collect();
        let out_names: Vec<*const c_char> = m.c_outputs.iter().map(|c| c.as_ptr()).collect();
        let in_ptrs: Vec<*const ort::sys::OrtValue> = values.iter().map(|v| v.ptr()).collect();
        let mut out_ptrs: Vec<*mut ort::sys::OrtValue> = vec![ptr::null_mut(); out_names.len()];
        // SAFETY: the session pointer lives as long as `m.inner`; name and
        // value arrays outlive the call and have the lengths passed; ONNX
        // Runtime allows concurrent `Run` on one session.
        let status = unsafe {
            (ort::api().Run)(
                m.inner.ptr().cast_mut(),
                ptr::null(),
                in_names.as_ptr(),
                in_ptrs.as_ptr(),
                in_ptrs.len(),
                out_names.as_ptr(),
                out_names.len(),
                out_ptrs.as_mut_ptr(),
            )
        };
        // SAFETY: `status` comes straight from the C API.
        unsafe { ort::error::status_to_result(status) }.map_err(ort_err("run"))?;
        // Take ownership of every output before any early return, so none leak.
        let outputs: Vec<Option<DynValue>> = out_ptrs
            .into_iter()
            .map(|p| {
                // SAFETY: a non-null pointer returned by `Run` is an owned
                // OrtValue allocated by this session.
                NonNull::new(p)
                    .map(|p| unsafe { DynValue::from_ptr(p, Some(Arc::clone(&m.inner))) })
            })
            .collect();
        // An output may alias an input buffer (a graph that returns an input
        // as an output), so `values` stays alive until every output is copied
        // and dropped.
        let result = m
            .output_names
            .iter()
            .zip(outputs)
            .map(|(name, v)| {
                let v = v.ok_or_else(|| {
                    Error::Inference(format!("ort: no value for output `{name}`"))
                })?;
                to_named(name, &v)
            })
            .collect();
        drop(values);
        result
    }
}
