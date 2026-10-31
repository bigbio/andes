//! Shared by the `rescore_pin` and `model_parity` examples: the library
//! predictors ([`onnx_infer::predictors`]) plus the model families only the
//! examples run (MS²PIP, andes-ms2 and the Koina-exported graphs).
//!
//! Peptidoforms a model cannot encode are rejected by
//! [`Predictor::supports_rt`] / [`Predictor::supports_ms2`] and become missing
//! features downstream, never a silently wrong tensor.

// Each example uses a different subset.
#![allow(dead_code, unused_imports)]

use std::path::Path;
use std::sync::OnceLock;

use onnx_infer::ms2_native::{Ms2Conditions, Ms2Model, Ms2Scratch, N_BASE_ION, N_ION};
use onnx_infer::prosit::{self, PrositPeptide, ResidueState};
use onnx_infer::rt_native::RtPeptide;
use onnx_infer::{ms2pip, BackendKind, BackendOptions, OnnxSession, Result};

pub use onnx_infer::predictors::form::{parse_peptide, Form, Mass, Mod};
pub use onnx_infer::predictors::spectral::Ion;
pub use onnx_infer::predictors::{DeepLc, Engine, PeptDeep, Predictor};

static BACKEND: OnceLock<(BackendKind, BackendOptions)> = OnceLock::new();

/// Parse a `--backend` value.
pub fn parse_backend(s: &str) -> BackendKind {
    match s {
        "tract" => BackendKind::Tract,
        "ort" | "onnxruntime" => BackendKind::Ort,
        other => panic!("unknown --backend {other} (tract | ort)"),
    }
}

/// Engine for every graph [`load`]ed afterwards; call once, before loading.
pub fn set_backend(kind: BackendKind, opts: BackendOptions) {
    BACKEND.set((kind, opts)).expect("set_backend called twice");
}

/// The engine set by [`set_backend`]; tract with default options otherwise.
pub fn backend() -> &'static (BackendKind, BackendOptions) {
    BACKEND.get_or_init(|| (BackendKind::Tract, BackendOptions::default()))
}

/// The engine set by [`set_backend`], for the library predictors.
pub fn engine() -> Engine {
    let (kind, opts) = backend();
    Engine {
        kind: *kind,
        opts: opts.clone(),
    }
}

pub fn load(models: &Path, manifest: &str) -> OnnxSession {
    engine()
        .load(models, manifest)
        .unwrap_or_else(|e| panic!("load {manifest}: {e}"))
}

// ---------------------------------------------------------------- MS²PIP

pub struct Ms2pip {
    pub dir: std::path::PathBuf,
    /// ms2pip model name (`HCD2021`, `CID`, `TMT`, …).
    pub model: String,
}

const MS2PIP_IONS: [Ion; 2] = [(true, 1), (false, 1)];

impl Predictor for Ms2pip {
    fn prefix(&self) -> &str {
        "ms2pip"
    }
    fn ions(&self) -> &[Ion] {
        &MS2PIP_IONS
    }
    fn has_rt(&self) -> bool {
        false
    }
    fn supports_rt(&self, _: &Form) -> bool {
        false
    }
    fn supports_ms2(&self, f: &Form, z: u8) -> bool {
        z > 0 && ms2pip::compute_features(&f.seq, z).is_ok()
    }
    fn batch_key(&self, _: &Form) -> usize {
        0
    }
    fn sessions(&self) -> Result<Vec<OnnxSession>> {
        Ok(["b", "y"]
            .iter()
            .map(|ion| load(&self.dir, &format!("ms2pip_{}_{ion}.json", self.model)))
            .collect())
    }
    fn rt(&self, _: &mut [OnnxSession], _: &[&Form]) -> Result<Vec<f32>> {
        unreachable!("MS²PIP has no RT model")
    }
    fn ms2(&self, s: &mut [OnnxSession], items: &[(&Form, u8)]) -> Result<Vec<Vec<f32>>> {
        let peps: Vec<(&str, u8)> = items.iter().map(|(f, z)| (f.seq.as_str(), *z)).collect();
        let batch = ms2pip::pack_batch(&peps).unwrap();
        let input = std::slice::from_ref(&batch.features);
        let b = ms2pip::log2_predictions(&s[0].run(input).expect("ms2pip b")[0]).unwrap();
        let y = ms2pip::log2_predictions(&s[1].run(input).expect("ms2pip y")[0]).unwrap();
        Ok(batch
            .offsets
            .windows(2)
            .map(|w| {
                (w[0]..w[1])
                    .flat_map(|r| [ms2pip::to_linear(b[r]), ms2pip::to_linear(y[r])])
                    .collect()
            })
            .collect())
    }
}

// ---------------------------------------------------------------- andes-ms2

/// b z1, b z2, y z1, y z2: the base ions andes-ms2 predicts per site.
const ANDES_MS2_IONS: [Ion; 4] = [(true, 1), (true, 2), (false, 1), (false, 2)];

/// The native andes-ms2 model (`onnx_infer::ms2_native`); no ONNX graphs.
pub struct AndesMs2 {
    pub model: Ms2Model,
    pub nce: Option<f32>,
    pub activation: String,
    pub instrument: String,
}

impl AndesMs2 {
    /// `path` is the `AMS20001` weight file, or a directory holding `andes_ms2.ams2`.
    pub fn new(path: &Path, nce: Option<f32>, activation: &str, instrument: &str) -> Self {
        let file = if path.is_dir() {
            path.join("andes_ms2.ams2")
        } else {
            path.to_path_buf()
        };
        let model =
            Ms2Model::load(&file).unwrap_or_else(|e| panic!("load {}: {e}", file.display()));
        Self {
            model,
            nce,
            activation: activation.into(),
            instrument: instrument.into(),
        }
    }

    /// The [`RtPeptide`] the form was parsed from (see [`parse_peptide`]).
    pub fn peptide(f: &Form) -> RtPeptide {
        f.rt_peptide()
    }

    fn conditions(&self, z: u8) -> Ms2Conditions {
        Ms2Conditions::new(z, self.nce, &self.activation, &self.instrument)
    }
}

impl Predictor for AndesMs2 {
    fn prefix(&self) -> &str {
        "ams2"
    }
    fn ions(&self) -> &[Ion] {
        &ANDES_MS2_IONS
    }
    fn has_rt(&self) -> bool {
        false
    }
    fn supports_rt(&self, _: &Form) -> bool {
        false
    }
    fn supports_ms2(&self, f: &Form, z: u8) -> bool {
        z > 0 && f.seq.len() >= 2 && Self::peptide(f).encode().is_some()
    }
    fn batch_key(&self, _: &Form) -> usize {
        0
    }
    fn sessions(&self) -> Result<Vec<OnnxSession>> {
        Ok(Vec::new())
    }
    fn rt(&self, _: &mut [OnnxSession], _: &[&Form]) -> Result<Vec<f32>> {
        unreachable!("andes-ms2 has no RT model")
    }
    fn ms2(&self, _: &mut [OnnxSession], items: &[(&Form, u8)]) -> Result<Vec<Vec<f32>>> {
        let mut scratch = Ms2Scratch::default();
        Ok(items
            .iter()
            .map(|(f, z)| {
                let full = self
                    .model
                    .predict_with(&Self::peptide(f), &self.conditions(*z), &mut scratch)
                    .expect("supports_ms2 checked the encoding");
                full.chunks_exact(N_ION)
                    .flat_map(|site| site[..N_BASE_ION].to_vec())
                    .collect()
            })
            .collect())
    }
}

// ---------------------------------------------------------------- Prosit

pub struct Prosit {
    pub dir: std::path::PathBuf,
    pub nce: f32,
}

const PROSIT_IONS: [Ion; 6] = [
    (true, 1),
    (true, 2),
    (true, 3),
    (false, 1),
    (false, 2),
    (false, 3),
];

impl Prosit {
    /// `None` when the peptidoform carries anything but CAM-C / oxidised M.
    pub fn peptide(f: &Form) -> Option<PrositPeptide<'_>> {
        let mut states = vec![ResidueState::Plain; f.seq.len()];
        for (site, m) in &f.mods {
            states[site.checked_sub(1)?] = match m {
                Mod::Carbamidomethyl => ResidueState::Carbamidomethyl,
                Mod::Oxidation => ResidueState::Oxidation,
                _ => return None,
            };
        }
        let p = PrositPeptide {
            sequence: &f.seq,
            states,
        };
        prosit::encode_sequence(p.sequence, &p.states).ok()?;
        Some(p)
    }
}

impl Predictor for Prosit {
    fn prefix(&self) -> &str {
        "prosit"
    }
    fn ions(&self) -> &[Ion] {
        &PROSIT_IONS
    }
    fn has_rt(&self) -> bool {
        true
    }
    fn supports_rt(&self, f: &Form) -> bool {
        Self::peptide(f).is_some()
    }
    fn supports_ms2(&self, f: &Form, z: u8) -> bool {
        (1..=prosit::MAX_CHARGE).contains(&z) && f.seq.len() >= 2 && Self::peptide(f).is_some()
    }
    fn batch_key(&self, _: &Form) -> usize {
        0
    }
    fn sessions(&self) -> Result<Vec<OnnxSession>> {
        Ok(vec![
            load(&self.dir, "prosit_irt.json"),
            load(&self.dir, "prosit_intensity_hcd.json"),
        ])
    }
    fn rt(&self, s: &mut [OnnxSession], forms: &[&Form]) -> Result<Vec<f32>> {
        let peps: Vec<PrositPeptide<'_>> =
            forms.iter().map(|f| Self::peptide(f).unwrap()).collect();
        let input = prosit::pack_irt_batch(&peps).unwrap();
        Ok(
            s[0].run(std::slice::from_ref(&input)).expect("prosit irt")[0]
                .data
                .clone(),
        )
    }
    fn ms2(&self, s: &mut [OnnxSession], items: &[(&Form, u8)]) -> Result<Vec<Vec<f32>>> {
        let peps: Vec<PrositPeptide<'_>> = items
            .iter()
            .map(|(f, _)| Self::peptide(f).unwrap())
            .collect();
        let charges: Vec<u8> = items.iter().map(|(_, z)| *z).collect();
        let inputs = prosit::pack_intensity_batch(&peps, &charges, self.nce).unwrap();
        let out = &s[1].run(&inputs).expect("prosit intensity")[0];
        Ok(items
            .iter()
            .enumerate()
            .map(|(b, (f, z))| {
                let raw = &out.data[b * prosit::N_OUT..(b + 1) * prosit::N_OUT];
                prosit::decode_intensity(raw, f.seq.len(), *z).unwrap()
            })
            .collect())
    }
}
