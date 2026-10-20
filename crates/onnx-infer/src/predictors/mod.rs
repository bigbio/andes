//! Retention-time and fragment-intensity predictors over andes peptidoforms.
//!
//! A [`Predictor`] turns parsed peptidoforms ([`form::Form`]) into retention
//! time and/or fragment-intensity predictions. Peptidoforms a model cannot
//! encode are rejected by [`Predictor::supports_rt`] /
//! [`Predictor::supports_ms2`] and become missing values downstream, never a
//! silently wrong tensor.
//!
//! - [`NativeRt`]: the andes-rt model ([`crate::rt_native`]), plain Rust.
//! - [`PeptDeep`]: AlphaPeptDeep RT and MS2 graphs exported to ONNX.
//! - [`DeepLc`]: the DeepLC RT graph exported to ONNX.
//!
//! [`spectral::ms2_features`] compares a prediction with an observed
//! spectrum; [`calibration::RtCalibration`] maps predicted RT to minutes.

pub mod calibration;
pub mod form;
pub mod spectral;

use std::collections::HashMap;
use std::fs::File;
use std::path::{Path, PathBuf};

use serde::Deserialize;

use crate::deeplc::{self, DeepLcPeptide};
use crate::error::{Error, Result};
use crate::peptdeep::{pack_ms2_batch, pack_rt_batch, PeptDeepPeptide};
use crate::rt_native::{RtModel, Scratch};
use crate::{BackendKind, BackendOptions, OnnxSession};

use form::{Form, Mod};
use spectral::Ion;

/// The engine ONNX graphs are loaded with.
#[derive(Clone, Debug)]
pub struct Engine {
    /// tract or ONNX Runtime.
    pub kind: BackendKind,
    /// Engine options.
    pub opts: BackendOptions,
}

impl Engine {
    /// Load the graph a sidecar manifest in `dir` names. The error names the
    /// manifest path.
    pub fn load(&self, dir: &Path, manifest: &str) -> Result<OnnxSession> {
        let path = dir.join(manifest);
        if !path.is_file() {
            return Err(Error::InvalidArg(format!(
                "model manifest {} not found",
                path.display()
            )));
        }
        OnnxSession::from_manifest_with(&path, self.kind, &self.opts)
            .map_err(|e| Error::InvalidArg(format!("loading {}: {e}", path.display())))
    }
}

/// A model family that predicts RT and/or fragment intensities.
pub trait Predictor: Sync {
    /// Feature-column prefix (`pdeep`, `deeplc`, `andesrt`).
    fn prefix(&self) -> &str;
    /// Per-site ion layout of [`Self::ms2`] output; empty when no MS2 model.
    fn ions(&self) -> &[Ion];
    /// Whether the family has an RT model.
    fn has_rt(&self) -> bool;
    /// Whether the RT model can encode `form`.
    fn supports_rt(&self, form: &Form) -> bool;
    /// Whether the MS2 model can encode `form` at `charge`.
    fn supports_ms2(&self, form: &Form, charge: u8) -> bool;
    /// Forms in one RT / MS2 batch must share this key.
    fn batch_key(&self, form: &Form) -> usize;
    /// Load the model graphs once; worker threads get [`OnnxSession::share`]d
    /// handles.
    fn sessions(&self) -> Result<Vec<OnnxSession>>;
    /// RT prediction (model units) per form.
    fn rt(&self, s: &mut [OnnxSession], forms: &[&Form]) -> Result<Vec<f32>>;
    /// Per form: (len − 1) sites × `ions().len()` linear intensities, site `i`
    /// = cleavage after residue `i + 1`; ions above the precursor charge are 0.
    fn ms2(&self, s: &mut [OnnxSession], items: &[(&Form, u8)]) -> Result<Vec<Vec<f32>>>;
}

// ---------------------------------------------------------------- andes-rt

/// The andes-rt retention-time model; no ONNX graphs.
pub struct NativeRt {
    model: RtModel,
}

impl NativeRt {
    /// The model from an `ARTN0001` weight buffer.
    pub fn from_bytes(buf: &[u8]) -> Result<Self> {
        let model = RtModel::from_bytes(buf)
            .map_err(|e| Error::InvalidArg(format!("andes-rt weights: {e}")))?;
        Ok(Self { model })
    }

    /// The model from an `ARTN0001` weight file.
    pub fn load(path: &Path) -> Result<Self> {
        let model = RtModel::load(path)
            .map_err(|e| Error::InvalidArg(format!("andes-rt weights {}: {e}", path.display())))?;
        Ok(Self { model })
    }
}

impl Predictor for NativeRt {
    fn prefix(&self) -> &str {
        "andesrt"
    }
    fn ions(&self) -> &[Ion] {
        &[]
    }
    fn has_rt(&self) -> bool {
        true
    }
    fn supports_rt(&self, f: &Form) -> bool {
        f.rt_peptide().encode().is_some()
    }
    fn supports_ms2(&self, _: &Form, _: u8) -> bool {
        false
    }
    fn batch_key(&self, _: &Form) -> usize {
        0
    }
    fn sessions(&self) -> Result<Vec<OnnxSession>> {
        Ok(Vec::new())
    }
    fn rt(&self, _: &mut [OnnxSession], forms: &[&Form]) -> Result<Vec<f32>> {
        let mut scratch = Scratch::default();
        Ok(forms
            .iter()
            .map(|f| {
                self.model
                    .predict_with(&f.rt_peptide(), &mut scratch)
                    .unwrap_or(f32::NAN)
            })
            .collect())
    }
    fn ms2(&self, _: &mut [OnnxSession], _: &[(&Form, u8)]) -> Result<Vec<Vec<f32>>> {
        Err(Error::InvalidArg("andes-rt has no MS2 model".into()))
    }
}

// ---------------------------------------------------------------- PeptDeep

/// `peptdeep_mods.json`: the alphabase modification vectors and the
/// instrument table the exported graphs were built with.
#[derive(Deserialize)]
pub struct PeptDeepMeta {
    /// Length of one modification vector.
    pub mod_dim: usize,
    /// alphabase modification name -> feature vector.
    pub mods: HashMap<String, Vec<f32>>,
    /// MS2 output columns.
    pub ms2_frag_types: Vec<String>,
    /// Upper-case instrument name -> index.
    pub instrument_index: HashMap<String, i64>,
    /// Index of an instrument outside `instrument_index`.
    pub unknown_instrument_index: i64,
}

/// Graph manifests and the metadata file of an exported PeptDeep model.
pub const PEPTDEEP_FILES: [&str; 3] = [
    "peptdeep_rt_dynamic.json",
    "peptdeep_ms2_dynamic.json",
    "peptdeep_mods.json",
];

/// AlphaPeptDeep RT and MS2 graphs.
pub struct PeptDeep {
    dir: PathBuf,
    meta: PeptDeepMeta,
    nce: f32,
    instrument: i64,
    /// Columns of b_z1, b_z2, y_z1, y_z2 in the MS2 output.
    cols: [usize; 4],
    engine: Engine,
    load_rt: bool,
    load_ms2: bool,
}

const PDEEP_IONS: [Ion; 4] = [(true, 1), (true, 2), (false, 1), (false, 2)];

impl PeptDeep {
    /// Read `peptdeep_mods.json` from `dir`. `instrument` is matched
    /// case-insensitively against the model's instrument table (an unknown name
    /// takes the model's unknown-instrument index). `load_rt` / `load_ms2`
    /// select which graphs [`Predictor::sessions`] loads; a missing requested
    /// graph manifest is an error here.
    pub fn new(
        dir: &Path,
        nce: f32,
        instrument: &str,
        engine: Engine,
        load_rt: bool,
        load_ms2: bool,
    ) -> Result<Self> {
        let need = [
            (load_rt, PEPTDEEP_FILES[0]),
            (load_ms2, PEPTDEEP_FILES[1]),
            (true, PEPTDEEP_FILES[2]),
        ];
        for (want, name) in need {
            if want && !dir.join(name).is_file() {
                return Err(Error::InvalidArg(format!(
                    "PeptDeep model file {} not found",
                    dir.join(name).display()
                )));
            }
        }
        let meta_path = dir.join(PEPTDEEP_FILES[2]);
        let file = File::open(&meta_path)
            .map_err(|e| Error::InvalidArg(format!("{}: {e}", meta_path.display())))?;
        let meta: PeptDeepMeta = serde_json::from_reader(std::io::BufReader::new(file))
            .map_err(|e| Error::InvalidArg(format!("{}: {e}", meta_path.display())))?;
        let col = |n: &str| {
            meta.ms2_frag_types
                .iter()
                .position(|f| f == n)
                .ok_or_else(|| {
                    Error::InvalidArg(format!(
                        "{}: ms2_frag_types has no {n}",
                        meta_path.display()
                    ))
                })
        };
        let cols = [col("b_z1")?, col("b_z2")?, col("y_z1")?, col("y_z2")?];
        let instrument = *meta
            .instrument_index
            .get(&instrument.to_uppercase())
            .unwrap_or(&meta.unknown_instrument_index);
        Ok(Self {
            dir: dir.to_path_buf(),
            meta,
            nce,
            instrument,
            cols,
            engine,
            load_rt,
            load_ms2,
        })
    }

    /// alphabase modification names the model's metadata carries.
    pub fn known_mods(&self) -> impl Iterator<Item = &str> {
        self.meta.mods.keys().map(String::as_str)
    }

    /// The MS2 graph's output columns, by fragment type.
    pub fn ms2_frag_types(&self) -> &[String] {
        &self.meta.ms2_frag_types
    }

    /// Output columns of b_z1, b_z2, y_z1, y_z2 (the order of [`Predictor::ms2`]).
    pub fn ms2_columns(&self) -> [usize; 4] {
        self.cols
    }

    fn supported(&self, f: &Form) -> bool {
        f.mods
            .iter()
            .all(|(s, m)| self.meta.mods.contains_key(&m.alphabase(*s, f.residue(*s))))
    }

    fn peptides<'f>(
        &'f self,
        forms: &[&'f Form],
        names: &'f [Vec<String>],
    ) -> Vec<PeptDeepPeptide<'f>> {
        forms
            .iter()
            .zip(names)
            .map(|(f, n)| PeptDeepPeptide {
                sequence: &f.seq,
                mods: f
                    .mods
                    .iter()
                    .zip(n)
                    .map(|((s, _), name)| (*s, self.meta.mods[name].as_slice()))
                    .collect(),
            })
            .collect()
    }

    fn names(forms: &[&Form]) -> Vec<Vec<String>> {
        forms
            .iter()
            .map(|f| {
                f.mods
                    .iter()
                    .map(|(s, m)| m.alphabase(*s, f.residue(*s)))
                    .collect()
            })
            .collect()
    }

    /// Session slots: RT graph first (when loaded), then MS2.
    fn slot(&self, ms2: bool) -> usize {
        usize::from(ms2 && self.load_rt)
    }
}

impl Predictor for PeptDeep {
    fn prefix(&self) -> &str {
        "pdeep"
    }
    fn ions(&self) -> &[Ion] {
        if self.load_ms2 {
            &PDEEP_IONS
        } else {
            &[]
        }
    }
    fn has_rt(&self) -> bool {
        self.load_rt
    }
    fn supports_rt(&self, f: &Form) -> bool {
        self.load_rt && self.supported(f)
    }
    fn supports_ms2(&self, f: &Form, z: u8) -> bool {
        self.load_ms2 && z > 0 && f.seq.len() >= 2 && self.supported(f)
    }
    fn batch_key(&self, f: &Form) -> usize {
        f.seq.len()
    }
    fn sessions(&self) -> Result<Vec<OnnxSession>> {
        let mut s = Vec::new();
        if self.load_rt {
            s.push(self.engine.load(&self.dir, PEPTDEEP_FILES[0])?);
        }
        if self.load_ms2 {
            s.push(self.engine.load(&self.dir, PEPTDEEP_FILES[1])?);
        }
        Ok(s)
    }
    fn rt(&self, s: &mut [OnnxSession], forms: &[&Form]) -> Result<Vec<f32>> {
        let names = Self::names(forms);
        let len = forms[0].seq.len();
        let batch = pack_rt_batch(&self.peptides(forms, &names), len, self.meta.mod_dim)?;
        let out = s[self.slot(false)].run(&batch.as_inputs())?;
        let first = out
            .into_iter()
            .next()
            .ok_or_else(|| Error::Inference("PeptDeep RT produced no output".into()))?;
        Ok(first.data)
    }
    fn ms2(&self, s: &mut [OnnxSession], items: &[(&Form, u8)]) -> Result<Vec<Vec<f32>>> {
        let len = items[0].0.seq.len();
        let forms: Vec<&Form> = items.iter().map(|(f, _)| *f).collect();
        let charges: Vec<u8> = items.iter().map(|(_, z)| *z).collect();
        let names = Self::names(&forms);
        let batch = pack_ms2_batch(
            &self.peptides(&forms, &names),
            &charges,
            len,
            self.meta.mod_dim,
            self.nce,
            self.instrument,
        )?;
        let outs = s[self.slot(true)].run(&batch.as_inputs())?;
        let out = outs
            .first()
            .ok_or_else(|| Error::Inference("PeptDeep MS2 produced no output".into()))?;
        if out.shape.len() != 3 || out.shape[1] != len - 1 {
            return Err(Error::Inference(format!(
                "PeptDeep MS2 output shape {:?}, expected [batch, {}, columns]",
                out.shape,
                len - 1
            )));
        }
        let (rows, ncol) = (out.shape[1], out.shape[2]);
        Ok((0..items.len())
            .map(|b| {
                let mut v = Vec::with_capacity(rows * 4);
                for r in 0..rows {
                    for (k, &c) in self.cols.iter().enumerate() {
                        let x = out.data[(b * rows + r) * ncol + c].max(0.0);
                        v.push(if PDEEP_IONS[k].1 > charges[b] { 0.0 } else { x });
                    }
                }
                v
            })
            .collect())
    }
}

// ---------------------------------------------------------------- DeepLC

/// Graph manifest of an exported DeepLC model.
pub const DEEPLC_FILE: &str = "deeplc_rt.json";

/// The DeepLC RT graph.
pub struct DeepLc {
    dir: PathBuf,
    engine: Engine,
}

impl DeepLc {
    /// DeepLC from `dir`; a missing graph manifest is an error.
    pub fn new(dir: &Path, engine: Engine) -> Result<Self> {
        if !dir.join(DEEPLC_FILE).is_file() {
            return Err(Error::InvalidArg(format!(
                "DeepLC model file {} not found",
                dir.join(DEEPLC_FILE).display()
            )));
        }
        Ok(Self {
            dir: dir.to_path_buf(),
            engine,
        })
    }

    /// pyro-Glu goes on the first residue (site 1), as `Q[Gln->pyro-Glu]` in
    /// the ProForma of `scripts/deeplc_reference.py`; every other mod keeps its site.
    pub fn peptide(f: &Form) -> DeepLcPeptide<'_> {
        let site = |s: usize, m: Mod| match m {
            Mod::GlnPyroGlu | Mod::GluPyroGlu => 1,
            _ => s,
        };
        DeepLcPeptide {
            sequence: &f.seq,
            mods: f
                .mods
                .iter()
                .map(|&(s, m)| (site(s, m), m.composition()))
                .collect(),
        }
    }
}

impl Predictor for DeepLc {
    fn prefix(&self) -> &str {
        "deeplc"
    }
    fn ions(&self) -> &[Ion] {
        &[]
    }
    fn has_rt(&self) -> bool {
        true
    }
    fn supports_rt(&self, f: &Form) -> bool {
        deeplc::encode(&Self::peptide(f)).is_ok()
    }
    fn supports_ms2(&self, _: &Form, _: u8) -> bool {
        false
    }
    fn batch_key(&self, _: &Form) -> usize {
        0
    }
    fn sessions(&self) -> Result<Vec<OnnxSession>> {
        Ok(vec![self.engine.load(&self.dir, DEEPLC_FILE)?])
    }
    fn rt(&self, s: &mut [OnnxSession], forms: &[&Form]) -> Result<Vec<f32>> {
        let peps: Vec<DeepLcPeptide<'_>> = forms.iter().map(|f| Self::peptide(f)).collect();
        let inputs = deeplc::pack_batch(&peps)?;
        let out = s[0].run(&inputs)?;
        let first = out
            .into_iter()
            .next()
            .ok_or_else(|| Error::Inference("DeepLC produced no output".into()))?;
        Ok(first.data)
    }
    fn ms2(&self, _: &mut [OnnxSession], _: &[(&Form, u8)]) -> Result<Vec<Vec<f32>>> {
        Err(Error::InvalidArg("DeepLC has no MS2 model".into()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn engine() -> Engine {
        Engine {
            kind: BackendKind::Tract,
            opts: BackendOptions::default(),
        }
    }

    #[test]
    fn missing_model_files_are_named_in_the_error() {
        let dir = tempfile::tempdir().unwrap();
        let e = PeptDeep::new(dir.path(), 30.0, "QE", engine(), false, true)
            .err()
            .unwrap()
            .to_string();
        assert!(e.contains("peptdeep_ms2_dynamic.json"), "{e}");
        let e = DeepLc::new(dir.path(), engine()).err().unwrap().to_string();
        assert!(e.contains("deeplc_rt.json"), "{e}");
        let e = engine()
            .load(dir.path(), "x.json")
            .err()
            .unwrap()
            .to_string();
        assert!(e.contains("x.json"), "{e}");
    }

    #[test]
    fn peptdeep_marks_forms_with_unknown_mods_unsupported() {
        let dir = tempfile::tempdir().unwrap();
        for f in &PEPTDEEP_FILES[..2] {
            std::fs::write(dir.path().join(f), "{}").unwrap();
        }
        let meta = r#"{"mod_dim": 2, "mods": {"Carbamidomethyl@C": [1.0, 0.0]},
            "ms2_frag_types": ["b_z1","b_z2","y_z1","y_z2"],
            "instrument_index": {"QE": 0}, "unknown_instrument_index": 7}"#;
        std::fs::write(dir.path().join(PEPTDEEP_FILES[2]), meta).unwrap();
        let p = PeptDeep::new(dir.path(), 30.0, "qe", engine(), true, true).unwrap();
        let (cam, _) = form::parse_peptide("K.PEPC+57.02146K.R").unwrap();
        let (ox, _) = form::parse_peptide("K.PEPM+15.99491K.R").unwrap();
        assert!(p.supports_rt(&cam) && p.supports_ms2(&cam, 2));
        assert!(!p.supports_rt(&ox) && !p.supports_ms2(&ox, 2));
        assert!(!p.supports_ms2(&cam, 0));
        let rt_only = PeptDeep::new(dir.path(), 30.0, "QE", engine(), true, false).unwrap();
        assert!(rt_only.ions().is_empty() && !rt_only.supports_ms2(&cam, 2));
    }
}
