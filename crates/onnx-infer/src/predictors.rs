pub mod calibration;
pub mod form;
pub mod spectral;

use std::collections::HashMap;
use std::fs::File;
use std::path::{Path, PathBuf};

use serde::Deserialize;

use crate::deeplc::{self, DeepLcPeptide};
use crate::error::{Error, Result};
use crate::ms2::{Ms2Conditions, Ms2Model, Ms2Scratch, N_BASE_ION, N_ION};
use crate::peptdeep::{pack_ms2_batch, pack_rt_batch, PeptDeepPeptide};
use crate::rt::{RtModel, Scratch};
use crate::{BackendKind, BackendOptions, OnnxSession};

use form::{Form, Mod};
use spectral::Ion;

#[derive(Clone, Debug)]
pub struct Engine {
    pub kind: BackendKind,
    pub opts: BackendOptions,
}

impl Engine {
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

pub trait Predictor: Sync {
    fn prefix(&self) -> &str;
    fn ions(&self) -> &[Ion];
    fn has_rt(&self) -> bool;
    fn supports_rt(&self, form: &Form) -> bool;
    fn supports_ms2(&self, form: &Form, charge: u8) -> bool;
    fn batch_key(&self, form: &Form) -> usize;
    fn sessions(&self) -> Result<Vec<OnnxSession>>;
    fn rt(&self, s: &mut [OnnxSession], forms: &[&Form]) -> Result<Vec<f32>>;
    fn ms2(&self, s: &mut [OnnxSession], items: &[(&Form, u8)]) -> Result<Vec<Vec<f32>>>;
}

pub struct NativeRt {
    model: RtModel,
}

impl NativeRt {
    pub fn from_bytes(buf: &[u8]) -> Result<Self> {
        let model = RtModel::from_bytes(buf)
            .map_err(|e| Error::InvalidArg(format!("andes-rt weights: {e}")))?;
        Ok(Self { model })
    }

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

const MS2_IONS: [Ion; 4] = [(true, 1), (true, 2), (false, 1), (false, 2)];

pub struct NativeMs2 {
    model: Ms2Model,
    nce: Option<f32>,
    activation: String,
    instrument: String,
}

impl NativeMs2 {
    pub fn from_bytes(
        buf: &[u8],
        nce: Option<f32>,
        activation: &str,
        instrument: &str,
    ) -> Result<Self> {
        let model = Ms2Model::from_bytes(buf)
            .map_err(|e| Error::InvalidArg(format!("andes-ms2 weights: {e}")))?;
        Ok(Self {
            model,
            nce,
            activation: activation.to_string(),
            instrument: instrument.to_string(),
        })
    }

    fn conditions(&self, charge: u8) -> Ms2Conditions {
        Ms2Conditions::new(charge, self.nce, &self.activation, &self.instrument)
    }
}

impl Predictor for NativeMs2 {
    fn prefix(&self) -> &str {
        "andesms2"
    }
    fn ions(&self) -> &[Ion] {
        &MS2_IONS
    }
    fn has_rt(&self) -> bool {
        false
    }
    fn supports_rt(&self, _: &Form) -> bool {
        false
    }
    fn supports_ms2(&self, f: &Form, z: u8) -> bool {
        z > 0 && f.seq.len() >= 2 && f.rt_peptide().encode().is_some()
    }
    fn batch_key(&self, _: &Form) -> usize {
        0
    }
    fn sessions(&self) -> Result<Vec<OnnxSession>> {
        Ok(Vec::new())
    }
    fn rt(&self, _: &mut [OnnxSession], _: &[&Form]) -> Result<Vec<f32>> {
        Err(Error::InvalidArg("andes-ms2 has no RT model".into()))
    }
    fn ms2(&self, _: &mut [OnnxSession], items: &[(&Form, u8)]) -> Result<Vec<Vec<f32>>> {
        let mut scratch = Ms2Scratch::default();
        items
            .iter()
            .map(|(f, z)| {
                let full = self
                    .model
                    .predict_with(&f.rt_peptide(), &self.conditions(*z), &mut scratch)
                    .ok_or_else(|| {
                        Error::Inference("andes-ms2 could not encode a peptidoform".into())
                    })?;
                Ok(full
                    .chunks_exact(N_ION)
                    .flat_map(|site| site[..N_BASE_ION].to_vec())
                    .collect())
            })
            .collect()
    }
}

#[derive(Deserialize)]
pub struct PeptDeepMeta {
    pub mod_dim: usize,
    pub mods: HashMap<String, Vec<f32>>,
    pub ms2_frag_types: Vec<String>,
    pub instrument_index: HashMap<String, i64>,
    pub unknown_instrument_index: i64,
}

pub const PEPTDEEP_FILES: [&str; 3] = [
    "peptdeep_rt_dynamic.json",
    "peptdeep_ms2_dynamic.json",
    "peptdeep_mods.json",
];

pub struct PeptDeep {
    dir: PathBuf,
    meta: PeptDeepMeta,
    nce: f32,
    instrument: i64,
    cols: [usize; 4],
    engine: Engine,
    load_rt: bool,
    load_ms2: bool,
}

const PDEEP_IONS: [Ion; 4] = [(true, 1), (true, 2), (false, 1), (false, 2)];

impl PeptDeep {
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

    pub fn known_mods(&self) -> impl Iterator<Item = &str> {
        self.meta.mods.keys().map(String::as_str)
    }

    pub fn ms2_frag_types(&self) -> &[String] {
        &self.meta.ms2_frag_types
    }

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

pub const DEEPLC_FILE: &str = "deeplc_rt.json";

pub struct DeepLc {
    dir: PathBuf,
    engine: Engine,
}

impl DeepLc {
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
    fn native_ms2_predicts_the_four_plain_channels_per_site() {
        let buf = std::fs::read(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/fixtures/ms2_native/tiny.ams2"
        ))
        .unwrap();
        let p = NativeMs2::from_bytes(&buf, Some(28.0), "HCD-FTMS", "Q Exactive").unwrap();
        assert_eq!(p.prefix(), "andesms2");
        assert_eq!(p.ions().len(), 4);
        let (f, _) = form::parse_peptide("K.PEPM+15.99491TIDEK.R").unwrap();
        assert!(p.supports_ms2(&f, 2));
        assert!(!p.supports_ms2(&f, 0));
        let out = p.ms2(&mut [], &[(&f, 2), (&f, 3)]).unwrap();
        let sites = f.seq.len() - 1;
        assert_eq!(out.len(), 2);
        assert!(out.iter().all(|v| v.len() == sites * 4));
        assert!(out.iter().flatten().all(|x| x.is_finite()));
        assert!(p.rt(&mut [], &[&f]).is_err());
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
