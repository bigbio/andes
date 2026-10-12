//! Shared by the `rescore_pin` and `model_parity` examples: andes peptide
//! parsing and one predictor per supported model family.
//!
//! A predictor turns peptidoforms into retention-time and/or fragment
//! intensity predictions with tract. Peptidoforms a model cannot encode are
//! rejected by [`Predictor::supports_rt`] / [`Predictor::supports_ms2`] and
//! become missing features downstream, never a silently wrong tensor.

// Each example uses a different subset.
#![allow(dead_code)]

use std::collections::HashMap;
use std::fs::File;
use std::path::Path;

use onnx_infer::peptdeep::{PeptDeepPeptide, pack_ms2_batch, pack_rt_batch};
use onnx_infer::{BackendKind, OnnxSession, ms2pip};
use serde::Deserialize;

pub const PROTON: f64 = 1.007_276_466_8;
pub const H2O: f64 = 18.010_564_684;

pub fn residue_mass(c: char) -> Option<f64> {
    Some(match c {
        'G' => 57.021_463_72,
        'A' => 71.037_113_79,
        'S' => 87.032_028_41,
        'P' => 97.052_763_85,
        'V' => 99.068_413_91,
        'T' => 101.047_678_5,
        'C' => 103.009_184_5,
        'L' | 'I' => 113.084_064,
        'N' => 114.042_927_4,
        'D' => 115.026_943,
        'Q' => 128.058_577_5,
        'K' => 128.094_963,
        'E' => 129.042_593_1,
        'M' => 131.040_484_6,
        'H' => 137.058_911_9,
        'F' => 147.068_413_9,
        'R' => 156.101_111,
        'Y' => 163.063_328_6,
        'W' => 186.079_312_9,
        _ => return None,
    })
}

/// Modifications the examples recognise in andes peptide strings.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Mod {
    Carbamidomethyl,
    Oxidation,
    Acetyl,
    Tmt6plex,
    Phospho,
    Deamidated,
    GlnPyroGlu,
    GluPyroGlu,
}

/// Where a modification may sit.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Placement {
    /// On one of the listed residues.
    Residue,
    /// Peptide N-terminus only (site 0).
    NTerm,
    /// On one of the listed residues, or the peptide N-terminus.
    ResidueOrNTerm,
    /// On the first residue only, which must be one of the listed residues.
    FirstResidue,
}

/// (mod, monoisotopic delta, residues, placement)
const MODS: &[(Mod, f64, &str, Placement)] = &[
    (Mod::Carbamidomethyl, 57.021_464, "C", Placement::Residue),
    (Mod::Oxidation, 15.994_915, "M", Placement::Residue),
    (Mod::Acetyl, 42.010_565, "", Placement::NTerm),
    (Mod::Tmt6plex, 229.162_932, "K", Placement::ResidueOrNTerm),
    (Mod::Phospho, 79.966_331, "STY", Placement::Residue),
    (Mod::Deamidated, 0.984_016, "NQ", Placement::Residue),
    (Mod::GlnPyroGlu, -17.026_549, "Q", Placement::FirstResidue),
    (Mod::GluPyroGlu, -18.010_565, "E", Placement::FirstResidue),
];

impl Mod {
    /// alphabase / peptdeep name at a site (0 = N-terminus).
    pub fn alphabase(self, site: usize, aa: char) -> String {
        match self {
            Mod::Carbamidomethyl => "Carbamidomethyl@C".into(),
            Mod::Oxidation => "Oxidation@M".into(),
            Mod::Acetyl => "Acetyl@Protein_N-term".into(),
            Mod::Tmt6plex if site == 0 => "TMT6plex@Any_N-term".into(),
            Mod::Tmt6plex => "TMT6plex@K".into(),
            Mod::Phospho => format!("Phospho@{aa}"),
            Mod::Deamidated => format!("Deamidated@{aa}"),
            Mod::GlnPyroGlu => "Gln->pyro-Glu@Q^Any_N-term".into(),
            Mod::GluPyroGlu => "Glu->pyro-Glu@E^Any_N-term".into(),
        }
    }

    /// Elemental composition delta (C, H, N, O, S, P), unimod; heavy
    /// isotopes counted as their element.
    pub fn composition(self) -> [i32; 6] {
        match self {
            Mod::Carbamidomethyl => [2, 3, 1, 1, 0, 0],
            Mod::Oxidation => [0, 0, 0, 1, 0, 0],
            Mod::Acetyl => [2, 2, 0, 1, 0, 0],
            Mod::Tmt6plex => [12, 20, 2, 2, 0, 0],
            Mod::Phospho => [0, 1, 0, 3, 0, 1],
            Mod::Deamidated => [0, -1, -1, 1, 0, 0],
            Mod::GlnPyroGlu => [0, -3, -1, 0, 0, 0],
            Mod::GluPyroGlu => [0, -2, 0, -1, 0, 0],
        }
    }
}

/// One parsed peptidoform.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Form {
    pub seq: String,
    /// (site, mod) sorted — site 0 = N-term, 1..=n residue.
    pub mods: Vec<(usize, Mod)>,
}

impl Form {
    pub fn residue(&self, site: usize) -> char {
        self.seq.as_bytes()[site.max(1) - 1] as char
    }
}

/// Fragment-mass inputs for one peptidoform.
pub struct Mass {
    /// Residue + residue-mod masses, len n.
    pub res: Vec<f64>,
    /// N-terminal mod mass.
    pub nterm: f64,
}

/// Parse `K.PEPM+15.99491K.D` into a peptidoform; `None` for anything outside
/// the residue alphabet or [`MODS`]. andes prints N-terminal deltas after the
/// first residue, so a delta there is placed on the N-terminus or the residue
/// by [`Placement`].
pub fn parse_peptide(s: &str) -> Option<(Form, Mass)> {
    let b = s.as_bytes();
    let core = if b.len() > 4 && b[1] == b'.' && b[b.len() - 2] == b'.' {
        &s[2..s.len() - 2]
    } else {
        s
    };
    let mut seq = String::new();
    let mut res = Vec::new();
    let mut mods: Vec<(usize, Mod)> = Vec::new();
    let mut nterm = 0.0;
    let chars: Vec<char> = core.chars().collect();
    let mut i = 0;
    while i < chars.len() {
        let c = chars[i];
        if c.is_ascii_uppercase() {
            seq.push(c);
            res.push(residue_mass(c)?);
            i += 1;
        } else if c == '+' || c == '-' {
            let start = i;
            i += 1;
            while i < chars.len() && (chars[i].is_ascii_digit() || chars[i] == '.') {
                i += 1;
            }
            let delta: f64 = chars[start..i].iter().collect::<String>().parse().ok()?;
            let pos = seq.len();
            if pos == 0 {
                return None;
            }
            let aa = seq.as_bytes()[pos - 1] as char;
            let free = |site: usize, mods: &[(usize, Mod)]| mods.iter().all(|(s, _)| *s != site);
            let cands: Vec<(Mod, &str, Placement)> = MODS
                .iter()
                .filter(|m| (delta - m.1).abs() < 2e-3)
                .map(|m| (m.0, m.2, m.3))
                .collect();
            let on_res = |p: Placement| matches!(p, Placement::Residue | Placement::ResidueOrNTerm);
            let pick = cands
                .iter()
                .find(|(_, _, p)| pos == 1 && *p == Placement::NTerm && free(0, &mods))
                .map(|(m, _, _)| (0, *m))
                .or_else(|| {
                    cands
                        .iter()
                        .find(|(_, r, p)| {
                            pos == 1 && *p == Placement::FirstResidue && r.contains(aa) && free(1, &mods)
                        })
                        .map(|(m, _, _)| (1, *m))
                })
                .or_else(|| {
                    cands
                        .iter()
                        .find(|(_, r, p)| on_res(*p) && r.contains(aa) && free(pos, &mods))
                        .map(|(m, _, _)| (pos, *m))
                })
                .or_else(|| {
                    cands
                        .iter()
                        .find(|(_, _, p)| pos == 1 && *p == Placement::ResidueOrNTerm && free(0, &mods))
                        .map(|(m, _, _)| (0, *m))
                })?;
            if pick.0 == 0 {
                nterm += delta;
            } else {
                res[pos - 1] += delta;
            }
            mods.push(pick);
        } else {
            return None;
        }
    }
    if seq.is_empty() {
        return None;
    }
    mods.sort();
    Some((Form { seq, mods }, Mass { res, nterm }))
}

/// One fragment ion type in a predictor's per-site output: (is b ion, charge).
pub type Ion = (bool, u8);

/// A model family that predicts RT and/or fragment intensities.
pub trait Predictor: Sync {
    /// Feature-column prefix (`pdeep`, `ms2pip`, …).
    fn prefix(&self) -> &str;
    /// Per-site ion layout of [`Self::ms2`] output; empty when no MS2 model.
    fn ions(&self) -> &[Ion];
    fn has_rt(&self) -> bool;
    fn supports_rt(&self, form: &Form) -> bool;
    fn supports_ms2(&self, form: &Form, charge: u8) -> bool;
    /// Forms in one RT / MS2 batch must share this key.
    fn batch_key(&self, form: &Form) -> usize;
    /// Load the model graphs once; worker threads get [`OnnxSession::share`]d handles.
    fn sessions(&self) -> Vec<OnnxSession>;
    fn rt(&self, s: &mut [OnnxSession], forms: &[&Form]) -> Vec<f32>;
    /// Per form: (len − 1) sites × `ions().len()` linear intensities, site `i`
    /// = cleavage after residue `i + 1`; ions above the precursor charge are 0.
    fn ms2(&self, s: &mut [OnnxSession], items: &[(&Form, u8)]) -> Vec<Vec<f32>>;
}

pub fn load(models: &Path, manifest: &str) -> OnnxSession {
    OnnxSession::from_manifest(models.join(manifest), BackendKind::Tract)
        .unwrap_or_else(|e| panic!("load {manifest}: {e}"))
}

// ---------------------------------------------------------------- PeptDeep

#[derive(Deserialize)]
pub struct PeptDeepMeta {
    pub mod_dim: usize,
    pub mods: HashMap<String, Vec<f32>>,
    pub ms2_frag_types: Vec<String>,
    pub instrument_index: HashMap<String, i64>,
    pub unknown_instrument_index: i64,
}

pub struct PeptDeep {
    pub dir: std::path::PathBuf,
    pub meta: PeptDeepMeta,
    pub nce: f32,
    pub instrument: i64,
    /// Columns of b_z1, b_z2, y_z1, y_z2 in the MS2 output.
    pub cols: [usize; 4],
}

const PDEEP_IONS: [Ion; 4] = [(true, 1), (true, 2), (false, 1), (false, 2)];

impl PeptDeep {
    pub fn new(dir: &Path, nce: f32, instrument: &str) -> Self {
        let meta: PeptDeepMeta =
            serde_json::from_reader(File::open(dir.join("peptdeep_mods.json")).unwrap()).unwrap();
        let col = |n: &str| meta.ms2_frag_types.iter().position(|f| f == n).unwrap();
        let cols = [col("b_z1"), col("b_z2"), col("y_z1"), col("y_z2")];
        let instrument = *meta
            .instrument_index
            .get(&instrument.to_uppercase())
            .unwrap_or(&meta.unknown_instrument_index);
        Self { dir: dir.to_path_buf(), meta, nce, instrument, cols }
    }

    fn supported(&self, f: &Form) -> bool {
        f.mods.iter().all(|(s, m)| self.meta.mods.contains_key(&m.alphabase(*s, f.residue(*s))))
    }

    fn peptides<'f>(&'f self, forms: &[&'f Form], names: &'f [Vec<String>]) -> Vec<PeptDeepPeptide<'f>> {
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
            .map(|f| f.mods.iter().map(|(s, m)| m.alphabase(*s, f.residue(*s))).collect())
            .collect()
    }
}

impl Predictor for PeptDeep {
    fn prefix(&self) -> &str {
        "pdeep"
    }
    fn ions(&self) -> &[Ion] {
        &PDEEP_IONS
    }
    fn has_rt(&self) -> bool {
        true
    }
    fn supports_rt(&self, f: &Form) -> bool {
        self.supported(f)
    }
    fn supports_ms2(&self, f: &Form, z: u8) -> bool {
        z > 0 && f.seq.len() >= 2 && self.supported(f)
    }
    fn batch_key(&self, f: &Form) -> usize {
        f.seq.len()
    }
    fn sessions(&self) -> Vec<OnnxSession> {
        vec![
            load(&self.dir, "peptdeep_rt_dynamic.json"),
            load(&self.dir, "peptdeep_ms2_dynamic.json"),
        ]
    }
    fn rt(&self, s: &mut [OnnxSession], forms: &[&Form]) -> Vec<f32> {
        let names = Self::names(forms);
        let len = forms[0].seq.len();
        let batch = pack_rt_batch(&self.peptides(forms, &names), len, self.meta.mod_dim).unwrap();
        s[0].run(&batch.as_inputs()).expect("rt run")[0].data.clone()
    }
    fn ms2(&self, s: &mut [OnnxSession], items: &[(&Form, u8)]) -> Vec<Vec<f32>> {
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
        )
        .unwrap();
        let out = &s[1].run(&batch.as_inputs()).expect("ms2 run")[0];
        let (rows, ncol) = (out.shape[1], out.shape[2]);
        assert_eq!(rows, len - 1, "MS2 output rows");
        (0..items.len())
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
            .collect()
    }
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
    fn sessions(&self) -> Vec<OnnxSession> {
        ["b", "y"]
            .iter()
            .map(|ion| load(&self.dir, &format!("ms2pip_{}_{ion}.json", self.model)))
            .collect()
    }
    fn rt(&self, _: &mut [OnnxSession], _: &[&Form]) -> Vec<f32> {
        unreachable!("MS²PIP has no RT model")
    }
    fn ms2(&self, s: &mut [OnnxSession], items: &[(&Form, u8)]) -> Vec<Vec<f32>> {
        let peps: Vec<(&str, u8)> = items.iter().map(|(f, z)| (f.seq.as_str(), *z)).collect();
        let batch = ms2pip::pack_batch(&peps).unwrap();
        let input = std::slice::from_ref(&batch.features);
        let b = ms2pip::log2_predictions(&s[0].run(input).expect("ms2pip b")[0]).unwrap();
        let y = ms2pip::log2_predictions(&s[1].run(input).expect("ms2pip y")[0]).unwrap();
        batch
            .offsets
            .windows(2)
            .map(|w| {
                (w[0]..w[1])
                    .flat_map(|r| [ms2pip::to_linear(b[r]), ms2pip::to_linear(y[r])])
                    .collect()
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_andes_peptides() {
        let (f, m) = parse_peptide("K.M+15.99491EEHILR.P").unwrap();
        assert_eq!(f.seq, "MEEHILR");
        assert_eq!(f.mods, vec![(1, Mod::Oxidation)]);
        assert!((m.res[0] - (131.0404846 + 15.99491)).abs() < 1e-6);
        let (f, m) = parse_peptide("M.R+42.01056AC+57.02146INNQQIR.H").unwrap();
        assert_eq!(f.seq, "RACINNQQIR");
        assert_eq!(f.mods, vec![(0, Mod::Acetyl), (3, Mod::Carbamidomethyl)]);
        assert!((m.nterm - 42.01056).abs() < 1e-9);
        let (f, _) = parse_peptide("_.M+42.01056+15.99491LRK.A").unwrap();
        assert_eq!(f.mods, vec![(0, Mod::Acetyl), (1, Mod::Oxidation)]);
        assert!(parse_peptide("K.PEPXIDE.K").is_none());
        assert!(parse_peptide("K.PEP+12.34TIDE.K").is_none());
    }

    #[test]
    fn places_tmt_and_phospho() {
        let (f, m) = parse_peptide("K.K+229.16293+229.16293PEPTK+229.16293.D").unwrap();
        assert_eq!(f.mods, vec![(0, Mod::Tmt6plex), (1, Mod::Tmt6plex), (6, Mod::Tmt6plex)]);
        assert!((m.nterm - 229.16293).abs() < 1e-9);
        let (f, _) = parse_peptide("K.A+229.16293PES+79.96633K+229.16293.D").unwrap();
        assert_eq!(f.mods, vec![(0, Mod::Tmt6plex), (4, Mod::Phospho), (5, Mod::Tmt6plex)]);
        let (f, _) = parse_peptide("K.Q-17.02655PEPTIDEK.D").unwrap();
        assert_eq!(f.mods, vec![(1, Mod::GlnPyroGlu)]);
        assert!(parse_peptide("K.AQ-17.02655PEPTIDEK.D").is_none());
    }
}
