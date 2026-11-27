//! Identification scores of the study's PSMs: loaded from the rescorer or Percolator output,
//! keyed by spectrum and peptidoform.

use std::collections::{BTreeMap, HashMap};
use std::path::{Path, PathBuf};

use model::isomers::{bundled_unimod_mass, mod_units_match};
use model::IsomerClass;
use output::PercolatorPsm;

#[derive(Clone, Default)]
pub(super) struct SpecIdKeys(std::collections::hash_map::RandomState);

impl SpecIdKeys {
    fn hash(&self, tag: u8, bytes: &[u8]) -> u64 {
        use std::hash::{BuildHasher, Hasher};
        let mut h = self.0.build_hasher();
        h.write_u8(tag);
        h.write(bytes);
        h.finish()
    }

    pub(super) fn key(&self, spec_id: &str) -> u128 {
        (u128::from(self.hash(0, spec_id.as_bytes())) << 64)
            | u128::from(self.hash(1, spec_id.as_bytes()))
    }

    fn residues(&self, residues: &[u8]) -> u64 {
        self.hash(2, residues)
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Scored {
    pub score: f64,
    pub q: f64,
    pub pooled_q: f64,
    pub pep: f64,
    tag: u64,
}

const RESIDUE_BITS: u32 = 40;
const MOD_SET_BITS: u32 = 64 - RESIDUE_BITS - 1;

impl Scored {
    pub(crate) fn is_decoy(&self) -> bool {
        self.tag & 1 == 1
    }

    fn residues(&self) -> u64 {
        self.tag >> (MOD_SET_BITS + 1)
    }

    fn mod_set(&self) -> usize {
        ((self.tag >> 1) & ((1 << MOD_SET_BITS) - 1)) as usize
    }
}

pub(crate) struct ScoreSet {
    pub source: String,
    pub(super) keys: SpecIdKeys,
    pub(super) psms: HashMap<u128, Scored, rustc_hash::FxBuildHasher>,
    mod_sets: Vec<Box<[i64]>>,
    mod_set_ids: HashMap<Box<[i64]>, u32, rustc_hash::FxBuildHasher>,
    pub(super) class: IsomerClass,
}

impl ScoreSet {
    pub(super) fn with_capacity(source: String, n: usize) -> Self {
        Self {
            source,
            keys: SpecIdKeys::default(),
            psms: HashMap::with_capacity_and_hasher(n, rustc_hash::FxBuildHasher),
            mod_sets: Vec::new(),
            mod_set_ids: HashMap::default(),
            class: IsomerClass::default(),
        }
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn insert(
        &mut self,
        file: &str,
        spec_id: &str,
        score: Option<f64>,
        q: f64,
        pooled_q: Option<f64>,
        pep: f64,
        peptide: &str,
        is_decoy: bool,
    ) -> Result<(), String> {
        let score = score.ok_or_else(|| {
            format!(
                "{file}: PSM {spec_id:?} has no finite `score` (missing column or value); the library ranks precursors by the rescorer's PSM score"
            )
        })?;
        self.class
            .parse_into(peptide, bundled_unimod_mass)
            .map_err(|e| format!("{file}: PSM {spec_id:?}: {e}"))?;
        let mods = match self.mod_set_ids.get(self.class.mod_units()) {
            Some(&i) => i,
            None => {
                let i = self.mod_sets.len() as u64;
                if i >> MOD_SET_BITS != 0 {
                    return Err(format!(
                        "{}: more than 2^{MOD_SET_BITS} distinct modification combinations",
                        self.source
                    ));
                }
                let set: Box<[i64]> = self.class.mod_units().into();
                self.mod_sets.push(set.clone());
                self.mod_set_ids.insert(set, i as u32);
                i as u32
            }
        };
        let tag = ((self.keys.residues(self.class.residues()) >> (64 - RESIDUE_BITS))
            << (MOD_SET_BITS + 1))
            | (u64::from(mods) << 1)
            | u64::from(is_decoy);
        let scored = Scored {
            score,
            q,
            pooled_q: pooled_q.unwrap_or(f64::NAN),
            pep,
            tag,
        };
        match self.psms.entry(self.keys.key(spec_id)) {
            std::collections::hash_map::Entry::Occupied(o) => {
                if o.get().is_decoy() == is_decoy {
                    Err(format!(
                        "{file}: {}",
                        output::PercolatorError::DuplicatePsmId {
                            psm_id: spec_id.to_string()
                        }
                    ))
                } else {
                    Err(format!(
                        "{}: PSM {spec_id:?} is both a target and a decoy",
                        self.source
                    ))
                }
            }
            std::collections::hash_map::Entry::Vacant(v) => {
                v.insert(scored);
                Ok(())
            }
        }
    }

    pub(crate) fn get(&self, spec_id: &str) -> Option<&Scored> {
        self.psms.get(&self.keys.key(spec_id))
    }

    #[cfg(test)]
    pub(super) fn len(&self) -> usize {
        self.psms.len()
    }

    pub(super) fn agrees(&self, s: &Scored, is_decoy: bool, class: &IsomerClass) -> bool {
        s.is_decoy() == is_decoy
            && s.residues() == self.keys.residues(class.residues()) >> (64 - RESIDUE_BITS)
            && mod_units_match(&self.mod_sets[s.mod_set()], class.mod_units())
    }
}

pub(crate) fn rescore_scores(
    targets: &HashMap<String, PercolatorPsm>,
    decoys: &HashMap<String, PercolatorPsm>,
    source: &str,
) -> Result<ScoreSet, String> {
    let mut set = ScoreSet::with_capacity(source.to_string(), targets.len() + decoys.len());
    for (map, is_decoy) in [(targets, false), (decoys, true)] {
        for (id, p) in map {
            set.insert(
                source,
                id,
                p.score,
                p.q_value,
                p.pooled_q_value,
                p.pep,
                &p.peptide,
                is_decoy,
            )?;
        }
    }
    Ok(set)
}

pub(super) const MAX_JOIN_MISMATCH: f64 = 0.01;

pub(crate) enum ScoreLoader<'a> {
    Merged(Box<dyn FnOnce(usize) -> Result<ScoreSet, String> + 'a>),
    PerRun(RunScores<'a>),
}

pub(crate) type RunScores<'a> = Box<dyn Fn(&str) -> Result<ScoreSet, String> + Sync + 'a>;

#[derive(Debug, PartialEq)]
pub(super) enum ScoreSource {
    Merged(PathBuf, PathBuf),
    PerRun(BTreeMap<String, (PathBuf, PathBuf)>),
}

pub(super) fn score_source(
    psms: &[String],
    decoys: &[String],
    stems: &[String],
) -> Result<ScoreSource, String> {
    let split = |v: &str| -> Option<(String, PathBuf)> {
        let (stem, path) = v.split_once('=')?;
        stems
            .iter()
            .any(|s| s == stem)
            .then(|| (stem.to_string(), PathBuf::from(path)))
    };
    let per_run = |flag: &str, values: &[String]| -> Result<BTreeMap<String, PathBuf>, String> {
        let mut m = BTreeMap::new();
        for v in values {
            let (stem, path) = split(v).ok_or_else(|| {
                format!("{flag} {v:?}: per-run files are <stem>=<file> with a stem of the library's runs; plain paths and per-run files cannot be mixed")
            })?;
            if m.insert(stem.clone(), path).is_some() {
                return Err(format!("{flag}: run {stem:?} is given twice"));
            }
        }
        let missing: Vec<&str> = stems
            .iter()
            .filter(|s| !m.contains_key(*s))
            .map(String::as_str)
            .collect();
        if !missing.is_empty() {
            return Err(format!("{flag}: no file for runs {}", missing.join(", ")));
        }
        Ok(m)
    };
    let any_per_run = psms.iter().chain(decoys).any(|v| split(v).is_some());
    if !any_per_run {
        return match (psms, decoys) {
            ([t], [d]) => Ok(ScoreSource::Merged(t.into(), d.into())),
            _ => Err(format!(
                "give one merged --psms and --decoy-psms file, or <stem>=<file> for every run (got {} and {})",
                psms.len(),
                decoys.len()
            )),
        };
    }
    let t = per_run("--psms", psms)?;
    let mut d = per_run("--decoy-psms", decoys)?;
    Ok(ScoreSource::PerRun(
        t.into_iter()
            .map(|(stem, tp)| {
                let dp = d.remove(&stem).expect("every run has a decoy file");
                (stem, (tp, dp))
            })
            .collect(),
    ))
}

pub(super) fn score_loader(source: &ScoreSource) -> ScoreLoader<'_> {
    match source {
        ScoreSource::Merged(t, d) => {
            ScoreLoader::Merged(Box::new(move |n: usize| percolator_scores(t, d, n)))
        }
        ScoreSource::PerRun(m) => ScoreLoader::PerRun(Box::new(move |stem: &str| {
            let (t, d) = m
                .get(stem)
                .ok_or_else(|| format!("no rescoring scores for run {stem:?}"))?;
            percolator_scores(t, d, 0)
        })),
    }
}

pub(super) fn percolator_scores(
    targets: &Path,
    decoys: &Path,
    capacity: usize,
) -> Result<ScoreSet, String> {
    let mut set = ScoreSet::with_capacity(
        format!("{} and {}", targets.display(), decoys.display()),
        capacity,
    );
    for (path, is_decoy) in [(targets, false), (decoys, true)] {
        let file = path.display().to_string();
        let err = |e: output::PercolatorError| format!("{file}: {e}");
        let f = std::fs::File::open(path).map_err(|e| format!("{file}: {e}"))?;
        let mut reader =
            output::PsmResultReader::new(std::io::BufReader::with_capacity(1 << 20, f))
                .map_err(err)?;
        while let Some(row) = reader.next_row().map_err(err)? {
            set.insert(
                &file,
                row.psm_id,
                row.score,
                row.q_value,
                row.pooled_q,
                row.pep,
                row.peptide,
                is_decoy,
            )?;
        }
    }
    Ok(set)
}
