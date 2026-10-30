//! The match-between-runs stages as commands: `andes quant library`,
//! `andes quant extract` and `andes quant finalize`, and the in-process run of
//! all three behind `andes search --lfq --mbr`. Both paths call the same
//! functions on the same files in the quant directory, so a scattered run and
//! an in-process run give the same `<base>.mbr.parquet`.
//!
//! Files in the quant directory, per run stem: `<stem>.quant.parquet`,
//! `<stem>.ms1.parquet` and `<stem>.ms2.parquet` (written by `andes search
//! --quant-dir`),
//! `<stem>.extract.parquet` (extract); once per experiment: `library.parquet`
//! and `alignment.parquet` (library).
//!
//! Run index: the position of the run's stem in the sorted stem list.
//! Precursor: a positional-isomer class of the run tables' precursor keys
//! (residues, modification masses and charge, [`model::IsomerClass`]), so the
//! site forms of one peptide that win in different spectra or runs are one
//! precursor; its key is the `<peptidoform>/<charge>` of the site form of its
//! best PSM. Precursor id: the position of its class in the sorted union of
//! every run's classes. The rescorer's PSM score, PSM q-value and PEP
//! (Percolator's `score`, `q-value` and `posterior_error_prob` columns, or the
//! native rescorer's) are joined to the run tables on (run, SpecId). Library
//! confidence ranks precursors by the score; the q-value, and the PEP when
//! `--pep` is given, pick the donor PSMs, as they pick the PSMs of the per-run
//! quantification.

use std::collections::{BTreeMap, HashMap};
use std::path::{Path, PathBuf};

use model::isomer_class::{bundled_unimod_mass, mod_units_match};
use model::{IsomerClass, Tolerance};
use output::extract_table::{read_extract_tags, write_extract_table_tagged};
use output::feature_parquet::{FeatureGroup, FeatureRecord, MbrColumns};
use output::library_table::{
    library_fingerprint, read_library, write_library, StoredLibrary, FINGERPRINT_KEY,
};
use output::mbr_out::EvidenceType;
use output::ms2_table::read_ms2_table;
use output::run_quant_table::{open_run_quant_table, read_run_quant_mods};
use output::{
    FeatureCounts, FeatureParquetWriter, GroupAt, Inference, ModRecord, PercolatorPsm,
    ProteinPosition,
};
use quant::align::AlignParams;
use quant::evidence::EvidenceParams;
use quant::library::{Library, LibraryBuilder, PrecursorMeta, PsmIn};
use quant::ms2::{fragment_ladder, library_spectrum, residue_masses, LibrarySpectrum, Ms2RunIndex};
use rayon::prelude::*;

use crate::cli::{QuantCommand, QuantExtractArgs, QuantFinalizeArgs, QuantLibraryArgs};
use crate::quant_finalize::{
    finalize_stream, ExtractFiles, ExtractSource, FinalizeLimits, FinalizeParams, MbrFeature,
};
use crate::quant_mbr::{
    attach_ms2, extract_run, mbr_paths, ExtractCounts, ExtractRow, Ms2Counts, TransferSelection,
};

/// Seed of the transfer FDR's folds and models, in every path.
pub(crate) const MBR_SEED: u64 = 42;

const TABLE_SUFFIX: &str = ".quant.parquet";
const MS1_SUFFIX: &str = ".ms1.parquet";
const MS2_SUFFIX: &str = ".ms2.parquet";
const EXTRACT_SUFFIX: &str = ".extract.parquet";
/// Extract-table metadata key of the template-fit m/z tolerance, ppm.
const TOL_KEY: &str = "andes.extract.tol_ppm";
/// Extract table metadata key of the transfer selection (absent when every
/// transfer pair was measured).
const SELECTION_KEY: &str = "andes.extract.transfer_selection";
/// Extract table tag: the run's MS1 scan spacing and elution-width scale
/// ([`crate::quant_mbr::RunScales::describe`]).
const RUN_SCALES_KEY: &str = "andes.extract.run_scales";
/// Isotopes per precursor envelope in the run tables and the library.
const ENVELOPE_LEN: usize = 4;

fn artifact(dir: &Path, stem: &str, suffix: &str) -> PathBuf {
    dir.join(format!("{stem}{suffix}"))
}

/// The run stems of the library, sorted: `runs` when given, else every
/// `<stem>.quant.parquet` in `dir`. Errors: no runs, a repeated stem, or a
/// run whose PSM, MS1 or MS2 table is missing (every such file is named).
pub(crate) fn library_runs(dir: &Path, runs: &[String]) -> Result<Vec<String>, String> {
    let mut stems: Vec<String> = if runs.is_empty() {
        let entries = std::fs::read_dir(dir).map_err(|e| format!("{}: {e}", dir.display()))?;
        let mut v = Vec::new();
        for e in entries {
            let name = e
                .map_err(|e| format!("{}: {e}", dir.display()))?
                .file_name()
                .to_string_lossy()
                .into_owned();
            if let Some(stem) = name.strip_suffix(TABLE_SUFFIX) {
                v.push(stem.to_string());
            }
        }
        v
    } else {
        runs.to_vec()
    };
    stems.sort();
    if stems.is_empty() {
        return Err(format!(
            "no runs: {} holds no *{TABLE_SUFFIX} and no --runs were given",
            dir.display()
        ));
    }
    if let Some(w) = stems.windows(2).find(|w| w[0] == w[1]) {
        return Err(format!(
            "run {:?} is given twice; run stems must be unique",
            w[0]
        ));
    }
    if stems.len() > u16::MAX as usize {
        return Err(format!("{} runs exceed the supported number", stems.len()));
    }
    let missing: Vec<String> = stems
        .iter()
        .flat_map(|s| {
            [
                artifact(dir, s, TABLE_SUFFIX),
                artifact(dir, s, MS1_SUFFIX),
                artifact(dir, s, MS2_SUFFIX),
            ]
        })
        .filter(|p| !p.is_file())
        .map(|p| p.display().to_string())
        .collect();
    if !missing.is_empty() {
        return Err(format!(
            "missing per-run artifacts (written by `andes search --quant-dir`): {}",
            missing.join(", ")
        ));
    }
    Ok(stems)
}

/// Hashes SpecIds and residue sequences under one random key (SipHash, the
/// standard library's keyed hash). A SpecId becomes 128 bits: two 64-bit
/// hashes of it, under distinct one-byte prefixes. Two distinct SpecIds of a
/// set of `n` share a key with probability about `n^2 / 2^129` (below
/// `1e-22` for 100 million SpecIds), so a map keyed by the hash joins as a
/// map keyed by the SpecId does, without holding the SpecId strings.
#[derive(Clone, Default)]
struct SpecIdKeys(std::collections::hash_map::RandomState);

impl SpecIdKeys {
    fn hash(&self, tag: u8, bytes: &[u8]) -> u64 {
        use std::hash::{BuildHasher, Hasher};
        let mut h = self.0.build_hasher();
        h.write_u8(tag);
        h.write(bytes);
        h.finish()
    }

    /// The 128-bit key of a SpecId.
    fn key(&self, spec_id: &str) -> u128 {
        (u128::from(self.hash(0, spec_id.as_bytes())) << 64)
            | u128::from(self.hash(1, spec_id.as_bytes()))
    }

    /// A 64-bit hash of a residue sequence.
    fn residues(&self, residues: &[u8]) -> u64 {
        self.hash(2, residues)
    }
}

/// One rescored PSM as the library joins it: the rescorer's score, PSM
/// q-value and PEP, the label (from the decoy result file or not) and the
/// positional-isomer class of the rescorer's peptide string, packed in `tag`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Scored {
    pub score: f64,
    /// The rescorer's PSM q-value.
    pub q: f64,
    /// The rescorer's PSM posterior error probability.
    pub pep: f64,
    /// From the high bit down: [`RESIDUE_BITS`] bits of the residue hash,
    /// [`MOD_SET_BITS`] bits of the modification multiset's index in
    /// [`ScoreSet::mod_sets`], and the decoy label.
    tag: u64,
}

/// Bits of [`Scored::tag`] that hold the residue hash: two different residue
/// sequences pass the join check together with probability 2^-40.
const RESIDUE_BITS: u32 = 40;
/// Bits of [`Scored::tag`] that hold the modification multiset index.
const MOD_SET_BITS: u32 = 64 - RESIDUE_BITS - 1;

impl Scored {
    /// From the decoy result file.
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

/// Rescored PSMs by SpecId and the file(s) they came from. A SpecId is held
/// as its 128-bit hash ([`SpecIdKeys`]): 48 bytes per PSM in the map.
pub(crate) struct ScoreSet {
    pub source: String,
    keys: SpecIdKeys,
    psms: HashMap<u128, Scored, rustc_hash::FxBuildHasher>,
    /// The distinct modification mass multisets of the PSMs' peptides
    /// ([`IsomerClass::mod_units`]), and the index of each.
    mod_sets: Vec<Box<[i64]>>,
    mod_set_ids: HashMap<Box<[i64]>, u32, rustc_hash::FxBuildHasher>,
    /// Scratch space of the peptide parse.
    class: IsomerClass,
}

impl ScoreSet {
    fn with_capacity(source: String, n: usize) -> Self {
        Self {
            source,
            keys: SpecIdKeys::default(),
            psms: HashMap::with_capacity_and_hasher(n, rustc_hash::FxBuildHasher),
            mod_sets: Vec::new(),
            mod_set_ids: HashMap::default(),
            class: IsomerClass::default(),
        }
    }

    /// Adds a PSM. Errors (the message names [`Self::source`]): a PSM without
    /// a score, a peptide whose modifications have no mass in the bundled
    /// Unimod table, or a SpecId already added; `file` names where the PSM is
    /// read from, for a SpecId repeated within one label or an unreadable
    /// peptide.
    #[allow(clippy::too_many_arguments)]
    fn insert(
        &mut self,
        file: &str,
        spec_id: &str,
        score: Option<f64>,
        q: f64,
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
        let scored = Scored { score, q, pep, tag };
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

    /// The rescored PSM of a SpecId.
    pub(crate) fn get(&self, spec_id: &str) -> Option<&Scored> {
        self.psms.get(&self.keys.key(spec_id))
    }

    #[cfg(test)]
    fn len(&self) -> usize {
        self.psms.len()
    }

    /// Whether `s` has label `is_decoy` and the positional-isomer class
    /// `class` ([`IsomerClass::matches`]): the same residues and modification
    /// masses, on any sites.
    fn agrees(&self, s: &Scored, is_decoy: bool, class: &IsomerClass) -> bool {
        s.is_decoy() == is_decoy
            && s.residues() == self.keys.residues(class.residues()) >> (64 - RESIDUE_BITS)
            && mod_units_match(&self.mod_sets[s.mod_set()], class.mod_units())
    }
}

/// The target and decoy PSMs of one rescoring result. Errors: a PSM without
/// a finite score or a SpecId among both (the message names `source`).
pub(crate) fn rescore_scores(
    targets: &HashMap<String, PercolatorPsm>,
    decoys: &HashMap<String, PercolatorPsm>,
    source: &str,
) -> Result<ScoreSet, String> {
    let mut set = ScoreSet::with_capacity(source.to_string(), targets.len() + decoys.len());
    for (map, is_decoy) in [(targets, false), (decoys, true)] {
        for (id, p) in map {
            set.insert(source, id, p.score, p.q_value, p.pep, &p.peptide, is_decoy)?;
        }
    }
    Ok(set)
}

/// Largest share of a run's joined PSMs that may disagree with the rescoring
/// result on the peptide's positional-isomer class (residues and modification
/// masses) or the target/decoy label.
const MAX_JOIN_MISMATCH: f64 = 0.01;

/// Where the rescorer's scores come from: one set for every run (SpecIds
/// must then be distinct across runs), loaded once; or one set per run stem,
/// loaded when that run is joined and dropped after it, so only one run's
/// scores are held at a time.
pub(crate) enum ScoreLoader<'a> {
    /// Called with the number of distinct SpecIds of the run tables, to size
    /// the set.
    Merged(Box<dyn FnOnce(usize) -> Result<ScoreSet, String> + 'a>),
    PerRun(RunScores<'a>),
}

/// The score set of one run stem.
pub(crate) type RunScores<'a> = Box<dyn Fn(&str) -> Result<ScoreSet, String> + 'a>;

/// Where the Percolator results of `quant library` come from.
#[derive(Debug, PartialEq)]
enum ScoreSource {
    /// One target and one decoy file over every run.
    Merged(PathBuf, PathBuf),
    /// Per run stem, its target and decoy file.
    PerRun(BTreeMap<String, (PathBuf, PathBuf)>),
}

/// Reads `--psms` / `--decoy-psms`: either one plain path each, or
/// `<stem>=<path>` once per run of `stems` in each.
fn score_source(
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

/// The loader of `quant library`'s Percolator results.
fn score_loader(source: &ScoreSource) -> ScoreLoader<'_> {
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

/// Scores from Percolator's target and decoy PSM files, read one line at a
/// time; `capacity` sizes the map (the expected number of PSMs).
fn percolator_scores(targets: &Path, decoys: &Path, capacity: usize) -> Result<ScoreSet, String> {
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
                row.pep,
                row.peptide,
                is_decoy,
            )?;
        }
    }
    Ok(set)
}

/// The columns of a run table that give a precursor's metadata.
const META_COLUMNS: [&str; 6] = [
    "precursor_key",
    "is_decoy",
    "mz",
    "charge",
    "envelope",
    "proteins",
];
/// The columns of a run table that join its PSMs to the rescoring result.
const JOIN_COLUMNS: [&str; 7] = [
    "spec_id",
    "precursor_key",
    "peptidoform",
    "is_decoy",
    "rt_s",
    "faims_cv",
    "ion_mobility",
];

/// A precursor key's metadata row: the key's first target row, else its first
/// row, and the run it came from.
struct KeyMeta {
    has_target: bool,
    run: usize,
    mz: f64,
    charge: u8,
    /// The isotope fractions, or the number of them when it is not
    /// [`ENVELOPE_LEN`].
    envelope: Result<[f32; ENVELOPE_LEN], usize>,
    proteins: String,
}

fn envelope_of(v: &[f32]) -> Result<[f32; ENVELOPE_LEN], usize> {
    v.try_into().map_err(|_| v.len())
}

/// Alignment settings with each run's residual-width floor: half the median
/// spacing of its MS1 scans of one FAIMS voltage, read from
/// `<stem>.ms1.parquet`. A run without a readable table, or with fewer than two
/// MS1 scans, keeps the default floor.
fn align_params(dir: &Path, stems: &[String]) -> AlignParams {
    let mut ap = AlignParams::default();
    ap.run_min_spacing_s = stems
        .iter()
        .map(|stem| {
            let path = artifact(dir, stem, MS1_SUFFIX);
            let spacing = output::ms1_table::read_ms1_scan_spacing(&path)
                .ok()
                .flatten();
            let floor = spacing.map_or(ap.min_spacing_s, |s| s / 2.0);
            match spacing {
                Some(s) => eprintln!(
                    "[mbr] run {stem}: MS1 scan spacing {s:.2} s, RT residual-width floor {floor:.2} s"
                ),
                None => eprintln!(
                    "[mbr] run {stem}: MS1 scan spacing unknown, RT residual-width floor {floor:.2} s"
                ),
            }
            floor
        })
        .collect();
    ap
}

/// Log, per aligned run, how many residual-width knots sit at the run's floor:
/// where the floor binds, the MS1 sampling and not the observed scatter sets
/// the transfer RT uncertainty.
fn log_floor_hits(library: &Library, ap: &AlignParams, stems: &[String]) {
    for (r, knots) in library.alignment.residual_sd.iter().enumerate() {
        let Some(knots) = knots else { continue };
        let floor = ap.floor_s(r);
        let at_floor = knots.iter().filter(|k| k.1 <= floor * (1.0 + 1e-9)).count();
        eprintln!(
            "[mbr] run {}: {at_floor} of {} RT residual-width knots at the {floor:.2} s floor",
            stems.get(r).map_or("?", String::as_str),
            knots.len()
        );
    }
}

/// Builds the library over the run tables of `stems` (sorted, as from
/// [`library_runs`]) with the rescorer's scores joined on (run, SpecId), and
/// writes it to `dir`. `max_q` is the entry gate, `donor_q` the PSM q-value
/// and `donor_pep` the optional PSM PEP cap of the donor PSMs. A PSM whose
/// rescoring result names another positional-isomer class (residues or
/// modification masses) or target/decoy label is skipped.
///
/// A precursor is a positional-isomer class with a charge: the precursor keys
/// of its site forms join it. Its key is the site form of its best PSM: a donor
/// PSM when it has one, then the highest score, earliest RT and lowest run, as
/// [`DonorPsms`] picks the PSM whose scan gives the library spectrum.
///
/// The run tables are read twice, one run and one record batch at a time:
/// first the precursor keys and their metadata (the union of their classes
/// numbers the precursors), then the PSMs, joined to their run's scores and folded into
/// per-precursor aggregates ([`quant::library::LibraryBuilder`] and the
/// donor scans). Memory follows the precursors and their donor runs, plus one
/// run's scores; with merged scores, the one score set (48 bytes per PSM: a
/// SpecId is held as its 128-bit hash, see [`ScoreSet`]) and, before it is
/// read, every SpecId's hash and first run, to name a SpecId shared by two
/// runs. The Percolator files are read one line at a time. Merged scores are
/// loaded after the tables' SpecIds are checked, so that a SpecId shared by
/// two runs fails with the runs named.
pub(crate) fn build_and_write_library(
    dir: &Path,
    stems: &[String],
    scores: ScoreLoader<'_>,
    max_q: f64,
    donor_q: f64,
    donor_pep: Option<f64>,
) -> Result<StoredLibrary, String> {
    if !(max_q > 0.0 && max_q <= 1.0) {
        return Err(format!("library entry gate {max_q} is outside (0, 1]"));
    }
    if !(donor_q > 0.0 && donor_q <= 1.0) {
        return Err(format!(
            "donor PSM q-value threshold {donor_q} is outside (0, 1]"
        ));
    }
    if let Some(c) = donor_pep.filter(|c| !(0.0..=1.0).contains(c)) {
        return Err(format!("donor PSM PEP cap {c} is outside [0, 1]"));
    }
    let paths: Vec<PathBuf> = stems
        .iter()
        .map(|s| artifact(dir, s, TABLE_SUFFIX))
        .collect();
    let open = |r: usize, columns: &[&str]| {
        let batches = open_run_quant_table(&paths[r], columns)?;
        if batches.run() != stems[r] {
            return Err(format!(
                "{}: the table is of run {:?}, not {:?}",
                paths[r].display(),
                batches.run(),
                stems[r]
            ));
        }
        Ok(batches)
    };
    let merged = matches!(scores, ScoreLoader::Merged(_));
    // Pass 1: every precursor key with its metadata row, and (merged scores
    // only) the first SpecId shared by two runs.
    let mut keys: BTreeMap<String, KeyMeta> = BTreeMap::new();
    // SpecId (as its 128-bit hash, see `SpecIdKeys`) -> first run.
    let spec_keys = SpecIdKeys::default();
    let mut seen: HashMap<u128, u32, rustc_hash::FxBuildHasher> = HashMap::default();
    let mut shared: Option<(String, usize, usize)> = None;
    let mut n_rows = 0usize;
    for r in 0..stems.len() {
        let columns: Vec<&str> = META_COLUMNS
            .iter()
            .copied()
            .chain(merged.then_some("spec_id"))
            .collect();
        for batch in open(r, &columns)? {
            for row in batch? {
                n_rows += 1;
                if merged {
                    let other = *seen.entry(spec_keys.key(&row.spec_id)).or_insert(r as u32);
                    if shared.is_none() && other as usize != r {
                        shared = Some((row.spec_id.clone(), other as usize, r));
                    }
                }
                let target = !row.is_decoy;
                match keys.get_mut(row.precursor_key.as_str()) {
                    Some(m) => {
                        if target && !m.has_target {
                            *m = KeyMeta {
                                has_target: true,
                                run: r,
                                mz: row.mz,
                                charge: row.charge,
                                envelope: envelope_of(&row.envelope),
                                proteins: row.proteins,
                            };
                        }
                    }
                    None => {
                        keys.insert(
                            row.precursor_key,
                            KeyMeta {
                                has_target: target,
                                run: r,
                                mz: row.mz,
                                charge: row.charge,
                                envelope: envelope_of(&row.envelope),
                                proteins: row.proteins,
                            },
                        );
                    }
                }
            }
        }
    }
    let n_spec_ids = seen.len();
    drop(seen);
    let (merged_scores, per_run) = match scores {
        ScoreLoader::Merged(load) => {
            let set = load(n_spec_ids)?;
            if let Some((id, a, b)) = &shared {
                return Err(format!(
                    "SpecId {id:?} is in runs {} and {}, so one merged Percolator result cannot be joined; pass per-run results: --psms <stem>=<target.tsv> --decoy-psms <stem>=<decoy.tsv> for every run",
                    stems[*a], stems[*b]
                ));
            }
            (Some(set), None)
        }
        ScoreLoader::PerRun(load) => (None, Some(load)),
    };
    // Every run's modification masses by name, to read a peptidoform's class.
    let mut mod_masses: HashMap<String, f64> = HashMap::new();
    for p in &paths {
        mod_masses.extend(read_run_quant_mods(p)?);
    }
    let resolve = |tag: &str| {
        mod_masses
            .get(tag)
            .copied()
            .or_else(|| bundled_unimod_mass(tag))
    };
    // The precursor keys of each positional-isomer class and charge.
    let site_keys: Vec<(String, KeyMeta)> = std::mem::take(&mut keys).into_iter().collect();
    let mut classes: BTreeMap<(IsomerClass, u8), Vec<u32>> = BTreeMap::new();
    for (i, (key, m)) in site_keys.iter().enumerate() {
        let class = IsomerClass::parse(split_key(key), resolve)
            .map_err(|e| format!("{}: precursor {key:?}: {e}", paths[m.run].display()))?;
        classes.entry((class, m.charge)).or_default().push(i as u32);
    }
    // Dense precursor ids over the sorted classes; the first target row of
    // the class's first site form with one (else of its first site form) gives
    // its metadata. A decoy row of a class that also has target rows is
    // dropped: a precursor is all-target or all-decoy.
    let mut meta = Vec::with_capacity(classes.len());
    let mut proteins = Vec::with_capacity(classes.len());
    let mut target_key = Vec::with_capacity(classes.len());
    let mut forms: Vec<Vec<u32>> = Vec::with_capacity(classes.len());
    for (_, sites) in std::mem::take(&mut classes) {
        let first = *sites
            .iter()
            .find(|&&i| site_keys[i as usize].1.has_target)
            .unwrap_or(&sites[0]);
        let (key, m) = &site_keys[first as usize];
        let envelope = m.envelope.map_err(|len| {
            format!(
                "{}: precursor {key:?} has an envelope of {len} isotopes, expected {ENVELOPE_LEN}",
                paths[m.run].display()
            )
        })?;
        meta.push(PrecursorMeta {
            key: key.clone(),
            mz: m.mz,
            charge: m.charge,
            envelope,
        });
        proteins.push(m.proteins.clone());
        target_key.push(m.has_target);
        forms.push(sites);
    }
    // Precursor key -> (precursor id, site form index in `site_keys`).
    let site_keys_ref = &site_keys;
    let id_of: HashMap<&str, (u32, u32)> = forms
        .iter()
        .enumerate()
        .flat_map(|(id, sites)| {
            sites
                .iter()
                .map(move |&i| (site_keys_ref[i as usize].0.as_str(), (id as u32, i)))
        })
        .collect();
    // Pass 2: the PSMs, joined to their run's scores.
    let mut builder = LibraryBuilder::new(meta.len(), stems.len(), donor_q, donor_pep);
    let mut donors = DonorPsms::new(meta.len());
    let mut best = vec![None::<BestForm>; meta.len()];
    let (mut unscored_t, mut unscored_d, mut shadowed, mut mismatched) = (0usize, 0, 0, 0);
    let (mut n_psms, mut donor_psms, mut n_joined) = (0usize, 0usize, 0usize);
    let mut class = IsomerClass::default();
    for r in 0..stems.len() {
        let own;
        let set = match (&merged_scores, &per_run) {
            (Some(m), _) => m,
            (None, Some(load)) => {
                own = load(&stems[r])?;
                &own
            }
            (None, None) => unreachable!("one score source"),
        };
        let (mut run_rows, mut run_joined, mut run_mismatched) = (0usize, 0usize, 0usize);
        for batch in open(r, &JOIN_COLUMNS)? {
            for row in batch? {
                run_rows += 1;
                let (id, site) = *id_of.get(row.precursor_key.as_str()).ok_or_else(|| {
                    format!(
                        "{}: precursor {:?} changed between two reads of the table",
                        paths[r].display(),
                        row.precursor_key
                    )
                })?;
                if row.is_decoy && target_key[id as usize] {
                    shadowed += 1;
                    continue;
                }
                let Some(&j) = set.get(&row.spec_id) else {
                    if row.is_decoy {
                        unscored_d += 1;
                    } else {
                        unscored_t += 1;
                    }
                    continue;
                };
                run_joined += 1;
                class
                    .parse_into(&row.peptidoform, resolve)
                    .map_err(|e| format!("{}: {e}", paths[r].display()))?;
                if !set.agrees(&j, row.is_decoy, &class) {
                    run_mismatched += 1;
                    continue;
                }
                let p = PsmIn {
                    run: r as u16,
                    precursor: id,
                    rt: row.rt_s,
                    score: j.score,
                    q: j.q,
                    pep: j.pep,
                    is_decoy: row.is_decoy,
                    faims_cv: row.faims_cv,
                    ion_mobility: row.ion_mobility,
                };
                n_psms += 1;
                let donor = builder.is_donor(&p);
                donor_psms += usize::from(donor && !p.is_decoy);
                builder.add(&p);
                let form = BestForm {
                    donor,
                    score: p.score,
                    rt: p.rt,
                    run: p.run,
                    site,
                };
                let slot = &mut best[id as usize];
                if slot.as_ref().is_none_or(|b| form.beats(b)) {
                    *slot = Some(form);
                }
                if donor {
                    donors.add(&p, &row.spec_id);
                }
            }
        }
        // Per-run results whose SpecIds all belong to another run (every SpecId
        // starts with its run's file stem) would join nothing and pass silently.
        if per_run.is_some() && run_joined == 0 && run_rows > 0 && !set.psms.is_empty() {
            return Err(format!(
                "run {}: none of the {} rescored PSMs in {} has a SpecId of this run; are these the rescoring results of this run?",
                stems[r],
                set.psms.len(),
                set.source
            ));
        }
        if run_mismatched as f64 > MAX_JOIN_MISMATCH * run_joined as f64 {
            return Err(format!(
                "run {}: {run_mismatched} of {run_joined} PSMs joined on SpecId disagree with {} on the peptide, its modifications or the target/decoy label; are these the rescoring results of this run?",
                stems[r], set.source
            ));
        }
        mismatched += run_mismatched;
        n_joined += run_joined;
    }
    if n_psms == 0 && n_rows > 0 {
        return Err(format!(
            "none of the {n_rows} PSMs of the run tables has a rescoring score; \
             the PSM files must come from the PINs of the same searches"
        ));
    }
    drop(id_of);
    drop(merged_scores);
    // Each precursor's key: the site form of its best PSM.
    let (mut multi, mut n_forms) = (0usize, 0usize);
    for (id, b) in best.iter().enumerate() {
        if let Some(b) = b {
            meta[id].key = site_keys[b.site as usize].0.clone();
        }
        if target_key[id] {
            n_forms += forms[id].len();
            multi += usize::from(forms[id].len() > 1);
        }
    }
    drop(forms);
    drop(site_keys);
    eprintln!(
        "[mbr] positional isomers: {} target precursors hold {n_forms} site forms; {multi} precursors join several site forms and are keyed by the site form of their best PSM",
        target_key.iter().filter(|&&t| t).count()
    );
    let donor_gate = match donor_pep {
        Some(c) => format!("PSM q<={donor_q} and PEP<={c}"),
        None => format!("PSM q<={donor_q}"),
    };
    let ap = align_params(dir, stems);
    let library = builder.finish(meta, max_q, &ap);
    log_floor_hits(&library, &ap, stems);
    let (spectra, donor_scans, sc) = donors.finish(dir, stems, &library)?;
    eprintln!(
        "[mbr] library: join check: {mismatched} of {n_joined} PSMs joined on SpecId ({:.4}%) disagree with the rescoring result on the peptide's residues, modification masses or target/decoy label (a run fails above {}%)",
        100.0 * mismatched as f64 / n_joined.max(1) as f64,
        100.0 * MAX_JOIN_MISMATCH
    );
    eprintln!(
        "[mbr] library: {} runs, {} PSMs scored, {donor_psms} target PSMs at {donor_gate} identify their precursor in their run; skipped: {} PSMs whose SpecId joins a rescoring PSM of another peptide, modification composition or target/decoy label, {} target and {} decoy PSMs without a rescoring score, {} decoy PSMs sharing a target's key; {} precursors, {} entries at the internal entry gate q<={max_q} (precursor-level, not an identification FDR) with a PSM at {donor_gate}, {} of them identity decoys (decoy precursors with decoy PSMs at {donor_gate}, measured as targets and never reported)",
        stems.len(),
        n_psms,
        mismatched,
        unscored_t,
        unscored_d,
        shadowed,
        library.meta.len(),
        library.entries.len(),
        library.entries.iter().filter(|e| e.identity_decoy).count()
    );
    eprintln!(
        "[mbr] library: fragment spectra (best donor PSM, at most {} b/y fragments of {}+ residues matched on its scan) for {} of {} entries; none for {} whose scan is not in the MS2 table, {} without fragment masses, {} with fewer than {} matched fragments",
        quant::ms2::LIBRARY_FRAGMENTS,
        quant::ms2::MIN_ION_LENGTH,
        sc.with_spectrum,
        library.entries.len(),
        sc.no_scan,
        sc.no_masses,
        sc.too_few,
        quant::ms2::MIN_LIBRARY_FRAGMENTS
    );
    let a = &library.alignment;
    let no_consensus = library
        .entries
        .iter()
        .filter(|e| !e.consensus_rt.is_finite())
        .count();
    if no_consensus > 0 {
        eprintln!(
            "[mbr] library: {no_consensus} entries have no consensus RT (no observation in an aligned run) and are not transferred"
        );
    }
    if a.n_components() > 1 {
        let mut groups: Vec<Vec<&str>> = vec![Vec::new(); a.n_components()];
        for (r, stem) in stems.iter().enumerate() {
            if let Some(c) = a.component(r as u16) {
                groups[c as usize].push(stem);
            }
        }
        eprintln!(
            "[mbr] alignment: the aligned runs form {} components (runs joined when they share at least {} precursors); transfers stay inside a component: {}",
            groups.len(),
            ap.min_shared,
            groups
                .iter()
                .enumerate()
                .map(|(c, g)| format!("component {c}: {}", g.join(", ")))
                .collect::<Vec<_>>()
                .join("; ")
        );
    }
    for (r, stem) in stems.iter().enumerate() {
        if a.maps[r].is_none() {
            let isolated = a.anchors[r] >= ap.min_anchors && a.coverage[r] >= ap.min_coverage;
            eprintln!(
                "[mbr] run {stem} is unaligned ({} anchors, {:.2} of the consensus RT range{}): it neither gives nor receives transfers and keeps its identified quantities",
                a.anchors[r],
                a.coverage[r],
                if isolated {
                    format!(
                        "; shares fewer than {} precursors with every other aligned run",
                        ap.min_shared
                    )
                } else {
                    String::new()
                }
            );
        }
    }
    let stored = StoredLibrary {
        runs: stems.to_vec(),
        max_q,
        donor_q,
        donor_pep,
        library,
        proteins,
        spectra,
        donor_scans,
    };
    write_library(dir, &stored)?;
    eprintln!(
        "[mbr] wrote {} and {}",
        dir.join(output::library_table::LIBRARY_FILE).display(),
        dir.join(output::library_table::ALIGNMENT_FILE).display()
    );
    Ok(stored)
}

/// What [`entry_spectra`] gave, for the log.
#[derive(Debug, Default, PartialEq)]
pub(crate) struct SpectrumCounts {
    pub with_spectrum: usize,
    /// The best donor PSM's scan is not in its run's MS2 table.
    pub no_scan: usize,
    /// The peptidoform holds a residue or modification without a mass.
    pub no_masses: usize,
    /// Fewer than [`quant::ms2::MIN_LIBRARY_FRAGMENTS`] fragments matched.
    pub too_few: usize,
}

/// Per precursor id: the library spectrum and the donor PSMs' `(run, RT)`;
/// and the counts.
pub(crate) type EntrySpectra = (
    Vec<Option<LibrarySpectrum>>,
    Vec<Vec<(u16, f64)>>,
    SpectrumCounts,
);

/// The MS2 table scan of a PSM SpecId (`<spec_key>_<rank>` or
/// `<spec_key>_<rank>_<row>`).
fn scan_of(by_key: &HashMap<&str, usize>, spec_id: &str) -> Option<usize> {
    let mut id = spec_id;
    for _ in 0..2 {
        id = id.rsplit_once('_')?.0;
        if let Some(&i) = by_key.get(id) {
            return Some(i);
        }
    }
    None
}

/// The best PSM of a precursor so far and its site form, which keys it: a
/// donor PSM beats any other PSM, then the higher score, earlier RT and lower
/// run, as [`DonorPsms`] picks the best donor.
#[derive(Clone, Copy)]
struct BestForm {
    donor: bool,
    score: f64,
    rt: f64,
    run: u16,
    /// Index of the site form's precursor key.
    site: u32,
}

impl BestForm {
    fn beats(&self, other: &Self) -> bool {
        if self.donor != other.donor {
            return self.donor;
        }
        self.score > other.score
            || (self.score == other.score
                && (self.rt.total_cmp(&other.rt).then(self.run.cmp(&other.run))).is_lt())
    }
}

/// The donor PSMs of every precursor as the library's spectra need them,
/// added one at a time in run order: per precursor id with a donor PSM, the
/// `(run, RT)` of each donor PSM and the best one (highest score; then
/// earliest RT, lowest run) with its SpecId.
pub(crate) struct DonorPsms {
    /// Per precursor id: its index in `scans` and `best`, or `u32::MAX`.
    slot: Vec<u32>,
    scans: Vec<Vec<(u16, f64)>>,
    best: Vec<(BestDonor, Box<str>)>,
}

/// What [`DonorPsms`] compares of a donor PSM.
#[derive(Clone, Copy)]
struct BestDonor {
    score: f64,
    rt: f64,
    run: u16,
}

impl DonorPsms {
    pub(crate) fn new(n_precursors: usize) -> Self {
        Self {
            slot: vec![u32::MAX; n_precursors],
            scans: Vec::new(),
            best: Vec::new(),
        }
    }

    /// Adds donor PSM `p` with SpecId `spec_id`; ignored when its id is out of
    /// range.
    pub(crate) fn add(&mut self, p: &PsmIn, spec_id: &str) {
        let id = p.precursor as usize;
        if id >= self.slot.len() {
            return;
        }
        let new = BestDonor {
            score: p.score,
            rt: p.rt,
            run: p.run,
        };
        let i = match self.slot[id] {
            u32::MAX => {
                self.slot[id] = u32::try_from(self.scans.len())
                    .ok()
                    .filter(|&i| i != u32::MAX)
                    .expect("fewer than 2^32 - 1 donor precursors");
                self.scans.push(Vec::new());
                self.best.push((new, spec_id.into()));
                self.scans.len() - 1
            }
            i => {
                let i = i as usize;
                let (b, s) = &mut self.best[i];
                if p.score > b.score
                    || (p.score == b.score
                        && (p.rt.total_cmp(&b.rt).then(p.run.cmp(&b.run))).is_lt())
                {
                    *b = new;
                    *s = spec_id.into();
                }
                i
            }
        };
        self.scans[i].push((p.run, p.rt));
    }

    /// Per precursor id of `lib`: the empirical fragment spectrum of an
    /// entry's best donor PSM, matched on its scan in the run's MS2 table at
    /// the table's fragment tolerance, and the `(run, RT)` of every donor PSM
    /// of the entry, sorted; empty for a precursor that is not an entry.
    /// The MS2 tables are read one run at a time. Errors: an unreadable MS2
    /// or PSM table (the message names it).
    pub(crate) fn finish(
        self,
        dir: &Path,
        stems: &[String],
        lib: &Library,
    ) -> Result<EntrySpectra, String> {
        let n = lib.meta.len();
        let Self {
            slot,
            mut scans,
            best: donor_best,
        } = self;
        // Per precursor id of an entry: its donor scans, sorted, and its best
        // donor PSM's run and SpecId.
        let mut entry_scans: Vec<Vec<(u16, f64)>> = vec![Vec::new(); n];
        let mut best: Vec<Option<(u16, Box<str>)>> = vec![None; n];
        let mut donor_best: Vec<Option<(BestDonor, Box<str>)>> =
            donor_best.into_iter().map(Some).collect();
        for e in &lib.entries {
            let id = e.precursor as usize;
            let Some(&i) = slot.get(id).filter(|&&i| i != u32::MAX) else {
                continue;
            };
            let mut d = std::mem::take(&mut scans[i as usize]);
            d.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.total_cmp(&b.1)));
            entry_scans[id] = d;
            best[id] = donor_best[i as usize].take().map(|(b, s)| (b.run, s));
        }
        drop(scans);
        drop(donor_best);
        drop(slot);
        let mut spectra: Vec<Option<LibrarySpectrum>> = vec![None; n];
        let mut c = SpectrumCounts::default();
        for (r, stem) in stems.iter().enumerate() {
            let wanted: Vec<(usize, &str)> = best
                .iter()
                .enumerate()
                .filter_map(|(id, b)| {
                    b.as_ref()
                        .filter(|(run, _)| *run as usize == r)
                        .map(|(_, s)| (id, &**s))
                })
                .collect();
            if wanted.is_empty() {
                continue;
            }
            let ms2_path = artifact(dir, stem, MS2_SUFFIX);
            let table = read_ms2_table(&ms2_path)?;
            if table.run != *stem {
                return Err(format!(
                    "{}: the table is of run {:?}, not {stem:?}",
                    ms2_path.display(),
                    table.run
                ));
            }
            let mods = read_run_quant_mods(&artifact(dir, stem, TABLE_SUFFIX))?;
            let by_key: HashMap<&str, usize> = table
                .scans
                .iter()
                .enumerate()
                .map(|(i, s)| (s.spec_key.as_str(), i))
                .collect();
            let made: Vec<(usize, Result<LibrarySpectrum, u8>)> = wanted
                .par_iter()
                .map(|&(id, spec_id)| {
                    let m = &lib.meta[id];
                    let Some(scan) = scan_of(&by_key, spec_id) else {
                        return (id, Err(0));
                    };
                    let Ok(masses) = residue_masses(split_key(&m.key), &mods) else {
                        return (id, Err(1));
                    };
                    let ladder = fragment_ladder(&masses, m.charge.max(1));
                    let peaks = &table.scans[scan].peaks;
                    (
                        id,
                        library_spectrum(&ladder, peaks, table.fragment_tol).ok_or(2),
                    )
                })
                .collect();
            for (id, s) in made {
                match s {
                    Ok(s) => {
                        c.with_spectrum += 1;
                        spectra[id] = Some(s);
                    }
                    Err(0) => c.no_scan += 1,
                    Err(1) => c.no_masses += 1,
                    Err(_) => c.too_few += 1,
                }
            }
        }
        Ok((spectra, entry_scans, c))
    }
}

fn log_extract(stem: &str, rows: usize, c: &ExtractCounts) {
    eprintln!(
        "[mbr] extract {stem}: {rows} rows{}; skipped: {} entries at a FAIMS voltage the run lacks, {} identified pairs with a non-finite best PSM RT, {} identified pairs without a prediction from other runs (no pseudo-transfer), {} transfers without a prediction, {} pseudo-transfers and {} transfers whose donor runs are all in other alignment components, {} transfers left out by the transfer selection; MS1 scan spacing {}, width scale {:.3} s ({} fits, {} rounds{}); width prior from {} fits{}",
        if c.unaligned {
            " (unaligned run: identified rows only)"
        } else {
            ""
        },
        c.foreign_cv,
        c.nonfinite_best_rt,
        c.no_pseudo_prediction,
        c.no_prediction,
        c.pseudo_other_component,
        c.other_component,
        c.not_selected,
        c.scales
            .spacing_s
            .map_or_else(|| "unknown".to_string(), |s| format!("{s:.3} s")),
        c.scales.width_s,
        c.scales.width_fits,
        c.scales.rounds,
        if c.scales.fallback {
            "; too few fits: the seed width is used"
        } else {
            ""
        },
        c.width_fits,
        if c.width_prior_fallback {
            " (too few: the width scale is used)"
        } else {
            ""
        }
    );
}

fn log_ms2(stem: &str, c: &Ms2Counts, secs: f64) {
    let part = |role: usize| {
        ["target", "mass-shift decoy", "RT-shift decoy"]
            .iter()
            .enumerate()
            .filter(|&(f, _)| c.rows[role][f] > 0)
            .map(|(f, name)| {
                let n = c.rows[role][f] as f64;
                format!(
                    "{name} {} (scan in window {:.3}, supported {:.4})",
                    c.rows[role][f],
                    c.with_scan[role][f] as f64 / n,
                    c.supported[role][f] as f64 / n
                )
            })
            .collect::<Vec<_>>()
            .join(", ")
    };
    eprintln!(
        "[mbr] extract {stem}: MS2 evidence ({secs:.1} s): transfers: {}; pseudo-transfers: {}; {} rows of entries without a library spectrum",
        part(0),
        part(1),
        c.no_spectrum
    );
}

/// Measures run `run` of the library on its MS1 table (the transfers that
/// `selection` keeps), checks its transfer and pseudo-transfer rows against
/// its MS2 table, and writes `<stem>.extract.parquet`.
pub(crate) fn extract_and_write(
    dir: &Path,
    stored: &StoredLibrary,
    fingerprint: &str,
    run: u16,
    ep: &EvidenceParams,
    selection: &TransferSelection,
) -> Result<Vec<ExtractRow>, String> {
    let stem = &stored.runs[run as usize];
    let ms1 = artifact(dir, stem, MS1_SUFFIX);
    let mut ex = extract_run(run, &ms1, &stored.library, ep, selection)?;
    log_extract(stem, ex.rows.len(), &ex.counts);
    let t = std::time::Instant::now();
    let ms2_path = artifact(dir, stem, MS2_SUFFIX);
    let table = read_ms2_table(&ms2_path)?;
    if table.run != *stem {
        return Err(format!(
            "{}: the table is of run {:?}, not {stem:?}",
            ms2_path.display(),
            table.run
        ));
    }
    let tol = table.fragment_tol;
    let index = Ms2RunIndex::new(table.scans);
    let c = attach_ms2(
        &mut ex.rows,
        &stored.library,
        &stored.spectra,
        &stored.donor_scans,
        &index,
        tol,
    );
    drop(index);
    log_ms2(stem, &c, t.elapsed().as_secs_f64());
    let path = artifact(dir, stem, EXTRACT_SUFFIX);
    let mut tags = vec![
        (FINGERPRINT_KEY.to_string(), fingerprint.to_string()),
        (TOL_KEY.to_string(), ep.tol_ppm.to_string()),
        (RUN_SCALES_KEY.to_string(), ex.counts.scales.describe()),
    ];
    if !selection.is_off() {
        let aligned = stored.library.alignment.maps.iter().flatten().count();
        tags.push((SELECTION_KEY.to_string(), selection.describe(aligned)));
    }
    write_extract_table_tagged(&path, &ex.rows, &tags)?;
    Ok(ex.rows)
}

/// Template-fit settings for an MS1 tolerance (`--lfq-tol`, ppm only;
/// default 10 ppm).
pub(crate) fn evidence_params(tol: Option<Tolerance>) -> Result<EvidenceParams, String> {
    match tol {
        None => Ok(EvidenceParams::default()),
        Some(Tolerance::Ppm(v)) => Ok(EvidenceParams {
            tol_ppm: v,
            ..EvidenceParams::default()
        }),
        Some(Tolerance::Da(v)) => Err(format!(
            "match-between-runs extraction needs a ppm --lfq-tol, got {v} Da"
        )),
    }
}

/// Peptidoform and charge of a `"<peptidoform>/<charge>"` key.
fn split_key(key: &str) -> &str {
    key.rsplit_once('/').map_or(key, |(p, _)| p)
}

/// Bare sequence and modification records of a peptidoform written as
/// residues with `[...]` tags (a UNIMOD accession or a signed mass delta).
fn parse_peptidoform(peptidoform: &str) -> (String, Vec<ModRecord>) {
    let mut seq = String::new();
    let mut mods = Vec::new();
    let mut chars = peptidoform.chars();
    while let Some(c) = chars.next() {
        if c != '[' {
            seq.push(c);
            continue;
        }
        let tag: String = chars.by_ref().take_while(|&c| c != ']').collect();
        if let Some(aa) = seq.chars().last() {
            mods.push(ModRecord {
                accession: tag.starts_with("UNIMOD:").then(|| tag.clone()),
                name: tag,
                position: seq.chars().count() as i32,
                amino_acid: Some(aa),
            });
        }
    }
    (seq, mods)
}

/// The protein groups of the QPX views: the inference over every run's PSMs
/// and the protein q-value cut.
pub(crate) struct ProteinGroups<'a> {
    pub inf: &'a Inference,
    pub max_q: f64,
}

/// What the feature writers did with the features of the protein groups.
#[derive(Default)]
pub(crate) struct GroupTally {
    /// By group id, the features written.
    pub counts: HashMap<usize, FeatureCounts>,
    /// Features whose peptide counts for no target group.
    pub no_group: usize,
    /// Features whose group is above the protein q-value cut.
    pub above_cut: usize,
}

impl GroupTally {
    /// The group accessions and the group of a feature of the peptide
    /// `sequence`, counted; `None` (counted) when it has no group at the cut.
    pub(crate) fn assign(
        &mut self,
        pg: &ProteinGroups<'_>,
        sequence: &str,
    ) -> Option<(Vec<ProteinPosition>, FeatureGroup)> {
        match pg.inf.group_at(sequence, pg.max_q) {
            GroupAt::Group { group, unique } => {
                let c = self.counts.entry(group.id).or_default();
                c.total += 1;
                c.unique += unique as usize;
                Some((
                    group
                        .accessions
                        .iter()
                        .map(|a| (a.clone(), None, None, None, None))
                        .collect(),
                    FeatureGroup {
                        anchor: group.lead.clone(),
                        unique,
                        global_q: group.q_value,
                    },
                ))
            }
            GroupAt::NoGroup => {
                self.no_group += 1;
                None
            }
            GroupAt::AboveCut => {
                self.above_cut += 1;
                None
            }
        }
    }

    /// Writes `quantms.pg.parquet` into `dir`: the groups of `pg` at its cut,
    /// with the features counted here, over `runs`.
    pub(crate) fn write_pg(
        &self,
        dir: &Path,
        pg: &ProteinGroups<'_>,
        runs: &[String],
    ) -> Result<(), String> {
        let path = dir.join("quantms.pg.parquet");
        output::write_pg_parquet(
            &path,
            pg.inf,
            runs,
            runs.len() == 1,
            pg.max_q,
            Some(&self.counts),
        )
        .map_err(|e| format!("write {}: {e}", path.display()))?;
        let kept = pg.inf.groups.iter().filter(|g| g.q_value <= pg.max_q);
        let decoys = pg
            .inf
            .decoy_groups
            .iter()
            .filter(|g| g.q_value <= pg.max_q)
            .count();
        eprintln!(
            "[proteins] wrote {} ({} target groups at protein q<={} over {} run(s), {} with features; {} decoy groups at the cut); features left out: {} of peptides in no protein group, {} of groups above the cut",
            path.display(),
            kept.count(),
            pg.max_q,
            runs.len(),
            self.counts.len(),
            decoys,
            self.no_group,
            self.above_cut
        );
        Ok(())
    }
}

/// The QPX feature rows of the MBR table: every row except the
/// below-threshold transfers and identified rows, which stay in
/// `<base>.mbr.parquet` only. With `groups`, each row carries its peptide's
/// protein group, and rows whose peptide has no group at or below the cut
/// are left out (counted in `tally`). Feature ids continue from `last_id`.
fn feature_records(
    stored: &StoredLibrary,
    features: &[MbrFeature],
    last_id: &mut i64,
    groups: Option<&ProteinGroups<'_>>,
    tally: &mut GroupTally,
) -> Vec<FeatureRecord> {
    features
        .iter()
        .filter(|f| {
            !matches!(
                f.evidence_type,
                EvidenceType::BelowThreshold | EvidenceType::IdentifiedBelowThreshold
            )
        })
        .filter_map(|f| {
            let m = &stored.library.meta[f.precursor as usize];
            let peptidoform = split_key(&m.key).to_string();
            let (sequence, modifications) = parse_peptidoform(&peptidoform);
            let (proteins, group): (Vec<ProteinPosition>, _) = match groups {
                Some(pg) => {
                    let (proteins, group) = tally.assign(pg, &sequence)?;
                    (proteins, Some(group))
                }
                None => (
                    stored.proteins[f.precursor as usize]
                        .split(';')
                        .filter(|a| !a.is_empty())
                        .map(|a| (a.to_string(), None, None, None, None))
                        .collect(),
                    None,
                ),
            };
            *last_id += 1;
            let run = stored.runs[f.run as usize].clone();
            Some(FeatureRecord {
                feature_id: *last_id,
                sequence,
                peptidoform,
                modifications,
                charge: m.charge as i16,
                pep: None,
                is_decoy: false,
                calculated_mz: m.mz as f32,
                observed_mz: f.mass_err_ppm.map_or(m.mz, |e| m.mz * (1.0 + e * 1e-6)) as f32,
                mass_error_ppm: f.mass_err_ppm.map(|e| e as f32),
                additional_scores: vec![("library_q".into(), f.library_q, false)],
                run_file_name: run.clone(),
                cv_params: Vec::new(),
                scan: Vec::new(),
                rt: f.apex_rt_s.map(|t| t as f32),
                missed_cleavages: None,
                intensities: f
                    .quantity
                    .map(|q| vec![("LFQ".to_string(), q as f32)])
                    .unwrap_or_default(),
                additional_intensities: Vec::new(),
                proteins,
                id_run_file_name: (f.evidence_type == EvidenceType::Identified).then_some(run),
                rt_start: None,
                rt_stop: None,
                psm_ids: Vec::new(),
                mbr: Some(MbrColumns {
                    evidence_type: f.evidence_type.as_str().to_string(),
                    transfer_q: f.transfer_q,
                    quant_q: f.quant_q,
                    quantity_se: f.quantity_se,
                }),
                group,
            })
        })
        .collect()
}

/// Transfer FDR, identified quality q-values and the MBR outputs over the
/// extract rows of `source` ([`finalize_stream`] within `limits`); with
/// `parquet_dir`, also `<parquet_dir>/quantms.feature.parquet`, written as the
/// rows come, and with `groups` the features' protein groups and
/// `<parquet_dir>/quantms.pg.parquet`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn finalize_and_write(
    stored: &StoredLibrary,
    source: &dyn ExtractSource,
    mbr_fdr: f64,
    quant_feature_fdr: Option<f64>,
    seed: u64,
    limits: FinalizeLimits,
    base: &Path,
    parquet_dir: Option<&Path>,
    groups: Option<&ProteinGroups<'_>>,
) -> Result<(), String> {
    let (table, validation) = mbr_paths(base);
    let mut feature_out = match parquet_dir {
        Some(dir) => {
            std::fs::create_dir_all(dir).map_err(|e| format!("create {}: {e}", dir.display()))?;
            let path = dir.join("quantms.feature.parquet");
            Some((FeatureParquetWriter::new_mbr(&path), path, 0i64))
        }
        None => None,
    };
    let params = FinalizeParams {
        mbr_fdr,
        quant_feature_fdr,
        seed,
        limits,
        mbr_path: &table,
        validation_path: &validation,
    };
    let mut tally = GroupTally::default();
    let summary = finalize_stream(&stored.runs, &stored.library, source, &params, &mut |f| {
        if let Some((w, path, last_id)) = feature_out.as_mut() {
            w.push(feature_records(stored, &f, last_id, groups, &mut tally))
                .map_err(|e| format!("write {}: {e}", path.display()))?;
        }
        Ok(())
    })?;
    eprintln!(
        "[mbr] wrote {} ({} rows; transfer_q and quant_q are internal quant-quality q-values, not identification FDR) and {}",
        table.display(),
        summary.rows_written,
        validation.display()
    );
    if let Some((w, path, n)) = feature_out {
        w.finish()
            .map_err(|e| format!("write {}: {e}", path.display()))?;
        eprintln!(
            "[mbr] wrote {} ({n} feature rows with evidence type; below-threshold transfers and identified rows are in the mbr.parquet only)",
            path.display(),
        );
    }
    if let (Some(dir), Some(pg)) = (parquet_dir, groups) {
        tally.write_pg(dir, pg, &stored.runs)?;
    }
    Ok(())
}

/// Library, extract (one run at a time) and finalize in this process, from
/// the per-run artifacts in `dir`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn run_in_process(
    dir: &Path,
    stems: &[String],
    scores: ScoreSet,
    ep: &EvidenceParams,
    max_q: f64,
    donor_q: f64,
    donor_pep: Option<f64>,
    mbr_fdr: f64,
    quant_feature_fdr: Option<f64>,
    base: &Path,
    parquet_dir: Option<&Path>,
    groups: Option<&ProteinGroups<'_>>,
) -> Result<(), String> {
    let t0 = std::time::Instant::now();
    let stems = library_runs(dir, stems)?;
    let stored = build_and_write_library(
        dir,
        &stems,
        ScoreLoader::Merged(Box::new(move |_| Ok(scores))),
        max_q,
        donor_q,
        donor_pep,
    )?;
    let fingerprint = library_fingerprint(&stored);
    // One run at a time: each extraction is parallel over its precursors and
    // holds that run's MS1; finalize reads the extract tables back.
    for r in 0..stems.len() {
        extract_and_write(
            dir,
            &stored,
            &fingerprint,
            r as u16,
            ep,
            &TransferSelection::default(),
        )?;
    }
    let files = ExtractFiles(
        stems
            .iter()
            .map(|s| artifact(dir, s, EXTRACT_SUFFIX))
            .collect(),
    );
    finalize_and_write(
        &stored,
        &files,
        mbr_fdr,
        quant_feature_fdr,
        MBR_SEED,
        FinalizeLimits::default(),
        base,
        parquet_dir,
        groups,
    )?;
    eprintln!("[PHASE mbr: {:.2}s]", t0.elapsed().as_secs_f64());
    Ok(())
}

fn run_library(a: QuantLibraryArgs) -> Result<(), String> {
    let stems = library_runs(&a.quant_dir, &a.runs)?;
    let source = score_source(&a.psms, &a.decoy_psms, &stems)?;
    build_and_write_library(
        &a.quant_dir,
        &stems,
        score_loader(&source),
        a.mbr_library_fdr,
        a.quant_fdr,
        a.pep,
    )
    .map(|_| ())
}

fn run_extract(a: QuantExtractArgs) -> Result<(), String> {
    let ep = evidence_params(a.lfq_tol)?;
    let stored = read_library(&a.quant_dir)?;
    let run = stored
        .runs
        .iter()
        .position(|s| *s == a.run)
        .ok_or_else(|| {
            format!(
                "run {:?} is not in the library of {} (runs: {})",
                a.run,
                a.quant_dir.display(),
                stored.runs.join(", ")
            )
        })?;
    let fingerprint = library_fingerprint(&stored);
    let selection = if a.mbr_preselect {
        TransferSelection {
            donor_gate: true,
            min_donor_fraction: a.mbr_min_donor_fraction,
        }
    } else {
        TransferSelection::OFF
    };
    extract_and_write(
        &a.quant_dir,
        &stored,
        &fingerprint,
        run as u16,
        &ep,
        &selection,
    )
    .map(|_| ())
}

/// Returns the allocator's free pages to the system. glibc keeps memory freed
/// by many small allocations (the tens of millions of PSM rows of the
/// protein inference) in its arenas, where the process's RSS still counts it.
fn release_free_memory() {
    #[cfg(all(target_os = "linux", target_env = "gnu"))]
    // SAFETY: malloc_trim only releases free heap memory; it is thread safe.
    unsafe {
        libc::malloc_trim(0);
    }
}

fn run_finalize(a: QuantFinalizeArgs) -> Result<(), String> {
    let stored = read_library(&a.quant_dir)?;
    let fingerprint = library_fingerprint(&stored);
    let mut paths = Vec::with_capacity(stored.runs.len());
    let mut tol: Option<(String, PathBuf)> = None;
    for stem in &stored.runs {
        let path = artifact(&a.quant_dir, stem, EXTRACT_SUFFIX);
        let tags = read_extract_tags(&path)?;
        let get = |k: &str| tags.iter().find(|t| t.0 == k).map(|t| t.1.clone());
        if get(FINGERPRINT_KEY).as_deref() != Some(fingerprint.as_str()) {
            return Err(format!(
                "{}: extracted with another library than {} (stale extract); run `andes quant extract` again",
                path.display(),
                a.quant_dir.join(output::library_table::LIBRARY_FILE).display()
            ));
        }
        let t = get(TOL_KEY).ok_or_else(|| format!("{}: no {TOL_KEY}", path.display()))?;
        match &tol {
            None => tol = Some((t, path.clone())),
            Some((t0, p0)) if *t0 != t => {
                return Err(format!(
                    "{}: extracted at {t} ppm, {} at {t0} ppm; every run needs the same --lfq-tol",
                    path.display(),
                    p0.display()
                ));
            }
            _ => {}
        }
        if let Some(scales) = get(RUN_SCALES_KEY) {
            eprintln!("[mbr] run {stem}: {scales}");
        }
        paths.push(path);
    }
    let limits = FinalizeLimits {
        max_train_rows: a.max_train_rows,
        transfer_models: a.transfer_models,
        max_norm_cells: a.max_norm_cells,
        max_donor_fits: a.max_donor_fits,
        max_held_scores: a.max_held_scores,
        temp_dir: a.temp_dir,
    };
    let inf = match (&a.output_parquet, &a.database) {
        (Some(_), Some(database)) => Some(
            crate::inference_cmd::infer_from_files(database, &a.proteins, &a.psms, &a.decoy_psms)
                .map_err(|e| e.to_string())?,
        ),
        _ => None,
    };
    // The inference's PSM maps and database are dropped by now; hand their
    // pages back before the finalize passes allocate their own.
    if inf.is_some() {
        release_free_memory();
    }
    let groups = inf.as_ref().map(|inf| ProteinGroups {
        inf,
        max_q: a.proteins.protein_fdr,
    });
    finalize_and_write(
        &stored,
        &ExtractFiles(paths),
        a.mbr_fdr,
        a.quant_feature_fdr.0,
        a.seed,
        limits,
        &a.output,
        a.output_parquet.as_deref(),
        groups.as_ref(),
    )
}

/// Runs one `andes quant` subcommand.
pub(crate) fn run_quant_command(cmd: QuantCommand) -> Result<(), Box<dyn std::error::Error>> {
    if let Some(n) = cmd.threads() {
        // Ignore the error of a pool that is already built (tests).
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global();
    }
    match cmd {
        QuantCommand::Library(a) => run_library(a),
        QuantCommand::Extract(a) => run_extract(a),
        QuantCommand::Finalize(a) => run_finalize(*a),
    }
    .map_err(Into::into)
}

#[cfg(test)]
mod tests {
    use super::*;
    use quant::library::build_library;

    #[test]
    fn peptidoform_tags_become_modifications() {
        let (seq, mods) = parse_peptidoform("LC[UNIMOD:4]AK[+8.01420]");
        assert_eq!(seq, "LCAK");
        assert_eq!(mods.len(), 2);
        assert_eq!(
            (
                mods[0].position,
                mods[0].amino_acid,
                mods[0].accession.as_deref()
            ),
            (2, Some('C'), Some("UNIMOD:4"))
        );
        assert_eq!((mods[1].position, mods[1].accession.clone()), (4, None));
        assert_eq!(split_key("LC[UNIMOD:4]AK/3"), "LC[UNIMOD:4]AK");
    }

    #[test]
    fn library_runs_checks_artifacts_and_duplicates() {
        let dir = tempfile::tempdir().unwrap();
        let d = dir.path();
        for f in [
            "b.quant.parquet",
            "b.ms1.parquet",
            "b.ms2.parquet",
            "a.quant.parquet",
            "a.ms2.parquet",
        ] {
            std::fs::write(d.join(f), b"x").unwrap();
        }
        let e = library_runs(d, &[]).unwrap_err();
        assert!(e.contains("a.ms1.parquet") && !e.contains("b.ms"), "{e}");
        std::fs::remove_file(d.join("a.ms2.parquet")).unwrap();
        std::fs::write(d.join("a.ms1.parquet"), b"x").unwrap();
        let e = library_runs(d, &[]).unwrap_err();
        assert!(e.contains("a.ms2.parquet"), "{e}");
        std::fs::write(d.join("a.ms2.parquet"), b"x").unwrap();
        assert_eq!(library_runs(d, &[]).unwrap(), vec!["a", "b"]);
        assert_eq!(
            library_runs(d, &["b".into(), "a".into()]).unwrap(),
            vec!["a", "b"]
        );
        let e = library_runs(d, &["a".into(), "a".into()]).unwrap_err();
        assert!(e.contains("twice"), "{e}");
        let e = library_runs(d, &["c".into()]).unwrap_err();
        assert!(
            e.contains("c.quant.parquet")
                && e.contains("c.ms1.parquet")
                && e.contains("c.ms2.parquet"),
            "{e}"
        );
    }

    #[test]
    fn score_files_are_merged_or_one_pair_per_run() {
        let stems: Vec<String> = vec!["a".into(), "b".into()];
        let v = |x: &[&str]| x.iter().map(|s| s.to_string()).collect::<Vec<_>>();
        assert_eq!(
            score_source(&v(&["t.tsv"]), &v(&["d.tsv"]), &stems).unwrap(),
            ScoreSource::Merged("t.tsv".into(), "d.tsv".into())
        );
        let s = score_source(&v(&["b=tb", "a=ta"]), &v(&["a=da", "b=db"]), &stems).unwrap();
        let ScoreSource::PerRun(m) = s else {
            panic!("{s:?}")
        };
        assert_eq!(m["a"], ("ta".into(), "da".into()));
        assert_eq!(m["b"], ("tb".into(), "db".into()));
        let e = score_source(&v(&["a=ta"]), &v(&["a=da"]), &stems).unwrap_err();
        assert!(e.contains("no file for runs b"), "{e}");
        let e = score_source(&v(&["a=ta", "b.tsv"]), &v(&["a=da", "b=db"]), &stems).unwrap_err();
        assert!(e.contains("cannot be mixed"), "{e}");
        let e = score_source(&v(&["a=ta", "a=tb"]), &v(&["a=da", "b=db"]), &stems).unwrap_err();
        assert!(e.contains("twice"), "{e}");
        assert!(score_source(&v(&["t1", "t2"]), &v(&["d"]), &stems).is_err());
    }

    #[test]
    fn extraction_tolerance_is_ppm() {
        assert_eq!(evidence_params(None).unwrap().tol_ppm, 10.0);
        assert_eq!(
            evidence_params(Some(Tolerance::Ppm(5.0))).unwrap().tol_ppm,
            5.0
        );
        assert!(evidence_params(Some(Tolerance::Da(0.01))).is_err());
    }

    #[test]
    fn scores_need_the_score_column() {
        let mut t = HashMap::new();
        t.insert(
            "s1".to_string(),
            PercolatorPsm {
                psm_id: "s1".into(),
                q_value: 0.0,
                pep: 0.0,
                score: None,
                peptide: "K.A.R".into(),
                proteins: "P".into(),
            },
        );
        let none = HashMap::new();
        let e = rescore_scores(&t, &none, "t.tsv").err().unwrap();
        assert!(e.contains("t.tsv") && e.contains("score"), "{e}");
        t.get_mut("s1").unwrap().score = Some(2.0);
        let e = rescore_scores(&t, &t.clone(), "t.tsv").err().unwrap();
        assert!(e.contains("both"), "{e}");
        let mut d = t.clone();
        let p = d.remove("s1").unwrap();
        d.insert(
            "s2".into(),
            PercolatorPsm {
                psm_id: "s2".into(),
                ..p
            },
        );
        let m = rescore_scores(&t, &d, "t.tsv").unwrap();
        assert_eq!(m.len(), 2);
        assert!(!m.get("s1").unwrap().is_decoy() && m.get("s2").unwrap().is_decoy());
        assert!(m.get("s3").is_none());
    }

    #[test]
    fn percolator_files_are_read_line_by_line() {
        let dir = tempfile::tempdir().unwrap();
        let write = |name: &str, rows: &[&str]| {
            let path = dir.path().join(name);
            let mut text =
                "PSMId\tscore\tq-value\tposterior_error_prob\tpeptide\tproteinIds\r\n".to_string();
            for r in rows {
                text.push_str(r);
                text.push_str("\r\n");
            }
            std::fs::write(&path, text).unwrap();
            path
        };
        let t = write(
            "t.tsv",
            &[
                "a_1_1\t2.5\t0.001\t0.01\tK.PEPC+57.02146TIDEK.A\tP1\tP2",
                "",
                "b_2_1\t1.5\t0.002\t0.02\t-.PEPM[+15.99491]K.-\tP3",
                "e_4_1\t1.0\t0.003\t0.03\tR.PEPS+79.96633TK.A\tP4",
            ],
        );
        let d = write("d.tsv", &["c_3_1\t-1.0\t0.5\t0.9\tK.KEDITPEP.A\tDECOY_P1"]);
        let set = percolator_scores(&t, &d, 0).unwrap();
        assert_eq!(set.len(), 4);
        let a = *set.get("a_1_1").unwrap();
        assert_eq!(
            (a.score, a.q, a.pep, a.is_decoy()),
            (2.5, 0.001, 0.01, false)
        );
        let c = *set.get("c_3_1").unwrap();
        assert!(c.is_decoy() && c.score == -1.0);
        assert!(set.get("a_1").is_none());
        let cls = |t: &str| IsomerClass::parse(t, bundled_unimod_mass).unwrap();
        assert!(set.agrees(&a, false, &cls("PEPC[UNIMOD:4]TIDEK")));
        assert!(set.agrees(&a, false, &cls("PEPC[+57.021464]TIDEK")));
        assert!(
            !set.agrees(&a, false, &cls("PEPCTIDEK")),
            "another modification composition"
        );
        assert!(!set.agrees(&a, true, &cls("PEPC[UNIMOD:4]TIDEK")));
        assert!(!set.agrees(&a, false, &cls("PEPC[UNIMOD:4]TIDEKK")));
        let b = *set.get("b_2_1").unwrap();
        assert!(set.agrees(&b, false, &cls("PEPM[UNIMOD:35]K")));
        assert!(!set.agrees(&b, false, &cls("PEPMK")));
        let e = *set.get("e_4_1").unwrap();
        assert!(
            set.agrees(&e, false, &cls("PEPST[UNIMOD:21]K")),
            "a positional isomer"
        );
        assert!(!set.agrees(&e, false, &cls("PEPS[UNIMOD:21]T[UNIMOD:21]K")));
        assert!(!set.agrees(&e, false, &cls("PEPST[+79.95682]K")));

        let dup = write(
            "dup.tsv",
            &[
                "a_1_1\t2.5\t0.001\t0.01\tK.PEPTIDEK.A\tP1",
                "a_1_1\t2.0\t0.001\t0.01\tK.PEPTIDEK.A\tP1",
            ],
        );
        let e = percolator_scores(&dup, &d, 0).err().unwrap();
        assert!(e.contains("dup.tsv") && e.contains("a_1_1"), "{e}");
        let e = percolator_scores(&t, &t, 0).err().unwrap();
        assert!(e.contains("both a target and a decoy"), "{e}");
        let no_score = write("ns.tsv", &["a_1_1\tnan\t0.001\t0.01\tK.PEPTIDEK.A\tP1"]);
        let e = percolator_scores(&no_score, &d, 0).err().unwrap();
        assert!(e.contains("ns.tsv") && e.contains("score"), "{e}");
        let bad = write("bad.tsv", &["a_1_1\t2.5\tx\t0.01\tK.PEPTIDEK.A\tP1"]);
        let e = percolator_scores(&bad, &d, 0).err().unwrap();
        assert!(e.contains("bad.tsv") && e.contains("q-value"), "{e}");
        let e = percolator_scores(&dir.path().join("missing.tsv"), &d, 0)
            .err()
            .unwrap();
        assert!(e.contains("missing.tsv"), "{e}");
    }

    #[test]
    fn feature_records_leave_out_below_threshold_rows() {
        let meta: Vec<PrecursorMeta> = (0..2)
            .map(|i| PrecursorMeta {
                key: format!("PEPT{}IDEK/2", ["", "M[UNIMOD:35]"][i]),
                mz: 500.0 + i as f64,
                charge: 2,
                envelope: [0.55, 0.3, 0.1, 0.05],
            })
            .collect();
        let library = build_library(&[], meta, 2, 0.01, 0.01, None, &AlignParams::default());
        let stored = StoredLibrary {
            runs: vec!["a".into(), "b".into()],
            max_q: 0.01,
            donor_q: 0.01,
            donor_pep: None,
            library,
            proteins: vec!["P1".into(), "P2".into()],
            spectra: vec![None; 2],
            donor_scans: vec![Vec::new(); 2],
        };
        let feature = |precursor: u32, run: u16, evidence_type, transfer_q| MbrFeature {
            precursor,
            run,
            evidence_type,
            quantity: Some(1e6),
            quantity_se: Some(1e4),
            quantity_normalized: Some(1e6),
            transfer_q,
            quant_q: transfer_q.or(Some(0.003)),
            apex_rt_s: Some(100.0),
            library_q: 0.001,
            mass_err_ppm: Some(0.5),
        };
        let features = [
            feature(0, 0, EvidenceType::Identified, None),
            feature(0, 1, EvidenceType::BelowThreshold, Some(0.2)),
            feature(1, 0, EvidenceType::TransferMs1, Some(0.004)),
            feature(1, 1, EvidenceType::IdentifiedBelowThreshold, None),
        ];
        let records = feature_records(&stored, &features, &mut 0, None, &mut GroupTally::default());
        let got: Vec<(i64, &str, &str)> = records
            .iter()
            .map(|r| {
                (
                    r.feature_id,
                    r.run_file_name.as_str(),
                    r.mbr.as_ref().unwrap().evidence_type.as_str(),
                )
            })
            .collect();
        assert_eq!(
            got,
            [(1, "a", "identified"), (2, "a", "transfer_ms1")],
            "below-threshold transfers and identified rows are not feature rows"
        );
        assert_eq!(records[1].peptidoform, "PEPTM[UNIMOD:35]IDEK");
        let quant_q: Vec<Option<f64>> = records
            .iter()
            .map(|r| r.mbr.as_ref().unwrap().quant_q)
            .collect();
        assert_eq!(quant_q, [Some(0.003), Some(0.004)]);
        assert_eq!(records[0].proteins[0].0, "P1");
        assert!(records[0].group.is_none());

        // With protein groups: PEPTIDEK's group passes the cut, PEPTMIDEK's
        // does not.
        let group = |id: usize, lead: &str, q: f64| output::ProteinGroup {
            id,
            lead: lead.into(),
            families: vec![lead.into()],
            accessions: vec![lead.into(), format!("{lead}-2")],
            peptides: 1,
            unique: 0,
            razor: 1,
            psms: 1,
            score: 3.0,
            q_value: q,
            subsumed: Vec::new(),
            peptides_per_accession: Vec::new(),
            is_decoy: false,
        };
        let inf = Inference {
            groups: vec![group(1, "Q1", 0.0), group(2, "Q2", 0.05)],
            peptide_group: [("PEPTIDEK".into(), 0), ("PEPTMIDEK".into(), 1)].into(),
            shared_peptides: ["PEPTIDEK".to_string()].into(),
            decoy_groups: Vec::new(),
        };
        let groups = ProteinGroups {
            inf: &inf,
            max_q: 0.01,
        };
        let mut tally = GroupTally::default();
        let records = feature_records(&stored, &features, &mut 0, Some(&groups), &mut tally);
        assert_eq!(records.len(), 1);
        let r = &records[0];
        assert_eq!(r.feature_id, 1);
        let accs: Vec<&str> = r.proteins.iter().map(|p| p.0.as_str()).collect();
        assert_eq!(accs, ["Q1", "Q1-2"]);
        assert_eq!(
            r.group,
            Some(FeatureGroup {
                anchor: "Q1".into(),
                unique: false,
                global_q: 0.0
            })
        );
        assert_eq!(
            tally.counts[&1],
            FeatureCounts {
                unique: 0,
                total: 1
            }
        );
        assert_eq!((tally.no_group, tally.above_cut), (0, 1));
    }

    /// 200 target PSMs of distinct peptides in run `a`, scored by a rescoring
    /// result at q = 0.001 and PEP = `pep`; the result names `wrong.0` (a
    /// SpecId) with peptide `wrong.1`. The MS2 table holds every PSM's scan:
    /// the even ones carry the peptide's b/y ladder, the odd ones noise.
    fn one_run_library(
        dir: &Path,
        pep: f64,
        wrong: Option<(&str, &str)>,
    ) -> (Vec<String>, HashMap<String, PercolatorPsm>) {
        use output::run_quant_table::{write_run_quant_table, RunPsmRow};
        let peptide = |i: usize| -> String {
            let aa = b"ACDEFGHIKLMNPQRSTVWY";
            let tail: String = [aa[i % 20], aa[(i / 20) % 20]]
                .iter()
                .map(|&c| c as char)
                .collect();
            format!("PEPTIDE{tail}K")
        };
        let rows: Vec<RunPsmRow> = (0..200)
            .map(|i| RunPsmRow {
                spec_id: format!("s{i}_{i}_1"),
                rt_s: 100.0 + i as f64,
                precursor_key: format!("{}/2", peptide(i)),
                peptidoform: peptide(i),
                charge: 2,
                mz: 500.0 + i as f64,
                is_decoy: false,
                search_score: 10.0,
                faims_cv: None,
                envelope: vec![0.55, 0.3, 0.1, 0.05],
                proteins: "P1".into(),
                ion_mobility: None,
            })
            .collect();
        write_run_quant_table(&artifact(dir, "a", TABLE_SUFFIX), "a", &rows, &[]).unwrap();
        let scans: Vec<quant::ms2::Ms2Scan> = rows
            .iter()
            .enumerate()
            .map(|(i, r)| {
                let ladder =
                    fragment_ladder(&residue_masses(&r.peptidoform, &HashMap::new()).unwrap(), 2);
                let peaks: Vec<(f32, f32)> = if i % 2 == 0 {
                    ladder
                        .iter()
                        .enumerate()
                        .map(|(k, f)| (f.mz as f32, 100.0 + k as f32))
                        .collect()
                } else {
                    vec![(150.05, 10.0), (250.5, 20.0)]
                };
                quant::ms2::Ms2Scan {
                    native_id: format!("scan={i}"),
                    spec_key: format!("s{i}_{i}"),
                    rt: r.rt_s,
                    precursor_mz: r.mz,
                    charge: 2,
                    isolation: Some((r.mz, 0.7, 0.7)),
                    faims_cv: None,
                    peaks,
                }
            })
            .collect();
        output::ms2_table::write_ms2_table(
            &artifact(dir, "a", MS2_SUFFIX),
            &output::ms2_table::Ms2Table {
                run: "a".into(),
                fragment_tol: Tolerance::Ppm(20.0),
                top_peaks: 100,
                scans,
            },
        )
        .unwrap();
        let scores = rows
            .iter()
            .map(|r| {
                let peptide = if let Some((_, p)) = wrong.filter(|w| w.0 == r.spec_id) {
                    p.to_string()
                } else {
                    format!("K.{}.A", r.peptidoform)
                };
                let psm = PercolatorPsm {
                    psm_id: r.spec_id.clone(),
                    q_value: 0.001,
                    pep,
                    score: Some(100.0 + r.rt_s),
                    peptide,
                    proteins: "P1".into(),
                };
                (r.spec_id.clone(), psm)
            })
            .collect();
        (vec!["a".to_string()], scores)
    }

    fn open_peptide(i: usize) -> String {
        let aa = b"ACDEFGHIKLMNPQRSTVWY";
        let tail: String = [aa[i % 20], aa[(i / 20) % 20]]
            .iter()
            .map(|&c| c as char)
            .collect();
        format!("PEPTIDE{tail}K")
    }

    #[test]
    fn scan_of_strips_rank_and_row() {
        let keys: HashMap<&str, usize> = [("a/scan=5_5", 0usize), ("x_1", 1)].into_iter().collect();
        assert_eq!(scan_of(&keys, "a/scan=5_5_1"), Some(0));
        assert_eq!(scan_of(&keys, "a/scan=5_5_1_2"), Some(0));
        assert_eq!(scan_of(&keys, "x_1_1"), Some(1));
        assert_eq!(scan_of(&keys, "a/scan=6_6_1"), None);
        assert_eq!(scan_of(&keys, "nounderscore"), None);
    }

    #[test]
    fn donors_honour_the_pep_cap() {
        // q = 0.001 and PEP = 0.5 for every PSM.
        let dir = tempfile::tempdir().unwrap();
        let (stems, t) = one_run_library(dir.path(), 0.5, None);
        let build = |pep: Option<f64>| {
            let set = rescore_scores(&t, &HashMap::new(), "t.tsv").unwrap();
            build_and_write_library(
                dir.path(),
                &stems,
                ScoreLoader::Merged(Box::new(move |_| Ok(set))),
                0.01,
                0.01,
                pep,
            )
            .unwrap()
        };
        let open = build(None);
        assert_eq!(open.library.entries.len(), 200);
        // The even PSMs' scans carry their ladder: a spectrum each; every
        // entry keeps its one donor scan.
        let with: Vec<usize> = (0..open.spectra.len())
            .filter(|&id| open.spectra[id].is_some())
            .collect();
        assert_eq!(with.len(), 100);
        for &id in &with {
            let key = &open.library.meta[id].key;
            let i = (0..200)
                .find(|&i| key == &format!("{}/2", open_peptide(i)))
                .unwrap();
            assert_eq!(i % 2, 0, "{key}");
            let s = open.spectra[id].as_ref().unwrap();
            assert_eq!(s.fragments.len(), quant::ms2::LIBRARY_FRAGMENTS);
        }
        assert!(open.donor_scans.iter().all(|d| d.len() == 1));
        assert_eq!(read_library(dir.path()).unwrap().spectra, open.spectra);
        let capped = build(Some(0.01));
        assert!(
            capped.library.entries.is_empty(),
            "donors above the PEP cap"
        );
        assert_eq!(capped.donor_pep, Some(0.01));
        assert_eq!(read_library(dir.path()).unwrap().donor_pep, Some(0.01));
        assert_ne!(library_fingerprint(&open), library_fingerprint(&capped));
    }

    #[test]
    fn mismatched_joins_are_not_donors() {
        // s7_7_1 is PEPTIDEIAK in the run table; the rescoring result names
        // another peptide, or the same residues with another modification
        // composition; the same peptide joins.
        for (wrong, joins) in [
            ("K.WRONGPEPTIDEK.A", false),
            ("K.PEPTIDEIAK+8.01420.A", false),
            ("K.PEPTIDEIAK.A", true),
        ] {
            let dir = tempfile::tempdir().unwrap();
            let (stems, t) = one_run_library(dir.path(), 0.001, Some(("s7_7_1", wrong)));
            let set = rescore_scores(&t, &HashMap::new(), "t.tsv").unwrap();
            let stored = build_and_write_library(
                dir.path(),
                &stems,
                ScoreLoader::Merged(Box::new(move |_| Ok(set))),
                0.01,
                0.01,
                None,
            )
            .unwrap();
            let key_of = |e: &quant::library::LibraryEntry| {
                stored.library.meta[e.precursor as usize].key.clone()
            };
            let keys: Vec<String> = stored.library.entries.iter().map(key_of).collect();
            assert_eq!(keys.len(), 199 + usize::from(joins), "{wrong}");
            assert_eq!(
                keys.iter().any(|k| k == "PEPTIDEIAK/2"),
                joins,
                "{wrong}: the PSM is a donor"
            );
        }
    }

    #[test]
    fn a_run_with_mismatched_modifications_fails_the_join_check() {
        // Every rescoring peptide carries an extra modification the run table
        // lacks: same residues, another composition.
        let dir = tempfile::tempdir().unwrap();
        let (stems, mut t) = one_run_library(dir.path(), 0.001, None);
        for p in t.values_mut() {
            p.peptide = p.peptide.replace("K.A", "K+8.01420.A");
        }
        let set = rescore_scores(&t, &HashMap::new(), "t.tsv").unwrap();
        let e = build_and_write_library(
            dir.path(),
            &stems,
            ScoreLoader::Merged(Box::new(move |_| Ok(set))),
            0.01,
            0.01,
            None,
        )
        .err()
        .unwrap();
        assert!(
            e.contains("200 of 200") && e.contains("modifications"),
            "{e}"
        );
    }

    /// Runs `a`, `b` and `c` share 150 peptides (the alignment anchors); one
    /// phosphopeptide wins as its S-form in `a` and `c` and as its T-form in
    /// `b`, where it scores best.
    #[test]
    fn site_forms_of_one_peptide_are_one_precursor() {
        use output::run_quant_table::{write_run_quant_table, RunPsmRow};
        let dir = tempfile::tempdir().unwrap();
        let s_form = "PEPS[UNIMOD:21]TIDEK";
        let t_form = "PEPST[UNIMOD:21]IDEK";
        let row = |spec_id: String, peptidoform: &str, rt_s: f64| RunPsmRow {
            spec_id,
            rt_s,
            precursor_key: format!("{peptidoform}/2"),
            peptidoform: peptidoform.to_string(),
            charge: 2,
            mz: 520.0,
            is_decoy: false,
            search_score: 10.0,
            faims_cv: None,
            envelope: vec![0.55, 0.3, 0.1, 0.05],
            proteins: "P1".into(),
            ion_mobility: None,
        };
        let mut scores = HashMap::new();
        let stems: Vec<String> = ["a", "b", "c"].map(String::from).to_vec();
        for (r, stem) in stems.iter().enumerate() {
            let shift = 5.0 * r as f64;
            let mut rows: Vec<RunPsmRow> = (0..150)
                .map(|i| {
                    row(
                        format!("{stem}{i}_{i}_1"),
                        &open_peptide(i),
                        100.0 + 10.0 * i as f64 + shift,
                    )
                })
                .collect();
            let (form, score) = match r {
                1 => (t_form, 600.0),
                _ => (s_form, 500.0 - 100.0 * r as f64),
            };
            rows.push(row(format!("{stem}p_900_1"), form, 900.0 + shift));
            for (i, x) in rows.iter().enumerate() {
                let pin = x.peptidoform.replace("[UNIMOD:21]", "+79.96633");
                let score = if i == 150 { score } else { 100.0 + x.rt_s };
                scores.insert(
                    x.spec_id.clone(),
                    PercolatorPsm {
                        psm_id: x.spec_id.clone(),
                        q_value: 0.001,
                        pep: 0.001,
                        score: Some(score),
                        peptide: format!("K.{pin}.A"),
                        proteins: "P1".into(),
                    },
                );
            }
            write_run_quant_table(
                &artifact(dir.path(), stem, TABLE_SUFFIX),
                stem,
                &rows,
                &[("UNIMOD:21".to_string(), 79.966331)],
            )
            .unwrap();
            output::ms2_table::write_ms2_table(
                &artifact(dir.path(), stem, MS2_SUFFIX),
                &output::ms2_table::Ms2Table {
                    run: stem.clone(),
                    fragment_tol: Tolerance::Ppm(20.0),
                    top_peaks: 100,
                    scans: Vec::new(),
                },
            )
            .unwrap();
        }
        let set = rescore_scores(&scores, &HashMap::new(), "t.tsv").unwrap();
        let stored = build_and_write_library(
            dir.path(),
            &stems,
            ScoreLoader::Merged(Box::new(move |_| Ok(set))),
            0.01,
            0.01,
            None,
        )
        .unwrap();
        let lib = &stored.library;
        assert_eq!(lib.meta.len(), 151, "150 peptides and one phosphopeptide");
        let phospho: Vec<usize> = (0..lib.meta.len())
            .filter(|&id| lib.meta[id].key.contains("UNIMOD:21"))
            .collect();
        assert_eq!(phospho.len(), 1, "one precursor for both site forms");
        let id = phospho[0];
        assert_eq!(
            lib.meta[id].key,
            format!("{t_form}/2"),
            "keyed by its best PSM"
        );
        let entries: Vec<&quant::library::LibraryEntry> = lib
            .entries
            .iter()
            .filter(|e| e.precursor as usize == id)
            .collect();
        assert_eq!(entries.len(), 1);
        assert_eq!(entries[0].donors, vec![0, 1, 2], "a donor in every run");
        let key_of = |id: usize| lib.meta[id].key.clone();
        let read = read_library(dir.path()).unwrap();
        assert_eq!(read.library.meta.len(), 151);
        assert_eq!(read.library.meta[id].key, key_of(id));
    }

    #[test]
    fn a_rescoring_peptide_without_a_mass_fails() {
        let mut t = HashMap::new();
        t.insert(
            "s1".to_string(),
            PercolatorPsm {
                psm_id: "s1".into(),
                q_value: 0.0,
                pep: 0.0,
                score: Some(1.0),
                peptide: "K.PEPC[NotAModification]K.A".into(),
                proteins: "P".into(),
            },
        );
        let e = rescore_scores(&t, &HashMap::new(), "t.tsv").err().unwrap();
        assert!(e.contains("NotAModification") && e.contains("s1"), "{e}");
    }
}
