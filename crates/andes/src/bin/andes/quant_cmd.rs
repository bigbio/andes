//! The match-between-runs stages as commands: `andes quant library`,
//! `andes quant extract` and `andes quant finalize`, and the in-process run of
//! all three behind `andes search --lfq --mbr`. Both paths call the same
//! functions on the same files in the quant directory, so a scattered run and
//! an in-process run give the same `<base>.mbr.tsv`.
//!
//! Files in the quant directory, per run stem: `<stem>.quant.parquet` and
//! `<stem>.cache` (written by `andes search --quant-dir`),
//! `<stem>.extract.parquet` (extract); once per experiment: `library.parquet`
//! and `alignment.parquet` (library).
//!
//! Run index: the position of the run's stem in the sorted stem list.
//! Precursor id: the position of its key in the sorted union of every run's
//! precursor keys. The rescorer's PSM score and PSM q-value (Percolator's
//! `score` and `q-value` columns, or the native rescorer's held-out score and
//! q-value) are joined to the run tables on (run, SpecId). Library confidence
//! ranks precursors by the score; the q-value picks the donor PSMs.

use std::collections::{BTreeMap, HashMap};
use std::path::{Path, PathBuf};

use model::Tolerance;
use output::extract_table::{read_extract_table, read_extract_tags, write_extract_table_tagged};
use output::feature_parquet::{FeatureRecord, MbrColumns};
use output::library_table::{
    library_fingerprint, read_library, write_library, StoredLibrary, FINGERPRINT_KEY,
};
use output::mbr_out::EvidenceType;
use output::run_quant_table::read_run_quant_table;
use output::{ModRecord, PercolatorPsm, ProteinPosition};
use quant::align::AlignParams;
use quant::evidence::EvidenceParams;
use quant::library::{build_library, PrecursorMeta, PsmIn};

use crate::cli::{QuantCommand, QuantExtractArgs, QuantFinalizeArgs, QuantLibraryArgs};
use crate::quant_mbr::{extract_run, finalize, mbr_paths, ExtractCounts, ExtractRow, MbrFeature};

/// Seed of the transfer FDR's folds and models, in every path.
pub(crate) const MBR_SEED: u64 = 42;

const TABLE_SUFFIX: &str = ".quant.parquet";
const CACHE_SUFFIX: &str = ".cache";
const EXTRACT_SUFFIX: &str = ".extract.parquet";
/// Extract-table metadata key of the template-fit m/z tolerance, ppm.
const TOL_KEY: &str = "andes.extract.tol_ppm";
/// Isotopes per precursor envelope in the run tables and the library.
const ENVELOPE_LEN: usize = 4;

fn artifact(dir: &Path, stem: &str, suffix: &str) -> PathBuf {
    dir.join(format!("{stem}{suffix}"))
}

/// The run stems of the library, sorted: `runs` when given, else every
/// `<stem>.quant.parquet` in `dir`. Errors: no runs, a repeated stem, or a
/// run whose PSM table or MS1 cache is missing (every such file is named).
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
                artifact(dir, s, CACHE_SUFFIX),
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

/// One rescored PSM as the library joins it.
pub(crate) struct Joined {
    pub score: f64,
    /// The rescorer's PSM q-value.
    pub q: f64,
    /// From the decoy result file.
    pub is_decoy: bool,
    /// The rescorer's peptide string (flanks and modification tags allowed).
    pub peptide: String,
}

/// Rescored PSMs by SpecId and the file(s) they came from.
pub(crate) struct ScoreSet {
    pub source: String,
    pub psms: HashMap<String, Joined>,
}

/// The target and decoy PSMs of one rescoring result. Errors: a PSM without
/// a finite score or a SpecId among both (the message names `source`).
pub(crate) fn rescore_scores(
    targets: &HashMap<String, PercolatorPsm>,
    decoys: &HashMap<String, PercolatorPsm>,
    source: &str,
) -> Result<ScoreSet, String> {
    let mut psms = HashMap::with_capacity(targets.len() + decoys.len());
    for (map, is_decoy) in [(targets, false), (decoys, true)] {
        for (id, p) in map {
            let score = p.score.ok_or_else(|| {
                format!(
                    "{source}: PSM {id:?} has no finite `score`; the library ranks precursors by the rescorer's PSM score"
                )
            })?;
            let joined = Joined {
                score,
                q: p.q_value,
                is_decoy,
                peptide: p.peptide.clone(),
            };
            if psms.insert(id.clone(), joined).is_some() {
                return Err(format!("{source}: PSM {id:?} is both a target and a decoy"));
            }
        }
    }
    Ok(ScoreSet {
        source: source.to_string(),
        psms,
    })
}

/// Residue letters of a peptide string: Percolator-style flanks (`K.` and
/// `.R`) removed, bracketed tags removed, then every character that is not
/// an uppercase letter dropped (mass deltas such as `+57.02146`).
fn residues(peptide: &str) -> String {
    let b = peptide.as_bytes();
    let mut p = peptide;
    if b.len() >= 4 && b[1] == b'.' && b[b.len() - 2] == b'.' {
        p = &peptide[2..peptide.len() - 2];
    }
    let mut out = String::with_capacity(p.len());
    let mut depth = 0usize;
    for c in p.chars() {
        match c {
            '[' => depth += 1,
            ']' => depth = depth.saturating_sub(1),
            c if depth == 0 && c.is_ascii_uppercase() => out.push(c),
            _ => {}
        }
    }
    out
}

/// Largest share of a run's joined PSMs that may disagree with the rescoring
/// result on peptide or target/decoy label.
const MAX_JOIN_MISMATCH: f64 = 0.01;

/// Reads a Percolator PSM result file; the error names it.
fn read_percolator(path: &Path) -> Result<HashMap<String, PercolatorPsm>, String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
    output::parse_psm_results(&text).map_err(|e| format!("{}: {e}", path.display()))
}

/// Rescorer scores by SpecId: one map for every run (SpecIds must then be
/// distinct across runs), or one map per run stem.
pub(crate) enum PsmScores {
    Merged(ScoreSet),
    PerRun(HashMap<String, ScoreSet>),
}

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

fn load_scores(source: &ScoreSource) -> Result<PsmScores, String> {
    Ok(match source {
        ScoreSource::Merged(t, d) => PsmScores::Merged(percolator_scores(t, d)?),
        ScoreSource::PerRun(m) => PsmScores::PerRun(
            m.iter()
                .map(|(stem, (t, d))| Ok((stem.clone(), percolator_scores(t, d)?)))
                .collect::<Result<_, String>>()?,
        ),
    })
}

/// Scores from Percolator's target and decoy PSM files.
fn percolator_scores(targets: &Path, decoys: &Path) -> Result<ScoreSet, String> {
    let t = read_percolator(targets)?;
    let d = read_percolator(decoys)?;
    for (path, m) in [(targets, &t), (decoys, &d)] {
        if let Some(p) = m.values().find(|p| p.score.is_none()) {
            return Err(format!(
                "{}: PSM {:?} has no finite `score` (missing column or value); the library ranks precursors by the rescorer's PSM score",
                path.display(),
                p.psm_id
            ));
        }
    }
    rescore_scores(
        &t,
        &d,
        &format!("{} and {}", targets.display(), decoys.display()),
    )
}

/// Builds the library over the run tables of `stems` (sorted, as from
/// [`library_runs`]) with the scores of `load_scores` joined on (run,
/// SpecId), and writes it to `dir`. `max_q` is the entry gate, `donor_q` the
/// PSM q-value of the donor PSMs. The scores are loaded after the tables, so
/// that merged scores over runs with a shared SpecId fail with the runs
/// named.
pub(crate) fn build_and_write_library(
    dir: &Path,
    stems: &[String],
    load_scores: impl FnOnce() -> Result<PsmScores, String>,
    max_q: f64,
    donor_q: f64,
) -> Result<StoredLibrary, String> {
    if !(max_q > 0.0 && max_q <= 1.0) {
        return Err(format!("library entry gate {max_q} is outside (0, 1]"));
    }
    if !(donor_q > 0.0 && donor_q <= 1.0) {
        return Err(format!(
            "donor PSM q-value threshold {donor_q} is outside (0, 1]"
        ));
    }
    // Every run's rows, and the first SpecId shared by two runs.
    let mut tables = Vec::with_capacity(stems.len());
    let mut seen: HashMap<String, usize> = HashMap::new();
    let mut shared: Option<(String, usize, usize)> = None;
    for (r, stem) in stems.iter().enumerate() {
        let path = artifact(dir, stem, TABLE_SUFFIX);
        let (run, rows) = read_run_quant_table(&path)?;
        if run != *stem {
            return Err(format!(
                "{}: the table is of run {run:?}, not {stem:?}",
                path.display()
            ));
        }
        for row in &rows {
            if let Some(&other) = seen.get(&row.spec_id) {
                if shared.is_none() && other != r {
                    shared = Some((row.spec_id.clone(), other, r));
                }
                continue;
            }
            seen.insert(row.spec_id.clone(), r);
        }
        tables.push((path, rows));
    }
    drop(seen);
    let scores = load_scores()?;
    let run_scores: Vec<&ScoreSet> = match &scores {
        PsmScores::Merged(m) => {
            if let Some((id, a, b)) = &shared {
                return Err(format!(
                    "SpecId {id:?} is in runs {} and {}, so one merged Percolator result cannot be joined; pass per-run results: --psms <stem>=<target.tsv> --decoy-psms <stem>=<decoy.tsv> for every run",
                    stems[*a], stems[*b]
                ));
            }
            vec![m; stems.len()]
        }
        PsmScores::PerRun(m) => stems
            .iter()
            .map(|s| {
                m.get(s)
                    .ok_or_else(|| format!("no rescoring scores for run {s:?}"))
            })
            .collect::<Result<_, _>>()?,
    };
    // Dense precursor ids over the sorted union of keys; the first target row
    // of a key (else its first row) gives its metadata. A decoy row whose key
    // is also a target's is dropped: a precursor is all-target or all-decoy.
    let mut keys: BTreeMap<&str, (bool, &output::run_quant_table::RunPsmRow, &Path)> =
        BTreeMap::new();
    for (path, rows) in &tables {
        for row in rows {
            let e = keys
                .entry(row.precursor_key.as_str())
                .or_insert((false, row, path.as_path()));
            if !row.is_decoy && !e.0 {
                *e = (true, row, path.as_path());
            }
        }
    }
    let mut meta = Vec::with_capacity(keys.len());
    let mut proteins = Vec::with_capacity(keys.len());
    let mut target_key = Vec::with_capacity(keys.len());
    let mut id_of: HashMap<&str, u32> = HashMap::with_capacity(keys.len());
    for (i, (&key, &(has_target, row, path))) in keys.iter().enumerate() {
        let envelope: [f32; ENVELOPE_LEN] = row.envelope.as_slice().try_into().map_err(|_| {
            format!(
                "{}: precursor {key:?} has an envelope of {} isotopes, expected {ENVELOPE_LEN}",
                path.display(),
                row.envelope.len()
            )
        })?;
        meta.push(PrecursorMeta {
            key: key.to_string(),
            mz: row.mz,
            charge: row.charge,
            envelope,
        });
        proteins.push(row.proteins.clone());
        target_key.push(has_target);
        id_of.insert(key, i as u32);
    }
    let (mut unscored_t, mut unscored_d, mut shadowed, mut mismatched) = (0usize, 0, 0, 0);
    let mut psms = Vec::new();
    for (r, (_, rows)) in tables.iter().enumerate() {
        let (mut run_joined, mut run_mismatched) = (0usize, 0usize);
        for row in rows {
            let id = id_of[row.precursor_key.as_str()];
            if row.is_decoy && target_key[id as usize] {
                shadowed += 1;
                continue;
            }
            let Some(j) = run_scores[r].psms.get(&row.spec_id) else {
                if row.is_decoy {
                    unscored_d += 1;
                } else {
                    unscored_t += 1;
                }
                continue;
            };
            run_joined += 1;
            if j.is_decoy != row.is_decoy || residues(&j.peptide) != residues(&row.peptidoform) {
                run_mismatched += 1;
            }
            psms.push(PsmIn {
                run: r as u16,
                precursor: id,
                rt: row.rt_s,
                score: j.score,
                q: j.q,
                is_decoy: row.is_decoy,
                faims_cv: row.faims_cv,
            });
        }
        if run_mismatched as f64 > MAX_JOIN_MISMATCH * run_joined as f64 {
            return Err(format!(
                "run {}: {run_mismatched} of {run_joined} PSMs joined on SpecId disagree with {} on the peptide or the target/decoy label; are these the rescoring results of this run?",
                stems[r], run_scores[r].source
            ));
        }
        mismatched += run_mismatched;
    }
    let n_rows: usize = tables.iter().map(|t| t.1.len()).sum();
    if psms.is_empty() && n_rows > 0 {
        return Err(format!(
            "none of the {n_rows} PSMs of the run tables has a rescoring score; \
             the PSM files must come from the PINs of the same searches"
        ));
    }
    drop(id_of);
    drop(keys);
    drop(tables);
    let donor_psms = psms
        .iter()
        .filter(|p| !p.is_decoy && p.q <= donor_q)
        .count();
    let library = build_library(
        &psms,
        meta,
        stems.len(),
        max_q,
        donor_q,
        &AlignParams::default(),
    );
    eprintln!(
        "[mbr] library: {} runs, {} PSMs scored ({} of them disagreeing with the rescoring result on peptide or label), {donor_psms} target PSMs at PSM q<={donor_q} identify their precursor in their run; skipped: {} target and {} decoy PSMs without a rescoring score, {} decoy PSMs sharing a target's key; {} precursors, {} entries at the internal entry gate q<={max_q} (precursor-level, not an identification FDR) with a PSM at PSM q<={donor_q}",
        stems.len(),
        psms.len(),
        mismatched,
        unscored_t,
        unscored_d,
        shadowed,
        library.meta.len(),
        library.entries.len()
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
    for (r, stem) in stems.iter().enumerate() {
        if a.maps[r].is_none() {
            eprintln!(
                "[mbr] run {stem} is unaligned ({} anchors, {:.2} of the consensus RT range): it neither gives nor receives transfers and keeps its identified quantities",
                a.anchors[r], a.coverage[r]
            );
        }
    }
    let stored = StoredLibrary {
        runs: stems.to_vec(),
        max_q,
        donor_q,
        library,
        proteins,
    };
    write_library(dir, &stored)?;
    eprintln!(
        "[mbr] wrote {} and {}",
        dir.join(output::library_table::LIBRARY_FILE).display(),
        dir.join(output::library_table::ALIGNMENT_FILE).display()
    );
    Ok(stored)
}

fn log_extract(stem: &str, rows: usize, c: &ExtractCounts) {
    eprintln!(
        "[mbr] extract {stem}: {rows} rows{}; skipped: {} entries at a FAIMS voltage the run lacks, {} identified pairs with a non-finite best PSM RT, {} identified pairs without a prediction from other runs (no pseudo-transfer), {} transfers without a prediction; width prior from {} fits{}",
        if c.unaligned {
            " (unaligned run: identified rows only)"
        } else {
            ""
        },
        c.foreign_cv,
        c.nonfinite_best_rt,
        c.no_pseudo_prediction,
        c.no_prediction,
        c.width_fits,
        if c.width_prior_fallback {
            " (too few: the default width is used)"
        } else {
            ""
        }
    );
}

/// Measures run `run` of the library on its cache and writes
/// `<stem>.extract.parquet`.
pub(crate) fn extract_and_write(
    dir: &Path,
    stored: &StoredLibrary,
    fingerprint: &str,
    run: u16,
    ep: &EvidenceParams,
) -> Result<Vec<ExtractRow>, String> {
    let stem = &stored.runs[run as usize];
    let cache = artifact(dir, stem, CACHE_SUFFIX);
    let ex = extract_run(run, &cache, &stored.library, ep)?;
    log_extract(stem, ex.rows.len(), &ex.counts);
    let path = artifact(dir, stem, EXTRACT_SUFFIX);
    write_extract_table_tagged(
        &path,
        &ex.rows,
        &[
            (FINGERPRINT_KEY.to_string(), fingerprint.to_string()),
            (TOL_KEY.to_string(), ep.tol_ppm.to_string()),
        ],
    )?;
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
                amino_acid: aa,
            });
        }
    }
    (seq, mods)
}

/// The QPX feature rows of the MBR table: every row except the
/// below-threshold transfers, which stay in `<base>.mbr.tsv` only.
fn feature_records(stored: &StoredLibrary, features: &[MbrFeature]) -> Vec<FeatureRecord> {
    features
        .iter()
        .filter(|f| f.evidence_type != EvidenceType::BelowThreshold)
        .enumerate()
        .map(|(i, f)| {
            let m = &stored.library.meta[f.precursor as usize];
            let peptidoform = split_key(&m.key).to_string();
            let (sequence, modifications) = parse_peptidoform(&peptidoform);
            let run = stored.runs[f.run as usize].clone();
            let proteins: Vec<ProteinPosition> = stored.proteins[f.precursor as usize]
                .split(';')
                .filter(|a| !a.is_empty())
                .map(|a| (a.to_string(), None, None, None, None))
                .collect();
            FeatureRecord {
                feature_id: i as i64 + 1,
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
                    quantity_se: f.quantity_se,
                }),
            }
        })
        .collect()
}

/// Transfer FDR and the MBR outputs over the extract rows of every run (in
/// run order); with `parquet_dir`, also `<parquet_dir>/quantms.feature.parquet`.
pub(crate) fn finalize_and_write(
    stored: &StoredLibrary,
    runs: &[Vec<ExtractRow>],
    mbr_fdr: f64,
    base: &Path,
    parquet_dir: Option<&Path>,
) -> Result<(), String> {
    let (tsv, validation) = mbr_paths(base);
    let summary = finalize(
        &stored.runs,
        &stored.library,
        runs,
        mbr_fdr,
        MBR_SEED,
        &tsv,
        &validation,
    )?;
    eprintln!(
        "[mbr] wrote {} ({} rows; transfer_q is an internal quant-quality q-value, not an identification FDR) and {}",
        tsv.display(),
        summary.rows_written,
        validation.display()
    );
    if let Some(dir) = parquet_dir {
        std::fs::create_dir_all(dir).map_err(|e| format!("create {}: {e}", dir.display()))?;
        let path = dir.join("quantms.feature.parquet");
        let records = feature_records(stored, &summary.features);
        output::write_feature_parquet(&path, &records)
            .map_err(|e| format!("write {}: {e}", path.display()))?;
        eprintln!(
            "[mbr] wrote {} ({} feature rows with evidence type; below-threshold transfers are in the mbr.tsv only)",
            path.display(),
            records.len()
        );
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
    mbr_fdr: f64,
    base: &Path,
    parquet_dir: Option<&Path>,
) -> Result<(), String> {
    let t0 = std::time::Instant::now();
    let stems = library_runs(dir, stems)?;
    let stored = build_and_write_library(
        dir,
        &stems,
        || Ok(PsmScores::Merged(scores)),
        max_q,
        donor_q,
    )?;
    let fingerprint = library_fingerprint(&stored);
    // One run at a time: each extraction is parallel over its precursors and
    // holds that run's MS1.
    let runs: Vec<Vec<ExtractRow>> = (0..stems.len())
        .map(|r| extract_and_write(dir, &stored, &fingerprint, r as u16, ep))
        .collect::<Result<_, _>>()?;
    finalize_and_write(&stored, &runs, mbr_fdr, base, parquet_dir)?;
    eprintln!("[PHASE mbr: {:.2}s]", t0.elapsed().as_secs_f64());
    Ok(())
}

fn run_library(a: QuantLibraryArgs) -> Result<(), String> {
    let stems = library_runs(&a.quant_dir, &a.runs)?;
    let source = score_source(&a.psms, &a.decoy_psms, &stems)?;
    build_and_write_library(
        &a.quant_dir,
        &stems,
        || load_scores(&source),
        a.mbr_library_fdr,
        a.quant_fdr,
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
    extract_and_write(&a.quant_dir, &stored, &fingerprint, run as u16, &ep).map(|_| ())
}

fn run_finalize(a: QuantFinalizeArgs) -> Result<(), String> {
    let stored = read_library(&a.quant_dir)?;
    let fingerprint = library_fingerprint(&stored);
    let mut runs = Vec::with_capacity(stored.runs.len());
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
        runs.push(read_extract_table(&path)?);
    }
    finalize_and_write(
        &stored,
        &runs,
        a.mbr_fdr,
        &a.output,
        a.output_parquet.as_deref(),
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
        QuantCommand::Finalize(a) => run_finalize(a),
    }
    .map_err(Into::into)
}

#[cfg(test)]
mod tests {
    use super::*;

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
            (2, 'C', Some("UNIMOD:4"))
        );
        assert_eq!((mods[1].position, mods[1].accession.clone()), (4, None));
        assert_eq!(split_key("LC[UNIMOD:4]AK/3"), "LC[UNIMOD:4]AK");
    }

    #[test]
    fn library_runs_checks_artifacts_and_duplicates() {
        let dir = tempfile::tempdir().unwrap();
        let d = dir.path();
        for f in ["b.quant.parquet", "b.cache", "a.quant.parquet"] {
            std::fs::write(d.join(f), b"x").unwrap();
        }
        let e = library_runs(d, &[]).unwrap_err();
        assert!(e.contains("a.cache"), "{e}");
        std::fs::write(d.join("a.cache"), b"x").unwrap();
        assert_eq!(library_runs(d, &[]).unwrap(), vec!["a", "b"]);
        assert_eq!(
            library_runs(d, &["b".into(), "a".into()]).unwrap(),
            vec!["a", "b"]
        );
        let e = library_runs(d, &["a".into(), "a".into()]).unwrap_err();
        assert!(e.contains("twice"), "{e}");
        let e = library_runs(d, &["c".into()]).unwrap_err();
        assert!(
            e.contains("c.quant.parquet") && e.contains("c.cache"),
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
        assert_eq!(m.psms.len(), 2);
        assert!(!m.psms["s1"].is_decoy && m.psms["s2"].is_decoy);
    }

    #[test]
    fn feature_records_leave_out_below_threshold_transfers() {
        let meta: Vec<PrecursorMeta> = (0..2)
            .map(|i| PrecursorMeta {
                key: format!("PEPT{}IDEK/2", ["", "M[UNIMOD:35]"][i]),
                mz: 500.0 + i as f64,
                charge: 2,
                envelope: [0.55, 0.3, 0.1, 0.05],
            })
            .collect();
        let library = build_library(&[], meta, 2, 0.01, 0.01, &AlignParams::default());
        let stored = StoredLibrary {
            runs: vec!["a".into(), "b".into()],
            max_q: 0.01,
            donor_q: 0.01,
            library,
            proteins: vec!["P1".into(), "P2".into()],
        };
        let feature = |precursor: u32, run: u16, evidence_type, transfer_q| MbrFeature {
            precursor,
            run,
            evidence_type,
            quantity: Some(1e6),
            quantity_se: Some(1e4),
            transfer_q,
            apex_rt_s: Some(100.0),
            library_q: 0.001,
            mass_err_ppm: Some(0.5),
        };
        let features = [
            feature(0, 0, EvidenceType::Identified, None),
            feature(0, 1, EvidenceType::BelowThreshold, Some(0.2)),
            feature(1, 0, EvidenceType::TransferMs1, Some(0.004)),
            feature(1, 1, EvidenceType::BelowThreshold, Some(0.5)),
        ];
        let records = feature_records(&stored, &features);
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
            "below-threshold transfers are not feature rows"
        );
        assert_eq!(records[1].peptidoform, "PEPTM[UNIMOD:35]IDEK");
    }

    #[test]
    fn residues_ignore_flanks_and_modification_notation() {
        assert_eq!(residues("K.GAC+57.02146LLPK.I"), "GACLLPK");
        assert_eq!(residues("GAC[UNIMOD:4]LLPK"), "GACLLPK");
        assert_eq!(residues("-.PEPM[+15.99491]K.-"), "PEPMK");
        assert_eq!(residues("PEPTIDEK"), "PEPTIDEK");
    }
}
