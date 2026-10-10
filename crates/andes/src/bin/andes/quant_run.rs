//! Quantification driver: collects what the search stream would otherwise
//! drop (reporter ions, MS1 scans, MS3 scans), then — after rescoring — turns
//! the confident rank-1 PSMs into quantified features and writes the tables.
//!
//! Position in the run: `spectra → search → PIN → rescore → quant → outputs`.
//! Everything above the feature level (normalization across samples, protein
//! roll-up, statistics) is left to quantms / OpenMS, which know the
//! experimental design.

use std::collections::HashMap;
use std::ops::Range;
use std::path::{Path, PathBuf};

use model::mass::{ISOTOPE, PROTON};
use model::scan::{Ms1Scan, ProductScan, ScanContext};
use model::spectrum::Spectrum;
use model::Tolerance;
use output::feature_parquet::FeatureRecord;
use output::{LfqRow, ModRecord, PercolatorPsm, QuantId, TmtRow};
use quant::formula::Formula;
use quant::isobaric::{extract_reporters, CorrectionMatrix, Plex};
use quant::lfq::{
    competition_score, feature_score, quantify, quantify_decoy, FeatureQuant, LfqParams, LfqTarget,
};
use quant::ms1_index::Ms1RunIndex;
use quant::tdc::{picked_qvalues, Pair};
use rayon::prelude::*;
use search::candidate_gen::Candidate;
use search::psm::{PsmMatch, TopNQueue};
use search::search_index::SearchIndex;

use crate::cli::Cli;

#[path = "ms1_cache.rs"]
mod ms1_cache;
use ms1_cache::Ms1Cache;

/// Isotope-peak tolerance of the precursor-purity ladder (OpenMS default).
const PURITY_ISOTOPE_PPM: f64 = 10.0;

/// MS1 resolution check for `--lfq`: isotope partners are matched at this
/// tolerance, and a run whose most intense MS1 centroids find one less often
/// than [`MS1_MIN_ISOTOPE_PARTNER_FRACTION`] is reported as low-resolution.
/// Orbitrap survey scans measured 0.70-0.93 (PXD000001, PXD001819, the OpenMS
/// BSA example); the same scans with ion-trap-sized centroid errors 0.09-0.32.
const MS1_RESOLUTION_CHECK_PPM: f64 = 10.0;
const MS1_MIN_ISOTOPE_PARTNER_FRACTION: f64 = 0.5;

/// Resolved quantification settings for a run.
#[derive(Debug, Clone)]
pub(crate) struct QuantSettings {
    pub plex: Option<Plex>,
    pub tmt_level: u8,
    pub tmt_tol: Tolerance,
    pub correction: Option<CorrectionMatrix>,
    pub min_purity: f32,
    pub lfq: bool,
    pub lfq_params: LfqParams,
    pub feature_fdr: f64,
    /// PSM q-value cut for quantification targets (applied when rescoring ran).
    pub quant_fdr: f64,
    pub pep_cap: Option<f64>,
}

impl QuantSettings {
    /// `None` when neither `--tmt` nor `--lfq` is set.
    pub(crate) fn from_cli(cli: &Cli, high_res: bool) -> Result<Option<Self>, String> {
        if cli.tmt.is_none() && !cli.lfq {
            return Ok(None);
        }
        if !(2..=3).contains(&cli.tmt_level) {
            return Err(format!("--tmt-level must be 2 or 3, got {}", cli.tmt_level));
        }
        // The default follows the scan the reporters are read from. SPS-MS3
        // reporters are read in the Orbitrap even when the identifying MS2 is
        // an ion-trap scan (Fusion/Lumos SPS methods), so the MS2 resolution
        // says nothing about them.
        let tmt_tol = cli.tmt_tol.unwrap_or(if high_res || cli.tmt_level >= 3 {
            Tolerance::Ppm(20.0)
        } else {
            Tolerance::Da(0.3)
        });
        if let Some(plex) = cli.tmt {
            plex.check_tolerance(tmt_tol).map_err(|e| {
                format!(
                    "--tmt {}: {e}. {}",
                    plex.name(),
                    if cli.tmt_tol.is_some() {
                        "Use a narrower --tmt-tol (e.g. 20ppm)."
                    } else {
                        "The MS2 is low-resolution, so its reporter ions cannot separate this \
                         kit's channels: use --tmt-level 3 for SPS-MS3 data, or set --tmt-tol."
                    }
                )
            })?;
        }
        let correction = match (&cli.tmt_correction, cli.tmt) {
            (Some(path), Some(plex)) => {
                let text = std::fs::read_to_string(path)
                    .map_err(|e| format!("--tmt-correction {}: {e}", path.display()))?;
                let m = CorrectionMatrix::parse(plex, &text)
                    .map_err(|e| format!("--tmt-correction {}: {e}", path.display()))?;
                if m.is_identity() {
                    eprintln!(
                        "WARN: --tmt-correction {} is an identity matrix; no correction applied",
                        path.display()
                    );
                    None
                } else {
                    Some(m)
                }
            }
            (Some(_), None) => return Err("--tmt-correction needs --tmt <plex>".into()),
            _ => None,
        };
        let lfq_params = LfqParams {
            tol: cli.lfq_tol.unwrap_or(Tolerance::Ppm(10.0)),
            rt_window_s: cli.lfq_rt_window,
            gaussian_width_s: cli.lfq_gaussian_width,
            min_cosine: cli.lfq_min_cosine as f32,
            ..LfqParams::default()
        };
        // Label-free quant reads the MS1 scans, whose analyzer the MS2 model
        // does not describe: Velos/Elite/Fusion/Lumos high-low methods pair an
        // ion-trap MS2 with Orbitrap survey scans. Each run's MS1 is checked
        // once it is captured (`finish_file`).
        if cli.lfq && !high_res {
            eprintln!(
                "WARN: --lfq: the MS2 scans are low-resolution; label-free quantification \
                 assumes high-resolution (Orbitrap/TOF) MS1 scans and checks them per run"
            );
        }
        Ok(Some(Self {
            plex: cli.tmt,
            tmt_level: cli.tmt_level,
            tmt_tol,
            correction,
            min_purity: cli.tmt_min_purity as f32,
            lfq: cli.lfq,
            lfq_params,
            feature_fdr: cli.lfq_feature_fdr,
            quant_fdr: cli.quant_fdr.or(cli.fdr).unwrap_or(0.01),
            pep_cap: cli.pep,
        }))
    }

    /// MS1 scans are needed for label-free quant and for the precursor purity
    /// of isobaric PSMs.
    pub(crate) fn needs_ms1(&self) -> bool {
        self.lfq || self.plex.is_some()
    }

    pub(crate) fn needs_ms3(&self) -> bool {
        self.plex.is_some() && self.tmt_level >= 3
    }

    pub(crate) fn describe(&self) -> String {
        let mut parts = Vec::new();
        if let Some(p) = self.plex {
            parts.push(format!(
                "isobaric {} at MS{} (reporter tol {}, correction {}, min purity {})",
                p.name(),
                self.tmt_level,
                match self.tmt_tol {
                    Tolerance::Ppm(v) => format!("{v} ppm"),
                    Tolerance::Da(v) => format!("{v} Da"),
                },
                if self.correction.is_some() {
                    "on"
                } else {
                    "off"
                },
                self.min_purity
            ));
        }
        if self.lfq {
            parts.push(format!(
                "label-free MS1 (tol {}, rt window ±{} s, {} isotopes, min cosine {}, feature FDR {}, smoothing {})",
                match self.lfq_params.tol {
                    Tolerance::Ppm(v) => format!("{v} ppm"),
                    Tolerance::Da(v) => format!("{v} Da"),
                },
                self.lfq_params.rt_window_s,
                self.lfq_params.n_isotopes,
                self.lfq_params.min_cosine,
                self.feature_fdr,
                match self.lfq_params.gaussian_width_s {
                    Some(width) => format!("Gaussian {width} s"),
                    None => "Savitzky–Golay 5 scans".into(),
                }
            ));
        }
        parts.join("; ")
    }
}

/// Scans a reader captured for one input file besides the MS2 stream.
#[derive(Debug, Default)]
pub(crate) struct RunScans {
    pub ms1: Vec<Ms1Scan>,
    pub ms3: Vec<ProductScan>,
    /// Acquisition context of each emitted MS2 (survey scan, m/z range, FAIMS).
    pub contexts: Vec<ScanContext>,
}

/// One input file's quantification state.
struct RunCapture {
    stem: String,
    /// `"<stem>/"` under multi-file searches (spectrum titles carry it).
    title_prefix: Option<String>,
    /// Indices of this file's spectra in the global spectrum list.
    span: Range<usize>,
    ms1: Option<Ms1Cache>,
    /// Parent MS2 native id → the reporter scan read for it.
    ms3_reporters: HashMap<String, ReporterScan>,
    /// MS2 native id → its acquisition context.
    contexts: HashMap<String, ScanContext>,
}

/// An MS3 reporter scan: native id, raw reporter vector, acquired m/z range.
struct ReporterScan {
    id: String,
    raw: Vec<f32>,
    scan_window: Option<(f64, f64)>,
}

/// Collects quantification inputs while the search streams the spectra.
pub(crate) struct QuantCollector {
    pub settings: QuantSettings,
    /// Raw MS2 reporter vectors, one per pushed spectrum (isobaric only).
    ms2_reporters: Vec<Vec<f32>>,
    runs: Vec<RunCapture>,
    next_span_start: usize,
    /// False when the read path does not capture MS1 scans for quantification
    /// (`--chimeric`); the missing purity was announced once at setup.
    ms1_captured: bool,
}

impl QuantCollector {
    pub(crate) fn new(settings: QuantSettings) -> Self {
        Self {
            settings,
            ms2_reporters: Vec::new(),
            runs: Vec::new(),
            next_span_start: 0,
            ms1_captured: true,
        }
    }

    /// The read path does not capture MS1 scans (`--chimeric`): skip the
    /// per-file warning about them.
    pub(crate) fn ms1_not_captured(&mut self) {
        self.ms1_captured = false;
    }

    /// Call for every spectrum right before it is pushed to the global list
    /// (while its peaks are still present).
    pub(crate) fn on_spectrum(&mut self, spec: &Spectrum) {
        if let Some(plex) = self.settings.plex {
            self.ms2_reporters
                .push(extract_reporters(&spec.peaks, plex, self.settings.tmt_tol));
        }
    }

    /// Call after a file's spectra have all been pushed (`end` = global count).
    pub(crate) fn finish_file(
        &mut self,
        path: &Path,
        title_prefix: Option<&str>,
        scans: RunScans,
        end: usize,
    ) -> Result<(), String> {
        let stem = run_stem(path);
        let ms1 = if self.settings.needs_ms1() && !scans.ms1.is_empty() {
            Some(Ms1RunIndex::new(scans.ms1))
        } else {
            None
        };
        let mut ms3_reporters: HashMap<String, ReporterScan> = HashMap::new();
        if let Some(plex) = self.settings.plex {
            if self.settings.tmt_level >= 3 {
                for p in scans.ms3 {
                    let Some(parent) = p.parent_id else { continue };
                    let raw = extract_reporters(&p.peaks, plex, self.settings.tmt_tol);
                    // Keep the most intense MS3 if several reference the same MS2.
                    let total: f32 = raw.iter().sum();
                    match ms3_reporters.get(&parent) {
                        Some(old) if old.raw.iter().sum::<f32>() >= total => {}
                        _ => {
                            ms3_reporters.insert(
                                parent,
                                ReporterScan {
                                    id: p.id,
                                    raw,
                                    scan_window: p.scan_window,
                                },
                            );
                        }
                    }
                }
            }
        }
        let contexts: HashMap<String, ScanContext> = scans
            .contexts
            .into_iter()
            .map(|c| (c.id.clone(), c))
            .collect();
        let span = self.next_span_start..end;
        self.next_span_start = end;
        if let (true, Some(idx)) = (self.settings.lfq, &ms1) {
            if let Some(frac) = idx.isotope_partner_fraction(MS1_RESOLUTION_CHECK_PPM) {
                if frac < MS1_MIN_ISOTOPE_PARTNER_FRACTION {
                    eprintln!(
                        "WARN: quant: {stem}: the MS1 scans do not look high-resolution \
                         ({:.0}% of their most intense centroids have an isotope partner within \
                         {MS1_RESOLUTION_CHECK_PPM} ppm; Orbitrap/TOF survey scans reach 70-95%). \
                         Label-free areas from low-resolution MS1 are not meaningful.",
                        100.0 * frac
                    );
                }
            }
        }
        if let Some(idx) = &ms1 {
            eprintln!(
                "quant: {}: {} MS1 scans ({} centroids) indexed{}",
                stem,
                idx.len(),
                idx.total_peaks(),
                if self.settings.tmt_level >= 3 && self.settings.plex.is_some() {
                    format!(", {} MS3 scans linked to an MS2", ms3_reporters.len())
                } else {
                    String::new()
                }
            );
        } else if self.settings.needs_ms1() && self.ms1_captured {
            eprintln!(
                "WARN: quant: {}: no MS1 scans captured ({}); {}",
                stem,
                path.display(),
                if self.settings.lfq {
                    "label-free quantification is skipped for this file"
                } else {
                    "precursor purity is unavailable for this file"
                }
            );
        }
        let ms1 = ms1
            .as_ref()
            .map(Ms1Cache::store)
            .transpose()
            .map_err(|e| format!("cache MS1 scans for {stem}: {e}"))?;
        self.runs.push(RunCapture {
            stem,
            title_prefix: title_prefix.map(str::to_string),
            span,
            ms1,
            ms3_reporters,
            contexts,
        });
        Ok(())
    }

    fn run_of(&self, spec_idx: usize) -> Option<usize> {
        self.runs.iter().position(|r| r.span.contains(&spec_idx))
    }
}

/// A rank-1 target PSM selected for quantification.
struct Hit<'a> {
    spec_idx: usize,
    run_idx: usize,
    psm: PsmMatch,
    cand: &'a Candidate,
    spec_id: String,
    q_value: Option<f64>,
    pep: Option<f64>,
    psm_index: i64,
}

/// 1-based ranks over a rank-sorted PSM list; ties share a rank (mirrors the
/// PIN writer, so SpecIds match the Percolator join key).
fn ranks_of(psms: &[PsmMatch]) -> Vec<u32> {
    let mut out = Vec::with_capacity(psms.len());
    let mut rank = 0u32;
    let mut prev: Option<f32> = None;
    for p in psms {
        let ties = prev.is_some_and(|q| p.rank_score == q || (p.rank_score.is_nan() && q.is_nan()));
        if !ties {
            rank += 1;
            prev = Some(p.rank_score);
        }
        out.push(rank);
    }
    out
}

/// Select the quantification targets: rank-1 target PSMs, at `q <= quant_fdr`
/// (and `pep <= --pep`) when a rescoring result is available. The running
/// `psm_index` reproduces the row order of `psms.parquet`.
fn select_hits<'a>(
    collector: &QuantCollector,
    spectra: &[Spectrum],
    queues: &[TopNQueue],
    candidates: &'a [Candidate],
    rescore: Option<&HashMap<String, PercolatorPsm>>,
) -> Vec<Hit<'a>> {
    let s = &collector.settings;
    let mut hits = Vec::new();
    let mut psm_index: i64 = 0;
    for (spec_idx, queue) in queues.iter().enumerate() {
        if queue.is_empty() {
            continue;
        }
        let spec = &spectra[spec_idx];
        let psms = queue.clone().into_rank_sorted_vec();
        let ranks = ranks_of(&psms);
        let multi_row = psms.len() > 1;
        let scan = spec.scan.unwrap_or(0);
        let spec_id_base = if spec.title.is_empty() {
            format!("scan={scan}")
        } else {
            spec.title.clone()
        };
        for (hit, psm) in psms.into_iter().enumerate() {
            let this_index = psm_index;
            psm_index += 1;
            if ranks[hit] != 1 {
                continue;
            }
            let cand = &candidates[psm.primary_candidate_idx() as usize];
            if cand.is_decoy {
                continue;
            }
            let spec_id = output::format_spec_id(&spec_id_base, scan, ranks[hit], hit, multi_row);
            let (q_value, pep) = match rescore {
                Some(m) => match m.get(&spec_id) {
                    Some(r) => {
                        if r.q_value > s.quant_fdr || s.pep_cap.is_some_and(|c| r.pep > c) {
                            continue;
                        }
                        (Some(r.q_value), Some(r.pep))
                    }
                    None => continue,
                },
                None => (None, None),
            };
            let Some(run_idx) = collector.run_of(spec_idx) else {
                continue;
            };
            hits.push(Hit {
                spec_idx,
                run_idx,
                psm,
                cand,
                spec_id,
                q_value,
                pep,
                psm_index: this_index,
            });
        }
    }
    hits
}

fn peptidoform_of(cand: &Candidate) -> String {
    output::qpx::peptidoform_string(cand)
}

fn modifications_of(cand: &Candidate) -> Vec<ModRecord> {
    cand.peptide
        .residues
        .iter()
        .enumerate()
        .filter_map(|(i, aa)| {
            aa.mod_.as_ref().map(|m| ModRecord {
                name: m.name.clone(),
                accession: m.accession.clone().filter(|a| !a.is_empty()),
                position: i as i32 + 1,
                amino_acid: aa.residue as char,
            })
        })
        .collect()
}

/// Accessions of every protein the PSM's candidates come from, with the
/// peptide's position in each (`start`, `end`, `pre`, `post`).
type ProteinPos = output::ProteinPosition;

fn proteins_of(psm: &PsmMatch, candidates: &[Candidate], index: &SearchIndex) -> Vec<ProteinPos> {
    let mut out: Vec<ProteinPos> = Vec::new();
    for &ci in &psm.candidate_idxs {
        let c = &candidates[ci as usize];
        let acc = match index.protein_at(c.protein_index) {
            Some(p) => p.accession.clone(),
            None => format!("PROT_{}", c.protein_index),
        };
        if out.iter().any(|(a, ..)| *a == acc) {
            continue;
        }
        let start = c.start_offset_in_protein as i32 + 1;
        let end = start + c.peptide.length() as i32 - 1;
        out.push((
            acc,
            Some(start),
            Some(end),
            Some((c.peptide.pre as char).to_string()),
            Some((c.peptide.post as char).to_string()),
        ));
    }
    out
}

fn quant_id(
    hit: &Hit<'_>,
    run: &RunCapture,
    spec: &Spectrum,
    candidates: &[Candidate],
    index: &SearchIndex,
) -> QuantId {
    let z = hit.psm.charge_used.max(1) as f64;
    let precursor_mz = hit.psm.precursor_mz_override.unwrap_or(spec.precursor_mz);
    let observed_mz = precursor_mz - ISOTOPE * (hit.psm.isotope_offset as f64) / z;
    let calculated_mz = hit.cand.peptide.mass() / z + PROTON;
    QuantId {
        file_idx: hit.run_idx,
        run: run.stem.clone(),
        spec_id: hit.spec_id.clone(),
        scan: spec.scan.unwrap_or(0),
        rt_seconds: spec.rt_seconds,
        sequence: output::qpx::bare_sequence(hit.cand),
        peptidoform: peptidoform_of(hit.cand),
        modifications: modifications_of(hit.cand),
        charge: hit.psm.charge_used,
        proteins: proteins_of(&hit.psm, candidates, index)
            .into_iter()
            .map(|(a, ..)| a)
            .collect(),
        is_decoy: hit.cand.is_decoy,
        calculated_mz,
        observed_mz,
        q_value: hit.q_value,
        pep: hit.pep,
        missed_cleavages: None,
        psm_index: hit.psm_index,
    }
}

fn base_record(id: &QuantId, proteins: Vec<ProteinPos>, feature_id: i64) -> FeatureRecord {
    let mass_error_ppm = if id.calculated_mz > 0.0 {
        Some(((id.observed_mz - id.calculated_mz) / id.calculated_mz * 1e6) as f32)
    } else {
        None
    };
    let mut additional_scores: Vec<(String, f64, bool)> = Vec::new();
    if let Some(q) = id.q_value {
        additional_scores.push(("psm_q_value".into(), q, false));
    }
    FeatureRecord {
        feature_id,
        sequence: id.sequence.clone(),
        peptidoform: id.peptidoform.clone(),
        modifications: id.modifications.clone(),
        charge: id.charge as i16,
        pep: id.pep,
        is_decoy: id.is_decoy,
        calculated_mz: id.calculated_mz as f32,
        observed_mz: id.observed_mz as f32,
        mass_error_ppm,
        additional_scores,
        run_file_name: id.run.clone(),
        cv_params: Vec::new(),
        // No scan number (e.g. `spectrum=N` native ids): an empty list, as
        // OpenMS writes it, not a made-up scan 0.
        scan: if id.scan > 0 {
            vec![id.scan]
        } else {
            Vec::new()
        },
        rt: id.rt_seconds.map(|t| t as f32),
        missed_cleavages: id.missed_cleavages,
        intensities: Vec::new(),
        additional_intensities: Vec::new(),
        proteins,
        id_run_file_name: Some(id.run.clone()),
        rt_start: None,
        rt_stop: None,
        psm_ids: vec![id.psm_index],
    }
}

/// A feature `lfq.tsv` reports: it won its decoy competition at the feature
/// FDR and its isotope envelope matches theory (`--lfq-min-cosine`).
/// `lfq_features.tsv` and the parquet keep every feature.
fn is_confident(row: &LfqRow, s: &QuantSettings) -> bool {
    row.feature_q_value.is_some_and(|q| q <= s.feature_fdr)
        && row.feature.cosine >= s.lfq_params.min_cosine
}

/// Run name of an input path: the file name without its extension, and
/// without a `.gz` before it (`run1.mzML.gz` → `run1`).
fn run_stem(path: &Path) -> String {
    let name = path
        .file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| path.display().to_string());
    let name = name.strip_suffix(".gz").map(str::to_string).unwrap_or(name);
    match name.rfind('.') {
        Some(i) if i > 0 => name[..i].to_string(),
        _ => name,
    }
}

/// Strip the multi-file title prefix to get the native id MS3 parents use.
fn native_id<'a>(run: &RunCapture, title: &'a str) -> &'a str {
    match &run.title_prefix {
        Some(p) => title.strip_prefix(p.as_str()).unwrap_or(title),
        None => title,
    }
}

/// The survey scan an MS2 was triggered from: the scan its precursor
/// references, else the last MS1 at or before it at the same FAIMS voltage,
/// else the last MS1 at or before it. On FAIMS runs the last MS1 in time is
/// usually a survey scan at another voltage, which does not hold the precursor.
fn survey_scan(ms1: &Ms1RunIndex, context: Option<&ScanContext>, rt: f64) -> Option<usize> {
    const MAX_BACKTRACK: usize = 64;
    if let Some(i) = context
        .and_then(|c| c.parent_id.as_deref())
        .and_then(|id| ms1.position_of(id))
    {
        return Some(i);
    }
    let last = ms1.scan_at_or_before(rt)?;
    let Some(cv) = context.and_then(|c| c.faims_cv) else {
        return Some(last);
    };
    (last.saturating_sub(MAX_BACKTRACK)..=last)
        .rev()
        .find(|&i| ms1.scan(i).faims_cv == Some(cv))
}

/// Whether a reporter scan acquired every channel of `plex`: a scan window
/// that starts above the lowest reporter (or ends below the highest) reads
/// zeros for channels that were never measured.
fn covers_reporters(scan_window: Option<(f64, f64)>, plex: Plex, tol: Tolerance) -> bool {
    let Some((lo, hi)) = scan_window else {
        return true;
    };
    let channels = plex.channels();
    let first = channels[0].mz;
    let last = channels[channels.len() - 1].mz;
    lo <= first - tol.as_da(first) && hi >= last + tol.as_da(last)
}

/// A target, its feature and its decoy twin's feature.
type Quantified = (usize, Option<FeatureQuant>, Option<FeatureQuant>);

/// One chromatographic series per FAIMS compensation voltage (the whole run
/// when it has at most one voltage).
fn faims_groups(ms1: Ms1RunIndex) -> Vec<(Option<f32>, Ms1RunIndex)> {
    if ms1.faims_voltages().len() > 1 {
        ms1.split_by_faims()
            .into_iter()
            .map(|(cv, idx)| (Some(cv), idx))
            .collect()
    } else {
        vec![(None, ms1)]
    }
}

/// The series a PSM's precursor is extracted from: its MS2's FAIMS voltage.
fn group_of(
    groups: &[(Option<f32>, Ms1RunIndex)],
    run: &RunCapture,
    spec: &Spectrum,
) -> Option<usize> {
    match groups.first()?.0 {
        None => Some(0),
        Some(_) => {
            let cv = run
                .contexts
                .get(native_id(run, &spec.title))
                .and_then(|c| c.faims_cv)?;
            groups.iter().position(|(g, _)| *g == Some(cv))
        }
    }
}

/// The label-free target of a PSM's peptidoform at `anchor_rt`.
fn lfq_target(h: &Hit<'_>, anchor_rt: f64, n_iso: usize) -> LfqTarget {
    let z = h.psm.charge_used.max(1);
    let formula = Formula::peptide(
        h.cand
            .peptide
            .residues
            .iter()
            .map(|aa| (aa.residue, aa.mod_.as_ref().map_or(0.0, |m| m.mass_delta))),
    );
    LfqTarget {
        mono_mz: h.cand.peptide.mass() / z as f64 + PROTON,
        charge: z,
        anchor_rt,
        envelope: formula.envelope(n_iso),
    }
}

fn quantify_targets(
    groups: &[(Option<f32>, Ms1RunIndex)],
    targets: &[(usize, usize, LfqTarget)],
    params: &LfqParams,
) -> Vec<Quantified> {
    targets
        .par_iter()
        .map(|(hi, g, t)| {
            let ms1 = &groups[*g].1;
            (
                *hi,
                quantify(ms1, t, params),
                quantify_decoy(ms1, t, params),
            )
        })
        .collect()
}

/// The height scale of the feature score: the highest apex over targets and
/// decoys alike.
fn max_apex_of(quantified: &[Quantified]) -> f32 {
    quantified
        .iter()
        .flat_map(|(_, f, d)| [f.as_ref(), d.as_ref()])
        .flatten()
        .map(|f| f.apex_intensity)
        .fold(0.0f32, f32::max)
}

fn competition_pairs(quantified: &[Quantified], params: &LfqParams, max_apex: f32) -> Vec<Pair> {
    quantified
        .iter()
        .map(|(_, f, d)| Pair {
            target: f
                .as_ref()
                .and_then(|f| competition_score(f, params, max_apex)),
            decoy: d
                .as_ref()
                .and_then(|f| competition_score(f, params, max_apex)),
        })
        .collect()
}

/// Run quantification over the search result and write the tables.
#[allow(clippy::too_many_arguments)]
pub(crate) fn run_quant(
    collector: &QuantCollector,
    spectra: &[Spectrum],
    queues: &[TopNQueue],
    candidates: &[Candidate],
    index: &SearchIndex,
    rescore: Option<&HashMap<String, PercolatorPsm>>,
    report_base: &Path,
    parquet_dir: Option<&Path>,
) -> Result<(), String> {
    let s = &collector.settings;
    let t0 = std::time::Instant::now();
    let hits = select_hits(collector, spectra, queues, candidates, rescore);
    eprintln!(
        "quant: {} rank-1 target PSMs selected{} ({})",
        hits.len(),
        match rescore {
            Some(_) => format!(" at q<={}", s.quant_fdr),
            None => " (no rescoring: every rank-1 target PSM; filter downstream)".to_string(),
        },
        s.describe()
    );
    let stem = report_base
        .file_stem()
        .map(|x| x.to_string_lossy().into_owned())
        .unwrap_or_else(|| "andes".to_string());
    let dir: PathBuf = report_base
        .parent()
        .map(Path::to_path_buf)
        .unwrap_or_else(|| ".".into());
    // Resolve output rows by their unique PSM row index, not by repeatedly
    // scanning every selected hit (quadratic, and SpecIds can be ambiguous).
    let hits_by_psm: HashMap<i64, &Hit<'_>> = if parquet_dir.is_some() {
        hits.iter().map(|h| (h.psm_index, h)).collect()
    } else {
        HashMap::new()
    };
    let mut records: Vec<FeatureRecord> = Vec::new();
    let mut next_feature_id: i64 = 1;

    // ── Isobaric ──────────────────────────────────────────────────────────────
    if let Some(plex) = s.plex {
        let mut rows: Vec<TmtRow> = Vec::new();
        let (mut no_scan, mut low_purity, mut out_of_range) = (0usize, 0usize, 0usize);
        for (run_idx, run) in collector.runs.iter().enumerate() {
            let ms1 = run
                .ms1
                .as_ref()
                .map(Ms1Cache::load)
                .transpose()
                .map_err(|e| format!("load MS1 scans for {}: {e}", run.stem))?;
            for hit in hits.iter().filter(|h| h.run_idx == run_idx) {
                let spec = &spectra[hit.spec_idx];
                let context = run.contexts.get(native_id(run, &spec.title));
                let (quant_scan_id, raw, scan_window) = if s.tmt_level >= 3 {
                    match run.ms3_reporters.get(native_id(run, &spec.title)) {
                        Some(r) => (r.id.clone(), r.raw.clone(), r.scan_window),
                        None => {
                            no_scan += 1;
                            continue;
                        }
                    }
                } else {
                    match collector.ms2_reporters.get(hit.spec_idx) {
                        Some(raw) => (
                            spec.title.clone(),
                            raw.clone(),
                            context.and_then(|c| c.scan_window),
                        ),
                        None => {
                            no_scan += 1;
                            continue;
                        }
                    }
                };
                if !covers_reporters(scan_window, plex, s.tmt_tol) {
                    out_of_range += 1;
                    continue;
                }
                let purity = ms1.as_ref().and_then(|ms1| {
                    let rt = spec.rt_seconds?;
                    let scan = survey_scan(ms1, context, rt)?;
                    let lo = spec.isolation_lower_offset?;
                    let hi = spec.isolation_upper_offset?;
                    let pmz = hit.psm.precursor_mz_override.unwrap_or(spec.precursor_mz);
                    quant::precursor_purity(
                        ms1.peaks(scan),
                        pmz,
                        hit.psm.charge_used,
                        lo,
                        hi,
                        PURITY_ISOTOPE_PPM,
                    )
                });
                if let Some(p) = purity {
                    if p < s.min_purity {
                        low_purity += 1;
                    }
                }
                let corrected = match &s.correction {
                    Some(m) => m.correct(&raw),
                    None => raw.clone(),
                };
                rows.push(TmtRow {
                    id: quant_id(hit, run, spec, candidates, index),
                    quant_scan_id,
                    quant_ms_level: s.tmt_level,
                    purity,
                    raw,
                    corrected,
                });
            }
        }
        let tsv_rows: Vec<TmtRow> = rows
            .iter()
            .filter(|r| r.purity.is_none_or(|p| p >= s.min_purity))
            .cloned()
            .collect();
        let tmt_path = dir.join(format!("{stem}.tmt.tsv"));
        output::write_tmt_tsv(&tmt_path, plex, &tsv_rows, s.correction.is_some())
            .map_err(|e| format!("write {}: {e}", tmt_path.display()))?;
        eprintln!(
            "quant: wrote {} ({} PSMs with reporter ions{}{}{})",
            tmt_path.display(),
            tsv_rows.len(),
            if no_scan > 0 {
                format!(", {no_scan} without a quantification scan")
            } else {
                String::new()
            },
            if out_of_range > 0 {
                format!(
                    ", {out_of_range} whose scan range does not cover the reporter ions \
                     (not quantified)"
                )
            } else {
                String::new()
            },
            if low_purity > 0 {
                format!(
                    ", {low_purity} below purity {} kept only in the parquet",
                    s.min_purity
                )
            } else {
                String::new()
            }
        );
        for r in rows.iter().filter(|_| parquet_dir.is_some()) {
            let hit = hits_by_psm[&r.id.psm_index];
            let proteins = proteins_of(&hit.psm, candidates, index);
            let mut rec = base_record(&r.id, proteins, next_feature_id);
            next_feature_id += 1;
            let labels: Vec<String> = plex
                .channels()
                .iter()
                .map(|c| format!("{}{}", plex.label_prefix(), c.name))
                .collect();
            rec.intensities = labels
                .iter()
                .cloned()
                .zip(r.corrected.iter().copied())
                .collect();
            if s.correction.is_some() {
                rec.additional_intensities = labels
                    .iter()
                    .cloned()
                    .zip(r.raw.iter().copied())
                    .map(|(l, v)| (l, vec![("raw".to_string(), v)]))
                    .collect();
            }
            if let Some(p) = r.purity {
                rec.additional_scores
                    .push(("precursor_purity".into(), p as f64, true));
            }
            rec.cv_params
                .push(("quant_scan".into(), r.quant_scan_id.clone()));
            rec.cv_params
                .push(("quant_ms_level".into(), r.quant_ms_level.to_string()));
            records.push(rec);
        }
    }

    // ── Label-free ────────────────────────────────────────────────────────────
    if s.lfq {
        let n_iso = s.lfq_params.n_isotopes;
        // Best PSM per (run, peptidoform, charge): lowest q, then highest score.
        let mut best: HashMap<(usize, String, u8), usize> = HashMap::new();
        for (i, h) in hits.iter().enumerate() {
            let key = (h.run_idx, peptidoform_of(h.cand), h.psm.charge_used);
            match best.get(&key) {
                Some(&j) => {
                    let other = &hits[j];
                    let better = match (h.q_value, other.q_value) {
                        (Some(a), Some(b)) if a != b => a < b,
                        _ => h.psm.rank_score > other.psm.rank_score,
                    };
                    if better {
                        best.insert(key, i);
                    }
                }
                None => {
                    best.insert(key, i);
                }
            }
        }
        let params = &s.lfq_params;
        let mut all_rows: Vec<LfqRow> = Vec::new();
        for (run_idx, run) in collector.runs.iter().enumerate() {
            let Some(cache) = run.ms1.as_ref() else {
                continue;
            };
            let ms1 = cache
                .load()
                .map_err(|e| format!("load MS1 scans for {}: {e}", run.stem))?;
            let groups = faims_groups(ms1);
            let mut targets: Vec<(usize, usize, LfqTarget)> = best
                .iter()
                .filter(|((r, _, _), _)| *r == run_idx)
                .filter_map(|(_, &hi)| {
                    let h = &hits[hi];
                    let spec = &spectra[h.spec_idx];
                    let group = group_of(&groups, run, spec)?;
                    Some((hi, group, lfq_target(h, spec.rt_seconds?, n_iso)))
                })
                .collect();
            // HashMap iteration must not determine feature row/ID order.
            targets.sort_unstable_by_key(|(hi, _, _)| *hi);
            let quantified = quantify_targets(&groups, &targets, params);
            let max_apex = max_apex_of(&quantified);
            let qvals = picked_qvalues(&competition_pairs(&quantified, params, max_apex));
            let mut n_found = 0usize;
            for (k, (hi, feature, decoy)) in quantified.into_iter().enumerate() {
                let Some(feature) = feature else { continue };
                n_found += 1;
                let hit = &hits[hi];
                let spec = &spectra[hit.spec_idx];
                all_rows.push(LfqRow {
                    id: quant_id(hit, run, spec, candidates, index),
                    score: feature_score(&feature, params.rt_window_s, max_apex),
                    feature_q_value: qvals[k],
                    decoy_score: decoy
                        .as_ref()
                        .map(|d| feature_score(d, params.rt_window_s, max_apex)),
                    feature,
                });
            }
            eprintln!(
                "quant: {}: {} of {} precursors quantified, {} at feature q<={} and envelope cosine>={}",
                run.stem,
                n_found,
                targets.len(),
                all_rows
                    .iter()
                    .filter(|r| r.id.file_idx == run_idx && is_confident(r, s))
                    .count(),
                s.feature_fdr,
                s.lfq_params.min_cosine
            );
        }
        let runs: Vec<String> = collector.runs.iter().map(|r| r.stem.clone()).collect();
        let confident: Vec<LfqRow> = all_rows
            .iter()
            .filter(|r| is_confident(r, s))
            .cloned()
            .collect();
        let lfq_path = dir.join(format!("{stem}.lfq.tsv"));
        output::write_lfq_tsv(&lfq_path, &runs, &confident)
            .map_err(|e| format!("write {}: {e}", lfq_path.display()))?;
        let feat_path = dir.join(format!("{stem}.lfq_features.tsv"));
        output::write_lfq_features_tsv(&feat_path, &all_rows)
            .map_err(|e| format!("write {}: {e}", feat_path.display()))?;
        eprintln!(
            "quant: wrote {} ({} accepted run-level features at feature q<={}, cosine>={}) and {} ({} features)",
            lfq_path.display(),
            confident.len(),
            s.feature_fdr,
            s.lfq_params.min_cosine,
            feat_path.display(),
            all_rows.len()
        );
        for r in all_rows.iter().filter(|_| parquet_dir.is_some()) {
            let hit = hits_by_psm[&r.id.psm_index];
            let proteins = proteins_of(&hit.psm, candidates, index);
            let mut rec = base_record(&r.id, proteins, next_feature_id);
            next_feature_id += 1;
            let f = &r.feature;
            rec.rt = Some(f.apex_rt as f32);
            rec.rt_start = Some(f.rt_start as f32);
            rec.rt_stop = Some(f.rt_stop as f32);
            rec.intensities = vec![("LFQ".to_string(), f.area as f32)];
            let mut extra = vec![("apex".to_string(), f.apex_intensity)];
            for (k, a) in f.isotope_areas.iter().enumerate() {
                extra.push((format!("isotope_{k}_area"), *a as f32));
            }
            rec.additional_intensities = vec![("LFQ".to_string(), extra)];
            rec.additional_scores
                .push(("envelope_cosine".into(), f.cosine as f64, true));
            rec.additional_scores
                .push(("rt_delta_s".into(), f.rt_delta_s, false));
            rec.additional_scores
                .push(("n_isotopes".into(), f.n_isotopes_found as f64, true));
            rec.additional_scores
                .push(("feature_score".into(), r.score as f64, true));
            if let Some(q) = r.feature_q_value {
                rec.additional_scores
                    .push(("feature_q_value".into(), q, false));
            }
            if let Some(d) = r.decoy_score {
                rec.additional_scores
                    .push(("decoy_feature_score".into(), d as f64, true));
            }
            rec.cv_params.push((
                "psm_rt".into(),
                format!("{:.3}", r.id.rt_seconds.unwrap_or(0.0)),
            ));
            records.push(rec);
        }
    }

    if let Some(pdir) = parquet_dir {
        std::fs::create_dir_all(pdir).map_err(|e| format!("create {}: {e}", pdir.display()))?;
        let path = pdir.join("quantms.feature.parquet");
        output::write_feature_parquet(&path, &records)
            .map_err(|e| format!("write {}: {e}", path.display()))?;
        eprintln!(
            "quant: wrote {} ({} feature rows)",
            path.display(),
            records.len()
        );
    }
    eprintln!("[PHASE quant: {:.2}s]", t0.elapsed().as_secs_f64());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn faims_run() -> Ms1RunIndex {
        // Survey scans alternate between two voltages, one per second.
        Ms1RunIndex::new(
            (0..6)
                .map(|i| Ms1Scan {
                    rt: i as f64,
                    peaks: Vec::new(),
                    id: format!("scan={}", i + 1),
                    faims_cv: Some(if i % 2 == 0 { -40.0 } else { -60.0 }),
                })
                .collect(),
        )
    }

    fn context(parent: Option<&str>, cv: Option<f32>) -> ScanContext {
        ScanContext {
            id: "scan=99".into(),
            parent_id: parent.map(str::to_string),
            scan_window: None,
            faims_cv: cv,
        }
    }

    #[test]
    fn survey_scan_follows_the_precursor_reference_then_the_voltage() {
        let ms1 = faims_run();
        // An MS2 at 3.5 s whose precursor references scan=3 (rt 2, -40 V).
        let c = context(Some("scan=3"), Some(-40.0));
        assert_eq!(survey_scan(&ms1, Some(&c), 3.5), Some(2));
        // No reference: the last scan at the MS2's voltage, not the last in time.
        let c = context(None, Some(-40.0));
        assert_eq!(survey_scan(&ms1, Some(&c), 3.5), Some(2));
        // An unknown reference falls back the same way.
        let c = context(Some("scan=1000"), Some(-60.0));
        assert_eq!(survey_scan(&ms1, Some(&c), 3.5), Some(3));
        // No context at all: the last scan in time.
        assert_eq!(survey_scan(&ms1, None, 4.5), Some(4));
        assert_eq!(survey_scan(&ms1, None, -1.0), None);
    }

    #[test]
    fn reporter_coverage_needs_the_whole_channel_range() {
        let tol = Tolerance::Ppm(20.0);
        assert!(covers_reporters(None, Plex::Tmt16, tol));
        assert!(covers_reporters(Some((100.0, 2000.0)), Plex::Tmt16, tol));
        assert!(covers_reporters(Some((126.0, 2000.0)), Plex::Tmt16, tol));
        // An auto scan range that starts above 126 never measured the low channels.
        assert!(!covers_reporters(Some((130.0, 2000.0)), Plex::Tmt16, tol));
        assert!(!covers_reporters(Some((148.0, 1866.0)), Plex::Tmt10, tol));
        assert!(!covers_reporters(Some((100.0, 134.0)), Plex::Tmt16, tol));
    }
}
