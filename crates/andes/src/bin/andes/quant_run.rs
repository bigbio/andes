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
use model::scan::{Ms1Scan, ProductScan};
use model::spectrum::Spectrum;
use model::Tolerance;
use output::feature_parquet::FeatureRecord;
use output::{LfqRow, ModRecord, PercolatorPsm, QuantId, TmtRow};
use quant::formula::Formula;
use quant::isobaric::{extract_reporters, CorrectionMatrix, Plex};
use quant::lfq::{feature_score, quantify, quantify_decoy, LfqParams, LfqTarget};
use quant::ms1_index::Ms1RunIndex;
use quant::tdc::{picked_qvalues, Pair};
use rayon::prelude::*;
use search::candidate_gen::Candidate;
use search::psm::{PsmMatch, TopNQueue};
use search::search_index::SearchIndex;

use crate::cli::Cli;

/// Isotope-peak tolerance of the precursor-purity ladder (OpenMS default).
const PURITY_ISOTOPE_PPM: f64 = 10.0;

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
        let tmt_tol = cli.tmt_tol.unwrap_or(if high_res {
            Tolerance::Ppm(20.0)
        } else {
            Tolerance::Da(0.3)
        });
        if !(2..=3).contains(&cli.tmt_level) {
            return Err(format!("--tmt-level must be 2 or 3, got {}", cli.tmt_level));
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
            min_cosine: cli.lfq_min_cosine as f32,
            ..LfqParams::default()
        };
        if cli.lfq && !high_res {
            return Err(
                "--lfq needs high-resolution MS1 (Orbitrap/TOF); the resolved model is low-resolution"
                    .into(),
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
                "label-free MS1 (tol {}, rt window ±{} s, {} isotopes, min cosine {}, feature FDR {})",
                match self.lfq_params.tol {
                    Tolerance::Ppm(v) => format!("{v} ppm"),
                    Tolerance::Da(v) => format!("{v} Da"),
                },
                self.lfq_params.rt_window_s,
                self.lfq_params.n_isotopes,
                self.lfq_params.min_cosine,
                self.feature_fdr
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
}

/// One input file's quantification state.
struct RunCapture {
    stem: String,
    /// `"<stem>/"` under multi-file searches (spectrum titles carry it).
    title_prefix: Option<String>,
    /// Indices of this file's spectra in the global spectrum list.
    span: Range<usize>,
    ms1: Option<Ms1RunIndex>,
    /// Parent MS2 native id → (MS3 native id, raw reporter vector).
    ms3_reporters: HashMap<String, (String, Vec<f32>)>,
}

/// Collects quantification inputs while the search streams the spectra.
pub(crate) struct QuantCollector {
    pub settings: QuantSettings,
    /// Raw MS2 reporter vectors, one per pushed spectrum (isobaric only).
    ms2_reporters: Vec<Vec<f32>>,
    runs: Vec<RunCapture>,
    next_span_start: usize,
}

impl QuantCollector {
    pub(crate) fn new(settings: QuantSettings) -> Self {
        Self {
            settings,
            ms2_reporters: Vec::new(),
            runs: Vec::new(),
            next_span_start: 0,
        }
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
    ) {
        let stem = run_stem(path);
        let ms1 = if self.settings.needs_ms1() && !scans.ms1.is_empty() {
            Some(Ms1RunIndex::new(scans.ms1))
        } else {
            None
        };
        let mut ms3_reporters: HashMap<String, (String, Vec<f32>)> = HashMap::new();
        if let Some(plex) = self.settings.plex {
            if self.settings.tmt_level >= 3 {
                for p in scans.ms3 {
                    let Some(parent) = p.parent_id else { continue };
                    let raw = extract_reporters(&p.peaks, plex, self.settings.tmt_tol);
                    // Keep the most intense MS3 if several reference the same MS2.
                    let total: f32 = raw.iter().sum();
                    match ms3_reporters.get(&parent) {
                        Some((_, old)) if old.iter().sum::<f32>() >= total => {}
                        _ => {
                            ms3_reporters.insert(parent, (p.id, raw));
                        }
                    }
                }
            }
        }
        let span = self.next_span_start..end;
        self.next_span_start = end;
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
        } else if self.settings.needs_ms1() {
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
        self.runs.push(RunCapture {
            stem,
            title_prefix: title_prefix.map(str::to_string),
            span,
            ms1,
            ms3_reporters,
        });
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
        scan: vec![id.scan],
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
    let mut records: Vec<FeatureRecord> = Vec::new();
    let mut next_feature_id: i64 = 1;

    // ── Isobaric ──────────────────────────────────────────────────────────────
    if let Some(plex) = s.plex {
        let mut rows: Vec<TmtRow> = Vec::new();
        let (mut no_scan, mut low_purity) = (0usize, 0usize);
        for hit in &hits {
            let spec = &spectra[hit.spec_idx];
            let run = &collector.runs[hit.run_idx];
            let (quant_scan_id, raw) = if s.tmt_level >= 3 {
                match run.ms3_reporters.get(native_id(run, &spec.title)) {
                    Some((id, raw)) => (id.clone(), raw.clone()),
                    None => {
                        no_scan += 1;
                        continue;
                    }
                }
            } else {
                match collector.ms2_reporters.get(hit.spec_idx) {
                    Some(raw) => (spec.title.clone(), raw.clone()),
                    None => {
                        no_scan += 1;
                        continue;
                    }
                }
            };
            let purity = run.ms1.as_ref().and_then(|ms1| {
                let rt = spec.rt_seconds?;
                let scan = ms1.scan_at_or_before(rt)?;
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
        let tsv_rows: Vec<TmtRow> = rows
            .iter()
            .filter(|r| r.purity.is_none_or(|p| p >= s.min_purity))
            .cloned()
            .collect();
        let tmt_path = dir.join(format!("{stem}.tmt.tsv"));
        output::write_tmt_tsv(&tmt_path, plex, &tsv_rows, s.correction.is_some())
            .map_err(|e| format!("write {}: {e}", tmt_path.display()))?;
        eprintln!(
            "quant: wrote {} ({} PSMs with reporter ions{}{})",
            tmt_path.display(),
            tsv_rows.len(),
            if no_scan > 0 {
                format!(", {no_scan} without a quantification scan")
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
        for r in &rows {
            let proteins = proteins_of(
                &hits
                    .iter()
                    .find(|h| h.spec_id == r.id.spec_id)
                    .expect("row from hit")
                    .psm,
                candidates,
                index,
            );
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
        let mut all_rows: Vec<LfqRow> = Vec::new();
        for (run_idx, run) in collector.runs.iter().enumerate() {
            let Some(ms1) = run.ms1.as_ref() else {
                continue;
            };
            let targets: Vec<(usize, LfqTarget)> =
                best.iter()
                    .filter(|((r, _, _), _)| *r == run_idx)
                    .filter_map(|(_, &hi)| {
                        let h = &hits[hi];
                        let spec = &spectra[h.spec_idx];
                        let rt = spec.rt_seconds?;
                        let z = h.psm.charge_used.max(1);
                        let formula = Formula::peptide(h.cand.peptide.residues.iter().map(|aa| {
                            (aa.residue, aa.mod_.as_ref().map_or(0.0, |m| m.mass_delta))
                        }));
                        Some((
                            hi,
                            LfqTarget {
                                mono_mz: h.cand.peptide.mass() / z as f64 + PROTON,
                                charge: z,
                                anchor_rt: rt,
                                envelope: formula.envelope(n_iso),
                            },
                        ))
                    })
                    .collect();
            let params = &s.lfq_params;
            let quantified: Vec<(usize, Option<_>, Option<_>)> = targets
                .par_iter()
                .map(|(hi, t)| {
                    (
                        *hi,
                        quantify(ms1, t, params),
                        quantify_decoy(ms1, t, params),
                    )
                })
                .collect();
            let max_apex = quantified
                .iter()
                .filter_map(|(_, f, _)| f.as_ref().map(|f| f.apex_intensity))
                .fold(0.0f32, f32::max);
            let pairs: Vec<Pair> = quantified
                .iter()
                .map(|(_, f, d)| Pair {
                    target: f
                        .as_ref()
                        .map(|f| feature_score(f, params.rt_window_s, max_apex)),
                    decoy: d
                        .as_ref()
                        .map(|f| feature_score(f, params.rt_window_s, max_apex)),
                })
                .collect();
            let qvals = picked_qvalues(&pairs);
            let mut n_found = 0usize;
            for (k, (hi, feature, decoy)) in quantified.into_iter().enumerate() {
                let Some(feature) = feature else { continue };
                n_found += 1;
                let hit = &hits[hi];
                let spec = &spectra[hit.spec_idx];
                all_rows.push(LfqRow {
                    id: quant_id(hit, run, spec, candidates, index),
                    score: pairs[k].target.unwrap_or(0.0),
                    feature_q_value: qvals[k],
                    decoy_score: pairs[k].decoy,
                    feature: feature.clone(),
                });
                let _ = decoy;
            }
            eprintln!(
                "quant: {}: {} of {} precursors quantified, {} at feature q<={}",
                run.stem,
                n_found,
                targets.len(),
                all_rows
                    .iter()
                    .filter(|r| r.id.file_idx == run_idx
                        && r.feature_q_value.is_some_and(|q| q <= s.feature_fdr))
                    .count(),
                s.feature_fdr
            );
        }
        let runs: Vec<String> = collector.runs.iter().map(|r| r.stem.clone()).collect();
        let confident: Vec<LfqRow> = all_rows
            .iter()
            .filter(|r| r.feature_q_value.is_some_and(|q| q <= s.feature_fdr))
            .cloned()
            .collect();
        let lfq_path = dir.join(format!("{stem}.lfq.tsv"));
        output::write_lfq_tsv(&lfq_path, &runs, &confident)
            .map_err(|e| format!("write {}: {e}", lfq_path.display()))?;
        let feat_path = dir.join(format!("{stem}.lfq_features.tsv"));
        output::write_lfq_features_tsv(&feat_path, &all_rows)
            .map_err(|e| format!("write {}: {e}", feat_path.display()))?;
        eprintln!(
            "quant: wrote {} ({} precursor-charge rows at feature q<={}) and {} ({} features)",
            lfq_path.display(),
            confident.len(),
            s.feature_fdr,
            feat_path.display(),
            all_rows.len()
        );
        for r in &all_rows {
            let hit = hits
                .iter()
                .find(|h| h.spec_id == r.id.spec_id)
                .expect("row from hit");
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
