use std::collections::{HashMap, HashSet};
use std::ops::Range;
use std::path::{Path, PathBuf};

use model::mass::{ISOTOPE, PROTON};
use model::scan::{Ms1Scan, ProductScan, ScanContext};
use model::spectrum::Spectrum;
use model::stats::{median_f32_mut, Even};
use model::{IsomerClass, Tolerance};
use output::features::{FeatureRecord, MbrColumns};
use output::ms2::{write_ms2_table, Ms2Table};
use output::quantities::runs::{write_run_quant_table, RunPsmRow};
use output::{LfqRow, ModRecord, PercolatorPsm, QuantId, TmtRow};
use quant::formula::Formula;
use quant::isobaric::{extract_reporters, reporter_noise, signal_to_noise, CorrectionMatrix, Plex};
use quant::lfq::{
    competition_score, feature_score, quantify, quantify_decoy, FeatureQuant, LfqParams, LfqTarget,
};
use quant::ms1::Ms1RunIndex;
use quant::ms2::{Ms2Scan, ScanReporters};
use quant::tdc::{picked_qvalues, Pair};
use quant::unmix::{unmix_run, window_record, ScanInput, Strength, Survey, UnmixRun, WindowRecord};
use rayon::prelude::*;
use search::candidates::Candidate;
use search::psm::{PsmMatch, TopNQueue};
use search::index::SearchIndex;

use crate::cli::{Cli, TmtUnmix};

pub(crate) mod labels;
mod ms1;
pub(crate) use labels::LabelSettings;
use ms1::Ms1Cache;

const PURITY_ISOTOPE_PPM: f64 = 10.0;

const MS1_RESOLUTION_CHECK_PPM: f64 = 10.0;
const MS1_MIN_ISOTOPE_PARTNER_FRACTION: f64 = 0.5;

const MS2_TOP_PEAKS: usize = 100;

const REPORTER_REGION_MARGIN: f64 = 1.0;

#[derive(Debug, Clone)]
pub(crate) struct QuantSettings {
    pub plex: Option<Plex>,
    pub tmt_level: u8,
    pub tmt_tol: Tolerance,
    pub correction: Option<CorrectionMatrix>,
    pub min_purity: f32,
    pub min_sn: f32,
    pub lfq: bool,
    pub lfq_params: LfqParams,
    pub feature_fdr: f64,
    pub psm_fdr: f64,
    pub pep_cap: Option<f64>,
    pub quant_dir: Option<PathBuf>,
    pub fragment_tol: Option<Tolerance>,
    pub unmix: Option<Strength>,
    pub tmt_mbr: bool,
    pub scratch_dir: PathBuf,
    pub label: Option<LabelSettings>,
}

impl QuantSettings {
    pub(crate) fn check_args(cli: &Cli) -> Result<(), String> {
        let on = cli.tmt_unmix == TmtUnmix::On;
        if on && cli.tmt.is_none() {
            return Err("--tmt-unmix on needs --tmt <plex>".into());
        }
        if on && cli.tmt_level.is_some_and(|level| level >= 3) {
            return Err(
                "--tmt-unmix on models the co-isolation interference of MS2 reporter ions; it does \
                 not apply at --tmt-level 3"
                    .into(),
            );
        }
        if let Some(text) = &cli.label {
            if cli.lfq || cli.tmt.is_some() {
                return Err(
                    "--label quantifies labelled channels on its own; it is not combined with \
                     --lfq or --tmt"
                        .into(),
                );
            }
            let spec = label_spec(cli, text)?;
            if cli.label_reference == 0 || cli.label_reference > spec.len() {
                return Err(format!(
                    "--label-reference {} is not a channel of --label {text} (1..={})",
                    cli.label_reference,
                    spec.len()
                ));
            }
            let unnamed: Vec<&str> = spec
                .channels
                .iter()
                .filter(|c| c.qpx_label.is_none())
                .map(|c| c.name.as_str())
                .collect();
            if !unnamed.is_empty() {
                log::warn!(
                    target: "quant",
                    "--label channel(s) {} have no SDRF/QPX label (SILAC light, \
                     medium, heavy; DIMETHYL0/4/6/8): their intensities are left out of the QPX \
                     features; label.parquet has every channel",
                    unnamed.join(", ")
                );
            }
        }
        if cli.tmt.is_none() && !cli.lfq && cli.quant_dir.is_none() {
            return Ok(());
        }
        if let Some(level) = cli.tmt_level.filter(|l| !(2..=3).contains(l)) {
            return Err(format!("--tmt-level must be 2 or 3, got {level}"));
        }
        let Some(plex) = cli.tmt else {
            if cli.tmt_correction.is_some() {
                return Err("--tmt-correction needs --tmt <plex>".into());
            }
            return Ok(());
        };
        if let Some(tol) = cli.tmt_tol {
            plex.check_tolerance(tol).map_err(|e| {
                format!(
                    "--tmt {}: {e}. Use a narrower --tmt-tol (e.g. {}ppm).",
                    plex.name(),
                    plex.default_reporter_ppm()
                )
            })?;
        }
        if let Some(path) = &cli.tmt_correction {
            CorrectionMatrix::check_supported(plex)
                .map_err(|e| format!("--tmt-correction {}: {e}", path.display()))?;
        }
        Ok(())
    }

    pub(crate) fn from_cli(
        cli: &Cli,
        high_res: bool,
        tmt_level: u8,
    ) -> Result<Option<Self>, String> {
        if cli.tmt.is_none() && !cli.lfq && cli.quant_dir.is_none() && cli.label.is_none() {
            return Ok(None);
        }
        Self::check_args(cli)?;
        if cli.label.is_some() && !high_res {
            log::warn!(
                target: "quant",
                "--label: the MS2 scans are low-resolution; labelled-pair quantification \
                 assumes high-resolution (Orbitrap/TOF) MS1 scans and checks them per run"
            );
        }
        if cli.mbr && cli.tmt.is_some() && (tmt_level >= 3 || !high_res) {
            return Err(if tmt_level >= 3 {
                "--mbr with --tmt transfers the reporter ions of MS2 scans and needs them in \
                 the MS2; it does not apply to SPS-MS3 reporters (--tmt-level 3), whose \
                 transferred MS2 scans have no MS3"
                    .into()
            } else {
                "--mbr with --tmt needs high-resolution MS2 scans: this run's MS2 is \
                 low-resolution (ion trap), whose transferred scans are not reliable enough to \
                 quantify"
                    .into()
            });
        }
        let tmt_tol = cli.tmt_tol.unwrap_or(if high_res || tmt_level >= 3 {
            Tolerance::Ppm(cli.tmt.map_or(20.0, |p| p.default_reporter_ppm()))
        } else {
            Tolerance::Da(0.3)
        });
        if let (Some(plex), None) = (cli.tmt, cli.tmt_tol) {
            plex.check_tolerance(tmt_tol).map_err(|e| {
                format!(
                    "--tmt {}: {e}. The MS2 is low-resolution, so its reporter ions cannot \
                     separate this kit's channels: use --tmt-level 3 for SPS-MS3 data, or set \
                     --tmt-tol.",
                    plex.name()
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
                    log::warn!(
                        target: "quant",
                        "--tmt-correction {} is an identity matrix; no correction applied",
                        path.display()
                    );
                    None
                } else {
                    Some(m)
                }
            }
            _ => None,
        };
        let lfq_params = LfqParams {
            tol: cli.lfq_tol.unwrap_or(Tolerance::Ppm(10.0)),
            rt_window_s: cli.lfq_rt_window.unwrap_or(60.0),
            gaussian_width_s: cli.lfq_gaussian_width,
            min_cosine: cli.lfq_min_cosine as f32,
            ..LfqParams::default()
        };
        if cli.tmt_unmix == TmtUnmix::On && !high_res {
            log::warn!(
                target: "quant",
                "--tmt-unmix on: the MS2 scans are low-resolution; interference unmixing \
                 needs high-resolution MS2 reporter ions and is not applied"
            );
        }
        if cli.lfq && !high_res {
            log::warn!(
                target: "quant",
                "--lfq: the MS2 scans are low-resolution; label-free quantification \
                 assumes high-resolution (Orbitrap/TOF) MS1 scans and checks them per run"
            );
        }
        Ok(Some(Self {
            plex: cli.tmt,
            tmt_level,
            tmt_tol,
            correction,
            min_purity: cli.tmt_min_purity as f32,
            min_sn: cli.tmt_min_sn as f32,
            lfq: cli.lfq,
            lfq_params,
            feature_fdr: cli.lfq_feature_fdr(),
            psm_fdr: cli.psm_fdr(),
            pep_cap: cli.pep,
            quant_dir: cli.quant_dir.clone(),
            fragment_tol: None,
            unmix: match cli.tmt_unmix {
                _ if !high_res || tmt_level >= 3 || cli.tmt.is_none() => None,
                TmtUnmix::Off => None,
                TmtUnmix::On => Some(Strength::Fixed(1.0)),
                TmtUnmix::Auto => Some(Strength::Learned),
            },
            tmt_mbr: cli.mbr && cli.tmt.is_some(),
            scratch_dir: cli.temp_dir.clone().unwrap_or_else(|| PathBuf::from(".")),
            label: match &cli.label {
                Some(text) => {
                    let mut spec = label_spec(cli, text)?;
                    spec.reference = cli.label_reference - 1;
                    Some(LabelSettings {
                        spec,
                        fdr: cli.quant_channel_fdr(),
                        tol_ppm: match cli.lfq_tol {
                            Some(Tolerance::Ppm(v)) => v,
                            Some(Tolerance::Da(_)) => {
                                return Err("--label needs --lfq-tol in ppm".into())
                            }
                            None => 10.0,
                        },
                        mbr_fdr: cli.mbr.then_some(cli.coarse_transfer_fdr()),
                        arg_pro: cli.label_arg_pro,
                    })
                }
                None => None,
            },
        }))
    }

    pub(crate) fn needs_ms1(&self) -> bool {
        self.lfq || self.plex.is_some() || self.quant_dir.is_some() || self.label.is_some()
    }

    pub(crate) fn needs_ms1_series(&self) -> bool {
        self.lfq || self.quant_dir.is_some() || self.label.is_some()
    }

    pub(crate) fn needs_ms3(&self) -> bool {
        self.plex.is_some() && self.tmt_level >= 3
    }

    pub(crate) fn describe(&self) -> String {
        let mut parts = Vec::new();
        if let Some(p) = self.plex {
            parts.push(format!(
                "isobaric {} at MS{} (reporter tol {}, correction {}, quantification gate: \
                 purity >= {}, summed S/N >= {}; weight purity^2 x (1 - PEP))",
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
                self.min_purity,
                self.min_sn
            ));
            if let Some(strength) = self.unmix {
                parts.push(format!(
                    "interference unmixing of the MS2 reporter ions ({})",
                    match strength {
                        Strength::Fixed(_) => "full strength",
                        Strength::Learned => "learned strength",
                    }
                ));
            }
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
        if let Some(l) = &self.label {
            parts.push(l.describe());
        }
        if let Some(dir) = &self.quant_dir {
            parts.push(format!(
                "per-run MS1, MS2 and PSM tables in {}",
                dir.display()
            ));
        }
        parts.join("; ")
    }
}

pub(crate) fn label_spec(cli: &Cli, text: &str) -> Result<quant::label::LabelSpec, String> {
    let mut spec = quant::label::LabelSpec::parse(text).map_err(|e| format!("--label: {e}"))?;
    spec.enrichment_2h_metabolic = cli.label_2h_enrichment_silac;
    spec.enrichment_2h_chemical = cli.label_2h_enrichment_chemical;
    spec.enrichment_heavy = cli.label_heavy_enrichment;
    Ok(spec)
}

#[derive(Debug, Default)]
pub(crate) struct RunScans {
    pub ms1: Vec<Ms1Scan>,
    pub ms3: Vec<ProductScan>,
    pub contexts: Vec<ScanContext>,
}

struct RunCapture {
    stem: String,
    title_prefix: Option<String>,
    span: Range<usize>,
    ms1: Option<Ms1Cache>,
    ms3_reporters: HashMap<String, ReporterScan>,
    contexts: HashMap<String, ScanContext>,
}

struct ReporterScan {
    id: String,
    raw: Vec<f32>,
    scan_window: Option<(f64, f64)>,
    noise: Option<Vec<f32>>,
}

pub(crate) struct QuantCollector {
    pub settings: QuantSettings,
    ms2_reporters: Vec<Vec<f32>>,
    ms2_pic: Vec<f32>,
    runs: Vec<RunCapture>,
    next_span_start: usize,
    ms2: Vec<Ms2Scan>,
}

impl QuantCollector {
    pub(crate) fn new(settings: QuantSettings) -> Self {
        Self {
            settings,
            ms2_reporters: Vec::new(),
            ms2_pic: Vec::new(),
            runs: Vec::new(),
            next_span_start: 0,
            ms2: Vec::new(),
        }
    }

    pub(crate) fn on_spectrum(&mut self, spec: &Spectrum) {
        if let Some(plex) = self.settings.plex {
            self.ms2_reporters
                .push(extract_reporters(&spec.peaks, plex, self.settings.tmt_tol));
            if self.settings.unmix.is_some() {
                self.ms2_pic.push(peptide_ion_current(&spec.peaks, plex));
            }
        }
        if self.settings.quant_dir.is_some() {
            self.ms2.push(reduced_ms2(spec));
        }
    }

    pub(crate) fn finish_file(
        &mut self,
        path: &Path,
        title_prefix: Option<&str>,
        scans: RunScans,
        end: usize,
    ) -> Result<(), String> {
        let stem = run_stem(path);
        if self.settings.quant_dir.is_some() && scans.ms1.is_empty() {
            return Err(format!(
                "--quant-dir: no MS1 scans captured from {}; multi-run quantification needs \
                 every input's MS1 scans",
                path.display()
            ));
        }
        if self.settings.needs_ms1()
            && scans.ms1.is_empty()
            && crate::spectra::input_format_flags(path).2
        {
            return Err(format!(
                "no MS1 frames could be read from {}; the run's MS1 data is needed for \
                 quantification",
                path.display()
            ));
        }
        let ms1 = if self.settings.needs_ms1() && !scans.ms1.is_empty() {
            Some(Ms1RunIndex::new(scans.ms1))
        } else {
            None
        };
        let mut ms3_reporters: HashMap<String, ReporterScan> = HashMap::new();
        let mut ms3_no_parent = 0usize;
        if let Some(plex) = self.settings.plex {
            if self.settings.tmt_level >= 3 {
                for p in scans.ms3 {
                    let Some(parent) = p.parent_id else {
                        ms3_no_parent += 1;
                        continue;
                    };
                    let raw = extract_reporters(&p.peaks, plex, self.settings.tmt_tol);
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
                                    noise: p.noise.as_ref().map(|n| reporter_noise(n, plex)),
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
        let ms1_quant = if self.settings.lfq {
            Some("Label-free areas")
        } else if self.settings.label.is_some() {
            Some("Labelled-pair ratios")
        } else {
            None
        };
        if let (Some(what), Some(idx)) = (ms1_quant, &ms1) {
            if let Some(frac) = idx.isotope_partner_fraction(MS1_RESOLUTION_CHECK_PPM) {
                if frac < MS1_MIN_ISOTOPE_PARTNER_FRACTION {
                    log::warn!(
                        target: "quant",
                        "{stem}: the MS1 scans do not look high-resolution \
                         ({:.0}% of their most intense centroids have an isotope partner within \
                         {MS1_RESOLUTION_CHECK_PPM} ppm; Orbitrap/TOF survey scans reach 70-95%). \
                         {what} from low-resolution MS1 are not meaningful.",
                        100.0 * frac
                    );
                }
            }
        }
        if ms3_no_parent > 0 {
            log::warn!(
                target: "quant",
                "{stem}: {ms3_no_parent} MS3 scans without a parent reference \
                 (not linked to an MS2, not quantified)"
            );
        }
        if let Some(idx) = &ms1 {
            log::info!(
                target: "quant",
                "{}: {} MS1 scans ({} centroids) indexed{}",
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
            log::warn!(
                target: "quant",
                "{}: no MS1 scans captured ({}); {}",
                stem,
                path.display(),
                if self.settings.needs_ms1_series() {
                    "label-free quantification is skipped for this file"
                } else {
                    "precursor purity is unavailable for this file"
                }
            );
        }
        let mut ms2 = std::mem::take(&mut self.ms2);
        if self.settings.tmt_mbr {
            attach_reporters(
                &self.settings,
                &self.ms2_reporters[span.clone()],
                &mut ms2,
                ms1.as_ref(),
                &contexts,
                title_prefix,
            )?;
        }
        let ms1 = ms1
            .as_ref()
            .map(|index| Ms1Cache::store(index, &stem, &self.settings.scratch_dir))
            .transpose()
            .map_err(|e| format!("cache MS1 scans for {stem}: {e}"))?;
        if let Some(dir) = &self.settings.quant_dir {
            let contexts_cv = |id: &str| contexts.get(id).and_then(|c| c.faims_cv);
            write_run_ms2(
                dir,
                &stem,
                title_prefix,
                ms2,
                self.settings.fragment_tol,
                contexts_cv,
            )?;
        }
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

    pub(crate) fn run_ms1(&self) -> Vec<RunMs1<'_>> {
        self.runs
            .iter()
            .filter_map(|r| {
                r.ms1.as_ref().map(|cache| RunMs1 {
                    stem: &r.stem,
                    span: r.span.clone(),
                    cache,
                })
            })
            .collect()
    }
}

pub(crate) struct RunMs1<'a> {
    pub(crate) stem: &'a str,
    pub(crate) span: Range<usize>,
    cache: &'a Ms1Cache,
}

impl RunMs1<'_> {
    pub(crate) fn load(&self) -> Result<Ms1RunIndex, String> {
        self.cache.load()
    }
}

fn attach_reporters(
    settings: &QuantSettings,
    reporters: &[Vec<f32>],
    scans: &mut [Ms2Scan],
    ms1: Option<&Ms1RunIndex>,
    contexts: &HashMap<String, ScanContext>,
    title_prefix: Option<&str>,
) -> Result<(), String> {
    let plex = settings.plex.ok_or("--mbr with --tmt needs --tmt <plex>")?;
    if reporters.len() != scans.len() {
        return Err(format!(
            "internal: {} reporter vectors for {} MS2 scans of one run",
            reporters.len(),
            scans.len()
        ));
    }
    for (scan, raw) in scans.iter_mut().zip(reporters) {
        let id = title_prefix
            .and_then(|p| scan.native_id.strip_prefix(p))
            .unwrap_or(scan.native_id.as_str());
        let context = contexts.get(id);
        let purity = ms1.and_then(|ms1| {
            let (target, lo, hi) = scan.isolation?;
            let survey = survey_scan(ms1, context, scan.rt)?;
            quant::precursor_purity(
                ms1.peaks(survey),
                scan.precursor_mz,
                Some(target),
                scan.charge,
                lo,
                hi,
                PURITY_ISOTOPE_PPM,
            )
        });
        let sn_sum = context
            .and_then(|c| c.noise.as_ref())
            .map(|n| reporter_noise(n, plex))
            .map(|noise| {
                let corrected = match &settings.correction {
                    Some(m) => m.correct(raw),
                    None => raw.clone(),
                };
                signal_to_noise(&corrected, &noise).iter().sum::<f32>()
            });
        scan.reporters = Some(ScanReporters {
            raw: raw.clone(),
            purity,
            sn_sum,
        });
    }
    Ok(())
}

pub(crate) struct TmtMbrInput<'a> {
    pub decoys: &'a HashMap<String, PercolatorPsm>,
    pub source: &'a str,
    pub library_fdr: f64,
    pub accept_q: f64,
    pub wrong_log2: f64,
}

#[derive(Debug, Default)]
struct GateCounts {
    purity: usize,
    sn: usize,
}

impl GateCounts {
    fn check(&mut self, s: &QuantSettings, purity: Option<f32>, sn_sum: Option<f32>) -> bool {
        if purity.is_some_and(|p| p < s.min_purity) {
            self.purity += 1;
        } else if sn_sum.is_some_and(|v| v < s.min_sn) {
            self.sn += 1;
        } else {
            return true;
        }
        false
    }

    fn describe(&self, s: &QuantSettings) -> String {
        let mut parts = Vec::new();
        if self.purity > 0 {
            parts.push(format!("{} below purity {}", self.purity, s.min_purity));
        }
        if self.sn > 0 {
            parts.push(format!("{} below summed S/N {}", self.sn, s.min_sn));
        }
        if parts.is_empty() {
            String::new()
        } else {
            format!(
                ", {} failed the quantification gate (kept with quant_gate_passed false and \
                 quant_weight 0)",
                parts.join(", ")
            )
        }
    }
}

/// Roll-up weight of a reporter row that passes the gate: its purity squared (1 when unknown)
/// times the probability that its identification is right, 1 - PEP (1 when no PEP is known).
pub(crate) fn quant_weight(purity: Option<f32>, pep: Option<f64>) -> f32 {
    let right = 1.0 - pep.unwrap_or(0.0).clamp(0.0, 1.0);
    purity.unwrap_or(1.0).clamp(0.0, 1.0).powi(2) * right as f32
}

fn reduced_ms2(spec: &Spectrum) -> Ms2Scan {
    let scan = spec.scan.unwrap_or(0);
    let base = if spec.title.is_empty() {
        format!("scan={scan}")
    } else {
        spec.title.clone()
    };
    let mut peaks: Vec<(f32, f32)> = spec
        .peaks
        .iter()
        .filter(|p| p.0.is_finite() && p.1.is_finite() && p.1 > 0.0)
        .map(|&(mz, i)| (mz as f32, i))
        .collect();
    if peaks.len() > MS2_TOP_PEAKS {
        peaks.select_nth_unstable_by(MS2_TOP_PEAKS - 1, |a, b| {
            b.1.total_cmp(&a.1).then(a.0.total_cmp(&b.0))
        });
        peaks.truncate(MS2_TOP_PEAKS);
    }
    peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
    let isolation = match (
        spec.isolation_target_mz.or(Some(spec.precursor_mz)),
        spec.isolation_lower_offset,
        spec.isolation_upper_offset,
    ) {
        (Some(t), Some(lo), Some(hi)) => Some((t, lo, hi)),
        _ => None,
    };
    Ms2Scan {
        native_id: spec.title.clone(),
        spec_key: format!("{base}_{scan}"),
        rt: spec.rt_seconds.unwrap_or(f64::NAN),
        precursor_mz: spec.precursor_mz,
        charge: spec
            .precursor_charge
            .and_then(|z| u8::try_from(z).ok())
            .unwrap_or(0),
        isolation,
        faims_cv: None,
        peaks,
        reporters: None,
    }
}

fn write_run_ms2(
    dir: &Path,
    stem: &str,
    title_prefix: Option<&str>,
    mut scans: Vec<Ms2Scan>,
    fragment_tol: Option<Tolerance>,
    faims_cv: impl Fn(&str) -> Option<f32>,
) -> Result<(), String> {
    let fragment_tol = fragment_tol
        .ok_or_else(|| format!("--quant-dir: no fragment tolerance for the MS2 table of {stem}"))?;
    std::fs::create_dir_all(dir)
        .map_err(|e| format!("create --quant-dir {}: {e}", dir.display()))?;
    for s in &mut scans {
        if let Some(p) = title_prefix {
            if let Some(id) = s.native_id.strip_prefix(p) {
                s.native_id = id.to_string();
            }
        }
        s.faims_cv = faims_cv(&s.native_id);
    }
    scans.retain(|s| s.rt.is_finite());
    scans.sort_by(|a, b| a.spec_key.cmp(&b.spec_key));
    scans.dedup_by(|a, b| a.spec_key == b.spec_key);
    scans.sort_by(|a, b| a.rt.total_cmp(&b.rt).then(a.spec_key.cmp(&b.spec_key)));
    let path = dir.join(format!("{stem}.ms2.parquet"));
    let n = scans.len();
    write_ms2_table(
        &path,
        &Ms2Table {
            run: stem.to_string(),
            fragment_tol,
            top_peaks: MS2_TOP_PEAKS,
            scans,
        },
    )
    .map_err(|e| format!("write MS2 table: {e}"))?;
    let bytes = std::fs::metadata(&path).map(|m| m.len()).unwrap_or(0);
    log::info!(
        target: "quant",
        "{stem}: wrote {} ({n} MS2 scans, {:.1} MB)",
        path.display(),
        bytes as f64 / 1e6
    );
    Ok(())
}

struct Hit<'a> {
    spec_idx: usize,
    co_isolated: bool,
    run_idx: usize,
    psm: PsmMatch,
    cand: &'a Candidate,
    spec_id: String,
    q_value: Option<f64>,
    pep: Option<f64>,
    psm_index: i64,
}

fn select_hits<'a>(
    collector: &QuantCollector,
    spectra: &[Spectrum],
    queues: &[TopNQueue],
    candidates: &'a [Candidate],
    rescore: Option<&HashMap<String, PercolatorPsm>>,
    omitted: &HashSet<String>,
    cutoffs: Option<&output::Cutoffs>,
) -> Vec<Hit<'a>> {
    select_rank1(
        collector, spectra, queues, candidates, rescore, Some(omitted), cutoffs, false,
    )
}

fn select_rank1_all<'a>(
    collector: &QuantCollector,
    spectra: &[Spectrum],
    queues: &[TopNQueue],
    candidates: &'a [Candidate],
) -> Vec<Hit<'a>> {
    select_rank1(collector, spectra, queues, candidates, None, None, None, true)
}

#[allow(clippy::too_many_arguments)]
fn select_rank1<'a>(
    collector: &QuantCollector,
    spectra: &[Spectrum],
    queues: &[TopNQueue],
    candidates: &'a [Candidate],
    rescore: Option<&HashMap<String, PercolatorPsm>>,
    omitted: Option<&HashSet<String>>,
    cutoffs: Option<&output::Cutoffs>,
    include_decoys: bool,
) -> Vec<Hit<'a>> {
    let s = &collector.settings;
    let mut hits = Vec::new();
    let mut psm_index: i64 = 0;
    for (spec_idx, queue) in queues.iter().enumerate() {
        if queue.is_empty() {
            continue;
        }
        let spec = &spectra[spec_idx];
        let sorted = queue.clone().into_rank_sorted_vec();
        let rows = output::emitted_rows(spec, &sorted, candidates);
        let multi_row = rows.len() > 1;
        let mut best_rank: Vec<(Option<u64>, u32)> = Vec::new();
        for (rank, psm) in &rows {
            let key = psm.precursor_mz_override.map(f64::to_bits);
            match best_rank.iter_mut().find(|(k, _)| *k == key) {
                Some((_, r)) => *r = (*r).min(*rank),
                None => best_rank.push((key, *rank)),
            }
        }
        let scan = spec.scan.unwrap_or(0);
        let spec_id_base = if spec.title.is_empty() {
            format!("scan={scan}")
        } else {
            spec.title.clone()
        };
        for (hit, (rank, psm)) in rows.into_iter().enumerate() {
            let psm = psm.into_owned();
            let cand = &candidates[psm.primary_candidate_idx() as usize];
            let spec_id = output::format_spec_id(&spec_id_base, scan, rank, hit, multi_row);
            let rescored = rescore.map(|m| m.get(&spec_id));
            if omitted.is_some_and(|o| o.contains(&spec_id))
                || (!cand.is_decoy && matches!(rescored, Some(None)))
            {
                continue;
            }
            let this_index = psm_index;
            psm_index += 1;
            let key = psm.precursor_mz_override.map(f64::to_bits);
            if !best_rank.iter().any(|&(k, r)| k == key && r == rank) {
                continue;
            }
            if cand.is_decoy && !include_decoys {
                continue;
            }
            let (q_value, pep) = match rescored {
                Some(Some(r)) => {
                    let cutoff = cutoffs.map_or(s.psm_fdr, |c| c.or(&spec_id, s.psm_fdr));
                    if r.q_value > cutoff || s.pep_cap.is_some_and(|c| r.pep > c) {
                        continue;
                    }
                    (Some(r.q_value), Some(r.pep))
                }
                Some(None) => continue,
                None => (None, None),
            };
            let Some(run_idx) = collector.run_of(spec_idx) else {
                continue;
            };
            hits.push(Hit {
                spec_idx,
                co_isolated: key.is_some(),
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
        .flat_map(|(i, aa)| ModRecord::of_residue(i, cand.peptide.residues.len(), aa))
        .collect()
}

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

fn sole_protein(hit: &Hit<'_>, candidates: &[Candidate]) -> Option<u32> {
    let mut proteins = hit
        .psm
        .candidate_idxs
        .iter()
        .map(|&ci| candidates[ci as usize].protein_index);
    let first = proteins.next()?;
    proteins.all(|p| p == first).then_some(first as u32)
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
        mbr: None,
        group: None,
    }
}

fn is_confident(row: &LfqRow, s: &QuantSettings) -> bool {
    row.feature_q_value.is_some_and(|q| q <= s.feature_fdr)
        && row.feature.cosine >= s.lfq_params.min_cosine
}

pub(crate) fn run_stem(path: &Path) -> String {
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

pub(crate) fn check_unique_stems(paths: &[PathBuf]) -> Result<(), String> {
    let mut seen: HashMap<String, &Path> = HashMap::new();
    for p in paths {
        if let Some(first) = seen.insert(run_stem(p), p) {
            return Err(format!(
                "--quant-dir: inputs {} and {} share the run name {:?}; the per-run artifacts \
                 are named by it, so rename one of the files",
                first.display(),
                p.display(),
                run_stem(p)
            ));
        }
    }
    Ok(())
}

fn native_id<'a>(run: &RunCapture, title: &'a str) -> &'a str {
    match &run.title_prefix {
        Some(p) => title.strip_prefix(p.as_str()).unwrap_or(title),
        None => title,
    }
}

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

fn covers_reporters(scan_window: Option<(f64, f64)>, plex: Plex, tol: Tolerance) -> bool {
    let Some((lo, hi)) = scan_window else {
        return true;
    };
    let (first, last) = plex
        .channels()
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), c| {
            (a.min(c.mz), b.max(c.mz))
        });
    lo <= first - tol.as_da(first) && hi >= last + tol.as_da(last)
}

fn next_survey(ms1: &Ms1RunIndex, after: usize, faims_cv: Option<f32>) -> Option<usize> {
    const MAX_FORWARD: usize = 64;
    (after + 1..ms1.len().min(after + 1 + MAX_FORWARD))
        .find(|&i| faims_cv.is_none() || ms1.scan(i).faims_cv == faims_cv)
}

fn isolation_window_record(
    ms1: &Ms1RunIndex,
    context: Option<&ScanContext>,
    spec: &Spectrum,
    mono_mz: f64,
    charge: u8,
) -> Option<WindowRecord> {
    let rt = spec.rt_seconds?;
    let centre = spec.isolation_target_mz.unwrap_or(spec.precursor_mz);
    let lo = centre - spec.isolation_lower_offset?;
    let hi = centre + spec.isolation_upper_offset?;
    let survey = |i: usize| Survey {
        rt: ms1.rt(i),
        peaks: ms1.peaks(i),
        detection_limit: ms1.detection_limit(i),
    };
    let prev = survey_scan(ms1, context, rt);
    let cv = context.and_then(|c| c.faims_cv);
    let next = match prev {
        Some(p) => next_survey(ms1, p, cv),
        None => (0..ms1.len().min(64)).find(|&i| cv.is_none() || ms1.scan(i).faims_cv == cv),
    };
    window_record(
        prev.map(survey),
        next.map(survey),
        rt,
        lo,
        hi,
        mono_mz,
        charge,
    )
}

fn peptide_ion_current(peaks: &[(f64, f32)], plex: Plex) -> f32 {
    let top = plex
        .channels()
        .iter()
        .fold(f64::NEG_INFINITY, |a, c| a.max(c.mz))
        + REPORTER_REGION_MARGIN;
    peaks
        .iter()
        .filter(|p| p.0 > top && p.1.is_finite())
        .map(|p| p.1.max(0.0))
        .sum()
}

fn median_f32(v: &[f32]) -> f32 {
    median_f32_mut(&mut v.to_vec(), Even::Upper).unwrap_or(f32::NAN)
}

type Quantified = (usize, Option<FeatureQuant>, Option<FeatureQuant>);

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

fn lfq_target(h: &Hit<'_>, anchor_rt: f64, n_iso: usize) -> LfqTarget {
    let z = h.psm.charge_used.max(1);
    let formula = Formula::peptide(
        h.cand
            .peptide
            .residues
            .iter()
            .map(|aa| (aa.residue, aa.mod_.as_ref().map_or(0.0, |m| m.mass_delta))),
    );
    let (envelope, isotope_share) = formula.envelope_and_share(n_iso);
    LfqTarget {
        mono_mz: h.cand.peptide.mass() / z as f64 + PROTON,
        charge: z,
        anchor_rt,
        envelope,
        isotope_share,
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

fn envelope_of(cand: &Candidate, n_iso: usize) -> (Vec<f32>, f32) {
    let formula = Formula::peptide(
        cand.peptide
            .residues
            .iter()
            .map(|aa| (aa.residue, aa.mod_.as_ref().map_or(0.0, |m| m.mass_delta))),
    );
    let (env, share) = formula.envelope_and_share(n_iso);
    let total: f64 = env.iter().sum();
    let env = if total > 0.0 {
        env.iter().map(|v| (v / total) as f32).collect()
    } else {
        env.iter().map(|v| *v as f32).collect()
    };
    (env, share as f32)
}

const TABLE_ENVELOPE_ISOTOPES: usize = 4;

fn write_run_artifacts(
    collector: &QuantCollector,
    spectra: &[Spectrum],
    queues: &[TopNQueue],
    candidates: &[Candidate],
    index: &SearchIndex,
    quant_dir: &Path,
) -> Result<(), String> {
    std::fs::create_dir_all(quant_dir)
        .map_err(|e| format!("create --quant-dir {}: {e}", quant_dir.display()))?;
    let hits = select_rank1_all(collector, spectra, queues, candidates);
    for (run_idx, run) in collector.runs.iter().enumerate() {
        let ms1 = run
            .ms1
            .as_ref()
            .ok_or_else(|| format!("--quant-dir: run {} has no MS1 scans", run.stem))?;
        let ms1_path = quant_dir.join(format!("{}.ms1.parquet", run.stem));
        ms1.persist(&ms1_path)
            .map_err(|e| format!("write MS1 table: {e}"))?;
        let mut rows: Vec<RunPsmRow> = Vec::new();
        let mut mods: std::collections::BTreeMap<String, f64> = std::collections::BTreeMap::new();
        let mut no_rt = 0usize;
        for hit in hits.iter().filter(|h| h.run_idx == run_idx) {
            let spec = &spectra[hit.spec_idx];
            let Some(rt_s) = spec.rt_seconds else {
                no_rt += 1;
                continue;
            };
            for m in hit.cand.peptide.residues.iter().filter_map(|aa| aa.mod_) {
                if let Some(acc) = m.accession.as_deref().filter(|a| !a.is_empty()) {
                    mods.entry(acc.to_string()).or_insert(m.mass_delta);
                }
            }
            let z = hit.psm.charge_used.max(1);
            let peptidoform = peptidoform_of(hit.cand);
            let context = run.contexts.get(native_id(run, &spec.title));
            let (envelope, isotope_share) = envelope_of(hit.cand, TABLE_ENVELOPE_ISOTOPES);
            rows.push(RunPsmRow {
                spec_id: hit.spec_id.clone(),
                rt_s,
                precursor_key: format!("{peptidoform}/{z}"),
                peptidoform,
                charge: z,
                mz: hit.cand.peptide.mass() / z as f64 + PROTON,
                is_decoy: hit.cand.is_decoy,
                search_score: hit.psm.features.strong_score,
                faims_cv: context.and_then(|c| c.faims_cv),
                envelope,
                proteins: proteins_of(&hit.psm, candidates, index)
                    .into_iter()
                    .map(|(a, ..)| a)
                    .collect::<Vec<_>>()
                    .join(";"),
                isotope_share,
                ion_mobility: context.and_then(|c| c.ion_mobility),
            });
        }
        let table_path = quant_dir.join(format!("{}.quant.parquet", run.stem));
        let mods: Vec<(String, f64)> = mods.into_iter().collect();
        write_run_quant_table(&table_path, &run.stem, &rows, &mods)
            .map_err(|e| format!("write {}: {e}", table_path.display()))?;
        log::info!(
            target: "quant",
            "{}: wrote {} and {} ({} PSMs, one per precursor of a spectrum, targets and decoys{})",
            run.stem,
            ms1_path.display(),
            table_path.display(),
            rows.len(),
            if no_rt > 0 {
                format!("; {no_rt} without a retention time left out")
            } else {
                String::new()
            }
        );
    }
    Ok(())
}

fn write_feature_outputs(
    pdir: &Path,
    groups: Option<&crate::quant::ProteinGroups<'_>>,
    mut records: Vec<FeatureRecord>,
) -> Result<(), String> {
    std::fs::create_dir_all(pdir).map_err(|e| format!("create {}: {e}", pdir.display()))?;
    let mut tally = crate::quant::GroupTally::default();
    if let Some(pg) = groups {
        records.retain_mut(|r| match tally.assign(pg, &r.sequence) {
            Some((proteins, group)) => {
                r.proteins = proteins;
                r.group = Some(group);
                true
            }
            None => false,
        });
    }
    for r in &records {
        tally.quantify(r, &r.run_file_name)?;
    }
    let path = pdir.join("quantms.feature.parquet");
    output::write_feature_parquet(&path, &records)
        .map_err(|e| format!("write {}: {e}", path.display()))?;
    log::info!(
        target: "quant",
        "wrote {} ({} feature rows)",
        path.display(),
        records.len()
    );
    if let Some(pg) = groups {
        let mut runs: Vec<String> = records.iter().map(|r| r.run_file_name.clone()).collect();
        runs.sort();
        runs.dedup();
        tally.write_pg(pdir, pg, &runs)?;
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn quantify_label(
    collector: &QuantCollector,
    ls: &LabelSettings,
    spectra: &[Spectrum],
    queues: &[TopNQueue],
    candidates: &[Candidate],
    index: &SearchIndex,
    hits: &[Hit<'_>],
    label_decoys: Option<&HashMap<String, PercolatorPsm>>,
    cutoffs: Option<&output::Cutoffs>,
    parquet_dir: Option<&Path>,
    groups: Option<&crate::quant::ProteinGroups<'_>>,
    next_feature_id: &mut i64,
    dir: &Path,
    stem: &str,
) -> Result<Vec<FeatureRecord>, String> {
    let decoy_hits: Vec<Hit<'_>> = match label_decoys {
        Some(d) => select_rank1(
            collector,
            spectra,
            queues,
            candidates,
            Some(d),
            None,
            cutoffs,
            true,
        )
        .into_iter()
        .filter(|h| h.cand.is_decoy)
        .collect(),
        None => Vec::new(),
    };
    let (rows, label_records, label_meta) = labels::quantify(
        collector,
        ls,
        spectra,
        hits,
        &decoy_hits,
        candidates,
        index,
        parquet_dir.is_some(),
        next_feature_id,
    )?;
    let runs: Vec<String> = collector.runs.iter().map(|r| r.stem.clone()).collect();
    let channels: Vec<String> = ls.spec.channels.iter().map(|c| c.name.clone()).collect();
    let path = dir.join(format!("{stem}.label.parquet"));
    output::write_label_parquet(
        &path,
        &rows,
        &output::LabelTableMeta {
            channels: &channels,
            reference: ls.spec.reference,
            runs: &runs,
            extra: &label_meta,
        },
    )?;
    log::info!(
        target: "quant",
        "wrote {} ({} channel rows)",
        path.display(),
        rows.len()
    );
    if parquet_dir.is_some() && ls.spec.channels.iter().any(|c| c.qpx_label.is_none()) {
        log::warn!(
            target: "quant",
            "--label: channels outside the SDRF/QPX label vocabulary (SILAC light / \
             medium / heavy, DIMETHYL0/4/6/8) are not written to quantms.feature.parquet"
        );
    }
    let ratios = labels::protein_ratios(
        &rows,
        |seq| {
            let pg = groups?;
            match pg.inf.group_at(seq, pg.max_q) {
                output::protein::inference::GroupAt::Group { group, .. } => {
                    Some((group.accessions.join(";"), group.lead.clone()))
                }
                _ => None,
            }
        },
        labels::MIN_RATIO_COUNT,
    );
    let ppath = dir.join(format!("{stem}.label_proteins.tsv"));
    std::fs::write(&ppath, labels::protein_ratio_tsv(&ratios))
        .map_err(|e| format!("write {}: {e}", ppath.display()))?;
    log::info!(
        target: "quant",
        "wrote {} ({} protein ratios from at least {} peptides)",
        ppath.display(),
        ratios.len(),
        labels::MIN_RATIO_COUNT
    );
    Ok(label_records)
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn run_quant(
    collector: &QuantCollector,
    spectra: &[Spectrum],
    queues: &[TopNQueue],
    candidates: &[Candidate],
    index: &SearchIndex,
    rescore: Option<&HashMap<String, PercolatorPsm>>,
    omitted: &HashSet<String>,
    cutoffs: Option<&output::Cutoffs>,
    report_base: &Path,
    parquet_dir: Option<&Path>,
    groups: Option<&crate::quant::ProteinGroups<'_>>,
    mbr: Option<&TmtMbrInput<'_>>,
    label_decoys: Option<&HashMap<String, PercolatorPsm>>,
) -> Result<output::QuantSummary, String> {
    let s = &collector.settings;
    let mut summary = output::QuantSummary::default();
    let t0 = std::time::Instant::now();
    if let Some(dir) = &s.quant_dir {
        write_run_artifacts(collector, spectra, queues, candidates, index, dir)?;
    }
    let hits = select_hits(collector, spectra, queues, candidates, rescore, omitted, cutoffs);
    log::info!(
        target: "quant",
        "{} rank-1 target PSMs selected{} ({})",
        hits.len(),
        match rescore {
            Some(_) => format!(" at q<={}", s.psm_fdr),
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
    let hits_by_psm: HashMap<i64, &Hit<'_>> = if parquet_dir.is_some() {
        hits.iter().map(|h| (h.psm_index, h)).collect()
    } else {
        HashMap::new()
    };
    let mut records: Vec<FeatureRecord> = Vec::new();
    let mut next_feature_id: i64 = 1;

    if let Some(plex) = s.plex {
        let mut rows: Vec<TmtRow> = Vec::new();
        let (mut no_scan, mut out_of_range) = (0usize, 0usize);
        let mut gated = GateCounts::default();
        let mut with_noise = 0usize;
        let mut unmix_strengths: Vec<(String, f64)> = Vec::new();
        for (run_idx, run) in collector.runs.iter().enumerate() {
            let ms1 = run
                .ms1
                .as_ref()
                .map(Ms1Cache::load)
                .transpose()
                .map_err(|e| format!("load MS1 scans for {}: {e}", run.stem))?;
            let run_start = rows.len();
            let mut unmix_inputs: Vec<ScanInput> = Vec::new();
            for hit in hits
                .iter()
                .filter(|h| h.run_idx == run_idx && !h.co_isolated)
            {
                let spec = &spectra[hit.spec_idx];
                let context = run.contexts.get(native_id(run, &spec.title));
                let (quant_scan_id, raw, scan_window, noise) = if s.tmt_level >= 3 {
                    match run.ms3_reporters.get(native_id(run, &spec.title)) {
                        Some(r) => (r.id.clone(), r.raw.clone(), r.scan_window, r.noise.clone()),
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
                            context
                                .and_then(|c| c.noise.as_ref())
                                .map(|n| reporter_noise(n, plex)),
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
                        spec.isolation_target_mz,
                        hit.psm.charge_used,
                        lo,
                        hi,
                        PURITY_ISOTOPE_PPM,
                    )
                });
                let corrected = match &s.correction {
                    Some(m) => m.correct(&raw),
                    None => raw.clone(),
                };
                let sn = noise.map(|n| signal_to_noise(&corrected, &n));
                with_noise += usize::from(sn.is_some());
                let sn_sum = sn.as_ref().map(|v| v.iter().sum::<f32>());
                let passes = gated.check(s, purity, sn_sum);
                if s.unmix.is_some() {
                    let z = hit.psm.charge_used.max(1);
                    let window = ms1.as_ref().and_then(|ms1| {
                        isolation_window_record(
                            ms1,
                            context,
                            spec,
                            hit.cand.peptide.mass() / z as f64 + PROTON,
                            z,
                        )
                    });
                    unmix_inputs.push(ScanInput {
                        charge: z,
                        reporters: corrected.clone(),
                        pic: collector
                            .ms2_pic
                            .get(hit.spec_idx)
                            .map_or(0.0, |&v| v as f64),
                        window,
                        protein: sole_protein(hit, candidates),
                    });
                }
                rows.push(TmtRow {
                    id: quant_id(hit, run, spec, candidates, index),
                    quant_scan_id,
                    quant_ms_level: s.tmt_level,
                    purity,
                    raw,
                    corrected,
                    sn,
                    quant_weight: if passes {
                        quant_weight(purity, hit.pep)
                    } else {
                        0.0
                    },
                    quant_gate_passed: passes,
                    interference: None,
                    unmixed: None,
                    transfer: None,
                });
            }
            if let Some(strength) = s.unmix {
                let UnmixRun {
                    model,
                    strength: applied,
                    scans: unmixed,
                } = unmix_run(&unmix_inputs, strength);
                unmix_strengths.push((run.stem.clone(), applied));
                match &model {
                    Some(m) => {
                        let fractions: Vec<f32> =
                            unmixed.iter().flatten().map(|u| u.interference).collect();
                        log::info!(
                            target: "quant",
                            "{}: interference model on {} scans (precursor yield {:.3}/{:.3}/{:.3} \
                             at charge <=2/3/>=4, co-isolated {:.3}, hidden {:.3}, R2(log) {:.2}); \
                             {} scans unmixed, median interference {:.2}, strength {:.2}",
                            run.stem,
                            m.n,
                            m.precursor[0],
                            m.precursor[1],
                            m.precursor[2],
                            m.co_isolated,
                            m.hidden,
                            m.r2_log,
                            fractions.len(),
                            median_f32(&fractions),
                            applied
                        );
                    }
                    None => log::warn!(
                        target: "quant",
                        "{}: too few scans with an isolation-window record and \
                         reporter ions to fit the interference model ({} scans, {} needed); \
                         its unmix_ columns are empty",
                        run.stem,
                        unmix_inputs.len(),
                        quant::unmix::MIN_FIT_SCANS
                    ),
                }
                for (row, u) in rows[run_start..].iter_mut().zip(unmixed) {
                    if let Some(u) = u {
                        row.interference = Some(u.interference);
                        row.unmixed = Some(u.reporters);
                    }
                }
            }
        }
        if s.tmt_mbr {
            let input = mbr.ok_or("--mbr with --tmt: the rescoring result is missing")?;
            let targets = rescore.ok_or("--mbr with --tmt needs --rescore")?;
            let quant_dir = s
                .quant_dir
                .as_deref()
                .ok_or("--mbr with --tmt needs the per-run tables of --quant-dir")?;
            let params = crate::quant::tmt::mbr::Params {
                plex,
                correction: s.correction.as_ref(),
                library_fdr: input.library_fdr,
                donor_q: s.psm_fdr,
                donor_pep: s.pep_cap,
                accept_q: input.accept_q,
                wrong_log2: input.wrong_log2,
            };
            let stems: Vec<String> = collector.runs.iter().map(|r| r.stem.clone()).collect();
            let make = || crate::quant::rescore_scores(targets, input.decoys, input.source);
            let (transfers, outcome) = crate::quant::tmt::mbr::run(quant_dir, &stems, &make, &params)?;
            log::info!(
                target: "quant",
                "--mbr transferred {} reporter rows, {} at quant_q<={}{}",
                outcome.candidates,
                outcome.accepted,
                input.accept_q,
                if outcome.model_fitted {
                    ""
                } else {
                    " (no risk model could be fitted: none accepted)"
                }
            );
            for t in transfers {
                let (sequence, modifications) = crate::quant::parse_peptidoform(&t.peptidoform);
                let file_idx = collector
                    .runs
                    .iter()
                    .position(|r| r.stem == t.run)
                    .unwrap_or(0);
                rows.push(TmtRow {
                    id: QuantId {
                        file_idx,
                        run: t.run.clone(),
                        spec_id: t.spec_key.clone(),
                        scan: t.scan,
                        rt_seconds: Some(t.rt),
                        sequence,
                        peptidoform: t.peptidoform.clone(),
                        modifications,
                        charge: t.charge,
                        proteins: t.proteins.clone(),
                        is_decoy: false,
                        calculated_mz: t.calc_mz,
                        observed_mz: t.precursor_mz,
                        q_value: None,
                        pep: None,
                        missed_cleavages: None,
                        psm_index: -1,
                    },
                    quant_scan_id: t.native_id.clone(),
                    quant_ms_level: s.tmt_level,
                    purity: t.purity,
                    raw: t.raw.clone(),
                    corrected: t.corrected.clone(),
                    sn: None,
                    quant_weight: if t.accepted {
                        quant_weight(t.purity, Some(f64::from(t.info.p_wrong)))
                    } else {
                        0.0
                    },
                    quant_gate_passed: t.accepted,
                    interference: None,
                    unmixed: None,
                    transfer: Some(t.info.clone()),
                });
            }
        }
        summary.columns = plex.channels().iter().map(|c| c.name.to_string()).collect();
        for r in rows.iter().filter(|r| r.quant_gate_passed) {
            for (col, &v) in r.corrected.iter().enumerate() {
                summary
                    .values
                    .push((r.id.sequence.clone(), r.id.proteins.clone(), col, v as f64));
            }
        }
        let tmt_path = dir.join(format!("{stem}.tmt.parquet"));
        let unmix = s.unmix.map(|strength| output::UnmixMeta {
            mode: match strength {
                Strength::Fixed(_) => "on",
                Strength::Learned => "auto",
            },
            strengths: unmix_strengths.clone(),
        });
        if s.tmt_mbr {
            output::write_tmt_parquet_mbr(
                &tmt_path,
                plex,
                &rows,
                s.correction.is_some(),
                unmix.as_ref(),
            )?;
        } else {
            output::write_tmt_parquet(
                &tmt_path,
                plex,
                &rows,
                s.correction.is_some(),
                unmix.as_ref(),
            )?;
        }
        log::info!(
            target: "quant",
            "wrote {} ({} PSMs with reporter ions passing the quantification gate{}{}{})",
            tmt_path.display(),
            rows.iter().filter(|r| r.quant_gate_passed).count(),
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
            gated.describe(s)
        );
        if with_noise == 0 && s.min_sn > 0.0 && !rows.is_empty() {
            log::warn!(
                target: "quant",
                "no reporter scan carries noise values (mzML converted without \
                 noise arrays, or ion-trap reporter scans); --tmt-min-sn was not applied and \
                 the sn_ columns are empty"
            );
        }
        for r in rows.iter().filter(|_| parquet_dir.is_some()) {
            let proteins = match &r.transfer {
                Some(_) if !r.quant_gate_passed => continue,
                Some(_) => {
                    r.id.proteins
                        .iter()
                        .map(|a| (a.clone(), None, None, None, None))
                        .collect()
                }
                None => proteins_of(&hits_by_psm[&r.id.psm_index].psm, candidates, index),
            };
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
            if s.correction.is_some() || r.sn.is_some() || r.unmixed.is_some() {
                rec.additional_intensities = labels
                    .iter()
                    .enumerate()
                    .map(|(i, l)| {
                        let mut extra = Vec::new();
                        if s.correction.is_some() {
                            extra.push(("raw".to_string(), r.raw[i]));
                        }
                        if let Some(sn) = &r.sn {
                            extra.push(("sn".to_string(), sn[i]));
                        }
                        if let Some(u) = &r.unmixed {
                            extra.push(("unmixed".to_string(), u[i]));
                        }
                        (l.clone(), extra)
                    })
                    .collect();
            }
            if let Some(p) = r.purity {
                rec.additional_scores
                    .push(("precursor_purity".into(), p as f64, true));
            }
            if let Some(v) = r.sn_sum() {
                rec.additional_scores
                    .push(("reporter_sn_sum".into(), v as f64, true));
            }
            if let Some(v) = r.interference {
                rec.additional_scores
                    .push(("interference".into(), v as f64, false));
            }
            rec.additional_scores
                .push(("quant_weight".into(), r.quant_weight as f64, true));
            rec.cv_params
                .push(("quant_scan".into(), r.quant_scan_id.clone()));
            rec.cv_params
                .push(("quant_ms_level".into(), r.quant_ms_level.to_string()));
            if s.tmt_mbr {
                match &r.transfer {
                    Some(t) => {
                        rec.psm_ids.clear();
                        rec.id_run_file_name = None;
                        rec.additional_scores
                            .push(("p_wrong".into(), f64::from(t.p_wrong), false));
                        rec.additional_scores.push((
                            "noise_floor".into(),
                            f64::from(t.noise_floor),
                            false,
                        ));
                        rec.additional_scores.push((
                            "transfer_angle".into(),
                            f64::from(t.angle),
                            true,
                        ));
                        rec.additional_scores.push((
                            "donor_batch_cos".into(),
                            f64::from(t.batch_cos),
                            true,
                        ));
                        rec.mbr = Some(MbrColumns {
                            evidence_type: "transfer_ms2".into(),
                            transfer_q: Some(t.quant_q),
                            quant_q: Some(t.quant_q),
                            quantity_se: None,
                            ..MbrColumns::default()
                        });
                    }
                    None => {
                        rec.mbr = Some(MbrColumns {
                            evidence_type: "identified".into(),
                            ..MbrColumns::default()
                        });
                    }
                }
            }
            records.push(rec);
        }
    }

    if s.lfq {
        let n_iso = s.lfq_params.n_isotopes;
        let better = |h: &Hit<'_>, other: &Hit<'_>| match (h.q_value, other.q_value) {
            (Some(a), Some(b)) if a != b => a < b,
            _ => h.psm.rank_score > other.psm.rank_score,
        };
        type Precursor = (IsomerClass, u8, bool);
        let mut best: HashMap<(usize, Precursor), usize> = HashMap::new();
        let mut form: HashMap<Precursor, usize> = HashMap::new();
        for (i, h) in hits.iter().enumerate() {
            let precursor = (
                IsomerClass::of_peptide(&h.cand.peptide),
                h.psm.charge_used,
                h.cand.is_decoy,
            );
            for slot in [
                best.entry((h.run_idx, precursor.clone())).or_insert(i),
                form.entry(precursor).or_insert(i),
            ] {
                if better(h, &hits[*slot]) {
                    *slot = i;
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
            let (mut no_series, mut no_rt) = (0usize, 0usize);
            let mut targets: Vec<(usize, usize, LfqTarget)> = Vec::new();
            let mut form_of: HashMap<usize, usize> = HashMap::new();
            for ((_, precursor), &hi) in best.iter().filter(|((r, _), _)| *r == run_idx) {
                form_of.insert(hi, form[precursor]);
                let h = &hits[hi];
                let spec = &spectra[h.spec_idx];
                let Some(group) = group_of(&groups, run, spec) else {
                    no_series += 1;
                    continue;
                };
                let Some(rt) = spec.rt_seconds else {
                    no_rt += 1;
                    continue;
                };
                targets.push((hi, group, lfq_target(h, rt, n_iso)));
            }
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
                let mut id = quant_id(hit, run, spec, candidates, index);
                let shown = hits[form_of[&hi]].cand;
                id.peptidoform = peptidoform_of(shown);
                id.modifications = modifications_of(shown);
                all_rows.push(LfqRow {
                    id,
                    score: feature_score(&feature, params.rt_window_s, max_apex),
                    feature_q_value: qvals[k],
                    decoy_score: decoy
                        .as_ref()
                        .map(|d| feature_score(d, params.rt_window_s, max_apex)),
                    feature,
                });
            }
            log::info!(
                target: "quant",
                "{}: {} of {} precursors quantified, {} at feature q<={} and envelope cosine>={}{}{}",
                run.stem,
                n_found,
                targets.len(),
                all_rows
                    .iter()
                    .filter(|r| r.id.file_idx == run_idx && is_confident(r, s))
                    .count(),
                s.feature_fdr,
                s.lfq_params.min_cosine,
                if no_series > 0 {
                    format!(
                        ", {no_series} whose MS2 has no FAIMS voltage with an MS1 series \
                         (not quantified)"
                    )
                } else {
                    String::new()
                },
                if no_rt > 0 {
                    format!(", {no_rt} whose MS2 has no retention time (not quantified)")
                } else {
                    String::new()
                }
            );
        }
        let runs: Vec<String> = collector.runs.iter().map(|r| r.stem.clone()).collect();
        let confident: Vec<LfqRow> = all_rows
            .iter()
            .filter(|r| is_confident(r, s))
            .cloned()
            .collect();
        summary.columns = runs.clone();
        for r in &confident {
            summary.values.push((
                r.id.sequence.clone(),
                r.id.proteins.clone(),
                r.id.file_idx,
                r.feature.area,
            ));
        }
        let lfq_path = dir.join(format!("{stem}.lfq.parquet"));
        let accepted: Vec<bool> = all_rows.iter().map(|r| is_confident(r, s)).collect();
        output::write_lfq_parquet(
            &lfq_path,
            &all_rows,
            &accepted,
            &output::LfqTableMeta {
                runs: &runs,
                feature_fdr: s.feature_fdr,
                min_cosine: s.lfq_params.min_cosine,
            },
        )?;
        log::info!(
            target: "quant",
            "wrote {} ({} features, {} accepted at feature q<={}, cosine>={})",
            lfq_path.display(),
            all_rows.len(),
            confident.len(),
            s.feature_fdr,
            s.lfq_params.min_cosine,
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

    if let Some(ls) = &s.label {
        records.extend(quantify_label(
            collector,
            ls,
            spectra,
            queues,
            candidates,
            index,
            &hits,
            label_decoys,
            cutoffs,
            parquet_dir,
            groups,
            &mut next_feature_id,
            &dir,
            &stem,
        )?);
    }

    let quantified = s.plex.is_some() || s.lfq || s.label.is_some() || !records.is_empty();
    if let (Some(pdir), true) = (parquet_dir, quantified) {
        write_feature_outputs(pdir, groups, records)?;
    }
    log::info!(target: "quant", "done in {:.2} s", t0.elapsed().as_secs_f64());
    Ok(summary)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn faims_run() -> Ms1RunIndex {
        Ms1RunIndex::new(
            (0..6)
                .map(|i| Ms1Scan {
                    rt: i as f64,
                    peaks: Vec::new(),
                    id: format!("scan={}", i + 1),
                    faims_cv: Some(if i % 2 == 0 { -40.0 } else { -60.0 }),
                    mobility: Vec::new(),
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
            noise: None,
            ion_mobility: None,
        }
    }

    #[test]
    fn survey_scan_follows_the_precursor_reference_then_the_voltage() {
        let ms1 = faims_run();
        let c = context(Some("scan=3"), Some(-40.0));
        assert_eq!(survey_scan(&ms1, Some(&c), 3.5), Some(2));
        let c = context(None, Some(-40.0));
        assert_eq!(survey_scan(&ms1, Some(&c), 3.5), Some(2));
        let c = context(Some("scan=1000"), Some(-60.0));
        assert_eq!(survey_scan(&ms1, Some(&c), 3.5), Some(3));
        assert_eq!(survey_scan(&ms1, None, 4.5), Some(4));
        assert_eq!(survey_scan(&ms1, None, -1.0), None);
    }

    fn candidate(seq: &[u8], is_decoy: bool) -> Candidate {
        let residues = seq
            .iter()
            .map(|&r| model::AminoAcid::standard(r).unwrap())
            .collect();
        Candidate {
            peptide: model::Peptide::new(residues, b'K', b'S'),
            protein_index: 0,
            start_offset_in_protein: 0,
            is_decoy,
            is_protein_n_term: false,
            is_protein_c_term: false,
        }
    }

    fn psm(rank_score: f32, candidate: u32, co_isolated_mz: Option<f64>) -> PsmMatch {
        PsmMatch {
            spectrum_idx: 0,
            candidate_idxs: vec![candidate].into(),
            charge_used: 2,
            mass_error_ppm: 0.0,
            score: rank_score,
            rank_score,
            edge_score: 0,
            activation_method: None,
            features: search::psm::PsmFeatures::default(),
            isotope_offset: 0,
            precursor_mz_override: co_isolated_mz,
            site_evidence: None,
            tie_key: 0,
        }
    }

    fn one_spectrum(psms: Vec<PsmMatch>) -> (QuantCollector, Vec<Spectrum>, Vec<TopNQueue>) {
        let mut c = QuantCollector::new(settings(true, Some("q")));
        c.runs.push(RunCapture {
            stem: "run".into(),
            title_prefix: None,
            span: 0..1,
            ms1: None,
            ms3_reporters: HashMap::new(),
            contexts: HashMap::new(),
        });
        let spec = Spectrum {
            title: "Scan 1".into(),
            precursor_mz: 500.0,
            precursor_intensity: None,
            precursor_charge: Some(2),
            rt_seconds: Some(60.0),
            scan: Some(1),
            peaks: vec![],
            activation_method: None,
            isolation_lower_offset: None,
            isolation_upper_offset: None,
            isolation_target_mz: None,
            mobility_window: None,
        };
        let mut queue = TopNQueue::new(10);
        for p in psms {
            queue.push(p);
        }
        (c, vec![spec], vec![queue])
    }

    fn selected(hits: &[Hit<'_>]) -> Vec<(String, bool, bool)> {
        hits.iter()
            .map(|h| (h.spec_id.clone(), h.cand.is_decoy, h.co_isolated))
            .collect()
    }

    #[test]
    fn plain_search_quantifies_rank1_rows_only() {
        let cands = vec![
            candidate(b"PEPTIDEK", false),
            candidate(b"SAMPLERK", false),
            candidate(b"KEDITPEP", true),
        ];
        let (c, spectra, queues) = one_spectrum(vec![
            psm(10.0, 0, None),
            psm(9.0, 1, None),
            psm(8.0, 2, None),
        ]);
        let hits = select_rank1_all(&c, &spectra, &queues, &cands);
        assert_eq!(
            selected(&hits),
            vec![("Scan 1_1_1_0".to_string(), false, false)]
        );
    }

    #[test]
    fn chimeric_search_quantifies_each_co_isolated_precursor() {
        let cands = vec![
            candidate(b"PEPTIDEK", false),
            candidate(b"SAMPLERK", false),
            candidate(b"LLLLLLLK", false),
            candidate(b"KEDITPEP", true),
        ];
        let (c, spectra, queues) = one_spectrum(vec![
            psm(12.0, 1, Some(600.3)),
            psm(10.0, 0, None),
            psm(9.0, 2, None),
            psm(7.0, 3, Some(450.2)),
        ]);
        let hits = select_rank1_all(&c, &spectra, &queues, &cands);
        assert_eq!(
            selected(&hits),
            vec![
                ("Scan 1_1_1_0".to_string(), false, true),
                ("Scan 1_1_2_1".to_string(), false, false),
                ("Scan 1_1_4_3".to_string(), true, true),
            ]
        );
        let hits = select_hits(&c, &spectra, &queues, &cands, None, &HashSet::new(), None);
        assert_eq!(hits.len(), 2);
        let mut rescore = HashMap::new();
        for (id, q) in [("Scan 1_1_1_0", 0.2), ("Scan 1_1_2_1", 0.001)] {
            rescore.insert(
                id.to_string(),
                PercolatorPsm {
                    pooled_q_value: None,
                    psm_id: id.to_string(),
                    q_value: q,
                    pep: 0.0,
                    score: Some(1.0),
                    peptide: String::new(),
                    proteins: String::new(),
                },
            );
        }
        let hits = select_hits(&c, &spectra, &queues, &cands, Some(&rescore), &HashSet::new(), None);
        assert_eq!(
            selected(&hits),
            vec![("Scan 1_1_2_1".to_string(), false, false)]
        );
        assert_eq!(hits[0].psm_index, 1);
    }

    #[test]
    fn a_row_psms_parquet_leaves_out_takes_no_psm_index() {
        let cands = vec![candidate(b"KEDITPEP", true), candidate(b"PEPTIDEK", false)];
        let (c, spectra, queues) = one_spectrum(vec![psm(12.0, 0, Some(600.3)), psm(10.0, 1, None)]);
        let mut rescore = HashMap::new();
        rescore.insert(
            "Scan 1_1_2_1".to_string(),
            PercolatorPsm {
                pooled_q_value: None,
                psm_id: "Scan 1_1_2_1".to_string(),
                q_value: 0.001,
                pep: 0.0,
                score: Some(1.0),
                peptide: String::new(),
                proteins: String::new(),
            },
        );
        let omitted: HashSet<String> = ["Scan 1_1_1_0".to_string()].into();
        let hits = select_hits(&c, &spectra, &queues, &cands, Some(&rescore), &omitted, None);
        assert_eq!(selected(&hits), vec![("Scan 1_1_2_1".to_string(), false, false)]);
        assert_eq!(hits[0].psm_index, 0);
        let hits = select_hits(&c, &spectra, &queues, &cands, Some(&rescore), &HashSet::new(), None);
        assert_eq!(hits[0].psm_index, 1);
    }

    #[test]
    fn omitted_rows_are_skipped_without_rescoring() {
        let cands = vec![candidate(b"PEPTIDEK", false), candidate(b"KEDITPEP", false)];
        let (c, spectra, queues) = one_spectrum(vec![psm(12.0, 0, None)]);
        let hits = select_hits(&c, &spectra, &queues, &cands, None, &HashSet::new(), None);
        assert_eq!(hits.len(), 1);
        let omitted: HashSet<String> = [hits[0].spec_id.clone()].into();
        assert!(select_hits(&c, &spectra, &queues, &cands, None, &omitted, None).is_empty());
    }

    fn settings(lfq: bool, quant_dir: Option<&str>) -> QuantSettings {
        QuantSettings {
            plex: None,
            tmt_level: 2,
            tmt_tol: Tolerance::Ppm(20.0),
            correction: None,
            min_purity: 0.0,
            min_sn: 0.0,
            lfq,
            lfq_params: LfqParams::default(),
            feature_fdr: 0.05,
            psm_fdr: 0.01,
            pep_cap: None,
            quant_dir: quant_dir.map(PathBuf::from),
            fragment_tol: Some(Tolerance::Ppm(20.0)),
            unmix: None,
            tmt_mbr: false,
            scratch_dir: PathBuf::from("."),
            label: None,
        }
    }

    #[test]
    fn quant_dir_needs_ms1_scans_from_every_input() {
        let mut c = QuantCollector::new(settings(true, Some("q")));
        let e = c
            .finish_file(Path::new("/in/run_a.mzML"), None, RunScans::default(), 0)
            .unwrap_err();
        assert!(e.contains("/in/run_a.mzML") && e.contains("MS1"), "{e}");
        let mut c = QuantCollector::new(settings(true, None));
        c.finish_file(Path::new("/in/run_a.mzML"), None, RunScans::default(), 0)
            .unwrap();
    }

    #[test]
    fn duplicate_run_stems_name_both_inputs() {
        let ok = [
            PathBuf::from("a/run1.mzML"),
            PathBuf::from("a/run2.mzML.gz"),
        ];
        check_unique_stems(&ok).unwrap();
        let dup = [
            PathBuf::from("a/run1.mzML"),
            PathBuf::from("b/run2.raw"),
            PathBuf::from("b/run1.mzML.gz"),
        ];
        let e = check_unique_stems(&dup).unwrap_err();
        assert!(
            e.contains("a/run1.mzML ") && e.contains("b/run1.mzML.gz") && e.contains("\"run1\""),
            "{e}"
        );
    }

    #[test]
    fn quantification_gate_counts_the_first_failed_criterion() {
        let mut s = settings(false, None);
        s.min_purity = 0.7;
        s.min_sn = 100.0;
        let mut g = GateCounts::default();
        assert!(g.check(&s, Some(0.9), Some(500.0)));
        assert!(g.check(&s, None, None));
        assert!(!g.check(&s, Some(0.5), Some(50.0)));
        assert!(!g.check(&s, Some(0.8), Some(50.0)));
        assert_eq!((g.purity, g.sn), (1, 1));
        let text = g.describe(&s);
        assert!(text.contains("1 below purity 0.7"), "{text}");
        assert!(text.contains("1 below summed S/N 100"), "{text}");
        let off = settings(false, None);
        assert!(GateCounts::default().check(&off, Some(0.0), Some(0.0)));
        assert!(GateCounts::default().describe(&off).is_empty());
    }

    #[test]
    fn roll_up_weight_is_purity_squared_times_the_chance_the_identification_is_right() {
        assert_eq!(quant_weight(None, None), 1.0);
        assert!((quant_weight(Some(0.8), None) - 0.64).abs() < 1e-6);
        assert_eq!(quant_weight(Some(1.2), None), 1.0);
        assert!((quant_weight(Some(0.8), Some(0.25)) - 0.48).abs() < 1e-6);
        assert_eq!(quant_weight(None, Some(1.0)), 0.0);
    }

    #[test]
    fn reporter_coverage_needs_the_whole_channel_range() {
        let tol = Tolerance::Ppm(20.0);
        assert!(covers_reporters(None, Plex::Tmt16, tol));
        assert!(covers_reporters(Some((100.0, 2000.0)), Plex::Tmt16, tol));
        assert!(covers_reporters(Some((126.0, 2000.0)), Plex::Tmt16, tol));
        assert!(!covers_reporters(Some((130.0, 2000.0)), Plex::Tmt16, tol));
        assert!(!covers_reporters(Some((148.0, 1866.0)), Plex::Tmt10, tol));
        assert!(!covers_reporters(Some((100.0, 134.0)), Plex::Tmt16, tol));
        let tol = Tolerance::Ppm(10.0);
        assert!(covers_reporters(Some((126.0, 135.1625)), Plex::Tmt35, tol));
        assert!(!covers_reporters(Some((126.0, 135.1590)), Plex::Tmt35, tol));
        assert!(!covers_reporters(Some((126.2, 2000.0)), Plex::Tmt32, tol));
    }

    fn search_args(extra: &[&str]) -> Cli {
        use clap::{CommandFactory, FromArgMatches};
        let mut argv = vec![
            "andes",
            "--spectrum",
            "x",
            "--database",
            "y",
            "--output-pin",
            "o",
        ];
        argv.extend_from_slice(extra);
        let matches = crate::TopCli::command()
            .try_get_matches_from(argv)
            .expect("parse args");
        crate::SearchArgs::from_arg_matches(&matches).unwrap()
    }

    #[test]
    fn tmt_mbr_is_refused_for_ms3_reporters_and_low_resolution_ms2() {
        let from = |extra: &[&str], high_res: bool, level: u8| {
            let mut a = vec!["--tmt", "tmt11", "--mbr", "--rescore"];
            a.extend_from_slice(extra);
            QuantSettings::from_cli(&search_args(&a), high_res, level)
        };
        let s = from(&[], true, 2).unwrap().unwrap();
        assert!(s.tmt_mbr);
        let e = from(&[], true, 3).unwrap_err();
        assert!(e.contains("SPS-MS3"), "{e}");
        let e = from(&["--tmt-level", "3"], false, 3).unwrap_err();
        assert!(e.contains("SPS-MS3"), "{e}");
        let e = from(&[], false, 2).unwrap_err();
        assert!(e.contains("low-resolution"), "{e}");
        let plain = |high_res: bool, level: u8| {
            QuantSettings::from_cli(&search_args(&["--tmt", "tmt11"]), high_res, level)
                .unwrap()
                .unwrap()
        };
        assert!(!plain(true, 2).tmt_mbr);
        assert!(!plain(true, 3).tmt_mbr);
    }

    #[test]
    fn command_line_checks_need_no_input() {
        let ok = |extra: &[&str]| QuantSettings::check_args(&search_args(extra));
        assert!(ok(&[]).is_ok());
        assert!(ok(&["--tmt", "tmt35"]).is_ok());
        assert!(ok(&["--tmt", "tmt35", "--tmt-tol", "10ppm"]).is_ok());
        assert!(ok(&["--tmt", "tmt10", "--tmt-correction", "lot.txt"]).is_ok());
        let e = ok(&["--tmt", "tmt35", "--tmt-correction", "lot.txt"]).unwrap_err();
        assert!(e.contains("not supported for tmt35plex"), "{e}");
        let e = ok(&["--tmt", "tmt32", "--tmt-tol", "20ppm"]).unwrap_err();
        assert!(e.contains("e.g. 10ppm"), "{e}");
        let e = ok(&["--lfq", "--tmt-correction", "lot.txt"]).unwrap_err();
        assert!(e.contains("needs --tmt"), "{e}");
        let e = ok(&["--tmt", "tmt10", "--tmt-level", "4"]).unwrap_err();
        assert!(e.contains("--tmt-level"), "{e}");
        assert!(ok(&["--tmt", "tmt16", "--tmt-unmix", "on"]).is_ok());
        assert!(ok(&["--tmt-unmix", "auto"]).is_ok());
        assert!(ok(&["--tmt", "tmt16", "--tmt-level", "3"]).is_ok());
        assert!(ok(&["--tmt-unmix", "off"]).is_ok());
        let e = ok(&["--tmt-unmix", "on"]).unwrap_err();
        assert!(e.contains("needs --tmt"), "{e}");
        let e = ok(&["--tmt", "tmt16", "--tmt-level", "3", "--tmt-unmix", "on"]).unwrap_err();
        assert!(e.contains("--tmt-level 3"), "{e}");
    }

    #[test]
    fn peptide_ion_current_skips_the_reporter_region() {
        let peaks = [
            (126.13, 1e6f32),
            (134.15, 1e6),
            (135.1, 5.0),
            (300.0, 10.0),
            (800.0, 20.0),
        ];
        assert_eq!(peptide_ion_current(&peaks, Plex::Tmt16), 30.0);
    }

    #[test]
    fn next_survey_stays_in_the_faims_series() {
        let idx = faims_run();
        let cv0 = idx.scan(0).faims_cv;
        assert_eq!(next_survey(&idx, 0, cv0), Some(2));
        assert_eq!(next_survey(&idx, 0, None), Some(1));
        assert_eq!(next_survey(&idx, 5, cv0), None);
    }
}
