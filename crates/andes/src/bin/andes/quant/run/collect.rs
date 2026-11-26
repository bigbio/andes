//! Collecting what a run needs while it is read: MS1 scans, reporter scans, product scans and
//! scan contexts, and the gate counts of reporter attachment.

use std::collections::HashMap;
use std::ops::Range;
use std::path::Path;

use model::scan::{Ms1Scan, ProductScan, ScanContext};
use model::spectrum::Spectrum;
use output::PercolatorPsm;
use quant::isobaric::{extract_reporters, reporter_noise, signal_to_noise};
use quant::ms1::Ms1RunIndex;
use quant::ms2::{Ms2Scan, ScanReporters};

use super::artifacts::{reduced_ms2, write_run_ms2};
use super::ms1::Ms1Cache;
use super::scans::{peptide_ion_current, run_stem, survey_scan};
use super::settings::QuantSettings;
use super::{MS1_MIN_ISOTOPE_PARTNER_FRACTION, MS1_RESOLUTION_CHECK_PPM, PURITY_ISOTOPE_PPM};

#[derive(Debug, Default)]
pub(crate) struct RunScans {
    pub ms1: Vec<Ms1Scan>,
    pub ms3: Vec<ProductScan>,
    pub contexts: Vec<ScanContext>,
    /// 1/K0 tables per DIA precursor scan ID and TIMS scan number (empty otherwise).
    pub scan_mobility: Vec<Option<Vec<f32>>>,
}

pub(super) struct RunCapture {
    pub(super) stem: String,
    pub(super) title_prefix: Option<String>,
    pub(super) span: Range<usize>,
    pub(super) ms1: Option<Ms1Cache>,
    pub(super) ms3_reporters: HashMap<String, ReporterScan>,
    pub(super) contexts: HashMap<String, ScanContext>,
}

pub(super) struct ReporterScan {
    pub(super) id: String,
    pub(super) raw: Vec<f32>,
    pub(super) scan_window: Option<(f64, f64)>,
    pub(super) noise: Option<Vec<f32>>,
}

pub(crate) struct QuantCollector {
    pub settings: QuantSettings,
    pub(super) ms2_reporters: Vec<Vec<f32>>,
    pub(super) ms2_pic: Vec<f32>,
    pub(super) runs: Vec<RunCapture>,
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
        let settings = &self.settings;
        let (ms1, ms2_written) = rayon::join(
            || {
                ms1.as_ref()
                    .map(|index| Ms1Cache::store(index, &stem, &settings.scratch_dir))
                    .transpose()
                    .map_err(|e| format!("cache MS1 scans for {stem}: {e}"))
            },
            || match &settings.quant_dir {
                Some(dir) => {
                    let contexts_cv = |id: &str| contexts.get(id).and_then(|c| c.faims_cv);
                    write_run_ms2(
                        dir,
                        &stem,
                        title_prefix,
                        ms2,
                        settings.fragment_tol,
                        contexts_cv,
                    )
                }
                None => Ok(()),
            },
        );
        let ms1 = ms1?;
        ms2_written?;
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

    pub(super) fn run_of(&self, spec_idx: usize) -> Option<usize> {
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
pub(super) struct GateCounts {
    pub(super) purity: usize,
    pub(super) sn: usize,
}

impl GateCounts {
    pub(super) fn check(
        &mut self,
        s: &QuantSettings,
        purity: Option<f32>,
        sn_sum: Option<f32>,
    ) -> bool {
        if purity.is_some_and(|p| p < s.min_purity) {
            self.purity += 1;
        } else if sn_sum.is_some_and(|v| v < s.min_sn) {
            self.sn += 1;
        } else {
            return true;
        }
        false
    }

    pub(super) fn describe(&self, s: &QuantSettings) -> String {
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
