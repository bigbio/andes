//! Quantification settings of a run, resolved from the flags and the detected label.

use std::path::PathBuf;

use model::Tolerance;
use quant::isobaric::{CorrectionMatrix, Plex};
use quant::lfq::LfqParams;
use quant::unmix::Strength;

use super::labels::LabelSettings;
use crate::cli::{Cli, TmtUnmix};

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
