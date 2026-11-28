//! Label-free MS1 quantification of cross-linked peptide pairs (`--xl --lfq`): the precursors
//! of the accepted cross-link spectrum matches, their measurement in the run that identified
//! them, their transfer into the other runs (`--mbr`) and the tables.

mod accept;
mod align;
mod library;
mod measure;
mod tables;

use std::path::Path;

use model::spectrum::Spectrum;
use model::Tolerance;
use quant::lfq::{competition_score, FeatureQuant, LfqParams};
use quant::tdc::{picked_qvalues, Pair};
use search::crosslinks::XlSpectrumResult;
use search::PreparedSearch;

use crate::cli::Cli;
use crate::quant::run::RunMs1;
use crate::xl::fdr::{Levels, XlRun};

/// Share of the isotope envelope the traced isotopes hold.
const ISOTOPE_COVERAGE: f64 = 0.95;
/// Fewest and most isotopes traced for one precursor.
const MIN_ISOTOPES: usize = 4;
const MAX_ISOTOPES: usize = 14;
/// Isotopes computed before the envelope is cut to its coverage.
const ENVELOPE_ISOTOPES: usize = 28;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Evidence {
    Identified,
    Transferred,
}

impl Evidence {
    pub(crate) fn as_str(self) -> &'static str {
        match self {
            Evidence::Identified => "identified",
            Evidence::Transferred => "transferred",
        }
    }
}

/// One precursor in one run: where it was looked for and what was found there.
#[derive(Debug, Clone)]
pub(crate) struct Row {
    pub(crate) precursor: usize,
    pub(crate) run: usize,
    pub(crate) evidence: Evidence,
    /// The retention time of the identification, or the predicted one of a transfer, in seconds.
    pub(crate) anchor_rt: f64,
    /// The standard deviation of the prior on the apex around the anchor, in seconds.
    pub(crate) anchor_sd_s: f64,
    /// The retention-time half-window that is extracted around the anchor, in seconds.
    pub(crate) window_s: f64,
    pub(crate) scan: i32,
    pub(crate) csm_score: f64,
    pub(crate) n_csms: u32,
    pub(crate) feature: Option<FeatureQuant>,
    pub(crate) decoy: Option<FeatureQuant>,
    pub(crate) score: Option<f32>,
    pub(crate) q: Option<f64>,
}

pub(crate) struct Settings {
    /// The extraction settings that do not depend on the run: `rt_window_s` and
    /// `max_apex_climb` are set per run from the width of its elution peaks, as is
    /// `gaussian_width_s` unless `--lfq-gaussian-width` gives it.
    pub(crate) params: LfqParams,
    /// `--lfq-rt-window`: the half-window around an identification, in seconds.
    pub(crate) rt_window_s: Option<f64>,
    pub(crate) xl_fdr: f64,
    pub(crate) transfer_fdr: f64,
    pub(crate) mbr: bool,
}

impl Settings {
    fn from_cli(cli: &Cli) -> Self {
        Self {
            params: LfqParams {
                tol: cli.lfq_tol.unwrap_or(Tolerance::Ppm(10.0)),
                gaussian_width_s: cli.lfq_gaussian_width,
                min_cosine: cli.lfq_min_cosine as f32,
                ..LfqParams::default()
            },
            rt_window_s: cli.lfq_rt_window,
            xl_fdr: cli.xl_fdr,
            transfer_fdr: cli.coarse_transfer_fdr(),
            mbr: cli.mbr,
        }
    }

    pub(crate) fn params_of(&self, row: &Row) -> LfqParams {
        LfqParams {
            rt_window_s: row.window_s,
            ..self.params.clone()
        }
    }
}

/// Whether a row's value is reported as a quantity: an identification with a feature whose
/// isotope pattern fits, or a transfer at or below the transfer q-value threshold.
pub(crate) fn accepted(row: &Row, settings: &Settings) -> bool {
    if row.feature.is_none() {
        return false;
    }
    match row.evidence {
        Evidence::Identified => row.score.is_some(),
        Evidence::Transferred => row.q.is_some_and(|q| q <= settings.transfer_fdr),
    }
}

/// The feature score of every row and, within `rows`, the q-values of a target-decoy
/// competition between each precursor's feature and the feature at its shifted decoy mass.
fn score_rows(rows: &mut [Row], settings: &Settings) {
    let max_apex = rows
        .iter()
        .flat_map(|r| [r.feature.as_ref(), r.decoy.as_ref()])
        .flatten()
        .map(|f| f.apex_intensity)
        .fold(0.0f32, f32::max);
    let pairs: Vec<Pair> = rows
        .iter()
        .map(|r| {
            let params = settings.params_of(r);
            let score = |f: &Option<FeatureQuant>| {
                f.as_ref()
                    .and_then(|f| competition_score(f, &params, max_apex))
            };
            Pair {
                target: score(&r.feature),
                decoy: score(&r.decoy),
            }
        })
        .collect();
    let q = picked_qvalues(&pairs);
    for ((row, pair), q) in rows.iter_mut().zip(&pairs).zip(q) {
        row.score = pair.target;
        row.q = q;
    }
}

/// What the cross-link search hands to quantification.
pub(crate) struct Search<'a> {
    pub(crate) spectra: &'a [Spectrum],
    pub(crate) results: &'a [XlSpectrumResult],
    pub(crate) score: &'a [Option<f64>],
    pub(crate) levels: &'a Levels,
    pub(crate) linkers: &'a [XlRun],
    pub(crate) prepared: &'a PreparedSearch<'a>,
}

pub(crate) fn run(
    report_base: &Path,
    search: &Search<'_>,
    runs: &[RunMs1<'_>],
    cli: &Cli,
) -> Result<(), String> {
    let settings = Settings::from_cli(cli);
    let cohort = match &cli.xl_cohort_fdr {
        Some(base) => Some(accept::CohortFdr::read(base)?),
        None => None,
    };
    log::info!(
        target: "xl-quant",
        "cross-linked precursors from MS1 (matches and residue pairs at q <= {}{}, tol {}, rt \
         window {}, {MIN_ISOTOPES}-{MAX_ISOTOPES} isotopes, min cosine {}, transfer FDR {})",
        settings.xl_fdr,
        if cohort.is_some() {
            " of --xl-cohort-fdr"
        } else {
            ""
        },
        match settings.params.tol {
            Tolerance::Ppm(v) => format!("{v} ppm"),
            Tolerance::Da(v) => format!("{v} Da"),
        },
        match settings.rt_window_s {
            Some(v) => format!("±{v} s"),
            None => "from the elution peak width of each run".to_string(),
        },
        settings.params.min_cosine,
        settings.transfer_fdr
    );
    let lib = library::build(search, cohort.as_ref(), runs, &settings);
    if lib.precursors.is_empty() {
        log::warn!(
            target: "xl-quant",
            "nothing to quantify: no target cross-link spectrum match with a retention time \
             has a q-value and a residue-pair q-value at or below {}. --xl-fdr sets the \
             threshold",
            settings.xl_fdr
        );
        return Ok(());
    }
    let transfer = settings.mbr && runs.len() > 1;
    let anchors: Vec<Vec<align::LinearAnchor>> = if transfer {
        runs.iter()
            .map(|run| align::linear_peptides(search, run))
            .collect()
    } else {
        Vec::new()
    };
    let measure::Identified {
        mut rows,
        linear,
        scales,
    } = measure::identified(&lib, &anchors, runs, &settings)?;
    drop(anchors);
    score_rows(&mut rows, &settings);
    let run_names: Vec<&str> = runs.iter().map(|r| r.stem).collect();
    if transfer {
        let alignment = align::build(&lib, &rows, &linear, runs.len(), &settings);
        let mut counts = alignment.counts;
        let targets = align::transfers(&lib, &rows, &alignment, runs.len());
        counts.attempted = targets.len();
        let mut transferred = measure::transferred(&lib, targets, &scales, runs, &settings)?;
        score_rows(&mut transferred, &settings);
        counts.measured = transferred.iter().filter(|r| r.feature.is_some()).count();
        rows.append(&mut transferred);
        log::info!(
            target: "xl-mbr",
            "{} of {} runs aligned on {} shared anchors ({} linear peptides, {} cross-links); \
             {} transfers attempted, {} measured, {} accepted at q <= {}",
            counts.aligned_runs,
            runs.len(),
            counts.anchors,
            counts.linear_anchors,
            counts.xl_anchors,
            counts.attempted,
            counts.measured,
            rows.iter()
                .filter(|r| r.evidence == Evidence::Transferred && accepted(r, &settings))
                .count(),
            settings.transfer_fdr
        );
    } else if settings.mbr {
        log::info!(target: "xl-mbr", "--mbr needs two or more runs; nothing to transfer");
    } else if runs.len() > 1 {
        log::info!(
            target: "xl-mbr",
            "match-between-runs is off; pass --mbr to measure each cross-link in the runs that \
             did not identify it"
        );
    }
    rows.sort_by(|a, b| {
        (a.precursor, a.run, a.evidence as u8).cmp(&(b.precursor, b.run, b.evidence as u8))
    });
    let written = tables::write(report_base, &lib, &rows, &run_names, &settings)?;
    let n = |evidence: Evidence| {
        rows.iter()
            .filter(|r| r.evidence == evidence && accepted(r, &settings))
            .count()
    };
    log::info!(
        target: "xl-quant",
        "{} precursors of {} residue pairs; {} values accepted ({} identified, {} transferred) \
         across {} runs; wrote {}",
        lib.precursors.len(),
        written.residue_pairs,
        n(Evidence::Identified) + n(Evidence::Transferred),
        n(Evidence::Identified),
        n(Evidence::Transferred),
        runs.len(),
        written.features.display()
    );
    Ok(())
}

#[cfg(test)]
#[path = "quantities/tests/quantities.rs"]
mod tests;
