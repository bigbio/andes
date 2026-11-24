//! Glycopeptide quantification: targets from the identified glycopeptides, their MS1
//! measurement, transfers between runs, q-values and the output tables.

mod mbr;
mod measure;
mod tables;

use std::path::{Path, PathBuf};

use glyco::glycan::first::GlycanCore;
use model::mass::{ISOTOPE, PROTON};
use model::spectrum::Spectrum;
use model::Tolerance;
use output::glycopeptides::pin::select_emitted_hits;
use quant::formula::Formula;
use quant::lfq::{competition_score, FeatureQuant, LfqParams, LfqTarget};
use quant::tdc::{picked_qvalues, Pair};
use search::candidates::Candidate;
use search::glycopeptides::GlycoSpectrumResult;
use search::index::SearchIndex;

use crate::cli::Cli;
use crate::quant::run::RunMs1;
use mbr::{build_transfers, dump_match_training, Alignment, TransferCounts};
use measure::{
    build_match_references, extracted_envelope, measure, measure_transfers, Measurement,
};
use tables::{
    composition_label, output_path, proteins_of, site_label, write_feature_table,
    write_rollup_tables,
};

const ISOTOPE_COVERAGE: f64 = 0.95;

const MIN_ISOTOPES: usize = 4;
const MAX_ISOTOPES: usize = 14;

const ENVELOPE_ISOTOPES: usize = 28;

const DEFAULT_RT_WINDOW_S: f64 = 120.0;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Evidence {
    Identified,
    Transferred,
}

impl Evidence {
    fn as_str(self) -> &'static str {
        match self {
            Evidence::Identified => "Identified",
            Evidence::Transferred => "Transferred",
        }
    }
}

#[derive(Debug, Clone)]
struct Target {
    run_idx: usize,
    run: String,
    spec_id: String,
    scan: i32,
    evidence: Evidence,
    window_s: f64,
    search_s: f64,
    donor: Option<usize>,
    peptide: String,
    glycan: String,
    topology: String,
    site: String,
    proteins: String,
    charge: u8,
    isotope_offset: i8,
    mono_mz: f64,
    theo_mz: f64,
    ppm: f64,
    anchor_rt: f64,
    envelope: Vec<f64>,
}

impl Target {
    fn lfq_target(&self) -> LfqTarget {
        LfqTarget {
            mono_mz: self.mono_mz,
            charge: self.charge,
            anchor_rt: self.anchor_rt,
            envelope: self.envelope.clone(),
            isotope_share: 1.0,
        }
    }
}

struct FeatureRow {
    target: Target,
    feature: Option<FeatureQuant>,
    decoy: Option<FeatureQuant>,
    score: Option<f32>,
    q: Option<f64>,
    match_score: Option<f32>,
    decoy_match: Option<f32>,
}

struct Settings {
    params: LfqParams,
    transfer_fdr: f64,
    core: GlycanCore,
    debug: bool,
    mbr: bool,
    match_train: Option<PathBuf>,
}

impl Settings {
    fn from_cli(cli: &Cli) -> Self {
        Self {
            params: LfqParams {
                tol: cli.lfq_tol.unwrap_or(Tolerance::Ppm(10.0)),
                rt_window_s: cli.lfq_rt_window.unwrap_or(DEFAULT_RT_WINDOW_S),
                gaussian_width_s: cli.lfq_gaussian_width,
                min_cosine: cli.lfq_min_cosine as f32,
                ..LfqParams::default()
            },
            transfer_fdr: cli.coarse_transfer_fdr(),
            core: cli.glyco_core.into(),
            debug: cli.debug_glyco,
            mbr: cli.mbr,
            match_train: cli.glyco_mbr_match_train.clone(),
        }
    }
}

pub(crate) fn run(
    glyco_pin_path: &Path,
    spectra: &[Spectrum],
    results: &[GlycoSpectrumResult],
    candidates: &[Candidate],
    index: &SearchIndex,
    runs: &[RunMs1<'_>],
    cli: &Cli,
) -> Result<(), String> {
    let settings = Settings::from_cli(cli);
    log::info!(
        target: "glyco-quant",
        "intact glycopeptides from MS1 (tol {}, rt window ±{} s, {MIN_ISOTOPES}-\
         {MAX_ISOTOPES} isotopes, min cosine {}, transfer FDR {}, monoisotopic intensity)",
        match settings.params.tol {
            Tolerance::Ppm(v) => format!("{v} ppm"),
            Tolerance::Da(v) => format!("{v} Da"),
        },
        settings.params.rt_window_s,
        settings.params.min_cosine,
        settings.transfer_fdr
    );
    let targets = build_targets(spectra, results, candidates, index, runs, &settings);
    if targets.is_empty() {
        log::warn!(
            target: "glyco-quant",
            "nothing to quantify: the PIN has no emitted hit with a resolved \
             glycan composition and a precursor retention time"
        );
        return Ok(());
    }
    let (features, quantified) = measure(&targets, runs, &settings, "identified")?;
    let mut rows = feature_rows(targets, quantified, &settings);
    assign_qvalues(&mut rows, &settings);

    let run_names: Vec<&str> = runs.iter().map(|r| r.stem).collect();
    let alignment = if (settings.mbr || settings.match_train.is_some()) && runs.len() > 1 {
        Some(Alignment::build(&rows, &run_names, &settings))
    } else {
        None
    };
    if let (Some(path), Some(alignment)) = (&settings.match_train, &alignment) {
        dump_match_training(&rows, runs, &settings, alignment, path)?;
    }
    let mut transfer_counts = TransferCounts::default();
    if settings.mbr && runs.len() > 1 {
        let alignment = alignment
            .as_ref()
            .expect("alignment is built for --glyco --mbr");
        let (transfers, counts) = build_transfers(&rows, &run_names, &settings, alignment);
        transfer_counts = counts;
        if !transfers.is_empty() {
            let refs = build_match_references(&rows, &transfers, runs, &settings)?;
            let (found, measured) = measure_transfers(&transfers, &refs, runs, &settings)?;
            transfer_counts.measured = found;
            let mut transfer_rows = feature_rows(transfers, measured, &settings);
            assign_qvalues(&mut transfer_rows, &settings);
            rows.append(&mut transfer_rows);
        }
        log::info!(
            target: "glyco-mbr",
            "{} run maps ({} anchors); {} transfers attempted, {} measured, \
             {} accepted",
            transfer_counts.maps,
            transfer_counts.anchors,
            transfer_counts.candidates,
            transfer_counts.measured,
            rows.iter()
                .filter(|r| r.target.evidence == Evidence::Transferred && accepted(r, &settings))
                .count()
        );
    } else if settings.mbr {
        log::info!(target: "glyco-mbr", "--mbr needs two or more runs; nothing to transfer");
    } else if runs.len() > 1 {
        log::info!(
            target: "glyco-mbr",
            "match-between-runs is off; pass --mbr to transfer identified \
             glycopeptides into runs that did not identify them"
        );
    }

    write_feature_table(&output_path(glyco_pin_path, "lfq"), &rows, &settings)?;
    let kept = write_rollup_tables(glyco_pin_path, &rows, &run_names, &settings)?;
    let count_evidence = |evidence: Evidence| {
        rows.iter()
            .filter(|r| r.target.evidence == evidence && accepted(r, &settings))
            .count()
    };
    log::info!(
        target: "glyco-quant",
        "quantified {} features from {} targets; {} accepted ({} identified, \
         quality-gated only; {} transferred at feature q <= {}); {} glycopeptides, {} \
         glycoform-site spectra across {} runs",
        features + transfer_counts.measured,
        rows.len(),
        count_evidence(Evidence::Identified) + count_evidence(Evidence::Transferred),
        count_evidence(Evidence::Identified),
        count_evidence(Evidence::Transferred),
        settings.transfer_fdr,
        kept.glycopeptides,
        kept.sites,
        runs.len()
    );
    Ok(())
}

fn build_targets(
    spectra: &[Spectrum],
    results: &[GlycoSpectrumResult],
    candidates: &[Candidate],
    index: &SearchIndex,
    runs: &[RunMs1<'_>],
    settings: &Settings,
) -> Vec<Target> {
    let collapse = !settings.debug;
    let enumerated_only = !settings.debug;
    let mut out = Vec::new();
    let (mut no_rt, mut no_glycan, mut decoy_rows) = (0usize, 0usize, 0usize);
    for result in results {
        let spec_idx = result.spectrum_idx;
        let Some(spec) = spectra.get(spec_idx) else {
            continue;
        };
        let Some(anchor_rt) = spec.rt_seconds else {
            no_rt += 1;
            continue;
        };
        let Some(run_idx) = runs.iter().position(|r| r.span.contains(&spec_idx)) else {
            continue;
        };
        let spec_id = if spec.title.is_empty() {
            format!("scan={}", spec.scan.unwrap_or(0))
        } else {
            spec.title.clone()
        };
        let scan = spec.scan.unwrap_or(0);
        for hit_idx in select_emitted_hits(&result.hits, collapse, enumerated_only) {
            let hit = &result.hits[hit_idx];
            let cand_idx = hit.psm.primary_candidate_idx() as usize;
            let Some(cand) = candidates.get(cand_idx) else {
                continue;
            };
            if cand.is_decoy {
                decoy_rows += 1;
                continue;
            }
            let Some(comp) = hit.glycan_key.glycan.as_ref() else {
                no_glycan += 1;
                continue;
            };
            let charge = hit.psm.charge_used.max(1);
            let z = charge as f64;
            let intact_mass = cand.peptide.mass() + comp.mass;
            let theo_mz = intact_mass / z + PROTON;
            let offset = hit.psm.isotope_offset.max(0) as usize;
            let mono_mz = theo_mz + offset as f64 * ISOTOPE / z;
            let mut formula = peptide_formula(cand);
            formula.add(Formula::glycan(
                comp.hex,
                comp.hexnac,
                comp.fuc,
                comp.neuac,
                comp.neugc,
            ));
            let envelope = extracted_envelope(&formula.envelope(ENVELOPE_ISOTOPES), offset);
            let observed_mono_mz = spec.precursor_mz - offset as f64 * ISOTOPE / z;
            let ppm = (theo_mz - observed_mono_mz) / theo_mz * 1e6;
            let residues: Vec<u8> = cand.peptide.residues.iter().map(|aa| aa.residue).collect();
            out.push(Target {
                run_idx,
                run: runs[run_idx].stem.to_string(),
                spec_id: format!("{}_glyco_{}_{}", spec_id, scan, hit_idx + 1),
                scan,
                evidence: Evidence::Identified,
                window_s: settings.params.rt_window_s,
                search_s: 0.0,
                donor: None,
                peptide: cand.peptide.to_string(),
                glycan: composition_label(comp),
                topology: comp.tree.canonical(),
                site: site_label(&residues, cand.peptide.post, settings.core),
                proteins: proteins_of(&hit.psm, candidates, index),
                charge,
                isotope_offset: hit.psm.isotope_offset,
                mono_mz,
                theo_mz,
                ppm,
                anchor_rt,
                envelope,
            });
        }
    }
    if no_rt > 0 {
        log::warn!(
            target: "glyco-quant",
            "{no_rt} spectra have no retention time and are not quantified \
             (their XIC anchor would be undefined)"
        );
    }
    if no_glycan > 0 {
        log::warn!(
            target: "glyco-quant",
            "{no_glycan} emitted hits carry no enumerated glycan composition \
             and are not quantified (their intact precursor mass is unknown)"
        );
    }
    if decoy_rows > 0 {
        log::info!(
            target: "glyco-quant",
            "{decoy_rows} emitted rows are decoy-backbone hits and are not \
             quantified (the PIN carries them for Percolator; a decoy precursor is not a \
             glycopeptide)"
        );
    }
    out
}

fn peptide_formula(cand: &Candidate) -> Formula {
    Formula::peptide(
        cand.peptide
            .residues
            .iter()
            .map(|aa| (aa.residue, aa.mod_.as_ref().map_or(0.0, |m| m.mass_delta))),
    )
}

fn params_of(row: &FeatureRow, settings: &Settings) -> LfqParams {
    LfqParams {
        rt_window_s: row.target.window_s,
        ..settings.params.clone()
    }
}

fn feature_rows(
    targets: Vec<Target>,
    quantified: Vec<Measurement>,
    settings: &Settings,
) -> Vec<FeatureRow> {
    debug_assert_eq!(targets.len(), quantified.len());
    let max_apex = quantified
        .iter()
        .flat_map(|m| [m.feature.as_ref(), m.decoy.as_ref()])
        .flatten()
        .map(|f| f.apex_intensity)
        .fold(0.0f32, f32::max);
    targets
        .into_iter()
        .zip(quantified)
        .map(|(target, m)| {
            let params = LfqParams {
                rt_window_s: target.window_s,
                ..settings.params.clone()
            };
            let score = match (m.match_score, m.decoy_match) {
                (Some(score), Some(_)) => Some(score),
                _ => m
                    .feature
                    .as_ref()
                    .and_then(|f| competition_score(f, &params, max_apex)),
            };
            FeatureRow {
                target,
                feature: m.feature,
                decoy: m.decoy,
                score,
                q: None,
                match_score: m.match_score,
                decoy_match: m.decoy_match,
            }
        })
        .collect()
}

fn assign_qvalues(rows: &mut [FeatureRow], settings: &Settings) {
    let max_apex = rows
        .iter()
        .flat_map(|r| [r.feature.as_ref(), r.decoy.as_ref()])
        .flatten()
        .map(|f| f.apex_intensity)
        .fold(0.0f32, f32::max);
    let pairs: Vec<Pair> = rows
        .iter()
        .map(|r| Pair {
            target: r.score,
            decoy: r.decoy_match.or_else(|| {
                r.decoy
                    .as_ref()
                    .and_then(|f| competition_score(f, &params_of(r, settings), max_apex))
            }),
        })
        .collect();
    for (row, q) in rows.iter_mut().zip(picked_qvalues(&pairs)) {
        row.q = q;
    }
}

fn accepted(row: &FeatureRow, settings: &Settings) -> bool {
    if row.feature.is_none() {
        return false;
    }
    match row.target.evidence {
        Evidence::Identified => row.score.is_some(),
        Evidence::Transferred => row.q.is_some_and(|q| q <= settings.transfer_fdr),
    }
}

#[cfg(test)]
#[path = "quantities/tests/quantities.rs"]
mod tests;
