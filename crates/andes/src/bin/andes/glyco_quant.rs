use std::collections::{BTreeSet, HashMap};
use std::io::Write;
use std::path::{Path, PathBuf};

use glyco::glycan_db::GlycanComp;
use glyco::glycan_first::GlycanCore;
use model::mass::{ISOTOPE, PROTON};
use model::spectrum::Spectrum;
use model::stats::{median_mad, Even, MAD_TO_SD};
use model::Tolerance;
use output::glycopeptides::pin::select_emitted_hits;
use quant::formula::Formula;
use quant::lfq::{competition_score, quantify, quantify_decoy, FeatureQuant, LfqParams, LfqTarget};
use quant::tdc::{picked_qvalues, Pair};
use rayon::prelude::*;
use search::candidates::Candidate;
use search::glycopeptides::GlycoSpectrumResult;
use search::index::SearchIndex;

use crate::cli::Cli;
use crate::glyco_match;
use crate::quant_run::RunMs1;

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

fn extracted_envelope(full: &[f64], offset: usize) -> Vec<f64> {
    let tail = full.get(offset.min(full.len())..).unwrap_or(&[]);
    if tail.is_empty() {
        return vec![1.0];
    }
    let total: f64 = tail.iter().sum();
    let mut k = tail.len();
    let mut cum = 0.0;
    for (i, &v) in tail.iter().enumerate() {
        cum += v;
        if i + 1 >= MIN_ISOTOPES && cum >= ISOTOPE_COVERAGE * total {
            k = i + 1;
            break;
        }
    }
    let k = k.clamp(MIN_ISOTOPES, MAX_ISOTOPES).min(tail.len());
    let mut env: Vec<f64> = tail[..k].to_vec();
    let sum: f64 = env.iter().sum();
    if sum > 0.0 {
        for v in &mut env {
            *v /= sum;
        }
    }
    env
}

#[derive(Debug, Clone, Default)]
struct Measurement {
    feature: Option<FeatureQuant>,
    decoy: Option<FeatureQuant>,
    match_score: Option<f32>,
    decoy_match: Option<f32>,
}

fn measure(
    targets: &[Target],
    runs: &[RunMs1<'_>],
    settings: &Settings,
    phase: &str,
) -> Result<(usize, Vec<Measurement>), String> {
    measure_with(targets, runs, settings, phase, measure_identified)
}

fn measure_identified(
    ms1: &quant::ms1_index::Ms1RunIndex,
    t: &Target,
    params: &LfqParams,
) -> Measurement {
    let lfq = t.lfq_target();
    let params = LfqParams {
        n_isotopes: t.envelope.len(),
        rt_window_s: t.window_s,
        ..params.clone()
    };
    let mut feature = quantify(ms1, &lfq, &params);
    if let Some(f) = feature.as_mut() {
        f.area = mono_area(f, &lfq.envelope);
    }
    Measurement {
        feature,
        decoy: quantify_decoy(ms1, &lfq, &params),
        match_score: None,
        decoy_match: None,
    }
}

fn measure_with(
    targets: &[Target],
    runs: &[RunMs1<'_>],
    settings: &Settings,
    phase: &str,
    f: impl Fn(&quant::ms1_index::Ms1RunIndex, &Target, &LfqParams) -> Measurement + Sync,
) -> Result<(usize, Vec<Measurement>), String> {
    let mut by_run: Vec<Vec<usize>> = vec![Vec::new(); runs.len()];
    for (i, t) in targets.iter().enumerate() {
        by_run[t.run_idx].push(i);
    }
    let mut slot: Vec<Option<Measurement>> = vec![None; targets.len()];
    let mut found = 0usize;
    for (ri, run) in runs.iter().enumerate() {
        if by_run[ri].is_empty() {
            continue;
        }
        let t0 = std::time::Instant::now();
        let ms1 = run.load().map_err(|e| format!("--glyco --lfq: {e}"))?;
        let part: Vec<(usize, Measurement)> = by_run[ri]
            .par_iter()
            .map(|&i| (i, f(&ms1, &targets[i], &settings.params)))
            .collect();
        let run_found = part.iter().filter(|(_, m)| m.feature.is_some()).count();
        found += run_found;
        log::info!(
            target: "glyco-quant",
            "{phase} {}: {} features in {} targets, {:.2}s",
            run.stem,
            run_found,
            part.len(),
            t0.elapsed().as_secs_f64()
        );
        for (i, m) in part {
            slot[i] = Some(m);
        }
    }
    Ok((
        found,
        slot.into_iter().map(Option::unwrap_or_default).collect(),
    ))
}

struct MatchReferences {
    target: glyco_match::Reference,
    decoy: Option<glyco_match::Reference>,
}

fn build_match_references(
    rows: &[FeatureRow],
    transfers: &[Target],
    runs: &[RunMs1<'_>],
    settings: &Settings,
) -> Result<HashMap<usize, MatchReferences>, String> {
    let mut by_run: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); runs.len()];
    for t in transfers {
        if let Some(d) = t.donor {
            by_run[rows[d].target.run_idx].insert(d);
        }
    }
    let mut refs: HashMap<usize, MatchReferences> = HashMap::new();
    for (ri, ids) in by_run.iter().enumerate() {
        if ids.is_empty() {
            continue;
        }
        let ms1 = runs[ri]
            .load()
            .map_err(|e| format!("--glyco --mbr match reference: {e}"))?;
        for &id in ids {
            let row = &rows[id];
            let Some(f) = &row.feature else {
                continue;
            };
            let t = &row.target;
            let build = |mono_mz: f64| {
                glyco_match::Reference::build(
                    &ms1,
                    mono_mz,
                    t.charge,
                    t.envelope.len(),
                    t.anchor_rt,
                    t.window_s,
                    f.rt_start,
                    f.rt_stop,
                    settings.params.tol,
                )
            };
            let Some(target) = build(t.mono_mz) else {
                continue;
            };
            let decoy = build(t.mono_mz + settings.params.decoy_mz_shift);
            refs.insert(id, MatchReferences { target, decoy });
        }
    }
    Ok(refs)
}

fn measure_transfers(
    targets: &[Target],
    refs: &HashMap<usize, MatchReferences>,
    runs: &[RunMs1<'_>],
    settings: &Settings,
) -> Result<(usize, Vec<Measurement>), String> {
    measure_with(targets, runs, settings, "transferred", |ms1, t, params| {
        let reference = t.donor.and_then(|d| refs.get(&d));
        let mut anchor = t.anchor_rt;
        let mut window = t.window_s;
        let mut decoy_anchor = anchor;
        let mut decoy_window = window;
        let mut match_score = None;
        let mut decoy_match = None;
        if let Some(r) = reference {
            let predicted = r.target.window_start(t.anchor_rt);
            if let Some(m) =
                glyco_match::best_match(ms1, &r.target, predicted, t.search_s, params.tol)
            {
                anchor = 0.5 * (m.start_rt + m.stop_rt);
                window = m.window_s;
                match_score = Some(m.score);
            }
            if let Some(d) = &r.decoy {
                let predicted = d.window_start(t.anchor_rt);
                if let Some(m) = glyco_match::best_match(ms1, d, predicted, t.search_s, params.tol)
                {
                    decoy_anchor = 0.5 * (m.start_rt + m.stop_rt);
                    decoy_window = m.window_s;
                    decoy_match = Some(m.score);
                }
            }
        }
        if match_score.is_none() || decoy_match.is_none() {
            match_score = None;
            decoy_match = None;
            anchor = t.anchor_rt;
            window = t.window_s;
        }
        let lfq = LfqTarget {
            anchor_rt: anchor,
            ..t.lfq_target()
        };
        let p = LfqParams {
            n_isotopes: t.envelope.len(),
            rt_window_s: window,
            ..params.clone()
        };
        let mut feature = quantify(ms1, &lfq, &p);
        if let Some(f) = feature.as_mut() {
            f.area = mono_area(f, &t.envelope);
        }
        let decoy = if match_score.is_some() {
            let twin = LfqTarget {
                mono_mz: t.mono_mz + params.decoy_mz_shift,
                anchor_rt: decoy_anchor,
                ..t.lfq_target()
            };
            quantify(
                ms1,
                &twin,
                &LfqParams {
                    n_isotopes: t.envelope.len(),
                    rt_window_s: decoy_window,
                    ..params.clone()
                },
            )
        } else {
            quantify_decoy(ms1, &lfq, &p)
        };
        Measurement {
            feature,
            decoy,
            match_score,
            decoy_match,
        }
    })
}

fn mono_area(f: &FeatureQuant, envelope: &[f64]) -> f64 {
    match (f.isotope_areas.first(), envelope.first()) {
        (Some(&mono), Some(&share)) if share > 0.0 => mono / share,
        _ => 0.0,
    }
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

const MBR_MIN_ANCHORS: usize = 50;
const MBR_MAX_RESIDUAL_S: f64 = 45.0;
const MBR_ANCHORS_PER_KNOT: usize = 300;
const MBR_MIN_KNOTS: usize = 3;
const MBR_MAX_KNOTS: usize = 10;
const MBR_MIN_WINDOW_S: f64 = 10.0;
const MBR_MAX_WINDOW_S: f64 = 60.0;
const MBR_MATCH_MIN_SEARCH_S: f64 = 60.0;
const MBR_MATCH_MAX_SEARCH_S: f64 = 120.0;

#[derive(Debug, Clone, Copy, Default, PartialEq)]
struct TransferCounts {
    maps: usize,
    anchors: usize,
    candidates: usize,
    measured: usize,
}

#[derive(Debug, Clone, PartialEq)]
struct RtMap {
    xs: Vec<f64>,
    ys: Vec<f64>,
    sigma_s: f64,
}

impl RtMap {
    fn identity() -> RtMap {
        RtMap {
            xs: vec![0.0, 1.0],
            ys: vec![0.0, 1.0],
            sigma_s: 0.0,
        }
    }

    fn forward(&self, rt: f64) -> f64 {
        interpolate(&self.xs, &self.ys, rt)
    }

    fn inverse(&self, rt_ref: f64) -> Option<f64> {
        let ys_are_increasing = self.ys.windows(2).all(|w| w[1] > w[0]);
        ys_are_increasing.then(|| interpolate(&self.ys, &self.xs, rt_ref))
    }

    fn inverse_scale(&self, rt_ref: f64) -> f64 {
        let i = self
            .ys
            .partition_point(|&k| k <= rt_ref)
            .saturating_sub(1)
            .min(self.ys.len() - 2);
        (self.xs[i + 1] - self.xs[i]) / (self.ys[i + 1] - self.ys[i])
    }
}

fn interpolate(a: &[f64], b: &[f64], v: f64) -> f64 {
    debug_assert!(a.len() >= 2 && a.len() == b.len());
    let i = a
        .partition_point(|&k| k <= v)
        .saturating_sub(1)
        .min(a.len() - 2);
    let slope = (b[i + 1] - b[i]) / (a[i + 1] - a[i]);
    b[i] + slope * (v - a[i])
}

fn fit_line(pairs: &[(f64, f64)]) -> Option<(f64, f64)> {
    if pairs.len() < 3 {
        return None;
    }
    let n = pairs.len() as f64;
    let mx = pairs.iter().map(|p| p.0).sum::<f64>() / n;
    let my = pairs.iter().map(|p| p.1).sum::<f64>() / n;
    let sxx: f64 = pairs.iter().map(|p| (p.0 - mx).powi(2)).sum();
    if sxx.is_nan() || sxx <= 0.0 {
        return None;
    }
    let sxy: f64 = pairs.iter().map(|p| (p.0 - mx) * (p.1 - my)).sum();
    let slope = sxy / sxx;
    if !slope.is_finite() || slope <= 0.0 {
        return None;
    }
    Some((my - slope * mx, slope))
}

fn fit_rt_map(pairs: &[(f64, f64)]) -> Option<RtMap> {
    let (intercept, slope) = fit_line(pairs)?;
    let residuals: Vec<f64> = pairs
        .iter()
        .map(|&(x, y)| y - (intercept + slope * x))
        .collect();
    let (center, mad) = median_mad(&residuals, Even::Mean).unwrap_or((f64::NAN, f64::NAN));
    let mad = MAD_TO_SD * mad;
    let mut kept: Vec<(f64, f64)> = pairs
        .iter()
        .zip(&residuals)
        .filter(|(_, r)| mad <= 0.0 || ((*r - center).abs()) <= 3.0 * mad)
        .map(|(p, _)| *p)
        .collect();
    if kept.len() < 3 {
        return None;
    }
    kept.sort_by(|a, b| a.0.total_cmp(&b.0));

    let bins = (kept.len() / MBR_ANCHORS_PER_KNOT).clamp(MBR_MIN_KNOTS, MBR_MAX_KNOTS);
    let (mut xs, mut ys) = (Vec::with_capacity(bins), Vec::with_capacity(bins));
    for k in 0..bins {
        let bin = &kept[k * kept.len() / bins..(k + 1) * kept.len() / bins];
        let x = crate::quant_mbr::median(bin.iter().map(|p| p.0).collect());
        let y = crate::quant_mbr::median(bin.iter().map(|p| p.1).collect());
        let advances = xs.last().is_none_or(|&px| x > px) && ys.last().is_none_or(|&py| y > py);
        if advances {
            xs.push(x);
            ys.push(y);
        }
    }
    if xs.len() < 2 {
        return None;
    }
    let map = RtMap {
        xs,
        ys,
        sigma_s: 0.0,
    };
    let dev: Vec<f64> = kept.iter().map(|&(x, y)| y - map.forward(x)).collect();
    let sigma = MAD_TO_SD * median_mad(&dev, Even::Mean).map_or(f64::NAN, |(_, mad)| mad);
    let map = RtMap {
        sigma_s: sigma,
        ..map
    };
    map.sigma_s.is_finite().then_some(map)
}

type GlycoKey = (String, String, u8);

fn glyco_key(t: &Target) -> GlycoKey {
    (t.peptide.clone(), t.glycan.clone(), t.charge)
}

struct Alignment {
    identified: Vec<std::collections::HashSet<GlycoKey>>,
    donor: Vec<HashMap<GlycoKey, usize>>,
    maps: Vec<Option<RtMap>>,
    counts: TransferCounts,
}

impl Alignment {
    fn build(rows: &[FeatureRow], run_names: &[&str], settings: &Settings) -> Self {
        let mut counts = TransferCounts {
            maps: 1,
            ..TransferCounts::default()
        };
        let mut identified: Vec<std::collections::HashSet<GlycoKey>> =
            vec![std::collections::HashSet::new(); run_names.len()];
        let mut donor: Vec<HashMap<GlycoKey, usize>> = vec![HashMap::new(); run_names.len()];
        for (i, row) in rows.iter().enumerate() {
            if row.target.evidence != Evidence::Identified {
                continue;
            }
            let key = glyco_key(&row.target);
            identified[row.target.run_idx].insert(key.clone());
            if !accepted(row, settings) {
                continue;
            }
            let score = row.score.unwrap_or(f32::NEG_INFINITY);
            match donor[row.target.run_idx].get(&key) {
                Some(&j) if rows[j].score.unwrap_or(f32::NEG_INFINITY) >= score => {}
                _ => {
                    donor[row.target.run_idx].insert(key, i);
                }
            }
        }

        let mut maps: Vec<Option<RtMap>> = vec![None; run_names.len()];
        maps[0] = Some(RtMap::identity());
        for (r, map) in maps.iter_mut().enumerate().skip(1) {
            let pairs: Vec<(f64, f64)> = donor[r]
                .iter()
                .filter_map(|(key, &j)| {
                    let other = *donor[0].get(key)?;
                    let rt = rows[j].feature.as_ref()?.apex_rt;
                    let rt_ref = rows[other].feature.as_ref()?.apex_rt;
                    Some((rt, rt_ref))
                })
                .collect();
            if pairs.len() < MBR_MIN_ANCHORS {
                log::warn!(
                    target: "glyco-mbr",
                    "{}: {} shared glycopeptides, fewer than {MBR_MIN_ANCHORS} \
                     needed for an RT map; no transfers into or out of it",
                    run_names[r],
                    pairs.len()
                );
                continue;
            }
            let Some(fitted) = fit_rt_map(&pairs) else {
                continue;
            };
            if fitted.sigma_s > MBR_MAX_RESIDUAL_S {
                log::warn!(
                    target: "glyco-mbr",
                    "{}: RT map residual {:.1}s is above {MBR_MAX_RESIDUAL_S}s; \
                     no transfers into or out of it",
                    run_names[r], fitted.sigma_s
                );
                continue;
            }
            counts.maps += 1;
            counts.anchors += pairs.len();
            *map = Some(fitted);
        }
        Self {
            identified,
            donor,
            maps,
            counts,
        }
    }

    fn predict(&self, from_idx: usize, to_idx: usize, rt: f64) -> Option<(f64, f64)> {
        let from = self.maps[from_idx].as_ref()?;
        let to = self.maps[to_idx].as_ref()?;
        let t_ref = from.forward(rt);
        let predicted = to.inverse(t_ref)?;
        let sigma = (from.sigma_s.powi(2) + (to.sigma_s * to.inverse_scale(t_ref)).powi(2)).sqrt();
        Some((predicted, sigma))
    }
}

fn build_transfers(
    rows: &[FeatureRow],
    run_names: &[&str],
    settings: &Settings,
    alignment: &Alignment,
) -> (Vec<Target>, TransferCounts) {
    let Alignment {
        identified,
        donor,
        maps,
        counts,
    } = alignment;
    let mut counts = *counts;
    let mut best: HashMap<(GlycoKey, usize), (f32, Target)> = HashMap::new();
    for (i, row) in rows.iter().enumerate() {
        if row.target.evidence != Evidence::Identified
            || !accepted(row, settings)
            || donor[row.target.run_idx].get(&glyco_key(&row.target)) != Some(&i)
        {
            continue;
        }
        let Some(feature) = &row.feature else {
            continue;
        };
        if maps[row.target.run_idx].is_none() {
            continue;
        }
        let key = glyco_key(&row.target);
        for (r, to) in maps.iter().enumerate() {
            if r == row.target.run_idx {
                continue;
            }
            if to.is_none() {
                continue;
            }
            if identified[r].contains(&key) {
                continue;
            }
            let Some((rt, sigma)) = alignment.predict(row.target.run_idx, r, feature.apex_rt)
            else {
                continue;
            };
            let window = (3.0 * sigma).clamp(MBR_MIN_WINDOW_S, MBR_MAX_WINDOW_S);
            let search = (3.0 * sigma).clamp(MBR_MATCH_MIN_SEARCH_S, MBR_MATCH_MAX_SEARCH_S);
            let score = row.score.unwrap_or(f32::NEG_INFINITY);
            let slot = (key.clone(), r);
            if best
                .get(&slot)
                .is_some_and(|(existing, _)| *existing >= score)
            {
                continue;
            }
            let mut target = row.target.clone();
            target.run_idx = r;
            target.run = run_names[r].to_string();
            target.spec_id = format!("{}_mbr_{}", target.spec_id, run_names[r]);
            target.scan = 0;
            target.evidence = Evidence::Transferred;
            target.anchor_rt = rt;
            target.window_s = window;
            target.search_s = search;
            target.donor = Some(i);
            best.insert(slot, (score, target));
        }
    }
    let mut transfers: Vec<Target> = best.into_values().map(|(_, t)| t).collect();
    transfers.sort_by(|a, b| {
        (
            a.run_idx,
            &a.peptide,
            &a.glycan,
            &a.topology,
            a.charge,
            &a.spec_id,
        )
            .cmp(&(
                b.run_idx,
                &b.peptide,
                &b.glycan,
                &b.topology,
                b.charge,
                &b.spec_id,
            ))
    });
    counts.candidates = transfers.len();
    (transfers, counts)
}

fn dump_match_training(
    rows: &[FeatureRow],
    runs: &[RunMs1<'_>],
    settings: &Settings,
    alignment: &Alignment,
    path: &Path,
) -> Result<(), String> {
    let mut out = std::io::BufWriter::new(
        std::fs::File::create(path).map_err(|e| format!("{}: {e}", path.display()))?,
    );
    writeln!(
        out,
        "# glyco-mbr match training table: label, match-model features \
         f0..f9, peptide, donor/target run"
    )
    .map_err(|e| format!("{}: {e}", path.display()))?;
    let mut pairs = 0usize;
    let mut positives = 0usize;
    let mut negatives = 0usize;
    for d in 0..runs.len() {
        if alignment.maps[d].is_none() || alignment.donor[d].is_empty() {
            continue;
        }
        let ms1_d = runs[d]
            .load()
            .map_err(|e| format!("--glyco-mbr-match-train: {e}"))?;
        let mut refs: Vec<(
            GlycoKey,
            usize,
            glyco_match::Reference,
            Option<glyco_match::Reference>,
        )> = Vec::with_capacity(alignment.donor[d].len());
        for (key, &j) in &alignment.donor[d] {
            let row = &rows[j];
            let Some(f) = &row.feature else { continue };
            let t = &row.target;
            let build = |mono_mz: f64| {
                glyco_match::Reference::build(
                    &ms1_d,
                    mono_mz,
                    t.charge,
                    t.envelope.len(),
                    t.anchor_rt,
                    t.window_s,
                    f.rt_start,
                    f.rt_stop,
                    settings.params.tol,
                )
            };
            let Some(target) = build(t.mono_mz) else {
                continue;
            };
            let decoy = build(t.mono_mz + settings.params.decoy_mz_shift);
            refs.push((key.clone(), j, target, decoy));
        }
        drop(ms1_d);
        for (t_run, run) in runs.iter().enumerate() {
            if t_run == d || alignment.maps[t_run].is_none() {
                continue;
            }
            let mut truth: HashMap<&GlycoKey, usize> = HashMap::new();
            for (key, &j) in &alignment.donor[t_run] {
                truth.insert(key, j);
            }
            if truth.is_empty() {
                continue;
            }
            let ms1 = run
                .load()
                .map_err(|e| format!("--glyco-mbr-match-train: {e}"))?;
            let (mut pair_pos, mut pair_neg) = (0usize, 0usize);
            for (key, donor_row, reference, decoy) in &refs {
                let Some(&truth_row) = truth.get(key) else {
                    continue;
                };
                let donor = &rows[*donor_row];
                let true_row = &rows[truth_row];
                let (Some(df), Some(tf)) = (&donor.feature, &true_row.feature) else {
                    continue;
                };
                let Some((predicted_rt, sigma)) = alignment.predict(d, t_run, df.apex_rt) else {
                    continue;
                };
                let search = (3.0 * sigma).clamp(MBR_MATCH_MIN_SEARCH_S, MBR_MATCH_MAX_SEARCH_S);
                let Some(candidates) = glyco_match::score_candidates(
                    &ms1,
                    reference,
                    reference.window_start(predicted_rt),
                    search,
                    settings.params.tol,
                ) else {
                    continue;
                };
                let scored = glyco_match::scored_top(&candidates);
                if scored.is_empty() {
                    continue;
                }
                let truth_width = (tf.rt_stop - tf.rt_start).max(0.0);
                let min_overlap = 0.5 * truth_width;
                let overlaps: Vec<f64> = scored
                    .iter()
                    .map(|m| (m.stop_rt.min(tf.rt_stop) - m.start_rt.max(tf.rt_start)).max(0.0))
                    .collect();
                if !overlaps.iter().any(|&o| o >= min_overlap && o > 0.0) {
                    continue;
                }
                let group = &key.0;
                for (m, &overlap) in scored.iter().zip(&overlaps) {
                    let label = u8::from(overlap >= min_overlap && overlap > 0.0);
                    write_training_row(&mut out, label, m, group, d, t_run)
                        .map_err(|e| format!("{}: {e}", path.display()))?;
                    if label == 1 {
                        pair_pos += 1;
                    } else {
                        pair_neg += 1;
                    }
                }
                if let Some(decoy) = decoy {
                    if let Some(candidates) = glyco_match::score_candidates(
                        &ms1,
                        decoy,
                        decoy.window_start(predicted_rt),
                        search,
                        settings.params.tol,
                    ) {
                        for m in glyco_match::scored_top(&candidates) {
                            write_training_row(&mut out, 0, &m, group, d, t_run)
                                .map_err(|e| format!("{}: {e}", path.display()))?;
                            pair_neg += 1;
                        }
                    }
                }
                pairs += 1;
            }
            positives += pair_pos;
            negatives += pair_neg;
            log::debug!(
                target: "glyco-mbr-match-train",
                "runs {d}->{t_run}: {pair_pos} positive, \
                 {pair_neg} negative candidate windows"
            );
        }
    }
    out.flush()
        .map_err(|e| format!("{}: {e}", path.display()))?;
    log::info!(
        target: "glyco-mbr-match-train",
        "{} glycopeptide pairs, {positives} positive and \
         {negatives} negative rows written to {}",
        pairs,
        path.display()
    );
    Ok(())
}

fn write_training_row(
    out: &mut impl Write,
    label: u8,
    m: &glyco_match::Match,
    group: &str,
    donor_run: usize,
    target_run: usize,
) -> std::io::Result<()> {
    write!(out, "{label}")?;
    for f in m.features {
        write!(out, "\t{f:.6}")?;
    }
    writeln!(out, "\t{group}\t{donor_run}>{target_run}")
}

fn output_path(pin: &Path, kind: &str) -> PathBuf {
    let s = pin.to_string_lossy();
    match s.strip_suffix(".glyco.pin") {
        Some(base) => PathBuf::from(format!("{base}.glyco.{kind}.tsv")),
        None => PathBuf::from(format!("{s}.glyco.{kind}.tsv")),
    }
}

fn write_feature_table(
    path: &Path,
    rows: &[FeatureRow],
    settings: &Settings,
) -> Result<(), String> {
    let mut out = std::io::BufWriter::new(
        std::fs::File::create(path).map_err(|e| format!("{}: {e}", path.display()))?,
    );
    writeln!(
        out,
        "Run\tSpecId\tScan\tPeptide\tGlycosite\tGlycan\tTopology\tProteins\tCharge\t\
         IsotopeOffset\tAnchorRT\tApexRT\tRTStart\tRTStop\tRTDelta_s\tMonoMz\tTheoMz\tPpm\t\
         Area\tMonoArea\tApexIntensity\tCosine\tNIsotopesFound\tNScans\tFeatureScore\t\
         FeatureQValue\tDecoyCosine\tAccepted\tEvidence\tMatchScore"
    )
    .map_err(|e| format!("{}: {e}", path.display()))?;
    for row in rows {
        let t = &row.target;
        let (apex, start, stop, area, mono_area, apex_int, cos, n_found, n_scans, delta) =
            match &row.feature {
                Some(f) => (
                    f.apex_rt,
                    f.rt_start,
                    f.rt_stop,
                    f.area,
                    f.isotope_areas.first().copied().unwrap_or(0.0),
                    f.apex_intensity as f64,
                    f.cosine as f64,
                    f.n_isotopes_found,
                    f.n_scans,
                    f.rt_delta_s,
                ),
                None => (
                    f64::NAN,
                    f64::NAN,
                    f64::NAN,
                    f64::NAN,
                    f64::NAN,
                    f64::NAN,
                    f64::NAN,
                    0,
                    0,
                    f64::NAN,
                ),
            };
        let decoy_cos = row.decoy.as_ref().map_or(f64::NAN, |f| f.cosine as f64);
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t\
             {}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            t.run,
            t.spec_id,
            t.scan,
            t.peptide,
            t.site,
            t.glycan,
            t.topology,
            t.proteins,
            t.charge,
            t.isotope_offset,
            fmt(t.anchor_rt),
            fmt(apex),
            fmt(start),
            fmt(stop),
            fmt(delta),
            fmt(t.mono_mz),
            fmt(t.theo_mz),
            fmt(t.ppm),
            fmt(area),
            fmt(mono_area),
            fmt(apex_int),
            fmt(cos),
            n_found,
            n_scans,
            fmt(row.score.map_or(f64::NAN, f64::from)),
            fmt(row.q.unwrap_or(f64::NAN)),
            fmt(decoy_cos),
            u8::from(accepted(row, settings)),
            t.evidence.as_str(),
            fmt(row.match_score.map_or(f64::NAN, f64::from)),
        )
        .map_err(|e| format!("{}: {e}", path.display()))?;
    }
    out.flush().map_err(|e| format!("{}: {e}", path.display()))
}

struct Kept {
    glycopeptides: usize,
    sites: usize,
}

fn write_rollup_tables(
    pin: &Path,
    rows: &[FeatureRow],
    run_names: &[&str],
    settings: &Settings,
) -> Result<Kept, String> {
    #[derive(PartialEq, Eq, Hash)]
    struct Key(String, String, String, u8);
    let mut best: HashMap<Key, usize> = HashMap::new();
    for (i, row) in rows.iter().enumerate() {
        if !accepted(row, settings) || row.feature.is_none() {
            continue;
        }
        let key = Key(
            row.target.run.clone(),
            row.target.peptide.clone(),
            row.target.glycan.clone(),
            row.target.charge,
        );
        let better = match best.get(&key) {
            None => true,
            Some(&j) => {
                let sj = rows[j].score.unwrap_or(f32::NEG_INFINITY);
                let si = row.score.unwrap_or(f32::NEG_INFINITY);
                si > sj
                    || (si == sj
                        && row.feature.as_ref().map_or(0.0, |f| f.area)
                            > rows[j].feature.as_ref().map_or(0.0, |f| f.area))
            }
        };
        if better {
            best.insert(key, i);
        }
    }

    struct Glyco {
        area: f64,
        best_cosine: f32,
        best_q: f64,
        charges: BTreeSet<u8>,
        scans: BTreeSet<i32>,
        transferred: usize,
        proteins: String,
        site: String,
    }
    let mut glyco: HashMap<(String, String, String), Glyco> = HashMap::new();
    for &i in best.values() {
        let row = &rows[i];
        let Some(f) = row.feature.as_ref() else {
            continue;
        };
        let key = (
            row.target.run.clone(),
            row.target.peptide.clone(),
            row.target.glycan.clone(),
        );
        let e = glyco.entry(key).or_insert_with(|| Glyco {
            area: 0.0,
            best_cosine: 0.0,
            best_q: 1.0,
            charges: BTreeSet::new(),
            scans: BTreeSet::new(),
            transferred: 0,
            proteins: row.target.proteins.clone(),
            site: row.target.site.clone(),
        });
        e.area += f.area;
        e.best_cosine = e.best_cosine.max(f.cosine);
        e.best_q = e.best_q.min(row.q.unwrap_or(1.0));
        e.charges.insert(row.target.charge);
        match row.target.evidence {
            Evidence::Identified => {
                e.scans.insert(row.target.scan);
            }
            Evidence::Transferred => e.transferred += 1,
        }
    }
    let quant_path = output_path(pin, "quant");
    let mut out = std::io::BufWriter::new(
        std::fs::File::create(&quant_path).map_err(|e| format!("{}: {e}", quant_path.display()))?,
    );
    writeln!(
        out,
        "Run\tPeptide\tGlycosite\tGlycan\tProteins\tCharges\tNScans\tNTransferred\tBestCosine\t\
         BestQValue\tArea"
    )
    .map_err(|e| format!("{}: {e}", quant_path.display()))?;
    let mut glyco_rows: Vec<_> = glyco.iter().collect();
    glyco_rows.sort_by(|a, b| {
        a.0 .0
            .cmp(&b.0 .0)
            .then_with(|| a.0 .1.cmp(&b.0 .1))
            .then_with(|| a.0 .2.cmp(&b.0 .2))
    });
    for ((run, peptide, glycan), g) in &glyco_rows {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            run,
            peptide,
            g.site,
            glycan,
            g.proteins,
            g.charges
                .iter()
                .map(|c| c.to_string())
                .collect::<Vec<_>>()
                .join(";"),
            g.scans.len(),
            g.transferred,
            fmt(g.best_cosine as f64),
            fmt(g.best_q),
            fmt(g.area),
        )
        .map_err(|e| format!("{}: {e}", quant_path.display()))?;
    }
    out.flush()
        .map_err(|e| format!("{}: {e}", quant_path.display()))?;

    let mut site_total: HashMap<(String, String), f64> = HashMap::new();
    for ((run, peptide, _), g) in &glyco_rows {
        *site_total
            .entry((run.clone(), peptide.clone()))
            .or_insert(0.0) += g.area;
    }
    let site_path = output_path(pin, "site");
    let mut out = std::io::BufWriter::new(
        std::fs::File::create(&site_path).map_err(|e| format!("{}: {e}", site_path.display()))?,
    );
    writeln!(
        out,
        "Run\tPeptide\tGlycosite\tProteins\tGlycan\tArea\tSiteGlycoArea\tGlycoformShare"
    )
    .map_err(|e| format!("{}: {e}", site_path.display()))?;
    for ((run, peptide, glycan), g) in &glyco_rows {
        let total = site_total
            .get(&(run.clone(), peptide.clone()))
            .copied()
            .unwrap_or(0.0);
        let share = if total > 0.0 {
            g.area / total
        } else {
            f64::NAN
        };
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            run,
            peptide,
            g.site,
            g.proteins,
            glycan,
            fmt(g.area),
            fmt(total),
            fmt(share),
        )
        .map_err(|e| format!("{}: {e}", site_path.display()))?;
    }
    out.flush()
        .map_err(|e| format!("{}: {e}", site_path.display()))?;

    struct MatrixRow {
        site: String,
        proteins: String,
        areas: HashMap<String, f64>,
    }
    let mut matrix: HashMap<(String, String), MatrixRow> = HashMap::new();
    for ((run, peptide, glycan), g) in &glyco_rows {
        let e = matrix
            .entry((peptide.clone(), glycan.clone()))
            .or_insert_with(|| MatrixRow {
                site: g.site.clone(),
                proteins: g.proteins.clone(),
                areas: HashMap::new(),
            });
        e.areas.insert(run.clone(), g.area);
    }
    let matrix_path = output_path(pin, "matrix");
    let mut out = std::io::BufWriter::new(
        std::fs::File::create(&matrix_path)
            .map_err(|e| format!("{}: {e}", matrix_path.display()))?,
    );
    write!(out, "Peptide\tGlycosite\tGlycan\tProteins\tNRuns")
        .map_err(|e| format!("{}: {e}", matrix_path.display()))?;
    for name in run_names {
        write!(out, "\t{name}").map_err(|e| format!("{}: {e}", matrix_path.display()))?;
    }
    writeln!(out).map_err(|e| format!("{}: {e}", matrix_path.display()))?;
    let mut matrix_rows: Vec<_> = matrix.iter().collect();
    matrix_rows.sort_by(|a, b| a.0 .0.cmp(&b.0 .0).then_with(|| a.0 .1.cmp(&b.0 .1)));
    for ((peptide, glycan), m) in matrix_rows {
        write!(
            out,
            "{}\t{}\t{}\t{}\t{}",
            peptide,
            m.site,
            glycan,
            m.proteins,
            m.areas.len()
        )
        .map_err(|e| format!("{}: {e}", matrix_path.display()))?;
        for name in run_names {
            match m.areas.get(*name) {
                Some(area) => write!(out, "\t{}", fmt(*area)),
                None => write!(out, "\t"),
            }
            .map_err(|e| format!("{}: {e}", matrix_path.display()))?;
        }
        writeln!(out).map_err(|e| format!("{}: {e}", matrix_path.display()))?;
    }
    out.flush()
        .map_err(|e| format!("{}: {e}", matrix_path.display()))?;

    let mut per_protein: HashMap<(String, String), Vec<f64>> = HashMap::new();
    for ((run, _, _), g) in &glyco_rows {
        for acc in g.proteins.split(';').filter(|s| !s.is_empty()) {
            per_protein
                .entry((run.clone(), acc.to_string()))
                .or_default()
                .push(g.area);
        }
    }
    let protein_path = output_path(pin, "protein");
    let mut out = std::io::BufWriter::new(
        std::fs::File::create(&protein_path)
            .map_err(|e| format!("{}: {e}", protein_path.display()))?,
    );
    writeln!(out, "Run\tProtein\tGlycopeptides\tSummedArea\tTop3")
        .map_err(|e| format!("{}: {e}", protein_path.display()))?;
    let mut protein_matrix: HashMap<String, HashMap<String, (f64, f64)>> = HashMap::new();
    for ((run, acc), areas) in &per_protein {
        let mut sorted = areas.clone();
        sorted.sort_by(|a, b| b.total_cmp(a));
        let top3: f64 = sorted.iter().take(3).sum();
        protein_matrix
            .entry(acc.clone())
            .or_default()
            .insert(run.clone(), (areas.iter().sum::<f64>(), top3));
    }
    let protein_matrix_path = output_path(pin, "protein_matrix");
    let mut pout = std::io::BufWriter::new(
        std::fs::File::create(&protein_matrix_path)
            .map_err(|e| format!("{}: {e}", protein_matrix_path.display()))?,
    );
    write!(pout, "Protein").map_err(|e| format!("{}: {e}", protein_matrix_path.display()))?;
    for name in run_names {
        write!(pout, "\t{name}_SummedArea\t{name}_Top3")
            .map_err(|e| format!("{}: {e}", protein_matrix_path.display()))?;
    }
    writeln!(pout).map_err(|e| format!("{}: {e}", protein_matrix_path.display()))?;
    let mut protein_matrix_rows: Vec<_> = protein_matrix.iter().collect();
    protein_matrix_rows.sort_by(|a, b| a.0.cmp(b.0));
    for (acc, per_run) in protein_matrix_rows {
        write!(pout, "{acc}").map_err(|e| format!("{}: {e}", protein_matrix_path.display()))?;
        for name in run_names {
            match per_run.get(*name) {
                Some((sum, top3)) => write!(pout, "\t{}\t{}", fmt(*sum), fmt(*top3)),
                None => write!(pout, "\t\t"),
            }
            .map_err(|e| format!("{}: {e}", protein_matrix_path.display()))?;
        }
        writeln!(pout).map_err(|e| format!("{}: {e}", protein_matrix_path.display()))?;
    }
    pout.flush()
        .map_err(|e| format!("{}: {e}", protein_matrix_path.display()))?;

    let mut protein_rows: Vec<_> = per_protein.into_iter().collect();
    protein_rows.sort_by(|a, b| a.0 .0.cmp(&b.0 .0).then_with(|| a.0 .1.cmp(&b.0 .1)));
    for ((run, acc), mut areas) in protein_rows {
        let total: f64 = areas.iter().sum();
        areas.sort_by(|a, b| b.total_cmp(a));
        let top3: f64 = areas.iter().take(3).sum();
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}",
            run,
            acc,
            areas.len(),
            fmt(total),
            fmt(top3),
        )
        .map_err(|e| format!("{}: {e}", protein_path.display()))?;
    }
    out.flush()
        .map_err(|e| format!("{}: {e}", protein_path.display()))?;

    Ok(Kept {
        glycopeptides: glyco_rows.len(),
        sites: site_total.len(),
    })
}

fn composition_label(comp: &GlycanComp) -> String {
    format!(
        "HexNAc{}Hex{}Fuc{}NeuAc{}NeuGc{}",
        comp.hexnac, comp.hex, comp.fuc, comp.neuac, comp.neugc
    )
}

fn site_label(residues: &[u8], post: u8, core: GlycanCore) -> String {
    match core {
        GlycanCore::NGlycan => match glyco::sequon::glycosites(residues, post).as_slice() {
            [only] => format!("N{}", only + 1),
            _ => "N?".to_string(),
        },
        GlycanCore::OGlycan => match glyco::sequon::all_ser_thr_sites(residues).as_slice() {
            [only] => format!("{}{}", residues[*only] as char, only + 1),
            _ => "S/T?".to_string(),
        },
    }
}

fn proteins_of(
    psm: &search::psm::PsmMatch,
    candidates: &[Candidate],
    index: &SearchIndex,
) -> String {
    psm.candidate_idxs
        .iter()
        .filter_map(|&c| candidates.get(c as usize))
        .map(|c| {
            index.protein_at(c.protein_index).map_or_else(
                || format!("PROT_{}", c.protein_index),
                |p| p.accession.clone(),
            )
        })
        .collect::<Vec<_>>()
        .join(";")
}

fn fmt(v: f64) -> String {
    if v.is_finite() {
        format!("{v:.6}")
    } else {
        String::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_extracted_envelope_starts_at_the_isotope_offset() {
        let full = Formula::averagine(3000.0).envelope(ENVELOPE_ISOTOPES);
        let mono = extracted_envelope(&full, 0);
        let m2 = extracted_envelope(&full, 2);
        assert!(m2[0] > mono[0], "{mono:?} {m2:?}");
        assert!((m2.iter().sum::<f64>() - 1.0).abs() < 1e-9);
        assert!(m2.len() <= MAX_ISOTOPES && m2.len() >= MIN_ISOTOPES);
    }

    #[test]
    fn an_offset_past_the_envelope_is_a_single_isotope() {
        assert_eq!(extracted_envelope(&[0.5, 0.5], 9), vec![1.0]);
    }

    #[test]
    fn a_large_glycopeptide_extracts_more_isotopes_than_a_peptide() {
        let peptide = Formula::peptide(b"PEPTIDEK".iter().map(|&r| (r, 0.0)));
        let mut glyco = Formula::peptide(b"IGNTTNA".iter().map(|&r| (r, 0.0)));
        glyco.add(Formula::glycan(8, 2, 0, 0, 0));
        let a = extracted_envelope(&peptide.envelope(ENVELOPE_ISOTOPES), 0);
        let b = extracted_envelope(&glyco.envelope(ENVELOPE_ISOTOPES), 0);
        assert!(b.len() >= a.len(), "{} vs {}", a.len(), b.len());
        assert!(b.len() > 4, "a 2 kDa glycopeptide needs more than 4: {b:?}");
    }

    #[test]
    fn composition_labels_match_the_pin_tag() {
        let comp = GlycanComp::from_tree_str("(N(N(H(H(H(H))(H)(H))(H(H)))))").unwrap();
        assert_eq!(
            composition_label(&comp),
            format!(
                "HexNAc{}Hex{}Fuc{}NeuAc{}NeuGc{}",
                comp.hexnac, comp.hex, comp.fuc, comp.neuac, comp.neugc
            )
        );
    }

    #[test]
    fn output_paths_sit_beside_the_pin() {
        let p = output_path(Path::new("/d/run.glyco.pin"), "quant");
        assert_eq!(p, PathBuf::from("/d/run.glyco.quant.tsv"));
    }

    fn settings() -> Settings {
        Settings {
            params: LfqParams::default(),
            transfer_fdr: 0.05,
            core: GlycanCore::NGlycan,
            debug: false,
            mbr: false,
            match_train: None,
        }
    }

    fn target(mono_mz: f64, charge: u8, anchor_rt: f64) -> Target {
        let mut f = Formula::peptide(b"IGNTTNA".iter().map(|&r| (r, 0.0)));
        f.add(Formula::glycan(8, 2, 0, 0, 0));
        Target {
            run_idx: 0,
            run: "r1".into(),
            spec_id: "scan=1_glyco_1_1".into(),
            scan: 1,
            evidence: Evidence::Identified,
            window_s: 60.0,
            peptide: "K.IGNTTNA.-".into(),
            glycan: "HexNAc2Hex8Fuc0NeuAc0NeuGc0".into(),
            topology: "(N(N(H(H(H(H))(H)(H))(H(H)))))".into(),
            site: "N3".into(),
            proteins: "P1".into(),
            charge,
            isotope_offset: 0,
            mono_mz,
            theo_mz: mono_mz,
            ppm: 0.0,
            anchor_rt,
            envelope: extracted_envelope(&f.envelope(ENVELOPE_ISOTOPES), 0),
            search_s: 0.0,
            donor: None,
        }
    }

    fn synthetic_run(t: &Target) -> quant::ms1_index::Ms1RunIndex {
        let z = t.charge as f64;
        let scans = (0..200)
            .map(|i| {
                let rt = i as f64;
                let g = (-0.5 * ((rt - 100.0) / 4.0).powi(2)).exp();
                let mut peaks: Vec<(f64, f32)> = Vec::new();
                for (k, e) in t.envelope.iter().enumerate() {
                    let inten = (1e6 * e * g) as f32;
                    if inten > 1.0 {
                        peaks.push((t.mono_mz + k as f64 * ISOTOPE / z, inten));
                    }
                }
                peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
                model::scan::Ms1Scan {
                    rt,
                    peaks,
                    ..model::scan::Ms1Scan::default()
                }
            })
            .collect();
        quant::ms1_index::Ms1RunIndex::new(scans)
    }

    #[test]
    fn measures_a_synthetic_glycopeptide_envelope() {
        let t = target(800.0, 2, 97.0);
        let index = synthetic_run(&t);
        let m = measure_identified(&index, &t, &LfqParams::default());
        let f = m.feature.as_ref().expect("feature");
        assert!((f.apex_rt - 100.0).abs() <= 1.0, "{f:?}");
        assert!(f.cosine > 0.99, "{f:?}");
        assert!(f.n_isotopes_found >= MIN_ISOTOPES as u32, "{f:?}");
        let expected = 1e6 * 4.0 * (2.0 * std::f64::consts::PI).sqrt();
        assert!((f.area / expected - 1.0).abs() < 0.03, "{f:?}");
        assert!(m.decoy.is_none(), "{:?}", m.decoy);
    }

    #[test]
    fn an_isotope_offset_extracts_the_observed_peak() {
        let mut t = target(800.0, 2, 97.0);
        t.mono_mz += ISOTOPE / 2.0;
        t.isotope_offset = 1;
        let mut f = Formula::peptide(b"IGNTTNA".iter().map(|&r| (r, 0.0)));
        f.add(Formula::glycan(8, 2, 0, 0, 0));
        t.envelope = extracted_envelope(&f.envelope(ENVELOPE_ISOTOPES), 1);
        let index = synthetic_run(&t);
        let m = measure_identified(&index, &t, &LfqParams::default());
        let f = m.feature.as_ref().expect("feature");
        assert!((f.apex_rt - 100.0).abs() <= 1.0, "{f:?}");
        assert!(f.cosine > 0.99, "{f:?}");
    }

    fn synthetic_run_mono_only(t: &Target) -> quant::ms1_index::Ms1RunIndex {
        let scans = (0..200)
            .map(|i| {
                let rt = i as f64;
                let g = (-0.5 * ((rt - 100.0) / 4.0).powi(2)).exp();
                let inten = (1e6 * t.envelope[0] * g) as f32;
                let peaks = if inten > 1.0 {
                    vec![(t.mono_mz, inten)]
                } else {
                    Vec::new()
                };
                model::scan::Ms1Scan {
                    rt,
                    peaks,
                    ..model::scan::Ms1Scan::default()
                }
            })
            .collect();
        quant::ms1_index::Ms1RunIndex::new(scans)
    }

    #[test]
    fn mono_intensity_survives_dropped_isotope_traces() {
        let t = target(800.0, 2, 97.0);
        let full = measure_identified(&synthetic_run(&t), &t, &LfqParams::default());
        let mono_only = measure_identified(&synthetic_run_mono_only(&t), &t, &LfqParams::default());
        let a = full.feature.as_ref().expect("full feature").area;
        let b = mono_only.feature.as_ref().expect("mono-only feature").area;
        assert!((a / b - 1.0).abs() < 0.03, "{a} vs {b}");
        let f = mono_only.feature.as_ref().unwrap();
        let expected = f.isotope_areas[0] / t.envelope[0];
        assert!((f.area / expected - 1.0).abs() < 1e-12, "{f:?}");
    }

    fn feature(area: f64, cosine: f32) -> FeatureQuant {
        FeatureQuant {
            apex_rt: 100.0,
            rt_start: 90.0,
            rt_stop: 110.0,
            area,
            apex_intensity: 1e6,
            isotope_areas: vec![area],
            cosine,
            n_isotopes_found: 4,
            rt_delta_s: 3.0,
            n_scans: 21,
        }
    }

    fn row(t: Target, area: f64, score: f32, q: f64) -> FeatureRow {
        FeatureRow {
            target: t,
            feature: Some(feature(area, 0.99)),
            decoy: None,
            score: Some(score),
            q: Some(q),
            match_score: None,
            decoy_match: None,
        }
    }

    fn read(path: &Path) -> Vec<Vec<String>> {
        let text = std::fs::read_to_string(path).unwrap();
        text.lines()
            .map(|l| l.split('\t').map(str::to_string).collect())
            .collect()
    }

    #[test]
    fn the_rollup_dedupes_scans_sums_charges_and_drops_rejected_features() {
        let dir = tempfile::tempdir().unwrap();
        let pin = dir.path().join("run.glyco.pin");
        let mut a = row(target(800.0, 2, 97.0), 100.0, 0.9, 0.001);
        a.target.scan = 1;
        let mut b = row(target(800.0, 2, 97.0), 100.0, 0.5, 0.001);
        b.target.scan = 2;
        let mut c = row(target(800.0, 3, 97.0), 50.0, 0.8, 0.002);
        c.target.scan = 3;
        c.target.charge = 3;
        let mut rejected = row(target(800.0, 2, 97.0), 999.0, 0.99, 0.5);
        rejected.target.evidence = Evidence::Transferred;
        rejected.target.scan = 0;
        let rows = vec![a, b, c, rejected];
        let kept = write_rollup_tables(&pin, &rows, &["r1"], &settings()).unwrap();
        assert_eq!(kept.glycopeptides, 1);
        assert_eq!(kept.sites, 1);

        let q = read(&dir.path().join("run.glyco.quant.tsv"));
        assert_eq!(q.len(), 2);
        assert_eq!(
            q[0],
            vec![
                "Run",
                "Peptide",
                "Glycosite",
                "Glycan",
                "Proteins",
                "Charges",
                "NScans",
                "NTransferred",
                "BestCosine",
                "BestQValue",
                "Area"
            ]
        );
        let r = &q[1];
        assert_eq!(r[0], "r1");
        assert_eq!(r[1], "K.IGNTTNA.-");
        assert_eq!(r[5], "2;3");
        assert_eq!(r[6], "2");
        assert_eq!(r[7], "0", "an identified-only row has no transfers");
        assert!((r[10].parse::<f64>().unwrap() - 150.0).abs() < 1e-6);

        let s = read(&dir.path().join("run.glyco.site.tsv"));
        assert_eq!(s.len(), 2);
        assert!((s[1][7].parse::<f64>().unwrap() - 1.0).abs() < 1e-9);

        let p = read(&dir.path().join("run.glyco.protein.tsv"));
        assert_eq!(p.len(), 2);
        assert!((p[1][3].parse::<f64>().unwrap() - 150.0).abs() < 1e-6);
        assert!((p[1][4].parse::<f64>().unwrap() - 150.0).abs() < 1e-6);
    }

    #[test]
    fn an_identified_feature_is_not_vetoed_by_its_decoy_q() {
        let dir = tempfile::tempdir().unwrap();
        let pin = dir.path().join("run.glyco.pin");
        let identified = row(target(800.0, 2, 97.0), 100.0, 0.9, 0.5);
        let mut transfer = row(target(800.0, 2, 98.0), 50.0, 0.9, 0.5);
        transfer.target.evidence = Evidence::Transferred;
        transfer.target.scan = 0;
        transfer.target.spec_id = "donor_mbr_1".into();
        let kept =
            write_rollup_tables(&pin, &[identified, transfer], &["r1"], &settings()).unwrap();
        assert_eq!(kept.glycopeptides, 1);

        let q = read(&dir.path().join("run.glyco.quant.tsv"));
        let r = &q[1];
        assert_eq!(r[7], "0", "the transfer is q-gated out: {r:?}");
        assert!((r[10].parse::<f64>().unwrap() - 100.0).abs() < 1e-6);
    }

    #[test]
    fn the_matrices_pivot_runs_and_leave_missing_cells_empty() {
        let dir = tempfile::tempdir().unwrap();
        let pin = dir.path().join("run.glyco.pin");
        let mut a = row(target(800.0, 2, 97.0), 100.0, 0.9, 0.001);
        a.target.run = "r1".into();
        let mut b = row(target(800.0, 2, 97.0), 40.0, 0.9, 0.001);
        b.target.run = "r2".into();
        let mut other = row(target(800.0, 2, 97.0), 7.0, 0.9, 0.001);
        other.target.run = "r2".into();
        other.target.glycan = "HexNAc3Hex4Fuc0NeuAc0NeuGc0".into();
        let kept = write_rollup_tables(&pin, &[a, b, other], &["r1", "r2"], &settings()).unwrap();
        assert_eq!((kept.glycopeptides, kept.sites), (3, 2));

        let m = read(&dir.path().join("run.glyco.matrix.tsv"));
        assert_eq!(
            m[0],
            vec![
                "Peptide",
                "Glycosite",
                "Glycan",
                "Proteins",
                "NRuns",
                "r1",
                "r2"
            ]
        );
        assert_eq!(m[1][2], "HexNAc2Hex8Fuc0NeuAc0NeuGc0");
        assert_eq!(m[1][4], "2");
        assert!((m[1][5].parse::<f64>().unwrap() - 100.0).abs() < 1e-6);
        assert!((m[1][6].parse::<f64>().unwrap() - 40.0).abs() < 1e-6);
        assert_eq!(m[2][2], "HexNAc3Hex4Fuc0NeuAc0NeuGc0");
        assert_eq!(m[2][4], "1");
        assert_eq!(m[2][5], "", "r1 did not quantify this glycoform: {m:?}");
        assert!((m[2][6].parse::<f64>().unwrap() - 7.0).abs() < 1e-6);

        let p = read(&dir.path().join("run.glyco.protein_matrix.tsv"));
        assert_eq!(
            p[0],
            vec![
                "Protein",
                "r1_SummedArea",
                "r1_Top3",
                "r2_SummedArea",
                "r2_Top3"
            ]
        );
        assert_eq!(p[1][0], "P1");
        assert!((p[1][1].parse::<f64>().unwrap() - 100.0).abs() < 1e-6);
        assert!((p[1][3].parse::<f64>().unwrap() - 47.0).abs() < 1e-6);
        assert!((p[1][4].parse::<f64>().unwrap() - 47.0).abs() < 1e-6);
    }

    #[test]
    fn the_rt_map_recovers_a_linear_shift() {
        let pairs: Vec<(f64, f64)> = (0..100)
            .map(|i| {
                let x = i as f64 * 3.0;
                (x, 4.0 + 1.1 * x)
            })
            .collect();
        let map = fit_rt_map(&pairs).expect("map");
        assert!(map.sigma_s < 1e-9, "{map:?}");
        for x in [0.0, 123.0, 297.0] {
            let want = 4.0 + 1.1 * x;
            assert!((map.forward(x) - want).abs() < 1e-6, "{x}: {map:?}");
        }
        let rt_ref = map.forward(123.0);
        assert!((map.inverse(rt_ref).unwrap() - 123.0).abs() < 1e-6);
    }

    #[test]
    fn the_rt_map_follows_a_non_linear_warp() {
        let pairs: Vec<(f64, f64)> = (0..900)
            .map(|i| {
                let x = i as f64 * 4.0;
                let phase = x / 3600.0 * std::f64::consts::PI;
                (x, x + 40.0 * (phase).sin())
            })
            .collect();
        let map = fit_rt_map(&pairs).expect("map");
        assert!(map.sigma_s < 10.0, "{map:?}");
        for x in [400.0, 1800.0, 3200.0] {
            let want = x + 40.0 * (x / 3600.0 * std::f64::consts::PI).sin();
            assert!(
                (map.forward(x) - want).abs() < 5.0,
                "{x}: {} vs {want}",
                map.forward(x)
            );
            assert!((map.inverse(want).unwrap() - x).abs() < 5.0, "{x}");
        }
    }

    #[test]
    fn the_rt_map_drops_gross_outliers() {
        let mut pairs: Vec<(f64, f64)> = (0..900)
            .map(|i| {
                let x = i as f64 * 3.0;
                (x, 4.0 + 1.1 * x)
            })
            .collect();
        for i in (0..900).step_by(180) {
            pairs[i].1 += 1000.0;
        }
        let map = fit_rt_map(&pairs).expect("map");
        for x in [100.0, 1400.0, 2600.0] {
            assert!(
                (map.forward(x) - (4.0 + 1.1 * x)).abs() < 1.0,
                "{x}: {map:?}"
            );
        }
        assert!(map.sigma_s < 0.5, "{map:?}");
    }

    fn identified_row(
        key: usize,
        run_idx: usize,
        run: &str,
        apex_rt: f64,
        q: f64,
        score: f32,
    ) -> FeatureRow {
        let mut t = target(800.0, 2, apex_rt);
        t.run_idx = run_idx;
        t.run = run.to_string();
        t.peptide = format!("K.PEPTIDE{key}K.-");
        t.spec_id = format!("{run}_scan{key}");
        let mut f = feature(100.0, 0.99);
        f.apex_rt = apex_rt;
        FeatureRow {
            target: t,
            feature: Some(f),
            decoy: None,
            score: Some(score),
            q: Some(q),
            match_score: None,
            decoy_match: None,
        }
    }

    #[test]
    fn transfers_fill_the_runs_that_missed_the_glycopeptide() {
        let n = MBR_MIN_ANCHORS + 10;
        let mut rows = Vec::new();
        for k in 0..n {
            let rt = 100.0 + 2.0 * k as f64;
            rows.push(identified_row(k, 0, "r1", rt, 0.01, 10.0));
            rows.push(identified_row(k, 1, "r2", 0.8 * rt + 50.0, 0.01, 10.0));
        }
        rows.push(identified_row(n, 0, "r1", 500.0, 0.01, 20.0));
        rows.push(identified_row(n + 1, 0, "r1", 520.0, 0.01, 20.0));
        let mut rejected = identified_row(n + 1, 1, "r2", 466.0, 0.5, 20.0);
        rejected.score = None;
        rows.push(rejected);

        let settings = settings();
        let alignment = Alignment::build(&rows, &["r1", "r2"], &settings);
        let (transfers, counts) = build_transfers(&rows, &["r1", "r2"], &settings, &alignment);
        assert_eq!(counts.maps, 2, "both runs have an RT map");
        assert_eq!(counts.anchors, n);
        assert_eq!(counts.candidates, 1);
        let t = &transfers[0];
        assert_eq!(t.run_idx, 1);
        assert_eq!(t.evidence, Evidence::Transferred);
        assert_eq!(t.scan, 0);
        assert!(t.spec_id.ends_with("_mbr_r2"), "{}", t.spec_id);
        let expected = 0.8 * 500.0 + 50.0;
        assert!(
            (t.anchor_rt - expected).abs() < 1e-6,
            "{} vs {expected}",
            t.anchor_rt
        );
        assert!(
            (t.window_s - MBR_MIN_WINDOW_S).abs() < 1e-9,
            "{:?}",
            t.window_s
        );
    }

    #[test]
    fn too_few_shared_glycopeptides_mean_no_map_and_no_transfers() {
        let n = MBR_MIN_ANCHORS - 1;
        let mut rows = Vec::new();
        for k in 0..n {
            let rt = 100.0 + 2.0 * k as f64;
            rows.push(identified_row(k, 0, "r1", rt, 0.01, 10.0));
            rows.push(identified_row(k, 1, "r2", rt + 10.0, 0.01, 10.0));
        }
        rows.push(identified_row(999, 0, "r1", 500.0, 0.01, 20.0));
        let settings = settings();
        let alignment = Alignment::build(&rows, &["r1", "r2"], &settings);
        let (transfers, counts) = build_transfers(&rows, &["r1", "r2"], &settings, &alignment);
        assert_eq!(counts.maps, 1, "only the reference run is mapped");
        assert_eq!(counts.candidates, 0);
        assert!(transfers.is_empty(), "{transfers:?}");
    }

    #[test]
    fn a_transferred_row_is_counted_and_not_given_a_scan_number() {
        let dir = tempfile::tempdir().unwrap();
        let pin = dir.path().join("run.glyco.pin");
        let mut a = row(target(800.0, 2, 97.0), 100.0, 0.9, 0.001);
        a.target.run = "r1".into();
        a.target.scan = 7;
        let mut b = row(target(800.0, 2, 97.0), 40.0, 0.9, 0.001);
        b.target.run = "r2".into();
        b.target.scan = 0;
        b.target.evidence = Evidence::Transferred;
        b.target.spec_id = "r1_scan7_glyco_1_1_mbr_r2".into();
        write_rollup_tables(&pin, &[a, b], &["r1", "r2"], &settings()).unwrap();

        let q = read(&dir.path().join("run.glyco.quant.tsv"));
        let r2 = q.iter().find(|r| r[0] == "r2").expect("r2 row");
        assert_eq!(r2[6], "0", "a transfer has no identifying scan: {r2:?}");
        assert_eq!(r2[7], "1");
        assert!((r2[10].parse::<f64>().unwrap() - 40.0).abs() < 1e-6);
        let r1 = q.iter().find(|r| r[0] == "r1").expect("r1 row");
        assert_eq!(r1[6], "1");
        assert_eq!(r1[7], "0");
    }
}
