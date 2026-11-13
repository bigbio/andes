//! Intact-glycopeptide label-free quantification from MS1, modelled on the
//! published method of Zeng et al., *Nat. Commun.* 2022.
//!
//! What is quantified: every TARGET row the glyco PIN emits (the top-1-per-scan
//! collapse, or the full dump under `--debug-glyco`), identified by its INTACT
//! precursor — peptide backbone + glycan composition + charge. The PIN also
//! carries decoy-backbone rows for Percolator; those are skipped here, since a
//! decoy precursor is not a glycopeptide and the decoy evidence this
//! quantification uses is its own exchangeable twin (+11 Th, extracted for
//! every target). A glycan topology is carried through from the search for
//! reporting, but it is NOT part of the target: two topologies of one
//! composition are the same mass and therefore one MS1 feature. That is a
//! property of the measurement, not a simplification.
//!
//! Where this deliberately differs from the peptide `--lfq` path, and why:
//!
//! * **The envelope is the real glycopeptide formula.** The peptide's
//!   elemental formula (from its residues and modification deltas) plus the
//!   glycan's residue composition (`Formula::glycan`). A glycan carries ~twice
//!   the oxygen and ~14 % fewer carbons per dalton than an averagine peptide,
//!   which moves the M+1/M+2 shares — the isotope-cosine gate reads exactly
//!   those shares, so an averagine envelope would score a true glycopeptide
//!   like a wrong one.
//! * **The envelope is longer.** A 3–4 kDa glycopeptide's monoisotopic peak is
//!   a minority of its envelope (M+2/M+3 can be taller). Isotopes are
//!   extracted until the extracted window holds [`ISOTOPE_COVERAGE`] of the
//!   theoretical distribution, between [`MIN_ISOTOPES`] and
//!   [`MAX_ISOTOPES`]. Extracting a fixed four, as for peptides, would read a
//!   truncated envelope and systematically distort the cosine and the mono
//!   share.
//! * **The intensity is the monoisotopic peak, not the envelope sum.** As in
//!   the published method, `Area` is the monoisotopic XIC integrated over the
//!   peak and divided by its theoretical share of the extracted envelope; the
//!   higher isotopes feed only the cosine. For a 3–4 kDa glycopeptide the
//!   monoisotope is 16–28 % of the envelope (a 3.5–6x multiplier), and the
//!   M+1..M+n traces are the ones most exposed to co-eluting interferences,
//!   which an envelope sum accumulates. The trade-off is that monoisotopic
//!   noise is amplified by the same factor. The raw observed mono area is
//!   still reported (`MonoArea`) next to the scaled `Area`.
//! * **The search's isotope offset is honoured.** When the search matched the
//!   precursor at M+i, the observed peak sits at `mono + i·ISOTOPE/z`; the
//!   extraction starts there and the theoretical envelope is rotated by i.
//!   Otherwise a true M+1 precursor would have its monoisotope extracted from
//!   the M−1 position and fail `min_mono_share`.
//!
//! Cross-run matching (`--glyco --mbr`) is deliberately simpler than the peptide
//! `--mbr`: a robust piecewise-linear RT map per run to the first input,
//! fitted on the glycopeptides both runs accepted (outlier-trimmed, one knot
//! per 300 anchors). The match itself is the published cross-run scorer
//! (`glyco_match`): the donor's isotope reference XIC is slid across the target
//! run over the alignment uncertainty (60–120 s) and scored by a small MLP, and
//! the chosen window — not the predicted apex — is what gets quantified; with
//! no fitted map, a transfer falls back to the mapped apex with a three-sigma
//! window clamped to 10–60 s. The transfers compete against their own decoy
//! twins, as a separate competition from the identified features, so a transfer
//! never borrows an identified feature's score distribution.
//!
//! What is NOT here (yet): FQR (a two-component GMM over a +10 Da decoy's
//! scores, as in the published method). The feature q-value below is a picked
//! target/decoy q-value against the same exchangeable twin the peptide LFQ
//! path uses
//! (`mono_mz + 11 Th`, no RT shift), not FQR. It gates **transfers** only:
//! a transfer's sole evidence is its MS1 feature, whereas an identification
//! is confirmed by its MS2 (the identification FDR is controlled upstream by
//! Percolator), so for identified features the q-value is reported but never
//! vetoes a feature that passed the cosine and mono-share gates.
//!
//! Outputs, next to `<base>.glyco.pin`:
//!
//! * `<base>.glyco.lfq.tsv` — one row per quantified (run, scan, glycopeptide,
//!   charge) feature, including the decoy evidence and the picked feature
//!   q-value. `Evidence` says `Identified` or `Transferred` (the latter has
//!   scan 0 and a `{donor SpecId}_mbr_{run}` SpecId). `SpecId` matches the PIN
//!   byte for byte, so a Percolator `q-value` column joins on it.
//! * `<base>.glyco.quant.tsv` — per run and glycopeptide (peptide + site +
//!   composition), summed over charge states. `NScans` counts the run's own MS2
//!   identifications; `NTransferred` counts the features transferred in by
//!   `--glyco --mbr` (which have no scan).
//! * `<base>.glyco.site.tsv` — per run, peptide and site: the glycoform shares
//!   (each composition's area over the site's total glyco area).
//! * `<base>.glyco.protein.tsv` — per run and protein: the summed glycopeptide
//!   area and Top3 (sum of the three largest glycopeptide areas).
//! * `<base>.glyco.matrix.tsv` — the same glycopeptide areas pivoted: one row
//!   per (peptide, site, composition), one column per run, empty where a run
//!   did not quantify it. This is the multi-fraction comparison view; the
//!   per-run tables above are its long form.
//! * `<base>.glyco.protein_matrix.tsv` — the protein pivot: per run, the
//!   protein's summed glycopeptide area and Top3.

use std::collections::{BTreeSet, HashMap};
use std::io::Write;
use std::path::{Path, PathBuf};

use andes_glyco::glycan_db::GlycanComp;
use andes_glyco::glycan_first::GlycanCore;
use model::mass::{ISOTOPE, PROTON};
use model::spectrum::Spectrum;
use model::Tolerance;
use output::glyco_pin::select_emitted_hits;
use quant::formula::Formula;
use quant::lfq::{competition_score, quantify, quantify_decoy, FeatureQuant, LfqParams, LfqTarget};
use quant::tdc::{picked_qvalues, Pair};
use rayon::prelude::*;
use search::candidate_gen::Candidate;
use search::glyco_search::GlycoSpectrumResult;
use search::search_index::SearchIndex;

use crate::cli::Cli;
use crate::glyco_match;
use crate::quant_run::RunMs1;

/// Fraction of the theoretical precursor envelope the extracted isotope
/// window must cover. The published method sums the whole theoretical
/// distribution; beyond ~95 % the remaining isotopes are below the MS1 noise
/// floor.
const ISOTOPE_COVERAGE: f64 = 0.95;

/// Isotopes extracted per feature, bounds. The floor matches the peptide
/// path's four; the ceiling keeps a very large glycopeptide's XIC count
/// bounded (24 isotope traces at 10 ppm is still cheap, but the marginal
/// signal past ~14 is noise).
const MIN_ISOTOPES: usize = 4;
const MAX_ISOTOPES: usize = 14;

/// Isotopes the theoretical envelope is computed over before trimming. An
/// M+3-selected precursor of a 5 kDa glycopeptide still needs its upper tail.
const ENVELOPE_ISOTOPES: usize = 28;

/// Default RT half-window around the identifying scan, s. Twice the peptide
/// `--lfq` default: an intact glycopeptide's XIC is wider, and a HILIC anchor
/// can sit a scan or two off the apex, so 60 s clips the peak on the late
/// side. `--lfq-rt-window` still overrides it.
const DEFAULT_RT_WINDOW_S: f64 = 120.0;

/// How a feature was measured: from its own identifying MS2 scan, or
/// transferred into a run that did not identify it (match-between-runs).
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

/// One (run, scan, glycopeptide, charge) target.
#[derive(Debug, Clone)]
struct Target {
    run_idx: usize,
    /// Run name (input file stem), repeated per row: the tables are per run.
    run: String,
    /// The PIN's `SpecId`, byte for byte, for an identified target;
    /// `{donor SpecId}_mbr_{run}` for a transfer, which no PIN row names.
    spec_id: String,
    /// 0 for a transfer: no MS2 scan in this run identified it.
    scan: i32,
    evidence: Evidence,
    /// Extraction half-window around `anchor_rt`, seconds. The identified
    /// window is `--lfq-rt-window`; a transfer's is set by the match scorer
    /// (its chosen peak window), or the alignment uncertainty when the
    /// scorer could not run.
    window_s: f64,
    /// Match search half-window (`--glyco --mbr` transfers), seconds: how far
    /// the donor's extraction window may slide from the RT prediction.
    search_s: f64,
    /// The identified row this transfer was built from (a `rows` index).
    donor: Option<usize>,
    /// `pre.SEQ_WITH_MODS.post` — the peptide identity across proteins and
    /// charge states.
    peptide: String,
    /// Peptide + glycans are grouped by this (peptide, composition) pair.
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
    /// Theoretical envelope from the observed peak, trimmed and normalised.
    envelope: Vec<f64>,
}

impl Target {
    fn lfq_target(&self) -> LfqTarget {
        LfqTarget {
            mono_mz: self.mono_mz,
            charge: self.charge,
            anchor_rt: self.anchor_rt,
            envelope: self.envelope.clone(),
            // 1: `FeatureQuant::area` stays the summed traced isotopes, as
            // before the whole-pattern scaling of the peptide path; the
            // reported `Area` is computed here from the monoisotopic trace.
            isotope_share: 1.0,
        }
    }
}

/// One target's measurement.
struct FeatureRow {
    target: Target,
    feature: Option<FeatureQuant>,
    decoy: Option<FeatureQuant>,
    score: Option<f32>,
    q: Option<f64>,
    /// The match score (transfers scored by the matcher; `score`
    /// is this value then).
    match_score: Option<f32>,
    /// The decoy twin's match score, the other half of that competition.
    decoy_match: Option<f32>,
}

/// Settings resolved from the `--lfq*` flags (the glyco path reuses them).
struct Settings {
    params: LfqParams,
    feature_fdr: f64,
    core: GlycanCore,
    debug: bool,
    /// `--glyco --mbr`: transfer accepted glycopeptides into runs that did not
    /// identify them.
    mbr: bool,
    /// `--glyco-mbr-match-train`: where to write the match-training table.
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
            feature_fdr: cli.lfq_feature_fdr,
            core: cli.glyco_core.into(),
            debug: cli.debug_glyco,
            mbr: cli.mbr,
            match_train: cli.glyco_mbr_match_train.clone(),
        }
    }
}

/// Quantify the emitted glyco PSMs and write the four tables.
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
    eprintln!(
        "Quantification: intact glycopeptides from MS1 (tol {}, rt window ±{} s, {MIN_ISOTOPES}-\
         {MAX_ISOTOPES} isotopes, min cosine {}, feature FDR {}, monoisotopic intensity)",
        match settings.params.tol {
            Tolerance::Ppm(v) => format!("{v} ppm"),
            Tolerance::Da(v) => format!("{v} Da"),
        },
        settings.params.rt_window_s,
        settings.params.min_cosine,
        settings.feature_fdr
    );
    let targets = build_targets(spectra, results, candidates, index, runs, &settings);
    if targets.is_empty() {
        eprintln!(
            "[glyco-quant] nothing to quantify: the PIN has no emitted hit with a resolved \
             glycan composition and a precursor retention time"
        );
        return Ok(());
    }
    let (features, quantified) = measure(&targets, runs, &settings, "identified")?;
    let mut rows = feature_rows(targets, quantified, &settings);
    assign_qvalues(&mut rows, &settings);

    // ── Match-between-runs (`--glyco --mbr`) ─────────────────────────────────
    // Transfers are built from the accepted identified features and quantified
    // in the runs that did not identify the glycopeptide, then scored as their
    // own target/decoy competition (a transfer has a predicted RT and a
    // narrow window, so its score distribution is not the identified one's).
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
            // The match scorer needs the donor's own XIC as well as the
            // target run's; references are read run by run, then the
            // transfers are measured.
            let refs = build_match_references(&rows, &transfers, runs, &settings)?;
            let (found, measured) = measure_transfers(&transfers, &refs, runs, &settings)?;
            transfer_counts.measured = found;
            let mut transfer_rows = feature_rows(transfers, measured, &settings);
            assign_qvalues(&mut transfer_rows, &settings);
            rows.append(&mut transfer_rows);
        }
        eprintln!(
            "[glyco-mbr] {} run maps ({} anchors); {} transfers attempted, {} measured, \
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
        eprintln!("[glyco-mbr] --mbr needs two or more runs; nothing to transfer");
    } else if runs.len() > 1 {
        eprintln!(
            "[glyco-mbr] match-between-runs is off; pass --mbr to transfer identified \
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
    eprintln!(
        "[glyco-quant] quantified {} features from {} targets; {} accepted ({} identified, \
         quality-gated only; {} transferred at feature q <= {}); {} glycopeptides, {} \
         glycoform-site spectra across {} runs",
        features + transfer_counts.measured,
        rows.len(),
        count_evidence(Evidence::Identified) + count_evidence(Evidence::Transferred),
        count_evidence(Evidence::Identified),
        count_evidence(Evidence::Transferred),
        settings.feature_fdr,
        kept.glycopeptides,
        kept.sites,
        runs.len()
    );
    Ok(())
}

/// Build one target per emitted glyco hit that carries a composition.
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
            // The glyco PIN is a Percolator input: it carries the decoy-backbone
            // rows the search scored. Those are not measurements of a
            // glycopeptide — their intact precursor is a random mass, so any
            // feature found there is noise that would enter the target side of
            // the q-value competition and the roll-up tables. The peptide LFQ
            // path also quantifies targets only. The decoy evidence this
            // quantification uses is its own exchangeable twin (the same
            // glycopeptide at `mono_mz + 11 Th`), extracted for every target.
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
            // Intact precursor: the peptide's exact neutral mass (the residue
            // masses the search scored) plus the glycan's residue mass. The
            // elemental formula is the same glycan composition; it is what the
            // isotope envelope is computed from.
            let intact_mass = cand.peptide.mass() + comp.mass;
            let theo_mz = intact_mass / z + PROTON;
            // The search matched this precursor at M+offset: the observed peak
            // is offset isotopes above the monoisotope. Glyco searches use the
            // 0..=2 window, so a negative offset cannot arise here; clamping it
            // keeps a user-widened window from extracting below the envelope.
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
            // The mass error of the MONOISOTOPIC hypothesis: the observed
            // precursor is offset isotopes above it, so remove that offset
            // before taking the difference (the PIN's `dm` does the same).
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
        eprintln!(
            "[glyco-quant] {no_rt} spectra have no retention time and are not quantified \
             (their XIC anchor would be undefined)"
        );
    }
    if no_glycan > 0 {
        eprintln!(
            "[glyco-quant] {no_glycan} emitted hits carry no enumerated glycan composition \
             and are not quantified (their intact precursor mass is unknown)"
        );
    }
    if decoy_rows > 0 {
        eprintln!(
            "[glyco-quant] {decoy_rows} emitted rows are decoy-backbone hits and are not \
             quantified (the PIN carries them for Percolator; a decoy precursor is not a \
             glycopeptide)"
        );
    }
    out
}

/// The peptide's elemental formula, from its residues and modification deltas
/// (the same construction as the peptide `--lfq` path).
fn peptide_formula(cand: &Candidate) -> Formula {
    Formula::peptide(
        cand.peptide
            .residues
            .iter()
            .map(|aa| (aa.residue, aa.mod_.as_ref().map_or(0.0, |m| m.mass_delta))),
    )
}

/// The theoretical envelope from isotope `offset` onward, trimmed to
/// [`ISOTOPE_COVERAGE`] of the distribution and normalised to sum 1.
///
/// The trimmed total is the coverage denominator: `envelope(n)` normalises
/// over the `n` computed isotopes, so the tail past `ENVELOPE_ISOTOPES` is
/// ignored — it is below 10⁻⁵ of the distribution at these masses.
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

/// One target's measurement: its feature and its decoy twin, both optional
/// (no peak in the window).
#[derive(Debug, Clone, Default)]
struct Measurement {
    feature: Option<FeatureQuant>,
    decoy: Option<FeatureQuant>,
    /// The match score of the transferred window. Set only when both
    /// this row's and its decoy's match scores exist, so the transfer's
    /// target/decoy competition stays symmetric.
    match_score: Option<f32>,
    /// The decoy twin's match score (set together with `match_score`).
    decoy_match: Option<f32>,
}

/// Quantify every target in its run. One run's MS1 scans are resident at a
/// time (as in the peptide LFQ path); the returned vector is indexed like
/// `targets`.
fn measure(
    targets: &[Target],
    runs: &[RunMs1<'_>],
    settings: &Settings,
    phase: &str,
) -> Result<(usize, Vec<Measurement>), String> {
    measure_with(targets, runs, settings, phase, measure_identified)
}

/// Quantify one identified target in its own run: the envelope over the search
/// RT, mono-scaled to the monoisotopic intensity, plus the +11 Th decoy twin.
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

/// Quantify every target in its run with `f`, one run's MS1 scans resident at
/// a time (as in the peptide LFQ path); the returned vector is indexed like
/// `targets`.
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
        eprintln!(
            "[glyco-quant] {phase} {}: {} features in {} targets, {:.2}s",
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

/// One donor feature's match references in its own run: the identified
/// feature's isotope XIC and the same XIC at the decoy twin's m/z.
struct MatchReferences {
    target: glyco_match::Reference,
    decoy: Option<glyco_match::Reference>,
}

/// Build the match references of every donor a transfer names, grouped by run
/// so one run's MS1 scans are resident at a time. References are keyed by the
/// donor's `rows` index; one donor may serve several target runs.
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

/// Measure the transfers with the match scorer: slide the donor's
/// extraction window across the target run, let the model choose the window
/// and score it, then quantify that window with the ordinary extractor. The
/// decoy twin gets its own reference and its own scored window, so the
/// transfer q-value stays a symmetric competition.
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
        // The competition is target against decoy: use the match score only
        // when both sides have one.
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

/// The published label-free intensity: the extracted isotope at the search's
/// offset divided by its theoretical share of the extracted envelope
/// (`Σ_i theo_i / theo_offset`). `envelope` is normalised over the extracted
/// isotopes, so its first element is that share; the higher isotopes enter
/// only through the envelope cosine, so interference in an M+1..M+n trace
/// does not move the reported quantity.
fn mono_area(f: &FeatureQuant, envelope: &[f64]) -> f64 {
    match (f.isotope_areas.first(), envelope.first()) {
        (Some(&mono), Some(&share)) if share > 0.0 => mono / share,
        _ => 0.0,
    }
}

/// The `--lfq` parameters of one row: its own extraction window (a transfer's
/// is the alignment uncertainty, an identified feature's `--lfq-rt-window`).
fn params_of(row: &FeatureRow, settings: &Settings) -> LfqParams {
    LfqParams {
        rt_window_s: row.target.window_s,
        ..settings.params.clone()
    }
}

/// Rows with their feature scores, before any q-value. The decoy twin is
/// extracted with the same procedure at a shifted m/z; the competition follows
/// the peptide LFQ path.
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
            // A matched transfer is scored by the model, not by the cosine
            // gate: eligibility is the model's (and `quantify`'s own)
            // judgement, applied identically to the decoy twin.
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

/// Picked target/decoy q-values over `rows`, which must be one evidence class
/// (identified features and transfers are separate competitions: their windows
/// and their search spaces differ). The height scale is the class's own.
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

/// Whether a row's feature enters the quant tables.
///
/// A transfer has no MS2 of its own, so its only FDR control is the picked
/// feature q against the +11 Th decoy twin. Its competition score is the
/// match score of the window the model chose (falling back to the
/// isotope-envelope competition when either side has no match). An
/// identification does not need that gate: the peptide, composition and site
/// are already confirmed by the MS2 evidence (Percolator controls their FDR),
/// and the quantification quality gates that do matter are the envelope cosine
/// (`score`, None when below `--lfq-min-cosine`) and `min_mono_share`, which
/// reject a feature extracted on the wrong isotope. The published method
/// likewise never re-filters an identification on MS1 decoy evidence; the
/// q-value is still computed and reported in `FeatureQValue` for both classes.
fn accepted(row: &FeatureRow, settings: &Settings) -> bool {
    if row.feature.is_none() {
        return false;
    }
    match row.target.evidence {
        Evidence::Identified => row.score.is_some(),
        Evidence::Transferred => row.q.is_some_and(|q| q <= settings.feature_fdr),
    }
}

// ── Match-between-runs ──────────────────────────────────────────────────────

/// Fewest shared glycopeptides (same peptide, composition and charge,
/// quantified in a run and in the reference run) that make a cross-run RT map.
const MBR_MIN_ANCHORS: usize = 50;
/// Largest RT-map residual (robust σ) a transfer may be predicted from, s.
/// Above this the alignment is not informative enough to place a predicted
/// apex (the transfer would carry no RT evidence the identified window does
/// not already carry), so the run is not used.
const MBR_MAX_RESIDUAL_S: f64 = 45.0;
/// Anchors per RT knot: the map is piecewise-linear between equal-count bins
/// of [`MBR_ANCHORS_PER_KNOT`] anchors, which follows a non-linear gradient
/// warp one straight line cannot (measured on two mouse-liver runs: 36 s
/// residual for one line, 28 s piecewise).
const MBR_ANCHORS_PER_KNOT: usize = 300;
const MBR_MIN_KNOTS: usize = 3;
const MBR_MAX_KNOTS: usize = 10;
/// Transfer search half-window bounds, s: never tighter than the minimum (an
/// apex can sit a scan or two from the prediction) nor wider than the
/// identified window.
const MBR_MIN_WINDOW_S: f64 = 10.0;
const MBR_MAX_WINDOW_S: f64 = 60.0;
/// Bounds of the match search half-window, s: how far the donor's whole
/// extraction window may slide from its RT prediction. The original searches
/// ±2 min; the RT map is a better prior, so andes searches 3σ within these.
const MBR_MATCH_MIN_SEARCH_S: f64 = 60.0;
const MBR_MATCH_MAX_SEARCH_S: f64 = 120.0;

/// What the transfer stage did, for the log.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
struct TransferCounts {
    /// Runs with a usable RT map to the reference (run 0 always has one).
    maps: usize,
    /// Shared glycopeptides the maps were fitted on, summed over runs.
    anchors: usize,
    /// (glycopeptide, run) transfers attempted.
    candidates: usize,
    /// Of those, how many had a feature in the target run.
    measured: usize,
}

/// A monotone piecewise-linear RT map between one run and the reference run
/// (run 0), defined by knots `(xs, ys)` that increase strictly in both
/// coordinates: `rt_ref = f(rt_run)`. Between knots the map is linear; past the
/// outermost anchors it extrapolates on the nearest segment's slope.
#[derive(Debug, Clone, PartialEq)]
struct RtMap {
    xs: Vec<f64>,
    ys: Vec<f64>,
    /// Robust residual σ of the map, seconds on the reference run's scale.
    sigma_s: f64,
}

impl RtMap {
    /// The identity map of the reference run itself.
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

    /// The run's RT for a reference-run RT; `None` if the knots are not
    /// strictly increasing (never built that way, but inversion would divide
    /// by zero).
    fn inverse(&self, rt_ref: f64) -> Option<f64> {
        let ys_are_increasing = self.ys.windows(2).all(|w| w[1] > w[0]);
        ys_are_increasing.then(|| interpolate(&self.ys, &self.xs, rt_ref))
    }

    /// d(run RT)/d(reference RT) at `rt_ref`: converts a residual on the
    /// reference scale onto the run's own scale.
    fn inverse_scale(&self, rt_ref: f64) -> f64 {
        let i = self
            .ys
            .partition_point(|&k| k <= rt_ref)
            .saturating_sub(1)
            .min(self.ys.len() - 2);
        (self.xs[i + 1] - self.xs[i]) / (self.ys[i + 1] - self.ys[i])
    }
}

/// Piecewise-linear `v -> w` through knots `(a, b)`, extrapolating on the
/// first/last segment.
fn interpolate(a: &[f64], b: &[f64], v: f64) -> f64 {
    debug_assert!(a.len() >= 2 && a.len() == b.len());
    let i = a
        .partition_point(|&k| k <= v)
        .saturating_sub(1)
        .min(a.len() - 2);
    let slope = (b[i + 1] - b[i]) / (a[i + 1] - a[i]);
    b[i] + slope * (v - a[i])
}

/// Least-squares line `y = a + b·x` of `pairs`, `None` for fewer than three
/// points, a degenerate x spread, or a non-increasing slope (an RT map must be
/// monotone).
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

/// A robust RT map from `(rt_run, rt_ref)` anchors. A least-squares line
/// trims obvious mis-pairs (beyond 3 MAD of its residuals); the kept anchors
/// are binned into equal counts and each bin's median becomes a knot, so the
/// map follows a non-linear gradient warp. `sigma_s` is the MAD-scaled
/// residual of the kept anchors around the finished map.
fn fit_rt_map(pairs: &[(f64, f64)]) -> Option<RtMap> {
    let (intercept, slope) = fit_line(pairs)?;
    let residuals: Vec<f64> = pairs
        .iter()
        .map(|&(x, y)| y - (intercept + slope * x))
        .collect();
    let center = crate::quant_mbr::median(residuals.clone());
    let mad =
        1.4826 * crate::quant_mbr::median(residuals.iter().map(|r| (r - center).abs()).collect());
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
        // A knot that does not advance would make the map non-monotone; drop
        // it (the neighbouring segments already cover that RT span).
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
    let mut dev: Vec<f64> = kept.iter().map(|&(x, y)| y - map.forward(x)).collect();
    let dev_center = crate::quant_mbr::median(dev.clone());
    let sigma =
        1.4826 * crate::quant_mbr::median(dev.drain(..).map(|d| (d - dev_center).abs()).collect());
    let map = RtMap {
        sigma_s: sigma,
        ..map
    };
    map.sigma_s.is_finite().then_some(map)
}

/// A glycopeptide's measurement identity across runs: peptide, composition,
/// charge (the same intact precursor).
type GlycoKey = (String, String, u8);

fn glyco_key(t: &Target) -> GlycoKey {
    (t.peptide.clone(), t.glycan.clone(), t.charge)
}

/// Build the transfers: for every accepted identified feature, every other
/// run with an RT map that did NOT identify its glycopeptide gets a target at
/// the predicted RT, searched in a window set by the alignment uncertainty.
///
/// A key's best-scoring accepted donor is used, so one (glycopeptide, run)
/// pair is measured once even when several runs identified the glycopeptide.
/// A run that identified the glycopeptide but rejected the feature (q above
/// the threshold) is not a transfer target: resurrecting a rejected
/// identification as a transfer would bypass the rejection.
/// The cross-run alignment `--glyco --mbr` works from: every glycopeptide each
/// run identified, its best accepted donor, and the run's RT map to the first
/// input. Shared by the transfer builder and the match-model training dump.
struct Alignment {
    /// Every key each run identified, accepted or not.
    identified: Vec<std::collections::HashSet<GlycoKey>>,
    /// The best accepted donor per key and run (a `rows` index).
    donor: Vec<HashMap<GlycoKey, usize>>,
    maps: Vec<Option<RtMap>>,
    counts: TransferCounts,
}

impl Alignment {
    fn build(rows: &[FeatureRow], run_names: &[&str], settings: &Settings) -> Self {
        let mut counts = TransferCounts {
            maps: 1, // the reference run is its own map
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

        // RT maps to the reference run.
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
                eprintln!(
                    "[glyco-mbr] {}: {} shared glycopeptides, fewer than {MBR_MIN_ANCHORS} \
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
                eprintln!(
                    "[glyco-mbr] {}: RT map residual {:.1}s is above {MBR_MAX_RESIDUAL_S}s; \
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

    /// The RT of `rt` in run `from_idx` expressed in run `to_idx`'s scale, and
    /// the combined alignment σ of the prediction, seconds.
    fn predict(&self, from_idx: usize, to_idx: usize, rt: f64) -> Option<(f64, f64)> {
        let from = self.maps[from_idx].as_ref()?;
        let to = self.maps[to_idx].as_ref()?;
        let t_ref = from.forward(rt);
        let predicted = to.inverse(t_ref)?;
        // Both ends contribute: the donor run's map error at the donor RT and
        // the target run's inverse map error at the prediction (the residual is
        // on the reference scale; the local inverse slope converts it to the
        // target run's scale).
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
    // One transfer per (key, target run), from the best-scoring donor. The
    // donor's score orders two runs that both identified the key.
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
            // The matcher slides the donor's whole extraction window across
            // the alignment uncertainty; 60 s is the floor (the map's residual
            // is a robust σ, so 3σ can be small on a clean gradient).
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
    // HashMap order is not deterministic; sort so the tables are.
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

/// Write the training table of the match scorer (`--glyco-mbr-match-train`).
///
/// For every ordered pair of runs that both accepted a glycopeptide, score the
/// donor's reference against the target run exactly as `--glyco --mbr` does, and
/// label each scored window: positive when it covers the peak the target run
/// identified on its own, negative otherwise. Every window of the decoy twin
/// (+11 Th) is negative. Positives therefore come from runs where the answer is
/// known, which is what lets the model learn the target-versus-decoy boundary
/// the transfer q-value needs.
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
        // The donor run's references for every glycopeptide it accepted.
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
            // This run's own accepted feature for a key, when it has one.
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
                    // The prediction put the true peak outside every scored
                    // window; the pair carries no positive and is skipped.
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
            eprintln!(
                "[glyco-mbr-match-train] runs {d}->{t_run}: {pair_pos} positive, \
                 {pair_neg} negative candidate windows"
            );
        }
    }
    out.flush()
        .map_err(|e| format!("{}: {e}", path.display()))?;
    eprintln!(
        "[glyco-mbr-match-train] {} glycopeptide pairs, {positives} positive and \
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

/// The path of an output table: `<base>.glyco.<kind>.tsv` beside the PIN.
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

/// Counts reported after the roll-up.
struct Kept {
    glycopeptides: usize,
    sites: usize,
}

/// The roll-up tables: per run and glycopeptide, per run and site, per run and
/// protein, and the two cross-run matrices (`run_names` is the column order).
/// Only accepted features enter them: a feature the exchangeable decoy
/// explains is not a measurement of that glycopeptide.
fn write_rollup_tables(
    pin: &Path,
    rows: &[FeatureRow],
    run_names: &[&str],
    settings: &Settings,
) -> Result<Kept, String> {
    // One feature per (run, peptidoform, composition, charge): the same
    // peptide+charge elutes once, so PSMs from different scans of one elution
    // describe the same XIC. Keep the best-scoring one; the others are
    // duplicates of one measurement, not replicates.
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

    // Glycopeptide rows: sum over charge states.
    struct Glyco {
        area: f64,
        best_cosine: f32,
        best_q: f64,
        charges: BTreeSet<u8>,
        scans: BTreeSet<i32>,
        /// Features transferred into this run (no MS2 of their own).
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
            // A transfer has no identifying scan in this run (scan 0 is a
            // placeholder); NScans counts MS2 evidence only.
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

    // Site rows: each glycoform's share of its site's total glyco signal.
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

    // Cross-run glycopeptide matrix: one row per (peptide, composition), one
    // area column per run in `run_names` order, empty where that run did not
    // quantify the glycopeptide. This is the view a multi-fraction experiment
    // compares across; the per-run tables above are its long form.
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

    // Protein rows: the summed glycopeptide area and Top3.
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
    // Cross-run protein matrix, the same pivot: per run, the protein's summed
    // glycopeptide area and Top3 (empty where the run quantified none).
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

/// `HexNAc2Hex8Fuc0NeuAc0NeuGc0` — the composition half of the PIN's glycan
/// tag, identical to `output::glyco_pin`'s rendering.
fn composition_label(comp: &GlycanComp) -> String {
    format!(
        "HexNAc{}Hex{}Fuc{}NeuAc{}NeuGc{}",
        comp.hexnac, comp.hex, comp.fuc, comp.neuac, comp.neugc
    )
}

/// The glycosite the PIN reports for this backbone: the sequon position when
/// the sequence determines it, `N?` / `S/T?` when it does not (the same rule
/// as `output::glyco_pin::glycan_tag_with_site`).
fn site_label(residues: &[u8], post: u8, core: GlycanCore) -> String {
    match core {
        GlycanCore::NGlycan => match andes_glyco::sequon::glycosites(residues, post).as_slice() {
            [only] => format!("N{}", only + 1),
            _ => "N?".to_string(),
        },
        GlycanCore::OGlycan => match andes_glyco::sequon::all_ser_thr_sites(residues).as_slice() {
            [only] => format!("{}{}", residues[*only] as char, only + 1),
            _ => "S/T?".to_string(),
        },
    }
}

/// Every protein accession this PSM's candidate list resolves to, `;`-joined.
/// The accession rule matches `output::row_context::resolve_accession` (the
/// PIN's `Proteins` columns), including its `PROT_<index>` fallback for a
/// candidate whose protein record is absent from the index.
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

/// A float cell: six decimals, or empty for a value with no measurement
/// (NaN). Empty keeps the column numeric in pandas rather than mixing in the
/// string `NaN`.
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
        // The M+2 extraction starts at the theoretical M+2 share, so its
        // first element is larger relative to the rest than the M+0 window's.
        assert!(m2[0] > mono[0], "{mono:?} {m2:?}");
        assert!((m2.iter().sum::<f64>() - 1.0).abs() < 1e-9);
        // And the M+2 window keeps the isotope count bounded.
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
            feature_fdr: 0.05,
            core: GlycanCore::NGlycan,
            debug: false,
            mbr: false,
            match_train: None,
        }
    }

    /// A glycopeptide target with a real envelope, at `mono_mz`, anchored at
    /// `anchor_rt`.
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

    /// One run with a single glycopeptide eluting as a Gaussian at 100 s
    /// (sigma 4 s, 1 Hz scans), its isotope pattern from the target's own
    /// theoretical envelope.
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
        // Mono-scaled area = mono area / mono share. With the observed
        // isotope shares equal to theory, that is the Gaussian area of the
        // whole envelope (height · sigma · √(2π), envelope normalised to 1).
        let expected = 1e6 * 4.0 * (2.0 * std::f64::consts::PI).sqrt();
        assert!((f.area / expected - 1.0).abs() < 0.03, "{f:?}");
        // The shifted decoy twin finds nothing in this clean run.
        assert!(m.decoy.is_none(), "{:?}", m.decoy);
    }

    #[test]
    fn an_isotope_offset_extracts_the_observed_peak() {
        let mut t = target(800.0, 2, 97.0);
        // The search matched M+1: the feature sits one isotope above the
        // monoisotope, and the envelope starts there.
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

    /// A run with only the monoisotopic trace: the higher isotope channels,
    /// which interference and low S/N mostly affect, are absent.
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
        // The point of the mono estimator: losing the M+1..M+n channels
        // changes the cosine, not the quantity. With the envelope sum it
        // would lose the fraction those channels carried.
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
        b.target.scan = 2; // same peptide+charge: a duplicate XIC, dropped
        let mut c = row(target(800.0, 3, 97.0), 50.0, 0.8, 0.002);
        c.target.scan = 3;
        c.target.charge = 3;
        // A transfer is the class the feature q still gates.
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
        // q = 0.5 is above the feature FDR; the MS2 identification is the
        // evidence, so the feature still enters the tables (the method never
        // re-filters an identification on MS1 decoy evidence).
        let identified = row(target(800.0, 2, 97.0), 100.0, 0.9, 0.5);
        // A transfer with the same q is rejected: its feature is its only
        // evidence, so the q gates it.
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
        // A second glycoform, quantified only in r2.
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
        // Sorted by peptide, then composition.
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

    // ── Match-between-runs ─────────────────────────────────────────────────

    #[test]
    fn the_rt_map_recovers_a_linear_shift() {
        // rt_ref = 4 + 1.1 · rt_run, exactly.
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
        // The inverse prediction round-trips.
        let rt_ref = map.forward(123.0);
        assert!((map.inverse(rt_ref).unwrap() - 123.0).abs() < 1e-6);
    }

    #[test]
    fn the_rt_map_follows_a_non_linear_warp() {
        // A gradient drift that no straight line follows: ±40 s in the
        // middle of the run, back to 0 at both ends.
        let pairs: Vec<(f64, f64)> = (0..900)
            .map(|i| {
                let x = i as f64 * 4.0;
                let phase = x / 3600.0 * std::f64::consts::PI;
                (x, x + 40.0 * (phase).sin())
            })
            .collect();
        // A straight line through these anchors would leave a ~25 s residual
        // (the sine peak), which is above the 3 MAD trim width.
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
        // Five mis-identified anchors: without trimming they would move the
        // map away from every true anchor.
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

    /// An accepted identified row for `key` in `run_idx` with its apex at
    /// `apex_rt`.
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
            // run 2 is 0.8x the reference scale plus 50 s.
            rows.push(identified_row(k, 1, "r2", 0.8 * rt + 50.0, 0.01, 10.0));
        }
        // Identified only in run 1: the one transfer.
        rows.push(identified_row(n, 0, "r1", 500.0, 0.01, 20.0));
        // Identified in run 2 but its feature failed the cosine gate (score
        // None): not an anchor donor and not a transfer target. The decoy q
        // alone would no longer reject an identified feature.
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
        // The prediction is the inverse map of the donor's apex: run 2 elutes
        // at 0.8 · rt_ref + 50 s.
        let expected = 0.8 * 500.0 + 50.0;
        assert!(
            (t.anchor_rt - expected).abs() < 1e-6,
            "{} vs {expected}",
            t.anchor_rt
        );
        // Clean anchors leave no residual, so the window is the floor.
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
