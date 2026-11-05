//! Match-between-runs for isobaric (TMT) reporter ions: `andes search --tmt
//! --mbr`. High-resolution MS2 reporters only.
//!
//! A donor is a precursor identified at the donor q-value in at least one
//! run (the library of `--lfq --mbr`, built by the same code). In a run that
//! lacks it, the candidates are the MS2 scans that carry no identification
//! at the donor q-value, whose selected precursor m/z lies within the
//! learned ppm window of the donor's theoretical m/z, with the donor's
//! charge, inside the aligned retention-time window. The candidate whose
//! fragments match the donor peptide's library spectrum best (spectral
//! angle) is the transfer; a scan claimed by several donors goes to the best
//! match. The transferred quantity is the candidate's impurity-corrected
//! reporter vector.
//!
//! Window calibration: the m/z window and the retention-time multiplier come
//! from the run's own identified pairs. Risk: three small logistic models
//! (see [`quant::tmt_transfer`]) are fitted on the same data without the
//! experimental design: the identified noise floor (an identified scan
//! against another identified scan of its precursor), the value risk (a free
//! sibling scan of an identified precursor against the identified scan, or an
//! empty-channel level above the identified P99) and the identity risk
//! (siblings against mass-shifted decoy scans). `quant_q` of a transfer is the
//! running mean of the predicted excess error over the transfers at least as
//! reliable.

use std::collections::HashMap;
use std::path::Path;

use model::Tolerance;
use output::library_table::StoredLibrary;
use output::ms2_table::read_ms2_table;
use output::run_quant_table::{open_run_quant_table, read_run_quant_mods};
use output::TransferInfo;
use quant::decoys::{probes, Hypothesis};
use quant::isobaric::{CorrectionMatrix, Plex};
use quant::ms2::{
    fragment_ladder, residue_masses, score_scan, Fragment, LibrarySpectrum, Ms2RunIndex, ScanMatch,
};
use quant::tmt_transfer::{
    calibration, fragment_evidence, learn_rt_multiplier, log2_deviation, quant_q, scan_features,
    PpmWindow, RunProfile, Training, TransferRisk, N_SCAN, N_XFER,
};
use rayon::prelude::*;

use crate::quant_cmd::{
    artifact, build_and_write_library, library_runs, scan_of, split_key, ScoreLoader, ScoreSet,
    MS2_SUFFIX, TABLE_SUFFIX,
};

/// Index of the donor-run loading cosine within the transfer covariates.
const BATCH_COS_AT: usize = 6;

/// Settings of the transfer.
pub(crate) struct Params<'a> {
    pub plex: Plex,
    pub correction: Option<&'a CorrectionMatrix>,
    /// Library entry gate and donor PSM thresholds (as `--lfq --mbr`).
    pub library_fdr: f64,
    pub donor_q: f64,
    pub donor_pep: Option<f64>,
    /// Transfers with `quant_q` at or below this are accepted.
    pub accept_q: f64,
    /// RMS log2 profile deviation that defines a wrong value.
    pub wrong_log2: f64,
}

/// One transferred reporter row.
#[derive(Debug, Clone)]
pub(crate) struct TransferRow {
    pub run: String,
    pub native_id: String,
    pub spec_key: String,
    pub scan: i32,
    pub rt: f64,
    pub precursor_mz: f64,
    pub peptidoform: String,
    pub charge: u8,
    pub calc_mz: f64,
    pub proteins: Vec<String>,
    pub raw: Vec<f32>,
    pub corrected: Vec<f32>,
    pub purity: Option<f32>,
    pub info: TransferInfo,
    pub accepted: bool,
    scan_features: [f64; N_SCAN],
    xfer: [f64; N_XFER],
    donors: Vec<u16>,
    run_index: u16,
    donor_pep: f64,
}

/// What a transfer run did.
#[derive(Debug, Default)]
pub(crate) struct Summary {
    pub candidates: usize,
    pub accepted: usize,
    pub model_fitted: bool,
}

/// A training row awaiting the cross-run covariate (donor-run loading cosine).
struct TrainRow {
    x: Vec<f64>,
    y: f64,
    run: u16,
    donors: Vec<u16>,
}

#[derive(Default)]
struct Counts {
    siblings: usize,
    targets_real: usize,
    decoys_real: usize,
}

struct RunOut {
    profile: RunProfile,
    floor: Vec<TrainRow>,
    value: Vec<TrainRow>,
    identity: Vec<TrainRow>,
    cands: Vec<TransferRow>,
    counts: Counts,
}

/// An identified PSM of a run joined to its rescoring result.
struct IdRow {
    key: String,
    peptidoform: String,
    charge: u8,
    mz: f64,
    rt: f64,
    spec_id: String,
    score: f64,
}

fn read_ident(
    dir: &Path,
    stem: &str,
    scores: &ScoreSet,
    p: &Params<'_>,
) -> Result<Vec<IdRow>, String> {
    let path = artifact(dir, stem, TABLE_SUFFIX);
    let mut out = Vec::new();
    let columns = [
        "spec_id",
        "rt_s",
        "precursor_key",
        "peptidoform",
        "charge",
        "mz",
        "is_decoy",
    ];
    for batch in open_run_quant_table(&path, &columns)? {
        for row in batch? {
            if row.is_decoy {
                continue;
            }
            let Some(s) = scores.get(&row.spec_id) else {
                continue;
            };
            if s.is_decoy() || s.q > p.donor_q || p.donor_pep.is_some_and(|c| s.pep > c) {
                continue;
            }
            out.push(IdRow {
                key: row.precursor_key,
                peptidoform: row.peptidoform,
                charge: row.charge,
                mz: row.mz,
                rt: row.rt_s,
                spec_id: row.spec_id,
                score: s.score,
            });
        }
    }
    Ok(out)
}

fn corrected_of(p: &Params<'_>, raw: &[f32]) -> Vec<f32> {
    match p.correction {
        Some(m) => m.correct(raw),
        None => raw.to_vec(),
    }
}

/// The scan number of a native id (`... scan=1234`), 0 when it has none.
fn scan_number(native_id: &str) -> i32 {
    native_id
        .rsplit_once("scan=")
        .map(|(_, tail)| {
            tail.chars()
                .take_while(char::is_ascii_digit)
                .collect::<String>()
        })
        .and_then(|d| d.parse().ok())
        .unwrap_or(0)
}

struct Found {
    pos: usize,
    m: ScanMatch,
    rt_z: f64,
    ppm: f64,
}

/// The free scan best matching `spec` among those whose retention time is in
/// `window`, whose selected m/z is within `ppm` of `mz` and whose charge is
/// `charge` (or unrecorded), at the entry's FAIMS voltage.
#[allow(clippy::too_many_arguments)]
fn best_candidate(
    index: &Ms2RunIndex,
    spec: &LibrarySpectrum,
    mz: f64,
    charge: u8,
    faims: Option<f32>,
    window: (f64, f64),
    t_hat: f64,
    s: f64,
    ppm: &PpmWindow,
    free: &[bool],
    tol: Tolerance,
) -> Option<Found> {
    let scans = index.scans();
    let start = scans.partition_point(|x| x.rt < window.0);
    let mut best: Option<Found> = None;
    for (pos, sc) in scans.iter().enumerate().skip(start) {
        if sc.rt > window.1 {
            break;
        }
        if !free[pos] || (sc.charge != 0 && sc.charge != charge) {
            continue;
        }
        if faims.is_some_and(|cv| sc.faims_cv.is_some_and(|v| v != cv)) {
            continue;
        }
        let err = (sc.precursor_mz - mz) / mz * 1e6;
        if !ppm.contains(err) {
            continue;
        }
        let m = score_scan(spec, &sc.peaks, tol);
        let rt_z = (sc.rt - t_hat) / s;
        let better = match &best {
            None => true,
            Some(b) => {
                m.angle > b.m.angle
                    || (m.angle == b.m.angle
                        && (m.n_matched > b.m.n_matched
                            || (m.n_matched == b.m.n_matched && rt_z.abs() < b.rt_z.abs())))
            }
        };
        if better {
            best = Some(Found {
                pos,
                m,
                rt_z,
                ppm: err,
            });
        }
    }
    best
}

struct Ctx<'a> {
    dir: &'a Path,
    stems: &'a [String],
    stored: &'a StoredLibrary,
    scores: &'a ScoreSet,
    p: &'a Params<'a>,
    mods: HashMap<String, f64>,
    /// Precursor key -> (score, run) of its best identified PSM.
    best_run: HashMap<String, (f64, u16)>,
    /// Precursor key -> entry position in the library.
    key_entry: HashMap<&'a str, usize>,
    min_mz: f64,
}

impl Ctx<'_> {
    fn ladder(&self, peptidoform: &str, charge: u8) -> Vec<Fragment> {
        match residue_masses(peptidoform, &self.mods) {
            Ok(m) => fragment_ladder(&m, charge.max(1)),
            Err(_) => Vec::new(),
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn covariates(
    ctx: &Ctx<'_>,
    scan: &quant::ms2::Ms2Scan,
    corrected: &[f32],
    profile: &RunProfile,
    ladder: &[Fragment],
    tol: Tolerance,
    charge: u8,
) -> [f64; N_SCAN] {
    let (current, top10) = fragment_evidence(ladder, &scan.peaks, tol, ctx.min_mz);
    let total: f64 = scan
        .peaks
        .iter()
        .filter(|p| p.0 as f64 >= ctx.min_mz)
        .map(|p| f64::from(p.1))
        .sum();
    let rep = scan.reporters.as_ref();
    scan_features(
        corrected,
        rep.and_then(|r| r.purity),
        rep.and_then(|r| r.sn_sum),
        profile,
        current,
        top10,
        charge,
        total,
    )
}

fn transfer_features(
    f: &Found,
    ppm: &PpmWindow,
    donor_q: f64,
    n_donors_other: usize,
) -> [f64; N_XFER] {
    let mut x = [0.0; N_XFER];
    x[0] = f64::from(f.m.angle);
    x[1] = f64::from(f.m.explained);
    x[2] = f64::from(f.m.n_matched);
    x[3] = f.rt_z.abs();
    x[4] = ppm.z(f.ppm);
    x[5] = donor_q.max(1e-6).log10();
    x[BATCH_COS_AT] = f64::NAN;
    x[7] = (1.0 + n_donors_other as f64).log2();
    x
}

fn process_run(ctx: &Ctx<'_>, r: usize) -> Result<RunOut, String> {
    let stem = &ctx.stems[r];
    let r16 = r as u16;
    let lib = &ctx.stored.library;
    let ms2_path = artifact(ctx.dir, stem, MS2_SUFFIX);
    let table = read_ms2_table(&ms2_path)?;
    if table.run != *stem {
        return Err(format!(
            "{}: the table is of run {:?}, not {stem:?}",
            ms2_path.display(),
            table.run
        ));
    }
    let tol = table.fragment_tol;
    let index = Ms2RunIndex::new(table.scans);
    let scans = index.scans();
    if scans.iter().any(|s| s.reporters.is_none()) {
        return Err(format!(
            "{}: the MS2 table holds no reporter ions; --tmt --mbr needs a table written by `andes search --tmt --mbr`",
            ms2_path.display()
        ));
    }
    let by_key = index.key_positions();
    let ids = read_ident(ctx.dir, stem, ctx.scores, ctx.p)?;
    let id_pos: Vec<Option<usize>> = ids.iter().map(|i| scan_of(&by_key, &i.spec_id)).collect();
    let mut free = vec![true; scans.len()];
    for &p in id_pos.iter().flatten() {
        free[p] = false;
    }
    let corrected: Vec<Vec<f32>> = scans
        .par_iter()
        .map(|s| corrected_of(ctx.p, s.reporters.as_ref().map_or(&[][..], |r| &r.raw)))
        .collect();
    let mut seen: Vec<usize> = id_pos.iter().flatten().copied().collect();
    seen.sort_unstable();
    seen.dedup();
    let ident_vectors: Vec<Vec<f32>> = seen.iter().map(|&p| corrected[p].clone()).collect();
    let profile = RunProfile::from_identified(&ident_vectors);
    drop(ident_vectors);

    // Window calibration on the identified pairs.
    let ppm_errors: Vec<f64> = ids
        .iter()
        .zip(&id_pos)
        .filter_map(|(i, p)| p.map(|p| (scans[p].precursor_mz - i.mz) / i.mz * 1e6))
        .collect();
    let ppm = PpmWindow::learn(&ppm_errors).unwrap_or_else(PpmWindow::fallback);
    let mut z = Vec::new();
    for (i, p) in ids.iter().zip(&id_pos) {
        if p.is_none() {
            continue;
        }
        let Some(&e) = ctx.key_entry.get(i.key.as_str()) else {
            continue;
        };
        let pid = lib.entries[e].precursor;
        if let Some((t, s)) = lib.alignment.predict_excluding(pid, r16) {
            if s > 0.0 && s.is_finite() {
                z.push((i.rt - t) / s);
            }
        }
    }
    let (k, k_learned) = learn_rt_multiplier(&z);

    let mut by_ident: HashMap<&str, Vec<usize>> = HashMap::new();
    for (i, p) in ids.iter().zip(&id_pos) {
        if let Some(p) = p {
            let v = by_ident.entry(i.key.as_str()).or_default();
            if !v.contains(p) {
                v.push(*p);
            }
        }
    }

    let mut out = RunOut {
        profile,
        floor: Vec::new(),
        value: Vec::new(),
        identity: Vec::new(),
        cands: Vec::new(),
        counts: Counts::default(),
    };
    // A pair of scans deviates by the root of the sum of both scans' squared
    // errors: against another noisy scan the wrong-value definition widens by
    // the square root of two.
    let tau = ctx.p.wrong_log2 * std::f64::consts::SQRT_2;

    // Floor rows: an identified scan against the nearest other identified
    // scan of its precursor.
    for (key, positions) in &by_ident {
        if positions.len() < 2 {
            continue;
        }
        let first = ids
            .iter()
            .find(|i| i.key == *key)
            .expect("a key of by_ident is in ids");
        let ladder = ctx.ladder(&first.peptidoform, first.charge);
        for &b in positions {
            let a = positions
                .iter()
                .copied()
                .filter(|&a| a != b)
                .min_by(|&x, &y| {
                    (scans[x].rt - scans[b].rt)
                        .abs()
                        .total_cmp(&(scans[y].rt - scans[b].rt).abs())
                })
                .expect("two or more positions");
            let Some(dev) = log2_deviation(&corrected[b], &corrected[a]) else {
                continue;
            };
            let wrong = dev > tau || out.profile.leak(&corrected[b]);
            let x = covariates(
                ctx,
                &scans[b],
                &corrected[b],
                &out.profile,
                &ladder,
                tol,
                first.charge,
            );
            out.floor.push(TrainRow {
                x: x.to_vec(),
                y: f64::from(u8::from(wrong)),
                run: r16,
                donors: Vec::new(),
            });
        }
    }

    // Entries: masked attempts (identified here, spectrum from another run)
    // train the models; the rest are transfers.
    let mut cands: Vec<TransferRow> = Vec::new();
    let mut half_windows: Vec<f64> = Vec::new();
    for e in &lib.entries {
        if e.identity_decoy {
            continue;
        }
        let pid = e.precursor;
        let m = &lib.meta[pid as usize];
        let Some(spec) = ctx.stored.spectra[pid as usize].as_ref() else {
            continue;
        };
        let ident_here = e.donors.binary_search(&r16).is_ok();
        if ident_here && !matches!(ctx.best_run.get(&m.key), Some(&(_, br)) if br != r16) {
            continue;
        }
        let pred = if ident_here {
            lib.alignment.predict_excluding(pid, r16)
        } else {
            lib.alignment.predict(pid, r16)
        };
        let Some((t_hat, s)) = pred else {
            continue;
        };
        if !(s > 0.0 && s.is_finite() && t_hat.is_finite()) {
            continue;
        }
        let window = (t_hat - k * s, t_hat + k * s);
        half_windows.push(k * s);
        let donors: Vec<u16> = e.donors.iter().copied().filter(|&d| d != r16).collect();
        let best = best_candidate(
            &index, spec, m.mz, m.charge, e.faims_cv, window, t_hat, s, &ppm, &free, tol,
        );
        let dec_mz = probes(pid, r16, m.mz, m.charge.max(1), t_hat, s, 0.0, (0.0, 1.0e9))
            .iter()
            .find(|p| matches!(p.kind, Hypothesis::MassShiftDecoy(_)))
            .map(|p| p.mz);
        let decoy = dec_mz.and_then(|mz| {
            best_candidate(
                &index, spec, mz, m.charge, e.faims_cv, window, t_hat, s, &ppm, &free, tol,
            )
        });
        let peptidoform = split_key(&m.key);
        if best.is_none() && decoy.is_none() {
            continue;
        }
        let ladder = ctx.ladder(peptidoform, m.charge);
        let donor_q = e.donor_psm.q;
        let scan_x = |f: &Found, ppm_w: &PpmWindow| {
            let sc = &scans[f.pos];
            (
                covariates(
                    ctx,
                    sc,
                    &corrected[f.pos],
                    &out.profile,
                    &ladder,
                    tol,
                    m.charge,
                ),
                transfer_features(f, ppm_w, donor_q, donors.len()),
            )
        };
        if ident_here {
            if let Some(f) = &best {
                let (sx, xx) = scan_x(f, &ppm);
                let mut full = sx.to_vec();
                full.extend_from_slice(&xx);
                out.counts.siblings += 1;
                if let Some(refs) = by_ident.get(m.key.as_str()) {
                    let reference = refs
                        .iter()
                        .copied()
                        .min_by(|&a, &b| {
                            (scans[a].rt - scans[f.pos].rt)
                                .abs()
                                .total_cmp(&(scans[b].rt - scans[f.pos].rt).abs())
                        })
                        .expect("a precursor identified here has a scan");
                    if let Some(dev) = log2_deviation(&corrected[f.pos], &corrected[reference]) {
                        let wrong = dev > tau || out.profile.leak(&corrected[f.pos]);
                        out.value.push(TrainRow {
                            x: full,
                            y: f64::from(u8::from(wrong)),
                            run: r16,
                            donors: donors.clone(),
                        });
                    }
                }
            }
        } else {
            // Target and decoy attempts of one entry have the same breadth:
            // together they train the identity model.
            if let Some(f) = &decoy {
                let (sx, xx) = scan_x(f, &ppm);
                let mut full = sx.to_vec();
                full.extend_from_slice(&xx);
                out.counts.decoys_real += 1;
                out.identity.push(TrainRow {
                    x: full,
                    y: 1.0,
                    run: r16,
                    donors: donors.clone(),
                });
            }
            let Some(f) = &best else {
                continue;
            };
            let (sx, xx) = scan_x(f, &ppm);
            {
                let mut full = sx.to_vec();
                full.extend_from_slice(&xx);
                out.counts.targets_real += 1;
                out.identity.push(TrainRow {
                    x: full,
                    y: 0.0,
                    run: r16,
                    donors: donors.clone(),
                });
            }
            let sc = &scans[f.pos];
            let rep = sc.reporters.as_ref().expect("checked above");
            let donor_run = ctx
                .best_run
                .get(&m.key)
                .map(|&(_, d)| ctx.stems[d as usize].clone())
                .unwrap_or_default();
            cands.push(TransferRow {
                run: stem.clone(),
                native_id: sc.native_id.clone(),
                spec_key: sc.spec_key.clone(),
                scan: scan_number(&sc.native_id),
                rt: sc.rt,
                precursor_mz: sc.precursor_mz,
                peptidoform: peptidoform.to_string(),
                charge: m.charge,
                calc_mz: m.mz,
                proteins: ctx.stored.proteins[pid as usize]
                    .split(';')
                    .filter(|a| !a.is_empty())
                    .map(str::to_string)
                    .collect(),
                raw: rep.raw.clone(),
                corrected: corrected[f.pos].clone(),
                purity: rep.purity,
                info: TransferInfo {
                    quant_q: 1.0,
                    p_wrong: 1.0,
                    noise_floor: 0.0,
                    p_identity: 1.0,
                    p_value: 1.0,
                    angle: f.m.angle,
                    n_matched: f.m.n_matched,
                    explained: f.m.explained,
                    rt_z: f.rt_z as f32,
                    ppm: f.ppm as f32,
                    donor_q,
                    n_donors: e.donors.len().min(u16::MAX as usize) as u16,
                    donor_run,
                    batch_cos: f32::NAN,
                },
                accepted: false,
                scan_features: sx,
                xfer: xx,
                donors: donors.clone(),
                run_index: r16,
                donor_pep: e.donor_psm.pep,
            });
        }
    }

    // One precursor per scan: the best match keeps it.
    let mut winner: HashMap<&str, usize> = HashMap::new();
    let key_of = |c: &TransferRow| (c.info.angle, c.info.n_matched);
    for (i, c) in cands.iter().enumerate() {
        match winner.get(c.spec_key.as_str()) {
            Some(&w) if key_of(&cands[w]) >= key_of(c) => {}
            _ => {
                winner.insert(c.spec_key.as_str(), i);
            }
        }
    }
    let keep: std::collections::HashSet<usize> = winner.values().copied().collect();
    drop(winner);
    out.cands = cands
        .into_iter()
        .enumerate()
        .filter(|(i, _)| keep.contains(i))
        .map(|(_, c)| c)
        .collect();
    eprintln!(
        "[tmt-mbr] {stem}: {} identified PSMs on {} scans; precursor window {:+.2} +- {:.2} ppm{}, RT window {:.1} predicted widths (median half-window {:.1} s){}; {} empty channel(s) detected; training rows: floor {}, value {}, identity {}; {} transfer candidates",
        ids.len(),
        seen.len(),
        ppm.centre,
        ppm.half_width,
        if ppm_errors.len() < quant::tmt_transfer::MIN_CALIBRATION { " (fallback)" } else { "" },
        k,
        quant::tmt_transfer::median_finite(&half_windows),
        if k_learned { "" } else { " (fallback)" },
        out.profile.empty.len(),
        out.floor.len(),
        out.value.len(),
        out.identity.len(),
        out.cands.len()
    );
    Ok(out)
}

/// Library, candidates, risk model and `quant_q`. `make_scores` gives the
/// rescoring result (called twice: the library takes one set, the identified
/// scans are joined to the other).
pub(crate) fn run(
    dir: &Path,
    stems: &[String],
    make_scores: &dyn Fn() -> Result<ScoreSet, String>,
    p: &Params<'_>,
) -> Result<(Vec<TransferRow>, Summary), String> {
    let t0 = std::time::Instant::now();
    let stems = library_runs(dir, stems)?;
    let stored = build_and_write_library(
        dir,
        &stems,
        ScoreLoader::Merged(Box::new(|_| make_scores())),
        p.library_fdr,
        p.donor_q,
        p.donor_pep,
    )?;
    let scores = make_scores()?;
    let mut mods: HashMap<String, f64> = HashMap::new();
    for s in &stems {
        mods.extend(read_run_quant_mods(&artifact(dir, s, TABLE_SUFFIX))?);
    }
    let mut best_run: HashMap<String, (f64, u16)> = HashMap::new();
    for (r, s) in stems.iter().enumerate() {
        for row in read_ident(dir, s, &scores, p)? {
            let r16 = r as u16;
            match best_run.get_mut(&row.key) {
                Some(b) if b.0 >= row.score => {}
                Some(b) => *b = (row.score, r16),
                None => {
                    best_run.insert(row.key, (row.score, r16));
                }
            }
        }
    }
    let key_entry: HashMap<&str, usize> = stored
        .library
        .entries
        .iter()
        .enumerate()
        .filter(|(_, e)| !e.identity_decoy)
        .map(|(i, e)| (stored.library.meta[e.precursor as usize].key.as_str(), i))
        .collect();
    let reporter_top = p
        .plex
        .channels()
        .iter()
        .fold(f64::NEG_INFINITY, |a, c| a.max(c.mz));
    let ctx = Ctx {
        dir,
        stems: &stems,
        stored: &stored,
        scores: &scores,
        p,
        mods,
        best_run,
        key_entry,
        min_mz: reporter_top + 1.0,
    };
    let mut outs = Vec::with_capacity(stems.len());
    for r in 0..stems.len() {
        outs.push(process_run(&ctx, r)?);
    }

    // Cross-run covariate and training set.
    let profiles: Vec<RunProfile> = outs.iter().map(|o| o.profile.clone()).collect();
    let batch_cos = |run: u16, donors: &[u16]| -> f64 {
        donors
            .iter()
            .map(|&d| profiles[run as usize].loading_cosine(&profiles[d as usize]))
            .filter(|c| c.is_finite())
            .fold(f64::NAN, f64::max)
    };
    let mut training = Training::default();
    let mut cands: Vec<TransferRow> = Vec::new();
    let mut counts = Counts::default();
    for o in outs {
        for row in o.floor {
            training.floor_x.push(row.x);
            training.floor_y.push(row.y);
        }
        for mut row in o.value {
            row.x[N_SCAN + BATCH_COS_AT] = batch_cos(row.run, &row.donors);
            training.value_x.push(row.x);
            training.value_y.push(row.y);
        }
        for mut row in o.identity {
            row.x[N_SCAN + BATCH_COS_AT] = batch_cos(row.run, &row.donors);
            training.id_x.push(row.x);
            training.id_y.push(row.y);
        }
        counts.siblings += o.counts.siblings;
        counts.targets_real += o.counts.targets_real;
        counts.decoys_real += o.counts.decoys_real;
        cands.extend(o.cands);
    }
    for c in &mut cands {
        let cos = batch_cos(c.run_index, &c.donors);
        c.xfer[BATCH_COS_AT] = cos;
        c.info.batch_cos = cos as f32;
    }
    let mut summary = Summary {
        candidates: cands.len(),
        ..Summary::default()
    };
    match TransferRisk::fit(&training) {
        Ok(risk) => {
            summary.model_fitted = true;
            if !risk.value_fitted() {
                eprintln!(
                    "WARN: [tmt-mbr] too few masked sibling scans ({}) to fit the value-risk model: the identified noise floor stands in for it, so only identity and floor-level risk are priced",
                    training.value_y.len()
                );
            }
            let mut pw = Vec::with_capacity(cands.len());
            let mut fl = Vec::with_capacity(cands.len());
            let mut parts = Vec::with_capacity(cands.len());
            for c in &cands {
                let (p_scan, g, f) = risk.parts(&c.scan_features, &c.xfer);
                // The donor's own identification can be false too.
                let p_id = p_scan + (1.0 - p_scan) * c.donor_pep.clamp(0.0, 1.0);
                pw.push(p_id + (1.0 - p_id) * g);
                // Identified rows carry no identity error beyond their FDR, so
                // the identity risk is all excess; the value risk is excess
                // only above the identified noise floor.
                fl.push(pw[pw.len() - 1] - (p_id + (1.0 - p_id) * (g - f).max(0.0)));
                parts.push((p_id, g, f));
            }
            let q = quant_q(&pw, &fl);
            for (i, c) in cands.iter_mut().enumerate() {
                c.info.p_wrong = pw[i] as f32;
                c.info.noise_floor = parts[i].2 as f32;
                c.info.p_identity = parts[i].0 as f32;
                c.info.p_value = parts[i].1 as f32;
                c.info.quant_q = q[i];
                c.accepted = q[i] <= p.accept_q;
            }
            summary.accepted = cands.iter().filter(|c| c.accepted).count();
            let masked_pred: Vec<f64> = training
                .value_x
                .iter()
                .map(|x| {
                    let mut s = [0.0; N_SCAN];
                    s.copy_from_slice(&x[..N_SCAN]);
                    let mut t = [0.0; N_XFER];
                    t.copy_from_slice(&x[N_SCAN..]);
                    risk.value_risk(&s, &t)
                })
                .collect();
            let table: Vec<String> = calibration(&masked_pred, &training.value_y, 5)
                .iter()
                .map(|(p, o, n)| format!("{p:.4}->{o:.4} (n={n})"))
                .collect();
            eprintln!(
                "[tmt-mbr] masked-pair calibration of the value risk, predicted -> observed wrong rate by quintile: {}",
                table.join(", ")
            );
            eprintln!(
                "[tmt-mbr] risk model: floor {} rows (wrong {:.4}), value {} rows (wrong {:.4}; {} masked siblings), identity {} rows ({} target and {} decoy candidates of the real attempts); mean predicted over the transfers: identity {:.4}, value {:.4}, floor {:.4}; {} of {} transfers at quant_q<={}",
                training.floor_y.len(),
                training.floor_y.iter().sum::<f64>() / training.floor_y.len().max(1) as f64,
                training.value_y.len(),
                training.value_y.iter().sum::<f64>() / training.value_y.len().max(1) as f64,
                counts.siblings,
                training.id_y.len(),
                counts.targets_real,
                counts.decoys_real,
                cands.iter().map(|c| f64::from(c.info.p_identity)).sum::<f64>() / cands.len().max(1) as f64,
                cands.iter().map(|c| f64::from(c.info.p_value)).sum::<f64>() / cands.len().max(1) as f64,
                cands.iter().map(|c| f64::from(c.info.noise_floor)).sum::<f64>() / cands.len().max(1) as f64,
                summary.accepted,
                cands.len(),
                p.accept_q
            );
        }
        Err(e) => eprintln!(
            "WARN: [tmt-mbr] {e}; the {} transfer candidates are reported with quant_q 1 and are not accepted",
            cands.len()
        ),
    }
    eprintln!("[PHASE tmt_mbr: {:.2}s]", t0.elapsed().as_secs_f64());
    Ok((cands, summary))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scan_numbers_are_read_from_native_ids() {
        assert_eq!(
            scan_number("controllerType=0 controllerNumber=1 scan=1234"),
            1234
        );
        assert_eq!(scan_number("scan=7"), 7);
        assert_eq!(scan_number("index=5"), 0);
    }

    type Peaks = Vec<(f32, f32)>;

    fn spectrum_and_peaks() -> (LibrarySpectrum, Peaks, Peaks) {
        let ladder = fragment_ladder(
            &residue_masses("LGEYGFQNALIVR", &HashMap::new()).unwrap(),
            2,
        );
        let peaks: Vec<(f32, f32)> = ladder
            .iter()
            .enumerate()
            .map(|(k, f)| (f.mz as f32, 100.0 + k as f32))
            .collect();
        let spec = quant::ms2::library_spectrum(&ladder, &peaks, Tolerance::Ppm(20.0)).unwrap();
        (spec, peaks, vec![(321.1, 10.0), (654.3, 20.0)])
    }

    fn scan(rt: f64, mz: f64, charge: u8, peaks: &[(f32, f32)]) -> quant::ms2::Ms2Scan {
        quant::ms2::Ms2Scan {
            native_id: format!("scan={}", (rt * 10.0) as i64),
            spec_key: format!("k{rt}"),
            rt,
            precursor_mz: mz,
            charge,
            isolation: Some((mz, 0.7, 0.7)),
            faims_cv: None,
            peaks: peaks.to_vec(),
            reporters: Some(quant::ms2::ScanReporters {
                raw: vec![1.0; 11],
                purity: Some(0.9),
                sn_sum: None,
            }),
        }
    }

    #[test]
    fn candidate_rule_is_mz_charge_rt_and_free_scan_and_picks_the_best_match() {
        let (spec, donor_like, noise) = spectrum_and_peaks();
        let mz = 700.0;
        let index = Ms2RunIndex::new(vec![
            scan(100.0, mz, 2, &donor_like),
            scan(101.0, mz * (1.0 + 4e-6), 2, &noise),
            scan(102.0, mz * (1.0 + 40e-6), 2, &donor_like),
            scan(103.0, mz, 3, &donor_like),
            scan(104.0, mz, 0, &donor_like),
            scan(160.0, mz, 2, &donor_like),
            scan(105.0, mz, 2, &donor_like),
        ]);
        let ppm = PpmWindow {
            centre: 0.0,
            half_width: 10.0,
        };
        let tol = Tolerance::Ppm(20.0);
        let mut free = vec![true; index.len()];
        let found = best_candidate(
            &index,
            &spec,
            mz,
            2,
            None,
            (90.0, 120.0),
            105.0,
            5.0,
            &ppm,
            &free,
            tol,
        )
        .unwrap();
        // The scan at 102 s is 40 ppm away, 103 s has another charge, 160 s is
        // outside the window; 100, 101, 104 (charge unrecorded) and 105 are
        // eligible and the donor-like ones tie on angle, so the one closest in
        // retention time wins.
        assert_eq!(index.scans()[found.pos].rt, 105.0);
        assert!(found.m.angle > 0.9 && found.rt_z.abs() < 1e-9);
        // An identified scan is not a candidate.
        let taken = index.scans().iter().position(|s| s.rt == 105.0).unwrap();
        free[taken] = false;
        let next = best_candidate(
            &index,
            &spec,
            mz,
            2,
            None,
            (90.0, 120.0),
            105.0,
            5.0,
            &ppm,
            &free,
            tol,
        )
        .unwrap();
        assert_eq!(index.scans()[next.pos].rt, 104.0);
        // A window with nothing in it, and a donor m/z nothing matches.
        assert!(best_candidate(
            &index,
            &spec,
            mz,
            2,
            None,
            (200.0, 210.0),
            205.0,
            5.0,
            &ppm,
            &free,
            tol
        )
        .is_none());
        assert!(best_candidate(
            &index,
            &spec,
            900.0,
            2,
            None,
            (90.0, 120.0),
            105.0,
            5.0,
            &ppm,
            &free,
            tol
        )
        .is_none());
        // A scan with a lone noise spectrum is still a candidate but a poor one.
        let only_noise = Ms2RunIndex::new(vec![scan(100.0, mz, 2, &noise)]);
        let weak = best_candidate(
            &only_noise,
            &spec,
            mz,
            2,
            None,
            (90.0, 110.0),
            100.0,
            5.0,
            &ppm,
            &[true],
            tol,
        )
        .unwrap();
        assert!(weak.m.angle < found.m.angle);
    }
}
