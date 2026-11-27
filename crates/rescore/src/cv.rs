use crate::budget::Reservation;
use crate::codec::RowCodec;
use crate::direction::{init_directions, pin_direction};
use crate::folds::folds_by_precursor;
use crate::normalize::{per_row_sorted, raw_weights, score, Norm};
use crate::pin::Pin;
use crate::scores::tdc_cut;
use crate::svm::{mfn_from, warm_scores, Examples, MfnOpts, Rows, GRAM_BYTES};
use rayon::prelude::*;
use std::sync::Arc;

const FOLDS: usize = 3;

const GRID_ITERATIONS: usize = 3;

const MIN_GROWTH: f64 = 0.001;

pub struct CvConfig {
    pub iterations: usize,
    pub test_fdr: f64,
    pub selection_fdr: f64,
    pub initial_selection_fdr: f64,
    pub cpos: Vec<f64>,
    pub cfrac: Vec<f64>,
    pub seed: u64,
    pub overrule: bool,
    pub verbose: u32,
}

impl Default for CvConfig {
    fn default() -> Self {
        CvConfig {
            iterations: 10,
            test_fdr: 0.01,
            selection_fdr: 0.01,
            initial_selection_fdr: 0.01,
            cpos: vec![10.0, 1.0, 0.1],
            cfrac: vec![1.0, 3.0, 10.0],
            seed: 1,
            overrule: false,
            verbose: 0,
        }
    }
}

pub struct CvResult {
    pub w_norm: [Vec<f64>; 3],
    pub w_raw: [Vec<f64>; 3],
    pub test_scores: Vec<f64>,
    pub fold: Vec<u8>,
}

struct Fold {
    train: Vec<usize>,
    test: Vec<usize>,
    train_decoy: Vec<bool>,
    test_decoy: Vec<bool>,
}

struct TrainSet {
    ex: Examples,
    negatives: usize,
}

/// Bytes one SVM fit holds per training example: `o`, `o_bar`, `subset`, a `deltas` entry, and
/// when screening the `act` flag and the row norm.
pub(crate) const TASK_BYTES_PER_EXAMPLE: u64 = (3 * std::mem::size_of::<f64>()
    + std::mem::size_of::<(f64, usize, f64)>()
    + std::mem::size_of::<bool>()
    + std::mem::size_of::<f32>()) as u64;

fn score_fits(pin: &Pin, norm: &Norm, rows: &[usize], ws: &[&[f64]]) -> Vec<Vec<f64>> {
    let nf = pin.rows.nfeat();
    if nf == 0 {
        return ws.iter().map(|w| vec![w[0]; rows.len()]).collect();
    }
    per_row_sorted(&pin.rows, rows, ws.len(), |row, out, normed| {
        normed.clear();
        normed.extend(
            row.iter()
                .zip(&norm.mean)
                .zip(&norm.inv_sd)
                .map(|((&x, &m), &sd)| (f64::from(x) - m) * sd),
        );
        for (o, w) in out.iter_mut().zip(ws) {
            *o = dot(normed, w);
        }
    })
}

fn dot(row: &[f64], w: &[f64]) -> f64 {
    let nf = w.len() - 1;
    let mut s = w[nf];
    for (&x, &wj) in row[..nf].iter().zip(&w[..nf]).rev() {
        s += x * wj;
    }
    s
}

const SCORE_BYTES_PER_ROW: u64 = 8;

fn fit_positives(
    pin: &Pin,
    norm: &Norm,
    fold: &Fold,
    fits: &[Option<Vec<f64>>],
    fdr: f64,
) -> Option<(usize, usize, Vec<f64>)> {
    let present: Vec<usize> = (0..fits.len()).filter(|&k| fits[k].is_some()).collect();
    let per_fit = (fold.train.len() as u64).max(1) * SCORE_BYTES_PER_ROW;
    let _best_held = Reservation::take(&pin.budget, per_fit);
    let group = ((pin.budget.available() / per_fit) as usize).clamp(1, present.len().max(1));
    let mut best: Option<(usize, usize, Vec<f64>)> = None;
    for ks in present.chunks(group) {
        let _held = Reservation::take(&pin.budget, ks.len() as u64 * per_fit);
        let ws: Vec<&[f64]> = ks
            .iter()
            .map(|&k| fits[k].as_deref().expect("present fit"))
            .collect();
        let scores = score_fits(pin, norm, &fold.train, &ws);
        let counted: Vec<(usize, Vec<f64>)> = scores
            .into_par_iter()
            .map(|s| (tdc_cut(&s, &fold.train_decoy, fdr, true).below, s))
            .collect();
        for (&k, (tp, s)) in ks.iter().zip(counted) {
            if best.as_ref().is_none_or(|b| tp >= b.1) {
                best = Some((k, tp, s));
            }
        }
    }
    best
}

fn train_set(fold: &Fold, s: &[f64], fdr: f64) -> TrainSet {
    let cut = tdc_cut(s, &fold.train_decoy, fdr, true).at_most;
    let mut ex = Examples {
        rows: Vec::new(),
        pos: Vec::new(),
    };
    let mut negatives = 0usize;
    for ((&row, &v), &d) in fold.train.iter().zip(s).zip(&fold.train_decoy) {
        if d {
            ex.rows.push(row);
            ex.pos.push(false);
            negatives += 1;
        } else if cut.is_some_and(|c| v >= c) {
            ex.rows.push(row);
            ex.pos.push(true);
        }
    }
    TrainSet { ex, negatives }
}

pub(crate) fn make_room(pin: &mut Pin) -> Result<(), String> {
    let n = pin.rows.len();
    let need = Rows::bytes_for(n, pin.rows.nfeat()) + (n as u64 * 2 / 3) * TASK_BYTES_PER_EXAMPLE;
    if pin.budget.available() < need {
        pin.texts.spill(&pin.budget)?;
    }
    Ok(())
}

fn fit_reserve(folds: &[Fold]) -> u64 {
    folds
        .iter()
        .map(|f| f.train.len() as u64 * TASK_BYTES_PER_EXAMPLE)
        .max()
        .unwrap_or(0)
}

fn example_union(sets: &[TrainSet]) -> Vec<usize> {
    let mut all: Vec<usize> = sets
        .iter()
        .flat_map(|t| t.ex.rows.iter().copied())
        .collect();
    all.par_sort_unstable();
    all.dedup();
    all
}

fn to_positions(rows: &mut [usize], union: &[usize]) {
    let mut j = 0usize;
    for r in rows.iter_mut() {
        while union[j] < *r {
            j += 1;
        }
        *r = j;
    }
}

fn positives(s: &[f64], decoy: &[bool], fdr: f64, skip_decoys_plus_one: bool) -> usize {
    tdc_cut(s, decoy, fdr, skip_decoys_plus_one).below
}

/// Trains the cross-validated semi-supervised rescorer on a PIN table, after the Percolator
/// recipe [cite:kall2007_percolator].
pub fn train(
    pin: &Pin,
    norm: &Norm,
    cfg: &CvConfig,
    pool: &rayon::ThreadPool,
) -> Result<CvResult, String> {
    pool.install(|| train_in_pool(pin, norm, cfg))
}

fn train_in_pool(pin: &Pin, norm: &Norm, cfg: &CvConfig) -> Result<CvResult, String> {
    let n = pin.meta.len();
    let fold = folds_by_precursor(&pin.meta, FOLDS, cfg.seed);
    let is_decoy: Vec<bool> = pin.meta.iter().map(|m| m.label < 0).collect();
    let folds = make_folds(&fold, &is_decoy);
    let targets = is_decoy.iter().filter(|&&d| !d).count();
    if targets == 0 {
        return Err("Error: no target PSMs were provided.".to_string());
    }
    if targets == n {
        return Err("Error: no decoy PSMs were provided.".to_string());
    }
    let ratio = targets as f64 / (n - targets) as f64;
    let pairs: Vec<(usize, f64, f64)> = (0..FOLDS)
        .flat_map(|f| {
            cfg.cpos
                .iter()
                .flat_map(move |&cpos| cfg.cfrac.iter().map(move |&c| (f, cpos, c * ratio)))
        })
        .collect();
    let per_fold = cfg.cpos.len() * cfg.cfrac.len();
    let opts = MfnOpts::default();

    let all_sets = (1u8 << FOLDS) - 1;
    let mask: Vec<u8> = fold.iter().map(|&f| all_sets & !(1u8 << f)).collect();
    let init = init_directions(
        pin,
        norm,
        &mask,
        FOLDS,
        cfg.initial_selection_fdr,
        cfg.overrule,
    )?;
    let init_positives = initial_positives(pin, norm, cfg, &folds, &init);
    let mut w = init.clone();
    let mut train_scores: Vec<Vec<f64>> = folds
        .iter()
        .zip(&w)
        .map(|(fd, w)| score(&pin.rows, &fd.train, w, Some(norm)))
        .collect();

    let keep = fit_reserve(&folds);
    let full = Rows::all_in_ram(&pin.rows, norm, &pin.budget, keep);
    if full.is_none() {
        log::debug!(target: "rescore", "training rows gathered per iteration (the budget does not hold every row)");
    }
    let mut codec: Option<Arc<RowCodec>> = None;
    let mut warm: Vec<Option<Vec<f64>>> = vec![None; pairs.len()];
    let mut fixed: Option<[usize; FOLDS]> = None;
    let (mut found_old_old, mut found_old) = (0i64, 0i64);
    let mut widest = 0usize;
    for it in 0..cfg.iterations {
        let fdr = if it == 0 {
            cfg.initial_selection_fdr
        } else {
            cfg.selection_fdr
        };
        let mut sets: Vec<TrainSet> = folds
            .par_iter()
            .zip(&train_scores)
            .map(|(fd, s)| train_set(fd, s, fdr))
            .collect();
        let compact;
        let x: &Rows = match &full {
            Some(x) => x,
            None => {
                let union = example_union(&sets);
                let codec = codec.get_or_insert_with(|| {
                    let stride = crate::svm::stride_for(pin.rows.nfeat() + 1);
                    Arc::new(RowCodec::scan(&pin.rows, norm, stride))
                });
                compact = Rows::gather_within(&pin.rows, norm, &union, &pin.budget, codec)?;
                if compact.is_spilled() {
                    log::debug!(
                        target: "rescore",
                        "iteration {}: training rows spilled, {} bytes a row ({} of {} features as one byte)",
                        it + 1,
                        codec.width(),
                        codec.coded(),
                        pin.rows.nfeat()
                    );
                }
                for t in sets.iter_mut() {
                    to_positions(&mut t.ex.rows, &union);
                }
                &compact
            }
        };
        if cfg.verbose > 1 {
            log::info!(target: "rescore", "Iteration {}:", it + 1);
        }
        if cfg.verbose > 2 {
            for (f, t) in sets.iter().enumerate() {
                log::info!(
                    target: "rescore",
                    "Split {}: Training with {} positives and {} negatives",
                    f + 1,
                    t.ex.rows.len() - t.negatives,
                    t.negatives
                );
            }
        }
        let task_bytes: Vec<u64> = pairs
            .iter()
            .map(|&(f, _, _)| {
                sets[f].ex.rows.len() as u64 * TASK_BYTES_PER_EXAMPLE + GRAM_BYTES as u64
            })
            .collect();
        let mut fits: Vec<Option<Vec<f64>>> = vec![None; pairs.len()];
        let fit_one = |k: usize,
                       start: Option<&[f64]>,
                       pre: Option<(Vec<f64>, &[f32])>|
         -> Option<Vec<f64>> {
            let (f, cpos, cfrac) = pairs[k];
            let t = &sets[f];
            if t.ex.rows.is_empty() {
                return None;
            }
            Some(mfn_from(x, &t.ex, cpos, cfrac * cpos, &opts, start, pre))
        };
        let ks: Vec<usize> = match fixed {
            Some(sel) => (0..FOLDS).map(|f| f * per_fold + sel[f]).collect(),
            None => (0..pairs.len()).collect(),
        };
        {
            let mut done = 0usize;
            while done < ks.len() {
                let mut held = Vec::new();
                while done + held.len() < ks.len() {
                    match Reservation::take(&pin.budget, task_bytes[ks[done + held.len()]]) {
                        Some(r) => held.push(r),
                        None => break,
                    }
                }
                let end = (done + held.len()).max(done + 1);
                widest = widest.max(end - done);
                let batch = &ks[done..end];
                let starts: Vec<Option<&[f64]>> =
                    batch.iter().map(|&k| warm[k].as_deref()).collect();
                let (pre, norms) = warm_starts(x, &sets, &pairs, batch, &starts);
                let norms = &norms;
                let got: Vec<Option<Vec<f64>>> = batch
                    .par_iter()
                    .zip(starts.par_iter())
                    .zip(pre.into_par_iter())
                    .map(|((&k, &start), o)| {
                        let nrm = norms[pairs[k].0].as_slice();
                        fit_one(k, start, o.map(|o| (o, nrm)))
                    })
                    .collect();
                for (&k, g) in ks[done..end].iter().zip(got) {
                    fits[k] = g;
                }
                done = end;
            }
        }
        for (wk, fit) in warm.iter_mut().zip(&fits) {
            if fit.is_some() {
                wk.clone_from(fit);
            }
        }
        drop(sets);
        let mut picked = [0usize; FOLDS];
        for f in 0..FOLDS {
            let fold_fits = &fits[f * per_fold..(f + 1) * per_fold];
            match fit_positives(pin, norm, &folds[f], fold_fits, cfg.test_fdr) {
                Some((k, best, s)) => {
                    train_scores[f] = s;
                    picked[f] = k;
                    w[f].clone_from(fold_fits[k].as_ref().expect("chosen fit"));
                    if cfg.verbose > 2 {
                        let (_, cpos, cfrac) = pairs[f * per_fold + k];
                        log::info!(
                            target: "rescore",
                            "Split {}: selected Cpos={cpos} Cneg={} with {best} training set positives",
                            f + 1,
                            cfrac * cpos
                        );
                    }
                }
                None => {
                    train_scores[f] =
                        score_fits(pin, norm, &folds[f].train, &[&w[f]]).swap_remove(0);
                }
            }
        }
        if fixed.is_none() && it + 1 >= GRID_ITERATIONS && per_fold > 1 {
            fixed = Some(picked);
            if cfg.verbose > 1 {
                for (f, &k) in picked.iter().enumerate() {
                    let (_, cpos, cfrac) = pairs[f * per_fold + k];
                    log::info!(
                        target: "rescore",
                        "Split {}: using Cpos={cpos}, Cneg={} in the remaining iterations",
                        f + 1,
                        cfrac * cpos
                    );
                }
            }
        }
        let found: usize = folds
            .par_iter()
            .zip(&train_scores)
            .map(|(fd, s)| positives(s, &fd.train_decoy, cfg.test_fdr, false))
            .sum();
        let estimated = (found as f64 / (FOLDS - 1) as f64) as i64;
        if cfg.verbose > 1 {
            log::info!(target: "rescore", "Estimated {estimated} PSMs with q<{}", cfg.test_fdr);
        }
        if estimated > 0
            && found_old_old > 0
            && ((estimated - found_old_old) as f64) <= found_old_old as f64 * MIN_GROWTH
        {
            if cfg.verbose > 0 {
                log::info!(
                    target: "rescore",
                    "Performance increase over the last two iterations indicate that the algorithm has converged ({estimated} vs {found_old_old})"
                );
            }
            break;
        }
        found_old_old = found_old;
        found_old = estimated;
    }

    log::debug!(target: "rescore", "at most {widest} SVM fits ran at once");
    let mut test = fold_test_scores(pin, norm, &folds, &w);
    let found: usize = folds
        .iter()
        .zip(&test)
        .map(|(fd, s)| positives(s, &fd.test_decoy, cfg.test_fdr, false))
        .sum();
    if cfg.verbose > 0 {
        log::info!(target: "rescore", "Found {found} test set PSMs with q<{}.", cfg.test_fdr);
    }
    if (found == 0 || init_positives > found) && !cfg.overrule {
        w = init;
        test = fold_test_scores(pin, norm, &folds, &w);
    }

    let test_scores = scatter_test_scores(&folds, &test, n);
    Ok(CvResult {
        w_raw: std::array::from_fn(|f| raw_weights(norm, &w[f])),
        w_norm: std::array::from_fn(|f| w[f].clone()),
        test_scores,
        fold,
    })
}

/// The `FOLDS` cross-validation folds: each fold's test rows are the rows assigned to it, its
/// training rows all the others, with their decoy labels.
fn make_folds(fold: &[u8], is_decoy: &[bool]) -> Vec<Fold> {
    let n = fold.len();
    (0..FOLDS)
        .into_par_iter()
        .map(|f| {
            let (test, train): (Vec<usize>, Vec<usize>) =
                (0..n).partition(|&i| fold[i] as usize == f);
            Fold {
                train_decoy: train.iter().map(|&i| is_decoy[i]).collect(),
                test_decoy: test.iter().map(|&i| is_decoy[i]).collect(),
                train,
                test,
            }
        })
        .collect()
}

/// For the fits of `batch` that start from a previous solution, their starting outputs on the
/// fold's training examples, and each fold's example norms.
fn warm_starts(
    x: &Rows,
    sets: &[TrainSet],
    pairs: &[(usize, f64, f64)],
    batch: &[usize],
    starts: &[Option<&[f64]>],
) -> (Vec<Option<Vec<f64>>>, Vec<Vec<f32>>) {
    let mut pre: Vec<Option<Vec<f64>>> = batch.iter().map(|_| None).collect();
    let mut norms: Vec<Vec<f32>> = vec![Vec::new(); FOLDS];
    for (f, nf) in norms.iter_mut().enumerate() {
        if sets[f].ex.rows.is_empty() {
            continue;
        }
        let js: Vec<usize> = (0..batch.len())
            .filter(|&j| pairs[batch[j]].0 == f && starts[j].is_some())
            .collect();
        if js.is_empty() {
            continue;
        }
        let ws: Vec<&[f64]> = js.iter().filter_map(|&j| starts[j]).collect();
        let (os, nrm) = warm_scores(x, &sets[f].ex, &ws);
        for (&j, o) in js.iter().zip(os) {
            pre[j] = Some(o);
        }
        *nf = nrm;
    }
    (pre, norms)
}

/// Each fold's scores of its test rows under its weights `w`.
fn fold_test_scores(pin: &Pin, norm: &Norm, folds: &[Fold], w: &[Vec<f64>]) -> Vec<Vec<f64>> {
    folds
        .iter()
        .zip(w)
        .map(|(fd, w)| score(&pin.rows, &fd.test, w, Some(norm)))
        .collect()
}

/// Every row's score from the fold that tests it.
fn scatter_test_scores(folds: &[Fold], test: &[Vec<f64>], n: usize) -> Vec<f64> {
    let mut test_scores = vec![0.0f64; n];
    for (fd, s) in folds.iter().zip(test) {
        for (&i, &v) in fd.test.iter().zip(s) {
            test_scores[i] = v;
        }
    }
    test_scores
}

/// The test-set positives the initial directions separate, summed over the folds; at verbose
/// above 1, without a fixed PIN direction, also logs each fold's chosen initial feature.
fn initial_positives(
    pin: &Pin,
    norm: &Norm,
    cfg: &CvConfig,
    folds: &[Fold],
    init: &[Vec<f64>],
) -> usize {
    let nf = pin.rows.nfeat();
    if cfg.verbose > 1 && pin_direction(pin).is_none() {
        for (f, (fd, w)) in folds.iter().zip(init).enumerate() {
            let s = score(&pin.rows, &fd.train, w, Some(norm));
            let p = positives(&s, &fd.train_decoy, cfg.initial_selection_fdr, true);
            let j = w[..nf].iter().position(|&v| v != 0.0).map_or(0, |j| j + 1);
            log::info!(
                target: "rescore",
                "Split {}:\tSelected feature {j} as initial direction. Could separate {p} training set positives with q<{} in that direction.",
                f + 1,
                cfg.initial_selection_fdr
            );
        }
    }
    let init_positives: usize = folds
        .iter()
        .zip(init)
        .map(|(fd, w)| {
            let s = score(&pin.rows, &fd.test, w, Some(norm));
            positives(&s, &fd.test_decoy, cfg.test_fdr, false)
        })
        .sum();
    if cfg.verbose > 0 {
        log::info!(
            target: "rescore",
            "Found {init_positives} test set positives with q<{} in initial direction",
            cfg.test_fdr
        );
    }
    init_positives
}
