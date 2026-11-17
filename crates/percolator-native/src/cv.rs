use crate::budget::Reservation;
use crate::folds::folds_by_spectrum;
use crate::normalize::{raw_weights, score, Norm};
use crate::pin::Pin;
use crate::scores::{count_at, tdc_q, training_sets};
use crate::svm::{fill, l2_svm_mfn, MfnOpts, Rows};
use rayon::prelude::*;
use std::cmp::Ordering;

const FOLDS: usize = 3;

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
    x: Rows,
    idx_neg: Vec<usize>,
    idx_pos: Vec<usize>,
}

fn pin_direction(pin: &Pin) -> Option<&[f64]> {
    pin.default_direction
        .as_deref()
        .filter(|d| d.iter().any(|&v| v != 0.0))
}

const NO_DIRECTION: &str = "Error in the input data: cannot find an initial direction with positive training examples. Consider setting/raising the initial training FDR threshold (--train-fdr-initial). Terminating.";

pub fn init_direction(
    pin: &Pin,
    norm: &Norm,
    rows: &[usize],
    fdr: f64,
    overrule: bool,
) -> Result<Vec<f64>, String> {
    let nf = pin.rows.nfeat();
    let mut w = vec![0.0f64; nf + 1];
    if let Some(d) = pin_direction(pin) {
        for (wi, &di) in w.iter_mut().zip(d) {
            *wi = di;
        }
        return Ok(w);
    }
    if nf == 0 {
        return Ok(w);
    }
    let decoy: Vec<bool> = rows.iter().map(|&i| pin.meta[i].label < 0).collect();
    let per_feature = (rows.len() as u64).max(1) * COLUMN_BYTES_PER_ROW;
    let chunk = ((pin.budget.available() / per_feature) as usize).clamp(1, nf);
    let mut counts: Vec<[usize; 2]> = Vec::with_capacity(nf);
    for j0 in (0..nf).step_by(chunk) {
        let j1 = (j0 + chunk).min(nf);
        let _held = Reservation::take(&pin.budget, (j1 - j0) as u64 * per_feature);
        let mut cols = vec![vec![0.0f64; rows.len()]; j1 - j0];
        fill(&pin.rows, norm, rows, |p, row| {
            for (c, &v) in cols.iter_mut().zip(&row[j0..j1]) {
                c[p] = v;
            }
        });
        counts.par_extend(cols.into_par_iter().map(|high| {
            let low: Vec<f64> = high.iter().map(|&v| -v).collect();
            [
                count_at(&tdc_q(&low, &decoy, true), &decoy, fdr),
                count_at(&tdc_q(&high, &decoy, true), &decoy, fdr),
            ]
        }));
    }
    let mut best: Option<(usize, usize, bool)> = None;
    for (j, c) in counts.iter().enumerate() {
        for (i, &p) in c.iter().enumerate() {
            let better = match best {
                None => true,
                Some((b, _, _)) => p > b,
            };
            if better {
                best = Some((p, j, i == 0));
            }
        }
    }
    if let Some((p, j, low)) = best {
        if p == 0 && !overrule {
            return Err(NO_DIRECTION.to_string());
        }
        w[j] = if low { -1.0 } else { 1.0 };
    }
    Ok(w)
}

const COLUMN_BYTES_PER_ROW: u64 = 48;
const TASK_BYTES_PER_EXAMPLE: u64 = 64;
const TASK_BYTES_PER_ROW: u64 = 40;

fn score_rows(pin: &Pin, norm: &Norm, rows: &[usize], w: &[f64]) -> Vec<f64> {
    let mut out = vec![0.0f64; rows.len()];
    fill(&pin.rows, norm, rows, |p, row| out[p] = dot(row, w));
    out
}

fn dot(row: &[f64], w: &[f64]) -> f64 {
    let nf = w.len() - 1;
    let mut s = w[nf];
    for (&x, &wj) in row[..nf].iter().zip(&w[..nf]).rev() {
        s += x * wj;
    }
    s
}

fn descending(pin: &Pin, fold: &Fold, s: &[f64], a: usize, b: usize) -> Ordering {
    let (ma, mb) = (&pin.meta[fold.train[a]], &pin.meta[fold.train[b]]);
    (s[b] + 0.0)
        .total_cmp(&(s[a] + 0.0))
        .then(mb.scan.cmp(&ma.scan))
        .then(mb.exp_mass.total_cmp(&ma.exp_mass))
        .then(a.cmp(&b))
}

fn train_set(pin: &Pin, norm: &Norm, fold: &Fold, s: &[f64], fdr: f64) -> Result<TrainSet, String> {
    let q = tdc_q(s, &fold.train_decoy, true);
    let (mut pos, mut neg) = training_sets(&q, &fold.train_decoy, fdr);
    neg.sort_by(|&a, &b| descending(pin, fold, s, a, b));
    pos.sort_by(|&a, &b| descending(pin, fold, s, a, b));
    let global: Vec<usize> = neg.iter().chain(&pos).map(|&l| fold.train[l]).collect();
    Ok(TrainSet {
        x: Rows::gather_within(&pin.rows, norm, &global, &pin.budget)?,
        idx_neg: (0..neg.len()).collect(),
        idx_pos: (neg.len()..neg.len() + pos.len()).collect(),
    })
}

fn positives(s: &[f64], decoy: &[bool], fdr: f64, skip_decoys_plus_one: bool) -> usize {
    count_at(&tdc_q(s, decoy, skip_decoys_plus_one), decoy, fdr)
}

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
    let nf = pin.rows.nfeat();
    let fold = folds_by_spectrum(&pin.meta, FOLDS, cfg.seed);
    let is_decoy: Vec<bool> = pin.meta.iter().map(|m| m.label < 0).collect();
    let folds: Vec<Fold> = (0..FOLDS)
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
        .collect();
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

    let init: Vec<Vec<f64>> = folds
        .iter()
        .map(|fd| {
            init_direction(
                pin,
                norm,
                &fd.train,
                cfg.initial_selection_fdr,
                cfg.overrule,
            )
        })
        .collect::<Result<_, _>>()?;
    if cfg.verbose > 1 && pin_direction(pin).is_none() {
        for (f, (fd, w)) in folds.iter().zip(&init).enumerate() {
            let s = score(&pin.rows, &fd.train, w, Some(norm));
            let p = positives(&s, &fd.train_decoy, cfg.initial_selection_fdr, true);
            let j = w[..nf].iter().position(|&v| v != 0.0).map_or(0, |j| j + 1);
            eprintln!(
                "Split {}:\tSelected feature {j} as initial direction. Could separate {p} training set positives with q<{} in that direction.",
                f + 1,
                cfg.initial_selection_fdr
            );
        }
    }
    let init_positives: usize = folds
        .iter()
        .zip(&init)
        .map(|(fd, w)| {
            let s = score(&pin.rows, &fd.test, w, Some(norm));
            positives(&s, &fd.test_decoy, cfg.test_fdr, false)
        })
        .sum();
    if cfg.verbose > 0 {
        eprintln!(
            "Found {init_positives} test set positives with q<{} in initial direction",
            cfg.test_fdr
        );
    }
    let mut w = init.clone();
    let mut train_scores: Vec<Vec<f64>> = folds
        .iter()
        .zip(&w)
        .map(|(fd, w)| score(&pin.rows, &fd.train, w, Some(norm)))
        .collect();

    for it in 0..cfg.iterations {
        let fdr = if it == 0 {
            cfg.initial_selection_fdr
        } else {
            cfg.selection_fdr
        };
        let sets: Vec<TrainSet> = folds
            .iter()
            .zip(&train_scores)
            .map(|(fd, s)| train_set(pin, norm, fd, s, fdr))
            .collect::<Result<_, _>>()?;
        if cfg.verbose > 1 {
            eprintln!("Iteration {}:", it + 1);
        }
        if cfg.verbose > 2 {
            for (f, t) in sets.iter().enumerate() {
                eprintln!(
                    "Split {}: Training with {} positives and {} negatives",
                    f + 1,
                    t.idx_pos.len(),
                    t.idx_neg.len()
                );
            }
        }
        let task_bytes: Vec<u64> = pairs
            .iter()
            .map(|&(f, _, _)| {
                let m = (sets[f].idx_neg.len() + sets[f].idx_pos.len()) as u64;
                m * TASK_BYTES_PER_EXAMPLE + folds[f].train.len() as u64 * TASK_BYTES_PER_ROW
            })
            .collect();
        let mut fits: Vec<Option<(Vec<f64>, usize)>> = Vec::with_capacity(pairs.len());
        while fits.len() < pairs.len() {
            let start = fits.len();
            let mut held = Vec::new();
            while start + held.len() < pairs.len() {
                match Reservation::take(&pin.budget, task_bytes[start + held.len()]) {
                    Some(r) => held.push(r),
                    None => break,
                }
            }
            let end = (start + held.len()).max(start + 1);
            fits.par_extend(pairs[start..end].par_iter().map(|&(f, cpos, cfrac)| {
                let t = &sets[f];
                if t.idx_neg.is_empty() && t.idx_pos.is_empty() {
                    return None;
                }
                let ww = l2_svm_mfn(&t.x, &t.idx_pos, &t.idx_neg, cpos, cfrac * cpos, &opts);
                let s = score_rows(pin, norm, &folds[f].train, &ww);
                let tp = positives(&s, &folds[f].train_decoy, cfg.test_fdr, true);
                Some((ww, tp))
            }));
        }
        let mut estimated = 0usize;
        for f in 0..sets.len() {
            let mut best = 0usize;
            let mut chosen = None;
            for (k, fit) in fits[f * per_fold..(f + 1) * per_fold].iter().enumerate() {
                if let Some((ww, tp)) = fit {
                    if *tp >= best {
                        best = *tp;
                        chosen = Some(k);
                        w[f].clone_from(ww);
                    }
                }
            }
            estimated += best;
            if cfg.verbose > 2 {
                if let Some(k) = chosen {
                    let (_, cpos, cfrac) = pairs[f * per_fold + k];
                    eprintln!(
                        "Split {}: selected Cpos={cpos} Cneg={} with {best} training set positives",
                        f + 1,
                        cfrac * cpos
                    );
                }
            }
            train_scores[f] = score_rows(pin, norm, &folds[f].train, &w[f]);
        }
        if cfg.verbose > 1 {
            eprintln!("Estimated {estimated} PSMs with q<{}", cfg.test_fdr);
        }
    }

    let mut test: Vec<Vec<f64>> = folds
        .iter()
        .zip(&w)
        .map(|(fd, w)| score(&pin.rows, &fd.test, w, Some(norm)))
        .collect();
    let found: usize = folds
        .iter()
        .zip(&test)
        .map(|(fd, s)| positives(s, &fd.test_decoy, cfg.test_fdr, false))
        .sum();
    if cfg.verbose > 0 {
        eprintln!("Found {found} test set PSMs with q<{}.", cfg.test_fdr);
    }
    if found == 0 || init_positives > found {
        if !cfg.overrule {
            let from_pin = pin_direction(pin).is_some();
            for (wf, i) in w.iter_mut().zip(&init) {
                if from_pin {
                    wf.clone_from(i);
                } else {
                    wf[..nf].copy_from_slice(&i[..nf]);
                }
            }
        }
        test = folds
            .iter()
            .map(|fd| score(&pin.rows, &fd.test, &w[0], Some(norm)))
            .collect();
    }

    let mut test_scores = vec![0.0f64; n];
    for (fd, s) in folds.iter().zip(&test) {
        for (&i, &v) in fd.test.iter().zip(s) {
            test_scores[i] = v;
        }
    }
    Ok(CvResult {
        w_raw: std::array::from_fn(|f| raw_weights(norm, &w[f])),
        w_norm: std::array::from_fn(|f| w[f].clone()),
        test_scores,
        fold,
    })
}
