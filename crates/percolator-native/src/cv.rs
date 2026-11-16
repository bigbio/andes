use crate::folds::folds_by_spectrum;
use crate::normalize::{raw_weights, score, Norm};
use crate::pin::Pin;
use crate::scores::{count_at, tdc_q, training_sets};
use crate::svm::{l2_svm_mfn, MfnOpts, Rows};
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
    perm: Vec<usize>,
    decoy: Vec<bool>,
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
    let x = Rows::gather(&pin.rows, norm, rows);
    let decoy: Vec<bool> = rows.iter().map(|&i| pin.meta[i].label < 0).collect();
    let counts: Vec<[usize; 2]> = (0..nf)
        .into_par_iter()
        .map(|j| {
            let high: Vec<f64> = (0..x.nrows).map(|r| x.row(r)[j]).collect();
            let low: Vec<f64> = high.iter().map(|&v| -v).collect();
            [
                count_at(&tdc_q(&low, &decoy, true), &decoy, fdr),
                count_at(&tdc_q(&high, &decoy, true), &decoy, fdr),
            ]
        })
        .collect();
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

fn train_set(pin: &Pin, norm: &Norm, fold: &Fold, s: &[f64], fdr: f64) -> TrainSet {
    let q = tdc_q(s, &fold.train_decoy, true);
    let (mut pos, mut neg) = training_sets(&q, &fold.train_decoy, fdr);
    neg.sort_by(|&a, &b| descending(pin, fold, s, a, b));
    pos.sort_by(|&a, &b| descending(pin, fold, s, a, b));
    let mut used = vec![false; s.len()];
    for &i in neg.iter().chain(&pos) {
        used[i] = true;
    }
    let mut perm = Vec::with_capacity(s.len());
    perm.extend_from_slice(&neg);
    perm.extend_from_slice(&pos);
    perm.extend((0..s.len()).filter(|&i| !used[i]));
    let global: Vec<usize> = perm.iter().map(|&l| fold.train[l]).collect();
    TrainSet {
        x: Rows::gather(&pin.rows, norm, &global),
        decoy: perm.iter().map(|&l| fold.train_decoy[l]).collect(),
        idx_neg: (0..neg.len()).collect(),
        idx_pos: (neg.len()..neg.len() + pos.len()).collect(),
        perm,
    }
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
    let init_positives: usize = folds
        .iter()
        .zip(&init)
        .map(|(fd, w)| {
            let s = score(&pin.rows, &fd.test, w, Some(norm));
            positives(&s, &fd.test_decoy, cfg.test_fdr, false)
        })
        .sum();
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
            .collect();
        let fits: Vec<Option<(Vec<f64>, usize)>> = pairs
            .par_iter()
            .map(|&(f, cpos, cfrac)| {
                let t = &sets[f];
                if t.idx_neg.is_empty() && t.idx_pos.is_empty() {
                    return None;
                }
                let ww = l2_svm_mfn(&t.x, &t.idx_pos, &t.idx_neg, cpos, cfrac * cpos, &opts);
                let s: Vec<f64> = (0..t.x.nrows).map(|r| dot(t.x.row(r), &ww)).collect();
                let tp = positives(&s, &t.decoy, cfg.test_fdr, true);
                Some((ww, tp))
            })
            .collect();
        for (f, t) in sets.iter().enumerate() {
            let mut best = 0usize;
            for (ww, tp) in fits[f * per_fold..(f + 1) * per_fold].iter().flatten() {
                if *tp >= best {
                    best = *tp;
                    w[f].clone_from(ww);
                }
            }
            for (r, &local) in t.perm.iter().enumerate() {
                train_scores[f][local] = dot(t.x.row(r), &w[f]);
            }
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
