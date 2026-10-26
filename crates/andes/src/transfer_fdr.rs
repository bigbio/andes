//! Transfer FDR of match-between-runs: a cross-validated GBDT separates
//! transfer targets from their decoys, and the K-decoy estimator turns the
//! held-out scores into q-values.
//!
//! The transfer q-value is an internal quant-quality score that decides which
//! transferred MS1 quantities are reported; it is not an identification FDR.
//! The identified rows' quality q-value (`quant_q`) is the same estimate over
//! the identified targets and their decoys, fitted separately with
//! [`IDENTIFIED_FEATURES`].
//!
//! Folds: a hash of the precursor id (splitmix64, seeded) gives each precursor
//! one of three folds, so its targets and decoys in every run share a fold.
//! For each fold a classifier (target = 1, decoy = 0, validation groups = the
//! precursor) is trained on the other two folds and scores the held-out fold.
//! Scores are the raw boosted margins and are compared within a fold only.
//!
//! Bagging: a fold's classifier can be the mean of several members, each
//! trained on its own row set (a [`TrainingSample`] of its own seed) with its
//! own trainer seed ([`member_seed`]); a row's score is the mean of the
//! members' margins. Every member of a fold is trained on the other folds
//! only, so the mean is held out like a single model's. The members average
//! out the variance of one boosted model trained on a sample of the rows,
//! which at hundreds of runs moves the scores of the transfers near the
//! threshold more than their evidence does.
//!
//! q-values, per fold: rank the fold's scored rows by score; at threshold `c`,
//! `FDR(c) = (D(c) + 1) / (K * max(T(c), 1))`, capped at 1, with `D`/`T` the
//! decoys/targets at or above `c` and `K` the fold's decoy-to-target row
//! ratio; q is the running minimum of `FDR` from the lowest score up. The `+1`
//! is the finite-sample correction: a target above every decoy is not given
//! q = 0, and no target reaches q <= 0.01 unless at least `1 / (0.01 K)`
//! targets score above every decoy. The accepted set at
//! any threshold is the union of the folds' accepted sets, so its estimated FDR
//! is at most the threshold. The same estimate is made with one decoy family
//! at a time (`K` = that family's ratio).
//!
//! Identity decoys (family [`IDENTITY`]): the transfers of decoy precursors
//! that passed the library's donor rules with decoy PSMs. The mass- and
//! RT-shift decoys model a transfer of the right peptide onto the wrong peak;
//! a false donor identification sits on a real MS1 feature that recurs in
//! every run, which they cannot model. Under target-decoy competition each
//! identity-decoy entry stands for one false target entry, so the accepted
//! identity-decoy transfers estimate the accepted target transfers of false
//! donors, and the combined estimate adds that share:
//! `FDR(c) = (D(c) + 1) / (K * max(T(c), 1)) + D_id(c) / max(T(c), 1)`,
//! with `D` and `K` over the wrong-peak families only. The two error classes
//! overlap (a false donor's transfer can also land on a wrong peak), so the
//! sum is conservative. Without identity-decoy rows it is the wrong-peak
//! estimate. Identity decoys are negatives of the classifier's training like
//! the other decoys (held out by the same folds), so it learns what sets a
//! false donor's transfer apart where it can.
//!
//! `K` is counted over the rows passed in, which are the rows that have
//! features: a target or decoy without evidence takes no part. When decoys are
//! missing (no RT-shift decoy near a run edge, a window too short to fit), `K`
//! falls below the nominal number of decoys per target and each remaining
//! decoy counts for more targets, rather than the missing decoys being counted
//! as absent false hits; this errs on the conservative side.
//!
//! A row gets no score when its fold's training set is single-class or fails
//! the trainer's quality gate, and no q-value when it has no score or its fold
//! (restricted to the estimate's families) lacks targets or decoys.

use model_train::gbdt::train::{train_gbdt, Dataset, TrainParams};
use rayon::prelude::*;
use scoring_crate::gbdt_eval::GbdtPeakModel;

/// Cross-validation folds.
pub const N_FOLDS: usize = 3;

/// Feature names of the transfer classifier, in [`FdrRow::x`] order.
pub const TRANSFER_FEATURES: &[&str] = &[
    "llr",
    "rt_z",
    "abs_rt_z",
    "width_ratio",
    "envelope_cos",
    "coherence",
    "m1_ratio",
    "abs_mass_err_ppm",
    "n_points",
    "n_isotopes_used",
    "log10_amplitude_snr",
    "n_donor_runs",
    "library_q",
    "neighbour_donor_rate",
    "shape_log_width",
    "shape_tailing_diff",
    "shape_iso_cos",
];

/// Feature names of the identified-row quality classifier: the template-fit
/// evidence of [`TRANSFER_FEATURES`] without `n_donor_runs`, `library_q` and
/// `neighbour_donor_rate`, which are the same for an identified target and
/// its decoys, and without the donor shape comparison (an identified row is
/// itself a donor).
pub const IDENTIFIED_FEATURES: &[&str] = &[
    "llr",
    "rt_z",
    "abs_rt_z",
    "width_ratio",
    "envelope_cos",
    "coherence",
    "m1_ratio",
    "abs_mass_err_ppm",
    "n_points",
    "n_isotopes_used",
    "log10_amplitude_snr",
];

/// [`FdrRow::family`] of a target.
pub const TARGET: u8 = 0;
/// [`FdrRow::family`] of a mass-shift decoy.
pub const MASS_SHIFT: u8 = 1;
/// [`FdrRow::family`] of an RT-shift decoy.
pub const RT_SHIFT: u8 = 2;
/// [`FdrRow::family`] of an identity decoy: the target probe of a decoy
/// precursor's transfer.
pub const IDENTITY: u8 = 3;

/// One target or decoy to be scored.
#[derive(Debug, Clone, PartialEq)]
pub struct FdrRow {
    pub precursor: u32,
    pub is_decoy: bool,
    /// [`TARGET`], [`MASS_SHIFT`], [`RT_SHIFT`] or [`IDENTITY`].
    pub family: u8,
    pub x: Vec<f32>,
}

/// Per-row output of [`transfer_qvalues`], in input order.
#[cfg(test)]
#[derive(Debug, Clone, PartialEq)]
pub struct FdrResult {
    /// The wrong-peak families and the identity decoys.
    pub q: Vec<Option<f64>>,
    /// Both wrong-peak families, without the identity decoys.
    pub q_wrong_peak: Vec<Option<f64>>,
    /// Mass-shift decoys only; `None` for RT-shift decoy rows.
    pub q_mass_only: Vec<Option<f64>>,
    /// RT-shift decoys only; `None` for mass-shift decoy rows.
    pub q_rt_only: Vec<Option<f64>>,
    /// Held-out raw margin, comparable within a fold; NaN when unscored.
    pub score: Vec<f32>,
}

/// What the cross-validation could not score, for the caller to log.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct FdrCounts {
    /// One message per fold without a model.
    pub failed_folds: Vec<String>,
    /// Rows without a score.
    pub unscored: usize,
    /// Scored rows without a pooled q-value.
    pub no_q: usize,
}

/// The fold of a precursor (splitmix64 of its id).
pub fn fold_of(precursor: u32, seed: u64) -> usize {
    let mut z = (seed ^ (precursor as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15))
        .wrapping_add(0x9e37_79b9_7f4a_7c15);
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^= z >> 31;
    (z % N_FOLDS as u64) as usize
}

/// Cross-validated q-values of transfer targets and decoys.
#[cfg(test)]
pub fn transfer_qvalues(rows: &[FdrRow], seed: u64) -> FdrResult {
    transfer_fdr(rows, seed).result
}

/// Trains the classifier of `fold` on the rows of the other folds.
fn train_fold(
    rows: &[FdrRow],
    folds: &[usize],
    nf: usize,
    fold: usize,
    seed: u64,
) -> Result<GbdtPeakModel, String> {
    let mut x = Vec::new();
    let mut y = Vec::new();
    let mut groups = Vec::new();
    for (r, &f) in rows.iter().zip(folds) {
        if f == fold {
            continue;
        }
        x.extend_from_slice(&r.x);
        y.push(u8::from(!r.is_decoy));
        groups.push(r.precursor);
    }
    let positives = y.iter().filter(|&&v| v == 1).count();
    if positives == 0 || positives == y.len() {
        return Err(format!(
            "fold {fold}: training rows are single-class ({positives} targets of {})",
            y.len()
        ));
    }
    let ds = Dataset {
        x,
        y,
        groups,
        n_features: nf,
        feature_names: scoring_crate::gbdt_eval::indexed_feature_names(nf as u32),
    };
    train_gbdt(
        &ds,
        &TrainParams::default(),
        seed.wrapping_add(fold as u64 + 1),
    )
    .map_err(|e| format!("fold {fold}: {e:?}"))
}

/// The fold models of one classifier: fold `f`'s members are trained on the
/// rows of the other folds and score the rows of fold `f`.
pub struct FoldModels {
    seed: u64,
    /// Per fold, the members that trained (empty when none did).
    models: Vec<Vec<GbdtPeakModel>>,
    /// One message per fold without a model.
    pub failed_folds: Vec<String>,
}

/// The seed of bagging member `m` (its training sample and its trainer):
/// `seed` itself for the first member, so that one member is the classifier
/// of a single model.
pub static SAMPLE_SEED: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(u64::MAX);
pub static TRAIN_SEED: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(u64::MAX);
fn ovr(a: &std::sync::atomic::AtomicU64, seed: u64) -> u64 {
    let v = a.load(std::sync::atomic::Ordering::Relaxed);
    if v == u64::MAX { seed } else { v }
}
pub fn sample_seed(seed: u64, m: usize) -> u64 {
    member_seed(ovr(&SAMPLE_SEED, seed), m)
}
pub fn train_seed(seed: u64, m: usize) -> u64 {
    member_seed(ovr(&TRAIN_SEED, seed), m)
}
pub fn member_seed(seed: u64, m: usize) -> u64 {
    if m == 0 {
        seed
    } else {
        sample_hash(m as u64, seed ^ 0x6d65_6d62_6572_7321)
    }
}

/// Trains the fold models on `rows` (all the rows, or a sample of whole
/// precursors: see [`TrainingSample`]).
///
/// # Panics
/// When the rows do not all have the same number of features.
pub fn train_folds(rows: &[FdrRow], seed: u64) -> FoldModels {
    train_bagged_folds(&[rows], seed)
}

/// Trains the fold models of a bagged classifier: member `m` is trained on
/// `sets[m]` with the trainer seed of [`member_seed`]`(seed, m)`; the folds
/// are those of `seed` for every member. A fold whose members all fail has
/// no model; a fold with some members left averages those (the others are
/// logged).
///
/// # Panics
/// When `sets` is empty or the rows do not all have the same number of
/// features.
pub fn train_bagged_folds(sets: &[&[FdrRow]], seed: u64) -> FoldModels {
    assert!(!sets.is_empty(), "a classifier needs at least one member");
    let nf = sets
        .iter()
        .find_map(|rows| rows.first())
        .map_or(0, |r| r.x.len());
    assert!(
        sets.iter()
            .flat_map(|rows| rows.iter())
            .all(|r| r.x.len() == nf),
        "transfer FDR rows have different feature counts"
    );
    let folds: Vec<Vec<usize>> = sets
        .iter()
        .map(|rows| rows.iter().map(|r| fold_of(r.precursor, seed)).collect())
        .collect();
    let jobs: Vec<(usize, usize)> = (0..N_FOLDS)
        .flat_map(|f| (0..sets.len()).map(move |m| (f, m)))
        .collect();
    let fitted: Vec<Result<GbdtPeakModel, String>> = jobs
        .par_iter()
        .map(|&(f, m)| train_fold(sets[m], &folds[m], nf, f, train_seed(seed, m)))
        .collect();
    let mut models: Vec<Vec<GbdtPeakModel>> = (0..N_FOLDS).map(|_| Vec::new()).collect();
    let mut errors: Vec<Vec<String>> = vec![Vec::new(); N_FOLDS];
    for (&(f, m), fit) in jobs.iter().zip(fitted) {
        match fit {
            Ok(model) => models[f].push(model),
            Err(e) if sets.len() == 1 => errors[f].push(e),
            Err(e) => errors[f].push(format!("member {m}: {e}")),
        }
    }
    let mut failed_folds = Vec::new();
    for (f, e) in errors.into_iter().enumerate() {
        if e.is_empty() {
            continue;
        }
        if models[f].is_empty() {
            failed_folds.push(e.join("; "));
        } else {
            eprintln!(
                "[mbr] fold {f}: {} of {} classifier members failed, the fold averages the others ({})",
                e.len(),
                sets.len(),
                e.join("; ")
            );
        }
    }
    FoldModels {
        seed,
        models,
        failed_folds,
    }
}

/// Rows per call of the batched GBDT evaluation; the batches of one fold
/// are scored in parallel (each row's margin does not depend on the batch).
const PREDICT_BATCH: usize = 4096;

impl FoldModels {
    /// The held-out raw margins of `xs` (rows of precursor `precursors[i]`)
    /// into `out`: `None` for a row whose fold has no model. A bagged fold
    /// scores the mean of its members' margins (summed in member order).
    pub fn score(&self, precursors: &[u32], xs: &[&[f32]], out: &mut Vec<Option<f32>>) {
        out.clear();
        out.resize(xs.len(), None);
        for (fold, members) in self.models.iter().enumerate() {
            if members.is_empty() {
                continue;
            }
            let at: Vec<usize> = (0..xs.len())
                .filter(|&i| fold_of(precursors[i], self.seed) == fold)
                .collect();
            if at.is_empty() {
                continue;
            }
            let batch: Vec<&[f32]> = at.iter().map(|&i| xs[i]).collect();
            let mut s = vec![0.0f32; at.len()];
            batch
                .par_chunks(PREDICT_BATCH)
                .zip(s.par_chunks_mut(PREDICT_BATCH))
                .for_each(|(b, o)| predict_mean(members, b, o));
            for (&i, v) in at.iter().zip(s) {
                out[i] = Some(v);
            }
        }
    }

    /// Members per fold (the most of any fold).
    pub fn members(&self) -> usize {
        self.models.iter().map(Vec::len).max().unwrap_or(0)
    }
}

/// The mean margin of `members` on `rows` into `out`; one member's margin
/// as it is.
fn predict_mean(members: &[GbdtPeakModel], rows: &[&[f32]], out: &mut [f32]) {
    let Some((first, rest)) = members.split_first() else {
        return;
    };
    first.predict_value_batch(rows, out);
    if rest.is_empty() {
        return;
    }
    let mut tmp = vec![0.0f32; out.len()];
    for m in rest {
        m.predict_value_batch(rows, &mut tmp);
        for (o, v) in out.iter_mut().zip(&tmp) {
            *o += v;
        }
    }
    let n = members.len() as f32;
    for o in out.iter_mut() {
        *o /= n;
    }
}

/// Items (training rows, normalisation cells) collected over a stream, at
/// most `cap` of them, by sampling units: every item carries a unit id (a
/// precursor, or a precursor in one run), a unit's items are kept when the
/// top `level` bits of a hash of its id are zero, and `level` rises by one
/// each time the kept items exceed `cap`. The kept set is thus every item of
/// the units at the final level, whatever the order the items came in, and
/// all items when they never exceed `cap`. Items carry a sort key so that
/// [`TrainingSample::into_items`] can return them in a fixed order.
pub struct TrainingSample<K, T = FdrRow> {
    cap: usize,
    seed: u64,
    level: u32,
    seen: usize,
    items: Vec<(K, u64, T)>,
}

impl<K: Ord, T> TrainingSample<K, T> {
    pub fn new(cap: usize, seed: u64) -> Self {
        Self {
            cap: cap.max(1),
            seed,
            level: 0,
            seen: 0,
            items: Vec::new(),
        }
    }

    fn kept(seed: u64, unit: u64, level: u32) -> bool {
        level == 0 || sample_hash(unit, seed) >> (64 - level) == 0
    }

    /// Whether an item of sampling unit `unit` offered now would be kept
    /// (so that an item need only be built when it is).
    pub fn wants(&self, unit: u64) -> bool {
        Self::kept(self.seed, unit, self.level)
    }

    /// An item of unit `unit` that [`TrainingSample::wants`] declined:
    /// counted as offered.
    pub fn skip(&mut self) {
        self.seen += 1;
    }

    /// Offers `item` of sampling unit `unit`.
    pub fn push(&mut self, key: K, unit: u64, item: T) {
        self.seen += 1;
        if !Self::kept(self.seed, unit, self.level) {
            return;
        }
        self.items.push((key, unit, item));
        while self.items.len() > self.cap && self.level < 63 {
            self.level += 1;
            let (seed, level) = (self.seed, self.level);
            self.items.retain(|(_, u, _)| Self::kept(seed, *u, level));
        }
    }

    /// Items offered, items kept and the share of units kept.
    pub fn stats(&self) -> (usize, usize, f64) {
        (self.seen, self.items.len(), 0.5f64.powi(self.level as i32))
    }

    /// The kept items, sorted by key (a stable sort: items of equal keys keep
    /// their order of arrival).
    pub fn into_items(self) -> Vec<T> {
        let mut items = self.items;
        items.sort_by(|a, b| a.0.cmp(&b.0));
        items.into_iter().map(|(_, _, t)| t).collect()
    }
}

/// The sampling unit of a training row: its precursor in its run, so that a
/// sample keeps every precursor (in some runs) and a target stays with its
/// decoys. Folds are by precursor, so a held-out precursor is never trained
/// on.
pub fn row_unit(precursor: u32, run: u16) -> u64 {
    ((precursor as u64) << 16) | run as u64
}

/// The unit hash of [`TrainingSample`] (splitmix64 of the id, on a seed
/// apart from the folds').
pub fn sample_hash(unit: u64, seed: u64) -> u64 {
    let mut z = (seed ^ 0x5bd1_e995_7a4f_3c21 ^ unit.wrapping_mul(0x9e37_79b9_7f4a_7c15))
        .wrapping_add(0x9e37_79b9_7f4a_7c15);
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

/// The q-value estimates of a classifier's rows.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Estimate {
    /// The wrong-peak families and the identity decoys.
    Pooled = 0,
    /// Both wrong-peak families.
    WrongPeak = 1,
    /// Mass-shift decoys only.
    MassOnly = 2,
    /// RT-shift decoys only.
    RtOnly = 3,
}

const ESTIMATES: [Estimate; 4] = [
    Estimate::Pooled,
    Estimate::WrongPeak,
    Estimate::MassOnly,
    Estimate::RtOnly,
];

impl Estimate {
    /// Whether rows of `family` take part in this estimate.
    fn keeps(self, family: u8) -> bool {
        match self {
            Estimate::Pooled => true,
            Estimate::WrongPeak => family != IDENTITY,
            Estimate::MassOnly => family == TARGET || family == MASS_SHIFT,
            Estimate::RtOnly => family == TARGET || family == RT_SHIFT,
        }
    }
}

/// `0.0` for `-0.0`: scores equal under `==` are one value.
fn canonical(s: f32) -> f32 {
    if s == 0.0 {
        0.0
    } else {
        s
    }
}

/// A fold's q-value as a step function of the score: the ascending scores at
/// which q drops, with the q from each up to the next.
#[derive(Debug, Clone, Default)]
pub(crate) struct Steps {
    at: Vec<f32>,
    q: Vec<f64>,
}

impl Steps {
    /// The q of a row of the estimate scoring `s` (a score that took part).
    fn of_member(&self, s: f32) -> f64 {
        let s = canonical(s);
        let k = self.at.partition_point(|&a| a.total_cmp(&s).is_le());
        self.q[k.max(1) - 1]
    }

    /// The q of the highest fitted score at or below `s`, or of the lowest
    /// when `s` is below them all ([`TransferFdr::pseudo_q`]).
    fn at_or_below(&self, s: f32) -> f64 {
        let k = self.at.partition_point(|&a| a <= s);
        self.q[k.max(1) - 1]
    }
}

/// The held-out scores of one classifier's rows, per fold and family, kept
/// to turn them into q-values: 4 bytes per scored row held.
///
/// A mergeable summary: [`FoldScores::merge`] of the scores of disjoint row
/// sets gives the scores of their union, whatever the order, so the rows can
/// be scored in parallel pieces. The scores held can be taken out
/// ([`FoldScores::take_held`], to be sorted and written elsewhere) and given
/// back to [`FoldScores::finish_with`] as ascending streams: the step
/// functions only need each fold's scores in ascending order.
pub struct FoldScores {
    seed: u64,
    scores: Vec<[Vec<f32>; 4]>,
    /// Rows with a finite score, per fold and family, held or taken out.
    counts: Vec<[usize; 4]>,
    members: Vec<usize>,
    unscored: usize,
}

impl FoldScores {
    pub fn new(seed: u64) -> Self {
        Self {
            seed,
            scores: (0..N_FOLDS).map(|_| Default::default()).collect(),
            counts: vec![[0; 4]; N_FOLDS],
            members: vec![0; N_FOLDS],
            unscored: 0,
        }
    }

    /// One row of `precursor` and `family` with its held-out score (`None`
    /// when its fold has no model).
    pub fn add(&mut self, precursor: u32, family: u8, score: Option<f32>) {
        let fold = fold_of(precursor, self.seed);
        match score {
            None => self.unscored += 1,
            Some(s) => {
                self.members[fold] += 1;
                if s.is_finite() {
                    let f = family.min(3) as usize;
                    self.scores[fold][f].push(canonical(s));
                    self.counts[fold][f] += 1;
                }
            }
        }
    }

    /// Adds the rows of `other` (scored with the same seed's folds).
    ///
    /// # Panics
    /// When `other` has another seed.
    pub fn merge(&mut self, other: FoldScores) {
        assert_eq!(self.seed, other.seed, "fold scores of different seeds");
        for (fold, by_family) in other.scores.into_iter().enumerate() {
            for (f, v) in by_family.into_iter().enumerate() {
                let mine = &mut self.scores[fold][f];
                if mine.is_empty() {
                    *mine = v;
                } else {
                    mine.extend(v);
                }
                self.counts[fold][f] += other.counts[fold][f];
            }
            self.members[fold] += other.members[fold];
        }
        self.unscored += other.unscored;
    }

    /// Scores held in memory.
    pub fn held(&self) -> usize {
        self.scores.iter().flatten().map(Vec::len).sum()
    }

    /// The scores held, per fold and family, in arrival order, taken out (the
    /// counts stay).
    pub fn take_held(&mut self) -> Vec<[Vec<f32>; 4]> {
        self.scores
            .iter_mut()
            .map(|by_family| by_family.each_mut().map(std::mem::take))
            .collect()
    }

    /// The q-value step functions; `failed_folds` are the folds without a
    /// model (for the counts).
    #[cfg(test)]
    pub fn finish(self, failed_folds: &[String]) -> QValues {
        self.finish_with(failed_folds, None, |_, _| Vec::new())
    }

    /// [`FoldScores::finish`] with, for each fold and family, the scores
    /// taken out by [`FoldScores::take_held`] given back by `taken(fold,
    /// family)` as ascending streams (any number of them). The result is the
    /// one of every score held in memory. With `thresholds`, only the pooled
    /// estimate keeps its whole step function; the others keep the steps
    /// where q first reaches each threshold, which answer
    /// [`QValues::accepts`] at those thresholds exactly (about a quarter of
    /// the memory).
    pub fn finish_with<'a>(
        self,
        failed_folds: &[String],
        thresholds: Option<&[f64]>,
        mut taken: impl FnMut(usize, usize) -> Vec<Box<dyn Iterator<Item = f32> + 'a>>,
    ) -> QValues {
        let mut steps: Vec<[Option<Steps>; 4]> = Vec::with_capacity(N_FOLDS);
        let mut no_q = 0usize;
        for (fold, by_family) in self.scores.into_iter().enumerate() {
            let mut k = 0usize;
            let streams: [SortedMerge<'a>; 4] = by_family.map(|mut v| {
                v.sort_unstable_by(f32::total_cmp);
                let mut runs = taken(fold, k);
                k += 1;
                runs.push(Box::new(v.into_iter()));
                SortedMerge::new(runs)
            });
            let totals = self.counts[fold];
            let mut builder = StepsBuilder::new(totals, thresholds);
            walk_scores(streams, |c, at| builder.visit(c, at));
            let fold_steps = builder.finish();
            let with_q = if fold_steps[Estimate::Pooled as usize].is_some() {
                totals.iter().sum()
            } else {
                0
            };
            no_q += self.members[fold] - with_q;
            steps.push(fold_steps);
        }
        QValues {
            seed: self.seed,
            steps,
            whole: ESTIMATES.map(|e| thresholds.is_none() || e == Estimate::Pooled),
            thresholds: thresholds.map(<[f64]>::to_vec).unwrap_or_default(),
            counts: FdrCounts {
                failed_folds: failed_folds.to_vec(),
                unscored: self.unscored,
                no_q,
            },
        }
    }
}

/// Ascending streams of scores merged into one ascending stream.
struct SortedMerge<'a> {
    runs: Vec<Box<dyn Iterator<Item = f32> + 'a>>,
    heads: Vec<Option<f32>>,
}

impl<'a> SortedMerge<'a> {
    fn new(mut runs: Vec<Box<dyn Iterator<Item = f32> + 'a>>) -> Self {
        let heads = runs.iter_mut().map(|r| r.next()).collect();
        Self { runs, heads }
    }

    /// The lowest score left.
    fn peek(&self) -> Option<f32> {
        self.heads.iter().flatten().copied().min_by(f32::total_cmp)
    }

    /// Takes every score equal to `c` (the lowest left) and counts them.
    fn take_equal(&mut self, c: f32) -> usize {
        let mut n = 0;
        for (head, run) in self.heads.iter_mut().zip(&mut self.runs) {
            while *head == Some(c) {
                n += 1;
                *head = run.next();
            }
        }
        n
    }
}

/// Calls `visit(c, at)` for each distinct score `c` of the four families'
/// ascending streams, lowest first, with `at[f]` the scores of family `f`
/// equal to `c`.
fn walk_scores(mut by_family: [SortedMerge<'_>; 4], mut visit: impl FnMut(f32, [usize; 4])) {
    loop {
        let low = by_family
            .iter()
            .filter_map(SortedMerge::peek)
            .min_by(f32::total_cmp);
        let Some(c) = low else { break };
        let mut at = [0usize; 4];
        for (n, s) in at.iter_mut().zip(&mut by_family) {
            *n = s.take_equal(c);
        }
        visit(c, at);
    }
}

/// The K-decoy step functions of one fold, every [`Estimate`] at once, built
/// from the fold's distinct scores visited in ascending order (see
/// [`kdecoy_steps`] for the estimate).
pub(crate) struct StepsBuilder {
    walks: [Option<KdecoyWalk>; 4],
}

/// The state of one estimate's walk up the scores.
struct KdecoyWalk {
    e: Estimate,
    k: f64,
    /// Targets, wrong-peak decoys and identity decoys at or above the
    /// current score.
    t: usize,
    d: usize,
    d_id: usize,
    running: f64,
    steps: Steps,
    /// Keep only the first step and those where q first reaches one of
    /// these; `None` keeps every step.
    only_at: Option<Vec<f64>>,
}

impl StepsBuilder {
    /// The walks over a fold whose finite scores number `totals[f]` per
    /// family. With `thresholds`, the estimates other than the pooled one
    /// keep only their first step and the steps where q first reaches each
    /// threshold.
    pub fn new(totals: [usize; 4], thresholds: Option<&[f64]>) -> Self {
        Self {
            walks: ESTIMATES.map(|e| {
                let total = |f: u8| if e.keeps(f) { totals[f as usize] } else { 0 };
                let targets = total(TARGET);
                let decoys = total(MASS_SHIFT) + total(RT_SHIFT);
                (targets > 0 && decoys > 0).then(|| KdecoyWalk {
                    e,
                    k: decoys as f64 / targets as f64,
                    t: targets,
                    d: decoys,
                    d_id: total(IDENTITY),
                    running: 1.0,
                    steps: Steps::default(),
                    only_at: thresholds
                        .filter(|_| e != Estimate::Pooled)
                        .map(<[f64]>::to_vec),
                })
            }),
        }
    }

    /// The next distinct score `c` (ascending, canonical) with `at[f]` rows of
    /// family `f` scoring it.
    pub(crate) fn visit(&mut self, c: f32, at: [usize; 4]) {
        for w in self.walks.iter_mut().flatten() {
            let kept = |f: usize| w.e.keeps(f as u8);
            if !(0..4).any(|f| kept(f) && at[f] > 0) {
                continue;
            }
            let tf = w.t.max(1) as f64;
            let v = ((w.d + 1) as f64 / (w.k * tf) + w.d_id as f64 / tf).min(1.0);
            if v < w.running || w.steps.at.is_empty() {
                let before = w.running;
                w.running = w.running.min(v);
                let kept = match &w.only_at {
                    None => true,
                    Some(ts) => {
                        w.steps.at.is_empty() || ts.iter().any(|&t| w.running <= t && before > t)
                    }
                };
                if kept {
                    w.steps.at.push(c);
                    w.steps.q.push(w.running);
                }
            }
            for (f, &n) in at.iter().enumerate() {
                if !kept(f) {
                    continue;
                }
                match f as u8 {
                    TARGET => w.t -= n,
                    IDENTITY => w.d_id -= n,
                    _ => w.d -= n,
                }
            }
        }
    }

    /// The step functions by [`Estimate`]; `None` for an estimate without
    /// targets or wrong-peak decoys.
    pub fn finish(self) -> [Option<Steps>; 4] {
        self.walks.map(|w| w.map(|w| w.steps))
    }
}

/// The K-decoy q-values (`FDR = (D + 1) / (K * max(T, 1)) + D_id / max(T, 1)`) over one fold's sorted scores by
/// family, as a step function; `None` without targets or wrong-peak decoys
/// among the estimate's rows. Walks the scores upwards: at each distinct
/// score `c`, `T`, `D` and `D_id` at or above it are the totals minus the
/// rows below, and q is the running minimum of the FDR from the bottom.
#[cfg(test)]
fn kdecoy_steps(by_family: &[Vec<f32>; 4], e: Estimate) -> Option<Steps> {
    let mut builder = StepsBuilder::new([0, 1, 2, 3].map(|f| by_family[f].len()), None);
    let streams = [0, 1, 2, 3].map(|f| {
        SortedMerge::new(vec![
            Box::new(by_family[f].iter().copied()) as Box<dyn Iterator<Item = f32>>
        ])
    });
    walk_scores(streams, |c, at| builder.visit(c, at));
    let [a, b, c, d] = builder.finish();
    [a, b, c, d].into_iter().nth(e as usize).flatten()
}

/// The q-values of one classifier's rows: per fold and [`Estimate`], a step
/// function of the held-out score.
pub struct QValues {
    seed: u64,
    steps: Vec<[Option<Steps>; 4]>,
    /// Per [`Estimate`], whether its step functions are whole (else they
    /// answer [`QValues::accepts`] at `thresholds` only).
    whole: [bool; 4],
    thresholds: Vec<f64>,
    pub counts: FdrCounts,
}

impl QValues {
    /// Steps of every step function, held in memory (12 bytes each).
    pub fn steps(&self) -> usize {
        self.steps
            .iter()
            .flatten()
            .flatten()
            .map(|s| s.at.len())
            .sum()
    }

    /// Steps of the pooled estimate's step functions.
    pub fn pooled_steps(&self) -> usize {
        self.steps
            .iter()
            .filter_map(|f| f[Estimate::Pooled as usize].as_ref())
            .map(|s| s.at.len())
            .sum()
    }

    /// The q-value under estimate `e` of a row of `precursor` and `family`
    /// scoring `score`; `None` when it has no score, its family takes no part
    /// in `e`, or its fold lacks targets or decoys.
    ///
    /// # Panics
    /// When the step functions of `e` were kept at thresholds only.
    pub fn q(&self, e: Estimate, precursor: u32, family: u8, score: Option<f32>) -> Option<f64> {
        assert!(
            self.whole[e as usize],
            "{e:?} q-values were kept at thresholds only"
        );
        self.step_q(e, precursor, family, score)
    }

    /// Whether the row's q-value under `e` (see [`QValues::q`]) is at most
    /// `t`; `false` without one.
    ///
    /// # Panics
    /// When the step functions of `e` were kept at thresholds only and `t`
    /// is not one of them.
    pub fn accepts(
        &self,
        e: Estimate,
        precursor: u32,
        family: u8,
        score: Option<f32>,
        t: f64,
    ) -> bool {
        assert!(
            self.whole[e as usize] || self.thresholds.contains(&t),
            "{e:?} q-values were not kept at threshold {t}"
        );
        self.step_q(e, precursor, family, score)
            .is_some_and(|q| q <= t)
    }

    fn step_q(&self, e: Estimate, precursor: u32, family: u8, score: Option<f32>) -> Option<f64> {
        let s = score.filter(|s| s.is_finite())?;
        if !e.keeps(family) {
            return None;
        }
        let steps = self.steps[fold_of(precursor, self.seed)][e as usize].as_ref()?;
        Some(steps.of_member(s))
    }

    /// The pooled q of a row outside the fit (a pseudo-transfer) scoring
    /// `score` under its fold's model: the q of the fold's scored row with the
    /// highest score at or below it, or of the fold's lowest row when it
    /// scores below all of them. `None` without a score or without q-values
    /// in the fold.
    pub fn pseudo_q(&self, precursor: u32, score: Option<f32>) -> Option<f64> {
        let s = score?;
        let steps =
            self.steps[fold_of(precursor, self.seed)][Estimate::Pooled as usize].as_ref()?;
        Some(steps.at_or_below(s))
    }
}

/// The fitted fold models with their q-values over a slice of rows, kept to
/// score rows that took no part in the fit.
#[cfg(test)]
pub struct TransferFdr {
    pub result: FdrResult,
    pub counts: FdrCounts,
    models: FoldModels,
    q: QValues,
}

/// [`transfer_qvalues`] keeping the fold models: [`train_folds`] on every
/// row, then each row scored by its fold's model and given its q-values.
///
/// # Panics
/// When the rows do not all have the same number of features.
#[cfg(test)]
pub fn transfer_fdr(rows: &[FdrRow], seed: u64) -> TransferFdr {
    let models = train_folds(rows, seed);
    let precursors: Vec<u32> = rows.iter().map(|r| r.precursor).collect();
    let xs: Vec<&[f32]> = rows.iter().map(|r| r.x.as_slice()).collect();
    let mut scores = Vec::new();
    models.score(&precursors, &xs, &mut scores);
    let mut fs = FoldScores::new(seed);
    for (r, s) in rows.iter().zip(&scores) {
        fs.add(r.precursor, r.family, *s);
    }
    let q = fs.finish(&models.failed_folds);
    let each = |e: Estimate| -> Vec<Option<f64>> {
        rows.iter()
            .zip(&scores)
            .map(|(r, s)| q.q(e, r.precursor, r.family, *s))
            .collect()
    };
    let result = FdrResult {
        q: each(Estimate::Pooled),
        q_wrong_peak: each(Estimate::WrongPeak),
        q_mass_only: each(Estimate::MassOnly),
        q_rt_only: each(Estimate::RtOnly),
        score: scores.iter().map(|s| s.unwrap_or(f32::NAN)).collect(),
    };
    TransferFdr {
        result,
        counts: q.counts.clone(),
        models,
        q,
    }
}

#[cfg(test)]
impl TransferFdr {
    /// The q-value of a row outside the fit (a pseudo-transfer), scored by
    /// the model of its precursor's fold: see [`QValues::pseudo_q`].
    pub fn pseudo_q(&self, precursor: u32, x: &[f32]) -> Option<f64> {
        let mut s = Vec::new();
        self.models.score(&[precursor], &[x], &mut s);
        self.q.pseudo_q(precursor, s[0])
    }
}

/// K-decoy q-values of `(row, score, family)` into `q[row]`, with `K` the
/// wrong-peak-decoy-to-target ratio of `items` and
/// `FDR = (D + 1) / (K * max(T, 1)) + D_id / max(T, 1)` (`D` wrong-peak
/// decoys, `D_id` identity decoys at or above the score). Leaves `q`
/// untouched when `items` has no targets or no wrong-peak decoys. The
/// reference of [`kdecoy_steps`].
#[cfg(test)]
fn kdecoy_q(items: &[(usize, f32, u8)], q: &mut [Option<f64>]) {
    let targets = items.iter().filter(|it| it.2 == TARGET).count();
    let decoys = items
        .iter()
        .filter(|it| it.2 == MASS_SHIFT || it.2 == RT_SHIFT)
        .count();
    if targets == 0 || decoys == 0 {
        return;
    }
    let k = decoys as f64 / targets as f64;
    let mut order: Vec<usize> = (0..items.len()).collect();
    order.sort_by(|&a, &b| items[b].1.total_cmp(&items[a].1));
    let mut fdr = vec![0.0f64; order.len()];
    let (mut t, mut d, mut d_id) = (0usize, 0usize, 0usize);
    let mut start = 0;
    while start < order.len() {
        let score = items[order[start]].1;
        let mut end = start;
        while end < order.len() && (end == start || items[order[end]].1 == score) {
            match items[order[end]].2 {
                TARGET => t += 1,
                IDENTITY => d_id += 1,
                _ => d += 1,
            }
            end += 1;
        }
        let tf = t.max(1) as f64;
        let v = ((d + 1) as f64 / (k * tf) + d_id as f64 / tf).min(1.0);
        fdr[start..end].fill(v);
        start = end;
    }
    let mut running = 1.0f64;
    for (pos, &j) in order.iter().enumerate().rev() {
        running = running.min(fdr[pos]);
        q[items[j].0] = Some(running);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Deterministic standard normal draws (splitmix64 + Box-Muller).
    struct Normal(u64);

    impl Normal {
        fn unit(&mut self) -> f64 {
            self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
            let mut z = self.0;
            z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
            z ^= z >> 31;
            ((z >> 11) as f64 + 0.5) / (1u64 << 53) as f64
        }
        fn draw(&mut self) -> f32 {
            let (u1, u2) = (self.unit(), self.unit());
            ((-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()) as f32
        }
        fn vector(&mut self, mean: f32) -> Vec<f32> {
            (0..TRANSFER_FEATURES.len())
                .map(|_| mean + self.draw())
                .collect()
        }
    }

    /// `n` precursors, each with one target (mean `target_mean`) and one
    /// decoy per listed family (mean 0).
    fn rows(n: u32, target_mean: f32, families: &[u8], seed: u64) -> Vec<FdrRow> {
        let mut g = Normal(seed);
        let mut out = Vec::new();
        for p in 0..n {
            out.push(FdrRow {
                precursor: p,
                is_decoy: false,
                family: 0,
                x: g.vector(target_mean),
            });
            for &f in families {
                out.push(FdrRow {
                    precursor: p,
                    is_decoy: true,
                    family: f,
                    x: g.vector(0.0),
                });
            }
        }
        out
    }

    fn accepted_targets(rows: &[FdrRow], q: &[Option<f64>]) -> usize {
        rows.iter()
            .zip(q)
            .filter(|(r, q)| !r.is_decoy && q.is_some_and(|q| q <= 0.01))
            .count()
    }

    #[test]
    fn separable_targets_get_low_q() {
        let rows = rows(500, 3.0, &[1, 2], 11);
        let r = transfer_qvalues(&rows, 42);
        let n = accepted_targets(&rows, &r.q);
        eprintln!(
            "separable: {n} of 500 targets at q <= 0.01 (mass-only {}, rt-only {})",
            accepted_targets(&rows, &r.q_mass_only),
            accepted_targets(&rows, &r.q_rt_only)
        );
        assert!(n >= 400, "{n} targets at q <= 0.01");
        assert!(r.score.iter().all(|s| s.is_finite()));
    }

    #[test]
    fn null_targets_are_not_accepted() {
        for seed in [3u64, 5, 9] {
            let rows = rows(500, 0.0, &[1, 2], seed);
            let r = transfer_qvalues(&rows, seed);
            let n = accepted_targets(&rows, &r.q);
            eprintln!("null seed {seed}: {n} targets at q <= 0.01");
            assert!(n <= 5, "seed {seed}: {n} null targets at q <= 0.01");
        }
    }

    #[test]
    fn null_experiments_rarely_accept_a_target() {
        // 10,000 experiments of 20 targets and 40 decoys, every score drawn
        // from one distribution: any accepted target is a false discovery,
        // so the share of experiments that accept one is their mean FDP.
        let mut g = Normal(77);
        let n_exp = 10_000;
        let mut with_false = 0usize;
        let mut q = vec![None; 60];
        for _ in 0..n_exp {
            let items: Vec<(usize, f32, u8)> = (0..60)
                .map(|i| (i, g.draw(), if i >= 20 { MASS_SHIFT } else { TARGET }))
                .collect();
            q.fill(None);
            kdecoy_q(&items, &mut q);
            if q[..20].iter().any(|q| q.is_some_and(|q| q <= 0.01)) {
                with_false += 1;
            }
        }
        let frac = with_false as f64 / n_exp as f64;
        eprintln!(
            "null experiments: {with_false} of {n_exp} ({frac:.4}) accept a target at q <= 0.01"
        );
        // 0.01 plus about five binomial standard errors at n = 10,000.
        assert!(
            frac <= 0.015,
            "{frac} of the null experiments accept a target"
        );
    }

    #[test]
    fn mixture_fdp_stays_near_the_estimate() {
        // 30% of the targets are present (shifted on four features), the rest
        // are drawn like the decoys.
        let mut g = Normal(21);
        let mut rows = Vec::new();
        let mut present = Vec::new();
        for p in 0..4000u32 {
            let real = p % 10 < 3;
            let mut x = g.vector(0.0);
            if real {
                for v in &mut x[..4] {
                    *v += 1.5;
                }
            }
            rows.push(FdrRow {
                precursor: p,
                is_decoy: false,
                family: 0,
                x,
            });
            present.push(real);
            for f in [1, 2] {
                rows.push(FdrRow {
                    precursor: p,
                    is_decoy: true,
                    family: f,
                    x: g.vector(0.0),
                });
                present.push(false);
            }
        }
        let r = transfer_qvalues(&rows, 17);
        for thr in [0.01, 0.05, 0.1] {
            let (mut acc, mut false_acc) = (0usize, 0usize);
            for (i, row) in rows.iter().enumerate() {
                if !row.is_decoy && r.q[i].is_some_and(|q| q <= thr) {
                    acc += 1;
                    false_acc += usize::from(!present[i]);
                }
            }
            let fdp = false_acc as f64 / acc.max(1) as f64;
            eprintln!("mixture: q <= {thr}: {acc} accepted, {false_acc} absent, FDP {fdp:.4}");
            assert!(acc > 100, "q <= {thr}: {acc} accepted");
            assert!(fdp <= 2.0 * thr + 0.01, "q <= {thr}: FDP {fdp}");
        }
    }

    #[test]
    fn folds_group_by_precursor() {
        let rows = rows(600, 3.0, &[1, 2], 1);
        let mut fold_of_precursor = std::collections::HashMap::new();
        let mut sizes = [0usize; N_FOLDS];
        for r in &rows {
            let f = fold_of(r.precursor, 7);
            assert!(f < N_FOLDS);
            let first = *fold_of_precursor.entry(r.precursor).or_insert(f);
            assert_eq!(first, f, "precursor {} spans folds", r.precursor);
        }
        for f in fold_of_precursor.values() {
            sizes[*f] += 1;
        }
        assert!(sizes.iter().all(|&s| s > 150), "fold sizes {sizes:?}");
        let moved = (0..600u32)
            .filter(|&p| fold_of(p, 7) != fold_of(p, 8))
            .count();
        assert!(moved > 200, "the seed changes {moved} of 600 folds");
    }

    #[test]
    fn fdr_degenerate_inputs_yield_none() {
        let empty = transfer_qvalues(&[], 1);
        assert!(empty.q.is_empty() && empty.score.is_empty());

        let no_decoys = rows(300, 3.0, &[], 2);
        let r = transfer_fdr(&no_decoys, 1);
        assert!(r.result.q.iter().all(Option::is_none));
        assert_eq!(r.counts.failed_folds.len(), N_FOLDS);
        assert_eq!(r.counts.unscored, no_decoys.len());

        let only_decoys: Vec<FdrRow> = rows(300, 3.0, &[1, 2], 3)
            .into_iter()
            .filter(|r| r.is_decoy)
            .collect();
        let r = transfer_qvalues(&only_decoys, 1);
        assert!(r.q.iter().all(Option::is_none));

        // Fold 0 holds targets only: its rows have no decoys to estimate from.
        let seed = 5;
        let mixed: Vec<FdrRow> = rows(600, 3.0, &[1, 2], 4)
            .into_iter()
            .filter(|r| !(r.is_decoy && fold_of(r.precursor, seed) == 0))
            .collect();
        let r = transfer_fdr(&mixed, seed);
        for (row, q) in mixed.iter().zip(&r.result.q) {
            if fold_of(row.precursor, seed) == 0 {
                assert!(q.is_none(), "one-class fold row got q {q:?}");
            }
        }
        assert!(accepted_targets(&mixed, &r.result.q) > 200);
        assert!(r.counts.no_q > 0);
    }

    #[test]
    fn per_family_estimates_reported() {
        let rows = rows(500, 3.0, &[1], 6);
        let r = transfer_qvalues(&rows, 42);
        assert!(r.q_rt_only.iter().all(Option::is_none));
        assert!(accepted_targets(&rows, &r.q_mass_only) >= 400);
        assert_eq!(r.q, r.q_mass_only);
    }

    #[test]
    fn q_values_are_monotone_in_score_within_a_fold() {
        let rows = rows(400, 1.0, &[1, 2], 8);
        let r = transfer_qvalues(&rows, 3);
        let mut by_fold: Vec<Vec<(f32, f64)>> = vec![Vec::new(); N_FOLDS];
        for (i, row) in rows.iter().enumerate() {
            if let Some(q) = r.q[i] {
                by_fold[fold_of(row.precursor, 3)].push((r.score[i], q));
            }
        }
        for f in by_fold.iter_mut() {
            assert!(!f.is_empty());
            f.sort_by(|a, b| a.0.total_cmp(&b.0));
            for w in f.windows(2) {
                assert!(w[1].1 <= w[0].1, "q rises with score: {w:?}");
            }
            assert!(f.iter().all(|&(_, q)| (0.0..=1.0).contains(&q)));
        }
    }

    #[test]
    fn k_follows_the_decoy_to_target_ratio() {
        // Four decoys per target at the same scores as the targets: FDR = 1
        // only when K is the observed ratio 4, not an assumed 2.
        let mut rows = Vec::new();
        for p in 0..400u32 {
            let v = vec![p as f32; TRANSFER_FEATURES.len()];
            rows.push(FdrRow {
                precursor: p,
                is_decoy: false,
                family: 0,
                x: v.clone(),
            });
            for _ in 0..4 {
                rows.push(FdrRow {
                    precursor: p,
                    is_decoy: true,
                    family: 1,
                    x: v.clone(),
                });
            }
        }
        let mut q = vec![None; rows.len()];
        let scored: Vec<(usize, f32, u8)> = rows
            .iter()
            .enumerate()
            .map(|(i, r)| (i, r.x[0], r.family))
            .collect();
        kdecoy_q(&scored, &mut q);
        assert!(
            q.iter().all(|q| (q.unwrap() - 1.0).abs() < 1e-12),
            "{:?}",
            &q[..5]
        );
    }

    #[test]
    fn pseudo_rows_take_the_q_of_the_nearest_lower_score() {
        let rows = rows(500, 3.0, &[1, 2], 12);
        let fdr = transfer_fdr(&rows, 4);
        // A pseudo row identical to a scored target gets that target's q.
        for (i, row) in rows
            .iter()
            .enumerate()
            .filter(|(_, r)| !r.is_decoy)
            .take(50)
        {
            assert_eq!(fdr.pseudo_q(row.precursor, &row.x), fdr.result.q[i]);
        }
        // Far below every row: the fold's largest q.
        let low = vec![-50.0f32; TRANSFER_FEATURES.len()];
        assert!(fdr.pseudo_q(0, &low).is_some_and(|q| q > 0.1));
    }

    /// `rows(n, 3.0, [1, 2])` plus `n_id` identity decoys drawn like the
    /// targets: false donors whose transfers look like real ones.
    fn with_identity(n: u32, n_id: u32, seed: u64) -> Vec<FdrRow> {
        let mut out = rows(n, 3.0, &[MASS_SHIFT, RT_SHIFT], seed);
        let mut g = Normal(seed ^ 0xabc);
        out.extend((0..n_id).map(|i| FdrRow {
            precursor: n + i,
            is_decoy: true,
            family: IDENTITY,
            x: g.vector(3.0),
        }));
        out
    }

    #[test]
    fn identity_decoys_add_their_share_to_the_estimate() {
        // 1,000 well separated targets and 30 identity decoys that look like
        // targets: the wrong-peak estimate accepts nearly every target at
        // 0.01, the combined one sees about 3 % false donors and keeps few.
        let rows = with_identity(1000, 30, 31);
        {
            let r = transfer_fdr(&rows, 5).result;
            let wp = accepted_targets(&rows, &r.q_wrong_peak);
            let all = accepted_targets(&rows, &r.q);
            eprintln!("identity: wrong-peak {wp}, combined {all} at q <= 0.01");
            assert!(wp >= 800, "{wp}");
            // Only targets above every identity decoy of their fold remain.
            assert!(5 * all < wp, "{all} of {wp}");
            let at10 = rows
                .iter()
                .zip(&r.q)
                .filter(|(row, q)| row.family == TARGET && q.is_some_and(|q| q <= 0.1))
                .count();
            assert!(at10 >= 800, "{at10}");
            // Identity decoys get q-values like the targets they mimic.
            assert!(rows
                .iter()
                .zip(&r.q)
                .filter(|(row, _)| row.family == IDENTITY)
                .all(|(_, q)| q.is_some()));
        }
    }

    #[test]
    fn kdecoy_counts_identity_decoys_once_each() {
        // 100 targets, then 50 mass-shift decoys and 2 identity decoys below
        // them: K = 0.5. Over all rows (D + 1) / (K T) + D_id / T =
        // 51 / 50 + 2 / 100, capped at 1; over the targets alone 1 / 50.
        let mut items: Vec<(usize, f32, u8)> =
            (0..100).map(|i| (i, 1000.0 - i as f32, TARGET)).collect();
        items.extend((100..150).map(|i| (i, 10.0 - i as f32 * 0.01, MASS_SHIFT)));
        items.extend([(150, 0.5, IDENTITY), (151, 0.4, IDENTITY)]);
        let mut q = vec![None; items.len()];
        kdecoy_q(&items, &mut q);
        assert!((q[0].unwrap() - 1.0 / 50.0).abs() < 1e-12, "{:?}", q[0]);
        assert_eq!(q[151], Some(1.0));
        // An identity decoy above the targets counts against every target.
        items[150].1 = 2000.0;
        kdecoy_q(&items, &mut q);
        assert!((q[0].unwrap() - (1.0 / 50.0 + 1.0 / 100.0)).abs() < 1e-12);
    }

    #[test]
    fn results_are_deterministic_for_a_seed() {
        let rows = rows(300, 1.5, &[1, 2], 13);
        assert_eq!(transfer_qvalues(&rows, 9), transfer_qvalues(&rows, 9));
    }

    #[test]
    fn step_functions_match_the_per_row_estimate() {
        // Coarse scores (many ties, both zeros) over every family.
        let mut g = Normal(99);
        for trial in 0..20u64 {
            let n = 50 + 37 * trial as usize;
            let items: Vec<(usize, f32, u8)> = (0..n)
                .map(|i| {
                    let fam = match i % 7 {
                        0..=2 => TARGET,
                        3 | 4 => MASS_SHIFT,
                        5 => RT_SHIFT,
                        _ => IDENTITY,
                    };
                    let shift = if fam == TARGET { 1.0 } else { 0.0 };
                    let mut s = ((g.draw() + shift) * 4.0).round() / 4.0;
                    if s == 0.0 && i % 2 == 0 {
                        s = -0.0;
                    }
                    (i, s, fam)
                })
                .collect();
            let mut by_family: [Vec<f32>; 4] = Default::default();
            for &(_, s, f) in &items {
                by_family[f as usize].push(canonical(s));
            }
            for v in &mut by_family {
                v.sort_by(f32::total_cmp);
            }
            for e in ESTIMATES {
                let kept: Vec<(usize, f32, u8)> =
                    items.iter().copied().filter(|it| e.keeps(it.2)).collect();
                let mut want = vec![None; n];
                kdecoy_q(&kept, &mut want);
                let steps = kdecoy_steps(&by_family, e);
                for &(i, s, f) in &items {
                    let got = steps
                        .as_ref()
                        .filter(|_| e.keeps(f))
                        .map(|st| st.of_member(s));
                    assert_eq!(got, want[i], "trial {trial} {e:?} row {i} score {s}");
                }
                // A row outside the fit takes the q of the highest score at
                // or below it, or the lowest row's.
                if let (Some(st), Estimate::Pooled) = (&steps, e) {
                    let mut table: Vec<(f32, f64)> = kept
                        .iter()
                        .filter_map(|it| want[it.0].map(|q| (it.1, q)))
                        .collect();
                    table.sort_by(|a, b| a.0.total_cmp(&b.0).then(b.1.total_cmp(&a.1)));
                    for probe in [-9.0f32, -0.3, -0.0, 0.0, 0.1, 0.25, 1.3, 9.0, f32::NAN] {
                        let at = table.partition_point(|&(s, _)| s <= probe);
                        let old = if at == 0 { table[0].1 } else { table[at - 1].1 };
                        assert_eq!(st.at_or_below(probe), old, "probe {probe}");
                    }
                }
            }
        }
    }

    #[test]
    fn threshold_steps_accept_as_the_whole_ones() {
        let mut g = Normal(7);
        let (mut whole, mut cut) = (FoldScores::new(3), FoldScores::new(3));
        let mut items = Vec::new();
        for i in 0..6000u32 {
            let fam = match i % 7 {
                0..=2 => TARGET,
                3 | 4 => MASS_SHIFT,
                5 => RT_SHIFT,
                _ => IDENTITY,
            };
            let shift = if fam == TARGET { 1.5 } else { 0.0 };
            let s = Some(((g.draw() + shift) * 8.0).round() / 8.0);
            items.push((i % 800, fam, s));
            whole.add(i % 800, fam, s);
            cut.add(i % 800, fam, s);
        }
        let ts = [0.01, 0.05, 0.3];
        let whole = whole.finish(&[]);
        let cut = cut.finish_with(&[], Some(&ts), |_, _| Vec::new());
        assert!(cut.steps() < whole.steps());
        assert_eq!(cut.pooled_steps(), whole.pooled_steps());
        for &(p, f, s) in &items {
            for e in ESTIMATES {
                for t in ts {
                    assert_eq!(whole.accepts(e, p, f, s, t), cut.accepts(e, p, f, s, t));
                }
            }
            assert_eq!(
                whole.q(Estimate::Pooled, p, f, s),
                cut.q(Estimate::Pooled, p, f, s)
            );
            assert_eq!(whole.pseudo_q(p, s), cut.pseudo_q(p, s));
        }
    }

    /// Pooled q-values of `rows` under a classifier of `members` members,
    /// each trained on all the rows.
    fn bagged_q(rows: &[FdrRow], members: usize, seed: u64) -> Vec<Option<f64>> {
        let sets: Vec<&[FdrRow]> = (0..members).map(|_| rows).collect();
        let models = train_bagged_folds(&sets, seed);
        assert_eq!(models.members(), members);
        let precursors: Vec<u32> = rows.iter().map(|r| r.precursor).collect();
        let xs: Vec<&[f32]> = rows.iter().map(|r| r.x.as_slice()).collect();
        let mut scores = Vec::new();
        models.score(&precursors, &xs, &mut scores);
        let mut fs = FoldScores::new(seed);
        for (r, s) in rows.iter().zip(&scores) {
            fs.add(r.precursor, r.family, *s);
        }
        let q = fs.finish(&models.failed_folds);
        rows.iter()
            .zip(&scores)
            .map(|(r, s)| q.q(Estimate::Pooled, r.precursor, r.family, *s))
            .collect()
    }

    #[test]
    fn one_member_is_the_single_model() {
        let rows = rows(400, 1.0, &[1, 2], 23);
        assert_eq!(bagged_q(&rows, 1, 6), transfer_qvalues(&rows, 6).q);
    }

    #[test]
    fn bagged_scores_are_the_mean_of_the_members() {
        let rows = rows(300, 1.0, &[1, 2], 29);
        let sets: Vec<&[FdrRow]> = vec![&rows; 3];
        let bagged = train_bagged_folds(&sets, 4);
        let xs: Vec<&[f32]> = rows.iter().map(|r| r.x.as_slice()).collect();
        let precursors: Vec<u32> = rows.iter().map(|r| r.precursor).collect();
        let mut got = Vec::new();
        bagged.score(&precursors, &xs, &mut got);
        for (fold, members) in bagged.models.iter().enumerate() {
            assert_eq!(members.len(), 3);
            for (i, r) in rows.iter().enumerate() {
                if fold_of(r.precursor, 4) != fold {
                    continue;
                }
                let sum: f32 = members.iter().map(|m| m.predict_value(&r.x)).sum();
                let want = sum / 3.0;
                let g = got[i].unwrap();
                assert!(
                    (g - want).abs() <= 1e-5 * want.abs().max(1.0),
                    "{g} vs {want}"
                );
            }
        }
        // Members differ (their trainer seeds), so the mean is not one of them.
        let m = &bagged.models[0];
        assert_ne!(
            m[0].predict_value(&rows[0].x),
            m[1].predict_value(&rows[0].x)
        );
    }

    #[test]
    fn bagging_keeps_the_mixture_calibrated_and_steadier_across_seeds() {
        // As mixture_fdp_stays_near_the_estimate, with weaker targets: the
        // accepted sets of two seeds overlap more with 8 members than with 1.
        let mut g = Normal(41);
        let mut rows = Vec::new();
        let mut present = Vec::new();
        for p in 0..3000u32 {
            let real = p % 10 < 3;
            let mut x = g.vector(0.0);
            if real {
                for v in &mut x[..4] {
                    *v += 1.2;
                }
            }
            rows.push(FdrRow {
                precursor: p,
                is_decoy: false,
                family: 0,
                x,
            });
            present.push(real);
            for f in [1, 2] {
                rows.push(FdrRow {
                    precursor: p,
                    is_decoy: true,
                    family: f,
                    x: g.vector(0.0),
                });
                present.push(false);
            }
        }
        let accepted = |q: &[Option<f64>]| -> std::collections::BTreeSet<usize> {
            (0..rows.len())
                .filter(|&i| !rows[i].is_decoy && q[i].is_some_and(|q| q <= 0.05))
                .collect()
        };
        let jaccard = |a: &std::collections::BTreeSet<usize>,
                       b: &std::collections::BTreeSet<usize>| {
            a.intersection(b).count() as f64 / a.union(b).count().max(1) as f64
        };
        let mut overlap = Vec::new();
        for members in [1, 8] {
            let (a, b) = (
                accepted(&bagged_q(&rows, members, 1)),
                accepted(&bagged_q(&rows, members, 2)),
            );
            for set in [&a, &b] {
                let false_acc = set.iter().filter(|&&i| !present[i]).count();
                let fdp = false_acc as f64 / set.len().max(1) as f64;
                assert!(set.len() > 100, "{members} members: {} accepted", set.len());
                assert!(fdp <= 0.11, "{members} members: FDP {fdp}");
            }
            let j = jaccard(&a, &b);
            eprintln!(
                "{members} members: accepted {} and {}, Jaccard {j:.4}",
                a.len(),
                b.len()
            );
            overlap.push(j);
        }
        assert!(overlap[1] > overlap[0], "{overlap:?}");
    }

    #[test]
    fn training_sample_wants_what_push_keeps() {
        let mut a: TrainingSample<usize, u32> = TrainingSample::new(500, 9);
        let mut b: TrainingSample<usize, u32> = TrainingSample::new(500, 9);
        for i in 0..5000usize {
            let unit = (i / 3) as u64;
            a.push(i, unit, i as u32);
            if b.wants(unit) {
                b.push(i, unit, i as u32);
            } else {
                b.skip();
            }
        }
        assert_eq!(a.stats(), b.stats());
        assert_eq!(a.into_items(), b.into_items());
    }

    #[test]
    fn training_sample_keeps_whole_precursors_whatever_the_order() {
        let rows = rows(3000, 1.0, &[1, 2], 5);
        let mut a: TrainingSample<usize> = TrainingSample::new(2500, 42);
        for (i, r) in rows.iter().enumerate() {
            a.push(i, r.precursor as u64, r.clone());
        }
        let mut b: TrainingSample<usize> = TrainingSample::new(2500, 42);
        for (i, r) in rows.iter().enumerate().rev() {
            b.push(i, r.precursor as u64, r.clone());
        }
        let (seen, kept, share) = a.stats();
        assert_eq!(seen, 9000);
        assert!(kept <= 2500 && kept > 1000, "{kept}");
        assert!(share < 1.0);
        let (ra, rb) = (a.into_items(), b.into_items());
        assert_eq!(
            ra, rb,
            "the kept set and its order do not depend on arrival"
        );
        // Whole precursors: three rows each.
        let mut per: std::collections::HashMap<u32, usize> = std::collections::HashMap::new();
        for r in &ra {
            *per.entry(r.precursor).or_default() += 1;
        }
        assert!(per.values().all(|&c| c == 3));
        // Under the cap every row is kept, in key order.
        let mut all: TrainingSample<usize> = TrainingSample::new(10_000, 42);
        for (i, r) in rows.iter().enumerate() {
            all.push(i, r.precursor as u64, r.clone());
        }
        assert_eq!(all.into_items(), rows);
    }
}
