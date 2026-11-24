use model::hash::{splitmix64, SPLITMIX_GOLDEN};
use model_train::gbdt::train::{train_gbdt, Dataset, TrainParams};
use rayon::prelude::*;
use scoring_crate::gbdt_eval::GbdtPeakModel;

pub const N_FOLDS: usize = 3;

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
    "log10_donor_pep",
];

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

pub const MAX_STRATA: usize = 16;
pub const EB_LEVEL: f64 = 0.01;
pub fn stratum_of(prior: f64, edges: &[f64]) -> u8 {
    if !prior.is_finite() {
        return 0;
    }
    edges.iter().filter(|&&e| prior > e).count() as u8
}

pub fn stratum_kappa(strata: &[(f64, f64)], k: f64) -> f64 {
    let live: Vec<&(f64, f64)> = strata.iter().filter(|(n, _)| *n > 0.0).collect();
    let s = live.len() as f64;
    let n_tot: f64 = live.iter().map(|(n, _)| n).sum();
    let x_tot: f64 = live.iter().map(|(_, x)| x).sum();
    if live.len() < 2 || x_tot <= 0.0 || k <= 0.0 {
        return f64::INFINITY;
    }
    let p = x_tot / n_tot;
    let v0 = p / k;
    let z: f64 = live
        .iter()
        .map(|(n, x)| n * (x / n - p).powi(2))
        .sum::<f64>()
        / v0;
    let sq: f64 = live.iter().map(|(n, _)| n * n).sum::<f64>() / n_tot;
    let denom = n_tot - sq - (s - 1.0);
    if denom <= 0.0 {
        return f64::INFINITY;
    }
    let rho = (z - (s - 1.0)) / denom;
    if rho <= 0.0 {
        f64::INFINITY
    } else {
        ((1.0 - rho) / rho).max(0.0)
    }
}

pub fn pack_family(family: u8, stratum: impl Into<u32>) -> u32 {
    u32::from(family & 0x0F) | (stratum.into() << 4)
}

pub const fn unpack_family(packed: u32) -> (u8, usize) {
    ((packed & 0x0F) as u8, (packed >> 4) as usize)
}

pub const TARGET: u8 = 0;
pub const MASS_SHIFT: u8 = 1;
pub const RT_SHIFT: u8 = 2;
pub const IDENTITY: u8 = 3;

#[derive(Debug, Clone, PartialEq)]
pub struct FdrRow {
    pub precursor: u32,
    pub is_decoy: bool,
    pub family: u8,
    pub x: Vec<f32>,
}

#[cfg(test)]
#[derive(Debug, Clone, PartialEq)]
pub struct FdrResult {
    pub q: Vec<Option<f64>>,
    pub q_wrong_peak: Vec<Option<f64>>,
    pub q_mass_only: Vec<Option<f64>>,
    pub q_rt_only: Vec<Option<f64>>,
    pub score: Vec<f32>,
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct FdrCounts {
    pub failed_folds: Vec<String>,
    pub unscored: usize,
    pub no_q: usize,
}

pub fn fold_of(precursor: u32, seed: u64) -> usize {
    let z = splitmix64(seed ^ (precursor as u64).wrapping_mul(SPLITMIX_GOLDEN));
    (z % N_FOLDS as u64) as usize
}

#[cfg(test)]
pub fn transfer_qvalues(rows: &[FdrRow], seed: u64) -> FdrResult {
    transfer_fdr(rows, seed).result
}

fn train_fold(
    rows: &[FdrRow],
    folds: &[usize],
    nf: usize,
    fold: usize,
    seed: u64,
    params: &TrainParams,
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
    train_gbdt(&ds, params, seed.wrapping_add(fold as u64 + 1))
        .map_err(|e| format!("fold {fold}: {e:?}"))
}

pub struct FoldModels {
    seed: u64,
    models: Vec<Vec<GbdtPeakModel>>,
    pub failed_folds: Vec<String>,
}

pub fn member_seed(seed: u64, m: usize) -> u64 {
    if m == 0 {
        seed
    } else {
        sample_hash(m as u64, seed ^ 0x6d65_6d62_6572_7321)
    }
}

pub fn train_folds(rows: &[FdrRow], seed: u64) -> FoldModels {
    train_bagged_folds(&[rows], seed, &TrainParams::default())
}

pub fn transfer_train_params() -> TrainParams {
    TrainParams {
        learning_rate: TRANSFER_LEARNING_RATE,
        ..TrainParams::default()
    }
}

pub const TRANSFER_LEARNING_RATE: f32 = 0.3;

pub fn train_bagged_folds(sets: &[&[FdrRow]], seed: u64, params: &TrainParams) -> FoldModels {
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
        .map(|&(f, m)| train_fold(sets[m], &folds[m], nf, f, member_seed(seed, m), params))
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
            log::warn!(
                target: "mbr",
                "fold {f}: {} of {} classifier members failed, the fold averages the others ({})",
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

const PREDICT_BATCH: usize = 4096;

impl FoldModels {
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

    pub fn members(&self) -> usize {
        self.models.iter().map(Vec::len).max().unwrap_or(0)
    }
}

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

    pub fn wants(&self, unit: u64) -> bool {
        Self::kept(self.seed, unit, self.level)
    }

    pub fn skip(&mut self) {
        self.seen += 1;
    }

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

    pub fn stats(&self) -> (usize, usize, f64) {
        (self.seen, self.items.len(), 0.5f64.powi(self.level as i32))
    }

    pub fn into_items(self) -> Vec<T> {
        let mut items = self.items;
        items.sort_by(|a, b| a.0.cmp(&b.0));
        items.into_iter().map(|(_, _, t)| t).collect()
    }
}

pub fn row_unit(precursor: u32, run: u16) -> u64 {
    ((precursor as u64) << 16) | run as u64
}

pub fn sample_hash(unit: u64, seed: u64) -> u64 {
    splitmix64(seed ^ 0x5bd1_e995_7a4f_3c21 ^ unit.wrapping_mul(SPLITMIX_GOLDEN))
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Estimate {
    Pooled = 0,
    WrongPeak = 1,
    MassOnly = 2,
    RtOnly = 3,
}

const ESTIMATES: [Estimate; 4] = [
    Estimate::Pooled,
    Estimate::WrongPeak,
    Estimate::MassOnly,
    Estimate::RtOnly,
];

impl Estimate {
    fn keeps(self, family: u8) -> bool {
        match self {
            Estimate::Pooled => true,
            Estimate::WrongPeak => family != IDENTITY,
            Estimate::MassOnly => family == TARGET || family == MASS_SHIFT,
            Estimate::RtOnly => family == TARGET || family == RT_SHIFT,
        }
    }
}

fn canonical(s: f32) -> f32 {
    if s == 0.0 {
        0.0
    } else {
        s
    }
}

#[derive(Debug, Clone, Default)]
pub(crate) struct Steps {
    at: Vec<f32>,
    q: Vec<f64>,
}

impl Steps {
    fn of_member(&self, s: f32) -> f64 {
        let s = canonical(s);
        let k = self.at.partition_point(|&a| a.total_cmp(&s).is_le());
        self.q[k.max(1) - 1]
    }

    fn at_or_below(&self, s: f32) -> f64 {
        let k = self.at.partition_point(|&a| a <= s);
        self.q[k.max(1) - 1]
    }
}

pub struct FoldScores {
    seed: u64,
    strata: usize,
    scores: Vec<[Vec<f32>; 4]>,
    counts: Vec<[usize; 4]>,
    members: Vec<usize>,
    unscored: usize,
}

impl FoldScores {
    pub fn new(seed: u64) -> Self {
        Self::with_strata(seed, 1)
    }

    pub fn with_strata(seed: u64, strata: usize) -> Self {
        let groups = N_FOLDS * strata.max(1);
        Self {
            seed,
            strata: strata.max(1),
            scores: (0..groups).map(|_| Default::default()).collect(),
            counts: vec![[0; 4]; groups],
            members: vec![0; groups],
            unscored: 0,
        }
    }

    pub fn add(&mut self, precursor: u32, family: u32, score: Option<f32>) {
        let (family, stratum) = unpack_family(family);
        let group = fold_of(precursor, self.seed) * self.strata + stratum.min(self.strata - 1);
        match score {
            None => self.unscored += 1,
            Some(s) => {
                self.members[group] += 1;
                if s.is_finite() {
                    let f = family.min(3) as usize;
                    self.scores[group][f].push(canonical(s));
                    self.counts[group][f] += 1;
                }
            }
        }
    }

    pub fn merge(&mut self, other: FoldScores) {
        assert_eq!(self.seed, other.seed, "fold scores of different seeds");
        assert_eq!(self.strata, other.strata, "fold scores of different strata");
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

    pub fn held(&self) -> usize {
        self.scores.iter().flatten().map(Vec::len).sum()
    }

    pub fn take_held(&mut self) -> Vec<[Vec<f32>; 4]> {
        self.scores
            .iter_mut()
            .map(|by_family| by_family.each_mut().map(std::mem::take))
            .collect()
    }

    #[cfg(test)]
    pub fn finish(self, failed_folds: &[String]) -> QValues {
        self.finish_with(failed_folds, None, |_, _| Vec::new())
    }

    pub fn finish_with<'a>(
        mut self,
        failed_folds: &[String],
        thresholds: Option<&[f64]>,
        mut taken: impl FnMut(usize, usize) -> Vec<Box<dyn Iterator<Item = f32> + 'a>>,
    ) -> QValues {
        let strata = self.strata;
        let mut steps: Vec<[Option<Steps>; 4]> = Vec::with_capacity(N_FOLDS);
        let mut stratum_steps: Vec<Vec<Option<Steps>>> = Vec::new();
        let mut no_q = 0usize;
        let mut peps: Vec<Option<PepSteps>> = Vec::with_capacity(N_FOLDS);
        let mut stratum_beta: Vec<Vec<f64>> = Vec::new();
        for v in self.scores.iter_mut().flatten() {
            v.sort_unstable_by(f32::total_cmp);
        }
        for fold in 0..N_FOLDS {
            let groups = fold * strata..(fold + 1) * strata;
            let mut totals = [0usize; 4];
            for g in groups.clone() {
                for (t, c) in totals.iter_mut().zip(self.counts[g]) {
                    *t += c;
                }
            }
            let streams: [SortedMerge<'_>; 4] = std::array::from_fn(|f| {
                let mut runs: Vec<Box<dyn Iterator<Item = f32> + '_>> = Vec::new();
                for g in groups.clone() {
                    runs.extend(taken(g, f));
                    runs.push(Box::new(self.scores[g][f].iter().copied()));
                }
                SortedMerge::new(runs)
            });
            let mut builder = StepsBuilder::new(totals, thresholds);
            let mut pep = PepBuilder::new(totals);
            walk_scores(streams, |c, at| {
                builder.visit(c, at);
                pep.visit(c, at);
                true
            });
            let fold_steps = builder.finish();
            peps.push(pep.finish());
            let with_q = if fold_steps[Estimate::Pooled as usize].is_some() {
                totals.iter().sum()
            } else {
                0
            };
            no_q += groups.clone().map(|g| self.members[g]).sum::<usize>() - with_q;
            if strata > 1 {
                let mut per_stratum: Vec<Option<Steps>> = Vec::with_capacity(strata);
                let k = (totals[MASS_SHIFT as usize] + totals[RT_SHIFT as usize]).max(1) as f64
                    / totals[TARGET as usize].max(1) as f64;
                let global = fold_steps[Estimate::Pooled as usize].as_ref();
                let mut kappa = f64::INFINITY;
                if let Some(c_op) = global.and_then(|g| {
                    g.at.iter()
                        .zip(&g.q)
                        .find(|(_, &q)| q <= EB_LEVEL)
                        .or_else(|| g.at.iter().zip(&g.q).next_back())
                        .map(|(a, _)| *a)
                }) {
                    let mut counts = Vec::with_capacity(strata);
                    for g in groups.clone() {
                        let local = self.counts[g];
                        let streams: [SortedMerge<'_>; 4] = std::array::from_fn(|f| {
                            let mut runs: Vec<Box<dyn Iterator<Item = f32> + '_>> = taken(g, f);
                            runs.push(Box::new(self.scores[g][f].iter().copied()));
                            SortedMerge::new(runs)
                        });
                        let mut below = [0usize; 4];
                        walk_scores(streams, |c, at| {
                            if c >= c_op {
                                return false;
                            }
                            for (b, a) in below.iter_mut().zip(at) {
                                *b += a;
                            }
                            true
                        });
                        let n = (local[TARGET as usize] - below[TARGET as usize]) as f64;
                        let d = (local[MASS_SHIFT as usize] + local[RT_SHIFT as usize]
                            - below[MASS_SHIFT as usize]
                            - below[RT_SHIFT as usize]) as f64;
                        let d_id = (local[IDENTITY as usize] - below[IDENTITY as usize]) as f64;
                        counts.push((n, d / k + d_id));
                    }
                    kappa = stratum_kappa(&counts, k);
                }
                log::debug!(target: "mbr", "transfer q-values, fold {fold}: {strata} strata, shrinkage {kappa:.1} pseudo-targets");
                let mut fold_betas = Vec::with_capacity(strata);
                if let Some(pool) = peps.last().and_then(Option::as_ref) {
                    for g in groups.clone() {
                        let mut z_of = |f: u8| -> Vec<f64> {
                            let mut runs = taken(g, f as usize);
                            runs.push(Box::new(self.scores[g][f as usize].iter().copied()));
                            runs.into_iter()
                                .flatten()
                                .filter(|s| s.is_finite())
                                .map(|s| logit(pool.of(s).0))
                                .collect()
                        };
                        let targets = z_of(TARGET);
                        let mut decoys = z_of(MASS_SHIFT);
                        decoys.extend(z_of(RT_SHIFT));
                        fold_betas.push(logit_shift(&targets, &decoys, k));
                    }
                    log::debug!(target: "mbr", "transfer PEP, fold {fold}: logit shifts by stratum {fold_betas:.2?}");
                }
                stratum_beta.push(fold_betas);
                for g in groups {
                    let local = self.counts[g];
                    let walk = global.filter(|_| local[TARGET as usize] > 0 && kappa.is_finite());
                    let Some(global) = walk else {
                        per_stratum.push(None);
                        continue;
                    };
                    let mut w = StratumWalk::new(k, local, global, kappa);
                    let streams: [SortedMerge<'_>; 4] = std::array::from_fn(|f| {
                        let mut runs: Vec<Box<dyn Iterator<Item = f32> + '_>> = taken(g, f);
                        runs.push(Box::new(self.scores[g][f].iter().copied()));
                        SortedMerge::new(runs)
                    });
                    walk_scores(streams, |c, at| {
                        w.visit(c, at);
                        true
                    });
                    per_stratum.push(w.finish());
                }
                stratum_steps.push(per_stratum);
            }
            steps.push(fold_steps);
        }
        QValues {
            seed: self.seed,
            strata,
            steps,
            stratum_steps,
            peps,
            stratum_beta,
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

fn logit(p: f64) -> f64 {
    let p = p.clamp(1e-9, 1.0 - 1e-9);
    (p / (1.0 - p)).ln()
}

fn sigmoid(z: f64) -> f64 {
    1.0 / (1.0 + (-z).exp())
}

fn logit_shift(targets: &[f64], decoys: &[f64], k: f64) -> f64 {
    if targets.is_empty() || decoys.is_empty() {
        return 0.0;
    }
    let score = |b: f64| {
        let term = |z: f64| {
            let pep = sigmoid(z + b);
            let p = k * pep / (k * pep + 1.0);
            (p, 1.0 - pep)
        };
        let up: f64 = decoys.iter().map(|&z| { let (p, w) = term(z); (1.0 - p) * w }).sum();
        let down: f64 = targets.iter().map(|&z| { let (p, w) = term(z); p * w }).sum();
        up - down
    };
    let (mut lo, mut hi) = (-30.0, 30.0);
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        if score(mid) > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

#[derive(Debug, Clone, Default)]
pub(crate) struct PepSteps {
    wrong_peak: Antitone,
    identity: Antitone,
}

impl PepSteps {
    pub(crate) fn of(&self, s: f32) -> (f64, f64) {
        (self.wrong_peak.of(s), self.identity.of(s))
    }
}

#[derive(Debug, Clone, Default)]
pub(crate) struct Antitone {
    lo: Vec<f32>,
    rate: Vec<f64>,
}

impl Antitone {
    fn of(&self, s: f32) -> f64 {
        if self.rate.is_empty() {
            return 0.0;
        }
        let k = self.lo.partition_point(|&a| a <= s);
        self.rate[k.saturating_sub(1)]
    }
}

#[derive(Default)]
struct AntitoneFit {
    blocks: Vec<(f32, f64, f64)>,
}

impl AntitoneFit {
    fn push(&mut self, c: f32, t: f64, f: f64) {
        if t == 0.0 && f == 0.0 {
            return;
        }
        self.blocks.push((c, t, f));
        let rate = |t: f64, f: f64| if t > 0.0 { f / t } else { f64::INFINITY };
        while self.blocks.len() >= 2 {
            let n = self.blocks.len();
            let (_, ta, fa) = self.blocks[n - 2];
            let (_, tb, fb) = self.blocks[n - 1];
            if rate(tb, fb) <= rate(ta, fa) {
                break;
            }
            self.blocks.pop();
            let a = self.blocks.last_mut().expect("two blocks");
            a.1 += tb;
            a.2 += fb;
        }
    }

    fn finish(self) -> Antitone {
        let (lo, rate) = self
            .blocks
            .iter()
            .map(|&(c, t, f)| (c, if t > 0.0 { (f / t).min(1.0) } else { 1.0 }))
            .unzip();
        Antitone { lo, rate }
    }
}

pub(crate) struct PepBuilder {
    k: f64,
    wrong_peak: AntitoneFit,
    identity: AntitoneFit,
}

impl PepBuilder {
    pub(crate) fn new(totals: [usize; 4]) -> Self {
        let decoys = totals[MASS_SHIFT as usize] + totals[RT_SHIFT as usize];
        let targets = totals[TARGET as usize];
        Self {
            k: if targets > 0 && decoys > 0 {
                decoys as f64 / targets as f64
            } else {
                0.0
            },
            wrong_peak: AntitoneFit::default(),
            identity: AntitoneFit::default(),
        }
    }

    pub(crate) fn visit(&mut self, c: f32, at: [usize; 4]) {
        if self.k <= 0.0 {
            return;
        }
        let t = at[TARGET as usize] as f64;
        let d = (at[MASS_SHIFT as usize] + at[RT_SHIFT as usize]) as f64 / self.k;
        self.wrong_peak.push(c, t, d);
        self.identity.push(c, t, at[IDENTITY as usize] as f64);
    }

    pub(crate) fn finish(self) -> Option<PepSteps> {
        if self.wrong_peak.blocks.is_empty() {
            return None;
        }
        Some(PepSteps {
            wrong_peak: self.wrong_peak.finish(),
            identity: self.identity.finish(),
        })
    }
}

struct StratumWalk<'g> {
    k: f64,
    t: usize,
    d: usize,
    d_id: usize,
    global: &'g Steps,
    kappa: f64,
    running: f64,
    steps: Steps,
}

impl<'g> StratumWalk<'g> {
    fn new(k: f64, local: [usize; 4], global: &'g Steps, kappa: f64) -> Self {
        Self {
            k,
            t: local[TARGET as usize],
            d: local[MASS_SHIFT as usize] + local[RT_SHIFT as usize],
            d_id: local[IDENTITY as usize],
            global,
            kappa,
            running: 1.0,
            steps: Steps::default(),
        }
    }

    fn visit(&mut self, c: f32, at: [usize; 4]) {
        if at.iter().all(|&n| n == 0) {
            return;
        }
        let kappa = self.kappa;
        let local =
            ((self.d + 1) as f64 / self.k + self.d_id as f64 + kappa * self.global.at_or_below(c))
                / (self.t as f64 + kappa).max(1.0);
        let v = local.min(1.0);
        if v < self.running || self.steps.at.is_empty() {
            self.running = self.running.min(v);
            self.steps.at.push(c);
            self.steps.q.push(self.running);
        }
        self.t -= at[TARGET as usize];
        self.d -= at[MASS_SHIFT as usize] + at[RT_SHIFT as usize];
        self.d_id -= at[IDENTITY as usize];
    }

    fn finish(self) -> Option<Steps> {
        (!self.steps.at.is_empty()).then_some(self.steps)
    }
}

struct SortedMerge<'a> {
    runs: Vec<Box<dyn Iterator<Item = f32> + 'a>>,
    heads: Vec<Option<f32>>,
}

impl<'a> SortedMerge<'a> {
    fn new(mut runs: Vec<Box<dyn Iterator<Item = f32> + 'a>>) -> Self {
        let heads = runs.iter_mut().map(|r| r.next()).collect();
        Self { runs, heads }
    }

    fn peek(&self) -> Option<f32> {
        self.heads.iter().flatten().copied().min_by(f32::total_cmp)
    }

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

fn walk_scores(
    mut by_family: [SortedMerge<'_>; 4],
    mut visit: impl FnMut(f32, [usize; 4]) -> bool,
) {
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
        if !visit(c, at) {
            break;
        }
    }
}

pub(crate) struct StepsBuilder {
    walks: [Option<KdecoyWalk>; 4],
}

struct KdecoyWalk {
    e: Estimate,
    k: f64,
    t: usize,
    d: usize,
    d_id: usize,
    running: f64,
    steps: Steps,
    only_at: Option<Vec<f64>>,
}

impl StepsBuilder {
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

    pub fn finish(self) -> [Option<Steps>; 4] {
        self.walks.map(|w| w.map(|w| w.steps))
    }
}

#[cfg(test)]
fn kdecoy_steps(by_family: &[Vec<f32>; 4], e: Estimate) -> Option<Steps> {
    let mut builder = StepsBuilder::new([0, 1, 2, 3].map(|f| by_family[f].len()), None);
    let streams = [0, 1, 2, 3].map(|f| {
        SortedMerge::new(vec![
            Box::new(by_family[f].iter().copied()) as Box<dyn Iterator<Item = f32>>
        ])
    });
    walk_scores(streams, |c, at| {
        builder.visit(c, at);
        true
    });
    let [a, b, c, d] = builder.finish();
    [a, b, c, d].into_iter().nth(e as usize).flatten()
}

pub struct QValues {
    seed: u64,
    strata: usize,
    steps: Vec<[Option<Steps>; 4]>,
    stratum_steps: Vec<Vec<Option<Steps>>>,
    peps: Vec<Option<PepSteps>>,
    stratum_beta: Vec<Vec<f64>>,
    whole: [bool; 4],
    thresholds: Vec<f64>,
    pub counts: FdrCounts,
}

impl QValues {
    pub fn pep(&self, precursor: u32, family: u32, score: Option<f32>) -> Option<f64> {
        let s = score.filter(|s| s.is_finite())?;
        let fold = fold_of(precursor, self.seed);
        let (wrong_peak, identity) = self.peps.get(fold)?.as_ref()?.of(s);
        let stratum = unpack_family(family).1.min(self.strata - 1);
        let beta = self
            .stratum_beta
            .get(fold)
            .and_then(|b| b.get(stratum))
            .copied()
            .unwrap_or(0.0);
        Some((sigmoid(logit(wrong_peak) + beta) + identity).min(1.0))
    }

    pub fn steps(&self) -> usize {
        self.steps
            .iter()
            .flatten()
            .flatten()
            .map(|s| s.at.len())
            .sum()
    }

    pub fn pooled_steps(&self) -> usize {
        self.steps
            .iter()
            .filter_map(|f| f[Estimate::Pooled as usize].as_ref())
            .map(|s| s.at.len())
            .sum()
    }

    pub fn q(&self, e: Estimate, precursor: u32, family: u32, score: Option<f32>) -> Option<f64> {
        assert!(
            self.whole[e as usize],
            "{e:?} q-values were kept at thresholds only"
        );
        self.step_q(e, precursor, family, score)
    }

    pub fn accepts(
        &self,
        e: Estimate,
        precursor: u32,
        family: u32,
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

    fn step_q(&self, e: Estimate, precursor: u32, family: u32, score: Option<f32>) -> Option<f64> {
        let s = score.filter(|s| s.is_finite())?;
        let (family, stratum) = unpack_family(family);
        if !e.keeps(family) {
            return None;
        }
        let fold = fold_of(precursor, self.seed);
        if e == Estimate::Pooled {
            if let Some(Some(steps)) = self.stratum_steps_of(fold, stratum) {
                return Some(steps.of_member(s));
            }
        }
        let steps = self.steps[fold][e as usize].as_ref()?;
        Some(steps.of_member(s))
    }

    fn stratum_steps_of(&self, fold: usize, stratum: usize) -> Option<&Option<Steps>> {
        self.stratum_steps
            .get(fold)
            .and_then(|v| v.get(stratum.min(self.strata - 1)))
    }

    pub fn pseudo_q(&self, precursor: u32, family: u32, score: Option<f32>) -> Option<f64> {
        let s = score?;
        let fold = fold_of(precursor, self.seed);
        if let Some(Some(steps)) = self.stratum_steps_of(fold, unpack_family(family).1) {
            return Some(steps.at_or_below(s));
        }
        let steps = self.steps[fold][Estimate::Pooled as usize].as_ref()?;
        Some(steps.at_or_below(s))
    }
}

#[cfg(test)]
pub struct TransferFdr {
    pub result: FdrResult,
    pub counts: FdrCounts,
    models: FoldModels,
    q: QValues,
}

#[cfg(test)]
pub fn transfer_fdr(rows: &[FdrRow], seed: u64) -> TransferFdr {
    let models = train_folds(rows, seed);
    let precursors: Vec<u32> = rows.iter().map(|r| r.precursor).collect();
    let xs: Vec<&[f32]> = rows.iter().map(|r| r.x.as_slice()).collect();
    let mut scores = Vec::new();
    models.score(&precursors, &xs, &mut scores);
    let mut fs = FoldScores::new(seed);
    for (r, s) in rows.iter().zip(&scores) {
        fs.add(r.precursor, r.family.into(), *s);
    }
    let q = fs.finish(&models.failed_folds);
    let each = |e: Estimate| -> Vec<Option<f64>> {
        rows.iter()
            .zip(&scores)
            .map(|(r, s)| q.q(e, r.precursor, r.family.into(), *s))
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
    pub fn pseudo_q(&self, precursor: u32, x: &[f32]) -> Option<f64> {
        let mut s = Vec::new();
        self.models.score(&[precursor], &[x], &mut s);
        self.q.pseudo_q(precursor, 0, s[0])
    }
}

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
    let ranked = model::fdr::ranked_q::<(usize, usize, usize)>(
        order.len(),
        |(t, d, d_id), rank| match items[order[rank]].2 {
            TARGET => *t += 1,
            IDENTITY => *d_id += 1,
            _ => *d += 1,
        },
        |start, rank| items[order[rank]].1 == items[order[start]].1,
        model::fdr::Eval::GroupEnd,
        |&(t, d, d_id)| {
            let tf = t.max(1) as f64;
            ((d + 1) as f64 / (k * tf) + d_id as f64 / tf).min(1.0)
        },
        1.0,
    );
    for (&j, v) in order.iter().zip(ranked) {
        q[items[j].0] = Some(v);
    }
}

#[cfg(test)]
#[path = "transfers/tests/transfers.rs"]
mod tests;
