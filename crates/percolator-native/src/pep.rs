use crate::scores::{descending_order, tdc_q};

const NO_INTERVALS: usize = 500;
const NUM_LAMBDA: usize = 100;
const MAX_LAMBDA: f64 = 0.5;
const NUM_BOOT: usize = 100;
const MAX_BOOT_DRAW: usize = 1000;
const G_RANGE: f64 = 35.0;
const STEP_EPSILON: f64 = 1e-8;
const WEIGHT_SLOPE: f64 = 1e1;
const SCALE_ALPHA: f64 = 1.0;
const ACCEPTED_Q: f64 = 0.01;
const LCG_MULT: u64 = 279_470_273;
const LCG_MOD: u64 = 4_294_967_291;

fn tao() -> f64 {
    2.0 / (1.0 + 5.0f64.sqrt())
}

struct Lcg(u64);

impl Lcg {
    fn uniform(&mut self) -> f64 {
        self.0 = (self.0 * LCG_MULT) % LCG_MOD;
        self.0 as f64 / (LCG_MOD as f64 + 1.0)
    }
}

fn combined_desc(scores: &[f64], is_decoy: &[bool]) -> (Vec<usize>, Vec<(f64, bool)>) {
    let order = descending_order(scores);
    let combined = order.iter().map(|&i| (scores[i], !is_decoy[i])).collect();
    (order, combined)
}

fn p_values(combined: &[(f64, bool)]) -> Vec<f64> {
    let mut p = Vec::new();
    let mut n_decoys = 1usize;
    let mut pos_same = 0usize;
    let mut neg_same = 0usize;
    for (i, &(s, target)) in combined.iter().enumerate() {
        if target {
            pos_same += 1;
        } else {
            neg_same += 1;
        }
        if i + 1 == combined.len() || s != combined[i + 1].0 {
            for ix in 0..pos_same {
                p.push(n_decoys as f64 + (neg_same * (ix + 1)) as f64 / (pos_same + 1) as f64);
            }
            n_decoys += neg_same;
            neg_same = 0;
            pos_same = 0;
        }
    }
    let n = n_decoys as f64;
    for v in p.iter_mut() {
        *v /= n;
    }
    p
}

fn count_ge(sorted: &[f64], lambda: f64) -> usize {
    sorted.len() - sorted.partition_point(|&v| v < lambda)
}

pub fn too_good_separation(scores: &[f64], is_decoy: &[bool]) -> bool {
    let (_, combined) = combined_desc(scores, is_decoy);
    let p = p_values(&combined);
    let min_lambda = (1.0 / NUM_LAMBDA as f64) * MAX_LAMBDA;
    count_ge(&p, min_lambda) == 0
}

pub fn estimate_pi0(scores: &[f64], is_decoy: &[bool]) -> f64 {
    let (_, combined) = combined_desc(scores, is_decoy);
    let p = p_values(&combined);
    let n = p.len() as f64;
    let mut lambdas = Vec::new();
    let mut pi0s = Vec::new();
    for ix in 0..=NUM_LAMBDA {
        let lambda = ((ix + 1) as f64 / NUM_LAMBDA as f64) * MAX_LAMBDA;
        let wl = count_ge(&p, lambda) as f64;
        let pi0 = wl / n / (1.0 - lambda);
        if pi0 > 0.0 {
            lambdas.push(lambda);
            pi0s.push(pi0);
        }
    }
    if pi0s.is_empty() {
        return 1.0;
    }
    let min_pi0 = pi0s.iter().copied().fold(f64::INFINITY, f64::min);
    let mut mse = vec![0.0f64; pi0s.len()];
    let mut rng = Lcg(1);
    let num_draw = p.len().min(MAX_BOOT_DRAW);
    let mut boot = Vec::with_capacity(num_draw);
    for _ in 0..NUM_BOOT {
        boot.clear();
        for _ in 0..num_draw {
            let draw = (rng.uniform() * n) as usize;
            boot.push(p[draw]);
        }
        boot.sort_by(f64::total_cmp);
        let nb = boot.len() as f64;
        for (m, &lambda) in mse.iter_mut().zip(&lambdas) {
            let wl = count_ge(&boot, lambda) as f64;
            let pi0_boot = wl / nb / (1.0 - lambda);
            *m += (pi0_boot - min_pi0) * (pi0_boot - min_pi0);
        }
    }
    let mut min_ix = 0usize;
    for (i, &m) in mse.iter().enumerate() {
        if m < mse[min_ix] {
            min_ix = i;
        }
    }
    pi0s[min_ix].clamp(0.0, 1.0)
}

pub(crate) fn mix_max_counts(combined: &[(f64, bool)]) -> (Vec<f64>, Vec<f64>) {
    let mut h_w = Vec::new();
    let mut h_z = Vec::new();
    let (mut cnt_w, mut cnt_z, mut queue) = (0usize, 0usize, 0usize);
    for i in (0..combined.len()).rev() {
        if combined[i].1 {
            cnt_w += 1;
        } else {
            cnt_z += 1;
            queue += 1;
        }
        if i == 0 || combined[i].0 != combined[i - 1].0 {
            for _ in 0..queue {
                h_w.push(cnt_w as f64);
                h_z.push(cnt_z as f64);
            }
            queue = 0;
        }
    }
    (h_w, h_z)
}

struct Bins {
    medians: Vec<f64>,
    negatives: Vec<f64>,
    sizes: Vec<f64>,
}

fn bin_data(combined: &[(f64, bool)], pi0: f64) -> Bins {
    let (h_w, h_z) = if pi0 < 1.0 {
        mix_max_counts(combined)
    } else {
        (Vec::new(), Vec::new())
    };
    let n = combined.len();
    let mut bins = Bins {
        medians: Vec::new(),
        negatives: Vec::new(),
        sizes: Vec::new(),
    };
    let mut e_f1 = 0.0f64;
    let mut bins_left = NO_INTERVALS as i64 - 1;
    let target_size = (n as f64 / NO_INTERVALS as f64).max(1.0);
    let (mut n_z_ge_w, mut sum_n_z_ge_w) = (0i64, 0i64);
    let (mut decoy_queue, mut psms_in_bin, mut bin_start) = (0i64, 0i64, 0i64);
    for i in 0..n {
        if !combined[i].1 {
            n_z_ge_w += 1;
            decoy_queue += 1;
        }
        psms_in_bin += 1;
        if i + 1 == n || combined[i].0 != combined[i + 1].0 {
            if pi0 < 1.0 && decoy_queue > 0 {
                let j = (h_w.len() as i64 - sum_n_z_ge_w - n_z_ge_w) as usize;
                let (cnt_w, cnt_z) = (h_w[j], h_z[j]);
                let mut est = (cnt_w - pi0 * cnt_z) / ((1.0 - pi0) * cnt_z);
                est = if est > 1.0 { 1.0 } else { est };
                est = if est < 0.0 { 0.0 } else { est };
                e_f1 += decoy_queue as f64 * est * (1.0 - pi0);
            }
            decoy_queue = 0;
            if ((n as i64 - bin_start - psms_in_bin) as f64) <= bins_left as f64 * target_size {
                let median = combined[(bin_start + psms_in_bin / 2) as usize].0;
                let negatives = n_z_ge_w as f64 * pi0 + e_f1;
                let size = (psms_in_bin - n_z_ge_w) as f64 + negatives;
                let len = bins.medians.len();
                if len > 0 && bins.medians[len - 1] == median {
                    bins.negatives[len - 1] += negatives;
                    bins.sizes[len - 1] += size;
                } else {
                    bins.medians.push(median);
                    bins.negatives.push(negatives);
                    bins.sizes.push(size);
                }
                bin_start += psms_in_bin;
                bins_left -= 1;
                psms_in_bin = 0;
                e_f1 = 0.0;
                sum_n_z_ge_w += n_z_ge_w;
                n_z_ge_w = 0;
            }
        }
    }
    bins
}

#[derive(Clone, Copy)]
enum Transform {
    Identity,
    Logit { low: f64, high: f64 },
    Log { low: f64 },
}

impl Transform {
    fn of(xs: &[f64]) -> Transform {
        let min_v = xs.iter().copied().fold(f64::INFINITY, f64::min);
        let max_v = xs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        if min_v >= 0.0 && max_v <= 1.0 {
            Transform::Logit {
                low: if min_v > 0.0 { 0.0 } else { 1e-20 },
                high: if max_v < 1.0 { 0.0 } else { 1e-10 },
            }
        } else if min_v >= 0.0 {
            Transform::Log {
                low: if min_v > 0.0 { 0.0 } else { 1e-20 },
            }
        } else {
            Transform::Identity
        }
    }

    fn apply(self, mut xx: f64) -> f64 {
        match self {
            Transform::Identity => xx,
            Transform::Logit { low, high } => {
                if low > 0.0 || high > 0.0 {
                    xx *= 1.0 - high - low;
                    xx += low;
                }
                (xx / (1.0 - xx)).ln()
            }
            Transform::Log { low } => {
                if low > 0.0 {
                    xx += low;
                }
                xx.ln()
            }
        }
    }
}

struct Spline {
    transf: Transform,
    x: Vec<f64>,
    y: Vec<f64>,
    m: Vec<f64>,
    dx: Vec<f64>,
    g: Vec<f64>,
    gnew: Vec<f64>,
    w: Vec<f64>,
    z: Vec<f64>,
    gamma: Vec<f64>,
}

fn qm(dx: &[f64], i: usize, c: usize) -> f64 {
    if i == c {
        1.0 / dx[c]
    } else if i == c + 1 {
        -1.0 / dx[c] - 1.0 / dx[c + 1]
    } else if i == c + 2 {
        1.0 / dx[c + 1]
    } else {
        0.0
    }
}

fn rm(dx: &[f64], c: usize, d: usize) -> f64 {
    let (lo, hi) = if c <= d { (c, d) } else { (d, c) };
    if lo == hi {
        (dx[lo] + dx[lo + 1]) / 3.0
    } else if hi == lo + 1 {
        dx[hi] / 6.0
    } else {
        0.0
    }
}

fn solve_penta(a0: &[f64], a1: &[f64], a2: &[f64], b: &mut [f64]) {
    let k = a0.len();
    let mut d = vec![0.0f64; k];
    let mut l1 = vec![0.0f64; k];
    let mut l2 = vec![0.0f64; k];
    for j in 0..k {
        let mut dj = a0[j];
        if j >= 1 {
            dj -= l1[j - 1] * l1[j - 1] * d[j - 1];
        }
        if j >= 2 {
            dj -= l2[j - 2] * l2[j - 2] * d[j - 2];
        }
        d[j] = dj;
        if j + 1 < k {
            let mut v = a1[j];
            if j >= 1 {
                v -= l2[j - 1] * l1[j - 1] * d[j - 1];
            }
            l1[j] = v / dj;
        }
        if j + 2 < k {
            l2[j] = a2[j] / dj;
        }
    }
    for j in 0..k {
        if j >= 1 {
            b[j] -= l1[j - 1] * b[j - 1];
        }
        if j >= 2 {
            b[j] -= l2[j - 2] * b[j - 2];
        }
    }
    for (bj, dj) in b.iter_mut().zip(&d) {
        *bj /= dj;
    }
    for j in (0..k).rev() {
        if j + 1 < k {
            b[j] -= l1[j] * b[j + 1];
        }
        if j + 2 < k {
            b[j] -= l2[j] * b[j + 2];
        }
    }
}

impl Spline {
    fn new(medians: &[f64], negatives: &[f64], sizes: &[f64]) -> Spline {
        let transf = Transform::of(medians);
        let x: Vec<f64> = medians.iter().map(|&v| transf.apply(v)).collect();
        let n = x.len();
        let dx = (0..n - 1).map(|i| x[i + 1] - x[i]).collect();
        let gnew = negatives
            .iter()
            .zip(sizes)
            .map(|(&y, &m)| {
                let p = (y + 0.05) / (m + 0.1);
                (p / (1.0 - p)).ln()
            })
            .collect();
        Spline {
            transf,
            x,
            y: negatives.to_vec(),
            m: sizes.to_vec(),
            dx,
            g: vec![0.0; n],
            gnew,
            w: vec![0.0; n],
            z: vec![0.5; n],
            gamma: vec![0.0; n - 2],
        }
    }

    fn calc_pzw(&mut self) {
        let epsilon = 1e-15;
        let obs = self.g.iter().zip(&self.y).zip(&self.m);
        for (((&g, &y), &m), (w, z)) in obs.zip(self.w.iter_mut().zip(self.z.iter_mut())) {
            let e = g.exp();
            let p = (e / (1.0 + e)).clamp(epsilon, 1.0 - epsilon);
            *w = (m * p * (1.0 - p)).max(epsilon);
            *z = (g + (y - p * m) / *w).clamp(-G_RANGE, G_RANGE);
        }
    }

    fn irls(&mut self, alpha: f64) {
        let n = self.x.len();
        let k = n - 2;
        let mut iter = 0;
        loop {
            self.g.clone_from(&self.gnew);
            self.calc_pzw();
            let a: Vec<f64> = self.w.iter().map(|&w| (1.0 / w) * alpha).collect();
            let dx = &self.dx;
            let mut a0 = vec![0.0f64; k];
            let mut a1 = vec![0.0f64; k];
            let mut a2 = vec![0.0f64; k];
            for c in 0..k {
                for (off, slot) in [(0usize, &mut a0), (1, &mut a1), (2, &mut a2)] {
                    let d = c + off;
                    if d >= k {
                        continue;
                    }
                    let mut v = 0.0;
                    for (i, &ai) in a.iter().enumerate().take(c + 3).skip(d) {
                        v += qm(dx, i, c) * ai * qm(dx, i, d);
                    }
                    slot[c] = rm(dx, c, d) + v;
                }
            }
            let mut gamma: Vec<f64> = (0..k)
                .map(|c| (c..=c + 2).map(|i| qm(dx, i, c) * self.z[i]).sum())
                .collect();
            solve_penta(&a0, &a1, &a2, &mut gamma);
            let mut ss = 0.0;
            for (i, &ai) in a.iter().enumerate() {
                let lo = i.saturating_sub(2);
                let hi = i.min(k - 1);
                let mut qg = 0.0;
                for (c, &gc) in gamma.iter().enumerate().take(hi + 1).skip(lo) {
                    qg += ai * qm(dx, i, c) * gc;
                }
                let v = (self.z[i] - qg).clamp(-G_RANGE, G_RANGE);
                self.gnew[i] = v;
                let diff = self.g[i] - v;
                ss += diff * diff;
            }
            self.gamma = gamma;
            let step = ss.sqrt() / n as f64;
            iter += 1;
            if !(step > STEP_EPSILON && iter < 20) {
                break;
            }
        }
        self.g.clone_from(&self.gnew);
    }

    fn evaluate_slope(&mut self, alpha: f64) -> f64 {
        self.irls(alpha);
        let n = self.g.len();
        let mut mixg = 1usize;
        let mut maxg = self.g[1];
        for (ix, &gi) in self.g.iter().enumerate().take(n - 1).skip(1) {
            if gi >= maxg {
                maxg = gi;
                mixg = ix;
            }
        }
        let mut max_slope = -10e6;
        for ix in mixg + 1..n.saturating_sub(2) {
            let slope = self.g[ix - 1] - self.g[ix];
            if slope > max_slope {
                max_slope = slope;
            }
        }
        max_slope * WEIGHT_SLOPE + alpha
    }

    fn alpha_search(&mut self, p1: f64, p2: f64) -> f64 {
        let t = tao();
        let (mut min_p, mut max_p, mut p1, mut p2) = (0.0f64, 1.0f64, p1, p2);
        let mut cv2 = self.evaluate_slope(-SCALE_ALPHA * p2.ln());
        let mut cv1 = self.evaluate_slope(-SCALE_ALPHA * p1.ln());
        loop {
            let old;
            if cv2 < cv1 {
                min_p = p1;
                p1 = p2;
                p2 = min_p + t * (max_p - min_p);
                old = cv1;
                cv1 = cv2;
                cv2 = self.evaluate_slope(-SCALE_ALPHA * p2.ln());
            } else {
                max_p = p2;
                p2 = p1;
                p1 = min_p + (1.0 - t) * (max_p - min_p);
                old = cv2;
                cv2 = cv1;
                cv1 = self.evaluate_slope(-SCALE_ALPHA * p1.ln());
            }
            let best = if cv2 < cv1 { cv2 } else { cv1 };
            if (old - best) / old < 1e-5 || (p2 - p1).abs() < 1e-10 {
                return if cv1 < cv2 {
                    -SCALE_ALPHA * p1.ln()
                } else {
                    -SCALE_ALPHA * p2.ln()
                };
            }
        }
    }

    fn fit(&mut self) {
        let t = tao();
        let alpha = self.alpha_search(1.0 - t, t);
        self.irls(alpha);
    }

    fn predict(&self, xx: f64) -> f64 {
        let xx = self.transf.apply(xx);
        let (x, g, gamma) = (&self.x, &self.g, &self.gamma);
        let n = x.len();
        let rix = x.partition_point(|&v| v < xx);
        if rix == n {
            let derl = (g[n - 1] - g[n - 2]) / (x[n - 1] - x[n - 2])
                + (x[n - 1] - x[n - 2]) / 6.0 * gamma[n - 3];
            return g[n - 1] + (xx - x[n - 1]) * derl;
        }
        if x[rix] == xx {
            return g[rix];
        }
        if rix > 0 {
            let dr = x[rix] - xx;
            let dl = xx - x[rix - 1];
            let gamr = if rix < n - 1 { gamma[rix - 1] } else { 0.0 };
            let gaml = if rix > 1 { gamma[rix - 2] } else { 0.0 };
            let h = x[rix] - x[rix - 1];
            return (dl * g[rix] + dr * g[rix - 1]) / h
                - dl * dr / 6.0 * ((1.0 + dl / h) * gamr + (1.0 + dr / h) * gaml);
        }
        let derr = (g[1] - g[0]) / (x[1] - x[0]) - (x[1] - x[0]) / 6.0 * gamma[0];
        g[0] - (x[0] - xx) * derr
    }
}

fn bins_for(combined: &[(f64, bool)], use_pi0: bool, pi0: f64) -> Bins {
    let mut bins = if use_pi0 {
        bin_data(combined, pi0)
    } else {
        let ascending: Vec<(f64, bool)> = combined.iter().rev().copied().collect();
        bin_data(&ascending, pi0)
    };
    if bins.medians.len() > 1 && bins.medians[0] > bins.medians[bins.medians.len() - 1] {
        bins.medians.reverse();
        bins.negatives.reverse();
        bins.sizes.reverse();
    }
    bins
}

fn finish(mut peps: Vec<f64>) -> Vec<f64> {
    if peps.is_empty() {
        return peps;
    }
    let max = peps.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let top = 1.0f64.min(max.exp());
    let mut crap = false;
    for pep in peps.iter_mut() {
        if crap {
            *pep = top;
            continue;
        }
        *pep = pep.exp();
        if *pep >= top {
            *pep = top;
            crap = true;
        }
    }
    let mut run = f64::INFINITY;
    for pep in peps.iter_mut().rev() {
        run = if run > *pep { *pep } else { run };
        *pep = run;
    }
    peps
}

fn estimate_pep(combined: &[(f64, bool)], use_pi0: bool, pi0: f64) -> Option<Vec<f64>> {
    let bins = bins_for(combined, use_pi0, pi0);
    if bins.medians.len() < 4 {
        return None;
    }
    let mut spline = Spline::new(&bins.medians, &bins.negatives, &bins.sizes);
    spline.fit();
    let raw: Vec<f64> = combined.iter().map(|&(s, _)| spline.predict(s)).collect();
    let peps = finish(raw);
    if peps.iter().all(|p| p.is_finite()) {
        Some(peps)
    } else {
        None
    }
}

fn pava_non_increasing(values: &[f64], weights: &[f64]) -> Vec<f64> {
    let mut blocks: Vec<(f64, f64, usize)> = Vec::new();
    for (&v, &w) in values.iter().zip(weights) {
        blocks.push((v, w.max(f64::MIN_POSITIVE), 1));
        while blocks.len() > 1 {
            let (v2, w2, n2) = blocks[blocks.len() - 1];
            let (v1, w1, n1) = blocks[blocks.len() - 2];
            if v1 >= v2 {
                break;
            }
            blocks.pop();
            blocks.pop();
            let w = w1 + w2;
            blocks.push(((v1 * w1 + v2 * w2) / w, w, n1 + n2));
        }
    }
    let mut out = Vec::with_capacity(values.len());
    for (v, _, n) in blocks {
        out.resize(out.len() + n, v);
    }
    out
}

fn isotonic_pep(combined: &[(f64, bool)], use_pi0: bool, pi0: f64) -> Vec<f64> {
    let bins = bins_for(combined, use_pi0, pi0);
    if bins.medians.is_empty() {
        return vec![1.0; combined.len()];
    }
    let odds: Vec<f64> = bins
        .negatives
        .iter()
        .zip(&bins.sizes)
        .map(|(&y, &m)| {
            let p = (y + 0.05) / (m + 0.1);
            p / (1.0 - p)
        })
        .collect();
    let fitted = pava_non_increasing(&odds, &bins.sizes);
    let x = &bins.medians;
    let n = x.len();
    let peps: Vec<f64> = combined
        .iter()
        .map(|&(s, _)| {
            let r = x.partition_point(|&v| v < s);
            let v = if r == 0 {
                fitted[0]
            } else if r == n {
                fitted[n - 1]
            } else {
                let t = (s - x[r - 1]) / (x[r] - x[r - 1]);
                fitted[r - 1] + t * (fitted[r] - fitted[r - 1])
            };
            v.clamp(0.0, 1.0)
        })
        .collect();
    let mut run = f64::INFINITY;
    let mut out = peps;
    for pep in out.iter_mut().rev() {
        run = run.min(*pep);
        *pep = run;
    }
    out
}

pub fn pep_is_degenerate(pep: &[f64], q: &[f64], is_decoy: &[bool]) -> bool {
    let top = pep.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mut accepted = 0usize;
    for ((&p, &qi), &d) in pep.iter().zip(q).zip(is_decoy) {
        if d || qi > ACCEPTED_Q {
            continue;
        }
        if p < top {
            return false;
        }
        accepted += 1;
    }
    accepted > 0
}

fn in_input_order(order: &[usize], desc: Vec<f64>) -> Vec<f64> {
    let mut out = vec![0.0f64; desc.len()];
    for (&i, v) in order.iter().zip(desc) {
        out[i] = v;
    }
    out
}

fn guarded(scores: &[f64], is_decoy: &[bool], use_pi0: bool, pi0: f64) -> Vec<f64> {
    let (order, combined) = combined_desc(scores, is_decoy);
    let q = tdc_q(scores, is_decoy, false);
    if let Some(desc) = estimate_pep(&combined, use_pi0, pi0) {
        let pep = in_input_order(&order, desc);
        if !pep_is_degenerate(&pep, &q, is_decoy) {
            return pep;
        }
    }
    in_input_order(&order, isotonic_pep(&combined, use_pi0, pi0))
}

pub fn pep_tdc_fit(scores: &[f64], is_decoy: &[bool]) -> Option<Vec<f64>> {
    let (order, combined) = combined_desc(scores, is_decoy);
    estimate_pep(&combined, false, 1.0).map(|desc| in_input_order(&order, desc))
}

pub fn pep_tdc(scores: &[f64], is_decoy: &[bool]) -> Vec<f64> {
    guarded(scores, is_decoy, false, 1.0)
}

pub fn pep_mix_max(scores: &[f64], is_decoy: &[bool], pi0: f64) -> Vec<f64> {
    guarded(scores, is_decoy, true, pi0)
}
