use crate::scores::descending_order;

const NUM_LAMBDA: usize = 100;
const MAX_LAMBDA: f64 = 0.5;
const NUM_BOOT: usize = 100;
const MAX_BOOT_DRAW: usize = 1000;
const LCG_MULT: u64 = 279_470_273;
const LCG_MOD: u64 = 4_294_967_291;

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

const DEGREE: usize = 3;
const MAX_KNOTS: usize = 200;
const RIDGE_LAMBDA: f64 = 1e-4;
const SMOOTH_LAMBDA: f64 = 1e-3;
const PEP_FLOOR: f64 = 1e-10;

struct Dense {
    n: usize,
    v: Vec<f64>,
}

impl Dense {
    fn new(n: usize) -> Dense {
        Dense {
            n,
            v: vec![0.0; n * n],
        }
    }

    fn at(&self, i: usize, j: usize) -> f64 {
        self.v[i * self.n + j]
    }

    fn add(&mut self, i: usize, j: usize, x: f64) {
        self.v[i * self.n + j] += x;
    }

    fn times(&self, x: &[f64]) -> Vec<f64> {
        (0..self.n)
            .map(|i| {
                let mut s = 0.0;
                for (g, xj) in self.v[i * self.n..(i + 1) * self.n].iter().zip(x) {
                    s += g * xj;
                }
                s
            })
            .collect()
    }
}

fn solve_on_subset(g: &Dense, free: &[usize], rhs: &[f64]) -> Option<Vec<f64>> {
    let m = free.len();
    let mut l = vec![0.0f64; m * m];
    for j in 0..m {
        let mut s = g.at(free[j], free[j]);
        for k in 0..j {
            s -= l[j * m + k] * l[j * m + k];
        }
        if s.partial_cmp(&0.0) != Some(std::cmp::Ordering::Greater) {
            return None;
        }
        let ljj = s.sqrt();
        l[j * m + j] = ljj;
        for i in j + 1..m {
            let mut t = g.at(free[i], free[j]);
            for k in 0..j {
                t -= l[i * m + k] * l[j * m + k];
            }
            l[i * m + j] = t / ljj;
        }
    }
    let mut y = vec![0.0f64; m];
    for i in 0..m {
        let mut t = rhs[i];
        for k in 0..i {
            t -= l[i * m + k] * y[k];
        }
        y[i] = t / l[i * m + i];
    }
    for i in (0..m).rev() {
        let mut t = y[i];
        for k in i + 1..m {
            t -= l[k * m + i] * y[k];
        }
        y[i] = t / l[i * m + i];
    }
    Some(y)
}

fn solve_box_qp(g: &Dense, h: &[f64], lb: &[f64], ub: &[f64]) -> (Vec<f64>, usize, bool) {
    let p = h.len();
    let mut x: Vec<f64> = (0..p).map(|i| 0.0f64.clamp(lb[i], ub[i])).collect();
    let mut at_bound: Vec<bool> = (0..p).map(|i| x[i] == lb[i] || x[i] == ub[i]).collect();
    let max_iters = 10 * p + 100;
    let mut it = 0usize;
    let mut converged = false;
    while it < max_iters {
        let free: Vec<usize> = (0..p).filter(|&i| !at_bound[i]).collect();
        let rhs: Vec<f64> = free
            .iter()
            .map(|&i| {
                let mut t = h[i];
                for j in 0..p {
                    if at_bound[j] {
                        t -= g.at(i, j) * x[j];
                    }
                }
                t
            })
            .collect();
        let Some(xf) = solve_on_subset(g, &free, &rhs) else {
            break;
        };
        let mut alpha = 1.0f64;
        let mut blocking = None;
        for (a, &i) in free.iter().enumerate() {
            let step = xf[a] - x[i];
            if step < 0.0 && xf[a] < lb[i] {
                let r = (lb[i] - x[i]) / step;
                if r < alpha {
                    alpha = r;
                    blocking = Some(i);
                }
            } else if step > 0.0 && xf[a] > ub[i] {
                let r = (ub[i] - x[i]) / step;
                if r < alpha {
                    alpha = r;
                    blocking = Some(i);
                }
            }
        }
        for (a, &i) in free.iter().enumerate() {
            x[i] = if alpha == 1.0 {
                xf[a]
            } else {
                x[i] + alpha * (xf[a] - x[i])
            };
        }
        if let Some(b) = blocking {
            x[b] = if x[b] - lb[b] < ub[b] - x[b] {
                lb[b]
            } else {
                ub[b]
            };
            at_bound[b] = true;
            it += 1;
            continue;
        }
        let gx = g.times(&x);
        let mut release = None;
        let mut most_inward = 1e-12;
        for i in (0..p).filter(|&i| at_bound[i]) {
            let grad = 2.0 * (gx[i] - h[i]);
            let inward = if x[i] == lb[i] { -grad } else { grad };
            if inward > most_inward {
                most_inward = inward;
                release = Some(i);
            }
        }
        match release {
            Some(r) => at_bound[r] = false,
            None => {
                converged = true;
                break;
            }
        }
        it += 1;
    }
    (x, it, converged)
}

fn span_values(x: f64, knots: &[f64], cum: &mut [f64; DEGREE + 1]) -> usize {
    let d = DEGREE;
    let n = knots.len();
    let t = knots;
    let k = if x >= t[n - d - 1] {
        n - d - 2
    } else if x <= t[d] {
        d
    } else {
        let (mut lo, mut hi) = (d, n - d - 1);
        while lo < hi {
            let mid = (lo + hi) / 2;
            if t[mid + 1] <= x {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        lo
    };
    let mut b = [0.0f64; DEGREE + 1];
    b[0] = 1.0;
    for j in 1..=d {
        let mut saved = 0.0;
        for r in 0..j {
            let left = x - t[k + 1 - j + r];
            let right = t[k + 1 + r] - x;
            let denom = left + right;
            let temp = if denom > 0.0 { b[r] / denom } else { 0.0 };
            b[r] = saved + right * temp;
            saved = left * temp;
        }
        b[j] = saved;
    }
    cum[d] = b[d];
    for r in (0..d).rev() {
        cum[r] = cum[r + 1] + b[r];
    }
    k
}

fn unit_interval(x: &[f64]) -> Result<Vec<f64>, String> {
    if x.is_empty() {
        return Ok(Vec::new());
    }
    let lo = x.iter().copied().fold(f64::INFINITY, f64::min);
    let hi = x.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let span = hi - lo;
    if span
        .partial_cmp(&f64::MAX)
        .is_none_or(|o| o == std::cmp::Ordering::Greater)
    {
        return Err(format!(
            "ERROR: Cannot estimate PEPs by I-spline regression: the scores range from {} to {}, and their span is not a finite number. Check the input for scores that are too large in magnitude.",
            crate::write::fmt_g(-hi, 6),
            crate::write::fmt_g(-lo, 6)
        ));
    }
    if span.partial_cmp(&0.0) != Some(std::cmp::Ordering::Greater) {
        return Ok(vec![0.5; x.len()]);
    }
    Ok(x.iter().map(|&v| (v - lo) / span).collect())
}

fn default_knots(sorted: &[f64]) -> Vec<f64> {
    let n = sorted.len();
    if n == 0 {
        return Vec::new();
    }
    let lo = sorted[0];
    let mut hi = sorted[n - 1];
    if hi <= lo {
        hi = lo + 1e-6;
    }
    let num_internal = MAX_KNOTS.min((n as f64).sqrt() as usize);
    let mut internal = Vec::with_capacity(num_internal);
    let mut prev = lo;
    for k in 1..=num_internal {
        let q = k as f64 / (num_internal + 1) as f64;
        let idx = ((q * (n - 1) as f64) as usize).min(n - 1);
        let v = sorted[idx];
        if v > lo + 1e-12 && v < hi - 1e-12 && v > prev + 1e-12 {
            internal.push(v);
            prev = v;
        }
    }
    let mut knots = vec![lo; DEGREE + 1];
    knots.extend(internal);
    knots.extend(std::iter::repeat_n(hi, DEGREE + 1));
    knots
}

fn ispline_fit(x: &[f64], y: &[f64], clip_lo: f64, clip_hi: f64) -> Result<Vec<f64>, String> {
    let n = x.len();
    let neg: Vec<f64> = x.iter().map(|&v| -v).collect();
    let xn = unit_interval(&neg)?;
    drop(neg);
    let knots = if xn
        .windows(2)
        .all(|p| p[1].partial_cmp(&p[0]) != Some(std::cmp::Ordering::Less))
    {
        default_knots(&xn)
    } else {
        let mut sorted = xn.clone();
        sorted.sort_by(f64::total_cmp);
        default_knots(&sorted)
    };
    let d = DEGREE;
    let n_bspline = knots.len() as i64 - d as i64 - 1;
    let splines = n_bspline > 1;
    let p = if splines { n_bspline as usize } else { 1 };
    let mut g = Dense::new(p);
    let mut h = vec![0.0f64; p];
    let mut cum = [0.0f64; DEGREE + 1];
    if !splines {
        let mut y_sum = 0.0;
        for &v in y {
            y_sum += v;
        }
        g.add(0, 0, n as f64);
        h[0] = y_sum;
    } else {
        let spans = knots.len();
        let dd = d + 1;
        let mut count = vec![0.0f64; spans];
        let mut y_sum = vec![0.0f64; spans];
        let mut cum_sum = vec![0.0f64; spans * dd];
        let mut y_cum_sum = vec![0.0f64; spans * dd];
        let mut cum_cum_sum = vec![0.0f64; spans * dd * dd];
        for (&xi, &yi) in xn.iter().zip(y) {
            let k = span_values(xi, &knots, &mut cum);
            count[k] += 1.0;
            y_sum[k] += yi;
            for t in 1..=d {
                cum_sum[k * dd + t] += cum[t];
                y_cum_sum[k * dd + t] += yi * cum[t];
                for u in t..=d {
                    cum_cum_sum[k * dd * dd + t * dd + u] += cum[t] * cum[u];
                }
            }
        }
        for k in d..spans {
            if count[k] == 0.0 {
                continue;
            }
            let a = k - d;
            let cs = &cum_sum[k * dd..(k + 1) * dd];
            let ycs = &y_cum_sum[k * dd..(k + 1) * dd];
            let ccs = &cum_cum_sum[k * dd * dd..(k + 1) * dd * dd];
            for (r, hr) in h.iter_mut().enumerate().take(a + 1) {
                *hr += y_sum[k];
                for c in 0..=a {
                    g.add(r, c, count[k]);
                }
                for (t, &cst) in cs.iter().enumerate().skip(1) {
                    let c = a + t;
                    g.add(r, c, cst);
                    g.add(c, r, cst);
                }
            }
            for t in 1..=d {
                let r = a + t;
                h[r] += ycs[t];
                for u in t..=d {
                    let c = a + u;
                    g.add(r, c, ccs[t * dd + u]);
                    if c != r {
                        g.add(c, r, ccs[t * dd + u]);
                    }
                }
            }
        }
    }
    let inv_n = 1.0 / n.max(1) as f64;
    for v in g.v.iter_mut() {
        *v *= inv_n;
    }
    for v in h.iter_mut() {
        *v *= inv_n;
    }
    for j in 0..p {
        g.add(j, j, RIDGE_LAMBDA);
    }
    let n_ispline = p - 1;
    if n_ispline > 2 {
        let coef = [1.0f64, -2.0, 1.0];
        for di in 0..n_ispline - 2 {
            for (r, &cr) in coef.iter().enumerate() {
                for (c, &cc) in coef.iter().enumerate() {
                    g.add(1 + di + r, 1 + di + c, SMOOTH_LAMBDA * cr * cc);
                }
            }
        }
    }
    let mut lb = vec![0.0f64; p];
    let ub = vec![f64::INFINITY; p];
    lb[0] = f64::NEG_INFINITY;
    let (sol, iters, converged) = solve_box_qp(&g, &h, &lb, &ub);
    if !converged {
        eprintln!(
            "Warning: the I-spline regression of the PEPs did not converge ({iters} iterations, {p} coefficients); the PEPs may be inaccurate."
        );
    }
    if !splines {
        return Ok(vec![sol[0].max(clip_lo).min(clip_hi); n]);
    }
    let nb = n_bspline as usize;
    let mut prefix = vec![0.0f64; nb];
    for j in 1..nb {
        prefix[j] = prefix[j - 1] + sol[j];
    }
    Ok(xn
        .iter()
        .map(|&xi| {
            let k = span_values(xi, &knots, &mut cum);
            let a = k - d;
            let mut v = sol[0] + prefix[a];
            for t in 1..=d {
                v += cum[t] * sol[a + t];
            }
            v.max(clip_lo).min(clip_hi)
        })
        .collect())
}

fn interpolate(q: f64, q1: f64, q2: f64, pep1: f64, pep2: f64) -> f64 {
    if (q2 - q1).abs() <= 1e-12 {
        return pep2.clamp(0.0, 1.0);
    }
    (pep1 + (q - q1) * (pep2 - pep1) / (q2 - q1)).clamp(0.0, 1.0)
}

pub fn pep_from_q(
    order: &[usize],
    scores: &[f64],
    q: &[f64],
    is_decoy: &[bool],
) -> Result<Vec<f64>, String> {
    let targets: Vec<usize> = order.iter().copied().filter(|&k| !is_decoy[k]).collect();
    let mut target_q: Vec<f64> = targets.iter().map(|&k| q[k]).collect();
    let sc: Vec<f64> = targets.iter().map(|&k| scores[k]).collect();
    let pseudo = 0.5 / target_q.len().max(1) as f64;
    let mut prev = 0.0f64;
    let raw: Vec<f64> = target_q
        .iter()
        .enumerate()
        .map(|(i, &qi)| {
            let qn = qi * (i + 1) as f64;
            let r = qn - prev;
            prev = qn;
            r + pseudo
        })
        .collect();
    let mut target_pep = ispline_fit(&sc, &raw, PEP_FLOOR, 1.0)?;
    target_pep.push(1.0);
    target_q.push(1.0);
    let mut pep = vec![0.0f64; scores.len()];
    let (mut l_q, mut l_pep) = (0.0f64, 0.0f64);
    let mut ti = 0usize;
    for &k in order {
        if is_decoy[k] {
            pep[k] = interpolate(q[k], l_q, target_q[ti], l_pep, target_pep[ti]);
        } else {
            pep[k] = target_pep[ti];
            l_pep = target_pep[ti];
            l_q = target_q[ti];
            ti += 1;
        }
    }
    Ok(pep)
}
