use crate::normalize::Norm;
use crate::store::Store;

pub struct Rows {
    pub data: Vec<f64>,
    pub nrows: usize,
    pub ncols: usize,
}

impl Rows {
    pub fn gather(store: &Store, norm: &Norm, idx: &[usize]) -> Rows {
        let nf = store.nfeat();
        let ncols = nf + 1;
        let mut data = vec![0.0f64; idx.len() * ncols];
        for r in 0..idx.len() {
            data[r * ncols + nf] = 1.0;
        }
        if !idx.is_empty() && nf > 0 {
            let mut order: Vec<usize> = (0..idx.len()).collect();
            order.sort_by_key(|&p| idx[p]);
            let mut next = 0usize;
            store.for_each_block(|start, block| {
                let n = if nf == 0 { 0 } else { block.len() / nf };
                let end = start + n;
                while next < order.len() && idx[order[next]] < end {
                    let p = order[next];
                    let r = idx[p] - start;
                    let src = &block[r * nf..(r + 1) * nf];
                    let dst = &mut data[p * ncols..p * ncols + nf];
                    for (((d, &x), &m), &sd) in
                        dst.iter_mut().zip(src).zip(&norm.mean).zip(&norm.inv_sd)
                    {
                        *d = (f64::from(x) - m) * sd;
                    }
                    next += 1;
                }
            });
            assert_eq!(
                next,
                order.len(),
                "Rows::gather: row {} not in a store of {} rows",
                idx[order[next.min(order.len() - 1)]],
                store.len()
            );
        }
        Rows {
            data,
            nrows: idx.len(),
            ncols,
        }
    }

    pub fn row(&self, i: usize) -> &[f64] {
        &self.data[i * self.ncols..(i + 1) * self.ncols]
    }
}

pub struct MfnOpts {
    pub lambda: f64,
    pub epsilon: f64,
    pub cgitermax: usize,
    pub mfnitermax: usize,
    pub big_epsilon: f64,
    pub small_cgitermax: usize,
    pub relative_stop_eps: f64,
}

impl Default for MfnOpts {
    fn default() -> Self {
        MfnOpts {
            lambda: 1.0,
            epsilon: 1e-7,
            cgitermax: 10000,
            mfnitermax: 50,
            big_epsilon: 0.01,
            small_cgitermax: 10,
            relative_stop_eps: 1e-9,
        }
    }
}

struct Train<'a> {
    x: &'a Rows,
    idx_pos: &'a [usize],
    idx_neg: &'a [usize],
    cpos: f64,
    cneg: f64,
}

impl Train<'_> {
    fn len(&self) -> usize {
        self.idx_neg.len() + self.idx_pos.len()
    }

    fn row(&self, i: usize) -> &[f64] {
        let n = self.idx_neg.len();
        if i < n {
            self.x.row(self.idx_neg[i])
        } else {
            self.x.row(self.idx_pos[i - n])
        }
    }

    fn y(&self, i: usize) -> f64 {
        if i < self.idx_neg.len() {
            -1.0
        } else {
            1.0
        }
    }

    fn c(&self, i: usize) -> f64 {
        if i < self.idx_neg.len() {
            self.cneg
        } else {
            self.cpos
        }
    }
}

fn ddot(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len();
    let mut s = 0.0;
    let mut i = 0usize;
    while i + 4 < n {
        s += x[i] * y[i]
            + x[i + 1] * y[i + 1]
            + x[i + 2] * y[i + 2]
            + x[i + 3] * y[i + 3]
            + x[i + 4] * y[i + 4];
        i += 5;
    }
    while i < n {
        s += x[i] * y[i];
        i += 1;
    }
    s
}

fn daxpy(a: f64, x: &[f64], y: &mut [f64]) {
    if a != 0.0 {
        for (yi, &xi) in y.iter_mut().zip(x) {
            *yi += a * xi;
        }
    }
}

fn dscal(a: f64, x: &mut [f64]) {
    for v in x.iter_mut() {
        *v *= a;
    }
}

fn gemv_t(row: &[f64], p: &[f64]) -> f64 {
    let mut t = 0.0;
    for (&a, &b) in row.iter().zip(p) {
        t += a * b;
    }
    t
}

struct CgWork {
    z: Vec<f64>,
    q: Vec<f64>,
    r: Vec<f64>,
    p: Vec<f64>,
}

#[allow(clippy::too_many_arguments)]
fn cgls(
    t: &Train,
    lambda: f64,
    cgitermax: usize,
    epsilon: f64,
    subset: &[usize],
    beta: &mut [f64],
    o: &mut [f64],
    ws: &mut CgWork,
) -> bool {
    let active = subset.len();
    let z = &mut ws.z[..active];
    let q = &mut ws.q[..active];
    let r = &mut ws.r;
    let p = &mut ws.p;
    r.fill(0.0);
    for (zi, &ii) in z.iter_mut().zip(subset) {
        *zi = t.c(ii) * (t.y(ii) - o[ii]);
        daxpy(*zi, t.row(ii), r);
    }
    daxpy(-lambda, beta, r);
    p.copy_from_slice(r);
    let mut omega1 = ddot(r, r);
    let mut omega_p = omega1;
    let epsilon2 = epsilon * epsilon;
    let mut cgiter = 0usize;
    while cgiter < cgitermax {
        cgiter += 1;
        let mut omega_q = 0.0;
        for (qi, &ii) in q.iter_mut().zip(subset) {
            *qi = gemv_t(t.row(ii), p);
            omega_q += t.c(ii) * *qi * *qi;
        }
        let gamma = omega1 / (lambda * omega_p + omega_q);
        let inv_omega2 = 1.0 / omega1;
        r.copy_from_slice(beta);
        dscal(-lambda, r);
        daxpy(gamma, p, beta);
        dscal(gamma, q);
        for ((zi, &qi), &ii) in z.iter_mut().zip(q.iter()).zip(subset) {
            o[ii] += qi;
            *zi -= t.c(ii) * qi;
            daxpy(*zi, t.row(ii), r);
        }
        let omega_z = ddot(z, z);
        omega1 = ddot(r, r);
        if omega1 < epsilon2 * omega_z {
            return true;
        }
        let scale = omega1 * inv_omega2;
        dscal(scale, p);
        daxpy(1.0, r, p);
        omega_p = ddot(p, p);
    }
    false
}

fn line_search(
    t: &Train,
    w: &[f64],
    w_bar: &[f64],
    lambda: f64,
    o: &[f64],
    o_bar: &[f64],
    deltas: &mut Vec<(f64, usize, f64)>,
) -> f64 {
    let mut omega_l = 0.0;
    let mut omega_r = 0.0;
    for i in (0..w.len()).rev() {
        let diff = w_bar[i] - w[i];
        omega_l += w[i] * diff;
        omega_r += w_bar[i] * diff;
    }
    let mut l = lambda * omega_l;
    let mut r = lambda * omega_r;
    deltas.clear();
    for i in 0..t.len() {
        let y = t.y(i);
        let diff = y * (o_bar[i] - o[i]);
        if y * o[i] < 1.0 {
            let d2 = t.c(i) * (o_bar[i] - o[i]);
            l += (o[i] - y) * d2;
            r += (o_bar[i] - y) * d2;
            if diff > 0.0 {
                deltas.push(((1.0 - y * o[i]) / diff, i, -1.0));
            }
        } else if diff < 0.0 {
            deltas.push(((1.0 - y * o[i]) / diff, i, 1.0));
        }
    }
    deltas.sort_by(|a, b| a.0.total_cmp(&b.0));
    for &(delta, ii, s) in deltas.iter() {
        let delta_prime = l + delta * (r - l);
        if delta_prime >= 0.0 {
            break;
        }
        let diff = s * t.c(ii) * (o_bar[ii] - o[ii]);
        l += diff * (o[ii] - t.y(ii));
        r += diff * (o_bar[ii] - t.y(ii));
    }
    -l / (r - l)
}

fn partition(t: &Train, o: &[f64], subset: &mut [usize], f_start: f64) -> (usize, f64) {
    let m = t.len();
    let mut active = 0usize;
    let mut inactive = m;
    let mut f = f_start;
    for (i, &oi) in o.iter().enumerate() {
        let diff = 1.0 - t.y(i) * oi;
        if diff > 0.0 {
            subset[active] = i;
            active += 1;
            f += 0.5 * t.c(i) * diff * diff;
        } else {
            inactive -= 1;
            subset[inactive] = i;
        }
    }
    (active, f)
}

pub fn l2_svm_mfn(
    x: &Rows,
    idx_pos: &[usize],
    idx_neg: &[usize],
    cpos: f64,
    cneg: f64,
    opts: &MfnOpts,
) -> Vec<f64> {
    let t = Train {
        x,
        idx_pos,
        idx_neg,
        cpos,
        cneg,
    };
    let n = x.ncols;
    let n0 = n - 1;
    let m = t.len();
    let lambda = opts.lambda;
    let mut epsilon = opts.big_epsilon;
    let mut cgitermax = opts.small_cgitermax;
    let mut w = vec![0.0f64; n];
    let mut o = vec![0.0f64; m];
    let mut subset = vec![0usize; m];
    let (mut active, mut f) = partition(&t, &o, &mut subset, 0.5 * lambda * ddot(&w, &w));
    let mut w_bar = vec![0.0f64; n];
    let mut o_bar = vec![0.0f64; m];
    let mut ws = CgWork {
        z: vec![0.0; m],
        q: vec![0.0; m],
        r: vec![0.0; n],
        p: vec![0.0; n],
    };
    let mut deltas: Vec<(f64, usize, f64)> = Vec::with_capacity(m);
    let mut ini = false;
    let mut iter = 0usize;
    while iter < opts.mfnitermax {
        iter += 1;
        w_bar.copy_from_slice(&w);
        o_bar.copy_from_slice(&o);
        let opt = cgls(
            &t,
            lambda,
            cgitermax,
            epsilon,
            &subset[..active],
            &mut w_bar,
            &mut o_bar,
            &mut ws,
        );
        for &ii in &subset[active..] {
            o_bar[ii] = ddot(&t.row(ii)[..n0], &w_bar[..n0]) + w_bar[n0];
        }
        if !ini {
            cgitermax = opts.cgitermax;
            ini = true;
        }
        let mut opt2 = true;
        for (k, &ii) in subset.iter().enumerate() {
            let yo = t.y(ii) * o_bar[ii];
            opt2 = if k < active {
                yo <= 1.0 + epsilon
            } else {
                yo >= 1.0 - epsilon
            };
            if !opt2 {
                break;
            }
        }
        if opt && opt2 {
            if epsilon == opts.big_epsilon {
                epsilon = opts.epsilon;
                continue;
            }
            return w_bar;
        }
        let delta = line_search(&t, &w, &w_bar, lambda, &o, &o_bar, &mut deltas);
        let f_old = f;
        dscal(1.0 - delta, &mut w);
        daxpy(delta, &w_bar, &mut w);
        for (oi, &ob) in o.iter_mut().zip(&o_bar) {
            *oi += delta * (ob - *oi);
        }
        let (a, fl) = partition(&t, &o, &mut subset, 0.5 * lambda * ddot(&w, &w));
        active = a;
        f = fl;
        if (f - f_old).abs() < opts.relative_stop_eps * f_old.abs() {
            return w;
        }
    }
    w
}
