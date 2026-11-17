use crate::budget::{Budget, Reservation};
use crate::normalize::Norm;
use crate::store::Store;
use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom, Write};
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;

pub struct Rows {
    data: RowData,
    pub nrows: usize,
    pub ncols: usize,
    stride: usize,
    _ram: Reservation,
}

enum RowData {
    Ram(Vec<f32>),
    Spill(PathBuf),
}

static ROWS_COUNTER: AtomicU64 = AtomicU64::new(0);

pub(crate) const LANES: usize = 8;

pub(crate) fn stride_for(ncols: usize) -> usize {
    ncols.div_ceil(LANES) * LANES
}

pub(crate) fn fill(store: &Store, norm: &Norm, idx: &[usize], mut put: impl FnMut(usize, &[f32])) {
    let nf = store.nfeat();
    let stride = stride_for(nf + 1);
    let mut buf = vec![0.0f32; stride];
    buf[nf] = 1.0;
    if idx.is_empty() {
        return;
    }
    if nf == 0 {
        for p in 0..idx.len() {
            put(p, &buf);
        }
        return;
    }
    let mut order: Vec<usize> = (0..idx.len()).collect();
    order.sort_by_key(|&p| idx[p]);
    let mut next = 0usize;
    store.for_each_block(|start, block| {
        let end = start + block.len() / nf;
        while next < order.len() && idx[order[next]] < end {
            let p = order[next];
            let r = idx[p] - start;
            let src = &block[r * nf..(r + 1) * nf];
            for (((d, &x), &m), &sd) in buf[..nf]
                .iter_mut()
                .zip(src)
                .zip(&norm.mean)
                .zip(&norm.inv_sd)
            {
                *d = ((f64::from(x) - m) * sd) as f32;
            }
            put(p, &buf);
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

impl Rows {
    pub fn gather(store: &Store, norm: &Norm, idx: &[usize]) -> Rows {
        let ncols = store.nfeat() + 1;
        let stride = stride_for(ncols);
        let mut data = vec![0.0f32; idx.len() * stride];
        fill(store, norm, idx, |p, row| {
            data[p * stride..(p + 1) * stride].copy_from_slice(row)
        });
        Rows {
            data: RowData::Ram(data),
            nrows: idx.len(),
            ncols,
            stride,
            _ram: Reservation::none(),
        }
    }

    pub fn all_within(store: &Store, norm: &Norm, budget: &Arc<Budget>) -> Result<Rows, String> {
        let nf = store.nfeat();
        let ncols = nf + 1;
        let stride = stride_for(ncols);
        let n = store.len();
        let bytes = (n * stride * 4) as u64;
        if nf == 0 || n == 0 {
            let all: Vec<usize> = (0..n).collect();
            return Rows::gather_within(store, norm, &all, budget);
        }
        let Some(ram) = Reservation::take(budget, bytes) else {
            let all: Vec<usize> = (0..n).collect();
            return Rows::gather_within(store, norm, &all, budget);
        };
        let mut data = vec![0.0f32; n * stride];
        let br = store.block_rows();
        let parts: Vec<std::sync::Mutex<&mut [f32]>> = data
            .chunks_mut(br * stride)
            .map(std::sync::Mutex::new)
            .collect();
        store.par_for_each_block(|start, block| {
            let mut dst = parts[start / br].lock().expect("rows block lock");
            for (out, src) in dst.chunks_exact_mut(stride).zip(block.chunks_exact(nf)) {
                for (((d, &x), &m), &sd) in out[..nf]
                    .iter_mut()
                    .zip(src)
                    .zip(&norm.mean)
                    .zip(&norm.inv_sd)
                {
                    *d = ((f64::from(x) - m) * sd) as f32;
                }
                out[nf] = 1.0;
            }
        });
        drop(parts);
        Ok(Rows {
            data: RowData::Ram(data),
            nrows: n,
            ncols,
            stride,
            _ram: ram,
        })
    }

    pub fn gather_within(
        store: &Store,
        norm: &Norm,
        idx: &[usize],
        budget: &Arc<Budget>,
    ) -> Result<Rows, String> {
        let ncols = store.nfeat() + 1;
        let stride = stride_for(ncols);
        let bytes = (idx.len() * stride * 4) as u64;
        if let Some(ram) = Reservation::take(budget, bytes) {
            let mut rows = Rows::gather(store, norm, idx);
            rows._ram = ram;
            return Ok(rows);
        }
        let k = ROWS_COUNTER.fetch_add(1, Ordering::Relaxed);
        let path = budget.spill_dir().join(format!(
            "percolator-native-{}-{}.rows",
            std::process::id(),
            k
        ));
        let rows = Rows {
            data: RowData::Spill(path.clone()),
            nrows: idx.len(),
            ncols,
            stride,
            _ram: Reservation::none(),
        };
        let file = File::create(&path)
            .map_err(|e| format!("create spill file {}: {e}", path.display()))?;
        file.set_len(bytes)
            .map_err(|e| format!("size spill file {}: {e}", path.display()))?;
        let mut file = file;
        let mut bytes_row = vec![0u8; stride * 4];
        let mut failed = None;
        fill(store, norm, idx, |p, row| {
            if failed.is_some() {
                return;
            }
            for (c, v) in bytes_row.chunks_exact_mut(4).zip(row) {
                c.copy_from_slice(&v.to_le_bytes());
            }
            let res = file
                .seek(SeekFrom::Start((p * stride * 4) as u64))
                .and_then(|_| file.write_all(&bytes_row));
            if let Err(e) = res {
                failed = Some(e);
            }
        });
        if let Some(e) = failed {
            return Err(format!("write spill file {}: {e}", path.display()));
        }
        Ok(rows)
    }

    pub fn from_dense(data: Vec<f64>, nrows: usize, ncols: usize) -> Rows {
        assert_eq!(data.len(), nrows * ncols);
        let stride = stride_for(ncols);
        let mut out = vec![0.0f32; nrows * stride];
        for (dst, src) in out
            .chunks_exact_mut(stride)
            .zip(data.chunks_exact(ncols.max(1)))
        {
            for (d, &v) in dst.iter_mut().zip(src) {
                *d = v as f32;
            }
        }
        Rows {
            data: RowData::Ram(out),
            nrows,
            ncols,
            stride,
            _ram: Reservation::none(),
        }
    }

    pub fn is_spilled(&self) -> bool {
        matches!(self.data, RowData::Spill(_))
    }

    pub fn row(&self, i: usize) -> &[f32] {
        match &self.data {
            RowData::Ram(v) => &v[i * self.stride..(i + 1) * self.stride],
            RowData::Spill(_) => panic!("Rows::row: the rows are spilled to disk; use a scan"),
        }
    }

    pub(crate) fn scan(&self, ids: &[usize], mut f: impl FnMut(usize, &[f32])) {
        match &self.data {
            RowData::Ram(_) => {
                for (k, &i) in ids.iter().enumerate() {
                    f(k, self.row(i));
                }
            }
            RowData::Spill(path) => {
                let file = File::open(path).expect("open rows spill file");
                let mut r = BufReader::with_capacity(8 << 20, file);
                let width = self.stride * 4;
                let mut bytes = vec![0u8; width];
                let mut vals = vec![0.0f32; self.stride];
                let mut cur = 0usize;
                for (k, &i) in ids.iter().enumerate() {
                    assert!(i >= cur, "Rows::scan: ids must increase");
                    if i > cur {
                        r.seek_relative(((i - cur) * width) as i64)
                            .expect("seek rows spill file");
                    }
                    r.read_exact(&mut bytes).expect("read rows spill file");
                    for (v, c) in vals.iter_mut().zip(bytes.chunks_exact(4)) {
                        *v = f32::from_le_bytes([c[0], c[1], c[2], c[3]]);
                    }
                    f(k, &vals);
                    cur = i + 1;
                }
            }
        }
    }
}

impl Drop for Rows {
    fn drop(&mut self) {
        if let RowData::Spill(p) = &self.data {
            let _ = std::fs::remove_file(p);
        }
    }
}

pub(crate) fn dot_row(x: &[f32], p: &[f64]) -> f64 {
    let mut acc = [0.0f64; LANES];
    for (xc, pc) in x.chunks_exact(LANES).zip(p.chunks_exact(LANES)) {
        for l in 0..LANES {
            acc[l] += f64::from(xc[l]) * pc[l];
        }
    }
    ((acc[0] + acc[4]) + (acc[2] + acc[6])) + ((acc[1] + acc[5]) + (acc[3] + acc[7]))
}

pub(crate) fn axpy_row(a: f64, x: &[f32], r: &mut [f64]) {
    if a != 0.0 {
        for (rc, xc) in r.chunks_exact_mut(LANES).zip(x.chunks_exact(LANES)) {
            for l in 0..LANES {
                rc[l] += a * f64::from(xc[l]);
            }
        }
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
    rows: Vec<usize>,
    pos: Vec<bool>,
    cpos: f64,
    cneg: f64,
}

impl<'a> Train<'a> {
    fn new(x: &'a Rows, idx_pos: &[usize], idx_neg: &[usize], cpos: f64, cneg: f64) -> Self {
        let mut ex: Vec<(usize, bool)> = idx_neg
            .iter()
            .map(|&r| (r, false))
            .chain(idx_pos.iter().map(|&r| (r, true)))
            .collect();
        ex.sort_unstable();
        Train {
            x,
            rows: ex.iter().map(|e| e.0).collect(),
            pos: ex.iter().map(|e| e.1).collect(),
            cpos,
            cneg,
        }
    }

    fn len(&self) -> usize {
        self.rows.len()
    }

    fn scan(&self, subset: &[usize], f: impl FnMut(usize, &[f32])) {
        let ids: Vec<usize> = subset.iter().map(|&i| self.rows[i]).collect();
        self.x.scan(&ids, f);
    }

    fn y(&self, i: usize) -> f64 {
        if self.pos[i] {
            1.0
        } else {
            -1.0
        }
    }

    fn c(&self, i: usize) -> f64 {
        if self.pos[i] {
            self.cpos
        } else {
            self.cneg
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

struct CgWork {
    z: Vec<f64>,
    q: Vec<f64>,
    r: Vec<f64>,
    p: Vec<f64>,
    u: Vec<f64>,
    s: Vec<f64>,
    beta_old: Vec<f64>,
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
    let n = t.x.ncols;
    let z = &mut ws.z[..active];
    let q = &mut ws.q[..active];
    let r = &mut ws.r;
    let p = &mut ws.p;
    let u = &mut ws.u;
    let sv = &mut ws.s;
    let beta_old = &mut ws.beta_old;
    u.fill(0.0);
    t.scan(subset, |k, row| {
        let ii = subset[k];
        z[k] = t.c(ii) * (t.y(ii) - o[ii]);
        axpy_row(z[k], row, u);
    });
    r.copy_from_slice(u);
    daxpy(-lambda, beta, r);
    p.copy_from_slice(r);
    let mut omega1 = ddot(&r[..n], &r[..n]);
    let mut omega_p = omega1;
    let epsilon2 = epsilon * epsilon;
    let mut cgiter = 0usize;
    while cgiter < cgitermax {
        cgiter += 1;
        sv.fill(0.0);
        let mut omega_q = 0.0;
        t.scan(subset, |k, row| {
            let tk = dot_row(row, p);
            q[k] = tk;
            let ct = t.c(subset[k]) * tk;
            omega_q += ct * tk;
            axpy_row(ct, row, sv);
        });
        let gamma = omega1 / (lambda * omega_p + omega_q);
        let inv_omega2 = 1.0 / omega1;
        beta_old.copy_from_slice(beta);
        daxpy(gamma, p, beta);
        for (k, &ii) in subset.iter().enumerate() {
            let g = gamma * q[k];
            o[ii] += g;
            z[k] -= t.c(ii) * g;
        }
        daxpy(-gamma, sv, u);
        r.copy_from_slice(u);
        daxpy(-lambda, beta_old, r);
        let omega_z = ddot(z, z);
        omega1 = ddot(&r[..n], &r[..n]);
        if omega1 < epsilon2 * omega_z {
            return true;
        }
        let scale = omega1 * inv_omega2;
        dscal(scale, p);
        daxpy(1.0, r, p);
        omega_p = ddot(&p[..n], &p[..n]);
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
    let t = Train::new(x, idx_pos, idx_neg, cpos, cneg);
    let n = x.ncols;
    let stride = x.stride;
    let m = t.len();
    let lambda = opts.lambda;
    let mut epsilon = opts.big_epsilon;
    let mut cgitermax = opts.small_cgitermax;
    let mut w = vec![0.0f64; stride];
    let mut o = vec![0.0f64; m];
    let mut subset = vec![0usize; m];
    let (mut active, mut f) = partition(&t, &o, &mut subset, 0.5 * lambda * ddot(&w[..n], &w[..n]));
    let mut w_bar = vec![0.0f64; stride];
    let mut o_bar = vec![0.0f64; m];
    let mut ws = CgWork {
        z: vec![0.0; m],
        q: vec![0.0; m],
        r: vec![0.0; stride],
        p: vec![0.0; stride],
        u: vec![0.0; stride],
        s: vec![0.0; stride],
        beta_old: vec![0.0; stride],
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
        let mut rest = subset[active..].to_vec();
        rest.sort_unstable();
        t.scan(&rest, |k, row| {
            o_bar[rest[k]] = dot_row(row, &w_bar);
        });
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
            w_bar.truncate(n);
            return w_bar;
        }
        let delta = line_search(&t, &w, &w_bar, lambda, &o, &o_bar, &mut deltas);
        let f_old = f;
        dscal(1.0 - delta, &mut w);
        daxpy(delta, &w_bar, &mut w);
        for (oi, &ob) in o.iter_mut().zip(&o_bar) {
            *oi += delta * (ob - *oi);
        }
        let (a, fl) = partition(&t, &o, &mut subset, 0.5 * lambda * ddot(&w[..n], &w[..n]));
        active = a;
        f = fl;
        if (f - f_old).abs() < opts.relative_stop_eps * f_old.abs() {
            w.truncate(n);
            return w;
        }
    }
    w.truncate(n);
    w
}
