use crate::budget::{Budget, Reservation};
use crate::normalize::Norm;
use crate::store::Store;
use rayon::prelude::*;
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

const CHUNK: usize = 4096;

pub(crate) const GRAM_BYTES: usize = 64 << 20;

const GROUP: usize = 4;

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

    pub fn bytes_for(nrows: usize, nfeat: usize) -> u64 {
        (nrows * stride_for(nfeat + 1) * 4) as u64
    }

    pub fn all_in_ram(store: &Store, norm: &Norm, budget: &Arc<Budget>, keep: u64) -> Option<Rows> {
        let nf = store.nfeat();
        let ncols = nf + 1;
        let stride = stride_for(ncols);
        let n = store.len();
        if nf == 0 || n == 0 {
            return None;
        }
        let bytes = Rows::bytes_for(n, nf);
        let mut ram = Reservation::none();
        if !ram.grow(budget, bytes, keep) {
            return None;
        }
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
        Some(Rows {
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

pub(crate) struct Examples {
    pub(crate) rows: Vec<usize>,
    pub(crate) pos: Vec<bool>,
}

impl Examples {
    pub(crate) fn from_lists(idx_pos: &[usize], idx_neg: &[usize]) -> Examples {
        let mut ex: Vec<(usize, bool)> = idx_neg
            .iter()
            .map(|&r| (r, false))
            .chain(idx_pos.iter().map(|&r| (r, true)))
            .collect();
        ex.sort_unstable();
        Examples {
            rows: ex.iter().map(|e| e.0).collect(),
            pos: ex.iter().map(|e| e.1).collect(),
        }
    }
}

struct Train<'a> {
    x: &'a Rows,
    rows: &'a [usize],
    pos: &'a [bool],
    cpos: f64,
    cneg: f64,
}

fn pass_ram<'r, P, I, B, R>(row: R, out: &mut [f64], init: I, body: B) -> Vec<P>
where
    P: Send,
    I: Fn() -> P + Sync,
    B: Fn(&mut P, usize, &[f32], &mut f64) + Sync,
    R: Fn(usize) -> &'r [f32] + Sync,
{
    out.par_chunks_mut(CHUNK)
        .enumerate()
        .map(|(c, outs)| {
            let mut part = init();
            for (j, o) in outs.iter_mut().enumerate() {
                let k = c * CHUNK + j;
                body(&mut part, k, row(k), o);
            }
            part
        })
        .collect()
}

impl Train<'_> {
    fn len(&self) -> usize {
        self.rows.len()
    }

    fn pass<P, I, B>(&self, subset: &[usize], out: &mut [f64], init: I, body: B) -> Vec<P>
    where
        P: Send,
        I: Fn() -> P + Sync,
        B: Fn(&mut P, usize, &[f32], &mut f64) + Sync,
    {
        match &self.x.data {
            RowData::Ram(_) => pass_ram(
                |k| self.x.row(self.rows[subset[k]]),
                &mut out[..subset.len()],
                init,
                body,
            ),
            RowData::Spill(_) => {
                let ids: Vec<usize> = subset.iter().map(|&i| self.rows[i]).collect();
                let mut parts: Vec<P> = Vec::with_capacity(subset.len().div_ceil(CHUNK));
                self.x.scan(&ids, |k, row| {
                    if k % CHUNK == 0 {
                        parts.push(init());
                    }
                    let part = parts.last_mut().expect("chunk partial");
                    body(part, k, row, &mut out[k]);
                });
                parts
            }
        }
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
    r: Vec<f64>,
    p: Vec<f64>,
    u: Vec<f64>,
    s: Vec<f64>,
    v: Vec<f64>,
    hd: Vec<f64>,
    d: Vec<f64>,
    beta0: Vec<f64>,
    beta_old: Vec<f64>,
    g: Vec<f64>,
    h: Vec<f64>,
}

struct GramPart {
    g: [Vec<f64>; 2],
    buf: [Vec<f64>; 2],
    fill: [usize; 2],
    u: Vec<f64>,
    v: Vec<f64>,
    zz: f64,
}

impl GramPart {
    fn new(stride: usize) -> GramPart {
        GramPart {
            g: [vec![0.0; stride * stride], vec![0.0; stride * stride]],
            buf: [vec![0.0; GROUP * stride], vec![0.0; GROUP * stride]],
            fill: [0, 0],
            u: vec![0.0; stride],
            v: vec![0.0; stride],
            zz: 0.0,
        }
    }

    fn push(&mut self, class: usize, row: &[f32], stride: usize, ncols: usize) {
        let k = self.fill[class];
        for (d, &x) in self.buf[class][k * stride..(k + 1) * stride]
            .iter_mut()
            .zip(row)
        {
            *d = f64::from(x);
        }
        self.fill[class] = k + 1;
        if k + 1 == GROUP {
            self.flush(class, stride, ncols);
        }
    }

    fn flush(&mut self, class: usize, stride: usize, ncols: usize) {
        let rows = self.fill[class];
        syrk(&mut self.g[class], &self.buf[class], rows, stride, ncols);
        self.fill[class] = 0;
    }
}

fn syrk(g: &mut [f64], buf: &[f64], rows: usize, stride: usize, ncols: usize) {
    for j in 0..ncols {
        let k0 = j / LANES * LANES;
        let gj = &mut g[j * stride + k0..(j + 1) * stride];
        if rows == GROUP {
            let (a0, a1, a2, a3) = (
                buf[j],
                buf[stride + j],
                buf[2 * stride + j],
                buf[3 * stride + j],
            );
            let x0 = &buf[k0..stride];
            let x1 = &buf[stride + k0..2 * stride];
            let x2 = &buf[2 * stride + k0..3 * stride];
            let x3 = &buf[3 * stride + k0..4 * stride];
            for ((((gv, &b0), &b1), &b2), &b3) in gj.iter_mut().zip(x0).zip(x1).zip(x2).zip(x3) {
                *gv += (a0 * b0 + a1 * b1) + (a2 * b2 + a3 * b3);
            }
        } else {
            for r in 0..rows {
                let a = buf[r * stride + j];
                for (gv, &b) in gj.iter_mut().zip(&buf[r * stride + k0..(r + 1) * stride]) {
                    *gv += a * b;
                }
            }
        }
    }
}

fn symv(g: &[f64], x: &[f64], out: &mut [f64], ncols: usize) {
    let stride = x.len();
    for (j, oj) in out.iter_mut().enumerate() {
        *oj = if j < ncols {
            let mut acc = [0.0f64; LANES];
            for (gc, xc) in g[j * stride..(j + 1) * stride]
                .chunks_exact(LANES)
                .zip(x.chunks_exact(LANES))
            {
                for l in 0..LANES {
                    acc[l] += gc[l] * xc[l];
                }
            }
            ((acc[0] + acc[4]) + (acc[2] + acc[6])) + ((acc[1] + acc[5]) + (acc[3] + acc[7]))
        } else {
            0.0
        };
    }
}

fn gram(t: &Train, subset: &[usize], o: &[f64], ws: &mut CgWork) -> f64 {
    let stride = t.x.stride;
    let ncols = t.x.ncols;
    let part_bytes = (2 * stride * stride + 2 * GROUP * stride + 2 * stride) * 8;
    let batch = (GRAM_BYTES / part_bytes / CHUNK).max(1) * CHUNK;
    let mut gs = [vec![0.0f64; stride * stride], vec![0.0f64; stride * stride]];
    ws.u.fill(0.0);
    ws.v.fill(0.0);
    let mut zz = 0.0;
    let mut scratch = vec![0.0f64; batch.min(subset.len())];
    for sub in subset.chunks(batch) {
        let parts = t.pass(
            sub,
            &mut scratch,
            || GramPart::new(stride),
            |acc, k, row, _| {
                let ii = sub[k];
                let c = t.c(ii);
                let zk = c * (t.y(ii) - o[ii]);
                axpy_row(zk, row, &mut acc.u);
                axpy_row(c * zk, row, &mut acc.v);
                acc.zz += zk * zk;
                acc.push(usize::from(t.pos[ii]), row, stride, ncols);
            },
        );
        for mut part in parts {
            for (class, gc) in gs.iter_mut().enumerate() {
                if part.fill[class] > 0 {
                    part.flush(class, stride, ncols);
                }
                for (a, &b) in gc.iter_mut().zip(&part.g[class]) {
                    *a += b;
                }
            }
            for (a, &b) in ws.u.iter_mut().zip(&part.u) {
                *a += b;
            }
            for (a, &b) in ws.v.iter_mut().zip(&part.v) {
                *a += b;
            }
            zz += part.zz;
        }
    }
    for j in 0..ncols {
        for k in j + 1..ncols {
            for gc in gs.iter_mut() {
                gc[k * stride + j] = gc[j * stride + k];
            }
        }
    }
    let (cn, cp) = (t.cneg, t.cpos);
    for (((g, h), &n), &p) in ws.g.iter_mut().zip(ws.h.iter_mut()).zip(&gs[0]).zip(&gs[1]) {
        *g = cn * n + cp * p;
        *h = (cn * cn) * n + (cp * cp) * p;
    }
    zz
}

#[allow(clippy::too_many_arguments)]
fn cgls(
    t: &Train,
    lambda: f64,
    cgitermax: usize,
    epsilon: f64,
    subset: &[usize],
    beta: &mut [f64],
    o: &[f64],
    ws: &mut CgWork,
) -> bool {
    let n = t.x.ncols;
    let zz0 = gram(t, subset, o, ws);
    let CgWork {
        r,
        p,
        u,
        s: sv,
        v,
        hd,
        d,
        beta0,
        beta_old,
        g,
        h,
    } = ws;
    beta0.copy_from_slice(beta);
    r.copy_from_slice(u);
    daxpy(-lambda, beta, r);
    p.copy_from_slice(r);
    let mut omega1 = ddot(&r[..n], &r[..n]);
    let mut omega_p = omega1;
    let epsilon2 = epsilon * epsilon;
    let mut cgiter = 0usize;
    while cgiter < cgitermax {
        cgiter += 1;
        symv(g, p, sv, n);
        let omega_q = ddot(&p[..n], &sv[..n]);
        let gamma = omega1 / (lambda * omega_p + omega_q);
        let inv_omega2 = 1.0 / omega1;
        beta_old.copy_from_slice(beta);
        daxpy(gamma, p, beta);
        daxpy(-gamma, sv, u);
        r.copy_from_slice(u);
        daxpy(-lambda, beta_old, r);
        for ((dk, &b), &b0) in d.iter_mut().zip(beta.iter()).zip(beta0.iter()) {
            *dk = b - b0;
        }
        symv(h, d, hd, n);
        let omega_z = zz0 - 2.0 * ddot(&d[..n], &v[..n]) + ddot(&d[..n], &hd[..n]);
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
    mfn_from(
        x,
        &Examples::from_lists(idx_pos, idx_neg),
        cpos,
        cneg,
        opts,
        None,
    )
}

pub(crate) fn mfn_from(
    x: &Rows,
    ex: &Examples,
    cpos: f64,
    cneg: f64,
    opts: &MfnOpts,
    start: Option<&[f64]>,
) -> Vec<f64> {
    let t = Train {
        x,
        rows: &ex.rows,
        pos: &ex.pos,
        cpos,
        cneg,
    };
    let n = x.ncols;
    let stride = x.stride;
    let m = t.len();
    let lambda = opts.lambda;
    let mut epsilon = opts.big_epsilon;
    let mut cgitermax = opts.small_cgitermax;
    let mut w = vec![0.0f64; stride];
    let mut o = vec![0.0f64; m];
    let all: Vec<usize> = (0..m).collect();
    if let Some(w0) = start {
        w[..n].copy_from_slice(&w0[..n]);
        let wv: &[f64] = &w;
        t.pass(&all, &mut o, || (), |_, _, row, ok| *ok = dot_row(row, wv));
    }
    let mut subset = vec![0usize; m];
    let (mut active, mut f) = partition(&t, &o, &mut subset, 0.5 * lambda * ddot(&w[..n], &w[..n]));
    let mut w_bar = vec![0.0f64; stride];
    let mut o_bar = vec![0.0f64; m];
    let mut ws = CgWork {
        r: vec![0.0; stride],
        p: vec![0.0; stride],
        u: vec![0.0; stride],
        s: vec![0.0; stride],
        v: vec![0.0; stride],
        hd: vec![0.0; stride],
        d: vec![0.0; stride],
        beta0: vec![0.0; stride],
        beta_old: vec![0.0; stride],
        g: vec![0.0; stride * stride],
        h: vec![0.0; stride * stride],
    };
    let mut act = vec![false; m];
    let mut step = vec![0.0f64; stride];
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
            &o,
            &mut ws,
        );
        act.fill(false);
        for &ii in &subset[..active] {
            act[ii] = true;
        }
        for ((sk, &b), &a) in step.iter_mut().zip(&w_bar).zip(&w) {
            *sk = b - a;
        }
        let wb: &[f64] = &w_bar;
        let dv: &[f64] = &step;
        let av: &[bool] = &act;
        t.pass(
            &all,
            &mut o_bar,
            || (),
            |_, k, row, ob| {
                if av[k] {
                    *ob += dot_row(row, dv);
                } else {
                    *ob = dot_row(row, wb);
                }
            },
        );
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
