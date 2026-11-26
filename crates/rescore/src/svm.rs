use crate::budget::{Budget, Reservation};
use crate::normalize::Norm;
use crate::store::Store;
use rayon::prelude::*;
use std::borrow::Cow;
use std::fs::File;
use std::io::{BufWriter, Seek, SeekFrom, Write};
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

const REACH_MIN: f64 = 1.0;

const REACH_START: f64 = 2.0;

const REACH_MAX: f64 = 2.0;

const FRESH: u32 = u32::MAX;

/// Breakpoints from which the line search sorts them on the thread pool.
const PAR_SORT_MIN: usize = 1 << 15;

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
        let path = budget
            .spill_dir()
            .join(format!("rescore-{}-{}.rows", std::process::id(), k));
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
        let mut file = BufWriter::with_capacity(8 << 20, file);
        let width = (stride * 4) as u64;
        let mut at = 0u64;
        let mut bytes_row = vec![0u8; stride * 4];
        let mut failed = None;
        fill(store, norm, idx, |p, row| {
            if failed.is_some() {
                return;
            }
            for (c, v) in bytes_row.chunks_exact_mut(4).zip(row) {
                c.copy_from_slice(&v.to_le_bytes());
            }
            let to = p as u64 * width;
            let res = if to == at {
                Ok(0)
            } else {
                file.seek(SeekFrom::Start(to))
            }
            .and_then(|_| file.write_all(&bytes_row));
            at = to + width;
            if let Err(e) = res {
                failed = Some(e);
            }
        });
        if let Some(e) = failed.or_else(|| file.flush().err()) {
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

    #[cfg(test)]
    pub(crate) fn scan(&self, ids: &[usize], mut f: impl FnMut(usize, &[f32])) {
        match &self.data {
            RowData::Ram(_) => {
                for (k, &i) in ids.iter().enumerate() {
                    f(k, self.row(i));
                }
            }
            RowData::Spill(path) => {
                let file = File::open(path).expect("open rows spill file");
                read_rows(&file, self.stride, ids, f);
            }
        }
    }
}

/// Rows at most this many bytes apart are read with one `pread`, the bytes between them
/// included: reading them costs less than another call.
const READ_GAP_BYTES: usize = 64 << 10;

/// Bytes one `pread` of spilled rows reads at most.
const READ_SPAN_BYTES: usize = 4 << 20;

/// Calls `f(j, row)` for each row `ids[j]` of a spilled row file, in order. `ids` must not
/// decrease; nearby rows are read together with positioned reads, so several threads can read
/// one file at once.
fn read_rows(file: &File, stride: usize, ids: &[usize], mut f: impl FnMut(usize, &[f32])) {
    use std::os::unix::fs::FileExt;
    let width = stride * 4;
    let span_rows = (READ_SPAN_BYTES / width).max(1);
    let gap_rows = (READ_GAP_BYTES / width).max(1);
    let mut bytes: Vec<u8> = Vec::new();
    let mut vals = vec![0.0f32; stride];
    let mut a = 0usize;
    while a < ids.len() {
        let lo = ids[a];
        let mut b = a + 1;
        while b < ids.len() {
            assert!(ids[b] >= ids[b - 1], "Rows::pass: ids must not decrease");
            if ids[b] - ids[b - 1] > gap_rows || ids[b] - lo >= span_rows {
                break;
            }
            b += 1;
        }
        let hi = ids[b - 1] + 1;
        bytes.resize((hi - lo) * width, 0);
        file.read_exact_at(&mut bytes, (lo * width) as u64)
            .expect("read rows spill file");
        for (j, &i) in ids.iter().enumerate().take(b).skip(a) {
            let src = &bytes[(i - lo) * width..(i - lo + 1) * width];
            for (v, c) in vals.iter_mut().zip(src.chunks_exact(4)) {
                *v = f32::from_le_bytes([c[0], c[1], c[2], c[3]]);
            }
            f(j, &vals);
        }
        a = b;
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
    pub mfnitermax: usize,
    pub big_epsilon: f64,
    pub relative_stop_eps: f64,
}

impl Default for MfnOpts {
    fn default() -> Self {
        MfnOpts {
            lambda: 1.0,
            epsilon: 1e-7,
            mfnitermax: 50,
            big_epsilon: 0.01,
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
            RowData::Spill(path) => {
                let ids: Vec<usize> = subset.iter().map(|&i| self.rows[i]).collect();
                let file = File::open(path).expect("open rows spill file");
                let stride = self.x.stride;
                out[..subset.len()]
                    .par_chunks_mut(CHUNK)
                    .enumerate()
                    .map(|(c, outs)| {
                        let mut part = init();
                        let first = c * CHUNK;
                        read_rows(&file, stride, &ids[first..first + outs.len()], |j, row| {
                            body(&mut part, first + j, row, &mut outs[j]);
                        });
                        part
                    })
                    .collect()
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

struct Solve {
    b: Vec<f64>,
    g: Vec<f64>,
}

struct GramPart {
    g: [Vec<f64>; 2],
    buf: [Vec<f64>; 2],
    fill: [usize; 2],
    b: Vec<f64>,
}

impl GramPart {
    fn new(stride: usize) -> GramPart {
        GramPart {
            g: [vec![0.0; stride * stride], vec![0.0; stride * stride]],
            buf: [vec![0.0; GROUP * stride], vec![0.0; GROUP * stride]],
            fill: [0, 0],
            b: vec![0.0; stride],
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

fn gram_batch(stride: usize) -> usize {
    let part_bytes = (2 * stride * stride + 2 * GROUP * stride + stride) * 8;
    (GRAM_BYTES / part_bytes).max(1) * CHUNK
}

fn gram(t: &Train, subset: &[usize], ws: &mut Solve) {
    let stride = t.x.stride;
    let ncols = t.x.ncols;
    let batch = gram_batch(stride);
    let mut gs = [vec![0.0f64; stride * stride], vec![0.0f64; stride * stride]];
    ws.b.fill(0.0);
    let mut scratch = vec![0.0f64; batch.min(subset.len())];
    for sub in subset.chunks(batch) {
        let parts = t.pass(
            sub,
            &mut scratch,
            || GramPart::new(stride),
            |acc, k, row, _| {
                let ii = sub[k];
                axpy_row(t.c(ii) * t.y(ii), row, &mut acc.b);
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
            for (a, &b) in ws.b.iter_mut().zip(&part.b) {
                *a += b;
            }
        }
    }
    let (cn, cp) = (t.cneg, t.cpos);
    for j in 0..ncols {
        for k in j..ncols {
            let v = cn * gs[0][j * stride + k] + cp * gs[1][j * stride + k];
            ws.g[j * stride + k] = v;
            ws.g[k * stride + j] = v;
        }
    }
}

fn direct(t: &Train, lambda: f64, subset: &[usize], beta: &mut [f64], ws: &mut Solve) -> bool {
    let n = t.x.ncols;
    let stride = t.x.stride;
    gram(t, subset, ws);
    let g = &mut ws.g;
    for j in 0..n {
        g[j * stride + j] += lambda;
    }
    for j in 0..n {
        let mut s = g[j * stride + j];
        for k in 0..j {
            s -= g[j * stride + k] * g[j * stride + k];
        }
        if s.partial_cmp(&0.0) != Some(std::cmp::Ordering::Greater) {
            return false;
        }
        let ljj = s.sqrt();
        g[j * stride + j] = ljj;
        for i in j + 1..n {
            let mut v = g[i * stride + j];
            for k in 0..j {
                v -= g[i * stride + k] * g[j * stride + k];
            }
            g[i * stride + j] = v / ljj;
        }
    }
    for i in 0..n {
        let mut v = ws.b[i];
        for k in 0..i {
            v -= g[i * stride + k] * beta[k];
        }
        beta[i] = v / g[i * stride + i];
    }
    for i in (0..n).rev() {
        let mut v = beta[i];
        for k in i + 1..n {
            v -= g[k * stride + i] * beta[k];
        }
        beta[i] = v / g[i * stride + i];
    }
    beta[n..].fill(0.0);
    true
}

fn norm_up(row: &[f32]) -> f32 {
    let s: f64 = row.iter().map(|&v| f64::from(v) * f64::from(v)).sum();
    (s.sqrt() * (1.0 + 1e-6)) as f32
}

fn l2(a: &[f64]) -> f64 {
    a.iter().map(|v| v * v).sum::<f64>().sqrt()
}

fn dist(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y) * (x - y))
        .sum::<f64>()
        .sqrt()
}

struct Screen<'a> {
    nrm: Cow<'a, [f32]>,
    since: Vec<u32>,
    live: Vec<u32>,
    w_at: Vec<Vec<f64>>,
    w_bar_at: Vec<Vec<f64>>,
    delta_at: Vec<f64>,
    reach: Vec<f64>,
    slack: f64,
    lim: f64,
    last: f64,
    redo: bool,
}

impl<'a> Screen<'a> {
    fn new(nrm: Cow<'a, [f32]>, w: &[f64]) -> Screen<'a> {
        Screen {
            since: vec![FRESH; nrm.len()],
            live: (0..nrm.len() as u32).collect(),
            nrm,
            w_at: vec![w.to_vec()],
            w_bar_at: Vec::new(),
            delta_at: Vec::new(),
            reach: Vec::new(),
            slack: 0.0,
            lim: REACH_START,
            last: REACH_START,
            redo: false,
        }
    }

    fn bound(&mut self, w_bar: &[f64]) {
        if !self.lim.is_finite() {
            self.reach.fill(f64::INFINITY);
            return;
        }
        let lim = self.lim;
        let w = self.w_at.last().expect("current w");
        let e: Vec<f64> = w
            .iter()
            .zip(w_bar)
            .map(|(&a, &b)| a + lim * (b - a))
            .collect();
        let big = self
            .w_at
            .iter()
            .chain(&self.w_bar_at)
            .map(|v| l2(v.as_slice()))
            .fold(l2(w_bar).max(l2(&e)), f64::max);
        self.reach = self
            .w_at
            .iter()
            .map(|wr| dist(w, wr).max(dist(&e, wr)))
            .collect();
        self.slack = 1e-10 * big;
    }

    fn any_stale(&self) -> bool {
        self.live.len() < self.since.len()
    }

    fn record(&mut self, w_bar: &[f64], delta: f64, w: &[f64]) {
        self.w_bar_at.push(w_bar.to_vec());
        self.delta_at.push(delta);
        self.w_at.push(w.to_vec());
    }
}

fn screened_pass(
    t: &Train,
    sc: &mut Screen<'_>,
    act: &[bool],
    o: &mut [f64],
    o_bar: &mut [f64],
    w_bar: &[f64],
) {
    let u = sc.w_at.len() - 1;
    let nrm: &[f32] = &sc.nrm;
    let w_bar_at: &[Vec<f64>] = &sc.w_bar_at;
    let delta_at: &[f64] = &sc.delta_at;
    let reach: &[f64] = &sc.reach;
    let slack = sc.slack;
    let redo = sc.redo;
    let since = &mut sc.since;
    let visit = |k: usize, ob: &mut f64, oi: &mut f64, s: &mut u32| {
        if redo && *s == FRESH {
            return;
        }
        let r = if *s == FRESH { u } else { *s as usize };
        if !act[k] {
            let lim = 1.0 + 1e-6 + 1e-10 * oi.abs() + f64::from(nrm[k]) * (reach[r] + slack);
            if t.y(k) * *oi > lim {
                if *s == FRESH {
                    *s = u as u32;
                }
                return;
            }
        }
        let row = t.x.row(t.rows[k]);
        if *s != FRESH {
            for (wb, &d) in w_bar_at[r..u].iter().zip(&delta_at[r..u]) {
                let v = dot_row(row, wb);
                *oi += d * (v - *oi);
            }
            *s = FRESH;
        }
        *ob = dot_row(row, w_bar);
    };
    let parts: Vec<Vec<u32>> = o_bar
        .par_chunks_mut(CHUNK)
        .zip(o.par_chunks_mut(CHUNK))
        .zip(since.par_chunks_mut(CHUNK))
        .enumerate()
        .map(|(c, ((obs, os), ss))| {
            let mut fresh = Vec::new();
            for (j, ((ob, oi), s)) in obs
                .iter_mut()
                .zip(os.iter_mut())
                .zip(ss.iter_mut())
                .enumerate()
            {
                let k = c * CHUNK + j;
                visit(k, ob, oi, s);
                if *s == FRESH {
                    fresh.push(k as u32);
                }
            }
            fresh
        })
        .collect();
    sc.live = parts.concat();
}

fn each_row(m: usize, sc: Option<&Screen<'_>>, mut f: impl FnMut(usize)) {
    match sc {
        Some(sc) => sc.live.iter().for_each(|&i| f(i as usize)),
        None => (0..m).for_each(f),
    }
}

/// A line-search breakpoint: the step at which a row's loss changes, the row, and the sign of
/// the change.
type Breakpoint = (f64, usize, f64);

/// One chunk of rows' line-search slope terms and breakpoints, in row order.
type RowTerms = (Vec<(f64, f64)>, Vec<Breakpoint>);

#[allow(clippy::too_many_arguments)]
fn line_search(
    t: &Train,
    w: &[f64],
    w_bar: &[f64],
    lambda: f64,
    o: &[f64],
    o_bar: &[f64],
    deltas: &mut Vec<(f64, usize, f64)>,
    sc: Option<&Screen<'_>>,
) -> (f64, f64) {
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
    let visit = |i: usize, terms: &mut Vec<(f64, f64)>, found: &mut Vec<Breakpoint>| {
        let y = t.y(i);
        let diff = y * (o_bar[i] - o[i]);
        if y * o[i] < 1.0 {
            let d2 = t.c(i) * (o_bar[i] - o[i]);
            terms.push(((o[i] - y) * d2, (o_bar[i] - y) * d2));
            if diff > 0.0 {
                found.push(((1.0 - y * o[i]) / diff, i, -1.0));
            }
        } else if diff < 0.0 {
            found.push(((1.0 - y * o[i]) / diff, i, 1.0));
        }
    };
    let n = sc.map_or(t.len(), |sc| sc.live.len());
    let row = |k: usize| sc.map_or(k, |sc| sc.live[k] as usize);
    let parts: Vec<RowTerms> = (0..n.div_ceil(CHUNK))
        .into_par_iter()
        .map(|c| {
            let (mut terms, mut found) = (Vec::new(), Vec::new());
            for k in c * CHUNK..((c + 1) * CHUNK).min(n) {
                visit(row(k), &mut terms, &mut found);
            }
            (terms, found)
        })
        .collect();
    for (terms, found) in parts {
        for (dl, dr) in terms {
            l += dl;
            r += dr;
        }
        deltas.extend(found);
    }
    if deltas.len() >= PAR_SORT_MIN {
        deltas.par_sort_by(|a, b| a.0.total_cmp(&b.0));
    } else {
        deltas.sort_by(|a, b| a.0.total_cmp(&b.0));
    }
    let mut stop = None;
    for &(delta, ii, s) in deltas.iter() {
        let delta_prime = l + delta * (r - l);
        if delta_prime >= 0.0 {
            stop = Some(delta);
            break;
        }
        let diff = s * t.c(ii) * (o_bar[ii] - o[ii]);
        l += diff * (o[ii] - t.y(ii));
        r += diff * (o_bar[ii] - t.y(ii));
    }
    let step = -l / (r - l);
    let reach = match stop {
        Some(d) => d,
        None if r - l > 0.0 => deltas.last().map_or(step, |d| d.0.max(step)),
        None => f64::INFINITY,
    };
    (step, reach)
}

fn partition(
    t: &Train,
    o: &[f64],
    subset: &mut [usize],
    f_start: f64,
    sc: Option<&Screen<'_>>,
) -> (usize, f64) {
    let m = t.len();
    let mut active = 0usize;
    let mut inactive = m;
    let mut f = f_start;
    each_row(m, sc, |i| {
        let diff = 1.0 - t.y(i) * o[i];
        if diff > 0.0 {
            subset[active] = i;
            active += 1;
            f += 0.5 * t.c(i) * diff * diff;
        } else {
            inactive -= 1;
            subset[inactive] = i;
        }
    });
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
        None,
    )
}

fn padded(w0: &[f64], stride: usize) -> Vec<f64> {
    let mut w = vec![0.0f64; stride];
    w[..w0.len()].copy_from_slice(w0);
    w
}

pub(crate) fn warm_scores(x: &Rows, ex: &Examples, starts: &[&[f64]]) -> (Vec<Vec<f64>>, Vec<f32>) {
    let m = ex.rows.len();
    let k = starts.len();
    if k == 0 {
        return (Vec::new(), Vec::new());
    }
    let n = x.ncols;
    let ws: Vec<Vec<f64>> = starts.iter().map(|w0| padded(&w0[..n], x.stride)).collect();
    let screen = !x.is_spilled();
    let t = Train {
        x,
        rows: &ex.rows,
        pos: &ex.pos,
        cpos: 0.0,
        cneg: 0.0,
    };
    let all: Vec<usize> = (0..m).collect();
    let mut first = vec![0.0f64; m];
    let parts = t.pass(
        &all,
        &mut first,
        || (Vec::new(), Vec::new()),
        |p: &mut (Vec<f64>, Vec<f32>), _, row, o0| {
            *o0 = dot_row(row, &ws[0]);
            for w in &ws[1..] {
                p.0.push(dot_row(row, w));
            }
            if screen {
                p.1.push(norm_up(row));
            }
        },
    );
    let mut os = vec![first];
    for j in 1..k {
        let mut o = Vec::with_capacity(m);
        for p in &parts {
            o.extend(p.0.iter().skip(j - 1).step_by(k - 1));
        }
        os.push(o);
    }
    let nrm = if screen {
        parts.iter().flat_map(|p| p.1.iter().copied()).collect()
    } else {
        Vec::new()
    };
    (os, nrm)
}

pub(crate) fn mfn_from(
    x: &Rows,
    ex: &Examples,
    cpos: f64,
    cneg: f64,
    opts: &MfnOpts,
    start: Option<&[f64]>,
    pre: Option<(Vec<f64>, &[f32])>,
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
    let mut w = vec![0.0f64; stride];
    let mut o = vec![0.0f64; m];
    let all: Vec<usize> = (0..m).collect();
    let mut o_bar = vec![0.0f64; m];
    let screen = !x.is_spilled() && m > 0;
    let mut nrm = Cow::Owned(Vec::new());
    if let (Some(w0), Some((o0, n0))) = (start, pre) {
        w[..n].copy_from_slice(&w0[..n]);
        o = o0;
        nrm = Cow::Borrowed(n0);
    } else if let Some(w0) = start {
        w[..n].copy_from_slice(&w0[..n]);
        let wv: &[f64] = &w;
        let parts = t.pass(
            &all,
            &mut o,
            || Vec::with_capacity(if screen { CHUNK } else { 0 }),
            |p, _, row, ok| {
                *ok = dot_row(row, wv);
                if screen {
                    p.push(norm_up(row));
                }
            },
        );
        nrm = Cow::Owned(parts.concat());
    } else if screen {
        nrm = Cow::Owned(
            t.pass(
                &all,
                &mut o_bar,
                || Vec::with_capacity(CHUNK),
                |p, _, row, _| p.push(norm_up(row)),
            )
            .concat(),
        );
    }
    let mut sc = screen.then(|| Screen::new(nrm, &w));
    let mut subset = vec![0usize; m];
    let (mut active, mut f) = partition(
        &t,
        &o,
        &mut subset,
        0.5 * lambda * ddot(&w[..n], &w[..n]),
        sc.as_ref(),
    );
    let mut w_bar = vec![0.0f64; stride];
    let mut ws = Solve {
        b: vec![0.0; stride],
        g: vec![0.0; stride * stride],
    };
    let mut act = vec![false; m];
    let mut deltas: Vec<(f64, usize, f64)> = Vec::with_capacity(m);
    let mut iter = 0usize;
    let mut opt = false;
    let mut recheck = false;
    while iter < opts.mfnitermax {
        iter += 1;
        if !recheck {
            w_bar.copy_from_slice(&w);
            opt = direct(&t, lambda, &subset[..active], &mut w_bar, &mut ws);
            if !opt {
                w_bar.copy_from_slice(&w);
            }
            act.fill(false);
            for &ii in &subset[..active] {
                act[ii] = true;
            }
            if let Some(sc) = sc.as_mut() {
                sc.lim = (1.25 * sc.last).clamp(REACH_MIN, REACH_MAX);
                sc.redo = false;
                sc.bound(&w_bar);
                screened_pass(&t, sc, &act, &mut o, &mut o_bar, &w_bar);
            } else {
                let wb: &[f64] = &w_bar;
                t.pass(
                    &all,
                    &mut o_bar,
                    || (),
                    |_, _, row, ob| {
                        *ob = dot_row(row, wb);
                    },
                );
            }
        }
        recheck = false;
        let mut opt2 = true;
        match sc.as_ref() {
            Some(sc) => {
                opt2 = subset[..active]
                    .iter()
                    .all(|&ii| t.y(ii) * o_bar[ii] <= 1.0 + epsilon)
                    && sc.live.iter().all(|&i| {
                        let i = i as usize;
                        act[i] || t.y(i) * o_bar[i] >= 1.0 - epsilon
                    });
            }
            None => {
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
            }
        }
        if opt && opt2 {
            if epsilon == opts.big_epsilon {
                epsilon = opts.epsilon;
                recheck = true;
                continue;
            }
            w_bar.truncate(n);
            return w_bar;
        }
        let (mut delta, mut need) =
            line_search(&t, &w, &w_bar, lambda, &o, &o_bar, &mut deltas, sc.as_ref());
        if let Some(sc) = sc.as_mut() {
            loop {
                let within = need <= sc.lim * (1.0 - 1e-6);
                if within || !sc.any_stale() {
                    break;
                }
                sc.lim = if need.is_finite() {
                    (1.25 * need).max(2.0 * sc.lim)
                } else {
                    f64::INFINITY
                };
                sc.redo = true;
                sc.bound(&w_bar);
                screened_pass(&t, sc, &act, &mut o, &mut o_bar, &w_bar);
                (delta, need) =
                    line_search(&t, &w, &w_bar, lambda, &o, &o_bar, &mut deltas, Some(&*sc));
            }
            sc.last = delta;
        }
        let f_old = f;
        dscal(1.0 - delta, &mut w);
        daxpy(delta, &w_bar, &mut w);
        match sc.as_mut() {
            Some(sc) => {
                for &i in &sc.live {
                    let i = i as usize;
                    o[i] += delta * (o_bar[i] - o[i]);
                }
                sc.record(&w_bar, delta, &w);
            }
            None => {
                for (oi, &ob) in o.iter_mut().zip(&o_bar) {
                    *oi += delta * (ob - *oi);
                }
            }
        }
        let (a, fl) = partition(
            &t,
            &o,
            &mut subset,
            0.5 * lambda * ddot(&w[..n], &w[..n]),
            sc.as_ref(),
        );
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn spilled_scan_repeats_a_row_as_the_ram_scan_does() {
        let dir = tempfile::tempdir().unwrap();
        let unlimited = Budget::new(u64::MAX, dir.path());
        let mut b = crate::store::StoreBuilder::new(3, &unlimited);
        for r in 0..6 {
            b.push_row(&[r as f32, 2.0 * r as f32, 0.5]).unwrap();
        }
        let store = b.finish().unwrap();
        let norm = Norm {
            mean: vec![0.0; 3],
            inv_sd: vec![1.0; 3],
        };
        let all: Vec<usize> = (0..6).collect();
        let ram = Rows::gather(&store, &norm, &all);
        let spilled =
            Rows::gather_within(&store, &norm, &all, &Budget::new(0, dir.path())).unwrap();
        assert!(spilled.is_spilled());
        let ids = [0, 0, 1, 3, 3, 3, 5];
        let collect = |rows: &Rows| {
            let mut v = Vec::new();
            rows.scan(&ids, |k, row| v.push((k, row.to_vec())));
            v
        };
        assert_eq!(collect(&spilled), collect(&ram));
        assert_eq!(collect(&ram).len(), ids.len());
    }

    #[test]
    fn a_spilled_pass_gives_the_ram_pass_partials_and_outputs() {
        let dir = tempfile::tempdir().unwrap();
        let unlimited = Budget::new(u64::MAX, dir.path());
        let nf = 5;
        let mut b = crate::store::StoreBuilder::new(nf, &unlimited);
        let n = 3 * CHUNK + 1234;
        for r in 0..n {
            let row: Vec<f32> = (0..nf)
                .map(|c| ((r * 7919 + c * 104_729) % 1000) as f32 * 0.37 - 150.0)
                .collect();
            b.push_row(&row).unwrap();
        }
        let store = b.finish().unwrap();
        let norm = Norm {
            mean: vec![1.5; nf],
            inv_sd: vec![0.25; nf],
        };
        let all: Vec<usize> = (0..n).collect();
        let ram = Rows::gather(&store, &norm, &all);
        let spilled =
            Rows::gather_within(&store, &norm, &all, &Budget::new(0, dir.path())).unwrap();
        assert!(spilled.is_spilled());
        let rows: Vec<usize> = (0..n).filter(|r| r % 3 != 1 || r % 41 == 0).collect();
        let pos: Vec<bool> = rows.iter().map(|r| r % 5 == 0).collect();
        let subset: Vec<usize> = (0..rows.len())
            .filter(|k| ((k / 2500) % 2 == 0 && k % 97 < 60) || k % 1777 == 0)
            .collect();
        let w: Vec<f64> = (0..ram.stride).map(|j| 0.1 * j as f64 - 0.3).collect();
        let run = |x: &Rows| {
            let t = Train {
                x,
                rows: &rows,
                pos: &pos,
                cpos: 1.0,
                cneg: 1.0,
            };
            let mut out = vec![0.0f64; subset.len()];
            let parts = t.pass(
                &subset,
                &mut out,
                || (0.0f64, Vec::new()),
                |p, k, row, o| {
                    *o = dot_row(row, &w);
                    p.0 += *o;
                    p.1.push(k);
                },
            );
            let parts: Vec<(u64, Vec<usize>)> =
                parts.into_iter().map(|(a, ks)| (a.to_bits(), ks)).collect();
            let out: Vec<u64> = out.iter().map(|v| v.to_bits()).collect();
            (parts, out)
        };
        let (a, b) = (run(&ram), run(&spilled));
        assert_eq!(a.0.len(), subset.len().div_ceil(CHUNK));
        assert!(a == b, "spilled pass differs from the RAM pass");
    }

    #[test]
    fn gram_batch_holds_many_chunk_partials_within_the_workspace() {
        for stride in [8, 24, 64, 104, 1024, 4096] {
            let batch = gram_batch(stride);
            let parts = batch / CHUNK;
            assert_eq!(batch % CHUNK, 0);
            let part_bytes = (2 * stride * stride + 2 * GROUP * stride + stride) * 8;
            assert!(
                parts == 1 || parts * part_bytes <= GRAM_BYTES,
                "stride {stride}"
            );
            if part_bytes * 2 <= GRAM_BYTES {
                assert!(parts >= 2, "stride {stride}: {parts} partials");
            }
        }
        assert_eq!(gram_batch(64) / CHUNK, 956);
    }
}
