use crate::pin::Pin;
use crate::store::Store;
use std::sync::Mutex;

pub struct Norm {
    pub mean: Vec<f64>,
    pub inv_sd: Vec<f64>,
}

pub fn stdv(pin: &Pin) -> Norm {
    let nf = pin.rows.nfeat();
    let n = pin.rows.len() as f64;
    let mut mean = vec![0.0f64; nf];
    let mut lo = vec![f32::INFINITY; nf];
    let mut hi = vec![f32::NEG_INFINITY; nf];
    pin.rows.for_each_block(|_, block| {
        for row in block.chunks_exact(nf) {
            for (((m, l), h), &x) in mean.iter_mut().zip(&mut lo).zip(&mut hi).zip(row) {
                *m += f64::from(x);
                *l = l.min(x);
                *h = h.max(x);
            }
        }
    });
    if n > 0.0 {
        for m in mean.iter_mut() {
            *m /= n;
        }
    }
    for ((m, &l), &h) in mean.iter_mut().zip(&lo).zip(&hi) {
        if l == h {
            *m = f64::from(l);
        }
    }
    let mut ss = vec![0.0f64; nf];
    pin.rows.for_each_block(|_, block| {
        for row in block.chunks_exact(nf) {
            for ((s, &m), &x) in ss.iter_mut().zip(&mean).zip(row) {
                let d = f64::from(x) - m;
                *s += d * d;
            }
        }
    });
    let inv_sd = ss
        .iter()
        .map(|&s| {
            if s <= 0.0 || n == 0.0 {
                1.0
            } else {
                1.0 / (s / n).sqrt()
            }
        })
        .collect();
    Norm { mean, inv_sd }
}

pub fn raw_weights(norm: &Norm, w: &[f64]) -> Vec<f64> {
    let nf = norm.mean.len();
    let mut out = Vec::with_capacity(nf + 1);
    let mut sum = 0.0;
    for ((&wi, &sd), &m) in w.iter().zip(&norm.inv_sd).zip(&norm.mean) {
        let r = wi * sd;
        sum += m * r;
        out.push(r);
    }
    out.push(w[nf] - sum);
    out
}

pub(crate) fn per_row_sorted<F>(rows: &Store, idx: &[usize], k: usize, f: F) -> Vec<Vec<f64>>
where
    F: Fn(&[f32], &mut [f64], &mut Vec<f64>) + Sync,
{
    let n = idx.len();
    let mut out: Vec<Vec<f64>> = (0..k).map(|_| vec![0.0f64; n]).collect();
    let nf = rows.nfeat();
    if n == 0 || rows.is_empty() || nf == 0 {
        return out;
    }
    let br = rows.block_rows();
    let nblocks = rows.len().div_ceil(br);
    let bounds: Vec<usize> = (0..=nblocks)
        .map(|b| idx.partition_point(|&r| r < b * br))
        .collect();
    let mut per_block: Vec<Vec<&mut [f64]>> = (0..nblocks).map(|_| Vec::with_capacity(k)).collect();
    for o in out.iter_mut() {
        let mut rest: &mut [f64] = o;
        for (b, slot) in per_block.iter_mut().enumerate() {
            let (head, tail) = rest.split_at_mut(bounds[b + 1] - bounds[b]);
            slot.push(head);
            rest = tail;
        }
    }
    let per_block: Vec<Mutex<Vec<&mut [f64]>>> = per_block.into_iter().map(Mutex::new).collect();
    rows.par_for_each_block(|start, block| {
        let b = start / br;
        let mut dst = per_block[b].lock().expect("score block lock");
        let mut vals = vec![0.0f64; k];
        let mut scratch = Vec::new();
        for (p, &r) in idx[bounds[b]..bounds[b + 1]].iter().enumerate() {
            let row = &block[(r - start) * nf..(r - start + 1) * nf];
            f(row, &mut vals, &mut scratch);
            for (o, &v) in dst.iter_mut().zip(&vals) {
                o[p] = v;
            }
        }
    });
    out
}

pub fn score(rows: &Store, idx: &[usize], w: &[f64], norm: Option<&Norm>) -> Vec<f64> {
    let nf = rows.nfeat();
    let one = |row: &[f32], out: &mut [f64], _: &mut Vec<f64>| {
        let mut s = w[nf];
        match norm {
            Some(nm) => {
                for (((&x, &wj), &m), &sd) in row.iter().zip(w).zip(&nm.mean).zip(&nm.inv_sd) {
                    s += wj * (f64::from(x) - m) * sd;
                }
            }
            None => {
                for (&x, &wj) in row.iter().zip(&w[..nf]).rev() {
                    s += wj * f64::from(x);
                }
            }
        }
        out[0] = s;
    };
    if idx.windows(2).all(|p| p[0] <= p[1]) {
        return per_row_sorted(rows, idx, 1, one).swap_remove(0);
    }
    let mut order: Vec<usize> = (0..idx.len()).collect();
    order.sort_by_key(|&p| idx[p]);
    let sorted: Vec<usize> = order.iter().map(|&p| idx[p]).collect();
    let got = per_row_sorted(rows, &sorted, 1, one).swap_remove(0);
    let mut out = vec![0.0f64; idx.len()];
    for (&p, v) in order.iter().zip(got) {
        out[p] = v;
    }
    out
}
