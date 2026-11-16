use crate::pin::Pin;
use crate::store::Store;

pub struct Norm {
    pub mean: Vec<f64>,
    pub inv_sd: Vec<f64>,
}

pub fn stdv(pin: &Pin) -> Norm {
    let nf = pin.rows.nfeat();
    let n = pin.rows.len() as f64;
    let mut mean = vec![0.0f64; nf];
    pin.rows.for_each_block(|_, block| {
        for row in block.chunks_exact(nf) {
            for (m, &x) in mean.iter_mut().zip(row) {
                *m += f64::from(x);
            }
        }
    });
    if n > 0.0 {
        for m in mean.iter_mut() {
            *m /= n;
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

pub fn score(rows: &Store, idx: &[usize], w: &[f64], norm: Option<&Norm>) -> Vec<f64> {
    let nf = rows.nfeat();
    let mut out = vec![0.0f64; idx.len()];
    if idx.is_empty() || rows.is_empty() {
        return out;
    }
    let mut order: Vec<usize> = (0..idx.len()).collect();
    order.sort_by_key(|&p| idx[p]);
    let mut next = 0usize;
    rows.for_each_block(|start, block| {
        let n = if nf == 0 { 0 } else { block.len() / nf };
        let end = start + n;
        while next < order.len() && idx[order[next]] < end {
            let p = order[next];
            let r = idx[p] - start;
            let row = &block[r * nf..(r + 1) * nf];
            let mut s = w[nf];
            match norm {
                Some(nm) => {
                    for (((&x, &wj), &m), &sd) in row.iter().zip(w).zip(&nm.mean).zip(&nm.inv_sd) {
                        s += wj * (f64::from(x) - m) * sd;
                    }
                }
                None => {
                    for (&x, &wj) in row.iter().zip(w) {
                        s += wj * f64::from(x);
                    }
                }
            }
            out[p] = s;
            next += 1;
        }
    });
    out
}
