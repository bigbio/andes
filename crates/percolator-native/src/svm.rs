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
        if !idx.is_empty() && !store.is_empty() {
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

pub fn l2_svm_mfn(
    x: &Rows,
    idx_pos: &[usize],
    idx_neg: &[usize],
    cpos: f64,
    cneg: f64,
    opts: &MfnOpts,
) -> Vec<f64> {
    let _ = (idx_pos, idx_neg, cpos, cneg, opts);
    vec![0.0; x.ncols]
}
