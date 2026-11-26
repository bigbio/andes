//! The layout of spilled training rows. A feature whose values are all integers within a span
//! of 256 is written as one byte, an index into the table of its normalised values; every other
//! feature is written as its normalised f32. The bias column and the padding are not written.
//! Decoding gives back, bit for bit, the normalised row that the rows in memory hold, with
//! fewer bytes to read in each pass over spilled rows.

use crate::normalize::Norm;
use crate::store::Store;
use std::sync::Mutex;

/// Values a feature written as one byte may take.
const CODES: usize = 256;

/// Integers from which an f32 no longer holds every integer.
const EXACT_INTEGERS: f32 = 16_777_216.0;

/// A feature's normalised value, as the training rows hold it.
pub(crate) fn normalised(x: f32, mean: f64, inv_sd: f64) -> f32 {
    ((f64::from(x) - mean) * inv_sd) as f32
}

fn integer(x: f32) -> bool {
    x.abs() < EXACT_INTEGERS && (x as i32 as f32).to_bits() == x.to_bits()
}

/// A feature's lowest and highest value, and whether every value is an integer.
#[derive(Clone, Copy)]
struct Span {
    lo: f32,
    hi: f32,
    integers: bool,
}

impl Span {
    fn new() -> Span {
        Span {
            lo: f32::INFINITY,
            hi: f32::NEG_INFINITY,
            integers: true,
        }
    }

    fn add(&mut self, x: f32) {
        if !self.integers {
            return;
        }
        if integer(x) {
            self.lo = self.lo.min(x);
            self.hi = self.hi.max(x);
        } else {
            self.integers = false;
        }
    }

    fn merge(&mut self, o: &Span) {
        self.integers &= o.integers;
        self.lo = self.lo.min(o.lo);
        self.hi = self.hi.max(o.hi);
    }

    fn coded(&self) -> bool {
        self.integers && self.lo <= self.hi && self.hi - self.lo < CODES as f32
    }
}

pub(crate) struct RowCodec {
    nfeat: usize,
    stride: usize,
    coded: Vec<usize>,
    lo: Vec<f32>,
    tables: Vec<Vec<f32>>,
    plain: Vec<usize>,
}

impl RowCodec {
    /// The layout of `store`'s rows normalised with `norm`, from one pass over the store.
    pub(crate) fn scan(store: &Store, norm: &Norm, stride: usize) -> RowCodec {
        let nf = store.nfeat();
        let all = Mutex::new(vec![Span::new(); nf]);
        store.par_for_each_block(|_, block| {
            let mut part = vec![Span::new(); nf];
            for row in block.chunks_exact(nf.max(1)) {
                for (s, &x) in part.iter_mut().zip(row) {
                    s.add(x);
                }
            }
            let mut a = all.lock().unwrap_or_else(|p| p.into_inner());
            for (s, p) in a.iter_mut().zip(&part) {
                s.merge(p);
            }
        });
        let spans = all.into_inner().unwrap_or_else(|p| p.into_inner());
        let mut c = RowCodec {
            nfeat: nf,
            stride,
            coded: Vec::new(),
            lo: Vec::new(),
            tables: Vec::new(),
            plain: Vec::new(),
        };
        for (j, s) in spans.iter().enumerate() {
            if !s.coded() {
                c.plain.push(j);
                continue;
            }
            let n = (s.hi - s.lo) as usize + 1;
            c.coded.push(j);
            c.lo.push(s.lo);
            c.tables.push(
                (0..n)
                    .map(|k| normalised(s.lo + k as f32, norm.mean[j], norm.inv_sd[j]))
                    .collect(),
            );
        }
        c
    }

    /// Bytes of one written row.
    pub(crate) fn width(&self) -> usize {
        self.coded.len() + 4 * self.plain.len()
    }

    /// Features written as one byte.
    pub(crate) fn coded(&self) -> usize {
        self.coded.len()
    }

    /// A row before decoding: zeros, the bias 1 after the features.
    pub(crate) fn blank(&self) -> Vec<f32> {
        let mut row = vec![0.0f32; self.stride];
        if self.nfeat < self.stride {
            row[self.nfeat] = 1.0;
        }
        row
    }

    /// Writes the row of raw features `raw`, normalised as `normed`, to `out` (`width` bytes).
    pub(crate) fn encode(&self, raw: &[f32], normed: &[f32], out: &mut [u8]) {
        for ((o, &j), &lo) in out.iter_mut().zip(&self.coded).zip(&self.lo) {
            *o = (raw[j] - lo) as u8;
        }
        let at = self.coded.len();
        for (c, &j) in out[at..].chunks_exact_mut(4).zip(&self.plain) {
            c.copy_from_slice(&normed[j].to_le_bytes());
        }
    }

    /// Sets the features of `row` (from `blank`) to those of the written row `bytes`.
    pub(crate) fn decode(&self, bytes: &[u8], row: &mut [f32]) {
        for ((&b, &j), t) in bytes.iter().zip(&self.coded).zip(&self.tables) {
            row[j] = t[usize::from(b)];
        }
        let at = self.coded.len();
        for (c, &j) in bytes[at..].chunks_exact(4).zip(&self.plain) {
            row[j] = f32::from_le_bytes([c[0], c[1], c[2], c[3]]);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::budget::Budget;

    #[test]
    fn decoding_gives_back_the_bits_of_every_normalised_row() {
        let dir = tempfile::tempdir().unwrap();
        let budget = Budget::new(u64::MAX, dir.path());
        let nf = 6;
        let mut b = crate::store::StoreBuilder::new(nf, &budget);
        let mut raw = Vec::new();
        for r in 0..5000usize {
            let row = [
                (r % 7) as f32,
                (r % 300) as f32,
                if r == 17 { -0.0 } else { (r % 3) as f32 },
                if r == 99 { f32::NAN } else { (r % 5) as f32 - 2.0 },
                3.0,
                (r % 11) as f32 * 0.5 - 1.0,
            ];
            raw.push(row);
            b.push_row(&row).unwrap();
        }
        let store = b.finish().unwrap();
        let norm = Norm {
            mean: vec![1.3, 140.0, 0.9, -0.1, 3.0, 0.7],
            inv_sd: vec![0.37, 0.011, 1.7, 0.66, 1.0, 2.3],
        };
        let stride = crate::svm::stride_for(nf + 1);
        let c = RowCodec::scan(&store, &norm, stride);
        assert_eq!(c.coded, vec![0, 4]);
        assert_eq!(c.width(), 2 + 4 * 4);
        let mut bytes = vec![0u8; c.width()];
        let mut got = c.blank();
        for row in &raw {
            let mut want = vec![0.0f32; stride];
            want[nf] = 1.0;
            for j in 0..nf {
                want[j] = normalised(row[j], norm.mean[j], norm.inv_sd[j]);
            }
            c.encode(row, &want, &mut bytes);
            c.decode(&bytes, &mut got);
            let bits = |v: &[f32]| v.iter().map(|x| x.to_bits()).collect::<Vec<u32>>();
            assert_eq!(bits(&got), bits(&want));
        }
    }
}
