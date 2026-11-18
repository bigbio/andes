use crate::pin::RowMeta;
use rayon::prelude::*;

const PARK_MILLER_MODULUS: u64 = 4_294_967_291;

struct ParkMiller(u64);

impl ParkMiller {
    fn new(seed: u64) -> ParkMiller {
        match seed % PARK_MILLER_MODULUS {
            0 => ParkMiller(1),
            s => ParkMiller(s),
        }
    }

    fn next(&mut self) -> u64 {
        self.0 = (self.0 * 279_470_273) % PARK_MILLER_MODULUS;
        self.0
    }
}

fn fast_uint_hash(x: u32) -> u32 {
    let x = ((x >> 16) ^ x).wrapping_mul(0x045d_9f3b);
    let x = ((x >> 16) ^ x).wrapping_mul(0x045d_9f3b);
    (x >> 16) ^ x
}

fn scan_hash(m: &RowMeta) -> u32 {
    fast_uint_hash(m.file) ^ fast_uint_hash(m.scan as u32).wrapping_shl(1)
}

pub fn folds_by_spectrum(meta: &[RowMeta], k: usize, seed: u64) -> Vec<u8> {
    let k = k.max(1);
    let n = meta.len();
    let mut fold = vec![0u8; n];
    if n == 0 {
        return fold;
    }
    let mut order: Vec<(u32, u32, u32, i8, u32)> = (0..n)
        .into_par_iter()
        .map(|i| {
            let m = &meta[i];
            (scan_hash(m), m.file, m.scan as u32, -m.label, i as u32)
        })
        .collect();
    order.par_sort_unstable();
    let mut remain = vec![0i64; k];
    let mut left = n as i64;
    for f in (0..k).rev() {
        remain[f] = left / (f as i64 + 1);
        left -= remain[f];
    }
    let mut spectra_left = 1 + order.windows(2).filter(|p| p[0].2 != p[1].2).count();
    let mut empty = vec![true; k];
    let mut empty_folds = k;
    let mut rng = ParkMiller::new(seed);
    let mut previous = order[0].2;
    let mut chosen = (rng.next() % k as u64) as usize;
    for (pos, &(_, _, current, _, i)) in order.iter().enumerate() {
        let new_spectrum = previous != current;
        if new_spectrum {
            chosen = (rng.next() % k as u64) as usize;
            while remain[chosen] <= 0 {
                chosen = (rng.next() % k as u64) as usize;
            }
        }
        if pos == 0 || new_spectrum {
            if spectra_left == empty_folds && !empty[chosen] {
                chosen = empty.iter().position(|&e| e).unwrap_or(chosen);
            }
            if empty[chosen] {
                empty[chosen] = false;
                empty_folds -= 1;
            }
            spectra_left -= 1;
        }
        fold[i as usize] = chosen as u8;
        remain[chosen] -= 1;
        previous = current;
    }
    fold
}
