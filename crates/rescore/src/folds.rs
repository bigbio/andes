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

/// The cross-validation unit of each row: the peptide sequence (`RowMeta::sequence`) of the
/// first row of its spectrum (file and scan), which is the spectrum's rank-1 match because the
/// PIN lists a spectrum's rows best first. Every row of a spectrum gets the same unit, and
/// spectra whose best match is the same peptide (neighbouring scans of one precursor, other
/// charges, other peptidoforms) share a unit. A spectrum's key never depends on its lower
/// ranks, so units do not chain into one large group.
pub fn fold_groups(meta: &[RowMeta]) -> Vec<u32> {
    let n = meta.len();
    let mut by_spectrum: Vec<(u32, u64, u32)> = (0..n)
        .into_par_iter()
        .map(|i| (meta[i].file, meta[i].scan, i as u32))
        .collect();
    by_spectrum.par_sort_unstable();
    let mut unit = vec![0u32; n];
    let mut key = 0u32;
    for (pos, &(file, scan, i)) in by_spectrum.iter().enumerate() {
        if pos == 0 || (by_spectrum[pos - 1].0, by_spectrum[pos - 1].1) != (file, scan) {
            key = meta[i as usize].sequence;
        }
        unit[i as usize] = key;
    }
    unit
}

/// Cross-validation folds in which no precursor is split: every unit of `fold_groups` goes
/// to one fold, so no spectrum is scored by a model trained on another spectrum whose best
/// match is the same peptide sequence. Units are placed as `folds_by_spectrum` places
/// spectra: a random fold per unit, sizes balanced, no fold left empty while units remain.
pub fn folds_by_precursor(meta: &[RowMeta], k: usize, seed: u64) -> Vec<u8> {
    let k = k.max(1);
    let n = meta.len();
    let mut fold = vec![0u8; n];
    if n == 0 {
        return fold;
    }
    let unit = fold_groups(meta);
    let mut order: Vec<(u32, u32, i8, u32)> = (0..n)
        .into_par_iter()
        .map(|i| (fast_uint_hash(unit[i]), unit[i], -meta[i].label, i as u32))
        .collect();
    drop(unit);
    order.par_sort_unstable();
    let mut remain = vec![0i64; k];
    let mut left = n as i64;
    for f in (0..k).rev() {
        remain[f] = left / (f as i64 + 1);
        left -= remain[f];
    }
    let mut units_left = 1 + order.windows(2).filter(|p| p[0].1 != p[1].1).count();
    let mut empty = vec![true; k];
    let mut empty_folds = k;
    let mut rng = ParkMiller::new(seed);
    let mut previous = order[0].1;
    let mut chosen = (rng.next() % k as u64) as usize;
    for (pos, &(_, current, _, i)) in order.iter().enumerate() {
        let new_unit = previous != current;
        if new_unit {
            chosen = (rng.next() % k as u64) as usize;
            if remain.iter().any(|&r| r > 0) {
                while remain[chosen] <= 0 {
                    chosen = (rng.next() % k as u64) as usize;
                }
            }
        }
        if pos == 0 || new_unit {
            if units_left == empty_folds && !empty[chosen] {
                chosen = empty.iter().position(|&e| e).unwrap_or(chosen);
            }
            if empty[chosen] {
                empty[chosen] = false;
                empty_folds -= 1;
            }
            units_left -= 1;
        }
        fold[i as usize] = chosen as u8;
        remain[chosen] -= 1;
        previous = current;
    }
    fold
}
