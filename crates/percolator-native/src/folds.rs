use crate::pin::RowMeta;
use std::collections::HashMap;

struct SplitMix(u64);

impl SplitMix {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
}

pub fn folds_by_spectrum(meta: &[RowMeta], k: usize, seed: u64) -> Vec<u8> {
    let k = k.max(1);
    let mut ids: HashMap<(u32, u64), usize> = HashMap::new();
    let mut sizes: Vec<usize> = Vec::new();
    let mut row_spec = Vec::with_capacity(meta.len());
    for m in meta {
        let next = ids.len();
        let id = *ids.entry((m.file, m.scan)).or_insert(next);
        if id == sizes.len() {
            sizes.push(0);
        }
        sizes[id] += 1;
        row_spec.push(id);
    }
    let mut order: Vec<usize> = (0..sizes.len()).collect();
    let mut rng = SplitMix(seed);
    for i in (1..order.len()).rev() {
        let j = (rng.next() % (i as u64 + 1)) as usize;
        order.swap(i, j);
    }
    let mut load = vec![0usize; k];
    let mut fold_of = vec![0u8; sizes.len()];
    for &s in &order {
        let f = (0..k).min_by_key(|&f| (load[f], f)).unwrap_or(0);
        load[f] += sizes[s];
        fold_of[s] = f as u8;
    }
    row_spec.iter().map(|&s| fold_of[s]).collect()
}
