use crate::pin::RowMeta;
use rayon::prelude::*;

pub(crate) fn order_key64(x: f64) -> u64 {
    let b = x.to_bits();
    if b >> 63 == 1 {
        !b
    } else {
        b | (1 << 63)
    }
}

fn from_key64(k: u64) -> f64 {
    f64::from_bits(if k >> 63 == 1 { k & !(1 << 63) } else { !k })
}

pub(crate) fn descending_order(scores: &[f64]) -> Vec<usize> {
    let mut keyed: Vec<(u64, usize)> = scores
        .par_iter()
        .enumerate()
        .map(|(i, &s)| (!order_key64(s), i))
        .collect();
    keyed.par_sort_unstable();
    keyed.into_par_iter().map(|(_, i)| i).collect()
}

pub(crate) fn radix_sort_u64(v: &mut Vec<u64>) {
    let n = v.len();
    if n < 2 {
        return;
    }
    let mut tmp = vec![0u64; n];
    for shift in (0..64).step_by(8) {
        let mut count = [0usize; 256];
        for &x in v.iter() {
            count[((x >> shift) & 0xff) as usize] += 1;
        }
        if count.contains(&n) {
            continue;
        }
        let mut sum = 0usize;
        for c in count.iter_mut() {
            let t = *c;
            *c = sum;
            sum += t;
        }
        for &x in v.iter() {
            let d = ((x >> shift) & 0xff) as usize;
            tmp[count[d]] = x;
            count[d] += 1;
        }
        std::mem::swap(v, &mut tmp);
    }
}

pub(crate) struct TdcCut {
    pub below: usize,
    pub at_most: Option<f64>,
}

/// Rows below which [`tdc_cut`] sorts all keys instead of selecting.
const SELECT_MIN: usize = 1 << 16;

/// Bits of the order key that one histogram level resolves.
const LEVEL_BITS: u32 = 16;

/// Rows per task of the parallel histogram and gather passes.
const SELECT_CHUNK: usize = 1 << 15;

/// Walks the training targets and decoys from the best score down, grouping
/// equal scores, and returns the target count of the last group whose FDR is
/// below `fdr` and the score of the last group whose FDR is at most `fdr`.
///
/// Large inputs are not sorted. A parallel histogram of the top 16 bits of
/// the order key gives, for every bucket, the counts above it; from those a
/// lower bound on the FDR of every group inside a bucket follows, because
/// correctly rounded division is monotone. Only buckets whose bound passes
/// are gathered and walked exactly, from the lowest score up, so the result
/// is the one the full sort gives, ties included.
pub(crate) fn tdc_cut(
    scores: &[f64],
    is_decoy: &[bool],
    fdr: f64,
    skip_decoys_plus_one: bool,
) -> TdcCut {
    let n = scores.len().min(is_decoy.len());
    if n < SELECT_MIN {
        return tdc_cut_sorted(scores, is_decoy, fdr, skip_decoys_plus_one);
    }
    let (scores, is_decoy) = (&scores[..n], &is_decoy[..n]);
    let offset: u64 = if skip_decoys_plus_one { 0 } else { 1 };
    let sel = Selector {
        scores,
        is_decoy,
        offset,
        hist: histogram(scores, is_decoy),
    };
    let below = sel.last_group(|f| f < fdr).map_or(0, |(t, _)| t);
    let at_most = sel.last_group(|f| f <= fdr).map(|(_, v)| v);
    TdcCut {
        below: below as usize,
        at_most,
    }
}

const BUCKETS: usize = 1 << LEVEL_BITS;

fn bucket_of(s: f64) -> usize {
    (order_key64(s + 0.0) >> (64 - LEVEL_BITS)) as usize
}

fn histogram(scores: &[f64], is_decoy: &[bool]) -> Vec<[u64; 2]> {
    let chunk = scores
        .len()
        .div_ceil(rayon::current_num_threads().max(1))
        .max(SELECT_CHUNK);
    let parts: Vec<Vec<[u32; 2]>> = scores
        .par_chunks(chunk)
        .zip(is_decoy.par_chunks(chunk))
        .map(|(s, d)| {
            let mut h = vec![[0u32; 2]; BUCKETS];
            for (&v, &dd) in s.iter().zip(d) {
                h[bucket_of(v)][usize::from(dd)] += 1;
            }
            h
        })
        .collect();
    let mut hist = vec![[0u64; 2]; BUCKETS];
    for p in &parts {
        for (x, y) in hist.iter_mut().zip(p) {
            x[0] += u64::from(y[0]);
            x[1] += u64::from(y[1]);
        }
    }
    hist
}

struct Selector<'a> {
    scores: &'a [f64],
    is_decoy: &'a [bool],
    offset: u64,
    hist: Vec<[u64; 2]>,
}

impl Selector<'_> {
    fn last_group(&self, pass: impl Fn(f64) -> bool) -> Option<(u64, f64)> {
        let mut t_above: u64 = self.hist.iter().map(|h| h[0]).sum();
        let mut d_above: u64 = self.hist.iter().map(|h| h[1]).sum();
        for (b, h) in self.hist.iter().enumerate() {
            if h[0] + h[1] == 0 {
                continue;
            }
            let t_through = t_above;
            t_above -= h[0];
            d_above -= h[1];
            let decoys = d_above + self.offset;
            let bound = (decoys as f64 / t_through.max(1) as f64).min(1.0);
            if !pass(bound) {
                continue;
            }
            if let Some(hit) = self.walk_bucket(b, t_above, decoys, &pass) {
                return Some(hit);
            }
        }
        None
    }

    fn walk_bucket(
        &self,
        b: usize,
        mut targets: u64,
        mut decoys: u64,
        pass: &impl Fn(f64) -> bool,
    ) -> Option<(u64, f64)> {
        let (mut t, mut d): (Vec<u64>, Vec<u64>) = self
            .scores
            .par_chunks(SELECT_CHUNK)
            .zip(self.is_decoy.par_chunks(SELECT_CHUNK))
            .map(|(s, dd)| {
                let mut tt = Vec::new();
                let mut ds = Vec::new();
                for (&v, &isd) in s.iter().zip(dd) {
                    if bucket_of(v) == b {
                        if isd {
                            ds.push(order_key64(v));
                        } else {
                            tt.push(order_key64(v));
                        }
                    }
                }
                (tt, ds)
            })
            .reduce(
                || (Vec::new(), Vec::new()),
                |mut a, c| {
                    a.0.extend(c.0);
                    a.1.extend(c.1);
                    a
                },
            );
        t.par_sort_unstable();
        d.par_sort_unstable();
        let mut hit = None;
        let (mut i, mut j) = (t.len(), d.len());
        while i > 0 || j > 0 {
            let from_t = j == 0 || (i > 0 && t[i - 1] >= d[j - 1]);
            let v = if from_t {
                from_key64(t[i - 1])
            } else {
                from_key64(d[j - 1])
            };
            if from_t {
                targets += 1;
                i -= 1;
            } else {
                decoys += 1;
                j -= 1;
            }
            while i > 0 && from_key64(t[i - 1]) == v {
                targets += 1;
                i -= 1;
            }
            while j > 0 && from_key64(d[j - 1]) == v {
                decoys += 1;
                j -= 1;
            }
            let f = (decoys as f64 / targets.max(1) as f64).min(1.0);
            if pass(f) {
                hit = Some((targets, v));
            }
        }
        hit
    }
}

fn tdc_cut_sorted(
    scores: &[f64],
    is_decoy: &[bool],
    fdr: f64,
    skip_decoys_plus_one: bool,
) -> TdcCut {
    let mut t: Vec<u64> = Vec::new();
    let mut d: Vec<u64> = Vec::new();
    for (&s, &dd) in scores.iter().zip(is_decoy) {
        if dd {
            d.push(order_key64(s));
        } else {
            t.push(order_key64(s));
        }
    }
    radix_sort_u64(&mut t);
    radix_sort_u64(&mut d);
    let mut decoys: u64 = if skip_decoys_plus_one { 0 } else { 1 };
    let mut targets: u64 = 0;
    let mut below = 0u64;
    let mut at_most = None;
    let (mut i, mut j) = (t.len(), d.len());
    while i > 0 || j > 0 {
        let from_t = j == 0 || (i > 0 && t[i - 1] >= d[j - 1]);
        let v = if from_t {
            from_key64(t[i - 1])
        } else {
            from_key64(d[j - 1])
        };
        if from_t {
            targets += 1;
            i -= 1;
        } else {
            decoys += 1;
            j -= 1;
        }
        while i > 0 && from_key64(t[i - 1]) == v {
            targets += 1;
            i -= 1;
        }
        while j > 0 && from_key64(d[j - 1]) == v {
            decoys += 1;
            j -= 1;
        }
        let f = (decoys as f64 / targets.max(1) as f64).min(1.0);
        if f < fdr {
            below = targets;
        }
        if f <= fdr {
            at_most = Some(v);
        }
    }
    TdcCut {
        below: below as usize,
        at_most,
    }
}

pub fn tdc_q(scores: &[f64], is_decoy: &[bool], skip_decoys_plus_one: bool) -> Vec<f64> {
    let n = scores.len();
    let order = descending_order(scores);
    let mut fdr = vec![0.0f64; n];
    let mut decoys: u64 = if skip_decoys_plus_one { 0 } else { 1 };
    let mut targets: u64 = 0;
    let mut group_start = 0usize;
    for pos in 0..n {
        if is_decoy[order[pos]] {
            decoys += 1;
        } else {
            targets += 1;
        }
        if pos + 1 == n || scores[order[pos]] != scores[order[pos + 1]] {
            let f = (decoys as f64 / targets.max(1) as f64).min(1.0);
            for &i in &order[group_start..=pos] {
                fdr[i] = f;
            }
            group_start = pos + 1;
        }
    }
    let mut run = f64::INFINITY;
    for &i in order.iter().rev() {
        run = run.min(fdr[i]);
        fdr[i] = run;
    }
    fdr
}

pub fn count_at(q: &[f64], is_decoy: &[bool], fdr: f64) -> usize {
    q.iter()
        .zip(is_decoy)
        .filter(|&(&q, &d)| !d && q < fdr)
        .count()
}

pub fn training_sets(q: &[f64], is_decoy: &[bool], fdr: f64) -> (Vec<usize>, Vec<usize>) {
    let mut positives = Vec::new();
    let mut negatives = Vec::new();
    for (i, (&qi, &d)) in q.iter().zip(is_decoy).enumerate() {
        if d {
            negatives.push(i);
        } else if qi <= fdr {
            positives.push(i);
        }
    }
    (positives, negatives)
}

pub fn tdc_winners<F>(meta: &[RowMeta], scores: &[f64], psm_id: F) -> Result<Vec<usize>, String>
where
    F: Fn(usize) -> Result<String, String>,
{
    let mut keyed: Vec<(u32, u64, u64, u64, i8, usize)> = meta
        .par_iter()
        .enumerate()
        .map(|(i, m)| {
            (
                m.file,
                m.scan,
                order_key64(m.exp_mass + 0.0),
                !order_key64(scores[i] + 0.0),
                m.label,
                i,
            )
        })
        .collect();
    keyed.par_sort_unstable();
    let mut out = Vec::new();
    for group in keyed.chunk_by(|a, b| (a.0, a.1, a.2) == (b.0, b.1, b.2)) {
        let tied = group
            .iter()
            .take_while(|k| (k.3, k.4) == (group[0].3, group[0].4))
            .count();
        let mut best = group[0].5;
        if tied > 1 {
            let mut best_id = psm_id(best)?;
            for k in &group[1..tied] {
                let id = psm_id(k.5)?;
                if id < best_id {
                    best = k.5;
                    best_id = id;
                }
            }
        }
        out.push(best);
    }
    out.par_sort_unstable();
    Ok(out)
}

pub(crate) fn psm_order<F>(
    meta: &[RowMeta],
    rows: &[usize],
    scores: &[f64],
    psm_id: F,
) -> Result<Vec<usize>, String>
where
    F: Fn(usize) -> Result<String, String>,
{
    let key = |k: usize| {
        let m = &meta[rows[k]];
        (
            !order_key64(scores[k] + 0.0),
            !m.scan,
            !order_key64(m.exp_mass + 0.0),
            -m.label,
            !m.file,
        )
    };
    let mut order: Vec<usize> = (0..rows.len()).collect();
    order.par_sort_unstable_by(|&a, &b| key(a).cmp(&key(b)).then(rows[b].cmp(&rows[a])));
    let mut start = 0usize;
    while start < order.len() {
        let mut end = start + 1;
        while end < order.len() && key(order[end]) == key(order[start]) {
            end += 1;
        }
        if end - start > 1 {
            let mut run: Vec<(String, usize)> = order[start..end]
                .iter()
                .map(|&k| psm_id(rows[k]).map(|id| (id, k)))
                .collect::<Result<_, _>>()?;
            run.sort_by(|a, b| b.0.cmp(&a.0).then(rows[b.1].cmp(&rows[a.1])));
            for (slot, (_, k)) in order[start..end].iter_mut().zip(run) {
                *slot = k;
            }
        }
        start = end;
    }
    Ok(order)
}

pub fn mix_max_q(scores: &[f64], is_decoy: &[bool], pi0: f64) -> Vec<f64> {
    if pi0 >= 1.0 {
        return tdc_q(scores, is_decoy, false);
    }
    let n = scores.len();
    let order = descending_order(scores);
    let combined: Vec<(f64, bool)> = order.iter().map(|&i| (scores[i], !is_decoy[i])).collect();
    let (h_w, h_z) = crate::pep::mix_max_counts(&combined);
    let mut q = vec![0.0f64; n];
    let mut n_z: i64 = 1;
    let mut n_w: i64 = 0;
    let mut decoy_queue: i64 = 0;
    let mut correction = 0.0f64;
    let mut group_start = 0usize;
    for pos in 0..n {
        if combined[pos].1 {
            n_w += 1;
        } else {
            n_z += 1;
            decoy_queue += 1;
        }
        if pos + 1 == n || combined[pos].0 != combined[pos + 1].0 {
            if decoy_queue > 0 {
                let j = h_w.len() - (n_z - 1) as usize;
                let cnt_w = h_w[j].trunc();
                let cnt_z = h_z[j].trunc();
                let est = ((cnt_w - pi0 * cnt_z) / ((1.0 - pi0) * cnt_z)).clamp(0.0, 1.0);
                correction += decoy_queue as f64 * est * (1.0 - pi0);
            }
            let fdr = ((n_z as f64 * pi0 + correction) / n_w.max(1) as f64).min(1.0);
            for &i in &order[group_start..=pos] {
                q[i] = fdr;
            }
            group_start = pos + 1;
            decoy_queue = 0;
        }
    }
    let mut run = f64::INFINITY;
    for &i in order.iter().rev() {
        run = run.min(q[i]);
        q[i] = run;
    }
    q
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lcg(state: &mut u64) -> u64 {
        *state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        *state >> 33
    }

    fn same(scores: &[f64], decoy: &[bool], fdr: f64, skip: bool) {
        let a = tdc_cut(scores, decoy, fdr, skip);
        let b = tdc_cut_sorted(scores, decoy, fdr, skip);
        assert_eq!(a.below, b.below, "below at fdr {fdr} skip {skip}");
        assert_eq!(
            a.at_most.map(f64::to_bits),
            b.at_most.map(f64::to_bits),
            "at_most at fdr {fdr} skip {skip}"
        );
    }

    #[test]
    fn selected_cut_equals_the_sorted_walk_ties_and_signed_zeros_included() {
        let mut st = 7u64;
        for (case, n) in [70_000usize, 150_000, 300_000].into_iter().enumerate() {
            let mut scores = Vec::with_capacity(n);
            let mut decoy = Vec::with_capacity(n);
            for _ in 0..n {
                let d = lcg(&mut st) % 3 == 0;
                let u = (lcg(&mut st) % 1_000_000) as f64 / 1_000_000.0;
                let v = match case {
                    0 => (u * 64.0).round() / 16.0 - if d { 2.5 } else { 0.5 },
                    1 => u * 6.0 - if d { 4.0 } else { 1.0 },
                    _ => {
                        if lcg(&mut st) % 50 == 0 {
                            if lcg(&mut st) % 2 == 0 {
                                0.0
                            } else {
                                -0.0
                            }
                        } else {
                            (u - 0.5) * 1e-3 + if d { -1e-4 } else { 1e-4 }
                        }
                    }
                };
                scores.push(v);
                decoy.push(d);
            }
            for fdr in [0.0, 1e-4, 0.001, 0.01, 0.05, 0.3, 0.5, 1.0] {
                same(&scores, &decoy, fdr, true);
                same(&scores, &decoy, fdr, false);
            }
        }
    }

    #[test]
    fn selected_cut_on_constant_and_one_sided_inputs() {
        let n = 80_000;
        let alt: Vec<bool> = (0..n).map(|i| i % 4 == 0).collect();
        for fdr in [0.01, 0.25, 0.3, 1.0] {
            same(&vec![1.5; n], &alt, fdr, true);
            same(&vec![1.5; n], &alt, fdr, false);
            same(&vec![-0.0; n], &alt, fdr, true);
            let s: Vec<f64> = (0..n).map(|i| i as f64).collect();
            same(&s, &vec![false; n], fdr, false);
            same(&s, &vec![true; n], fdr, true);
            same(&s, &alt, fdr, true);
            let r: Vec<f64> = (0..n).map(|i| -(i as f64)).collect();
            same(&r, &alt, fdr, false);
        }
    }
}
