use crate::pin::RowMeta;
use std::collections::HashMap;

pub fn tdc_q(scores: &[f64], is_decoy: &[bool], skip_decoys_plus_one: bool) -> Vec<f64> {
    let n = scores.len();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| scores[b].total_cmp(&scores[a]));
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

pub fn tdc_winners(meta: &[RowMeta], scores: &[f64]) -> Vec<usize> {
    let mut best: HashMap<(u32, u64, u64), usize> = HashMap::new();
    for (i, m) in meta.iter().enumerate() {
        let key = (m.file, m.scan, m.exp_mass.to_bits());
        match best.get(&key) {
            Some(&b)
                if scores[b]
                    .total_cmp(&scores[i])
                    .then(meta[b].label.cmp(&m.label))
                    .is_ge() => {}
            _ => {
                best.insert(key, i);
            }
        }
    }
    let mut out: Vec<usize> = best.into_values().collect();
    out.sort_unstable();
    out
}

pub fn mix_max_q(scores: &[f64], is_decoy: &[bool], pi0: f64) -> Vec<f64> {
    if pi0 >= 1.0 {
        return tdc_q(scores, is_decoy, false);
    }
    let n = scores.len();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| scores[b].total_cmp(&scores[a]).then(a.cmp(&b)));
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
