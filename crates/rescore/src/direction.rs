use crate::budget::Reservation;
use crate::normalize::Norm;
use crate::pin::Pin;
use rayon::prelude::*;
use std::sync::Mutex;

const NO_DIRECTION: &str = "Error in the input data: cannot find an initial direction with positive training examples. Consider setting/raising the initial training FDR threshold (--train-fdr-initial). Terminating.";

/// Bytes per row of one feature's order keys, sorted in place.
const SORT_BYTES_PER_ROW: u64 = 8;

pub(crate) fn pin_direction(pin: &Pin) -> Option<&[f64]> {
    pin.default_direction
        .as_deref()
        .filter(|d| d.iter().any(|&v| v != 0.0))
}

pub fn init_direction(
    pin: &Pin,
    norm: &Norm,
    rows: &[usize],
    fdr: f64,
    overrule: bool,
) -> Result<Vec<f64>, String> {
    let mut mask = vec![0u8; pin.meta.len()];
    for &r in rows {
        mask[r] = 1;
    }
    let mut w = init_directions(pin, norm, &mask, 1, fdr, overrule)?;
    Ok(w.remove(0))
}

pub(crate) fn init_directions(
    pin: &Pin,
    norm: &Norm,
    mask: &[u8],
    sets: usize,
    fdr: f64,
    overrule: bool,
) -> Result<Vec<Vec<f64>>, String> {
    let nf = pin.rows.nfeat();
    if let Some(d) = pin_direction(pin) {
        let mut w = vec![0.0f64; nf + 1];
        for (wi, &di) in w.iter_mut().zip(d) {
            *wi = di;
        }
        return Ok(vec![w; sets]);
    }
    if nf == 0 {
        return Ok(vec![vec![0.0f64; 1]; sets]);
    }
    let counts = direction_counts(pin, norm, mask, sets, fdr)?;
    counts
        .iter()
        .map(|c| {
            let mut w = vec![0.0f64; nf + 1];
            let mut best: Option<(usize, usize, bool)> = None;
            for (j, pair) in c.iter().enumerate() {
                for (i, &p) in pair.iter().enumerate() {
                    let better = match best {
                        None => true,
                        Some((b, _, _)) => p > b,
                    };
                    if better {
                        best = Some((p, j, i == 0));
                    }
                }
            }
            if let Some((p, j, low)) = best {
                if p == 0 && !overrule {
                    return Err(NO_DIRECTION.to_string());
                }
                w[j] = if low { -1.0 } else { 1.0 };
            }
            Ok(w)
        })
        .collect()
}

fn order_key(x: f32) -> u32 {
    let b = x.to_bits();
    if b >> 31 == 1 {
        !b
    } else {
        b | 0x8000_0000
    }
}

fn from_key(k: u32) -> f32 {
    f32::from_bits(if k >> 31 == 1 { k & 0x7fff_ffff } else { !k })
}

fn close_group(decoys: &[u64], targets: &[u64], best: &mut [u64], fdr: f64) {
    for ((&d, &t), b) in decoys.iter().zip(targets).zip(best.iter_mut()) {
        let f = (d as f64 / t.max(1) as f64).min(1.0);
        if f < fdr {
            *b = t;
        }
    }
}

fn walk<'a>(
    col: impl Iterator<Item = &'a u64>,
    mean: f64,
    inv_sd: f64,
    sets: usize,
    fdr: f64,
) -> Vec<u64> {
    let mut decoys = vec![0u64; sets];
    let mut targets = vec![0u64; sets];
    let mut best = vec![0u64; sets];
    let mut prev: Option<f64> = None;
    for &e in col {
        let v = (f64::from(from_key((e >> 32) as u32)) - mean) * inv_sd;
        if prev.is_some_and(|p| p != v) {
            close_group(&decoys, &targets, &mut best, fdr);
        }
        prev = Some(v);
        let member = (e >> 1) & 0x7f;
        let decoy = e & 1 == 1;
        for s in 0..sets {
            if (member >> s) & 1 == 1 {
                if decoy {
                    decoys[s] += 1;
                } else {
                    targets[s] += 1;
                }
            }
        }
    }
    if prev.is_some() {
        close_group(&decoys, &targets, &mut best, fdr);
    }
    best
}

/// One feature's (low, high) direction counts per training set, from its order keys.
fn column_counts(
    mut col: Vec<u64>,
    mean: f64,
    inv_sd: f64,
    sets: usize,
    fdr: f64,
) -> Vec<[usize; 2]> {
    col.par_sort_unstable();
    let high = walk(col.iter().rev(), mean, inv_sd, sets, fdr);
    let low = walk(col.iter(), mean, inv_sd, sets, fdr);
    (0..sets)
        .map(|s| [low[s] as usize, high[s] as usize])
        .collect()
}

/// For each training set and feature, the positives the feature separates at `fdr` read low
/// and high. Features are sorted as many at a time as the budget holds; a spilled store is
/// first written as one file per column, so each further pass reads only its features.
fn direction_counts(
    pin: &Pin,
    norm: &Norm,
    mask: &[u8],
    sets: usize,
    fdr: f64,
) -> Result<Vec<Vec<[usize; 2]>>, String> {
    let nf = pin.rows.nfeat();
    let n = pin.meta.len();
    let flags: Vec<u8> = mask
        .iter()
        .zip(&pin.meta)
        .map(|(&m, meta)| {
            if m == 0 {
                0
            } else {
                (m << 1) | u8::from(meta.label < 0)
            }
        })
        .collect();
    let per_feature = (n as u64).max(1) * SORT_BYTES_PER_ROW;
    let chunk = ((pin.budget.available() / per_feature) as usize).clamp(1, nf);
    let columns = if chunk < nf && pin.rows.is_spilled() {
        Some(pin.rows.columns()?)
    } else {
        None
    };
    let br = pin.rows.block_rows();
    let mut out = vec![vec![[0usize; 2]; nf]; sets];
    for j0 in (0..nf).step_by(chunk) {
        let j1 = (j0 + chunk).min(nf);
        let _held = Reservation::take(&pin.budget, (j1 - j0) as u64 * per_feature);
        let res: Vec<Vec<[usize; 2]>> = match &columns {
            Some(columns) => (j0..j1)
                .into_par_iter()
                .map(|j| {
                    let mut col = Vec::with_capacity(n);
                    columns.for_each_value(j, |x| {
                        let fl = flags.get(col.len()).copied().unwrap_or(0);
                        col.push((u64::from(order_key(x)) << 32) | u64::from(fl));
                    })?;
                    if col.len() != n {
                        return Err(format!(
                            "column file of feature {}: {} values for {n} rows",
                            j + 1,
                            col.len()
                        ));
                    }
                    Ok(column_counts(col, norm.mean[j], norm.inv_sd[j], sets, fdr))
                })
                .collect::<Result<_, String>>()?,
            None => {
                let mut cols: Vec<Vec<u64>> = (j0..j1).map(|_| vec![0u64; n]).collect();
                {
                    let mut per_block: Vec<Vec<&mut [u64]>> = (0..n.div_ceil(br))
                        .map(|_| Vec::with_capacity(j1 - j0))
                        .collect();
                    for col in cols.iter_mut() {
                        for (b, c) in col.chunks_mut(br).enumerate() {
                            per_block[b].push(c);
                        }
                    }
                    let per_block: Vec<Mutex<Vec<&mut [u64]>>> =
                        per_block.into_iter().map(Mutex::new).collect();
                    pin.rows.par_for_each_block(|start, block| {
                        let mut dst = per_block[start / br]
                            .lock()
                            .expect("direction column block lock");
                        for (r, row) in block.chunks_exact(nf).enumerate() {
                            let fl = u64::from(flags[start + r]);
                            for (d, &x) in dst.iter_mut().zip(&row[j0..j1]) {
                                d[r] = (u64::from(order_key(x)) << 32) | fl;
                            }
                        }
                    });
                }
                cols.into_par_iter()
                    .enumerate()
                    .map(|(k, col)| {
                        column_counts(col, norm.mean[j0 + k], norm.inv_sd[j0 + k], sets, fdr)
                    })
                    .collect()
            }
        };
        for (k, r) in res.into_iter().enumerate() {
            for (s, pair) in r.into_iter().enumerate() {
                out[s][j0 + k] = pair;
            }
        }
    }
    Ok(out)
}
