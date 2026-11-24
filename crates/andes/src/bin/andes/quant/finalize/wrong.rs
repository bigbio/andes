//! Calibration of wrong-transfer probabilities from pseudo-transfers of identified precursors.

use quant::decoys::Hypothesis;
use quant::evidence::Evidence;

use super::extract::BlockContext;
use super::index::{EntryIndex, Table};
use crate::quant::mbr::{amplitude_snr, ExtractRow, Role};
use crate::transfers::{FoldModels, TARGET};

const WRONG_FEATURES: usize = 11;
pub(super) const WRONG_CALIBRATION_ROWS: usize = 4_000_000;

pub(crate) fn wrong_features(e: &Evidence) -> Option<Vec<f32>> {
    if e.quantity.is_nan() || e.quantity <= 0.0 {
        return None;
    }
    let x = [
        e.llr,
        e.rt_z.abs(),
        e.width_ratio.ln().abs(),
        e.envelope_cos,
        e.coherence,
        e.m1_ratio,
        e.mass_err_ppm.abs(),
        e.n_points as f64,
        (amplitude_snr(e) + 1.0).log10(),
        (e.quantity_se / e.quantity).ln(),
        e.quantity.log2(),
    ];
    debug_assert_eq!(x.len(), WRONG_FEATURES);
    let x: Vec<f32> = x.iter().map(|&v| v as f32).collect();
    x.iter().all(|v| v.is_finite()).then_some(x)
}

pub(crate) fn running_mean_q(pool: &mut [(f64, usize)], mut set: impl FnMut(usize, f64)) {
    pool.sort_unstable_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
    let mut sum = 0.0;
    let mut start = 0;
    while start < pool.len() {
        let p = pool[start].0;
        let end = start + pool[start..].partition_point(|r| r.0.total_cmp(&p).is_eq());
        sum += pool[start..end].iter().map(|r| r.0).sum::<f64>();
        let q = sum / end as f64;
        for &(_, i) in &pool[start..end] {
            set(i, q);
        }
        start = end;
    }
}

pub(crate) fn transfer_quant_q(transfer_q: Option<f64>, measurement_q: f64) -> f64 {
    transfer_q.map_or(measurement_q, |t| t.max(measurement_q))
}

fn same_fit(a: &Evidence, b: &Evidence) -> bool {
    [
        (a.t0, b.t0),
        (a.w, b.w),
        (a.amplitude, b.amplitude),
        (a.llr, b.llr),
        (a.quantity, b.quantity),
    ]
    .iter()
    .all(|(x, y)| x.to_bits() == y.to_bits())
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Mismatch {
    pub(crate) apex: bool,
    pub(crate) quantity: Option<bool>,
}

impl Mismatch {
    pub(crate) fn wrong(self) -> bool {
        self.apex || self.quantity == Some(true)
    }
}

pub(crate) fn mismatch(pseudo: &Evidence, id: &Evidence) -> Option<Mismatch> {
    if same_fit(pseudo, id) {
        return None;
    }
    let quantity = |e: &Evidence| Some(e.quantity).filter(|q| *q > 0.0 && q.is_finite());
    (pseudo.t0.is_finite() && id.t0.is_finite() && id.w.is_finite()).then(|| Mismatch {
        apex: (pseudo.t0 - id.t0).abs() > id.w,
        quantity: quantity(pseudo)
            .zip(quantity(id))
            .map(|(p, i)| (p / i).log2().abs() > 1.0),
    })
}

pub(crate) fn wrong_label(pseudo: &Evidence, id: &Evidence) -> Option<bool> {
    let m = mismatch(pseudo, id)?;
    m.quantity.map(|_| m.wrong())
}

pub(crate) struct WrongCalibration {
    pub(super) lo: Vec<f32>,
    pub(super) p: Vec<f64>,
}

impl WrongCalibration {
    pub(crate) fn fit(mut pairs: Vec<(f32, bool)>) -> Option<Self> {
        if pairs.is_empty() {
            return None;
        }
        pairs.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));
        let mut pts: Vec<(f32, f64, f64)> = Vec::new();
        for &(s, w) in pairs.iter().rev() {
            match pts.last_mut() {
                Some(l) if l.0 == s => {
                    l.1 += f64::from(u8::from(w));
                    l.2 += 1.0;
                }
                _ => pts.push((s, f64::from(u8::from(w)), 1.0)),
            }
        }
        let mut blocks: Vec<(f32, f64, f64)> = Vec::new();
        for pt in pts {
            blocks.push(pt);
            while blocks.len() > 1 {
                let n = blocks.len();
                let (a, b) = (blocks[n - 2], blocks[n - 1]);
                if a.1 / a.2 <= b.1 / b.2 {
                    break;
                }
                blocks.truncate(n - 2);
                blocks.push((b.0, a.1 + b.1, a.2 + b.2));
            }
        }
        blocks.reverse();
        Some(Self {
            lo: blocks.iter().map(|b| b.0).collect(),
            p: blocks.iter().map(|b| b.1 / b.2).collect(),
        })
    }

    pub(crate) fn at(&self, s: f32) -> f64 {
        let k = self.lo.partition_point(|&b| b <= s);
        self.p[k.saturating_sub(1)]
    }
}

pub(crate) fn pseudo_example(
    r: &ExtractRow,
    ctx: &BlockContext<'_>,
    idx: &EntryIndex<'_>,
) -> Option<(bool, Vec<f32>)> {
    if r.role != Role::PseudoTransfer
        || r.kind != Hypothesis::Target
        || idx.is_identity(r.precursor)
    {
        return None;
    }
    let e = r.evidence.as_ref()?;
    let id = ctx.identified.get(&(r.precursor, r.run))?;
    Some((wrong_label(e, id)?, wrong_features(e)?))
}

pub(super) fn wrong_peak_probabilities(
    rows: &[ExtractRow],
    idx: &EntryIndex<'_>,
    w_models: &FoldModels,
    wrong_cal: Option<&WrongCalibration>,
) -> Vec<f64> {
    let mut p_wrong: Vec<f64> = vec![1.0; rows.len()];
    let at: Vec<usize> = (0..rows.len())
        .filter(|&i| {
            rows[i].role == Role::Transfer
                && matches!(idx.table(&rows[i]), Some((Table::Transfer, TARGET)))
        })
        .collect();
    let xs: Vec<Option<Vec<f32>>> = at
        .iter()
        .map(|&i| rows[i].evidence.as_ref().and_then(wrong_features))
        .collect();
    let ok: Vec<usize> = (0..at.len()).filter(|&k| xs[k].is_some()).collect();
    let precursors: Vec<u32> = ok.iter().map(|&k| rows[at[k]].precursor).collect();
    let x_refs: Vec<&[f32]> = ok
        .iter()
        .map(|&k| xs[k].as_deref().expect("filtered on features"))
        .collect();
    let mut margins = Vec::new();
    w_models.score(&precursors, &x_refs, &mut margins);
    for (&k, m) in ok.iter().zip(margins) {
        if let (Some(m), Some(c)) = (m, wrong_cal) {
            p_wrong[at[k]] = c.at(m);
        }
    }
    p_wrong
}
