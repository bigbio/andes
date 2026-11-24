//! Per-precursor peak-shape summaries: the median window, donor fits, identified SNR and the
//! quantity rule, and the shape reference the transfer features compare against.

use std::collections::HashMap;

use quant::align::RunMap;
use quant::decoys::Hypothesis;

use crate::quant::mbr::{
    amplitude_snr, cosine, iso_profile, local_slope, median, DonorFit, ExtractRow, Role,
    ShapeReference, NO_SHAPE,
};

#[derive(Debug, Clone, Copy)]
pub(super) struct MedianWindow {
    pub(super) n: u32,
    pub(super) first: u32,
    pub(super) v: [f64; 3],
}

impl MedianWindow {
    pub(super) fn new(sorted: &[f64]) -> Self {
        let n = sorted.len();
        let first = (n / 2).saturating_sub(1);
        let mut v = [f64::NAN; 3];
        for (slot, x) in v.iter_mut().zip(sorted.iter().skip(first)) {
            *slot = *x;
        }
        Self {
            n: n as u32,
            first: first as u32,
            v,
        }
    }

    pub(super) fn at(&self, i: usize) -> f64 {
        self.v[i - self.first as usize]
    }

    pub(super) fn median(&self) -> Option<f64> {
        let n = self.n as usize;
        match n {
            0 => None,
            _ if n % 2 == 1 => Some(self.at(n / 2)),
            _ => Some(0.5 * (self.at(n / 2 - 1) + self.at(n / 2))),
        }
    }

    pub(super) fn median_without(&self, x: f64) -> Option<f64> {
        let m = (self.n as usize).saturating_sub(1);
        if m == 0 {
            return None;
        }
        let left = |j: usize| {
            let v = self.at(j);
            if v.total_cmp(&x).is_lt() {
                v
            } else {
                self.at(j + 1)
            }
        };
        Some(if m % 2 == 1 {
            left(m / 2)
        } else {
            0.5 * (left(m / 2 - 1) + left(m / 2))
        })
    }

    fn median_of(&self, x: Option<f64>) -> Option<f64> {
        match x {
            Some(x) => self.median_without(x),
            None => self.median(),
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub(super) struct ShapeSummary {
    width: MedianWindow,
    tailing: MedianWindow,
    iso: [MedianWindow; 4],
}

#[derive(Debug, Clone, Copy, Default)]
pub(super) struct PrecursorSummary {
    pub(super) relative_width: Option<f64>,
    pub(super) shape: Option<ShapeSummary>,
    pub(super) snr: Option<f64>,
    pub(super) trace: bool,
}

pub(super) type Summaries = HashMap<u32, PrecursorSummary>;

pub(super) type SummaryItem = (u32, Option<f64>, Option<DonorFit>, Option<f64>);

pub(super) fn summarise(mut items: Vec<SummaryItem>) -> Summaries {
    items.sort_unstable_by_key(|it| it.0);
    let mut out = Summaries::new();
    for g in items.chunk_by(|a, b| a.0 == b.0) {
        let ratios: Vec<f64> = g.iter().filter_map(|it| it.1).collect();
        let snrs: Vec<f64> = g.iter().filter_map(|it| it.3).collect();
        let fits: Vec<&DonorFit> = g.iter().filter_map(|it| it.2.as_ref()).collect();
        let window = |value: &dyn Fn(&DonorFit) -> f64| {
            let mut v: Vec<f64> = fits.iter().map(|f| value(f)).collect();
            v.sort_by(f64::total_cmp);
            MedianWindow::new(&v)
        };
        let shape = (!fits.is_empty()).then(|| ShapeSummary {
            width: window(&|f| f.1),
            tailing: window(&|f| f.2),
            iso: [0, 1, 2, 3].map(|i| window(&|f: &DonorFit| f.3[i])),
        });
        out.insert(
            g[0].0,
            PrecursorSummary {
                relative_width: (!ratios.is_empty()).then(|| median(ratios)),
                shape,
                snr: (!snrs.is_empty()).then(|| median(snrs)),
                trace: false,
            },
        );
    }
    out
}

pub(super) fn identified_snr(r: &ExtractRow) -> Option<f64> {
    if (r.role, r.kind) != (Role::Identified, Hypothesis::Target) {
        return None;
    }
    let e = r.evidence.as_ref()?;
    Some(amplitude_snr(e)).filter(|v| *v > 0.0 && v.is_finite())
}

pub(super) fn quantity_rule(summaries: &Summaries, precursor: u32) -> (Option<f64>, bool) {
    let s = summaries.get(&precursor);
    (s.and_then(|s| s.relative_width), s.is_some_and(|s| s.trace))
}

pub(super) fn mark_trace_precursors(summaries: &mut Summaries, fraction: f64) {
    let mut snrs: Vec<f64> = summaries.values().filter_map(|s| s.snr).collect();
    if snrs.is_empty() {
        log::info!(target: "finalize", "trace area for 0 precursors (no identified fit with a signal-to-noise)");
        return;
    }
    snrs.sort_unstable_by(f64::total_cmp);
    let n = snrs.len();
    let k = ((fraction * n as f64).ceil() as usize).clamp(1, n);
    let cut = snrs[n - k];
    let mut marked = 0usize;
    for s in summaries.values_mut() {
        s.trace = s.snr.is_some_and(|v| v >= cut);
        marked += usize::from(s.trace);
    }
    let at = |q: f64| snrs[((q * (n - 1) as f64).round() as usize).min(n - 1)];
    log::info!(
        target: "finalize",
        "trace area for {marked} of {n} precursors with an identified fit (share {fraction}): median amplitude signal-to-noise >= {cut:.1}; distribution p50 {:.1}, p90 {:.1}, p99 {:.1}",
        at(0.5),
        at(0.9),
        at(0.99)
    );
}

pub(super) struct BlockShapes<'a> {
    pub(super) summaries: &'a Summaries,
    pub(super) own: &'a HashMap<(u32, u16), DonorFit>,
}

impl ShapeReference for BlockShapes<'_> {
    fn features(&self, row: &ExtractRow, maps: &[Option<RunMap>]) -> [f64; 3] {
        let Some(e) = row.evidence.as_ref().filter(|e| e.amplitude > 0.0) else {
            return NO_SHAPE;
        };
        let shape = self
            .summaries
            .get(&row.precursor)
            .and_then(|s| s.shape.as_ref());
        let own = self.own.get(&(row.precursor, row.run));
        let refs = shape
            .map_or(0, |s| s.width.n as usize)
            .saturating_sub(usize::from(own.is_some()));
        let mine = local_slope(maps, row.run, e.t0).zip(iso_profile(&e.iso_heights));
        let (Some((slope, iso)), Some(s)) = (mine, shape.filter(|_| refs > 0)) else {
            return NO_SHAPE;
        };
        let reference = |w: &MedianWindow, x: Option<f64>| {
            w.median_of(x)
                .expect("a precursor with references has a median")
        };
        let w_ref = reference(&s.width, own.map(|f| f.1));
        let tail_ref = reference(&s.tailing, own.map(|f| f.2));
        let mut iso_ref = [0.0; 4];
        for (i, v) in iso_ref.iter_mut().enumerate() {
            *v = reference(&s.iso[i], own.map(|f| f.3[i]));
        }
        [
            ((e.w / slope) / w_ref).ln(),
            e.tailing - tail_ref,
            cosine(&iso, &iso_ref),
        ]
    }
}
