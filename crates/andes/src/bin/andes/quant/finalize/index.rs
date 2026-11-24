//! The library entry index: neighbour runs, score strata, the feature table of a transfer and
//! the walk over scored rows.

use quant::align::RunMap;
use quant::decoys::Hypothesis;
use quant::library::{Library, LibraryEntry};
use rayon::prelude::*;

use super::SCORE_CHUNK;
use crate::quant::mbr::{
    family, identified_features, transfer_features, ExtractRow, Role, ShapeReference,
};
use crate::transfers::{pack_family, stratum_of, FoldModels, IDENTITY, MAX_STRATA};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Table {
    Transfer,
    Identified,
    Pseudo,
}

pub(crate) const PRIOR_NEIGHBOURS: usize = 20;
const PRIOR_SAMPLE: usize = 20_000;

pub(crate) const NEIGHBOUR_RUNS: usize = 50;

pub(crate) struct EntryIndex<'a> {
    pub(super) entry_of: Vec<Option<&'a LibraryEntry>>,
    pub(super) neighbours: Vec<Vec<u16>>,
    edges: Vec<f64>,
}

impl<'a> EntryIndex<'a> {
    pub(crate) fn new(lib: &'a Library, runs: usize) -> Self {
        let mut entry_of: Vec<Option<&LibraryEntry>> = vec![None; lib.meta.len()];
        for e in &lib.entries {
            if let Some(slot) = entry_of.get_mut(e.precursor as usize) {
                *slot = Some(e);
            }
        }
        let neighbours = run_neighbours(lib, runs, NEIGHBOUR_RUNS);
        let edges = strata_edges(&entry_of, &neighbours, runs, PRIOR_NEIGHBOURS);
        Self {
            entry_of,
            neighbours,
            edges,
        }
    }

    pub(crate) fn strata(&self) -> usize {
        if self.by_run() {
            self.neighbours.len()
        } else {
            self.edges.len() + 1
        }
    }

    fn by_run(&self) -> bool {
        (2..=MAX_STRATA).contains(&self.neighbours.len())
    }

    pub(crate) fn presence_prior(&self, r: &ExtractRow) -> f64 {
        let (Some(Some(e)), Some(nn)) = (
            self.entry_of.get(r.precursor as usize),
            self.neighbours.get(r.run as usize),
        ) else {
            return 0.0;
        };
        let nn = &nn[..nn.len().min(PRIOR_NEIGHBOURS)];
        if nn.is_empty() {
            return 0.0;
        }
        let hits = nn
            .iter()
            .filter(|n| e.donors.binary_search(n).is_ok())
            .count();
        hits as f64 / nn.len() as f64
    }

    pub(crate) fn packed(&self, r: &ExtractRow, t: Table, family: u8) -> u32 {
        match t {
            Table::Transfer | Table::Pseudo if self.by_run() => pack_family(family, r.run),
            Table::Transfer | Table::Pseudo => {
                pack_family(family, stratum_of(self.presence_prior(r), &self.edges))
            }
            Table::Identified => family.into(),
        }
    }

    pub(crate) fn neighbour_rate(&self, r: &ExtractRow) -> f64 {
        let (Some(Some(e)), Some(nn)) = (
            self.entry_of.get(r.precursor as usize),
            self.neighbours.get(r.run as usize),
        ) else {
            return 0.0;
        };
        if nn.is_empty() {
            return 0.0;
        }
        let hits = nn
            .iter()
            .filter(|n| e.donors.binary_search(n).is_ok())
            .count();
        hits as f64 / nn.len() as f64
    }

    pub(crate) fn is_identity(&self, p: u32) -> bool {
        self.entry_of[p as usize].is_some_and(|e| e.identity_decoy)
    }

    pub(crate) fn donor_pep(&self, p: u32) -> f64 {
        self.entry_of
            .get(p as usize)
            .copied()
            .flatten()
            .map_or(1.0, |e| {
                if e.donor_psm.pep.is_finite() {
                    e.donor_psm.pep
                } else {
                    1.0
                }
            })
    }

    pub(crate) fn table(&self, r: &ExtractRow) -> Option<(Table, u8)> {
        let identity = self.is_identity(r.precursor);
        let t = match r.role {
            Role::Identified if identity => return None,
            Role::Identified => Table::Identified,
            Role::Transfer if identity && r.kind != Hypothesis::Target => return None,
            Role::Transfer => Table::Transfer,
            Role::PseudoTransfer => Table::Pseudo,
        };
        Some((t, if identity { IDENTITY } else { family(r.kind) }))
    }
}

pub(crate) fn strata_edges(
    entry_of: &[Option<&LibraryEntry>],
    neighbours: &[Vec<u16>],
    runs: usize,
    k: usize,
) -> Vec<f64> {
    let entries: Vec<&LibraryEntry> = entry_of
        .iter()
        .flatten()
        .copied()
        .filter(|e| !e.identity_decoy)
        .collect();
    let m = neighbours.iter().map(Vec::len).min().unwrap_or(0).min(k);
    if m == 0 {
        return Vec::new();
    }
    let stride = entries.len().div_ceil(PRIOR_SAMPLE).max(1);
    let hist = entries
        .iter()
        .step_by(stride)
        .par_bridge()
        .fold(
            || vec![0u64; m + 1],
            |mut h, e| {
                for (r, nn) in neighbours.iter().enumerate().take(runs) {
                    if e.donors.binary_search(&(r as u16)).is_ok() || nn.is_empty() {
                        continue;
                    }
                    let nn = &nn[..m];
                    let hits = nn
                        .iter()
                        .filter(|n| e.donors.binary_search(n).is_ok())
                        .count();
                    h[hits] += 1;
                }
                h
            },
        )
        .reduce(
            || vec![0u64; m + 1],
            |mut a, b| {
                a.iter_mut().zip(b).for_each(|(x, y)| *x += y);
                a
            },
        );
    let edges = edges_of_hist(&hist);
    log::debug!(target: "mbr", "presence prior strata edges {edges:?}");
    edges
}

pub(super) fn edges_of_hist(hist: &[u64]) -> Vec<f64> {
    let m = hist.len().saturating_sub(1);
    let total: u64 = hist.iter().sum();
    let mut edges = Vec::new();
    if total > 0 && m > 0 {
        let share = total as f64 / MAX_STRATA as f64;
        let mut open = 0.0;
        for (i, &c) in hist.iter().enumerate().take(m) {
            open += c as f64;
            if open >= share && edges.len() + 1 < MAX_STRATA {
                edges.push(i as f64 / m as f64);
                open = 0.0;
            }
        }
        if open + (hist[m] as f64) < share && !edges.is_empty() {
            edges.pop();
        }
    }
    edges
}

pub(crate) fn run_neighbours(lib: &Library, runs: usize, k: usize) -> Vec<Vec<u16>> {
    let words = lib.entries.len().div_ceil(64);
    let mut bits = vec![vec![0u64; words]; runs];
    for (i, e) in lib.entries.iter().enumerate() {
        if e.identity_decoy {
            continue;
        }
        for &r in &e.donors {
            if let Some(b) = bits.get_mut(r as usize) {
                b[i / 64] |= 1u64 << (i % 64);
            }
        }
    }
    let ones = |b: &[u64]| -> u64 { b.iter().map(|w| u64::from(w.count_ones())).sum() };
    let count: Vec<u64> = bits.iter().map(|b| ones(b)).collect();
    (0..runs)
        .into_par_iter()
        .map(|a| {
            let mut sim: Vec<(f64, u16)> = (0..runs)
                .filter(|&b| b != a)
                .map(|b| {
                    let inter: u64 = bits[a]
                        .iter()
                        .zip(&bits[b])
                        .map(|(x, y)| u64::from((x & y).count_ones()))
                        .sum();
                    let union = count[a] + count[b] - inter;
                    let j = if union == 0 {
                        0.0
                    } else {
                        inter as f64 / union as f64
                    };
                    (j, b as u16)
                })
                .collect();
            sim.sort_by(|x, y| y.0.total_cmp(&x.0).then(x.1.cmp(&y.1)));
            sim.truncate(k);
            sim.into_iter().map(|(_, b)| b).collect()
        })
        .collect()
}

pub(crate) fn table_features(
    rows: &[ExtractRow],
    idx: &EntryIndex<'_>,
    shapes: &(impl ShapeReference + Sync),
    maps: &[Option<RunMap>],
    want: impl Fn(Table) -> bool + Sync,
) -> Vec<(usize, Table, u8, Option<Vec<f32>>)> {
    rows.par_iter()
        .enumerate()
        .filter_map(|(i, r)| {
            let (t, fam) = idx.table(r).filter(|(t, _)| want(*t))?;
            let x = match t {
                Table::Identified => identified_features(r),
                Table::Transfer | Table::Pseudo => transfer_features(
                    r,
                    shapes,
                    maps,
                    idx.neighbour_rate(r),
                    idx.donor_pep(r.precursor),
                ),
            };
            Some((i, t, fam, x))
        })
        .collect()
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn for_each_scored(
    rows: &[ExtractRow],
    idx: &EntryIndex<'_>,
    shapes: &(impl ShapeReference + Sync),
    maps: &[Option<RunMap>],
    transfer: &FoldModels,
    identified: &FoldModels,
    want: impl Fn(Table) -> bool + Sync + Copy,
    mut f: impl FnMut(usize, Table, u32, Option<f32>),
) {
    let mut out = Vec::new();
    for (c, chunk) in rows.chunks(SCORE_CHUNK).enumerate() {
        let base = c * SCORE_CHUNK;
        let feats = table_features(chunk, idx, shapes, maps, want);
        let mut scores: Vec<Option<f32>> = vec![None; feats.len()];
        for (models, is_id) in [(transfer, false), (identified, true)] {
            let at: Vec<usize> = (0..feats.len())
                .filter(|&k| feats[k].3.is_some() && ((feats[k].1 == Table::Identified) == is_id))
                .collect();
            if at.is_empty() {
                continue;
            }
            let precursors: Vec<u32> = at.iter().map(|&k| chunk[feats[k].0].precursor).collect();
            let xs: Vec<&[f32]> = at
                .iter()
                .map(|&k| feats[k].3.as_deref().expect("filtered on features"))
                .collect();
            models.score(&precursors, &xs, &mut out);
            for (&k, s) in at.iter().zip(&out) {
                scores[k] = *s;
            }
        }
        for ((i, t, fam, x), s) in feats.iter().zip(scores) {
            if x.is_some() {
                f(base + i, *t, idx.packed(&chunk[*i], *t, *fam), s);
            }
        }
    }
}
