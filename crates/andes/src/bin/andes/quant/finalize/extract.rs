//! Reading the per-run extract tables block by block (from runs in memory or from files), with
//! the bounds and context each block needs.

use std::collections::{HashMap, HashSet};
use std::path::PathBuf;

use output::extract::ExtractTableReader;
use quant::align::RunMap;
use quant::decoys::Hypothesis;
use quant::evidence::Evidence;
use rayon::prelude::*;

use super::READ_BATCH_ROWS;
use crate::quant::mbr::{donor_fit, DonorFit, ExtractRow, Role};

pub(crate) trait ExtractSource: Sync {
    fn chunks(&self) -> usize;
    fn name(&self, chunk: usize) -> String;
    fn for_each_batch(
        &self,
        chunk: usize,
        f: &mut dyn FnMut(Vec<ExtractRow>) -> Result<(), String>,
    ) -> Result<(), String>;

    fn for_each_fit_batch(
        &self,
        chunk: usize,
        f: &mut dyn FnMut(Vec<ExtractRow>) -> Result<(), String>,
    ) -> Result<(), String> {
        self.for_each_batch(chunk, &mut |rows| {
            f(rows
                .into_iter()
                .filter(|r| (r.role, r.kind) == (Role::Identified, Hypothesis::Target))
                .collect())
        })
    }
}

#[cfg(test)]
pub(crate) struct ExtractRuns<'a> {
    pub runs: &'a [Vec<ExtractRow>],
    pub batch_rows: usize,
}

#[cfg(test)]
impl ExtractSource for ExtractRuns<'_> {
    fn chunks(&self) -> usize {
        self.runs.len()
    }

    fn name(&self, chunk: usize) -> String {
        format!("extract rows {chunk}")
    }

    fn for_each_batch(
        &self,
        chunk: usize,
        f: &mut dyn FnMut(Vec<ExtractRow>) -> Result<(), String>,
    ) -> Result<(), String> {
        self.runs[chunk]
            .chunks(self.batch_rows.max(1))
            .try_for_each(|b| f(b.to_vec()))
    }
}

pub(crate) struct ExtractFiles(pub Vec<PathBuf>);

impl ExtractSource for ExtractFiles {
    fn chunks(&self) -> usize {
        self.0.len()
    }

    fn name(&self, chunk: usize) -> String {
        self.0[chunk].display().to_string()
    }

    fn for_each_batch(
        &self,
        chunk: usize,
        f: &mut dyn FnMut(Vec<ExtractRow>) -> Result<(), String>,
    ) -> Result<(), String> {
        for b in ExtractTableReader::open(&self.0[chunk], READ_BATCH_ROWS)? {
            f(b?)?;
        }
        Ok(())
    }

    fn for_each_fit_batch(
        &self,
        chunk: usize,
        f: &mut dyn FnMut(Vec<ExtractRow>) -> Result<(), String>,
    ) -> Result<(), String> {
        for b in ExtractTableReader::open_identified_fits(&self.0[chunk], READ_BATCH_ROWS)? {
            f(b?)?;
        }
        Ok(())
    }
}

#[derive(Debug, Clone, Copy)]
pub(super) struct Bounds {
    pub(super) runs: usize,
    pub(super) precursors: usize,
}

pub(super) fn for_each_block<S>(
    source: &dyn ExtractSource,
    fits_only: bool,
    bounds: Bounds,
    init: impl Fn(usize) -> S + Sync,
    f: impl Fn(&mut S, usize, u64, &[ExtractRow]) -> Result<(), String> + Sync,
    done: impl Fn(usize, S) -> Result<(), String> + Sync,
) -> Result<(), String> {
    (0..source.chunks()).into_par_iter().try_for_each(|chunk| {
        let mut state = init(chunk);
        let mut closed: HashSet<u32> = HashSet::new();
        let mut emit = |state: &mut S, rows: &[ExtractRow], first: u64| {
            if let Some(r) = rows
                .iter()
                .find(|r| r.run as usize >= bounds.runs || r.precursor as usize >= bounds.precursors)
            {
                return Err(format!(
                    "extract row of precursor {} in run {} is outside the {} precursors and {} runs",
                    r.precursor, r.run, bounds.precursors, bounds.runs
                ));
            }
            for g in rows.chunk_by(|a, b| a.precursor == b.precursor) {
                if !closed.insert(g[0].precursor) {
                    return Err(format!(
                        "{}: the rows of precursor {} are not contiguous",
                        source.name(chunk),
                        g[0].precursor
                    ));
                }
            }
            f(state, chunk, first, rows)
        };
        let mut buf: Vec<ExtractRow> = Vec::new();
        let mut first = 0u64;
        let mut each = |batch: Vec<ExtractRow>| {
            buf.extend(batch);
            let Some(last) = buf.last().map(|r| r.precursor) else {
                return Ok(());
            };
            let cut = buf.len() - buf.iter().rev().take_while(|r| r.precursor == last).count();
            if cut > 0 {
                emit(&mut state, &buf[..cut], first)?;
                first += cut as u64;
                buf.drain(..cut);
            }
            Ok(())
        };
        if fits_only {
            source.for_each_fit_batch(chunk, &mut each)?;
        } else {
            source.for_each_batch(chunk, &mut each)?;
        }
        if !buf.is_empty() {
            emit(&mut state, &buf, first)?;
        }
        done(chunk, state)
    })
}

pub(crate) struct BlockContext<'a> {
    pub(super) identified: HashMap<(u32, u16), &'a Evidence>,
    pub(super) own: HashMap<(u32, u16), DonorFit>,
}

impl<'a> BlockContext<'a> {
    pub(crate) fn new(rows: &'a [ExtractRow], maps: &[Option<RunMap>]) -> Result<Self, String> {
        let mut identified = HashMap::new();
        let mut own = HashMap::new();
        for r in rows {
            if let (Role::Identified, Hypothesis::Target, Some(e)) = (r.role, r.kind, &r.evidence) {
                if identified.insert((r.precursor, r.run), e).is_some() {
                    return Err(format!(
                        "two identified rows of precursor {} in run {}",
                        r.precursor, r.run
                    ));
                }
                if let Some(f) = donor_fit(r, maps) {
                    own.insert((r.precursor, r.run), f);
                }
            }
        }
        Ok(Self { identified, own })
    }
}
