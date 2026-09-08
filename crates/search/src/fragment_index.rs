//! Per-chunk fragment-ion index for PTM-rich out-of-core searches (issue #76).
//!
//! The enumeration path materialises every peptidoform of every base record in
//! a spectrum's precursor windows and runs the full scorer on each. On a
//! phospho search that is tens of thousands of forms per spectrum, and the
//! generation alone was most of the CPU. This index inverts the work: the
//! chunk's peptidoforms are enumerated ONCE, each form's singly-charged b and y
//! ions are binned, and a spectrum's peaks vote for the forms whose ions they
//! hit. Only forms with at least `min_matched` votes, and at most `top_k` of
//! them, are materialised and handed to the unchanged scoring path.
//!
//! What changes when it is on: the candidate SET entering the scorer. Forms
//! with fewer than `min_matched` fragment hits are never scored, and the
//! enumeration multiplicity copies (the PIN `Proteins` repeats) are not
//! reproduced. Everything downstream of candidate assembly is untouched.
//! With it off every path is byte-identical.

use rayon::prelude::*;
use rustc_hash::FxHashMap;

use crate::candidate_gen::{
    base_record_key, expand_base_record, for_each_record_form_masses, BaseRecordKey, Candidate,
};
use crate::candidate_index::IndexRecord;
use crate::precursor_cal::adjusted_observed_neutral_mass;
use crate::search_index::SearchIndex;
use crate::search_params::SearchParams;
use model::mass::{H2O, ISOTOPE, PROTON};
use model::spectrum::Spectrum;
use model::tolerance::Tolerance;

pub struct ChunkFragmentIndex {
    /// The chunk's distinct base records, in `enumerate_candidates` order.
    records: Vec<IndexRecord>,
    form_record: Vec<u32>,
    form_k: Vec<u16>,
    form_mass: Vec<f64>,
    bin_width: f64,
    /// CSR over fragment bins: entries of bin `b` are
    /// `entries[bin_start[b]..bin_start[b + 1]]`, sorted by the form's mass.
    bin_start: Vec<u64>,
    /// (form id, ion m/z).
    entries: Vec<(u32, f32)>,
}

/// Singly-charged b (prefix) and y (suffix) ion m/z of one peptidoform from its
/// per-residue masses (modification deltas folded in), appended to `out`.
fn by_ions_from_masses(masses: &[f64], out: &mut Vec<f32>) {
    let n = masses.len();
    if n < 2 {
        return;
    }
    let mut prefix = 0.0;
    for m in &masses[..n - 1] {
        prefix += m;
        out.push((prefix + PROTON) as f32);
    }
    let mut suffix = 0.0;
    for m in masses[1..].iter().rev() {
        suffix += m;
        out.push((suffix + H2O + PROTON) as f32);
    }
}

#[inline]
fn neutral_mass(masses: &[f64]) -> f64 {
    masses.iter().sum::<f64>() + H2O
}

impl ChunkFragmentIndex {
    /// Enumerate every peptidoform of `records` once and bin its ions. The bin
    /// width is the fragment tolerance at the top of the fragment m/z range
    /// (2,500), so a peak checking its own bin and both neighbours finds every
    /// ion within tolerance whatever its m/z.
    /// Only forms whose neutral mass lies in `[mass_lo, mass_hi]` (the union
    /// of the chunk's precursor windows) are indexed: a record reachable
    /// through one modification offset carries dozens of forms at other
    /// masses, and indexing them all made one chunk's index exceed 2^32 ions.
    ///
    /// Two passes over the expansion so nothing per form is held beyond its
    /// packed entries: pass 1 counts forms per record and ions per bin, pass 2
    /// re-expands and writes entries straight into the packed table.
    pub fn build(
        records: Vec<IndexRecord>,
        db: &SearchIndex,
        params: &SearchParams,
        fragment_tol: Tolerance,
        mass_lo: f64,
        mass_hi: f64,
    ) -> Self {
        use std::sync::atomic::{AtomicU32, AtomicU64, Ordering};
        let bin_width = fragment_tol.as_da(2500.0).max(0.001);
        // Pass 1: forms per record, ions per bin. The lean expansion visits
        // forms as residue-mass slices; nothing per form is allocated.
        let (forms_per_record, counts): (Vec<u32>, Vec<u64>) = {
            let per: Vec<(u32, Vec<(u32, u32)>)> = records
                .par_iter()
                .map(|rec| {
                    let mut n = 0u32;
                    let mut bins: Vec<(u32, u32)> = Vec::new();
                    let mut ions: Vec<f32> = Vec::new();
                    for_each_record_form_masses(db, params, rec, |_, masses| {
                        let m = neutral_mass(masses);
                        if m < mass_lo || m > mass_hi {
                            return;
                        }
                        n += 1;
                        ions.clear();
                        by_ions_from_masses(masses, &mut ions);
                        for &mz in &ions {
                            bins.push(((mz as f64 / bin_width) as u32, 1));
                        }
                    });
                    bins.sort_unstable();
                    bins.dedup_by(|a, b| {
                        if a.0 == b.0 {
                            b.1 += a.1;
                            true
                        } else {
                            false
                        }
                    });
                    (n, bins)
                })
                .collect();
            let max_bin = per
                .iter()
                .flat_map(|(_, b)| b.iter().map(|x| x.0 as usize))
                .max()
                .unwrap_or(0);
            let mut counts = vec![0u64; max_bin + 3];
            for (_, bins) in &per {
                for &(b, c) in bins {
                    counts[b as usize + 1] += c as u64;
                }
            }
            (per.into_iter().map(|(n, _)| n).collect(), counts)
        };
        let n_bins = counts.len() - 1;
        let mut bin_start = counts;
        for b in 0..n_bins {
            bin_start[b + 1] += bin_start[b];
        }
        let n_entries = bin_start[n_bins] as usize;
        let n_forms: usize = forms_per_record.iter().map(|&n| n as usize).sum();
        let mut form_base: Vec<u32> = Vec::with_capacity(records.len());
        let mut acc = 0u32;
        for &n in &forms_per_record {
            form_base.push(acc);
            acc += n;
        }
        eprintln!(
            "fragment-index: chunk window {:.1}-{:.1} Da, {} records, {} forms, {} ion entries, {} bins",
            mass_lo,
            mass_hi,
            records.len(),
            n_forms,
            n_entries,
            n_bins
        );
        // Pass 2: fill the packed tables; per-bin write cursors are atomics so
        // records fill in parallel. Bins are sorted by form mass afterwards, so
        // the result does not depend on scheduling.
        let form_record: Vec<AtomicU64> = (0..n_forms).map(|_| AtomicU64::new(0)).collect();
        let form_mass_atomic: Vec<AtomicU64> = (0..n_forms).map(|_| AtomicU64::new(0)).collect();
        let fill: Vec<AtomicU64> = bin_start.iter().map(|&s| AtomicU64::new(s)).collect();
        let entries: Vec<(AtomicU32, AtomicU32)> = (0..n_entries)
            .map(|_| (AtomicU32::new(0), AtomicU32::new(0)))
            .collect();
        records.par_iter().enumerate().for_each(|(r, rec)| {
            let mut form_id = form_base[r];
            let mut ions: Vec<f32> = Vec::new();
            for_each_record_form_masses(db, params, rec, |k, masses| {
                let m = neutral_mass(masses);
                if m < mass_lo || m > mass_hi {
                    return;
                }
                form_record[form_id as usize]
                    .store(((r as u64) << 16) | (k as u64 & 0xffff), Ordering::Relaxed);
                form_mass_atomic[form_id as usize].store(m.to_bits(), Ordering::Relaxed);
                ions.clear();
                by_ions_from_masses(masses, &mut ions);
                for &mz in &ions {
                    let b = (mz as f64 / bin_width) as usize;
                    let pos = fill[b].fetch_add(1, Ordering::Relaxed) as usize;
                    entries[pos].0.store(form_id, Ordering::Relaxed);
                    entries[pos].1.store(mz.to_bits(), Ordering::Relaxed);
                }
                form_id += 1;
            });
        });
        let form_record_k: Vec<u64> = form_record.into_iter().map(|a| a.into_inner()).collect();
        let form_record: Vec<u32> = form_record_k.iter().map(|v| (v >> 16) as u32).collect();
        let form_k: Vec<u16> = form_record_k.iter().map(|v| (v & 0xffff) as u16).collect();
        let form_mass: Vec<f64> = form_mass_atomic
            .into_iter()
            .map(|a| f64::from_bits(a.into_inner()))
            .collect();
        let mut entries: Vec<(u32, f32)> = entries
            .into_iter()
            .map(|(a, b)| (a.into_inner(), f32::from_bits(b.into_inner())))
            .collect();
        for b in 0..n_bins {
            let (lo, hi) = (bin_start[b] as usize, bin_start[b + 1] as usize);
            entries[lo..hi].sort_unstable_by(|x, y| {
                form_mass[x.0 as usize]
                    .partial_cmp(&form_mass[y.0 as usize])
                    .unwrap()
                    .then(x.0.cmp(&y.0))
            });
        }
        Self {
            records,
            form_record,
            form_k,
            form_mass,
            bin_width,
            bin_start,
            entries,
        }
    }

    pub fn n_forms(&self) -> usize {
        self.form_mass.len()
    }

    /// Forms with at least `min_matched` singly-charged b/y ions within
    /// `fragment_tol` of a peak (ppm tolerances are evaluated per peak),
    /// restricted to the spectrum's precursor windows over every charge in
    /// `charges` and every isotope offset in `params`. At most `top_k`, best
    /// votes first.
    pub fn query(
        &self,
        spec: &Spectrum,
        charges: &[u8],
        params: &SearchParams,
        fragment_tol: Tolerance,
        top_k: usize,
        min_matched: u16,
    ) -> Vec<(u32, u16)> {
        // One over-inclusive precursor mass interval; the scoring loop applies
        // the exact per-offset test afterwards.
        let shift_ppm = params.precursor_mass_shift_ppm;
        let mut lo = f64::MAX;
        let mut hi = f64::MIN;
        for &z in charges {
            let zf = z as f64;
            let obs =
                adjusted_observed_neutral_mass(spec.precursor_mz * zf - zf * PROTON, shift_ppm);
            for o in params.isotope_error_range.clone() {
                let c = obs - (o as f64) * ISOTOPE;
                lo = lo.min(c - params.precursor_tolerance.left.as_da(c));
                hi = hi.max(c + params.precursor_tolerance.right.as_da(c));
            }
        }
        if lo > hi {
            return Vec::new();
        }
        let mut votes: FxHashMap<u32, u16> = FxHashMap::default();
        let n_bins = self.bin_start.len() - 1;
        for &(mz, _) in &spec.peaks {
            let tol_da = fragment_tol.as_da(mz);
            let b = (mz / self.bin_width) as usize;
            for bb in b.saturating_sub(1)..=(b + 1).min(n_bins - 1) {
                let seg =
                    &self.entries[self.bin_start[bb] as usize..self.bin_start[bb + 1] as usize];
                let start = seg.partition_point(|e| self.form_mass[e.0 as usize] < lo);
                for e in &seg[start..] {
                    if self.form_mass[e.0 as usize] > hi {
                        break;
                    }
                    if ((e.1 as f64) - mz).abs() <= tol_da {
                        *votes.entry(e.0).or_insert(0) += 1;
                    }
                }
            }
        }
        let mut sel: Vec<(u32, u16)> = votes
            .into_iter()
            .filter(|&(_, v)| v >= min_matched)
            .collect();
        sel.sort_unstable_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
        sel.truncate(top_k);
        sel
    }

    /// Materialise selected forms as candidates, in (record, k) order. Cached
    /// expansions are used when present; otherwise the record is re-expanded.
    pub fn materialise(
        &self,
        selected: &[(u32, u16)],
        db: &SearchIndex,
        params: &SearchParams,
        cache: Option<&FxHashMap<BaseRecordKey, Vec<Candidate>>>,
        relabel: &dyn Fn(&mut [Candidate]),
    ) -> Vec<Candidate> {
        let mut by_record: Vec<(u32, u16)> = selected
            .iter()
            .map(|&(f, _)| (self.form_record[f as usize], self.form_k[f as usize]))
            .collect();
        by_record.sort_unstable();
        let mut out = Vec::with_capacity(by_record.len());
        let mut i = 0;
        while i < by_record.len() {
            let r = by_record[i].0;
            let mut j = i;
            while j < by_record.len() && by_record[j].0 == r {
                j += 1;
            }
            let rec = &self.records[r as usize];
            let key = base_record_key(rec);
            let owned: Option<Vec<Candidate>> = match cache.and_then(|c| c.get(&key)) {
                Some(_) => None,
                None => {
                    let mut v = expand_base_record(db, params, rec);
                    relabel(&mut v);
                    Some(v)
                }
            };
            let forms: &[Candidate] = match &owned {
                Some(v) => v.as_slice(),
                None => cache.unwrap().get(&key).unwrap().as_slice(),
            };
            for &(_, k) in &by_record[i..j] {
                if let Some(c) = forms.get(k as usize) {
                    out.push(c.clone());
                }
            }
            i = j;
        }
        out
    }
}
