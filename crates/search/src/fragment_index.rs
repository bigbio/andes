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

use crate::candidate_gen::{base_record_key, expand_base_record, BaseRecordKey, Candidate};
use crate::candidate_index::IndexRecord;
use crate::precursor_cal::adjusted_observed_neutral_mass;
use crate::search_index::SearchIndex;
use crate::search_params::SearchParams;
use model::mass::{H2O, ISOTOPE, PROTON};
use model::spectrum::Spectrum;
use model::tolerance::Tolerance;

/// One indexed peptidoform: which base record it came from, its position in
/// that record's `expand_mod_combinations` order, and its neutral mass.
struct FormIons {
    record: u32,
    k: u16,
    mass: f64,
    ions: Vec<f32>,
}

pub struct ChunkFragmentIndex {
    /// The chunk's distinct base records, in `enumerate_candidates` order.
    records: Vec<IndexRecord>,
    form_record: Vec<u32>,
    form_k: Vec<u16>,
    form_mass: Vec<f64>,
    bin_width: f64,
    /// CSR over fragment bins: entries of bin `b` are
    /// `entries[bin_start[b]..bin_start[b + 1]]`, sorted by the form's mass.
    bin_start: Vec<u32>,
    /// (form id, ion m/z).
    entries: Vec<(u32, f32)>,
}

/// Singly-charged b (prefix) and y (suffix) ion m/z of one peptidoform, from
/// its residue masses with the modification deltas folded in.
fn by_ions(cand: &Candidate) -> Vec<f32> {
    let res = &cand.peptide.residues;
    let n = res.len();
    if n < 2 {
        return Vec::new();
    }
    let masses: Vec<f64> = res
        .iter()
        .map(|aa| aa.mass + aa.mod_.as_ref().map_or(0.0, |m| m.mass_delta))
        .collect();
    let mut out = Vec::with_capacity(2 * (n - 1));
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
    out
}

impl ChunkFragmentIndex {
    /// Enumerate every peptidoform of `records` once and bin its ions. The bin
    /// width is the fragment tolerance at the top of the fragment m/z range
    /// (2,500), so a peak checking its own bin and both neighbours finds every
    /// ion within tolerance whatever its m/z.
    pub fn build(
        records: Vec<IndexRecord>,
        db: &SearchIndex,
        params: &SearchParams,
        fragment_tol: Tolerance,
    ) -> Self {
        let bin_width = fragment_tol.as_da(2500.0).max(0.001);
        let per_record: Vec<Vec<FormIons>> = records
            .par_iter()
            .enumerate()
            .map(|(r, rec)| {
                expand_base_record(db, params, rec)
                    .iter()
                    .enumerate()
                    .map(|(k, cand)| FormIons {
                        record: r as u32,
                        k: k as u16,
                        mass: cand.peptide.mass(),
                        ions: by_ions(cand),
                    })
                    .collect()
            })
            .collect();
        let n_forms: usize = per_record.iter().map(|v| v.len()).sum();
        let mut form_record = Vec::with_capacity(n_forms);
        let mut form_k = Vec::with_capacity(n_forms);
        let mut form_mass = Vec::with_capacity(n_forms);
        let mut max_bin = 0usize;
        let mut n_entries = 0usize;
        for f in per_record.iter().flatten() {
            form_record.push(f.record);
            form_k.push(f.k);
            form_mass.push(f.mass);
            n_entries += f.ions.len();
            for &mz in &f.ions {
                max_bin = max_bin.max((mz as f64 / bin_width) as usize);
            }
        }
        let n_bins = max_bin + 2;
        // Counting sort of (bin, form, ion) into CSR, then sort each bin by
        // form mass so a query binary-searches its precursor window.
        let mut counts = vec![0u32; n_bins + 1];
        for f in per_record.iter().flatten() {
            for &mz in &f.ions {
                counts[(mz as f64 / bin_width) as usize + 1] += 1;
            }
        }
        for b in 0..n_bins {
            counts[b + 1] += counts[b];
        }
        let bin_start = counts.clone();
        let mut fill = counts;
        let mut entries = vec![(0u32, 0f32); n_entries];
        for (form_id, f) in per_record.iter().flatten().enumerate() {
            for &mz in &f.ions {
                let b = (mz as f64 / bin_width) as usize;
                entries[fill[b] as usize] = (form_id as u32, mz);
                fill[b] += 1;
            }
        }
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
