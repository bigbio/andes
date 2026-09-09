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
//!
//! Sizing, measured on the PXD007653 phospho search space: ~88k–120k
//! peptidoforms per Da, so a 150 Da slice is ~13M forms and ~1G ion entries.
//! Forms are enumerated with a mass-bounded walk (only forms inside the
//! slice's window are visited), form ids are assigned in mass order so each
//! fragment bin sorts by a plain integer, and nothing per form is allocated.

use rayon::prelude::*;
use rustc_hash::FxHashMap;

use crate::candidate_gen::{
    expand_base_record_bounded, for_each_record_form_masses_bounded, BaseRecordKey, Candidate,
};
use crate::candidate_index::IndexRecord;
use crate::precursor_cal::adjusted_observed_neutral_mass;
use crate::search_index::SearchIndex;
use crate::search_params::SearchParams;
use model::mass::{H2O, ISOTOPE, PROTON};
use model::spectrum::Spectrum;
use model::tolerance::Tolerance;

/// Pass-1 output per record: (mass, pruned k) per in-window form, and the
/// record's (bin, count) pairs.
type RecordPass1 = (Vec<(f64, u32)>, Vec<(u32, u32)>);

pub struct ChunkFragmentIndex {
    /// The chunk's distinct base records, in `enumerate_candidates` order.
    records: Vec<IndexRecord>,
    /// The slice's precursor mass window; forms outside it are not indexed,
    /// and materialisation walks the same bounded enumeration.
    mass_lo: f64,
    mass_hi: f64,
    /// Per form id (ids are in ascending mass order): its record and its
    /// pruned index within that record's bounded enumeration.
    form_record: Vec<u32>,
    form_k: Vec<u32>,
    form_mass: Vec<f64>,
    bin_width: f64,
    /// CSR over fragment bins: entries of bin `b` are
    /// `entries[bin_start[b]..bin_start[b + 1]]`, sorted by form id (= mass).
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
    /// Enumerate every peptidoform of `records` with neutral mass in
    /// `[mass_lo, mass_hi]` once and bin its ions. The bin width is the
    /// fragment tolerance at the top of the fragment m/z range (2,500), so a
    /// peak checking its own bin and both neighbours finds every ion within
    /// tolerance whatever its m/z.
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

        // Pass 1: every in-window form's (mass, record, k) plus per-record bin
        // counts. The walk visits only in-window subtrees.
        let per: Vec<RecordPass1> = records
            .par_iter()
            .map(|rec| {
                let mut forms: Vec<(f64, u32)> = Vec::new();
                let mut bins: Vec<(u32, u32)> = Vec::new();
                let mut ions: Vec<f32> = Vec::new();
                for_each_record_form_masses_bounded(
                    db,
                    params,
                    rec,
                    mass_lo,
                    mass_hi,
                    |k, masses| {
                        forms.push((neutral_mass(masses), k as u32));
                        ions.clear();
                        by_ions_from_masses(masses, &mut ions);
                        for &mz in &ions {
                            bins.push(((mz as f64 / bin_width) as u32, 1));
                        }
                    },
                );
                bins.sort_unstable();
                bins.dedup_by(|a, b| {
                    if a.0 == b.0 {
                        b.1 += a.1;
                        true
                    } else {
                        false
                    }
                });
                (forms, bins)
            })
            .collect();
        let max_bin = per
            .iter()
            .flat_map(|(_, b)| b.iter().map(|x| x.0 as usize))
            .max()
            .unwrap_or(0);
        let n_bins = max_bin + 2;
        let mut bin_start = vec![0u64; n_bins + 1];
        for (_, bins) in &per {
            for &(b, c) in bins {
                bin_start[b as usize + 1] += c as u64;
            }
        }
        for b in 0..n_bins {
            bin_start[b + 1] += bin_start[b];
        }
        let n_entries = bin_start[n_bins] as usize;

        // Form ids in ascending mass order, so a bin sorted by id is sorted by
        // mass and the per-bin sort is a plain integer sort.
        let mut all: Vec<(f64, u32, u32)> = Vec::new();
        for (r, (forms, _)) in per.iter().enumerate() {
            all.extend(forms.iter().map(|&(m, k)| (m, r as u32, k)));
        }
        all.par_sort_unstable_by(|a, b| {
            a.0.partial_cmp(&b.0)
                .unwrap()
                .then(a.1.cmp(&b.1))
                .then(a.2.cmp(&b.2))
        });
        let n_forms = all.len();
        let mut form_record = Vec::with_capacity(n_forms);
        let mut form_k = Vec::with_capacity(n_forms);
        let mut form_mass = Vec::with_capacity(n_forms);
        // ids_by_record[r][k] = form id, for pass 2.
        let mut ids_by_record: Vec<Vec<u32>> = per
            .iter()
            .map(|(forms, _)| vec![0u32; forms.len()])
            .collect();
        for (id, &(m, r, k)) in all.iter().enumerate() {
            form_record.push(r);
            form_k.push(k);
            form_mass.push(m);
            ids_by_record[r as usize][k as usize] = id as u32;
        }
        drop(all);
        eprintln!(
            "fragment-index: chunk window {:.1}-{:.1} Da, {} records, {} forms, {} ion entries, {} bins",
            mass_lo,
            mass_hi,
            records.len(),
            n_forms,
            n_entries,
            n_bins
        );

        // Pass 2: fill the packed entries; per-bin write cursors are atomics so
        // records fill in parallel, then each bin sorts by form id.
        let fill: Vec<AtomicU64> = bin_start.iter().map(|&s| AtomicU64::new(s)).collect();
        let entries: Vec<(AtomicU32, AtomicU32)> = (0..n_entries)
            .map(|_| (AtomicU32::new(0), AtomicU32::new(0)))
            .collect();
        records.par_iter().enumerate().for_each(|(r, rec)| {
            let ids = &ids_by_record[r];
            let mut ions: Vec<f32> = Vec::new();
            for_each_record_form_masses_bounded(db, params, rec, mass_lo, mass_hi, |k, masses| {
                let form_id = ids[k];
                ions.clear();
                by_ions_from_masses(masses, &mut ions);
                for &mz in &ions {
                    let b = (mz as f64 / bin_width) as usize;
                    let pos = fill[b].fetch_add(1, Ordering::Relaxed) as usize;
                    entries[pos].0.store(form_id, Ordering::Relaxed);
                    entries[pos].1.store(mz.to_bits(), Ordering::Relaxed);
                }
            });
        });
        let mut entries: Vec<(u32, f32)> = entries
            .into_iter()
            .map(|(a, b)| (a.into_inner(), f32::from_bits(b.into_inner())))
            .collect();
        {
            // Sort each bin by form id, in parallel over bins.
            let starts = &bin_start;
            let base = entries.as_mut_ptr() as usize;
            (0..n_bins).into_par_iter().for_each(|b| {
                let (lo, hi) = (starts[b] as usize, starts[b + 1] as usize);
                if hi > lo + 1 {
                    // SAFETY: bins are disjoint ranges of `entries`, and no
                    // other reference to `entries` is live during this loop.
                    let seg = unsafe {
                        std::slice::from_raw_parts_mut((base as *mut (u32, f32)).add(lo), hi - lo)
                    };
                    seg.sort_unstable_by_key(|e| e.0);
                }
            });
        }
        Self {
            records,
            mass_lo,
            mass_hi,
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
        // Form ids are mass-ordered, so the window is an id range.
        let id_lo = self.form_mass.partition_point(|&m| m < lo) as u32;
        let id_hi = self.form_mass.partition_point(|&m| m <= hi) as u32;
        if id_lo >= id_hi {
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
                let start = seg.partition_point(|e| e.0 < id_lo);
                for e in &seg[start..] {
                    if e.0 >= id_hi {
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

    /// Materialise selected forms as candidates, in (record, k) order, by
    /// re-walking each record's bounded enumeration (the same walk the build
    /// used, so `k` selects the same form).
    pub fn materialise(
        &self,
        selected: &[(u32, u16)],
        db: &SearchIndex,
        params: &SearchParams,
        _cache: Option<&FxHashMap<BaseRecordKey, Vec<Candidate>>>,
        relabel: &dyn Fn(&mut [Candidate]),
    ) -> Vec<Candidate> {
        let mut by_record: Vec<(u32, u32)> = selected
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
            let mut forms = expand_base_record_bounded(db, params, rec, self.mass_lo, self.mass_hi);
            relabel(&mut forms);
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
