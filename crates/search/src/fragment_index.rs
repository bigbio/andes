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
    expand_base_record_bounded, for_each_record_form_masses_bounded, record_may_have_forms_in,
    BaseRecordKey, Candidate, FormWalkTables,
};
use crate::candidate_index::IndexRecord;
use crate::precursor_cal::adjusted_observed_neutral_mass;
use crate::search_index::SearchIndex;
use crate::search_params::SearchParams;
use model::mass::{H2O, ISOTOPE, PROTON};
use model::spectrum::Spectrum;
use model::tolerance::Tolerance;

/// Pass-1 output of one group of records.
struct GroupPass1 {
    /// (mass, record, pruned k) of every in-window form, in record then k order.
    forms: Vec<(f64, u32, u32)>,
    /// In-window form count of each record in the group.
    forms_per_record: Vec<u32>,
    /// The group's ion count per fragment bin.
    counts: Vec<u32>,
}

/// `query`'s reusable vote buffers: per-id (count, summed intensity) over the
/// current window, zeroed again after each query, and the ids touched.
type VoteBuffers = (Vec<(u16, f32)>, Vec<u32>);

thread_local! {
    static VOTE_BUF: std::cell::RefCell<VoteBuffers> =
        const { std::cell::RefCell::new((Vec::new(), Vec::new())) };
}

pub struct ChunkFragmentIndex {
    /// The chunk's distinct base records, in `enumerate_candidates` order.
    records: Vec<IndexRecord>,
    /// Variant tables for the bounded walk, shared by build and materialisation.
    tables: FormWalkTables,
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
    /// Form ids are cut into `n_blocks` precursor-mass blocks of equal width;
    /// block `k` holds ids `block_first_id[k]..block_first_id[k + 1]`.
    n_blocks: usize,
    block_first_id: Vec<u32>,
    /// Per bin, `n_blocks + 1` offsets into the bin's entries: block `k`'s
    /// entries start at `block_off[b * (n_blocks + 1) + k]`. A query reads its
    /// window's place in a bin from here instead of searching the bin.
    block_off: Vec<u32>,
}

/// Average entries per (bin, precursor-mass block) cell the block count aims
/// for; the offset table then costs about 4 bytes per this many entries.
const ENTRIES_PER_CELL: usize = 16;

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
        // NOTE: `fragment_tol` here is `RankScorer::feature_match_tolerance()`, a
        // CONSTANT 20 ppm on high-resolution data — not the model's `mme`, which
        // scoring matches at and which is 0.5 Da in every bundled model. Retrieval
        // is therefore far tighter than scoring, in both directions, and the effect
        // on identifications has never been measured. See "Known gaps" in
        // docs/benchmarks/README.md before changing this.
        let bin_width = fragment_tol.as_da(2500.0).max(0.001);
        let tables = FormWalkTables::new(params);
        // Records whose walk cannot reach the window contribute nothing; drop
        // them once here instead of in both passes. The order of the rest is
        // kept, so form ids and materialisation order are unchanged.
        let records: Vec<IndexRecord> = records
            .into_par_iter()
            .filter(|rec| record_may_have_forms_in(db, &tables, rec, mass_lo, mass_hi))
            .collect();
        // Every singly-charged b/y ion of an in-window form lies below
        // `mass_hi + PROTON`; one Da of margin covers rounding.
        let n_bins = if mass_hi > 0.0 {
            ((mass_hi + PROTON + 1.0) / bin_width) as usize + 2
        } else {
            2
        };

        // Records are split into contiguous groups; each group is walked by one
        // task with its own dense per-bin counters, so neither pass shares a
        // write cursor with another thread.
        let n_groups = (rayon::current_num_threads() * 4).clamp(1, records.len().max(1));
        let group_len = records.len().div_ceil(n_groups).max(1);
        let groups: Vec<std::ops::Range<usize>> = (0..records.len())
            .step_by(group_len)
            .map(|a| a..(a + group_len).min(records.len()))
            .collect();

        // Pass 1, per group: every in-window form's (mass, record, k), each
        // record's form count, and the group's ion count per bin. The walk
        // visits only in-window subtrees.
        let pass1: Vec<GroupPass1> = groups
            .par_iter()
            .map(|range| {
                let mut forms: Vec<(f64, u32, u32)> = Vec::new();
                let mut forms_per_record: Vec<u32> = Vec::with_capacity(range.len());
                let mut counts = vec![0u32; n_bins];
                let mut ions: Vec<f32> = Vec::new();
                for r in range.clone() {
                    let before = forms.len();
                    for_each_record_form_masses_bounded(
                        db,
                        params,
                        &tables,
                        &records[r],
                        mass_lo,
                        mass_hi,
                        |k, masses| {
                            forms.push((neutral_mass(masses), r as u32, k as u32));
                            ions.clear();
                            by_ions_from_masses(masses, &mut ions);
                            for &mz in &ions {
                                counts[(mz as f64 / bin_width) as usize] += 1;
                            }
                        },
                    );
                    forms_per_record.push((forms.len() - before) as u32);
                }
                GroupPass1 {
                    forms,
                    forms_per_record,
                    counts,
                }
            })
            .collect();

        // CSR bin starts, and each group's first write offset within every bin
        // (the ions of earlier groups in that bin), stored in place of its counts.
        let mut cursors: Vec<Vec<u32>> = Vec::with_capacity(pass1.len());
        let mut bin_start = vec![0u64; n_bins + 1];
        let mut running = vec![0u32; n_bins];
        let mut forms_by_group: Vec<Vec<(f64, u32, u32)>> = Vec::with_capacity(pass1.len());
        let mut record_form_start: Vec<u64> = Vec::with_capacity(records.len() + 1);
        record_form_start.push(0);
        for g in pass1 {
            let mut c = g.counts;
            for (slot, run) in c.iter_mut().zip(running.iter_mut()) {
                let n = *slot;
                *slot = *run;
                *run += n;
            }
            cursors.push(c);
            for n in g.forms_per_record {
                let last = *record_form_start.last().unwrap();
                record_form_start.push(last + n as u64);
            }
            forms_by_group.push(g.forms);
        }
        for b in 0..n_bins {
            bin_start[b + 1] = bin_start[b] + running[b] as u64;
        }
        drop(running);
        let n_entries = bin_start[n_bins] as usize;

        // Form ids in ascending mass order, so a bin sorted by id is sorted by
        // mass and the per-bin sort is a plain integer sort.
        let mut all: Vec<(f64, u32, u32)> = forms_by_group.into_iter().flatten().collect();
        all.par_sort_unstable_by(|a, b| {
            a.0.partial_cmp(&b.0)
                .unwrap()
                .then(a.1.cmp(&b.1))
                .then(a.2.cmp(&b.2))
        });
        let n_forms = all.len();
        let form_record: Vec<u32> = all.par_iter().map(|f| f.1).collect();
        let form_k: Vec<u32> = all.par_iter().map(|f| f.2).collect();
        let form_mass: Vec<f64> = all.par_iter().map(|f| f.0).collect();
        // form_ids[record_form_start[r] + k] = form id, for pass 2.
        let mut form_ids = vec![0u32; n_forms];
        for (id, &(_, r, k)) in all.iter().enumerate() {
            form_ids[(record_form_start[r as usize] + k as u64) as usize] = id as u32;
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

        // Pass 2: each group re-walks its records and writes its ions at its
        // own cursors, then each bin sorts by form id.
        let mut entries: Vec<(u32, f32)> = vec![(0, 0.0); n_entries];
        {
            let base = entries.as_mut_ptr() as usize;
            groups
                .par_iter()
                .zip(cursors.into_par_iter())
                .for_each(|(range, mut cursor)| {
                    let mut ions: Vec<f32> = Vec::new();
                    for r in range.clone() {
                        let first = record_form_start[r] as usize;
                        for_each_record_form_masses_bounded(
                            db,
                            params,
                            &tables,
                            &records[r],
                            mass_lo,
                            mass_hi,
                            |k, masses| {
                                let form_id = form_ids[first + k];
                                ions.clear();
                                by_ions_from_masses(masses, &mut ions);
                                for &mz in &ions {
                                    let b = (mz as f64 / bin_width) as usize;
                                    let pos = bin_start[b] as usize + cursor[b] as usize;
                                    cursor[b] += 1;
                                    debug_assert!(pos < bin_start[b + 1] as usize);
                                    // SAFETY: pass 1 counted exactly these ions
                                    // for this group, so `pos` lies in this
                                    // group's own slots of bin `b`, which no
                                    // other task writes.
                                    unsafe {
                                        (base as *mut (u32, f32)).add(pos).write((form_id, mz));
                                    }
                                }
                            },
                        );
                    }
                });
        }
        // Cut the form ids into precursor-mass blocks, sort each bin by form
        // id, and record where each block starts in it. Parallel over bins.
        let n_blocks = (n_entries / (n_bins * ENTRIES_PER_CELL)).clamp(1, u16::MAX as usize);
        let block_first_id: Vec<u32> = (0..=n_blocks)
            .map(|k| {
                if k == 0 {
                    0
                } else if k == n_blocks {
                    n_forms as u32
                } else {
                    let edge = mass_lo + (mass_hi - mass_lo) * (k as f64) / (n_blocks as f64);
                    form_mass.partition_point(|&m| m < edge) as u32
                }
            })
            .collect();
        let mut block_off = vec![0u32; n_bins * (n_blocks + 1)];
        {
            let starts = &bin_start;
            let first = &block_first_id;
            let base = entries.as_mut_ptr() as usize;
            block_off
                .par_chunks_mut(n_blocks + 1)
                .enumerate()
                .for_each(|(b, offs)| {
                    let (lo, hi) = (starts[b] as usize, starts[b + 1] as usize);
                    // SAFETY: bins are disjoint ranges of `entries`, and no
                    // other reference to `entries` is live during this loop.
                    let seg = unsafe {
                        std::slice::from_raw_parts_mut((base as *mut (u32, f32)).add(lo), hi - lo)
                    };
                    sort_bin_by_block(seg, first, offs);
                });
        }
        Self {
            records,
            tables,
            mass_lo,
            mass_hi,
            form_record,
            form_k,
            form_mass,
            bin_width,
            bin_start,
            entries,
            n_blocks,
            block_first_id,
            block_off,
        }
    }

    pub fn n_forms(&self) -> usize {
        self.form_mass.len()
    }

    /// Forms with at least `min_matched` singly-charged b/y ions within
    /// `fragment_tol` of a peak (ppm tolerances are evaluated per peak),
    /// restricted to the spectrum's precursor windows over every charge in
    /// `charges` and every isotope offset in `params`. At most `top_k`, best
    /// votes first. Only the `vote_peaks` most intense peaks vote (0 = every peak).
    ///
    /// Votes are counted over one interval spanning every (charge, isotope
    /// offset) window, and `top_k` is taken over that interval. Of the kept
    /// forms, only those inside at least one exact window are returned, since
    /// the scorer's precursor test rejects every other one; the test here is
    /// widened by [`PRECURSOR_GATE_MARGIN_DA`] so it never drops a form the
    /// scorer would accept.
    #[allow(clippy::too_many_arguments)]
    pub fn query(
        &self,
        spec: &Spectrum,
        charges: &[u8],
        params: &SearchParams,
        fragment_tol: Tolerance,
        top_k: usize,
        min_matched: u16,
        vote_peaks: usize,
    ) -> Vec<(u32, u16)> {
        // One over-inclusive precursor mass interval; the scoring loop applies
        // the exact per-offset test afterwards.
        let shift_ppm = params.precursor_mass_shift_ppm;
        let mut lo = f64::MAX;
        let mut hi = f64::MIN;
        let mut exact: smallvec::SmallVec<[(f64, f64); 16]> = smallvec::SmallVec::new();
        for &z in charges {
            let zf = z as f64;
            let obs =
                adjusted_observed_neutral_mass(spec.precursor_mz * zf - zf * PROTON, shift_ppm);
            for o in params.isotope_error_range.clone() {
                let c = obs - (o as f64) * ISOTOPE;
                let (w_lo, w_hi) = (
                    c - params.precursor_tolerance.left.as_da(c),
                    c + params.precursor_tolerance.right.as_da(c),
                );
                lo = lo.min(w_lo);
                hi = hi.max(w_hi);
                exact.push((w_lo - PRECURSOR_GATE_MARGIN_DA, w_hi + PRECURSOR_GATE_MARGIN_DA));
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
        // Vote count, plus the summed intensity of the peaks that produced it.
        // The count stays the primary key — it is what `min_matched` means — and the
        // intensity is only ever used to break ties among equal counts.
        // Votes go into a dense per-thread array over the window's id range, plus
        // the list of ids touched, instead of a hash map grown per spectrum. Each
        // id's sums accumulate in peak order as before, and `select_top_k` fully
        // orders the result, so the selection is identical.
        let n_bins = self.bin_start.len() - 1;
        let width = (id_hi - id_lo) as usize;
        // The blocks holding the first and the last id of the window; in every
        // bin the window's entries lie between their offsets.
        let block_of = |id: u32| self.block_first_id.partition_point(|&f| f <= id) - 1;
        let (blk_lo, blk_hi) = (block_of(id_lo), block_of(id_hi - 1));
        let row_len = self.n_blocks + 1;
        let mut sel: Vec<(u32, u16, f32)> = VOTE_BUF.with(|cell| {
            let mut buf = cell.borrow_mut();
            let (counts, touched) = &mut *buf;
            if counts.len() < width {
                counts.resize(width, (0, 0.0));
            }
            touched.clear();
            let floor = vote_intensity_floor(&spec.peaks, vote_peaks);
            for &(mz, intensity) in &spec.peaks {
                if intensity < floor {
                    continue;
                }
                let tol_da = fragment_tol.as_da(mz);
                let b = (mz / self.bin_width) as usize;
                for bb in b.saturating_sub(1)..=(b + 1).min(n_bins - 1) {
                    let base = self.bin_start[bb] as usize;
                    let row = &self.block_off[bb * row_len..(bb + 1) * row_len];
                    let seg = &self.entries
                        [base + row[blk_lo] as usize..base + row[blk_hi + 1] as usize];
                    for e in seg {
                        if e.0 < id_lo {
                            continue;
                        }
                        if e.0 >= id_hi {
                            break;
                        }
                        if ((e.1 as f64) - mz).abs() <= tol_da {
                            let slot = &mut counts[(e.0 - id_lo) as usize];
                            if slot.0 == 0 {
                                touched.push(e.0);
                            }
                            slot.0 += 1;
                            slot.1 += intensity;
                        }
                    }
                }
            }
            let mut sel = Vec::new();
            for &id in touched.iter() {
                let slot = &mut counts[(id - id_lo) as usize];
                if slot.0 >= min_matched {
                    sel.push((id, slot.0, slot.1));
                }
                *slot = (0, 0.0);
            }
            sel
        });
        select_top_k(&mut sel, top_k);
        sel.into_iter()
            .filter(|&(id, _, _)| {
                let m = self.form_mass[id as usize];
                exact.iter().any(|&(a, b)| m >= a && m <= b)
            })
            .map(|(id, v, _)| (id, v))
            .collect()
    }

    /// Materialise selected forms as candidates, in (record, k) order, by
    /// re-walking each record's bounded enumeration (the same walk the build
    /// used, so `k` selects the same form). Only the selected forms are built.
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
        let mut ks: Vec<u32> = Vec::new();
        let mut i = 0;
        while i < by_record.len() {
            let r = by_record[i].0;
            let mut j = i;
            while j < by_record.len() && by_record[j].0 == r {
                j += 1;
            }
            ks.clear();
            ks.extend(by_record[i..j].iter().map(|&(_, k)| k));
            let rec = &self.records[r as usize];
            let mut forms = expand_base_record_bounded(
                db,
                params,
                &self.tables,
                rec,
                self.mass_lo,
                self.mass_hi,
                &ks,
            );
            relabel(&mut forms);
            out.extend(forms);
            i = j;
        }
        out
    }
}

/// Slack on the exact precursor windows `query` filters by. A form's indexed
/// mass and its candidate's `Peptide::mass` are the same residue masses summed
/// in possibly different order, so they differ by float rounding only; the
/// scorer's precursor test stays the authority on which forms match.
const PRECURSOR_GATE_MARGIN_DA: f64 = 1e-6;

/// Sort one bin's entries by form id and write the bin's block offsets to
/// `offs`: `offs[k]` is the first entry whose id is at or past
/// `block_first_id[k]`, and `offs[n_blocks]` is the bin's length. Entries with
/// equal ids are one form's ions and may land in either order.
fn sort_bin_by_block(seg: &mut [(u32, f32)], block_first_id: &[u32], offs: &mut [u32]) {
    seg.sort_unstable_by_key(|e| e.0);
    let n_blocks = offs.len() - 1;
    let mut i = 0usize;
    for k in 0..n_blocks {
        while i < seg.len() && seg[i].0 < block_first_id[k] {
            i += 1;
        }
        offs[k] = i as u32;
    }
    offs[n_blocks] = seg.len() as u32;
}

/// Mass a spectrum is ordered by for fragment-index chunking: precursor m/z
/// times the reported charge (2 when unknown).
pub fn index_order_mass(s: &Spectrum) -> f64 {
    let z = s.precursor_charge.filter(|z| *z > 0).unwrap_or(2) as f64;
    s.precursor_mz * z
}

/// Most spectra in one fragment-index chunk.
pub const INDEX_CHUNK_SIZE: usize = 20_000;

/// Length of the next fragment-index chunk at the head of `rest` (sorted by
/// `mass`): at most [`INDEX_CHUNK_SIZE`] items within `span_da` of the first.
pub fn index_chunk_len<T>(rest: &[T], mass: impl Fn(&T) -> f64, span_da: f64) -> usize {
    let Some(head) = rest.first() else {
        return 0;
    };
    let first = mass(head);
    let mut take = 1;
    while take < rest.len() && take < INDEX_CHUNK_SIZE && mass(&rest[take]) - first <= span_da {
        take += 1;
    }
    take
}

/// Whether a fragment tolerance is tight enough for the index to bin usefully.
pub fn is_high_resolution(fragment_tol: Tolerance) -> bool {
    match fragment_tol {
        Tolerance::Ppm(_) => true,
        Tolerance::Da(d) => d <= 0.05,
    }
}

/// How many of a spectrum's most intense peaks vote in `query` (0 = every peak).
///
/// At a low-resolution tolerance a bin spans most of a Dalton, so every noise
/// peak lands on some form's ion; with all peaks voting the noise swamps the
/// count and the index loses most identifications. Letting only the
/// `LOWRES_VOTE_PEAKS` most intense peaks vote restores them. High-resolution
/// spectra are not capped: on dense Astral spectra the real fragments reach
/// well below the top 150 peaks, and a cap there costs identifications.
pub fn vote_peaks_for(fragment_tol: Tolerance) -> usize {
    if is_high_resolution(fragment_tol) {
        0
    } else {
        crate::search_params::FRAGMENT_INDEX_LOWRES_VOTE_PEAKS
    }
}

/// The intensity a peak needs to be among the `vote_peaks` most intense
/// (0.0, letting every peak vote, when `vote_peaks` is 0 or covers them all).
/// Peaks tied at the cut all vote, so the result does not depend on peak order.
pub(crate) fn vote_intensity_floor(peaks: &[(f64, f32)], vote_peaks: usize) -> f32 {
    if vote_peaks == 0 || peaks.len() <= vote_peaks {
        return 0.0;
    }
    let mut intensities: Vec<f32> = peaks.iter().map(|p| p.1).collect();
    let (_, nth, _) =
        intensities.select_nth_unstable_by(vote_peaks - 1, |a, b| b.total_cmp(a));
    *nth
}

/// Order `(form_id, votes, matched_intensity)` and keep the best `top_k`.
///
/// Vote counts are small integers over a candidate set that can run to thousands,
/// so the cut at `top_k` lands inside a large tie. That tie is ordered by the summed
/// intensity of the matched peaks, so the better-supported forms survive the cut.
/// The form id is the FINAL key, because the order must be total and reproducible —
/// this repo has had an FDR swing from a non-deterministic sort in the candidate path.
pub(crate) fn select_top_k(sel: &mut Vec<(u32, u16, f32)>, top_k: usize) {
    let order = |a: &(u32, u16, f32), b: &(u32, u16, f32)| {
        b.1.cmp(&a.1)
            .then_with(|| b.2.total_cmp(&a.2))
            .then_with(|| a.0.cmp(&b.0))
    };
    if top_k == 0 {
        sel.clear();
        return;
    }
    // The order is total, so partitioning at `top_k` and sorting only the kept
    // prefix selects and orders exactly what a full sort would.
    if sel.len() > top_k {
        sel.select_nth_unstable_by(top_k - 1, order);
        sel.truncate(top_k);
    }
    sel.sort_unstable_by(order);
}

#[cfg(test)]
mod select_top_k_tests {
    use super::select_top_k;

    /// Equal vote counts: the matched intensity decides, not the form id (mass order).
    #[test]
    fn a_tie_is_broken_by_evidence_not_by_position_in_the_index() {
        // Form 7 carries the stronger evidence; form 3 sits earlier in the index,
        // so the two keys disagree.
        let mut by_evidence = vec![(7u32, 5u16, 900.0f32), (3u32, 5u16, 10.0f32)];
        select_top_k(&mut by_evidence, 1);
        assert_eq!(
            by_evidence[0].0, 7,
            "with the intensity tie-break the better-supported form survives"
        );
    }

    /// A real vote difference must still win: intensity is a TIE-break only.
    #[test]
    fn a_higher_vote_count_still_wins_outright() {
        let mut sel = vec![(1u32, 9u16, 1.0f32), (2u32, 4u16, 5000.0f32)];
        select_top_k(&mut sel, 1);
        assert_eq!(sel[0].0, 1, "9 votes must beat 4 regardless of intensity");
    }

    /// The order is total, so a rerun selects the same set.
    #[test]
    fn ordering_is_deterministic_under_a_full_tie() {
        let mut a = vec![
            (5u32, 3u16, 1.0f32),
            (2u32, 3u16, 1.0f32),
            (9u32, 3u16, 1.0f32),
        ];
        let mut b = vec![
            (9u32, 3u16, 1.0f32),
            (5u32, 3u16, 1.0f32),
            (2u32, 3u16, 1.0f32),
        ];
        select_top_k(&mut a, 2);
        select_top_k(&mut b, 2);
        assert_eq!(a, b, "input order must not survive into the selection");
    }
}

#[cfg(test)]
mod sort_bin_by_block_tests {
    use super::sort_bin_by_block;

    /// The bin comes out sorted by id and each block's offset is the first
    /// entry whose id is at or past the block's first id.
    #[test]
    fn sorts_by_id_and_records_block_starts() {
        let first = vec![0u32, 10, 10, 25, 40];
        let ids = [33u32, 2, 39, 10, 0, 24, 11, 2, 30, 9];
        let mut seg: Vec<(u32, f32)> = ids.iter().map(|&i| (i, i as f32)).collect();
        let mut offs = vec![0u32; first.len()];
        sort_bin_by_block(&mut seg, &first, &mut offs);
        let got: Vec<u32> = seg.iter().map(|e| e.0).collect();
        let mut want = ids.to_vec();
        want.sort_unstable();
        assert_eq!(got, want);
        for (k, &off) in offs.iter().enumerate().take(first.len() - 1) {
            assert_eq!(off as usize, got.partition_point(|&id| id < first[k]), "block {k}");
        }
        assert_eq!(*offs.last().unwrap() as usize, ids.len());
        assert!(seg.iter().all(|e| e.1 == e.0 as f32), "entries keep their m/z");
    }

    #[test]
    fn a_single_block_is_a_plain_sort() {
        let mut seg = vec![(3u32, 0.0f32), (1, 0.0), (2, 0.0)];
        let mut offs = vec![0u32; 2];
        sort_bin_by_block(&mut seg, &[0, 4], &mut offs);
        assert_eq!(seg.iter().map(|e| e.0).collect::<Vec<_>>(), vec![1, 2, 3]);
        assert_eq!(offs, vec![0, 3]);
    }
}

#[cfg(test)]
mod vote_peaks_tests {
    use super::{vote_intensity_floor, vote_peaks_for};
    use model::tolerance::Tolerance;

    #[test]
    fn high_resolution_is_never_capped() {
        assert_eq!(vote_peaks_for(Tolerance::Ppm(20.0)), 0);
        assert_eq!(vote_peaks_for(Tolerance::Da(0.02)), 0);
    }

    #[test]
    fn low_resolution_is_capped() {
        assert_eq!(
            vote_peaks_for(Tolerance::Da(0.5)),
            crate::search_params::FRAGMENT_INDEX_LOWRES_VOTE_PEAKS
        );
    }

    #[test]
    fn floor_admits_exactly_the_most_intense_peaks() {
        let peaks = vec![(100.0, 5.0), (200.0, 1.0), (300.0, 9.0), (400.0, 3.0)];
        let floor = vote_intensity_floor(&peaks, 2);
        let voting: Vec<f64> = peaks.iter().filter(|p| p.1 >= floor).map(|p| p.0).collect();
        assert_eq!(voting, vec![100.0, 300.0]);
    }

    #[test]
    fn no_cap_or_a_cap_above_the_peak_count_lets_every_peak_vote() {
        let peaks = vec![(100.0, 5.0), (200.0, 1.0)];
        assert_eq!(vote_intensity_floor(&peaks, 0), 0.0);
        assert_eq!(vote_intensity_floor(&peaks, 2), 0.0);
        assert_eq!(vote_intensity_floor(&peaks, 10), 0.0);
    }

    #[test]
    fn peaks_tied_at_the_cut_all_vote() {
        let peaks = vec![(100.0, 4.0), (200.0, 2.0), (300.0, 2.0), (400.0, 1.0)];
        let floor = vote_intensity_floor(&peaks, 2);
        assert_eq!(peaks.iter().filter(|p| p.1 >= floor).count(), 3);
    }
}
