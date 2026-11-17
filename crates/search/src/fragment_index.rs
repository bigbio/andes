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

struct GroupPass1 {
    forms: Vec<(f64, u32, u32)>,
    forms_per_record: Vec<u32>,
    counts: Vec<u32>,
}

type VoteBuffers = (Vec<(u16, f32)>, Vec<u32>);

thread_local! {
    static VOTE_BUF: std::cell::RefCell<VoteBuffers> =
        const { std::cell::RefCell::new((Vec::new(), Vec::new())) };
}

pub struct ChunkFragmentIndex {
    records: Vec<IndexRecord>,
    tables: FormWalkTables,
    mass_lo: f64,
    mass_hi: f64,
    form_record: Vec<u32>,
    form_k: Vec<u32>,
    form_mass: Vec<f64>,
    bin_width: f64,
    bin_start: Vec<u64>,
    entries: Vec<(u32, f32)>,
    n_blocks: usize,
    block_first_id: Vec<u32>,
    block_off: Vec<u32>,
}

const ENTRIES_PER_CELL: usize = 16;

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
    pub fn build(
        records: Vec<IndexRecord>,
        db: &SearchIndex,
        params: &SearchParams,
        fragment_tol: Tolerance,
        mass_lo: f64,
        mass_hi: f64,
    ) -> Self {
        let bin_width = fragment_tol.as_da(2500.0).max(0.001);
        let tables = FormWalkTables::new(params);
        let records: Vec<IndexRecord> = records
            .into_par_iter()
            .filter(|rec| record_may_have_forms_in(db, &tables, rec, mass_lo, mass_hi))
            .collect();
        let n_bins = if mass_hi > 0.0 {
            ((mass_hi + PROTON + 1.0) / bin_width) as usize + 2
        } else {
            2
        };

        let n_groups = (rayon::current_num_threads() * 4).clamp(1, records.len().max(1));
        let group_len = records.len().div_ceil(n_groups).max(1);
        let groups: Vec<std::ops::Range<usize>> = (0..records.len())
            .step_by(group_len)
            .map(|a| a..(a + group_len).min(records.len()))
            .collect();

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
                                    unsafe {
                                        (base as *mut (u32, f32)).add(pos).write((form_id, mz));
                                    }
                                }
                            },
                        );
                    }
                });
        }
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
                exact.push((
                    w_lo - PRECURSOR_GATE_MARGIN_DA,
                    w_hi + PRECURSOR_GATE_MARGIN_DA,
                ));
            }
        }
        if lo > hi {
            return Vec::new();
        }
        let complement_masses: Vec<f64> = match (params.open_search.is_some(), charges.len()) {
            (true, 1) => {
                let z = f64::from(charges[0]);
                let obs = adjusted_observed_neutral_mass(
                    spec.precursor_mz * z - z * PROTON,
                    shift_ppm,
                );
                vec![obs]
            }
            _ => Vec::new(),
        };
        let mut sel = self.vote(
            spec,
            lo,
            hi,
            fragment_tol,
            min_matched,
            vote_peaks,
            &complement_masses,
        );
        select_top_k(&mut sel, top_k);
        sel.into_iter()
            .filter(|&(id, _, _)| {
                let m = self.form_mass[id as usize];
                exact.iter().any(|&(a, b)| m >= a && m <= b)
            })
            .map(|(id, v, _)| (id, v))
            .collect()
    }

    #[allow(clippy::too_many_arguments)]
    pub fn query_mass_range(
        &self,
        spec: &Spectrum,
        lo: f64,
        hi: f64,
        fragment_tol: Tolerance,
        top_k: usize,
        min_matched: u16,
        vote_peaks: usize,
    ) -> Vec<(u32, u16)> {
        if lo > hi {
            return Vec::new();
        }
        let mut sel = self.vote(
            spec,
            lo,
            hi,
            fragment_tol,
            min_matched,
            vote_peaks,
            &[],
        );
        select_top_k(&mut sel, top_k);
        sel.into_iter().map(|(id, v, _)| (id, v)).collect()
    }

    pub fn form_mass(&self, id: u32) -> f64 {
        self.form_mass[id as usize]
    }

    #[allow(clippy::too_many_arguments)]
    fn vote(
        &self,
        spec: &Spectrum,
        lo: f64,
        hi: f64,
        fragment_tol: Tolerance,
        min_matched: u16,
        vote_peaks: usize,
        complement_masses: &[f64],
    ) -> Vec<(u32, u16, f32)> {
        let id_lo = self.form_mass.partition_point(|&m| m < lo) as u32;
        let id_hi = self.form_mass.partition_point(|&m| m <= hi) as u32;
        if id_lo >= id_hi {
            return Vec::new();
        }
        let n_bins = self.bin_start.len() - 1;
        let width = (id_hi - id_lo) as usize;
        let block_of = |id: u32| self.block_first_id.partition_point(|&f| f <= id) - 1;
        let (blk_lo, blk_hi) = (block_of(id_lo), block_of(id_hi - 1));
        let row_len = self.n_blocks + 1;
        VOTE_BUF.with(|cell| {
            let mut buf = cell.borrow_mut();
            let (counts, touched) = &mut *buf;
            if counts.len() < width {
                counts.resize(width, (0, 0.0));
            }
            touched.clear();
            let floor = vote_intensity_floor(&spec.peaks, vote_peaks);
            let vote = |mz: f64,
                            intensity: f32,
                            counts: &mut Vec<(u16, f32)>,
                            touched: &mut Vec<u32>,
                            clip: bool| {
                let tol_da = fragment_tol.as_da(mz);
                let b = (mz / self.bin_width) as usize;
                if clip && b >= n_bins {
                    return;
                }
                let span = ((tol_da / self.bin_width).ceil() as usize).max(1);
                for bb in b.saturating_sub(span)..=(b + span).min(n_bins - 1) {
                    let base = self.bin_start[bb] as usize;
                    let row = &self.block_off[bb * row_len..(bb + 1) * row_len];
                    let seg =
                        &self.entries[base + row[blk_lo] as usize..base + row[blk_hi + 1] as usize];
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
            };
            for &(mz, intensity) in &spec.peaks {
                if intensity < floor {
                    continue;
                }
                vote(mz, intensity, counts, touched, false);
            }
            for &base in complement_masses {
                for &(mz, intensity) in &spec.peaks {
                    if intensity < floor {
                        continue;
                    }
                    let q = base + 2.0 * PROTON - mz;
                    if q <= 0.0 {
                        continue;
                    }
                    vote(q, intensity, counts, touched, true);
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
        })
    }

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

const PRECURSOR_GATE_MARGIN_DA: f64 = 1e-6;

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

pub fn index_order_mass(s: &Spectrum) -> f64 {
    let z = s.precursor_charge.filter(|z| *z > 0).unwrap_or(2) as f64;
    s.precursor_mz * z
}

pub const INDEX_CHUNK_SIZE: usize = 20_000;

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

pub fn covers_activation(activation: model::activation::ActivationMethod) -> bool {
    use model::activation::ActivationMethod;
    matches!(
        activation,
        ActivationMethod::CID
            | ActivationMethod::HCD
            | ActivationMethod::PQD
            | ActivationMethod::EThcD
    )
}

pub fn is_high_resolution(fragment_tol: Tolerance) -> bool {
    match fragment_tol {
        Tolerance::Ppm(_) => true,
        Tolerance::Da(d) => d <= 0.05,
    }
}

pub fn vote_peaks_for(fragment_tol: Tolerance) -> usize {
    if is_high_resolution(fragment_tol) {
        0
    } else {
        crate::search_params::FRAGMENT_INDEX_LOWRES_VOTE_PEAKS
    }
}

pub(crate) fn vote_intensity_floor(peaks: &[(f64, f32)], vote_peaks: usize) -> f32 {
    if vote_peaks == 0 || peaks.len() <= vote_peaks {
        return 0.0;
    }
    let mut intensities: Vec<f32> = peaks.iter().map(|p| p.1).collect();
    let (_, nth, _) = intensities.select_nth_unstable_by(vote_peaks - 1, |a, b| b.total_cmp(a));
    *nth
}

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
    if sel.len() > top_k {
        sel.select_nth_unstable_by(top_k - 1, order);
        sel.truncate(top_k);
    }
    sel.sort_unstable_by(order);
}

#[cfg(test)]
mod select_top_k_tests {
    use super::select_top_k;

    #[test]
    fn a_tie_is_broken_by_evidence_not_by_position_in_the_index() {
        let mut by_evidence = vec![(7u32, 5u16, 900.0f32), (3u32, 5u16, 10.0f32)];
        select_top_k(&mut by_evidence, 1);
        assert_eq!(
            by_evidence[0].0, 7,
            "with the intensity tie-break the better-supported form survives"
        );
    }

    #[test]
    fn a_higher_vote_count_still_wins_outright() {
        let mut sel = vec![(1u32, 9u16, 1.0f32), (2u32, 4u16, 5000.0f32)];
        select_top_k(&mut sel, 1);
        assert_eq!(sel[0].0, 1, "9 votes must beat 4 regardless of intensity");
    }

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
            assert_eq!(
                off as usize,
                got.partition_point(|&id| id < first[k]),
                "block {k}"
            );
        }
        assert_eq!(*offs.last().unwrap() as usize, ids.len());
        assert!(
            seg.iter().all(|e| e.1 == e.0 as f32),
            "entries keep their m/z"
        );
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
