//! Match-between-runs library confidence.
//!
//! The functions here decide which precursors enter the transfer library. The
//! values are internal quant-quality gates, not reported identification FDR;
//! identification FDR is only ever reported by the PSM-level rescorer.

use crate::align::{align, AlignParams, Alignment, Observation};

/// Precursor-level library entry gate (internal; not an identification FDR).
///
/// `rows` holds `(precursor_id, is_decoy, psm_score)` with higher scores
/// better. Each precursor keeps its best PSM score across all runs, then
/// precursors are ranked by `(score desc, decoys before targets, id asc)` (ties
/// are conservative: a target tied with a decoy never gets a lower q than it
/// would with the decoy first) and run through
/// target/decoy competition: `(D + 1) / max(T, 1)` at each rank, a running
/// minimum from the worst rank upwards, clamped to at most 1.
///
/// A NaN score is treated as the worst possible score. A precursor id is
/// expected to be all-target or all-decoy (checked by `debug_assert!`); in release
/// builds a mixed id is decided by the first row seen. Ids are assumed dense:
/// the result is indexed by precursor id and has length `max id + 1`, `None`
/// for ids without rows.
pub fn precursor_qvalues(rows: &[(u32, bool, f64)]) -> Vec<Option<f64>> {
    let Some(max_id) = rows.iter().map(|r| r.0).max() else {
        return Vec::new();
    };
    let n = max_id as usize + 1;
    // (is_decoy, best score) per id; NaN is mapped to -inf so it is the worst.
    let mut best: Vec<Option<(bool, f64)>> = vec![None; n];
    for &(id, is_decoy, score) in rows {
        let score = if score.is_nan() {
            f64::NEG_INFINITY
        } else {
            score
        };
        let slot = &mut best[id as usize];
        match slot {
            Some((d, s)) => {
                debug_assert_eq!(
                    *d, is_decoy,
                    "precursor id {id} mixes target and decoy rows"
                );
                if score > *s {
                    *s = score;
                }
            }
            None => *slot = Some((is_decoy, score)),
        }
    }
    let mut ranked: Vec<(u32, bool, f64)> = best
        .iter()
        .enumerate()
        .filter_map(|(id, b)| b.map(|(d, s)| (id as u32, d, s)))
        .collect();
    // Score desc; on ties decoys first (conservative), then id for determinism.
    ranked.sort_by(|a, b| b.2.total_cmp(&a.2).then(b.1.cmp(&a.1)).then(a.0.cmp(&b.0)));

    let mut fdr = Vec::with_capacity(ranked.len());
    let (mut t, mut d) = (0u64, 0u64);
    for &(_, is_decoy, _) in &ranked {
        if is_decoy {
            d += 1;
        } else {
            t += 1;
        }
        fdr.push((d + 1) as f64 / t.max(1) as f64);
    }
    let mut running = f64::INFINITY;
    let mut out = vec![None; n];
    for (i, &(id, _, _)) in ranked.iter().enumerate().rev() {
        running = running.min(fdr[i]);
        out[id as usize] = Some(running.min(1.0));
    }
    out
}

/// One PSM of one run, as the library sees it. `precursor` is a dense id in
/// `0..meta.len()` shared by every run.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PsmIn {
    pub run: u16,
    pub precursor: u32,
    pub rt: f64,
    /// The rescorer's PSM score (higher is better).
    pub score: f64,
    /// The rescorer's PSM q-value.
    pub q: f64,
    /// The rescorer's PSM posterior error probability.
    pub pep: f64,
    pub is_decoy: bool,
    pub faims_cv: Option<f32>,
}

/// Whether `p` is a donor PSM under the PSM q-value threshold `donor_q` and the
/// optional PEP cap `donor_pep`.
pub fn is_donor_psm(p: &PsmIn, donor_q: f64, donor_pep: Option<f64>) -> bool {
    p.q <= donor_q && donor_pep.is_none_or(|c| p.pep <= c)
}

/// Identity of a precursor: peptidoform and charge key, monoisotopic m/z and
/// the theoretical isotope fractions.
#[derive(Debug, Clone, PartialEq)]
pub struct PrecursorMeta {
    pub key: String,
    pub mz: f64,
    pub charge: u8,
    pub envelope: [f32; 4],
}

/// The best donor PSM of an entry (highest score; then earliest RT, lowest
/// run): its rescorer score, PSM q-value and PEP.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DonorPsm {
    pub score: f64,
    pub q: f64,
    pub pep: f64,
}

/// A confident precursor that is measured (identified or transferred) in
/// every run. Its donor runs are the runs with a PSM of it at PSM q-value at
/// most the library's donor threshold; only those PSMs count below.
///
/// An identity decoy is a decoy precursor that passes the same entry gate and
/// donor rules with its decoy PSMs. It is measured exactly like a target
/// entry, so its transfers show how often the transfer of a false donor
/// identification is accepted; it is never reported.
#[derive(Debug, Clone, PartialEq)]
pub struct LibraryEntry {
    pub precursor: u32,
    /// True for an identity decoy (a decoy precursor's entry).
    pub identity_decoy: bool,
    /// Library confidence (internal quant gate, not an identification FDR).
    pub q: f64,
    /// The best donor PSM.
    pub donor_psm: DonorPsm,
    /// Consensus RT from the joint alignment (NaN without an observation in
    /// an aligned run).
    pub consensus_rt: f64,
    /// Most frequent FAIMS compensation voltage among its donor PSMs; the
    /// smallest value wins a tie; `None` when no donor PSM has one.
    pub faims_cv: Option<f32>,
    /// Donor runs, ascending.
    pub donors: Vec<u16>,
    /// Per donor run (ascending), the RT of its best-scoring donor PSM
    /// (earliest RT on equal scores).
    pub best_rt: Vec<(u16, f64)>,
}

/// The transfer library: entries sorted by m/z, the joint alignment of the
/// entries' observations, and the precursor table indexed by id.
#[derive(Debug, Clone)]
pub struct Library {
    pub entries: Vec<LibraryEntry>,
    pub alignment: Alignment,
    pub meta: Vec<PrecursorMeta>,
}

/// Builds the library from the PSMs of all runs.
///
/// Precursors are independent: confidence is [`precursor_qvalues`] over the
/// best score per precursor id (every PSM, whatever its PSM q-value), with no
/// per-scan competition. Donor PSMs are the PSMs with PSM `q <= donor_q` and,
/// when `donor_pep` is given, PSM `pep <= donor_pep`. Entries are the
/// precursors with `q <= max_q` and at least one donor PSM, sorted by
/// `(mz, id)`: target precursors, and decoy precursors as identity decoys
/// (the same rules on their decoy PSMs). The alignment is fitted on one
/// observation per `(target entry, donor run)`, the RT of the best donor PSM
/// there; an identity decoy's observations are added after the fit
/// ([`Alignment::add_observations`]), so they never move a run map but give
/// the decoy its consensus RT and predictions as a target's would. PSMs with
/// an id outside `0..meta.len()` or a run outside `0..n_runs` are ignored.
/// The PSMs are added to a [`LibraryBuilder`] in run order, each run's in
/// slice order.
pub fn build_library(
    psms: &[PsmIn],
    meta: Vec<PrecursorMeta>,
    n_runs: usize,
    max_q: f64,
    donor_q: f64,
    donor_pep: Option<f64>,
    ap: &AlignParams,
) -> Library {
    let mut b = LibraryBuilder::new(meta.len(), n_runs, donor_q, donor_pep);
    let mut order: Vec<&PsmIn> = psms.iter().collect();
    order.sort_by_key(|p| p.run);
    for p in order {
        b.add(p);
    }
    b.finish(meta, max_q, ap)
}

/// NaN scores rank worst, as in the confidence gate.
fn score_of(score: f64) -> f64 {
    if score.is_nan() {
        f64::NEG_INFINITY
    } else {
        score
    }
}

/// A donor PSM reduced to what an entry keeps of it.
#[derive(Debug, Clone, Copy)]
struct Donor {
    score: f64,
    rt: f64,
    q: f64,
    pep: f64,
    run: u16,
    is_decoy: bool,
}

impl Donor {
    fn of(p: &PsmIn) -> Self {
        Self {
            score: p.score,
            rt: p.rt,
            q: p.q,
            pep: p.pep,
            run: p.run,
            is_decoy: p.is_decoy,
        }
    }

    /// Whether `self` replaces `best`: a higher score (NaN as the worst),
    /// or an equal score and an earlier RT.
    fn beats(&self, best: &Donor) -> bool {
        let (a, b) = (score_of(self.score), score_of(best.score));
        a > b || (a == b && self.rt < best.rt)
    }
}

/// The donor PSMs of one precursor, reduced to what an entry needs. Runs
/// arrive in ascending order: `current` is the best donor PSM of the latest
/// run; each earlier run left its best donor PSM's `(run, RT)` in `best_rt`
/// and was folded into `best`, the precursor's best donor PSM over those
/// runs. `cvs` counts each FAIMS voltage over all its donor PSMs.
#[derive(Debug, Clone)]
struct DonorRuns {
    current: Donor,
    best: Option<Donor>,
    /// The label of the best donor PSM of the lowest donor run.
    first_is_decoy: bool,
    best_rt: Vec<(u16, f64)>,
    cvs: Vec<(f32, u32)>,
}

impl DonorRuns {
    fn new(d: Donor) -> Self {
        Self {
            current: d,
            best: None,
            first_is_decoy: false,
            best_rt: Vec::new(),
            cvs: Vec::new(),
        }
    }

    /// Ends the run of `current`: its best donor PSM is final.
    fn close_run(&mut self) {
        let c = self.current;
        if self.best_rt.is_empty() {
            self.first_is_decoy = c.is_decoy;
        }
        self.best_rt.push((c.run, c.rt));
        // Runs ascend, so a tie on score and RT keeps the lowest run.
        if self.best.is_none_or(|b| c.beats(&b)) {
            self.best = Some(c);
        }
    }
}

/// No donor slot.
const NO_SLOT: u32 = u32::MAX;

/// [`build_library`] over a stream of PSMs: [`LibraryBuilder::add`] takes
/// them one at a time and keeps, per precursor, its best score and, per
/// donor run, the RT of its best donor PSM, plus its best donor PSM and
/// voltage counts. Memory follows the precursors and their donor runs (16
/// bytes per donor run) rather than the PSMs. Each precursor's PSMs must
/// arrive in ascending run order; adding the PSMs of [`build_library`]'s
/// slice sorted by run (stably) gives the same library.
#[derive(Debug, Clone)]
pub struct LibraryBuilder {
    n_runs: usize,
    donor_q: f64,
    donor_pep: Option<f64>,
    /// Per precursor id: the label of its first PSM and its best score (NaN
    /// as the worst), over every PSM.
    best: Vec<Option<(bool, f64)>>,
    /// Per precursor id: its index in `donors`, or [`NO_SLOT`].
    slot: Vec<u32>,
    /// The donor PSMs of the precursors with one, in order of their first.
    donors: Vec<DonorRuns>,
}

impl LibraryBuilder {
    /// A builder for `n_precursors` precursor ids and `n_runs` runs.
    pub fn new(n_precursors: usize, n_runs: usize, donor_q: f64, donor_pep: Option<f64>) -> Self {
        Self {
            n_runs,
            donor_q,
            donor_pep,
            best: vec![None; n_precursors],
            slot: vec![NO_SLOT; n_precursors],
            donors: Vec::new(),
        }
    }

    /// Whether `p` is a donor PSM: PSM q-value and PEP within the donor
    /// thresholds.
    pub fn is_donor(&self, p: &PsmIn) -> bool {
        is_donor_psm(p, self.donor_q, self.donor_pep)
    }

    /// Adds one PSM; ignored when its id or run is out of range.
    ///
    /// # Panics
    /// When a donor PSM's run is lower than that of an earlier donor PSM of
    /// the same precursor.
    pub fn add(&mut self, p: &PsmIn) {
        let id = p.precursor as usize;
        if id >= self.best.len() || p.run as usize >= self.n_runs {
            return;
        }
        let score = score_of(p.score);
        match &mut self.best[id] {
            Some((d, s)) => {
                debug_assert_eq!(
                    *d, p.is_decoy,
                    "precursor id {id} mixes target and decoy rows"
                );
                if score > *s {
                    *s = score;
                }
            }
            slot @ None => *slot = Some((p.is_decoy, score)),
        }
        if !self.is_donor(p) {
            return;
        }
        let new = Donor::of(p);
        let d = match self.slot[id] {
            NO_SLOT => {
                self.slot[id] = u32::try_from(self.donors.len())
                    .ok()
                    .filter(|&i| i != NO_SLOT)
                    .expect("fewer than 2^32 - 1 donor precursors");
                self.donors.push(DonorRuns::new(new));
                self.donors.last_mut().expect("just pushed")
            }
            i => {
                let d = &mut self.donors[i as usize];
                let run = d.current.run;
                assert!(
                    p.run >= run,
                    "precursor {id}: a donor PSM of run {} after one of run {run}; add the PSMs in run order",
                    p.run
                );
                if p.run > run {
                    d.close_run();
                    d.current = new;
                } else if new.beats(&d.current) {
                    d.current = new;
                }
                d
            }
        };
        if let Some(cv) = p.faims_cv {
            // Voltages equal under `==` count together, under the smallest.
            match d.cvs.iter_mut().find(|c| c.0 == cv) {
                Some(c) => {
                    c.1 += 1;
                    if cv.total_cmp(&c.0).is_lt() {
                        c.0 = cv;
                    }
                }
                None => d.cvs.push((cv, 1)),
            }
        }
    }

    /// The library over the PSMs added: see [`build_library`].
    pub fn finish(self, meta: Vec<PrecursorMeta>, max_q: f64, ap: &AlignParams) -> Library {
        let n_prec = meta.len();
        let n_runs = self.n_runs;
        let rows: Vec<(u32, bool, f64)> = self
            .best
            .iter()
            .enumerate()
            .filter_map(|(id, b)| b.map(|(d, s)| (id as u32, d, s)))
            .collect();
        let qs = precursor_qvalues(&rows);
        drop(rows);

        let Self {
            slot, mut donors, ..
        } = self;
        let mut entries: Vec<LibraryEntry> = Vec::new();
        let mut obs: Vec<Observation> = Vec::new();
        let mut decoy_obs: Vec<Observation> = Vec::new();
        for (id, &i) in slot.iter().enumerate() {
            if i == NO_SLOT {
                continue;
            }
            let Some(q) = qs.get(id).copied().flatten().filter(|&q| q <= max_q) else {
                continue;
            };
            let id = id as u32;
            let d = &mut donors[i as usize];
            d.close_run();
            let identity_decoy = d.first_is_decoy;
            let best_rt = std::mem::take(&mut d.best_rt);
            let donor_runs: Vec<u16> = best_rt.iter().map(|&(r, _)| r).collect();
            for &(run, rt) in &best_rt {
                if rt.is_finite() {
                    let o = Observation {
                        precursor: id,
                        run,
                        rt,
                    };
                    if identity_decoy {
                        decoy_obs.push(o);
                    } else {
                        obs.push(o);
                    }
                }
            }
            // Majority CV over all donor PSMs; smallest value on ties.
            d.cvs.sort_by(|a, b| a.0.total_cmp(&b.0));
            let mut faims_cv: Option<f32> = None;
            let mut best_n = 0u32;
            for &(cv, n) in &d.cvs {
                if n > best_n {
                    best_n = n;
                    faims_cv = Some(cv);
                }
            }
            let b = d.best.expect("an entry has a donor PSM");
            entries.push(LibraryEntry {
                precursor: id,
                identity_decoy,
                q,
                donor_psm: DonorPsm {
                    score: b.score,
                    q: b.q,
                    pep: b.pep,
                },
                consensus_rt: f64::NAN,
                faims_cv,
                donors: donor_runs,
                best_rt,
            });
        }
        drop(donors);
        drop(slot);

        let mut alignment = align(&obs, n_prec, n_runs, ap);
        drop(obs);
        alignment
            .add_observations(&decoy_obs)
            .expect("one observation per identity decoy and run");
        for e in &mut entries {
            e.consensus_rt = alignment
                .consensus
                .get(e.precursor as usize)
                .copied()
                .unwrap_or(f64::NAN);
        }
        entries.sort_by(|a, b| {
            meta[a.precursor as usize]
                .mz
                .total_cmp(&meta[b.precursor as usize].mz)
                .then(a.precursor.cmp(&b.precursor))
        });
        Library {
            entries,
            alignment,
            meta,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn meta(n: usize) -> Vec<PrecursorMeta> {
        (0..n)
            .map(|i| PrecursorMeta {
                key: format!("P{i}/2"),
                mz: 500.0 + i as f64,
                charge: 2,
                envelope: [0.5, 0.3, 0.15, 0.05],
            })
            .collect()
    }

    fn psm(run: u16, precursor: u32, rt: f64, score: f64, is_decoy: bool) -> PsmIn {
        PsmIn {
            run,
            precursor,
            rt,
            score,
            q: 0.0,
            pep: 0.0,
            is_decoy,
            faims_cv: None,
        }
    }

    /// 200 strong targets over two runs plus the supplied extra PSMs.
    fn base(extra: Vec<PsmIn>) -> Vec<PsmIn> {
        let mut v = Vec::new();
        for p in 0..200u32 {
            for r in 0..2u16 {
                v.push(psm(r, p, 10.0 + p as f64, 100.0 + p as f64, false));
            }
        }
        v.extend(extra);
        v
    }

    #[test]
    fn decoys_enter_as_identity_decoys_without_moving_the_alignment() {
        // Decoy 200 outscores every target and is seen in both runs, at RTs
        // that would bend the maps if it were an anchor.
        let psms = base(vec![
            psm(0, 200, 5.0, 1000.0, true),
            psm(1, 200, 900.0, 990.0, true),
        ]);
        let ap = AlignParams::default();
        let lib = build_library(&psms, meta(201), 2, 1.0, 1.0, None, &ap);
        let targets = build_library(&base(vec![]), meta(201), 2, 1.0, 1.0, None, &ap);
        assert_eq!(lib.entries.len(), 201);
        let d = lib.entries.iter().find(|e| e.precursor == 200).unwrap();
        assert!(d.identity_decoy);
        assert_eq!(d.donors, vec![0, 1]);
        assert_eq!(
            d.donor_psm,
            DonorPsm {
                score: 1000.0,
                q: 0.0,
                pep: 0.0
            }
        );
        assert!(d.consensus_rt.is_finite());
        assert!(lib
            .entries
            .iter()
            .all(|e| e.identity_decoy == (e.precursor == 200)));
        // The target entries and their predictions are those of a library
        // built without the decoy (the entry gate q-values aside, which the
        // decoy's score moves).
        for (a, b) in lib
            .entries
            .iter()
            .filter(|e| !e.identity_decoy)
            .zip(&targets.entries)
        {
            assert_eq!(
                LibraryEntry {
                    q: b.q,
                    ..a.clone()
                },
                *b
            );
            for r in 0..2 {
                assert_eq!(
                    lib.alignment.predict(a.precursor, r),
                    targets.alignment.predict(b.precursor, r)
                );
            }
        }
        assert_eq!(lib.alignment.anchors, targets.alignment.anchors);
    }

    #[test]
    fn identity_decoys_follow_the_target_entry_rules() {
        // Decoy 200: no PSM at the donor threshold. Decoy 201: one donor
        // PSM, but its precursor fails the entry gate (a weak best score).
        let weak = |run, precursor, rt, score| PsmIn {
            q: 0.2,
            ..psm(run, precursor, rt, score, true)
        };
        let psms = base(vec![
            weak(0, 200, 20.0, 900.0),
            psm(1, 201, 30.0, -5.0, true),
        ]);
        let lib = build_library(
            &psms,
            meta(202),
            2,
            0.01,
            0.01,
            None,
            &AlignParams::default(),
        );
        assert!(lib.entries.iter().all(|e| !e.identity_decoy));
        let open = build_library(&psms, meta(202), 2, 1.0, 1.0, None, &AlignParams::default());
        let ids: Vec<u32> = open
            .entries
            .iter()
            .filter(|e| e.identity_decoy)
            .map(|e| e.precursor)
            .collect();
        assert_eq!(ids, vec![200, 201]);
    }

    #[test]
    fn donor_psm_is_the_best_scoring_donor_psm() {
        let with = |run, rt, score, q, pep| PsmIn {
            q,
            pep,
            ..psm(run, 0, rt, score, false)
        };
        let mut psms: Vec<PsmIn> = base(vec![])
            .into_iter()
            .filter(|p| p.precursor != 0)
            .collect();
        psms.extend([
            with(0, 10.0, 50.0, 0.001, 0.01),
            with(1, 11.0, 70.0, 0.0005, 0.002),
            with(1, 12.0, 90.0, 0.2, 0.5),
        ]);
        let lib = build_library(
            &psms,
            meta(200),
            2,
            1.0,
            0.01,
            None,
            &AlignParams::default(),
        );
        let e0 = lib.entries.iter().find(|e| e.precursor == 0).unwrap();
        // The PSM scoring 90 is above the donor threshold and does not count.
        assert_eq!(
            e0.donor_psm,
            DonorPsm {
                score: 70.0,
                q: 0.0005,
                pep: 0.002
            }
        );
    }

    #[test]
    fn entries_sorted_by_mz() {
        let mut m = meta(200);
        // Reverse the m/z order relative to ids.
        for (i, x) in m.iter_mut().enumerate() {
            x.mz = 900.0 - i as f64;
        }
        let lib = build_library(
            &base(vec![]),
            m,
            2,
            0.01,
            1.0,
            None,
            &AlignParams::default(),
        );
        assert_eq!(lib.entries.len(), 200);
        let mzs: Vec<f64> = lib
            .entries
            .iter()
            .map(|e| lib.meta[e.precursor as usize].mz)
            .collect();
        assert!(mzs.windows(2).all(|w| w[0] <= w[1]));
        assert_eq!(lib.entries[0].precursor, 199);
    }

    #[test]
    fn faims_cv_is_majority_of_donor_ids() {
        let cv = |run, rt, c: Option<f32>| PsmIn {
            faims_cv: c,
            ..psm(run, 0, rt, 50.0, false)
        };
        let mut psms: Vec<PsmIn> = base(vec![])
            .into_iter()
            .filter(|p| p.precursor != 0)
            .collect();
        psms.extend([
            cv(0, 10.0, Some(-50.0)),
            cv(1, 10.0, Some(-65.0)),
            cv(1, 11.0, Some(-65.0)),
            cv(0, 12.0, None),
        ]);
        let lib = build_library(&psms, meta(200), 2, 1.0, 1.0, None, &AlignParams::default());
        let e0 = lib.entries.iter().find(|e| e.precursor == 0).unwrap();
        assert_eq!(e0.faims_cv, Some(-65.0));
        // A tie goes to the smallest value; no CV gives None.
        let mut psms2: Vec<PsmIn> = psms.iter().copied().filter(|p| p.precursor != 0).collect();
        psms2.extend([cv(0, 10.0, Some(-50.0)), cv(1, 10.0, Some(-65.0))]);
        let lib = build_library(
            &psms2,
            meta(200),
            2,
            1.0,
            1.0,
            None,
            &AlignParams::default(),
        );
        let e0 = lib.entries.iter().find(|e| e.precursor == 0).unwrap();
        assert_eq!(e0.faims_cv, Some(-65.0));
        let e1 = lib.entries.iter().find(|e| e.precursor == 1).unwrap();
        assert_eq!(e1.faims_cv, None);
    }

    #[test]
    fn best_rt_uses_best_scoring_psm_per_run() {
        let mut psms: Vec<PsmIn> = base(vec![])
            .into_iter()
            .filter(|p| p.precursor != 0)
            .collect();
        psms.extend([
            psm(0, 0, 10.0, 40.0, false),
            psm(0, 0, 12.5, 90.0, false),
            psm(2, 0, 30.0, 60.0, false),
            psm(2, 0, 20.0, 60.0, false),
        ]);
        let lib = build_library(&psms, meta(200), 3, 1.0, 1.0, None, &AlignParams::default());
        let e0 = lib.entries.iter().find(|e| e.precursor == 0).unwrap();
        assert_eq!(e0.donors, vec![0, 2]);
        // Later PSM scores higher in run 0; equal scores take the earliest RT.
        assert_eq!(e0.best_rt, vec![(0, 12.5), (2, 20.0)]);
    }

    #[test]
    fn donors_are_runs_with_a_psm_at_the_psm_threshold() {
        let weak = |run, precursor, rt, score| PsmIn {
            q: 0.2,
            ..psm(run, precursor, rt, score, false)
        };
        let mut psms: Vec<PsmIn> = base(vec![])
            .into_iter()
            .filter(|p| p.precursor > 1)
            .collect();
        // Precursor 0: confident in run 0; in run 1 only a PSM above the PSM
        // threshold, which scores best of all and sits at another RT.
        psms.extend([psm(0, 0, 10.0, 90.0, false), weak(1, 0, 55.0, 500.0)]);
        // Precursor 1: only PSMs above the PSM threshold.
        psms.extend([weak(0, 1, 11.0, 400.0), weak(1, 1, 11.0, 400.0)]);
        let lib = build_library(
            &psms,
            meta(200),
            2,
            0.01,
            0.01,
            None,
            &AlignParams::default(),
        );
        let e0 = lib.entries.iter().find(|e| e.precursor == 0).unwrap();
        assert_eq!(e0.donors, vec![0]);
        assert_eq!(e0.best_rt, vec![(0, 10.0)]);
        // The entry gate still ranks by the best score of every PSM.
        let all = build_library(
            &psms,
            meta(200),
            2,
            0.01,
            1.0,
            None,
            &AlignParams::default(),
        );
        let a0 = all.entries.iter().find(|e| e.precursor == 0).unwrap();
        assert_eq!(e0.q, a0.q);
        assert_eq!(a0.donors, vec![0, 1]);
        assert_eq!(a0.best_rt, vec![(0, 10.0), (1, 55.0)]);
        // No donor PSM: no entry.
        assert!(lib.entries.iter().all(|e| e.precursor != 1));
        assert!(all.entries.iter().any(|e| e.precursor == 1));
    }

    #[test]
    fn donors_need_the_pep_cap_when_given() {
        // Precursor 0: PEP 0.5 in run 0, PEP 0.001 in run 1; precursor 1:
        // PEP 0.5 in both runs. Every PSM is at q = 0.
        let with_pep = |run, precursor, rt, pep| PsmIn {
            pep,
            ..psm(run, precursor, rt, 300.0, false)
        };
        let mut psms: Vec<PsmIn> = base(vec![])
            .into_iter()
            .filter(|p| p.precursor > 1)
            .collect();
        psms.extend([
            with_pep(0, 0, 10.0, 0.5),
            with_pep(1, 0, 12.0, 0.001),
            with_pep(0, 1, 20.0, 0.5),
            with_pep(1, 1, 22.0, 0.5),
        ]);
        let ap = AlignParams::default();
        let capped = build_library(&psms, meta(200), 2, 0.01, 0.01, Some(0.01), &ap);
        let e0 = capped.entries.iter().find(|e| e.precursor == 0).unwrap();
        assert_eq!(e0.donors, vec![1]);
        assert!(capped.entries.iter().all(|e| e.precursor != 1));
        let open = build_library(&psms, meta(200), 2, 0.01, 0.01, None, &ap);
        let o0 = open.entries.iter().find(|e| e.precursor == 0).unwrap();
        assert_eq!(o0.donors, vec![0, 1]);
        assert!(open.entries.iter().any(|e| e.precursor == 1));
    }

    #[test]
    fn best_psm_per_precursor_then_tdc() {
        // p0 target scores 10 and 3 -> best 10; p1 decoy 9; p2 target 8; p3 target 1
        let rows = [
            (0, false, 10.0),
            (0, false, 3.0),
            (1, true, 9.0),
            (2, false, 8.0),
            (3, false, 1.0),
        ];
        let q = precursor_qvalues(&rows);
        // raw (D+1)/T down the list: p0 1/1, p1 2/1 -> 1 (clamp), p2 2/2, p3 2/3;
        // the running minimum from the bottom makes every one 2/3
        assert_eq!(q.len(), 4);
        for (i, qi) in q.iter().enumerate() {
            assert!((qi.unwrap() - 2.0 / 3.0).abs() < 1e-12, "{i}: {qi:?}");
        }
    }

    #[test]
    fn many_targets_reach_low_q() {
        let mut rows: Vec<(u32, bool, f64)> =
            (0..1000).map(|i| (i, false, 1000.0 - i as f64)).collect();
        rows.push((1000, true, 0.5));
        let q = precursor_qvalues(&rows);
        assert!(q[0].unwrap() <= 0.01 && q[999].unwrap() <= 0.01);
    }

    #[test]
    fn missing_ids_are_none_and_empty_input_is_empty() {
        assert!(precursor_qvalues(&[]).is_empty());
        let q = precursor_qvalues(&[(2, false, 1.0)]);
        assert_eq!(q.len(), 3);
        assert!(q[0].is_none() && q[1].is_none() && q[2].is_some());
    }

    #[test]
    fn nan_scores_rank_worst_and_do_not_panic() {
        // 100 strong targets, a decoy at 1.0, and a NaN-scored target. If NaN
        // were ranked best the NaN target would get the low q of the strong
        // targets; ranked worst it sits below the decoy.
        let mut rows: Vec<(u32, bool, f64)> =
            (0..100).map(|i| (i, false, 100.0 + i as f64)).collect();
        rows.push((100, true, 1.0));
        rows.push((101, false, f64::NAN));
        let q = precursor_qvalues(&rows);
        let (strong, decoy, nan) = (q[0].unwrap(), q[100].unwrap(), q[101].unwrap());
        assert!(nan >= decoy, "{nan} < {decoy}");
        assert!(nan > strong, "{nan} <= {strong}");
        // A NaN in one PSM does not hide a finite score of the same precursor.
        let q = precursor_qvalues(&[(0, false, f64::NAN), (0, false, 4.0)]);
        assert_eq!(q[0], Some(1.0));
    }

    #[test]
    fn target_decoy_tie_is_conservative() {
        // 9 strong targets, then a target and a decoy tied at 5.0. Decoy-first
        // gives the tied target (D+1)/T = 2/10; target-first would give 1/10.
        let strong: Vec<(u32, bool, f64)> = (0..9).map(|i| (i, false, 100.0 + i as f64)).collect();
        let mut a = strong.clone();
        a.extend([(9, false, 5.0), (10, true, 5.0)]);
        let mut b = strong;
        b.extend([(9, true, 5.0), (10, false, 5.0)]);
        let ta = precursor_qvalues(&a)[9].unwrap();
        let tb = precursor_qvalues(&b)[10].unwrap();
        assert!((ta - 0.2).abs() < 1e-12, "{ta}");
        assert_eq!(ta, tb);
    }

    #[test]
    #[should_panic(expected = "add the PSMs in run order")]
    fn builder_needs_donor_psms_in_run_order() {
        let mut b = LibraryBuilder::new(1, 3, 1.0, None);
        b.add(&psm(2, 0, 10.0, 5.0, false));
        b.add(&psm(1, 0, 10.0, 5.0, false));
    }

    #[test]
    fn identity_label_is_the_lowest_runs_best_donor() {
        // Precursor 0 mixes labels (a violated precondition, so release
        // builds only): its label is that of run 0's best donor PSM.
        if cfg!(debug_assertions) {
            return;
        }
        let mut b = LibraryBuilder::new(1, 2, 1.0, None);
        b.add(&psm(0, 0, 10.0, 1.0, false));
        b.add(&psm(0, 0, 11.0, 9.0, true));
        b.add(&psm(1, 0, 12.0, 5.0, false));
        let lib = b.finish(meta(1), 1.0, &AlignParams::default());
        assert!(lib.entries[0].identity_decoy);
        assert_eq!(lib.entries[0].best_rt, vec![(0, 11.0), (1, 12.0)]);
    }

    /// The pre-streaming `build_library`, kept to check [`LibraryBuilder`].
    fn reference_library(
        psms: &[PsmIn],
        meta: Vec<PrecursorMeta>,
        n_runs: usize,
        max_q: f64,
        donor_q: f64,
        donor_pep: Option<f64>,
        ap: &AlignParams,
    ) -> Library {
        let n_prec = meta.len();
        let valid: Vec<&PsmIn> = psms
            .iter()
            .filter(|p| (p.precursor as usize) < n_prec && (p.run as usize) < n_runs)
            .collect();
        let rows: Vec<(u32, bool, f64)> = valid
            .iter()
            .map(|p| (p.precursor, p.is_decoy, p.score))
            .collect();
        let qs = precursor_qvalues(&rows);

        // Group the donor PSMs of confident precursors by precursor id.
        let mut by_prec: Vec<&PsmIn> = valid
            .iter()
            .copied()
            .filter(|p| {
                p.q <= donor_q
                    && donor_pep.is_none_or(|c| p.pep <= c)
                    && qs
                        .get(p.precursor as usize)
                        .copied()
                        .flatten()
                        .is_some_and(|q| q <= max_q)
            })
            .collect();
        by_prec.sort_by_key(|p| (p.precursor, p.run));

        let mut entries: Vec<LibraryEntry> = Vec::new();
        let mut obs: Vec<Observation> = Vec::new();
        let mut decoy_obs: Vec<Observation> = Vec::new();
        // NaN scores rank worst, as in the confidence gate.
        let score_of = |p: &PsmIn| {
            if p.score.is_nan() {
                f64::NEG_INFINITY
            } else {
                p.score
            }
        };
        let mut i = 0;
        while i < by_prec.len() {
            let id = by_prec[i].precursor;
            let mut j = i;
            while j < by_prec.len() && by_prec[j].precursor == id {
                j += 1;
            }
            let group = &by_prec[i..j];
            i = j;
            let identity_decoy = group[0].is_decoy;
            debug_assert!(
                group.iter().all(|p| p.is_decoy == identity_decoy),
                "precursor {id} mixes target and decoy PSMs"
            );
            let mut best_psm: Option<&PsmIn> = None;

            let mut donors: Vec<u16> = Vec::new();
            let mut best_rt: Vec<(u16, f64)> = Vec::new();
            let mut k = 0;
            while k < group.len() {
                let run = group[k].run;
                let mut best: Option<&PsmIn> = None;
                while k < group.len() && group[k].run == run {
                    let c = group[k];
                    let better = match best {
                        None => true,
                        Some(b) => {
                            score_of(c) > score_of(b) || (score_of(c) == score_of(b) && c.rt < b.rt)
                        }
                    };
                    if better {
                        best = Some(c);
                    }
                    k += 1;
                }
                let b = best.expect("run group is non-empty");
                // Runs ascend, so a tie on score and RT keeps the lowest run.
                if best_psm.is_none_or(|a| {
                    score_of(b) > score_of(a) || (score_of(b) == score_of(a) && b.rt < a.rt)
                }) {
                    best_psm = Some(b);
                }
                donors.push(run);
                best_rt.push((run, b.rt));
                if b.rt.is_finite() {
                    let o = Observation {
                        precursor: id,
                        run,
                        rt: b.rt,
                    };
                    if identity_decoy {
                        decoy_obs.push(o);
                    } else {
                        obs.push(o);
                    }
                }
            }

            // Majority CV over all PSMs of the precursor; smallest value on ties.
            let mut cvs: Vec<f32> = group.iter().filter_map(|p| p.faims_cv).collect();
            cvs.sort_by(f32::total_cmp);
            let mut faims_cv: Option<f32> = None;
            let mut best_n = 0usize;
            let mut a = 0;
            while a < cvs.len() {
                let mut b = a;
                while b < cvs.len() && cvs[b] == cvs[a] {
                    b += 1;
                }
                if b - a > best_n {
                    best_n = b - a;
                    faims_cv = Some(cvs[a]);
                }
                a = b;
            }

            entries.push(LibraryEntry {
                precursor: id,
                identity_decoy,
                q: qs[id as usize].expect("entry has a q-value"),
                donor_psm: best_psm
                    .map(|b| DonorPsm {
                        score: b.score,
                        q: b.q,
                        pep: b.pep,
                    })
                    .expect("an entry has a donor PSM"),
                consensus_rt: f64::NAN,
                faims_cv,
                donors,
                best_rt,
            });
        }

        let mut alignment = align(&obs, n_prec, n_runs, ap);
        alignment
            .add_observations(&decoy_obs)
            .expect("one observation per identity decoy and run");
        for e in &mut entries {
            e.consensus_rt = alignment
                .consensus
                .get(e.precursor as usize)
                .copied()
                .unwrap_or(f64::NAN);
        }
        entries.sort_by(|a, b| {
            meta[a.precursor as usize]
                .mz
                .total_cmp(&meta[b.precursor as usize].mz)
                .then(a.precursor.cmp(&b.precursor))
        });
        Library {
            entries,
            alignment,
            meta,
        }
    }

    #[test]
    fn streaming_builder_matches_the_slice_build() {
        // Deterministic pseudo-random PSMs over 5 runs: targets and decoys,
        // tied scores and RTs, NaN scores, FAIMS voltages and PSMs above the
        // donor thresholds.
        let mut x = 0x1234_5678_u64;
        let mut next = || {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (x >> 33) as u32
        };
        let n_prec = 400u32;
        let mut psms = Vec::new();
        for _ in 0..6000 {
            let precursor = next() % n_prec;
            let run = (next() % 5) as u16;
            let score = match next() % 50 {
                0 => f64::NAN,
                k => (k % 25) as f64 + (precursor % 7) as f64,
            };
            psms.push(PsmIn {
                run,
                precursor,
                rt: 100.0 + precursor as f64 + (next() % 3) as f64 + run as f64,
                score,
                q: (next() % 4) as f64 * 0.01,
                pep: (next() % 5) as f64 * 0.05,
                is_decoy: precursor % 9 == 0,
                faims_cv: match next() % 4 {
                    0 => None,
                    1 => Some(-0.0),
                    2 => Some(0.0),
                    _ => Some(-45.0),
                },
            });
        }
        // Interleaved runs (build_library adds them in run order), and runs
        // in order with the PSMs of one run in slice order (the stream order).
        let interleaved = psms.clone();
        psms.sort_by_key(|p| p.run);
        let ap = AlignParams {
            min_anchors: 5,
            ..AlignParams::default()
        };
        for (max_q, donor_q, pep, psms) in [
            (1.0, 1.0, None, &psms),
            (0.3, 0.02, Some(0.1), &psms),
            (0.05, 0.01, None, &psms),
            (1.0, 1.0, None, &interleaved),
            (0.3, 0.02, Some(0.1), &interleaved),
        ] {
            let a = reference_library(psms, meta(n_prec as usize), 5, max_q, donor_q, pep, &ap);
            let b = build_library(psms, meta(n_prec as usize), 5, max_q, donor_q, pep, &ap);
            if max_q == 1.0 {
                assert!(!a.entries.is_empty());
            }
            assert_eq!(format!("{:?}", a.entries), format!("{:?}", b.entries));
            assert_eq!(a.alignment.anchors, b.alignment.anchors);
            for p in 0..n_prec {
                for r in 0..5u16 {
                    assert_eq!(
                        format!("{:?}", a.alignment.predict(p, r)),
                        format!("{:?}", b.alignment.predict(p, r))
                    );
                }
            }
        }
    }
}
