use crate::align::{align, AlignParams, Alignment, Observation};
use model::stats::{median_f32_mut, Even};

pub fn precursor_qvalues(rows: &[(u32, bool, f64)]) -> Vec<Option<f64>> {
    let Some(max_id) = rows.iter().map(|r| r.0).max() else {
        return Vec::new();
    };
    let n = max_id as usize + 1;
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

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PsmIn {
    pub run: u16,
    pub precursor: u32,
    pub rt: f64,
    pub score: f64,
    pub q: f64,
    pub q_gate: f64,
    pub pep: f64,
    pub is_decoy: bool,
    pub faims_cv: Option<f32>,
    pub ion_mobility: Option<f32>,
}

impl PsmIn {
    pub fn gate_q(&self) -> f64 {
        if self.q_gate.is_nan() {
            self.q
        } else {
            self.q_gate
        }
    }
}

pub fn is_donor_psm(p: &PsmIn, donor_q: f64, donor_pep: Option<f64>) -> bool {
    p.q <= donor_q && donor_pep.is_none_or(|c| p.pep <= c)
}

#[derive(Debug, Clone, PartialEq)]
pub struct PrecursorMeta {
    pub key: String,
    pub mz: f64,
    pub charge: u8,
    pub envelope: [f32; 4],
    pub isotope_share: f32,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DonorPsm {
    pub score: f64,
    pub q: f64,
    pub pep: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct LibraryEntry {
    pub precursor: u32,
    pub identity_decoy: bool,
    pub q: f64,
    pub donor_psm: DonorPsm,
    pub consensus_rt: f64,
    pub faims_cv: Option<f32>,
    pub ion_mobility: Option<f32>,
    pub donors: Vec<u16>,
    pub best_rt: Vec<(u16, f64)>,
}

#[derive(Debug, Clone)]
pub struct Library {
    pub entries: Vec<LibraryEntry>,
    pub alignment: Alignment,
    pub meta: Vec<PrecursorMeta>,
}

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

fn score_of(score: f64) -> f64 {
    if score.is_nan() {
        f64::NEG_INFINITY
    } else {
        score
    }
}

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
            q: p.gate_q(),
            pep: p.pep,
            run: p.run,
            is_decoy: p.is_decoy,
        }
    }

    fn beats(&self, best: &Donor) -> bool {
        let (a, b) = (score_of(self.score), score_of(best.score));
        a > b || (a == b && self.rt < best.rt)
    }
}

#[derive(Debug, Clone)]
struct DonorRuns {
    current: Donor,
    best: Option<Donor>,
    first_is_decoy: bool,
    best_rt: Vec<(u16, f64)>,
    cvs: Vec<(f32, u32)>,
    k0s: Vec<f32>,
}

impl DonorRuns {
    fn new(d: Donor) -> Self {
        Self {
            current: d,
            best: None,
            first_is_decoy: false,
            best_rt: Vec::new(),
            cvs: Vec::new(),
            k0s: Vec::new(),
        }
    }

    fn close_run(&mut self) {
        let c = self.current;
        if self.best_rt.is_empty() {
            self.first_is_decoy = c.is_decoy;
        }
        self.best_rt.push((c.run, c.rt));
        if self.best.is_none_or(|b| c.beats(&b)) {
            self.best = Some(c);
        }
    }
}

const NO_SLOT: u32 = u32::MAX;

#[derive(Debug, Clone)]
pub struct LibraryBuilder {
    n_runs: usize,
    donor_q: f64,
    donor_pep: Option<f64>,
    best: Vec<Option<(bool, f64)>>,
    slot: Vec<u32>,
    donors: Vec<DonorRuns>,
}

impl LibraryBuilder {
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

    pub fn is_donor(&self, p: &PsmIn) -> bool {
        is_donor_psm(p, self.donor_q, self.donor_pep)
    }

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
        if let Some(k0) = p.ion_mobility.filter(|v| v.is_finite()) {
            d.k0s.push(k0);
        }
        if let Some(cv) = p.faims_cv {
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
            d.cvs.sort_by(|a, b| a.0.total_cmp(&b.0));
            let mut faims_cv: Option<f32> = None;
            let mut best_n = 0u32;
            for &(cv, n) in &d.cvs {
                if n > best_n {
                    best_n = n;
                    faims_cv = Some(cv);
                }
            }
            let ion_mobility = median_f32_mut(&mut d.k0s, Even::Mean);
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
                ion_mobility,
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
                isotope_share: 1.0,
            })
            .collect()
    }

    fn psm(run: u16, precursor: u32, rt: f64, score: f64, is_decoy: bool) -> PsmIn {
        PsmIn {
            q_gate: f64::NAN,
            run,
            precursor,
            rt,
            score,
            q: 0.0,
            pep: 0.0,
            is_decoy,
            faims_cv: None,
            ion_mobility: None,
        }
    }

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
        let weak = |run, precursor, rt, score| PsmIn {
            q_gate: f64::NAN,
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
            q_gate: f64::NAN,
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
    fn best_donor_is_gated_by_its_pooled_q_while_donors_follow_the_run_q() {
        let with = |run, rt, score, q, q_gate| PsmIn {
            q_gate,
            q,
            pep: 0.001,
            ..psm(run, 0, rt, score, false)
        };
        let mut psms: Vec<PsmIn> = base(vec![])
            .into_iter()
            .filter(|p| p.precursor != 0)
            .collect();
        psms.extend([
            with(0, 10.0, 50.0, 0.004, 0.0002),
            with(1, 11.0, 70.0, 0.0005, 0.03),
            with(1, 12.0, 90.0, 0.2, 0.0001),
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
        assert_eq!(e0.donors, vec![0, 1]);
        assert_eq!(e0.donor_psm.score, 70.0);
        assert_eq!(e0.donor_psm.q, 0.03);
        let psms_nan: Vec<PsmIn> = psms
            .iter()
            .map(|p| PsmIn {
                q_gate: f64::NAN,
                ..*p
            })
            .collect();
        let lib = build_library(
            &psms_nan,
            meta(200),
            2,
            1.0,
            0.01,
            None,
            &AlignParams::default(),
        );
        let e0 = lib.entries.iter().find(|e| e.precursor == 0).unwrap();
        assert_eq!(e0.donor_psm.q, 0.0005);
    }

    #[test]
    fn an_identity_decoy_is_gated_like_a_target() {
        let decoy = |run, rt, score, q, q_gate| PsmIn {
            q,
            q_gate,
            ..psm(run, 200, rt, score, true)
        };
        let pooled = base(vec![
            decoy(0, 5.0, 1000.0, 0.3, 0.0002),
            decoy(1, 900.0, 990.0, 0.004, 0.02),
        ]);
        let ap = AlignParams::default();
        let lib = build_library(&pooled, meta(201), 2, 1.0, 0.01, None, &ap);
        let d = lib.entries.iter().find(|e| e.precursor == 200).unwrap();
        assert!(d.identity_decoy);
        assert_eq!(d.donors, vec![1]);
        assert_eq!(d.donor_psm.score, 990.0);
        assert_eq!(d.donor_psm.q, 0.02);
        let plain: Vec<PsmIn> = pooled
            .iter()
            .map(|p| PsmIn {
                q_gate: f64::NAN,
                ..*p
            })
            .collect();
        let lib = build_library(&plain, meta(201), 2, 1.0, 0.01, None, &ap);
        let d = lib.entries.iter().find(|e| e.precursor == 200).unwrap();
        assert_eq!(d.donor_psm.q, 0.004);
    }

    #[test]
    fn entries_sorted_by_mz() {
        let mut m = meta(200);
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
            q_gate: f64::NAN,
            faims_cv: c,
            ion_mobility: None,
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
        assert_eq!(e0.best_rt, vec![(0, 12.5), (2, 20.0)]);
    }

    #[test]
    fn donors_are_runs_with_a_psm_at_the_psm_threshold() {
        let weak = |run, precursor, rt, score| PsmIn {
            q_gate: f64::NAN,
            q: 0.2,
            ..psm(run, precursor, rt, score, false)
        };
        let mut psms: Vec<PsmIn> = base(vec![])
            .into_iter()
            .filter(|p| p.precursor > 1)
            .collect();
        psms.extend([psm(0, 0, 10.0, 90.0, false), weak(1, 0, 55.0, 500.0)]);
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
        assert!(lib.entries.iter().all(|e| e.precursor != 1));
        assert!(all.entries.iter().any(|e| e.precursor == 1));
    }

    #[test]
    fn donors_need_the_pep_cap_when_given() {
        let with_pep = |run, precursor, rt, pep| PsmIn {
            q_gate: f64::NAN,
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
        let rows = [
            (0, false, 10.0),
            (0, false, 3.0),
            (1, true, 9.0),
            (2, false, 8.0),
            (3, false, 1.0),
        ];
        let q = precursor_qvalues(&rows);
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
        let mut rows: Vec<(u32, bool, f64)> =
            (0..100).map(|i| (i, false, 100.0 + i as f64)).collect();
        rows.push((100, true, 1.0));
        rows.push((101, false, f64::NAN));
        let q = precursor_qvalues(&rows);
        let (strong, decoy, nan) = (q[0].unwrap(), q[100].unwrap(), q[101].unwrap());
        assert!(nan >= decoy, "{nan} < {decoy}");
        assert!(nan > strong, "{nan} <= {strong}");
        let q = precursor_qvalues(&[(0, false, f64::NAN), (0, false, 4.0)]);
        assert_eq!(q[0], Some(1.0));
    }

    #[test]
    fn target_decoy_tie_is_conservative() {
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

            let mut k0s: Vec<f32> = group
                .iter()
                .filter_map(|p| p.ion_mobility.filter(|v| v.is_finite()))
                .collect();
            let ion_mobility = median_f32_mut(&mut k0s, Even::Mean);
            entries.push(LibraryEntry {
                precursor: id,
                identity_decoy,
                q: qs[id as usize].expect("entry has a q-value"),
                donor_psm: best_psm
                    .map(|b| DonorPsm {
                        score: b.score,
                        q: b.gate_q(),
                        pep: b.pep,
                    })
                    .expect("an entry has a donor PSM"),
                consensus_rt: f64::NAN,
                faims_cv,
                ion_mobility,
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
                q_gate: f64::NAN,
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
                ion_mobility: None,
            });
        }
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
