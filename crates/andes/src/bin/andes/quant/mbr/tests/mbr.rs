use std::path::PathBuf;

use model::mass::ISOTOPE;
use model::scan::Ms1Scan;
use output::extract::{read_extract_table, write_extract_table};
use output::ms1::write_ms1_table;
use quant::align::AlignParams;
use quant::decoys::Hypothesis;
use quant::library::{build_library, PrecursorMeta, PsmIn};
use quant::ms1_index::Ms1RunIndex;

use super::*;
use crate::quant::finalize::finalize;

const ENVELOPE: [f32; 4] = [0.55, 0.30, 0.11, 0.04];
const CHARGE: u8 = 2;
const WIDTH: f64 = 4.2;
const AMPLITUDE: f64 = 1e6;
const BOTH: u32 = 30;
const ONLY_RUN0: u32 = 30;
const ONLY_RUN1: u32 = 31;
const N_PREC: u32 = 32;
const N_RUNS: usize = 3;
const FAIMS_CVS: [f32; 2] = [-45.0, -60.0];

fn mz_of(p: u32) -> f64 {
    450.0 + 13.7 * p as f64
}

fn apex(p: u32, run: u16) -> f64 {
    200.0 + 30.0 * p as f64 + 4.0 * run as f64
}

struct Lcg(u64);
impl Lcg {
    fn uniform(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
}

fn run_scans(run: u16, faims: bool, silent: &dyn Fn(u32) -> bool) -> Vec<Ms1Scan> {
    run_scans_shaped(run, faims, silent, WIDTH, 1.0, 0.0)
}

const ECHO_AFTER: f64 = 19.0;

fn run_scans_shaped(
    run: u16,
    faims: bool,
    silent: &dyn Fn(u32) -> bool,
    width: f64,
    step: f64,
    echo: f64,
) -> Vec<Ms1Scan> {
    let mut rng = Lcg(17 + run as u64);
    (0..(1300.0 / step) as usize)
        .map(|k| {
            let t = k as f64 * step;
            let cv = faims.then_some(FAIMS_CVS[k % 2]);
            let mut peaks = Vec::new();
            for _ in 0..40 {
                peaks.push((
                    300.0 + 1200.0 * rng.uniform(),
                    (500.0 + 1000.0 * rng.uniform()) as f32,
                ));
            }
            for p in 0..N_PREC {
                let elutes = !silent(p) && cv.is_none_or(|cv| cv == FAIMS_CVS[p as usize % 2]);
                let profile = |c: f64| (-(t - c).powi(2) / (2.0 * width * width)).exp();
                let g = profile(apex(p, run)) + echo * profile(apex(p, run) + ECHO_AFTER);
                for (i, &e) in ENVELOPE.iter().enumerate() {
                    let signal = if elutes {
                        AMPLITUDE * e as f64 * g
                    } else {
                        0.0
                    };
                    let inten = 200.0 + signal + 50.0 * (rng.uniform() - 0.5);
                    peaks.push((mz_of(p) + i as f64 * ISOTOPE / CHARGE as f64, inten as f32));
                }
            }
            Ms1Scan {
                rt: t,
                peaks,
                faims_cv: cv,
                ..Ms1Scan::default()
            }
        })
        .collect()
}

struct Fixture {
    _dir: tempfile::TempDir,
    lib: Library,
    ms1_tables: Vec<PathBuf>,
}

fn fixture(faims: bool) -> Fixture {
    fixture_with(faims, &|_| false, &|_| false)
}

const DONOR_Q: f64 = 0.01;

fn fixture_with(
    faims: bool,
    silent_in_run1: &dyn Fn(u32) -> bool,
    weak_in_run1: &dyn Fn(u32) -> bool,
) -> Fixture {
    fixture_decoys(faims, silent_in_run1, weak_in_run1, &|_| false)
}

fn fixture_decoys(
    faims: bool,
    silent_in_run1: &dyn Fn(u32) -> bool,
    weak_in_run1: &dyn Fn(u32) -> bool,
    decoy: &dyn Fn(u32) -> bool,
) -> Fixture {
    fixture_psm_rt(faims, silent_in_run1, weak_in_run1, decoy, &|p, run| {
        apex(p, run) + 2.0
    })
}

fn fixture_psm_rt(
    faims: bool,
    silent_in_run1: &dyn Fn(u32) -> bool,
    weak_in_run1: &dyn Fn(u32) -> bool,
    decoy: &dyn Fn(u32) -> bool,
    psm_rt: &dyn Fn(u32, u16) -> f64,
) -> Fixture {
    let dir = tempfile::tempdir().unwrap();
    let meta: Vec<PrecursorMeta> = (0..N_PREC)
        .map(|p| PrecursorMeta {
            key: format!("PEP{p}/2"),
            mz: mz_of(p),
            charge: CHARGE,
            envelope: ENVELOPE,
            isotope_share: 1.0,
        })
        .collect();
    let mut psms = Vec::new();
    let mut psm = |p: u32, run: u16| {
        psms.push(PsmIn {
            q_gate: f64::NAN,
            run,
            precursor: p,
            rt: psm_rt(p, run),
            score: 10.0,
            q: if run == 1 && weak_in_run1(p) {
                0.2
            } else {
                0.001
            },
            pep: 0.0,
            is_decoy: decoy(p),
            faims_cv: (faims && run == 1).then_some(FAIMS_CVS[p as usize % 2]),
            ion_mobility: None,
        })
    };
    for p in 0..BOTH {
        psm(p, 0);
        psm(p, 1);
    }
    psm(ONLY_RUN0, 0);
    psm(ONLY_RUN1, 1);
    psm(0, 2);
    let ap = AlignParams {
        min_anchors: 10,
        ..AlignParams::default()
    };
    let lib = build_library(&psms, meta, N_RUNS, 1.0, DONOR_Q, None, &ap);
    let ms1_tables = (0..N_RUNS as u16)
        .map(|r| {
            let path = dir.path().join(format!("run{r}.ms1.parquet"));
            let silent = |p: u32| r == 1 && silent_in_run1(p);
            let index = Ms1RunIndex::new(run_scans(r, faims && r == 1, &silent));
            write_ms1_table(&path, &format!("run{r}"), &index).unwrap();
            path
        })
        .collect();
    Fixture {
        _dir: dir,
        lib,
        ms1_tables,
    }
}

fn params() -> EvidenceParams {
    EvidenceParams::default()
}

fn extract_all(
    run: u16,
    ms1_table: &Path,
    lib: &Library,
    ep: &EvidenceParams,
) -> Result<RunExtract, String> {
    extract_run(run, ms1_table, lib, ep, &TransferSelection::OFF)
}

fn rows_of(rows: &[ExtractRow], p: u32) -> Vec<&ExtractRow> {
    rows.iter().filter(|r| r.precursor == p).collect()
}

fn is_identified_target(r: &ExtractRow) -> bool {
    (r.role, r.kind) == (Role::Identified, Hypothesis::Target)
}

#[test]
fn identified_pairs_get_identified_and_pseudo_rows() {
    let fx = fixture(false);
    assert!(fx.lib.alignment.maps[1].is_some());
    let out = extract_all(1, &fx.ms1_tables[1], &fx.lib, &params()).unwrap();
    for p in 0..BOTH {
        let rows = rows_of(&out.rows, p);
        assert_eq!(rows.len(), 4, "precursor {p}: {rows:?}");
        let (id, pseudo) = (rows[0], rows[3]);
        assert_eq!((id.role, id.kind), (Role::Identified, Hypothesis::Target));
        assert_eq!(rows[1].role, Role::Identified);
        assert!(matches!(rows[1].kind, Hypothesis::MassShiftDecoy(_)));
        assert_eq!(
            (rows[2].role, rows[2].kind),
            (Role::Identified, Hypothesis::RtShiftDecoy)
        );
        for d in &rows[1..3] {
            assert_eq!(
                (d.n_donor_runs, d.library_q),
                (id.n_donor_runs, id.library_q)
            );
            let e = d.evidence.as_ref().unwrap();
            let target_llr = id.evidence.as_ref().unwrap().llr;
            assert!(e.llr < target_llr / 10.0, "precursor {p}: {d:?}");
        }
        assert_eq!(
            (pseudo.role, pseudo.kind),
            (Role::PseudoTransfer, Hypothesis::Target)
        );
        assert_eq!((id.n_donor_runs, pseudo.n_donor_runs), (2, 1));
        assert_eq!((id.run, pseudo.run), (1, 1));
        for r in [id, pseudo] {
            let e = r.evidence.as_ref().unwrap();
            assert!(e.llr > 100.0, "precursor {p}: {e:?}");
            assert!((e.t0 - apex(p, 1)).abs() <= 1.0, "precursor {p}: {e:?}");
            assert!((e.amplitude / AMPLITUDE - 1.0).abs() < 0.1, "{e:?}");
        }
    }
    let rows = rows_of(&out.rows, ONLY_RUN1);
    assert_eq!(rows.len(), 3);
    assert!(rows.iter().all(|r| r.role == Role::Identified));
    assert_eq!(out.counts.no_pseudo_prediction, 1);
    assert!(!out.counts.unaligned);
}

#[test]
fn identified_decoys_are_the_transfer_probes_at_the_identified_prior() {
    let fx = fixture(false);
    let ms1 = Ms1RunIndex::new(run_scans(1, false, &|_| false));
    let scales = RunScales {
        spacing_s: ms1.median_scan_spacing(),
        width_s: WIDTH,
        ..RunScales::default()
    };
    let wp = WidthPrior::from_fits(Vec::new(), &scales);
    let entry = fx.lib.entries.iter().find(|e| e.precursor == 4).unwrap();
    let off = PairGate::new(&TransferSelection::OFF, &fx.lib);
    let (rows, _, _) = entry_rows(entry, 1, true, &fx.lib, &off, &ms1, &wp, &rt_pool(&fx.lib, 1));
    let id: Vec<&(ExtractRow, Probe)> = rows
        .iter()
        .filter(|(r, _)| r.role == Role::Identified)
        .collect();
    assert_eq!(id.len(), 3, "{rows:?}");
    let rt = identified_rt(entry, 1).unwrap();
    let w = wp.at(rt);
    let want = probes(
        4,
        1,
        mz_of(4),
        CHARGE,
        rt,
        w.max(1.0),
        w,
        ms1.rt_range().unwrap(),
    );
    let got: Vec<Probe> = id.iter().map(|(_, p)| *p).collect();
    assert_eq!(got, want);
    let target = got[0];
    assert_eq!(target.kind, Hypothesis::Target);
    assert_eq!((target.t_hat, target.mz), (rt, mz_of(4)));
    for (row, probe) in &id {
        assert_eq!(row.kind, probe.kind);
        assert_eq!((probe.s, probe.w_prior), (target.s, target.w_prior));
    }
    let mass = got[1];
    let Hypothesis::MassShiftDecoy(n) = mass.kind else {
        panic!("{mass:?}")
    };
    assert_eq!(mass.t_hat, target.t_hat);
    let shift = n as f64 * quant::decoys::PEPTIDE_MASS_SLOPE / CHARGE as f64;
    assert!((mass.mz - target.mz - shift).abs() < 1e-9);
    let rt_decoy = got[2];
    assert_eq!(rt_decoy.mz, target.mz);
    assert!((rt_decoy.t_hat - rt).abs() >= 8.0 * target.s + 4.0 * target.w_prior);
}

#[test]
fn identity_decoys_are_measured_as_targets() {
    let target = fixture(false);
    let decoy = fixture_decoys(false, &|_| false, &|_| false, &|p| p == ONLY_RUN0);
    let e = decoy
        .lib
        .entries
        .iter()
        .find(|e| e.precursor == ONLY_RUN0)
        .unwrap();
    assert!(e.identity_decoy);
    assert_eq!(e.donors, vec![0]);
    let same = |rows: Vec<&ExtractRow>| -> Vec<ExtractRow> {
        rows.into_iter()
            .map(|r| ExtractRow {
                library_q: 0.0,
                ..r.clone()
            })
            .collect()
    };
    for run in [0u16, 1] {
        let a = extract_all(
            run,
            &target.ms1_tables[run as usize],
            &target.lib,
            &params(),
        )
        .unwrap();
        let b =
            extract_all(run, &decoy.ms1_tables[run as usize], &decoy.lib, &params()).unwrap();
        let rows = rows_of(&b.rows, ONLY_RUN0);
        if run == 0 {
            assert_eq!(rows.len(), 3, "run {run}: {rows:?}");
        } else {
            assert_transfer_layout(&rows);
        }
        let role = if run == 0 {
            Role::Identified
        } else {
            Role::Transfer
        };
        assert!(rows.iter().all(|r| r.role == role));
        if run == 1 {
            let ev = rows[0].evidence.as_ref().unwrap();
            assert!(ev.llr > 100.0, "{ev:?}");
            assert_eq!(same(rows_of(&a.rows, ONLY_RUN0)), same(rows));
            assert_eq!(same(a.rows.iter().collect()), same(b.rows.iter().collect()));
        } else {
            assert_eq!(a.counts.width_fits, b.counts.width_fits + 1);
        }
    }
}

#[test]
fn missing_pairs_get_target_and_decoys() {
    let fx = fixture(false);
    let out = extract_all(1, &fx.ms1_tables[1], &fx.lib, &params()).unwrap();
    let rows = rows_of(&out.rows, ONLY_RUN0);
    assert_transfer_layout(&rows);
    assert!(rows
        .iter()
        .all(|r| r.role == Role::Transfer && r.n_donor_runs == 1));
    let target = rows[0].evidence.as_ref().unwrap();
    assert!((target.t0 - apex(ONLY_RUN0, 1)).abs() <= 1.0, "{target:?}");
    assert!(target.llr > 100.0, "{target:?}");
    for d in &rows[1..] {
        let e = d.evidence.as_ref().unwrap();
        assert!(e.llr < target.llr / 10.0, "{d:?}");
    }
    let out0 = extract_all(0, &fx.ms1_tables[0], &fx.lib, &params()).unwrap();
    let rows = rows_of(&out0.rows, ONLY_RUN1);
    assert_transfer_layout(&rows);
    assert!(rows.iter().all(|r| r.role == Role::Transfer));
    let lq = fx
        .lib
        .entries
        .iter()
        .find(|e| e.precursor == ONLY_RUN1)
        .unwrap()
        .q;
    assert!(rows.iter().all(|r| r.library_q == lq));
}

fn assert_transfer_layout(rows: &[&ExtractRow]) {
    assert!(rows.len() >= 2, "{rows:?}");
    assert!(rows.len() <= 2 + quant::decoys::RANDOM_RT_DECOYS, "{rows:?}");
    assert_eq!(rows[0].kind, Hypothesis::Target);
    assert!(matches!(rows[1].kind, Hypothesis::MassShiftDecoy(_)));
    assert!(
        rows[2..]
            .iter()
            .all(|r| matches!(r.kind, Hypothesis::RandomRtDecoy(_))),
        "{rows:?}"
    );
}

fn with_entry(lib: &Library, p: u32, f: impl Fn(&mut LibraryEntry)) -> Library {
    let mut lib = lib.clone();
    f(lib.entries.iter_mut().find(|e| e.precursor == p).unwrap());
    lib
}

fn without(all: &[ExtractRow], skipped: &[(u32, u16)]) -> Vec<ExtractRow> {
    all.iter()
        .filter(|r| !(r.role == Role::Transfer && skipped.contains(&(r.precursor, r.run))))
        .cloned()
        .collect()
}

const SELECT: TransferSelection = TransferSelection {
    min_donor_fraction: 0.5,
};

#[test]
fn a_marginal_best_donor_psm_does_not_skip_the_transfer() {
    let fx = fixture(false);
    let lib = with_entry(&fx.lib, ONLY_RUN0, |e| e.donor_psm.q = 0.002);
    let all = extract_all(1, &fx.ms1_tables[1], &lib, &params()).unwrap();
    let off = TransferSelection { min_donor_fraction: 0.0 };
    let kept = extract_run(1, &fx.ms1_tables[1], &lib, &params(), &off).unwrap();
    assert!(!rows_of(&kept.rows, ONLY_RUN0).is_empty());
    assert_eq!(kept.rows, all.rows);
}

#[test]
fn selection_skips_transfers_with_too_few_donor_runs() {
    let fx = fixture(false);
    let lib = with_entry(&fx.lib, ONLY_RUN1, |e| e.donors.clear());
    let all = extract_all(0, &fx.ms1_tables[0], &lib, &params()).unwrap();
    assert_transfer_layout(&rows_of(&all.rows, ONLY_RUN1));
    let sel = extract_run(0, &fx.ms1_tables[0], &lib, &params(), &SELECT).unwrap();
    assert_eq!(sel.rows, without(&all.rows, &[(ONLY_RUN1, 0)]));
    assert_eq!(sel.counts.not_selected, 1);
    let roles =
        |rows: &[ExtractRow], role: Role| rows.iter().filter(|r| r.role == role).count();
    for role in [Role::Identified, Role::PseudoTransfer] {
        assert_eq!(roles(&sel.rows, role), roles(&all.rows, role));
    }
    let sel = extract_run(0, &fx.ms1_tables[0], &fx.lib, &params(), &SELECT).unwrap();
    assert_eq!(
        sel,
        extract_all(0, &fx.ms1_tables[0], &fx.lib, &params()).unwrap()
    );
}

#[test]
fn min_donors_scale_with_the_aligned_runs() {
    let sel = |f: f64| TransferSelection {
        min_donor_fraction: f,
    };
    assert_eq!(sel(0.0).min_donors(500), 0);
    assert_eq!(sel(0.004).min_donors(6), 1);
    assert_eq!(sel(0.004).min_donors(250), 1);
    assert_eq!(sel(0.004).min_donors(500), 2);
    assert_eq!(sel(0.004).min_donors(1000), 4);
    assert_eq!(sel(1.0).min_donors(1), 0);
    assert!(TransferSelection::OFF.is_off());
    assert!(TransferSelection::default().is_off());
}

#[test]
fn transfers_stay_inside_an_alignment_component() {
    let dir = tempfile::tempdir().unwrap();
    let n = 122u32;
    let meta: Vec<PrecursorMeta> = (0..n)
        .map(|p| PrecursorMeta {
            key: format!("PEP{p}/2"),
            mz: 400.0 + 7.3 * p as f64,
            charge: CHARGE,
            envelope: ENVELOPE,
            isotope_share: 1.0,
        })
        .collect();
    let mut psms = Vec::new();
    let mut psm = |p: u32, run: u16| {
        psms.push(PsmIn {
            q_gate: f64::NAN,
            run,
            precursor: p,
            rt: 100.0 + 8.0 * p as f64 + if run >= 2 { 120.0 } else { 0.0 },
            score: 10.0,
            q: 0.001,
            pep: 0.0,
            is_decoy: false,
            faims_cv: None,
            ion_mobility: None,
        })
    };
    for p in 0..120 {
        let pair = if p < 60 { [0, 1] } else { [2, 3] };
        pair.into_iter().for_each(|r| psm(p, r));
    }
    psm(120, 0);
    psm(121, 0);
    psm(121, 2);
    let ap = AlignParams {
        min_anchors: 10,
        ..AlignParams::default()
    };
    let lib = build_library(&psms, meta, 4, 1.0, DONOR_Q, None, &ap);
    assert_eq!(lib.alignment.n_components(), 2);
    let path = dir.path().join("run3.ms1.parquet");
    let index = Ms1RunIndex::new(run_scans(3, false, &|_| false));
    write_ms1_table(&path, "run3", &index).unwrap();
    let out = extract_all(3, &path, &lib, &params()).unwrap();
    assert_eq!(out.counts.other_component, 61);
    assert_eq!(out.counts.no_prediction, 0);
    for p in (0..60).chain([120]) {
        assert!(rows_of(&out.rows, p).is_empty(), "{p}");
    }
    let rows = rows_of(&out.rows, 121);
    assert!(!rows.is_empty());
    assert!(rows.iter().all(|r| r.role == Role::Transfer), "{rows:?}");
    assert_eq!(rows[0].kind, Hypothesis::Target);
    let path = dir.path().join("run2.ms1.parquet");
    let index = Ms1RunIndex::new(run_scans(2, false, &|_| false));
    write_ms1_table(&path, "run2", &index).unwrap();
    let out = extract_all(2, &path, &lib, &params()).unwrap();
    assert_eq!(out.counts.pseudo_other_component, 1);
    let rows = rows_of(&out.rows, 121);
    assert!(rows.iter().all(|r| r.role == Role::Identified), "{rows:?}");
    assert_eq!(rows.iter().filter(|r| is_identified_target(r)).count(), 1);
}

#[test]
fn unaligned_run_keeps_its_identified_rows() {
    let fx = fixture(false);
    assert!(fx.lib.alignment.maps[2].is_none());
    let out = extract_all(2, &fx.ms1_tables[2], &fx.lib, &params()).unwrap();
    assert!(out.counts.unaligned);
    assert_eq!(out.rows.len(), 3, "{:?}", out.rows);
    assert!(out.rows.iter().all(|r| r.role == Role::Identified));
    let r = &out.rows[0];
    assert_eq!(
        (r.precursor, r.run, r.role, r.kind),
        (0, 2, Role::Identified, Hypothesis::Target)
    );
    let e = r.evidence.as_ref().unwrap();
    assert!((e.t0 - apex(0, 2)).abs() <= 1.0, "{e:?}");
    assert!((e.amplitude / AMPLITUDE - 1.0).abs() < 0.1, "{e:?}");
    assert_eq!(
        (out.counts.no_prediction, out.counts.no_pseudo_prediction),
        (0, 0)
    );
    assert!(extract_all(N_RUNS as u16, &fx.ms1_tables[2], &fx.lib, &params()).is_err());

    let runs: Vec<Vec<ExtractRow>> = (0..N_RUNS)
        .map(|r| {
            extract_all(r as u16, &fx.ms1_tables[r], &fx.lib, &params())
                .unwrap()
                .rows
        })
        .collect();
    let stems: Vec<String> = ["r0", "r1", "r2"].map(String::from).to_vec();
    let (mbr, val) = mbr_paths(&fx._dir.path().join("x"));
    finalize(&stems, &fx.lib, &runs, 0.01, None, 7, &mbr, &val).unwrap();
    let rows = output::mbr::read_mbr_parquet(&mbr).unwrap();
    let in_r2: Vec<_> = rows.iter().filter(|r| r.run == "r2").collect();
    assert_eq!(in_r2.len(), 1, "{rows:?}");
    assert_eq!(
        (
            in_r2[0].precursor_key.as_str(),
            in_r2[0].evidence_type.as_str()
        ),
        ("PEP0/2", "identified"),
        "{in_r2:?}"
    );
    assert!(
        in_r2[0].quantity.is_some() && in_r2[0].quantity_se.is_some(),
        "{in_r2:?}"
    );
}

#[test]
fn psms_above_the_donor_threshold_make_a_transfer() {
    let fx = fixture_with(false, &|_| false, &|p| p == 7);
    let e = fx.lib.entries.iter().find(|e| e.precursor == 7).unwrap();
    assert_eq!(e.donors, vec![0]);
    let out = extract_all(1, &fx.ms1_tables[1], &fx.lib, &params()).unwrap();
    let rows = rows_of(&out.rows, 7);
    assert_transfer_layout(&rows);
    assert!(rows.iter().all(|r| r.role == Role::Transfer), "{rows:?}");
    assert_eq!(rows[0].kind, Hypothesis::Target);
    let target = rows[0].evidence.as_ref().unwrap();
    assert!((target.t0 - apex(7, 1)).abs() <= 1.0, "{target:?}");
    let rows = rows_of(&out.rows, 8);
    assert_eq!(rows[0].role, Role::Identified);
}

#[test]
fn foreign_voltage_entries_are_skipped_and_counted() {
    let mut fx = fixture(true);
    for e in fx.lib.entries.iter_mut() {
        if e.precursor == ONLY_RUN0 || e.precursor == 3 {
            e.faims_cv = Some(-30.0);
        }
    }
    let out = extract_all(1, &fx.ms1_tables[1], &fx.lib, &params()).unwrap();
    assert!(rows_of(&out.rows, ONLY_RUN0).is_empty());
    assert!(rows_of(&out.rows, 3).is_empty());
    assert_eq!(out.counts.foreign_cv, 2);
    assert!(!rows_of(&out.rows, 4).is_empty());
    let out0 = extract_all(0, &fx.ms1_tables[0], &fx.lib, &params()).unwrap();
    assert_eq!(out0.counts.foreign_cv, 0);
    assert_eq!(rows_of(&out0.rows, ONLY_RUN0)[0].role, Role::Identified);
    assert_eq!(rows_of(&out0.rows, 3)[0].role, Role::Identified);
}

#[test]
fn output_order_is_deterministic() {
    let fx = fixture(false);
    let a = extract_all(1, &fx.ms1_tables[1], &fx.lib, &params()).unwrap();
    let b = extract_all(1, &fx.ms1_tables[1], &fx.lib, &params()).unwrap();
    assert_eq!(a, b);
    let order: Vec<u32> = fx.lib.entries.iter().map(|e| e.precursor).collect();
    let mut seen: Vec<u32> = a.rows.iter().map(|r| r.precursor).collect();
    seen.dedup();
    let expected: Vec<u32> = order.into_iter().filter(|p| seen.contains(p)).collect();
    assert_eq!(seen, expected);
    assert_eq!(seen.len(), N_PREC as usize);
}

#[test]
fn nonfinite_best_rt_is_skipped_and_counted() {
    let mut fx = fixture(false);
    let e = fx
        .lib
        .entries
        .iter_mut()
        .find(|e| e.precursor == 5)
        .unwrap();
    for (r, rt) in e.best_rt.iter_mut() {
        if *r == 1 {
            *rt = f64::NAN;
        }
    }
    let out = extract_all(1, &fx.ms1_tables[1], &fx.lib, &params()).unwrap();
    assert!(rows_of(&out.rows, 5).is_empty());
    assert_eq!(out.counts.nonfinite_best_rt, 1);
}

#[test]
fn width_prior_follows_identified_fits() {
    let fx = fixture(false);
    let out = extract_all(1, &fx.ms1_tables[1], &fx.lib, &params()).unwrap();
    assert!(out.counts.width_fits >= 20, "{:?}", out.counts);
    assert!(!out.counts.width_prior_fallback);
    for r in out.rows.iter().filter(|r| is_identified_target(r)) {
        let e = r.evidence.as_ref().unwrap();
        let w_prior = e.w / e.width_ratio;
        assert!(
            (w_prior / WIDTH - 1.0).abs() < 0.1,
            "w_prior {w_prior}: {e:?}"
        );
    }
}

#[test]
fn run_scales_come_from_the_run() {
    let fx = fixture(false);
    let out = extract_all(1, &fx.ms1_tables[1], &fx.lib, &params()).unwrap();
    let sc = out.counts.scales;
    assert_eq!(sc.spacing_s, Some(1.0));
    assert!(!sc.fallback, "{sc:?}");
    assert!(sc.width_fits >= 20, "{sc:?}");
    assert!((sc.width_s / WIDTH - 1.0).abs() < 0.1, "{sc:?}");
    assert!(
        sc.describe().contains("ms1_spacing_s=1.0000;"),
        "{}",
        sc.describe()
    );
}

fn narrow_run(fx: &Fixture, width: f64) -> PathBuf {
    let path = fx._dir.path().join("narrow.ms1.parquet");
    let index = Ms1RunIndex::new(run_scans_shaped(1, false, &|_| false, width, 0.25, 5.0));
    write_ms1_table(&path, "run1", &index).unwrap();
    path
}

fn echo_run(fx: &Fixture, echo: f64) -> PathBuf {
    let path = fx._dir.path().join("echo.ms1.parquet");
    let index = Ms1RunIndex::new(run_scans_shaped(1, false, &|_| false, 1.0, 0.25, echo));
    write_ms1_table(&path, "run1", &index).unwrap();
    path
}

fn late(p: u32) -> bool {
    matches!(p, 8 | 15 | 22)
}

fn fixture_late_psms() -> Fixture {
    fixture_psm_rt(false, &|_| false, &|_| false, &|_| false, &|p, run| {
        let after = if run == 1 && late(p) { ECHO_AFTER } else { 0.0 };
        apex(p, run) + 2.0 + after
    })
}

#[test]
fn identified_pair_with_a_psm_on_the_tail_is_measured_on_the_elution_peak() {
    let fx = fixture_late_psms();
    let path = echo_run(&fx, 0.05);
    let out = extract_all(1, &path, &fx.lib, &params()).unwrap();
    assert_eq!(out.counts.consensus_peak, 3);
    for p in 0..BOTH {
        let rows = rows_of(&out.rows, p);
        let (id, pseudo) = (rows[0], rows[3]);
        assert!(is_identified_target(id));
        assert_eq!(pseudo.role, Role::PseudoTransfer);
        let e = id.evidence.as_ref().unwrap();
        assert!((e.t0 - apex(p, 1)).abs() <= 0.5, "precursor {p}: {e:?}");
        assert!((e.amplitude / AMPLITUDE - 1.0).abs() < 0.1, "{p}: {e:?}");
        if late(p) {
            assert_eq!(id.evidence, pseudo.evidence, "precursor {p}");
        }
        for d in &rows[1..3] {
            assert_eq!(d.role, Role::Identified);
            assert_ne!(d.evidence, id.evidence);
        }
    }
}

#[test]
fn identified_pair_keeps_a_psm_peak_that_is_not_much_lower() {
    let fx = fixture_late_psms();
    let path = echo_run(&fx, 0.6);
    let out = extract_all(1, &path, &fx.lib, &params()).unwrap();
    assert_eq!(out.counts.consensus_peak, 0);
    for p in (0..BOTH).filter(|&p| late(p)) {
        let rows = rows_of(&out.rows, p);
        let e = rows[0].evidence.as_ref().unwrap();
        assert!(
            (e.t0 - apex(p, 1) - ECHO_AFTER).abs() <= 0.5,
            "precursor {p}: {e:?}"
        );
        assert!((e.amplitude / (0.6 * AMPLITUDE) - 1.0).abs() < 0.1, "{e:?}");
    }
}

#[test]
fn consensus_peak_needs_another_larger_peak_and_a_signal_at_the_psm() {
    let fx = fixture(false);
    let out = extract_all(1, &fx.ms1_tables[1], &fx.lib, &params()).unwrap();
    let base: Vec<ExtractRow> = rows_of(&out.rows, 3).into_iter().cloned().collect();
    let fit = base[0].evidence.clone().unwrap();
    let taken = |id: Option<&Evidence>, pseudo: Option<&Evidence>| {
        let mut rows = base.clone();
        rows[0].evidence = id.cloned();
        rows[3].evidence = pseudo.cloned();
        let before = rows.clone();
        let (took, _) = prefer_consensus_peak(&mut rows, &[]);
        if took {
            assert_eq!(rows[0].evidence, rows[3].evidence);
            assert_eq!(rows[1..], before[1..]);
        } else {
            assert_eq!(rows, before);
        }
        took
    };
    let at = |t0: f64, amplitude: f64| Evidence {
        t0,
        amplitude,
        w: 4.0,
        ..fit.clone()
    };
    let id = at(500.0, 1e4);
    assert!(taken(Some(&id), Some(&at(440.0, 2.1e4))));
    assert!(taken(Some(&id), Some(&at(508.1, 1e6))));
    assert!(!taken(Some(&id), Some(&at(440.0, 2e4))));
    assert!(!taken(Some(&id), Some(&at(440.0, 5e3))));
    assert!(!taken(Some(&id), Some(&at(508.0, 1e6))));
    assert!(!taken(Some(&id), Some(&at(500.0, 1e6))));
    assert!(!taken(Some(&at(500.0, 0.0)), Some(&at(440.0, 1e6))));
    assert!(!taken(None, Some(&at(440.0, 1e6))));
    assert!(!taken(Some(&id), None));
    let mut rows = base[..3].to_vec();
    assert_eq!(prefer_consensus_peak(&mut rows, &[]), (false, 0));
    let mut rows = base.clone();
    rows[0].evidence = Some(id.clone());
    rows[1].evidence = Some(id.clone());
    rows[2].evidence = Some(id.clone());
    rows[3].evidence = Some(at(500.0, 1.5e4));
    let far = at(440.0, 3e4);
    let (took, decoys) = prefer_consensus_peak(
        &mut rows,
        &[(1, Some(far.clone())), (2, Some(at(440.0, 1.5e4)))],
    );
    assert_eq!((took, decoys), (false, 1));
    assert_eq!(rows[1].evidence.as_ref(), Some(&far));
    assert_eq!(rows[2].evidence.as_ref(), Some(&id));
}

#[test]
fn narrow_peaks_set_a_narrow_width_scale() {
    let fx = fixture(false);
    let path = narrow_run(&fx, 1.0);
    let out = extract_all(1, &path, &fx.lib, &params()).unwrap();
    let sc = out.counts.scales;
    assert_eq!(sc.spacing_s, Some(0.25));
    assert!(!sc.fallback && sc.rounds > 1, "{sc:?}");
    assert!((sc.width_s - 1.0).abs() < 0.1, "{sc:?}");
    for r in out.rows.iter().filter(|r| is_identified_target(r)) {
        let e = r.evidence.as_ref().unwrap();
        let w_prior = e.w / e.width_ratio;
        assert!(
            (w_prior - 1.0).abs() < 0.15,
            "precursor {}: {e:?}",
            r.precursor
        );
    }
}

#[test]
fn identified_fit_stays_on_its_own_peak() {
    let fx = fixture(false);
    let path = narrow_run(&fx, 1.0);
    let out = extract_all(1, &path, &fx.lib, &params()).unwrap();
    let mut n = 0;
    for r in out.rows.iter().filter(|r| is_identified_target(r)) {
        let e = r.evidence.as_ref().unwrap();
        assert!(
            (e.t0 - apex(r.precursor, 1)).abs() <= 0.5,
            "precursor {}: {e:?}",
            r.precursor
        );
        n += 1;
    }
    assert!(n >= 30, "{n}");
}

#[test]
fn too_few_identified_fits_fall_back_to_the_seed_width() {
    let fx = fixture_with(false, &|_| true, &|_| false);
    let out = extract_all(1, &fx.ms1_tables[1], &fx.lib, &params()).unwrap();
    let sc = out.counts.scales;
    assert!(sc.fallback, "{sc:?}");
    assert_eq!((sc.rounds, sc.width_s), (1, SEED_WIDTH_S));
    assert!(out.counts.width_prior_fallback);
}

#[test]
fn faims_entry_is_fitted_on_its_voltage_series() {
    let fx = fixture(true);
    let out = extract_all(1, &fx.ms1_tables[1], &fx.lib, &params()).unwrap();
    for r in out.rows.iter().filter(|r| is_identified_target(r)) {
        let e = r.evidence.as_ref().unwrap();
        assert!(
            (e.amplitude / AMPLITUDE - 1.0).abs() < 0.1,
            "precursor {}: {e:?}",
            r.precursor
        );
    }
}

#[test]
fn rows_round_trip_through_the_extract_table() {
    let fx = fixture(false);
    let out = extract_all(1, &fx.ms1_tables[1], &fx.lib, &params()).unwrap();
    let path = fx._dir.path().join("run1.extract.parquet");
    write_extract_table(&path, &out.rows).unwrap();
    assert_eq!(read_extract_table(&path).unwrap(), out.rows);
}

#[test]
fn faims_run_with_whole_run_entry_holds_one_ms1_copy() {
    let fx = fixture(true);
    let whole = fx
        .lib
        .entries
        .iter()
        .find(|e| e.precursor == ONLY_RUN0)
        .unwrap();
    assert_eq!(whole.faims_cv, None);
    tracking::reset();
    let out = extract_all(1, &fx.ms1_tables[1], &fx.lib, &params()).unwrap();
    let (peak, loads) = tracking::stats();
    assert_eq!(peak, 1, "at most one MS1 copy is live");
    assert_eq!(loads, 2, "the whole run is read again for its entries");
    assert_transfer_layout(&rows_of(&out.rows, ONLY_RUN0));
    let mut seen: Vec<u32> = out.rows.iter().map(|r| r.precursor).collect();
    seen.dedup();
    let order: Vec<u32> = fx.lib.entries.iter().map(|e| e.precursor).collect();
    assert_eq!(seen, order);

    let fx = fixture(false);
    tracking::reset();
    extract_all(1, &fx.ms1_tables[1], &fx.lib, &params()).unwrap();
    assert_eq!(tracking::stats(), (1, 1));
}

#[test]
fn width_prior_ignores_identified_pairs_without_signal() {
    let silent = |p: u32| p % 10 < 3;
    let fx = fixture_with(false, &silent, &|_| false);
    let out = extract_all(1, &fx.ms1_tables[1], &fx.lib, &params()).unwrap();
    assert!(out.counts.width_fits >= 20, "{:?}", out.counts);
    assert!(!out.counts.width_prior_fallback);
    for r in out.rows.iter().filter(|r| is_identified_target(r)) {
        let e = r.evidence.as_ref().unwrap();
        let w_prior = e.w / e.width_ratio;
        assert!(
            (w_prior / WIDTH - 1.0).abs() < 0.1,
            "precursor {}: w_prior {w_prior}",
            r.precursor
        );
    }
}

#[test]
fn faims_transfer_decoys_use_the_target_voltage_series() {
    let mut fx = fixture(true);
    let e = fx
        .lib
        .entries
        .iter_mut()
        .find(|e| e.precursor == ONLY_RUN0)
        .unwrap();
    e.faims_cv = Some(FAIMS_CVS[ONLY_RUN0 as usize % 2]);
    let out = extract_all(1, &fx.ms1_tables[1], &fx.lib, &params()).unwrap();
    let rows = rows_of(&out.rows, ONLY_RUN0);
    assert_transfer_layout(&rows);
    let target = rows[0].evidence.as_ref().unwrap();
    assert!(
        (target.amplitude / AMPLITUDE - 1.0).abs() < 0.1,
        "{target:?}"
    );
    for r in &rows {
        let e = r.evidence.as_ref().unwrap();
        assert!(
            (e.n_points as f64) <= 2.0 * e.w + 1.5,
            "{:?}: {e:?}",
            r.kind
        );
    }
}
