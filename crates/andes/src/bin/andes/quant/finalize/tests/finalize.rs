//! Tests of `quant finalize`: block reading, summaries, the entry index, score files, output rows
//! and wrong-transfer calibration.

use super::index::{edges_of_hist, strata_edges, EntryIndex};
use super::out::{merge_out, OutRow, OutWriter, EVIDENCE_CODES};
use super::scores::ScoreStore;
use super::shapes::MedianWindow;
use super::wrong::{running_mean_q, WrongCalibration};
use crate::quant::mbr::median;
use crate::transfers::{pack_family, stratum_of, MAX_STRATA, TARGET};

#[test]
fn zero_amplitude_transfer_is_not_accepted() {
    use super::transfer_accepted;
    assert!(transfer_accepted(Some(0.001), 0.01, true));
    assert!(!transfer_accepted(Some(0.001), 0.01, false));
    assert!(!transfer_accepted(Some(0.5), 0.01, true));
    assert!(!transfer_accepted(None, 0.01, true));
}

#[test]
fn the_brightest_share_of_precursors_reports_trace_areas() {
    use super::shapes::PrecursorSummary;
    use super::{mark_trace_precursors, Summaries};
    let summaries = |snrs: &[Option<f64>]| -> Summaries {
        snrs.iter()
            .enumerate()
            .map(|(p, &snr)| {
                (
                    p as u32,
                    PrecursorSummary {
                        snr,
                        ..PrecursorSummary::default()
                    },
                )
            })
            .collect()
    };
    let marked = |s: &Summaries| {
        let mut v: Vec<u32> = s.iter().filter(|(_, s)| s.trace).map(|(p, _)| *p).collect();
        v.sort_unstable();
        v
    };
    let mut snrs: Vec<Option<f64>> = (1..=10).map(|v| Some(v as f64)).collect();
    snrs.push(None);
    let mut s = summaries(&snrs);
    mark_trace_precursors(&mut s, 0.3);
    assert_eq!(marked(&s), vec![7, 8, 9]);
    snrs[6] = Some(8.0);
    let mut s = summaries(&snrs);
    mark_trace_precursors(&mut s, 0.25);
    assert_eq!(marked(&s), vec![6, 7, 8, 9]);
    let mut s = summaries(&snrs);
    mark_trace_precursors(&mut s, 1.0);
    assert_eq!(marked(&s), (0..10).collect::<Vec<u32>>());
    let mut s = summaries(&snrs);
    mark_trace_precursors(&mut s, 1e-9);
    assert_eq!(marked(&s), vec![9]);
    let mut s = summaries(&[None, None]);
    mark_trace_precursors(&mut s, 1.0);
    assert!(marked(&s).is_empty());
}

#[test]
fn strata_edges_are_ascending_and_bounded() {
    use quant::library::{DonorPsm, LibraryEntry};
    let psm = DonorPsm {
        score: 1.0,
        q: 0.001,
        pep: 0.001,
    };
    let entries: Vec<LibraryEntry> = (0..400u32)
        .map(|p| LibraryEntry {
            precursor: p,
            identity_decoy: false,
            q: 0.001,
            donor_psm: psm,
            consensus_rt: 10.0,
            faims_cv: None,
            ion_mobility: None,
            donors: (0..6u16)
                .map(|r| r + 6 * (p % 2) as u16)
                .filter(|r| (u32::from(*r) + p) % 7 != 0)
                .collect(),
            best_rt: Vec::new(),
        })
        .collect();
    let entry_of: Vec<Option<&LibraryEntry>> = entries.iter().map(Some).collect();
    let neighbours: Vec<Vec<u16>> = (0..12u16)
        .map(|r| {
            let c = r / 6;
            let mut n: Vec<u16> = (6 * c..6 * c + 6).filter(|&x| x != r).collect();
            n.extend((0..12u16).filter(|x| x / 6 != c));
            n
        })
        .collect();
    let edges = strata_edges(&entry_of, &neighbours, 12, 5);
    assert!(!edges.is_empty());
    assert!(edges.iter().all(|&e| (0.0..1.0).contains(&e)), "{edges:?}");
    assert!(edges.windows(2).all(|w| w[0] < w[1]), "{edges:?}");
    assert!(edges.len() < MAX_STRATA);
}

#[test]
fn the_donor_pep_feature_reads_the_best_donor_psm() {
    use quant::align::AlignParams;
    use quant::library::{build_library, PrecursorMeta, PsmIn};
    let meta: Vec<PrecursorMeta> = (0..200)
        .map(|i| PrecursorMeta {
            key: format!("P{i}/2"),
            mz: 500.0 + i as f64,
            charge: 2,
            envelope: [0.5, 0.3, 0.15, 0.05],
            isotope_share: 1.0,
        })
        .collect();
    let psms = |q: f64, q_gate: f64| -> Vec<PsmIn> {
        (0..200u32)
            .flat_map(|p| {
                (0..2u16).map(move |r| PsmIn {
                    run: r,
                    precursor: p,
                    rt: 10.0 + p as f64,
                    score: 100.0 + p as f64,
                    q,
                    q_gate,
                    pep: 0.001,
                    is_decoy: false,
                    faims_cv: None,
                    ion_mobility: None,
                })
            })
            .collect()
    };
    let ap = AlignParams::default();
    for (q, q_gate) in [(0.0005, 0.005), (0.005, 0.0005), (0.005, f64::NAN)] {
        let lib = build_library(&psms(q, q_gate), meta.clone(), 2, 1.0, 0.01, None, &ap);
        assert_eq!(EntryIndex::new(&lib, 2).donor_pep(0), 0.001);
    }
}

#[test]
fn an_even_prior_histogram_gives_max_strata() {
    let edges = edges_of_hist(&[1; 32]);
    assert_eq!(edges.len() + 1, MAX_STRATA, "{edges:?}");
    assert!(edges.windows(2).all(|w| w[0] < w[1]), "{edges:?}");
    let last = stratum_of(1.0, &edges);
    assert_eq!(usize::from(last), MAX_STRATA - 1);
    let packed = crate::transfers::unpack_family(pack_family(TARGET, last));
    assert_eq!(packed, (TARGET, MAX_STRATA - 1));
    assert!(edges_of_hist(&[]).is_empty());
    assert!(edges_of_hist(&[0, 0, 0]).is_empty());
}

use super::*;

struct Lcg(u64);
impl Lcg {
    fn below(&mut self, n: u64) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 33) % n
    }
}

#[test]
fn median_window_gives_every_median_and_leave_one_out_median() {
    let mut g = Lcg(3);
    for n in 0..40usize {
        for _ in 0..20 {
            let v: Vec<f64> = (0..n)
                .map(|_| match g.below(9) {
                    0 => 0.0,
                    1 => -0.0,
                    k => (k as f64 - 4.0) * 0.25,
                })
                .collect();
            let mut sorted = v.clone();
            sorted.sort_by(f64::total_cmp);
            let w = MedianWindow::new(&sorted);
            let full = (!v.is_empty()).then(|| median(v.clone()));
            assert_eq!(
                w.median().map(f64::to_bits),
                full.map(f64::to_bits),
                "{v:?}"
            );
            for i in 0..n {
                let mut rest = v.clone();
                let x = rest.remove(i);
                let want = (!rest.is_empty()).then(|| median(rest));
                assert_eq!(
                    w.median_without(x).map(f64::to_bits),
                    want.map(f64::to_bits),
                    "{v:?} without {x}"
                );
            }
        }
    }
}

fn out_row(precursor: u32, run: u16) -> OutRow {
    OutRow {
        precursor,
        run,
        evidence_type: EVIDENCE_CODES[(precursor as usize + run as usize) % 5],
        quantity: Some(precursor as f64 * 1.5),
        quantity_se: (run % 2 == 0).then_some(f64::NAN),
        quantity_normalized: None,
        transfer_q: Some(0.25),
        quant_q: Some(-0.0),
        apex_rt_s: Some(f64::INFINITY),
        measurement_q: (run % 3 == 0).then_some(0.125),
        transfer_pep: None,
        library_q: 1e-3,
        mass_err_ppm: Some(-1.25),
        claim: Claim {
            claimant: (precursor % 2 == 0).then_some(precursor + 1),
            dt_s: (precursor % 2 == 0).then_some(-2.5),
            identified: precursor % 4 == 0,
            co_apex: precursor % 2 == 0,
            n: 2,
            n_co_apex: 1,
            n_identified_co_apex: u16::from(precursor % 4 == 0),
            ..Claim::default()
        },
    }
}

#[test]
fn output_tables_merge_in_precursor_and_run_order() {
    let dir = tempfile::tempdir().unwrap();
    let runs = MAX_OPEN_PARTS as u16 + 9;
    let mut paths = Vec::new();
    let mut want = Vec::new();
    for run in 0..runs {
        let rows: Vec<OutRow> = (0..7u32)
            .filter(|p| (p + run as u32) % 3 != 0)
            .map(|p| out_row(p * 2 + u32::from(run % 2), run))
            .collect();
        want.extend(rows.iter().copied());
        let path = dir.path().join(format!("rows-{run:06}.parquet"));
        let mut w = OutWriter::create(&path).unwrap();
        w.write(&rows).unwrap();
        w.finish().unwrap();
        paths.push(path);
    }
    want.sort_by_key(|r| (r.precursor, r.run));
    let mut got = Vec::new();
    merge_out(&paths, dir.path(), 0, &mut |rows| {
        got.extend_from_slice(rows);
        Ok(())
    })
    .unwrap();
    assert_eq!(got.len(), want.len());
    for (a, b) in got.iter().zip(&want) {
        assert_eq!(
            (a.precursor, a.run, a.evidence_type),
            (b.precursor, b.run, b.evidence_type)
        );
        let bits = |r: &OutRow| {
            [
                r.quantity,
                r.quantity_se,
                r.quantity_normalized,
                r.transfer_q,
                r.quant_q,
                r.apex_rt_s,
                r.measurement_q,
                Some(r.library_q),
                r.mass_err_ppm,
            ]
            .map(|v| v.map(f64::to_bits))
        };
        assert_eq!(bits(a), bits(b));
        assert_eq!(a.claim, b.claim);
    }
}

#[test]
fn spilled_scores_give_the_in_memory_q_values() {
    spilled_scores_match(1);
    spilled_scores_match(3);
}

fn spilled_scores_match(strata: usize) {
    let dir = tempfile::tempdir().unwrap();
    let mut g = Lcg(17);
    let rows: Vec<(u32, u32, Option<f32>)> = (0..5000u32)
        .map(|i| {
            let fam = g.below(4) as u32 | ((g.below(strata as u64) as u32) << 4);
            let s = (g.below(41) as f32 - 20.0) / 4.0 + if fam & 0x0F == 0 { 1.0 } else { 0.0 };
            (
                i % 700,
                fam,
                (i % 97 != 0).then_some(if s == 0.0 && i % 2 == 0 { -0.0 } else { s }),
            )
        })
        .collect();
    let mut one = FoldScores::with_strata(5, strata);
    for &(p, f, s) in &rows {
        one.add(p, f, s);
    }
    let store = ScoreStore::new("t", dir.path(), 300, 5, strata);
    for piece in rows.chunks(211) {
        let mut s = FoldScores::with_strata(5, strata);
        for &(p, f, sc) in piece {
            s.add(p, f, sc);
        }
        store.add(s).unwrap();
    }
    let a = one.finish(&[]);
    let (b, spills) = store.finish(&[], &[0.05]).unwrap();
    assert!(spills > 3, "{spills}");
    assert_eq!(a.counts, b.counts);
    assert_eq!(a.pooled_steps(), b.pooled_steps());
    for &(p, f, s) in &rows {
        for e in [
            Estimate::Pooled,
            Estimate::WrongPeak,
            Estimate::MassOnly,
            Estimate::RtOnly,
        ] {
            assert_eq!(a.accepts(e, p, f, s, 0.05), b.accepts(e, p, f, s, 0.05));
        }
        assert_eq!(
            a.q(Estimate::Pooled, p, f, s).map(f64::to_bits),
            b.q(Estimate::Pooled, p, f, s).map(f64::to_bits)
        );
        assert_eq!(
            a.pseudo_q(p, f, s).map(f64::to_bits),
            b.pseudo_q(p, f, s).map(f64::to_bits)
        );
    }
    let left = std::fs::read_dir(dir.path()).unwrap().count();
    assert_eq!(left, 0, "spilled score tables are removed");
}

#[test]
fn wrong_calibration_is_a_non_increasing_step_function_of_the_margin() {
    let mut g = Lcg(5);
    let pairs: Vec<(f32, bool)> = (0..4000)
        .map(|_| {
            let s = g.below(2000) as f32 / 100.0 - 10.0;
            let p = 1.0 / (1.0 + f64::from(s).exp());
            (s, (g.below(1000) as f64) / 1000.0 < p)
        })
        .collect();
    let c = WrongCalibration::fit(pairs).unwrap();
    assert!(c.p.windows(2).all(|w| w[0] >= w[1]));
    assert!(c.lo.windows(2).all(|w| w[0] < w[1]));
    assert!(c.at(-100.0) >= c.at(0.0) && c.at(0.0) >= c.at(100.0));
    assert!(c.at(100.0) < 0.05 && c.at(-100.0) > 0.8);
    assert!(WrongCalibration::fit(Vec::new()).is_none());
    let c =
        WrongCalibration::fit(vec![(1.0, true), (1.0, false), (1.0, false), (1.0, false)]).unwrap();
    assert_eq!(c.p, vec![0.25]);
}

fn quant_qs(rows: &[f64]) -> Vec<f64> {
    let mut pool: Vec<(f64, usize)> = rows.iter().copied().zip(0..).collect();
    let mut q = vec![f64::NAN; rows.len()];
    running_mean_q(&mut pool, |i, v| q[i] = v);
    q
}

#[test]
fn tied_wrong_peak_probabilities_share_one_quant_q() {
    let rows = [0.2, 0.5, 0.1, 0.2, 0.5, 0.2];
    let q = quant_qs(&rows);
    let at = |p: f64| -> Vec<f64> {
        rows.iter()
            .zip(&q)
            .filter(|r| *r.0 == p)
            .map(|r| *r.1)
            .collect()
    };
    assert_eq!(at(0.1), vec![0.1]);
    let tie = (0.1 + 3.0 * 0.2) / 4.0;
    assert!(at(0.2).iter().all(|&v| (v - tie).abs() < 1e-12), "{q:?}");
    let all = (0.1 + 3.0 * 0.2 + 2.0 * 0.5) / 6.0;
    assert!(at(0.5).iter().all(|&v| (v - all).abs() < 1e-12), "{q:?}");
    let mut g = Lcg(11);
    let probs: Vec<f64> = (0..500)
        .map(|_| f64::from(g.below(7) as u32) / 10.0)
        .collect();
    let q = quant_qs(&probs);
    for _ in 0..5 {
        let mut order: Vec<usize> = (0..probs.len()).collect();
        for i in (1..order.len()).rev() {
            order.swap(i, g.below(i as u64 + 1) as usize);
        }
        let permuted: Vec<f64> = order.iter().map(|&i| probs[i]).collect();
        let qp = quant_qs(&permuted);
        for (j, &i) in order.iter().enumerate() {
            assert_eq!(qp[j].to_bits(), q[i].to_bits());
        }
    }
}
