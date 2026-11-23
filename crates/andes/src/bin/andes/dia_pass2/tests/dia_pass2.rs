use super::*;

fn nominee(decoy: bool, llr: f64, ions: Vec<(f64, f32)>) -> Nominee {
    Nominee {
        run: 0,
        charge: 2,
        decoy,
        prec_mz: 500.0,
        pred_rt: 10.0,
        rt_missing: false,
        raw_rt: f64::NAN,
        apex_rt: 10.0,
        dia_llr: llr,
        mobility: f64::NAN,
        ions,
        seq: "PEPTIDEK".into(),
        peptide: "-.PEPTIDEK.-".into(),
        proteins: "P1".into(),
        spec_id: "s_1_1".into(),
        xfer: false,
    }
}

fn window(n_cycles: usize, ions: &[f64], on: std::ops::Range<usize>) -> (Vec<Vec<(f64, f32)>>, Vec<f64>) {
    let mut scans = Vec::new();
    let mut rt = Vec::new();
    for c in 0..n_cycles {
        let mut p: Vec<(f64, f32)> = vec![(150.0, 10.0), (999.0, 10.0)];
        if on.contains(&c) {
            let h = 1.0 + (5.0 - (c as f64 - 33.0).abs()).max(0.0);
            for (j, &mz) in ions.iter().enumerate() {
                p.push((mz, (1000.0 * h * (1.0 + j as f64)) as f32));
            }
        }
        p.sort_by(|a, b| a.0.total_cmp(&b.0));
        scans.push(p);
        rt.push(c as f64 * 0.02);
    }
    (scans, rt)
}

fn place(centre: f64, half: f64, mobility: Option<(u32, u32)>, first: u32) -> Place {
    Place {
        centre,
        iso: (centre - half, centre + half),
        mobility,
        first,
    }
}

#[test]
fn window_of_picks_by_mz_then_mobility_and_breaks_ties_by_first_spectrum() {
    let slices = [
        place(500.0, 1.6, Some((600, 700)), 40),
        place(500.0, 1.6, Some((650, 750)), 10),
        place(500.0, 1.6, Some((700, 800)), 20),
        place(520.0, 1.6, Some((650, 750)), 5),
    ];
    assert_eq!(window_of(&slices, 500.2, 700.0), Some(1));
    assert_eq!(window_of(&slices, 500.2, 750.0), Some(2));
    assert_eq!(window_of(&slices, 500.2, 640.0), Some(0));
    assert_eq!(window_of(&slices, 500.2, 900.0), Some(2), "nearest slice when none holds it");
    assert_eq!(window_of(&slices, 500.2, f64::NAN), Some(1), "earliest of equal slices");
    let fenced = [place(500.0, 1.6, Some((600, 640)), 1), place(501.0, 1.6, Some((650, 700)), 2)];
    assert_eq!(window_of(&fenced, 500.2, 660.0), Some(1), "the window that holds it");
    assert_eq!(window_of(&fenced, 503.0, 660.0), None, "outside both isolation ranges");
    let astral = [place(512.0, 1.6, None, 30), place(500.0, 1.6, None, 7), place(503.0, 1.6, None, 3)];
    assert_eq!(window_of(&astral, 501.0, f64::NAN), Some(1));
    assert_eq!(window_of(&astral, 501.5, f64::NAN), Some(2), "a tie goes to the first spectrum");
    assert_eq!(window_of(&astral, 520.0, f64::NAN), None);
}

#[test]
fn a_precursor_inside_a_wide_isolation_window_is_placed_in_it() {
    let wide = [place(500.0, 10.0, None, 1), place(520.0, 10.0, None, 2)];
    assert_eq!(window_of(&wide, 503.0, f64::NAN), Some(0), "3 Th from the centre of a 20 Th window");
    assert_eq!(window_of(&wide, 509.9, f64::NAN), Some(0));
    assert_eq!(window_of(&wide, 512.5, f64::NAN), Some(1));
    assert_eq!(window_of(&wide, 531.0, f64::NAN), None, "above every window");
    let pasef = [
        place(612.5, 12.5, Some((300, 420)), 4),
        place(612.5, 12.5, Some((420, 540)), 9),
        place(637.5, 12.5, Some((300, 420)), 2),
    ];
    assert_eq!(window_of(&pasef, 621.0, 500.0), Some(1), "the m/z window, then the mobility slice");
    assert_eq!(window_of(&pasef, 630.0, 350.0), Some(2));
    let asym = [Place { centre: 700.0, iso: (699.0, 712.0), mobility: None, first: 0 }];
    assert_eq!(window_of(&asym, 711.0, f64::NAN), Some(0), "asymmetric offsets");
    assert_eq!(window_of(&asym, 698.5, f64::NAN), None);
    let gap = [place(381.42325, 1.0, None, 1), place(383.42416, 1.0, None, 2)];
    assert_eq!(window_of(&gap, 382.4237, f64::NAN), Some(0), "inside the 0.9 mTh gap");
}

#[test]
fn isolation_range_uses_the_offsets_and_falls_back_to_the_reach() {
    let mut s = Spectrum {
        isolation_target_mz: Some(600.0),
        isolation_lower_offset: Some(12.5),
        isolation_upper_offset: Some(12.5),
        ..Default::default()
    };
    assert_eq!(isolation_range(&s), (587.5, 612.5));
    s.isolation_upper_offset = None;
    assert_eq!(isolation_range(&s), (600.0 - WINDOW_REACH_MZ, 600.0 + WINDOW_REACH_MZ));
}

fn edge_features(on: std::ops::Range<usize>, rt_cal_cycle: usize, n_cycles: usize) -> RowFeatures {
    let mzs: Vec<f64> = (0..12).map(|j| 300.0 + 50.0 * j as f64).collect();
    let (scans, rt) = window(n_cycles, &mzs, on);
    let w = Window {
        centre: 500.0,
        iso: (499.0, 501.0),
        mobility: None,
        rt: rt.clone(),
        scans: scans.iter().map(|s| s.as_slice()).collect(),
    };
    let ids: Vec<u32> = (0..n_cycles as u32).collect();
    let n = nominee(false, 20.0, mzs.iter().map(|&m| (m, 1.0)).collect());
    let rt_cal = rt[rt_cal_cycle];
    let t = trace(&w, &n, rt_cal, 1.0);
    let terms = [(2.0, -1.0); TOP_IONS];
    row_features(&w, &ids, &n, &t, &terms, rt_cal)
}

#[test]
fn lsum3_stays_bounded_when_the_best_cycle_is_the_first_of_the_run() {
    let f = edge_features(0..1, 0, 120);
    let v = |c: &str| f.vals[col_ix(c)];
    assert_eq!(v("p2_llr"), 24.0);
    assert_eq!(v("p2_lsum3"), (LLR_FLOOR + 24.0 - 12.0) as f32);
    let f = edge_features(119..120, 119, 120);
    assert_eq!(f.vals[col_ix("p2_lsum3")], (LLR_FLOOR + 24.0 - 12.0) as f32, "last cycle of the run");
}

#[test]
fn lsum3_counts_the_best_cycle_once_at_the_trace_edge() {
    let f = edge_features(0..120, 60, 200);
    let v = |c: &str| f.vals[col_ix(c)];
    assert_eq!(v("p2_llr"), 24.0);
    assert_eq!(v("p2_lsum3"), (LLR_FLOOR + 24.0 + 24.0) as f32, "best at trace cycle 0");
    let short = edge_features(0..3, 1, 3);
    let v = |c: &str| short.vals[col_ix(c)];
    assert_eq!(v("p2_llr2"), LLR_FLOOR as f32, "no far cycle in a three-cycle run");
    assert_eq!(v("p2_lsum3"), (LLR_FLOOR + 24.0 + 24.0) as f32);
}

#[test]
fn rows_without_an_rt_prediction_stay_out_of_the_ridge() {
    let mut noms: Vec<Nominee> = (0..200)
        .map(|i| {
            let mut n = nominee(false, 20.0, vec![(300.0, 1.0)]);
            n.pred_rt = i as f64 * 0.1;
            n.apex_rt = 1.0 + 0.9 * n.pred_rt;
            n
        })
        .collect();
    let clean: Vec<&Nominee> = noms.iter().collect();
    let (cal_clean, sigma_clean) = calibrate_rt(&clean);
    noms.extend((0..200).map(|i| {
        let mut n = nominee(false, 20.0, vec![(300.0, 1.0)]);
        n.pred_rt = 10.0;
        n.apex_rt = i as f64 * 0.3;
        n.rt_missing = true;
        n
    }));
    let all: Vec<&Nominee> = noms.iter().collect();
    let (cal_all, sigma_all) = calibrate_rt(&all);
    assert_eq!(&cal_all[..200], &cal_clean[..]);
    assert_eq!(sigma_all, sigma_clean, "and out of the RT prior's spread");
}

#[test]
fn pass2_llr_peaks_at_the_coeluting_cycle() {
    let mzs: Vec<f64> = (0..12).map(|j| 300.0 + 50.0 * j as f64).collect();
    let (scans, rt) = window(120, &mzs, 32..36);
    let w = Window {
        centre: 500.0,
        iso: (499.0, 501.0),
        mobility: None,
        rt: rt.clone(),
        scans: scans.iter().map(|s| s.as_slice()).collect(),
    };
    let ids: Vec<u32> = (0..120).collect();
    let n = nominee(false, 20.0, mzs.iter().map(|&m| (m, 1.0)).collect());
    let rt_cal = rt[33];
    let t = trace(&w, &n, rt_cal, 1.0);
    let mut terms = [(2.0, -1.0); TOP_IONS];
    terms[11] = (1.0, -0.5);
    let f = row_features(&w, &ids, &n, &t, &terms, rt_cal);
    let v = |c: &str| f.vals[col_ix(c)];
    assert!(v("p2_rtstar") >= rt[30] as f32 && v("p2_rtstar") <= rt[36] as f32);
    assert_eq!(v("p2_n6"), 6.0);
    assert_eq!(v("p2_nall"), 12.0);
    assert!(v("p2_llr") > v("p2_llr2"));
    assert!(v("p2_cor") > 0.9);
}

#[test]
fn pass2_shared_peak_goes_to_the_higher_llr() {
    let mk = |llr: f32, key: PeakKey| {
        let mut vals = [0f32; COLUMNS.len()];
        vals[col_ix("p2_llr")] = llr;
        RowFeatures {
            vals,
            tr_full: [0.0; CYCLES],
            tt: [0.0; CYCLES],
            claims: vec![(key, 2.0)],
            apex: 0,
            quant: RowQuant::default(),
            quant_window: None,
        }
    };
    let mut rows = vec![(0, mk(10.0, (5, 7))), (1, mk(4.0, (5, 7)))];
    share_peaks(&mut rows);
    assert_eq!(rows[0].1.vals[col_ix("p2_shared")], 0.0);
    assert_eq!(rows[1].1.vals[col_ix("p2_shared")], 1.0);
    assert_eq!(rows[1].1.vals[col_ix("p2_lone")], 2.0);
}

#[test]
fn transferred_rows_stay_out_of_the_rt_ridge_and_the_rates() {
    let mut noms: Vec<Nominee> = (0..200)
        .map(|i| {
            let mut n = nominee(false, 20.0, vec![(300.0, 1.0)]);
            n.pred_rt = i as f64 * 0.1;
            n.apex_rt = 1.0 + 0.9 * n.pred_rt;
            n
        })
        .collect();
    let own: Vec<&Nominee> = noms.iter().collect();
    let (cal_own, sigma_own) = calibrate_rt(&own);
    noms.extend((0..200).map(|i| {
        let mut n = nominee(i % 2 == 0, 20.0, vec![(300.0, 1.0)]);
        n.pred_rt = i as f64 * 0.1;
        n.apex_rt = 50.0 - n.pred_rt;
        n.xfer = true;
        n
    }));
    let all: Vec<&Nominee> = noms.iter().collect();
    let (cal_all, sigma_all) = calibrate_rt(&all);
    assert_eq!(&cal_all[..200], &cal_own[..]);
    assert_eq!(sigma_all, sigma_own);

    let mzs: Vec<f64> = (0..12).map(|j| 300.0 + 50.0 * j as f64).collect();
    let (scans, rt) = window(120, &mzs, 32..36);
    let w = Window {
        centre: 500.0,
        iso: (499.0, 501.0),
        mobility: None,
        rt: rt.clone(),
        scans: scans.iter().map(|s| s.as_slice()).collect(),
    };
    let ions: Vec<(f64, f32)> = mzs.iter().map(|&m| (m, 1.0)).collect();
    for decoy in [false, true] {
        let mut n = nominee(decoy, 20.0, ions.clone());
        n.apex_rt = rt[33];
        let t = trace(&w, &n, rt[33], 1.0);
        let r = collect_rates(&w, &n, &t);
        assert!(r.qn.iter().chain(&r.pn).sum::<f64>() > 0.0);
        n.xfer = true;
        let r = collect_rates(&w, &n, &t);
        assert_eq!(r.qn.iter().chain(&r.pn).sum::<f64>(), 0.0);
        assert!(!confident(&n) && !null_member(&n));
    }
}

#[test]
fn xfer_column_marks_transferred_rows() {
    let mut noms = vec![
        nominee(false, 20.0, vec![(300.0, 1.0)]),
        nominee(true, 20.0, vec![(300.0, 1.0)]),
        nominee(false, TRANSFER_DIA_LLR, vec![(300.0, 1.0)]),
    ];
    let p2 = Pass2 {
        cols: (0..COLUMNS.len())
            .map(|c| if COLUMNS[c] == "p2_ok" { vec![1.0; 3] } else { vec![0.0; 3] })
            .collect(),
        joint: vec![vec![0.0; 3]; crate::dia_joint::COLUMNS.len()],
        quant: vec![
            RowQuant::default(),
            RowQuant {
                frag: Some(5.0),
                n_frag: 2,
                ms1: Some(7.0),
            },
            RowQuant::default(),
        ],
    };
    assert!(to_table(&noms, &p2).feature(XFER_COLUMN).is_none());
    noms[2].xfer = true;
    let t = to_table(&noms, &p2);
    assert_eq!(t.feature(XFER_COLUMN).unwrap(), &[0.0, 0.0, 1.0]);
    assert_eq!(t.features.last().unwrap().0, XFER_COLUMN);
    assert_eq!(t.quant_frag, vec![None, Some(5.0), None]);
    assert_eq!(t.quant_ms1, vec![None, Some(7.0), None]);
    assert_eq!(t.quant_n_frag, vec![0, 2, 0]);
    assert!(t.features.iter().all(|(n, _)| !n.starts_with("quant")));
}

/// Gaussian elution profile over `cycles` cycles, apex at `apex`, sd `sd` cycles.
fn gauss(cycles: usize, apex: f64, sd: f64) -> Vec<f32> {
    (0..cycles)
        .map(|k| (-((k as f64 - apex).powi(2)) / (2.0 * sd * sd)).exp() as f32)
        .collect()
}

fn area(t: &[f32], w: Range<usize>) -> f64 {
    t[w].iter().map(|&v| f64::from(v)).sum()
}

#[test]
fn fragment_quantity_sums_the_six_most_intense_predicted_fragments() {
    let g = gauss(41, 20.0, 2.0);
    let valid = [true; 41];
    let mut traces: Vec<Vec<f32>> = (0..6)
        .map(|j| g.iter().map(|v| v * 1000.0 * (j + 1) as f32).collect())
        .collect();
    traces.push(g.iter().map(|v| v * 1e6).collect());
    let w = 17..24;
    let (q, n) = fragment_quantity(&traces, &valid, w.clone(), 3);
    assert_eq!(n, 6);
    let want: f64 = traces[..6].iter().map(|t| area(t, w.clone())).sum();
    let got = q.unwrap();
    assert!((got - want).abs() < 1e-6 * want, "{got} vs {want}");
}

#[test]
fn fragment_quantity_subtracts_each_fragment_background() {
    let g = gauss(41, 20.0, 2.0);
    let valid = [true; 41];
    let floor = [50.0f32, 0.0, 400.0];
    let traces: Vec<Vec<f32>> = [300.0f32, 900.0, 150.0]
        .iter()
        .zip(floor)
        .map(|(a, b)| g.iter().map(|v| v * a + b).collect())
        .collect();
    let w = 17..24;
    let (q, n) = fragment_quantity(&traces, &valid, w.clone(), 3);
    assert_eq!(n, 3);
    let want: f64 = [300.0f64, 900.0, 150.0]
        .iter()
        .map(|a| a * area(&g, w.clone()))
        .sum();
    assert!((q.unwrap() - want).abs() < 1e-3 * want, "{q:?} vs {want}");
    let mut few = [false; 41];
    few[14..27].fill(true);
    let (q, _) = fragment_quantity(&traces, &few, w.clone(), 3);
    let raw: f64 = traces.iter().map(|t| area(t, w.clone())).sum();
    assert!((q.unwrap() - raw).abs() < 1e-6 * raw, "without background cycles nothing is subtracted");
}

#[test]
fn quantity_missing_without_signal() {
    let valid = [true; 41];
    let zero = vec![vec![0f32; 41]; 6];
    assert_eq!(fragment_quantity(&zero, &valid, 17..24, 3), (None, 0));
    let flat = vec![vec![70f32; 41]; 6];
    assert_eq!(fragment_quantity(&flat, &valid, 17..24, 3), (None, 0), "background only");
    let g = gauss(41, 20.0, 2.0);
    let one = vec![g.clone()];
    assert_eq!(fragment_quantity(&one, &valid, 20..20, 3), (None, 0), "empty window");
    assert_eq!(fragment_quantity(&[], &valid, 17..24, 3), (None, 0), "no fragments");
    let g = gauss(21, 10.0, 2.0);
    assert_eq!(ms1_quantity(&[0.0; 9], &[0.0; 9], 0..9), None);
    assert_eq!(ms1_quantity(&[0.0; 9], &[5.0; 9], 0..9), None, "M+1 alone is not a quantity");
    let mono: Vec<f32> = g.iter().map(|v| v * 1e4).collect();
    let m1: Vec<f32> = g.iter().map(|v| v * 5e3).collect();
    let q = ms1_quantity(&mono, &m1, 6..15).unwrap();
    let want = area(&mono, 6..15) + area(&m1, 6..15);
    assert!((q - want).abs() < 1e-6 * want);
    assert_eq!(ms1_quantity(&mono, &m1, 30..40), None, "window past the trace");
}

#[test]
fn transferred_row_at_window_edge_is_missing() {
    let mzs: Vec<f64> = (0..12).map(|j| 300.0 + 50.0 * j as f64).collect();
    let rt: Vec<f64> = (0..40).map(|c| c as f64 * 0.02).collect();
    let scans: Vec<Vec<(f64, f32)>> = (0..40)
        .map(|c| {
            let h = (c as usize).saturating_sub(29) as f32 * 1000.0;
            mzs.iter().filter(|_| h > 0.0).map(|&m| (m, h)).collect()
        })
        .collect();
    let w = Window {
        centre: 500.0,
        iso: (499.0, 501.0),
        mobility: None,
        rt: rt.clone(),
        scans: scans.iter().map(|s| s.as_slice()).collect(),
    };
    let ions: Vec<(f64, f32)> = mzs.iter().map(|&m| (m, 1.0)).collect();
    let mut n = nominee(false, TRANSFER_DIA_LLR, ions);
    n.xfer = true;
    let terms = [(2.0, -1.0); TOP_IONS];
    let ids: Vec<u32> = (0..40).collect();
    for rt_cal in [rt[39], rt[39] + 100.0] {
        let t = trace(&w, &n, rt_cal, 1.0);
        let f = row_features(&w, &ids, &n, &t, &terms, rt_cal);
        let (q, win) = row_quantity(&t, n.ions.len(), f.apex, 3);
        assert_eq!(q, RowQuant::default(), "rt_cal {rt_cal}");
        assert!(win.is_none());
    }
    let (scans, rt) = window(120, &mzs, 28..39);
    let w = Window {
        centre: 500.0,
        iso: (499.0, 501.0),
        mobility: None,
        rt: rt.clone(),
        scans: scans.iter().map(|s| s.as_slice()).collect(),
    };
    let ids: Vec<u32> = (0..120).collect();
    let t = trace(&w, &n, rt[33], 1.0);
    let f = row_features(&w, &ids, &n, &t, &terms, rt[33]);
    let (q, win) = row_quantity(&t, n.ions.len(), f.apex, 3);
    let win = win.unwrap();
    assert_eq!(win.len(), 7);
    assert_eq!(t.rt(&w, win.start + 3), rt[33], "the window is centred on the trace maximum");
    assert_eq!(q.n_frag, 6);
    assert!(q.frag.is_some_and(|v| v > 0.0));
}

#[test]
fn run_peak_width_is_the_median_half_maximum_width() {
    let mut few = vec![9usize; QUANT_WIDTH_MIN_TARGETS - 1];
    assert_eq!(quant_half_cycles(&mut few), QUANT_HALF_CYCLES_FALLBACK);
    let mut w: Vec<usize> = (0..101).map(|i| 1 + i % 9).collect();
    assert_eq!(quant_half_cycles(&mut w), 5);
    let mut narrow = vec![1usize; 60];
    assert_eq!(quant_half_cycles(&mut narrow), QUANT_MIN_HALF_CYCLES);
    let mut wide = vec![500usize; 60];
    assert_eq!(quant_half_cycles(&mut wide), HALF_CYCLES - 1);

    let mzs: Vec<f64> = (0..12).map(|j| 300.0 + 50.0 * j as f64).collect();
    let (scans, rt) = window(120, &mzs, 28..39);
    let win = Window {
        centre: 500.0,
        iso: (499.0, 501.0),
        mobility: None,
        rt: rt.clone(),
        scans: scans.iter().map(|s| s.as_slice()).collect(),
    };
    let mut n = nominee(false, 20.0, mzs.iter().map(|&m| (m, 1.0)).collect());
    n.apex_rt = rt[33];
    let t = trace(&win, &n, rt[33], 1.0);
    assert_eq!(half_max_width(&win, &n, &t), Some(7));
}

#[test]
fn ridge_recovers_a_linear_rt_map() {
    let noms: Vec<Nominee> = (0..200)
        .map(|i| {
            let mut n = nominee(false, 20.0, vec![(300.0, 1.0)]);
            n.pred_rt = i as f64 * 0.1;
            n.apex_rt = 1.0 + 0.9 * n.pred_rt;
            n
        })
        .collect();
    let refs: Vec<&Nominee> = noms.iter().collect();
    let (cal, _) = calibrate_rt(&refs);
    for (n, c) in noms.iter().zip(&cal) {
        assert!((n.apex_rt - c).abs() < 0.05, "{} vs {}", n.apex_rt, c);
    }
}
