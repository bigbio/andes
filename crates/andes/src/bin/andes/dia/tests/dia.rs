use super::*;
use onnx_infer::predictors::form::parse_peptide;

#[test]
fn classwise_q_competes_within_each_class() {
    let rows = vec![
        row(0, 10.0, false),
        row(1, 9.0, true),
        row(2, 8.0, false),
        row(3, 7.0, false),
        row(4, 6.0, true),
    ];
    let class = [false, true, false, true, false];
    let q = classwise_q_values(5, &rows, &class);
    assert!((q[0] - 0.5).abs() < 1e-12, "{}", q[0]);
    assert!((q[2] - 0.5).abs() < 1e-12, "{}", q[2]);
    assert_eq!(q[3], 1.0);
    let pooled = classwise_q_values(5, &rows, &[false; 5]);
    for i in [0, 2, 3] {
        assert!((pooled[i] - 2.0 / 3.0).abs() < 1e-12, "pooled q[{i}] = {}", pooled[i]);
    }
}

fn row(k: u32, score: f64, decoy: bool) -> (u32, f64, bool, String) {
    (k, score, decoy, format!("r{k}_{score}"))
}

#[test]
fn precursor_q_counts_each_precursor_once() {
    let rows = vec![
        row(0, 10.0, false),
        row(0, 9.0, false),
        row(0, 8.0, false),
        row(1, 7.0, false),
        row(2, 6.0, true),
        row(3, 5.0, false),
    ];
    let q = precursor_q_values(4, &rows);
    assert_eq!(q[0], Some((10.0, false, 0.5)));
    assert_eq!(q[1], Some((7.0, false, 0.5)));
    assert_eq!(q[2], Some((6.0, true, 2.0 / 3.0)));
    assert_eq!(q[3], Some((5.0, false, 2.0 / 3.0)));
}

#[test]
fn precursor_q_takes_the_best_row_and_counts_tied_decoys_first() {
    let rows = vec![row(0, 3.0, false), row(0, 4.0, false), row(1, 4.0, true)];
    let q = precursor_q_values(2, &rows);
    assert_eq!(q[0], Some((4.0, false, 1.0)));
    assert_eq!(q[1], Some((4.0, true, 1.0)));
}

#[test]
fn row_q_is_the_precursor_threshold_at_its_score() {
    let q = vec![
        Some((10.0, false, 0.1)),
        Some((7.0, false, 0.2)),
        Some((6.0, true, 0.3)),
    ];
    let curve = q_curve(&q);
    assert_eq!(curve_at(&curve, 11.0), 1.0);
    assert_eq!(curve_at(&curve, 10.0), 0.1);
    assert_eq!(curve_at(&curve, 8.0), 0.1);
    assert_eq!(curve_at(&curve, 7.0), 0.2);
    assert_eq!(curve_at(&curve, 1.0), 0.3);
}

#[test]
fn coeluting_fragments_and_precursor_score_high() {
    let (_, mass) = parse_peptide("K.PEPTIDEKLR.A").unwrap();
    let z = 2u8;
    let total: f64 = mass.res.iter().sum::<f64>() + mass.nterm;
    let mono = (total + H2O) / z as f64 + PROTON;
    let mut prefix = mass.nterm;
    let mut ys = Vec::new();
    for r in &mass.res[..mass.res.len() - 3] {
        prefix += r;
        ys.push(total - prefix + H2O + PROTON);
    }
    let shape = [0.0, 0.1, 0.4, 0.9, 1.0, 0.7, 0.3, 0.1, 0.0, 0.0, 0.0];
    let rts: Vec<f64> = (0..shape.len()).map(|k| k as f64 * 1.6).collect();
    let mut scans: Vec<Vec<(f64, f32)>> = Vec::new();
    let mut survey: Vec<Ms1Scan> = Vec::new();
    for (k, &h) in shape.iter().enumerate() {
        let mut p: Vec<(f64, f32)> = ys
            .iter()
            .enumerate()
            .map(|(j, &mz)| (mz, (1e5 * h * (1.0 + j as f64)) as f32 + 10.0))
            .collect();
        p.sort_by(|a, b| a.0.total_cmp(&b.0));
        scans.push(p);
        let env = averagine_isotope_envelope(total + H2O, 3);
        survey.push(Ms1Scan {
            rt: rts[k],
            peaks: (0..3)
                .map(|i| {
                    (
                        mono + i as f64 * ISOTOPE / 2.0,
                        (1e6 * h * env[i]) as f32 + 1.0,
                    )
                })
                .collect(),
            ..Default::default()
        });
    }
    let refs: Vec<&[(f64, f32)]> = scans.iter().map(Vec::as_slice).collect();
    let ms1 = Ms1RunIndex::new(survey);
    let f = row_features(&refs, 4, &rts, &mass, z, Some(&ms1));
    assert!(f[0] > 0.95, "frag corr {}", f[0]);
    assert_eq!(f[2], 0.0, "apex offset");
    assert!(f[4] > 0.95, "ms1 corr {}", f[4]);
    assert_eq!(f[5], 1.0);
    assert!(f[8] > 0.99, "iso cos {}", f[8]);

    let (_, other) = parse_peptide("K.GGGGWWWWK.A").unwrap();
    let g = row_features(&refs, 4, &rts, &other, z, Some(&ms1));
    assert_eq!((g[0], g[5], g[10]), (0.0, 0.0, 0.0));

    let n = ys.len() as f64;
    let (n2, _, gain) = apex_sum(refs[4], Some(refs[3]), Some(refs[5]), &mass, z);
    assert_eq!((n2, gain), (n, 0.0));
    assert_eq!(apex_sum(refs[4], None, None, &mass, z).0, 0.0);
    assert_eq!(apex_sum(refs[4], Some(refs[3]), None, &mass, z).0, n);
    assert_eq!(apex_sum(refs[4], None, None, &other, z).0, 0.0);
}

#[test]
fn windows_come_out_in_key_order() {
    let targets: Vec<f64> = (0..64).map(|i| 400.0 + ((i * 37) % 64) as f64 * 2.5).collect();
    let spectra: Vec<Spectrum> = targets
        .iter()
        .enumerate()
        .map(|(i, &t)| Spectrum {
            rt_seconds: Some(i as f64),
            isolation_target_mz: Some(t),
            ..Default::default()
        })
        .collect();
    let (_, windows) = windows(&spectra, &[spectra.len()]);
    let centre = |w: &Vec<u32>| spectra[w[0] as usize].isolation_target_mz.unwrap();
    let got: Vec<f64> = windows.iter().map(centre).collect();
    let mut want = targets.clone();
    want.sort_by(f64::total_cmp);
    assert_eq!(got, want);
}

#[test]
fn windows_keep_mobility_windows_with_one_centre_apart() {
    let spec = |rt: f64, mobility_window: Option<(u32, u32)>| Spectrum {
        rt_seconds: Some(rt),
        isolation_target_mz: Some(600.0),
        mobility_window,
        ..Default::default()
    };
    let spectra = [
        spec(2.0, Some((34, 500))),
        spec(1.0, Some((500, 935))),
        spec(1.0, Some((34, 500))),
        spec(2.0, Some((500, 935))),
        spec(3.0, None),
    ];
    let (slot, windows) = windows(&spectra, &[spectra.len()]);
    assert_eq!(windows.len(), 3);
    let at = |i: usize| &windows[slot[i].0 as usize];
    assert_eq!(at(0), &vec![2, 0]);
    assert_eq!(at(1), &vec![1, 3]);
    assert_eq!(at(4), &vec![4]);
    assert_eq!((slot[0].1, slot[2].1), (1, 0));
}

#[test]
fn refine_rt_learns_a_run_offset() {
    let dir = tempfile::tempdir().unwrap();
    let pin = dir.path().join("x.pin");
    let mut text = String::from(
        "SpecId\tLabel\tScanNr\tFileName\tRawScore\tcharge2\tandesrt_predicted_rt\tPeptide\tProteins\n",
    );
    let mut spectra = Vec::new();
    for i in 0..1200usize {
        let core: String = (0..4)
            .map(|k| b"AELKSVTG"[(i >> (3 * k)) & 7] as char)
            .collect();
        let w = i % 3 == 0;
        let pep = format!("K.{core}{}R.A", if w { "W" } else { "F" });
        let pred = 10.0 + (i % 100) as f64 * 0.5;
        let obs = pred + if w { 2.0 } else { 0.0 };
        let scan = i as i32 + 1;
        spectra.push(Spectrum {
            title: format!("t.{scan}"),
            precursor_mz: 500.0,
            precursor_intensity: None,
            precursor_charge: Some(2),
            rt_seconds: Some(obs * 60.0),
            scan: Some(scan),
            peaks: vec![(100.0, 1.0)],
            activation_method: None,
            isolation_lower_offset: None,
            isolation_upper_offset: None,
            isolation_target_mz: None,
            mobility_window: None,
        });
        text.push_str(&format!(
            "t.{scan}_{scan}_1\t1\t{scan}\tr\t{}\t1\t{pred}\t{pep}\tP1\n",
            100 - i % 50
        ));
    }
    std::fs::write(&pin, &text).unwrap();
    assert!(refine_rt(&pin, &spectra, &|_| true).unwrap());
    let out = std::fs::read_to_string(&pin).unwrap();
    let h: Vec<&str> = out.lines().next().unwrap().split('\t').collect();
    let c = h.iter().position(|c| *c == RT_REFINED).unwrap();
    assert_eq!(h[c + 1], "Peptide");
    let mut e: Vec<f64> = out
        .lines()
        .skip(1)
        .map(|l| l.split('\t').nth(c).unwrap().parse().unwrap())
        .collect();
    e.sort_by(f64::total_cmp);
    assert!(
        e[e.len() * 9 / 10] < 0.5,
        "90th percentile {}",
        e[e.len() * 9 / 10]
    );

    let bare = text.replace("\tandesrt_predicted_rt", "\tother");
    std::fs::write(&pin, &bare).unwrap();
    assert!(!refine_rt(&pin, &spectra, &|_| true).unwrap());
    assert_eq!(std::fs::read_to_string(&pin).unwrap(), bare);

    for (i, s) in spectra.iter_mut().enumerate() {
        let pred = 10.0 + (i % 100) as f64 * 0.5;
        s.rt_seconds = Some((pred + ((i * 7919) % 13) as f64 * 0.1) * 60.0);
    }
    std::fs::write(&pin, &text).unwrap();
    assert!(refine_rt(&pin, &spectra, &|_| true).unwrap());
}

#[test]
fn train_fdr_defaults_unless_set() {
    let mut a = vec!["--testFDR".to_string(), "0.05".to_string()];
    default_train_fdr(&mut a);
    assert_eq!(a, ["--testFDR", "0.05", "--trainFDR", "0.05"]);
    for set in [vec!["-F", "0.02"], vec!["--trainFDR=0.1"]] {
        let mut b: Vec<String> = set.iter().map(|s| s.to_string()).collect();
        default_train_fdr(&mut b);
        assert_eq!(b, set);
    }
}

#[test]
fn rerank_keeps_top_rows_per_scan_by_score_and_angle() {
    let dir = tempfile::tempdir().unwrap();
    let pin = dir.path().join("x.pin");
    let rows = [
        ("r/c=1 scan=1_1_1_0", 1, 10.0, 0.1),
        ("r/c=1 scan=1_1_2_1", 1, 9.0, 0.9),
        ("r/c=1 scan=1_1_3_2", 1, 8.0, 0.2),
        ("r/c=2 scan=1_1_1_0", 1, 1.0, 0.1),
        ("r/c=2 scan=1_1_2_1", 1, 0.5, 0.1),
        ("r/c=1 scan=2_2_1", 2, 1.0, 0.1),
    ];
    let mut s = String::from("SpecId\tLabel\tScanNr\tRawScore\tandesms2_spectral_angle\n");
    for (id, scan, raw, ms2) in rows {
        s.push_str(&format!("{id}\t1\t{scan}\t{raw}\t{ms2}\n"));
    }
    let quiet = dir.path().join("q.pin");
    std::fs::write(&quiet, &s).unwrap();
    std::fs::write(&pin, s).unwrap();
    let dropped = rerank_pin(&pin, 1, true).unwrap();
    assert!(rerank_pin(&quiet, 1, false).unwrap().is_empty());
    assert_eq!(
        std::fs::read(&quiet).unwrap(),
        std::fs::read(&pin).unwrap(),
        "the kept rows do not depend on collecting the dropped ids"
    );
    let want: std::collections::HashSet<String> = [
        "r/c=1 scan=1_1_1_0",
        "r/c=1 scan=1_1_3_2",
        "r/c=2 scan=1_1_2_1",
    ]
    .map(String::from)
    .into();
    assert_eq!(dropped, want);
    let ids: Vec<String> = std::fs::read_to_string(&pin)
        .unwrap()
        .lines()
        .skip(1)
        .map(|l| l.split('\t').next().unwrap().to_string())
        .collect();
    assert_eq!(
        ids,
        [
            "r/c=1 scan=1_1_2_1",
            "r/c=2 scan=1_1_1_0",
            "r/c=1 scan=2_2_1"
        ]
    );
}

#[test]
fn rerank_drops_a_row_far_from_its_predicted_rt() {
    let dir = tempfile::tempdir().unwrap();
    let pin = dir.path().join("x.pin");
    let mut s = String::from(
        "SpecId\tLabel\tScanNr\tFileName\tRawScore\tandesms2_spectral_angle\t\
         andesrt_rt_diff\tandesrt_missing\n",
    );
    let id = |scan: u32, k: u32| format!("r/c=1 scan={scan}_{scan}_{k}_{}", k - 1);
    for scan in 1..=120u32 {
        let d = 0.3 * f64::from(scan % 5) / 4.0;
        s.push_str(&format!("{}\t1\t{scan}\tr\t10\t0.5\t{d}\t0\n", id(scan, 2)));
    }
    s.push_str(&format!("{}\t1\t1\tr\t10.5\t0.5\t5.0\t0\n", id(1, 3)));
    s.push_str(&format!("{}\t1\t1\tr\t10.2\t0.5\t0.0\t0\n", id(1, 4)));
    std::fs::write(&pin, &s).unwrap();
    let dropped = rerank_pin(&pin, 2, true).unwrap();
    assert_eq!(dropped, [id(1, 3)].into());

    let bare = s.replace("andesrt_rt_diff", "other");
    std::fs::write(&pin, &bare).unwrap();
    let dropped = rerank_pin(&pin, 2, true).unwrap();
    assert_eq!(dropped, [id(1, 2)].into(), "no RT column, no RT term");
}

#[test]
fn rerank_keeps_the_earlier_of_tied_rows_of_a_scan_spread_over_the_pin() {
    let dir = tempfile::tempdir().unwrap();
    let pin = dir.path().join("x.pin");
    let s = "SpecId\tLabel\tScanNr\tRawScore\tandesms2_spectral_angle\r\n\
             r_7_7_1\t1\t7\t2.0\t0.5\n\
             r_8_8_1\t1\t8\t5.0\t0.5\n\
             \n\
             r_7_7_2\t-1\t7\t2.0\t0.5\n\
             r_8_8_2\t1\t8\t1.0\tx\n\
             r_7_7_3\t1\t7\t1.0\t0.1";
    std::fs::write(&pin, s).unwrap();
    let dropped = rerank_pin(&pin, 1, true).unwrap();
    let want: std::collections::HashSet<String> =
        ["r_7_7_2", "r_8_8_2", "r_7_7_3"].map(String::from).into();
    assert_eq!(dropped, want);
    assert_eq!(
        std::fs::read_to_string(&pin).unwrap(),
        "SpecId\tLabel\tScanNr\tRawScore\tandesms2_spectral_angle\n\
         r_7_7_1\t1\t7\t2.0\t0.5\n\
         r_8_8_1\t1\t8\t5.0\t0.5\n"
    );
}

#[test]
fn evidence_ignores_an_interfered_fragment() {
    let (_, mass) = parse_peptide("K.PEPTIDEKLR.A").unwrap();
    let total: f64 = mass.res.iter().sum::<f64>() + mass.nterm;
    let ions = [(false, 1u8)];
    let n = mass.res.len() - 1;
    let pred: Vec<f32> = (0..n).map(|i| 1.0 + i as f32).collect();
    let mut prefix = mass.nterm;
    let ys: Vec<f64> = mass.res[..n]
        .iter()
        .map(|r| {
            prefix += r;
            total - prefix + H2O + PROTON
        })
        .collect();
    let shape = [0.0, 0.2, 0.6, 1.0, 0.6, 0.2, 0.0];
    let build = |bad: u8| -> Vec<Vec<(f64, f32)>> {
        shape
            .iter()
            .enumerate()
            .map(|(k, &h)| {
                let mut p: Vec<(f64, f32)> = ys
                    .iter()
                    .zip(&pred)
                    .enumerate()
                    .map(|(j, (&mz, &w))| {
                        let top = j == n - 1;
                        let v = match bad {
                            1 if top => h + 3.0,
                            2 if top => 5.0 * (k as f64 / 6.0),
                            _ => h,
                        };
                        (mz, (1e4 * w as f64 * v) as f32 + 1.0)
                    })
                    .collect();
                p.sort_by(|a, b| a.0.total_cmp(&b.0));
                p
            })
            .collect()
    };
    let score = |s: &Vec<Vec<(f64, f32)>>| {
        let v: Vec<&[(f64, f32)]> = s.iter().map(|x| x.as_slice()).collect();
        evidence(&pred, &ions, &mass, 2, &v, 3, FragmentTolerance::Ppm(7.0))
    };
    let clean = score(&build(0));
    assert!(clean[0] > 0.99 && clean[4] > 0.99, "{clean:?}");
    let off = build(0);
    let v: Vec<&[(f64, f32)]> = off[1..].iter().map(|x| x.as_slice()).collect();
    let from_side = evidence(&pred, &ions, &mass, 2, &v, 1, FragmentTolerance::Ppm(7.0));
    let at_apex = evidence(&pred, &ions, &mass, 2, &v, 2, FragmentTolerance::Ppm(7.0));
    assert_eq!(from_side[0], at_apex[0], "the sum is centred on the apex");
    let shifted: Vec<Vec<(f64, f32)>> = off
        .iter()
        .map(|s| s.iter().map(|&(mz, i)| (mz + 0.02, i)).collect())
        .collect();
    let v: Vec<&[(f64, f32)]> = shifted.iter().map(|x| x.as_slice()).collect();
    assert_eq!(
        evidence(&pred, &ions, &mass, 2, &v, 3, FragmentTolerance::Ppm(7.0))[0],
        0.0
    );
    assert!(evidence(&pred, &ions, &mass, 2, &v, 3, FragmentTolerance::Da(0.05))[0] > 0.99);
    let flat = score(&build(1));
    assert!(
        flat[0] < 0.95,
        "background lowers the summed angle: {flat:?}"
    );
    assert!(flat[4] > 0.99, "the baseline fit removes it: {flat:?}");
    let ramp = score(&build(2));
    assert!(
        ramp[3] < 1.0,
        "an interfered trace is not counted clean: {ramp:?}"
    );
}

#[test]
fn fragment_traces_match_a_binary_search_per_fragment() {
    let mut x = 0x2545_f491_4f6c_dd1du64;
    let mut next = || {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        (x >> 11) as f64 / (1u64 << 53) as f64
    };
    for tol in [FragmentTolerance::Ppm(20.0), FragmentTolerance::Da(0.02)] {
        for _ in 0..50 {
            let scans: Vec<Vec<(f64, f32)>> = (0..7)
                .map(|_| {
                    let mut p: Vec<(f64, f32)> = (0..300)
                        .map(|_| (100.0 + 400.0 * next(), (1000.0 * next()) as f32))
                        .collect();
                    p.sort_by(|a, b| a.0.total_cmp(&b.0));
                    p
                })
                .collect();
            let refs: Vec<&[(f64, f32)]> = scans.iter().map(Vec::as_slice).collect();
            let mut frags: Vec<(f64, f64)> = (0..120)
                .map(|_| (100.0 + 400.0 * next(), next()))
                .collect();
            frags.push(frags[3]);
            frags.push((refs[2][10].0, 1.0));
            let got = fragment_traces(&frags, &refs, tol);
            for (k, &(mz, _)) in frags.iter().enumerate() {
                for (s, peaks) in refs.iter().enumerate() {
                    let w = tol.window(mz);
                    let lo = peaks.partition_point(|p| p.0 < mz - w);
                    let want = peaks[lo..]
                        .iter()
                        .take_while(|p| p.0 <= mz + w)
                        .fold(0.0, |m, p| f64::max(m, p.1 as f64));
                    assert_eq!(got[k * refs.len() + s].to_bits(), want.to_bits());
                }
            }
        }
    }
}
