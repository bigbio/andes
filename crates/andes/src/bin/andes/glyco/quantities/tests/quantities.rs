//! Tests of glycopeptide quantification.

use super::mbr::{build_transfers, fit_rt_map, Alignment, MBR_MIN_ANCHORS, MBR_MIN_WINDOW_S};
use super::measure::{extracted_envelope, measure_identified};
use super::tables::{composition_label, output_path, write_rollup_tables};
use glyco::glycan::db::GlycanComp;

use super::*;

#[test]
fn the_extracted_envelope_starts_at_the_isotope_offset() {
    let full = Formula::averagine(3000.0).envelope(ENVELOPE_ISOTOPES);
    let mono = extracted_envelope(&full, 0);
    let m2 = extracted_envelope(&full, 2);
    assert!(m2[0] > mono[0], "{mono:?} {m2:?}");
    assert!((m2.iter().sum::<f64>() - 1.0).abs() < 1e-9);
    assert!(m2.len() <= MAX_ISOTOPES && m2.len() >= MIN_ISOTOPES);
}

#[test]
fn an_offset_past_the_envelope_is_a_single_isotope() {
    assert_eq!(extracted_envelope(&[0.5, 0.5], 9), vec![1.0]);
}

#[test]
fn a_large_glycopeptide_extracts_more_isotopes_than_a_peptide() {
    let peptide = Formula::peptide(b"PEPTIDEK".iter().map(|&r| (r, 0.0)));
    let mut glyco = Formula::peptide(b"IGNTTNA".iter().map(|&r| (r, 0.0)));
    glyco.add(Formula::glycan(8, 2, 0, 0, 0));
    let a = extracted_envelope(&peptide.envelope(ENVELOPE_ISOTOPES), 0);
    let b = extracted_envelope(&glyco.envelope(ENVELOPE_ISOTOPES), 0);
    assert!(b.len() >= a.len(), "{} vs {}", a.len(), b.len());
    assert!(b.len() > 4, "a 2 kDa glycopeptide needs more than 4: {b:?}");
}

#[test]
fn composition_labels_match_the_pin_tag() {
    let comp = GlycanComp::from_tree_str("(N(N(H(H(H(H))(H)(H))(H(H)))))").unwrap();
    assert_eq!(
        composition_label(&comp),
        format!(
            "HexNAc{}Hex{}Fuc{}NeuAc{}NeuGc{}",
            comp.hexnac, comp.hex, comp.fuc, comp.neuac, comp.neugc
        )
    );
}

#[test]
fn output_paths_sit_beside_the_pin() {
    let p = output_path(Path::new("/d/run.glyco.pin"), "quant");
    assert_eq!(p, PathBuf::from("/d/run.glyco.quant.tsv"));
}

fn settings() -> Settings {
    Settings {
        params: LfqParams::default(),
        transfer_fdr: 0.05,
        core: GlycanCore::NGlycan,
        debug: false,
        mbr: false,
        match_train: None,
    }
}

fn target(mono_mz: f64, charge: u8, anchor_rt: f64) -> Target {
    let mut f = Formula::peptide(b"IGNTTNA".iter().map(|&r| (r, 0.0)));
    f.add(Formula::glycan(8, 2, 0, 0, 0));
    Target {
        run_idx: 0,
        run: "r1".into(),
        spec_id: "scan=1_glyco_1_1".into(),
        scan: 1,
        evidence: Evidence::Identified,
        window_s: 60.0,
        peptide: "K.IGNTTNA.-".into(),
        glycan: "HexNAc2Hex8Fuc0NeuAc0NeuGc0".into(),
        topology: "(N(N(H(H(H(H))(H)(H))(H(H)))))".into(),
        site: "N3".into(),
        proteins: "P1".into(),
        charge,
        isotope_offset: 0,
        mono_mz,
        theo_mz: mono_mz,
        ppm: 0.0,
        anchor_rt,
        envelope: extracted_envelope(&f.envelope(ENVELOPE_ISOTOPES), 0),
        search_s: 0.0,
        donor: None,
    }
}

fn synthetic_run(t: &Target) -> quant::ms1::Ms1RunIndex {
    let z = t.charge as f64;
    let scans = (0..200)
        .map(|i| {
            let rt = i as f64;
            let g = (-0.5 * ((rt - 100.0) / 4.0).powi(2)).exp();
            let mut peaks: Vec<(f64, f32)> = Vec::new();
            for (k, e) in t.envelope.iter().enumerate() {
                let inten = (1e6 * e * g) as f32;
                if inten > 1.0 {
                    peaks.push((t.mono_mz + k as f64 * ISOTOPE / z, inten));
                }
            }
            peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
            model::scan::Ms1Scan {
                rt,
                peaks,
                ..model::scan::Ms1Scan::default()
            }
        })
        .collect();
    quant::ms1::Ms1RunIndex::new(scans)
}

#[test]
fn measures_a_synthetic_glycopeptide_envelope() {
    let t = target(800.0, 2, 97.0);
    let index = synthetic_run(&t);
    let m = measure_identified(&index, &t, &LfqParams::default());
    let f = m.feature.as_ref().expect("feature");
    assert!((f.apex_rt - 100.0).abs() <= 1.0, "{f:?}");
    assert!(f.cosine > 0.99, "{f:?}");
    assert!(f.n_isotopes_found >= MIN_ISOTOPES as u32, "{f:?}");
    let expected = 1e6 * 4.0 * (2.0 * std::f64::consts::PI).sqrt();
    assert!((f.area / expected - 1.0).abs() < 0.03, "{f:?}");
    assert!(m.decoy.is_none(), "{:?}", m.decoy);
}

#[test]
fn an_isotope_offset_extracts_the_observed_peak() {
    let mut t = target(800.0, 2, 97.0);
    t.mono_mz += ISOTOPE / 2.0;
    t.isotope_offset = 1;
    let mut f = Formula::peptide(b"IGNTTNA".iter().map(|&r| (r, 0.0)));
    f.add(Formula::glycan(8, 2, 0, 0, 0));
    t.envelope = extracted_envelope(&f.envelope(ENVELOPE_ISOTOPES), 1);
    let index = synthetic_run(&t);
    let m = measure_identified(&index, &t, &LfqParams::default());
    let f = m.feature.as_ref().expect("feature");
    assert!((f.apex_rt - 100.0).abs() <= 1.0, "{f:?}");
    assert!(f.cosine > 0.99, "{f:?}");
}

fn synthetic_run_mono_only(t: &Target) -> quant::ms1::Ms1RunIndex {
    let scans = (0..200)
        .map(|i| {
            let rt = i as f64;
            let g = (-0.5 * ((rt - 100.0) / 4.0).powi(2)).exp();
            let inten = (1e6 * t.envelope[0] * g) as f32;
            let peaks = if inten > 1.0 {
                vec![(t.mono_mz, inten)]
            } else {
                Vec::new()
            };
            model::scan::Ms1Scan {
                rt,
                peaks,
                ..model::scan::Ms1Scan::default()
            }
        })
        .collect();
    quant::ms1::Ms1RunIndex::new(scans)
}

#[test]
fn mono_intensity_survives_dropped_isotope_traces() {
    let t = target(800.0, 2, 97.0);
    let full = measure_identified(&synthetic_run(&t), &t, &LfqParams::default());
    let mono_only = measure_identified(&synthetic_run_mono_only(&t), &t, &LfqParams::default());
    let a = full.feature.as_ref().expect("full feature").area;
    let b = mono_only.feature.as_ref().expect("mono-only feature").area;
    assert!((a / b - 1.0).abs() < 0.03, "{a} vs {b}");
    let f = mono_only.feature.as_ref().unwrap();
    let expected = f.isotope_areas[0] / t.envelope[0];
    assert!((f.area / expected - 1.0).abs() < 1e-12, "{f:?}");
}

fn feature(area: f64, cosine: f32) -> FeatureQuant {
    FeatureQuant {
        apex_rt: 100.0,
        rt_start: 90.0,
        rt_stop: 110.0,
        area,
        apex_intensity: 1e6,
        isotope_areas: vec![area],
        cosine,
        n_isotopes_found: 4,
        rt_delta_s: 3.0,
        n_scans: 21,
    }
}

fn row(t: Target, area: f64, score: f32, q: f64) -> FeatureRow {
    FeatureRow {
        target: t,
        feature: Some(feature(area, 0.99)),
        decoy: None,
        score: Some(score),
        q: Some(q),
        match_score: None,
        decoy_match: None,
    }
}

fn read(path: &Path) -> Vec<Vec<String>> {
    let text = std::fs::read_to_string(path).unwrap();
    text.lines()
        .map(|l| l.split('\t').map(str::to_string).collect())
        .collect()
}

#[test]
fn the_rollup_dedupes_scans_sums_charges_and_drops_rejected_features() {
    let dir = tempfile::tempdir().unwrap();
    let pin = dir.path().join("run.glyco.pin");
    let mut a = row(target(800.0, 2, 97.0), 100.0, 0.9, 0.001);
    a.target.scan = 1;
    let mut b = row(target(800.0, 2, 97.0), 100.0, 0.5, 0.001);
    b.target.scan = 2;
    let mut c = row(target(800.0, 3, 97.0), 50.0, 0.8, 0.002);
    c.target.scan = 3;
    c.target.charge = 3;
    let mut rejected = row(target(800.0, 2, 97.0), 999.0, 0.99, 0.5);
    rejected.target.evidence = Evidence::Transferred;
    rejected.target.scan = 0;
    let rows = vec![a, b, c, rejected];
    let kept = write_rollup_tables(&pin, &rows, &["r1"], &settings()).unwrap();
    assert_eq!(kept.glycopeptides, 1);
    assert_eq!(kept.sites, 1);

    let q = read(&dir.path().join("run.glyco.quant.tsv"));
    assert_eq!(q.len(), 2);
    assert_eq!(
        q[0],
        vec![
            "Run",
            "Peptide",
            "Glycosite",
            "Glycan",
            "Proteins",
            "Charges",
            "NScans",
            "NTransferred",
            "BestCosine",
            "BestQValue",
            "Area"
        ]
    );
    let r = &q[1];
    assert_eq!(r[0], "r1");
    assert_eq!(r[1], "K.IGNTTNA.-");
    assert_eq!(r[5], "2;3");
    assert_eq!(r[6], "2");
    assert_eq!(r[7], "0", "an identified-only row has no transfers");
    assert!((r[10].parse::<f64>().unwrap() - 150.0).abs() < 1e-6);

    let s = read(&dir.path().join("run.glyco.site.tsv"));
    assert_eq!(s.len(), 2);
    assert!((s[1][7].parse::<f64>().unwrap() - 1.0).abs() < 1e-9);

    let p = read(&dir.path().join("run.glyco.protein.tsv"));
    assert_eq!(p.len(), 2);
    assert!((p[1][3].parse::<f64>().unwrap() - 150.0).abs() < 1e-6);
    assert!((p[1][4].parse::<f64>().unwrap() - 150.0).abs() < 1e-6);
}

#[test]
fn an_identified_feature_is_not_vetoed_by_its_decoy_q() {
    let dir = tempfile::tempdir().unwrap();
    let pin = dir.path().join("run.glyco.pin");
    let identified = row(target(800.0, 2, 97.0), 100.0, 0.9, 0.5);
    let mut transfer = row(target(800.0, 2, 98.0), 50.0, 0.9, 0.5);
    transfer.target.evidence = Evidence::Transferred;
    transfer.target.scan = 0;
    transfer.target.spec_id = "donor_mbr_1".into();
    let kept = write_rollup_tables(&pin, &[identified, transfer], &["r1"], &settings()).unwrap();
    assert_eq!(kept.glycopeptides, 1);

    let q = read(&dir.path().join("run.glyco.quant.tsv"));
    let r = &q[1];
    assert_eq!(r[7], "0", "the transfer is q-gated out: {r:?}");
    assert!((r[10].parse::<f64>().unwrap() - 100.0).abs() < 1e-6);
}

#[test]
fn the_matrices_pivot_runs_and_leave_missing_cells_empty() {
    let dir = tempfile::tempdir().unwrap();
    let pin = dir.path().join("run.glyco.pin");
    let mut a = row(target(800.0, 2, 97.0), 100.0, 0.9, 0.001);
    a.target.run = "r1".into();
    let mut b = row(target(800.0, 2, 97.0), 40.0, 0.9, 0.001);
    b.target.run = "r2".into();
    let mut other = row(target(800.0, 2, 97.0), 7.0, 0.9, 0.001);
    other.target.run = "r2".into();
    other.target.glycan = "HexNAc3Hex4Fuc0NeuAc0NeuGc0".into();
    let kept = write_rollup_tables(&pin, &[a, b, other], &["r1", "r2"], &settings()).unwrap();
    assert_eq!((kept.glycopeptides, kept.sites), (3, 2));

    let m = read(&dir.path().join("run.glyco.matrix.tsv"));
    assert_eq!(
        m[0],
        vec![
            "Peptide",
            "Glycosite",
            "Glycan",
            "Proteins",
            "NRuns",
            "r1",
            "r2"
        ]
    );
    assert_eq!(m[1][2], "HexNAc2Hex8Fuc0NeuAc0NeuGc0");
    assert_eq!(m[1][4], "2");
    assert!((m[1][5].parse::<f64>().unwrap() - 100.0).abs() < 1e-6);
    assert!((m[1][6].parse::<f64>().unwrap() - 40.0).abs() < 1e-6);
    assert_eq!(m[2][2], "HexNAc3Hex4Fuc0NeuAc0NeuGc0");
    assert_eq!(m[2][4], "1");
    assert_eq!(m[2][5], "", "r1 did not quantify this glycoform: {m:?}");
    assert!((m[2][6].parse::<f64>().unwrap() - 7.0).abs() < 1e-6);

    let p = read(&dir.path().join("run.glyco.protein_matrix.tsv"));
    assert_eq!(
        p[0],
        vec![
            "Protein",
            "r1_SummedArea",
            "r1_Top3",
            "r2_SummedArea",
            "r2_Top3"
        ]
    );
    assert_eq!(p[1][0], "P1");
    assert!((p[1][1].parse::<f64>().unwrap() - 100.0).abs() < 1e-6);
    assert!((p[1][3].parse::<f64>().unwrap() - 47.0).abs() < 1e-6);
    assert!((p[1][4].parse::<f64>().unwrap() - 47.0).abs() < 1e-6);
}

#[test]
fn the_rt_map_recovers_a_linear_shift() {
    let pairs: Vec<(f64, f64)> = (0..100)
        .map(|i| {
            let x = i as f64 * 3.0;
            (x, 4.0 + 1.1 * x)
        })
        .collect();
    let map = fit_rt_map(&pairs).expect("map");
    assert!(map.sigma_s < 1e-9, "{map:?}");
    for x in [0.0, 123.0, 297.0] {
        let want = 4.0 + 1.1 * x;
        assert!((map.forward(x) - want).abs() < 1e-6, "{x}: {map:?}");
    }
    let rt_ref = map.forward(123.0);
    assert!((map.inverse(rt_ref).unwrap() - 123.0).abs() < 1e-6);
}

#[test]
fn the_rt_map_follows_a_non_linear_warp() {
    let pairs: Vec<(f64, f64)> = (0..900)
        .map(|i| {
            let x = i as f64 * 4.0;
            let phase = x / 3600.0 * std::f64::consts::PI;
            (x, x + 40.0 * (phase).sin())
        })
        .collect();
    let map = fit_rt_map(&pairs).expect("map");
    assert!(map.sigma_s < 10.0, "{map:?}");
    for x in [400.0, 1800.0, 3200.0] {
        let want = x + 40.0 * (x / 3600.0 * std::f64::consts::PI).sin();
        assert!(
            (map.forward(x) - want).abs() < 5.0,
            "{x}: {} vs {want}",
            map.forward(x)
        );
        assert!((map.inverse(want).unwrap() - x).abs() < 5.0, "{x}");
    }
}

#[test]
fn the_rt_map_drops_gross_outliers() {
    let mut pairs: Vec<(f64, f64)> = (0..900)
        .map(|i| {
            let x = i as f64 * 3.0;
            (x, 4.0 + 1.1 * x)
        })
        .collect();
    for i in (0..900).step_by(180) {
        pairs[i].1 += 1000.0;
    }
    let map = fit_rt_map(&pairs).expect("map");
    for x in [100.0, 1400.0, 2600.0] {
        assert!(
            (map.forward(x) - (4.0 + 1.1 * x)).abs() < 1.0,
            "{x}: {map:?}"
        );
    }
    assert!(map.sigma_s < 0.5, "{map:?}");
}

fn identified_row(
    key: usize,
    run_idx: usize,
    run: &str,
    apex_rt: f64,
    q: f64,
    score: f32,
) -> FeatureRow {
    let mut t = target(800.0, 2, apex_rt);
    t.run_idx = run_idx;
    t.run = run.to_string();
    t.peptide = format!("K.PEPTIDE{key}K.-");
    t.spec_id = format!("{run}_scan{key}");
    let mut f = feature(100.0, 0.99);
    f.apex_rt = apex_rt;
    FeatureRow {
        target: t,
        feature: Some(f),
        decoy: None,
        score: Some(score),
        q: Some(q),
        match_score: None,
        decoy_match: None,
    }
}

#[test]
fn transfers_fill_the_runs_that_missed_the_glycopeptide() {
    let n = MBR_MIN_ANCHORS + 10;
    let mut rows = Vec::new();
    for k in 0..n {
        let rt = 100.0 + 2.0 * k as f64;
        rows.push(identified_row(k, 0, "r1", rt, 0.01, 10.0));
        rows.push(identified_row(k, 1, "r2", 0.8 * rt + 50.0, 0.01, 10.0));
    }
    rows.push(identified_row(n, 0, "r1", 500.0, 0.01, 20.0));
    rows.push(identified_row(n + 1, 0, "r1", 520.0, 0.01, 20.0));
    let mut rejected = identified_row(n + 1, 1, "r2", 466.0, 0.5, 20.0);
    rejected.score = None;
    rows.push(rejected);

    let settings = settings();
    let alignment = Alignment::build(&rows, &["r1", "r2"], &settings);
    let (transfers, counts) = build_transfers(&rows, &["r1", "r2"], &settings, &alignment);
    assert_eq!(counts.maps, 2, "both runs have an RT map");
    assert_eq!(counts.anchors, n);
    assert_eq!(counts.candidates, 1);
    let t = &transfers[0];
    assert_eq!(t.run_idx, 1);
    assert_eq!(t.evidence, Evidence::Transferred);
    assert_eq!(t.scan, 0);
    assert!(t.spec_id.ends_with("_mbr_r2"), "{}", t.spec_id);
    let expected = 0.8 * 500.0 + 50.0;
    assert!(
        (t.anchor_rt - expected).abs() < 1e-6,
        "{} vs {expected}",
        t.anchor_rt
    );
    assert!(
        (t.window_s - MBR_MIN_WINDOW_S).abs() < 1e-9,
        "{:?}",
        t.window_s
    );
}

#[test]
fn too_few_shared_glycopeptides_mean_no_map_and_no_transfers() {
    let n = MBR_MIN_ANCHORS - 1;
    let mut rows = Vec::new();
    for k in 0..n {
        let rt = 100.0 + 2.0 * k as f64;
        rows.push(identified_row(k, 0, "r1", rt, 0.01, 10.0));
        rows.push(identified_row(k, 1, "r2", rt + 10.0, 0.01, 10.0));
    }
    rows.push(identified_row(999, 0, "r1", 500.0, 0.01, 20.0));
    let settings = settings();
    let alignment = Alignment::build(&rows, &["r1", "r2"], &settings);
    let (transfers, counts) = build_transfers(&rows, &["r1", "r2"], &settings, &alignment);
    assert_eq!(counts.maps, 1, "only the reference run is mapped");
    assert_eq!(counts.candidates, 0);
    assert!(transfers.is_empty(), "{transfers:?}");
}

#[test]
fn a_transferred_row_is_counted_and_not_given_a_scan_number() {
    let dir = tempfile::tempdir().unwrap();
    let pin = dir.path().join("run.glyco.pin");
    let mut a = row(target(800.0, 2, 97.0), 100.0, 0.9, 0.001);
    a.target.run = "r1".into();
    a.target.scan = 7;
    let mut b = row(target(800.0, 2, 97.0), 40.0, 0.9, 0.001);
    b.target.run = "r2".into();
    b.target.scan = 0;
    b.target.evidence = Evidence::Transferred;
    b.target.spec_id = "r1_scan7_glyco_1_1_mbr_r2".into();
    write_rollup_tables(&pin, &[a, b], &["r1", "r2"], &settings()).unwrap();

    let q = read(&dir.path().join("run.glyco.quant.tsv"));
    let r2 = q.iter().find(|r| r[0] == "r2").expect("r2 row");
    assert_eq!(r2[6], "0", "a transfer has no identifying scan: {r2:?}");
    assert_eq!(r2[7], "1");
    assert!((r2[10].parse::<f64>().unwrap() - 40.0).abs() < 1e-6);
    let r1 = q.iter().find(|r| r[0] == "r1").expect("r1 row");
    assert_eq!(r1[6], "1");
    assert_eq!(r1[7], "0");
}
