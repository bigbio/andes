use super::*;

fn row(id: &str, decoy: bool, delta: f64, site: usize, seq: &str) -> PinOffset {
    PinOffset {
        psm_id: id.to_string(),
        decoy,
        delta,
        calc_mass: 1500.0,
        refined: false,
        site,
        sequence: seq.to_string(),
        run: "r1".to_string(),
    }
}

fn report(rows: &[PinOffset]) -> OffsetReport {
    let accepted: HashSet<String> = rows
        .iter()
        .filter(|r| !r.decoy)
        .map(|r| r.psm_id.clone())
        .collect();
    let target_scores: HashMap<String, f64> = rows
        .iter()
        .filter(|r| !r.decoy)
        .map(|r| (r.psm_id.clone(), 1.0))
        .collect();
    let decoy_scores: HashMap<String, f64> = rows
        .iter()
        .filter(|r| r.decoy)
        .map(|r| (r.psm_id.clone(), 2.0))
        .collect();
    build_report(
        rows,
        &accepted,
        &target_scores,
        &decoy_scores,
        model::unimod::bundled(),
    )
}

fn zero_rows(n: usize, ppm: f64) -> Vec<PinOffset> {
    (0..n)
        .map(|i| {
            let jitter = 0.0001 * (i % 5) as f64 - 0.0002;
            row(&format!("z{i}"), false, ppm * 1e-6 * 1500.0 + jitter, 0, "PEPTIDEK")
        })
        .collect()
}

#[test]
fn stripped_sequence_drops_flanks_and_mods() {
    assert_eq!(stripped_sequence("K.PEM[+15.995]TIDEK.R"), "PEMTIDEK");
    assert_eq!(stripped_sequence("-.ACDK.-"), "ACDK");
    assert_eq!(stripped_sequence("PEPTIDE"), "PEPTIDE");
}

#[test]
fn unimod_specificities_keep_hidden_flags_positions_and_label_classes() {
    let db = model::unimod::bundled();
    let ox = db.by_accession("UNIMOD:35").expect("Oxidation");
    assert!(ox.specificities.iter().any(|s| s.site == "M" && !s.hidden), "{:?}", ox.specificities);
    assert!(ox.specificities.iter().any(|s| s.hidden), "{:?}", ox.specificities);
    assert!(!ox.isotopic_label);
    let pyro = db.by_accession("UNIMOD:28").expect("Gln->pyro-Glu");
    assert!(
        pyro.specificities
            .iter()
            .any(|s| s.site == "Q" && s.position.needs_n_terminus()),
        "{:?}",
        pyro.specificities
    );
    assert!(db.by_accession("UNIMOD:994").expect("Label:15N(1)").isotopic_label);
}

#[test]
fn an_n_terminal_specificity_supports_only_an_n_terminal_residue() {
    let pyro = single_candidates(model::unimod::bundled())
        .into_iter()
        .find(|c| c.label.starts_with("Gln->pyro-Glu (UNIMOD:28)"))
        .expect("pyro-Glu candidate");
    let first_q = site_support(&pyro, &[("Q".to_string(), true, false)]);
    let middle_q = site_support(&pyro, &[("Q".to_string(), false, false)]);
    assert!(first_q > middle_q, "{first_q} {middle_q}");
    assert!((middle_q - 0.25).abs() < 1e-12, "{middle_q}");
}

#[test]
fn each_run_is_calibrated_on_its_own_zero_offset_psms() {
    let true_delta = 79.96633;
    let mut rows: Vec<PinOffset> = Vec::new();
    for (run, ppm) in [("a", 5.0), ("b", -5.0)] {
        for i in 0..30 {
            let jitter = 0.0001 * (i % 5) as f64 - 0.0002;
            rows.push(PinOffset {
                run: run.to_string(),
                calc_mass: 2000.0,
                ..row(&format!("{run}z{i}"), false, ppm * 1e-6 * 2000.0 + jitter, 0, "PEPTIDEK")
            });
        }
        for i in 0..10 {
            let observed = true_delta + ppm * 1e-6 * (2000.0 + true_delta) + 0.0001 * (i % 3) as f64;
            rows.push(PinOffset {
                run: run.to_string(),
                calc_mass: 2000.0,
                ..row(&format!("{run}p{i}"), false, observed, 2, "ASPEK")
            });
        }
    }
    rows.push(PinOffset {
        run: "b".to_string(),
        refined: true,
        ..row("bref", false, 15.9949, 3, "PEMTIDEK")
    });
    let r = report(&rows);
    assert_eq!(r.run_ppm_shifts.len(), 2);
    assert!((r.run_ppm_shifts[0].1 - 5.0).abs() < 0.2, "{:?}", r.run_ppm_shifts);
    assert!((r.run_ppm_shifts[1].1 + 5.0).abs() < 0.2, "{:?}", r.run_ppm_shifts);
    let phospho: Vec<&OffsetPeak> = r
        .peaks
        .iter()
        .filter(|p| (p.mass - 79.966).abs() < 0.05)
        .collect();
    assert_eq!(phospho.len(), 1, "{:?}", r.peaks);
    assert_eq!(phospho[0].psms, 20);
    assert!((phospho[0].mass - true_delta).abs() < 0.001, "{}", phospho[0].mass);
    let refined = r.assignments.iter().find(|a| a.psm_id == "bref").expect("refined row");
    assert!((refined.delta - 15.9949).abs() < 1e-12, "{}", refined.delta);
}

#[test]
fn methionine_peak_at_plus_16_is_oxidation() {
    let mut rows = zero_rows(30, 0.0);
    rows.extend((0..12).map(|i| {
        row(&format!("o{i}"), false, 15.9949 + 0.0003 * (i % 3) as f64 - 0.0003, 3, "PEMTIDEK")
    }));
    let r = report(&rows);
    let ox = r
        .peaks
        .iter()
        .find(|p| (p.mass - 15.9949).abs() < 0.01)
        .expect("+16 peak");
    assert_eq!(ox.psms, 12);
    assert_eq!(ox.residues, vec![("M".to_string(), 12)]);
    let top = &ox.explanations[0];
    assert!(top.label.starts_with("Oxidation (UNIMOD:35)"), "{:?}", ox.explanations);
    assert!(top.posterior > 0.5, "{:?}", ox.explanations);
    let zero = r.peaks.iter().find(|p| p.mass.abs() < 0.01).expect("zero peak");
    assert_eq!(zero.explanations[0].label, "unmodified");
}

#[test]
fn a_ppm_calibration_offset_is_removed_before_the_peaks_are_found() {
    let ppm = 4.0;
    let mut rows = zero_rows(40, ppm);
    rows.extend((0..12).map(|i| {
        let true_delta = 79.96633;
        let observed = true_delta + ppm * 1e-6 * (1500.0 + true_delta) + 0.0001 * (i % 3) as f64;
        row(&format!("p{i}"), false, observed, 2, "ASPEK")
    }));
    let r = report(&rows);
    assert!((r.ppm_shift - ppm).abs() < 0.2, "{}", r.ppm_shift);
    let zero = r.peaks.iter().find(|p| p.mass.abs() < 0.01).expect("zero peak");
    assert!(zero.mass.abs() < 0.0005, "{}", zero.mass);
    let phospho = r.peaks.iter().find(|p| (p.mass - 79.966).abs() < 0.01).expect("phospho");
    assert!((phospho.mass - 79.96633).abs() < 0.001, "{}", phospho.mass);
    assert!(phospho.explanations[0].label.starts_with("Phospho (UNIMOD:21)"), "{:?}", phospho.explanations);
}

#[test]
fn fragment_refined_offsets_keep_their_value() {
    let ppm = 4.0;
    let mut rows = zero_rows(40, ppm);
    rows.extend((0..12).map(|i| PinOffset {
        refined: true,
        ..row(&format!("o{i}"), false, 15.9949 + 0.0001 * (i % 3) as f64, 3, "PEMTIDEK")
    }));
    let r = report(&rows);
    let ox = r.peaks.iter().find(|p| (p.mass - 15.99).abs() < 0.05).expect("+16 peak");
    assert!((ox.mass - 15.995).abs() < 0.0005, "{}", ox.mass);
}

#[test]
fn deamidation_and_an_isotope_step_split_at_the_density_valley() {
    let mut sorted: Vec<f64> = Vec::new();
    for i in 0..60 {
        sorted.push(0.9840 + 0.0004 * ((i % 7) as f64 - 3.0));
        sorted.push(1.0034 + 0.0004 * ((i % 7) as f64 - 3.0));
    }
    for i in 0..8 {
        sorted.push(0.9900 + 0.0014 * i as f64);
    }
    sorted.sort_by(|a, b| a.total_cmp(b));
    let groups = offset_peaks(&sorted, 0.002);
    let big: Vec<_> = groups.iter().filter(|g| g.len() >= 40).collect();
    assert_eq!(big.len(), 2, "{groups:?}");
    let centre = |g: &std::ops::Range<usize>| median(&mut sorted[g.clone()].to_vec());
    assert!((centre(big[0]) - 0.9840).abs() < 0.002);
    assert!((centre(big[1]) - 1.0034).abs() < 0.002);
}

#[test]
fn a_mixed_peak_reports_both_explanations_in_proportion() {
    let mut rows = zero_rows(30, 0.0);
    rows.extend((0..30).map(|i| {
        row(&format!("n{i}"), false, 0.9840 + 0.0005 * (i % 5) as f64 - 0.001, 3, "PENTIDEK")
    }));
    rows.extend((0..20).map(|i| {
        row(&format!("j{i}"), false, 1.0034 + 0.0005 * (i % 5) as f64 - 0.001, 0, "PEPTIDEK")
    }));
    rows.extend((0..12).map(|i| row(&format!("b{i}"), false, 0.9870 + 0.0013 * i as f64, 0, "PEPTIDEK")));
    let r = report(&rows);
    let n0_peak = r
        .assignments
        .iter()
        .find(|a| a.psm_id == "n0")
        .and_then(|a| a.peak)
        .expect("n0 is in a peak");
    let peak = r.peaks.iter().find(|p| p.id == n0_peak).expect("its peak");
    let label_of = |e: &Explanation| e.label.clone();
    let labels: Vec<String> = peak.explanations.iter().map(label_of).collect();
    assert!(labels.iter().any(|l| l.starts_with("Deamidated (UNIMOD:7)")), "{:?}", peak.explanations);
    let deamidated: Vec<&PsmAssignment> = r
        .assignments
        .iter()
        .filter(|a| a.psm_id.starts_with('n'))
        .collect();
    assert!(deamidated.iter().all(|a| a.explanation.starts_with("Deamidated (UNIMOD:7)")), "{deamidated:?}");
}

#[test]
fn refined_offsets_are_judged_by_the_spread_of_refined_peaks() {
    let mut rows = zero_rows(30, 0.0);
    rows.extend((0..40).map(|i| PinOffset {
        refined: true,
        ..row(&format!("p{i}"), false, 79.9663 + 0.0012 * ((i % 9) as f64 - 4.0), 2, "ASPEK")
    }));
    let r = report(&rows);
    assert!(r.refined_precision > r.precision, "{} {}", r.refined_precision, r.precision);
    let phospho = r.peaks.iter().find(|p| (p.mass - 79.966).abs() < 0.01).expect("phospho");
    assert!(phospho.explanations[0].label.starts_with("Phospho (UNIMOD:21)"), "{:?}", phospho.explanations);
    assert!(phospho.explanations[0].posterior > 0.8, "{:?}", phospho.explanations);
}

#[test]
fn a_peak_one_neutron_up_is_an_isotope_error() {
    let rows: Vec<PinOffset> = (0..10)
        .map(|i| row(&format!("i{i}"), false, 1.0034 + 0.0002 * (i % 2) as f64, 0, "PEPTIDEK"))
        .collect();
    let r = report(&rows);
    assert_eq!(r.peaks.len(), 1);
    assert!(
        r.peaks[0]
            .explanations
            .iter()
            .any(|e| e.label == "isotope error +1"),
        "{:?}",
        r.peaks[0].explanations
    );
}

#[test]
fn decoys_above_the_threshold_are_counted_in_their_peak() {
    let mut rows: Vec<PinOffset> = (0..5)
        .map(|i| row(&format!("t{i}"), false, 42.0106, 1, "AKDE"))
        .collect();
    rows.push(row("d0", true, 42.0107, 1, "EDKA"));
    rows.push(row("d1", true, 80.0, 1, "EDKA"));
    let r = report(&rows);
    assert_eq!(r.peaks.len(), 1);
    assert_eq!(r.peaks[0].decoys, 1);
}

#[test]
fn neighbouring_decoy_windows_do_not_overlap() {
    let windows = decoy_windows(&[(10.0, 10.0055), (10.0065, 10.012)], 0.002);
    assert!(windows[0].1 <= windows[1].0, "{windows:?}");
    assert!((windows[0].1 - 10.006).abs() < 1e-9, "{windows:?}");
    assert!((windows[0].0 - 9.998).abs() < 1e-9, "{windows:?}");
    assert!((windows[1].1 - 10.014).abs() < 1e-9, "{windows:?}");
}

#[test]
fn singletons_are_not_peaks_but_keep_an_assignment_row() {
    let rows = vec![
        row("a", false, 10.0, 0, "PEPTIDE"),
        row("b", false, 55.5, 0, "PEPTIDE"),
    ];
    let r = report(&rows);
    assert!(r.peaks.is_empty());
    assert_eq!(r.assignments.len(), 2);
    assert!(r.assignments.iter().all(|a| a.peak.is_none()));
}

#[test]
fn report_tables_round_trip_through_parquet() {
    let mut rows = zero_rows(30, 0.0);
    rows.extend((0..12).map(|i| {
        row(&format!("o{i}"), false, 15.9949 + 0.0003 * (i % 3) as f64 - 0.0003, 3, "PEMTIDEK")
    }));
    rows.push(row("lone", false, 55.5, 0, "PEPTIDE"));
    let r = report(&rows);
    let dir = tempfile::tempdir().unwrap();
    let peaks = dir.path().join("x.open_offsets.parquet");
    let psms = dir.path().join("x.open_psms.parquet");
    write_peaks_table(&peaks, &r).unwrap();
    write_psms_table(&psms, &r).unwrap();
    let t = read_peaks_table(&peaks).unwrap();
    let n: usize = t.batches.iter().map(|b| b.num_rows()).sum();
    assert_eq!(n, r.peaks.len());
    assert_eq!(t.meta.get(PRECISION_KEY).map(String::as_str), Some(r.precision.to_string().as_str()));
    let t = read_psms_table(&psms).unwrap();
    let n: usize = t.batches.iter().map(|b| b.num_rows()).sum();
    assert_eq!(n, r.assignments.len());
    assert!(read_peaks_table(&psms).is_err());
}

#[test]
fn summary_block_counts_each_psm_once_and_folds_isotope_picks() {
    let mut rows = zero_rows(30, 0.0);
    rows.extend((0..12).map(|i| {
        row(&format!("o{i}"), false, 15.9949 + 0.0003 * (i % 3) as f64 - 0.0003, 3, "PEMTIDEK")
    }));
    rows.push(row("lone", false, 55.5, 0, "PEPTIDE"));
    let r = report(&rows);
    let s = summary_block(&r);
    assert!(s.contains("PSMs                         : 43"), "{s}");
    assert!(s.contains("unmodified (zero offset)   : 30"), "{s}");
    assert!(s.contains("outside every peak         : 1"), "{s}");
    assert!(s.contains("Oxidation (UNIMOD:35)"), "{s}");
    assert!(s.contains("Largest peaks:"), "{s}");
    assert_eq!(modification_of("Phospho (UNIMOD:21) +1 isotope"), "Phospho (UNIMOD:21)");
    assert_eq!(modification_of("isotope error -1"), ISOTOPE_PICK);
}

#[test]
fn pin_columns_are_read_by_name() {
    let text = "SpecId\tLabel\tScanNr\tExpMass\tCalcMass\tisotope_error\tDeltaMass\tAbsDeltaMass\tOpenSite\tPeptide\tProteins\n\
                s1\t1\t5\t1015.99\t1000.0\t0\t15.9949\t15.9949\t3\tK.PEMTIDEK.R\tP1\n\
                s2\t-1\t6\t1000.001\t1000.0\t0\t0.0010\t0.0010\t0\tK.KEDITMEP.R\tdecoy_P1\n";
    let rows = parse_pin_offsets(std::io::Cursor::new(text)).unwrap();
    assert_eq!(rows.len(), 2);
    assert_eq!(rows[0].sequence, "PEMTIDEK");
    assert_eq!(rows[0].site, 3);
    assert_eq!(rows[0].calc_mass, 1000.0);
    assert!(rows[0].refined);
    assert!(rows[1].decoy);
    assert!(!rows[1].refined);
}
