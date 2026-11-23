use super::*;

#[test]
fn statistics_log_is_written_through_a_partial_file() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("run.statistics.log");
    let stats = RunStatistics {
        precursor_tolerance: "10 ppm".to_string(),
        precursor_calibration: "Auto".to_string(),
        fragment_tolerance: "0.02 Da".to_string(),
        node_score_bin: "0.5 Da".to_string(),
        spectra_with_match: 3,
        target_psms: 2,
        decoy_psms: 1,
        unmodified_target_psms: 2,
        ptm_counts: Vec::new(),
    };
    write_statistics_log(&stats, None, &path).unwrap();
    assert_eq!(std::fs::read_to_string(&path).unwrap(), stats.render());
    assert!(!crate::atomic::partial(&path).exists());
}

#[test]
fn render_lists_ptms_desc_and_unmodified() {
    let stats = RunStatistics {
        precursor_tolerance: "Symmetric(10.0 ppm)".to_string(),
        precursor_calibration: "Auto".to_string(),
        fragment_tolerance: "0.5 Da".to_string(),
        node_score_bin: "0.5 Da".to_string(),
        spectra_with_match: 100,
        target_psms: 60,
        decoy_psms: 40,
        unmodified_target_psms: 25,
        ptm_counts: vec![
            ("Carbamidomethyl".to_string(), 30),
            ("Oxidation".to_string(), 10),
        ],
    };
    let out = stats.render();
    assert!(out.contains("Final fragment tolerance  : 0.5 Da"));
    assert!(out.contains("Carbamidomethyl : 30"));
    assert!(out.contains("Oxidation       : 10"));
    assert!(out.contains("(unmodified)    : 25"));
    let cam = out.find("Carbamidomethyl").unwrap();
    let ox = out.find("Oxidation").unwrap();
    assert!(cam < ox, "PTMs should be sorted count-descending");
}

#[test]
fn high_res_run_reports_matching_tolerance_apart_from_node_score_bin() {
    let params = SearchParams::default_tryptic(
        model::aa_set::AminoAcidSetBuilder::new_standard()
            .build()
            .unwrap(),
    );
    let stats = RunStatistics::compute(
        &[],
        &[],
        &params,
        &Tolerance::Ppm(20.0),
        &Tolerance::Ppm(9.87654),
        &Tolerance::Da(0.5),
    );
    let out = stats.render();
    assert!(out.contains(
        "Final fragment tolerance  : 20 ppm (retrieval and model features); \
         9.877 ppm (PIN match features, measured)"
    ));
    assert!(out.contains("Node-score bin width      : 0.5 Da"));
    assert!(!out.contains("Final fragment tolerance  : 0.5 Da"));

    let low = RunStatistics::compute(
        &[],
        &[],
        &params,
        &Tolerance::Da(0.5),
        &Tolerance::Da(0.5),
        &Tolerance::Da(0.5),
    );
    assert!(low
        .render()
        .contains("Final fragment tolerance  : 0.5 Da (retrieval and match features)"));
}

fn fact(run: usize, charge: u8, pep: &str, ppm: f64, missed: u32, mods: &[&str]) -> PsmFacts {
    PsmFacts {
        run,
        charge,
        peptidoform: pep.to_string(),
        sequence: pep.chars().filter(|c| c.is_ascii_uppercase()).collect(),
        ppm,
        ppm_calibrated: ppm,
        spectrum_charge: Some(charge),
        secondary: false,
        missed_cleavages: missed,
        variable_mods: mods.iter().map(|m| m.to_string()).collect(),
        fixed_carried: Vec::new(),
        fixed_eligible: Vec::new(),
    }
}

fn names(n: usize) -> Vec<String> {
    (1..=n).map(|i| format!("run{i}")).collect()
}

fn tally(per_run: Vec<usize>, charge2: usize, unknown: usize) -> SpectraTally {
    SpectraTally {
        per_run,
        by_charge: BTreeMap::from([(2, charge2)]),
        unknown_charge: unknown,
    }
}

fn ctx() -> SummaryContext<'static> {
    SummaryContext {
        rescorer: "Percolator",
        cut: "q ≤ 0.01",
        decoys_at_cut: 0,
        protein_groups: None,
        fixed: &[],
        calibration_shift_ppm: 0.0,
    }
}

fn summarise(
    tally: &SpectraTally,
    runs: usize,
    facts: &[PsmFacts],
    unlinked: usize,
) -> IdentificationSummary {
    IdentificationSummary::compute(tally, &names(runs), facts, unlinked, &ctx())
}

fn single_run_summary() -> IdentificationSummary {
    let tally = SpectraTally {
        per_run: vec![20],
        by_charge: BTreeMap::from([(2, 12), (3, 8)]),
        unknown_charge: 0,
    };
    let facts = vec![
        fact(0, 2, "AAK", 1.0, 0, &["Carbamidomethyl"]),
        fact(0, 2, "AAK", 3.0, 0, &["Carbamidomethyl"]),
        fact(0, 3, "AAK", -1.0, 1, &["Carbamidomethyl"]),
        fact(0, 2, "BBR", 2.0, 0, &[]),
    ];
    let ctx = SummaryContext {
        decoys_at_cut: 1,
        protein_groups: Some((0.01, 5)),
        ..ctx()
    };
    IdentificationSummary::compute(&tally, &names(1), &facts, 0, &ctx)
}

#[test]
fn summary_counts_rates_and_distinct_levels() {
    let s = single_run_summary();
    assert_eq!(s.spectra, 20);
    assert_eq!(s.psms, 4);
    assert_eq!((s.precursors, s.peptidoforms, s.sequences), (3, 2, 2));
    assert_eq!(s.decoys_at_cut, 1);
    assert_eq!(s.missed_cleavages, vec![(0, 3), (1, 1)]);
    assert_eq!(s.ptm_counts, vec![("Carbamidomethyl".to_string(), 3)]);
    assert_eq!(s.unmodified_psms, 1);
    let z2 = &s.by_charge[0];
    assert_eq!((z2.label.as_str(), z2.spectra, z2.psms), ("2+", 12, 3));
    assert!((z2.rate.unwrap() - 0.25).abs() < 1e-12);
    let (m, q1, q3) = s.mass_error_ppm.unwrap();
    assert!((m - 1.5).abs() < 1e-12 && (q1 - 0.5).abs() < 1e-12 && (q3 - 2.25).abs() < 1e-12);
}

#[test]
fn single_run_render_has_no_run_table() {
    let out = single_run_summary().render();
    assert!(out.contains("andes identification summary (Percolator, q ≤ 0.01)"));
    assert!(out.contains("  Spectra searched           : 20\n"));
    assert!(!out.contains("MS2"));
    assert!(out.contains("  Identified PSMs            : 4  (20.0 % of spectra)\n"));
    assert!(out.contains("  Precursors / peptidoforms  : 3 / 2   sequences: 2\n"));
    assert!(out.contains("  Protein groups (q ≤ 0.01)  : 5\n"));
    assert!(out.contains("  Decoys at the cut          : 1\n"));
    assert!(out.contains("Precursor mass error (raw)"));
    assert!(out.contains("median +1.5 ppm (IQR +0.5 … +2.3)"));
    assert!(!out.contains("after calibration shift"));
    assert!(!out.contains("Fixed:"));
    assert!(out.contains("  Missed cleavages           : 0: 75.0 %  1: 25.0 %\n"));
    assert!(out.contains("Carbamidomethyl : 3  (75.0 %)"));
    assert!(out.contains("(unmodified)    : 1  (25.0 %)"));
    assert!(!out.contains("By run"));
    assert!(!out.contains("per-charge rates are not shown"));
    assert!(out.contains("25.0 %"));
}

#[test]
fn heading_names_the_rescoring_backend() {
    let native = SummaryContext {
        rescorer: "native rescorer",
        ..ctx()
    };
    let s = IdentificationSummary::compute(&tally(vec![5], 5, 0), &names(1), &[], 0, &native);
    let out = s.render();
    assert!(out.contains("andes identification summary (native rescorer, q ≤ 0.01)"));
    assert!(!out.contains("Percolator"));
}

#[test]
fn unknown_charge_spectra_get_a_row_and_a_note_and_no_rates() {
    let mut psm3 = fact(0, 3, "AAK", 0.0, 0, &[]);
    psm3.spectrum_charge = None;
    let facts = vec![fact(0, 2, "AAK", 0.0, 0, &[]), psm3];
    let s = summarise(&tally(vec![10], 6, 4), 1, &facts, 0);
    assert_eq!(s.spectra, 10);
    assert!(!s.charge_rates_comparable);
    let out = s.render();
    assert!(out.contains("  unknown"));
    assert!(out.contains("per-charge rates are not shown"));
    assert!(out.contains("  3+"));
    assert!(!out.contains("16.7 %"), "per-charge rates are withheld");
    assert!(!out.contains("Protein groups"));
    assert!(out.contains("median +0.0 ppm"));
}

#[test]
fn a_psm_at_another_charge_or_a_secondary_withholds_the_charge_rates() {
    let base = vec![
        fact(0, 2, "AAK", 0.0, 0, &[]),
        fact(0, 2, "BBK", 0.0, 0, &[]),
    ];
    let ok = summarise(&tally(vec![4], 4, 0), 1, &base, 0);
    assert!(ok.charge_rates_comparable);
    assert_eq!(ok.by_charge[0].rate, Some(0.5));

    let mut other = base.clone();
    other[1].charge = 3;
    let s = summarise(&tally(vec![4], 4, 0), 1, &other, 0);
    assert!(!s.charge_rates_comparable);
    assert!(s.by_charge.iter().all(|r| r.rate.is_none()));
    assert!(s.by_charge.iter().any(|r| r.psms > 0));

    let mut secondary = base.clone();
    secondary[1].secondary = true;
    let s = summarise(&tally(vec![4], 4, 0), 1, &secondary, 0);
    assert!(!s.charge_rates_comparable);
    assert!(s.render().contains("per-charge rates are not shown"));

    assert!(!summarise(&tally(vec![4], 4, 0), 1, &base, 1).charge_rates_comparable);
}

#[test]
fn no_rate_exceeds_one_hundred_percent() {
    let facts = vec![
        fact(0, 2, "AAK", 0.0, 0, &[]),
        fact(0, 2, "BBK", 0.0, 0, &[]),
        fact(0, 2, "CCK", 0.0, 0, &[]),
    ];
    let one = summarise(&tally(vec![2], 2, 0), 1, &facts, 0);
    assert_eq!(one.by_charge[0].rate, None);
    let out = one.render();
    assert!(out.contains("Identified PSMs            : 3  (- of spectra)"));
    assert!(!out.contains("150.0 %"));
    let two = summarise(&tally(vec![2, 5], 7, 0), 2, &facts, 0);
    assert_eq!(two.by_run[0].rate, None);
    assert!(two.by_run[1].rate.is_none_or(|v| v <= 1.0));
    assert!(two
        .by_charge
        .iter()
        .all(|r| r.rate.is_none_or(|v| v <= 1.0)));
}

#[test]
fn run_table_sums_and_flags_a_run_far_below_the_others() {
    let mut facts = Vec::new();
    for (run, n) in [(0usize, 50usize), (1, 48), (2, 52), (3, 5)] {
        for i in 0..n {
            facts.push(fact(run, 2, &format!("PEP{i}K"), 0.0, 0, &[]));
        }
    }
    let s = summarise(&tally(vec![100, 100, 100, 100], 400, 0), 4, &facts, 0);
    assert_eq!(s.by_run.iter().map(|r| r.psms).sum::<usize>(), s.psms);
    assert_eq!(s.by_run.iter().map(|r| r.spectra).sum::<usize>(), s.spectra);
    let flagged: Vec<&str> = s
        .by_run
        .iter()
        .filter(|r| r.z.is_some_and(|z| z <= LOW_RATE_Z))
        .map(|r| r.label.as_str())
        .collect();
    assert_eq!(flagged, vec!["run4"]);
    let out = s.render();
    assert!(out.contains("By run"));
    assert!(out.contains("run4") && out.contains("  <\n"));
    assert!(out.contains("changes no result"));
}

#[test]
fn similar_runs_and_short_run_lists_are_not_flagged() {
    let facts: Vec<PsmFacts> = [(0usize, 50usize), (1, 49), (2, 47)]
        .iter()
        .flat_map(|&(r, n)| (0..n).map(move |i| fact(r, 2, &format!("P{i}K"), 0.0, 0, &[])))
        .collect();
    let s = summarise(&tally(vec![100, 100, 100], 300, 0), 3, &facts, 0);
    assert!(s.by_run.iter().all(|r| r.z.is_none_or(|z| z > LOW_RATE_Z)));
    assert!(!s.render().contains("changes no result"));
    let two = summarise(
        &tally(vec![100, 100], 200, 0),
        2,
        &[fact(0, 2, "AAK", 0.0, 0, &[])],
        0,
    );
    assert!(two.by_run.iter().all(|r| r.z.is_none()));
}

#[test]
fn empty_accepted_set_renders_without_panicking() {
    let s = summarise(&tally(vec![5], 5, 0), 1, &[], 0);
    assert_eq!(s.psms, 0);
    assert!(s.mass_error_ppm.is_none());
    let out = s.render();
    assert!(out.contains("Identified PSMs            : 0  (0.0 % of spectra)"));
    assert!(!out.contains("Precursor mass error"));
    assert!(!out.contains("Missed cleavages"));
    assert!(out.contains("(unmodified)"));
}

#[test]
fn ptms_sort_by_count_then_name_and_numbers_get_separators() {
    let facts = vec![
        fact(0, 2, "AK", 0.0, 0, &["Oxidation"]),
        fact(0, 2, "BK", 0.0, 0, &["Acetyl"]),
        fact(0, 2, "CK", 0.0, 0, &["Acetyl", "Oxidation"]),
        fact(0, 2, "DK", 0.0, 0, &["Zeta"]),
    ];
    let s = summarise(&tally(vec![2_000_000], 2_000_000, 0), 1, &facts, 0);
    let order: Vec<&str> = s.ptm_counts.iter().map(|(n, _)| n.as_str()).collect();
    assert_eq!(order, vec!["Acetyl", "Oxidation", "Zeta"]);
    assert!(s.render().contains("2,000,000"));
    assert_eq!(thousands(0), "0");
    assert_eq!(thousands(999), "999");
    assert_eq!(thousands(1000), "1,000");
    assert_eq!(thousands(1234567), "1,234,567");
}

#[test]
fn unlinked_psms_are_counted_and_called_out() {
    let facts = vec![fact(0, 2, "AAK", 0.0, 0, &[])];
    let s = summarise(&tally(vec![10], 10, 0), 1, &facts, 2);
    assert_eq!(s.psms, 3);
    assert!(s
        .render()
        .contains("2 identified PSMs could not be matched"));
}

#[test]
fn spectra_tally_splits_runs_and_charges() {
    let spec = |z: Option<i32>| Spectrum {
        title: String::new(),
        precursor_mz: 500.0,
        precursor_intensity: None,
        precursor_charge: z,
        rt_seconds: None,
        scan: None,
        peaks: Vec::new(),
        activation_method: None,
        isolation_lower_offset: None,
        isolation_upper_offset: None,
        isolation_target_mz: None,
        mobility_window: None,
    };
    let spectra = vec![spec(Some(2)), spec(None), spec(Some(3)), spec(Some(2))];
    let t = SpectraTally::of(&spectra, &[2, 4]);
    assert_eq!(t.per_run, vec![2, 2]);
    assert_eq!(t.by_charge, BTreeMap::from([(2, 2), (3, 1)]));
    assert_eq!(t.unknown_charge, 1);
    assert_eq!(t.total(), 4);
}

#[test]
fn rescored_report_replaces_the_pre_fdr_ptm_list_and_appends_the_summary() {
    let stats = RunStatistics {
        precursor_tolerance: "10 ppm".to_string(),
        precursor_calibration: "Auto".to_string(),
        fragment_tolerance: "0.02 Da".to_string(),
        node_score_bin: "0.5 Da".to_string(),
        spectra_with_match: 3,
        target_psms: 2,
        decoy_psms: 1,
        unmodified_target_psms: 2,
        ptm_counts: Vec::new(),
    };
    let plain = stats.render();
    assert_eq!(plain, stats.render_with(None));
    assert!(plain.contains("PTM report (rank-1 target PSMs"));
    let summary = single_run_summary();
    let out = stats.render_with(Some(&summary));
    assert!(!out.contains("PTM report"));
    assert!(out.contains("Rank-1 PSMs (pre-FDR)     : 2 target, 1 decoy"));
    assert!(out.ends_with(&summary.render()));
}

fn spec(name: &str, residue: Option<char>, location: ModLocation) -> FixedModSpec {
    FixedModSpec {
        name: name.to_string(),
        residue,
        location,
    }
}

fn peptide<'a>(seq: &str, mods: &[(usize, &'a str)]) -> Vec<(u8, Option<&'a str>)> {
    seq.bytes()
        .enumerate()
        .map(|(i, r)| (r, mods.iter().find(|(p, _)| *p == i).map(|(_, m)| *m)))
        .collect()
}

#[test]
fn components_keep_declared_names_with_a_plus_whole() {
    let declared = ["Oxidation", "Acetyl", "Label+Heavy", "Label"];
    assert_eq!(
        components("Oxidation+Acetyl", &declared),
        vec!["Oxidation", "Acetyl"]
    );
    assert_eq!(components("Label+Heavy", &declared), vec!["Label+Heavy"]);
    assert_eq!(
        components("Label+Heavy+Oxidation", &declared),
        vec!["Label+Heavy", "Oxidation"]
    );
    assert_eq!(
        components("Found+Elsewhere", &declared),
        vec!["Found+Elsewhere"]
    );
}

#[test]
fn classify_splits_fixed_from_variable_by_declaration() {
    let fixed = vec![spec("Carbamidomethyl", Some('C'), ModLocation::Anywhere)];
    let declared = ["Carbamidomethyl", "Oxidation", "Acetyl"];
    let residues = peptide("ACMK", &[(1, "Carbamidomethyl+Acetyl"), (2, "Oxidation")]);
    let m = classify_mods(&residues, false, false, &fixed, &declared);
    assert_eq!(m.variable, vec!["Acetyl", "Oxidation"]);
    assert_eq!(m.fixed_carried, vec![0]);
    assert_eq!(m.fixed_eligible, vec![0]);
    let none = classify_mods(&peptide("AMK", &[]), false, false, &fixed, &declared);
    assert!(none.fixed_eligible.is_empty() && none.fixed_carried.is_empty());
}

#[test]
fn a_name_with_a_plus_is_one_modification_and_one_declaration() {
    let fixed = vec![spec("Label+Heavy", Some('K'), ModLocation::Anywhere)];
    let declared = ["Label+Heavy", "Label", "Heavy"];
    let m = classify_mods(
        &peptide("AKR", &[(1, "Label+Heavy")]),
        false,
        false,
        &fixed,
        &declared,
    );
    assert_eq!(m.fixed_carried, vec![0]);
    assert!(m.variable.is_empty());
}

#[test]
fn terminal_fixed_declarations_are_separate_rows_with_their_own_denominator() {
    let fixed = vec![
        spec("TMT", Some('K'), ModLocation::Anywhere),
        spec("TMT", None, ModLocation::NTerm),
    ];
    let declared = ["TMT"];
    let both = classify_mods(
        &peptide("AKR", &[(0, "TMT"), (1, "TMT")]),
        false,
        false,
        &fixed,
        &declared,
    );
    assert_eq!(both.fixed_carried, vec![0, 1]);
    assert_eq!(both.fixed_eligible, vec![0, 1]);
    let k_first = classify_mods(
        &peptide("KA", &[(0, "TMT+TMT")]),
        false,
        false,
        &fixed,
        &declared,
    );
    assert_eq!(k_first.fixed_carried, vec![0, 1]);
    let prot = classify_mods(
        &peptide("KA", &[(0, "TMT")]),
        true,
        false,
        &fixed,
        &declared,
    );
    assert_eq!(prot.fixed_carried, vec![0]);
    assert_eq!(prot.fixed_eligible, vec![0]);
}

#[test]
fn fixed_rows_count_declarations_against_their_eligible_psms_and_stay_within_100() {
    let fixed = vec![
        spec("Carbamidomethyl", Some('C'), ModLocation::Anywhere),
        spec("TMT", None, ModLocation::NTerm),
    ];
    let mut a = fact(0, 2, "ACK", 0.0, 0, &["Oxidation"]);
    a.fixed_carried = vec![0, 1];
    a.fixed_eligible = vec![0, 1];
    let mut b = fact(0, 2, "AAK", 0.0, 0, &[]);
    b.fixed_carried = vec![1];
    b.fixed_eligible = vec![1];
    let mut c = fact(0, 2, "CCK", 0.0, 0, &[]);
    c.fixed_carried = vec![0, 1];
    c.fixed_eligible = vec![0, 1];
    let ctx = SummaryContext {
        fixed: &fixed,
        ..ctx()
    };
    let s =
        IdentificationSummary::compute(&tally(vec![10], 10, 0), &names(1), &[a, b, c], 0, &ctx);
    assert_eq!(
        s.fixed_mods,
        vec![
            FixedModRow {
                name: "Carbamidomethyl".to_string(),
                target: "C".to_string(),
                site: "a C".to_string(),
                psms: 2,
                eligible: 2,
            },
            FixedModRow {
                name: "TMT".to_string(),
                target: "peptide N-terminus".to_string(),
                site: "the peptide N-terminus".to_string(),
                psms: 3,
                eligible: 3,
            },
        ]
    );
    assert!(s.fixed_mods.iter().all(|f| f.psms <= f.eligible));
    let out = s.render();
    assert!(out.contains("  Fixed: Carbamidomethyl (C) on 2 PSMs (100.0 % of PSMs with a C)\n"));
    assert!(out.contains(
        "  Fixed: TMT (peptide N-terminus) on 3 PSMs (100.0 % of PSMs with the peptide N-terminus)\n"
    ));
    assert_eq!(s.ptm_counts, vec![("Oxidation".to_string(), 1)]);
    assert_eq!(s.unmodified_psms, 2);
}

#[test]
fn calibration_shift_adds_the_error_after_the_shift() {
    let mut a = fact(0, 2, "AAK", 2.0, 0, &[]);
    a.ppm_calibrated = -1.0;
    let mut b = fact(0, 2, "BBK", 4.0, 0, &[]);
    b.ppm_calibrated = 1.0;
    let with_shift = SummaryContext {
        calibration_shift_ppm: 3.0,
        ..ctx()
    };
    let s = IdentificationSummary::compute(
        &tally(vec![10], 10, 0),
        &names(1),
        &[a.clone(), b.clone()],
        0,
        &with_shift,
    );
    assert_eq!(s.mass_error_ppm.unwrap().0, 3.0);
    assert_eq!(s.mass_error_after_shift_ppm.unwrap().0, 0.0);
    let out = s.render();
    assert!(out.contains("Precursor mass error (raw) : median +3.0 ppm"));
    assert!(out.contains("after calibration shift"));
    assert!(out.contains("median +0.0 ppm (IQR −0.5 … +0.5)  [shift +3.0 ppm]"));
    let none = summarise(&tally(vec![10], 10, 0), 1, &[a, b], 0);
    assert!(none.mass_error_after_shift_ppm.is_none());
    assert!(!none.render().contains("after calibration shift"));
}

#[test]
fn signed_ppm_uses_a_true_minus_and_no_negative_zero() {
    assert_eq!(signed(-1.14), "−1.1");
    assert_eq!(signed(0.4), "+0.4");
    assert_eq!(signed(-0.04), "+0.0");
}
