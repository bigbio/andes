//! Tests of per-run quantification.

use super::artifacts::quant_weight;
use super::collect::{GateCounts, QuantCollector, RunCapture, RunScans};
use super::hits::{select_hits, select_rank1_all, Hit};
use super::scans::{
    check_unique_stems, covers_reporters, next_survey, peptide_ion_current, survey_scan,
};
use super::settings::QuantSettings;
use crate::cli::Cli;
use model::scan::{Ms1Scan, ScanContext};
use model::Tolerance;
use quant::isobaric::Plex;
use quant::lfq::LfqParams;
use quant::ms1::Ms1RunIndex;
use search::psm::PsmMatch;

use super::*;

fn faims_run() -> Ms1RunIndex {
    Ms1RunIndex::new(
        (0..6)
            .map(|i| Ms1Scan {
                rt: i as f64,
                peaks: Vec::new(),
                id: format!("scan={}", i + 1),
                faims_cv: Some(if i % 2 == 0 { -40.0 } else { -60.0 }),
                mobility: Vec::new(),
            })
            .collect(),
    )
}

fn context(parent: Option<&str>, cv: Option<f32>) -> ScanContext {
    ScanContext {
        id: "scan=99".into(),
        parent_id: parent.map(str::to_string),
        scan_window: None,
        faims_cv: cv,
        noise: None,
        ion_mobility: None,
    }
}

#[test]
fn survey_scan_follows_the_precursor_reference_then_the_voltage() {
    let ms1 = faims_run();
    let c = context(Some("scan=3"), Some(-40.0));
    assert_eq!(survey_scan(&ms1, Some(&c), 3.5), Some(2));
    let c = context(None, Some(-40.0));
    assert_eq!(survey_scan(&ms1, Some(&c), 3.5), Some(2));
    let c = context(Some("scan=1000"), Some(-60.0));
    assert_eq!(survey_scan(&ms1, Some(&c), 3.5), Some(3));
    assert_eq!(survey_scan(&ms1, None, 4.5), Some(4));
    assert_eq!(survey_scan(&ms1, None, -1.0), None);
}

fn candidate(seq: &[u8], is_decoy: bool) -> Candidate {
    let residues = seq
        .iter()
        .map(|&r| model::AminoAcid::standard(r).unwrap())
        .collect();
    Candidate {
        peptide: model::Peptide::new(residues, b'K', b'S'),
        protein_index: 0,
        start_offset_in_protein: 0,
        is_decoy,
        is_protein_n_term: false,
        is_protein_c_term: false,
    }
}

fn psm(rank_score: f32, candidate: u32, co_isolated_mz: Option<f64>) -> PsmMatch {
    PsmMatch {
        spectrum_idx: 0,
        candidate_idxs: vec![candidate].into(),
        charge_used: 2,
        mass_error_ppm: 0.0,
        score: rank_score,
        rank_score,
        edge_score: 0,
        activation_method: None,
        features: search::psm::PsmFeatures::default(),
        isotope_offset: 0,
        precursor_mz_override: co_isolated_mz,
        site_evidence: None,
        tie_key: 0,
    }
}

fn one_spectrum(psms: Vec<PsmMatch>) -> (QuantCollector, Vec<Spectrum>, Vec<TopNQueue>) {
    let mut c = QuantCollector::new(settings(true, Some("q")));
    c.runs.push(RunCapture {
        stem: "run".into(),
        title_prefix: None,
        span: 0..1,
        ms1: None,
        ms3_reporters: HashMap::new(),
        contexts: HashMap::new(),
    });
    let spec = Spectrum {
        title: "Scan 1".into(),
        precursor_mz: 500.0,
        precursor_intensity: None,
        precursor_charge: Some(2),
        rt_seconds: Some(60.0),
        scan: Some(1),
        peaks: vec![],
        activation_method: None,
        isolation_lower_offset: None,
        isolation_upper_offset: None,
        isolation_target_mz: None,
        mobility_window: None,
    };
    let mut queue = TopNQueue::new(10);
    for p in psms {
        queue.push(p);
    }
    (c, vec![spec], vec![queue])
}

fn selected(hits: &[Hit<'_>]) -> Vec<(String, bool, bool)> {
    hits.iter()
        .map(|h| (h.spec_id.clone(), h.cand.is_decoy, h.co_isolated))
        .collect()
}

#[test]
fn plain_search_quantifies_rank1_rows_only() {
    let cands = vec![
        candidate(b"PEPTIDEK", false),
        candidate(b"SAMPLERK", false),
        candidate(b"KEDITPEP", true),
    ];
    let (c, spectra, queues) = one_spectrum(vec![
        psm(10.0, 0, None),
        psm(9.0, 1, None),
        psm(8.0, 2, None),
    ]);
    let hits = select_rank1_all(&c, &spectra, &queues, &cands);
    assert_eq!(
        selected(&hits),
        vec![("Scan 1_1_1_0".to_string(), false, false)]
    );
}

#[test]
fn chimeric_search_quantifies_each_co_isolated_precursor() {
    let cands = vec![
        candidate(b"PEPTIDEK", false),
        candidate(b"SAMPLERK", false),
        candidate(b"LLLLLLLK", false),
        candidate(b"KEDITPEP", true),
    ];
    let (c, spectra, queues) = one_spectrum(vec![
        psm(12.0, 1, Some(600.3)),
        psm(10.0, 0, None),
        psm(9.0, 2, None),
        psm(7.0, 3, Some(450.2)),
    ]);
    let hits = select_rank1_all(&c, &spectra, &queues, &cands);
    assert_eq!(
        selected(&hits),
        vec![
            ("Scan 1_1_1_0".to_string(), false, true),
            ("Scan 1_1_2_1".to_string(), false, false),
            ("Scan 1_1_4_3".to_string(), true, true),
        ]
    );
    let hits = select_hits(&c, &spectra, &queues, &cands, None, &HashSet::new(), None);
    assert_eq!(hits.len(), 2);
    let mut rescore = HashMap::new();
    for (id, q) in [("Scan 1_1_1_0", 0.2), ("Scan 1_1_2_1", 0.001)] {
        rescore.insert(
            id.to_string(),
            PercolatorPsm {
                pooled_q_value: None,
                psm_id: id.to_string(),
                q_value: q,
                pep: 0.0,
                score: Some(1.0),
                peptide: String::new(),
                proteins: String::new(),
            },
        );
    }
    let hits = select_hits(
        &c,
        &spectra,
        &queues,
        &cands,
        Some(&rescore),
        &HashSet::new(),
        None,
    );
    assert_eq!(
        selected(&hits),
        vec![("Scan 1_1_2_1".to_string(), false, false)]
    );
    assert_eq!(hits[0].psm_index, 1);
}

#[test]
fn a_row_psms_parquet_leaves_out_takes_no_psm_index() {
    let cands = vec![candidate(b"KEDITPEP", true), candidate(b"PEPTIDEK", false)];
    let (c, spectra, queues) = one_spectrum(vec![psm(12.0, 0, Some(600.3)), psm(10.0, 1, None)]);
    let mut rescore = HashMap::new();
    rescore.insert(
        "Scan 1_1_2_1".to_string(),
        PercolatorPsm {
            pooled_q_value: None,
            psm_id: "Scan 1_1_2_1".to_string(),
            q_value: 0.001,
            pep: 0.0,
            score: Some(1.0),
            peptide: String::new(),
            proteins: String::new(),
        },
    );
    let omitted: HashSet<String> = ["Scan 1_1_1_0".to_string()].into();
    let hits = select_hits(
        &c,
        &spectra,
        &queues,
        &cands,
        Some(&rescore),
        &omitted,
        None,
    );
    assert_eq!(
        selected(&hits),
        vec![("Scan 1_1_2_1".to_string(), false, false)]
    );
    assert_eq!(hits[0].psm_index, 0);
    let hits = select_hits(
        &c,
        &spectra,
        &queues,
        &cands,
        Some(&rescore),
        &HashSet::new(),
        None,
    );
    assert_eq!(hits[0].psm_index, 1);
}

#[test]
fn omitted_rows_are_skipped_without_rescoring() {
    let cands = vec![candidate(b"PEPTIDEK", false), candidate(b"KEDITPEP", false)];
    let (c, spectra, queues) = one_spectrum(vec![psm(12.0, 0, None)]);
    let hits = select_hits(&c, &spectra, &queues, &cands, None, &HashSet::new(), None);
    assert_eq!(hits.len(), 1);
    let omitted: HashSet<String> = [hits[0].spec_id.clone()].into();
    assert!(select_hits(&c, &spectra, &queues, &cands, None, &omitted, None).is_empty());
}

fn settings(lfq: bool, quant_dir: Option<&str>) -> QuantSettings {
    QuantSettings {
        plex: None,
        tmt_level: 2,
        tmt_tol: Tolerance::Ppm(20.0),
        correction: None,
        min_purity: 0.0,
        min_sn: 0.0,
        lfq,
        lfq_params: LfqParams::default(),
        feature_fdr: 0.05,
        psm_fdr: 0.01,
        pep_cap: None,
        quant_dir: quant_dir.map(PathBuf::from),
        fragment_tol: Some(Tolerance::Ppm(20.0)),
        unmix: None,
        tmt_mbr: false,
        scratch_dir: PathBuf::from("."),
        label: None,
    }
}

#[test]
fn quant_dir_needs_ms1_scans_from_every_input() {
    let mut c = QuantCollector::new(settings(true, Some("q")));
    let e = c
        .finish_file(Path::new("/in/run_a.mzML"), None, RunScans::default(), 0)
        .unwrap_err();
    assert!(e.contains("/in/run_a.mzML") && e.contains("MS1"), "{e}");
    let mut c = QuantCollector::new(settings(true, None));
    c.finish_file(Path::new("/in/run_a.mzML"), None, RunScans::default(), 0)
        .unwrap();
}

#[test]
fn duplicate_run_stems_name_both_inputs() {
    let ok = [
        PathBuf::from("a/run1.mzML"),
        PathBuf::from("a/run2.mzML.gz"),
    ];
    check_unique_stems(&ok).unwrap();
    let dup = [
        PathBuf::from("a/run1.mzML"),
        PathBuf::from("b/run2.raw"),
        PathBuf::from("b/run1.mzML.gz"),
    ];
    let e = check_unique_stems(&dup).unwrap_err();
    assert!(
        e.contains("a/run1.mzML ") && e.contains("b/run1.mzML.gz") && e.contains("\"run1\""),
        "{e}"
    );
}

#[test]
fn quantification_gate_counts_the_first_failed_criterion() {
    let mut s = settings(false, None);
    s.min_purity = 0.7;
    s.min_sn = 100.0;
    let mut g = GateCounts::default();
    assert!(g.check(&s, Some(0.9), Some(500.0)));
    assert!(g.check(&s, None, None));
    assert!(!g.check(&s, Some(0.5), Some(50.0)));
    assert!(!g.check(&s, Some(0.8), Some(50.0)));
    assert_eq!((g.purity, g.sn), (1, 1));
    let text = g.describe(&s);
    assert!(text.contains("1 below purity 0.7"), "{text}");
    assert!(text.contains("1 below summed S/N 100"), "{text}");
    let off = settings(false, None);
    assert!(GateCounts::default().check(&off, Some(0.0), Some(0.0)));
    assert!(GateCounts::default().describe(&off).is_empty());
}

#[test]
fn roll_up_weight_is_purity_squared_times_the_chance_the_identification_is_right() {
    assert_eq!(quant_weight(None, None), 1.0);
    assert!((quant_weight(Some(0.8), None) - 0.64).abs() < 1e-6);
    assert_eq!(quant_weight(Some(1.2), None), 1.0);
    assert!((quant_weight(Some(0.8), Some(0.25)) - 0.48).abs() < 1e-6);
    assert_eq!(quant_weight(None, Some(1.0)), 0.0);
}

#[test]
fn reporter_coverage_needs_the_whole_channel_range() {
    let tol = Tolerance::Ppm(20.0);
    assert!(covers_reporters(None, Plex::Tmt16, tol));
    assert!(covers_reporters(Some((100.0, 2000.0)), Plex::Tmt16, tol));
    assert!(covers_reporters(Some((126.0, 2000.0)), Plex::Tmt16, tol));
    assert!(!covers_reporters(Some((130.0, 2000.0)), Plex::Tmt16, tol));
    assert!(!covers_reporters(Some((148.0, 1866.0)), Plex::Tmt10, tol));
    assert!(!covers_reporters(Some((100.0, 134.0)), Plex::Tmt16, tol));
    let tol = Tolerance::Ppm(10.0);
    assert!(covers_reporters(Some((126.0, 135.1625)), Plex::Tmt35, tol));
    assert!(!covers_reporters(Some((126.0, 135.1590)), Plex::Tmt35, tol));
    assert!(!covers_reporters(Some((126.2, 2000.0)), Plex::Tmt32, tol));
}

fn search_args(extra: &[&str]) -> Cli {
    use clap::{CommandFactory, FromArgMatches};
    let mut argv = vec![
        "andes",
        "--spectrum",
        "x",
        "--database",
        "y",
        "--output-pin",
        "o",
    ];
    argv.extend_from_slice(extra);
    let matches = crate::TopCli::command()
        .try_get_matches_from(argv)
        .expect("parse args");
    crate::SearchArgs::from_arg_matches(&matches).unwrap()
}

#[test]
fn tmt_mbr_is_refused_for_ms3_reporters_and_low_resolution_ms2() {
    let from = |extra: &[&str], high_res: bool, level: u8| {
        let mut a = vec!["--tmt", "tmt11", "--mbr", "--rescore"];
        a.extend_from_slice(extra);
        QuantSettings::from_cli(&search_args(&a), high_res, level)
    };
    let s = from(&[], true, 2).unwrap().unwrap();
    assert!(s.tmt_mbr);
    let e = from(&[], true, 3).unwrap_err();
    assert!(e.contains("SPS-MS3"), "{e}");
    let e = from(&["--tmt-level", "3"], false, 3).unwrap_err();
    assert!(e.contains("SPS-MS3"), "{e}");
    let e = from(&[], false, 2).unwrap_err();
    assert!(e.contains("low-resolution"), "{e}");
    let plain = |high_res: bool, level: u8| {
        QuantSettings::from_cli(&search_args(&["--tmt", "tmt11"]), high_res, level)
            .unwrap()
            .unwrap()
    };
    assert!(!plain(true, 2).tmt_mbr);
    assert!(!plain(true, 3).tmt_mbr);
}

#[test]
fn command_line_checks_need_no_input() {
    let ok = |extra: &[&str]| QuantSettings::check_args(&search_args(extra));
    assert!(ok(&[]).is_ok());
    assert!(ok(&["--tmt", "tmt35"]).is_ok());
    assert!(ok(&["--tmt", "tmt35", "--tmt-tol", "10ppm"]).is_ok());
    assert!(ok(&["--tmt", "tmt10", "--tmt-correction", "lot.txt"]).is_ok());
    let e = ok(&["--tmt", "tmt35", "--tmt-correction", "lot.txt"]).unwrap_err();
    assert!(e.contains("not supported for tmt35plex"), "{e}");
    let e = ok(&["--tmt", "tmt32", "--tmt-tol", "20ppm"]).unwrap_err();
    assert!(e.contains("e.g. 10ppm"), "{e}");
    let e = ok(&["--lfq", "--tmt-correction", "lot.txt"]).unwrap_err();
    assert!(e.contains("needs --tmt"), "{e}");
    let e = ok(&["--tmt", "tmt10", "--tmt-level", "4"]).unwrap_err();
    assert!(e.contains("--tmt-level"), "{e}");
    assert!(ok(&["--tmt", "tmt16", "--tmt-unmix", "on"]).is_ok());
    assert!(ok(&["--tmt-unmix", "auto"]).is_ok());
    assert!(ok(&["--tmt", "tmt16", "--tmt-level", "3"]).is_ok());
    assert!(ok(&["--tmt-unmix", "off"]).is_ok());
    let e = ok(&["--tmt-unmix", "on"]).unwrap_err();
    assert!(e.contains("needs --tmt"), "{e}");
    let e = ok(&["--tmt", "tmt16", "--tmt-level", "3", "--tmt-unmix", "on"]).unwrap_err();
    assert!(e.contains("--tmt-level 3"), "{e}");
}

#[test]
fn peptide_ion_current_skips_the_reporter_region() {
    let peaks = [
        (126.13, 1e6f32),
        (134.15, 1e6),
        (135.1, 5.0),
        (300.0, 10.0),
        (800.0, 20.0),
    ];
    assert_eq!(peptide_ion_current(&peaks, Plex::Tmt16), 30.0);
}

#[test]
fn next_survey_stays_in_the_faims_series() {
    let idx = faims_run();
    let cv0 = idx.scan(0).faims_cv;
    assert_eq!(next_survey(&idx, 0, cv0), Some(2));
    assert_eq!(next_survey(&idx, 0, None), Some(1));
    assert_eq!(next_survey(&idx, 5, cv0), None);
}
