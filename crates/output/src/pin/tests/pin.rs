use super::*;
use model::amino_acid::AminoAcid;
use model::peptide::Peptide;
use model::protein::{Protein, ProteinDb};
use model::tolerance::PrecursorTolerance;
use model::tolerance::Tolerance;
use search::candidate_gen::Candidate;
use search::search_index::SearchIndex;

#[test]
fn write_mass_keeps_micro_dalton_precision() {
    let fmt = |v: f64| {
        let mut out = Vec::new();
        write_mass(&mut out, v).unwrap();
        String::from_utf8(out).unwrap()
    };
    assert_eq!(fmt(1674.791234), "1674.791234");
    assert_eq!(fmt(1674.801234), "1674.801234");
    assert_eq!(fmt(790.3345), "790.3345");
    assert_eq!(fmt(1500.0), "1500");
    assert_eq!(fmt(0.0), "0");
    assert_eq!(fmt(f64::NAN), "0");
}

fn make_search_index(accession: &str) -> SearchIndex {
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: accession.to_string(),
            description: String::new(),
            sequence: b"MKWVTFISLL".to_vec(),
            annotation: None,
        }],
    };
    SearchIndex::from_target_db(&target, "XXX_")
}

fn make_empty_search_index() -> SearchIndex {
    let target = ProteinDb { proteins: vec![] };
    SearchIndex::from_target_db(&target, "XXX_")
}

fn make_spectrum(title: &str, scan: i32, precursor_mz: f64) -> Spectrum {
    Spectrum {
        title: title.to_string(),
        precursor_mz,
        precursor_intensity: None,
        precursor_charge: Some(2),
        rt_seconds: None,
        scan: Some(scan),
        peaks: vec![],
        activation_method: None,
        isolation_lower_offset: None,
        isolation_upper_offset: None,
        isolation_target_mz: None,
        mobility_window: None,
    }
}

fn make_candidate(protein_index: usize, is_decoy: bool) -> Candidate {
    let aa = AminoAcid::standard(b'A').unwrap();
    let peptide = Peptide::new(vec![aa], b'K', b'S');
    Candidate {
        peptide,
        protein_index,
        start_offset_in_protein: 0,
        is_decoy,
        is_protein_n_term: false,
        is_protein_c_term: false,
    }
}

fn make_candidate_seq(seq: &[u8], mods: &[(usize, f64)], protein_index: usize) -> Candidate {
    use model::modification::{ModLocation, Modification, ResidueSpec};
    let residues = seq
        .iter()
        .enumerate()
        .map(|(i, &r)| {
            let mut aa = AminoAcid::standard(r).unwrap();
            if let Some(&(_, delta)) = mods.iter().find(|&&(p, _)| p == i) {
                let m: &'static Modification = Box::leak(Box::new(Modification {
                    name: format!("mod{delta}"),
                    mass_delta: delta,
                    residue: ResidueSpec::Specific(r),
                    location: ModLocation::Anywhere,
                    fixed: false,
                    accession: None,
                    neutral_losses: Vec::new(),
                    loss_class: 0,
                    annotated_site_mass: None,
                    blocks_cleavage: false,
                    parts: Vec::new(),
                    diagnostic_ions: Vec::new(),
                }));
                aa.mod_ = Some(m);
            }
            aa
        })
        .collect();
    Candidate {
        peptide: Peptide::new(residues, b'K', b'S'),
        protein_index,
        start_offset_in_protein: 0,
        is_decoy: false,
        is_protein_n_term: false,
        is_protein_c_term: false,
    }
}

fn emitted(psms: Vec<PsmMatch>, cands: &[Candidate]) -> Vec<(String, String)> {
    let params = make_params(2..=3);
    let spectra = vec![make_spectrum("Scan 1", 1, 500.0)];
    let mut queue = TopNQueue::new(10);
    for p in psms {
        queue.push(p);
    }
    let mut buf = Vec::<u8>::new();
    write_pin_to(
        &mut buf,
        &spectra,
        &[queue],
        cands,
        &params,
        &make_empty_search_index(),
    )
    .unwrap();
    let header = parse_header(&buf);
    let pep_col = header.iter().position(|h| h == "Peptide").unwrap();
    parse_rows(&buf)
        .into_iter()
        .map(|r| (r[0].clone(), r[pep_col].clone()))
        .collect()
}

fn make_psm(
    spectrum_idx: usize,
    score: f32,
    rank_score: f32,
    candidate_idx: u32,
    charge: u8,
) -> PsmMatch {
    PsmMatch {
        spectrum_idx,
        candidate_idxs: vec![candidate_idx].into(),
        charge_used: charge,
        mass_error_ppm: 1.5,
        score,
        rank_score,
        edge_score: 0,
        activation_method: Some(model::activation::ActivationMethod::HCD),
        features: search::psm::PsmFeatures::default(),
        isotope_offset: 0,
        precursor_mz_override: None,
        site_evidence: None,
        tie_key: 0,
    }
}

fn make_params(charge_range: std::ops::RangeInclusive<u8>) -> SearchParams {
    use model::aa_set::AminoAcidSetBuilder;
    let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
    SearchParams {
        aa_set,
        enzyme: model::enzyme::Enzyme::Trypsin,
        extra_enzymes: Vec::new(),
        min_length: 6,
        max_length: 40,
        max_missed_cleavages: 1,
        max_variable_mods_per_peptide: 3,
        max_label_sites: 0,
        precursor_tolerance: PrecursorTolerance::symmetric(Tolerance::Ppm(20.0)),
        charge_range,
        isotope_error_range: -1..=2,
        top_n_psms_per_spectrum: 10,
        num_tolerable_termini: 2,
        min_peaks: 10,
        precursor_cal_mode: search::PrecursorCalMode::Auto,
        cal_min_spec_keys: search::precursor_cal::constants::MIN_SPECKEYS_FOR_PREPASS,
        fragment_index_top_k: 0,
        precursor_mass_shift_ppm: 0.0,
        chimeric: false,
        dia: false,
        dia_ms1_tol_ppm: 5.0,
        dia_vote_filter: None,
        chimeric_isolation_halfwidth_da: 1.5,
        chimeric_max_coisolated: 2,
        chimeric_max_kl: 0.3,
        score_mode: search::ScoreMode::Rank,
        refine_select_psm_fdr: 0.01,
        require_nxst_sequon: false,
        open_search: None,
    }
}

fn parse_header(output: &[u8]) -> Vec<String> {
    let text = std::str::from_utf8(output).unwrap();
    let first_line = text.lines().next().unwrap_or("");
    first_line.split('\t').map(|s| s.to_string()).collect()
}

fn parse_rows(output: &[u8]) -> Vec<Vec<String>> {
    let text = std::str::from_utf8(output).unwrap();
    text.lines()
        .skip(1)
        .filter(|l| !l.is_empty())
        .map(|l| l.split('\t').map(|s| s.to_string()).collect())
        .collect()
}

#[test]
fn pin_header_columns_are_gf_free_schema() {
    let expected: Vec<&str> = vec![
        "SpecId",
        "Label",
        "ScanNr",
        "ExpMass",
        "CalcMass",
        "FileName",
        "mass",
        "RankScore",
        "isotope_error",
        "peplen",
        "dm",
        "absdm",
        "charge2",
        "charge3",
        "enzN",
        "enzC",
        "enzInt",
        "NumMatchedMainIons",
        "longest_b",
        "longest_y",
        "longest_y_pct",
        "ExplainedIonCurrentRatio",
        "NTermIonCurrentRatio",
        "CTermIonCurrentRatio",
        "MS2IonCurrent",
        "IsolationWindowEfficiency",
        "MeanErrorTop7",
        "StdevErrorTop7",
        "MeanRelErrorTop7",
        "StdevRelErrorTop7",
        "matchedIonRatio",
        "EdgeScore",
        "PrecursorIsotopeKL",
        "PrecursorSNR",
        "DeltaRankScore",
        "TailorScore",
        "PpmGaussianScore",
        "NeutralLossIonCount",
        "LongestComplementaryLadder",
        "ComplementaryIonBalance",
        "MeanMatchedIntensityRank",
        "DoublyChargedMatchedIonCount",
        "UniqueMatchFraction",
        "ChanceMatchSurprise",
        "IntensitySignal",
        "FragPredExplained",
        "FragPredChanceLLR",
        "FragTopKObserved",
        "RichIonLLR",
        "IsRefinement",
        "NumMods",
        "RefinementModClass",
        "ModSiteShiftedMatched",
        "ModSiteShiftedFrac",
        "ModSiteIntensFrac",
        "ModSiteLocalized",
        "ModSiteDetCount",
        "MassCompetitionEvidence",
        "CandidateRankEntropy",
        "ListwiseScoreGap",
        "RawScore",
        "RawScoreCal",
        "RankScoreFloat",
        "DeltaRT",
        "AbsDeltaRT",
        "DeltaRTNorm",
        "IsVariant",
        "IsIsoform",
        "IsChain",
        "RefCompetitor",
        "DeltaRefScore",
        "IsotopeMimic",
        "IsotopeErrorAbs",
        "IsotopeErrorNeg",
        "Peptide",
        "Proteins",
    ];

    let params = make_params(2..=3);
    let spectra: Vec<Spectrum> = vec![];
    let queues: Vec<TopNQueue> = vec![];
    let idx = make_empty_search_index();

    let mut buf = Vec::<u8>::new();
    let cands: Vec<Candidate> = vec![];
    write_pin_to(&mut buf, &spectra, &queues, &cands, &params, &idx).unwrap();

    let cols = parse_header(&buf);
    assert_eq!(
        cols, expected,
        "PIN header columns must match the reference fixture column order exactly"
    );
}

#[test]
fn pin_header_has_refinement_columns() {
    let params = make_params(2..=3);
    let spectra: Vec<Spectrum> = vec![];
    let queues: Vec<TopNQueue> = vec![];
    let idx = make_empty_search_index();

    let mut buf = Vec::<u8>::new();
    let cands: Vec<Candidate> = vec![];
    write_pin_to(&mut buf, &spectra, &queues, &cands, &params, &idx).unwrap();

    let cols = parse_header(&buf);
    for col in ["IsRefinement", "NumMods", "RefinementModClass"] {
        assert!(cols.iter().any(|c| c == col), "header missing column {col}");
    }
}

#[test]
fn pin_writes_label_minus_one_for_decoy() {
    let params = make_params(2..=3);
    let spectra = vec![make_spectrum("Scan 1", 1, 500.0)];

    let mut queue = TopNQueue::new(10);
    queue.push(make_psm(0, 10.0, 10.0, 0, 2));
    let queues = vec![queue];
    let idx = make_empty_search_index();

    let mut buf = Vec::<u8>::new();
    let cands = vec![make_candidate(0, true)];
    write_pin_to(&mut buf, &spectra, &queues, &cands, &params, &idx).unwrap();

    let rows = parse_rows(&buf);
    assert_eq!(rows.len(), 1, "should have 1 data row");

    assert_eq!(rows[0][1], "-1", "decoy PSM should have Label = -1");
}

#[test]
fn pin_writes_the_isotope_offset_class_columns() {
    let params = make_params(2..=3);
    let spectra = vec![
        make_spectrum("a", 1, 500.0),
        make_spectrum("b", 2, 500.0),
        make_spectrum("c", 3, 500.0),
        make_spectrum("d", 4, 500.0),
    ];
    let queues: Vec<TopNQueue> = [-1i8, 0, 1, 2]
        .iter()
        .enumerate()
        .map(|(i, &offset)| {
            let mut psm = make_psm(i, 10.0, 10.0, 0, 2);
            psm.isotope_offset = offset;
            let mut q = TopNQueue::new(10);
            q.push(psm);
            q
        })
        .collect();
    let mut buf = Vec::<u8>::new();
    let cands = vec![make_candidate(0, false)];
    write_pin_to(
        &mut buf,
        &spectra,
        &queues,
        &cands,
        &params,
        &make_empty_search_index(),
    )
    .unwrap();
    let header = parse_header(&buf);
    let col = |name: &str| header.iter().position(|h| h == name).unwrap();
    let (iso, abs, neg) = (
        col("isotope_error"),
        col("IsotopeErrorAbs"),
        col("IsotopeErrorNeg"),
    );
    let got: Vec<(String, String, String)> = parse_rows(&buf)
        .into_iter()
        .map(|r| (r[iso].clone(), r[abs].clone(), r[neg].clone()))
        .collect();
    let want = [
        ("-1", "1", "1"),
        ("0", "0", "0"),
        ("1", "1", "0"),
        ("2", "2", "0"),
    ];
    assert_eq!(got.len(), want.len());
    for (g, w) in got.iter().zip(want) {
        assert_eq!((g.0.as_str(), g.1.as_str(), g.2.as_str()), w);
    }
}

#[test]
fn write_pin_errors_on_duplicate_spec_id_base() {
    let params = make_params(2..=3);
    let spectra = vec![
        make_spectrum("dup", 7, 500.0),
        make_spectrum("dup", 7, 600.0),
    ];
    let mut q0 = TopNQueue::new(10);
    q0.push(make_psm(0, 10.0, 10.0, 0, 2));
    let mut q1 = TopNQueue::new(10);
    q1.push(make_psm(1, 9.0, 9.0, 0, 2));
    let queues = vec![q0, q1];
    let idx = make_empty_search_index();
    let cands = vec![make_candidate(0, false)];

    let mut buf = Vec::<u8>::new();
    let err = write_pin_to(&mut buf, &spectra, &queues, &cands, &params, &idx)
        .expect_err("duplicate SpecId base must error");
    assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
    assert!(
        err.to_string().contains("duplicate PIN SpecId"),
        "got: {err}"
    );
}

#[test]
fn write_pin_is_atomic() {
    let params = make_params(2..=3);
    let idx = make_empty_search_index();
    let cands = vec![make_candidate(0, false)];
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("run.pin");
    std::fs::write(&path, "previous").unwrap();

    let dup = vec![
        make_spectrum("dup", 7, 500.0),
        make_spectrum("dup", 7, 600.0),
    ];
    let mut q0 = TopNQueue::new(10);
    q0.push(make_psm(0, 10.0, 10.0, 0, 2));
    let mut q1 = TopNQueue::new(10);
    q1.push(make_psm(1, 9.0, 9.0, 0, 2));
    assert!(write_pin(&path, &dup, &[q0, q1], &cands, &params, &idx).is_err());
    assert_eq!(std::fs::read_to_string(&path).unwrap(), "previous");
    assert!(!crate::atomic::partial(&path).exists());

    let one = vec![make_spectrum("Scan 1", 1, 500.0)];
    let mut q = TopNQueue::new(10);
    q.push(make_psm(0, 10.0, 10.0, 0, 2));
    write_pin(&path, &one, &[q], &cands, &params, &idx).unwrap();
    assert!(std::fs::read_to_string(&path)
        .unwrap()
        .starts_with("SpecId\t"));
    assert!(!crate::atomic::partial(&path).exists());
}

#[test]
fn pin_writes_charge_one_hot_correctly() {
    let params = make_params(2..=3);
    let spectra = vec![make_spectrum("Scan 1", 1, 500.0)];

    let mut queue = TopNQueue::new(10);
    queue.push(make_psm(0, 10.0, 10.0, 0, 2));
    let queues = vec![queue];
    let idx = make_empty_search_index();

    let mut buf = Vec::<u8>::new();
    let cands = vec![make_candidate(0, false)];
    write_pin_to(&mut buf, &spectra, &queues, &cands, &params, &idx).unwrap();

    let cols = parse_header(&buf);
    let rows = parse_rows(&buf);
    assert_eq!(rows.len(), 1);

    let charge2_idx = cols
        .iter()
        .position(|c| c == "charge2")
        .expect("charge2 column missing");
    let charge3_idx = cols
        .iter()
        .position(|c| c == "charge3")
        .expect("charge3 column missing");

    assert_eq!(
        rows[0][charge2_idx], "1",
        "charge2 should be 1 for a charge-2 PSM"
    );
    assert_eq!(
        rows[0][charge3_idx], "0",
        "charge3 should be 0 for a charge-2 PSM"
    );
}

#[test]
fn pin_enzyme_features_use_every_enzyme_of_a_multi_protease_digest() {
    let cands = vec![make_candidate_seq(b"AEAKAE", &[], 0)];
    let spectra = vec![make_spectrum("Scan 1", 1, 500.0)];
    let idx = make_empty_search_index();
    let enz = |extra: Vec<model::enzyme::Enzyme>| {
        let mut params = make_params(2..=3);
        params.extra_enzymes = extra;
        let mut queue = TopNQueue::new(10);
        queue.push(make_psm(0, 10.0, 10.0, 0, 2));
        let mut buf = Vec::<u8>::new();
        write_pin_to(&mut buf, &spectra, &[queue], &cands, &params, &idx).unwrap();
        let cols = parse_header(&buf);
        let row = parse_rows(&buf).remove(0);
        ["enzN", "enzC", "enzInt"].map(|name| {
            let i = cols.iter().position(|c| c == name).unwrap();
            row[i].clone()
        })
    };
    assert_eq!(enz(Vec::new()), ["1", "0", "1"]);
    assert_eq!(enz(vec![model::enzyme::Enzyme::GluC]), ["1", "1", "2"]);
}

#[test]
fn pin_enz_n_is_one_for_met_clipped_targets_and_decoys() {
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".to_string(),
            description: String::new(),
            sequence: b"MASLVTEGLKPLLSEAVQKWEDNGKSSLLR".to_vec(),
            annotation: None,
        }],
    };
    let idx = SearchIndex::from_target_db(&target, "XXX_");
    let params = make_params(2..=3);
    let cands: Vec<Candidate> =
        search::candidate_gen::enumerate_candidates(&idx, &params, &idx.decoy_prefix).collect();
    let met_clipped = |decoy: bool| {
        cands
            .iter()
            .filter(|c| {
                c.is_decoy == decoy && c.is_protein_n_term && c.start_offset_in_protein == 1
            })
            .count()
    };
    assert!(met_clipped(false) > 0 && met_clipped(true) > 0);

    let spectra: Vec<Spectrum> = (0..cands.len())
        .map(|i| make_spectrum(&format!("Scan {i}"), i as i32, 500.0))
        .collect();
    let queues: Vec<TopNQueue> = (0..cands.len())
        .map(|i| {
            let mut q = TopNQueue::new(10);
            q.push(make_psm(i, 10.0, 10.0, i as u32, 2));
            q
        })
        .collect();
    let mut buf = Vec::<u8>::new();
    write_pin_to(&mut buf, &spectra, &queues, &cands, &params, &idx).unwrap();
    let cols = parse_header(&buf);
    let enz_n = cols.iter().position(|c| c == "enzN").unwrap();
    let peptide = cols.iter().position(|c| c == "Peptide").unwrap();
    let rows = parse_rows(&buf);
    assert_eq!(rows.len(), cands.len());
    for (row, cand) in rows.iter().zip(&cands) {
        assert_eq!(row[peptide], cand.peptide.to_string());
        let want = if cand.is_protein_n_term {
            "1"
        } else if cand.peptide.residues[0].residue == b'P' {
            "0"
        } else {
            "1"
        };
        assert_eq!(row[enz_n], want, "{}", row[peptide]);
    }
    assert!(rows.iter().any(|r| r[peptide] == "K.PLLSEAVQK.W"));
}

#[test]
fn pin_charge_columns_are_the_charge_range_for_any_psms() {
    let params = make_params(2..=3);
    let spectra = vec![
        make_spectrum("Scan 1", 1, 500.0),
        make_spectrum("Scan 2", 2, 500.0),
        make_spectrum("Scan 3", 3, 500.0),
    ];
    let mut queues = Vec::new();
    for (i, z) in [(0, 1u8), (1, 2), (2, 6)] {
        let mut q = TopNQueue::new(10);
        q.push(make_psm(i, 10.0, 10.0, 0, z));
        queues.push(q);
    }
    let idx = make_empty_search_index();
    let cands = vec![make_candidate(0, false)];
    let mut buf = Vec::<u8>::new();
    write_pin_to(&mut buf, &spectra, &queues, &cands, &params, &idx).unwrap();

    let cols = parse_header(&buf);
    let charge_cols: Vec<&String> = cols.iter().filter(|c| c.starts_with("charge")).collect();
    assert_eq!(charge_cols, ["charge2", "charge3"]);
    let rows = parse_rows(&buf);
    assert_eq!(rows.len(), 3);
    let first = cols.iter().position(|c| c == "charge2").unwrap();
    for (row, expected) in rows.iter().zip([["1", "0"], ["1", "0"], ["0", "1"]]) {
        let hot: Vec<&str> = row[first..first + 2].iter().map(String::as_str).collect();
        assert_eq!(hot, expected);
    }

    let mut buf = Vec::<u8>::new();
    write_pin_to(
        &mut buf,
        &spectra[1..2],
        &queues[1..2],
        &cands,
        &params,
        &idx,
    )
    .unwrap();
    let cols = parse_header(&buf);
    let charge_cols: Vec<&String> = cols.iter().filter(|c| c.starts_with("charge")).collect();
    assert_eq!(charge_cols, ["charge2", "charge3"]);
}

#[test]
fn chimeric_specids_unique_for_cofragmented_peptides_same_scan() {
    let mut params = make_params(2..=3);
    params.chimeric = true;
    let spectra = vec![make_spectrum("Scan 1", 1, 500.0)];

    let mut queue = TopNQueue::new(10);
    queue.push(make_psm(0, 10.0, 10.0, 0, 2));
    queue.push(make_psm(0, 9.0, 10.0, 1, 2));
    let queues = vec![queue];
    let idx = make_empty_search_index();
    let cands = vec![make_candidate(0, false), make_candidate_seq(b"G", &[], 1)];

    let mut buf = Vec::<u8>::new();
    write_pin_to(&mut buf, &spectra, &queues, &cands, &params, &idx).unwrap();

    let rows = parse_rows(&buf);
    assert_eq!(rows.len(), 2, "both co-fragmented PSMs should be emitted");
    assert_ne!(
        rows[0][0], rows[1][0],
        "chimeric SpecIds must be unique per row, got {:?} and {:?}",
        rows[0][0], rows[1][0]
    );
}

#[test]
fn non_chimeric_tied_psms_same_scan_get_distinct_specids() {
    let params = make_params(2..=3);
    let spectra = vec![make_spectrum("Scan 1", 1, 500.0)];

    let mut queue = TopNQueue::new(10);
    queue.push(make_psm(0, 10.0, 10.0, 0, 2));
    queue.push(make_psm(0, 10.0, 10.0, 1, 2));
    let queues = vec![queue];
    let idx = make_empty_search_index();
    let cands = vec![make_candidate(0, false), make_candidate_seq(b"G", &[], 1)];

    let mut buf = Vec::<u8>::new();
    write_pin_to(&mut buf, &spectra, &queues, &cands, &params, &idx).unwrap();

    let rows = parse_rows(&buf);
    assert_eq!(rows.len(), 2, "both tied PSMs should be emitted");
    assert_ne!(
        rows[0][0], rows[1][0],
        "non-chimeric tied SpecIds must be unique, got {:?} and {:?}",
        rows[0][0], rows[1][0]
    );
}

const PHOSPHO: f64 = 79.96633;
const OX: f64 = 15.99491;

#[test]
fn positional_isomers_keep_only_the_best_ranked_site() {
    let cands = vec![
        make_candidate_seq(b"SATSK", &[(0, PHOSPHO)], 0),
        make_candidate_seq(b"SATSK", &[(2, PHOSPHO)], 0),
        make_candidate_seq(b"SATSK", &[(3, PHOSPHO)], 0),
    ];
    let rows = emitted(
        vec![
            make_psm(0, 9.0, 9.0, 1, 2),
            make_psm(0, 10.0, 10.0, 2, 2),
            make_psm(0, 8.0, 8.0, 0, 2),
        ],
        &cands,
    );
    assert_eq!(rows.len(), 1, "{rows:?}");
    assert_eq!(rows[0].1, "K.SATS+79.96633K.S");
    assert_eq!(rows[0].0, "Scan 1_1_1");
}

#[test]
fn non_isomers_are_kept_and_ranks_stay_contiguous() {
    let cands = vec![
        make_candidate_seq(b"SATSK", &[(0, PHOSPHO)], 0),
        make_candidate_seq(b"SATSK", &[(3, PHOSPHO)], 0),
        make_candidate_seq(b"SATSK", &[], 0),
        make_candidate_seq(b"SATTK", &[(0, PHOSPHO)], 0),
        make_candidate_seq(b"MSATK", &[(0, OX), (1, PHOSPHO)], 0),
        make_candidate_seq(b"MSATK", &[(0, OX), (3, PHOSPHO)], 0),
        make_candidate_seq(b"MSATK", &[(1, PHOSPHO)], 0),
    ];
    let rows = emitted(
        vec![
            make_psm(0, 10.0, 10.0, 0, 2),
            make_psm(0, 9.0, 9.0, 1, 2),
            make_psm(0, 8.0, 8.0, 2, 2),
            make_psm(0, 7.0, 7.0, 3, 2),
            make_psm(0, 6.0, 6.0, 4, 2),
            make_psm(0, 5.0, 5.0, 5, 2),
            make_psm(0, 4.0, 4.0, 6, 2),
        ],
        &cands,
    );
    let peps: Vec<&str> = rows.iter().map(|r| r.1.as_str()).collect();
    assert_eq!(
        peps,
        [
            "K.S+79.96633ATSK.S",
            "K.SATSK.S",
            "K.S+79.96633ATTK.S",
            "K.M+15.99491S+79.96633ATK.S",
            "K.MS+79.96633ATK.S",
        ]
    );
    let ids: Vec<&str> = rows.iter().map(|r| r.0.as_str()).collect();
    assert_eq!(
        ids,
        [
            "Scan 1_1_1_0",
            "Scan 1_1_2_1",
            "Scan 1_1_3_2",
            "Scan 1_1_4_3",
            "Scan 1_1_5_4"
        ]
    );
}

#[test]
fn isomers_at_different_charge_or_precursor_are_kept() {
    let cands = vec![
        make_candidate_seq(b"SATSK", &[(0, PHOSPHO)], 0),
        make_candidate_seq(b"SATSK", &[(3, PHOSPHO)], 0),
    ];
    let rows = emitted(
        vec![make_psm(0, 10.0, 10.0, 0, 2), make_psm(0, 9.0, 9.0, 1, 3)],
        &cands,
    );
    assert_eq!(rows.len(), 2, "{rows:?}");
    let mut secondary = make_psm(0, 9.0, 9.0, 1, 2);
    secondary.precursor_mz_override = Some(501.0);
    let rows = emitted(vec![make_psm(0, 10.0, 10.0, 0, 2), secondary], &cands);
    assert_eq!(rows.len(), 2, "{rows:?}");
    let mut a = make_psm(0, 9.0, 9.0, 0, 2);
    a.precursor_mz_override = Some(501.0);
    let mut b = make_psm(0, 8.0, 8.0, 1, 2);
    b.precursor_mz_override = Some(501.0);
    let rows = emitted(vec![make_psm(0, 10.0, 10.0, 0, 2), a, b], &cands);
    assert_eq!(rows.len(), 2, "{rows:?}");
}

#[test]
fn terminal_mod_isomers_collapse() {
    const ACETYL: f64 = 42.01057;
    let cands = vec![
        make_candidate_seq(b"AKGK", &[(0, ACETYL)], 0),
        make_candidate_seq(b"AKGK", &[(1, ACETYL)], 0),
    ];
    let rows = emitted(
        vec![make_psm(0, 9.0, 9.0, 1, 2), make_psm(0, 8.0, 8.0, 0, 2)],
        &cands,
    );
    assert_eq!(rows.len(), 1, "{rows:?}");
    assert_eq!(rows[0].1, "K.AK+42.01057GK.S");
}

#[test]
fn tied_isomers_keep_the_first_in_emission_order() {
    let cands = vec![
        make_candidate_seq(b"SATSK", &[(0, PHOSPHO)], 0),
        make_candidate_seq(b"SATSK", &[(3, PHOSPHO)], 0),
        make_candidate_seq(b"GGGGK", &[], 0),
    ];
    let psms = vec![
        make_psm(0, 10.0, 10.0, 0, 2),
        make_psm(0, 10.0, 10.0, 1, 2),
        make_psm(0, 10.0, 10.0, 2, 2),
    ];
    let mut queue = TopNQueue::new(10);
    for p in psms.clone() {
        queue.push(p);
    }
    let order: Vec<u32> = queue
        .into_rank_sorted_vec()
        .iter()
        .map(|p| p.primary_candidate_idx())
        .collect();
    let first_isomer = *order.iter().find(|&&c| c != 2).unwrap();
    let rows = emitted(psms, &cands);
    assert_eq!(rows.len(), 2, "{rows:?}");
    let kept_isomer = rows.iter().find(|r| r.1.starts_with("K.S")).unwrap();
    assert_eq!(
        kept_isomer.1,
        cands[first_isomer as usize].peptide.to_string()
    );
    let ids: Vec<&str> = rows.iter().map(|r| r.0.as_str()).collect();
    assert_eq!(ids, ["Scan 1_1_1_0", "Scan 1_1_1_1"]);
}

#[test]
fn target_and_decoy_isomers_are_kept_apart() {
    let mut decoy = make_candidate_seq(b"SATSK", &[(3, PHOSPHO)], 1);
    decoy.is_decoy = true;
    let cands = vec![make_candidate_seq(b"SATSK", &[(0, PHOSPHO)], 0), decoy];
    let rows = emitted(
        vec![make_psm(0, 10.0, 10.0, 0, 2), make_psm(0, 9.0, 9.0, 1, 2)],
        &cands,
    );
    assert_eq!(rows.len(), 2, "{rows:?}");
}

#[test]
fn collapsed_isomers_add_their_proteins_to_the_kept_row() {
    let cands = vec![
        make_candidate_seq(b"SATSK", &[(0, PHOSPHO)], 0),
        make_candidate_seq(b"SATSK", &[(3, PHOSPHO)], 0),
        make_candidate_seq(b"SATSK", &[(3, PHOSPHO)], 1),
    ];
    let params = make_params(2..=3);
    let spectra = vec![make_spectrum("Scan 1", 1, 500.0)];
    let mut queue = TopNQueue::new(10);
    queue.push(make_psm(0, 10.0, 10.0, 0, 2));
    queue.push(make_psm(0, 9.0, 9.0, 1, 2));
    queue.push(make_psm(0, 8.0, 8.0, 2, 2));
    let mut buf = Vec::<u8>::new();
    write_pin_to(
        &mut buf,
        &spectra,
        &[queue],
        &cands,
        &params,
        &make_empty_search_index(),
    )
    .unwrap();
    let header = parse_header(&buf);
    let pep_col = header.iter().position(|h| h == "Peptide").unwrap();
    let rows = parse_rows(&buf);
    assert_eq!(rows.len(), 1, "{rows:?}");
    assert_eq!(rows[0][pep_col], "K.S+79.96633ATSK.S");
    assert_eq!(rows[0][pep_col + 1..], ["PROT_0", "PROT_1"]);
}

#[test]
fn pin_handles_empty_queue() {
    let params = make_params(2..=3);
    let spectra = vec![make_spectrum("Scan 1", 1, 500.0)];
    let queues = vec![TopNQueue::new(10)];
    let idx = make_empty_search_index();

    let mut buf = Vec::<u8>::new();
    let cands: Vec<Candidate> = vec![];
    write_pin_to(&mut buf, &spectra, &queues, &cands, &params, &idx).unwrap();

    let rows = parse_rows(&buf);
    assert!(rows.is_empty(), "empty queue should produce no data rows");
}

#[test]
fn pin_writes_real_accession_when_search_index_provided() {
    let accession = "sp|P02769|ALBU_BOVIN";
    let idx = make_search_index(accession);

    let params = make_params(2..=3);
    let spectra = vec![make_spectrum("Scan 1", 1, 500.0)];

    let psm = make_psm(0, 10.0, 10.0, 0, 2);

    let mut queue = TopNQueue::new(10);
    queue.push(psm);
    let queues = vec![queue];

    let mut buf = Vec::<u8>::new();
    let cands = vec![make_candidate(0, false)];
    write_pin_to(&mut buf, &spectra, &queues, &cands, &params, &idx).unwrap();

    let cols = parse_header(&buf);
    let rows = parse_rows(&buf);
    assert_eq!(rows.len(), 1);

    let prot_idx = cols
        .iter()
        .position(|c| c == "Proteins")
        .expect("Proteins column missing");
    assert_eq!(
        rows[0][prot_idx], accession,
        "Proteins column should contain the real accession, not a PROT_N placeholder"
    );
}

#[test]
fn pin_writes_decoy_prefix_for_decoy_protein() {
    let accession = "sp|P02769|ALBU_BOVIN";
    let idx = make_search_index(accession);

    let params = make_params(2..=3);
    let spectra = vec![make_spectrum("Scan 1", 1, 500.0)];

    let psm = make_psm(0, 10.0, 10.0, 0, 2);

    let mut queue = TopNQueue::new(10);
    queue.push(psm);
    let queues = vec![queue];

    let mut buf = Vec::<u8>::new();
    let cands = vec![make_candidate(1, true)];
    write_pin_to(&mut buf, &spectra, &queues, &cands, &params, &idx).unwrap();

    let cols = parse_header(&buf);
    let rows = parse_rows(&buf);
    assert_eq!(rows.len(), 1);

    let prot_idx = cols
        .iter()
        .position(|c| c == "Proteins")
        .expect("Proteins column missing");
    let expected_decoy = format!("XXX_{}", accession);
    assert_eq!(
        rows[0][prot_idx], expected_decoy,
        "Proteins column should carry decoy prefix for decoy PSM"
    );
}

#[test]
fn pin_emits_real_num_matched_main_ions_when_features_populated() {
    let params = make_params(2..=3);
    let spectra = vec![make_spectrum("Scan 1", 1, 500.0)];

    let mut psm = make_psm(0, 10.0, 10.0, 0, 2);
    psm.features.num_matched_main_ions = 5;

    let mut queue = TopNQueue::new(10);
    queue.push(psm);
    let queues = vec![queue];
    let idx = make_empty_search_index();

    let mut buf = Vec::<u8>::new();
    let cands = vec![make_candidate(0, false)];
    write_pin_to(&mut buf, &spectra, &queues, &cands, &params, &idx).unwrap();

    let cols = parse_header(&buf);
    let rows = parse_rows(&buf);
    assert_eq!(rows.len(), 1);

    let col_idx = cols
        .iter()
        .position(|c| c == "NumMatchedMainIons")
        .expect("NumMatchedMainIons column missing");
    assert_eq!(
        rows[0][col_idx], "5",
        "NumMatchedMainIons should be 5, not zero-stubbed"
    );
}

#[test]
fn pin_emits_delta_raw_score_on_rank1_only() {
    let params = make_params(2..=3);
    let spectra = vec![make_spectrum("Scan 1", 1, 500.0)];

    let mut psm1 = make_psm(0, 12.0, 12.0, 0, 2);
    psm1.features.delta_raw_score = 7.0;
    let mut psm2 = make_psm(0, 5.0, 5.0, 1, 2);
    psm2.features.delta_raw_score = 7.0;

    let mut queue = TopNQueue::new(10);
    queue.push(psm1);
    queue.push(psm2);
    let queues = vec![queue];
    let idx = make_empty_search_index();

    let mut buf = Vec::<u8>::new();
    let cands = vec![make_candidate(0, false), make_candidate_seq(b"G", &[], 1)];
    write_pin_to(&mut buf, &spectra, &queues, &cands, &params, &idx).unwrap();

    let cols = parse_header(&buf);
    let rows = parse_rows(&buf);
    assert_eq!(rows.len(), 2, "two PSMs → two rows");

    let col_idx = cols
        .iter()
        .position(|c| c == "DeltaRankScore")
        .expect("DeltaRankScore column missing");

    let r1: f64 = rows[0][col_idx]
        .parse()
        .expect("rank-1 DeltaRankScore numeric");
    let r2: f64 = rows[1][col_idx]
        .parse()
        .expect("rank-2 DeltaRankScore numeric");
    assert!(
        (r1 - 7.0).abs() < 1e-6,
        "rank-1 DeltaRankScore should be 7.0, got {r1}"
    );
    assert_eq!(
        r2, 0.0,
        "rank-2 DeltaRankScore should be gated to 0.0, got {r2}"
    );
}

#[test]
fn pin_emits_longest_y_pct_with_six_decimals() {
    let params = make_params(2..=3);
    let spectra = vec![make_spectrum("Scan 1", 1, 500.0)];

    let mut psm = make_psm(0, 10.0, 10.0, 0, 2);
    psm.features.longest_y = 1;
    psm.features.longest_y_pct = 0.5;

    let mut queue = TopNQueue::new(10);
    queue.push(psm);
    let queues = vec![queue];
    let idx = make_empty_search_index();

    let mut buf = Vec::<u8>::new();
    let cands = vec![make_candidate(0, false)];
    write_pin_to(&mut buf, &spectra, &queues, &cands, &params, &idx).unwrap();

    let cols = parse_header(&buf);
    let rows = parse_rows(&buf);
    assert_eq!(rows.len(), 1);

    let col_idx = cols
        .iter()
        .position(|c| c == "longest_y_pct")
        .expect("longest_y_pct column missing");
    assert_eq!(
        rows[0][col_idx], "0.500000",
        "longest_y_pct should be formatted with 6 decimal places"
    );
}
