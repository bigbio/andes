//! The out-of-core candidate index against independent references.
//!
//! * Candidate lists: the index's full candidate list must equal
//!   `enumerate_candidates` over the whole database with collision decoys
//!   relabelled globally — same peptidoforms, flags and ORDER (the queue's
//!   tie-breaks follow candidate order).
//! * PSMs: small searches (real peaks + a small FASTA; a variable mod, a
//!   semi-tryptic digest, an asymmetric precursor tolerance with a fixed mod,
//!   and a palindromic collision decoy) must reproduce golden PSM lists:
//!   peptide (residues + per-residue mod mass), charge, score, rank_score,
//!   isotope offset, decoy flag and, for the bit-identity golden, the resolved
//!   proteins (the PIN/TSV resolution of `candidate_idxs`).
//!
//! The PSM goldens were recorded on 2026-10-11 from the in-memory candidate list
//! (whole-database enumeration + a mass-bucket lookup), the path these searches
//! took before the out-of-core index became the only one; at that time the
//! out-of-core path was asserted equal to it on every one of these fixtures.

use std::collections::HashSet;

use rustc_hash::FxHashMap;

use model::{
    activation::ActivationMethod, instrument::InstrumentType, protocol::Protocol, AminoAcid,
    AminoAcidSetBuilder, ModLocation, Modification, Peptide, Protein, ProteinDb, ResidueSpec,
    Spectrum, Tolerance, H2O, PROTON,
};
use scoring_crate::param_model::{FragmentOffsetFrequency, IonType, Partition, SpecDataType};
use scoring_crate::scoring::fragment_ions::predict_by_ions;
use scoring_crate::{Param, RankScorer};
use search::candidate_gen::{enumerate_candidates, Candidate};
use search::psm::TopNQueue;
use search::{PreparedSearch, SearchIndex, SearchParams};

/// Realistic-enough scorer so candidates score on placed b/y peaks.
fn make_scorer(tol_da: f64) -> RankScorer {
    let part = Partition {
        charge: 2,
        parent_mass: 0.0,
        seg_num: 0,
    };
    let prefix1 = IonType::Prefix {
        charge: 1,
        offset_bits: (PROTON as f32).to_bits(),
        loss_class: 0,
    };
    let suffix1 = IonType::Suffix {
        charge: 1,
        offset_bits: ((H2O + PROTON) as f32).to_bits(),
        loss_class: 0,
    };
    let noise = IonType::Noise;
    let mut ion_table = FxHashMap::default();
    ion_table.insert(prefix1, vec![0.6_f32, 0.3, 0.05, 0.001]);
    ion_table.insert(suffix1, vec![0.6_f32, 0.3, 0.05, 0.001]);
    ion_table.insert(noise, vec![0.1_f32, 0.2, 0.3, 0.4]);
    let mut rank_dist_table = FxHashMap::default();
    rank_dist_table.insert(part, ion_table);
    let mut frag_off_table = FxHashMap::default();
    frag_off_table.insert(
        part,
        vec![
            FragmentOffsetFrequency {
                ion_type: prefix1,
                frequency: 0.7,
            },
            FragmentOffsetFrequency {
                ion_type: suffix1,
                frequency: 0.7,
            },
        ],
    );
    let mut param = Param {
        version: 10001,
        data_type: SpecDataType {
            activation: ActivationMethod::HCD,
            instrument: InstrumentType::QExactive,
            enzyme: None,
            protocol: Protocol::Automatic,
        },
        mme: Tolerance::Da(tol_da),
        apply_deconvolution: false,
        deconvolution_error_tolerance: 0.0,
        charge_hist: vec![(2, 100)],
        min_charge: 2,
        max_charge: 2,
        num_segments: 1,
        partitions: vec![part],
        num_precursor_off: 0,
        precursor_off_map: FxHashMap::default(),
        frag_off_table,
        max_rank: 3,
        rank_dist_table,
        error_scaling_factor: 0,
        ion_err_dist_table: FxHashMap::default(),
        noise_err_dist_table: FxHashMap::default(),
        ion_existence_table: FxHashMap::default(),
        partition_ion_types_cache: FxHashMap::default(),
        gbdt_peak_model: None,
        frag_intensity_model: None,
        rich_ion_model: None,
    };
    param.rebuild_cache();
    RankScorer::new(&param)
}

/// Build a peptide directly from residue bytes (no flanks of consequence here).
fn residues(bytes: &[u8]) -> Vec<AminoAcid> {
    bytes
        .iter()
        .map(|&b| AminoAcid::standard(b).unwrap())
        .collect()
}

/// Place charge-1 b/y peaks for `pep` so the candidate scores; add background.
fn peaks_for(pep: &Peptide) -> Vec<(f64, f32)> {
    let mut peaks: Vec<(f64, f32)> = predict_by_ions(pep, 1..=1)
        .iter()
        .enumerate()
        .map(|(i, p)| (p.mz, 100.0 - i as f32))
        .collect();
    peaks.push((37.5, 3.0)); // background
    peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    peaks
}

fn spectrum_for(pep: &Peptide, charge: u8, title: &str) -> Spectrum {
    let mz = (pep.mass() + charge as f64 * PROTON) / charge as f64;
    Spectrum {
        title: title.into(),
        precursor_mz: mz,
        precursor_intensity: None,
        precursor_charge: Some(charge as i32),
        rt_seconds: None,
        scan: None,
        peaks: peaks_for(pep),
        activation_method: None,
        isolation_lower_offset: None,
        isolation_upper_offset: None,
        isolation_target_mz: None,
    }
}

/// Small search: two proteins; spectra targeting an unmodified peptide, an
/// Oxidation-M modified peptide, and a peptide shared across both proteins.
fn small_search_fixture() -> (Vec<Spectrum>, SearchIndex, SearchParams) {
    // PEPTMIDEK contains an M (Ox-M), WVTFISLLR is a clean tryptic peptide,
    // SHAREDPEPK appears in BOTH proteins (multi-protein aggregation).
    let target = ProteinDb {
        proteins: vec![
            Protein {
                accession: "P1".into(),
                description: "fixture one".into(),
                sequence: b"MKPEPTMIDEKWVTFISLLRSHAREDPEPK".to_vec(),
                annotation: None,
            },
            Protein {
                accession: "P2".into(),
                description: "fixture two".into(),
                sequence: b"AAGTLLNRSHAREDPEPKGGGR".to_vec(),
                annotation: None,
            },
        ],
    };
    let idx = SearchIndex::from_target_db(&target, "XXX");

    let ox_m = Modification {
        name: "Oxidation".to_string(),
        mass_delta: 15.99491,
        residue: ResidueSpec::Specific(b'M'),
        location: ModLocation::Anywhere,
        fixed: false,
        accession: None,
        neutral_losses: Vec::new(),
        loss_class: 0,
        annotated_site_mass: None,
        blocks_cleavage: false,
    };
    let aa_set = AminoAcidSetBuilder::new_standard()
        .add_variable_mod(ox_m.clone())
        .build()
        .unwrap();
    let mut params = SearchParams::default_tryptic(aa_set);
    params.min_length = 3;
    params.max_variable_mods_per_peptide = 1;
    params.top_n_psms_per_spectrum = 5;
    params.min_peaks = 0;

    // Spectrum targeting the Ox-M modified PEPTMIDEK.
    let ox_residues = {
        let mut r = residues(b"PEPTMIDEK");
        // place Ox on the M (index 4)
        r[4].mod_ = Some(model::modification::leak_mod(ox_m));
        r
    };
    let ox_pep = Peptide::new(ox_residues, b'K', b'W');
    // Unmodified WVTFISLLR.
    let clean_pep = Peptide::new(residues(b"WVTFISLLR"), b'K', b'S');
    // Shared peptide (both proteins).
    let shared_pep = Peptide::new(residues(b"SHAREDPEPK"), b'R', b'-');

    let spectra = vec![
        spectrum_for(&ox_pep, 2, "scan=ox"),
        spectrum_for(&clean_pep, 2, "scan=clean"),
        spectrum_for(&shared_pep, 2, "scan=shared"),
    ];
    (spectra, idx, params)
}

/// Asymmetric precursor tolerance (wide on the `left` = "peptide lighter than
/// precursor" side, tight on the `right`) plus a FIXED carbamidomethyl-C and NO
/// variable mods. Several near-isobaric tryptic peptides sit on the WIDE side of
/// the precursor mass, where a symmetric `max(left,right)` window would admit a
/// different candidate set than the directional `matches_precursor` gate.
fn asymmetric_cam_only_fixture() -> (Vec<Spectrum>, SearchIndex, SearchParams) {
    // One protein with several tryptic peptides whose neutral masses are close
    // together, so the asymmetric window's two sides admit different sets.
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "cam-only asymmetric".into(),
            // Tryptic peptides (cut after K/R): ELVISCLIVEK, SAMPLERPEPTIDEK,
            // DAVIDCMENGEK, GASLYCDEFGHIK, ... several Cys-bearing peptides so the
            // fixed CAM mass shifts each, and several near-mass neighbours.
            sequence: b"ELVISCLIVEKSAMPLERPEPTIDEKDAVIDCMENGEKGASLYCDEFGHIKWANDACEDFGHIK".to_vec(),
            annotation: None,
        }],
    };
    let idx = SearchIndex::from_target_db(&target, "XXX");

    // Fixed carbamidomethyl on C (cam_only): NOT a variable mod.
    let cam = Modification {
        name: "Carbamidomethyl".to_string(),
        mass_delta: 57.02146,
        residue: ResidueSpec::Specific(b'C'),
        location: ModLocation::Anywhere,
        fixed: true,
        accession: None,
        neutral_losses: Vec::new(),
        loss_class: 0,
        annotated_site_mass: None,
        blocks_cleavage: false,
    };
    let aa_set = AminoAcidSetBuilder::new_standard()
        .add_fixed_mod(cam)
        .build()
        .unwrap();
    let mut params = SearchParams::default_tryptic(aa_set);
    params.min_length = 3;
    params.max_variable_mods_per_peptide = 0; // cam_only: no variable mods
    params.top_n_psms_per_spectrum = 10;
    params.min_peaks = 0;
    // STRONGLY ASYMMETRIC precursor tolerance: wide left, tight right. Da units so
    // the window is mass-independent and the directional asymmetry is stark. The
    // wide left side reaches peptides up to ~0.6 Da LIGHTER than the precursor;
    // the tight right side admits almost nothing heavier.
    params.precursor_tolerance =
        model::PrecursorTolerance::asymmetric(Tolerance::Da(0.6), Tolerance::Da(0.02));
    // Multiple isotope offsets exercise the per-offset directional window.
    params.isotope_error_range = -1..=2;

    // Build spectra: for EACH tryptic Cys-bearing peptide, place a precursor a
    // little HEAVIER than the peptide so the peptide falls on the WIDE (left)
    // side, and nearby lighter neighbours land in/out of the asymmetric window
    // differently than a symmetric one would.
    let build_pep = |seq: &[u8], pre: u8, post: u8| -> Peptide {
        let mut r = residues(seq);
        // Apply the fixed CAM to every C so the spectrum's precursor matches the
        // CAM-modified candidate the search will score.
        let cam_mod = Modification {
            name: "Carbamidomethyl".to_string(),
            mass_delta: 57.02146,
            residue: ResidueSpec::Specific(b'C'),
            location: ModLocation::Anywhere,
            fixed: true,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
            annotated_site_mass: None,
            blocks_cleavage: false,
        };
        for aa in r.iter_mut() {
            if aa.residue == b'C' {
                aa.mod_ = Some(model::modification::leak_mod(cam_mod.clone()));
            }
        }
        Peptide::new(r, pre, post)
    };

    let peps = [
        build_pep(b"ELVISCLIVEK", b'_', b'S'),
        build_pep(b"DAVIDCMENGEK", b'R', b'G'),
        build_pep(b"GASLYCDEFGHIK", b'K', b'W'),
        build_pep(b"WANDACEDFGHIK", b'K', b'-'),
    ];
    let mut spectra = Vec::new();
    for (i, pep) in peps.iter().enumerate() {
        // Nudge the precursor +0.3 Da heavier (within wide left, outside tight
        // right) so the matching peptide sits on the asymmetric WIDE side.
        let charge = 2u8;
        let mz = (pep.mass() + 0.3 + charge as f64 * PROTON) / charge as f64;
        let mut s = spectrum_for(pep, charge, &format!("scan=asym{i}"));
        s.precursor_mz = mz;
        spectra.push(s);
    }
    (spectra, idx, params)
}

/// A palindromic protein: its tryptic peptide `MAEK` reverses to ITSELF, so the
/// reverse decoy's bare sequence equals a real target's — a collision decoy that
/// must be relabelled as target using the GLOBAL target set.
fn palindrome_fixture() -> (Vec<Spectrum>, SearchIndex, SearchParams) {
    // Protein "MAEKKEAM" is a palindrome (reverse == itself), so its tryptic
    // peptide "MAEK" (cut after K) collides with its own reverse-decoy "MAEK".
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "palindrome collision".into(),
            sequence: b"MAEKKEAM".to_vec(),
            annotation: None,
        }],
    };
    let idx = SearchIndex::from_target_db(&target, "XXX");

    let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
    let mut params = SearchParams::default_tryptic(aa_set);
    params.min_length = 3;
    params.max_variable_mods_per_peptide = 0;
    params.top_n_psms_per_spectrum = 10;
    params.min_peaks = 0;

    // Spectrum targeting MAEK (the collision peptide).
    let maek = Peptide::new(residues(b"MAEK"), b'_', b'K');
    let spectra = vec![spectrum_for(&maek, 2, "scan=collision")];
    (spectra, idx, params)
}

/// Search `spectra` on the out-of-core path with per-spectrum enumeration
/// (`fragment_index_top_k` is 0 in these fixtures). Returns the per-scan queues
/// PLUS the materialized candidates the PSM `candidate_idxs` resolve against
/// (just like the binary holds `prepared.candidates`).
fn run_search(
    idx: &SearchIndex,
    params: &SearchParams,
    spectra: &[Spectrum],
    scorer: &RankScorer,
) -> (Vec<TopNQueue>, Vec<Candidate>) {
    let tmp = tempfile::NamedTempFile::new().unwrap();
    let mut prepared = PreparedSearch::prepare_mmap(idx, params, scorer, 0.05, "XXX", tmp.path())
        .expect("prepare_mmap");
    let queues = prepared.run_chunk(spectra, 0);
    prepared.sync_materialized_candidates();
    (queues, std::mem::take(&mut prepared.candidates))
}

fn bare(c: &Candidate) -> Vec<u8> {
    c.peptide.residues.iter().map(|a| a.residue).collect()
}

/// The whole database's candidate list as `enumerate_candidates` yields it, with
/// every decoy whose bare sequence is also a target's relabelled as target.
fn enumeration_oracle(idx: &SearchIndex, params: &SearchParams) -> Vec<Candidate> {
    let mut all: Vec<Candidate> = enumerate_candidates(idx, params, "XXX").collect();
    let targets: HashSet<Vec<u8>> = all.iter().filter(|c| !c.is_decoy).map(bare).collect();
    for c in all.iter_mut().filter(|c| c.is_decoy) {
        if targets.contains(&bare(c)) {
            c.is_decoy = false;
        }
    }
    all
}

/// Everything a candidate carries into scoring and the PIN.
type CandSig = (usize, usize, Vec<(u8, i64)>, u8, u8, bool, bool, bool);

fn cand_sig(c: &Candidate) -> CandSig {
    (
        c.protein_index,
        c.start_offset_in_protein,
        c.peptide
            .residues
            .iter()
            .map(|a| {
                let d = a.mod_.as_ref().map_or(0.0, |m| m.mass_delta);
                (a.residue, (d * 1e5).round() as i64)
            })
            .collect(),
        c.peptide.pre,
        c.peptide.post,
        c.is_decoy,
        c.is_protein_n_term,
        c.is_protein_c_term,
    )
}

/// The out-of-core index lists exactly the enumeration's candidates, with the
/// same collision-decoy labels, on a variable-mod search with a shared peptide, a
/// semi-tryptic digest, a fixed mod, and a palindrome whose reverse decoy
/// collides with a target. The order matches too, except semi-tryptic: the
/// enumeration lists a protein's strict spans before its semi-specific ones,
/// the index lists spans by start and length.
#[test]
fn out_of_core_candidate_list_equals_the_enumeration() {
    let scorer = make_scorer(0.05);
    let (_, idx, params) = small_search_fixture();
    let mut semi = params.clone();
    semi.num_tolerable_termini = 1;
    let (_, asym_idx, asym_params) = asymmetric_cam_only_fixture();
    let (_, pal_idx, pal_params) = palindrome_fixture();
    let mut relabelled = 0;
    for (label, idx, params, ordered) in [
        ("variable mod", &idx, &params, true),
        ("semi-tryptic", &idx, &semi, false),
        ("fixed mod", &asym_idx, &asym_params, true),
        ("palindrome", &pal_idx, &pal_params, true),
    ] {
        let tmp = tempfile::NamedTempFile::new().unwrap();
        let prepared = PreparedSearch::prepare_mmap(idx, params, &scorer, 0.05, "XXX", tmp.path())
            .expect("prepare_mmap");
        let mut got: Vec<CandSig> = prepared
            .out_of_core_candidates_where(|_, _| true)
            .iter()
            .map(cand_sig)
            .collect();
        let oracle = enumeration_oracle(idx, params);
        let raw: Vec<Candidate> = enumerate_candidates(idx, params, "XXX").collect();
        relabelled += raw
            .iter()
            .zip(&oracle)
            .filter(|(r, o)| r.is_decoy != o.is_decoy)
            .count();
        let mut want: Vec<CandSig> = oracle.iter().map(cand_sig).collect();
        assert!(!want.is_empty(), "[{label}] empty fixture");
        if ordered {
            assert_eq!(
                got, want,
                "[{label}] candidate list differs from the enumeration"
            );
        }
        got.sort();
        want.sort();
        assert_eq!(
            got, want,
            "[{label}] candidate list differs from the enumeration"
        );
    }
    assert!(relabelled > 0, "no fixture exercises the collision relabel");
}

/// Protein 0 starts with M; each of its peptides is enumerated once, so it is
/// listed once in `prots`.
const FIXTURE_SIGNATURE: &[&[&str]] = &[
    &[
        "pep=E[0]D[0]I[0]M[1599491]T[0]P[0]E[0]P[0]K[0] z=2 score=-87.000000 rank=-87.000000 iso=0 decoy=true prots=[2]",
        "pep=P[0]E[0]P[0]T[0]M[1599491]I[0]D[0]E[0]K[0] z=2 score=-43.000000 rank=-43.000000 iso=0 decoy=false prots=[0]",
    ],
    &[
        "pep=W[0]V[0]T[0]F[0]I[0]S[0]L[0]L[0]R[0] z=2 score=-53.000000 rank=-53.000000 iso=0 decoy=false prots=[0]",
    ],
    &[
        "pep=S[0]H[0]A[0]R[0]E[0]D[0]P[0]E[0]P[0]K[0] z=2 score=-43.000000 rank=-43.000000 iso=0 decoy=false prots=[0,1]",
    ],
];

const SEMITRYPTIC: &[&[&str]] = &[
    &[
        "pep=E[0]D[0]I[0]M[1599491]T[0]P[0]E[0]P[0]K[0] z=2 score=-87.0000 rank=-87.0000 iso=0 decoy=true",
        "pep=K[0]P[0]E[0]P[0]T[0]M[1599491]I[0]D[0]E[0] z=2 score=-97.0000 rank=-97.0000 iso=0 decoy=false",
        "pep=P[0]E[0]P[0]T[0]M[1599491]I[0]D[0]E[0]K[0] z=2 score=-43.0000 rank=-43.0000 iso=0 decoy=false",
    ],
    &[
        "pep=W[0]V[0]T[0]F[0]I[0]S[0]L[0]L[0]R[0] z=2 score=-53.0000 rank=-53.0000 iso=0 decoy=false",
    ],
    &[
        "pep=S[0]H[0]A[0]R[0]E[0]D[0]P[0]E[0]P[0]K[0] z=2 score=-43.0000 rank=-43.0000 iso=0 decoy=false",
    ],
];

const ASYMMETRIC_TOL_CAM_ONLY: &[&[&str]] = &[
    &[
        "pep=E[0]L[0]V[0]I[0]S[0]C[5702146]L[0]I[0]V[0]E[0]K[0] z=2 score=-61.0000 rank=-61.0000 iso=0 decoy=false",
    ],
    &[
        "pep=D[0]A[0]V[0]I[0]D[0]C[5702146]M[0]E[0]N[0]G[0]E[0]K[0] z=2 score=-86.0000 rank=-86.0000 iso=0 decoy=false",
        "pep=E[0]G[0]N[0]E[0]M[0]C[5702146]D[0]I[0]V[0]A[0]D[0]K[0] z=2 score=-115.0000 rank=-115.0000 iso=0 decoy=true",
    ],
    &[
        "pep=G[0]A[0]S[0]L[0]Y[0]C[5702146]D[0]E[0]F[0]G[0]H[0]I[0]K[0] z=2 score=-94.0000 rank=-94.0000 iso=0 decoy=false",
        "pep=I[0]H[0]G[0]F[0]E[0]D[0]C[5702146]Y[0]L[0]S[0]A[0]G[0]K[0] z=2 score=-135.0000 rank=-135.0000 iso=0 decoy=true",
    ],
    &[
        "pep=I[0]H[0]G[0]F[0]D[0]E[0]C[5702146]A[0]D[0]N[0]A[0]W[0]K[0] z=2 score=-135.0000 rank=-135.0000 iso=0 decoy=true",
        "pep=W[0]A[0]N[0]D[0]A[0]C[5702146]E[0]D[0]F[0]G[0]H[0]I[0]K[0] z=2 score=-103.0000 rank=-103.0000 iso=0 decoy=false",
    ],
];

const COLLISION_DECOY_RELABEL: &[&[&str]] = &[&[
    "pep=K[0]E[0]A[0]M[0] z=2 score=-37.0000 rank=-37.0000 iso=0 decoy=false",
    "pep=M[0]A[0]E[0]K[0] z=2 score=-5.0000 rank=-5.0000 iso=0 decoy=false",
]];

fn assert_golden(label: &str, got: &[Vec<String>], golden: &[&[&str]]) {
    assert!(
        got.iter().any(|scan| !scan.is_empty()),
        "[{label}] the fixture produced no PSMs, so the golden would be vacuous"
    );
    let want: Vec<Vec<String>> = golden
        .iter()
        .map(|scan| scan.iter().map(|r| r.to_string()).collect())
        .collect();
    assert_eq!(
        got, want,
        "[{label}] PSMs differ from the golden\nGOT : {got:#?}\nWANT: {want:#?}"
    );
}

/// Variable mod + a peptide shared by two proteins: full PSM rows including the
/// resolved proteins.
#[test]
fn psms_match_golden_with_variable_mod_and_shared_peptide() {
    let (spectra, idx, params) = small_search_fixture();
    let scorer = make_scorer(0.05);
    let (q, c) = run_search(&idx, &params, &spectra, &scorer);
    assert_golden("variable mod", &psm_signature(&q, &c), FIXTURE_SIGNATURE);
}

/// Semi-tryptic digest (`num_tolerable_termini = 1`).
#[test]
fn psms_match_golden_semitryptic() {
    let (spectra, idx, mut params) = small_search_fixture();
    params.num_tolerable_termini = 1;
    let scorer = make_scorer(0.05);
    let (q, c) = run_search(&idx, &params, &spectra, &scorer);
    assert_golden("semi-tryptic", &psm_result_set(&q, &c), SEMITRYPTIC);
}

/// Directional (asymmetric) precursor tolerance with a fixed mod: the scored
/// set must follow `matches_precursor`'s sides.
#[test]
fn psms_match_golden_asymmetric_tol_cam_only() {
    let (spectra, idx, params) = asymmetric_cam_only_fixture();
    let scorer = make_scorer(0.05);
    let (q, c) = run_search(&idx, &params, &spectra, &scorer);
    assert_golden(
        "asymmetric tolerance",
        &psm_result_set(&q, &c),
        ASYMMETRIC_TOL_CAM_ONLY,
    );
}

/// The palindrome's collision decoy is labelled as target in the PSMs.
#[test]
fn psms_match_golden_collision_decoy_relabel() {
    let (spectra, idx, params) = palindrome_fixture();
    let scorer = make_scorer(0.05);
    let (q, c) = run_search(&idx, &params, &spectra, &scorer);
    let got = psm_result_set(&q, &c);
    assert!(
        got.iter()
            .flatten()
            .any(|r| r.starts_with("pep=M[0]A[0]E[0]K[0] ") && r.ends_with("decoy=false")),
        "the collision MAEK must be scored as a target: {got:#?}"
    );
    assert_golden("collision decoy", &got, COLLISION_DECOY_RELABEL);
}

/// Canonical, resolution-aware signature of every scan's top-N PSM list.
/// Includes peptide string (residues + per-residue mod mass), charge, score,
/// rank_score, isotope_offset, and the resolved protein accessions (the PIN/TSV
/// `candidate_idxs` resolution). Sorted within scan so heap/order noise can't
/// hide a real divergence — but the PSM SET and resolution must match.
fn psm_signature(queues: &[TopNQueue], candidates: &[Candidate]) -> Vec<Vec<String>> {
    queues
        .iter()
        .map(|q| {
            let mut rows: Vec<String> = q
                .iter_psms()
                .map(|psm| {
                    let cand = &candidates[psm.primary_candidate_idx() as usize];
                    let pep: String = cand
                        .peptide
                        .residues
                        .iter()
                        .map(|aa| {
                            let m = aa
                                .mod_
                                .as_ref()
                                .map(|m| (m.mass_delta * 1e5).round() as i64)
                                .unwrap_or(0);
                            format!("{}[{}]", aa.residue as char, m)
                        })
                        .collect();
                    // Resolve EVERY candidate_idx into a protein accession (the
                    // PIN `Proteins` column), in the stored order.
                    let prots: Vec<String> = psm
                        .candidate_idxs
                        .iter()
                        .map(|&ci| candidates[ci as usize].protein_index.to_string())
                        .collect();
                    format!(
                        "pep={pep} z={} score={:.6} rank={:.6} iso={} decoy={} prots=[{}]",
                        psm.charge_used,
                        psm.score,
                        psm.rank_score,
                        psm.isotope_offset,
                        cand.is_decoy,
                        prots.join(","),
                    )
                })
                .collect();
            rows.sort();
            rows
        })
        .collect()
}

/// Canonical, order-independent set of accepted PSMs per scan.
///
/// Unlike `psm_signature` (which compares sorted rows — still sensitive to ordering
/// within the heap since `candidate_idxs` order affects the signature), this
/// function builds a per-scan **multiset** of `(peptide-residues+mods, charge,
/// score-rounded, rank-score-rounded, isotope_offset, is_decoy)` tuples and
/// sorts them. The protein column is deliberately excluded: for semi-tryptic
/// searches the multi-protein `Proteins` order depends on candidate visitation
/// order (cosmetic, FDR-neutral), and these tests assert result identity.
fn psm_result_set(queues: &[TopNQueue], candidates: &[Candidate]) -> Vec<Vec<String>> {
    queues
        .iter()
        .map(|q| {
            let mut rows: Vec<String> = q
                .iter_psms()
                .map(|psm| {
                    let cand = &candidates[psm.primary_candidate_idx() as usize];
                    let pep: String = cand
                        .peptide
                        .residues
                        .iter()
                        .map(|aa| {
                            let m = aa
                                .mod_
                                .as_ref()
                                .map(|m| (m.mass_delta * 1e5).round() as i64)
                                .unwrap_or(0);
                            format!("{}[{}]", aa.residue as char, m)
                        })
                        .collect();
                    // Round scores to 4 decimal places to tolerate any
                    // float-formatting noise while still catching real divergence.
                    format!(
                        "pep={pep} z={} score={:.4} rank={:.4} iso={} decoy={}",
                        psm.charge_used,
                        psm.score,
                        psm.rank_score,
                        psm.isotope_offset,
                        cand.is_decoy,
                    )
                })
                .collect();
            rows.sort();
            rows
        })
        .collect()
}

/// Issue #76: on spectra built from a peptide's exact b/y ions, the per-chunk
/// fragment-ion index must give that peptide the most votes and materialise it.
#[test]
fn fragment_index_ranks_the_true_peptide_first_on_exact_spectra() {
    use model::tolerance::Tolerance;
    use search::fragment_index::ChunkFragmentIndex;
    let (spectra, idx, params) = small_search_fixture();
    let scorer = make_scorer(0.05);
    let tmp = tempfile::NamedTempFile::new().unwrap();
    let prepared = PreparedSearch::prepare_mmap(&idx, &params, &scorer, 0.05, "XXX", tmp.path())
        .expect("prepare_mmap");
    let mi = &prepared.mmap_index;
    let records = mi.records();
    let fi = ChunkFragmentIndex::build(records, &idx, &params, Tolerance::Da(0.05), 0.0, 5000.0);
    assert!(fi.n_forms() > 0);
    // Reference: what the enumeration path scores as top-1 for each spectrum.
    let (queues, cands) = run_search(&idx, &params, &spectra, &scorer);
    for (spec, q) in spectra.iter().zip(&queues) {
        let best = q
            .iter_psms()
            .max_by(|a, b| a.rank_score.partial_cmp(&b.rank_score).unwrap())
            .expect("psm");
        // Modification-aware signature (residue, mod delta): the oxidised
        // PEPTMIDEK form must be distinguished from the unmodified one.
        let sig = |c: &Candidate| -> Vec<(u8, i64)> {
            c.peptide
                .residues
                .iter()
                .map(|a| {
                    let d = a.mod_.as_ref().map_or(0.0, |m| m.mass_delta);
                    (a.residue, (d * 1e4).round() as i64)
                })
                .collect()
        };
        let truth = sig(&cands[best.primary_candidate_idx() as usize]);
        let z = spec.precursor_charge.unwrap() as u8;
        let sel = fi.query(spec, &[z], &params, Tolerance::Da(0.05), 5, 3, 0);
        assert!(!sel.is_empty(), "{}: no votes", spec.title);
        let mats = fi.materialise(&sel, &idx, &params, None, &|_| {});
        let picks: Vec<Vec<(u8, i64)>> = mats.iter().map(sig).collect();
        assert!(
            picks.contains(&truth),
            "{}: enumeration top-1 peptidoform not among index picks (votes {:?})",
            spec.title,
            sel
        );
    }
}

/// Same as above with a tight mass window, so the bounded walk actually prunes
/// and the pruned index `k` must still select the same form on both sides.
#[test]
fn fragment_index_prunes_to_the_window_and_still_finds_the_true_peptide() {
    use model::tolerance::Tolerance;
    use search::fragment_index::ChunkFragmentIndex;
    let (spectra, idx, params) = small_search_fixture();
    let scorer = make_scorer(0.05);
    let tmp = tempfile::NamedTempFile::new().unwrap();
    let prepared = PreparedSearch::prepare_mmap(&idx, &params, &scorer, 0.05, "XXX", tmp.path())
        .expect("prepare_mmap");
    let mi = &prepared.mmap_index;
    let (queues, cands) = run_search(&idx, &params, &spectra, &scorer);
    for (spec, q) in spectra.iter().zip(&queues) {
        let best = q
            .iter_psms()
            .max_by(|a, b| a.rank_score.partial_cmp(&b.rank_score).unwrap())
            .expect("psm");
        let truth = &cands[best.primary_candidate_idx() as usize].peptide;
        let truth_seq: Vec<u8> = truth.residues.iter().map(|a| a.residue).collect();
        let m = truth.mass();
        let fi = ChunkFragmentIndex::build(
            mi.records(),
            &idx,
            &params,
            Tolerance::Da(0.05),
            m - 1.0,
            m + 1.0,
        );
        let z = spec.precursor_charge.unwrap() as u8;
        let sel = fi.query(spec, &[z], &params, Tolerance::Da(0.05), 5, 3, 0);
        assert!(!sel.is_empty(), "{}: no votes", spec.title);
        let mats = fi.materialise(&sel, &idx, &params, None, &|_| {});
        let seqs: Vec<Vec<u8>> = mats
            .iter()
            .map(|c| c.peptide.residues.iter().map(|a| a.residue).collect())
            .collect();
        assert!(
            seqs.contains(&truth_seq),
            "{}: truth {:?} not among picks {:?}",
            spec.title,
            String::from_utf8_lossy(&truth_seq),
            seqs.len()
        );
        for c in &mats {
            assert!(
                (c.peptide.mass() - m).abs() <= 1.0 + 1e-6,
                "materialised form outside the window"
            );
        }
    }
}

/// The single-shot out-of-core search returns one queue per input spectrum in
/// input order, each PSM's `spectrum_idx` is its input position, and its top
/// PSM is the one per-spectrum enumeration ranks first. The spectra are passed
/// in descending precursor mass, so the search has to reorder them.
#[test]
fn match_spectra_out_of_core_keeps_input_order_and_finds_the_enumeration_top_hit() {
    let (mut spectra, idx, params) = small_search_fixture();
    spectra.sort_by(|a, b| {
        search::fragment_index::index_order_mass(b)
            .total_cmp(&search::fragment_index::index_order_mass(a))
    });
    let scorer = make_scorer(0.05);
    let masses: Vec<f64> = spectra
        .iter()
        .map(search::fragment_index::index_order_mass)
        .collect();
    assert!(
        masses.windows(2).any(|w| w[0] > w[1]),
        "the fixture must exercise the mass reordering"
    );
    let (queues, cands) =
        search::match_spectra_out_of_core(&spectra, &idx, &params, &scorer, 0.05, "XXX")
            .expect("out-of-core search");
    let (enum_queues, enum_cands) = run_search(&idx, &params, &spectra, &scorer);
    assert_eq!(queues.len(), spectra.len());
    let top = |q: &TopNQueue| {
        q.iter_psms()
            .max_by(|a, b| a.rank_score.partial_cmp(&b.rank_score).unwrap())
            .cloned()
            .expect("psm")
    };
    for (i, (q, rq)) in queues.iter().zip(&enum_queues).enumerate() {
        assert!(q.iter_psms().all(|p| p.spectrum_idx == i));
        let (a, b) = (top(q), top(rq));
        assert_eq!(
            cands[a.primary_candidate_idx() as usize]
                .peptide
                .to_string(),
            enum_cands[b.primary_candidate_idx() as usize]
                .peptide
                .to_string(),
            "{}: top peptide differs from the enumeration",
            spectra[i].title
        );
    }
}

/// The glyco candidate table built from the out-of-core index must be the
/// whole database's enumeration restricted to N-X-S/T-bearing candidates: the same
/// candidates (peptidoform, flanks, protein, offset, terminal and decoy flags),
/// in the same relative order.
#[test]
fn glyco_table_from_the_index_is_the_enumeration_sequon_subset_in_order() {
    use search::candidate_gen::candidate_has_nxst_sequon;
    use search::glyco_search::GlycoCandidates;
    // Internal sequons (NGS, NLT), a sequon completed by the next residue
    // (...KAEN | TK: N-X at the peptide end, S/T after it), Met-initial proteins
    // (two enumeration walks), an oxidisable M next to a sequon, and a shared
    // peptide in two proteins.
    let target = ProteinDb {
        proteins: vec![
            Protein {
                accession: "G1".into(),
                description: "glyco one".into(),
                sequence: b"MNGSPEKLLNLTMRAENTKWVNASGRPEPNKTGK".to_vec(),
                annotation: None,
            },
            Protein {
                accession: "G2".into(),
                description: "glyco two".into(),
                sequence: b"MKNMSYKAANVTRLLNLTMRSTNGK".to_vec(),
                annotation: None,
            },
            Protein {
                accession: "G3".into(),
                description: "no sequon".into(),
                sequence: b"AAGTLLDRSHAREDPEPKGGGR".to_vec(),
                annotation: None,
            },
        ],
    };
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let ox_m = Modification {
        name: "Oxidation".to_string(),
        mass_delta: 15.99491,
        residue: ResidueSpec::Specific(b'M'),
        location: ModLocation::Anywhere,
        fixed: false,
        accession: None,
        neutral_losses: Vec::new(),
        loss_class: 0,
        annotated_site_mass: None,
        blocks_cleavage: false,
    };
    let aa_set = AminoAcidSetBuilder::new_standard()
        .add_variable_mod(ox_m)
        .build()
        .unwrap();
    let mut params = SearchParams::default_tryptic(aa_set);
    params.min_length = 3;
    params.max_missed_cleavages = 3;
    params.max_variable_mods_per_peptide = 2;
    let scorer = make_scorer(0.05);

    let all = enumeration_oracle(&idx, &params);
    let expected: Vec<&Candidate> = all
        .iter()
        .filter(|c| candidate_has_nxst_sequon(c))
        .collect();
    let tmp = tempfile::NamedTempFile::new().unwrap();
    let mmap = PreparedSearch::prepare_mmap(&idx, &params, &scorer, 0.05, "XXX", tmp.path())
        .expect("prepare_mmap");
    let table = GlycoCandidates::from_out_of_core(&mmap);

    type Sig = (usize, usize, Vec<(u8, i64)>, u8, u8, bool, bool, bool);
    let sig = |c: &Candidate| -> Sig {
        (
            c.protein_index,
            c.start_offset_in_protein,
            c.peptide
                .residues
                .iter()
                .map(|a| {
                    let d = a.mod_.as_ref().map_or(0.0, |m| m.mass_delta);
                    (a.residue, (d * 1e5).round() as i64)
                })
                .collect(),
            c.peptide.pre,
            c.peptide.post,
            c.is_decoy,
            c.is_protein_n_term,
            c.is_protein_c_term,
        )
    };
    let want: Vec<Sig> = expected.iter().map(|c| sig(c)).collect();
    let got: Vec<Sig> = table.candidates().iter().map(sig).collect();

    // Not vacuous: targets and decoys, a modified form, a flank-completed sequon,
    // and the Met-cleaved walk are all represented.
    assert!(want.iter().any(|s| s.5) && want.iter().any(|s| !s.5));
    assert!(want.iter().any(|s| s.2.iter().any(|&(_, d)| d != 0)));
    assert!(expected.iter().any(|c| {
        let r: Vec<u8> = c.peptide.residues.iter().map(|a| a.residue).collect();
        !andes_glyco::sequon::has_nxst_sequon(&r)
    }));
    assert!(expected
        .iter()
        .any(|c| c.is_protein_n_term && c.start_offset_in_protein == 1));
    assert!(want.len() < all.len());
    assert_eq!(got, want);
}
