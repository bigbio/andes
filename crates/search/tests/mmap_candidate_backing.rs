//! Task 4 decisive gate: the out-of-core (`Mmap`) candidate backing must yield
//! PSMs byte-identical to the in-RAM (`Ram`) path for the same search.
//!
//! Runs the SAME small search (real peaks + a small fasta + Oxidation-M variable
//! mod, so both unmodified AND modified candidates are scored) twice through a
//! shared `run_prepared` helper — once `Ram`, once `Mmap` — and asserts the
//! per-scan top-N PSM lists are IDENTICAL: peptide string (residues + per-residue
//! mod mass) + charge + score + rank_score + isotope_offset + the resolved
//! protein accession (i.e. the PIN/TSV resolution of `candidate_idxs`).

use rustc_hash::FxHashMap;

use model::{
    activation::ActivationMethod, instrument::InstrumentType, protocol::Protocol, AminoAcid,
    AminoAcidSetBuilder, ModLocation, Modification, Peptide, Protein, ProteinDb, ResidueSpec,
    Spectrum, Tolerance, H2O, PROTON,
};
use scoring_crate::param_model::{FragmentOffsetFrequency, IonType, Partition, SpecDataType};
use scoring_crate::scoring::fragment_ions::predict_by_ions;
use scoring_crate::{Param, RankScorer};
use search::candidate_gen::Candidate;
use search::match_engine::CandidateBacking;
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
            },
            Protein {
                accession: "P2".into(),
                description: "fixture two".into(),
                sequence: b"AAGTLLNRSHAREDPEPKGGGR".to_vec(),
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
        r[4].mod_ = Some(std::sync::Arc::new(ox_m));
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

/// Run the shared search once with the requested candidate backing. Returns the
/// per-scan queues PLUS the materialized candidates the PSM `candidate_idxs`
/// resolve against (just like the binary holds `prepared.candidates`).
fn run_prepared(
    idx: &SearchIndex,
    params: &SearchParams,
    spectra: &[Spectrum],
    backing: CandidateBacking,
    scorer: &RankScorer,
) -> (Vec<TopNQueue>, Vec<Candidate>) {
    match backing {
        CandidateBacking::Ram => {
            let prepared = PreparedSearch::prepare(idx, params, scorer, 0.05, "XXX");
            let queues = prepared.run_chunk(spectra, 0);
            (queues, prepared.candidates)
        }
        CandidateBacking::Mmap => {
            let tmp = tempfile::NamedTempFile::new().unwrap();
            let mut prepared =
                PreparedSearch::prepare_mmap(idx, params, scorer, 0.05, "XXX", tmp.path())
                    .expect("prepare_mmap");
            let queues = prepared.run_chunk(spectra, 0);
            prepared.sync_materialized_candidates();
            (queues, prepared.candidates)
        }
    }
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
/// sorts them.  The protein column is deliberately excluded: for semi-tryptic
/// searches the multi-protein `Proteins` order may differ between `Ram` and
/// `Mmap` (cosmetic, FDR-neutral), and this test asserts only result-identity.
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

/// Semi-tryptic result-identity test.
///
/// For `num_tolerable_termini = 1` the `Mmap` per-spectrum candidate order
/// (coordinate/mod-mass sort) can differ from `Ram` (strict spans → free-C →
/// free-N), so byte-identity of the full PIN is not guaranteed. This test asserts
/// the weaker but essential property: the **order-independent set** of accepted
/// PSMs (peptide residues+mods, charge, score, isotope offset, decoy flag) is
/// identical between `Ram` and `Mmap`. Protein column order is excluded (see
/// `psm_result_set`).
#[test]
fn mmap_result_identical_semitryptic() {
    let (spectra, idx, mut params) = small_search_fixture();
    // Switch to semi-tryptic: at least one terminus must be a cleavage site.
    params.num_tolerable_termini = 1;
    let scorer = make_scorer(0.05);

    let (ram_q, ram_c) = run_prepared(&idx, &params, &spectra, CandidateBacking::Ram, &scorer);
    let (mmap_q, mmap_c) = run_prepared(&idx, &params, &spectra, CandidateBacking::Mmap, &scorer);

    let ram_set = psm_result_set(&ram_q, &ram_c);
    let mmap_set = psm_result_set(&mmap_q, &mmap_c);

    // Sanity: the fixture must actually produce PSMs (else the test is vacuous).
    assert!(
        ram_set.iter().any(|scan| !scan.is_empty()),
        "semi-tryptic fixture produced no PSMs on the Ram path — test would be vacuous"
    );

    assert_eq!(
        ram_set, mmap_set,
        "mmap semi-tryptic: accepted PSM sets must be result-identical to the in-RAM path\n\
         (order-independent: candidate SET, scores, and accepted-PSM identities must match;\n\
         Proteins-column order and candidate_idxs[0] for shared peptides may legitimately differ)\n\
         RAM : {ram_set:#?}\nMMAP: {mmap_set:#?}"
    );
}

/// RC2 (directional precursor tolerance) result-identity gate.
///
/// Builds a search with a STRONGLY ASYMMETRIC precursor tolerance (wide on the
/// `left` = "peptide lighter than precursor" side, tight on the `right`) plus a
/// `cam_only` style FIXED carbamidomethyl-C mod and NO variable mods — the exact
/// configuration that diverged in the b1931 benchmark. Several near-isobaric
/// tryptic peptides sit on the WIDE (left) side of the precursor mass, where a
/// symmetric `max(left,right)` lazy fetch + symmetric final filter would admit a
/// candidate set DIFFERENT from the directional `matches_precursor` gate the RAM
/// path uses, perturbing the per-spectrum scored set (and, when a wrongly-admitted
/// candidate scores into the top-N, the accepted-PSM set itself).
///
/// Asserts the order-independent accepted-PSM set is identical RAM vs Mmap.
#[test]
fn mmap_result_identical_asymmetric_tol_cam_only() {
    use std::sync::Arc;

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

    let scorer = make_scorer(0.05);

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
        };
        for aa in r.iter_mut() {
            if aa.residue == b'C' {
                aa.mod_ = Some(Arc::new(cam_mod.clone()));
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

    let (ram_q, ram_c) = run_prepared(&idx, &params, &spectra, CandidateBacking::Ram, &scorer);
    let (mmap_q, mmap_c) = run_prepared(&idx, &params, &spectra, CandidateBacking::Mmap, &scorer);

    let ram_set = psm_result_set(&ram_q, &ram_c);
    let mmap_set = psm_result_set(&mmap_q, &mmap_c);

    assert!(
        ram_set.iter().any(|scan| !scan.is_empty()),
        "asymmetric cam_only fixture produced no PSMs on the Ram path — test would be vacuous"
    );
    assert_eq!(
        ram_set, mmap_set,
        "mmap asymmetric-tolerance cam_only: accepted PSM sets must be result-identical to RAM\n\
         (directional precursor tolerance must admit the SAME per-spectrum candidate set)\n\
         RAM : {ram_set:#?}\nMMAP: {mmap_set:#?}"
    );
}

/// RC1 (global collision-decoy relabel) result-identity gate.
///
/// A palindromic tryptic peptide reverses to ITSELF, so the reverse-decoy's bare
/// sequence equals a real target's bare sequence — a collision decoy. The in-RAM
/// path relabels it to target ONCE GLOBALLY at prepare; the Mmap path must reach
/// the SAME decision using the GLOBAL target bare-sequence set built once from the
/// full index (NOT a per-spectrum relabel from only the window candidates). This
/// test asserts the accepted-PSM set — including the decoy flag — is identical, so
/// a collision decoy is labeled the same way in both backings.
#[test]
fn mmap_result_identical_collision_decoy_relabel() {
    // Protein "MAEKKEAM" is a palindrome (reverse == itself), so its tryptic
    // peptide "MAEK" (cut after K) collides with its own reverse-decoy "MAEK".
    let target = ProteinDb {
        proteins: vec![Protein {
            accession: "P1".into(),
            description: "palindrome collision".into(),
            sequence: b"MAEKKEAM".to_vec(),
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
    let scorer = make_scorer(0.05);

    let (ram_q, ram_c) = run_prepared(&idx, &params, &spectra, CandidateBacking::Ram, &scorer);
    let (mmap_q, mmap_c) = run_prepared(&idx, &params, &spectra, CandidateBacking::Mmap, &scorer);

    let ram_set = psm_result_set(&ram_q, &ram_c);
    let mmap_set = psm_result_set(&mmap_q, &mmap_c);

    assert!(
        ram_set.iter().any(|scan| !scan.is_empty()),
        "collision fixture produced no PSMs on the Ram path — test would be vacuous"
    );
    // The collision MAEK candidate must carry the SAME decoy flag in both paths.
    assert_eq!(
        ram_set, mmap_set,
        "mmap collision-decoy relabel: accepted PSM sets (incl. decoy flag) must match RAM\n\
         RAM : {ram_set:#?}\nMMAP: {mmap_set:#?}"
    );
}

#[test]
fn mmap_path_bit_identical_to_ram_on_fixture() {
    let (spectra, idx, params) = small_search_fixture();
    let scorer = make_scorer(0.05);

    let (ram_q, ram_c) = run_prepared(&idx, &params, &spectra, CandidateBacking::Ram, &scorer);
    let (mmap_q, mmap_c) = run_prepared(&idx, &params, &spectra, CandidateBacking::Mmap, &scorer);

    let ram_sig = psm_signature(&ram_q, &ram_c);
    let mmap_sig = psm_signature(&mmap_q, &mmap_c);

    // Sanity: the fixture must actually produce PSMs (else the test is vacuous).
    assert!(
        ram_sig.iter().any(|scan| !scan.is_empty()),
        "fixture produced no PSMs on the Ram path — test would be vacuous"
    );

    assert_eq!(
        ram_sig, mmap_sig,
        "mmap candidate path must yield identical PSMs to the in-RAM path\nRAM : {ram_sig:#?}\nMMAP: {mmap_sig:#?}"
    );
}
