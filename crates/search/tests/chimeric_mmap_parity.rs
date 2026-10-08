//! The chimeric Pass 2 (`run_pass2_coisolation`) must emit the SAME secondary
//! PSMs on the out-of-core (`Mmap`) candidate backing as on the in-RAM one:
//! same peptidoform, protein, decoy flag, scores and features.
//!
//! Each spectrum carries a primary peptide (its precursor) plus the fragment
//! ions of a co-isolated peptide whose averagine envelope sits in a synthetic
//! linked MS1 scan inside the isolation window. Two cases:
//!
//!  * `DAFLGSFLYEYSR` — an ordinary target secondary;
//!  * protein-N-term-acetylated `VDLLAGEK` — a COLLISION DECOY: it is the
//!    N-terminal peptide of a reversed decoy protein, its bare sequence is also
//!    an internal target peptide, and the target twin has no acetylated form, so
//!    the twin lies 42 Da outside the secondary's window. Both backings must
//!    relabel it as target (the global relabel), not only when the twin happens
//!    to be among the window's candidates.

use rustc_hash::FxHashMap;

use input::Ms1Link;
use model::isotope::averagine_isotope_envelope;
use model::mass::ISOTOPE;
use model::{
    activation::ActivationMethod, instrument::InstrumentType, protocol::Protocol, AminoAcid,
    AminoAcidSetBuilder, ModLocation, Modification, Peptide, Protein, ProteinDb, ResidueSpec,
    Spectrum, Tolerance, H2O, PROTON,
};
use scoring_crate::param_model::{FragmentOffsetFrequency, IonType, Partition, SpecDataType};
use scoring_crate::scoring::fragment_ions::predict_by_ions;
use scoring_crate::{Param, RankScorer};
use search::candidate_gen::Candidate;
use search::psm::TopNQueue;
use search::{run_pass2_coisolation, PreparedSearch, SearchIndex, SearchParams};

const FRAG_TOL_DA: f64 = 0.05;

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

fn residues(bytes: &[u8]) -> Vec<AminoAcid> {
    bytes
        .iter()
        .map(|&b| AminoAcid::standard(b).unwrap())
        .collect()
}

fn acetyl() -> Modification {
    Modification {
        name: "Acetyl".to_string(),
        mass_delta: 42.010565,
        residue: ResidueSpec::Wildcard,
        location: ModLocation::ProtNTerm,
        fixed: false,
        accession: None,
        neutral_losses: Vec::new(),
        loss_class: 0,
    }
}

/// Target proteins. `P2` reversed starts with `VDLLAGEK|R...`, so the decoy
/// protein's N-terminal peptide is `VDLLAGEK`, which `P1` also contains
/// internally (never at a protein N-terminus).
fn fixture() -> (SearchIndex, SearchParams) {
    let target = ProteinDb {
        proteins: vec![
            Protein {
                accession: "P1".into(),
                description: "".into(),
                sequence: b"MKWVTFISLLRKDAFLGSFLYEYSRVDLLAGEKGR".to_vec(),
            },
            Protein {
                accession: "P2".into(),
                description: "".into(),
                sequence: b"MPSTRKEGALLDV".to_vec(),
            },
        ],
    };
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let aa_set = AminoAcidSetBuilder::new_standard()
        .add_variable_mod(acetyl())
        .build()
        .unwrap();
    let mut params = SearchParams::default_tryptic(aa_set);
    params.min_length = 4;
    params.max_variable_mods_per_peptide = 1;
    params.min_peaks = 0;
    params.chimeric = true;
    params.top_n_psms_per_spectrum = 1;
    (idx, params)
}

/// Averagine envelope (4 peaks) of a `charge` precursor of neutral `mass`.
fn envelope(mass: f64, charge: u8, scale: f32) -> Vec<(f64, f32)> {
    let mono_mz = (mass + charge as f64 * PROTON) / charge as f64;
    averagine_isotope_envelope(mass, 4)
        .iter()
        .enumerate()
        .map(|(k, &p)| {
            (
                mono_mz + k as f64 * ISOTOPE / charge as f64,
                p as f32 * scale,
            )
        })
        .collect()
}

/// A spectrum selected on `primary` (z=2) that also carries `secondary`'s b/y
/// ions, and its linked MS1 scan holding both precursors' envelopes. The
/// isolation window is widened so the z=3 secondary falls inside it.
fn chimeric_scan(
    primary: &Peptide,
    secondary: &Peptide,
    title: &str,
) -> (Spectrum, Vec<(f64, f32)>) {
    let mut peaks: Vec<(f64, f32)> = Vec::new();
    for (i, p) in predict_by_ions(primary, 1..=1).iter().enumerate() {
        peaks.push((p.mz, 200.0 - i as f32));
    }
    for (i, p) in predict_by_ions(secondary, 1..=1).iter().enumerate() {
        peaks.push((p.mz, 100.0 - i as f32));
    }
    peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
    let spec = Spectrum {
        title: title.into(),
        precursor_mz: (primary.mass() + 2.0 * PROTON) / 2.0,
        precursor_intensity: None,
        precursor_charge: Some(2),
        rt_seconds: None,
        scan: None,
        peaks,
        activation_method: None,
        isolation_lower_offset: Some(300.0),
        isolation_upper_offset: Some(300.0),
    };
    let mut ms1 = envelope(primary.mass(), 2, 1000.0);
    ms1.extend(envelope(secondary.mass(), 3, 500.0));
    ms1.sort_by(|a, b| a.0.total_cmp(&b.0));
    (spec, ms1)
}

fn spectra_and_link() -> (Vec<Spectrum>, Ms1Link) {
    let primary = Peptide::new(residues(b"WVTFISLLR"), b'K', b'K');
    let ordinary = Peptide::new(residues(b"DAFLGSFLYEYSR"), b'K', b'V');
    let mut ac = residues(b"VDLLAGEK");
    ac[0].mod_ = Some(model::modification::leak_mod(acetyl()));
    let collision = Peptide::new(ac, b'-', b'R');
    let (s0, m0) = chimeric_scan(&primary, &ordinary, "scan=ordinary");
    let (s1, m1) = chimeric_scan(&primary, &collision, "scan=collision");
    let link = Ms1Link {
        ms1_peaks: vec![m0, m1],
        ms2_to_ms1: vec![Some(0), Some(1)],
    };
    (vec![s0, s1], link)
}

/// One emitted row: (secondary?, peptidoform, protein, decoy, score, rank_score,
/// charge, mass error, features).
type Row = (bool, String, usize, bool, u32, u32, u8, u64, String);

fn rows(queues: &[TopNQueue], candidates: &[Candidate]) -> Vec<Vec<Row>> {
    queues
        .iter()
        .map(|q| {
            let mut out: Vec<Row> = q
                .iter_psms()
                .map(|psm| {
                    let c = &candidates[psm.primary_candidate_idx() as usize];
                    let pep: String = c
                        .peptide
                        .residues
                        .iter()
                        .map(|a| match &a.mod_ {
                            Some(m) => format!("{}[{:.4}]", a.residue as char, m.mass_delta),
                            None => (a.residue as char).to_string(),
                        })
                        .collect();
                    (
                        psm.precursor_mz_override.is_some(),
                        pep,
                        c.protein_index,
                        c.is_decoy,
                        psm.score.to_bits(),
                        psm.rank_score.to_bits(),
                        psm.charge_used,
                        psm.mass_error_ppm.to_bits(),
                        format!("{:?}", psm.features),
                    )
                })
                .collect();
            out.sort_by(|a, b| format!("{a:?}").cmp(&format!("{b:?}")));
            out
        })
        .collect()
}

fn run(mmap: bool, fragment_index: bool) -> Vec<Vec<Row>> {
    let (idx, mut params) = fixture();
    if fragment_index {
        params.fragment_index_top_k = search::search_params::FRAGMENT_INDEX_TOP_K;
    }
    let scorer = make_scorer(FRAG_TOL_DA);
    let (spectra, link) = spectra_and_link();
    let tmp = tempfile::NamedTempFile::new().unwrap();
    let mut prepared = if mmap {
        PreparedSearch::prepare_mmap(&idx, &params, &scorer, FRAG_TOL_DA, "XXX", tmp.path())
            .expect("prepare_mmap accepts --chimeric")
    } else {
        PreparedSearch::prepare(&idx, &params, &scorer, FRAG_TOL_DA, "XXX")
    };
    let mut queues = prepared.run_chunk(&spectra, 0);
    run_pass2_coisolation(&prepared, &spectra, &mut queues, &params, &link, 0);
    prepared.sync_materialized_candidates();
    rows(&queues, &prepared.candidates)
}

#[test]
fn mmap_pass2_secondaries_match_ram() {
    let ram = run(false, false);
    let mmap = run(true, false);

    // Non-vacuous: every scan has its primary AND one secondary on RAM.
    for (scan, r) in ram.iter().enumerate() {
        assert_eq!(
            r.iter().filter(|row| row.0).count(),
            1,
            "scan {scan}: expected one secondary on RAM, got {r:#?}"
        );
        assert_eq!(r.iter().filter(|row| !row.0).count(), 1);
    }
    let sec = |rs: &Vec<Row>| rs.iter().find(|r| r.0).cloned().unwrap();
    assert_eq!(sec(&ram[0]).1, "DAFLGSFLYEYSR");
    let coll = sec(&ram[1]);
    assert_eq!(coll.1, "V[42.0106]DLLAGEK", "collision secondary: {coll:?}");
    assert!(
        !coll.3,
        "the acetylated decoy N-terminal peptide collides with a target sequence and \
         must be relabeled as target"
    );

    assert_eq!(ram, mmap, "Pass-2 rows differ between RAM and mmap");
}

/// Fragment-index retrieval changes only Pass 1's candidate shortlist; Pass 2
/// still enumerates each secondary window, so on this fixture (where Pass 1
/// finds the same primary) its rows equal the RAM rows.
#[test]
fn mmap_fragment_index_pass2_secondaries_match_ram() {
    let ram = run(false, false);
    let mmap_index = run(true, true);
    let secondaries = |rs: &[Vec<Row>]| -> Vec<Vec<Row>> {
        rs.iter()
            .map(|r| r.iter().filter(|row| row.0).cloned().collect())
            .collect()
    };
    assert_eq!(secondaries(&ram), secondaries(&mmap_index));
}
