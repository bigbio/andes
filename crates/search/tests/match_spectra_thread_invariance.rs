//! Thread-count invariance: a search returns the same PSMs, in the same order,
//! whatever the Rayon thread count, with and without fragment-index retrieval.
//!
//! Each spectrum's pipeline runs on one worker, but PSMs reach a spectrum's queue
//! in candidate-enumeration order and the candidate store is filled in completion
//! order, so exact score ties used to be settled by push order. They are settled
//! by the `PsmMatch` total order now (`rank_score`, `score`, then peptide content),
//! and this test pins that down on real spectra whose scores tie often: the toy
//! scorer below gives integral scores.

mod common;
use common::*;

use std::fs::File;
use std::io::BufReader;

use rustc_hash::FxHashMap;

use input::{FastaReader, MgfReader};
use model::{
    activation::ActivationMethod, instrument::InstrumentType, protocol::Protocol, Spectrum,
    Tolerance, H2O, PROTON,
};
use scoring_crate::param_model::{FragmentOffsetFrequency, IonType, Partition, SpecDataType};
use scoring_crate::{Param, RankScorer};
use search::candidate_gen::Candidate;
use search::psm::TopNQueue;
use search::{PreparedSearch, SearchIndex, SearchParams};

const FRAG_TOL_DA: f64 = 0.5;

/// A b/y rank scorer with coarse rank probabilities, so many candidates tie.
fn scorer() -> RankScorer {
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
    let mut ion_table = FxHashMap::default();
    ion_table.insert(prefix1, vec![0.6_f32, 0.3, 0.05, 0.001]);
    ion_table.insert(suffix1, vec![0.6_f32, 0.3, 0.05, 0.001]);
    ion_table.insert(IonType::Noise, vec![0.1_f32, 0.2, 0.3, 0.4]);
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
            activation: ActivationMethod::CID,
            instrument: InstrumentType::LowRes,
            enzyme: None,
            protocol: Protocol::Automatic,
        },
        mme: Tolerance::Da(FRAG_TOL_DA),
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

fn fixture_inputs() -> (Vec<Spectrum>, SearchIndex) {
    let target = FastaReader::load_all(BufReader::new(
        File::open(fixture("test-fixtures/BSA.fasta")).unwrap(),
    ))
    .unwrap();
    let idx = SearchIndex::from_target_db(&target, "XXX");
    let reader = input::open_buf_maybe_gz(&fixture("test-fixtures/test.mgf.gz")).unwrap();
    let spectra: Vec<Spectrum> = MgfReader::new(reader).filter_map(|r| r.ok()).collect();
    assert!(
        spectra.len() > 1000,
        "fixture has {} spectra",
        spectra.len()
    );
    (spectra, idx)
}

/// One line per retained PSM, in output order, with everything the PIN writer
/// reads from it: peptide, decoy flag, proteins, charge, scores and errors.
fn signature(queues: &[TopNQueue], cands: &[Candidate]) -> Vec<String> {
    let mut out = Vec::new();
    for (i, q) in queues.iter().enumerate() {
        for (rank, p) in q.clone().into_rank_sorted_vec().iter().enumerate() {
            let c = &cands[p.primary_candidate_idx() as usize];
            let proteins: Vec<usize> = p
                .candidate_idxs
                .iter()
                .map(|&k| cands[k as usize].protein_index)
                .collect();
            out.push(format!(
                "{i}\t{rank}\t{}\t{}\t{proteins:?}\t{}\t{:08x}\t{:08x}\t{}\t{:016x}",
                c.peptide,
                c.is_decoy,
                p.charge_used,
                p.score.to_bits(),
                p.rank_score.to_bits(),
                p.isotope_offset,
                p.mass_error_ppm.to_bits(),
            ));
        }
    }
    out
}

fn run(
    threads: usize,
    spectra: &[Spectrum],
    idx: &SearchIndex,
    params: &SearchParams,
    scorer: &RankScorer,
) -> (Vec<String>, usize) {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("build pool");
    pool.install(|| {
        let mut prepared =
            PreparedSearch::prepare_in_memory(idx, params, scorer, FRAG_TOL_DA, "XXX")
                .expect("build the in-memory candidate index");
        let queues = prepared.run_chunk(spectra, 0);
        prepared.sync_materialized_candidates();
        let ties = queues
            .iter()
            .filter(|q| {
                let v = q.clone().into_rank_sorted_vec();
                v.len() > 1 && v[0].rank_score == v[1].rank_score
            })
            .count();
        (signature(&queues, &prepared.candidates), ties)
    })
}

fn assert_thread_invariant(top_k: u32) {
    let (spectra, idx) = fixture_inputs();
    let scorer = scorer();
    let mut params = SearchParams::default_tryptic(aa_set());
    params.precursor_tolerance =
        model::tolerance::PrecursorTolerance::symmetric(Tolerance::Ppm(20.0));
    params.charge_range = 2..=3;
    params.fragment_index_top_k = top_k;

    let (one, ties) = run(1, &spectra, &idx, &params, &scorer);
    assert!(!one.is_empty(), "no PSMs: fixture problem");
    assert!(
        ties > 0,
        "no spectrum has a top-1 tie; the fixture no longer exercises tie-breaking"
    );
    for threads in [4, 8] {
        let (many, _) = run(threads, &spectra, &idx, &params, &scorer);
        assert_eq!(
            one.len(),
            many.len(),
            "{threads} threads: PSM count differs"
        );
        if let Some((a, b)) = one.iter().zip(&many).find(|(a, b)| a != b) {
            panic!("{threads} threads differ from 1 thread:\n  1: {a}\n  {threads}: {b}");
        }
    }
    eprintln!(
        "top_k {top_k}: {} PSM rows identical across 1/4/8 threads; {ties} spectra with a top-1 tie",
        one.len()
    );
}

#[test]
fn match_spectra_output_invariant_across_thread_counts() {
    assert_thread_invariant(0);
}

#[test]
fn fragment_index_output_invariant_across_thread_counts() {
    assert_thread_invariant(search::search_params::FRAGMENT_INDEX_TOP_K);
}
