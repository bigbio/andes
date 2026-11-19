//! The calibrated strong score (PIN `RawScoreCal`) of a spectrum without a recorded charge is a
//! leave-one-out z-score over the candidates of every charge the spectrum is searched with.
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
use scoring_crate::scoring::strong_score_calibrated_loo;
use scoring_crate::{Param, RankScorer};
use search::{PreparedSearch, SearchIndex, SearchParams};

const FRAG_TOL_DA: f64 = 0.5;

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

#[test]
fn a_charge_less_spectrum_calibrates_its_strong_score_over_every_searched_charge() {
    let (mut spectra, idx) = fixture_inputs();
    for s in &mut spectra {
        s.precursor_charge = None;
    }
    let scorer = scorer();
    let mut params = SearchParams::default_tryptic(aa_set());
    params.precursor_tolerance =
        model::tolerance::PrecursorTolerance::symmetric(Tolerance::Ppm(20.0));
    params.charge_range = 1..=4;
    params.fragment_index_top_k = search::search_params::FRAGMENT_INDEX_TOP_K;
    params.top_n_psms_per_spectrum = 1000;

    let prepared = PreparedSearch::prepare_in_memory(&idx, &params, &scorer, FRAG_TOL_DA, "XXX")
        .expect("build the in-memory candidate index");
    let queues = prepared.run_chunk(&spectra, 0);

    let mut checked = 0usize;
    for q in &queues {
        let psms = q.clone().into_rank_sorted_vec();
        let charges: std::collections::BTreeSet<u8> = psms.iter().map(|p| p.charge_used).collect();
        if charges.len() < 2 {
            continue;
        }
        let pool: Vec<f32> = psms.iter().map(|p| p.features.strong_score).collect();
        for p in &psms {
            let want = strong_score_calibrated_loo(&pool, p.features.strong_score);
            let got = p.features.strong_score_cal;
            assert!(
                (got - want).abs() <= 1e-3 * want.abs().max(1.0),
                "charge {}: RawScoreCal {got} against {want} over the {} candidates of charges {charges:?}",
                p.charge_used,
                pool.len()
            );
        }
        checked += 1;
    }
    assert!(
        checked >= 10,
        "only {checked} spectra matched at two or more charges"
    );
}
