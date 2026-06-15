//! SP4 mechanism proof (synthetic): the full neutral-loss train→score round-trip.
//!
//! With no real glyco corpus available, this exercises the machinery end-to-end
//! on a synthetic spectrum whose fragments spanning a loss-bearing residue carry
//! loss-shifted peaks:
//!
//!   1. SP3: accumulate confident glyco PSMs → the estimator learns a per-class
//!      neutral-loss rank table (a `loss_class != 0` IonType key).
//!   2. Task 7: scoring a glyco peptide against that model adds positive
//!      loss-ion node score, so its RawScore exceeds the same peptide scored
//!      with its losses NOT declared (the loss peaks then go unprobed).
//!
//! The mod's `mass_delta` is 0 here so the with/without-losses peptides share
//! identical residue masses (hence identical intact ions) — this isolates the
//! loss machinery, which is what the test validates. The real biochemistry
//! (Unimod 393 Hex/Hex2 on a glyco residue) is exercised by the real-corpus
//! benchmark, not this mechanism test.

use std::sync::Arc;

use model::amino_acid::AminoAcid;
use model::mass::nominal_from;
use model::modification::{ModLocation, Modification, ResidueSpec};
use model::peptide::Peptide;
use model::spectrum::Spectrum;
use scoring_crate::param_model::{
    FragmentOffsetFrequency, IonType, Param, Partition, SpecDataType,
};
use scoring_crate::scoring::psm_score::score_psm;
use scoring_crate::scoring::rank_scorer::RankScorer;
use scoring_crate::scoring::scored_spectrum::ScoredSpectrum;

use model::activation::ActivationMethod;
use model::instrument::InstrumentType;
use model::protocol::Protocol;
use model::tolerance::Tolerance;

use model_train::accumulate::StatsAccumulator;
use model_train::counts::CountStats;
use model_train::estimate::{Estimator, EstimatorConfig};
use model_train::store::{write_models, ModelStore};

use rustc_hash::FxHashMap;

const HEX_LOSS: f64 = 162.0528; // Unimod 393 -Hex
const PROTON: f64 = 1.007_276_49;

/// CID / low-res seed template with a single intact b-ion type (charge 1,
/// proton offset). `parent_mass = 0.0` so the partition floor-matches any
/// peptide. CID predicts neutral losses, so the training matcher emits loss
/// facts for loss-bearing peptides.
fn seed_template() -> Param {
    let part = Partition { charge: 2, parent_mass: 0.0, seg_num: 0 };
    let bion = IonType::Prefix { charge: 1, offset_bits: (PROTON as f32).to_bits(), loss_class: 0 };
    let mut frag_off_table = FxHashMap::default();
    frag_off_table.insert(
        part,
        vec![FragmentOffsetFrequency { ion_type: bion, frequency: 0.9 }],
    );
    // Seed rank tables for the b-ion + Noise so the seed scorer's
    // partition_ion_logs is non-empty (the training matcher iterates that intact
    // vocabulary to derive loss variants). Values are placeholders — the real
    // tables are re-estimated from the accumulated counts.
    let max_rank = 10usize;
    let n_slots = max_rank + 1;
    let bion_freqs: Vec<f32> = (0..n_slots).map(|r| 0.5_f32 * 0.7_f32.powi(r as i32)).collect();
    let noise_freqs: Vec<f32> = vec![1.0_f32 / n_slots as f32; n_slots];
    let mut ion_table: FxHashMap<IonType, Vec<f32>> = FxHashMap::default();
    ion_table.insert(bion, bion_freqs);
    ion_table.insert(IonType::Noise, noise_freqs);
    let mut rank_dist_table: FxHashMap<Partition, FxHashMap<IonType, Vec<f32>>> = FxHashMap::default();
    rank_dist_table.insert(part, ion_table);
    let mut p = Param {
        version: 10001,
        data_type: SpecDataType {
            activation: ActivationMethod::CID,
            instrument: InstrumentType::LowRes,
            enzyme: None,
            protocol: Protocol::Automatic,
        },
        mme: Tolerance::Da(0.5),
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
        max_rank: max_rank as i32,
        rank_dist_table,
        error_scaling_factor: 0,
        ion_err_dist_table: FxHashMap::default(),
        noise_err_dist_table: FxHashMap::default(),
        ion_existence_table: FxHashMap::default(),
        partition_ion_types_cache: FxHashMap::default(),
    };
    p.rebuild_cache();
    p
}

/// `AAAT*AAA` with a Hex neutral-loss (class 1, mass_delta 0) on the T at index
/// 3. If `declare_losses` is false the same residues are built without the loss
/// attributes (identical masses, no loss ions).
fn glyco_peptide(declare_losses: bool) -> Peptide {
    let m = Modification {
        name: "Glyco".into(),
        mass_delta: 0.0,
        residue: ResidueSpec::Specific(b'T'),
        location: ModLocation::Anywhere,
        fixed: false,
        accession: Some("UNIMOD:393".into()),
        neutral_losses: if declare_losses { vec![HEX_LOSS] } else { vec![] },
        loss_class: if declare_losses { 1 } else { 0 },
    };
    let arc = Arc::new(m);
    let residues: Vec<AminoAcid> = b"AAATAAA"
        .iter()
        .enumerate()
        .map(|(i, &r)| {
            let aa = AminoAcid::standard(r).unwrap();
            if i == 3 { aa.with_mod(arc.clone()) } else { aa }
        })
        .collect();
    Peptide::new(residues, b'_', b'-')
}

/// Synthetic spectrum: a b-ion peak at every prefix split, plus a loss-shifted
/// (-Hex) peak for every prefix that spans the loss residue (index 3). Peaks
/// are intensity-ordered so split-1 ions rank high.
fn glyco_spectrum(peptide: &Peptide, charge: u8) -> Spectrum {
    let bion = IonType::Prefix { charge: 1, offset_bits: (PROTON as f32).to_bits(), loss_class: 0 };
    let n = peptide.length();
    let loss_index = 3usize;
    let mut peaks: Vec<(f64, f32)> = Vec::new();
    let mut prefix_acc = 0.0_f64;
    let mut intensity = 100_000.0_f32;
    for s in 1..n {
        prefix_acc += peptide.residues[s - 1].mass; // mass_delta is 0
        let prefix_nominal = nominal_from(prefix_acc);
        let intact_mz = bion.mz(prefix_nominal as f64);
        peaks.push((intact_mz, intensity));
        intensity *= 0.95;
        // Prefix spans the loss residue once we've passed index 3.
        if loss_index < s {
            let loss_mz = intact_mz - HEX_LOSS;
            if loss_mz > 0.0 {
                peaks.push((loss_mz, intensity));
                intensity *= 0.95;
            }
        }
    }
    peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    let precursor_mz = (peptide.mass() + charge as f64 * PROTON) / charge as f64;
    Spectrum {
        title: "synthetic-glyco".into(),
        precursor_mz,
        precursor_intensity: None,
        precursor_charge: Some(charge as i32),
        rt_seconds: None,
        scan: None,
        peaks,
        activation_method: None,
        isolation_lower_offset: None,
        isolation_upper_offset: None,
    }
}

#[test]
fn glyco_loss_table_is_learned_and_lifts_the_score() {
    let charge: u8 = 2;
    let template = seed_template();
    let seed_scorer = RankScorer::new(&template);
    let pep = glyco_peptide(true);
    let spec = glyco_spectrum(&pep, charge);

    // ── SP3: accumulate confident glyco PSMs and estimate a model ────────────
    let acc = StatsAccumulator::new(&seed_scorer);
    let mut counts = CountStats::new();
    for _ in 0..200 {
        acc.accumulate(&mut counts, &spec, &pep, charge);
    }
    let trained = Estimator::new(EstimatorConfig::default()).estimate(&counts, &template);

    // The estimator must have created a per-class loss rank table.
    let learned_loss_key = trained
        .rank_dist_table
        .values()
        .any(|ion_map| ion_map.keys().any(|ion| ion.loss_class() == 1));
    assert!(
        learned_loss_key,
        "SP3: estimator must learn a loss_class=1 rank table from glyco PSMs"
    );

    let trained_scorer = RankScorer::new(&trained);
    assert!(
        trained_scorer.has_loss_tables(),
        "trained model must expose loss tables to the scorer"
    );

    // ── Task 7: the loss ions lift the RawScore ──────────────────────────────
    // Same spectrum + same trained scorer; the only difference is whether the
    // peptide declares its losses (and is thus probed for loss ions).
    let ss = ScoredSpectrum::new(&spec, &trained_scorer, charge);
    let with_losses = score_psm(&ss, &glyco_peptide(true), &trained_scorer, charge, 0.5);
    let without_losses = score_psm(&ss, &glyco_peptide(false), &trained_scorer, charge, 0.5);

    assert!(
        with_losses > without_losses,
        "loss ions must lift the RawScore: with={with_losses} without={without_losses}"
    );

    // ── SP4: the learned loss table survives a parquet store round-trip ──────
    // (Closes the full pipeline: train → serialize → load → score.)
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("glyco_models.parquet");
    write_models(&path, &[("glyco".to_string(), &trained)]).unwrap();
    let store = ModelStore::open(&path).unwrap();
    let reloaded = store.load_param("glyco").unwrap();
    let reloaded_scorer = RankScorer::new(&reloaded);
    assert!(
        reloaded_scorer.has_loss_tables(),
        "loss tables must survive the parquet round-trip"
    );
    let ss2 = ScoredSpectrum::new(&spec, &reloaded_scorer, charge);
    let with_losses_reloaded =
        score_psm(&ss2, &glyco_peptide(true), &reloaded_scorer, charge, 0.5);
    assert_eq!(
        with_losses_reloaded, with_losses,
        "reloaded model must score the glyco peptide identically to the in-memory model"
    );
}
