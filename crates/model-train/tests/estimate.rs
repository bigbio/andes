//! TDD tests for [`model_train::estimate::Estimator`].
//!
//! The three tests correspond directly to the specification:
//! 1. An estimated `Param` from real counts is scorable (finite node scores, no panic).
//! 2. A partition with counts below `min_count` backs off rather than producing
//!    all-zero / -inf entries.
//! 3. Dense counts recover the empirical rank-1 fraction within tolerance.

use std::path::Path;

use model::amino_acid::AminoAcid;
use model::mass::nominal_from;
use model::peptide::Peptide;
use model::spectrum::Spectrum;
use scoring_crate::param_model::{IonType, Param, Partition};
use scoring_crate::scoring::rank_scorer::RankScorer;

use model_train::accumulate::StatsAccumulator;
use model_train::counts::CountStats;
use model_train::estimate::{Estimator, EstimatorConfig};

use rustc_hash::FxHashMap;
use scoring_crate::param_model::{FragmentOffsetFrequency, SpecDataType};
use model::activation::ActivationMethod;
use model::instrument::InstrumentType;
use model::protocol::Protocol;
use model::tolerance::Tolerance;

/// Minimal single-partition template with a configurable `max_rank` and one
/// prefix ion, used to probe rank-distribution smoothing in isolation.
fn one_partition_template(max_rank: i32) -> Param {
    let part = Partition { charge: 2, parent_mass: 1000.0, seg_num: 0 };
    let prefix1 = IonType::Prefix { charge: 1, offset_bits: 0.0_f32.to_bits() };
    let mut frag_off_table = FxHashMap::default();
    frag_off_table.insert(
        part,
        vec![FragmentOffsetFrequency { ion_type: prefix1, frequency: 0.7 }],
    );
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
        max_rank,
        rank_dist_table: FxHashMap::default(),
        error_scaling_factor: 0,
        ion_err_dist_table: FxHashMap::default(),
        noise_err_dist_table: FxHashMap::default(),
        ion_existence_table: FxHashMap::default(),
        partition_ion_types_cache: FxHashMap::default(),
    };
    p.rebuild_cache();
    p
}

/// Estimator-dilution regression: a strongly-peaked empirical NOISE rank
/// distribution must NOT be flattened by Laplace smoothing. With `max_rank=150`
/// (151 slots) the legacy add-1 over noise collapsed the peak to ~0.67, which
/// inflated `noise_freq[r]` at signal ranks and compressed `ln(ion/noise)`
/// node scores (the diagnosed −4.3% dilution). The noise model must stay sharp.
#[test]
fn noise_rank_dist_stays_sharp_not_flattened_by_smoothing() {
    let max_rank = 150;
    let part = Partition { charge: 2, parent_mass: 1000.0, seg_num: 0 };
    let prefix1 = IonType::Prefix { charge: 1, offset_bits: 0.0_f32.to_bits() };
    let template = one_partition_template(max_rank);

    let mut counts = CountStats::new();
    for _ in 0..300 {
        counts.bump_rank(part, prefix1, 0);
        counts.bump_rank(part, IonType::Noise, 0); // all noise mass on slot 0
    }

    let param = Estimator::new(EstimatorConfig::default()).estimate(&counts, &template);
    let noise = &param.rank_dist_table[&part][&IonType::Noise];
    let peak = noise.iter().cloned().fold(0.0_f32, f32::max);
    assert!(peak > 0.9, "noise peak should stay sharp (>0.9), got {peak}");
}

// ---------------------------------------------------------------------------
// Shared fixtures
// ---------------------------------------------------------------------------

fn make_peptide(seq: &[u8]) -> Peptide {
    let residues: Vec<AminoAcid> = seq
        .iter()
        .map(|&r| AminoAcid::standard(r).unwrap())
        .collect();
    Peptide::new(residues, b'_', b'-')
}

fn load_hcd_scorer() -> (Param, RankScorer) {
    let param_path = Path::new(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/HCD_QExactive_Tryp.param"
    ));
    let param = Param::load_from_file(param_path).expect("load HCD_QExactive_Tryp.param");
    let scorer = RankScorer::new(&param);
    (param, scorer)
}

/// Build a synthetic spectrum with peaks at every theoretical ion position for
/// `peptide`, with strictly decreasing intensities so that split-1 ions get
/// rank 1.
fn build_synthetic_spectrum(
    peptide: &Peptide,
    scorer: &RankScorer,
    precursor_mz: f64,
    charge: u8,
) -> Spectrum {
    let param = scorer.param();
    let n = peptide.length();
    assert!(n >= 2);

    let peptide_nominal = peptide.nominal_residue_mass();
    let mut prefix_acc = 0.0_f64;
    let mut prefix_nominals: Vec<i32> = vec![0];
    for s in 1..n {
        let aa = &peptide.residues[s - 1];
        prefix_acc += aa.mass + aa.mod_.as_ref().map_or(0.0, |m| m.mass_delta);
        prefix_nominals.push(nominal_from(prefix_acc));
    }

    const PROTON: f64 = 1.007_276_49;
    let parent_mass = (precursor_mz - PROTON) * (charge as f64);

    let num_segs = param.num_segments as usize;
    let mut all_ions: Vec<(f64, f32)> = Vec::new();
    let mut global_index = 0usize;

    #[allow(clippy::needless_range_loop)] // `split` indexes parallel arrays
    for split in 1..n {
        let prefix_nom = prefix_nominals[split] as f64;
        let suffix_nom = (peptide_nominal - prefix_nominals[split]) as f64;
        for (is_prefix, nominal_mass) in [(true, prefix_nom), (false, suffix_nom)] {
            for seg in 0..num_segs {
                let part = param.partition_for(charge, parent_mass, seg);
                let ions = param.ion_types_for_partition_slice(charge, parent_mass, seg);
                for &ion in ions {
                    let theo_mz = match (is_prefix, ion) {
                        (true, IonType::Prefix { .. }) => ion.mz(nominal_mass),
                        (false, IonType::Suffix { .. }) => ion.mz(nominal_mass),
                        _ => continue,
                    };
                    if param.segment_num(theo_mz, parent_mass) != seg { continue; }
                    if theo_mz <= 0.0 { continue; }
                    let _ = part; // partition is used implicitly via param lookups
                    let intensity = 100_000.0_f32 / (1.0 + global_index as f32);
                    all_ions.push((theo_mz, intensity));
                    global_index += 1;
                }
            }
        }
    }

    let mut peaks: Vec<(f64, f32)> = all_ions;
    peaks.push((3000.1, 0.001));
    peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    Spectrum {
        title: "synthetic".into(),
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

/// Accumulate N identical PSMs into a `CountStats` using the real HCD scorer.
fn accumulate_n(n: usize) -> (CountStats, Param) {
    let (template, scorer) = load_hcd_scorer();
    let peptide = make_peptide(b"PEPTIDE");
    let charge: u8 = 2;
    const PROTON: f64 = 1.007_276_49;
    let precursor_mz = (peptide.mass() + charge as f64 * PROTON) / charge as f64;
    let spectrum = build_synthetic_spectrum(&peptide, &scorer, precursor_mz, charge);

    let acc = StatsAccumulator::new(&scorer);
    let mut stats = CountStats::new();
    for _ in 0..n {
        acc.accumulate(&mut stats, &spectrum, &peptide, charge);
    }
    (stats, template)
}

// ---------------------------------------------------------------------------
// Test 1: estimated Param is scorable and produces finite node scores
// ---------------------------------------------------------------------------

#[test]
fn estimated_param_is_scorable_and_finite() {
    // Accumulate 200 identical synthetic PSMs so every partition sees counts.
    let (counts, template) = accumulate_n(200);

    let estimator = Estimator::new(EstimatorConfig::default());
    let param = estimator.estimate(&counts, &template);

    // RankScorer::new must not panic.
    let scorer = RankScorer::new(&param);

    // For every populated partition × ion, node_score at rank 1 must be finite
    // (not NaN or +/-inf).  The binary template has a rich ion set; the estimated
    // param should cover at least the partitions we observed.
    let mut checked = 0usize;
    for (&part, ion_table) in &param.rank_dist_table {
        for &ion in ion_table.keys() {
            if ion.is_noise() { continue; }
            let s = scorer.node_score(part, ion, 1);
            assert!(
                s.is_finite(),
                "node_score({:?}, {:?}, rank=1) = {s} is not finite",
                part, ion
            );
            let m = scorer.missing_ion_score(part, ion);
            assert!(
                m.is_finite(),
                "missing_ion_score({:?}, {:?}) = {m} is not finite",
                part, ion
            );
            checked += 1;
        }
    }
    assert!(checked > 0, "no (partition, ion) pairs were checked — rank_dist_table is empty");
}

// ---------------------------------------------------------------------------
// Test 2: empty (below min_count) partition backs off to a non-degenerate
//         distribution (never all-zero, never -inf log)
// ---------------------------------------------------------------------------

#[test]
fn empty_partition_backs_off_not_zero() {
    // Use an *empty* CountStats so every partition is below min_count.
    let (template, _scorer) = load_hcd_scorer();
    let empty_counts = CountStats::new();

    let estimator = Estimator::new(EstimatorConfig::default());
    let param = estimator.estimate(&empty_counts, &template);

    // Even with zero observed counts, every rank_dist_table entry must be
    // populated (Laplace + global-pool backoff) and the resulting log scores
    // must be finite.
    let scorer = RankScorer::new(&param);

    let mut checked = 0usize;
    for (&part, ion_table) in &param.rank_dist_table {
        // Noise entry must be present.
        assert!(
            ion_table.contains_key(&IonType::Noise),
            "partition {:?} missing Noise entry in rank_dist_table",
            part
        );
        for (&ion, freqs) in ion_table {
            // Every frequency must be positive (no zeros — Laplace guarantees this).
            for (i, &f) in freqs.iter().enumerate() {
                assert!(
                    f > 0.0,
                    "freq[{i}] = {f} <= 0 for partition {:?} ion {:?}",
                    part, ion
                );
            }
            if ion.is_noise() { continue; }
            let s = scorer.node_score(part, ion, 1);
            assert!(
                s.is_finite(),
                "node_score for {:?} {:?} = {s} is not finite after empty-count backoff",
                part, ion
            );
            checked += 1;
        }
    }
    assert!(checked > 0, "no entries checked — rank_dist_table unexpectedly empty");
}

// ---------------------------------------------------------------------------
// Test 3: dense counts recover empirical rank-1 fraction within tolerance
// ---------------------------------------------------------------------------

#[test]
fn dense_counts_recover_empirical_within_tolerance() {
    // With a very small pseudo (effectively zero relative to many counts) and
    // many PSMs, the normalised rank-1 frequency for the most-observed
    // (partition, ion) should be close to the empirical fraction.

    let (template, scorer_seed) = load_hcd_scorer();
    let peptide = make_peptide(b"PEPTIDE");
    let charge: u8 = 2;
    const PROTON: f64 = 1.007_276_49;
    let precursor_mz = (peptide.mass() + charge as f64 * PROTON) / charge as f64;
    let spectrum = build_synthetic_spectrum(&peptide, &scorer_seed, precursor_mz, charge);

    // 2000 PSMs: Laplace pseudo (1.0) is negligible relative to counts.
    let acc = StatsAccumulator::new(&scorer_seed);
    let mut counts = CountStats::new();
    for _ in 0..2000 {
        acc.accumulate(&mut counts, &spectrum, &peptide, charge);
    }

    // Estimator with tiny pseudo to minimise smoothing bias.
    let cfg = EstimatorConfig {
        pseudo: 0.001,
        min_count: 1,    // disable backoff (all partitions will exceed 1 count)
        backoff_weight: 0.0,
        ..Default::default()
    };
    let estimator = Estimator::new(cfg);
    let param = estimator.estimate(&counts, &template);

    // Find the (partition, ion) with the most counts at rank index 0.
    let mut best: Option<(Partition, IonType, u64, u64)> = None; // (part, ion, rank0, total)
    for (&(part, ion), v) in &counts.rank {
        if ion.is_noise() { continue; }
        let rank0 = v.first().copied().unwrap_or(0);
        let total: u64 = v.iter().sum();
        if rank0 > 0 && total > 0 && best.as_ref().is_none_or(|&(_, _, b_r, _)| rank0 > b_r) {
            best = Some((part, ion, rank0, total));
        }
    }
    let (best_part, best_ion, rank0_count, total_count) =
        best.expect("expected at least one (partition, ion) with rank-0 counts");

    // Empirical fraction for rank-1 (index 0).
    let empirical = rank0_count as f32 / total_count as f32;

    // Estimated frequency at index 0 from the rank_dist_table.
    let estimated = param
        .rank_dist_table
        .get(&best_part)
        .and_then(|t| t.get(&best_ion))
        .and_then(|v| v.first().copied())
        .expect("rank_dist_table entry missing for best (partition, ion)");

    // With pseudo=0.001 and n=2000, the smoothed value should be within 1%
    // of the empirical fraction.
    let tol = 0.01_f32;
    assert!(
        (estimated - empirical).abs() <= tol,
        "estimated rank-1 freq {estimated:.4} deviates from empirical {empirical:.4} by more than {tol} \
         for partition {:?} ion {:?}",
        best_part, best_ion
    );
}

/// A sparse partition (n < min_count) must shrink toward the INDEPENDENT PRIOR's
/// distribution, not the corpus-internal pool. Here the corpus empirical mass is
/// all on slot 0, but the prior is peaked on slot 5; the blended result must
/// carry materially more mass on slot 5 than the no-prior estimate does.
#[test]
fn sparse_partition_shrinks_toward_independent_prior() {
    let max_rank = 150;
    let n_slots = (max_rank + 1) as usize;
    let part = Partition { charge: 2, parent_mass: 1000.0, seg_num: 0 };
    let prefix1 = IonType::Prefix { charge: 1, offset_bits: 0.0_f32.to_bits() };
    let template = one_partition_template(max_rank);

    let mut prior = one_partition_template(max_rank);
    let mut prior_dist = vec![0.0_f32; n_slots];
    prior_dist[5] = 1.0;
    let mut ion_map = FxHashMap::default();
    ion_map.insert(prefix1, prior_dist);
    let mut noise_dist = vec![0.0_f32; n_slots];
    noise_dist[0] = 1.0;
    ion_map.insert(IonType::Noise, noise_dist);
    prior.rank_dist_table.insert(part, ion_map);

    let mut counts = CountStats::new();
    for _ in 0..10 {
        counts.bump_rank(part, prefix1, 0);
        counts.bump_rank(part, IonType::Noise, 0);
    }

    let est = Estimator::new(EstimatorConfig::default());
    let with_prior = est.estimate_with_prior(&counts, &template, Some(&prior));
    let no_prior = est.estimate_with_prior(&counts, &template, None);

    let p5_with = with_prior.rank_dist_table[&part][&prefix1][5];
    let p5_without = no_prior.rank_dist_table[&part][&prefix1][5];
    assert!(
        p5_with > p5_without + 0.05,
        "prior must pull mass toward slot 5: with={p5_with} without={p5_without}"
    );
}
