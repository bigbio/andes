//! Integration tests for [`model_train::accumulate::StatsAccumulator`].
//!
//! The TDD test builds a synthetic PSM where the spectrum peaks are placed
//! at the exact theoretical ion m/z values so we can assert deterministic
//! rank counts.

use std::path::Path;

use model::amino_acid::AminoAcid;
use model::mass::nominal_from;
use model::peptide::Peptide;
use model::spectrum::Spectrum;
use scoring_crate::param_model::{IonType, Param, Partition};
use scoring_crate::scoring::rank_scorer::RankScorer;

use model_train::accumulate::{merge, StatsAccumulator};
use model_train::counts::CountStats;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn make_peptide(seq: &[u8]) -> Peptide {
    let residues: Vec<AminoAcid> = seq
        .iter()
        .map(|&r| AminoAcid::standard(r).unwrap())
        .collect();
    Peptide::new(residues, b'_', b'-')
}

/// Load the bundled HCD_QExactive_Tryp param file and build a RankScorer.
fn load_hcd_scorer() -> (Param, RankScorer) {
    let param_path = Path::new(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/HCD_QExactive_Tryp.param"
    ));
    let param = Param::load_from_file(param_path).expect("load HCD_QExactive_Tryp.param");
    let scorer = RankScorer::new(&param);
    (param, scorer)
}

// ---------------------------------------------------------------------------
// Helper: build a synthetic spectrum for a known peptide
//
// The strategy:
//   - For each internal split of `peptide`, compute the theoretical m/z of
//     each ion type in the scorer's trained partition ion list.
//   - Place a peak at that exact m/z with a strictly decreasing intensity
//     (most-intense = split 1 prefix ion gets intensity N*1000, etc.).
//   - Add a few "decoy" peaks at odd m/z values with low intensity.
//
// With this construction, the peak placed for the *first* split's first ion
// gets rank 1.
//
// Returns (spectrum, expected_partition, expected_ion_type, expected_rank=1)
// for the single most-intense matched ion.
// ---------------------------------------------------------------------------

fn build_synthetic_spectrum(
    peptide: &Peptide,
    scorer: &RankScorer,
    precursor_mz: f64,
    charge: u8,
) -> (Spectrum, Partition, IonType) {
    let param = scorer.param();
    let n = peptide.length();
    assert!(n >= 2, "peptide must have at least 2 residues");

    // Compute nominal masses the same way score_psm does.
    let peptide_nominal = peptide.nominal_residue_mass();
    let mut prefix_acc = 0.0_f64;
    let mut prefix_nominals: Vec<i32> = vec![0];
    for s in 1..n {
        let aa = &peptide.residues[s - 1];
        prefix_acc += aa.mass + aa.mod_.as_ref().map_or(0.0, |m| m.mass_delta);
        prefix_nominals.push(nominal_from(prefix_acc));
    }

    // Derive the parent_mass from the precursor.
    const PROTON: f64 = 1.007_276_49;
    let parent_mass = (precursor_mz - PROTON) * (charge as f64);

    // We'll enumerate ions for the first split position, pick the first
    // prefix ion we find (HCD has many y-ions, but the partition will have
    // prefix AND suffix — we want the one we can control most easily).
    // In practice for HCD_QExactive we expect both b-ions (prefix) and y-ions (suffix).
    //
    // Strategy: for split=1, compute ALL ion m/z values; assign intensity
    // proportional to 10000 / (1 + ion_list_index) so the first one has rank 1.
    let num_segs = param.num_segments as usize;

    // Collect ALL (theo_mz, partition, ion_type, intensity) for all splits.
    // Splits processed in order 1..n, intensity descends so split=1 ions are rank-1.
    let mut all_ions: Vec<(f64, Partition, IonType, f32)> = Vec::new();
    let mut global_index = 0usize;

    #[allow(clippy::needless_range_loop)] // `split` indexes parallel arrays
    for split in 1..n {
        let prefix_nom = prefix_nominals[split] as f64;
        let suffix_nom = (peptide_nominal - prefix_nominals[split]) as f64;

        for (is_prefix, nominal_mass) in [(true, prefix_nom), (false, suffix_nom)] {
            for seg in 0..num_segs {
                let partition = param.partition_for(charge, parent_mass, seg);
                let ions = param.ion_types_for_partition_slice(charge, parent_mass, seg);
                for &ion in ions {
                    let theo_mz = match (is_prefix, ion) {
                        (true, IonType::Prefix { .. }) => ion.mz(nominal_mass),
                        (false, IonType::Suffix { .. }) => ion.mz(nominal_mass),
                        _ => continue,
                    };
                    if param.segment_num(theo_mz, parent_mass) != seg {
                        continue;
                    }
                    if theo_mz <= 0.0 {
                        continue;
                    }
                    // Strictly decreasing intensity: first ion gets highest.
                    let intensity = 100_000.0_f32 / (1.0 + global_index as f32);
                    all_ions.push((theo_mz, partition, ion, intensity));
                    global_index += 1;
                }
            }
        }
    }

    assert!(!all_ions.is_empty(), "expected at least one ion from the peptide/scorer combo");

    // Sort peaks by m/z (required by scoring) — this does NOT change relative
    // intensity ordering because we want rank determined by our intensity
    // assignment, which is independent of m/z ordering.
    let mut peaks: Vec<(f64, f32)> = all_ions
        .iter()
        .map(|&(mz, _, _, intensity)| (mz, intensity))
        .collect();
    // Add a few low-intensity decoys at arbitrary m/z far from any ion.
    peaks.push((2000.1, 0.001));
    peaks.push((2001.1, 0.001));
    peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    // Find expected: the first ion in our enumeration is the most intense,
    // so after ranking it will be rank 1.
    let (_, expected_partition, expected_ion, _) = all_ions[0];

    let spectrum = Spectrum {
        title: "synthetic_test".into(),
        precursor_mz,
        precursor_intensity: None,
        precursor_charge: Some(charge as i32),
        rt_seconds: None,
        scan: None,
        peaks,
        activation_method: None,
        isolation_lower_offset: None,
        isolation_upper_offset: None,
    };

    (spectrum, expected_partition, expected_ion)
}

// ---------------------------------------------------------------------------
// TDD test: write this first (it will fail until StatsAccumulator is correct),
// then implement until it passes.
// ---------------------------------------------------------------------------

/// Core known-rank synthetic test.
///
/// Builds a synthetic spectrum where the first ion gets rank 1 by design.
/// Verifies that `StatsAccumulator::accumulate` records `rank_index == 0`
/// (rank 1, index = rank - 1 = 0) for the (partition, ion_type) of that ion
/// in `CountStats::rank`, AND that the charge histogram is bumped.
#[test]
fn accumulator_records_rank_one_for_most_intense_matched_ion() {
    let (_param, scorer) = load_hcd_scorer();

    // Use a peptide long enough to have multiple splits so we exercise both
    // prefix and suffix ions.  "PEPTIDE" is 7 residues; tryptic (ends K/R not
    // required here — we only need the scoring path).
    let peptide = make_peptide(b"PEPTIDE");

    // Precursor: choose mz/charge so parent_mass is mid-range for the fixture.
    // Peptide mass(PEPTIDE) ≈ 799.36 Da; charge=2 → precursor_mz ≈ 400.69.
    let charge: u8 = 2;
    let peptide_mass = peptide.mass();
    const PROTON: f64 = 1.007_276_49;
    let precursor_mz = (peptide_mass + charge as f64 * PROTON) / charge as f64;

    let (spectrum, expected_partition, expected_ion) =
        build_synthetic_spectrum(&peptide, &scorer, precursor_mz, charge);

    // Run the accumulator.
    let accumulator = StatsAccumulator::new(&scorer);
    let mut stats = CountStats::new();
    accumulator.accumulate(&mut stats, &spectrum, &peptide, charge);

    // The first ion (highest intensity) must have been matched at rank 1.
    // In CountStats, rank_index = rank - 1 = 0.
    let rank_index = 0u32;
    let count = stats.rank_count(&expected_partition, expected_ion, rank_index);
    assert!(
        count >= 1,
        "expected rank_count({:?}, {:?}, rank_index={}) >= 1, got {}. \
         Ion list may be empty or partition mismatch.",
        expected_partition, expected_ion, rank_index, count
    );

    // Charge histogram must be bumped exactly once (one PSM).
    let charge_count = stats.charge.get(&(charge as i32)).copied().unwrap_or(0);
    assert_eq!(
        charge_count, 1,
        "expected charge histogram count = 1, got {charge_count}"
    );
}

/// Missing-ion slot test.
///
/// An empty spectrum means every ion is unmatched → the "missing ion" slot
/// (`rank_index = max_rank`) must accumulate counts.
#[test]
fn accumulator_records_missing_ions_on_empty_spectrum() {
    let (_param, scorer) = load_hcd_scorer();
    let peptide = make_peptide(b"PEPTIDE");
    let charge: u8 = 2;

    const PROTON: f64 = 1.007_276_49;
    let precursor_mz = (peptide.mass() + charge as f64 * PROTON) / charge as f64;

    let empty_spec = Spectrum {
        title: "empty".into(),
        precursor_mz,
        precursor_intensity: None,
        precursor_charge: Some(charge as i32),
        rt_seconds: None,
        scan: None,
        peaks: vec![],
        activation_method: None,
        isolation_lower_offset: None,
        isolation_upper_offset: None,
    };

    let accumulator = StatsAccumulator::new(&scorer);
    let mut stats = CountStats::new();
    accumulator.accumulate(&mut stats, &empty_spec, &peptide, charge);

    // With no peaks, every ion is unmatched.  The rank map must have at least
    // one (partition, ion) key with a non-zero count at max_rank index.
    let max_rank = scorer.max_rank();
    let has_missing = stats.rank.iter().any(|(_, vec)| {
        vec.get(max_rank as usize).copied().unwrap_or(0) > 0
    });
    assert!(
        has_missing,
        "expected missing-ion counts at rank_index={max_rank} on empty spectrum"
    );

    // Charge must still be bumped.
    assert_eq!(
        stats.charge.get(&(charge as i32)).copied().unwrap_or(0), 1,
        "charge must be bumped even for empty spectrum"
    );
}

/// Parallel merge test.
///
/// Accumulate into two separate CountStats, then merge.  The merged result
/// must equal the sum of the two individual accumulators (each run on the
/// same PSM, so merged counts = 2 × single).
#[test]
fn merge_doubles_counts_for_same_psm() {
    let (_param, scorer) = load_hcd_scorer();
    let peptide = make_peptide(b"PEPTIDE");
    let charge: u8 = 2;

    const PROTON: f64 = 1.007_276_49;
    let precursor_mz = (peptide.mass() + charge as f64 * PROTON) / charge as f64;
    let (spectrum, _, _) = build_synthetic_spectrum(&peptide, &scorer, precursor_mz, charge);

    let accumulator = StatsAccumulator::new(&scorer);

    let mut part1 = CountStats::new();
    accumulator.accumulate(&mut part1, &spectrum, &peptide, charge);

    let mut part2 = CountStats::new();
    accumulator.accumulate(&mut part2, &spectrum, &peptide, charge);

    // Expected: same as a single accumulation scaled by 2.
    let mut expected = CountStats::new();
    accumulator.accumulate(&mut expected, &spectrum, &peptide, charge);
    accumulator.accumulate(&mut expected, &spectrum, &peptide, charge);

    let merged = merge(vec![part1, part2]);
    assert_eq!(
        merged, expected,
        "merged(part1, part2) must equal 2× single-accumulation"
    );
}

/// `merge` on an empty Vec returns an empty CountStats.
#[test]
fn merge_empty_is_default() {
    let merged = merge(vec![]);
    assert_eq!(merged, CountStats::new(), "merge([]) must equal empty CountStats");
}

/// Verify that the scorer's `ion_match_facts` method (used internally) produces
/// non-empty results for our synthetic peptide — this is a sanity check that the
/// scoring crate's public surface is working correctly before testing the accumulator.
#[test]
fn ion_match_facts_non_empty_for_synthetic_psm() {
    use scoring_crate::ScoredSpectrum;

    let (_param, scorer) = load_hcd_scorer();
    let peptide = make_peptide(b"PEPTIDE");
    let charge: u8 = 2;

    const PROTON: f64 = 1.007_276_49;
    let precursor_mz = (peptide.mass() + charge as f64 * PROTON) / charge as f64;
    let (spectrum, _, _) = build_synthetic_spectrum(&peptide, &scorer, precursor_mz, charge);

    let scored_spec = ScoredSpectrum::new(&spectrum, &scorer, charge);
    let facts = scored_spec.ion_match_facts(&peptide, &scorer);

    assert!(
        !facts.is_empty(),
        "ion_match_facts should return at least one fact for a 7-residue peptide"
    );

    // At least one fact must be matched (rank = Some(_)) since we placed peaks
    // exactly at the theoretical ion positions.
    let matched = facts.iter().filter(|f| f.rank.is_some()).count();
    assert!(
        matched > 0,
        "at least one ion should be matched in the synthetic spectrum (got 0 matched out of {})",
        facts.len()
    );
}
