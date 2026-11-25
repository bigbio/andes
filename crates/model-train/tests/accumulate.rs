use std::path::Path;

use model::amino_acid::AminoAcid;
use model::mass::nominal_from;
use model::peptide::Peptide;
use model::spectrum::Spectrum;
use scoring_crate::params::{IonType, Param, Partition};
use scoring_crate::psm::rank::RankScorer;

use model_train::accumulate::{merge, StatsAccumulator};
use model_train::counts::CountStats;

fn make_peptide(seq: &[u8]) -> Peptide {
    let residues: Vec<AminoAcid> = seq
        .iter()
        .map(|&r| AminoAcid::standard(r).unwrap())
        .collect();
    Peptide::new(residues, b'_', b'-')
}

fn load_hcd_scorer() -> (Param, RankScorer) {
    let bundled = Path::new(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../resources/models"
    ));
    let store = model_train::store::ModelStore::open(bundled).expect("open bundled model store");
    let param = store
        .load_param("hcd_qexactive_tryp")
        .expect("load hcd_qexactive_tryp from store");
    let scorer = RankScorer::new(&param);
    (param, scorer)
}

fn build_synthetic_spectrum(
    peptide: &Peptide,
    scorer: &RankScorer,
    precursor_mz: f64,
    charge: u8,
) -> (Spectrum, Partition, IonType) {
    let param = scorer.param();
    let n = peptide.length();
    assert!(n >= 2, "peptide must have at least 2 residues");

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

    let mut all_ions: Vec<(f64, Partition, IonType, f32)> = Vec::new();
    let mut global_index = 0usize;

    #[allow(clippy::needless_range_loop)]
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
                    let intensity = 100_000.0_f32 / (1.0 + global_index as f32);
                    all_ions.push((theo_mz, partition, ion, intensity));
                    global_index += 1;
                }
            }
        }
    }

    assert!(
        !all_ions.is_empty(),
        "expected at least one ion from the peptide/scorer combo"
    );

    let mut peaks: Vec<(f64, f32)> = all_ions
        .iter()
        .map(|&(mz, _, _, intensity)| (mz, intensity))
        .collect();
    peaks.push((2000.1, 0.001));
    peaks.push((2001.1, 0.001));
    peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

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
        isolation_target_mz: None,
        mobility_window: None,
    };

    (spectrum, expected_partition, expected_ion)
}

#[test]
fn accumulator_records_rank_one_for_most_intense_matched_ion() {
    let (_param, scorer) = load_hcd_scorer();

    let peptide = make_peptide(b"PEPTIDE");

    let charge: u8 = 2;
    let peptide_mass = peptide.mass();
    const PROTON: f64 = 1.007_276_49;
    let precursor_mz = (peptide_mass + charge as f64 * PROTON) / charge as f64;

    let (spectrum, expected_partition, expected_ion) =
        build_synthetic_spectrum(&peptide, &scorer, precursor_mz, charge);

    let accumulator = StatsAccumulator::new(&scorer);
    let mut stats = CountStats::new();
    accumulator.accumulate(&mut stats, &spectrum, &peptide, charge);

    let rank_index = 0u32;
    let count = stats.rank_count(&expected_partition, expected_ion, rank_index);
    assert!(
        count >= 1,
        "expected rank_count({expected_partition:?}, {expected_ion:?}, rank_index={rank_index}) >= 1, got {count}. \
         Ion list may be empty or partition mismatch.",
    );

    let charge_count = stats.charge.get(&(charge as i32)).copied().unwrap_or(0);
    assert_eq!(
        charge_count, 1,
        "expected charge histogram count = 1, got {charge_count}"
    );
}

#[test]
fn accumulator_records_ion_existence_both_present() {
    let (_param, scorer) = load_hcd_scorer();
    let peptide = make_peptide(b"PEPTIDE");
    let charge: u8 = 2;

    const PROTON: f64 = 1.007_276_49;
    let precursor_mz = (peptide.mass() + charge as f64 * PROTON) / charge as f64;
    let (spectrum, _, _) = build_synthetic_spectrum(&peptide, &scorer, precursor_mz, charge);

    let accumulator = StatsAccumulator::new(&scorer);
    let mut stats = CountStats::new();
    accumulator.accumulate(&mut stats, &spectrum, &peptide, charge);

    assert!(
        !stats.existence.is_empty(),
        "existence counts must be recorded (was empty → uniform-0.25 bug)"
    );
    let both_present: u64 = stats
        .existence
        .iter()
        .filter(|((_, idx), _)| *idx == 3)
        .map(|(_, &c)| c)
        .sum();
    assert!(
        both_present >= 1,
        "expected >=1 existence count at index 3 (both nodes observed), got {both_present}"
    );
}

#[test]
fn accumulator_records_ion_existence_both_absent_on_empty_spectrum() {
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
        isolation_target_mz: None,
        mobility_window: None,
    };

    let accumulator = StatsAccumulator::new(&scorer);
    let mut stats = CountStats::new();
    accumulator.accumulate(&mut stats, &empty_spec, &peptide, charge);

    assert!(
        !stats.existence.is_empty(),
        "existence counts must be recorded"
    );
    let nonzero_idx: Vec<u32> = stats
        .existence
        .iter()
        .filter(|(_, &c)| c > 0)
        .map(|((_, idx), _)| *idx)
        .collect();
    assert!(
        nonzero_idx.iter().all(|&idx| idx == 0 || idx == 2),
        "empty spectrum → no real peak matched, so the 'cur observed' bit must \
         never be set; expected indices ⊆ {{0, 2}}, got {nonzero_idx:?}"
    );
}

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
        isolation_target_mz: None,
        mobility_window: None,
    };

    let accumulator = StatsAccumulator::new(&scorer);
    let mut stats = CountStats::new();
    accumulator.accumulate(&mut stats, &empty_spec, &peptide, charge);

    let max_rank = scorer.max_rank();
    let has_missing = stats
        .rank
        .iter()
        .any(|(_, vec)| vec.get(max_rank as usize).copied().unwrap_or(0) > 0);
    assert!(
        has_missing,
        "expected missing-ion counts at rank_index={max_rank} on empty spectrum"
    );

    assert_eq!(
        stats.charge.get(&(charge as i32)).copied().unwrap_or(0),
        1,
        "charge must be bumped even for empty spectrum"
    );
}

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

    let mut expected = CountStats::new();
    accumulator.accumulate(&mut expected, &spectrum, &peptide, charge);
    accumulator.accumulate(&mut expected, &spectrum, &peptide, charge);

    let merged = merge(vec![part1, part2]);
    assert_eq!(
        merged, expected,
        "merged(part1, part2) must equal 2× single-accumulation"
    );
}

#[test]
fn merge_empty_is_default() {
    let merged = merge(vec![]);
    assert_eq!(
        merged,
        CountStats::new(),
        "merge([]) must equal empty CountStats"
    );
}

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

    let matched = facts.iter().filter(|f| f.rank.is_some()).count();
    assert!(
        matched > 0,
        "at least one ion should be matched in the synthetic spectrum (got 0 matched out of {})",
        facts.len()
    );
}
