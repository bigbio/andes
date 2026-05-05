//! PSM scoring integration: ties together ScoredSpectrum (Task 1),
//! predict_by_ions (Task 2), and RankScorer (Task 3) into a single
//! `score_psm` entry point.

use crate::param_model::{IonType, Partition};
use crate::peptide::Peptide;
use crate::scoring::fragment_ions::{predict_by_ions, IonKind};
use crate::scoring::rank_scorer::RankScorer;
use crate::scoring::scored_spectrum::ScoredSpectrum;

/// Score a PSM by matching predicted b/y ions against the spectrum's
/// ranked peaks. Sum of log(ion / noise) per matched ion, plus
/// missing-ion penalties for unmatched.
///
/// `tolerance_da` is the m/z window for nearest-peak lookup. Typically
/// derived from the user's fragment-ion tolerance setting.
///
/// `charge` is the precursor charge; we predict ions at charges
/// `1..=max(1, charge-1)` (matching MS-GF+ default).
pub fn score_psm(
    scored_spec: &ScoredSpectrum,
    peptide: &Peptide,
    scorer: &RankScorer,
    charge: u8,
    tolerance_da: f64,
) -> f32 {
    if charge == 0 {
        return 0.0;
    }
    let partition = pick_partition(peptide, charge, scorer);

    // Predict b/y ions at charges 1..=max(1, charge-1).
    let max_ion_charge = charge.saturating_sub(1).max(1);
    let predicted = predict_by_ions(peptide, 1..=max_ion_charge);

    let mut total = 0.0_f32;
    for p in &predicted {
        let ion_type = ion_kind_to_param_ion_type(p.kind, p.charge);
        match scored_spec.nearest_peak_rank(p.mz, tolerance_da) {
            Some(rank) => total += scorer.node_score(partition, ion_type, rank),
            None => total += scorer.missing_ion_score(partition, ion_type),
        }
    }
    total
}

fn ion_kind_to_param_ion_type(kind: IonKind, charge: u8) -> IonType {
    let charge = charge as i32;
    let offset_bits = 0.0_f32.to_bits();
    match kind {
        IonKind::B => IonType::Prefix { charge, offset_bits },
        IonKind::Y => IonType::Suffix { charge, offset_bits },
    }
}

/// Pick the partition for this (peptide, charge). Scans the scorer's
/// log_table for partitions matching the requested charge, picking the
/// one whose `parent_mass` is closest to the peptide mass. Falls back to
/// a synthesised partition if nothing matches.
fn pick_partition(peptide: &Peptide, charge: u8, scorer: &RankScorer) -> Partition {
    let target_mass = peptide.mass() as f32;
    let z = charge as i32;
    let mut best: Option<(Partition, f32)> = None;
    for ((part, _ion), _) in scorer.log_table.iter() {
        if part.charge != z {
            continue;
        }
        let dist = (part.parent_mass - target_mass).abs();
        if best.as_ref().map_or(true, |(_, d)| dist < *d) {
            best = Some((*part, dist));
        }
    }
    best.map(|(p, _)| p).unwrap_or(Partition {
        charge: z,
        parent_mass: target_mass,
        seg_num: 0,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::amino_acid::AminoAcid;
    use crate::param_model::{IonType, Param, Partition};
    use crate::peptide::Peptide;
    use crate::scoring::fragment_ions::predict_by_ions;
    use crate::scoring::rank_scorer::RankScorer;
    use crate::scoring::scored_spectrum::ScoredSpectrum;
    use crate::spectrum::Spectrum;
    use std::collections::HashMap;

    /// Construct a minimal Param with one partition + Prefix(1) + Suffix(1) + Noise.
    /// Noise has uniform low frequency, ions have higher frequency at rank 1
    /// so a perfect match scores positively.
    fn tiny_param() -> Param {
        use crate::activation::ActivationMethod;
        use crate::instrument::InstrumentType;
        use crate::param_model::SpecDataType;
        use crate::protocol::Protocol;
        use crate::tolerance::Tolerance;

        let part = Partition { charge: 2, parent_mass: 500.0, seg_num: 0 };
        let prefix1 = IonType::Prefix { charge: 1, offset_bits: 0.0_f32.to_bits() };
        let suffix1 = IonType::Suffix { charge: 1, offset_bits: 0.0_f32.to_bits() };
        let noise = IonType::Noise;

        // max_rank=3 → arrays of length 4. ion freqs high at rank 1, low at missing.
        let prefix_freqs = vec![0.5_f32, 0.1, 0.05, 0.01];
        let suffix_freqs = vec![0.5_f32, 0.1, 0.05, 0.01];
        let noise_freqs = vec![0.05_f32, 0.05, 0.05, 0.05];

        let mut ion_table = HashMap::new();
        ion_table.insert(prefix1, prefix_freqs);
        ion_table.insert(suffix1, suffix_freqs);
        ion_table.insert(noise, noise_freqs);

        let mut rank_dist_table = HashMap::new();
        rank_dist_table.insert(part, ion_table);

        let mut frag_off_table = HashMap::new();
        frag_off_table.insert(part, vec![]);

        Param {
            version: 10001,
            data_type: SpecDataType {
                activation: ActivationMethod::HCD,
                instrument: InstrumentType::QExactive,
                enzyme: None,
                protocol: Protocol::Automatic,
            },
            mme: Tolerance::Ppm(20.0),
            apply_deconvolution: false,
            deconvolution_error_tolerance: 0.0,
            charge_hist: vec![(2, 100)],
            min_charge: 2,
            max_charge: 2,
            num_segments: 1,
            partitions: vec![part],
            num_precursor_off: 0,
            precursor_off_map: HashMap::new(),
            frag_off_table,
            max_rank: 3,
            rank_dist_table,
            error_scaling_factor: 0,
            ion_err_dist_table: HashMap::new(),
            noise_err_dist_table: HashMap::new(),
            ion_existence_table: HashMap::new(),
        }
    }

    fn pep(seq: &[u8]) -> Peptide {
        let residues: Vec<AminoAcid> = seq
            .iter()
            .map(|&r| AminoAcid::standard(r).unwrap())
            .collect();
        Peptide::new(residues, b'_', b'-')
    }

    #[test]
    fn empty_spectrum_returns_all_missing_score() {
        // No peaks → every predicted ion is missing → score = sum of missing_ion_score.
        let peptide = pep(b"AGR");
        let spec = Spectrum {
            title: "empty".into(),
            precursor_mz: 0.0,
            precursor_intensity: None,
            precursor_charge: Some(2),
            rt_seconds: None,
            scan: None,
            peaks: vec![],
        };
        let scored = ScoredSpectrum::new(&spec);
        let param = tiny_param();
        let scorer = RankScorer::new(&param);
        let s = score_psm(&scored, &peptide, &scorer, 2, 0.1);
        // No peak matches; all ions missing. Score = (2 b-ions + 2 y-ions) * missing_score.
        // For prefix1: missing = log(0.01 / (0.05 * min(1, 1))) = log(0.2) ≈ -1.609
        // Same for suffix1.
        // 4 ions × log(0.2) = 4 × -1.609 ≈ -6.44. Score should be NEGATIVE.
        assert!(s < 0.0, "score should be negative on empty spectrum, got {s}");
    }

    #[test]
    fn perfect_match_yields_positive_score() {
        // Build a spectrum where the highest-intensity peaks exactly match
        // every predicted b/y ion of "AGR" at charge 1.
        let peptide = pep(b"AGR");
        let predicted = predict_by_ions(&peptide, 1..=1);

        // Make each predicted ion the highest-intensity peak.
        let mut peaks: Vec<(f64, f32)> = predicted.iter().map(|p| (p.mz, 1000.0)).collect();
        peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        let spec = Spectrum {
            title: "match".into(),
            precursor_mz: 500.0,
            precursor_intensity: None,
            precursor_charge: Some(2),
            rt_seconds: None,
            scan: None,
            peaks,
        };
        let scored = ScoredSpectrum::new(&spec);
        let param = tiny_param();
        let scorer = RankScorer::new(&param);
        let s = score_psm(&scored, &peptide, &scorer, 2, 0.1);
        // Each ion finds a high-rank peak (rank 1..N out of N). Score should be positive.
        // For prefix1 rank 1: log(0.5 / (0.05 * 1)) = log(10) ≈ 2.30. Same for suffix1.
        assert!(s > 0.0, "perfect match score should be positive, got {s}");
    }

    #[test]
    fn perfect_match_outscores_empty_spectrum() {
        let peptide = pep(b"AGR");
        let predicted = predict_by_ions(&peptide, 1..=1);

        let mut match_peaks: Vec<(f64, f32)> = predicted.iter().map(|p| (p.mz, 1000.0)).collect();
        match_peaks.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        let match_spec = Spectrum {
            title: "match".into(),
            precursor_mz: 500.0,
            precursor_intensity: None,
            precursor_charge: Some(2),
            rt_seconds: None,
            scan: None,
            peaks: match_peaks,
        };
        let empty_spec = Spectrum {
            title: "empty".into(),
            precursor_mz: 500.0,
            precursor_intensity: None,
            precursor_charge: Some(2),
            rt_seconds: None,
            scan: None,
            peaks: vec![],
        };

        let param = tiny_param();
        let scorer = RankScorer::new(&param);
        let s_match = score_psm(&ScoredSpectrum::new(&match_spec), &peptide, &scorer, 2, 0.1);
        let s_empty = score_psm(&ScoredSpectrum::new(&empty_spec), &peptide, &scorer, 2, 0.1);
        assert!(s_match > s_empty);
    }
}
