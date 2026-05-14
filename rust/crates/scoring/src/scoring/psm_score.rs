//! PSM scoring integration.
//!
//! `score_psm` sums `ScoredSpectrum::node_score(prefix, suffix)` across each
//! peptide split position. The result is on the same score scale used by the
//! GF DP, so `GeneratingFunctionGroup::spectral_probability(psm.score)` is
//! calibrated.
//!
//! Per-split node score: `round(getNodeScore(prm, true) + getNodeScore(srm, false))`
//! where `prm` is the nominal prefix mass and `srm = peptideMass - prm`.

use model::mass::nominal_from;
use model::peptide::Peptide;
use crate::scoring::rank_scorer::RankScorer;
use crate::scoring::scored_spectrum::ScoredSpectrum;

/// Score a PSM as the sum of `ScoredSpectrum::node_score(prefix, suffix)`
/// across each peptide split position.  This produces a raw score on the
/// same scale as the GF distribution so that `GeneratingFunctionGroup::
/// spectral_probability(psm.score.round() as i32)` is calibrated.
///
/// For each split `i` in `1..n`:
/// - `nominal_prefix_mass[i] = nominal_from(sum of residues 0..i)`
/// - `peptide_mass = nominal_prefix_mass[n-1]` = nominal AA-only sum
/// - `score += round(prefix_score[prm] + suffix_score[srm])`
///
/// `fragment_tolerance_da` is forwarded to `ScoredSpectrum::node_score` for
/// peak-lookup.  The `charge` selects the partition; `parent_mass` is the
/// peptide neutral mass (residue_sum + H₂O), used for segment selection.
pub fn score_psm(
    scored_spec: &ScoredSpectrum,
    peptide: &Peptide,
    scorer: &RankScorer,
    charge: u8,
    fragment_tolerance_da: f64,
) -> f32 {
    if charge == 0 {
        return 0.0;
    }
    let n = peptide.length();
    if n < 2 {
        return 0.0;
    }

    // Two distinct masses with different roles:
    //  - `peptide_nominal`: candidate peptide's total nominal residue mass.
    //    Drives suffix lookup, built from the candidate's residues.
    //  - `spectrum_parent_mass`: spectrum's OBSERVED neutral mass.
    //    Drives partition + segment selection across all candidates,
    //    regardless of iso_off. Using `peptide.mass()` here would mismatch
    //    iso_off≥1 candidates and cause systematic top-1 flips.
    let spectrum_parent_mass = scored_spec.parent_mass();

    // Total nominal peptide mass = nominal(residue_sum) = nominal(mass - H2O).
    // Used to compute suffix_nominal = peptide_nominal - prefix_nominal.
    let peptide_nominal = peptide.nominal_residue_mass();

    // ── Score-traceability instrumentation ─────────────────────────────────
    // Gated by the `MSGF_TRACE_PEP` env var: if the peptide's unmodified
    // residue sequence contains the filter string, emit per-split trace
    // lines on stderr. Mirrors `FastScorer.getScoreWithTrace`, so the two
    // dumps line up split-by-split.
    let trace_pep = std::env::var("MSGF_TRACE_PEP").ok();
    let pep_seq_string: String = peptide.residues.iter().map(|aa| aa.residue as char).collect();
    let trace = trace_pep
        .as_ref()
        .map(|p| !p.is_empty() && pep_seq_string.contains(p.as_str()))
        .unwrap_or(false);
    if trace {
        eprintln!(
            "TRACE_RUST_HEADER\tpep={}\tcharge={}\tparent_mass={:.4}\tpeptide_nominal={}\tn={}\tfragment_tol_da={}",
            pep_seq_string, charge, spectrum_parent_mass, peptide_nominal, n, fragment_tolerance_da
        );
    }

    let mut total: i32 = 0;
    let mut prefix_mass_acc = 0.0_f64;
    // Split positions 1..n: after split s, prefix = residues[0..s], suffix = residues[s..n].
    for s in 1..n {
        // Accumulate exact float mass for residue s-1 (0-indexed).
        let aa = &peptide.residues[s - 1];
        let residue_mass = aa.mass + aa.mod_.as_ref().map_or(0.0, |m| m.mass_delta);
        prefix_mass_acc += residue_mass;

        // Nominal masses at the split position.
        let prefix_nominal = nominal_from(prefix_mass_acc);
        let suffix_nominal = peptide_nominal - prefix_nominal;

        // Sample cached prefix/suffix scores for trace diagnostics (no
        // change to the hot-path semantics): the cache exposes the same
        // FastScorer-style per-index tables Java prints from.
        let cached_pref = scored_spec.cached_prefix_score(prefix_nominal);
        let cached_suff = scored_spec.cached_suffix_score(suffix_nominal);

        let contribution = scored_spec
            .cached_split_score(prefix_nominal, suffix_nominal)
            .unwrap_or_else(|| {
                scored_spec.node_score(
                    prefix_nominal as f64,
                    suffix_nominal as f64,
                    scorer,
                    charge,
                    spectrum_parent_mass,
                    fragment_tolerance_da,
                )
            });
        total += contribution;

        if trace {
            let pref_str = cached_pref
                .map(|v| format!("{v}"))
                .unwrap_or_else(|| "NA".to_string());
            let suff_str = cached_suff
                .map(|v| format!("{v}"))
                .unwrap_or_else(|| "NA".to_string());
            eprintln!(
                "TRACE_RUST\tpep={}\tsplit={}\tprefMass={}\tsuffMass={}\tprefScore={}\tsuffScore={}\tcontribution={}\tcumulative={}\tprefAccF64={:.6}",
                pep_seq_string, s, prefix_nominal, suffix_nominal,
                pref_str, suff_str, contribution, total, prefix_mass_acc
            );
        }
    }
    if trace {
        eprintln!(
            "TRACE_RUST_FINAL\tpep={}\trawScore={}",
            pep_seq_string, total
        );
    }
    total as f32
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::amino_acid::AminoAcid;
    use crate::param_model::{FragmentOffsetFrequency, IonType, Param, Partition, SpecDataType};
    use model::peptide::Peptide;
    use crate::scoring::rank_scorer::RankScorer;
    use crate::scoring::scored_spectrum::ScoredSpectrum;
    use model::spectrum::Spectrum;
    use crate::testutil::tiny_param;
    use std::collections::HashMap;

    fn pep(seq: &[u8]) -> Peptide {
        let residues: Vec<AminoAcid> = seq
            .iter()
            .map(|&r| AminoAcid::standard(r).unwrap())
            .collect();
        Peptide::new(residues, b'_', b'-')
    }

    fn empty_spectrum(title: &str) -> Spectrum {
        Spectrum {
            title: title.into(),
            precursor_mz: 0.0,
            precursor_intensity: None,
            precursor_charge: Some(2),
            rt_seconds: None,
            scan: None,
            peaks: vec![],
        }
    }

    /// A param whose single partition has `parent_mass = 0.0`, so the floor-
    /// matching in `find_partition` returns it for *any* peptide mass.
    /// The prefix-ion frequencies are tuned so that rank-1 hits score positive.
    fn any_mass_param() -> Param {
        use model::activation::ActivationMethod;
        use model::instrument::InstrumentType;
        use model::protocol::Protocol;

        let part = Partition { charge: 2, parent_mass: 0.0, seg_num: 0 };
        let prefix_ion = IonType::Prefix { charge: 1, offset_bits: 0.0_f32.to_bits() };
        let noise_ion = IonType::Noise;

        let ion_freqs = vec![0.6_f32, 0.3, 0.05, 0.001];
        let noise_freqs = vec![0.1_f32, 0.2, 0.3, 0.4];

        let mut ion_table: HashMap<IonType, Vec<f32>> = HashMap::new();
        ion_table.insert(prefix_ion, ion_freqs);
        ion_table.insert(noise_ion, noise_freqs);

        let mut rank_dist_table: HashMap<Partition, HashMap<IonType, Vec<f32>>> = HashMap::new();
        rank_dist_table.insert(part, ion_table);

        let mut frag_off_table = HashMap::new();
        frag_off_table.insert(part, vec![FragmentOffsetFrequency { ion_type: prefix_ion, frequency: 0.7 }]);

        let mut p = Param {
            version: 10001,
            data_type: SpecDataType {
                activation: ActivationMethod::HCD,
                instrument: InstrumentType::QExactive,
                enzyme: None,
                protocol: Protocol::Automatic,
            },
            mme: model::tolerance::Tolerance::Da(0.2),
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
            partition_ion_types_cache: HashMap::new(),
        };
        p.rebuild_cache();
        p
    }

    #[test]
    fn empty_spectrum_returns_non_positive_score() {
        // No peaks → every node lookup is missing → score ≤ 0.
        // (With node_score iterating all ion types, missing_ion_score is
        // negative for all configured ions; the sum is non-positive.)
        let peptide = pep(b"AGR");
        let spec = empty_spectrum("empty");
        let scored = ScoredSpectrum::new_without_filtering(&spec);
        let param = any_mass_param();
        let scorer = RankScorer::new(&param);
        let s = score_psm(&scored, &peptide, &scorer, 2, 0.2);
        assert!(s <= 0.0, "score should be ≤ 0 on empty spectrum, got {s}");
    }

    #[test]
    fn perfect_match_yields_positive_score() {
        // Build a spectrum whose peaks fall exactly at the b-ion m/z of each
        // split position.  Uses `any_mass_param` so the partition lookup
        // succeeds for the small AGR peptide mass.
        let peptide = pep(b"AGR");
        let param = any_mass_param();

        // Compute b-ion m/z for each split position of AGR.
        let b_ion = IonType::Prefix { charge: 1, offset_bits: 0.0_f32.to_bits() };
        let mut prefix_acc = 0.0_f64;
        let mut peaks = Vec::new();
        for s in 1..peptide.length() {
            let aa = &peptide.residues[s - 1];
            prefix_acc += aa.mass;
            let nom = model::mass::nominal_from(prefix_acc) as f64;
            let mz = b_ion.mz(nom);
            peaks.push((mz, 1000.0_f32 / s as f32));  // rank-1 intensity
        }
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
        let scored = ScoredSpectrum::new_without_filtering(&spec);
        let scorer = RankScorer::new(&param);
        let s = score_psm(&scored, &peptide, &scorer, 2, 0.2);
        assert!(s > 0.0, "score with matched b-ions should be positive, got {s}");
    }

    #[test]
    fn perfect_match_outscores_empty_spectrum() {
        // A spectrum with matched peaks must outscore an empty spectrum.
        let peptide = pep(b"AGR");
        let param = any_mass_param();

        let b_ion = IonType::Prefix { charge: 1, offset_bits: 0.0_f32.to_bits() };
        let mut prefix_acc = 0.0_f64;
        let mut match_peaks = Vec::new();
        for s in 1..peptide.length() {
            let aa = &peptide.residues[s - 1];
            prefix_acc += aa.mass;
            let nom = model::mass::nominal_from(prefix_acc) as f64;
            let mz = b_ion.mz(nom);
            match_peaks.push((mz, 1000.0_f32));
        }
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

        let scorer = RankScorer::new(&param);
        let scored_match = ScoredSpectrum::new_without_filtering(&match_spec);
        let empty_spec = empty_spectrum("empty");
        let scored_empty = ScoredSpectrum::new_without_filtering(&empty_spec);
        let s_match = score_psm(&scored_match, &peptide, &scorer, 2, 0.2);
        let s_empty = score_psm(&scored_empty, &peptide, &scorer, 2, 0.2);
        assert!(s_match > s_empty, "matched spectrum ({s_match}) should outscore empty ({s_empty})");
    }

    /// Verify that `score_psm` equals the manually summed `node_score` calls
    /// across each split position (this is the definition of the new formula).
    #[test]
    fn score_psm_matches_sum_of_node_scores_across_splits() {
        use model::amino_acid::AminoAcid;
        use model::mass::nominal_from;

        let peptide = pep(b"AGR");
        let param = tiny_param();
        let scorer = RankScorer::new(&param);

        // Empty spectrum — all node scores are missing, but the sum should still match.
        let empty_spec = empty_spectrum("empty");
        let scored = ScoredSpectrum::new_without_filtering(&empty_spec);

        let parent_mass = peptide.mass();
        let peptide_nominal = peptide.nominal_residue_mass();
        let charge = 2u8;
        let tolerance_da = 0.05;

        let mut manual_total: i32 = 0;
        let mut prefix_acc = 0.0_f64;
        for s in 1..peptide.length() {
            let aa: &AminoAcid = &peptide.residues[s - 1];
            prefix_acc += aa.mass + aa.mod_.as_ref().map_or(0.0, |m| m.mass_delta);
            let pref = nominal_from(prefix_acc);
            let suf = peptide_nominal - pref;
            manual_total += scored.node_score(pref as f64, suf as f64, &scorer, charge, parent_mass, tolerance_da);
        }

        let computed = score_psm(&scored, &peptide, &scorer, charge, tolerance_da);
        assert_eq!(
            computed as i32, manual_total,
            "score_psm ({computed}) should equal manual split sum ({manual_total})"
        );
    }

    #[test]
    fn charge_zero_returns_zero() {
        let peptide = pep(b"AGR");
        let param = tiny_param();
        let scorer = RankScorer::new(&param);
        let spec = empty_spectrum("empty");
        let scored = ScoredSpectrum::new_without_filtering(&spec);
        assert_eq!(score_psm(&scored, &peptide, &scorer, 0, 0.1), 0.0);
    }

    #[test]
    fn single_residue_peptide_returns_zero() {
        let peptide = pep(b"A");
        let param = tiny_param();
        let scorer = RankScorer::new(&param);
        let spec = empty_spectrum("empty");
        let scored = ScoredSpectrum::new_without_filtering(&spec);
        assert_eq!(score_psm(&scored, &peptide, &scorer, 2, 0.1), 0.0);
    }
}
