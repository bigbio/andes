use model::mass::{nominal_from, H2O, INTEGER_MASS_SCALER, PROTON};
use model::peptide::Peptide;
use model::spectrum::Spectrum;
use model::stats::{median_mad, Even, MAD_TO_SD};
use scoring_crate::scoring::fragment_ions::{predict_by_ions, IonKind, PredictedIon};
use scoring_crate::scoring::FragmentErrorModel;

use crate::match_engine::PreparedSearch;
use crate::tdc::{confident_target_indices, ScoredLabel};

pub const MEASURE_WINDOW_PPM: f64 = 50.0;
pub const TOP_IONS_PER_PSM: usize = 7;
pub const MIN_CONFIDENT_PSMS: usize = 100;
pub const MEASURE_Q: f64 = 0.01;
pub const SAMPLE_SPECTRA: usize = 3000;
/// Spectra a `--dia` run measures on, spread over the run; with several candidates per scan
/// they hold many times [`MIN_CONFIDENT_PSMS`] confident PSMs.
pub const DIA_SAMPLE_SPECTRA: usize = 3_000;
/// Candidates per `--dia` scan that enter the measurement, best first: a DIA scan holds
/// several co-isolated peptides, and every confident one carries its own fragment errors.
pub const DIA_RANKS_PER_SPECTRUM: usize = 10;

/// The spectra a run's fragment-error measurement searches and the candidates it keeps.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MeasurePlan {
    /// Spectra searched, a contiguous block in precursor-mass order.
    pub sample_spectra: usize,
    /// Best candidates per spectrum that compete for 1 % q.
    pub ranks_per_spectrum: usize,
    /// Spread the sample evenly over the precursor-mass order instead of taking a contiguous
    /// block from its middle, so every isolation window and the whole gradient are sampled.
    pub strided: bool,
}

impl MeasurePlan {
    /// Positions, in precursor-mass order, of the spectra of a run of `len` spectra that the
    /// measurement searches, ascending.
    pub fn sample_positions(&self, len: usize) -> Vec<usize> {
        let n = self.sample_spectra.min(len);
        if !self.strided {
            let lo = (len / 2).saturating_sub(n / 2);
            return (lo..lo + n).collect();
        }
        (0..n).map(|i| i * len / n).collect()
    }
}

impl MeasurePlan {
    /// One candidate per spectrum from a small sample, for DDA.
    pub const DDA: MeasurePlan = MeasurePlan {
        sample_spectra: SAMPLE_SPECTRA,
        ranks_per_spectrum: 1,
        strided: false,
    };
    /// Several candidates per spectrum from a larger sample, for `--dia`.
    pub const DIA: MeasurePlan = MeasurePlan {
        sample_spectra: DIA_SAMPLE_SPECTRA,
        ranks_per_spectrum: DIA_RANKS_PER_SPECTRUM,
        strided: true,
    };
}

fn nominal_reference_mz(ion: &PredictedIon, peptide: &Peptide) -> f64 {
    let scaler = INTEGER_MASS_SCALER as f64;
    match ion.kind {
        IonKind::Y => {
            let suffix = ion.mz - H2O - PROTON;
            let prefix_nominal = nominal_from(peptide.residue_mass() - suffix);
            (peptide.nominal_residue_mass() - prefix_nominal) as f64 / scaler + H2O + PROTON
        }
        _ => nominal_from(ion.mz - PROTON) as f64 / scaler + PROTON,
    }
}

pub fn top_ion_errors(spec: &Spectrum, peptide: &Peptide, k: usize) -> Vec<(f64, f64)> {
    let peaks = &spec.peaks;
    let mut matched: Vec<(f32, (f64, f64))> = Vec::new();
    for ion in predict_by_ions(peptide, 1..=1) {
        let tol = ion.mz * MEASURE_WINDOW_PPM * 1e-6;
        let lo = peaks.partition_point(|p| p.0 < ion.mz - tol);
        let best = peaks[lo..]
            .iter()
            .take_while(|p| p.0 <= ion.mz + tol)
            .min_by(|a, b| (a.0 - ion.mz).abs().total_cmp(&(b.0 - ion.mz).abs()));
        if let Some(&(mz, intensity)) = best {
            let ppm = (mz - ion.mz) / ion.mz * 1e6;
            matched.push((intensity, (ppm, mz - nominal_reference_mz(&ion, peptide))));
        }
    }
    matched.sort_by(|a, b| b.0.total_cmp(&a.0));
    matched.truncate(k);
    matched.into_iter().map(|(_, e)| e).collect()
}

fn centre_and_spread(errors: &[f64]) -> (f64, f64) {
    median_mad(errors, Even::Mean).map_or((f64::NAN, f64::NAN), |(c, mad)| (c, MAD_TO_SD * mad))
}

pub fn estimate(errors: &[(f64, f64)], psms: usize) -> Option<FragmentErrorModel> {
    if errors.len() < MIN_CONFIDENT_PSMS {
        return None;
    }
    let ppm: Vec<f64> = errors.iter().map(|e| e.0).collect();
    let (c0, s0) = centre_and_spread(&ppm);
    let half = (4.0 * s0).max(1.0);
    let kept: Vec<(f64, f64)> = errors
        .iter()
        .copied()
        .filter(|e| (e.0 - c0).abs() <= half)
        .collect();
    if kept.len() < MIN_CONFIDENT_PSMS {
        return None;
    }
    let ppm: Vec<f64> = kept.iter().map(|e| e.0).collect();
    let nominal: Vec<f64> = kept.iter().map(|e| e.1).collect();
    let (offset_ppm, sigma_ppm) = centre_and_spread(&ppm);
    let (node_centre_da, node_sigma_da) = centre_and_spread(&nominal);
    Some(
        FragmentErrorModel::new(offset_ppm, sigma_ppm, kept.len(), psms)
            .with_node_error_da(node_centre_da, node_sigma_da),
    )
}

#[derive(Debug, Clone, PartialEq)]
pub enum NotMeasured {
    TooFewPsms { confident: usize, scored: usize },
    TooFewIons { confident: usize },
}

impl std::fmt::Display for NotMeasured {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            NotMeasured::TooFewPsms { confident, scored } => write!(
                f,
                "{confident} of {scored} sampled PSMs are targets at {MEASURE_Q} q-value \
                 (need {MIN_CONFIDENT_PSMS})"
            ),
            NotMeasured::TooFewIons { confident } => {
                write!(f, "too few matched ions in {confident} confident PSMs")
            }
        }
    }
}

/// The best `ranks` candidates of each queue as (spectrum, candidate, rank score), best first
/// within a spectrum.
pub fn ranked_candidates(queues: &[crate::TopNQueue], ranks: usize) -> Vec<(usize, u32, f32)> {
    let mut out = Vec::new();
    for (i, q) in queues.iter().enumerate() {
        if ranks <= 1 {
            if let Some(p) = q.peek_top() {
                out.push((i, p.primary_candidate_idx(), p.rank_score));
            }
            continue;
        }
        let mut v: Vec<&crate::PsmMatch> = q.iter_psms().collect();
        v.sort_by(|a, b| b.cmp(a));
        out.extend(
            v.into_iter()
                .take(ranks)
                .map(|p| (i, p.primary_candidate_idx(), p.rank_score)),
        );
    }
    out
}

pub fn measure_on_queues(
    prepared: &PreparedSearch<'_>,
    sample: &[Spectrum],
    queues: &[crate::TopNQueue],
    ranks: usize,
) -> Result<FragmentErrorModel, NotMeasured> {
    let mut tops: Vec<(usize, u32)> = Vec::new();
    let mut labels: Vec<ScoredLabel> = Vec::new();
    for (i, idx, score) in ranked_candidates(queues, ranks) {
        tops.push((i, idx));
        labels.push(ScoredLabel {
            score,
            is_decoy: prepared.candidate_is_decoy(idx),
        });
    }
    let confident = confident_target_indices(&labels, MEASURE_Q);
    if confident.len() < MIN_CONFIDENT_PSMS {
        return Err(NotMeasured::TooFewPsms {
            confident: confident.len(),
            scored: labels.len(),
        });
    }
    let mut errors: Vec<(f64, f64)> = Vec::new();
    for &k in &confident {
        let (si, idx) = tops[k];
        let Some(spec) = sample.get(si) else {
            continue;
        };
        let peptide = prepared.candidate_peptide(idx);
        errors.extend(top_ion_errors(spec, &peptide, TOP_IONS_PER_PSM));
    }
    estimate(&errors, confident.len()).ok_or(NotMeasured::TooFewIons {
        confident: confident.len(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::AminoAcid;

    fn peptide(residues: &[u8]) -> Peptide {
        let aas: Vec<AminoAcid> = residues
            .iter()
            .map(|&r| AminoAcid::standard(r).unwrap())
            .collect();
        Peptide::new(aas, b'K', b'K')
    }

    #[test]
    fn top_ion_errors_take_the_most_intense_matches() {
        let p = peptide(b"PEPTIDEK");
        let ions = predict_by_ions(&p, 1..=1);
        let mut peaks: Vec<(f64, f32)> = ions
            .iter()
            .enumerate()
            .map(|(i, ion)| (ion.mz * (1.0 + 3.0e-6), 100.0 + i as f32))
            .collect();
        peaks[0].0 = ions[0].mz * (1.0 + 80.0e-6);
        peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
        let spec = Spectrum {
            peaks,
            ..Default::default()
        };
        let errs = top_ion_errors(&spec, &p, 7);
        assert_eq!(errs.len(), 7);
        assert!(errs.iter().all(|e| (e.0 - 3.0).abs() < 1e-6), "{errs:?}");
        let all = top_ion_errors(&spec, &p, 100);
        assert_eq!(all.len(), ions.len() - 1);
        assert!(all.iter().all(|e| e.1.abs() < 0.2), "{all:?}");
        assert!(all.iter().any(|e| e.1.abs() > 0.005), "{all:?}");
    }

    #[test]
    fn estimate_is_robust_to_chance_matches() {
        let mut state = 12345u64;
        let mut uniform = || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        let mut errs = Vec::new();
        for i in 0..5000 {
            if i % 10 == 0 {
                errs.push(-50.0 + 100.0 * uniform());
            } else {
                let (u1, u2) = (uniform().max(1e-12), uniform());
                let n = (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();
                errs.push(-2.5 + 1.6 * n);
            }
        }
        let pairs: Vec<(f64, f64)> = errs
            .iter()
            .enumerate()
            .map(|(i, &e)| (e, (i % 21) as f64 * 0.01 - 0.1))
            .collect();
        let m = estimate(&pairs, 700).expect("estimate");
        let (nc, ns) = m.node_error_da.expect("nominal-reference error");
        assert!(nc.abs() < 0.011 && (ns - 0.074).abs() < 0.01, "{m:?}");
        let w = m.node_window_da().unwrap();
        let six = scoring_crate::scoring::fragment_error::NODE_WINDOW_SIGMAS;
        assert!((w - (nc.abs() + six * ns)).abs() < 1e-12, "{w}");
        assert!((m.offset_ppm + 2.5).abs() < 0.15, "{m:?}");
        assert!((m.sigma_ppm - 1.6).abs() < 0.25, "{m:?}");
        assert_eq!(m.psms, 700);
        assert!(m.feature_tolerance_ppm() > 5.0 && m.feature_tolerance_ppm() < 8.0);
    }

    fn queue_of(scores: &[(u32, f32)]) -> crate::TopNQueue {
        let mut q = crate::TopNQueue::new(scores.len() as u32);
        for &(idx, s) in scores {
            q.push(crate::PsmMatch {
                spectrum_idx: 0,
                candidate_idxs: vec![idx],
                charge_used: 2,
                mass_error_ppm: 0.0,
                score: s,
                rank_score: s,
                edge_score: 0,
                activation_method: None,
                features: Default::default(),
                isotope_offset: 0,
                precursor_mz_override: None,
                site_evidence: None,
                tie_key: idx as u64,
            });
        }
        q
    }

    #[test]
    fn a_dia_plan_keeps_the_best_ranks_of_every_spectrum() {
        let queues = vec![
            queue_of(&[(1, 3.0), (2, 9.0), (3, 5.0)]),
            queue_of(&[]),
            queue_of(&[(7, 1.0)]),
        ];
        assert_eq!(
            ranked_candidates(&queues, 1),
            vec![(0, 2, 9.0), (2, 7, 1.0)]
        );
        assert_eq!(
            ranked_candidates(&queues, 2),
            vec![(0, 2, 9.0), (0, 3, 5.0), (2, 7, 1.0)]
        );
        assert_eq!(ranked_candidates(&queues, 10).len(), 4);
    }

    #[test]
    fn a_strided_sample_spans_the_run_and_a_block_sits_in_its_middle() {
        let plan = |strided| MeasurePlan {
            sample_spectra: 4,
            ranks_per_spectrum: 1,
            strided,
        };
        assert_eq!(plan(true).sample_positions(10), vec![0, 2, 5, 7]);
        assert_eq!(plan(false).sample_positions(10), vec![3, 4, 5, 6]);
        assert_eq!(plan(true).sample_positions(3), vec![0, 1, 2]);
        assert_eq!(plan(false).sample_positions(3), vec![0, 1, 2]);
        assert!(plan(true).sample_positions(0).is_empty());
    }

    #[test]
    fn too_few_errors_give_no_estimate() {
        assert_eq!(estimate(&[(1.0, 0.01); 10], 10), None);
    }
}
