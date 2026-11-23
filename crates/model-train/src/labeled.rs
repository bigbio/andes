use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use input::FastaReader;
use model::fdr::{ranked_q, Counts, Eval};
use model::peptide::Peptide;
use model::spectrum::Spectrum;
use model::stats::nan_low_f32;
use scoring_crate::scoring::rank_scorer::RankScorer;
use search::decoy::{DecoyStrategy, Proteases, DEFAULT_DECOY_SEED};
use search::psm::PsmMatch;
use search::{match_spectra, SearchIndex, SearchParams};

use crate::TrainError;

pub const BOOTSTRAP_DECOY_PREFIX: &str = "XXX_";

const FRAGMENT_TOL_DA: f64 = 0.5;

#[derive(Debug, Clone)]
pub struct LabeledMatch {
    pub spectrum_index: usize,
    pub peptide: Peptide,
    pub charge: u8,
    pub confidence: f64,
}

pub fn bootstrap_labels(
    spectra: &[Spectrum],
    database: &Path,
    seed: &RankScorer,
    search_params: &SearchParams,
    train_fdr: f64,
    cache_dir: &Path,
) -> Result<Vec<LabeledMatch>, TrainError> {
    let file = File::open(database)?;
    let target_db = FastaReader::load_all(BufReader::new(file))
        .map_err(|e| TrainError::Other(format!("FASTA parse error: {e}")))?;
    let idx = SearchIndex::from_target_db_for(
        &target_db,
        BOOTSTRAP_DECOY_PREFIX,
        None,
        DecoyStrategy::default(),
        DEFAULT_DECOY_SEED,
        &Proteases::of(search_params),
    );

    let (queues, candidates) = match_spectra(
        spectra,
        &idx,
        search_params,
        seed,
        FRAGMENT_TOL_DA,
        BOOTSTRAP_DECOY_PREFIX,
        cache_dir,
    )?;

    struct BestPsm {
        spectrum_index: usize,
        peptide: Peptide,
        charge: u8,
        score: f32,
        is_decoy: bool,
    }

    let mut best_psms: Vec<BestPsm> = Vec::new();
    for (spec_idx, queue) in queues.iter().enumerate() {
        if queue.is_empty() {
            continue;
        }
        if let Some(psm) = queue.peek_top() {
            let cand_idx = psm.primary_candidate_idx() as usize;
            let cand = &candidates[cand_idx];
            best_psms.push(BestPsm {
                spectrum_index: spec_idx,
                peptide: cand.peptide.clone(),
                charge: psm.charge_used,
                score: tdc_score(psm),
                is_decoy: cand.is_decoy,
            });
        }
    }

    let mono_q = tdc_q_values(&mut best_psms, |p| p.score, |p| p.is_decoy);

    let mut labels: Vec<LabeledMatch> = Vec::new();
    for (psm, q) in best_psms.into_iter().zip(mono_q.into_iter()) {
        if !psm.is_decoy && q <= train_fdr {
            labels.push(LabeledMatch {
                spectrum_index: psm.spectrum_index,
                peptide: psm.peptide,
                charge: psm.charge,
                confidence: q,
            });
        }
    }

    Ok(labels)
}

pub(crate) fn tdc_score(psm: &PsmMatch) -> f32 {
    psm.score
}

pub(crate) fn tdc_q_values<T>(
    psms: &mut [T],
    score_of: impl Fn(&T) -> f32,
    is_decoy: impl Fn(&T) -> bool,
) -> Vec<f64> {
    let key = |p: &T| nan_low_f32(score_of(p));
    psms.sort_by(|a, b| {
        key(b)
            .partial_cmp(&key(a))
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| is_decoy(a).cmp(&is_decoy(b)))
    });
    ranked_q::<Counts>(
        psms.len(),
        |c, rank| c.add(is_decoy(&psms[rank])),
        |start, rank| scores_tie(score_of(&psms[rank]), score_of(&psms[start])),
        Eval::PerRowGroupMax,
        Counts::ratio,
        1.0,
    )
}

fn scores_tie(a: f32, b: f32) -> bool {
    if a.is_nan() && b.is_nan() {
        return true;
    }
    a.to_bits() == b.to_bits()
}

#[cfg(test)]
mod tests {
    use super::*;

    struct P {
        rank_score: f32,
        is_decoy: bool,
    }

    fn fdr(mut psms: Vec<P>, fdr: f64) -> (Vec<f64>, usize) {
        let q = tdc_q_values(&mut psms, |p| p.rank_score, |p| p.is_decoy);
        let accepted = psms
            .iter()
            .zip(q.iter())
            .filter(|(p, &qi)| !p.is_decoy && qi <= fdr)
            .count();
        (q, accepted)
    }

    #[test]
    fn tied_target_and_decoy_bucket_is_counted_conservatively() {
        let psms = vec![
            P {
                rank_score: 10.0,
                is_decoy: false,
            },
            P {
                rank_score: 5.0,
                is_decoy: false,
            },
            P {
                rank_score: 5.0,
                is_decoy: true,
            },
        ];
        let (q, accepted) = fdr(psms, 0.01);
        assert_eq!(
            accepted, 1,
            "tied target must not be accepted when its bucket has a decoy"
        );
        assert_eq!(q[1], q[2], "tied PSMs must share a q-value");
        assert!(q[1] > 0.01, "tied bucket q must reflect the decoy");
    }

    #[test]
    fn untied_scores_are_unaffected() {
        let psms = vec![
            P {
                rank_score: 10.0,
                is_decoy: false,
            },
            P {
                rank_score: 9.0,
                is_decoy: false,
            },
            P {
                rank_score: 8.0,
                is_decoy: true,
            },
        ];
        let (q, accepted) = fdr(psms, 0.01);
        assert_eq!(accepted, 2, "both distinct-score targets pass at 1% FDR");
        assert!(q[0] <= 0.01 && q[1] <= 0.01);
    }

    fn queued(score: f32, edge: i32) -> PsmMatch {
        PsmMatch {
            spectrum_idx: 0,
            candidate_idxs: vec![0].into(),
            charge_used: 2,
            mass_error_ppm: 0.0,
            score,
            rank_score: score + edge as f32,
            edge_score: edge,
            activation_method: None,
            features: Default::default(),
            isotope_offset: 0,
            precursor_mz_override: None,
            site_evidence: None,
            tie_key: 0,
        }
    }

    #[test]
    fn edge_inflated_decoys_do_not_outrank_targets() {
        let mut hits: Vec<(PsmMatch, bool)> = Vec::new();
        for i in 0..200 {
            hits.push((queued(40.0 + (i % 50) as f32, 0), false));
        }
        for i in 0..3 {
            hits.push((queued(15.0 + i as f32, 300), true));
        }
        for i in 0..60 {
            hits.push((queued(10.0 + (i % 20) as f32 * 0.5, 0), true));
        }
        let accepted = |score_of: fn(&PsmMatch) -> f32| {
            let mut v = hits.clone();
            let q = tdc_q_values(&mut v, |h| score_of(&h.0), |h| h.1);
            v.iter()
                .zip(q.iter())
                .filter(|(h, &qi)| !h.1 && qi <= 0.01)
                .count()
        };
        assert_eq!(
            accepted(|p| p.rank_score),
            0,
            "the edge key buries the targets"
        );
        assert_eq!(accepted(tdc_score), 200);
    }
}
