use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use input::FastaReader;
use model::spectrum::Spectrum;
use scoring_crate::scoring::rank_scorer::RankScorer;
use search::decoy::{DecoyStrategy, Proteases, DEFAULT_DECOY_SEED};
use search::{match_spectra, SearchIndex, SearchParams};

use crate::labeled::BOOTSTRAP_DECOY_PREFIX;
use crate::TrainError;

const FRAGMENT_TOL_DA: f64 = 0.5;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct YieldDelta {
    pub current_count: usize,
    pub candidate_count: usize,
}

impl YieldDelta {
    pub fn is_accepted(&self) -> bool {
        self.candidate_count >= self.current_count
    }
}

pub fn evaluate_candidate(
    spectra: &[Spectrum],
    database: &Path,
    current: &RankScorer,
    candidate: &RankScorer,
    search_params: &SearchParams,
    fdr: f64,
    cache_dir: &Path,
) -> Result<YieldDelta, TrainError> {
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

    let current_count = count_target_psms(spectra, &idx, search_params, current, fdr, cache_dir)?;
    let candidate_count =
        count_target_psms(spectra, &idx, search_params, candidate, fdr, cache_dir)?;

    Ok(YieldDelta {
        current_count,
        candidate_count,
    })
}

pub(crate) fn count_target_psms(
    spectra: &[Spectrum],
    idx: &SearchIndex,
    search_params: &SearchParams,
    scorer: &RankScorer,
    fdr: f64,
    cache_dir: &Path,
) -> Result<usize, TrainError> {
    let (queues, candidates) = match_spectra(
        spectra,
        idx,
        search_params,
        scorer,
        FRAGMENT_TOL_DA,
        BOOTSTRAP_DECOY_PREFIX,
        cache_dir,
    )?;

    struct BestPsm {
        score: f32,
        is_decoy: bool,
    }

    let mut best_psms: Vec<BestPsm> = Vec::new();
    for queue in queues.iter() {
        if queue.is_empty() {
            continue;
        }
        if let Some(psm) = queue.peek_top() {
            let cand_idx = psm.primary_candidate_idx() as usize;
            let cand = &candidates[cand_idx];
            best_psms.push(BestPsm {
                score: crate::labeled::tdc_score(psm),
                is_decoy: cand.is_decoy,
            });
        }
    }

    let mono_q = crate::labeled::tdc_q_values(&mut best_psms, |p| p.score, |p| p.is_decoy);

    let count = best_psms
        .iter()
        .zip(mono_q.iter())
        .filter(|(psm, &q)| !psm.is_decoy && q <= fdr)
        .count();

    Ok(count)
}
