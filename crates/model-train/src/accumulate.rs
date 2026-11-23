use model::peptide::Peptide;
use model::spectrum::Spectrum;
use scoring_crate::params::IonType;
use scoring_crate::psm::psm_edge_existence_facts;
use scoring_crate::psm::rank::RankScorer;
use scoring_crate::psm::spectrum::ScoredSpectrum;

use crate::counts::CountStats;

static DENSE_NOISE: std::sync::OnceLock<Option<usize>> = std::sync::OnceLock::new();

pub fn init_dense_noise(n: Option<usize>) {
    let _ = DENSE_NOISE.set(n.filter(|&n| n > 0));
}

#[inline]
fn dense_noise_sampling() -> Option<usize> {
    DENSE_NOISE.get().copied().flatten()
}

pub struct StatsAccumulator<'a> {
    scorer: &'a RankScorer,
}

impl<'a> StatsAccumulator<'a> {
    pub fn new(scorer: &'a RankScorer) -> Self {
        Self { scorer }
    }

    pub fn accumulate(
        &self,
        stats: &mut CountStats,
        spec: &Spectrum,
        peptide: &Peptide,
        charge: u8,
    ) {
        let scored_spec = ScoredSpectrum::new(spec, self.scorer, charge);

        let max_rank = self.scorer.max_rank();

        let facts = scored_spec.ion_match_facts(peptide, self.scorer);

        for fact in &facts {
            let rank_idx = match fact.rank {
                Some(r) => {
                    let clamped = r.min(max_rank).max(1);
                    clamped - 1
                }
                None => max_rank,
            };
            stats.bump_rank(fact.partition, fact.ion_type, rank_idx);

            if let (Some(_), Some(bin)) = (fact.rank, fact.error_bin) {
                stats.bump_error(fact.partition, bin);
            }
        }

        let noise_facts = match dense_noise_sampling() {
            Some(n) => scored_spec.dense_noise_facts(peptide, self.scorer, n),
            None => scored_spec.noise_match_facts(peptide, self.scorer),
        };
        for (partition, rank, error_bin) in noise_facts {
            let rank_idx = match rank {
                Some(r) => r.min(max_rank).max(1) - 1,
                None => max_rank,
            };
            stats.bump_rank(partition, IonType::Noise, rank_idx);

            if let (Some(_), Some(bin)) = (rank, error_bin) {
                stats.bump_noise_error(partition, bin);
            }
        }

        for (partition, idx) in psm_edge_existence_facts(&scored_spec, peptide, self.scorer, charge)
        {
            stats.bump_existence(partition, idx);
        }

        stats.bump_charge(charge as i32);
    }
}

pub fn merge(parts: Vec<CountStats>) -> CountStats {
    parts.into_iter().fold(CountStats::new(), |mut acc, part| {
        acc.add(&part);
        acc
    })
}
