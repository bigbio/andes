//! Configuring a `PreparedSearch` after it is built: scorers per scan, fragment-error model, RT
//! models, candidate filter, MS1 links, and the counters it reports.

use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, RwLock};

use input::Ms1Link;
use model::activation::ActivationMethod;
use model::spectrum::Spectrum;
use scoring_crate::psm::RankScorer;
use smallvec::SmallVec;

use super::dedup::*;
use super::{
    ActivationRoute, CandidateFilter, IpfScoring, PerScanScorers, PreparedSearch, ScorerLoader,
    ScorerSet,
};
use crate::candidates::Candidate;

impl<'a> PreparedSearch<'a> {
    /// Installs the per-form RT predictor of the DIA vote cut.
    pub fn set_form_rt_model(&mut self, model: Arc<crate::dia::FormRtModel>) {
        self.form_rt_model = Some(model);
    }

    pub fn form_rt_model(&self) -> Option<&crate::dia::FormRtModel> {
        self.form_rt_model.as_deref()
    }

    /// Sets (or with `None` clears) the run's RT prior; while set, every DIA chunk index
    /// carries calibrated per-form RTs and its vote cut uses the RT term.
    /// Clearing also drops the cached per-form RT predictions.
    pub fn set_dia_rt_prior(&self, prior: Option<crate::dia::DiaRtPrior>) {
        if prior.is_none() {
            self.form_rt_cache.clear();
        }
        *self.dia_rt_prior.write().expect("RT prior lock") = prior.map(Arc::new);
    }

    pub fn with_per_scan_scorers(
        mut self,
        run_activation: ActivationMethod,
        loader: ScorerLoader,
    ) -> Self {
        self.per_scan = Some(PerScanScorers {
            run_activation,
            loader,
            loaded: RwLock::new(Vec::new()),
        });
        self
    }

    pub fn set_fragment_error(&self, m: Option<scoring_crate::psm::FragmentErrorModel>) {
        self.scorer.set_fragment_error(m);
        if let Some(ps) = &self.per_scan {
            for route in ps.loaded.read().unwrap().iter() {
                if let Some(s) = &route.scorer {
                    s.set_fragment_error(m);
                }
            }
        }
    }

    pub fn per_scan_routes(&self) -> Vec<(ActivationMethod, bool, u64)> {
        self.per_scan.as_ref().map_or_else(Vec::new, |ps| {
            ps.loaded
                .read()
                .unwrap()
                .iter()
                .map(|r| {
                    (
                        r.activation,
                        r.scorer.is_some(),
                        r.scans.load(Ordering::Relaxed),
                    )
                })
                .collect()
        })
    }

    pub(crate) fn scorer_set(&self, spectra: &[Spectrum]) -> ScorerSet<'a> {
        let Some(ps) = &self.per_scan else {
            return ScorerSet {
                primary: self.scorer,
                run_activation: None,
                alts: Vec::new(),
            };
        };
        let mut seen: SmallVec<[(ActivationMethod, u64); 4]> = SmallVec::new();
        for a in spectra.iter().filter_map(|s| s.activation_method) {
            if a == ps.run_activation {
                continue;
            }
            match seen.iter_mut().find(|(m, _)| *m == a) {
                Some((_, n)) => *n += 1,
                None => seen.push((a, 1)),
            }
        }
        for &(a, n) in &seen {
            let known = ps.loaded.read().unwrap().iter().any(|r| r.activation == a);
            if !known {
                let mut loaded = ps.loaded.write().unwrap();
                if !loaded.iter().any(|r| r.activation == a) {
                    let scorer = (ps.loader)(a).map(|s| {
                        s.set_fragment_error(self.scorer.fragment_error());
                        Arc::new(s)
                    });
                    loaded.push(ActivationRoute {
                        activation: a,
                        scorer,
                        scans: AtomicU64::new(0),
                    });
                }
            }
            if let Some(r) = ps.loaded.read().unwrap().iter().find(|r| r.activation == a) {
                r.scans.fetch_add(n, Ordering::Relaxed);
            }
        }
        let alts = ps
            .loaded
            .read()
            .unwrap()
            .iter()
            .filter_map(|r| r.scorer.as_ref().map(|s| (r.activation, Arc::clone(s))))
            .collect();
        ScorerSet {
            primary: self.scorer,
            run_activation: Some(ps.run_activation),
            alts,
        }
    }

    pub fn ipf_scoring<'s>(&self, scorer: &'s RankScorer) -> IpfScoring<'s> {
        assert!(
            std::ptr::eq(scorer, self.scorer),
            "ipf_scoring needs the scorer the search was prepared with"
        );
        let set = self.scorer_set(&[]);
        IpfScoring {
            scorers: ScorerSet {
                primary: scorer,
                run_activation: set.run_activation,
                alts: set.alts,
            },
            fragment_tolerance_da: self.fragment_tolerance_da,
            declared: self.params.aa_set.distinct_mods(),
        }
    }

    pub fn chunk_prep_seconds(&self) -> f64 {
        self.chunk_prep_nanos.load(Ordering::Relaxed) as f64 * 1e-9
    }

    pub fn open_search_stats(&self) -> (usize, usize) {
        (
            self.open_rescore_attempts.load(Ordering::Relaxed),
            self.open_rescore_wins.load(Ordering::Relaxed),
        )
    }
    pub fn sync_materialized_candidates(&mut self) {
        let mut accum = self.mmap_accum.lock().unwrap();
        let materialized = std::mem::take(&mut accum.candidates);
        accum.index_of = Default::default();
        drop(accum);
        self.candidates = materialized;
    }

    /// Number of candidates materialised so far; the global index the next new one gets.
    pub fn materialized_len(&self) -> usize {
        self.mmap_accum
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .candidates
            .len()
    }

    /// Forgets the candidates materialised after the first `len`, and returns their memory.
    /// Only for a calibration search whose queues are dropped: a PSM that still holds one of
    /// these global indices would read another candidate.
    pub fn truncate_materialized(&self, len: usize) {
        let mut accum = self
            .mmap_accum
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        if accum.candidates.len() <= len {
            return;
        }
        accum.candidates.truncate(len);
        accum.candidates.shrink_to_fit();
        accum.index_of.retain(|_, idx| (*idx as usize) < len);
        accum.index_of.shrink_to_fit();
    }

    pub fn with_candidate_filter(mut self, keep: CandidateFilter) -> Self {
        self.target_bare_seqs =
            build_target_bare_seqs_where(&self.mmap_index, self.idx, self.params, &keep);
        self.candidate_filter = Some(keep);
        self
    }

    pub(super) fn keeps(&self, cand: &Candidate) -> bool {
        self.candidate_filter.as_ref().is_none_or(|keep| keep(cand))
    }

    pub fn with_ms1_link(mut self, ms1_link: Option<Ms1Link>) -> Self {
        self.ms1_link = ms1_link;
        self
    }
}
