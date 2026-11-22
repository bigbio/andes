//! Feature computation and learned scoring of candidate PSMs: fragment and ion features,
//! peak features, the retention-time model and the gradient-boosted tree evaluator. It sits
//! between candidate matching in `search` and the rescoring step, producing the features each
//! PSM carries.

pub mod frag_features;
pub mod gbdt_eval;
pub mod ion_features;
pub mod mod_site_features;
pub mod param_model;
pub mod peak_features;
pub mod rt_model;
pub mod scoring;

#[cfg(test)]
pub(crate) mod testutil;

pub use param_model::Param;
pub use scoring::{IonMatchFact, RankScorer, ScoredSpectrum};
