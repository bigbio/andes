//! Feature computation and learned scoring of candidate PSMs: fragment and ion features,
//! peak features, the retention-time model and the gradient-boosted tree evaluator. It sits
//! between candidate matching in `search` and the rescoring step, producing the features each
//! PSM carries.

pub mod features;
pub mod gbdt;
pub mod params;
pub mod psm;
pub mod rt;

#[cfg(test)]
pub(crate) mod testutil;

pub use params::Param;
pub use psm::{IonMatchFact, RankScorer, ScoredSpectrum};
