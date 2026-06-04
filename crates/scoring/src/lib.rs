//! Scoring model and ion-type prediction.
//!
//! Contains the parameter model, rank-based node scoring, and fragment ion
//! prediction. Depends only on the `model` crate.

pub mod param_model;
pub mod scoring;

#[cfg(test)]
pub(crate) mod testutil;

// Convenience re-exports.
pub use param_model::{Param, ParamParseError};
pub use scoring::{IonMatchFact, RankScorer, ScoredSpectrum};
