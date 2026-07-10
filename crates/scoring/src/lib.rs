//! Scoring model and ion-type prediction.
//!
//! Contains the parameter model, rank-based node scoring, and fragment ion
//! prediction. Depends only on the `model` crate.

pub mod frag_features;
pub mod ion_features;
pub mod gbdt_eval;
pub mod intensity_model;
pub mod mod_site_features;
pub mod param_model;
pub mod peak_features;
pub mod rt_model;
pub mod scoring;

#[cfg(test)]
pub(crate) mod testutil;

pub use intensity_model::{IntensityIonType, IntensityModel, IntensityModelError};
pub use param_model::Param;
pub use scoring::{IonMatchFact, RankScorer, ScoredSpectrum};
