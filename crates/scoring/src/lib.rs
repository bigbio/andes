//! Scoring sub-system for MS-GF+ Rust port.
//!
//! Contains the parameter model, rank-based scoring, fragment ion
//! prediction, and the generating-function DP for SpecEValue.
//! Depends only on the `model` crate.

pub mod gf;
pub mod param_model;
pub mod scoring;

#[cfg(test)]
pub(crate) mod testutil;

// Convenience re-exports.
pub use param_model::{Param, ParamParseError};
pub use scoring::{RankScorer, ScoredSpectrum};
