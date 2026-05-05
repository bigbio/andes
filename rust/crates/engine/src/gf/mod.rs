//! Phase 6 generating-function (GF) DP for SpecEValue computation.
//! Highest-correctness-risk phase per the parent design.
//!
//! Phase 6 Task 1: port Java's ScoreBound + ScoreDist (pure data
//! wrappers with indexed access; no algorithm logic).
//!
//! Tasks 2+: GF DP itself, spectral-probability computation, etc.

pub mod score_dist;
pub mod generating_function;
pub mod primitive_graph;

pub use score_dist::{ScoreBound, ScoreDist};
pub use generating_function::{GeneratingFunction, GfError};
pub use primitive_graph::PrimitiveAaGraph;
