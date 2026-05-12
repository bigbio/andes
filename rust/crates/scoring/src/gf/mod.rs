//! Generating-function (GF) DP for SpecEValue computation. Mirrors Java
//! `ScoreBound`, `ScoreDist`, `GeneratingFunction`, and `PrimitiveAaGraph`.

pub mod score_dist;
pub mod generating_function;
pub mod primitive_graph;
pub mod group;

pub use score_dist::{ScoreBound, ScoreDist};
pub use generating_function::{GeneratingFunction, GfError};
pub use primitive_graph::PrimitiveAaGraph;
pub use group::GeneratingFunctionGroup;
