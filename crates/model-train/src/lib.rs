//! Training of the learned models: accumulating labelled PSMs, selecting and gating training
//! data, fitting gradient-boosted trees and storing the resulting model artifacts. It runs
//! offline, outside a search, and produces the models that `scoring` and `search` load.

pub mod accumulate;
pub mod catalog;
pub mod counts;
pub mod estimate;
pub mod gate;
pub mod gbdt;
pub mod geometry;
pub mod labeled;
pub mod select;
pub mod store;

pub use select::{
    parse_experiment_class, select, select_nearest, select_nearest_for, SelectionEntry,
    SelectionKey,
};
pub use store::{protocol_to_experiment_class, ModelStore, RawManifestEntry};

#[derive(thiserror::Error, Debug)]
pub enum TrainError {
    #[error("io: {0}")]
    Io(#[from] std::io::Error),
    #[error("parquet: {0}")]
    Parquet(String),
    #[error("model not found for selection key: {0}")]
    NoModel(String),
    #[error("gbdt quality gate failed: {0}")]
    QualityGate(String),
    #[error("{0}")]
    Other(String),
}
