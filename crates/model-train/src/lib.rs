//! Scoring-model training and the Parquet model store.
pub mod gbdt;
pub mod accumulate;
pub mod counts;
pub mod estimate;
pub mod geometry;
pub mod catalog;
pub mod select;
pub mod labeled;
pub mod store;
pub mod gate;

// Re-export the most commonly used types at the crate root.
pub use store::{ModelStore, RawManifestEntry, protocol_to_experiment_class};
pub use select::{SelectionEntry, SelectionKey, select, parse_experiment_class};

#[derive(thiserror::Error, Debug)]
pub enum TrainError {
    #[error("io: {0}")]
    Io(#[from] std::io::Error),
    #[error("parquet: {0}")]
    Parquet(String),
    #[error("model not found for selection key: {0}")]
    NoModel(String),
    /// A GBDT trainer's hard quality gate failed (finding 3.6): too few rows,
    /// a single-class classification target, sub-threshold held-out
    /// AUC/Pearson/R², or an empty tree ensemble. The caller decides whether to
    /// hard-error or fall back, but the trainer never silently returns a
    /// non-deployable model.
    #[error("gbdt quality gate failed: {0}")]
    QualityGate(String),
    #[error("{0}")]
    Other(String),
}
