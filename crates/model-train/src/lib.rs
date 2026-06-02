//! Scoring-model training and the Parquet model store.
pub mod counts;
pub mod estimate;
pub mod catalog;
pub mod select;
pub mod labeled;
pub mod store;

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
    #[error("{0}")]
    Other(String),
}
