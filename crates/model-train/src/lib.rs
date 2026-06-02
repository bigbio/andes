//! Scoring-model training and the Parquet model store.
pub mod counts;
pub mod estimate;
pub mod catalog;
pub mod select;
pub mod labeled;
pub mod store;

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
