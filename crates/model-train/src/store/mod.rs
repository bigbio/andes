//! Parquet model store: schema definitions, reading, and writing.
pub mod schema;
pub mod write;
pub mod read;

pub use write::write_models;
pub use read::ModelStore;
