//! Rust-only GBDT training subsystem: labels → dataset → trainer → GbdtPeakModel.
pub mod dataset;
pub mod frag_dataset;
pub mod isotonic;
pub mod labels;
pub mod train;
pub mod tree;
