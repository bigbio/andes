//! quantms-rescoring–compatible feature contracts and top-N batch helpers.
//!
//! This module is the Rust-side counterpart to
//! [quantms-rescoring](https://github.com/bigbio/quantms-rescoring): feature
//! *names* and *roles* match the Python package (AlphaPeptDeep / DeepLC /
//! spectrum SNR columns) so andes PIN / idparquet writers and a future Rust
//! rescoring path can share one vocabulary with quantms / Percolator.
//!
//! Inference itself stays in [`crate::OnnxSession`]; this module only defines
//! the data contracts and batching discipline (top-N only).

mod features;
mod psm;

pub use features::{
    alphapeptdeep_feature_names, deeplc_feature_names, spectrum_feature_names,
    ALPHAPEPTDEEP_CORE_FEATURES, DEEPLC_CORE_FEATURES, SPECTRUM_CORE_FEATURES,
};
pub use psm::{CandidateKey, FeatureRow, RescoringBatch, TopNPolicy};
