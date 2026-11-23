//! Quantification: chromatogram extraction, retention-time alignment, label-free and
//! isobaric (TMT/iTRAQ) intensities, match-between-runs transfers and normalisation. It runs
//! after identification and turns the confident PSMs into per-run precursor intensities.

pub mod align;
pub mod claims;
pub mod decoys;
pub mod evidence;
pub mod formula;
pub mod gaussian;
pub mod isobaric;
pub mod label;
pub mod lfq;
pub mod library;
pub mod mobility;
pub mod ms1;
pub mod ms2;
pub mod nnls;
pub mod normalize;
pub mod purity;
pub mod shared;
pub mod tdc;
pub mod transfer;
pub mod unmix;
pub mod xic;

pub use isobaric::{extract_reporters, Channel, CorrectionMatrix, Plex};
pub use lfq::{FeatureQuant, LfqParams, LfqTarget};
pub use ms1::{MobilityGate, Ms1RunIndex, Ms1Scan};
pub use purity::precursor_purity;
