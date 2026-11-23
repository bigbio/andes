//! Precursor masses: matching a candidate to the precursor or isolation window, the
//! monoisotopic correction, the per-file calibration and the pre-pass that learns it.

pub mod calibration;
pub mod matching;
pub mod monoisotopic;
pub mod prepass;
