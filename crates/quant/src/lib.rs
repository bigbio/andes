//! Quantification: isobaric reporter ions (TMT / iTRAQ) and label-free MS1
//! precursor features.
//!
//! The crate is pure algorithms over plain data: peak lists, retention times
//! and peptide formulas go in, intensities and scores come out. It knows
//! nothing about the search (PSMs, candidates, FDR) — the `andes` binary
//! builds the targets and the `output` crate writes the tables.
//!
//! * [`isobaric`]: plex tables, reporter extraction, impurity correction (NNLS).
//! * [`purity`]: precursor isolation purity from the preceding MS1 scan.
//! * [`ms1_index`]: a run's MS1 scans indexed for chromatogram extraction.
//! * [`xic`]: extracted-ion chromatograms, apex/boundary detection, integration.
//! * [`formula`]: peptide elemental formula → theoretical isotope envelope.
//! * [`lfq`]: per-target feature quantification and the decoy-controlled
//!   feature q-value.

pub mod formula;
pub mod isobaric;
pub mod lfq;
pub mod ms1_index;
pub mod nnls;
pub mod purity;
pub mod tdc;
pub mod xic;

pub use isobaric::{extract_reporters, Channel, CorrectionMatrix, Plex};
pub use lfq::{FeatureQuant, LfqParams, LfqTarget};
pub use ms1_index::{Ms1RunIndex, Ms1Scan};
pub use purity::precursor_purity;
