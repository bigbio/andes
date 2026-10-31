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
//! * [`unmix`]: co-isolation interference of MS2 reporter ions, modelled per
//!   run from the isolation-window MS1 contents, and its removal.
//! * [`ms1_index`]: a run's MS1 scans indexed for chromatogram extraction.
//! * [`xic`]: extracted-ion chromatograms, apex/boundary detection, integration.
//! * [`formula`]: peptide elemental formula → theoretical isotope envelope.
//! * [`lfq`]: per-target feature quantification and the decoy-controlled
//!   feature q-value.
//! * [`library`]: precursor-level match-between-runs library entry gate.
//! * [`align`]: joint retention-time alignment across runs, without a
//!   reference run, and the predicted RT and its uncertainty per run.
//! * [`decoys`]: mass-defect-preserving and RT-shifted decoy hypotheses for
//!   the match-between-runs transfer FDR.
//! * [`normalize`]: run normalisation of MS1 quantities across runs (a level
//!   per run and an RT-local term).
//! * [`evidence`]: the likelihood template fit that measures one (precursor,
//!   run) hypothesis: amplitude, quantity and the evidence statistics.
//! * [`ms2`]: empirical fragment spectra of library entries and the MS2
//!   evidence of a transfer probe in the receiving run's MS2 scans.

pub mod align;
pub mod decoys;
pub mod evidence;
pub mod formula;
pub mod gaussian;
pub mod isobaric;
pub mod lfq;
pub mod library;
pub mod mobility;
pub mod ms1_index;
pub mod ms2;
pub mod nnls;
pub mod normalize;
pub mod purity;
pub mod tdc;
pub mod unmix;
pub mod xic;

pub use isobaric::{extract_reporters, Channel, CorrectionMatrix, Plex};
pub use lfq::{FeatureQuant, LfqParams, LfqTarget};
pub use ms1_index::{MobilityGate, Ms1RunIndex, Ms1Scan};
pub use purity::precursor_purity;
