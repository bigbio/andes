//! Output writers: Percolator PIN, TSV.
//!
//! # Known column behaviors
//!
//! * **FragMethod**: emitted via `ActivationMethod::name()` (e.g. `"HCD"`,
//!   `"CID"`). Unknown activation is written as `"UNKNOWN"`.
//!
//! * **IsotopeError**: the precursor-matching loop tries multiple isotope
//!   offsets but does not record *which* offset produced the match. The TSV
//!   column is always written as `0`. Will be fixed once the winning
//!   isotope offset is threaded into `PsmMatch`.
//!
//! * **Decoy filtering**: this writer emits decoy PSMs with the decoy
//!   prefix preserved in the Protein column; downstream Percolator handles
//!   decoy labelling.
//!
//! * **QValue / PepQValue**: Not emitted; TDA columns are not currently
//!   produced.

pub mod tsv;
pub use tsv::{write_tsv, write_tsv_to};

pub mod pin;
pub use pin::{format_spec_id, write_pin, write_pin_to};

pub mod stats;
pub use stats::{write_statistics_log, RunStatistics};

pub mod qpx;
pub use qpx::write_qpx;

pub mod rt_wiring;
pub use rt_wiring::populate_rt_features;

pub mod percolator;
pub use percolator::{
    parse_psm_results, resolve_backend, run_percolator, PercolatorBackend, PercolatorError,
    PercolatorPsm, DEFAULT_PERCOLATOR_IMAGE,
};

pub(crate) mod row_context;
pub(crate) mod percolator_enz;

pub mod glyco_pin;
pub use glyco_pin::write_glyco_pin;
