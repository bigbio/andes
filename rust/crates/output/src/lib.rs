//! Output writers for MS-GF+ search results.
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
pub use pin::{write_pin, write_pin_to};

pub(crate) mod row_context;
