//! Output writers for MS-GF+ search results.
//!
//! # Java divergences (Phase 7 MVP)
//!
//! * **FragMethod**: Java emits the string from `ActivationMethod`'s static
//!   table (e.g. `"HCD"`, `"CID"`). Rust uses `ActivationMethod::name()` which
//!   produces the same strings for the five canonical variants. Unknown
//!   activation is written as `"UNKNOWN"` (Java writes `"Unknown"`).
//!
//! * **IsotopeError**: Phase 4e's precursor-matching loop tries multiple
//!   isotope offsets but does not record *which* offset produced the match.
//!   The TSV column is always written as `0`. Fix in a later phase once the
//!   winning isotope offset is threaded into `PsmMatch`.
//!
//! * **Decoy filtering**: Java skips PSMs where the protein string is empty
//!   (all-decoy with no target alternates). This writer emits decoy PSMs with
//!   the decoy prefix preserved in the Protein column; downstream Percolator
//!   handles decoy labelling. Intentional simplification.
//!
//! * **QValue / PepQValue**: Not emitted. Phase 7 MVP omits TDA columns.
//!   Will be added in Task 6 or a later phase.

pub mod tsv;
pub use tsv::write_tsv;
