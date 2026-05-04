//! Input-side readers for MS-GF+ Rust port: spectrum file formats
//! (MGF in Phase 3a; mzML in Phase 3b) and `.fasta` (Phase 4).

pub mod mgf;
pub mod spectrum;

pub use mgf::{MgfParseError, MgfReader};
pub use spectrum::Spectrum;
