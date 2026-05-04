//! Input-side readers for MS-GF+ Rust port: spectrum file formats
//! (MGF in Phase 3a; mzML in Phase 3b) and `.fasta` (Phase 4a).

pub mod mgf;

pub use engine::Spectrum;
pub use mgf::{MgfParseError, MgfReader};
