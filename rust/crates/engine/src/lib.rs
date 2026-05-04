//! Domain model for MS-GF+ Rust port.
//!
//! Phase 1 milestone: amino acids, modifications, peptides, enzymes,
//! tolerances. Pure CPU + types. No I/O. See
//! `docs/superpowers/2026-05-03-phase1-engine-domain-model-design.md`.

pub mod mass;
pub mod tolerance;
pub mod enzyme;
pub mod modification;
pub mod amino_acid;
pub mod peptide;
pub mod aa_set;
