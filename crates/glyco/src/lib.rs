//! Glycopeptide identification: glycan compositions and masses, oxonium-ion evidence, sequon
//! handling and the glycan-first matching that pairs a peptide backbone with a glycan. It sits
//! between spectrum search and rescoring for glyco runs, supplying the glyco PSM candidates and
//! their features.

pub mod backbone;
pub mod psm;
pub mod hybrid;
pub mod oxonium;
pub mod rt;
pub mod sequon;
pub mod glycan;
