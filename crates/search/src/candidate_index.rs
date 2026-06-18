//! Out-of-core candidate index: binary flat-file of mass-sorted `IndexRecord`s.
//!
//! Each `IndexRecord` encodes a base peptide (no variable mods) compactly in 20
//! bytes using manual little-endian serialisation (no repr(C), no padding).
//!
//! ## Layout (little-endian)
//! ```text
//! 0..8   mass_milli: u64   — neutral monoisotopic mass × 1000, rounded
//! 8..12  protein_index: u32
//! 12..16 start_offset: u32
//! 16..18 length: u16
//! 18..20 flags: u16        — bit0=is_decoy, bit1=is_protein_n_term, bit2=is_protein_c_term
//! ```
//!
//! `build_base_peptide_index()` collects all base peptides from a `SearchIndex`,
//! sorts them by `mass_milli`, and writes the packed binary file. The caller then
//! uses it as an mmapped / seekable lookup for precursor-mass-windowed retrieval.

use std::io::{self, Write};

use crate::candidate_gen::enumerate_candidates;
use crate::search_index::SearchIndex;
use crate::search_params::SearchParams;

/// Serialised byte size of one `IndexRecord`.
pub const INDEX_RECORD_SIZE: usize = 20;

/// Flags bits inside `IndexRecord::flags`.
pub mod flags {
    pub const IS_DECOY: u16 = 1 << 0;
    pub const IS_PROTEIN_N_TERM: u16 = 1 << 1;
    pub const IS_PROTEIN_C_TERM: u16 = 1 << 2;
}

/// A single base-peptide entry in the flat binary index.
///
/// 20 bytes when serialised; sorted by `mass_milli` in the file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IndexRecord {
    /// Neutral monoisotopic mass × 1000, rounded (u64 avoids floating-point
    /// comparison issues in sort and range queries).
    pub mass_milli: u64,
    /// Index into the combined target+decoy `ProteinDb`.
    pub protein_index: u32,
    /// Byte offset of the peptide's first residue within the protein sequence.
    pub start_offset: u32,
    /// Peptide length in residues.
    pub length: u16,
    /// Packed flags: see [`flags`] constants.
    pub flags: u16,
}

impl IndexRecord {
    /// Is this a decoy peptide?
    pub fn is_decoy(&self) -> bool {
        self.flags & flags::IS_DECOY != 0
    }

    /// Serialise to 20 little-endian bytes.
    pub fn to_bytes(&self) -> [u8; INDEX_RECORD_SIZE] {
        let mut buf = [0u8; INDEX_RECORD_SIZE];
        buf[0..8].copy_from_slice(&self.mass_milli.to_le_bytes());
        buf[8..12].copy_from_slice(&self.protein_index.to_le_bytes());
        buf[12..16].copy_from_slice(&self.start_offset.to_le_bytes());
        buf[16..18].copy_from_slice(&self.length.to_le_bytes());
        buf[18..20].copy_from_slice(&self.flags.to_le_bytes());
        buf
    }

    /// Deserialise from 20 little-endian bytes.
    pub fn from_bytes(buf: &[u8; INDEX_RECORD_SIZE]) -> Self {
        let mass_milli = u64::from_le_bytes(buf[0..8].try_into().unwrap());
        let protein_index = u32::from_le_bytes(buf[8..12].try_into().unwrap());
        let start_offset = u32::from_le_bytes(buf[12..16].try_into().unwrap());
        let length = u16::from_le_bytes(buf[16..18].try_into().unwrap());
        let flags = u16::from_le_bytes(buf[18..20].try_into().unwrap());
        Self { mass_milli, protein_index, start_offset, length, flags }
    }
}

/// Enumerate all base peptides (no variable mods) from `idx` matching `params`.
///
/// Clones `params` and forces `max_variable_mods_per_peptide = 0` so only the
/// unmodified backbone is produced. Fixed mods in `aa_set` remain untouched —
/// they fold into per-residue mass exactly as in the normal search path.
///
/// The `decoy_prefix` is forwarded verbatim to [`enumerate_candidates`] for
/// decoy-membership detection (e.g. `"XXX"` when the index was built with that
/// prefix; `decoy_accession_needle` appends the `_` internally).
///
/// Collects into a `Vec` to satisfy lifetime requirements before converting
/// to `IndexRecord`s.
pub fn base_peptide_records(
    idx: &SearchIndex,
    params: &SearchParams,
    decoy_prefix: &str,
) -> Vec<IndexRecord> {
    // Build base params: zero variable mods, keep everything else.
    let base_params = SearchParams {
        max_variable_mods_per_peptide: 0,
        ..params.clone()
    };

    // Collect first (iterator borrows idx and base_params with the same lifetime).
    let candidates: Vec<_> = enumerate_candidates(idx, &base_params, decoy_prefix).collect();

    candidates
        .into_iter()
        .map(|c| {
            let mass_milli = (c.peptide.mass() * 1000.0).round() as u64;
            let mut f: u16 = 0;
            if c.is_decoy {
                f |= flags::IS_DECOY;
            }
            if c.is_protein_n_term {
                f |= flags::IS_PROTEIN_N_TERM;
            }
            if c.is_protein_c_term {
                f |= flags::IS_PROTEIN_C_TERM;
            }
            IndexRecord {
                mass_milli,
                protein_index: c.protein_index as u32,
                start_offset: c.start_offset_in_protein as u32,
                length: c.peptide.residues.len() as u16,
                flags: f,
            }
        })
        .collect()
}

/// Build a mass-sorted binary index of base peptides and write it to `out`.
///
/// Steps:
/// 1. Enumerate all base peptides (no variable mods) via [`base_peptide_records`].
/// 2. Sort by `mass_milli` ascending.
/// 3. Write each record as 20 packed little-endian bytes.
///
/// The `decoy_prefix` must match the prefix used when building the
/// `SearchIndex` (e.g. `"XXX"` — the `_` delimiter is appended internally by
/// `decoy_accession_needle`).
///
/// Returns the number of records written.
pub fn build_base_peptide_index<W: Write>(
    idx: &SearchIndex,
    params: &SearchParams,
    decoy_prefix: &str,
    out: &mut W,
) -> io::Result<usize> {
    let mut records = base_peptide_records(idx, params, decoy_prefix);
    records.sort_by_key(|r| r.mass_milli);

    let n = records.len();
    for rec in &records {
        out.write_all(&rec.to_bytes())?;
    }
    Ok(n)
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::aa_set::AminoAcidSetBuilder;
    use model::protein::{Protein, ProteinDb};
    use crate::search_index::SearchIndex;
    use crate::search_params::SearchParams;

    /// Read all records back from a binary blob.
    fn read_records(data: &[u8]) -> Vec<IndexRecord> {
        assert_eq!(
            data.len() % INDEX_RECORD_SIZE,
            0,
            "binary blob is not a multiple of {INDEX_RECORD_SIZE} bytes"
        );
        data.chunks_exact(INDEX_RECORD_SIZE)
            .map(|chunk| {
                let arr: &[u8; INDEX_RECORD_SIZE] = chunk.try_into().unwrap();
                IndexRecord::from_bytes(arr)
            })
            .collect()
    }

    /// Small toy database: MKWVRPK contains tryptic peptides WVR (len=3),
    /// MK (len=2), PK (len=2), plus missed-cleavage combinations. With
    /// min_length=3, WVR is the shortest included peptide.
    ///
    /// `SearchIndex::from_target_db` with prefix `"XXX"` builds decoys as
    /// `"XXX_<orig>"`, so we pass `"XXX"` to enumerate_candidates (it
    /// appends `"_"` internally via `decoy_accession_needle`).
    fn make_toy_index() -> SearchIndex {
        let target = ProteinDb {
            proteins: vec![Protein {
                accession: "P1".into(),
                description: "toy protein".into(),
                sequence: b"MKWVRPK".to_vec(),
            }],
        };
        SearchIndex::from_target_db(&target, "XXX")
    }

    #[test]
    fn builder_writes_mass_sorted_base_peptides() {
        let idx = make_toy_index();
        let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
        let mut params = SearchParams::default_tryptic(aa_set);
        // Lower min_length so WVR (len=3) is included alongside longer peptides.
        params.min_length = 3;

        let mut buf: Vec<u8> = Vec::new();
        let n = build_base_peptide_index(&idx, &params, "XXX", &mut buf)
            .expect("write failed");

        assert!(n > 0, "expected at least one record, got 0");
        assert_eq!(buf.len(), n * INDEX_RECORD_SIZE, "buffer size mismatch");

        let records = read_records(&buf);
        assert_eq!(records.len(), n);

        // Records must be sorted by mass_milli (ascending).
        for window in records.windows(2) {
            assert!(
                window[0].mass_milli <= window[1].mass_milli,
                "records not sorted: {} > {}",
                window[0].mass_milli,
                window[1].mass_milli
            );
        }

        // At least one non-decoy record must exist (targets are generated).
        assert!(
            records.iter().any(|r| !r.is_decoy()),
            "no target records found"
        );

        // At least one decoy record must exist (reversed decoys are generated by default).
        assert!(
            records.iter().any(|r| r.is_decoy()),
            "no decoy records found"
        );

        // There should be a record with length 3 (WVR from MKWVRPK).
        assert!(
            records.iter().any(|r| r.length == 3),
            "expected at least one length-3 record (WVR from MKWVRPK)"
        );
    }

    #[test]
    fn index_record_round_trips_through_bytes() {
        let rec = IndexRecord {
            mass_milli: 1_234_567_890,
            protein_index: 42,
            start_offset: 7,
            length: 9,
            flags: flags::IS_DECOY | flags::IS_PROTEIN_C_TERM,
        };
        let bytes = rec.to_bytes();
        assert_eq!(bytes.len(), INDEX_RECORD_SIZE);
        let back = IndexRecord::from_bytes(&bytes);
        assert_eq!(back, rec);
        assert!(back.is_decoy());
    }

    #[test]
    fn base_params_forces_zero_variable_mods() {
        // Verify that base_peptide_records does not blow up with max_variable_mods=0
        // and that it actually produces records.
        let idx = make_toy_index();
        let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
        let mut params = SearchParams::default_tryptic(aa_set);
        params.min_length = 3;
        // Even when caller sets 3 mods, base index should produce only base peptides.
        params.max_variable_mods_per_peptide = 3;

        let records = base_peptide_records(&idx, &params, "XXX");
        assert!(!records.is_empty(), "base records should not be empty");
    }
}
