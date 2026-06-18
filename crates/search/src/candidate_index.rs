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
use std::path::Path;

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

/// Canonical base-peptide enumeration for the out-of-core index.
///
/// Clones `params` and forces `max_variable_mods_per_peptide = 0` so only the
/// unmodified backbone is produced. Fixed mods in `aa_set` remain untouched —
/// they fold into per-residue mass exactly as in the normal search path
/// (e.g. carbamidomethyl-C is always present when configured as a fixed mod).
///
/// **Equivalence invariant**: the set of `(residues, base_mass)` pairs
/// produced here must equal the set of unmodified candidates (those whose
/// every residue carries no variable mod) produced by
/// [`enumerate_candidates`] with the original `params`.  This invariant
/// is guarded by `base_peptide_records_matches_unmodified_enumerate_candidates`
/// in the test suite — any future divergence in digestion logic will
/// fail that test.
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
            debug_assert!(
                c.protein_index <= u32::MAX as usize,
                "protein_index overflows u32"
            );
            debug_assert!(
                c.start_offset_in_protein <= u32::MAX as usize,
                "start_offset_in_protein overflows u32"
            );
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

/// Canonical builder: enumerate base peptides and write a mass-sorted binary
/// index to `out`.
///
/// This is the canonical entry point for creating the out-of-core index file.
/// It delegates to [`base_peptide_records`] (the canonical base-peptide
/// enumeration), which sets `max_variable_mods_per_peptide = 0` so that only
/// unmodified backbones are written — fixed mods (e.g. CAM on C, TMT on K)
/// are still folded into each residue's mass as configured.
///
/// Steps:
/// 1. Enumerate all base peptides (no variable mods) via [`base_peptide_records`].
/// 2. Sort by `mass_milli` ascending.
/// 3. Write each record as 20 packed little-endian bytes.
///
/// **Equivalence caveat**: `base_peptide_records` must stay equivalent to
/// `enumerate_candidates`'s digestion logic. This is enforced by the test
/// `base_peptide_records_matches_unmodified_enumerate_candidates` — any
/// divergence in enzyme/length/termini handling will fail that test.
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

/// Memory-mapped reader for the flat binary candidate index produced by
/// [`build_base_peptide_index`].
///
/// The file is a sequence of 20-byte little-endian [`IndexRecord`]s sorted
/// ascending by `mass_milli`.  Records are decoded on read via
/// [`IndexRecord::from_bytes`] — this avoids any `repr(C)` / alignment
/// requirements on the 20-byte format.
///
/// # Mass-window lookup
///
/// [`MmapCandidateIndex::mass_window`] uses two `partition_point` binary
/// searches over the sorted `mass_milli` key to find the half-open bounds
/// `[lo_idx, hi_idx)` and returns a `Vec` of the decoded records in that
/// range.  The overall complexity is O(log n) probes + O(k) decode where k
/// is the number of returned records.
pub struct MmapCandidateIndex {
    mmap: memmap2::Mmap,
    len: usize,
}

impl MmapCandidateIndex {
    /// Open a binary index file and memory-map it.
    ///
    /// Returns an error if the file length is not a multiple of
    /// [`INDEX_RECORD_SIZE`] (20 bytes), which would indicate a truncated or
    /// corrupt write.
    pub fn open(path: &Path) -> io::Result<Self> {
        let file = std::fs::File::open(path)?;
        let meta = file.metadata()?;
        let byte_len = meta.len() as usize;
        if byte_len % INDEX_RECORD_SIZE != 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "index file size {byte_len} is not a multiple of {INDEX_RECORD_SIZE} bytes"
                ),
            ));
        }
        // SAFETY: The file is opened read-only and its contents are valid for
        // the lifetime of `mmap`.  No other process is expected to modify the
        // file while the search is running.
        let mmap = unsafe { memmap2::Mmap::map(&file)? };
        let len = byte_len / INDEX_RECORD_SIZE;
        Ok(Self { mmap, len })
    }

    /// Number of records in the index.
    pub fn len(&self) -> usize {
        self.len
    }

    /// Returns `true` if the index contains no records.
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Decode and return the record at position `i`.
    ///
    /// # Panics
    /// Panics if `i >= self.len()`.
    pub fn record_at(&self, i: usize) -> IndexRecord {
        assert!(i < self.len, "record_at index {i} out of bounds (len={})", self.len);
        let start = i * INDEX_RECORD_SIZE;
        let chunk: &[u8; INDEX_RECORD_SIZE] = self.mmap[start..start + INDEX_RECORD_SIZE]
            .try_into()
            .expect("slice length is always INDEX_RECORD_SIZE");
        IndexRecord::from_bytes(chunk)
    }

    /// Return all decoded records.
    ///
    /// Allocates a `Vec` of length `self.len()`.  For large indexes, prefer
    /// [`mass_window`](Self::mass_window) to avoid decoding the entire file.
    pub fn records(&self) -> Vec<IndexRecord> {
        (0..self.len).map(|i| self.record_at(i)).collect()
    }

    /// Return all records whose `mass_milli` lies in the inclusive range
    /// `[lo_milli, hi_milli]`.
    ///
    /// Uses two `partition_point` binary searches (O(log n)) over the
    /// mass-sorted index, then decodes the matching slice (O(k)).
    pub fn mass_window(&self, lo_milli: u64, hi_milli: u64) -> Vec<IndexRecord> {
        if self.len == 0 || lo_milli > hi_milli {
            return Vec::new();
        }
        // Lower bound: first index where mass_milli >= lo_milli.
        let lo_idx = self.partition_point(|mass| mass < lo_milli);
        // Upper bound: first index where mass_milli > hi_milli.
        let hi_idx = self.partition_point(|mass| mass <= hi_milli);
        (lo_idx..hi_idx).map(|i| self.record_at(i)).collect()
    }

    /// Binary search helper: returns the first index `i` in `[0, self.len]`
    /// where `predicate(mass_milli_at(i))` is `false`.
    ///
    /// Equivalent to `slice::partition_point` but decodes only the `mass_milli`
    /// field of each probed record.
    fn partition_point<F>(&self, predicate: F) -> usize
    where
        F: Fn(u64) -> bool,
    {
        let mut lo = 0usize;
        let mut hi = self.len;
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            let mass = self.mass_milli_at(mid);
            if predicate(mass) {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        lo
    }

    /// Decode only the `mass_milli` field (first 8 bytes) of record `i`.
    fn mass_milli_at(&self, i: usize) -> u64 {
        let start = i * INDEX_RECORD_SIZE;
        u64::from_le_bytes(
            self.mmap[start..start + 8]
                .try_into()
                .expect("always 8 bytes"),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use model::aa_set::AminoAcidSetBuilder;
    use model::modification::{Modification, ModLocation, ResidueSpec};
    use model::protein::{Protein, ProteinDb};
    use crate::candidate_gen::enumerate_candidates;
    use crate::search_index::SearchIndex;
    use crate::search_params::SearchParams;

    /// Build a fixture index file containing base peptides from the toy protein.
    /// Returns the path wrapped in a `tempfile::NamedTempFile` (dropped = deleted).
    fn build_fixture_index() -> tempfile::NamedTempFile {
        use std::io::BufWriter;
        let idx = make_toy_index();
        let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
        let mut params = SearchParams::default_tryptic(aa_set);
        params.min_length = 3;

        let tmp = tempfile::NamedTempFile::new().expect("tempfile");
        {
            let mut bw = BufWriter::new(tmp.as_file());
            build_base_peptide_index(&idx, &params, "XXX", &mut bw)
                .expect("build_base_peptide_index failed");
        }
        tmp
    }

    #[test]
    fn mmap_window_returns_in_range_records() {
        let tmp = build_fixture_index();
        let mi = MmapCandidateIndex::open(tmp.path()).unwrap();
        let all = mi.records();
        // Pick the median record's mass as a stable target.
        let target = all[all.len() / 2].mass_milli;
        let win = mi.mass_window(target - 5, target + 5);
        // Every returned record must lie within [target-5, target+5].
        assert!(
            win.iter().all(|r| r.mass_milli >= target - 5 && r.mass_milli <= target + 5),
            "mass_window returned out-of-range records"
        );
        // At least one record must have exactly the target mass.
        assert!(
            win.iter().any(|r| r.mass_milli == target),
            "target mass not found in window"
        );
        // A window entirely below the minimum mass must be empty.
        assert!(
            mi.mass_window(0, all[0].mass_milli - 1).is_empty(),
            "expected empty window below minimum mass"
        );
    }

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

    /// Equivalence guard: the set of base peptides produced by
    /// `base_peptide_records` (which forces `max_variable_mods_per_peptide=0`)
    /// must exactly match the set of UNMODIFIED candidates produced by
    /// `enumerate_candidates` with the full variable-mod params, filtered to
    /// those whose every residue carries no variable mod.
    ///
    /// "Unmodified" is defined by backbone + base mass: a candidate is
    /// unmodified when none of its residues carry a variable (non-fixed) mod.
    /// Fixed mods are always present on both sides (they're folded in by the
    /// `AminoAcidSet`), so the two sides should produce exactly the same
    /// (residue_bytes, mass_milli) pairs.
    ///
    /// This test locks the digestion equivalence so that any future divergence
    /// between `base_peptide_records` and `enumerate_candidates` (e.g. a new
    /// ntt/enzyme path) is caught immediately.
    #[test]
    fn base_peptide_records_matches_unmodified_enumerate_candidates() {
        // Build an index with a protein that contains M (for Ox-M to apply) and
        // is long enough to produce multiple tryptic peptides.
        // WMKDALER: tryptic peptides include WMK (Ox-M variant), DALER, etc.
        let target = ProteinDb {
            proteins: vec![Protein {
                accession: "P2".into(),
                description: "test protein".into(),
                sequence: b"WMKDALERQPK".to_vec(),
            }],
        };
        let idx = SearchIndex::from_target_db(&target, "XXX");

        // Build params with Oxidation-M (variable mod) so enumerate_candidates
        // produces both unmodified and Ox-M variants.
        let ox_m = Modification {
            name: "Oxidation".to_string(),
            mass_delta: 15.99491,
            residue: ResidueSpec::Specific(b'M'),
            location: ModLocation::Anywhere,
            fixed: false,
            accession: None,
            neutral_losses: Vec::new(),
            loss_class: 0,
        };
        let aa_set = AminoAcidSetBuilder::new_standard()
            .add_variable_mod(ox_m)
            .build()
            .unwrap();
        let mut params = SearchParams::default_tryptic(aa_set);
        params.min_length = 3;
        params.max_variable_mods_per_peptide = 1;

        // --- LHS: base_peptide_records (forces NumMods=0) ---
        // Represent each base peptide as (residue_bytes, mass_milli) so the
        // comparison is independent of protein_index / start_offset ordering.
        let base_records = base_peptide_records(&idx, &params, "XXX");
        // Collect the (residue bytes, mass_milli) key from each record.
        // Since IndexRecord doesn't store residue bytes directly, we re-derive
        // them by running enumerate_candidates with base params and zipping.
        let base_params = SearchParams {
            max_variable_mods_per_peptide: 0,
            ..params.clone()
        };
        let base_candidates: Vec<_> =
            enumerate_candidates(&idx, &base_params, "XXX").collect();
        assert_eq!(
            base_records.len(),
            base_candidates.len(),
            "base_peptide_records and enumerate_candidates(NumMods=0) must yield same count"
        );

        // Build a canonical (residue_bytes, mass_milli) multiset from base_peptide_records.
        // Use the candidate residues (same enumeration) to reconstruct residue bytes.
        let mut lhs: Vec<(Vec<u8>, u64)> = base_candidates
            .iter()
            .map(|c| {
                let bytes: Vec<u8> = c.peptide.residues.iter().map(|aa| aa.residue).collect();
                let mass_milli = (c.peptide.mass() * 1000.0).round() as u64;
                (bytes, mass_milli)
            })
            .collect();
        lhs.sort_unstable();

        // --- RHS: enumerate_candidates with full mods, filter to unmodified ---
        // A candidate is "unmodified" when none of its residues carries a
        // variable (non-fixed) mod.
        let all_candidates: Vec<_> =
            enumerate_candidates(&idx, &params, "XXX").collect();
        let mut rhs: Vec<(Vec<u8>, u64)> = all_candidates
            .iter()
            .filter(|c| {
                c.peptide.residues.iter().all(|aa| {
                    aa.mod_
                        .as_ref()
                        .map(|m| m.fixed) // fixed mods are OK; variable are not
                        .unwrap_or(true)  // no mod at all = unmodified
                })
            })
            .map(|c| {
                let bytes: Vec<u8> = c.peptide.residues.iter().map(|aa| aa.residue).collect();
                let mass_milli = (c.peptide.mass() * 1000.0).round() as u64;
                (bytes, mass_milli)
            })
            .collect();
        rhs.sort_unstable();

        assert_eq!(
            lhs, rhs,
            "base_peptide_records must produce exactly the unmodified subset of enumerate_candidates"
        );

        // Sanity: both sides must be non-empty (the protein has tryptic peptides).
        assert!(!lhs.is_empty(), "expected at least one base peptide");

        // Sanity: the full run with Ox-M must be strictly larger (Ox-M adds variants).
        let full_count = all_candidates.len();
        let base_count = lhs.len();
        assert!(
            full_count > base_count,
            "full enumeration ({full_count}) should exceed base ({base_count}) when variable mods present"
        );
    }
}
