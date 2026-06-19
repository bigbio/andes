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
//!
//! ## Build-once cache
//!
//! [`index_cache_path`] derives a stable file path under the system temp dir from
//! a content-addressed key: the FASTA content (all target accessions + sequences),
//! the enzyme, missed cleavages, length bounds, and decoy prefix.  Callers may
//! pass this path to `prepare_mmap` (via [`MmapCandidateIndex::open_or_build`]);
//! when the file already exists at the key path the build step is skipped entirely
//! and the existing file is re-opened.  This makes repeated searches over the same
//! database pay the build cost only once.

use std::hash::{Hash, Hasher};
use std::collections::hash_map::DefaultHasher;
use std::io::{self, Write};
use std::path::{Path, PathBuf};

use crate::candidate_gen::enumerate_candidates;
use crate::search_index::SearchIndex;
use crate::search_params::SearchParams;

/// Serialised byte size of one `IndexRecord`.
pub const INDEX_RECORD_SIZE: usize = 20;

/// Magic bytes at the start of a candidate-index cache file. Lets the reader
/// reject foreign / corrupt files instead of trusting any 20-byte-aligned blob.
const INDEX_MAGIC: &[u8; 8] = b"ANDESCI\x00";
/// On-disk format version. BUMP whenever `IndexRecord`'s layout or the file
/// structure changes so a stale-format cache is rebuilt rather than misread.
const INDEX_FORMAT_VERSION: u32 = 1;
/// Header = 8 magic bytes + 4 version bytes (LE); records follow at this offset.
const INDEX_HEADER_SIZE: usize = INDEX_MAGIC.len() + 4;

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

/// Estimate whether the in-RAM candidate index (the `Vec<Candidate>` the default
/// `Ram` backing materializes) would fit within `budget_bytes`.
///
/// Used by `--candidate-index auto` to choose RAM vs out-of-core mmap WITHOUT
/// risking an OOM: the in-RAM path `.collect()`s every modified candidate (each
/// owning a `Peptide` residue `Vec`), which for large mod searches is the term
/// that blows up. We accumulate a per-candidate byte estimate over
/// [`enumerate_candidates`] and **bail early** the moment it crosses the budget,
/// so the cost is bounded by the budget (not by the — possibly enormous — total
/// candidate count). Returns `true` if the full set fits, `false` if it would
/// exceed `budget_bytes`.
///
/// The per-candidate cost is intentionally a slight OVER-estimate (so we err
/// toward mmap near the limit and never under-predict into an OOM): the
/// `Candidate` struct + its heap `Peptide` residues (the dominant term; each
/// `AminoAcid` carries an `Arc<Modification>`) + the mass-bucket index entry.
pub fn ram_candidate_index_fits(
    idx: &SearchIndex,
    params: &SearchParams,
    decoy_prefix: &str,
    budget_bytes: u64,
) -> bool {
    // Conservative per-residue / per-candidate constants (bytes).
    const PER_RESIDUE: u64 = 32; // AminoAcid + Arc refcount slot + slack
    const PER_CANDIDATE_OVERHEAD: u64 = 120; // struct + Vec header + bucket-idx entry + alloc overhead
    let mut acc: u64 = 0;
    for c in enumerate_candidates(idx, params, decoy_prefix) {
        acc = acc.saturating_add(
            PER_CANDIDATE_OVERHEAD + c.peptide.residues.len() as u64 * PER_RESIDUE,
        );
        if acc > budget_bytes {
            return false;
        }
    }
    true
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

    // Header: magic + format version, so a reader can reject a foreign/corrupt or
    // stale-format file (see `MmapCandidateIndex::open`).
    out.write_all(INDEX_MAGIC)?;
    out.write_all(&INDEX_FORMAT_VERSION.to_le_bytes())?;

    let n = records.len();
    for rec in &records {
        out.write_all(&rec.to_bytes())?;
    }
    Ok(n)
}

/// Derive a stable, content-addressed cache path for the base-peptide index.
///
/// The path is placed in the system temp directory and named
/// `andes-candidx-<hex16>.bin`, where `<hex16>` is a 16-character hex string
/// derived by hashing:
///
/// 1. Every **target** protein's accession and sequence (in iteration order from
///    [`SearchIndex::iter_target_proteins`]).  Decoy sequences are derived from
///    targets deterministically, so they need not be hashed separately.
/// 2. The search parameters that affect which base peptides are produced:
///    `enzyme`, `extra_enzymes`, `max_missed_cleavages`, `min_length`,
///    `max_length`, and `decoy_prefix`.
///
/// The hash uses Rust's stable [`DefaultHasher`].  It is **not** a
/// cryptographic hash — collisions are possible (though extremely unlikely in
/// practice).  The intent is cache correctness for a single machine / session,
/// not cross-machine portability.
///
/// # When to pass this path to `PreparedSearch::prepare_mmap`
///
/// Always.  `prepare_mmap` calls [`MmapCandidateIndex::open_or_build`], which
/// skips the (potentially slow) build step when the cache file already exists
/// and re-opens it directly.  The first call builds and caches; every
/// subsequent call with the same parameters opens the cached file.
pub fn index_cache_path(idx: &SearchIndex, params: &SearchParams) -> PathBuf {
    let mut h = DefaultHasher::new();

    // Hash ALL protein content (targets AND decoys): accession + sequence bytes.
    // Hashing the decoy sequences directly (not just `decoy_prefix`) captures the
    // decoy strategy and shuffle seed — a prefix-only key would let a
    // reverse-vs-shuffle (or reseeded shuffle) run silently reuse a stale decoy
    // index and corrupt FDR.
    for protein in &idx.db.proteins {
        protein.accession.hash(&mut h);
        protein.sequence.hash(&mut h);
    }
    idx.decoy_prefix.hash(&mut h);

    // Hash the FIXED modifications. Base-peptide masses fold fixed mods in (e.g.
    // Carbamidomethyl-C), so two runs with the same FASTA/enzyme/lengths but a
    // different fixed-mod set (e.g. label-free vs +TMT) MUST get different cache
    // keys — otherwise the second run reuses the first run's mass-shifted index and
    // silently loses IDs. Use the LOCATION-AWARE fingerprint, not `fixed_mod_deltas`
    // (residue+delta only): a fixed `*` +229 at N-term vs Anywhere fold different
    // base masses but share the same (residue, delta) set, so a delta-only key would
    // collide and reuse a stale index. (Variable mods are applied lazily at serve
    // time and do NOT affect the cached base index, so they are excluded.) The
    // fingerprint is already sorted + deduped → deterministic key.
    for (res, loc_rank, delta_bits) in params.aa_set.fixed_mod_fingerprint() {
        res.hash(&mut h);
        loc_rank.hash(&mut h);
        delta_bits.hash(&mut h);
    }

    // Hash all search params that affect base-peptide enumeration.
    // `enzyme` and `Enzyme` must impl `Hash` — confirmed by the model crate.
    params.enzyme.hash(&mut h);
    for extra in &params.extra_enzymes {
        extra.hash(&mut h);
    }
    params.max_missed_cleavages.hash(&mut h);
    params.min_length.hash(&mut h);
    params.max_length.hash(&mut h);
    params.num_tolerable_termini.hash(&mut h);

    let digest = h.finish();
    std::env::temp_dir().join(format!("andes-candidx-{digest:016x}.bin"))
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
    /// Validates the magic + format-version header and that the remaining bytes
    /// are a whole multiple of [`INDEX_RECORD_SIZE`]. Returns an error (treated
    /// upstream as a cache miss → rebuild) for a missing/foreign/truncated file
    /// or a format-version mismatch, rather than trusting any 20-byte-aligned
    /// blob as a valid index.
    pub fn open(path: &Path) -> io::Result<Self> {
        let file = std::fs::File::open(path)?;
        let meta = file.metadata()?;
        let byte_len = meta.len() as usize;
        if byte_len < INDEX_HEADER_SIZE {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("index file size {byte_len} is smaller than the {INDEX_HEADER_SIZE}-byte header"),
            ));
        }
        // SAFETY: The file is opened read-only and its contents are valid for
        // the lifetime of `mmap`.  No other process is expected to modify the
        // file while the search is running.
        let mmap = unsafe { memmap2::Mmap::map(&file)? };
        if &mmap[..INDEX_MAGIC.len()] != INDEX_MAGIC {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "index file is not an andes candidate index (bad magic)",
            ));
        }
        let version = u32::from_le_bytes(
            mmap[INDEX_MAGIC.len()..INDEX_HEADER_SIZE]
                .try_into()
                .expect("4 version bytes"),
        );
        if version != INDEX_FORMAT_VERSION {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("index file format version {version} != expected {INDEX_FORMAT_VERSION}"),
            ));
        }
        let data_len = byte_len - INDEX_HEADER_SIZE;
        if data_len % INDEX_RECORD_SIZE != 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("index payload {data_len} is not a multiple of {INDEX_RECORD_SIZE} bytes"),
            ));
        }
        let len = data_len / INDEX_RECORD_SIZE;
        Ok(Self { mmap, len })
    }

    /// Open the index from `path` if it already exists, or build it and then
    /// open it.
    ///
    /// This is the build-once-cache entry point.  Pass the path returned by
    /// [`index_cache_path`] to skip the build step on repeated invocations with
    /// the same parameters:
    ///
    /// - **Cache hit** (`path` exists and has a valid size): opens and mmaps
    ///   the existing file.  The build step is skipped entirely.
    /// - **Cache miss** (`path` does not exist): builds the index via
    ///   [`build_base_peptide_index`] and then opens the resulting file.
    ///
    /// Returns `(mmap_index, was_built)` where `was_built` is `true` when the
    /// index was freshly constructed and `false` when the existing cache was
    /// reused.  Callers may use `was_built` for diagnostic logging.
    pub fn open_or_build(
        path: &Path,
        idx: &SearchIndex,
        params: &SearchParams,
        decoy_prefix: &str,
    ) -> io::Result<(Self, bool)> {
        // Cache hit iff the file opens AND passes full validation (magic +
        // format version + payload alignment) — not just "exists and is
        // 20-byte-aligned". A foreign/corrupt/stale-format file fails `open` and
        // falls through to a rebuild rather than being trusted as an index.
        if let Ok(mmap_index) = Self::open(path) {
            return Ok((mmap_index, false));
        }

        // Build the index to a temp file alongside the target path, then rename
        // atomically so a partial write never leaves a corrupt cache entry.
        // The temp file lives in the same directory as `path` so the rename is
        // always within the same filesystem (avoiding cross-device link errors).
        let parent = path.parent().unwrap_or_else(|| Path::new("."));
        // Unique temp name per call: PID alone collides when several searches (or
        // tests) build distinct indexes concurrently in the SAME directory (e.g.
        // multiple `NamedTempFile`s under `/tmp`), racing on one tmp path and
        // corrupting each other's partial writes. Append a process-global atomic
        // counter so concurrent builds never share a tmp file.
        static BUILD_SEQ: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
        let seq = BUILD_SEQ.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let tmp_path = parent.join(format!(
            ".andes-candidx-build-{}-{}.tmp",
            std::process::id(),
            seq,
        ));
        {
            use std::io::BufWriter;
            let file = std::fs::File::create(&tmp_path)?;
            let mut bw = BufWriter::new(file);
            build_base_peptide_index(idx, params, decoy_prefix, &mut bw)?;
        }
        std::fs::rename(&tmp_path, path)?;

        let mmap_index = Self::open(path)?;
        Ok((mmap_index, true))
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
        let start = INDEX_HEADER_SIZE + i * INDEX_RECORD_SIZE;
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
        let start = INDEX_HEADER_SIZE + i * INDEX_RECORD_SIZE;
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
        // Skip + validate the magic/version header the builder now writes.
        assert!(
            data.len() >= INDEX_HEADER_SIZE && &data[..INDEX_MAGIC.len()] == INDEX_MAGIC,
            "binary blob missing the index header"
        );
        let payload = &data[INDEX_HEADER_SIZE..];
        assert_eq!(
            payload.len() % INDEX_RECORD_SIZE,
            0,
            "payload is not a multiple of {INDEX_RECORD_SIZE} bytes"
        );
        payload
            .chunks_exact(INDEX_RECORD_SIZE)
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
    fn ram_candidate_index_fits_respects_budget() {
        let idx = make_toy_index();
        let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
        let mut params = SearchParams::default_tryptic(aa_set);
        params.min_length = 3;
        // A generous budget fits the toy index; a zero budget never does.
        assert!(
            ram_candidate_index_fits(&idx, &params, "XXX", u64::MAX),
            "toy index must fit an unbounded budget"
        );
        assert!(
            !ram_candidate_index_fits(&idx, &params, "XXX", 0),
            "no candidate set fits a zero budget"
        );
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
        assert_eq!(buf.len(), INDEX_HEADER_SIZE + n * INDEX_RECORD_SIZE, "buffer size mismatch");

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

    // -----------------------------------------------------------------------
    // Build-once cache: `open_or_build` / `index_cache_path`
    // -----------------------------------------------------------------------

    /// `open_or_build` reuses the existing file on the second call (cache hit).
    ///
    /// Strategy: call `open_or_build` with an explicit temp path twice with
    /// the same inputs.  The first call must report `was_built = true`; the
    /// second must report `was_built = false` and must return a record set that
    /// is byte-identical to the first call's result.  The file's modification
    /// time must not change on the second call (proving the builder was not
    /// re-invoked).
    #[test]
    fn open_or_build_reuses_existing_cache_file() {
        let idx = make_toy_index();
        let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
        let mut params = SearchParams::default_tryptic(aa_set);
        params.min_length = 3;

        // Use a fresh named temp file as the cache path.  We keep the handle so
        // the file is not deleted between calls.
        let tmp = tempfile::NamedTempFile::new().expect("NamedTempFile");
        let cache_path = tmp.path();

        // First call: cache miss — must build.
        // We delete the content but keep the path alive by truncating, then
        // let open_or_build create a fresh file.  Actually `NamedTempFile::new`
        // creates a zero-byte file; `open_or_build` treats 0-byte files as a
        // cache miss (size % 20 == 0 but len == 0 fails the `len > 0` guard).
        let (mi1, was_built1) =
            MmapCandidateIndex::open_or_build(cache_path, &idx, &params, "XXX")
                .expect("first open_or_build failed");
        assert!(was_built1, "first call must build the index (cache miss)");
        assert!(mi1.len() > 0, "built index must be non-empty");

        // Record the modification time after the first build.
        let mtime_after_build = std::fs::metadata(cache_path)
            .expect("metadata after build")
            .modified()
            .expect("mtime after build");

        // Second call with identical inputs: cache hit — must NOT rebuild.
        let (mi2, was_built2) =
            MmapCandidateIndex::open_or_build(cache_path, &idx, &params, "XXX")
                .expect("second open_or_build failed");
        assert!(!was_built2, "second call must reuse the cache (cache hit)");

        // File must not have been rewritten (mtime unchanged).
        let mtime_after_reuse = std::fs::metadata(cache_path)
            .expect("metadata after reuse")
            .modified()
            .expect("mtime after reuse");
        assert_eq!(
            mtime_after_build, mtime_after_reuse,
            "cache file mtime changed on reuse — index was rebuilt unnecessarily"
        );

        // Record set must be byte-identical between build and reuse.
        let records1 = mi1.records();
        let records2 = mi2.records();
        assert_eq!(
            records1, records2,
            "cached records must be identical to freshly-built records"
        );
    }

    /// A foreign / header-less file (even one that is 20-byte-aligned) must be
    /// REJECTED by `open` and treated as a cache miss by `open_or_build`, not
    /// trusted as a valid index.
    #[test]
    fn open_rejects_header_less_file_and_rebuilds() {
        use std::io::Write;
        let idx = make_toy_index();
        let aa_set = AminoAcidSetBuilder::new_standard().build().unwrap();
        let mut params = SearchParams::default_tryptic(aa_set);
        params.min_length = 3;

        // A 40-byte blob (two 20-byte "records") with NO magic header.
        let tmp = tempfile::NamedTempFile::new().expect("NamedTempFile");
        tmp.as_file().write_all(&[0u8; 2 * INDEX_RECORD_SIZE]).unwrap();
        tmp.as_file().sync_all().unwrap();

        // `open` rejects it outright.
        assert!(
            MmapCandidateIndex::open(tmp.path()).is_err(),
            "open must reject a file with no magic header"
        );
        // `open_or_build` treats it as a miss and rebuilds a valid index.
        let (mi, was_built) =
            MmapCandidateIndex::open_or_build(tmp.path(), &idx, &params, "XXX")
                .expect("open_or_build");
        assert!(was_built, "header-less file must force a rebuild");
        assert!(mi.len() > 0, "rebuilt index must be non-empty");
    }

    /// `index_cache_path` returns distinct paths for different parameter sets.
    ///
    /// Two parameter sets that differ in a meaningful way (here: min_length)
    /// must produce different cache paths so a cached index built with the
    /// narrower set is not mistakenly reused for the wider set.
    #[test]
    fn index_cache_path_distinguishes_different_params() {
        let idx = make_toy_index();
        let aa_set_a = AminoAcidSetBuilder::new_standard().build().unwrap();
        let aa_set_b = AminoAcidSetBuilder::new_standard().build().unwrap();

        let mut params_a = SearchParams::default_tryptic(aa_set_a);
        params_a.min_length = 3;

        let mut params_b = SearchParams::default_tryptic(aa_set_b);
        params_b.min_length = 6; // different min_length

        let path_a = index_cache_path(&idx, &params_a);
        let path_b = index_cache_path(&idx, &params_b);

        assert_ne!(
            path_a, path_b,
            "cache paths must differ when min_length differs ({} vs {})",
            path_a.display(),
            path_b.display()
        );
    }

    /// `index_cache_path` distinguishes different FIXED-mod sets. Base masses fold
    /// fixed mods in, so a label-free vs +CAM-C run on the same FASTA/enzyme must
    /// NOT collide on one cache file (else the second run reuses a mass-shifted
    /// index and silently loses IDs).
    #[test]
    fn index_cache_path_distinguishes_fixed_mods() {
        let idx = make_toy_index();
        let mut p_nomod =
            SearchParams::default_tryptic(AminoAcidSetBuilder::new_standard().build().unwrap());
        p_nomod.min_length = 3;
        let mut p_cam = SearchParams::default_tryptic(
            AminoAcidSetBuilder::new_standard_with_carbamidomethyl_c().build().unwrap(),
        );
        p_cam.min_length = 3;
        assert_ne!(
            index_cache_path(&idx, &p_nomod),
            index_cache_path(&idx, &p_cam),
            "cache key must differ when the fixed-mod set differs"
        );
    }

    /// `index_cache_path` returns the same path when params are identical.
    #[test]
    fn index_cache_path_is_stable_for_identical_params() {
        let idx = make_toy_index();
        let aa_set_a = AminoAcidSetBuilder::new_standard().build().unwrap();
        let aa_set_b = AminoAcidSetBuilder::new_standard().build().unwrap();

        let mut params_a = SearchParams::default_tryptic(aa_set_a);
        params_a.min_length = 4;
        let mut params_b = SearchParams::default_tryptic(aa_set_b);
        params_b.min_length = 4; // same

        let path_a = index_cache_path(&idx, &params_a);
        let path_b = index_cache_path(&idx, &params_b);

        assert_eq!(
            path_a, path_b,
            "cache paths must be identical for identical params"
        );
    }
}
