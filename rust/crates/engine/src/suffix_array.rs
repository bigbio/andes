//! Suffix array + LCP over a CompactFastaSequence. Built via the
//! `suffix` crate (SA-IS algorithm). LCP via Kasai's algorithm.
//! Mirrors Java `edu.ucsd.msjava.msdbsearch.CompactSuffixArray` for
//! shape (`indices`, `nlcps`); algorithm differs (Java has its own
//! implementation), so byte-bit parity with Java's `.csarr` is NOT
//! required by Phase 4c (only candidate-set parity downstream in 4d).
//!
//! ## Wire format (`.csarr` / `.cnlcp`)
//!
//! Determined by reading `CompactSuffixArray.createSuffixArrayFiles()` in Java:
//!
//! ```text
//! .csarr:  i32 size  |  i32 id  |  i32[size] indices  |  i64 lastModified  |  i32 formatId
//! .cnlcp:  i32 size  |  i32 id  |  byte[size] nlcps   |  i64 lastModified  |  i32 formatId
//! ```
//!
//! `formatId` = 8294 (`COMPACT_SUFFIX_ARRAY_FILE_FORMAT_ID` in Java).
//! All multibyte integers are big-endian (Java `DataOutputStream`).
//! The `id` and `lastModified` fields are used by Java for cache validation;
//! Rust writes zeros for these fields (round-trip fidelity, not cache linking).

use std::io::{Read, Write};

use crate::compact_fasta::CompactFastaSequence;

#[derive(Debug, Clone)]
pub struct SuffixArray {
    /// Sorted suffix start positions over `compact.sequence`.
    pub indices: Vec<i32>,
    /// Nearest-LCP array. `nlcps[i]` = LCP between suffixes at
    /// `indices[i-1]` and `indices[i]`. `nlcps[0]` is conventionally 0.
    pub nlcps: Vec<i32>,
}

impl SuffixArray {
    /// Build a SA + LCP from a CompactFastaSequence.
    ///
    /// The `suffix` crate works on UTF-8 strings; the CompactFastaSequence
    /// guarantees ASCII content (residues + SEPARATOR/TERMINATOR) so the
    /// transmute through `from_utf8_unchecked` is safe.
    pub fn build(compact: &CompactFastaSequence) -> Self {
        if compact.sequence.is_empty() {
            return Self { indices: Vec::new(), nlcps: Vec::new() };
        }

        // SAFETY: CompactFastaSequence::from_protein_db only emits
        // ASCII bytes (uppercase residues + SEPARATOR=b'_' + TERMINATOR=0).
        // All ASCII bytes are valid single-byte UTF-8 codepoints.
        let s: &str = unsafe { std::str::from_utf8_unchecked(&compact.sequence) };

        let suffix_table = suffix::SuffixTable::new(s);
        // SuffixTable.table() returns &[u32] of byte positions into the
        // original string. We assume length fits in i32 for protein
        // databases (~600M bytes max). This is consistent with Java's i32 indices.
        let raw_indices = suffix_table.table();
        let indices: Vec<i32> = raw_indices.iter().map(|&i| i as i32).collect();

        let nlcps = compute_lcp(&compact.sequence, &indices);

        Self { indices, nlcps }
    }
}

/// Kasai's algorithm. Returns nearest-LCP array aligned with `indices`.
fn compute_lcp(text: &[u8], indices: &[i32]) -> Vec<i32> {
    let n = text.len();
    if n == 0 {
        return Vec::new();
    }
    // rank[i] = position of suffix starting at text[i..] in the sorted SA.
    let mut rank = vec![0i32; n];
    for (i, &sa_i) in indices.iter().enumerate() {
        rank[sa_i as usize] = i as i32;
    }
    let mut lcp = vec![0i32; n];
    let mut h: i32 = 0;
    for i in 0..n {
        if rank[i] > 0 {
            let j = indices[(rank[i] - 1) as usize] as usize;
            while i + (h as usize) < n
                && j + (h as usize) < n
                && text[i + h as usize] == text[j + h as usize]
            {
                h += 1;
            }
            lcp[rank[i] as usize] = h;
            if h > 0 {
                h -= 1;
            }
        } else {
            h = 0;
        }
    }
    lcp
}

/// Java `COMPACT_SUFFIX_ARRAY_FILE_FORMAT_ID`.
const FORMAT_ID: i32 = 8294;

impl SuffixArray {
    /// Serialize to `.csarr` and `.cnlcp` streams in Java-compatible wire format.
    ///
    /// Writes placeholder zeros for the `id` and `lastModified` header/footer
    /// fields (Java uses these for cache linking; Rust does not need them for
    /// round-trip or search purposes).
    pub fn write_to<W1: Write, W2: Write>(
        &self,
        csarr: &mut W1,
        cnlcp: &mut W2,
    ) -> Result<(), SuffixArrayError> {
        write_csarr(csarr, &self.indices)?;
        write_cnlcp(cnlcp, &self.nlcps)?;
        Ok(())
    }

    /// Deserialize from `.csarr` and `.cnlcp` streams in Java-compatible wire format.
    pub fn read_from<R1: Read, R2: Read>(
        csarr: &mut R1,
        cnlcp: &mut R2,
    ) -> Result<Self, SuffixArrayError> {
        let indices = read_csarr(csarr)?;
        let nlcps = read_cnlcp(cnlcp)?;
        if indices.len() != nlcps.len() {
            return Err(SuffixArrayError::LengthMismatch {
                indices: indices.len(),
                nlcps: nlcps.len(),
            });
        }
        Ok(Self { indices, nlcps })
    }
}

/// Write `.csarr`: `i32 size | i32 id=0 | i32[size] indices | i64 lastModified=0 | i32 formatId`.
fn write_csarr<W: Write>(w: &mut W, indices: &[i32]) -> Result<(), SuffixArrayError> {
    let size = indices.len() as i32;
    w.write_all(&size.to_be_bytes())?;
    w.write_all(&0_i32.to_be_bytes())?; // id placeholder
    for &v in indices {
        w.write_all(&v.to_be_bytes())?;
    }
    w.write_all(&0_i64.to_be_bytes())?; // lastModified placeholder
    w.write_all(&FORMAT_ID.to_be_bytes())?;
    Ok(())
}

/// Write `.cnlcp`: `i32 size | i32 id=0 | byte[size] nlcps | i64 lastModified=0 | i32 formatId`.
///
/// Java's LCP values are written as single bytes (via `writeByte`), capped at
/// [`i8::MAX`] (127). Values that exceed 127 are clamped before writing.
fn write_cnlcp<W: Write>(w: &mut W, nlcps: &[i32]) -> Result<(), SuffixArrayError> {
    let size = nlcps.len() as i32;
    w.write_all(&size.to_be_bytes())?;
    w.write_all(&0_i32.to_be_bytes())?; // id placeholder
    for &v in nlcps {
        let b = v.clamp(0, i8::MAX as i32) as u8;
        w.write_all(&[b])?;
    }
    w.write_all(&0_i64.to_be_bytes())?; // lastModified placeholder
    w.write_all(&FORMAT_ID.to_be_bytes())?;
    Ok(())
}

/// Read `.csarr`: parse size, skip id, read `size` i32 values, skip footer.
fn read_csarr<R: Read>(r: &mut R) -> Result<Vec<i32>, SuffixArrayError> {
    let mut buf4 = [0u8; 4];

    r.read_exact(&mut buf4)?;
    let size = i32::from_be_bytes(buf4) as usize;

    // skip id (4 bytes)
    r.read_exact(&mut buf4)?;

    let mut out = Vec::with_capacity(size);
    for _ in 0..size {
        r.read_exact(&mut buf4)?;
        out.push(i32::from_be_bytes(buf4));
    }

    // skip footer: i64 lastModified (8 bytes) + i32 formatId (4 bytes) = 12 bytes
    let mut footer = [0u8; 12];
    r.read_exact(&mut footer)?;

    Ok(out)
}

/// Read `.cnlcp`: parse size, skip id, read `size` bytes as i32 (sign-extended), skip footer.
fn read_cnlcp<R: Read>(r: &mut R) -> Result<Vec<i32>, SuffixArrayError> {
    let mut buf4 = [0u8; 4];

    r.read_exact(&mut buf4)?;
    let size = i32::from_be_bytes(buf4) as usize;

    // skip id (4 bytes)
    r.read_exact(&mut buf4)?;

    let mut out = Vec::with_capacity(size);
    let mut byte_buf = [0u8; 1];
    for _ in 0..size {
        r.read_exact(&mut byte_buf)?;
        // Java writeByte / readByte is signed; extend to i32 matching Java semantics.
        out.push(byte_buf[0] as i8 as i32);
    }

    // skip footer: i64 lastModified (8 bytes) + i32 formatId (4 bytes) = 12 bytes
    let mut footer = [0u8; 12];
    r.read_exact(&mut footer)?;

    Ok(out)
}

#[derive(thiserror::Error, Debug)]
pub enum SuffixArrayError {
    #[error("I/O error: {source}")]
    Io {
        #[from]
        source: std::io::Error,
    },
    #[error(".csarr length {indices} != .cnlcp length {nlcps}")]
    LengthMismatch { indices: usize, nlcps: usize },
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::protein::{Protein, ProteinDb};

    fn make_db(proteins: &[(&str, &[u8])]) -> ProteinDb {
        ProteinDb {
            proteins: proteins
                .iter()
                .map(|(acc, seq)| Protein {
                    accession: acc.to_string(),
                    description: String::new(),
                    sequence: seq.to_vec(),
                })
                .collect(),
        }
    }

    #[test]
    fn small_sa_has_expected_length() {
        let db = make_db(&[("P1", b"AB")]);
        let cf = CompactFastaSequence::from_protein_db(&db);
        let sa = SuffixArray::build(&cf);
        assert_eq!(sa.indices.len(), cf.sequence.len());
        assert_eq!(sa.nlcps.len(), cf.sequence.len());
    }

    #[test]
    fn sa_indices_are_a_permutation_of_positions() {
        let db = make_db(&[("P1", b"BANANA")]);
        let cf = CompactFastaSequence::from_protein_db(&db);
        let sa = SuffixArray::build(&cf);
        let n = cf.sequence.len();
        let mut seen = vec![false; n];
        for &i in &sa.indices {
            assert!((i as usize) < n, "index {i} out of bounds for len {n}");
            assert!(!seen[i as usize], "index {i} repeated");
            seen[i as usize] = true;
        }
        assert!(seen.iter().all(|&x| x), "not all positions covered");
    }

    #[test]
    fn sa_orders_suffixes_lexicographically() {
        let db = make_db(&[("P1", b"BANANA")]);
        let cf = CompactFastaSequence::from_protein_db(&db);
        let sa = SuffixArray::build(&cf);
        for i in 0..sa.indices.len() - 1 {
            let a = &cf.sequence[sa.indices[i] as usize..];
            let b = &cf.sequence[sa.indices[i + 1] as usize..];
            assert!(
                a <= b,
                "suffix order broken at i={}: {:?} vs {:?}",
                i,
                a,
                b
            );
        }
    }

    #[test]
    fn lcp_values_are_correct() {
        let db = make_db(&[("P1", b"ABAB")]);
        let cf = CompactFastaSequence::from_protein_db(&db);
        let sa = SuffixArray::build(&cf);
        for i in 1..sa.indices.len() {
            let a = &cf.sequence[sa.indices[i - 1] as usize..];
            let b = &cf.sequence[sa.indices[i] as usize..];
            let actual_lcp = a
                .iter()
                .zip(b.iter())
                .take_while(|(x, y)| x == y)
                .count();
            assert_eq!(
                sa.nlcps[i] as usize,
                actual_lcp,
                "LCP mismatch at i={}: indices[{}]={}, indices[{}]={}, suffixes={:?} vs {:?}",
                i,
                i - 1,
                sa.indices[i - 1],
                i,
                sa.indices[i],
                a,
                b
            );
        }
    }
}
