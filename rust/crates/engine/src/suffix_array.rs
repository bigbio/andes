//! Suffix array + LCP over a CompactFastaSequence. Built via the
//! `suffix` crate (SA-IS algorithm). LCP via Kasai's algorithm.
//! Mirrors Java `edu.ucsd.msjava.msdbsearch.CompactSuffixArray` for
//! shape (`indices`, `nlcps`); algorithm differs (Java has its own
//! implementation), so byte-bit parity with Java's `.csarr` is NOT
//! required by Phase 4c (only candidate-set parity downstream in 4d).

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

#[derive(thiserror::Error, Debug)]
pub enum SuffixArrayError {
    #[error("I/O error: {source}")]
    Io {
        #[from]
        source: std::io::Error,
    },
    #[error("file I/O not yet implemented (Phase 4c/Task 5 stub)")]
    NotYetImplemented,
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
