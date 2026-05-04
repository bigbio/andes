//! Concatenated-byte representation of a ProteinDb.
//! Mirrors Java `edu.ucsd.msjava.msdbsearch.CompactFastaSequence`.
//! Used as input to suffix-array construction.
//!
//! # Wire format details (verified against Java CompactFastaSequence + BSA fixtures)
//!
//! ## `.cseq` binary layout (big-endian)
//! ```text
//! i32   size          — number of body bytes (= total sequence length)
//! i32   formatId      — always 9873 (COMPACT_FASTA_SEQUENCE_FILE_FORMAT_ID)
//! i32   id            — UUID.randomUUID().hashCode() written at creation time
//! i64   lastModified  — milliseconds since epoch of source FASTA
//! u8[size]            — encoded residue body
//! ```
//! Total file size = 20 + size bytes. Verified: BSA.cseq is 629 bytes = 20 + 609.
//!
//! ## `.canno` text layout (line-based)
//! ```text
//! Line 1: formatId           e.g. "9873"
//! Line 2: id                 e.g. "816949726"
//! Line 3: lastModified ms    e.g. "1777316603419"
//! Line 4: alphabet           e.g. "A:B:C:D:E:F:G:H:I:J:K:L:M:N:O:P:Q:R:S:T:U:V:W:X:Y:Z"
//! Line 5+: <endOffset>:<annotation>   one per protein
//! ```
//! `endOffset` is the position of the TERMINATOR byte that follows the protein
//! (i.e., one past the last residue byte, same as Java's `offset` var after the terminator is written).
//! Verified: BSA.canno has "609:sp|P02769|ALBU_BOVIN ..." and BSA.cseq body[609-1] == TERMINATOR.
//!
//! ## Residue encoding (alphabet-indexed)
//! Java's `initializeAlphabet` assigns:
//! - byte 0  → TERMINATOR ('_')
//! - byte 1  → INVALID_CHAR_CODE ('?')
//! - byte 2  → first group in alphabet (= 'A' for CAPITAL_LETTERS_26)
//! - byte 3  → 'B', ..., byte 27 → 'Z'
//!
//! So `residue_to_byte('M') = ord('M') - ord('A') + 2 = 14`. Verified: BSA.cseq body[1] = 0x0e = 14,
//! and BSA starts with 'M' (Methionine).
//!
//! ## Sequence layout
//! `[TERM] <protein0 residues> [TERM] <protein1 residues> [TERM]`
//! The leading TERMINATOR is written before the first protein (Java writes a TERMINATOR each time
//! it sees a `>` header line, including the first one). The trailing TERMINATOR closes the last protein.
//! Each annotation's `endOffset` points to the TERMINATOR at the end of that protein (exclusive of residues).
//!
//! ## Rust representation difference from Java
//! Java's annotation TreeMap is keyed by `endOffset` (terminator position).
//! In this Rust struct, `ProteinAnnotation.start` stores the offset of the FIRST residue
//! (= Java endOffset of previous protein, which is the position after the leading terminator of this protein).
//! On write, we compute `end_offset = start + sequence_len + 1` (+ 1 for the trailing terminator).

use std::io::{Read, Write};

use crate::protein::ProteinDb;

/// Java's `COMPACT_FASTA_SEQUENCE_FILE_FORMAT_ID`.
pub const FORMAT_ID: i32 = 9873;

/// End-of-sequence / protein-delimiter terminator byte. Java Constants.TERMINATOR = 0.
pub const TERMINATOR: u8 = 0;

/// Invalid character code (byte 1). Java Constants.INVALID_CHAR_CODE = 1.
pub const INVALID_CHAR_CODE: u8 = 1;

/// Fixed alphabet matching Java's default `Constants.CAPITAL_LETTERS_26`.
/// Index 0 = TERMINATOR placeholder ('_'); indices 1+ are unused in this table.
/// Encoding: byte 0 = TERMINATOR, byte 1 = INVALID, byte 2 = 'A', ..., byte 27 = 'Z'.
/// Verified against BSA.cseq + BSA.canno fixtures.
pub const ALPHABET: &[u8] = b"_ABCDEFGHIJKLMNOPQRSTUVWXYZ";

/// Encode an ASCII uppercase residue to its storage byte.
/// Non-uppercase or unknown residues encode to INVALID_CHAR_CODE (1), matching Java's fallback
/// (Java uses `INVALID_CHAR_CODE` when `alpha2byte.get(c)` returns null).
#[inline]
pub fn residue_to_byte(residue: u8) -> u8 {
    if residue.is_ascii_uppercase() {
        residue - b'A' + 2
    } else {
        INVALID_CHAR_CODE
    }
}

/// Decode a storage byte back to its ASCII residue character.
/// Byte 0 → '_' (TERMINATOR), byte 1 → '?' (INVALID), bytes 2-27 → 'A'-'Z'.
#[inline]
pub fn byte_to_residue(b: u8) -> u8 {
    match b {
        0 => b'_',
        1 => b'?',
        2..=27 => b'A' + b - 2,
        _ => b'?',
    }
}

#[derive(Debug, Clone)]
pub struct CompactFastaSequence {
    /// Encoded sequence body: `[TERM] <protein0> [TERM] <protein1> [TERM] ...`
    /// Body bytes are alphabet indices, not raw ASCII.
    pub sequence: Vec<u8>,
    pub annotations: Vec<ProteinAnnotation>,
    /// Number of body bytes (= sequence.len()).
    pub size: u64,
}

#[derive(Debug, Clone)]
pub struct ProteinAnnotation {
    /// Offset into `sequence` of this protein's FIRST residue byte.
    /// (One past the leading TERMINATOR for this protein.)
    pub start: u64,
    pub accession: String,
    pub description: String,
}

impl CompactFastaSequence {
    /// Build an in-memory `CompactFastaSequence` from a `ProteinDb`.
    ///
    /// Layout: `[TERM] <encoded protein0> [TERM] <encoded protein1> [TERM]`
    pub fn from_protein_db(db: &ProteinDb) -> Self {
        if db.proteins.is_empty() {
            return Self {
                sequence: Vec::new(),
                annotations: Vec::new(),
                size: 0,
            };
        }

        let mut sequence = Vec::with_capacity(
            db.proteins.iter().map(|p| p.sequence.len() + 1).sum::<usize>() + 1,
        );
        let mut annotations = Vec::with_capacity(db.proteins.len());

        // Lead with TERMINATOR (matches Java: writeByte(TERMINATOR) on each '>' line).
        sequence.push(TERMINATOR);
        for p in &db.proteins {
            let start = sequence.len() as u64;
            for &residue in &p.sequence {
                sequence.push(residue_to_byte(residue));
            }
            sequence.push(TERMINATOR);
            annotations.push(ProteinAnnotation {
                start,
                accession: p.accession.clone(),
                description: p.description.clone(),
            });
        }

        let size = sequence.len() as u64;
        Self {
            sequence,
            annotations,
            size,
        }
    }

    pub fn protein_count(&self) -> usize {
        self.annotations.len()
    }

    /// Binary-search the annotation array for the protein containing
    /// position `pos`. Returns `None` for positions before the first protein.
    pub fn protein_index_at(&self, pos: u64) -> Option<usize> {
        if self.annotations.is_empty() {
            return None;
        }
        match self.annotations.binary_search_by(|a| a.start.cmp(&pos)) {
            Ok(idx) => Some(idx),
            Err(0) => None,
            Err(idx) => Some(idx - 1),
        }
    }

    /// Write `(.cseq, .canno)` byte streams in Java-compatible format.
    ///
    /// The `formatId` is written as 9873. `id` and `lastModified` are written as 0
    /// (placeholder values; Java regenerates the index on mismatch anyway).
    pub fn write_to<W1: Write, W2: Write>(
        &self,
        cseq: &mut W1,
        canno: &mut W2,
    ) -> Result<(), CompactFastaError> {
        // .cseq header: i32 size | i32 formatId | i32 id | i64 lastModified
        cseq.write_all(&(self.size as i32).to_be_bytes())?;
        cseq.write_all(&FORMAT_ID.to_be_bytes())?;
        cseq.write_all(&0_i32.to_be_bytes())?; // id placeholder
        cseq.write_all(&0_i64.to_be_bytes())?; // lastModified placeholder
        cseq.write_all(&self.sequence)?;

        // .canno: text format
        writeln!(canno, "{FORMAT_ID}")?; // formatId
        writeln!(canno, "0")?; // id placeholder
        writeln!(canno, "0")?; // lastModified placeholder
        // Alphabet: "A:B:C:...:Z"  (ALPHABET[1..] strips the leading '_' placeholder)
        let alpha_str: String = ALPHABET[1..]
            .iter()
            .map(|&c| (c as char).to_string())
            .collect::<Vec<_>>()
            .join(":");
        writeln!(canno, "{alpha_str}")?;

        // Annotation lines: <endOffset>:<accession> <description>
        //
        // Java emits the offset with an inconsistency between non-last and last proteins:
        // - Non-last protein: endOffset = position of the inter-protein TERMINATOR byte (0-indexed).
        //   Emitted BEFORE offset++ in Java's loop, so offset = TERM position.
        // - Last protein: endOffset = size (= TERM position + 1), emitted AFTER offset++.
        //
        // This means: on read, start_of_protein_N = canno_offset_of_(N-1) + 1.
        // We replicate this exactly so files are Java-compatible.
        let n = self.annotations.len();
        for (i, ann) in self.annotations.iter().enumerate() {
            let protein_len = self
                .sequence
                .get(ann.start as usize..)
                .map(|s| s.iter().position(|&b| b == TERMINATOR).unwrap_or(s.len()))
                .unwrap_or(0);
            // Non-last: TERM position = start + protein_len.
            // Last: size = start + protein_len + 1.
            let end_offset = if i + 1 < n {
                ann.start + protein_len as u64 // TERM position (0-indexed)
            } else {
                self.size // = start + protein_len + 1
            };
            if ann.description.is_empty() {
                writeln!(canno, "{}:{}", end_offset, ann.accession)?;
            } else {
                writeln!(
                    canno,
                    "{}:{} {}",
                    end_offset, ann.accession, ann.description
                )?;
            }
        }
        Ok(())
    }

    /// Read `(.cseq, .canno)` byte streams written in Java format.
    pub fn read_from<R1: Read, R2: Read>(
        cseq: &mut R1,
        canno: &mut R2,
    ) -> Result<Self, CompactFastaError> {
        // Parse .cseq header: i32 size | i32 formatId | i32 id | i64 lastModified
        let mut size_buf = [0u8; 4];
        cseq.read_exact(&mut size_buf)?;
        let size = i32::from_be_bytes(size_buf) as u64;

        // Skip formatId (i32), id (i32), lastModified (i64) = 16 bytes
        let mut skip_buf = [0u8; 16];
        cseq.read_exact(&mut skip_buf)?;

        // Read body
        let mut sequence = vec![0u8; size as usize];
        cseq.read_exact(&mut sequence)?;

        // Parse .canno text
        let mut canno_text = String::new();
        canno.read_to_string(&mut canno_text)?;
        let mut lines = canno_text.lines();

        let _format_id = lines.next().ok_or_else(|| CompactFastaError::MalformedCanno {
            line: 1,
            message: "missing line 1 (formatId)".to_string(),
        })?;
        let _id = lines.next().ok_or_else(|| CompactFastaError::MalformedCanno {
            line: 2,
            message: "missing line 2 (id)".to_string(),
        })?;
        let _last_modified = lines.next().ok_or_else(|| CompactFastaError::MalformedCanno {
            line: 3,
            message: "missing line 3 (lastModified)".to_string(),
        })?;
        let _alphabet = lines.next().ok_or_else(|| CompactFastaError::MalformedCanno {
            line: 4,
            message: "missing line 4 (alphabet)".to_string(),
        })?;

        // Parse annotation lines: <endOffset>:<annotation>
        // endOffset is the position of the trailing TERMINATOR (one past last residue).
        // We derive start = endOffset of previous protein (or 1 for the first protein,
        // because layout is [TERM=0] <protein0> [TERM=end0] <protein1> [TERM=end1] ...)
        let mut annotations = Vec::new();
        let mut prev_end: u64 = 1; // first protein starts at offset 1 (after leading TERM)

        for (i, line) in lines.enumerate() {
            let line_no = 5 + i;
            let (offset_str, ann_str) =
                line.split_once(':').ok_or_else(|| CompactFastaError::MalformedCanno {
                    line: line_no,
                    message: format!("expected `offset:annotation`, got {line:?}"),
                })?;
            let end_offset: u64 =
                offset_str
                    .parse()
                    .map_err(|e: std::num::ParseIntError| CompactFastaError::MalformedCanno {
                        line: line_no,
                        message: format!("bad offset {offset_str:?}: {e}"),
                    })?;
            let (accession, description) = match ann_str.split_once(' ') {
                Some((a, d)) => (a.to_string(), d.to_string()),
                None => (ann_str.to_string(), String::new()),
            };
            annotations.push(ProteinAnnotation {
                start: prev_end,
                accession,
                description,
            });
            // Next protein starts one byte after this protein's TERMINATOR.
            prev_end = end_offset + 1;
        }

        Ok(Self {
            sequence,
            annotations,
            size,
        })
    }
}

#[derive(thiserror::Error, Debug)]
pub enum CompactFastaError {
    #[error("I/O error: {source}")]
    Io {
        #[from]
        source: std::io::Error,
    },
    #[error("malformed .canno line {line}: {message}")]
    MalformedCanno { line: usize, message: String },
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
    fn empty_db_produces_zero_proteins() {
        let db = ProteinDb::new();
        let cf = CompactFastaSequence::from_protein_db(&db);
        assert_eq!(cf.protein_count(), 0);
        assert_eq!(cf.annotations.len(), 0);
    }

    #[test]
    fn single_protein_sequence_is_preserved() {
        let db = make_db(&[("P1", b"MKWV")]);
        let cf = CompactFastaSequence::from_protein_db(&db);
        assert_eq!(cf.protein_count(), 1);
        assert_eq!(cf.annotations[0].accession, "P1");
        let start = cf.annotations[0].start as usize;
        let expected_bytes: Vec<u8> = b"MKWV".iter().map(|&r| residue_to_byte(r)).collect();
        assert_eq!(&cf.sequence[start..start + 4], &expected_bytes[..]);
    }

    #[test]
    fn two_proteins_have_separator_between() {
        let db = make_db(&[("P1", b"AB"), ("P2", b"CD")]);
        let cf = CompactFastaSequence::from_protein_db(&db);
        assert_eq!(cf.protein_count(), 2);
        let start1 = cf.annotations[0].start as usize;
        let start2 = cf.annotations[1].start as usize;
        // Each protein 2 bytes; at least one separator byte between them.
        assert!(
            start2 > start1 + 2,
            "expected separator between proteins; start1={start1}, start2={start2}"
        );
        // The byte between protein 1's end and protein 2's start should be TERMINATOR.
        assert_eq!(cf.sequence[start1 + 2], TERMINATOR);
    }

    #[test]
    fn protein_index_at_returns_correct_index() {
        let db = make_db(&[("P1", b"ABC"), ("P2", b"DEF"), ("P3", b"GHI")]);
        let cf = CompactFastaSequence::from_protein_db(&db);
        let p1_start = cf.annotations[0].start;
        assert_eq!(cf.protein_index_at(p1_start), Some(0));
        let p2_start = cf.annotations[1].start;
        assert_eq!(cf.protein_index_at(p2_start), Some(1));
        let p3_start = cf.annotations[2].start;
        assert_eq!(cf.protein_index_at(p3_start), Some(2));
    }

    #[test]
    fn description_preserved() {
        let mut db = make_db(&[("P1", b"AB")]);
        db.proteins[0].description = "test description".into();
        let cf = CompactFastaSequence::from_protein_db(&db);
        assert_eq!(cf.annotations[0].description, "test description");
    }

    #[test]
    fn size_matches_sequence_length() {
        let db = make_db(&[("P1", b"AB"), ("P2", b"CD")]);
        let cf = CompactFastaSequence::from_protein_db(&db);
        assert_eq!(cf.size, cf.sequence.len() as u64);
    }
}
