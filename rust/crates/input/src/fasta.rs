//! Streaming FASTA reader. Sync I/O — FASTA is line-oriented text, no
//! async benefit. Handcrafted parser (no regex) — FASTA is simple
//! enough that hand-rolling is clearer than pulling in a dep.

use std::io::BufRead;

use model::{Protein, ProteinDb};

pub struct FastaReader<R: BufRead> {
    reader: R,
    line_no: usize,
    buf: String,
    /// Lookahead — when we read a `>` line that starts the NEXT protein
    /// while finishing the current one, stash it here.
    pending_header: Option<String>,
}

impl<R: BufRead> FastaReader<R> {
    pub fn new(reader: R) -> Self {
        Self { reader, line_no: 0, buf: String::new(), pending_header: None }
    }

    /// Eager-load all proteins into a `ProteinDb`.
    pub fn load_all(reader: R) -> Result<ProteinDb, FastaParseError> {
        let mut proteins = Vec::new();
        for result in FastaReader::new(reader) {
            proteins.push(result?);
        }
        Ok(ProteinDb { proteins })
    }

    /// Read one line into `self.buf`. Returns `Ok(None)` at EOF.
    /// Advances `line_no`.
    fn read_one_line(&mut self) -> Result<Option<()>, FastaParseError> {
        self.buf.clear();
        let n = self.reader.read_line(&mut self.buf)
            .map_err(|source| FastaParseError::Io { line: self.line_no + 1, source })?;
        if n == 0 {
            Ok(None)
        } else {
            self.line_no += 1;
            Ok(Some(()))
        }
    }
}

impl<R: BufRead> Iterator for FastaReader<R> {
    type Item = Result<Protein, FastaParseError>;

    fn next(&mut self) -> Option<Self::Item> {
        let header_line = match self.pending_header.take() {
            Some(h) => h,
            None => loop {
                match self.read_one_line() {
                    Ok(None) => return None,
                    Ok(Some(())) => {}
                    Err(e) => return Some(Err(e)),
                }
                let trimmed = self.buf.trim();
                if trimmed.is_empty() || trimmed.starts_with(';') {
                    continue;
                }
                if !trimmed.starts_with('>') {
                    return Some(Err(FastaParseError::OrphanSequence {
                        line: self.line_no, got: trimmed.to_string(),
                    }));
                }
                break trimmed.to_string();
            },
        };

        let header_line_no = self.line_no;
        let body = &header_line[1..];
        let (accession, description) = split_header(body);
        if accession.is_empty() {
            return Some(Err(FastaParseError::EmptyAccession { line: header_line_no }));
        }

        let mut sequence = Vec::with_capacity(64);
        loop {
            match self.read_one_line() {
                Ok(None) => break,
                Ok(Some(())) => {}
                Err(e) => return Some(Err(e)),
            }
            let trimmed = self.buf.trim();
            if trimmed.is_empty() || trimmed.starts_with(';') {
                continue;
            }
            if trimmed.starts_with('>') {
                self.pending_header = Some(trimmed.to_string());
                break;
            }
            for ch in trimmed.bytes() {
                if !ch.is_ascii_whitespace() {
                    sequence.push(ch.to_ascii_uppercase());
                }
            }
        }

        Some(Ok(Protein { accession, description, sequence }))
    }
}

fn split_header(s: &str) -> (String, String) {
    let s = s.trim_start();
    if let Some(idx) = s.find(char::is_whitespace) {
        let acc = s[..idx].to_string();
        let desc = s[idx..].trim().to_string();
        (acc, desc)
    } else {
        (s.to_string(), String::new())
    }
}

#[derive(thiserror::Error, Debug)]
pub enum FastaParseError {
    #[error("I/O error at line {line}: {source}")]
    Io { line: usize, #[source] source: std::io::Error },
    #[error("malformed FASTA at line {line}: expected `>` at start of header, got {got:?}")]
    NotAHeader { line: usize, got: String },
    #[error("FASTA header at line {line} has empty accession")]
    EmptyAccession { line: usize },
    #[error("sequence data at line {line} appears before any `>` header: {got:?}")]
    OrphanSequence { line: usize, got: String },
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn split_header_with_description() {
        let (a, d) = split_header("P1 some description here");
        assert_eq!(a, "P1");
        assert_eq!(d, "some description here");
    }

    #[test]
    fn split_header_no_description() {
        let (a, d) = split_header("P1");
        assert_eq!(a, "P1");
        assert_eq!(d, "");
    }

    #[test]
    fn split_header_empty() {
        let (a, d) = split_header("");
        assert_eq!(a, "");
        assert_eq!(d, "");
    }

    #[test]
    fn split_header_leading_whitespace_trimmed() {
        let (a, d) = split_header("  P1 desc");
        assert_eq!(a, "P1");
        assert_eq!(d, "desc");
    }
}
