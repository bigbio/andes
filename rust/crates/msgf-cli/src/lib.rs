//! msgf-diff library — schema parsing + comparison for .pin files.

use std::fs;
use std::path::Path;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum DiffError {
    #[error("cannot read {path}: {source}")]
    Read {
        path: String,
        source: std::io::Error,
    },
    #[error("file is empty or missing header: {path}")]
    EmptyHeader { path: String },
}

/// One parsed `.pin` file: header columns + a raw byte buffer for content
/// past the header.
pub struct PinFile {
    pub columns: Vec<String>,
    pub body: String,
}

impl PinFile {
    pub fn read(path: &Path) -> Result<Self, DiffError> {
        let text = fs::read_to_string(path).map_err(|source| DiffError::Read {
            path: path.display().to_string(),
            source,
        })?;
        let mut lines = text.splitn(2, '\n');
        let header = lines.next().ok_or_else(|| DiffError::EmptyHeader {
            path: path.display().to_string(),
        })?;
        if header.is_empty() {
            return Err(DiffError::EmptyHeader {
                path: path.display().to_string(),
            });
        }
        let columns = header.split('\t').map(|s| s.to_string()).collect();
        let body = lines.next().unwrap_or("").to_string();
        Ok(PinFile { columns, body })
    }
}

/// Compare two header schemas. Returns `Ok(())` if they match, otherwise an
/// error message naming the differing columns.
pub fn compare_schemas(a: &PinFile, b: &PinFile) -> Result<(), String> {
    if a.columns == b.columns {
        return Ok(());
    }
    let only_a: Vec<&str> = a
        .columns
        .iter()
        .filter(|c| !b.columns.contains(c))
        .map(|s| s.as_str())
        .collect();
    let only_b: Vec<&str> = b
        .columns
        .iter()
        .filter(|c| !a.columns.contains(c))
        .map(|s| s.as_str())
        .collect();
    Err(format!(
        "schema differs: only-in-A=[{}] only-in-B=[{}]",
        only_a.join(","),
        only_b.join(",")
    ))
}

pub mod compare;
pub use compare::{compare_with_tolerance, Tolerance};
