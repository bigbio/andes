//! Protein database loading: FASTA or UniProt XML, plain or gzipped, chosen
//! by file extension.

use std::path::Path;

use model::ProteinDb;

use crate::fasta::{FastaParseError, FastaReader};
use crate::uniprot::{load_uniprot_xml, UniProtLoadStats, UniProtOptions, UniProtParseError};

/// Format of a protein database file.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DatabaseFormat {
    Fasta,
    UniProtXml,
}

impl DatabaseFormat {
    /// `.xml` (or `.xml.gz`, case-insensitive) is UniProt XML; anything else
    /// is read as FASTA.
    pub fn of_path(path: &Path) -> Self {
        let name = path
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or_default()
            .to_ascii_lowercase();
        let name = name.strip_suffix(".gz").unwrap_or(&name);
        if name.ends_with(".xml") {
            Self::UniProtXml
        } else {
            Self::Fasta
        }
    }
}

/// A loaded target database.
#[derive(Debug, Clone)]
pub struct LoadedDatabase {
    pub db: ProteinDb,
    pub format: DatabaseFormat,
    /// Present for a UniProt XML database.
    pub uniprot: Option<UniProtLoadStats>,
}

#[derive(thiserror::Error, Debug)]
pub enum DatabaseError {
    #[error("opening {path}: {source}")]
    Open {
        path: String,
        #[source]
        source: std::io::Error,
    },
    #[error(transparent)]
    Fasta(#[from] FastaParseError),
    #[error(transparent)]
    UniProt(#[from] UniProtParseError),
}

/// Load the target database at `path` (gzip is decompressed transparently).
/// `uniprot` only applies to a UniProt XML database.
pub fn load_protein_database(
    path: &Path,
    uniprot: &UniProtOptions,
) -> Result<LoadedDatabase, DatabaseError> {
    let reader = crate::open_buf_maybe_gz(path).map_err(|source| DatabaseError::Open {
        path: path.display().to_string(),
        source,
    })?;
    let format = DatabaseFormat::of_path(path);
    Ok(match format {
        DatabaseFormat::Fasta => LoadedDatabase {
            db: FastaReader::load_all(reader)?,
            format,
            uniprot: None,
        },
        DatabaseFormat::UniProtXml => {
            let (db, stats) = load_uniprot_xml(reader, uniprot)?;
            LoadedDatabase {
                db,
                format,
                uniprot: Some(stats),
            }
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn format_follows_the_extension_through_gz() {
        let f = |p: &str| DatabaseFormat::of_path(Path::new(p));
        assert_eq!(f("human.fasta"), DatabaseFormat::Fasta);
        assert_eq!(f("human.fasta.gz"), DatabaseFormat::Fasta);
        assert_eq!(f("uniprot_sprot.xml.gz"), DatabaseFormat::UniProtXml);
        assert_eq!(f("dir.xml/UP000005640.XML"), DatabaseFormat::UniProtXml);
        assert_eq!(f("proteome.fa"), DatabaseFormat::Fasta);
    }
}
