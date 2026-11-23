//! The two tables `andes database` writes beside a target + decoy FASTA `<stem>.fasta`.
//!
//! - `<stem>.database.parquet` (`andes.database.schema`), the settings sidecar: one row per
//!   digest the decoys were checked for (`role` is `build` for the digest they were built for,
//!   `check` for a further digest of the cohort), with columns `role`, `enzymes`,
//!   `specificity`, `min_length`, `max_length` and `max_missed_cleavages`. The key-value
//!   metadata holds the rest of the decoy key (`decoy_strategy`, `decoy_seed`, `decoy_prefix`,
//!   `decoy_suffix`, `targets_sha256`), the digest of the FASTA bytes (`fasta_sha256`), the key
//!   itself (`decoy_key`), the entry counts (`targets`, `decoys`) and `andes_version`. A search
//!   given `<stem>.fasta` reads it to check that it searches with the same key.
//! - `<stem>.overlap.parquet` (`andes.decoy_overlap.schema`), the overlap report: one row per
//!   (digest, peptide length) with the digest columns, `length`, and the distinct target
//!   peptides, decoy peptides and decoy peptides equal to a target peptide (I and L folded,
//!   `target_peptides`, `decoy_peptides`, `shared_peptides`). Metadata: the decoy peptides the
//!   check redrew (`redrawn`) and the ones still equal to a target (`left`).

use std::path::{Path, PathBuf};
use std::sync::Arc;

use arrow::array::{ArrayRef, StringArray, UInt32Array, UInt64Array};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;

use crate::versioned::{column, read_table, write_table, TableVersion};

const SETTINGS_VERSION: &str = "1.0";
const SETTINGS: TableVersion = TableVersion {
    key: "andes.database.schema",
    version: SETTINGS_VERSION,
    major: 1,
    what: "decoy database settings table",
};

const OVERLAP_VERSION: &str = "1.0";
const OVERLAP: TableVersion = TableVersion {
    key: "andes.decoy_overlap.schema",
    version: OVERLAP_VERSION,
    major: 1,
    what: "decoy overlap table",
};

/// `role` of the digest the decoys were built for.
pub const ROLE_BUILD: &str = "build";
/// `role` of a further digest the decoys were checked against.
pub const ROLE_CHECK: &str = "check";

/// One digest: the enzymes (primary first, then the others as the search orders them), the
/// specificity (`fully`, `semi` or `non-specific`), the peptide lengths and the missed
/// cleavages.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Digest {
    pub enzymes: String,
    pub specificity: String,
    pub min_length: u32,
    pub max_length: u32,
    pub max_missed_cleavages: u32,
}

/// The settings sidecar of a decoy database.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DatabaseSettings {
    /// The build digest first, then the further digests.
    pub digests: Vec<Digest>,
    pub decoy_strategy: String,
    pub decoy_seed: u64,
    pub decoy_prefix: String,
    pub decoy_suffix: String,
    pub targets_sha256: String,
    pub fasta_sha256: String,
    pub decoy_key: String,
    pub targets: u64,
    pub decoys: u64,
}

/// One row of the overlap report.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OverlapRow {
    pub digest: Digest,
    pub length: u32,
    pub target_peptides: u64,
    pub decoy_peptides: u64,
    pub shared_peptides: u64,
}

fn sibling(fasta: &Path, suffix: &str) -> PathBuf {
    let name = fasta
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    let name = name.strip_suffix(".gz").unwrap_or(&name);
    let stem = [".fasta", ".fa", ".faa"]
        .iter()
        .find_map(|ext| name.strip_suffix(ext))
        .unwrap_or(name);
    fasta.with_file_name(format!("{stem}{suffix}"))
}

/// `<stem>.database.parquet` beside the FASTA `<stem>.fasta` (also `.fa`, `.faa`, gzipped).
pub fn settings_path(fasta: &Path) -> PathBuf {
    sibling(fasta, ".database.parquet")
}

/// `<stem>.overlap.parquet` beside the FASTA `<stem>.fasta`.
pub fn overlap_path(fasta: &Path) -> PathBuf {
    sibling(fasta, ".overlap.parquet")
}

fn digest_fields() -> Vec<Field> {
    vec![
        Field::new("enzymes", DataType::Utf8, false),
        Field::new("specificity", DataType::Utf8, false),
        Field::new("min_length", DataType::UInt32, false),
        Field::new("max_length", DataType::UInt32, false),
        Field::new("max_missed_cleavages", DataType::UInt32, false),
    ]
}

fn digest_columns<'a>(digests: impl Iterator<Item = &'a Digest> + Clone) -> Vec<ArrayRef> {
    vec![
        Arc::new(StringArray::from_iter_values(
            digests.clone().map(|d| d.enzymes.as_str()),
        )),
        Arc::new(StringArray::from_iter_values(
            digests.clone().map(|d| d.specificity.as_str()),
        )),
        Arc::new(UInt32Array::from_iter_values(
            digests.clone().map(|d| d.min_length),
        )),
        Arc::new(UInt32Array::from_iter_values(
            digests.clone().map(|d| d.max_length),
        )),
        Arc::new(UInt32Array::from_iter_values(
            digests.map(|d| d.max_missed_cleavages),
        )),
    ]
}

fn read_digests(b: &RecordBatch, path: &Path) -> Result<Vec<Digest>, String> {
    let enzymes = column::<StringArray>(b, "enzymes", path)?;
    let specificity = column::<StringArray>(b, "specificity", path)?;
    let min = column::<UInt32Array>(b, "min_length", path)?;
    let max = column::<UInt32Array>(b, "max_length", path)?;
    let missed = column::<UInt32Array>(b, "max_missed_cleavages", path)?;
    Ok((0..b.num_rows())
        .map(|i| Digest {
            enzymes: enzymes.value(i).to_string(),
            specificity: specificity.value(i).to_string(),
            min_length: min.value(i),
            max_length: max.value(i),
            max_missed_cleavages: missed.value(i),
        })
        .collect())
}

/// Writes the settings sidecar.
pub fn write_settings(path: &Path, s: &DatabaseSettings) -> Result<(), String> {
    let mut fields = vec![Field::new("role", DataType::Utf8, false)];
    fields.extend(digest_fields());
    let schema = Arc::new(Schema::new(fields));
    let roles = (0..s.digests.len()).map(|i| if i == 0 { ROLE_BUILD } else { ROLE_CHECK });
    let mut cols: Vec<ArrayRef> = vec![Arc::new(StringArray::from_iter_values(roles))];
    cols.extend(digest_columns(s.digests.iter()));
    let batch = RecordBatch::try_new(schema.clone(), cols).map_err(|e| e.to_string());
    let meta = [
        ("decoy_strategy", s.decoy_strategy.clone()),
        ("decoy_seed", s.decoy_seed.to_string()),
        ("decoy_prefix", s.decoy_prefix.clone()),
        ("decoy_suffix", s.decoy_suffix.clone()),
        ("targets_sha256", s.targets_sha256.clone()),
        ("fasta_sha256", s.fasta_sha256.clone()),
        ("decoy_key", s.decoy_key.clone()),
        ("targets", s.targets.to_string()),
        ("decoys", s.decoys.to_string()),
        ("andes_version", env!("CARGO_PKG_VERSION").to_string()),
    ];
    write_table(
        path,
        &SETTINGS,
        SETTINGS_VERSION,
        schema,
        std::iter::once(batch),
        &meta,
    )
}

/// Reads a settings sidecar written by [`write_settings`], checking its schema version.
pub fn read_settings(path: &Path) -> Result<DatabaseSettings, String> {
    let t = read_table(path, &SETTINGS)?;
    let mut digests = Vec::new();
    let mut roles = Vec::new();
    for b in &t.batches {
        let role = column::<StringArray>(b, "role", path)?;
        roles.extend((0..b.num_rows()).map(|i| role.value(i).to_string()));
        digests.extend(read_digests(b, path)?);
    }
    let first_is_build = roles.first().map(String::as_str) == Some(ROLE_BUILD);
    let rest_are_checks = roles.iter().skip(1).all(|r| r == ROLE_CHECK);
    if !first_is_build || !rest_are_checks {
        return Err(format!(
            "{}: the first digest must be the build digest and the others checks, got roles {roles:?}",
            path.display()
        ));
    }
    let get = |k: &str| {
        t.meta
            .get(k)
            .cloned()
            .ok_or_else(|| format!("{}: missing metadata key {k}", path.display()))
    };
    let number = |k: &str| -> Result<u64, String> {
        get(k)?
            .parse()
            .map_err(|e| format!("{}: metadata {k}: {e}", path.display()))
    };
    Ok(DatabaseSettings {
        digests,
        decoy_strategy: get("decoy_strategy")?,
        decoy_seed: number("decoy_seed")?,
        decoy_prefix: get("decoy_prefix")?,
        decoy_suffix: get("decoy_suffix")?,
        targets_sha256: get("targets_sha256")?,
        fasta_sha256: get("fasta_sha256")?,
        decoy_key: get("decoy_key")?,
        targets: number("targets")?,
        decoys: number("decoys")?,
    })
}

/// Writes the overlap report; `redrawn` and `left` are the totals of the peptide-level check.
pub fn write_overlap(
    path: &Path,
    rows: &[OverlapRow],
    redrawn: usize,
    left: usize,
) -> Result<(), String> {
    let mut fields = digest_fields();
    fields.extend([
        Field::new("length", DataType::UInt32, false),
        Field::new("target_peptides", DataType::UInt64, false),
        Field::new("decoy_peptides", DataType::UInt64, false),
        Field::new("shared_peptides", DataType::UInt64, false),
    ]);
    let schema = Arc::new(Schema::new(fields));
    let mut cols = digest_columns(rows.iter().map(|r| &r.digest));
    cols.push(Arc::new(UInt32Array::from_iter_values(
        rows.iter().map(|r| r.length),
    )));
    let counts: [fn(&OverlapRow) -> u64; 3] = [
        |r| r.target_peptides,
        |r| r.decoy_peptides,
        |r| r.shared_peptides,
    ];
    for f in counts {
        cols.push(Arc::new(UInt64Array::from_iter_values(rows.iter().map(f))));
    }
    let batch = RecordBatch::try_new(schema.clone(), cols).map_err(|e| e.to_string());
    let meta = [("redrawn", redrawn.to_string()), ("left", left.to_string())];
    write_table(
        path,
        &OVERLAP,
        OVERLAP_VERSION,
        schema,
        std::iter::once(batch),
        &meta,
    )
}

/// Reads an overlap report written by [`write_overlap`].
pub fn read_overlap(path: &Path) -> Result<Vec<OverlapRow>, String> {
    let t = read_table(path, &OVERLAP)?;
    let mut out = Vec::new();
    for b in &t.batches {
        let digests = read_digests(b, path)?;
        let length = column::<UInt32Array>(b, "length", path)?;
        let tp = column::<UInt64Array>(b, "target_peptides", path)?;
        let dp = column::<UInt64Array>(b, "decoy_peptides", path)?;
        let sp = column::<UInt64Array>(b, "shared_peptides", path)?;
        for (i, digest) in digests.into_iter().enumerate() {
            out.push(OverlapRow {
                digest,
                length: length.value(i),
                target_peptides: tp.value(i),
                decoy_peptides: dp.value(i),
                shared_peptides: sp.value(i),
            });
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn digest(enzymes: &str, missed: u32) -> Digest {
        Digest {
            enzymes: enzymes.into(),
            specificity: "fully".into(),
            min_length: 7,
            max_length: 30,
            max_missed_cleavages: missed,
        }
    }

    #[test]
    fn sidecar_names_follow_the_fasta_stem() {
        let p = Path::new("/d/td.fasta");
        assert_eq!(settings_path(p), Path::new("/d/td.database.parquet"));
        assert_eq!(overlap_path(p), Path::new("/d/td.overlap.parquet"));
        assert_eq!(
            settings_path(Path::new("x.fa.gz")),
            Path::new("x.database.parquet")
        );
        assert_eq!(
            settings_path(Path::new("x.db")),
            Path::new("x.db.database.parquet")
        );
    }

    #[test]
    fn settings_round_trip() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("td.database.parquet");
        let s = DatabaseSettings {
            digests: vec![digest("Trypsin", 2), digest("LysC", 1)],
            decoy_strategy: "reverse".into(),
            decoy_seed: 42,
            decoy_prefix: "XXX".into(),
            decoy_suffix: String::new(),
            targets_sha256: "ab".into(),
            fasta_sha256: "cd".into(),
            decoy_key: "ef".into(),
            targets: 3,
            decoys: 3,
        };
        write_settings(&path, &s).unwrap();
        assert_eq!(read_settings(&path).unwrap(), s);
    }

    #[test]
    fn overlap_round_trip() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("td.overlap.parquet");
        let rows = vec![
            OverlapRow {
                digest: digest("Trypsin", 2),
                length: 7,
                target_peptides: 10,
                decoy_peptides: 9,
                shared_peptides: 1,
            },
            OverlapRow {
                digest: digest("LysC", 1),
                length: 8,
                target_peptides: 5,
                decoy_peptides: 5,
                shared_peptides: 0,
            },
        ];
        write_overlap(&path, &rows, 4, 1).unwrap();
        assert_eq!(read_overlap(&path).unwrap(), rows);
    }

    #[test]
    fn another_table_is_refused() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("td.overlap.parquet");
        write_overlap(&path, &[], 0, 0).unwrap();
        let e = read_settings(&path).unwrap_err();
        assert!(e.contains("andes.database.schema"), "{e}");
    }
}
