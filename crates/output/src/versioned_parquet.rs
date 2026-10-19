//! Writer and reader shared by the final quantification tables
//! (`<stem>.tmt.parquet`, `<stem>.lfq.parquet`, `<stem>.mbr.parquet`,
//! `<stem>.mbr_validation.parquet`).
//!
//! Each table records its schema version under its own key-value metadata key
//! (`andes.<table>.schema`, `<major>.<minor>`); a reader rejects any other
//! major. A minor version only adds columns, so readers look columns up by
//! name. Files are written to `<path>.partial` and renamed into place, ZSTD
//! compressed.

use std::collections::HashMap;
use std::fs::File;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use arrow::array::Array;
use arrow::datatypes::Schema;
use arrow::record_batch::RecordBatch;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::ArrowWriter;
use parquet::basic::{Compression, ZstdLevel};
use parquet::file::properties::WriterProperties;
use parquet::format::KeyValue;

/// Rows per Parquet row group.
pub(crate) const BATCH_ROWS: usize = 65_536;

/// Identity of one versioned table: its metadata key, the version this build
/// writes and the major it reads, and a name for error messages.
pub(crate) struct TableVersion {
    pub key: &'static str,
    pub version: &'static str,
    pub major: u32,
    pub what: &'static str,
}

/// Writes `batches` (all with `schema`) to `path` through `<path>.partial`,
/// with the table's schema version and `extra` as key-value metadata. Every
/// error names the file.
pub(crate) fn write_table(
    path: &Path,
    table: &TableVersion,
    version: &str,
    schema: Arc<Schema>,
    batches: impl IntoIterator<Item = Result<RecordBatch, String>>,
    extra: &[(&str, String)],
) -> Result<(), String> {
    let mut partial = path.as_os_str().to_owned();
    partial.push(".partial");
    let partial = PathBuf::from(partial);
    let err = |p: &Path, e: &dyn std::fmt::Display| format!("{}: {e}", p.display());
    let meta: Vec<KeyValue> =
        std::iter::once(KeyValue::new(table.key.to_string(), version.to_string()))
            .chain(
                extra
                    .iter()
                    .map(|(k, v)| KeyValue::new(k.to_string(), v.clone())),
            )
            .collect();
    let props = WriterProperties::builder()
        .set_compression(Compression::ZSTD(ZstdLevel::default()))
        .set_max_row_group_size(BATCH_ROWS)
        .set_key_value_metadata(Some(meta))
        .build();
    let written = (|| {
        let file = File::create(&partial).map_err(|e| err(&partial, &e))?;
        let mut writer =
            ArrowWriter::try_new(file, schema, Some(props)).map_err(|e| err(&partial, &e))?;
        for b in batches {
            let b = b.map_err(|e| err(path, &e))?;
            writer.write(&b).map_err(|e| err(&partial, &e))?;
        }
        writer.close().map_err(|e| err(&partial, &e))?;
        std::fs::rename(&partial, path).map_err(|e| err(path, &e))
    })();
    if written.is_err() {
        let _ = std::fs::remove_file(&partial);
    }
    written
}

/// A table read back: its record batches and key-value metadata.
pub(crate) struct ReadTable {
    pub batches: Vec<RecordBatch>,
    pub meta: HashMap<String, String>,
}

/// Reads every batch of `path` after checking its schema version against
/// `table`. Every error names the file.
pub(crate) fn read_table(path: &Path, table: &TableVersion) -> Result<ReadTable, String> {
    let err = |e: &dyn std::fmt::Display| format!("{}: {e}", path.display());
    let file = File::open(path).map_err(|e| err(&e))?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file).map_err(|e| err(&e))?;
    let meta: HashMap<String, String> = builder
        .metadata()
        .file_metadata()
        .key_value_metadata()
        .cloned()
        .unwrap_or_default()
        .into_iter()
        .filter_map(|kv| kv.value.map(|v| (kv.key, v)))
        .collect();
    let version = meta.get(table.key).ok_or_else(|| {
        err(&format!(
            "missing metadata key {}; not an andes {}",
            table.key, table.what
        ))
    })?;
    let major: u32 = version
        .split('.')
        .next()
        .and_then(|m| m.parse().ok())
        .ok_or_else(|| err(&format!("unparseable schema version {version:?}")))?;
    if major != table.major {
        return Err(err(&format!(
            "unsupported {} schema {version} (this build reads major {})",
            table.what, table.major
        )));
    }
    let batches = builder
        .build()
        .map_err(|e| err(&e))?
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| err(&e))?;
    Ok(ReadTable { batches, meta })
}

/// Column `name` of `b` as array type `T`. The error names the file.
pub(crate) fn column<'a, T: 'static>(
    b: &'a RecordBatch,
    name: &str,
    path: &Path,
) -> Result<&'a T, String> {
    let c = b
        .column_by_name(name)
        .ok_or_else(|| format!("{}: missing column {name}", path.display()))?;
    c.as_any().downcast_ref::<T>().ok_or_else(|| {
        format!(
            "{}: column {name} has unexpected type {}",
            path.display(),
            c.data_type()
        )
    })
}

/// Value `i` of a nullable primitive column.
pub(crate) fn opt<T: arrow::datatypes::ArrowPrimitiveType>(
    a: &arrow::array::PrimitiveArray<T>,
    i: usize,
) -> Option<T::Native> {
    (!a.is_null(i)).then(|| a.value(i))
}

#[cfg(test)]
mod tests {
    use super::*;
    use arrow::array::Float64Array;
    use arrow::datatypes::{DataType, Field};

    const T: TableVersion = TableVersion {
        key: "andes.test.schema",
        version: "1.0",
        major: 1,
        what: "test table",
    };

    fn one(path: &Path, version: &str) -> Result<(), String> {
        let schema = Arc::new(Schema::new(vec![Field::new("x", DataType::Float64, true)]));
        let b = RecordBatch::try_new(
            schema.clone(),
            vec![Arc::new(Float64Array::from(vec![Some(1.5), None]))],
        )
        .map_err(|e| e.to_string());
        write_table(
            path,
            &T,
            version,
            schema,
            [b],
            &[("andes.test.tag", "v".into())],
        )
    }

    #[test]
    fn round_trip_keeps_rows_nulls_and_metadata() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("t.parquet");
        one(&p, T.version).unwrap();
        assert!(!dir.path().join("t.parquet.partial").exists());
        let t = read_table(&p, &T).unwrap();
        assert_eq!(t.meta.get("andes.test.tag").map(String::as_str), Some("v"));
        let x: &Float64Array = column(&t.batches[0], "x", &p).unwrap();
        assert_eq!((opt(x, 0), opt(x, 1)), (Some(1.5), None));
        assert!(column::<Float64Array>(&t.batches[0], "y", &p)
            .unwrap_err()
            .contains("missing column y"));
    }

    #[test]
    fn rejects_other_majors_and_foreign_files_naming_the_file() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("t.parquet");
        one(&p, "1.7").unwrap();
        assert!(read_table(&p, &T).is_ok());
        one(&p, "2.0").unwrap();
        let e = read_table(&p, &T).err().unwrap();
        assert!(e.contains("t.parquet") && e.contains("schema 2.0"), "{e}");
        let other = TableVersion {
            key: "andes.other.schema",
            ..T
        };
        let e = read_table(&p, &other).err().unwrap();
        assert!(e.contains("andes.other.schema"), "{e}");
    }
}
