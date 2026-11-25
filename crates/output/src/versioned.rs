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
use parquet::file::metadata::KeyValue;
use parquet::file::properties::{WriterProperties, WriterPropertiesBuilder};

use crate::atomic::partial_path;

pub(crate) const BATCH_ROWS: usize = 65_536;

pub(crate) struct TableVersion {
    pub key: &'static str,
    pub version: &'static str,
    pub major: u32,
    pub what: &'static str,
}

/// Writer properties shared by every ZSTD table: compression, row-group size and the
/// key-value footer. Returns the builder so a caller can add column encodings.
pub(crate) fn zstd_props(row_group_rows: usize, meta: Vec<KeyValue>) -> WriterPropertiesBuilder {
    WriterProperties::builder()
        .set_compression(Compression::ZSTD(ZstdLevel::default()))
        .set_max_row_group_row_count(Some(row_group_rows))
        .set_key_value_metadata(Some(meta))
}

pub(crate) struct TableWriter {
    path: PathBuf,
    partial: PathBuf,
    writer: Option<ArrowWriter<File>>,
}

impl TableWriter {
    pub(crate) fn create(
        path: &Path,
        table: &TableVersion,
        version: &str,
        schema: Arc<Schema>,
        extra: &[(&str, String)],
        row_group_rows: usize,
    ) -> Result<Self, String> {
        let partial = partial_path(path);
        let meta: Vec<KeyValue> =
            std::iter::once(KeyValue::new(table.key.to_string(), version.to_string()))
                .chain(
                    extra
                        .iter()
                        .map(|(k, v)| KeyValue::new(k.to_string(), v.clone())),
                )
                .collect();
        let props = zstd_props(row_group_rows, meta).build();
        let err = |e: &dyn std::fmt::Display| format!("{}: {e}", partial.display());
        let file = File::create(&partial).map_err(|e| err(&e))?;
        let writer = ArrowWriter::try_new(file, schema, Some(props));
        let writer = match writer {
            Ok(w) => w,
            Err(e) => {
                let _ = std::fs::remove_file(&partial);
                return Err(err(&e));
            }
        };
        Ok(Self {
            path: path.to_path_buf(),
            partial,
            writer: Some(writer),
        })
    }

    pub(crate) fn partial(&self) -> &Path {
        &self.partial
    }

    pub(crate) fn write(&mut self, b: &RecordBatch) -> Result<(), String> {
        let w = self.writer.as_mut().expect("writer is open until finish");
        w.write(b)
            .map_err(|e| format!("{}: {e}", self.partial.display()))
    }

    pub(crate) fn finish(mut self) -> Result<(), String> {
        let w = self.writer.take().expect("writer is open until finish");
        let done = w
            .close()
            .map_err(|e| format!("{}: {e}", self.partial.display()))
            .and_then(|_| {
                std::fs::rename(&self.partial, &self.path)
                    .map_err(|e| format!("{}: {e}", self.path.display()))
            });
        if done.is_err() {
            let _ = std::fs::remove_file(&self.partial);
        }
        done
    }
}

impl Drop for TableWriter {
    fn drop(&mut self) {
        if self.writer.take().is_some() {
            let _ = std::fs::remove_file(&self.partial);
        }
    }
}

pub(crate) fn write_table(
    path: &Path,
    table: &TableVersion,
    version: &str,
    schema: Arc<Schema>,
    batches: impl IntoIterator<Item = Result<RecordBatch, String>>,
    extra: &[(&str, String)],
) -> Result<(), String> {
    let mut w = TableWriter::create(path, table, version, schema, extra, BATCH_ROWS)?;
    for b in batches {
        let b = b.map_err(|e| format!("{}: {e}", path.display()))?;
        w.write(&b)?;
    }
    w.finish()
}

pub struct ReadTable {
    pub batches: Vec<RecordBatch>,
    pub meta: HashMap<String, String>,
}

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

/// Typed column by position. The error names the file when a path is given and is
/// otherwise the bare column message used by tables read from a stream.
fn downcast_column<'a, T: 'static>(
    b: &'a RecordBatch,
    i: usize,
    path: Option<&Path>,
) -> Result<&'a T, String> {
    b.column(i).as_any().downcast_ref::<T>().ok_or_else(|| {
        let what = format!(
            "column {} has unexpected type {}",
            b.schema().field(i).name(),
            b.column(i).data_type()
        );
        match path {
            Some(p) => format!("{}: {what}", p.display()),
            None => what,
        }
    })
}

pub(crate) fn column_at<'a, T: 'static>(
    b: &'a RecordBatch,
    i: usize,
    path: &Path,
) -> Result<&'a T, String> {
    downcast_column(b, i, Some(path))
}

pub(crate) fn column_at_unpathed<T: 'static>(b: &RecordBatch, i: usize) -> Result<&T, String> {
    downcast_column(b, i, None)
}

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
