//! Per-run MS2 table `<stem>.ms2.parquet`: one row per MS2 scan of a run,
//! reduced to what match-between-runs needs to check a transfer against the
//! run's fragment spectra.
//!
//! Columns: `native_id`, `spec_key` (the PSM SpecId of the scan's rank-1
//! match without its `_<rank>[_<row>]` suffix), `rt_s`, `precursor_mz`,
//! `charge` (0 when unknown), the isolation window as `isolation_target`,
//! `isolation_lower` and `isolation_upper` (null when the source records
//! none), `faims_cv` (null without FAIMS), and the scan's most intense
//! centroids sorted by m/z as two `float32` lists of equal length, `mz` and
//! `intensity`. Rows are in retention-time order.
//!
//! The key-value metadata carries the schema version (`andes.ms2.schema`; a
//! reader rejects another major), the run stem (`andes.ms2.run`), the
//! fragment tolerance of the search (`andes.ms2.fragment_tol`, e.g. `20ppm`
//! or `0.5Da`) and the number of centroids kept per scan
//! (`andes.ms2.top_peaks`). [`write_ms2_table`] writes `<path>.partial` and
//! renames it into place once complete.

use std::fs::File;
use std::path::Path;
use std::sync::Arc;

use arrow::array::{
    Array, ArrayRef, Float32Array, Float64Array, ListArray, StringArray, UInt8Array,
};
use arrow::buffer::{OffsetBuffer, ScalarBuffer};
use arrow::datatypes::{DataType, Field, FieldRef, Schema};
use arrow::record_batch::RecordBatch;
use model::Tolerance;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::ArrowWriter;
use parquet::basic::{Compression, Encoding, ZstdLevel};
use parquet::file::properties::WriterProperties;
use parquet::format::KeyValue;
use parquet::schema::types::ColumnPath;
use quant::ms2::Ms2Scan;

use crate::ms1_table::partial_path;

/// Major version of the table layout. A reader rejects any other major.
pub const MS2_SCHEMA_MAJOR: u32 = 1;

const SCHEMA_KEY: &str = "andes.ms2.schema";
const SCHEMA_VERSION: &str = "1.0";
const RUN_KEY: &str = "andes.ms2.run";
const TOL_KEY: &str = "andes.ms2.fragment_tol";
const TOP_KEY: &str = "andes.ms2.top_peaks";
/// Scans per row group / record batch.
const SCANS_PER_GROUP: usize = 4096;

const COLUMNS: [&str; 11] = [
    "native_id",
    "spec_key",
    "rt_s",
    "precursor_mz",
    "charge",
    "isolation_target",
    "isolation_lower",
    "isolation_upper",
    "faims_cv",
    "mz",
    "intensity",
];

/// A run's MS2 table.
#[derive(Debug, Clone, PartialEq)]
pub struct Ms2Table {
    pub run: String,
    pub fragment_tol: Tolerance,
    pub top_peaks: usize,
    pub scans: Vec<Ms2Scan>,
}

fn item(data_type: DataType) -> FieldRef {
    Arc::new(Field::new("item", data_type, false))
}

fn schema() -> Schema {
    Schema::new(vec![
        Field::new(COLUMNS[0], DataType::Utf8, false),
        Field::new(COLUMNS[1], DataType::Utf8, false),
        Field::new(COLUMNS[2], DataType::Float64, false),
        Field::new(COLUMNS[3], DataType::Float64, false),
        Field::new(COLUMNS[4], DataType::UInt8, false),
        Field::new(COLUMNS[5], DataType::Float64, true),
        Field::new(COLUMNS[6], DataType::Float64, true),
        Field::new(COLUMNS[7], DataType::Float64, true),
        Field::new(COLUMNS[8], DataType::Float32, true),
        Field::new(COLUMNS[9], DataType::List(item(DataType::Float32)), false),
        Field::new(COLUMNS[10], DataType::List(item(DataType::Float32)), false),
    ])
}

/// `20ppm` or `0.5Da`.
pub fn format_tolerance(t: Tolerance) -> String {
    match t {
        Tolerance::Ppm(v) => format!("{v}ppm"),
        Tolerance::Da(v) => format!("{v}Da"),
    }
}

fn parse_tolerance(s: &str) -> Option<Tolerance> {
    let s = s.trim();
    let (v, ppm) = if let Some(v) = s.strip_suffix("ppm") {
        (v, true)
    } else {
        (s.strip_suffix("Da")?, false)
    };
    let v: f64 = v.parse().ok()?;
    (v.is_finite() && v > 0.0).then_some(if ppm {
        Tolerance::Ppm(v)
    } else {
        Tolerance::Da(v)
    })
}

fn properties(version: &str, run: &str, tol: Tolerance, top_peaks: usize) -> WriterProperties {
    let list_values =
        |name: &str| ColumnPath::new(vec![name.to_string(), "list".into(), "item".into()]);
    WriterProperties::builder()
        .set_compression(Compression::ZSTD(ZstdLevel::default()))
        .set_dictionary_enabled(false)
        .set_column_encoding(ColumnPath::from(COLUMNS[2]), Encoding::BYTE_STREAM_SPLIT)
        .set_column_encoding(list_values(COLUMNS[9]), Encoding::BYTE_STREAM_SPLIT)
        .set_column_encoding(list_values(COLUMNS[10]), Encoding::BYTE_STREAM_SPLIT)
        .set_max_row_group_size(SCANS_PER_GROUP)
        .set_key_value_metadata(Some(vec![
            KeyValue::new(SCHEMA_KEY.to_string(), version.to_string()),
            KeyValue::new(RUN_KEY.to_string(), run.to_string()),
            KeyValue::new(TOL_KEY.to_string(), format_tolerance(tol)),
            KeyValue::new(TOP_KEY.to_string(), top_peaks.to_string()),
        ]))
        .build()
}

fn batch(schema: &Arc<Schema>, scans: &[Ms2Scan]) -> Result<RecordBatch, String> {
    let mut offs = vec![0i32];
    let mut end = 0usize;
    for s in scans {
        end += s.peaks.len();
        offs.push(i32::try_from(end).map_err(|_| "too many centroids in one row group")?);
    }
    let offs = OffsetBuffer::new(ScalarBuffer::from(offs));
    let mz = Float32Array::from_iter_values(scans.iter().flat_map(|s| s.peaks.iter().map(|p| p.0)));
    let int =
        Float32Array::from_iter_values(scans.iter().flat_map(|s| s.peaks.iter().map(|p| p.1)));
    let iso = |k: usize| -> Float64Array {
        scans
            .iter()
            .map(|s| s.isolation.map(|w| [w.0, w.1, w.2][k]))
            .collect()
    };
    let cols: Vec<ArrayRef> = vec![
        Arc::new(StringArray::from_iter_values(
            scans.iter().map(|s| s.native_id.as_str()),
        )),
        Arc::new(StringArray::from_iter_values(
            scans.iter().map(|s| s.spec_key.as_str()),
        )),
        Arc::new(Float64Array::from_iter_values(scans.iter().map(|s| s.rt))),
        Arc::new(Float64Array::from_iter_values(
            scans.iter().map(|s| s.precursor_mz),
        )),
        Arc::new(UInt8Array::from_iter_values(scans.iter().map(|s| s.charge))),
        Arc::new(iso(0)),
        Arc::new(iso(1)),
        Arc::new(iso(2)),
        Arc::new(scans.iter().map(|s| s.faims_cv).collect::<Float32Array>()),
        Arc::new(
            ListArray::try_new(item(DataType::Float32), offs.clone(), Arc::new(mz), None)
                .map_err(|e| e.to_string())?,
        ),
        Arc::new(
            ListArray::try_new(item(DataType::Float32), offs, Arc::new(int), None)
                .map_err(|e| e.to_string())?,
        ),
    ];
    RecordBatch::try_new(schema.clone(), cols).map_err(|e| e.to_string())
}

fn write_with_version(path: &Path, t: &Ms2Table, version: &str) -> Result<(), String> {
    let partial = partial_path(path);
    let written = (|| {
        let schema = Arc::new(schema());
        let file = File::create(&partial).map_err(|e| e.to_string())?;
        let props = properties(version, &t.run, t.fragment_tol, t.top_peaks);
        let mut w =
            ArrowWriter::try_new(file, schema.clone(), Some(props)).map_err(|e| e.to_string())?;
        for chunk in t.scans.chunks(SCANS_PER_GROUP) {
            w.write(&batch(&schema, chunk)?)
                .map_err(|e| e.to_string())?;
        }
        w.close().map_err(|e| e.to_string())?;
        std::fs::rename(&partial, path).map_err(|e| e.to_string())
    })()
    .map_err(|e| format!("{}: {e}", path.display()));
    if written.is_err() {
        let _ = std::fs::remove_file(&partial);
    }
    written
}

/// Writes `table` to `path` (scans in the order given). Every error names
/// `path`.
pub fn write_ms2_table(path: &Path, table: &Ms2Table) -> Result<(), String> {
    write_with_version(path, table, SCHEMA_VERSION)
}

fn column<T: 'static>(b: &RecordBatch, i: usize) -> Result<&T, String> {
    b.column(i).as_any().downcast_ref::<T>().ok_or_else(|| {
        format!(
            "column {} has unexpected type {}",
            b.schema().field(i).name(),
            b.column(i).data_type()
        )
    })
}

fn f32_list(b: &RecordBatch, i: usize) -> Result<(&Float32Array, &[i32]), String> {
    let list = column::<ListArray>(b, i)?;
    let values = list
        .values()
        .as_any()
        .downcast_ref::<Float32Array>()
        .ok_or_else(|| format!("column {} items are not float32", COLUMNS[i]))?;
    if list.null_count() != 0 || values.null_count() != 0 {
        return Err(format!("column {} holds nulls", COLUMNS[i]));
    }
    Ok((values, list.value_offsets()))
}

fn opt_f64(a: &Float64Array, i: usize) -> Option<f64> {
    a.is_valid(i).then(|| a.value(i))
}

fn decode(b: &RecordBatch, out: &mut Vec<Ms2Scan>) -> Result<(), String> {
    if b.num_columns() != COLUMNS.len() {
        return Err(format!(
            "expected {} columns, found {}",
            COLUMNS.len(),
            b.num_columns()
        ));
    }
    for (got, want) in b.schema().fields().iter().zip(COLUMNS) {
        if got.name() != want {
            return Err(format!(
                "column {:?} where {want:?} was expected",
                got.name()
            ));
        }
    }
    let ids = column::<StringArray>(b, 0)?;
    let keys = column::<StringArray>(b, 1)?;
    let rts = column::<Float64Array>(b, 2)?;
    let pmz = column::<Float64Array>(b, 3)?;
    let charge = column::<UInt8Array>(b, 4)?;
    let (it, il, iu) = (
        column::<Float64Array>(b, 5)?,
        column::<Float64Array>(b, 6)?,
        column::<Float64Array>(b, 7)?,
    );
    let cv = column::<Float32Array>(b, 8)?;
    if ids.null_count()
        + keys.null_count()
        + rts.null_count()
        + pmz.null_count()
        + charge.null_count()
        != 0
    {
        return Err("a required column holds nulls".to_string());
    }
    let (mz, mz_off) = f32_list(b, 9)?;
    let (int, int_off) = f32_list(b, 10)?;
    for i in 0..b.num_rows() {
        let (lo, hi) = (mz_off[i] as usize, mz_off[i + 1] as usize);
        let (ilo, ihi) = (int_off[i] as usize, int_off[i + 1] as usize);
        if hi - lo != ihi - ilo {
            return Err(format!(
                "scan {:?}: {} m/z values but {} intensities",
                ids.value(i),
                hi - lo,
                ihi - ilo
            ));
        }
        let isolation = match (opt_f64(it, i), opt_f64(il, i), opt_f64(iu, i)) {
            (Some(t), Some(l), Some(u)) => Some((t, l, u)),
            _ => None,
        };
        out.push(Ms2Scan {
            native_id: ids.value(i).to_string(),
            spec_key: keys.value(i).to_string(),
            rt: rts.value(i),
            precursor_mz: pmz.value(i),
            charge: charge.value(i),
            isolation,
            faims_cv: cv.is_valid(i).then(|| cv.value(i)),
            peaks: mz.values()[lo..hi]
                .iter()
                .copied()
                .zip(int.values()[ilo..ihi].iter().copied())
                .collect(),
        });
    }
    Ok(())
}

fn read(path: &Path) -> Result<Ms2Table, String> {
    let file = File::open(path).map_err(|e| e.to_string())?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file).map_err(|e| e.to_string())?;
    let meta = builder
        .metadata()
        .file_metadata()
        .key_value_metadata()
        .cloned()
        .unwrap_or_default();
    let get = |k: &str| {
        meta.iter()
            .find(|kv| kv.key == k)
            .and_then(|kv| kv.value.clone())
    };
    let version = get(SCHEMA_KEY)
        .ok_or_else(|| format!("missing metadata key {SCHEMA_KEY}; not an andes MS2 table"))?;
    let major: u32 = version
        .split('.')
        .next()
        .and_then(|m| m.parse().ok())
        .ok_or_else(|| format!("unparseable schema version {version:?}"))?;
    if major != MS2_SCHEMA_MAJOR {
        return Err(format!(
            "unsupported MS2 table schema {version} (this build reads major {MS2_SCHEMA_MAJOR})"
        ));
    }
    let run = get(RUN_KEY).ok_or_else(|| format!("missing metadata key {RUN_KEY}"))?;
    let tol_text = get(TOL_KEY).ok_or_else(|| format!("missing metadata key {TOL_KEY}"))?;
    let fragment_tol = parse_tolerance(&tol_text)
        .ok_or_else(|| format!("unreadable fragment tolerance {tol_text:?}"))?;
    let top_peaks = get(TOP_KEY)
        .and_then(|v| v.parse().ok())
        .ok_or_else(|| format!("missing or unreadable metadata key {TOP_KEY}"))?;
    let n = usize::try_from(builder.metadata().file_metadata().num_rows()).unwrap_or(0);
    let reader = builder
        .with_batch_size(SCANS_PER_GROUP)
        .build()
        .map_err(|e| e.to_string())?;
    let mut scans = Vec::with_capacity(n.min(1 << 21));
    for b in reader {
        decode(&b.map_err(|e| e.to_string())?, &mut scans)?;
    }
    Ok(Ms2Table {
        run,
        fragment_tol,
        top_peaks,
        scans,
    })
}

/// Reads the MS2 table at `path`. Every error names `path`.
pub fn read_ms2_table(path: &Path) -> Result<Ms2Table, String> {
    read(path).map_err(|e| format!("{}: {e}", path.display()))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn scan(i: usize) -> Ms2Scan {
        Ms2Scan {
            native_id: format!("controllerType=0 controllerNumber=1 scan={i}"),
            spec_key: format!("run/scan={i}_{i}"),
            rt: 10.0 + i as f64 * 0.1234567,
            precursor_mz: 500.0 + i as f64 * 1.000_123_456,
            charge: (i % 4) as u8,
            isolation: (i % 3 != 0).then_some((500.0 + i as f64, 0.7, 0.65)),
            faims_cv: (i % 2 == 0).then_some(-45.0),
            peaks: (0..(i % 7))
                .map(|k| (100.0 + k as f32 * 37.123_4, (k * i) as f32 + 0.5))
                .collect(),
        }
    }

    fn table(n: usize) -> Ms2Table {
        Ms2Table {
            run: "run".into(),
            fragment_tol: Tolerance::Ppm(20.0),
            top_peaks: 100,
            scans: (0..n).map(scan).collect(),
        }
    }

    #[test]
    fn round_trip_is_exact_over_row_groups() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("run.ms2.parquet");
        let t = table(SCANS_PER_GROUP + 13);
        write_ms2_table(&p, &t).unwrap();
        assert!(!partial_path(&p).exists());
        let back = read_ms2_table(&p).unwrap();
        assert_eq!(back, t);
        let groups = ParquetRecordBatchReaderBuilder::try_new(File::open(&p).unwrap())
            .unwrap()
            .metadata()
            .num_row_groups();
        assert_eq!(groups, 2);

        let e = dir.path().join("empty.ms2.parquet");
        let empty = Ms2Table {
            fragment_tol: Tolerance::Da(0.5),
            ..table(0)
        };
        write_ms2_table(&e, &empty).unwrap();
        assert_eq!(read_ms2_table(&e).unwrap(), empty);
    }

    #[test]
    fn tolerances_round_trip() {
        for t in [
            Tolerance::Ppm(20.0),
            Tolerance::Ppm(7.5),
            Tolerance::Da(0.02),
        ] {
            assert_eq!(parse_tolerance(&format_tolerance(t)), Some(t));
        }
        for bad in ["", "20", "ppm", "-1ppm", "0Da", "xDa"] {
            assert_eq!(parse_tolerance(bad), None, "{bad}");
        }
    }

    #[test]
    fn rejects_other_majors_foreign_and_renamed_tables_naming_the_file() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("v2.ms2.parquet");
        write_with_version(&p, &table(3), "2.0").unwrap();
        let e = read_ms2_table(&p).unwrap_err();
        assert!(
            e.contains("v2.ms2.parquet") && e.contains("schema 2.0"),
            "{e}"
        );

        let ms1 = dir.path().join("x.ms1.parquet");
        crate::ms1_table::write_ms1_table(&ms1, "x", &quant::ms1_index::Ms1RunIndex::default())
            .unwrap();
        let e = read_ms2_table(&ms1).unwrap_err();
        assert!(e.contains("x.ms1.parquet") && e.contains(SCHEMA_KEY), "{e}");

        let junk = dir.path().join("junk.ms2.parquet");
        std::fs::write(&junk, b"PAR1").unwrap();
        assert!(read_ms2_table(&junk)
            .unwrap_err()
            .contains("junk.ms2.parquet"));
        let e = read_ms2_table(&dir.path().join("nope.ms2.parquet")).unwrap_err();
        assert!(e.contains("nope.ms2.parquet"), "{e}");

        // A renamed column.
        let r = dir.path().join("r.ms2.parquet");
        let t = table(4);
        let good = batch(&Arc::new(schema()), &t.scans).unwrap();
        let mut fields: Vec<Field> = schema()
            .fields()
            .iter()
            .map(|f| f.as_ref().clone())
            .collect();
        fields[2] = Field::new("rt_min", DataType::Float64, false);
        let renamed = Arc::new(Schema::new(fields));
        let b = RecordBatch::try_new(renamed.clone(), good.columns().to_vec()).unwrap();
        let mut w = ArrowWriter::try_new(
            File::create(&r).unwrap(),
            renamed,
            Some(properties("1.0", "r", Tolerance::Ppm(20.0), 100)),
        )
        .unwrap();
        w.write(&b).unwrap();
        w.close().unwrap();
        let e = read_ms2_table(&r).unwrap_err();
        assert!(e.contains("r.ms2.parquet") && e.contains("rt_min"), "{e}");
    }

    #[test]
    fn failed_write_leaves_no_file() {
        let dir = tempfile::tempdir().unwrap();
        let bad = dir.path().join("missing").join("d.ms2.parquet");
        let e = write_ms2_table(&bad, &table(2)).unwrap_err();
        assert!(e.contains("d.ms2.parquet"), "{e}");
        assert!(!bad.exists() && !partial_path(&bad).exists());
    }
}
