//! Per-run MS1 table `<stem>.ms1.parquet`: one row per MS1 survey scan of a
//! run, in the order of its [`Ms1RunIndex`] (retention time). It is the run's
//! MS1 for the quantification stages, both the persistent `--quant-dir`
//! artifact and the in-process scratch copy kept while the search runs.
//!
//! Columns: `native_id` (empty when the source has none), `rt_s`, `faims_cv`
//! (null without FAIMS), and the scan's m/z-sorted centroids as two lists of
//! equal length, `mz` and `intensity`, and, for a run whose survey scans are
//! resolved by ion mobility, a last column `mobility` (the 1/K0 of each
//! centroid; schema 1.1). The values are stored bit for bit, so a
//! round trip rebuilds the same index. The key-value metadata carries the
//! schema version (`andes.ms1.schema`; a reader rejects another major) and the
//! run stem (`andes.ms1.run`). A row group holds a bounded number of scans, so
//! writing and reading stream; [`write_ms1_table`] writes `<path>.partial` and
//! renames it into place once complete.

use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use arrow::array::{Array, ArrayRef, Float32Array, Float64Array, ListArray, StringArray};
use arrow::buffer::{OffsetBuffer, ScalarBuffer};
use arrow::datatypes::{DataType, Field, FieldRef, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::{ArrowWriter, ProjectionMask};
use parquet::basic::{Compression, Encoding, ZstdLevel};
use parquet::file::properties::WriterProperties;
use parquet::format::KeyValue;
use parquet::schema::types::ColumnPath;
use quant::ms1_index::{Ms1RunIndex, Ms1Scan};

/// Major version of the table layout. A reader rejects any other major.
pub const MS1_SCHEMA_MAJOR: u32 = 1;

const SCHEMA_KEY: &str = "andes.ms1.schema";
const SCHEMA_VERSION: &str = "1.0";
/// Schema version of a table with the `mobility` column.
const SCHEMA_VERSION_MOBILITY: &str = "1.1";
/// Metadata key of the run stem.
pub const RUN_KEY: &str = "andes.ms1.run";
/// Scans per row group / record batch.
const SCANS_PER_GROUP: usize = 256;

const COLUMNS: [&str; 6] = [
    "native_id",
    "rt_s",
    "faims_cv",
    "mz",
    "intensity",
    "mobility",
];

fn item(data_type: DataType) -> FieldRef {
    Arc::new(Field::new("item", data_type, false))
}

fn schema(with_mobility: bool) -> Schema {
    let mut fields = vec![
        Field::new(COLUMNS[0], DataType::Utf8, false),
        Field::new(COLUMNS[1], DataType::Float64, false),
        Field::new(COLUMNS[2], DataType::Float32, true),
        Field::new(COLUMNS[3], DataType::List(item(DataType::Float64)), false),
        Field::new(COLUMNS[4], DataType::List(item(DataType::Float32)), false),
    ];
    if with_mobility {
        fields.push(Field::new(
            COLUMNS[5],
            DataType::List(item(DataType::Float32)),
            false,
        ));
    }
    Schema::new(fields)
}

/// ZSTD, no dictionaries (scan ids and centroid values rarely repeat), and
/// byte-stream-split floats, which compress better than plain ones.
fn properties(version: &str, run: &str, with_mobility: bool) -> WriterProperties {
    let list_values =
        |name: &str| ColumnPath::new(vec![name.to_string(), "list".into(), "item".into()]);
    let mut b = WriterProperties::builder()
        .set_compression(Compression::ZSTD(ZstdLevel::default()))
        .set_dictionary_enabled(false)
        .set_column_encoding(ColumnPath::from(COLUMNS[1]), Encoding::BYTE_STREAM_SPLIT)
        .set_column_encoding(list_values(COLUMNS[3]), Encoding::BYTE_STREAM_SPLIT)
        .set_column_encoding(list_values(COLUMNS[4]), Encoding::BYTE_STREAM_SPLIT);
    if with_mobility {
        b = b.set_column_encoding(list_values(COLUMNS[5]), Encoding::BYTE_STREAM_SPLIT);
    }
    b.set_max_row_group_size(SCANS_PER_GROUP)
        .set_key_value_metadata(Some(vec![
            KeyValue::new(SCHEMA_KEY.to_string(), version.to_string()),
            KeyValue::new(RUN_KEY.to_string(), run.to_string()),
        ]))
        .build()
}

fn offsets(lens: impl Iterator<Item = usize>) -> Result<OffsetBuffer<i32>, String> {
    let mut v = vec![0i32];
    let mut end = 0usize;
    for n in lens {
        end += n;
        v.push(i32::try_from(end).map_err(|_| "too many centroids in one row group".to_string())?);
    }
    Ok(OffsetBuffer::new(ScalarBuffer::from(v)))
}

fn batch(
    schema: &Arc<Schema>,
    index: &Ms1RunIndex,
    scans: std::ops::Range<usize>,
) -> Result<RecordBatch, String> {
    let ids = StringArray::from_iter_values(scans.clone().map(|i| index.scan(i).id.as_str()));
    let rts = Float64Array::from_iter_values(scans.clone().map(|i| index.rt(i)));
    let cvs: Float32Array = scans.clone().map(|i| index.scan(i).faims_cv).collect();
    let peaks = || scans.clone().flat_map(|i| index.peaks(i).iter());
    let offs = offsets(scans.clone().map(|i| index.peaks(i).len()))?;
    let mz = Float64Array::from_iter_values(peaks().map(|p| p.0));
    let intensity = Float32Array::from_iter_values(peaks().map(|p| p.1));
    let mz = ListArray::try_new(item(DataType::Float64), offs.clone(), Arc::new(mz), None)
        .map_err(|e| e.to_string())?;
    let intensity = ListArray::try_new(
        item(DataType::Float32),
        offs.clone(),
        Arc::new(intensity),
        None,
    )
    .map_err(|e| e.to_string())?;
    let mut cols: Vec<ArrayRef> = vec![
        Arc::new(ids),
        Arc::new(rts),
        Arc::new(cvs),
        Arc::new(mz),
        Arc::new(intensity),
    ];
    if schema.fields().len() == COLUMNS.len() {
        // A scan without mobility (none of its own) stores an empty list.
        let k0_offs = offsets(scans.clone().map(|i| index.mobility(i).len()))?;
        let k0 = Float32Array::from_iter_values(
            scans
                .clone()
                .flat_map(|i| index.mobility(i).iter().copied()),
        );
        cols.push(Arc::new(
            ListArray::try_new(item(DataType::Float32), k0_offs, Arc::new(k0), None)
                .map_err(|e| e.to_string())?,
        ));
    }
    RecordBatch::try_new(schema.clone(), cols).map_err(|e| e.to_string())
}

fn write_with_version<W: Write + Send>(
    out: W,
    run: &str,
    index: &Ms1RunIndex,
    version: &str,
) -> Result<(), String> {
    let with_mobility = index.has_mobility();
    let version = if with_mobility && version == SCHEMA_VERSION {
        SCHEMA_VERSION_MOBILITY
    } else {
        version
    };
    let schema = Arc::new(schema(with_mobility));
    let mut writer = ArrowWriter::try_new(
        out,
        schema.clone(),
        Some(properties(version, run, with_mobility)),
    )
    .map_err(|e| e.to_string())?;
    let mut start = 0;
    while start < index.len() {
        let end = (start + SCANS_PER_GROUP).min(index.len());
        writer
            .write(&batch(&schema, index, start..end)?)
            .map_err(|e| e.to_string())?;
        writer.flush().map_err(|e| e.to_string())?;
        start = end;
    }
    writer.close().map_err(|e| e.to_string())?;
    Ok(())
}

/// Write the MS1 table of run `run` to `out`.
pub fn write_ms1<W: Write + Send>(out: W, run: &str, index: &Ms1RunIndex) -> Result<(), String> {
    write_with_version(out, run, index, SCHEMA_VERSION)
}

/// `<path>.partial`, where a table is written before it is renamed to `path`.
pub fn partial_path(path: &Path) -> PathBuf {
    let mut partial = path.as_os_str().to_owned();
    partial.push(".partial");
    PathBuf::from(partial)
}

fn write_file_with_version(
    path: &Path,
    run: &str,
    index: &Ms1RunIndex,
    version: &str,
) -> Result<(), String> {
    let partial = partial_path(path);
    let written = File::create(&partial)
        .map_err(|e| e.to_string())
        .and_then(|f| write_with_version(f, run, index, version))
        .and_then(|()| std::fs::rename(&partial, path).map_err(|e| e.to_string()))
        .map_err(|e| format!("{}: {e}", path.display()));
    if written.is_err() {
        let _ = std::fs::remove_file(&partial);
    }
    written
}

/// Write the MS1 table of run `run` to `path`. Every error names `path`.
pub fn write_ms1_table(path: &Path, run: &str, index: &Ms1RunIndex) -> Result<(), String> {
    write_file_with_version(path, run, index, SCHEMA_VERSION)
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

/// The values of a list column and its offsets, checked to be null-free.
fn list<T: 'static + Array>(b: &RecordBatch, i: usize) -> Result<(&T, &[i32]), String> {
    let list = column::<ListArray>(b, i)?;
    let values = list.values().as_any().downcast_ref::<T>().ok_or_else(|| {
        format!(
            "column {} has unexpected item type {}",
            b.schema().field(i).name(),
            list.value_type()
        )
    })?;
    if list.null_count() != 0 || values.null_count() != 0 {
        return Err(format!("column {} holds nulls", b.schema().field(i).name()));
    }
    Ok((values, list.value_offsets()))
}

fn decode(b: &RecordBatch, scans: &mut Vec<Ms1Scan>) -> Result<(), String> {
    if b.num_columns() != COLUMNS.len() && b.num_columns() + 1 != COLUMNS.len() {
        return Err(format!(
            "expected {} or {} columns, found {}",
            COLUMNS.len() - 1,
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
    let rts = column::<Float64Array>(b, 1)?;
    let cvs = column::<Float32Array>(b, 2)?;
    if ids.null_count() != 0 || rts.null_count() != 0 {
        return Err("native_id or rt_s holds nulls".to_string());
    }
    let (mz, mz_off) = list::<Float64Array>(b, 3)?;
    let (intensity, int_off) = list::<Float32Array>(b, 4)?;
    let k0 = (b.num_columns() == COLUMNS.len())
        .then(|| list::<Float32Array>(b, 5))
        .transpose()?;
    let (mz, intensity) = (mz.values(), intensity.values());
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
        let mobility = match k0 {
            Some((values, off)) => {
                let (klo, khi) = (off[i] as usize, off[i + 1] as usize);
                if khi - klo != 0 && khi - klo != hi - lo {
                    return Err(format!(
                        "scan {:?}: {} m/z values but {} mobilities",
                        ids.value(i),
                        hi - lo,
                        khi - klo
                    ));
                }
                values.values()[klo..khi].to_vec()
            }
            None => Vec::new(),
        };
        scans.push(Ms1Scan {
            rt: rts.value(i),
            peaks: mz[lo..hi]
                .iter()
                .copied()
                .zip(intensity[ilo..ihi].iter().copied())
                .collect(),
            id: ids.value(i).to_string(),
            faims_cv: cvs.is_valid(i).then(|| cvs.value(i)),
            mobility,
        });
    }
    Ok(())
}

fn read(file: File) -> Result<Ms1RunIndex, String> {
    let builder = ParquetRecordBatchReaderBuilder::try_new(file).map_err(|e| e.to_string())?;
    let meta = builder.metadata().file_metadata();
    let version = meta
        .key_value_metadata()
        .and_then(|kv| kv.iter().find(|kv| kv.key == SCHEMA_KEY))
        .and_then(|kv| kv.value.clone())
        .ok_or_else(|| format!("missing metadata key {SCHEMA_KEY}; not an andes MS1 table"))?;
    let major: u32 = version
        .split('.')
        .next()
        .and_then(|m| m.parse().ok())
        .ok_or_else(|| format!("unparseable schema version {version:?}"))?;
    if major != MS1_SCHEMA_MAJOR {
        return Err(format!(
            "unsupported MS1 table schema {version} (this build reads major {MS1_SCHEMA_MAJOR})"
        ));
    }
    let n_rows = usize::try_from(meta.num_rows()).unwrap_or(0);
    let reader = builder
        .with_batch_size(SCANS_PER_GROUP)
        .build()
        .map_err(|e| e.to_string())?;
    let mut scans = Vec::with_capacity(n_rows.min(1 << 20));
    for b in reader {
        decode(&b.map_err(|e| e.to_string())?, &mut scans)?;
    }
    Ok(Ms1RunIndex::new(scans))
}

/// Read the MS1 table in `file`, opened from `path`. Every error names `path`.
pub fn read_ms1(file: File, path: &Path) -> Result<Ms1RunIndex, String> {
    read(file).map_err(|e| format!("{}: {e}", path.display()))
}

/// Median spacing in seconds between consecutive MS1 scans of the same FAIMS
/// compensation voltage (all scans when there is none), read from the `rt_s`
/// and `faims_cv` columns only. A precursor is sampled at its own voltage, so
/// this is the time between two of its MS1 points. `None` with fewer than two
/// scans of any voltage. Every error names `path`.
pub fn read_ms1_scan_spacing(path: &Path) -> Result<Option<f64>, String> {
    let named = |e: String| format!("{}: {e}", path.display());
    let file = File::open(path).map_err(|e| named(e.to_string()))?;
    let builder =
        ParquetRecordBatchReaderBuilder::try_new(file).map_err(|e| named(e.to_string()))?;
    let mask = ProjectionMask::roots(builder.parquet_schema(), [1, 2]);
    let reader = builder
        .with_projection(mask)
        .build()
        .map_err(|e| named(e.to_string()))?;
    // Scans are stored in retention-time order; the last RT per voltage.
    let mut last: std::collections::HashMap<Option<u32>, f64> = std::collections::HashMap::new();
    let mut gaps: Vec<f64> = Vec::new();
    for b in reader {
        let b = b.map_err(|e| named(e.to_string()))?;
        let rt = b
            .column_by_name(COLUMNS[1])
            .and_then(|c| c.as_any().downcast_ref::<Float64Array>())
            .ok_or_else(|| named(format!("column {} is not float64", COLUMNS[1])))?;
        let cv = b
            .column_by_name(COLUMNS[2])
            .and_then(|c| c.as_any().downcast_ref::<Float32Array>())
            .ok_or_else(|| named(format!("column {} is not float32", COLUMNS[2])))?;
        for i in 0..b.num_rows() {
            let key = cv.is_valid(i).then(|| cv.value(i).to_bits());
            let t = rt.value(i);
            if let Some(prev) = last.insert(key, t) {
                if t > prev {
                    gaps.push(t - prev);
                }
            }
        }
    }
    if gaps.is_empty() {
        return Ok(None);
    }
    gaps.sort_by(f64::total_cmp);
    Ok(Some(gaps[gaps.len() / 2]))
}

/// Read the MS1 table at `path`. Every error names `path`.
pub fn read_ms1_table(path: &Path) -> Result<Ms1RunIndex, String> {
    let file = File::open(path).map_err(|e| format!("{}: {e}", path.display()))?;
    read_ms1(file, path)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn scan(rt: f64, peaks: Vec<(f64, f32)>, id: &str, faims_cv: Option<f32>) -> Ms1Scan {
        Ms1Scan {
            rt,
            peaks,
            id: id.into(),
            faims_cv,
            mobility: Vec::new(),
        }
    }

    /// Bit-exact equality of two indexes: order, RT, ids, FAIMS voltage,
    /// centroids and detection limits.
    fn assert_same(a: &Ms1RunIndex, b: &Ms1RunIndex) {
        assert_eq!(a.len(), b.len());
        for i in 0..a.len() {
            assert_eq!(a.rt(i).to_bits(), b.rt(i).to_bits(), "rt of scan {i}");
            assert_eq!(a.scan(i).id, b.scan(i).id);
            assert_eq!(
                a.scan(i).faims_cv.map(f32::to_bits),
                b.scan(i).faims_cv.map(f32::to_bits)
            );
            assert_eq!(a.peaks(i).len(), b.peaks(i).len());
            for (p, q) in a.peaks(i).iter().zip(b.peaks(i)) {
                assert_eq!(p.0.to_bits(), q.0.to_bits());
                assert_eq!(p.1.to_bits(), q.1.to_bits());
            }
            assert_eq!(
                a.detection_limit(i).to_bits(),
                b.detection_limit(i).to_bits()
            );
            let bits = |v: &[f32]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
            assert_eq!(bits(a.mobility(i)), bits(b.mobility(i)));
        }
    }

    fn index() -> Ms1RunIndex {
        let nan = f32::from_bits(0x7fc0_1234);
        Ms1RunIndex::new(vec![
            scan(
                3.123456789,
                vec![(600.123456789, 1234.5678), (700.0, -0.0), (800.0, nan)],
                "scan=2",
                Some(-45.0),
            ),
            scan(1.0, vec![], "", None),
            // Equal RT: the stored order is kept.
            scan(1.0, vec![(-0.0, 5.0)], "scan=9", Some(-0.0)),
            scan(2.5, vec![(500.0, 1.0)], "contr\u{e9}l\u{e9}=3", Some(nan)),
        ])
    }

    #[test]
    fn scan_spacing_is_per_faims_voltage() {
        let dir = tempfile::tempdir().unwrap();
        // Two voltages alternating every 1 s: each voltage is sampled every 2 s.
        let faims: Vec<Ms1Scan> = (0..10)
            .map(|i| {
                let cv = if i % 2 == 0 { -45.0 } else { -60.0 };
                scan(i as f64, vec![(500.0, 1.0)], "", Some(cv))
            })
            .collect();
        let p = dir.path().join("faims.ms1.parquet");
        write_ms1_table(&p, "faims", &Ms1RunIndex::new(faims)).unwrap();
        assert_eq!(read_ms1_scan_spacing(&p).unwrap(), Some(2.0));
        // No FAIMS: every MS1 scan, 1.5 s apart.
        let plain: Vec<Ms1Scan> = (0..5)
            .map(|i| scan(1.5 * i as f64, vec![(500.0, 1.0)], "", None))
            .collect();
        let p = dir.path().join("plain.ms1.parquet");
        write_ms1_table(&p, "plain", &Ms1RunIndex::new(plain)).unwrap();
        assert_eq!(read_ms1_scan_spacing(&p).unwrap(), Some(1.5));
        // One scan: no spacing.
        let p = dir.path().join("one.ms1.parquet");
        let one = Ms1RunIndex::new(vec![scan(3.0, vec![], "", None)]);
        write_ms1_table(&p, "one", &one).unwrap();
        assert_eq!(read_ms1_scan_spacing(&p).unwrap(), None);
    }

    #[test]
    fn round_trip_is_bit_exact_and_keeps_scan_order() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("run.ms1.parquet");
        let idx = index();
        write_ms1_table(&p, "run", &idx).unwrap();
        assert!(!partial_path(&p).exists());
        let back = read_ms1_table(&p).unwrap();
        assert_same(&back, &idx);
        assert_eq!(back.scan(0).id, "");
        assert_eq!(back.scan(1).id, "scan=9");
        // Reading again gives the same index.
        assert_same(&read_ms1_table(&p).unwrap(), &idx);
    }

    #[test]
    fn mobility_round_trips_and_is_written_only_when_present() {
        let dir = tempfile::tempdir().unwrap();
        let mut scans = vec![
            scan(
                1.0,
                vec![(500.0, 2.0), (500.0, 3.0), (600.0, 1.0)],
                "frame=1",
                None,
            ),
            // A scan without centroids in a mobility-resolved run.
            scan(2.0, vec![], "frame=3", None),
        ];
        scans[0].mobility = vec![0.91, 1.02, 0.85];
        let idx = Ms1RunIndex::new(scans);
        assert!(idx.has_mobility());
        let p = dir.path().join("tims.ms1.parquet");
        write_ms1_table(&p, "tims", &idx).unwrap();
        let back = read_ms1_table(&p).unwrap();
        assert_same(&back, &idx);
        assert_eq!(back.mobility(0), &[0.91f32, 1.02, 0.85]);
        assert!(back.mobility(1).is_empty());
        let columns = |p: &Path| {
            let file = File::open(p).unwrap();
            ParquetRecordBatchReaderBuilder::try_new(file)
                .unwrap()
                .schema()
                .fields()
                .len()
        };
        assert_eq!(columns(&p), 6);
        // A run without mobility keeps the five-column layout.
        let q = dir.path().join("plain.ms1.parquet");
        write_ms1_table(&q, "plain", &index()).unwrap();
        assert_eq!(columns(&q), 5);
    }

    #[test]
    fn many_row_groups_and_empty_runs_round_trip() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("big.ms1.parquet");
        let scans: Vec<Ms1Scan> = (0..(2 * SCANS_PER_GROUP + 7))
            .map(|i| {
                let peaks = (0..(i % 5))
                    .map(|k| (400.0 + k as f64 * 1.000_1 + i as f64, (i * k) as f32))
                    .collect();
                scan(i as f64 * 0.5, peaks, &format!("scan={i}"), None)
            })
            .collect();
        let idx = Ms1RunIndex::new(scans);
        write_ms1_table(&p, "big", &idx).unwrap();
        let file = File::open(&p).unwrap();
        let groups = ParquetRecordBatchReaderBuilder::try_new(file)
            .unwrap()
            .metadata()
            .num_row_groups();
        assert_eq!(groups, 3);
        assert_same(&read_ms1_table(&p).unwrap(), &idx);

        let e = dir.path().join("empty.ms1.parquet");
        write_ms1_table(&e, "empty", &Ms1RunIndex::default()).unwrap();
        assert!(read_ms1_table(&e).unwrap().is_empty());
    }

    #[test]
    fn metadata_carries_schema_and_run() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("r1.ms1.parquet");
        write_ms1_table(&p, "r1", &index()).unwrap();
        let builder = ParquetRecordBatchReaderBuilder::try_new(File::open(&p).unwrap()).unwrap();
        let kv: Vec<(String, Option<String>)> = builder
            .metadata()
            .file_metadata()
            .key_value_metadata()
            .unwrap()
            .iter()
            .map(|kv| (kv.key.clone(), kv.value.clone()))
            .collect();
        assert!(
            kv.contains(&(SCHEMA_KEY.into(), Some("1.0".into()))),
            "{kv:?}"
        );
        assert!(kv.contains(&(RUN_KEY.into(), Some("r1".into()))), "{kv:?}");
    }

    #[test]
    fn rejects_other_majors_and_foreign_files_naming_the_file() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("v2.ms1.parquet");
        write_file_with_version(&p, "v2", &index(), "2.0").unwrap();
        let e = read_ms1_table(&p).unwrap_err();
        assert!(
            e.contains("v2.ms1.parquet") && e.contains("schema 2.0"),
            "{e}"
        );

        // A Parquet file without the schema key (another andes table).
        let other = dir.path().join("other.parquet");
        let schema = Arc::new(Schema::new(vec![Field::new("x", DataType::Float64, false)]));
        let mut w =
            ArrowWriter::try_new(File::create(&other).unwrap(), schema.clone(), None).unwrap();
        let col: ArrayRef = Arc::new(Float64Array::from(vec![1.0]));
        w.write(&RecordBatch::try_new(schema, vec![col]).unwrap())
            .unwrap();
        w.close().unwrap();
        let e = read_ms1_table(&other).unwrap_err();
        assert!(e.contains("other.parquet") && e.contains(SCHEMA_KEY), "{e}");

        let junk = dir.path().join("junk.ms1.parquet");
        std::fs::write(&junk, b"ANDESRC1").unwrap();
        let e = read_ms1_table(&junk).unwrap_err();
        assert!(e.contains("junk.ms1.parquet"), "{e}");

        let e = read_ms1_table(&dir.path().join("nope.ms1.parquet")).unwrap_err();
        assert!(e.contains("nope.ms1.parquet"), "{e}");
    }

    #[test]
    fn truncated_table_names_the_file() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("run.ms1.parquet");
        let idx = Ms1RunIndex::new(vec![scan(1.0, vec![(500.0, 1.0); 100], "", None)]);
        write_ms1_table(&p, "run", &idx).unwrap();
        let bytes = std::fs::read(&p).unwrap();
        for cut in [6, bytes.len() / 2] {
            std::fs::write(&p, &bytes[..bytes.len() - cut]).unwrap();
            let e = read_ms1_table(&p).unwrap_err();
            assert!(e.contains("run.ms1.parquet"), "{e}");
        }
    }

    #[test]
    fn renamed_column_is_rejected_naming_the_file() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("c.ms1.parquet");
        let idx = index();
        let good = batch(&Arc::new(schema(false)), &idx, 0..idx.len()).unwrap();
        let mut fields: Vec<Field> = schema(false)
            .fields()
            .iter()
            .map(|f| f.as_ref().clone())
            .collect();
        fields[1] = Field::new("rt_min", DataType::Float64, false);
        let renamed = Arc::new(Schema::new(fields));
        let b = RecordBatch::try_new(renamed.clone(), good.columns().to_vec()).unwrap();
        let mut w = ArrowWriter::try_new(
            File::create(&p).unwrap(),
            renamed,
            Some(properties("1.0", "c", false)),
        )
        .unwrap();
        w.write(&b).unwrap();
        w.close().unwrap();
        let e = read_ms1_table(&p).unwrap_err();
        assert!(e.contains("c.ms1.parquet") && e.contains("rt_min"), "{e}");
    }

    #[test]
    fn failed_write_leaves_no_file() {
        let dir = tempfile::tempdir().unwrap();
        let bad = dir.path().join("missing-dir").join("d.ms1.parquet");
        let e = write_ms1_table(&bad, "d", &index()).unwrap_err();
        assert!(e.contains("d.ms1.parquet"), "{e}");
        assert!(!bad.exists() && !partial_path(&bad).exists());
    }
}
