//! Per-run PSM table `<stem>.quant.parquet`: every rank-1 PSM of one run
//! (targets and decoys) with the fields the multi-run quantification stage
//! needs. The file carries a schema version in its key-value metadata; a
//! reader rejects a different major version. The metadata also maps each
//! modification accession used in the peptidoforms to its mass delta
//! (`andes.run_quant.mods`, `<accession>=<delta>` joined by `;`), so a
//! peptidoform's fragment masses can be computed from the table alone.

use std::collections::HashMap;
use std::fs::File;
use std::path::Path;
use std::sync::Arc;

use arrow::array::{
    Array, ArrayRef, BooleanArray, BooleanBuilder, Float32Array, Float32Builder, Float64Array,
    Float64Builder, ListArray, ListBuilder, StringArray, StringBuilder, UInt8Array, UInt8Builder,
};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::ArrowWriter;
use parquet::basic::{Compression, ZstdLevel};
use parquet::file::properties::WriterProperties;
use parquet::format::KeyValue;

/// Major version of the table layout. A reader rejects any other major.
pub const RUN_QUANT_SCHEMA_MAJOR: u32 = 1;

const SCHEMA_KEY: &str = "andes.run_quant.schema";
const RUN_KEY: &str = "andes.run_quant.run";
const MODS_KEY: &str = "andes.run_quant.mods";
const SCHEMA_VERSION: &str = "1.1";
/// Schema version of a table with the `ion_mobility` column, which is
/// written only when some PSM records a precursor 1/K0.
const SCHEMA_VERSION_MOBILITY: &str = "1.2";
/// Rows per Parquet row group / record batch.
const BATCH_ROWS: usize = 65_536;

/// One rank-1 PSM of a run, target or decoy.
#[derive(Debug, Clone, PartialEq)]
pub struct RunPsmRow {
    /// PIN SpecId, joins the Percolator PSMId.
    pub spec_id: String,
    /// MS2 retention time, seconds.
    pub rt_s: f64,
    /// `"<peptidoform>/<charge>"`.
    pub precursor_key: String,
    pub peptidoform: String,
    pub charge: u8,
    /// Theoretical monoisotopic m/z.
    pub mz: f64,
    pub is_decoy: bool,
    /// andes RawScore of the rank-1 match.
    pub search_score: f32,
    pub faims_cv: Option<f32>,
    /// Theoretical isotope fractions (sum 1).
    pub envelope: Vec<f32>,
    /// `;`-joined accessions.
    pub proteins: String,
    /// Precursor 1/K0 of the spectrum, when the source records it (Bruker
    /// timsTOF `.d`).
    pub ion_mobility: Option<f32>,
}

fn schema(with_mobility: bool) -> Schema {
    let mut fields = vec![
        Field::new("spec_id", DataType::Utf8, false),
        Field::new("rt_s", DataType::Float64, false),
        Field::new("precursor_key", DataType::Utf8, false),
        Field::new("peptidoform", DataType::Utf8, false),
        Field::new("charge", DataType::UInt8, false),
        Field::new("mz", DataType::Float64, false),
        Field::new("is_decoy", DataType::Boolean, false),
        Field::new("search_score", DataType::Float32, false),
        Field::new("faims_cv", DataType::Float32, true),
        Field::new(
            "envelope",
            DataType::List(Arc::new(Field::new("item", DataType::Float32, true))),
            false,
        ),
        Field::new("proteins", DataType::Utf8, false),
    ];
    if with_mobility {
        fields.push(Field::new("ion_mobility", DataType::Float32, true));
    }
    Schema::new(fields)
}

fn batch(
    schema: &Arc<Schema>,
    rows: &[RunPsmRow],
    with_mobility: bool,
) -> Result<RecordBatch, String> {
    let mut spec_id = StringBuilder::new();
    let mut rt_s = Float64Builder::new();
    let mut key = StringBuilder::new();
    let mut pepf = StringBuilder::new();
    let mut charge = UInt8Builder::new();
    let mut mz = Float64Builder::new();
    let mut decoy = BooleanBuilder::new();
    let mut score = Float32Builder::new();
    let mut cv = Float32Builder::new();
    let mut env = ListBuilder::new(Float32Builder::new());
    let mut prot = StringBuilder::new();
    let mut k0 = Float32Builder::new();
    for r in rows {
        spec_id.append_value(&r.spec_id);
        rt_s.append_value(r.rt_s);
        key.append_value(&r.precursor_key);
        pepf.append_value(&r.peptidoform);
        charge.append_value(r.charge);
        mz.append_value(r.mz);
        decoy.append_value(r.is_decoy);
        score.append_value(r.search_score);
        cv.append_option(r.faims_cv);
        env.values().append_slice(&r.envelope);
        env.append(true);
        prot.append_value(&r.proteins);
        k0.append_option(r.ion_mobility);
    }
    let mut cols: Vec<ArrayRef> = vec![
        Arc::new(spec_id.finish()),
        Arc::new(rt_s.finish()),
        Arc::new(key.finish()),
        Arc::new(pepf.finish()),
        Arc::new(charge.finish()),
        Arc::new(mz.finish()),
        Arc::new(decoy.finish()),
        Arc::new(score.finish()),
        Arc::new(cv.finish()),
        Arc::new(env.finish()),
        Arc::new(prot.finish()),
    ];
    if with_mobility {
        cols.push(Arc::new(k0.finish()));
    }
    RecordBatch::try_new(schema.clone(), cols).map_err(|e| e.to_string())
}

fn encode_mods(mods: &[(String, f64)]) -> Result<String, String> {
    let mut parts = Vec::with_capacity(mods.len());
    for (acc, delta) in mods {
        if acc.is_empty() || acc.contains([';', '=']) || !delta.is_finite() {
            return Err(format!("modification {acc:?} = {delta} cannot be stored"));
        }
        parts.push(format!("{acc}={delta:?}"));
    }
    Ok(parts.join(";"))
}

fn write_versioned(
    path: &Path,
    run_stem: &str,
    rows: &[RunPsmRow],
    mods: &[(String, f64)],
    version: &str,
) -> Result<(), String> {
    let err = |e: &dyn std::fmt::Display| format!("{}: {e}", path.display());
    let with_mobility = rows.iter().any(|r| r.ion_mobility.is_some());
    let version = if with_mobility && version == SCHEMA_VERSION {
        SCHEMA_VERSION_MOBILITY
    } else {
        version
    };
    let schema = Arc::new(schema(with_mobility));
    let props = WriterProperties::builder()
        .set_compression(Compression::ZSTD(ZstdLevel::default()))
        .set_max_row_group_size(BATCH_ROWS)
        .set_key_value_metadata(Some(vec![
            KeyValue::new(SCHEMA_KEY.to_string(), version.to_string()),
            KeyValue::new(RUN_KEY.to_string(), run_stem.to_string()),
            KeyValue::new(
                MODS_KEY.to_string(),
                encode_mods(mods).map_err(|e| err(&e))?,
            ),
        ]))
        .build();
    crate::atomic::write_atomic(path, |partial| {
        let io = |e: &dyn std::fmt::Display| std::io::Error::other(e.to_string());
        let file = File::create(partial)?;
        let mut writer =
            ArrowWriter::try_new(file, schema.clone(), Some(props)).map_err(|e| io(&e))?;
        for chunk in rows.chunks(BATCH_ROWS) {
            let b = batch(&schema, chunk, with_mobility).map_err(|e| io(&e))?;
            writer.write(&b).map_err(|e| io(&e))?;
        }
        writer.close().map_err(|e| io(&e))?;
        Ok(())
    })
    .map_err(|e| err(&e))
}

/// Write the per-run PSM table for run `run_stem`, with the mass delta of
/// every modification accession its peptidoforms use.
pub fn write_run_quant_table(
    path: &Path,
    run_stem: &str,
    rows: &[RunPsmRow],
    mods: &[(String, f64)],
) -> Result<(), String> {
    write_versioned(path, run_stem, rows, mods, SCHEMA_VERSION)
}

#[cfg(test)]
pub(crate) fn write_with_schema_version(
    path: &Path,
    run_stem: &str,
    rows: &[RunPsmRow],
    version: &str,
) -> Result<(), String> {
    write_versioned(path, run_stem, rows, &[], version)
}

/// The modification accessions of a per-run PSM table and their mass deltas
/// (empty for a table written before they were recorded). Every error names
/// `path`.
pub fn read_run_quant_mods(path: &Path) -> Result<HashMap<String, f64>, String> {
    let err = |e: &dyn std::fmt::Display| format!("{}: {e}", path.display());
    let file = File::open(path).map_err(|e| err(&e))?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file).map_err(|e| err(&e))?;
    let text = builder
        .metadata()
        .file_metadata()
        .key_value_metadata()
        .and_then(|kv| kv.iter().find(|kv| kv.key == MODS_KEY))
        .and_then(|kv| kv.value.clone())
        .unwrap_or_default();
    let mut out = HashMap::new();
    for part in text.split(';').filter(|p| !p.is_empty()) {
        let (acc, delta) = part
            .split_once('=')
            .ok_or_else(|| err(&format!("unreadable {MODS_KEY} entry {part:?}")))?;
        let delta: f64 = delta
            .parse()
            .map_err(|_| err(&format!("unreadable {MODS_KEY} entry {part:?}")))?;
        out.insert(acc.to_string(), delta);
    }
    Ok(out)
}

fn column<'a, T: 'static>(b: &'a RecordBatch, i: usize, path: &Path) -> Result<&'a T, String> {
    b.column(i).as_any().downcast_ref::<T>().ok_or_else(|| {
        format!(
            "{}: column {} has unexpected type {}",
            path.display(),
            b.schema().field(i).name(),
            b.column(i).data_type()
        )
    })
}

/// The columns of a per-run PSM table, in file order. The last,
/// `ion_mobility`, is present only in tables where some PSM records a
/// precursor 1/K0; a table without it reads `None`.
pub const RUN_QUANT_COLUMNS: [&str; 12] = [
    "spec_id",
    "rt_s",
    "precursor_key",
    "peptidoform",
    "charge",
    "mz",
    "is_decoy",
    "search_score",
    "faims_cv",
    "envelope",
    "proteins",
    "ion_mobility",
];

/// A per-run PSM table opened for reading one record batch at a time, with
/// only the requested columns decoded. A column left out reads as its empty
/// value (empty string or list, 0, `false`, `None`).
pub struct RunQuantBatches {
    path: std::path::PathBuf,
    run: String,
    wanted: [bool; 12],
    reader: parquet::arrow::arrow_reader::ParquetRecordBatchReader,
}

/// Opens a per-run PSM table for [`RunQuantBatches`], decoding only
/// `columns` (names of [`RUN_QUANT_COLUMNS`]; all of them when empty). Every
/// error names `path`.
pub fn open_run_quant_table(path: &Path, columns: &[&str]) -> Result<RunQuantBatches, String> {
    let err = |e: &dyn std::fmt::Display| format!("{}: {e}", path.display());
    if let Some(c) = columns.iter().find(|c| !RUN_QUANT_COLUMNS.contains(c)) {
        return Err(err(&format!("unknown run-quant column {c:?}")));
    }
    let file = File::open(path).map_err(|e| err(&e))?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file).map_err(|e| err(&e))?;
    let meta = builder
        .metadata()
        .file_metadata()
        .key_value_metadata()
        .cloned()
        .unwrap_or_default();
    let get = |key: &str| {
        meta.iter()
            .find(|kv| kv.key == key)
            .and_then(|kv| kv.value.clone())
    };
    let version =
        get(SCHEMA_KEY).ok_or_else(|| err(&format!("missing metadata key {SCHEMA_KEY}")))?;
    let major: u32 = version
        .split('.')
        .next()
        .and_then(|m| m.parse().ok())
        .ok_or_else(|| err(&format!("unparseable schema version {version:?}")))?;
    if major != RUN_QUANT_SCHEMA_MAJOR {
        return Err(err(&format!(
            "unsupported run-quant table schema {version} (this build reads major {RUN_QUANT_SCHEMA_MAJOR})"
        )));
    }
    let run = get(RUN_KEY).ok_or_else(|| err(&format!("missing metadata key {RUN_KEY}")))?;
    let fields = builder.schema().fields().clone();
    if fields.len() + 1 != RUN_QUANT_COLUMNS.len() && fields.len() != RUN_QUANT_COLUMNS.len() {
        return Err(err(&format!(
            "expected {} or {} columns, found {}",
            RUN_QUANT_COLUMNS.len() - 1,
            RUN_QUANT_COLUMNS.len(),
            fields.len()
        )));
    }
    let mut wanted = [false; 12];
    for (i, name) in RUN_QUANT_COLUMNS.iter().enumerate().take(fields.len()) {
        if fields[i].name() != name {
            return Err(err(&format!(
                "column {:?} where {name:?} was expected",
                fields[i].name()
            )));
        }
        wanted[i] = columns.is_empty() || columns.contains(name);
    }
    let mask = parquet::arrow::ProjectionMask::roots(
        builder.parquet_schema(),
        (0..wanted.len()).filter(|&i| wanted[i]),
    );
    let reader = builder
        .with_projection(mask)
        .with_batch_size(BATCH_ROWS)
        .build()
        .map_err(|e| err(&e))?;
    Ok(RunQuantBatches {
        path: path.to_path_buf(),
        run,
        wanted,
        reader,
    })
}

impl RunQuantBatches {
    /// The run stem recorded in the table.
    pub fn run(&self) -> &str {
        &self.run
    }

    fn rows(&self, b: &RecordBatch) -> Result<Vec<RunPsmRow>, String> {
        let path = self.path.as_path();
        let err = |e: &dyn std::fmt::Display| format!("{}: {e}", path.display());
        // Projected batches hold the wanted columns in file order.
        let mut at = [usize::MAX; 12];
        let mut k = 0;
        for (i, w) in self.wanted.iter().enumerate() {
            if *w {
                at[i] = k;
                k += 1;
            }
        }
        if b.num_columns() != k {
            return Err(err(&format!(
                "expected {k} columns, found {}",
                b.num_columns()
            )));
        }
        let col = |i: usize| (at[i] != usize::MAX).then_some(at[i]);
        let string = |i: usize| {
            col(i)
                .map(|c| column::<StringArray>(b, c, path))
                .transpose()
        };
        let (spec_id, key, pepf, prot) = (string(0)?, string(2)?, string(3)?, string(10)?);
        let rt_s = col(1)
            .map(|c| column::<Float64Array>(b, c, path))
            .transpose()?;
        let charge = col(4)
            .map(|c| column::<UInt8Array>(b, c, path))
            .transpose()?;
        let mz = col(5)
            .map(|c| column::<Float64Array>(b, c, path))
            .transpose()?;
        let decoy = col(6)
            .map(|c| column::<BooleanArray>(b, c, path))
            .transpose()?;
        let score = col(7)
            .map(|c| column::<Float32Array>(b, c, path))
            .transpose()?;
        let cv = col(8)
            .map(|c| column::<Float32Array>(b, c, path))
            .transpose()?;
        let env = col(9)
            .map(|c| column::<ListArray>(b, c, path))
            .transpose()?;
        let k0 = col(11)
            .map(|c| column::<Float32Array>(b, c, path))
            .transpose()?;
        let text =
            |a: Option<&StringArray>, i: usize| a.map_or(String::new(), |a| a.value(i).to_string());
        let mut rows = Vec::with_capacity(b.num_rows());
        for i in 0..b.num_rows() {
            let envelope = match env {
                Some(env) => {
                    let e = env.value(i);
                    e.as_any()
                        .downcast_ref::<Float32Array>()
                        .ok_or_else(|| err(&"envelope items are not Float32"))?
                        .values()
                        .to_vec()
                }
                None => Vec::new(),
            };
            rows.push(RunPsmRow {
                spec_id: text(spec_id, i),
                rt_s: rt_s.map_or(0.0, |a| a.value(i)),
                precursor_key: text(key, i),
                peptidoform: text(pepf, i),
                charge: charge.map_or(0, |a| a.value(i)),
                mz: mz.map_or(0.0, |a| a.value(i)),
                is_decoy: decoy.is_some_and(|a| a.value(i)),
                search_score: score.map_or(0.0, |a| a.value(i)),
                faims_cv: cv.and_then(|a| (!a.is_null(i)).then(|| a.value(i))),
                envelope,
                proteins: text(prot, i),
                ion_mobility: k0.and_then(|a| (!a.is_null(i)).then(|| a.value(i))),
            });
        }
        Ok(rows)
    }
}

impl Iterator for RunQuantBatches {
    type Item = Result<Vec<RunPsmRow>, String>;

    fn next(&mut self) -> Option<Self::Item> {
        let b = self.reader.next()?;
        Some(
            b.map_err(|e| format!("{}: {e}", self.path.display()))
                .and_then(|b| self.rows(&b)),
        )
    }
}

/// Read a per-run PSM table; returns `(run_stem, rows)`. Every error names
/// `path`.
pub fn read_run_quant_table(path: &Path) -> Result<(String, Vec<RunPsmRow>), String> {
    let mut batches = open_run_quant_table(path, &[])?;
    let run = batches.run().to_string();
    let mut rows = Vec::new();
    for b in batches.by_ref() {
        rows.extend(b?);
    }
    Ok((run, rows))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn row(i: usize) -> RunPsmRow {
        RunPsmRow {
            spec_id: format!("a/scan={i}_2_1"),
            rt_s: 1234.5 + i as f64,
            precursor_key: "PEPTIDEK/2".into(),
            peptidoform: "PEPTIDEK".into(),
            charge: 2,
            mz: 465.7380,
            is_decoy: i % 2 == 1,
            search_score: 41.0,
            faims_cv: if i % 2 == 0 { Some(-45.0) } else { None },
            envelope: vec![0.6, 0.3, 0.08, 0.02],
            proteins: "P1;P2".into(),
            ion_mobility: None,
        }
    }

    #[test]
    fn round_trip_preserves_rows_and_run() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("a.quant.parquet");
        let rows = vec![RunPsmRow {
            spec_id: "a/scan=7_2_1".into(),
            rt_s: 1234.5,
            precursor_key: "PEPTIDEK/2".into(),
            peptidoform: "PEPTIDEK".into(),
            charge: 2,
            mz: 465.7380,
            is_decoy: false,
            search_score: 41.0,
            faims_cv: Some(-45.0),
            envelope: vec![0.6, 0.3, 0.08, 0.02],
            proteins: "P1;P2".into(),
            ion_mobility: None,
        }];
        let mods = [
            ("UNIMOD:4".to_string(), 57.021464),
            ("UNIMOD:35".to_string(), 15.994915),
        ];
        write_run_quant_table(&p, "a", &rows, &mods).unwrap();
        let (stem, back) = read_run_quant_table(&p).unwrap();
        assert_eq!(stem, "a");
        assert_eq!(back, rows);
        let m = read_run_quant_mods(&p).unwrap();
        assert_eq!(m.len(), 2);
        assert_eq!(m["UNIMOD:4"].to_bits(), 57.021464f64.to_bits());
        assert_eq!(m["UNIMOD:35"].to_bits(), 15.994915f64.to_bits());
        let bad = [("A;B".to_string(), 1.0)];
        assert!(write_run_quant_table(&p, "a", &rows, &bad).is_err());
    }

    #[test]
    fn round_trip_mixed_nulls_decoys_and_empty() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("m.quant.parquet");
        let rows: Vec<RunPsmRow> = (0..5).map(row).collect();
        write_run_quant_table(&p, "m", &rows, &[]).unwrap();
        assert!(
            !crate::atomic::partial(&p).exists(),
            "written through <path>.partial and renamed"
        );
        assert!(read_run_quant_mods(&p).unwrap().is_empty());
        assert_eq!(read_run_quant_table(&p).unwrap(), ("m".to_string(), rows));

        let e = dir.path().join("e.quant.parquet");
        write_run_quant_table(&e, "e", &[], &[]).unwrap();
        assert_eq!(read_run_quant_table(&e).unwrap(), ("e".to_string(), vec![]));
    }

    #[test]
    fn projected_batches_read_only_the_requested_columns() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("m.quant.parquet");
        let rows: Vec<RunPsmRow> = (0..7).map(row).collect();
        write_run_quant_table(&p, "m", &rows, &[]).unwrap();
        let mut b = open_run_quant_table(&p, &["precursor_key", "is_decoy", "faims_cv"]).unwrap();
        assert_eq!(b.run(), "m");
        let got: Vec<RunPsmRow> = b.by_ref().flat_map(|b| b.unwrap()).collect();
        let want: Vec<RunPsmRow> = rows
            .iter()
            .map(|r| RunPsmRow {
                spec_id: String::new(),
                rt_s: 0.0,
                precursor_key: r.precursor_key.clone(),
                peptidoform: String::new(),
                charge: 0,
                mz: 0.0,
                is_decoy: r.is_decoy,
                search_score: 0.0,
                faims_cv: r.faims_cv,
                envelope: Vec::new(),
                proteins: String::new(),
                ion_mobility: None,
            })
            .collect();
        assert_eq!(got, want);
        let all: Vec<RunPsmRow> = open_run_quant_table(&p, &[])
            .unwrap()
            .flat_map(|b| b.unwrap())
            .collect();
        assert_eq!(all, rows);
        assert!(open_run_quant_table(&p, &["nope"]).is_err());
    }

    #[test]
    fn quant_table_rejects_unknown_major_version() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("b.quant.parquet");
        write_with_schema_version(&p, "b", &[], "2.0").unwrap();
        let e = read_run_quant_table(&p).unwrap_err();
        assert!(
            e.contains("b.quant.parquet") && e.contains("schema 2.0"),
            "{e}"
        );
    }

    #[test]
    fn missing_and_non_parquet_files_name_the_path() {
        let dir = tempfile::tempdir().unwrap();
        let missing = dir.path().join("nope.quant.parquet");
        let e = read_run_quant_table(&missing).unwrap_err();
        assert!(e.contains("nope.quant.parquet"), "{e}");

        let junk = dir.path().join("junk.quant.parquet");
        std::fs::write(&junk, b"not parquet").unwrap();
        let e = read_run_quant_table(&junk).unwrap_err();
        assert!(e.contains("junk.quant.parquet"), "{e}");
    }
}
