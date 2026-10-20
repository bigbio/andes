//! Per-run extraction table `<stem>.extract.parquet`: one row per fitted
//! (precursor, run, hypothesis) of the match-between-runs extract stage, with
//! the template-fit evidence. The file carries a schema version in its
//! key-value metadata; a reader rejects a different major version.
//!
//! Encoding: `kind` is 0 for a target, 1 for a mass-shift decoy (its nominal
//! shift in `mass_shift`) and 2 for an RT-shift decoy; `mass_shift` is 0 for
//! the other kinds. `role` is 0 identified, 1 transfer, 2 pseudo-transfer;
//! identified and transfer targets are each followed by their decoys (same
//! role), pseudo-transfers have none. Schema 2 added the identified rows'
//! decoys; a schema-1 table has none, so it is rejected. Schema 3 added the
//! MS2 evidence of the transfer and pseudo-transfer rows (`ms2_n_scans`,
//! `ms2_n_matched`, `ms2_angle`, `ms2_explained`; null for identified rows
//! and rows without a fit); a schema-2 table is rejected. Schema 4 added
//! the per-isotope heights `iso_h0`..`iso_h3` and the `tailing` of the fit;
//! a schema-3 table is rejected.
//! The evidence columns are null when the fit was not possible.
//! `n_donor_runs` counts the aligned runs with a PSM of the precursor, this
//! run included for an identified row and excluded otherwise: the runs that
//! informed the prediction. `library_q` is the internal library entry gate,
//! not an identification FDR. The file is written to `<path>.partial` and
//! renamed into place once complete.

use std::fs::File;
use std::path::Path;
use std::sync::Arc;

use arrow::array::{
    Array, ArrayRef, Float32Array, Float32Builder, Float64Array, Float64Builder, Int8Array,
    Int8Builder, UInt16Array, UInt16Builder, UInt32Array, UInt32Builder, UInt8Array, UInt8Builder,
};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::ArrowWriter;
use parquet::basic::{Compression, ZstdLevel};
use parquet::file::properties::WriterProperties;
use parquet::format::KeyValue;
use quant::decoys::Hypothesis;
use quant::evidence::Evidence;
use quant::ms2::{Ms2Evidence, ScanMatch};

/// Major version of the table layout. A reader rejects any other major.
pub const EXTRACT_SCHEMA_MAJOR: u32 = 4;

const SCHEMA_KEY: &str = "andes.extract.schema";
const SCHEMA_VERSION: &str = "4.0";
/// Rows per Parquet row group / record batch.
const BATCH_ROWS: usize = 65_536;

/// Why a (precursor, run) pair was measured.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExtractRole {
    /// Identified in this run; the prior is centred on its best PSM's RT.
    /// The target and its decoys carry this role.
    Identified,
    /// Not identified in this run; the prior comes from the alignment.
    Transfer,
    /// Identified in this run, measured as if it were a transfer (alignment
    /// prior without this run's observation).
    PseudoTransfer,
}

/// One fitted hypothesis of one library precursor in one run.
#[derive(Debug, Clone, PartialEq)]
pub struct ExtractRow {
    pub precursor: u32,
    pub run: u16,
    pub kind: Hypothesis,
    pub role: ExtractRole,
    /// `None` when the extraction window held too few scans.
    pub evidence: Option<Evidence>,
    /// Aligned runs with a PSM of the precursor, excluding this run except
    /// for an identified row.
    pub n_donor_runs: u16,
    /// Library confidence (internal quant gate, not an identification FDR).
    pub library_q: f64,
    /// The MS2 evidence of a transfer or pseudo-transfer row with a fit;
    /// `None` otherwise.
    pub ms2: Option<Ms2Evidence>,
}

const F64_EVIDENCE: [&str; 18] = [
    "amplitude",
    "amplitude_se",
    "quantity",
    "quantity_se",
    "llr",
    "rt_z",
    "width_ratio",
    "envelope_cos",
    "coherence",
    "m1_ratio",
    "mass_err_ppm",
    "t0",
    "w",
    "iso_h0",
    "iso_h1",
    "iso_h2",
    "iso_h3",
    "tailing",
];

fn f64_fields(e: &Evidence) -> [f64; 18] {
    [
        e.amplitude,
        e.amplitude_se,
        e.quantity,
        e.quantity_se,
        e.llr,
        e.rt_z,
        e.width_ratio,
        e.envelope_cos,
        e.coherence,
        e.m1_ratio,
        e.mass_err_ppm,
        e.t0,
        e.w,
        e.iso_heights[0],
        e.iso_heights[1],
        e.iso_heights[2],
        e.iso_heights[3],
        e.tailing,
    ]
}

/// Columns before the evidence block.
const N_KEY_COLUMNS: usize = 7;
/// The MS2 evidence columns, after the fit evidence.
const MS2_COLUMNS: [&str; 4] = ["ms2_n_scans", "ms2_n_matched", "ms2_angle", "ms2_explained"];
const N_COLUMNS: usize = N_KEY_COLUMNS + F64_EVIDENCE.len() + 2 + MS2_COLUMNS.len();

fn schema() -> Schema {
    let mut fields = vec![
        Field::new("precursor", DataType::UInt32, false),
        Field::new("run", DataType::UInt16, false),
        Field::new("kind", DataType::UInt8, false),
        Field::new("mass_shift", DataType::Int8, false),
        Field::new("role", DataType::UInt8, false),
        Field::new("n_donor_runs", DataType::UInt16, false),
        Field::new("library_q", DataType::Float64, false),
    ];
    fields.extend(
        F64_EVIDENCE
            .iter()
            .map(|n| Field::new(*n, DataType::Float64, true)),
    );
    fields.push(Field::new("n_points", DataType::UInt32, true));
    fields.push(Field::new("isotopes_used", DataType::UInt8, true));
    fields.push(Field::new(MS2_COLUMNS[0], DataType::UInt16, true));
    fields.push(Field::new(MS2_COLUMNS[1], DataType::UInt8, true));
    fields.push(Field::new(MS2_COLUMNS[2], DataType::Float32, true));
    fields.push(Field::new(MS2_COLUMNS[3], DataType::Float32, true));
    Schema::new(fields)
}

fn encode_kind(k: Hypothesis) -> (u8, i8) {
    match k {
        Hypothesis::Target => (0, 0),
        Hypothesis::MassShiftDecoy(n) => (1, n),
        Hypothesis::RtShiftDecoy => (2, 0),
    }
}

/// The hypothesis of a `(kind, mass_shift)` pair; `None` for an unknown kind
/// or a shift inconsistent with it (non-zero for a target or an RT-shift
/// decoy, zero for a mass-shift decoy).
fn decode_kind(code: u8, shift: i8) -> Option<Hypothesis> {
    match (code, shift) {
        (0, 0) => Some(Hypothesis::Target),
        (1, n) if n != 0 => Some(Hypothesis::MassShiftDecoy(n)),
        (2, 0) => Some(Hypothesis::RtShiftDecoy),
        _ => None,
    }
}

fn encode_role(r: ExtractRole) -> u8 {
    match r {
        ExtractRole::Identified => 0,
        ExtractRole::Transfer => 1,
        ExtractRole::PseudoTransfer => 2,
    }
}

fn decode_role(code: u8) -> Option<ExtractRole> {
    match code {
        0 => Some(ExtractRole::Identified),
        1 => Some(ExtractRole::Transfer),
        2 => Some(ExtractRole::PseudoTransfer),
        _ => None,
    }
}

fn batch(schema: &Arc<Schema>, rows: &[ExtractRow]) -> Result<RecordBatch, String> {
    let mut precursor = UInt32Builder::new();
    let mut run = UInt16Builder::new();
    let mut kind = UInt8Builder::new();
    let mut shift = Int8Builder::new();
    let mut role = UInt8Builder::new();
    let mut donors = UInt16Builder::new();
    let mut library_q = Float64Builder::new();
    let mut ev: Vec<Float64Builder> = F64_EVIDENCE.iter().map(|_| Float64Builder::new()).collect();
    let mut n_points = UInt32Builder::new();
    let mut isotopes = UInt8Builder::new();
    let mut ms2_scans = UInt16Builder::new();
    let mut ms2_matched = UInt8Builder::new();
    let mut ms2_angle = Float32Builder::new();
    let mut ms2_explained = Float32Builder::new();
    for r in rows {
        ms2_scans.append_option(r.ms2.map(|m| m.n_scans));
        ms2_matched.append_option(r.ms2.map(|m| m.best.n_matched));
        ms2_angle.append_option(r.ms2.map(|m| m.best.angle));
        ms2_explained.append_option(r.ms2.map(|m| m.best.explained));
        let (k, n) = encode_kind(r.kind);
        precursor.append_value(r.precursor);
        run.append_value(r.run);
        kind.append_value(k);
        shift.append_value(n);
        role.append_value(encode_role(r.role));
        donors.append_value(r.n_donor_runs);
        library_q.append_value(r.library_q);
        match &r.evidence {
            Some(e) => {
                for (b, v) in ev.iter_mut().zip(f64_fields(e)) {
                    b.append_value(v);
                }
                n_points.append_value(e.n_points);
                isotopes.append_value(e.isotopes_used);
            }
            None => {
                for b in ev.iter_mut() {
                    b.append_null();
                }
                n_points.append_null();
                isotopes.append_null();
            }
        }
    }
    let mut cols: Vec<ArrayRef> = vec![
        Arc::new(precursor.finish()),
        Arc::new(run.finish()),
        Arc::new(kind.finish()),
        Arc::new(shift.finish()),
        Arc::new(role.finish()),
        Arc::new(donors.finish()),
        Arc::new(library_q.finish()),
    ];
    cols.extend(ev.iter_mut().map(|b| Arc::new(b.finish()) as ArrayRef));
    cols.push(Arc::new(n_points.finish()));
    cols.push(Arc::new(isotopes.finish()));
    cols.push(Arc::new(ms2_scans.finish()));
    cols.push(Arc::new(ms2_matched.finish()));
    cols.push(Arc::new(ms2_angle.finish()));
    cols.push(Arc::new(ms2_explained.finish()));
    RecordBatch::try_new(schema.clone(), cols).map_err(|e| e.to_string())
}

fn write_versioned(
    path: &Path,
    rows: &[ExtractRow],
    version: &str,
    schema: Schema,
    tags: &[(String, String)],
) -> Result<(), String> {
    let mut partial = path.as_os_str().to_owned();
    partial.push(".partial");
    let partial = std::path::PathBuf::from(partial);
    let written = write_file(&partial, rows, version, schema, tags).and_then(|()| {
        std::fs::rename(&partial, path).map_err(|e| format!("{}: {e}", path.display()))
    });
    if written.is_err() {
        let _ = std::fs::remove_file(&partial);
    }
    written
}

fn write_file(
    path: &Path,
    rows: &[ExtractRow],
    version: &str,
    schema: Schema,
    tags: &[(String, String)],
) -> Result<(), String> {
    let err = |e: &dyn std::fmt::Display| format!("{}: {e}", path.display());
    let schema = Arc::new(schema);
    let props = WriterProperties::builder()
        .set_compression(Compression::ZSTD(ZstdLevel::default()))
        .set_max_row_group_size(BATCH_ROWS)
        .set_key_value_metadata(Some(
            std::iter::once(KeyValue::new(SCHEMA_KEY.to_string(), version.to_string()))
                .chain(
                    tags.iter()
                        .map(|(k, v)| KeyValue::new(k.clone(), v.clone())),
                )
                .collect(),
        ))
        .build();
    let file = File::create(path).map_err(|e| err(&e))?;
    let mut writer =
        ArrowWriter::try_new(file, schema.clone(), Some(props)).map_err(|e| err(&e))?;
    for chunk in rows.chunks(BATCH_ROWS) {
        let b = batch(&schema, chunk).map_err(|e| err(&e))?;
        writer.write(&b).map_err(|e| err(&e))?;
    }
    writer.close().map_err(|e| err(&e))?;
    Ok(())
}

/// Write the extraction table of one run.
pub fn write_extract_table(path: &Path, rows: &[ExtractRow]) -> Result<(), String> {
    write_versioned(path, rows, SCHEMA_VERSION, schema(), &[])
}

/// Write the extraction table of one run with extra key-value metadata
/// `tags` (read back by [`read_extract_tags`]).
pub fn write_extract_table_tagged(
    path: &Path,
    rows: &[ExtractRow],
    tags: &[(String, String)],
) -> Result<(), String> {
    write_versioned(path, rows, SCHEMA_VERSION, schema(), tags)
}

/// An extraction table written in pieces (e.g. a shard of several runs'
/// rows), with row groups of at most `row_group_rows` rows, to
/// `<path>.partial` and renamed into place by [`ExtractTableWriter::finish`].
/// Read back by [`read_extract_table`]. Dropped unfinished, the partial file
/// is removed.
pub struct ExtractTableWriter {
    path: std::path::PathBuf,
    partial: std::path::PathBuf,
    schema: Arc<Schema>,
    rows_per_group: usize,
    writer: Option<ArrowWriter<File>>,
}

impl ExtractTableWriter {
    pub fn create(
        path: &Path,
        tags: &[(String, String)],
        row_group_rows: usize,
    ) -> Result<Self, String> {
        let mut partial = path.as_os_str().to_owned();
        partial.push(".partial");
        let partial = std::path::PathBuf::from(partial);
        let err = |e: &dyn std::fmt::Display| format!("{}: {e}", partial.display());
        let rows_per_group = row_group_rows.max(1);
        let props = WriterProperties::builder()
            .set_compression(Compression::ZSTD(ZstdLevel::default()))
            .set_max_row_group_size(rows_per_group)
            .set_key_value_metadata(Some(
                std::iter::once(KeyValue::new(
                    SCHEMA_KEY.to_string(),
                    SCHEMA_VERSION.to_string(),
                ))
                .chain(
                    tags.iter()
                        .map(|(k, v)| KeyValue::new(k.clone(), v.clone())),
                )
                .collect(),
            ))
            .build();
        let schema = Arc::new(schema());
        let file = File::create(&partial).map_err(|e| err(&e))?;
        let writer = match ArrowWriter::try_new(file, schema.clone(), Some(props)) {
            Ok(w) => w,
            Err(e) => {
                let _ = std::fs::remove_file(&partial);
                return Err(err(&e));
            }
        };
        Ok(Self {
            path: path.to_path_buf(),
            partial,
            schema,
            rows_per_group,
            writer: Some(writer),
        })
    }

    pub fn write(&mut self, rows: &[ExtractRow]) -> Result<(), String> {
        let err = |e: &dyn std::fmt::Display| format!("{}: {e}", self.partial.display());
        let w = self.writer.as_mut().expect("writer is open until finish");
        for chunk in rows.chunks(self.rows_per_group) {
            let b = batch(&self.schema, chunk).map_err(|e| err(&e))?;
            w.write(&b).map_err(|e| err(&e))?;
        }
        Ok(())
    }

    pub fn finish(mut self) -> Result<(), String> {
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

impl Drop for ExtractTableWriter {
    fn drop(&mut self) {
        if self.writer.take().is_some() {
            let _ = std::fs::remove_file(&self.partial);
        }
    }
}

/// The number of rows of an extraction table, from its footer. Every error
/// names `path`.
pub fn extract_table_rows(path: &Path) -> Result<usize, String> {
    let err = |e: &dyn std::fmt::Display| format!("{}: {e}", path.display());
    let file = File::open(path).map_err(|e| err(&e))?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file).map_err(|e| err(&e))?;
    Ok(builder.metadata().file_metadata().num_rows().max(0) as usize)
}

/// The key-value metadata of an extraction table, without reading its rows.
/// Every error names `path`.
pub fn read_extract_tags(path: &Path) -> Result<Vec<(String, String)>, String> {
    let err = |e: &dyn std::fmt::Display| format!("{}: {e}", path.display());
    let file = File::open(path).map_err(|e| err(&e))?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file).map_err(|e| err(&e))?;
    Ok(builder
        .metadata()
        .file_metadata()
        .key_value_metadata()
        .cloned()
        .unwrap_or_default()
        .into_iter()
        .filter_map(|kv| kv.value.map(|v| (kv.key, v)))
        .collect())
}

#[cfg(test)]
pub(crate) fn write_with_schema_version(
    path: &Path,
    rows: &[ExtractRow],
    version: &str,
) -> Result<(), String> {
    write_versioned(path, rows, version, schema(), &[])
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

/// Read an extraction table. Every error names `path`.
pub fn read_extract_table(path: &Path) -> Result<Vec<ExtractRow>, String> {
    let mut rows = Vec::new();
    for b in ExtractTableReader::open(path, BATCH_ROWS)? {
        rows.extend(b?);
    }
    Ok(rows)
}

/// An extraction table read in record batches of at most `batch_rows` rows,
/// in file order (see [`read_extract_table`]). Every error names the file.
pub struct ExtractTableReader {
    path: std::path::PathBuf,
    reader: parquet::arrow::arrow_reader::ParquetRecordBatchReader,
    /// Read by [`ExtractTableReader::open_identified_fits`].
    fits_only: bool,
}

/// The evidence columns of [`ExtractTableReader::open_identified_fits`].
const FIT_EVIDENCE: [&str; 8] = [
    "amplitude",
    "t0",
    "w",
    "iso_h0",
    "iso_h1",
    "iso_h2",
    "iso_h3",
    "tailing",
];

impl ExtractTableReader {
    /// Opens `path` and checks its schema version.
    pub fn open(path: &Path, batch_rows: usize) -> Result<Self, String> {
        Self::open_with(path, batch_rows, false)
    }

    /// Opens `path` to read its identified target rows only, with the
    /// columns of their peak shape: the keys, and of the evidence the
    /// amplitude, apex, width, isotope heights and tailing. The other
    /// evidence values are NaN (counts 0) and the MS2 evidence is `None`.
    pub fn open_identified_fits(path: &Path, batch_rows: usize) -> Result<Self, String> {
        Self::open_with(path, batch_rows, true)
    }

    fn open_with(path: &Path, batch_rows: usize, fits_only: bool) -> Result<Self, String> {
        let err = |e: &dyn std::fmt::Display| format!("{}: {e}", path.display());
        let file = File::open(path).map_err(|e| err(&e))?;
        let builder = ParquetRecordBatchReaderBuilder::try_new(file).map_err(|e| err(&e))?;
        let meta = builder
            .metadata()
            .file_metadata()
            .key_value_metadata()
            .cloned()
            .unwrap_or_default();
        let version = meta
            .iter()
            .find(|kv| kv.key == SCHEMA_KEY)
            .and_then(|kv| kv.value.clone())
            .ok_or_else(|| err(&format!("missing metadata key {SCHEMA_KEY}")))?;
        let major: u32 = version
            .split('.')
            .next()
            .and_then(|m| m.parse().ok())
            .ok_or_else(|| err(&format!("unparseable schema version {version:?}")))?;
        if major != EXTRACT_SCHEMA_MAJOR {
            return Err(err(&format!(
                "unsupported extract table schema {version} (this build reads major {EXTRACT_SCHEMA_MAJOR}); run `andes quant extract` again"
            )));
        }
        let builder = builder.with_batch_size(batch_rows.max(1));
        let builder = if fits_only {
            let names: Vec<String> = schema()
                .fields()
                .iter()
                .take(N_KEY_COLUMNS)
                .map(|f| f.name().clone())
                .chain(FIT_EVIDENCE.iter().map(|n| n.to_string()))
                .collect();
            let leaves: Vec<usize> = names
                .iter()
                .map(|n| {
                    builder
                        .parquet_schema()
                        .columns()
                        .iter()
                        .position(|c| c.name() == n)
                        .ok_or_else(|| err(&format!("no column {n:?}")))
                })
                .collect::<Result<_, _>>()?;
            let mask = parquet::arrow::ProjectionMask::leaves(builder.parquet_schema(), leaves);
            builder.with_projection(mask)
        } else {
            builder
        };
        let reader = builder.build().map_err(|e| err(&e))?;
        Ok(Self {
            path: path.to_path_buf(),
            reader,
            fits_only,
        })
    }
}

impl Iterator for ExtractTableReader {
    type Item = Result<Vec<ExtractRow>, String>;

    fn next(&mut self) -> Option<Self::Item> {
        let b = self.reader.next()?;
        let path = self.path.as_path();
        let err = |e: &dyn std::fmt::Display| format!("{}: {e}", path.display());
        let fits_only = self.fits_only;
        Some(b.map_err(|e| err(&e)).and_then(|b| {
            if fits_only {
                decode_fit_batch(&b, path)
            } else {
                decode_batch(&b, path)
            }
        }))
    }
}

/// The identified target rows of one record batch read by
/// [`ExtractTableReader::open_identified_fits`].
fn decode_fit_batch(b: &RecordBatch, path: &Path) -> Result<Vec<ExtractRow>, String> {
    let err = |e: &dyn std::fmt::Display| format!("{}: {e}", path.display());
    let at = |name: &str| {
        b.schema()
            .index_of(name)
            .map_err(|_| err(&format!("no column {name:?}")))
    };
    let precursor = column::<UInt32Array>(b, at("precursor")?, path)?;
    let run = column::<UInt16Array>(b, at("run")?, path)?;
    let kind = column::<UInt8Array>(b, at("kind")?, path)?;
    let shift = column::<Int8Array>(b, at("mass_shift")?, path)?;
    let role = column::<UInt8Array>(b, at("role")?, path)?;
    let donors = column::<UInt16Array>(b, at("n_donor_runs")?, path)?;
    let library_q = column::<Float64Array>(b, at("library_q")?, path)?;
    let ev: Vec<&Float64Array> = FIT_EVIDENCE
        .iter()
        .map(|n| column::<Float64Array>(b, at(n)?, path))
        .collect::<Result<_, _>>()?;
    let mut rows = Vec::new();
    for i in 0..b.num_rows() {
        let kind = decode_kind(kind.value(i), shift.value(i)).ok_or_else(|| {
            err(&format!(
                "corrupt row: kind {} with mass_shift {}",
                kind.value(i),
                shift.value(i)
            ))
        })?;
        let role = decode_role(role.value(i))
            .ok_or_else(|| err(&format!("unknown role code {}", role.value(i))))?;
        if (role, kind) != (ExtractRole::Identified, Hypothesis::Target) {
            continue;
        }
        let evidence = (!ev[0].is_null(i)).then(|| {
            let v = |k: usize| ev[k].value(i);
            Evidence {
                amplitude: v(0),
                amplitude_se: f64::NAN,
                quantity: f64::NAN,
                quantity_se: f64::NAN,
                llr: f64::NAN,
                rt_z: f64::NAN,
                width_ratio: f64::NAN,
                envelope_cos: f64::NAN,
                coherence: f64::NAN,
                m1_ratio: f64::NAN,
                mass_err_ppm: f64::NAN,
                n_points: 0,
                isotopes_used: 0,
                t0: v(1),
                w: v(2),
                iso_heights: [v(3), v(4), v(5), v(6)],
                tailing: v(7),
            }
        });
        rows.push(ExtractRow {
            precursor: precursor.value(i),
            run: run.value(i),
            kind,
            role,
            evidence,
            n_donor_runs: donors.value(i),
            library_q: library_q.value(i),
            ms2: None,
        });
    }
    Ok(rows)
}

/// The rows of one record batch of an extraction table at `path`.
fn decode_batch(b: &RecordBatch, path: &Path) -> Result<Vec<ExtractRow>, String> {
    let err = |e: &dyn std::fmt::Display| format!("{}: {e}", path.display());
    let mut rows = Vec::with_capacity(b.num_rows());
    if b.num_columns() != N_COLUMNS {
        return Err(err(&format!(
            "expected {N_COLUMNS} columns, found {}",
            b.num_columns()
        )));
    }
    let expected = schema();
    for (got, want) in b.schema().fields().iter().zip(expected.fields()) {
        if got.name() != want.name() {
            return Err(err(&format!(
                "column {:?} where {:?} was expected",
                got.name(),
                want.name()
            )));
        }
    }
    let precursor = column::<UInt32Array>(b, 0, path)?;
    let run = column::<UInt16Array>(b, 1, path)?;
    let kind = column::<UInt8Array>(b, 2, path)?;
    let shift = column::<Int8Array>(b, 3, path)?;
    let role = column::<UInt8Array>(b, 4, path)?;
    let donors = column::<UInt16Array>(b, 5, path)?;
    let library_q = column::<Float64Array>(b, 6, path)?;
    let ev: Vec<&Float64Array> = (0..F64_EVIDENCE.len())
        .map(|k| column::<Float64Array>(b, N_KEY_COLUMNS + k, path))
        .collect::<Result<_, _>>()?;
    let after = N_KEY_COLUMNS + F64_EVIDENCE.len();
    let n_points = column::<UInt32Array>(b, after, path)?;
    let isotopes = column::<UInt8Array>(b, after + 1, path)?;
    let ms2_scans = column::<UInt16Array>(b, after + 2, path)?;
    let ms2_matched = column::<UInt8Array>(b, after + 3, path)?;
    let ms2_angle = column::<Float32Array>(b, after + 4, path)?;
    let ms2_explained = column::<Float32Array>(b, after + 5, path)?;

    for i in 0..b.num_rows() {
        let kind = decode_kind(kind.value(i), shift.value(i)).ok_or_else(|| {
            err(&format!(
                "corrupt row: kind {} with mass_shift {}",
                kind.value(i),
                shift.value(i)
            ))
        })?;
        let role = decode_role(role.value(i))
            .ok_or_else(|| err(&format!("unknown role code {}", role.value(i))))?;
        let evidence = if ev[0].is_null(i) {
            None
        } else {
            let v = |k: usize| ev[k].value(i);
            Some(Evidence {
                amplitude: v(0),
                amplitude_se: v(1),
                quantity: v(2),
                quantity_se: v(3),
                llr: v(4),
                rt_z: v(5),
                width_ratio: v(6),
                envelope_cos: v(7),
                coherence: v(8),
                m1_ratio: v(9),
                mass_err_ppm: v(10),
                n_points: n_points.value(i),
                isotopes_used: isotopes.value(i),
                t0: v(11),
                w: v(12),
                iso_heights: [v(13), v(14), v(15), v(16)],
                tailing: v(17),
            })
        };
        let ms2 = if ms2_scans.is_null(i) {
            None
        } else {
            if ms2_matched.is_null(i) || ms2_angle.is_null(i) || ms2_explained.is_null(i) {
                return Err(err(&"corrupt row: partial MS2 evidence"));
            }
            Some(Ms2Evidence {
                n_scans: ms2_scans.value(i),
                best: ScanMatch {
                    n_matched: ms2_matched.value(i),
                    angle: ms2_angle.value(i),
                    explained: ms2_explained.value(i),
                },
            })
        };
        rows.push(ExtractRow {
            precursor: precursor.value(i),
            run: run.value(i),
            kind,
            role,
            evidence,
            n_donor_runs: donors.value(i),
            library_q: library_q.value(i),
            ms2,
        });
    }
    Ok(rows)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn evidence(i: usize) -> Evidence {
        let f = i as f64;
        Evidence {
            amplitude: 1e5 + f,
            amplitude_se: 12.5,
            quantity: 1.2e6,
            quantity_se: 150.0,
            llr: 88.0 + f,
            rt_z: -0.4,
            width_ratio: 1.4,
            envelope_cos: 0.97,
            coherence: 0.91,
            m1_ratio: 0.02,
            mass_err_ppm: -1.5,
            n_points: 17 + i as u32,
            isotopes_used: 0b0111,
            t0: 1234.5,
            w: 4.2,
            iso_heights: [5.0 + f, 3.0, 1.0, -0.5],
            tailing: 0.1 * f,
        }
    }

    fn rows() -> Vec<ExtractRow> {
        let kinds = [
            Hypothesis::Target,
            Hypothesis::MassShiftDecoy(-11),
            Hypothesis::RtShiftDecoy,
            Hypothesis::MassShiftDecoy(7),
        ];
        let roles = [
            ExtractRole::Identified,
            ExtractRole::Transfer,
            ExtractRole::PseudoTransfer,
        ];
        (0..9)
            .map(|i| ExtractRow {
                precursor: 1000 + i as u32,
                run: 3,
                kind: kinds[i % kinds.len()],
                role: roles[i % roles.len()],
                evidence: (i % 4 != 2).then(|| evidence(i)),
                n_donor_runs: i as u16,
                library_q: 0.001 * i as f64,
                ms2: (i % 3 == 1).then_some(Ms2Evidence {
                    n_scans: i as u16,
                    best: ScanMatch {
                        n_matched: (i % 5) as u8,
                        angle: 0.123 * i as f32,
                        explained: 0.5,
                    },
                }),
            })
            .collect()
    }

    #[test]
    fn round_trip_preserves_rows_kinds_roles_and_nulls() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("a.extract.parquet");
        let rows = rows();
        write_extract_table(&p, &rows).unwrap();
        assert_eq!(read_extract_table(&p).unwrap(), rows);

        let e = dir.path().join("e.extract.parquet");
        write_extract_table(&e, &[]).unwrap();
        assert_eq!(read_extract_table(&e).unwrap(), vec![]);
    }

    #[test]
    fn identified_fits_reader_keeps_the_identified_targets_shape_columns() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("f.extract.parquet");
        let mut rows = rows();
        for (i, r) in rows.iter_mut().enumerate() {
            if i % 3 == 0 {
                (r.kind, r.role) = (Hypothesis::Target, ExtractRole::Identified);
            }
        }
        write_extract_table(&p, &rows).unwrap();
        let mut got = Vec::new();
        for b in ExtractTableReader::open_identified_fits(&p, 2).unwrap() {
            got.extend(b.unwrap());
        }
        let want: Vec<&ExtractRow> = rows
            .iter()
            .filter(|r| (r.role, r.kind) == (ExtractRole::Identified, Hypothesis::Target))
            .collect();
        assert_eq!(got.len(), want.len());
        assert!(want.iter().any(|r| r.evidence.is_none()));
        let keys = |r: &ExtractRow| (r.precursor, r.run, r.kind, r.role, r.n_donor_runs);
        for (g, w) in got.iter().zip(want) {
            assert_eq!(keys(g), keys(w));
            assert_eq!(g.library_q, w.library_q);
            assert_eq!(g.ms2, None);
            match (&g.evidence, &w.evidence) {
                (None, None) => {}
                (Some(a), Some(b)) => {
                    assert_eq!(
                        (a.amplitude, a.t0, a.w, a.iso_heights, a.tailing),
                        (b.amplitude, b.t0, b.w, b.iso_heights, b.tailing)
                    );
                    assert!(a.llr.is_nan() && a.quantity.is_nan());
                }
                other => panic!("{other:?}"),
            }
        }
    }

    #[test]
    fn extract_table_rejects_unknown_major_version() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("b.extract.parquet");
        for version in ["5.0", "3.0"] {
            write_with_schema_version(&p, &rows(), version).unwrap();
            let e = read_extract_table(&p).unwrap_err();
            assert!(
                e.contains("b.extract.parquet") && e.contains(&format!("schema {version}")),
                "{e}"
            );
        }
    }

    #[test]
    fn missing_and_non_parquet_files_name_the_path() {
        let dir = tempfile::tempdir().unwrap();
        let missing = dir.path().join("nope.extract.parquet");
        let e = read_extract_table(&missing).unwrap_err();
        assert!(e.contains("nope.extract.parquet"), "{e}");

        let junk = dir.path().join("junk.extract.parquet");
        std::fs::write(&junk, b"not parquet").unwrap();
        let e = read_extract_table(&junk).unwrap_err();
        assert!(e.contains("junk.extract.parquet"), "{e}");
    }

    #[test]
    fn write_is_atomic_and_leaves_no_partial_file() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("c.extract.parquet");
        write_extract_table(&p, &rows()).unwrap();
        assert!(p.exists());
        assert!(!dir.path().join("c.extract.parquet.partial").exists());
        // A failed write leaves neither file.
        let bad = dir.path().join("missing-dir").join("d.extract.parquet");
        let e = write_extract_table(&bad, &rows()).unwrap_err();
        assert!(e.contains("d.extract.parquet"), "{e}");
        assert!(!bad.exists());
    }

    #[test]
    fn renamed_column_is_rejected_naming_the_file() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("r.extract.parquet");
        let mut fields: Vec<Field> = schema()
            .fields()
            .iter()
            .map(|f| f.as_ref().clone())
            .collect();
        fields[4] = Field::new("rolle", DataType::UInt8, false);
        write_versioned(&p, &rows(), SCHEMA_VERSION, Schema::new(fields), &[]).unwrap();
        let e = read_extract_table(&p).unwrap_err();
        assert!(
            e.contains("r.extract.parquet") && e.contains("rolle"),
            "{e}"
        );
    }

    #[test]
    fn inconsistent_kind_and_mass_shift_is_corrupt() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("k.extract.parquet");
        let mut bad = rows();
        bad[1].kind = Hypothesis::MassShiftDecoy(0);
        write_extract_table(&p, &bad).unwrap();
        let e = read_extract_table(&p).unwrap_err();
        assert!(
            e.contains("k.extract.parquet") && e.contains("corrupt"),
            "{e}"
        );
    }
}
