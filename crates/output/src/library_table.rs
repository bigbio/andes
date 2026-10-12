//! The match-between-runs library on disk: `library.parquet` (one row per
//! precursor id: identity, proteins and, for library entries, the entry
//! fields) and `alignment.parquet` (one row per run: the run map, residual
//! widths, anchor statistics and the run's observations). Together they
//! rebuild the [`Library`] exactly, so an extract in another process gives the
//! same rows as one in the process that built the library.
//!
//! Each file carries a schema version in its key-value metadata; a reader
//! rejects a different major version. `library.parquet` also records the run
//! stems in run-index order and the library entry gate. Both files are written
//! to `<path>.partial` and renamed into place once complete.

use std::fs::File;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use arrow::array::{
    Array, ArrayRef, BooleanArray, BooleanBuilder, Float32Array, Float32Builder, Float64Array,
    Float64Builder, ListArray, ListBuilder, StringArray, StringBuilder, UInt16Array, UInt16Builder,
    UInt32Array, UInt32Builder, UInt64Array, UInt64Builder, UInt8Array, UInt8Builder,
};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::ArrowWriter;
use parquet::basic::{Compression, ZstdLevel};
use parquet::file::properties::WriterProperties;
use parquet::format::KeyValue;
use quant::align::{Alignment, RunAlignment, RunMap};
use quant::library::{Library, LibraryEntry, PrecursorMeta};

/// Major version of both layouts. A reader rejects any other major.
pub const LIBRARY_SCHEMA_MAJOR: u32 = 1;
/// File name of the precursor table inside the quant directory.
pub const LIBRARY_FILE: &str = "library.parquet";
/// File name of the alignment table inside the quant directory.
pub const ALIGNMENT_FILE: &str = "alignment.parquet";

const LIBRARY_KEY: &str = "andes.library.schema";
const ALIGNMENT_KEY: &str = "andes.alignment.schema";
const RUNS_KEY: &str = "andes.library.runs";
const MAX_Q_KEY: &str = "andes.library.max_q";
/// Key of [`library_fingerprint`] in `library.parquet`; extract tables copy it.
pub const FINGERPRINT_KEY: &str = "andes.library.fingerprint";
const ALIGNMENT_FINGERPRINT_KEY: &str = "andes.alignment.fingerprint";
const SCHEMA_VERSION: &str = "1.0";
/// Rows per Parquet row group / record batch of the precursor table.
const BATCH_ROWS: usize = 65_536;
const ENVELOPE_LEN: usize = 4;

/// A library with what the stages around it need: the run stems in
/// run-index order, the entry gate it was built with, and the `;`-joined
/// protein accessions per precursor id.
#[derive(Debug, Clone)]
pub struct StoredLibrary {
    pub runs: Vec<String>,
    /// Library entry gate (internal precursor-level gate, not an
    /// identification FDR).
    pub max_q: f64,
    pub library: Library,
    /// Per precursor id; same length as `library.meta`.
    pub proteins: Vec<String>,
}

fn list(item: DataType) -> DataType {
    DataType::List(Arc::new(Field::new("item", item, true)))
}

fn library_schema() -> Schema {
    Schema::new(vec![
        Field::new("key", DataType::Utf8, false),
        Field::new("mz", DataType::Float64, false),
        Field::new("charge", DataType::UInt8, false),
        Field::new("envelope", list(DataType::Float32), false),
        Field::new("proteins", DataType::Utf8, false),
        Field::new("entry_q", DataType::Float64, true),
        Field::new("faims_cv", DataType::Float32, true),
        Field::new("donors", list(DataType::UInt16), false),
        Field::new("best_rt_run", list(DataType::UInt16), false),
        Field::new("best_rt", list(DataType::Float64), false),
    ])
}

fn alignment_schema() -> Schema {
    Schema::new(vec![
        Field::new("run", DataType::UInt16, false),
        Field::new("stem", DataType::Utf8, false),
        Field::new("aligned", DataType::Boolean, false),
        Field::new("knots_x", list(DataType::Float64), false),
        Field::new("knots_y", list(DataType::Float64), false),
        Field::new("slope_lo", DataType::Float64, true),
        Field::new("slope_hi", DataType::Float64, true),
        Field::new("sd_t", list(DataType::Float64), false),
        Field::new("sd", list(DataType::Float64), false),
        Field::new("anchors", DataType::UInt64, false),
        Field::new("coverage", DataType::Float64, false),
        Field::new("obs_precursor", list(DataType::UInt32), false),
        Field::new("obs_rt", list(DataType::Float64), false),
    ])
}

/// 64-bit FNV-1a over the bytes fed to it.
struct Fnv(u64);

impl Fnv {
    fn new() -> Self {
        Fnv(0xcbf2_9ce4_8422_2325)
    }
    fn bytes(&mut self, b: &[u8]) {
        for &x in b {
            self.0 ^= x as u64;
            self.0 = self.0.wrapping_mul(0x0000_0100_0000_01b3);
        }
    }
    fn u64(&mut self, v: u64) {
        self.bytes(&v.to_le_bytes());
    }
    fn f64(&mut self, v: f64) {
        self.u64(v.to_bits());
    }
    fn str(&mut self, v: &str) {
        self.u64(v.len() as u64);
        self.bytes(v.as_bytes());
    }
    fn hex(&self) -> String {
        format!("{:016x}", self.0)
    }
}

/// A stable hash of an alignment's stored state (maps, residual widths,
/// anchor statistics and observations of every run).
pub fn alignment_fingerprint(a: &Alignment) -> String {
    let mut h = Fnv::new();
    let runs = a.runs();
    h.u64(runs.len() as u64);
    for r in &runs {
        match &r.map {
            Some(m) => {
                let (x, y) = m.knots();
                let (lo, hi) = m.end_slopes();
                h.u64(x.len() as u64);
                x.iter().chain(y).for_each(|&v| h.f64(v));
                h.f64(lo);
                h.f64(hi);
            }
            None => h.u64(u64::MAX),
        }
        let sd = r.residual_sd.as_deref().unwrap_or(&[]);
        h.u64(sd.len() as u64);
        for &(t, v) in sd {
            h.f64(t);
            h.f64(v);
        }
        h.u64(r.anchors as u64);
        h.f64(r.coverage);
        h.u64(r.observations.len() as u64);
        for &(p, t) in &r.observations {
            h.u64(p as u64);
            h.f64(t);
        }
    }
    h.hex()
}

/// A stable hash of a stored library: run stems, entry gate, precursor
/// metadata and proteins, entries and the alignment.
pub fn library_fingerprint(lib: &StoredLibrary) -> String {
    let mut h = Fnv::new();
    h.u64(lib.runs.len() as u64);
    lib.runs.iter().for_each(|r| h.str(r));
    h.f64(lib.max_q);
    h.u64(lib.library.meta.len() as u64);
    for (m, p) in lib.library.meta.iter().zip(&lib.proteins) {
        h.str(&m.key);
        h.f64(m.mz);
        h.u64(m.charge as u64);
        m.envelope.iter().for_each(|&v| h.u64(v.to_bits() as u64));
        h.str(p);
    }
    h.u64(lib.library.entries.len() as u64);
    for e in &lib.library.entries {
        h.u64(e.precursor as u64);
        h.f64(e.q);
        h.u64(e.faims_cv.map_or(u64::MAX, |v| v.to_bits() as u64));
        h.u64(e.donors.len() as u64);
        e.donors.iter().for_each(|&d| h.u64(d as u64));
        for &(r, t) in &e.best_rt {
            h.u64(r as u64);
            h.f64(t);
        }
    }
    h.str(&alignment_fingerprint(&lib.library.alignment));
    h.hex()
}

/// `<dir>/library.parquet` and `<dir>/alignment.parquet`.
pub fn library_paths(dir: &Path) -> (PathBuf, PathBuf) {
    (dir.join(LIBRARY_FILE), dir.join(ALIGNMENT_FILE))
}

fn library_batch(
    schema: &Arc<Schema>,
    lib: &StoredLibrary,
    entry_of: &[Option<&LibraryEntry>],
    ids: std::ops::Range<usize>,
) -> Result<RecordBatch, String> {
    let mut key = StringBuilder::new();
    let mut mz = Float64Builder::new();
    let mut charge = UInt8Builder::new();
    let mut env = ListBuilder::new(Float32Builder::new());
    let mut prot = StringBuilder::new();
    let mut q = Float64Builder::new();
    let mut cv = Float32Builder::new();
    let mut donors = ListBuilder::new(UInt16Builder::new());
    let mut best_run = ListBuilder::new(UInt16Builder::new());
    let mut best_rt = ListBuilder::new(Float64Builder::new());
    for id in ids {
        let m = &lib.library.meta[id];
        key.append_value(&m.key);
        mz.append_value(m.mz);
        charge.append_value(m.charge);
        env.values().append_slice(&m.envelope);
        env.append(true);
        prot.append_value(&lib.proteins[id]);
        let e = entry_of[id];
        q.append_option(e.map(|e| e.q));
        cv.append_option(e.and_then(|e| e.faims_cv));
        if let Some(e) = e {
            donors.values().append_slice(&e.donors);
            for &(r, t) in &e.best_rt {
                best_run.values().append_value(r);
                best_rt.values().append_value(t);
            }
        }
        donors.append(true);
        best_run.append(true);
        best_rt.append(true);
    }
    let cols: Vec<ArrayRef> = vec![
        Arc::new(key.finish()),
        Arc::new(mz.finish()),
        Arc::new(charge.finish()),
        Arc::new(env.finish()),
        Arc::new(prot.finish()),
        Arc::new(q.finish()),
        Arc::new(cv.finish()),
        Arc::new(donors.finish()),
        Arc::new(best_run.finish()),
        Arc::new(best_rt.finish()),
    ];
    RecordBatch::try_new(schema.clone(), cols).map_err(|e| e.to_string())
}

fn alignment_batch(
    schema: &Arc<Schema>,
    stems: &[String],
    runs: &[RunAlignment],
) -> Result<RecordBatch, String> {
    let mut run = UInt16Builder::new();
    let mut stem = StringBuilder::new();
    let mut aligned = BooleanBuilder::new();
    let mut kx = ListBuilder::new(Float64Builder::new());
    let mut ky = ListBuilder::new(Float64Builder::new());
    let mut slo = Float64Builder::new();
    let mut shi = Float64Builder::new();
    let mut sd_t = ListBuilder::new(Float64Builder::new());
    let mut sd = ListBuilder::new(Float64Builder::new());
    let mut anchors = UInt64Builder::new();
    let mut coverage = Float64Builder::new();
    let mut obs_p = ListBuilder::new(UInt32Builder::new());
    let mut obs_rt = ListBuilder::new(Float64Builder::new());
    for (r, ra) in runs.iter().enumerate() {
        run.append_value(r as u16);
        stem.append_value(&stems[r]);
        aligned.append_value(ra.map.is_some());
        match &ra.map {
            Some(m) => {
                let (x, y) = m.knots();
                kx.values().append_slice(x);
                ky.values().append_slice(y);
                let (lo, hi) = m.end_slopes();
                slo.append_value(lo);
                shi.append_value(hi);
            }
            None => {
                slo.append_null();
                shi.append_null();
            }
        }
        kx.append(true);
        ky.append(true);
        for &(t, s) in ra.residual_sd.iter().flatten() {
            sd_t.values().append_value(t);
            sd.values().append_value(s);
        }
        sd_t.append(true);
        sd.append(true);
        anchors.append_value(ra.anchors as u64);
        coverage.append_value(ra.coverage);
        for &(p, t) in &ra.observations {
            obs_p.values().append_value(p);
            obs_rt.values().append_value(t);
        }
        obs_p.append(true);
        obs_rt.append(true);
    }
    let cols: Vec<ArrayRef> = vec![
        Arc::new(run.finish()),
        Arc::new(stem.finish()),
        Arc::new(aligned.finish()),
        Arc::new(kx.finish()),
        Arc::new(ky.finish()),
        Arc::new(slo.finish()),
        Arc::new(shi.finish()),
        Arc::new(sd_t.finish()),
        Arc::new(sd.finish()),
        Arc::new(anchors.finish()),
        Arc::new(coverage.finish()),
        Arc::new(obs_p.finish()),
        Arc::new(obs_rt.finish()),
    ];
    RecordBatch::try_new(schema.clone(), cols).map_err(|e| e.to_string())
}

/// Writes one Parquet file through `<path>.partial`; `batches` produces the
/// record batches in order.
fn write_parquet(
    path: &Path,
    schema: Arc<Schema>,
    meta: Vec<KeyValue>,
    batches: impl Iterator<Item = Result<RecordBatch, String>>,
) -> Result<(), String> {
    let err = |e: &dyn std::fmt::Display| format!("{}: {e}", path.display());
    let mut partial = path.as_os_str().to_owned();
    partial.push(".partial");
    let partial = PathBuf::from(partial);
    let written = (|| {
        let props = WriterProperties::builder()
            .set_compression(Compression::ZSTD(ZstdLevel::default()))
            .set_max_row_group_size(BATCH_ROWS)
            .set_key_value_metadata(Some(meta))
            .build();
        let file = File::create(&partial).map_err(|e| err(&e))?;
        let mut writer = ArrowWriter::try_new(file, schema, Some(props)).map_err(|e| err(&e))?;
        for b in batches {
            writer
                .write(&b.map_err(|e| err(&e))?)
                .map_err(|e| err(&e))?;
        }
        writer.close().map_err(|e| err(&e))?;
        std::fs::rename(&partial, path).map_err(|e| err(&e))
    })();
    if written.is_err() {
        let _ = std::fs::remove_file(&partial);
    }
    written
}

fn write_versioned(dir: &Path, lib: &StoredLibrary, version: &str) -> Result<(), String> {
    let (lib_path, align_path) = library_paths(dir);
    let n = lib.library.meta.len();
    if lib.proteins.len() != n {
        return Err(format!(
            "{}: {} protein lists for {n} precursors",
            lib_path.display(),
            lib.proteins.len()
        ));
    }
    if lib.library.alignment.maps.len() != lib.runs.len() {
        return Err(format!(
            "{}: alignment of {} runs for {} run stems",
            align_path.display(),
            lib.library.alignment.maps.len(),
            lib.runs.len()
        ));
    }
    if let Some(s) = lib.runs.iter().find(|s| s.is_empty() || s.contains('\n')) {
        return Err(format!(
            "{}: run stem {s:?} is empty or holds a line break",
            lib_path.display()
        ));
    }
    if lib
        .library
        .meta
        .iter()
        .any(|m| m.envelope.len() != ENVELOPE_LEN)
    {
        return Err(format!(
            "{}: every envelope must have {ENVELOPE_LEN} isotopes",
            lib_path.display()
        ));
    }
    let mut entry_of: Vec<Option<&LibraryEntry>> = vec![None; n];
    for e in &lib.library.entries {
        let slot = entry_of.get_mut(e.precursor as usize).ok_or_else(|| {
            format!(
                "{}: entry {} outside the {n} precursors",
                lib_path.display(),
                e.precursor
            )
        })?;
        if slot.is_some() {
            return Err(format!(
                "{}: precursor {} has two entries",
                lib_path.display(),
                e.precursor
            ));
        }
        *slot = Some(e);
    }

    let runs = lib.library.alignment.runs();
    let align_schema = Arc::new(alignment_schema());
    let batch = alignment_batch(&align_schema, &lib.runs, &runs);
    write_parquet(
        &align_path,
        align_schema,
        vec![KeyValue::new(
            ALIGNMENT_KEY.to_string(),
            version.to_string(),
        )],
        std::iter::once(batch),
    )?;

    let lib_schema = Arc::new(library_schema());
    let chunks = (0..n.div_ceil(BATCH_ROWS)).map(|c| c * BATCH_ROWS..((c + 1) * BATCH_ROWS).min(n));
    write_parquet(
        &lib_path,
        lib_schema.clone(),
        vec![
            KeyValue::new(LIBRARY_KEY.to_string(), version.to_string()),
            KeyValue::new(RUNS_KEY.to_string(), lib.runs.join("\n")),
            KeyValue::new(MAX_Q_KEY.to_string(), lib.max_q.to_string()),
            KeyValue::new(FINGERPRINT_KEY.to_string(), library_fingerprint(lib)),
            KeyValue::new(
                ALIGNMENT_FINGERPRINT_KEY.to_string(),
                alignment_fingerprint(&lib.library.alignment),
            ),
        ],
        chunks.map(|ids| library_batch(&lib_schema, lib, &entry_of, ids)),
    )
}

/// Writes `<dir>/library.parquet` and `<dir>/alignment.parquet`.
pub fn write_library(dir: &Path, lib: &StoredLibrary) -> Result<(), String> {
    write_versioned(dir, lib, SCHEMA_VERSION)
}

#[cfg(test)]
pub(crate) fn write_with_schema_version(
    dir: &Path,
    lib: &StoredLibrary,
    version: &str,
) -> Result<(), String> {
    write_versioned(dir, lib, version)
}

/// Opens a versioned table: checks the schema version under `key` and the
/// column names and types. Returns the reader and the key-value metadata.
fn open(
    path: &Path,
    key: &str,
    what: &str,
    expected: &Schema,
) -> Result<
    (
        parquet::arrow::arrow_reader::ParquetRecordBatchReader,
        Vec<KeyValue>,
    ),
    String,
> {
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
        .find(|kv| kv.key == key)
        .and_then(|kv| kv.value.clone())
        .ok_or_else(|| err(&format!("missing metadata key {key}")))?;
    let major: u32 = version
        .split('.')
        .next()
        .and_then(|m| m.parse().ok())
        .ok_or_else(|| err(&format!("unparseable schema version {version:?}")))?;
    if major != LIBRARY_SCHEMA_MAJOR {
        return Err(err(&format!(
            "unsupported {what} schema {version} (this build reads major {LIBRARY_SCHEMA_MAJOR})"
        )));
    }
    let got = builder.schema();
    let names = |s: &Schema| -> Vec<(String, DataType)> {
        s.fields()
            .iter()
            .map(|f| (f.name().clone(), f.data_type().clone()))
            .collect()
    };
    if names(got) != names(expected) {
        return Err(err(&format!(
            "columns {:?} differ from the expected {:?}",
            got.fields().iter().map(|f| f.name()).collect::<Vec<_>>(),
            expected
                .fields()
                .iter()
                .map(|f| f.name())
                .collect::<Vec<_>>()
        )));
    }
    let reader = builder.build().map_err(|e| err(&e))?;
    Ok((reader, meta))
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

/// Values of row `i` of a list column.
fn items<T: Clone + 'static, A: Array + 'static>(
    l: &ListArray,
    i: usize,
    path: &Path,
    values: impl Fn(&A) -> &[T],
) -> Result<Vec<T>, String> {
    let v = l.value(i);
    let a = v.as_any().downcast_ref::<A>().ok_or_else(|| {
        format!(
            "{}: list items have unexpected type {}",
            path.display(),
            v.data_type()
        )
    })?;
    if a.null_count() > 0 {
        return Err(format!("{}: list items hold nulls", path.display()));
    }
    Ok(values(a).to_vec())
}

fn f64s(l: &ListArray, i: usize, path: &Path) -> Result<Vec<f64>, String> {
    items::<f64, Float64Array>(l, i, path, |a| a.values())
}

fn u16s(l: &ListArray, i: usize, path: &Path) -> Result<Vec<u16>, String> {
    items::<u16, UInt16Array>(l, i, path, |a| a.values())
}

/// Reads `<dir>/library.parquet` and `<dir>/alignment.parquet`. Every error
/// names the file.
pub fn read_library(dir: &Path) -> Result<StoredLibrary, String> {
    let (lib_path, align_path) = library_paths(dir);
    let lerr = |e: &dyn std::fmt::Display| format!("{}: {e}", lib_path.display());
    let aerr = |e: &dyn std::fmt::Display| format!("{}: {e}", align_path.display());

    let (reader, meta) = open(&lib_path, LIBRARY_KEY, "library", &library_schema())?;
    let get = |key: &str| {
        meta.iter()
            .find(|kv| kv.key == key)
            .and_then(|kv| kv.value.clone())
            .ok_or_else(|| lerr(&format!("missing metadata key {key}")))
    };
    let runs_text = get(RUNS_KEY)?;
    let runs: Vec<String> = if runs_text.is_empty() {
        Vec::new()
    } else {
        runs_text.split('\n').map(str::to_string).collect()
    };
    let max_q: f64 = get(MAX_Q_KEY)?
        .parse()
        .map_err(|e| lerr(&format!("{MAX_Q_KEY}: {e}")))?;
    let fingerprint = get(FINGERPRINT_KEY)?;
    let alignment_fp = get(ALIGNMENT_FINGERPRINT_KEY)?;

    let mut meta_rows: Vec<PrecursorMeta> = Vec::new();
    let mut proteins: Vec<String> = Vec::new();
    let mut entries: Vec<LibraryEntry> = Vec::new();
    for b in reader {
        let b = b.map_err(|e| lerr(&e))?;
        let key = column::<StringArray>(&b, 0, &lib_path)?;
        let mz = column::<Float64Array>(&b, 1, &lib_path)?;
        let charge = column::<UInt8Array>(&b, 2, &lib_path)?;
        let env = column::<ListArray>(&b, 3, &lib_path)?;
        let prot = column::<StringArray>(&b, 4, &lib_path)?;
        let q = column::<Float64Array>(&b, 5, &lib_path)?;
        let cv = column::<Float32Array>(&b, 6, &lib_path)?;
        let donors = column::<ListArray>(&b, 7, &lib_path)?;
        let best_run = column::<ListArray>(&b, 8, &lib_path)?;
        let best_rt = column::<ListArray>(&b, 9, &lib_path)?;
        for i in 0..b.num_rows() {
            let id = meta_rows.len();
            let e: Vec<f32> =
                items::<f32, Float32Array>(env, i, &lib_path, |a: &Float32Array| a.values())?;
            let envelope: [f32; ENVELOPE_LEN] = e.try_into().map_err(|e: Vec<f32>| {
                lerr(&format!(
                    "precursor {id}: envelope of {} isotopes, expected {ENVELOPE_LEN}",
                    e.len()
                ))
            })?;
            meta_rows.push(PrecursorMeta {
                key: key.value(i).to_string(),
                mz: mz.value(i),
                charge: charge.value(i),
                envelope,
            });
            proteins.push(prot.value(i).to_string());
            if q.is_null(i) {
                continue;
            }
            let runs_b = u16s(best_run, i, &lib_path)?;
            let rts_b = f64s(best_rt, i, &lib_path)?;
            if runs_b.len() != rts_b.len() {
                return Err(lerr(&format!(
                    "precursor {id}: {} best-RT runs for {} RTs",
                    runs_b.len(),
                    rts_b.len()
                )));
            }
            entries.push(LibraryEntry {
                precursor: u32::try_from(id)
                    .map_err(|_| lerr(&"more precursors than the id range"))?,
                q: q.value(i),
                consensus_rt: f64::NAN,
                faims_cv: (!cv.is_null(i)).then(|| cv.value(i)),
                donors: u16s(donors, i, &lib_path)?,
                best_rt: runs_b.into_iter().zip(rts_b).collect(),
            });
        }
    }
    let n = meta_rows.len();

    let (reader, _) = open(&align_path, ALIGNMENT_KEY, "alignment", &alignment_schema())?;
    let mut run_states: Vec<RunAlignment> = Vec::new();
    for b in reader {
        let b = b.map_err(|e| aerr(&e))?;
        let run = column::<UInt16Array>(&b, 0, &align_path)?;
        let stem = column::<StringArray>(&b, 1, &align_path)?;
        let aligned = column::<BooleanArray>(&b, 2, &align_path)?;
        let kx = column::<ListArray>(&b, 3, &align_path)?;
        let ky = column::<ListArray>(&b, 4, &align_path)?;
        let slo = column::<Float64Array>(&b, 5, &align_path)?;
        let shi = column::<Float64Array>(&b, 6, &align_path)?;
        let sd_t = column::<ListArray>(&b, 7, &align_path)?;
        let sd = column::<ListArray>(&b, 8, &align_path)?;
        let anchors = column::<UInt64Array>(&b, 9, &align_path)?;
        let coverage = column::<Float64Array>(&b, 10, &align_path)?;
        let obs_p = column::<ListArray>(&b, 11, &align_path)?;
        let obs_rt = column::<ListArray>(&b, 12, &align_path)?;
        for i in 0..b.num_rows() {
            let r = run_states.len();
            if run.value(i) as usize != r {
                return Err(aerr(&format!("row {r} holds run {}", run.value(i))));
            }
            if runs.get(r).map(String::as_str) != Some(stem.value(i)) {
                return Err(aerr(&format!(
                    "run {r} is {:?}, the library lists {:?}",
                    stem.value(i),
                    runs.get(r)
                )));
            }
            let (map, residual_sd) = if aligned.value(i) {
                if slo.is_null(i) || shi.is_null(i) {
                    return Err(aerr(&format!("run {r}: aligned without end slopes")));
                }
                let map = RunMap::from_parts(
                    f64s(kx, i, &align_path)?,
                    f64s(ky, i, &align_path)?,
                    slo.value(i),
                    shi.value(i),
                )
                .map_err(|e| aerr(&format!("run {r}: {e}")))?;
                let t = f64s(sd_t, i, &align_path)?;
                let s = f64s(sd, i, &align_path)?;
                if t.len() != s.len() {
                    return Err(aerr(&format!(
                        "run {r}: {} residual knots for {} widths",
                        t.len(),
                        s.len()
                    )));
                }
                (Some(map), Some(t.into_iter().zip(s).collect()))
            } else {
                (None, None)
            };
            let ps =
                items::<u32, UInt32Array>(obs_p, i, &align_path, |a: &UInt32Array| a.values())?;
            let ts = f64s(obs_rt, i, &align_path)?;
            if ps.len() != ts.len() {
                return Err(aerr(&format!(
                    "run {r}: {} observed precursors for {} RTs",
                    ps.len(),
                    ts.len()
                )));
            }
            run_states.push(RunAlignment {
                map,
                residual_sd,
                anchors: anchors.value(i) as usize,
                coverage: coverage.value(i),
                observations: ps.into_iter().zip(ts).collect(),
            });
        }
    }
    if run_states.len() != runs.len() {
        return Err(aerr(&format!(
            "{} runs, the library lists {}",
            run_states.len(),
            runs.len()
        )));
    }
    let alignment = Alignment::from_runs(run_states, n).map_err(|e| aerr(&e))?;
    if alignment_fingerprint(&alignment) != alignment_fp {
        return Err(aerr(&format!(
            "does not belong to {} (alignment fingerprint differs); rebuild the library",
            lib_path.display()
        )));
    }
    for e in &mut entries {
        if let Some(&(r, _)) = e.best_rt.iter().find(|&&(r, _)| r as usize >= runs.len()) {
            return Err(lerr(&format!(
                "precursor {}: best RT in run {r} outside the {} runs",
                e.precursor,
                runs.len()
            )));
        }
        e.consensus_rt = alignment
            .consensus
            .get(e.precursor as usize)
            .copied()
            .unwrap_or(f64::NAN);
    }
    entries.sort_by(|a, b| {
        meta_rows[a.precursor as usize]
            .mz
            .total_cmp(&meta_rows[b.precursor as usize].mz)
            .then(a.precursor.cmp(&b.precursor))
    });
    let stored = StoredLibrary {
        runs,
        max_q,
        library: Library {
            entries,
            alignment,
            meta: meta_rows,
        },
        proteins,
    };
    if library_fingerprint(&stored) != fingerprint {
        return Err(lerr(
            &"content does not match its fingerprint; rebuild the library",
        ));
    }
    Ok(stored)
}

#[cfg(test)]
mod tests {
    use super::*;
    use quant::align::AlignParams;
    use quant::library::{build_library, PsmIn};

    /// 300 precursors over three warped runs, 20 decoys, a FAIMS voltage on
    /// some, precursor 0 seen in run 0 only.
    fn stored() -> StoredLibrary {
        let n = 320usize;
        let meta: Vec<PrecursorMeta> = (0..n)
            .map(|i| PrecursorMeta {
                key: format!("PEP{i}/2"),
                mz: 400.0 + ((i * 37) % n) as f64 * 1.5,
                charge: 2,
                envelope: [0.55, 0.3, 0.1, 0.05],
            })
            .collect();
        let mut psms = Vec::new();
        for p in 0..n as u32 {
            let decoy = p >= 300;
            for r in 0..3u16 {
                if p == 0 && r > 0 {
                    continue;
                }
                let t = 100.0 + 5.0 * p as f64;
                let rt = match r {
                    0 => t,
                    1 => 1.05 * t + 20.0,
                    _ => t + 15.0 * (t / 600.0).sin(),
                } + ((p * 7 + r as u32) % 5) as f64 * 0.3;
                psms.push(PsmIn {
                    run: r,
                    precursor: p,
                    rt,
                    score: if decoy { 1.0 } else { 100.0 + p as f64 },
                    is_decoy: decoy,
                    faims_cv: (p % 3 == 0).then_some(-45.0),
                });
            }
        }
        let library = build_library(&psms, meta, 3, 0.01, &AlignParams::default());
        StoredLibrary {
            runs: vec!["a".into(), "b".into(), "c".into()],
            max_q: 0.01,
            library,
            proteins: (0..n).map(|i| format!("P{i};Q{i}")).collect(),
        }
    }

    fn bits(v: Option<(f64, f64)>) -> Option<(u64, u64)> {
        v.map(|(a, b)| (a.to_bits(), b.to_bits()))
    }

    #[test]
    fn library_round_trip_is_exact() {
        let dir = tempfile::tempdir().unwrap();
        let s = stored();
        assert!(s.library.alignment.maps.iter().all(Option::is_some));
        write_library(dir.path(), &s).unwrap();
        let back = read_library(dir.path()).unwrap();
        assert_eq!(back.runs, s.runs);
        assert_eq!(back.max_q, s.max_q);
        assert_eq!(back.proteins, s.proteins);
        assert_eq!(back.library.meta, s.library.meta);
        assert_eq!(back.library.entries.len(), s.library.entries.len());
        for (a, b) in s.library.entries.iter().zip(&back.library.entries) {
            assert_eq!(a.precursor, b.precursor);
            assert_eq!(a.q.to_bits(), b.q.to_bits());
            assert_eq!(a.consensus_rt.to_bits(), b.consensus_rt.to_bits());
            assert_eq!(
                (&a.faims_cv, &a.donors, &a.best_rt),
                (&b.faims_cv, &b.donors, &b.best_rt)
            );
        }
        let (a, b) = (&s.library.alignment, &back.library.alignment);
        assert_eq!(a.runs(), b.runs());
        for p in 0..s.library.meta.len() as u32 {
            for r in 0..3u16 {
                assert_eq!(bits(a.predict(p, r)), bits(b.predict(p, r)));
                assert_eq!(
                    bits(a.predict_excluding(p, r)),
                    bits(b.predict_excluding(p, r))
                );
            }
        }
        assert!(!dir.path().join("library.parquet.partial").exists());
    }

    #[test]
    fn library_rejects_unknown_major_version_naming_the_file() {
        let dir = tempfile::tempdir().unwrap();
        write_with_schema_version(dir.path(), &stored(), "2.0").unwrap();
        let e = read_library(dir.path()).unwrap_err();
        assert!(
            e.contains("alignment.parquet") || e.contains("library.parquet"),
            "{e}"
        );
        assert!(e.contains("2.0"), "{e}");
    }

    #[test]
    fn fingerprints_are_stable_and_detect_a_foreign_alignment() {
        let s = stored();
        assert_eq!(library_fingerprint(&s), library_fingerprint(&s.clone()));
        let mut other = s.clone();
        other.max_q = 0.05;
        assert_ne!(library_fingerprint(&s), library_fingerprint(&other));

        // An alignment.parquet from another library is rejected by name.
        let a = tempfile::tempdir().unwrap();
        let b = tempfile::tempdir().unwrap();
        write_library(a.path(), &s).unwrap();
        let mut moved = s.clone();
        let mut runs = moved.library.alignment.runs();
        runs[0].coverage += 0.01;
        moved.library.alignment =
            quant::align::Alignment::from_runs(runs, moved.library.meta.len()).unwrap();
        write_library(b.path(), &moved).unwrap();
        std::fs::copy(b.path().join(ALIGNMENT_FILE), a.path().join(ALIGNMENT_FILE)).unwrap();
        let e = read_library(a.path()).unwrap_err();
        assert!(
            e.contains("alignment.parquet") && e.contains("fingerprint"),
            "{e}"
        );
        assert_eq!(
            library_fingerprint(&read_library(b.path()).unwrap()),
            library_fingerprint(&moved)
        );
    }

    #[test]
    fn missing_library_names_the_file() {
        let dir = tempfile::tempdir().unwrap();
        let e = read_library(dir.path()).unwrap_err();
        assert!(e.contains("library.parquet"), "{e}");
    }
}
