//! Split a single-file model store into per-protocol Parquet partitions.
//!
//! The bundled store ships as one `resources/models.parquet` file. This module
//! rewrites it as a Hive-style partitioned directory:
//!
//! ```text
//! <dir>/protocol=Automatic/models.parquet
//! <dir>/protocol=TMT/models.parquet
//! <dir>/protocol=Phosphorylation/models.parquet
//! ...
//! ```
//!
//! Each model belongs to exactly one protocol (its manifest `protocol`
//! column). Every row that references a model — manifest, table, source and
//! stat rows — is routed into that model's protocol partition. The split is
//! **lossless**: rows are copied verbatim (via the Arrow `filter` kernel) so
//! the union of all partitions is row-for-row identical to the input file,
//! and [`crate::store::ModelStore::open`] on the directory yields the same
//! in-memory store as opening the single file.
//!
//! # On-disk size note — this is the DEV/PARITY splitter, not the shipping one
//!
//! The per-model GBDT blob columns (`rich_ion_model_bytes`,
//! `frag_intensity_model_bytes`) are high-entropy and sparse (one blob per
//! manifest row, NULL on every other row). The parquet-rs 53 writer encodes
//! these columns several times larger than parquet-cpp (pyarrow) — even with
//! ZSTD and dictionary disabled it produces a partition set of ~51 MB versus
//! the ~13.6 MB single-file source. This is a writer-efficiency gap, not a data
//! difference; the reader is codec/encoding-agnostic, so a pyarrow-written
//! partition set loads byte-for-byte identically (proven by the
//! `partition_parity` tests).
//!
//! Therefore the **shipped** `resources/models/protocol=<P>/models.parquet`
//! files are written by `scripts/split_store_by_protocol.py` (pyarrow/zstd),
//! which produces ~8.7 MB total (e.g. TMT ~2.1 MB, Automatic ~5.7 MB) — smaller
//! than the single file AND per-protocol granular. This Rust function remains
//! the dev/parity tool (it backs the parity test gate); do NOT use its output
//! as the shipping bundle.

use std::collections::BTreeMap;
use std::path::Path;

use arrow::array::{Array, BooleanArray, StringArray};
use arrow::compute::filter_record_batch;
use arrow::record_batch::RecordBatch;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::ArrowWriter;
use parquet::file::properties::WriterProperties;

use crate::TrainError;

/// Read the single-file store at `src` and write it as a per-protocol
/// partitioned directory rooted at `dir` (created if absent).
///
/// Returns the list of `(protocol, partition_path)` written, sorted by
/// protocol name.
pub fn split_store_by_protocol(
    src: &Path,
    dir: &Path,
) -> Result<Vec<(String, std::path::PathBuf)>, TrainError> {
    // 1. First pass: map every model_id → its protocol (from manifest rows).
    let model_protocol = read_model_protocols(src)?;

    // 2. Group batches by protocol. We accumulate, per protocol, the filtered
    //    record batches so each partition can be written in one shot.
    let mut by_protocol: BTreeMap<String, Vec<RecordBatch>> = BTreeMap::new();

    let file = std::fs::File::open(src)?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file)
        .map_err(|e| TrainError::Parquet(e.to_string()))?;
    let schema = builder.schema().clone();
    let reader = builder.build().map_err(|e| TrainError::Parquet(e.to_string()))?;

    for batch in reader {
        let batch = batch.map_err(|e| TrainError::Parquet(e.to_string()))?;
        let mid_col = batch
            .column_by_name("model_id")
            .ok_or_else(|| TrainError::Other("missing model_id column".into()))?
            .as_any()
            .downcast_ref::<StringArray>()
            .ok_or_else(|| TrainError::Other("model_id not StringArray".into()))?;

        // For each protocol present in this batch, build a boolean mask of the
        // rows whose model_id belongs to that protocol, then filter.
        // Determine the set of protocols touched by this batch first.
        let mut row_protocol: Vec<Option<&str>> = Vec::with_capacity(batch.num_rows());
        for i in 0..batch.num_rows() {
            if mid_col.is_null(i) {
                row_protocol.push(None);
                continue;
            }
            let mid = mid_col.value(i);
            match model_protocol.get(mid) {
                Some(p) => row_protocol.push(Some(p.as_str())),
                None => {
                    return Err(TrainError::Other(format!(
                        "model_id '{mid}' has no manifest row (cannot assign a protocol partition)"
                    )));
                }
            }
        }

        let protocols: std::collections::BTreeSet<&str> =
            row_protocol.iter().filter_map(|p| *p).collect();

        for proto in protocols {
            let mask: BooleanArray = row_protocol
                .iter()
                .map(|p| Some(*p == Some(proto)))
                .collect();
            let filtered = filter_record_batch(&batch, &mask)
                .map_err(|e| TrainError::Other(format!("filter_record_batch: {e}")))?;
            if filtered.num_rows() > 0 {
                by_protocol.entry(proto.to_string()).or_default().push(filtered);
            }
        }
    }

    // 3. Write one partition file per protocol.
    // Prune any pre-existing `protocol=*` partitions first so a rerun against a
    // smaller src can't leave stale partitions behind (mirrors the python
    // splitter's rmtree-then-recreate).
    if dir.exists() {
        for entry in std::fs::read_dir(dir)? {
            let entry = entry?;
            let p = entry.path();
            let name = p.file_name().and_then(|n| n.to_str()).unwrap_or("");
            if p.is_dir() && name.starts_with("protocol=") {
                std::fs::remove_dir_all(&p)?;
            }
        }
    }
    std::fs::create_dir_all(dir)?;
    let mut written: Vec<(String, std::path::PathBuf)> = Vec::new();
    for (proto, batches) in &by_protocol {
        let part_dir = dir.join(format!("protocol={proto}"));
        std::fs::create_dir_all(&part_dir)?;
        let part_path = part_dir.join("models.parquet");

        // Storage tuning (the reader is codec/encoding-agnostic, so this only
        // affects on-disk size):
        //  * ZSTD compresses the large per-model GBDT blob columns
        //    (rich_ion / frag_intensity model bytes) better than SNAPPY.
        //  * Dictionary encoding is DISABLED: every model's GBDT blob is
        //    distinct, so RLE_DICTIONARY only adds an index page on top of the
        //    plain values.
        //
        // NOTE: parquet-rs encodes the high-entropy/sparse `Binary` blob
        // columns several times larger than parquet-cpp, so this dev tool's
        // output is NOT the shipping bundle — the shipped partitions are
        // written by `scripts/split_store_by_protocol.py` (pyarrow). See the
        // module docs.
        let props = WriterProperties::builder()
            .set_compression(parquet::basic::Compression::ZSTD(
                parquet::basic::ZstdLevel::default(),
            ))
            .set_dictionary_enabled(false)
            .build();
        // Concatenate this protocol's filtered batches into one contiguous row
        // group (matches the single-file layout).
        let merged = arrow::compute::concat_batches(&schema, batches.iter())
            .map_err(|e| TrainError::Other(format!("concat_batches: {e}")))?;

        let out = std::fs::File::create(&part_path)?;
        let mut writer = ArrowWriter::try_new(out, schema.clone(), Some(props))
            .map_err(|e| TrainError::Parquet(e.to_string()))?;
        writer.write(&merged).map_err(|e| TrainError::Parquet(e.to_string()))?;
        writer.close().map_err(|e| TrainError::Parquet(e.to_string()))?;

        written.push((proto.clone(), part_path));
    }

    Ok(written)
}

/// First pass: build `model_id → protocol` from the manifest rows of `src`.
fn read_model_protocols(
    src: &Path,
) -> Result<std::collections::HashMap<String, String>, TrainError> {
    let file = std::fs::File::open(src)?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file)
        .map_err(|e| TrainError::Parquet(e.to_string()))?;
    let reader = builder.build().map_err(|e| TrainError::Parquet(e.to_string()))?;

    let mut map = std::collections::HashMap::new();
    for batch in reader {
        let batch = batch.map_err(|e| TrainError::Parquet(e.to_string()))?;
        let rk = str_col(&batch, "record_kind")?;
        let mid = str_col(&batch, "model_id")?;
        let proto = str_col(&batch, "protocol")?;
        for i in 0..batch.num_rows() {
            if !rk.is_null(i) && rk.value(i) == "manifest" {
                let model_id = mid.value(i).to_string();
                let protocol = if proto.is_null(i) {
                    String::new()
                } else {
                    proto.value(i).to_string()
                };
                map.insert(model_id, protocol);
            }
        }
    }
    Ok(map)
}

fn str_col<'a>(
    batch: &'a RecordBatch,
    name: &str,
) -> Result<&'a StringArray, TrainError> {
    batch
        .column_by_name(name)
        .ok_or_else(|| TrainError::Other(format!("missing column {name}")))?
        .as_any()
        .downcast_ref::<StringArray>()
        .ok_or_else(|| TrainError::Other(format!("column {name} not StringArray")))
}
