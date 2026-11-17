use std::collections::BTreeMap;
use std::path::Path;

use arrow::array::{Array, BooleanArray, StringArray};
use arrow::compute::filter_record_batch;
use arrow::record_batch::RecordBatch;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::ArrowWriter;
use parquet::file::properties::WriterProperties;

use crate::TrainError;

pub fn split_store_by_protocol(
    src: &Path,
    dir: &Path,
) -> Result<Vec<(String, std::path::PathBuf)>, TrainError> {
    let model_protocol = read_model_protocols(src)?;

    let mut by_protocol: BTreeMap<String, Vec<RecordBatch>> = BTreeMap::new();

    let file = std::fs::File::open(src)?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file)
        .map_err(|e| TrainError::Parquet(e.to_string()))?;
    let schema = builder.schema().clone();
    let reader = builder
        .build()
        .map_err(|e| TrainError::Parquet(e.to_string()))?;

    for batch in reader {
        let batch = batch.map_err(|e| TrainError::Parquet(e.to_string()))?;
        let mid_col = batch
            .column_by_name("model_id")
            .ok_or_else(|| TrainError::Other("missing model_id column".into()))?
            .as_any()
            .downcast_ref::<StringArray>()
            .ok_or_else(|| TrainError::Other("model_id not StringArray".into()))?;

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
                by_protocol
                    .entry(proto.to_string())
                    .or_default()
                    .push(filtered);
            }
        }
    }

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

        let props = WriterProperties::builder()
            .set_compression(parquet::basic::Compression::ZSTD(
                parquet::basic::ZstdLevel::default(),
            ))
            .set_dictionary_enabled(false)
            .build();
        let merged = arrow::compute::concat_batches(&schema, batches.iter())
            .map_err(|e| TrainError::Other(format!("concat_batches: {e}")))?;

        let out = std::fs::File::create(&part_path)?;
        let mut writer = ArrowWriter::try_new(out, schema.clone(), Some(props))
            .map_err(|e| TrainError::Parquet(e.to_string()))?;
        writer
            .write(&merged)
            .map_err(|e| TrainError::Parquet(e.to_string()))?;
        writer
            .close()
            .map_err(|e| TrainError::Parquet(e.to_string()))?;

        written.push((proto.clone(), part_path));
    }

    Ok(written)
}

fn read_model_protocols(
    src: &Path,
) -> Result<std::collections::HashMap<String, String>, TrainError> {
    let file = std::fs::File::open(src)?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file)
        .map_err(|e| TrainError::Parquet(e.to_string()))?;
    let reader = builder
        .build()
        .map_err(|e| TrainError::Parquet(e.to_string()))?;

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

fn str_col<'a>(batch: &'a RecordBatch, name: &str) -> Result<&'a StringArray, TrainError> {
    batch
        .column_by_name(name)
        .ok_or_else(|| TrainError::Other(format!("missing column {name}")))?
        .as_any()
        .downcast_ref::<StringArray>()
        .ok_or_else(|| TrainError::Other(format!("column {name} not StringArray")))
}
