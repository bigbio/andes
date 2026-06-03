//! Parquet writer: serialize one or more [`Param`] models to a single file.
//!
//! Call [`write_models`] with a slice of `(model_id, &Param)` pairs.
//! All models are written to **one** Parquet file using `combined_schema()`.
//! Manifest rows (scalar fields + charge_hist) come first, then table rows
//! (rank_dist, err_dist, frag_off, precursor_off, partition list).
//!
//! Use [`write_model_with_sources`] to also persist per-source sufficient
//! statistics ([`crate::counts::CountStats`]) and a sources ledger.

use std::path::Path;
use std::sync::Arc;

use arrow::array::{
    ArrayRef, BooleanBuilder, Float32Builder, Int32Builder, Int64Builder,
    ListBuilder, StringArray, StringBuilder,
    StructBuilder,
};
use arrow::datatypes::{DataType, Field, Fields, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use parquet::file::properties::WriterProperties;

use scoring_crate::param_model::{IonType, Param, Partition};

use crate::counts::CountStats;
use crate::store::schema::combined_schema;
use crate::TrainError;

// ── public types ─────────────────────────────────────────────────────────────

/// Metadata for a single training source stored in the Parquet ledger.
///
/// One [`SourceLedger`] row is written per `(model_id, source_id)` pair and
/// accompanies the matching [`CountStats`] sufficient-statistics rows.
#[derive(Debug, Clone, PartialEq)]
pub struct SourceLedger {
    /// Unique identifier for this source within the model (e.g. a dataset accession).
    pub source_id: String,
    /// Human-readable dataset name (e.g. `"PXD001819"`).
    pub dataset: String,
    /// Number of confident PSMs contributed by this source.
    pub n_psms: i64,
    /// ISO 8601 date string (e.g. `"2026-01-01"`).
    pub date: String,
    /// Relative weight applied when mixing this source into training.
    pub weight: f32,
    /// FDR threshold used to select confident PSMs from this source.
    pub train_fdr: f32,
    /// Instrument label (e.g. `"QExactive"`).
    pub instrument: String,
    /// Experiment class slug(s) (e.g. `"standard"`, `"tmt"`).
    pub experiment_class: String,
}

// ── type aliases ─────────────────────────────────────────────────────────────

/// One encoded precursor-offset entry: (reduced_charge, offset, tol_is_ppm, tol_val, frequency).
type PrecursorEntry = (i32, f32, bool, f32, f32);

// ── helpers ──────────────────────────────────────────────────────────────────

fn f32_to_i32_bits(v: f32) -> i32 {
    v.to_bits() as i32
}

// ── public entry point ───────────────────────────────────────────────────────

/// Write `models` to a Parquet file at `path`.
///
/// Rows are sorted by `model_id` within each record kind (manifest first,
/// then table rows).
pub fn write_models(path: &Path, models: &[(String, &Param)]) -> Result<(), TrainError> {
    // Delegate to write_model_with_sources with no source data.
    if models.is_empty() {
        let schema = combined_schema();
        let props = WriterProperties::builder().build();
        let file = std::fs::File::create(path)?;
        let writer = ArrowWriter::try_new(file, schema, Some(props))
            .map_err(|e| TrainError::Parquet(e.to_string()))?;
        writer.close().map_err(|e| TrainError::Parquet(e.to_string()))?;
        return Ok(());
    }

    let schema = combined_schema();

    // Sort by model_id for deterministic output.
    let mut sorted: Vec<(&str, &Param)> = models.iter().map(|(id, p)| (id.as_str(), *p)).collect();
    sorted.sort_by_key(|(id, _)| *id);

    // Build manifest rows + table rows.
    let manifest_batch = build_manifest_batch(&schema, &sorted)?;
    let table_batch = build_table_batch(&schema, &sorted)?;

    let props = WriterProperties::builder().build();
    let file = std::fs::File::create(path)?;
    let mut writer = ArrowWriter::try_new(file, schema.clone(), Some(props))
        .map_err(|e| TrainError::Parquet(e.to_string()))?;

    if manifest_batch.num_rows() > 0 {
        writer
            .write(&manifest_batch)
            .map_err(|e| TrainError::Parquet(e.to_string()))?;
    }
    if table_batch.num_rows() > 0 {
        writer
            .write(&table_batch)
            .map_err(|e| TrainError::Parquet(e.to_string()))?;
    }

    writer.close().map_err(|e| TrainError::Parquet(e.to_string()))?;
    Ok(())
}

/// Write a single model with per-source sufficient statistics to `path`.
///
/// Writes manifest + table rows for `param` (same as [`write_models`] for one
/// model) plus `"source"` (ledger) and `"stat"` (CountStats) rows for each
/// entry in `sources`.  `sources` may be empty.
pub fn write_model_with_sources(
    path: &Path,
    model_id: &str,
    param: &Param,
    sources: &[(SourceLedger, CountStats)],
) -> Result<(), TrainError> {
    let schema = combined_schema();
    let sorted: Vec<(&str, &Param)> = vec![(model_id, param)];

    let manifest_batch = build_manifest_batch(&schema, &sorted)?;
    let table_batch = build_table_batch(&schema, &sorted)?;

    let props = WriterProperties::builder().build();
    let file = std::fs::File::create(path)?;
    let mut writer = ArrowWriter::try_new(file, schema.clone(), Some(props))
        .map_err(|e| TrainError::Parquet(e.to_string()))?;

    if manifest_batch.num_rows() > 0 {
        writer.write(&manifest_batch).map_err(|e| TrainError::Parquet(e.to_string()))?;
    }
    if table_batch.num_rows() > 0 {
        writer.write(&table_batch).map_err(|e| TrainError::Parquet(e.to_string()))?;
    }

    if !sources.is_empty() {
        let source_batch = build_source_batch(&schema, model_id, sources)?;
        if source_batch.num_rows() > 0 {
            writer.write(&source_batch).map_err(|e| TrainError::Parquet(e.to_string()))?;
        }
        let stat_batch = build_stat_batch(&schema, model_id, sources)?;
        if stat_batch.num_rows() > 0 {
            writer.write(&stat_batch).map_err(|e| TrainError::Parquet(e.to_string()))?;
        }
    }

    writer.close().map_err(|e| TrainError::Parquet(e.to_string()))?;
    Ok(())
}

// ── manifest batch ───────────────────────────────────────────────────────────

fn build_manifest_batch(
    schema: &Arc<Schema>,
    models: &[(&str, &Param)],
) -> Result<RecordBatch, TrainError> {
    let n = models.len();

    // Shared columns (non-null).
    let record_kind: ArrayRef =
        Arc::new(StringArray::from(vec!["manifest"; n]));
    let model_id: ArrayRef =
        Arc::new(StringArray::from(models.iter().map(|(id, _)| *id).collect::<Vec<_>>()));

    // Nullable manifest scalars.
    let mut activation = StringBuilder::new();
    let mut instrument = StringBuilder::new();
    let mut enzyme = StringBuilder::new();
    let mut protocol = StringBuilder::new();
    let mut version = Int32Builder::new();
    let mut mme_val = Float32Builder::new();
    let mut mme_is_ppm = BooleanBuilder::new();
    let mut apply_deconv = BooleanBuilder::new();
    let mut deconv_tol = Float32Builder::new();
    let mut num_segments = Int32Builder::new();
    let mut max_rank = Int32Builder::new();
    let mut error_scaling_factor = Int32Builder::new();
    let mut min_charge = Int32Builder::new();
    let mut max_charge = Int32Builder::new();
    let mut num_precursor_off = Int32Builder::new();

    // charge_hist: List<Struct<charge: Int32, count: Int32>>
    // We use with_field() so the inner item is non-null (matching combined_schema's
    // list_of_struct definition which uses false for the struct item nullability).
    let charge_struct_fields: Fields = vec![
        Field::new("charge", DataType::Int32, false),
        Field::new("count", DataType::Int32, false),
    ]
    .into();
    let charge_item_field = Arc::new(Field::new(
        "item",
        DataType::Struct(charge_struct_fields.clone()),
        false, // non-null items, matching schema
    ));
    let charge_hist_builder = StructBuilder::new(
        charge_struct_fields.clone(),
        vec![
            Box::new(Int32Builder::new()) as Box<dyn arrow::array::ArrayBuilder>,
            Box::new(Int32Builder::new()) as Box<dyn arrow::array::ArrayBuilder>,
        ],
    );
    let mut charge_hist = ListBuilder::new(charge_hist_builder)
        .with_field(charge_item_field);

    for (_, param) in models {
        activation.append_value(param.data_type.activation.name());
        instrument.append_value(param.data_type.instrument.name());
        match param.data_type.enzyme {
            Some(e) => enzyme.append_value(e.name()),
            None => enzyme.append_null(),
        }
        protocol.append_value(param.data_type.protocol.name());
        version.append_value(param.version);
        let (is_ppm, raw) = match param.mme {
            model::tolerance::Tolerance::Ppm(v) => (true, v as f32),
            model::tolerance::Tolerance::Da(v) => (false, v as f32),
        };
        mme_val.append_value(raw);
        mme_is_ppm.append_value(is_ppm);
        apply_deconv.append_value(param.apply_deconvolution);
        deconv_tol.append_value(param.deconvolution_error_tolerance);
        num_segments.append_value(param.num_segments);
        max_rank.append_value(param.max_rank);
        error_scaling_factor.append_value(param.error_scaling_factor);
        min_charge.append_value(param.min_charge);
        max_charge.append_value(param.max_charge);
        num_precursor_off.append_value(param.num_precursor_off);

        // charge_hist entries
        let sb = charge_hist.values();
        for &(charge, count) in &param.charge_hist {
            sb.field_builder::<Int32Builder>(0).unwrap().append_value(charge);
            sb.field_builder::<Int32Builder>(1).unwrap().append_value(count);
            sb.append(true);
        }
        charge_hist.append(true);
    }

    // Table-only and source/stat-only columns: null for manifest rows.
    let null_i32 = null_i32_array(n);
    let null_i64 = null_i64_array(n);
    let null_utf8 = null_utf8_array(n);
    let null_f32 = null_f32_array(n);
    let null_float_list = null_float_list_array(n);
    let null_int64_list = null_int64_list_array(n);
    let null_int32_list = null_int32_list_array(n);
    let null_precursor_list = null_struct_list_array(
        n,
        vec![
            Field::new("reduced_charge", DataType::Int32, false),
            Field::new("offset", DataType::Float32, false),
            Field::new("tol_is_ppm", DataType::Boolean, false),
            Field::new("tol_val", DataType::Float32, false),
            Field::new("frequency", DataType::Float32, false),
        ],
    );

    let columns: Vec<ArrayRef> = vec![
        record_kind,
        model_id,
        Arc::new(activation.finish()),
        Arc::new(instrument.finish()),
        Arc::new(enzyme.finish()),
        Arc::new(protocol.finish()),
        Arc::new(version.finish()),
        Arc::new(mme_val.finish()),
        Arc::new(mme_is_ppm.finish()),
        Arc::new(apply_deconv.finish()),
        Arc::new(deconv_tol.finish()),
        Arc::new(num_segments.finish()),
        Arc::new(max_rank.finish()),
        Arc::new(error_scaling_factor.finish()),
        Arc::new(min_charge.finish()),
        Arc::new(max_charge.finish()),
        Arc::new(num_precursor_off.finish()),
        Arc::new(charge_hist.finish()),
        // table-only → null
        null_i32.clone(),     // part_charge
        null_i32.clone(),     // part_mass_bits
        null_i32.clone(),     // part_seg
        null_utf8.clone(),    // ion_kind
        null_i32.clone(),     // ion_charge
        null_i32.clone(),     // ion_offset_bits
        null_utf8.clone(),    // table_kind
        null_float_list,      // values
        null_precursor_list,  // precursor_offsets
        // source/stat-only → null
        null_utf8.clone(),    // source_id
        null_utf8.clone(),    // dataset
        null_i64,             // n_psms
        null_utf8.clone(),    // date
        null_f32.clone(),     // weight
        null_f32.clone(),     // train_fdr
        null_utf8.clone(),    // src_instrument
        null_utf8.clone(),    // src_experiment_class
        null_int64_list,      // counts
        null_int32_list,      // charge_keys
    ];

    RecordBatch::try_new(schema.clone(), columns).map_err(|e| TrainError::Parquet(e.to_string()))
}

// ── table batch ──────────────────────────────────────────────────────────────

fn build_table_batch(
    schema: &Arc<Schema>,
    models: &[(&str, &Param)],
) -> Result<RecordBatch, TrainError> {
    // Columns parallel to combined_schema table section.
    let mut record_kinds: Vec<&str> = Vec::new();
    let mut model_ids: Vec<&str> = Vec::new();
    let mut part_charges: Vec<Option<i32>> = Vec::new();
    let mut part_mass_bits: Vec<Option<i32>> = Vec::new();
    let mut part_segs: Vec<Option<i32>> = Vec::new();
    let mut ion_kinds: Vec<Option<&str>> = Vec::new();
    let mut ion_charges: Vec<Option<i32>> = Vec::new();
    let mut ion_offset_bits: Vec<Option<i32>> = Vec::new();
    let mut table_kinds: Vec<Option<&str>> = Vec::new();

    // Per-row payloads.
    let mut all_values: Vec<Option<Vec<f32>>> = Vec::new();
    let mut all_precursor_offsets: Vec<Option<Vec<PrecursorEntry>>> = Vec::new();

    for &(model_id, param) in models {
        // ── partition list (record_kind="partition") ──────────────────────────
        // Store the sorted partition list so the reader knows the exact order
        // used for rank_dist / frag_off assignment. We store one row per partition
        // with table_kind="partition" and no payload.
        for part in &param.partitions {
            record_kinds.push("table");
            model_ids.push(model_id);
            part_charges.push(Some(part.charge));
            part_mass_bits.push(Some(f32_to_i32_bits(part.parent_mass)));
            part_segs.push(Some(part.seg_num));
            ion_kinds.push(Some("-"));
            ion_charges.push(Some(0));
            ion_offset_bits.push(Some(0));
            table_kinds.push(Some("partition"));
            all_values.push(None);
            all_precursor_offsets.push(None);
        }

        // ── precursor_off_map ─────────────────────────────────────────────────
        // One row per charge key. Sort by charge for determinism.
        let mut precursor_charges: Vec<i32> = param.precursor_off_map.keys().copied().collect();
        precursor_charges.sort_unstable();
        for charge in precursor_charges {
            let entries = &param.precursor_off_map[&charge];
            let payload: Vec<(i32, f32, bool, f32, f32)> = entries
                .iter()
                .map(|e| {
                    let (is_ppm, tol_val) = match e.tolerance {
                        model::tolerance::Tolerance::Ppm(v) => (true, v as f32),
                        model::tolerance::Tolerance::Da(v) => (false, v as f32),
                    };
                    (e.reduced_charge, e.offset, is_ppm, tol_val, e.frequency)
                })
                .collect();
            record_kinds.push("table");
            model_ids.push(model_id);
            part_charges.push(Some(charge));
            part_mass_bits.push(Some(0));
            part_segs.push(Some(0));
            ion_kinds.push(Some("-"));
            ion_charges.push(Some(0));
            ion_offset_bits.push(Some(0));
            table_kinds.push(Some("precursor_off"));
            all_values.push(None);
            all_precursor_offsets.push(Some(payload));
        }

        // ── frag_off_table ────────────────────────────────────────────────────
        // Iterate partitions in sorted order (same order the binary reader stores them).
        for part in &param.partitions {
            if let Some(frags) = param.frag_off_table.get(part) {
                // Encode as flat f32 list: groups of 4 per entry.
                // offset_bits is already f32::to_bits() so we recover the exact f32
                // via f32::from_bits() — this round-trips without precision loss.
                let flat: Vec<f32> = frags
                    .iter()
                    .flat_map(|f| {
                        let (is_prefix, charge, off_f32) = match f.ion_type {
                            IonType::Prefix { charge, offset_bits } => (1.0f32, charge as f32, f32::from_bits(offset_bits)),
                            IonType::Suffix { charge, offset_bits } => (0.0f32, charge as f32, f32::from_bits(offset_bits)),
                            IonType::Noise => (-1.0f32, 0.0f32, 0.0f32),
                        };
                        [is_prefix, charge, off_f32, f.frequency]
                    })
                    .collect();

                record_kinds.push("table");
                model_ids.push(model_id);
                part_charges.push(Some(part.charge));
                part_mass_bits.push(Some(f32_to_i32_bits(part.parent_mass)));
                part_segs.push(Some(part.seg_num));
                ion_kinds.push(Some("-"));
                ion_charges.push(Some(0));
                ion_offset_bits.push(Some(0));
                table_kinds.push(Some("frag_off"));
                all_values.push(Some(flat));
                all_precursor_offsets.push(None);
            }
        }

        // ── rank_dist_table ───────────────────────────────────────────────────
        // One row per (partition, ion_type): values = Vec<f32> of length max_rank+1.
        // Iterate partitions in sorted order to be deterministic.
        for part in &param.partitions {
            if let Some(ion_map) = param.rank_dist_table.get(part) {
                // Iterate in a stable order: frag_off entries first (in insertion order),
                // then Noise. This mirrors the binary writer which writes ion types in the
                // order they appear in frag_off_table, then appends Noise.
                // We need to reproduce the same order on read. Store the ion_type info
                // in the row's ion_kind/ion_charge/ion_offset_bits columns.
                for (ion, freqs) in ion_map {
                    let (kind_str, ic, iob) = encode_ion_type(ion);
                    record_kinds.push("table");
                    model_ids.push(model_id);
                    part_charges.push(Some(part.charge));
                    part_mass_bits.push(Some(f32_to_i32_bits(part.parent_mass)));
                    part_segs.push(Some(part.seg_num));
                    ion_kinds.push(Some(kind_str));
                    ion_charges.push(Some(ic));
                    ion_offset_bits.push(Some(iob));
                    table_kinds.push(Some("rank_dist"));
                    all_values.push(Some(freqs.clone()));
                    all_precursor_offsets.push(None);
                }
            }
        }

        // ── error distribution tables ─────────────────────────────────────────
        for part in &param.partitions {
            if let Some(v) = param.ion_err_dist_table.get(part) {
                emit_dist_row(
                    &mut record_kinds, &mut model_ids, &mut part_charges,
                    &mut part_mass_bits, &mut part_segs, &mut ion_kinds,
                    &mut ion_charges, &mut ion_offset_bits, &mut table_kinds,
                    &mut all_values, &mut all_precursor_offsets,
                    model_id, part.charge, part.parent_mass, part.seg_num,
                    "ion_err", v,
                );
            }
            if let Some(v) = param.noise_err_dist_table.get(part) {
                emit_dist_row(
                    &mut record_kinds, &mut model_ids, &mut part_charges,
                    &mut part_mass_bits, &mut part_segs, &mut ion_kinds,
                    &mut ion_charges, &mut ion_offset_bits, &mut table_kinds,
                    &mut all_values, &mut all_precursor_offsets,
                    model_id, part.charge, part.parent_mass, part.seg_num,
                    "noise_err", v,
                );
            }
            if let Some(v) = param.ion_existence_table.get(part) {
                emit_dist_row(
                    &mut record_kinds, &mut model_ids, &mut part_charges,
                    &mut part_mass_bits, &mut part_segs, &mut ion_kinds,
                    &mut ion_charges, &mut ion_offset_bits, &mut table_kinds,
                    &mut all_values, &mut all_precursor_offsets,
                    model_id, part.charge, part.parent_mass, part.seg_num,
                    "ion_existence", v,
                );
            }
        }
    }

    let nrows = record_kinds.len();

    // ── build Arrow arrays ───────────────────────────────────────────────────
    let record_kind_arr: ArrayRef = Arc::new(StringArray::from(record_kinds));
    let model_id_arr: ArrayRef = Arc::new(StringArray::from(model_ids));
    let part_charge_arr: ArrayRef = Arc::new(arrow::array::Int32Array::from(part_charges));
    let part_mass_bits_arr: ArrayRef = Arc::new(arrow::array::Int32Array::from(part_mass_bits));
    let part_seg_arr: ArrayRef = Arc::new(arrow::array::Int32Array::from(part_segs));

    // ion_kind: Option<&str> → StringArray
    let ion_kind_arr: ArrayRef = Arc::new(StringArray::from(
        ion_kinds.into_iter().map(|s| s.map(|x| x.to_string())).collect::<Vec<_>>(),
    ));
    let ion_charge_arr: ArrayRef = Arc::new(arrow::array::Int32Array::from(ion_charges));
    let ion_offset_bits_arr: ArrayRef = Arc::new(arrow::array::Int32Array::from(ion_offset_bits));
    let table_kind_arr: ArrayRef = Arc::new(StringArray::from(
        table_kinds.into_iter().map(|s| s.map(|x| x.to_string())).collect::<Vec<_>>(),
    ));

    // values: List<Float32>
    let values_arr = build_float_list(all_values);

    // precursor_offsets: List<Struct<...>>
    let precursor_arr = build_precursor_list(all_precursor_offsets);

    // manifest-only and source/stat-only columns: null for table rows.
    let null_i32 = null_i32_array(nrows);
    let null_i64 = null_i64_array(nrows);
    let null_utf8 = null_utf8_array(nrows);
    let null_bool = null_bool_array(nrows);
    let null_f32 = null_f32_array(nrows);
    let null_int64_list = null_int64_list_array(nrows);
    let null_int32_list = null_int32_list_array(nrows);
    let null_charge_hist = null_struct_list_array(
        nrows,
        vec![
            Field::new("charge", DataType::Int32, false),
            Field::new("count", DataType::Int32, false),
        ],
    );

    let columns: Vec<ArrayRef> = vec![
        record_kind_arr,
        model_id_arr,
        // manifest-only → null
        null_utf8.clone(),   // activation
        null_utf8.clone(),   // instrument
        null_utf8.clone(),   // enzyme
        null_utf8.clone(),   // protocol
        null_i32.clone(),    // version
        null_f32.clone(),    // mme_val
        null_bool.clone(),   // mme_is_ppm
        null_bool.clone(),   // apply_deconv
        null_f32.clone(),    // deconv_tol
        null_i32.clone(),    // num_segments
        null_i32.clone(),    // max_rank
        null_i32.clone(),    // error_scaling_factor
        null_i32.clone(),    // min_charge
        null_i32.clone(),    // max_charge
        null_i32.clone(),    // num_precursor_off
        null_charge_hist,    // charge_hist
        // table columns
        part_charge_arr,
        part_mass_bits_arr,
        part_seg_arr,
        ion_kind_arr,
        ion_charge_arr,
        ion_offset_bits_arr,
        table_kind_arr,
        values_arr,
        precursor_arr,
        // source/stat-only → null
        null_utf8.clone(),   // source_id
        null_utf8.clone(),   // dataset
        null_i64,            // n_psms
        null_utf8.clone(),   // date
        null_f32.clone(),    // weight
        null_f32.clone(),    // train_fdr
        null_utf8.clone(),   // src_instrument
        null_utf8.clone(),   // src_experiment_class
        null_int64_list,     // counts
        null_int32_list,     // charge_keys
    ];

    RecordBatch::try_new(schema.clone(), columns).map_err(|e| TrainError::Parquet(e.to_string()))
}

// ── row emission helpers ─────────────────────────────────────────────────────

#[allow(clippy::too_many_arguments)]
fn emit_dist_row<'a>(
    record_kinds: &mut Vec<&'a str>,
    model_ids: &mut Vec<&'a str>,
    part_charges: &mut Vec<Option<i32>>,
    part_mass_bits: &mut Vec<Option<i32>>,
    part_segs: &mut Vec<Option<i32>>,
    ion_kinds: &mut Vec<Option<&'a str>>,
    ion_charges: &mut Vec<Option<i32>>,
    ion_offset_bits: &mut Vec<Option<i32>>,
    table_kinds: &mut Vec<Option<&'a str>>,
    all_values: &mut Vec<Option<Vec<f32>>>,
    all_precursor_offsets: &mut Vec<Option<Vec<PrecursorEntry>>>,
    model_id: &'a str,
    charge: i32,
    parent_mass: f32,
    seg_num: i32,
    kind: &'a str,
    values: &[f32],
) {
    record_kinds.push("table");
    model_ids.push(model_id);
    part_charges.push(Some(charge));
    part_mass_bits.push(Some(f32_to_i32_bits(parent_mass)));
    part_segs.push(Some(seg_num));
    ion_kinds.push(Some("-"));
    ion_charges.push(Some(0));
    ion_offset_bits.push(Some(0));
    table_kinds.push(Some(kind));
    all_values.push(Some(values.to_vec()));
    all_precursor_offsets.push(None);
}

fn encode_ion_type(ion: &IonType) -> (&'static str, i32, i32) {
    match ion {
        IonType::Prefix { charge, offset_bits } => ("prefix", *charge, *offset_bits as i32),
        IonType::Suffix { charge, offset_bits } => ("suffix", *charge, *offset_bits as i32),
        IonType::Noise => ("noise", 0, 0),
    }
}

// ── array builders ───────────────────────────────────────────────────────────

fn build_float_list(rows: Vec<Option<Vec<f32>>>) -> ArrayRef {
    let mut b = ListBuilder::new(Float32Builder::new());
    for row in rows {
        match row {
            Some(vals) => {
                for v in vals {
                    b.values().append_value(v);
                }
                b.append(true);
            }
            None => b.append(false),
        }
    }
    Arc::new(b.finish())
}

fn build_precursor_list(rows: Vec<Option<Vec<PrecursorEntry>>>) -> ArrayRef {
    let fields: Fields = vec![
        Field::new("reduced_charge", DataType::Int32, false),
        Field::new("offset", DataType::Float32, false),
        Field::new("tol_is_ppm", DataType::Boolean, false),
        Field::new("tol_val", DataType::Float32, false),
        Field::new("frequency", DataType::Float32, false),
    ]
    .into();

    let item_field = Arc::new(Field::new(
        "item",
        DataType::Struct(fields.clone()),
        false, // non-null items, matching schema
    ));
    let sb = StructBuilder::new(
        fields.clone(),
        vec![
            Box::new(Int32Builder::new()) as Box<dyn arrow::array::ArrayBuilder>,
            Box::new(Float32Builder::new()),
            Box::new(BooleanBuilder::new()),
            Box::new(Float32Builder::new()),
            Box::new(Float32Builder::new()),
        ],
    );
    let mut lb = ListBuilder::new(sb).with_field(item_field);

    for row in rows {
        match row {
            Some(entries) => {
                for (rc, off, is_ppm, tol, freq) in entries {
                    let sb = lb.values();
                    sb.field_builder::<Int32Builder>(0).unwrap().append_value(rc);
                    sb.field_builder::<Float32Builder>(1).unwrap().append_value(off);
                    sb.field_builder::<BooleanBuilder>(2).unwrap().append_value(is_ppm);
                    sb.field_builder::<Float32Builder>(3).unwrap().append_value(tol);
                    sb.field_builder::<Float32Builder>(4).unwrap().append_value(freq);
                    sb.append(true);
                }
                lb.append(true);
            }
            None => lb.append(false),
        }
    }
    Arc::new(lb.finish())
}

// ── null array helpers ───────────────────────────────────────────────────────

fn null_i32_array(n: usize) -> ArrayRef {
    let mut b = Int32Builder::new();
    for _ in 0..n { b.append_null(); }
    Arc::new(b.finish())
}

fn null_utf8_array(n: usize) -> ArrayRef {
    let mut b = StringBuilder::new();
    for _ in 0..n { b.append_null(); }
    Arc::new(b.finish())
}

fn null_bool_array(n: usize) -> ArrayRef {
    let mut b = BooleanBuilder::new();
    for _ in 0..n { b.append_null(); }
    Arc::new(b.finish())
}

fn null_f32_array(n: usize) -> ArrayRef {
    let mut b = Float32Builder::new();
    for _ in 0..n { b.append_null(); }
    Arc::new(b.finish())
}

fn null_float_list_array(n: usize) -> ArrayRef {
    let mut b = ListBuilder::new(Float32Builder::new());
    for _ in 0..n { b.append(false); }
    Arc::new(b.finish())
}

fn null_struct_list_array(n: usize, fields: Vec<Field>) -> ArrayRef {
    let fs: Fields = fields.into();
    let child_builders: Vec<Box<dyn arrow::array::ArrayBuilder>> = fs
        .iter()
        .map(|f| -> Box<dyn arrow::array::ArrayBuilder> {
            match f.data_type() {
                DataType::Int32 => Box::new(Int32Builder::new()),
                DataType::Float32 => Box::new(Float32Builder::new()),
                DataType::Boolean => Box::new(BooleanBuilder::new()),
                _ => Box::new(StringBuilder::new()),
            }
        })
        .collect();
    let item_field = Arc::new(Field::new(
        "item",
        DataType::Struct(fs.clone()),
        false, // non-null items, matching combined_schema's list_of_struct
    ));
    let sb = StructBuilder::new(fs.clone(), child_builders);
    let mut lb = ListBuilder::new(sb).with_field(item_field);
    for _ in 0..n { lb.append(false); }
    Arc::new(lb.finish())
}

fn null_i64_array(n: usize) -> ArrayRef {
    let mut b = Int64Builder::new();
    for _ in 0..n { b.append_null(); }
    Arc::new(b.finish())
}

fn null_int64_list_array(n: usize) -> ArrayRef {
    let mut b = ListBuilder::new(Int64Builder::new());
    for _ in 0..n { b.append(false); }
    Arc::new(b.finish())
}

fn null_int32_list_array(n: usize) -> ArrayRef {
    let mut b = ListBuilder::new(Int32Builder::new());
    for _ in 0..n { b.append(false); }
    Arc::new(b.finish())
}

// ── source batch (record_kind="source") ──────────────────────────────────────

/// Build the source-ledger batch: one row per (model_id, source).
fn build_source_batch(
    schema: &Arc<Schema>,
    model_id: &str,
    sources: &[(SourceLedger, CountStats)],
) -> Result<RecordBatch, TrainError> {
    let n = sources.len();

    // Shared columns.
    let record_kind: ArrayRef = Arc::new(StringArray::from(vec!["source"; n]));
    let model_id_col: ArrayRef = Arc::new(StringArray::from(vec![model_id; n]));

    // Source-ledger columns.
    let mut source_id_b = StringBuilder::new();
    let mut dataset_b = StringBuilder::new();
    let mut n_psms_b = Int64Builder::new();
    let mut date_b = StringBuilder::new();
    let mut weight_b = Float32Builder::new();
    let mut train_fdr_b = Float32Builder::new();
    let mut src_instrument_b = StringBuilder::new();
    let mut src_experiment_class_b = StringBuilder::new();

    for (ledger, _) in sources {
        source_id_b.append_value(&ledger.source_id);
        dataset_b.append_value(&ledger.dataset);
        n_psms_b.append_value(ledger.n_psms);
        date_b.append_value(&ledger.date);
        weight_b.append_value(ledger.weight);
        train_fdr_b.append_value(ledger.train_fdr);
        src_instrument_b.append_value(&ledger.instrument);
        src_experiment_class_b.append_value(&ledger.experiment_class);
    }

    // Columns that belong to other row kinds are null here.
    let null_i32 = null_i32_array(n);
    let null_utf8 = null_utf8_array(n);
    let null_bool = null_bool_array(n);
    let null_f32 = null_f32_array(n);
    let null_i32s = null_i32_array(n);
    let null_float_list = null_float_list_array(n);
    let null_int64_list = null_int64_list_array(n);
    let null_int32_list = null_int32_list_array(n);
    let null_charge_hist = null_struct_list_array(
        n,
        vec![
            Field::new("charge", DataType::Int32, false),
            Field::new("count", DataType::Int32, false),
        ],
    );
    let null_precursor_list = null_struct_list_array(
        n,
        vec![
            Field::new("reduced_charge", DataType::Int32, false),
            Field::new("offset", DataType::Float32, false),
            Field::new("tol_is_ppm", DataType::Boolean, false),
            Field::new("tol_val", DataType::Float32, false),
            Field::new("frequency", DataType::Float32, false),
        ],
    );

    let columns: Vec<ArrayRef> = vec![
        record_kind,
        model_id_col,
        // manifest-only → null
        null_utf8.clone(),   // activation
        null_utf8.clone(),   // instrument
        null_utf8.clone(),   // enzyme
        null_utf8.clone(),   // protocol
        null_i32.clone(),    // version
        null_f32.clone(),    // mme_val
        null_bool.clone(),   // mme_is_ppm
        null_bool.clone(),   // apply_deconv
        null_f32.clone(),    // deconv_tol
        null_i32.clone(),    // num_segments
        null_i32.clone(),    // max_rank
        null_i32.clone(),    // error_scaling_factor
        null_i32.clone(),    // min_charge
        null_i32.clone(),    // max_charge
        null_i32.clone(),    // num_precursor_off
        null_charge_hist,    // charge_hist
        // table-only → null
        null_i32s.clone(),   // part_charge
        null_i32s.clone(),   // part_mass_bits
        null_i32s.clone(),   // part_seg
        null_utf8.clone(),   // ion_kind
        null_i32s.clone(),   // ion_charge
        null_i32s.clone(),   // ion_offset_bits
        null_utf8.clone(),   // table_kind
        null_float_list,     // values
        null_precursor_list, // precursor_offsets
        // source-only → populated
        Arc::new(source_id_b.finish()),
        Arc::new(dataset_b.finish()),
        Arc::new(n_psms_b.finish()),
        Arc::new(date_b.finish()),
        Arc::new(weight_b.finish()),
        Arc::new(train_fdr_b.finish()),
        Arc::new(src_instrument_b.finish()),
        Arc::new(src_experiment_class_b.finish()),
        // stat-only → null
        null_int64_list,     // counts
        null_int32_list,     // charge_keys
    ];

    RecordBatch::try_new(schema.clone(), columns).map_err(|e| TrainError::Parquet(e.to_string()))
}

// ── stat batch (record_kind="stat") ──────────────────────────────────────────

/// Build the per-source sufficient-statistics batch.
///
/// Row encoding:
/// - `"rank"` — one row per `(source_id, Partition, IonType)`; counts = rank histogram Vec<u64>.
/// - `"error"` / `"noise_error"` — one row per `(source_id, Partition)`; counts = error histogram.
/// - `"existence"` — one row per `(source_id, Partition)`; counts[idx] = existence count for idx.
/// - `"charge"` — one row per source_id (partition zeroed); counts parallel to charge_keys.
fn build_stat_batch(
    schema: &Arc<Schema>,
    model_id: &str,
    sources: &[(SourceLedger, CountStats)],
) -> Result<RecordBatch, TrainError> {
    // Accumulate rows.
    let mut record_kinds: Vec<String> = Vec::new();
    let mut model_ids: Vec<String> = Vec::new();
    let mut source_ids: Vec<String> = Vec::new();
    let mut part_charges: Vec<Option<i32>> = Vec::new();
    let mut part_mass_bits_col: Vec<Option<i32>> = Vec::new();
    let mut part_segs: Vec<Option<i32>> = Vec::new();
    let mut ion_kinds: Vec<Option<String>> = Vec::new();
    let mut ion_charges: Vec<Option<i32>> = Vec::new();
    let mut ion_offset_bits: Vec<Option<i32>> = Vec::new();
    let mut table_kinds: Vec<Option<String>> = Vec::new();
    let mut all_counts: Vec<Option<Vec<i64>>> = Vec::new();
    let mut all_charge_keys: Vec<Option<Vec<i32>>> = Vec::new();

    for (ledger, stats) in sources {
        let sid = &ledger.source_id;

        // ── rank ──────────────────────────────────────────────────────────────
        // Sort for determinism: by (part, ion_kind, ion_charge, ion_offset_bits).
        let mut rank_keys: Vec<_> = stats.rank.keys().collect();
        rank_keys.sort_by_key(|(part, ion)| {
            let (ks, ic, iob) = encode_ion_type_str(ion);
            (*part, ks, ic, iob)
        });
        for (part, ion) in rank_keys {
            let counts_vec = &stats.rank[&(*part, *ion)];
            let (kind_str, ic, iob) = encode_ion_type_str(ion);
            record_kinds.push("stat".to_string());
            model_ids.push(model_id.to_string());
            source_ids.push(sid.clone());
            part_charges.push(Some(part.charge));
            part_mass_bits_col.push(Some(f32_to_i32_bits(part.parent_mass)));
            part_segs.push(Some(part.seg_num));
            ion_kinds.push(Some(kind_str.to_string()));
            ion_charges.push(Some(ic));
            ion_offset_bits.push(Some(iob));
            table_kinds.push(Some("rank".to_string()));
            all_counts.push(Some(counts_vec.iter().map(|&c| c as i64).collect()));
            all_charge_keys.push(None);
        }

        // ── error ─────────────────────────────────────────────────────────────
        let mut error_keys: Vec<_> = stats.error.keys().collect();
        error_keys.sort();
        for part in error_keys {
            let counts_vec = &stats.error[part];
            push_partition_stat_row(
                &mut record_kinds, &mut model_ids, &mut source_ids,
                &mut part_charges, &mut part_mass_bits_col, &mut part_segs,
                &mut ion_kinds, &mut ion_charges, &mut ion_offset_bits,
                &mut table_kinds, &mut all_counts, &mut all_charge_keys,
                model_id, sid, part, "error",
                counts_vec.iter().map(|&c| c as i64).collect(),
            );
        }

        // ── noise_error ───────────────────────────────────────────────────────
        let mut ne_keys: Vec<_> = stats.noise_error.keys().collect();
        ne_keys.sort();
        for part in ne_keys {
            let counts_vec = &stats.noise_error[part];
            push_partition_stat_row(
                &mut record_kinds, &mut model_ids, &mut source_ids,
                &mut part_charges, &mut part_mass_bits_col, &mut part_segs,
                &mut ion_kinds, &mut ion_charges, &mut ion_offset_bits,
                &mut table_kinds, &mut all_counts, &mut all_charge_keys,
                model_id, sid, part, "noise_error",
                counts_vec.iter().map(|&c| c as i64).collect(),
            );
        }

        // ── existence ─────────────────────────────────────────────────────────
        // Group by partition; within each partition, build counts[idx] = count for idx u32.
        // Collect all (partition, idx, count) tuples, group by partition.
        let mut existence_by_part: rustc_hash::FxHashMap<Partition, Vec<(u32, u64)>> =
            rustc_hash::FxHashMap::default();
        for (&(part, idx), &count) in &stats.existence {
            existence_by_part.entry(part).or_default().push((idx, count));
        }
        let mut ex_part_keys: Vec<Partition> = existence_by_part.keys().copied().collect();
        ex_part_keys.sort();
        for part in ex_part_keys {
            let entries = &existence_by_part[&part];
            // Build flat list: counts[idx] = count (fill gaps with 0).
            let max_idx = entries.iter().map(|(i, _)| *i).max().unwrap_or(0) as usize;
            let mut flat: Vec<i64> = vec![0i64; max_idx + 1];
            for &(idx, count) in entries {
                flat[idx as usize] = count as i64;
            }
            push_partition_stat_row(
                &mut record_kinds, &mut model_ids, &mut source_ids,
                &mut part_charges, &mut part_mass_bits_col, &mut part_segs,
                &mut ion_kinds, &mut ion_charges, &mut ion_offset_bits,
                &mut table_kinds, &mut all_counts, &mut all_charge_keys,
                model_id, sid, &part, "existence", flat,
            );
        }

        // ── charge ────────────────────────────────────────────────────────────
        // One row per source with part fields zeroed; charge map → two parallel lists.
        if !stats.charge.is_empty() {
            let mut charge_kv: Vec<(i32, u64)> = stats.charge.iter().map(|(&k, &v)| (k, v)).collect();
            charge_kv.sort_by_key(|(k, _)| *k);
            let keys: Vec<i32> = charge_kv.iter().map(|(k, _)| *k).collect();
            let counts: Vec<i64> = charge_kv.iter().map(|(_, v)| *v as i64).collect();

            record_kinds.push("stat".to_string());
            model_ids.push(model_id.to_string());
            source_ids.push(sid.clone());
            part_charges.push(Some(0));
            part_mass_bits_col.push(Some(0));
            part_segs.push(Some(0));
            ion_kinds.push(Some("-".to_string()));
            ion_charges.push(Some(0));
            ion_offset_bits.push(Some(0));
            table_kinds.push(Some("charge".to_string()));
            all_counts.push(Some(counts));
            all_charge_keys.push(Some(keys));
        }
    }

    let nrows = record_kinds.len();
    if nrows == 0 {
        return empty_stat_batch(schema);
    }

    // ── build Arrow arrays ───────────────────────────────────────────────────
    let record_kind_arr: ArrayRef = Arc::new(StringArray::from_iter_values(record_kinds.iter().map(|s| s.as_str())));
    let model_id_arr: ArrayRef = Arc::new(StringArray::from_iter_values(model_ids.iter().map(|s| s.as_str())));

    let source_id_arr: ArrayRef = Arc::new(StringArray::from(
        source_ids.iter().map(|s| s.as_str()).collect::<Vec<_>>(),
    ));

    let part_charge_arr: ArrayRef = Arc::new(arrow::array::Int32Array::from(part_charges));
    let part_mass_bits_arr: ArrayRef = Arc::new(arrow::array::Int32Array::from(part_mass_bits_col));
    let part_seg_arr: ArrayRef = Arc::new(arrow::array::Int32Array::from(part_segs));
    let ion_kind_arr: ArrayRef = Arc::new(StringArray::from(
        ion_kinds.into_iter().collect::<Vec<_>>(),
    ));
    let ion_charge_arr: ArrayRef = Arc::new(arrow::array::Int32Array::from(ion_charges));
    let ion_offset_bits_arr: ArrayRef = Arc::new(arrow::array::Int32Array::from(ion_offset_bits));
    let table_kind_arr: ArrayRef = Arc::new(StringArray::from(
        table_kinds.into_iter().collect::<Vec<_>>(),
    ));

    // counts: List<Int64>
    let counts_arr = build_int64_list(all_counts);
    // charge_keys: List<Int32>
    let charge_keys_arr = build_int32_list(all_charge_keys);

    // All other columns are null for stat rows.
    let null_i32 = null_i32_array(nrows);
    let null_i64 = null_i64_array(nrows);
    let null_utf8 = null_utf8_array(nrows);
    let null_bool = null_bool_array(nrows);
    let null_f32 = null_f32_array(nrows);
    let null_float_list = null_float_list_array(nrows);
    let null_charge_hist = null_struct_list_array(
        nrows,
        vec![
            Field::new("charge", DataType::Int32, false),
            Field::new("count", DataType::Int32, false),
        ],
    );
    let null_precursor_list = null_struct_list_array(
        nrows,
        vec![
            Field::new("reduced_charge", DataType::Int32, false),
            Field::new("offset", DataType::Float32, false),
            Field::new("tol_is_ppm", DataType::Boolean, false),
            Field::new("tol_val", DataType::Float32, false),
            Field::new("frequency", DataType::Float32, false),
        ],
    );

    let columns: Vec<ArrayRef> = vec![
        record_kind_arr,
        model_id_arr,
        // manifest-only → null
        null_utf8.clone(),   // activation
        null_utf8.clone(),   // instrument
        null_utf8.clone(),   // enzyme
        null_utf8.clone(),   // protocol
        null_i32.clone(),    // version
        null_f32.clone(),    // mme_val
        null_bool.clone(),   // mme_is_ppm
        null_bool.clone(),   // apply_deconv
        null_f32.clone(),    // deconv_tol
        null_i32.clone(),    // num_segments
        null_i32.clone(),    // max_rank
        null_i32.clone(),    // error_scaling_factor
        null_i32.clone(),    // min_charge
        null_i32.clone(),    // max_charge
        null_i32.clone(),    // num_precursor_off
        null_charge_hist,    // charge_hist
        // table-only → use the stat row partition/ion/table columns (populated above)
        part_charge_arr,
        part_mass_bits_arr,
        part_seg_arr,
        ion_kind_arr,
        ion_charge_arr,
        ion_offset_bits_arr,
        table_kind_arr,
        null_float_list,     // values (table payload, not used by stat)
        null_precursor_list, // precursor_offsets (table payload, not used by stat)
        // source/stat-only
        source_id_arr,
        null_utf8.clone(),   // dataset (stat rows don't repeat ledger fields)
        null_i64,            // n_psms
        null_utf8.clone(),   // date
        null_f32.clone(),    // weight
        null_f32.clone(),    // train_fdr
        null_utf8.clone(),   // src_instrument
        null_utf8.clone(),   // src_experiment_class
        counts_arr,
        charge_keys_arr,
    ];

    RecordBatch::try_new(schema.clone(), columns).map_err(|e| TrainError::Parquet(e.to_string()))
}

/// Build an empty stat batch (zero rows) — used when all sources have empty CountStats.
fn empty_stat_batch(schema: &Arc<Schema>) -> Result<RecordBatch, TrainError> {
    build_stat_batch_empty(schema)
}

fn build_stat_batch_empty(schema: &Arc<Schema>) -> Result<RecordBatch, TrainError> {
    let n = 0usize;
    let columns: Vec<ArrayRef> = schema
        .fields()
        .iter()
        .map(|f| -> ArrayRef {
            match f.data_type() {
                DataType::Utf8 => null_utf8_array(n),
                DataType::Int32 => null_i32_array(n),
                DataType::Int64 => null_i64_array(n),
                DataType::Float32 => null_f32_array(n),
                DataType::Boolean => null_bool_array(n),
                DataType::List(item) => match item.data_type() {
                    DataType::Float32 => null_float_list_array(n),
                    DataType::Int64 => null_int64_list_array(n),
                    DataType::Int32 => null_int32_list_array(n),
                    DataType::Struct(fields) => null_struct_list_array(n, fields.iter().map(|f| f.as_ref().clone()).collect()),
                    _ => null_float_list_array(n),
                },
                _ => null_utf8_array(n),
            }
        })
        .collect();
    RecordBatch::try_new(schema.clone(), columns).map_err(|e| TrainError::Parquet(e.to_string()))
}

#[allow(clippy::too_many_arguments)]
fn push_partition_stat_row(
    record_kinds: &mut Vec<String>,
    model_ids: &mut Vec<String>,
    source_ids: &mut Vec<String>,
    part_charges: &mut Vec<Option<i32>>,
    part_mass_bits_col: &mut Vec<Option<i32>>,
    part_segs: &mut Vec<Option<i32>>,
    ion_kinds: &mut Vec<Option<String>>,
    ion_charges: &mut Vec<Option<i32>>,
    ion_offset_bits: &mut Vec<Option<i32>>,
    table_kinds: &mut Vec<Option<String>>,
    all_counts: &mut Vec<Option<Vec<i64>>>,
    all_charge_keys: &mut Vec<Option<Vec<i32>>>,
    model_id: &str,
    source_id: &str,
    part: &Partition,
    table_kind: &str,
    counts: Vec<i64>,
) {
    record_kinds.push("stat".to_string());
    model_ids.push(model_id.to_string());
    source_ids.push(source_id.to_string());
    part_charges.push(Some(part.charge));
    part_mass_bits_col.push(Some(f32_to_i32_bits(part.parent_mass)));
    part_segs.push(Some(part.seg_num));
    ion_kinds.push(Some("-".to_string()));
    ion_charges.push(Some(0));
    ion_offset_bits.push(Some(0));
    table_kinds.push(Some(table_kind.to_string()));
    all_counts.push(Some(counts));
    all_charge_keys.push(None);
}

fn encode_ion_type_str(ion: &IonType) -> (&'static str, i32, i32) {
    match ion {
        IonType::Prefix { charge, offset_bits } => ("prefix", *charge, *offset_bits as i32),
        IonType::Suffix { charge, offset_bits } => ("suffix", *charge, *offset_bits as i32),
        IonType::Noise => ("noise", 0, 0),
    }
}

fn build_int64_list(rows: Vec<Option<Vec<i64>>>) -> ArrayRef {
    let mut b = ListBuilder::new(Int64Builder::new());
    for row in rows {
        match row {
            Some(vals) => {
                for v in vals {
                    b.values().append_value(v);
                }
                b.append(true);
            }
            None => b.append(false),
        }
    }
    Arc::new(b.finish())
}

fn build_int32_list(rows: Vec<Option<Vec<i32>>>) -> ArrayRef {
    let mut b = ListBuilder::new(Int32Builder::new());
    for row in rows {
        match row {
            Some(vals) => {
                for v in vals {
                    b.values().append_value(v);
                }
                b.append(true);
            }
            None => b.append(false),
        }
    }
    Arc::new(b.finish())
}
