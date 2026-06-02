//! Parquet writer: serialize one or more [`Param`] models to a single file.
//!
//! Call [`write_models`] with a slice of `(model_id, &Param)` pairs.
//! All models are written to **one** Parquet file using `combined_schema()`.
//! Manifest rows (scalar fields + charge_hist) come first, then table rows
//! (rank_dist, err_dist, frag_off, precursor_off, partition list).

use std::path::Path;
use std::sync::Arc;

use arrow::array::{
    ArrayRef, BooleanBuilder, Float32Builder, Int32Builder,
    ListBuilder, StringArray, StringBuilder,
    StructBuilder,
};
use arrow::datatypes::{DataType, Field, Fields, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use parquet::file::properties::WriterProperties;

use scoring_crate::param_model::{IonType, Param};

use crate::store::schema::combined_schema;
use crate::TrainError;

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

    // Table-only columns: null for manifest rows.
    let null_i32 = null_i32_array(n);
    let null_utf8 = null_utf8_array(n);
    let null_float_list = null_float_list_array(n);
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

    // manifest-only columns: null for table rows.
    let null_i32 = null_i32_array(nrows);
    let null_utf8 = null_utf8_array(nrows);
    let null_bool = null_bool_array(nrows);
    let null_f32 = null_f32_array(nrows);
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
